use crate::{cli::Reads, io::fastq, kmer};
use anyhow::{ensure, Result};
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};

type Code = kmer::Code;
pub(super) struct Graph {
    pub k: usize,
    pub counts: HashMap<Code, u32>,
    next: HashMap<Code, Vec<(Code, u8, u32)>>,
    previous: HashMap<Code, Vec<Code>>,
    pub pairs: usize,
    pub max_read: usize,
    pub weak_edges_removed: usize,
}

pub(super) fn code(sequence: &[u8]) -> Option<Code> {
    let mut result = None;
    if !sequence.is_empty() {
        kmer::visit(sequence, None, sequence.len(), 0, |_, c| result = Some(c));
    }
    result
}

impl Graph {
    pub fn build(reads: &Reads, k: usize, minq: u8, max_pairs: usize) -> Result<Self> {
        let mut counts = HashMap::<Code, u32>::new();
        let mut names = HashSet::new();
        let mut max_read = 0;
        fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
            if !names.insert(a.name(1)?.to_owned()) {
                return Ok(());
            }
            ensure!(
                names.len() <= max_pairs,
                "local assembly exceeds --max-records named FASTQ pairs"
            );
            let mut observed = HashSet::new();
            for read in [a, b] {
                max_read = max_read.max(read.sequence.len());
                kmer::visit(&read.sequence, Some(&read.quality), k, minq, |_, c| {
                    observed.insert(kmer::canonical(c, k));
                });
            }
            for c in observed {
                ensure!(
                    counts.len() < 2_000_000 || counts.contains_key(&c),
                    "local graph exceeds 2,000,000 distinct kmers"
                );
                *counts.entry(c).or_default() += 1;
            }
            Ok(())
        })?;
        let mut next = HashMap::<_, Vec<_>>::new();
        let mut previous = HashMap::<_, Vec<_>>::new();
        let mask = (1u128 << (2 * (k - 1))) - 1;
        for (&canonical, &count) in &counts {
            if count < 3 {
                continue;
            }
            for c in [canonical, kmer::reverse_code(canonical, k)] {
                let (from, to) = (c >> 2, c & mask);
                next.entry(from)
                    .or_default()
                    .push((to, b"ACGT"[(c & 3) as usize], count));
            }
        }
        let mut weak_edges_removed = 0;
        for edges in next.values_mut() {
            edges.sort_unstable();
            edges.dedup();
            // This is a dominant-consensus graph. A relative noise threshold
            // must still apply at high depth: ten observations out of thousands
            // must not evade it and exhaust paths through error combinations.
            // Retain competing outgoing edges at or above the 2% threshold.
            let maximum = edges.iter().map(|e| e.2).max().unwrap_or(0);
            let before = edges.len();
            edges.retain(|e| u64::from(e.2) * 50 >= u64::from(maximum));
            weak_edges_removed += before - edges.len();
        }
        for (&from, edges) in &next {
            for &(to, _, _) in edges {
                previous.entry(to).or_default().push(from);
            }
        }
        Ok(Self {
            k,
            counts,
            next,
            previous,
            pairs: names.len(),
            max_read,
            weak_edges_removed,
        })
    }

    pub fn support(&self, sequence: &[u8]) -> u32 {
        code(sequence)
            .and_then(|c| self.counts.get(&kmer::canonical(c, self.k)).copied())
            .unwrap_or(0)
    }

    pub fn paths(
        &self,
        left: &[u8],
        right: &[u8],
        max_span: usize,
        max_paths: usize,
        min_support: u32,
    ) -> Search {
        let start = code(&left[1..]).unwrap();
        let target = code(&right[1..]).unwrap();
        if start == target {
            return Search {
                limits: BTreeSet::from(["coincident_flank_nodes"]),
                ..Search::default()
            };
        }
        // Exclude excursions around the rest of the mitochondrial circle that
        // would revisit the unique left anchor before reaching the right one.
        let mut distance = HashMap::from([(target, 0usize)]);
        let mut queue = VecDeque::from([target]);
        while let Some(node) = queue.pop_front() {
            for &p in self.previous.get(&node).into_iter().flatten() {
                let edge = (p << 2) | (node & 3);
                if self
                    .counts
                    .get(&kmer::canonical(edge, self.k))
                    .copied()
                    .unwrap_or(0)
                    < min_support
                {
                    continue;
                }
                if p != start && !distance.contains_key(&p) {
                    distance.insert(p, distance[&node] + 1);
                    queue.push_back(p);
                }
            }
        }
        let mut search = Search::default();
        let mut low_support = HashSet::new();
        let mut sequence = left.to_vec();
        let mut visits = HashMap::from([(start, 1usize)]);
        let mut stack = vec![(start, 0usize, self.support(left))];
        while let Some(&(node, edge, minimum)) = stack.last() {
            if search.states >= 200_000 || search.paths.len() >= max_paths {
                search.limits.insert("path_or_state_limit");
                break;
            }
            let edges = self.next.get(&node).map(Vec::as_slice).unwrap_or(&[]);
            if edge == edges.len() || (node == target && stack.len() > 1) {
                if edge == 0 && node == target && sequence.ends_with(right) {
                    search.paths.push((sequence.clone(), minimum));
                }
                stack.pop();
                *visits.get_mut(&node).unwrap() -= 1;
                if !stack.is_empty() {
                    sequence.pop();
                }
                continue;
            }
            stack.last_mut().unwrap().1 += 1;
            let (to, base, count) = edges[edge];
            search.states += 1;
            if count < min_support {
                low_support.insert((node, to));
                continue;
            }
            if to == start || !distance.contains_key(&to) {
                continue;
            }
            if sequence.len() + 1 + distance[&to] > max_span {
                search.limits.insert("sequence_span_limit");
                continue;
            }
            if visits.get(&to).copied().unwrap_or(0) >= 4 {
                search.limits.insert("repeat_visit_limit");
                continue;
            }
            sequence.push(base);
            *visits.entry(to).or_default() += 1;
            stack.push((to, 0, minimum.min(count)));
        }
        search.paths.sort_by(|a, b| a.0.cmp(&b.0));
        search.paths.dedup_by(|a, b| a.0 == b.0);
        search.low_support_edges_excluded = low_support.len();
        search
    }
}

#[derive(Default)]
pub(super) struct Search {
    pub paths: Vec<(Vec<u8>, u32)>,
    pub limits: BTreeSet<&'static str>,
    pub states: usize,
    pub low_support_edges_excluded: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn depth_scaling_does_not_revive_rare_error_paths_or_remove_minor_paths() {
        use std::{fs, io::Write};
        let dir = std::env::temp_dir().join(format!("mapping-depth-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let reads = Reads {
            r1: dir.join("r1.fq"),
            r2: dir.join("r2.fq"),
        };
        let mut state = 891714u64;
        let sequence: Vec<_> = (0..90)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 32) as usize & 3]
            })
            .collect();
        let mut rare = sequence.clone();
        rare[30] = if rare[30] == b'A' { b'C' } else { b'A' };
        let mut minor = sequence.clone();
        minor[50] = if minor[50] == b'A' { b'C' } else { b'A' };
        let mut previous = None;
        for scale in [1, 2] {
            for mate in [1, 2] {
                let mut out =
                    fs::File::create(if mate == 1 { &reads.r1 } else { &reads.r2 }).unwrap();
                let mut id = 0;
                for (s, count) in [(&sequence, 1000), (&rare, 5), (&minor, 30)] {
                    for _ in 0..count * scale {
                        writeln!(
                            out,
                            "@pair{id}/{mate}\n{}\n+\n{}",
                            String::from_utf8_lossy(s),
                            "I".repeat(s.len())
                        )
                        .unwrap();
                        id += 1;
                    }
                }
            }
            let graph = Graph::build(&reads, 15, 20, 10000).unwrap();
            let paths = graph.paths(&sequence[..15], &sequence[75..], 200, 128, 3);
            assert!(paths.limits.is_empty());
            let actual: BTreeSet<_> = paths.paths.into_iter().map(|p| p.0).collect();
            assert!(actual.contains(&sequence) && actual.contains(&minor));
            assert!(!actual.contains(&rare));
            if let Some(old) = previous {
                assert_eq!(old, actual);
            }
            previous = Some(actual);
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn below_consensus_depth_cycle_does_not_exhaust_primary_path_search() {
        let ac = code(b"AC").unwrap();
        let cc = code(b"CC").unwrap();
        let cg = code(b"CG").unwrap();
        let gt = code(b"GT").unwrap();
        let counts = [
            (b"GAC", 20),
            (b"ACG", 20),
            (b"CGT", 20),
            (b"ACC", 5),
            (b"CCC", 5),
            (b"CCG", 5),
        ]
        .into_iter()
        .map(|(s, n)| (kmer::canonical(code(s).unwrap(), 3), n))
        .collect();
        let graph = Graph {
            k: 3,
            counts,
            next: HashMap::from([
                (ac, vec![(cc, b'C', 5), (cg, b'G', 20)]),
                (cc, vec![(cc, b'C', 5), (cg, b'G', 5)]),
                (cg, vec![(gt, b'T', 20)]),
            ]),
            previous: HashMap::from([(cc, vec![ac, cc]), (cg, vec![ac, cc]), (gt, vec![cg])]),
            pairs: 20,
            max_read: 100,
            weak_edges_removed: 0,
        };
        assert!(graph
            .paths(b"GAC", b"CGT", 100, 128, 3)
            .limits
            .contains("repeat_visit_limit"));
        let result = graph.paths(b"GAC", b"CGT", 100, 128, 10);
        assert!(result.limits.is_empty());
        assert_eq!(result.paths, vec![(b"GACGT".to_vec(), 20)]);
        assert_eq!(result.low_support_edges_excluded, 1);
        assert!(graph
            .paths(b"GAC", b"AAC", 100, 128, 10)
            .limits
            .contains("coincident_flank_nodes"));
    }
}
