//! Reference-span-guided proposal search, restricted to the mapping workflow.
//! The target is a search heuristic, never evidence for a copy number.
use super::{code, Code, Graph};
use crate::kmer;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};

const BEAM: usize = 256;
const STATES: usize = 50_000;
const REACHABLE: usize = 100_000;
const RETAINED: usize = 8;

#[derive(Clone)]
struct State {
    sequence: Vec<u8>,
    node: Code,
    visits: BTreeMap<Code, u8>,
    weight: u64,
    support: u32,
}

#[derive(Clone)]
struct Step {
    edge: Code,
    end: Code,
    bases: Vec<u8>,
    weight: u64,
    support: u32,
    visits: u8,
}

impl Graph {
    /// Keep the original complete flanks even when the graph k is very small.
    /// Raw counted edges allow the legacy threshold relaxation down to two
    /// named pairs; the strict graph and its relative-error cleanup are intact.
    pub fn legacy_paths(
        &self,
        left: &[u8],
        right: &[u8],
        target_len: usize,
        max_span: usize,
        minimum: u32,
    ) -> (Vec<(Vec<u8>, u32)>, Value) {
        let mut report = json!({"k":self.k,"minimum_pair_support":minimum,
            "target_length":target_len,"search_complete":false,
            "method":"legacy_reference_span_beam","max_states":STATES,
            "beam_width":BEAM,"max_reachable_nodes":REACHABLE,
            "retained_paths":RETAINED,"reference_bases_inserted":false});
        let count = |c| {
            self.counts
                .get(&kmer::canonical(c, self.k))
                .copied()
                .unwrap_or(0)
        };
        let flank_support = |s: &[u8]| {
            s.windows(self.k)
                .map(|s| code(s).map_or(0, count))
                .min()
                .unwrap_or(0)
        };
        if left.len() < self.k
            || right.len() < self.k
            || flank_support(left) < minimum
            || flank_support(right) < minimum
        {
            report["status"] = json!("unsupported_complete_flank");
            return (vec![], report);
        }
        let overlap = self.k - 1;
        let mask = (1u128 << (2 * overlap)) - 1;
        let start = code(&left[left.len() - overlap..]).unwrap();
        let target = code(&right[right.len() - overlap..]).unwrap();
        let ceiling =
            max_span.min(target_len.saturating_add(target_len.clamp(250, 2000) / 8 + 2 * overlap));
        report["maximum_sequence_length"] = json!(ceiling);
        let mut distances = HashMap::from([(target, 0usize)]);
        let mut queue = VecDeque::from([target]);
        let mut reachability_complete = true;
        while let Some(node) = queue.pop_front() {
            if distances.len() >= REACHABLE {
                reachability_complete = false;
                break;
            }
            for base in 0..4u128 {
                let previous = (base << (2 * (self.k - 2))) | (node >> 2);
                if count((previous << 2) | (node & 3)) >= minimum
                    && !distances.contains_key(&previous)
                {
                    distances.insert(previous, distances[&node] + 1);
                    queue.push_back(previous);
                }
            }
        }
        report["reachability_complete"] = json!(reachability_complete);
        report["reverse_reachable_nodes"] = json!(distances.len());
        if reachability_complete && !distances.contains_key(&start) {
            report["status"] = json!("no_graph_connection");
            return (vec![], report);
        }
        let outgoing = |node: Code| -> Vec<(Code, u32)> {
            (0..4u128)
                .filter_map(|base| {
                    let edge = (node << 2) | base;
                    let n = count(edge);
                    (n >= minimum
                        && (!reachability_complete || distances.contains_key(&(edge & mask))))
                    .then_some((edge, n))
                })
                .collect()
        };
        let incoming = |node: Code| {
            (0..4u128)
                .filter(|&base| count((base << (2 * overlap)) | node) >= minimum)
                .count()
        };
        let mut flank_counts: Vec<_> = left
            .windows(self.k)
            .chain(right.windows(self.k))
            .map(|s| count(code(s).unwrap()))
            .collect();
        flank_counts.sort_unstable();
        let baseline = flank_counts[flank_counts.len() / 2].max(1) as f64;
        let mut cache = HashMap::<Code, Vec<Step>>::new();
        let mut beam = vec![State {
            sequence: left.to_vec(),
            node: start,
            visits: BTreeMap::new(),
            weight: 0,
            support: flank_support(left),
        }];
        let mut found = BTreeMap::<Vec<u8>, (u32, u64)>::new();
        let mut limits = BTreeSet::new();
        let mut states = 0;
        // Old beam ordering: prefer progress toward the target, then support.
        // Complete paths rank by target distance, then support, with stable ties.
        while !beam.is_empty() {
            let mut next = Vec::new();
            for state in beam {
                if states >= STATES {
                    limits.insert("state_limit");
                    break;
                }
                states += 1;
                if state.sequence.len() >= left.len() + right.len()
                    && state.sequence.ends_with(right)
                {
                    found
                        .entry(state.sequence)
                        .or_insert((state.support, state.weight));
                    continue;
                }
                let steps = cache.entry(state.node).or_insert_with(|| {
                    outgoing(state.node)
                        .into_iter()
                        .map(|(edge, support)| {
                            let mut step = Step {
                                edge,
                                end: edge & mask,
                                bases: vec![b"ACGT"[(edge & 3) as usize]],
                                weight: u64::from(support),
                                support,
                                visits: 1,
                            };
                            let mut seen = HashSet::from([state.node]);
                            while step.end != target
                                && step.end != start
                                && seen.insert(step.end)
                                && incoming(step.end) == 1
                                && step.bases.len() < ceiling
                            {
                                let exits = outgoing(step.end);
                                if exits.len() != 1 {
                                    break;
                                }
                                let (edge, n) = exits[0];
                                step.end = edge & mask;
                                step.bases.push(b"ACGT"[(edge & 3) as usize]);
                                step.weight += u64::from(n);
                                step.support = step.support.min(n);
                            }
                            let ratio = step.weight as f64 / step.bases.len() as f64 / baseline;
                            let copies = if ratio >= 1.75
                                && (step.end == state.node
                                    || incoming(state.node) > 1
                                    || outgoing(step.end).len() > 1)
                            {
                                ((ratio + 0.25).floor() as u8).clamp(1, 4)
                            } else {
                                1
                            };
                            step.visits = if step.bases.len() + overlap <= (self.k * 12).max(96) {
                                (copies + 2).min(8)
                            } else {
                                (copies + 1).min(6)
                            };
                            step
                        })
                        .collect()
                });
                for step in steps.iter() {
                    if state.sequence.len() + step.bases.len() > ceiling {
                        limits.insert("sequence_span_limit");
                        continue;
                    }
                    if state.visits.get(&step.edge).copied().unwrap_or(0) >= step.visits {
                        limits.insert("repeat_visit_limit");
                        continue;
                    }
                    let mut child = state.clone();
                    child.sequence.extend_from_slice(&step.bases);
                    child.node = step.end;
                    child.weight += step.weight;
                    child.support = child.support.min(step.support);
                    *child.visits.entry(step.edge).or_default() += 1;
                    if child
                        .sequence
                        .windows(left.len())
                        .skip(1)
                        .any(|s| s == left)
                    {
                        limits.insert("complete_left_flank_revisit");
                        continue;
                    }
                    next.push(child);
                }
            }
            if states >= STATES {
                break;
            }
            next.sort_by(|a, b| {
                b.sequence
                    .len()
                    .min(target_len)
                    .cmp(&a.sequence.len().min(target_len))
                    .then(
                        a.sequence
                            .len()
                            .abs_diff(target_len)
                            .cmp(&b.sequence.len().abs_diff(target_len)),
                    )
                    .then(b.weight.cmp(&a.weight))
                    .then(a.sequence.cmp(&b.sequence))
            });
            next.dedup_by(|a, b| a.sequence == b.sequence);
            if next.len() > BEAM {
                limits.insert("beam_width_limit");
                next.truncate(BEAM);
            }
            beam = next;
        }
        let found_count = found.len();
        let mut found: Vec<_> = found.into_iter().collect();
        found.sort_by(|(a, (_, aw)), (b, (_, bw))| {
            a.len()
                .abs_diff(target_len)
                .cmp(&b.len().abs_diff(target_len))
                .then(bw.cmp(aw))
                .then(a.cmp(b))
        });
        if found.len() > RETAINED {
            limits.insert("retained_path_limit");
            found.truncate(RETAINED);
        }
        report["states"] = json!(states);
        report["limits"] = json!(limits);
        report["discovered_paths"] = json!(found_count);
        report["candidate_count"] = json!(found.len());
        report["status"] = json!(if found.is_empty() {
            "no_path_within_search_bounds"
        } else {
            "inferred_candidates"
        });
        (
            found.into_iter().map(|(s, (n, _))| (s, n)).collect(),
            report,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(n: usize, mut seed: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed >> 32) as usize & 3]
            })
            .collect()
    }
    fn graph(reads: &[Vec<u8>], k: usize) -> Graph {
        let mut counts = HashMap::new();
        for read in reads {
            kmer::visit(read, None, k, 0, |_, c| {
                *counts.entry(kmer::canonical(c, k)).or_default() += 20
            });
        }
        Graph {
            k,
            counts,
            r1_occurrences: HashMap::new(),
            next: HashMap::new(),
            previous: HashMap::new(),
            pairs: reads.len() * 20,
            max_read: reads.iter().map(Vec::len).max().unwrap_or(0),
            weak_edges_removed: 0,
        }
    }

    #[test]
    fn k11_connects_eleven_base_read_overlaps_but_never_invents_a_missing_edge() {
        let sequence = dna(311, 913);
        let reads: Vec<_> = (0..=270)
            .step_by(30)
            .map(|p| sequence[p..p + 41].to_vec())
            .collect();
        // These fragments overlap by 11 bases: k=11 works, k=13 does not.
        let left = &sequence[..31];
        let right = &sequence[280..];
        let (paths, report) = graph(&reads, 11).legacy_paths(left, right, 311, 500, 3);
        assert_eq!(paths[0].0, sequence, "{report}");
        assert!(!report["search_complete"].as_bool().unwrap());
        assert!(graph(&reads, 13)
            .legacy_paths(left, right, 311, 500, 3)
            .0
            .is_empty());
        let mut reversed = reads.clone();
        reversed.reverse();
        assert_eq!(
            paths,
            graph(&reversed, 11)
                .legacy_paths(left, right, 311, 500, 3)
                .0
        );
        let mut low_support = graph(&reads, 11);
        for count in low_support.counts.values_mut() {
            *count = 2;
        }
        assert!(low_support
            .legacy_paths(left, right, 311, 500, 3)
            .0
            .is_empty());
        assert_eq!(
            low_support.legacy_paths(left, right, 311, 500, 2).0[0].0,
            sequence
        );
        let mut missing = reads;
        missing.remove(5);
        assert!(graph(&missing, 11)
            .legacy_paths(left, right, 311, 500, 3)
            .0
            .is_empty());
    }

    #[test]
    fn length_prior_proposes_repeat_copies_without_claiming_they_are_true() {
        let left = dna(41, 812);
        let right = dna(41, 613);
        let unit = dna(37, 711);
        let truth = [left.clone(), unit.repeat(3), right.clone()].concat();
        let graph = graph(std::slice::from_ref(&truth), 11);
        let (paths, report) = graph.legacy_paths(&left, &right, truth.len(), 500, 3);
        assert_eq!(paths[0].0, truth, "{report}");
        assert!(paths
            .iter()
            .all(|p| p.0.starts_with(&left) && p.0.ends_with(&right)));
        let shorter = [left.clone(), unit.repeat(2), right.clone()].concat();
        let (biased, _) = graph.legacy_paths(&left, &right, shorter.len(), 500, 3);
        assert_eq!(
            biased[0].0, shorter,
            "a different length prior can prefer a wrong copy number"
        );
        assert_eq!(report["search_complete"], false);
    }

    #[test]
    fn internal_short_target_match_does_not_finish_before_the_complete_right_flank() {
        let left = dna(41, 235);
        let right = dna(41, 756);
        let sequence = [
            left.clone(),
            dna(60, 712),
            right[31..].to_vec(),
            dna(60, 771),
            right.clone(),
        ]
        .concat();
        let (paths, report) = graph(std::slice::from_ref(&sequence), 11).legacy_paths(
            &left,
            &right,
            sequence.len(),
            500,
            3,
        );
        assert_eq!(paths[0].0, sequence, "{report}");
        assert!(paths.iter().all(|p| p.0.ends_with(&right)));
    }
}
