//! Parameterized simple cycles, with occurrence-depth diagnostics kept separate
//! from fragment likelihoods. Graph topology alone never fixes the copy count.
use super::*;

#[derive(Clone)]
pub(super) struct Cycle {
    pub nodes: Vec<usize>,
    pub period: usize,
}

/// Iterative SCC decomposition, restricted to non-anchor segments. A simple
/// cycle has exactly one internal successor and predecessor at every node;
/// entry/exit edges may occur at any node. Complex SCCs use bounded cost-ordered traversal.
pub(super) fn cycles(graph: &PathGraph, anchors: &HashSet<usize>, overlap: usize) -> Vec<Cycle> {
    let n = graph.segments.len();
    let mut incoming = vec![Vec::new(); n];
    for (a, next) in graph.outgoing.iter().enumerate() {
        for &b in next {
            incoming[b].push(a);
        }
    }
    let excluded: HashSet<_> = anchors
        .iter()
        .copied()
        .chain(
            graph
                .segments
                .iter()
                .enumerate()
                .filter_map(|(id, s)| s.circular.then_some(id)),
        )
        .collect();
    let mut seen = excluded.clone();
    let mut order = Vec::new();
    for start in 0..n {
        if !seen.insert(start) {
            continue;
        }
        let mut stack = vec![(start, 0)];
        while let Some((node, offset)) = stack.last_mut() {
            if *offset == graph.outgoing[*node].len() {
                order.push(*node);
                stack.pop();
                continue;
            }
            let next = graph.outgoing[*node][*offset];
            *offset += 1;
            if seen.insert(next) {
                stack.push((next, 0));
            }
        }
    }
    seen = excluded;
    let mut result = Vec::new();
    for start in order.into_iter().rev() {
        if !seen.insert(start) {
            continue;
        }
        let mut component = HashSet::from([start]);
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            for &next in &incoming[node] {
                if seen.insert(next) {
                    component.insert(next);
                    stack.push(next);
                }
            }
        }
        if !component.iter().all(|&id| {
            graph.outgoing[id]
                .iter()
                .filter(|to| component.contains(to))
                .count()
                == 1
                && incoming[id]
                    .iter()
                    .filter(|from| component.contains(from))
                    .count()
                    == 1
        }) {
            continue;
        }
        let first = *component.iter().min().unwrap();
        let mut nodes = Vec::new();
        let mut node = first;
        loop {
            nodes.push(node);
            node = *graph.outgoing[node]
                .iter()
                .find(|id| component.contains(id))
                .unwrap();
            if node == first {
                break;
            }
        }
        let period = nodes
            .iter()
            .map(|&id| graph.segments[id].sequence.len().saturating_sub(overlap))
            .sum();
        if period > 0 {
            result.push(Cycle { nodes, period });
        }
    }
    result.sort_by_key(|c| c.nodes[0]);
    result
}

pub(super) struct Depth {
    counts: HashMap<kmer::Code, u64>,
}
impl Depth {
    pub fn new(loci: &[Locus], graph: &PathGraph) -> Self {
        let mut counts = HashMap::new();
        for l in loci.iter().filter(|l| !l.repeat_cycles.is_empty()) {
            for id in std::iter::once(l.from)
                .chain(l.routes.iter().map(|r| r.to))
                .chain(l.repeat_cycles.iter().flat_map(|c| c.nodes.iter().copied()))
            {
                kmer::visit(&graph.segments[id].sequence, None, 31, 0, |_, c| {
                    counts.insert(kmer::canonical(c, 31), 0);
                });
            }
        }
        Self { counts }
    }
    pub fn observe_r1(&mut self, read: &fastq::Record, minq: u8) {
        if self.counts.is_empty() {
            return;
        }
        // Count each observed window, including repeated offsets in R1. R2 is
        // absent from this statistic, so overlapping mates cannot inflate it.
        kmer::visit(&read.sequence, Some(&read.quality), 31, minq, |_, c| {
            if let Some(n) = self.counts.get_mut(&kmer::canonical(c, 31)) {
                *n = n.saturating_add(1);
            }
        });
    }
    pub fn report(
        &self,
        locus: &Locus,
        graph: &PathGraph,
        occurrence: &HashMap<kmer::Code, usize>,
    ) -> Value {
        if locus.repeat_cycles.is_empty() {
            return Value::Null;
        }
        let mut flanks = BTreeSet::new();
        for id in std::iter::once(locus.from).chain(locus.routes.iter().map(|r| r.to)) {
            let s = &graph.segments[id].sequence;
            let window = if id == locus.from {
                &s[s.len().saturating_sub(500)..]
            } else {
                &s[..s.len().min(500)]
            };
            kmer::visit(window, None, 31, 0, |_, c| {
                let c = kmer::canonical(c, 31);
                if occurrence.get(&c) == Some(&2) {
                    flanks.insert(c);
                }
            });
        }
        let mut depths: Vec<_> = flanks
            .iter()
            .map(|c| self.counts.get(c).copied().unwrap_or(0))
            .collect();
        depths.sort_unstable();
        let median = depths.get(depths.len() / 2).copied().unwrap_or(0) as f64;
        let models: Vec<_>=locus.repeat_cycles.iter().map(|cycle| {
            let mut repeat_codes=BTreeSet::new();
            for &id in &cycle.nodes {
                kmer::visit(&graph.segments[id].sequence,None,31,0,|_,c| {repeat_codes.insert(kmer::canonical(c,31));});
            }
            // A model predicts the number of genomic copies of every motif
            // 31-mer in this complete local route, including partial turns and
            // the graph-overlap bases. Do not round a raw graph dp ratio.
            let observed: u64=repeat_codes.iter().map(|c|self.counts.get(c).copied().unwrap_or(0)).sum();
            let candidates: Vec<_>=locus.routes.iter().enumerate().map(|(i,r)| {
                let mut instances=0;
                kmer::visit(&r.sequence,None,31,0,|_,c| {instances+=usize::from(repeat_codes.contains(&kmer::canonical(c,31)));});
                let turns=cycle.nodes.iter().map(|id|r.nodes.iter().filter(|n|*n==id).count()).min().unwrap_or(0);
                let expected=median*instances as f64;
                json!({"path":i+1,"full_cycle_turns":turns,"connecting_bases":r.right_start-r.left_end,
                    "predicted_motif_31mer_instances":instances,"observed_r1_occurrences":observed,
                    "expected_r1_occurrences_at_flank_depth":expected,
                    "observed_expected_ratio":if expected>0.0 {Some(observed as f64/expected)}else{None},
                    "fragment_score":locus.totals[i],"spanning_fragments":locus.spanning[i]})
            }).collect();
            json!({"segments":cycle.nodes.iter().map(|n|n+1).collect::<Vec<_>>(),"period_bases":cycle.period,
                "motif_31mers":repeat_codes.len(),"candidates":candidates,
                "higher_copy_tail":"unsearched_and_not_excluded","copy_number_identified":false})
        }).collect();
        json!({"simple_cycles":models,"flank_unique_31mers":flanks.len(),"flank_r1_occurrence_depth_median":median,
            "coverage_statistic":"high-quality R1 31-mer window occurrences; repeated offsets retained, R2 excluded to avoid overlap double counting",
            "coverage_use":"diagnostic prediction for each copy candidate; not an independent likelihood, hard copy bound or rounded dp ratio",
            "coverage_assumptions":"approximately uniform local sampling and no other genomic source of motif k-mers; local coverage bias and NUMTs may violate this"})
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::denovo::threading::Segment;
    #[test]
    fn two_segment_cycle_is_parameterized_but_complex_scc_is_not() {
        let mut graph = PathGraph {
            segments: (0..4)
                .map(|_| Segment {
                    sequence: vec![b'A'; 83],
                    circular: false,
                })
                .collect(),
            outgoing: vec![vec![1], vec![2, 3], vec![1], vec![]],
        };
        let cycles = cycles(&graph, &HashSet::from([0, 3]), 62);
        assert_eq!(cycles.len(), 1);
        assert_eq!(cycles[0].nodes, vec![1, 2]);
        assert_eq!(cycles[0].period, 42);
        graph.outgoing[2].push(2);
        assert!(super::cycles(&graph, &HashSet::from([0, 3]), 62).is_empty());
        let closed = PathGraph {
            segments: vec![Segment {
                sequence: vec![b'A'; 10],
                circular: true,
            }],
            outgoing: vec![vec![0]],
        };
        assert!(super::cycles(&closed, &HashSet::new(), 62).is_empty());
    }
    #[test]
    fn occurrence_depth_retains_repeated_windows_in_one_mate() {
        let mut depth = Depth {
            counts: HashMap::from([(0, 0)]),
        };
        let read = fastq::Record {
            header: "@test/1".into(),
            sequence: vec![b'A'; 100],
            quality: vec![b'I'; 100],
            plus: "+".into(),
        };
        depth.observe_r1(&read, 20);
        assert_eq!(depth.counts[&0], 70);
    }
}
