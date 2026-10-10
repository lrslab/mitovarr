//! Fill unknown positions at the existing sample-backbone span. The span is a
//! prior, not a repeat-copy measurement; every emitted k-mer comes from reads.
use super::{code, Code, Graph};
use crate::kmer;
use serde_json::{json, Value};
use std::collections::HashMap;

const WIDTH: usize = 1024;

#[derive(Clone)]
struct Walk {
    node: Code,
    sequence: Vec<u8>,
    weight: f64,
    support: u32,
}

fn preference(a: &Walk, b: &Walk) -> std::cmp::Ordering {
    b.weight
        .total_cmp(&a.weight)
        .then(b.support.cmp(&a.support))
        .then(a.sequence.cmp(&b.sequence))
}

impl Graph {
    /// A length-layered maximum-support walk with fixed called positions.
    /// Same-node histories have identical future constraints, so only their
    /// highest-scoring spelling is kept. Repeat visits are bounded by the span,
    /// not a fixed copy ceiling. This produces one representative, not phasing.
    pub fn fill_masked_span(
        &self,
        backbone: &[u8],
        minimum: u32,
    ) -> (Option<(Vec<u8>, u32)>, Value) {
        let mut report = json!({"method":"maximum_support_masked_span","k":self.k,
            "minimum_pair_support":minimum,"target_length":backbone.len(),
            "target_source":"current sample/pileup backbone span; reference-guided length prior",
            "called_bases_preserved":true,"reference_bases_inserted":false,
            "beam_width":WIDTH,"repeat_limit":"sequence span only",
            "search_complete":false,"structure_verified":false});
        let overlap = self.k - 1;
        if !backbone.contains(&b'N') || backbone.len() < 2 * self.k {
            report["status"] = json!("not_needed_or_span_too_short");
            return (None, report);
        }
        let Some(start) = code(&backbone[..overlap]) else {
            report["status"] = json!("unknown_left_flank");
            return (None, report);
        };
        let mask = (1u128 << (2 * overlap)) - 1;
        let mut beam = vec![Walk {
            node: start,
            sequence: backbone[..overlap].to_vec(),
            weight: 0.0,
            support: u32::MAX,
        }];
        let mut states = 0usize;
        let mut merged = 0usize;
        let mut pruned = 0usize;
        for (position, &known) in backbone.iter().enumerate().skip(overlap) {
            let mut next = HashMap::<Code, Walk>::new();
            for previous in &beam {
                for (base, &letter) in b"ACGT".iter().enumerate() {
                    if known != b'N' && known != letter {
                        continue;
                    }
                    let edge = (previous.node << 2) | base as Code;
                    let support = self
                        .counts
                        .get(&kmer::canonical(edge, self.k))
                        .copied()
                        .unwrap_or(0);
                    if support < minimum {
                        continue;
                    }
                    states += 1;
                    let mut walk = previous.clone();
                    walk.node = edge & mask;
                    walk.sequence.push(letter);
                    walk.weight += f64::from(support).ln_1p();
                    walk.support = walk.support.min(support);
                    match next.entry(walk.node) {
                        std::collections::hash_map::Entry::Vacant(e) => {
                            e.insert(walk);
                        }
                        std::collections::hash_map::Entry::Occupied(mut e) => {
                            merged += 1;
                            if preference(&walk, e.get()).is_lt() {
                                e.insert(walk);
                            }
                        }
                    }
                }
            }
            if next.is_empty() {
                report["status"] = json!("no_walk_matching_called_bases_and_span");
                report["stopped_position_1based"] = json!(position + 1);
                report["states"] = json!(states);
                report["merged_histories"] = json!(merged);
                report["beam_pruned_states"] = json!(pruned);
                return (None, report);
            }
            beam = next.into_values().collect();
            beam.sort_by(preference);
            pruned += beam.len().saturating_sub(WIDTH);
            beam.truncate(WIDTH);
        }
        let best = beam.remove(0);
        report["status"] = json!("inferred_fill");
        report["states"] = json!(states);
        report["merged_histories"] = json!(merged);
        report["beam_pruned_states"] = json!(pruned);
        report["sum_log_pair_support"] = json!(best.weight);
        report["minimum_path_pair_support"] = json!(best.support);
        report["filled_bases"] = json!(backbone.iter().filter(|&&b| b == b'N').count());
        (Some((best.sequence, best.support)), report)
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
                *counts.entry(kmer::canonical(c, k)).or_default() += 2;
            });
        }
        Graph {
            k,
            counts,
            r1_occurrences: HashMap::new(),
            next: HashMap::new(),
            previous: HashMap::new(),
            pairs: reads.len(),
            max_read: 41,
            weak_edges_removed: 0,
        }
    }

    #[test]
    fn fills_sample_overlap_without_deleting_called_bases_or_inventing_a_gap() {
        let truth = dna(311, 913);
        let reads: Vec<_> = (0..=270)
            .step_by(30)
            .map(|p| truth[p..p + 41].to_vec())
            .collect();
        let mut masked = truth.clone();
        masked[50..261].fill(b'N');
        let (found, _) = graph(&reads, 11).fill_masked_span(&masked, 2);
        assert_eq!(found.unwrap().0, truth);
        assert!(graph(&reads, 13).fill_masked_span(&masked, 2).0.is_none());
        let mut missing = reads.clone();
        missing.remove(5);
        assert!(graph(&missing, 11).fill_masked_span(&masked, 2).0.is_none());
        let reverse: Vec<_> = reads
            .iter()
            .rev()
            .map(|r| kmer::reverse_complement(r))
            .collect();
        assert_eq!(
            graph(&reverse, 11)
                .fill_masked_span(&masked, 2)
                .0
                .unwrap()
                .0,
            truth
        );
    }

    #[test]
    fn masked_length_allows_more_than_four_repeat_visits_and_preserves_known_islands() {
        let left = dna(41, 812);
        let right = dna(41, 613);
        let unit = dna(37, 711);
        let truth = [left.clone(), unit.repeat(6), right.clone()].concat();
        let mut masked = truth.clone();
        masked[41..263].fill(b'N');
        masked[147..159].copy_from_slice(&truth[147..159]);
        let g = graph(std::slice::from_ref(&truth), 17);
        let (found, report) = g.fill_masked_span(&masked, 2);
        assert_eq!(found.unwrap().0, truth, "{report}");
        assert_eq!(report["structure_verified"], false);
        assert_eq!(report["called_bases_preserved"], true);
    }
}
