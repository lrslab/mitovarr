//! Preserve alternatives removed only by permissive primary-path cleanup.
use super::{graph_kmer, kmer, Candidate, Unitig};
use crate::io::fasta::write_fasta_records;
use anyhow::Result;
use serde_json::{json, Value};
use std::{
    collections::{BTreeSet, HashMap},
    path::Path,
};

pub(super) struct Alternative {
    pub sequence: Vec<u8>,
    pub mean_count: f64,
    pub min_count: u32,
}

pub(super) struct Collapsed {
    pub removed: Vec<u8>,
    retained: Vec<u8>,
    removed_mean: f64,
    retained_mean: f64,
}

impl Collapsed {
    pub fn new(weak: &Unitig, alternative: Alternative) -> Self {
        let reverse = kmer::reverse_complement(&weak.sequence);
        let flipped = reverse < weak.sequence;
        Self {
            removed: if flipped {
                reverse
            } else {
                weak.sequence.clone()
            },
            retained: if flipped {
                kmer::reverse_complement(&alternative.sequence)
            } else {
                alternative.sequence
            },
            removed_mean: weak.mean_count,
            retained_mean: alternative.mean_count,
        }
    }
}

pub(super) fn reports(bubbles: &[Collapsed]) -> Vec<Value> {
    bubbles.iter().enumerate().map(|(i,b)| json!({
        "bubble":i+1,"removed_length":b.removed.len(),"retained_length":b.retained.len(),
        "removed_mean_kmer_support":b.removed_mean,"retained_mean_kmer_support":b.retained_mean,
        "weak_strong_ratio":b.removed_mean/b.retained_mean,
        "reason":"outside_conservative_length_edit_or_per_kmer_support_rule",
        "status":"retained_for_primary_only; requires_independent_resolution_before_final",
        "ambiguity_scope":"original supported graph component; not a localized variant interval"
    })).collect()
}

pub(super) fn write(out: &Path, bubbles: &[Collapsed]) -> Result<()> {
    if !bubbles.is_empty() {
        let records: Vec<_> = bubbles
            .iter()
            .enumerate()
            .flat_map(|(i, b)| {
                [
                    (format!("bubble_{}_removed", i + 1), b.removed.clone()),
                    (format!("bubble_{}_retained", i + 1), b.retained.clone()),
                ]
            })
            .collect();
        write_fasta_records(&out.join("bubble_alternatives.fasta"), &records)?;
    }
    Ok(())
}

/// Membership is measured before cleanup so later overlapping collapses cannot
/// erase the flag. Disconnected small components do not veto the main circle.
/// Canonical (k-1)-mer nodes unite the two orientations of each physical graph.
#[derive(Default)]
pub(super) struct Components {
    nodes: graph_kmer::Map<usize>,
    bubbles: HashMap<usize, Vec<usize>>,
}

fn root(parents: &mut [usize], mut node: usize) -> usize {
    while parents[node] != node {
        let next = parents[node];
        parents[node] = parents[next];
        node = next;
    }
    node
}

impl Components {
    pub fn new(
        counts: &graph_kmer::Map<u32>,
        bubbles: &[Collapsed],
        k: usize,
        minimum: u32,
    ) -> Self {
        if bubbles.is_empty() {
            return Self::default();
        }
        let mut ids = HashMap::new();
        let mut parents = Vec::new();
        let mut sizes = Vec::new();
        for (&code, _) in counts.iter().filter(|(_, n)| **n >= minimum) {
            let mut pair = [0; 2];
            for (i, node) in [code.prefix(), code.suffix(k - 1)].into_iter().enumerate() {
                let node = graph_kmer::canonical(node, k - 1);
                pair[i] = *ids.entry(node).or_insert_with(|| {
                    let id = parents.len();
                    parents.push(id);
                    sizes.push(1usize);
                    id
                });
            }
            let mut a = root(&mut parents, pair[0]);
            let mut b = root(&mut parents, pair[1]);
            if a != b {
                if sizes[a] < sizes[b] {
                    std::mem::swap(&mut a, &mut b);
                }
                parents[b] = a;
                sizes[a] += sizes[b];
            }
        }
        let mut affected = HashMap::<usize, Vec<usize>>::new();
        for (i, b) in bubbles.iter().enumerate() {
            graph_kmer::visit(&b.removed[..k], None, k, 0, |_, code| {
                let node = graph_kmer::canonical(code.prefix(), k - 1);
                affected
                    .entry(root(&mut parents, ids[&node]))
                    .or_default()
                    .push(i + 1);
            });
        }
        let nodes = ids
            .into_iter()
            .filter_map(|(code, id)| {
                let component = root(&mut parents, id);
                affected
                    .contains_key(&component)
                    .then_some((code, component))
            })
            .collect();
        Self {
            nodes,
            bubbles: affected,
        }
    }

    pub fn affected(&self, candidate: &Candidate, k: usize) -> Vec<usize> {
        if self.nodes.is_empty() {
            return Vec::new();
        }
        let mut sequence = candidate.sequence.clone();
        if candidate.circular {
            sequence.extend(candidate.sequence.iter().copied().cycle().take(k - 1));
        }
        let mut components = BTreeSet::new();
        graph_kmer::visit(&sequence, None, k, 0, |_, code| {
            let node = graph_kmer::canonical(code.prefix(), k - 1);
            if let Some(component) = self.nodes.get(&node) {
                components.insert(*component);
            }
        });
        let mut bubbles: Vec<_> = components
            .into_iter()
            .flat_map(|id| self.bubbles[&id].iter().copied())
            .collect();
        bubbles.sort_unstable();
        bubbles.dedup();
        bubbles
    }
}
