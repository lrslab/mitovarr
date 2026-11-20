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

const HOMOPOLYMER: &str = "homopolymer_minor_length_or_flank_variant_above_error_fraction";

pub(super) struct Collapsed {
    pub removed: Vec<u8>,
    retained: Vec<u8>,
    removed_mean: f64,
    retained_mean: f64,
    reason: &'static str,
}

impl Collapsed {
    pub fn new(weak: &Unitig, alternative: Alternative) -> Self {
        Self::oriented(
            weak.sequence.clone(),
            alternative.sequence,
            weak.mean_count,
            alternative.mean_count,
            "outside_conservative_length_edit_or_per_kmer_support_rule",
        )
    }

    /// A homopolymer tract's strongest minor length or flank variant against its
    /// selected spelling, each as one observed k-mer holding the whole run.
    pub fn homopolymer(
        removed: Vec<u8>,
        retained: Vec<u8>,
        removed_mean: f64,
        retained_mean: f64,
    ) -> Self {
        Self::oriented(removed, retained, removed_mean, retained_mean, HOMOPOLYMER)
    }

    /// Homopolymer tracts are marked by their own local primary interval.
    pub fn homopolymer_tract(&self) -> bool {
        self.reason == HOMOPOLYMER
    }

    /// Where `candidate` carries this bubble, as half-open intervals of the
    /// bases that differ from the other side. Both sides share their first and
    /// last k-1 bases (the bubble's junction nodes). First both whole flanks
    /// must match. Later cleanup of a nearby error can change a flank base
    /// further out, so when that finds nothing and both sides are at least
    /// 2(k-1) long, only the CORE_ANCHOR flank bases next to the differing
    /// bases must match (`Placement::CoreAnchors`). An empty result means the
    /// bubble was not placed on this sequence, not that it is absent from it.
    pub fn place(&self, candidate: &Candidate, k: usize) -> (Vec<(usize, usize)>, Placement) {
        let n = candidate.sequence.len();
        let flank = k - 1;
        if n < k || self.removed.len() < k || self.retained.len() < k {
            return (Vec::new(), Placement::NotPlaced);
        }
        let limit = self.retained.len().max(self.removed.len()) + k;
        let mut text = candidate.sequence.clone();
        if candidate.circular {
            text.extend(candidate.sequence.iter().copied().cycle().take(limit));
        }
        let site = Site {
            text: &text,
            length: n,
            circular: candidate.circular,
            limit,
        };
        let sides = [
            (self.retained.clone(), self.removed.clone()),
            (
                kmer::reverse_complement(&self.retained),
                kmer::reverse_complement(&self.removed),
            ),
        ];
        let mut intervals: Vec<_> = sides
            .iter()
            .flat_map(|(retained, removed)| site.find(retained, removed, flank))
            .collect();
        let mut placement = Placement::Flanks;
        let anchor = flank.min(CORE_ANCHOR);
        let trim = flank - anchor;
        if intervals.is_empty()
            && trim > 0
            && self.removed.len() >= 2 * flank
            && self.retained.len() >= 2 * flank
        {
            intervals = sides
                .iter()
                .flat_map(|(retained, removed)| {
                    site.find(
                        &retained[trim..retained.len() - trim],
                        &removed[trim..removed.len() - trim],
                        anchor,
                    )
                })
                .collect();
            placement = Placement::CoreAnchors;
        }
        if intervals.is_empty() {
            placement = Placement::NotPlaced;
        }
        intervals.sort_unstable();
        intervals.dedup();
        (intervals, placement)
    }

    fn oriented(
        removed: Vec<u8>,
        retained: Vec<u8>,
        removed_mean: f64,
        retained_mean: f64,
        reason: &'static str,
    ) -> Self {
        let reverse = kmer::reverse_complement(&removed);
        let flipped = reverse < removed;
        Self {
            removed: if flipped { reverse } else { removed },
            retained: if flipped {
                kmer::reverse_complement(&retained)
            } else {
                retained
            },
            removed_mean,
            retained_mean,
            reason,
        }
    }
}

/// Flank bases next to the differing bases that must still match when a base
/// further out in a k-1 junction flank differs from the sequence.
const CORE_ANCHOR: usize = 31;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Placement {
    /// Both whole junction flanks match.
    Flanks,
    /// Only the CORE_ANCHOR flank bases next to the differing bases match.
    CoreAnchors,
    NotPlaced,
}

/// A sequence searched for bubble sites. A circle's text runs `limit` bases
/// past its origin, so a site across the origin is found whole.
struct Site<'a> {
    text: &'a [u8],
    length: usize,
    circular: bool,
    limit: usize,
}

impl Site<'_> {
    /// Every start where the first `anchor` bases of `removed` match, closed by
    /// the exact `retained` spelling or else by the nearest copy of its last
    /// `anchor` bases within `limit`. Intervals cover the bases that differ.
    fn find(&self, retained: &[u8], removed: &[u8], anchor: usize) -> Vec<(usize, usize)> {
        let head = &removed[..anchor];
        let tail = &removed[removed.len() - anchor..];
        let n = self.length;
        let mut intervals = Vec::new();
        for start in 0..n {
            let rest = &self.text[start..];
            if !rest.starts_with(head) {
                continue;
            }
            let length = if rest.starts_with(retained) {
                retained.len()
            } else {
                let last = rest.len().min(self.limit).saturating_sub(anchor);
                match (1..=last).find(|&offset| rest[offset..].starts_with(tail)) {
                    Some(offset) => offset + anchor,
                    None => continue,
                }
            };
            let carried = &rest[..length];
            // The removed spelling itself can stay on a sequence through
            // another graph route; the same bases are then the open ones.
            let other = if carried == removed {
                retained
            } else {
                removed
            };
            let Some((from, to)) = differing_block(carried, other) else {
                continue;
            };
            // An empty allele marks one adjacent base.
            let from = from.min(length - 1);
            let (from, to) = (start + from, start + to.max(from + 1).min(length));
            if !self.circular || to <= n {
                intervals.push((from, to));
            } else if from >= n {
                intervals.push((from - n, to - n));
            } else {
                intervals.push((from, n));
                intervals.push((0, to - n));
            }
        }
        intervals
    }
}

/// The block of `a` where it differs from `b`, given shared flanks. Every
/// equivalent placement of an indel inside a run lies within it, so the block
/// can be empty only for an insertion into `b` with no repeated neighbour.
/// `None` when the two spellings are equal.
pub(super) fn differing_block(a: &[u8], b: &[u8]) -> Option<(usize, usize)> {
    if a == b {
        return None;
    }
    let prefix = a.iter().zip(b).take_while(|(x, y)| x == y).count();
    let suffix = a
        .iter()
        .rev()
        .zip(b.iter().rev())
        .take_while(|(x, y)| x == y)
        .count();
    let shortest = a.len().min(b.len());
    let start = prefix.min(shortest - suffix);
    let overlap = prefix.saturating_sub(shortest - suffix);
    Some((start, a.len() - suffix + overlap))
}

pub(super) fn reports(bubbles: &[Collapsed]) -> Vec<Value> {
    bubbles.iter().enumerate().map(|(i,b)| json!({
        "bubble":i+1,"removed_length":b.removed.len(),"retained_length":b.retained.len(),
        "removed_mean_kmer_support":b.removed_mean,"retained_mean_kmer_support":b.retained_mean,
        "weak_strong_ratio":b.removed_mean/b.retained_mean,
        "reason":b.reason,
        "status":"retained_for_primary_only; requires_independent_resolution_before_final",
        "ambiguity_scope":if b.homopolymer_tract() {
            "final publication: original supported graph component; primary: local homopolymer_length_family interval"
        } else {
            "final publication: original supported graph component; primary: local collapsed_bubble_alternative interval where its junction flanks are found on the sequence, otherwise listed as not placed"
        }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(n: usize, mut state: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 32) as usize & 3]
            })
            .collect()
    }

    fn bubble(removed: &[u8], retained: &[u8]) -> Collapsed {
        Collapsed::oriented(
            removed.to_vec(),
            retained.to_vec(),
            10.0,
            30.0,
            "outside_conservative_length_edit_or_per_kmer_support_rule",
        )
    }

    fn linear(sequence: Vec<u8>) -> Candidate {
        Candidate {
            sequence,
            circular: false,
        }
    }

    const K: usize = 15;

    /// A site between two 14-base junction flanks, inside a longer sequence.
    fn site(core: &[u8]) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let (left, right) = (dna(40, 77), dna(40, 91));
        let spelling = [&left[26..], core, &right[..14]].concat();
        let genome = [&left[..], core, &right[..]].concat();
        (spelling, genome, left)
    }

    #[test]
    fn differing_block_covers_every_placement_of_a_run_length_change() {
        // A substitution is its own base.
        assert_eq!(differing_block(b"ACGTACGT", b"ACGAACGT"), Some((3, 4)));
        // Ten against thirteen T: the whole run on either side.
        let short = [vec![b'C'; 20], vec![b'T'; 10], vec![b'G'; 20]].concat();
        let long = [vec![b'C'; 20], vec![b'T'; 13], vec![b'G'; 20]].concat();
        assert_eq!(differing_block(&short, &long), Some((20, 30)));
        assert_eq!(differing_block(&long, &short), Some((20, 33)));
        // An insertion with no repeated neighbour leaves an empty block on the
        // shorter side, at the insertion point.
        assert_eq!(differing_block(b"ACGTTGCA", b"ACGTCCTGCA"), Some((4, 4)));
        assert_eq!(differing_block(b"ACGTCCTGCA", b"ACGTTGCA"), Some((4, 6)));
        assert_eq!(differing_block(b"ACGT", b"ACGT"), None);
    }

    #[test]
    fn a_bubble_marks_only_its_own_site_on_either_strand() {
        let (retained, genome, _) = site(b"A");
        let (removed, _, _) = site(b"G");
        let b = bubble(&removed, &retained);
        assert_eq!(b.place(&linear(genome.clone()), K).0, vec![(40, 41)]);
        // The reverse strand carries the same site at the mirrored base.
        let reverse = kmer::reverse_complement(&genome);
        let n = genome.len();
        assert_eq!(b.place(&linear(reverse), K).0, vec![(n - 41, n - 40)]);
        // A sequence that lacks either junction flank is not marked at all.
        assert!(b.place(&linear(genome[..45].to_vec()), K).0.is_empty());
        assert!(b.place(&linear(dna(200, 5)), K).0.is_empty());
        // Sequences shorter than one k-mer cannot hold a bubble.
        assert!(b.place(&linear(genome[..10].to_vec()), K).0.is_empty());
    }

    #[test]
    fn a_site_across_a_circular_origin_is_split_not_widened() {
        let (retained, genome, _) = site(b"ACA");
        let (removed, _, _) = site(b"GTG");
        let b = bubble(&removed, &retained);
        let n = genome.len();
        // Rotate so that the three differing bases (40..43) straddle the origin.
        let mut rotated = genome.clone();
        rotated.rotate_left(41);
        let circle = Candidate {
            sequence: rotated,
            circular: true,
        };
        assert_eq!(b.place(&circle, K).0, vec![(0, 2), (n - 1, n)]);
        // The same sequence read as linear has no intact site.
        let cut = Candidate {
            sequence: circle.sequence.clone(),
            circular: false,
        };
        assert!(b.place(&cut, K).0.is_empty());
    }

    #[test]
    fn later_cleanup_of_the_retained_side_still_places_the_bubble() {
        // The sequence carries a third spelling between the two junction
        // flanks: neither saved side matches exactly.
        let (retained, _, _) = site(b"ACGGA");
        let (removed, _, _) = site(b"ACTTGA");
        let (_, genome, _) = site(b"ACGCA");
        let b = bubble(&removed, &retained);
        // Bases 42..45 of the carried spelling differ from the removed side.
        assert_eq!(b.place(&linear(genome), K).0, vec![(42, 44)]);
    }

    #[test]
    fn a_sequence_that_keeps_the_removed_spelling_marks_the_same_bases() {
        let (retained, _, _) = site(b"A");
        let (removed, genome, _) = site(b"G");
        let b = bubble(&removed, &retained);
        assert_eq!(b.place(&linear(genome), K).0, vec![(40, 41)]);
    }

    fn change(sequence: &mut [u8], at: usize) {
        sequence[at] = if sequence[at] == b'A' { b'C' } else { b'A' };
    }

    #[test]
    fn a_changed_outer_flank_base_falls_back_to_the_core_anchors() {
        // k = 63: 62-base junction flanks around a single-base site at 120.
        let (left, right) = (dna(120, 401), dna(120, 409));
        let side = |core: &[u8]| [&left[58..], core, &right[..62]].concat();
        let genome = [&left[..], b"A".as_slice(), &right[..]].concat();
        let n = genome.len();
        let b = bubble(&side(b"G"), &side(b"A"));
        assert_eq!(
            b.place(&linear(genome.clone()), 63),
            (vec![(120, 121)], Placement::Flanks)
        );
        // A changed base in the outer 31 bases of either flank, either strand.
        for at in [60, 85, 152, 180] {
            let mut changed = genome.clone();
            change(&mut changed, at);
            assert_eq!(
                b.place(&linear(changed.clone()), 63),
                (vec![(120, 121)], Placement::CoreAnchors),
                "{at}"
            );
            assert_eq!(
                b.place(&linear(kmer::reverse_complement(&changed)), 63),
                (vec![(n - 121, n - 120)], Placement::CoreAnchors),
                "{at}"
            );
        }
        // A changed base within 31 bases of the site is not placed.
        for at in [95, 119, 121, 151] {
            let mut changed = genome.clone();
            change(&mut changed, at);
            assert_eq!(
                b.place(&linear(changed), 63),
                (Vec::new(), Placement::NotPlaced),
                "{at}"
            );
        }
    }

    #[test]
    fn a_side_shorter_than_two_flanks_keeps_the_whole_flank_rule() {
        // One of two tandem copies of a 20-base unit deleted: the short side's
        // junction flanks overlap by one unit (104 < 2 x 62 bases).
        let (left, unit, right) = (dna(120, 501), dna(20, 503), dna(120, 509));
        let long = [&left[..], &unit[..], &unit[..], &right[..]].concat();
        let short = [&left[..], &unit[..], &right[..]].concat();
        let b = bubble(&short[78..182], &long[78..202]);
        assert_eq!(b.removed.len(), 104);
        // Either copy can hold the deletion, so both copies are open (plus any
        // flank base that happens to equal the unit's edge base).
        let (intervals, placement) = b.place(&linear(long.clone()), 63);
        assert_eq!(placement, Placement::Flanks);
        assert_eq!(intervals.len(), 1);
        let (from, to) = intervals[0];
        assert!(
            from <= 120 && to >= 160 && from + 3 >= 120 && to <= 163,
            "{intervals:?}"
        );
        let mut changed = long;
        change(&mut changed, 80);
        assert_eq!(
            b.place(&linear(changed), 63),
            (Vec::new(), Placement::NotPlaced)
        );
    }

    #[test]
    fn every_copy_of_a_repeated_site_is_marked() {
        let (retained, genome, _) = site(b"A");
        let (removed, _, _) = site(b"G");
        let b = bubble(&removed, &retained);
        let spacer = dna(30, 1234);
        let twice = [&genome[..], &spacer[..], &genome[..]].concat();
        let second = genome.len() + spacer.len() + 40;
        assert_eq!(
            b.place(&linear(twice), K).0,
            vec![(40, 41), (second, second + 1)]
        );
    }
}
