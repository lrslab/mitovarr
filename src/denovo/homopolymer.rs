//! Homopolymer tracts. When k exceeds a run's length, each run length has its
//! own spanning k-mers, so stutter and length heteroplasmy at one long run
//! interleave into bubbles that the pairwise bubble rule cannot simplify. Reads
//! crossing a long run also often miscall the next few bases, which adds flank
//! variants on the side each read reaches after the run. Each tract collapses to
//! its best-supported run length and flank spelling.
use super::{bubbles, graph_kmer, paths::primary::Region, Candidate, Edge};
use crate::kmer;
use serde_json::{json, Value};
use std::{
    cmp::Ordering,
    collections::{BTreeMap, HashMap},
};

/// A selected run at least this long defines a family.
const MIN_RUN: usize = 10;
/// Shorter runs between the same flanks still count as its length variants.
const MIN_VARIANT_RUN: usize = 6;
/// Flank bases identifying one tract on each side of the run.
const FLANK: usize = 12;
/// A flank variant keeps one flank exact and differs at most this much in the other.
const FLANK_MISMATCHES: usize = 3;
/// A flank variant's run length is within this many bases of the selected length.
const MAX_LENGTH_SHIFT: usize = 6;

/// Run base and flanks, spelled on the strand where the run is A or C.
type Key = (u8, [u8; FLANK], [u8; FLANK]);
type Side = (u8, [u8; FLANK]);

/// Maximal A or C runs of at least `MIN_VARIANT_RUN` bases: (base, start, length).
/// The reverse strand reports each T or G run once as an A or C run.
fn runs(sequence: &[u8]) -> Vec<(u8, usize, usize)> {
    let mut runs = Vec::new();
    let mut start = 0;
    while start < sequence.len() {
        let base = sequence[start];
        let length = sequence[start..].iter().take_while(|&&b| b == base).count();
        if matches!(base, b'A' | b'C') && length >= MIN_VARIANT_RUN {
            runs.push((base, start, length));
        }
        start += length;
    }
    runs
}

fn flank(sequence: &[u8]) -> [u8; FLANK] {
    sequence.try_into().expect("flank length")
}

fn mismatches(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

#[derive(Default)]
struct Observed {
    support: u64,
    kmers: usize,
    example: (u32, Vec<u8>),
    before: Vec<u8>,
    after: Vec<u8>,
}

struct Length {
    length: usize,
    mean: f64,
    kmers: usize,
    example: Vec<u8>,
}

/// One exact flank pair and the run lengths seen between it.
pub(super) struct Family {
    key: Key,
    lengths: Vec<Length>,
    selected: usize,
    strength: f64,
    /// Longest observed sequence before and after the selected run.
    before: Vec<u8>,
    after: Vec<u8>,
    /// The consensus family of this family's tract, unless it is that family.
    variant_of: Option<usize>,
    supported: bool,
}

impl Family {
    fn selected(&self) -> &Length {
        self.lengths
            .iter()
            .find(|l| l.length == self.selected)
            .expect("selected length")
    }

    fn strongest_minor(&self) -> Option<&Length> {
        self.lengths
            .iter()
            .filter(|l| l.length != self.selected)
            .max_by(|a, b| a.mean.total_cmp(&b.mean).then(b.length.cmp(&a.length)))
    }

    /// The tract's own spelling: exact flanks and its selected run length.
    fn accepts(&self, before: &[u8], after: &[u8], length: usize) -> bool {
        let (left, right, _) = self.compare(before, after);
        // A k-mer ending inside the run holds at least `length` bases of it.
        let length_kept = if before.is_empty() || after.is_empty() {
            length <= self.selected
        } else {
            length == self.selected
        };
        left + right == 0 && length_kept
    }

    /// Bases on each side of a run, as far as a k-mer holds them, against this
    /// family's flanks: (left mismatches, right mismatches, either side complete).
    fn compare(&self, before: &[u8], after: &[u8]) -> (usize, usize, bool) {
        let (l, r) = (before.len().min(FLANK), after.len().min(FLANK));
        let left = mismatches(&before[before.len() - l..], &self.key.1[FLANK - l..]);
        let right = mismatches(&after[..r], &self.key.2[..r]);
        let exact_side = (l == FLANK && left == 0) || (r == FLANK && right == 0);
        (left, right, exact_side)
    }

    /// The same tract seen in a k-mer: one complete flank exact and few
    /// mismatches in the other, within the flank bases held. A k-mer ending
    /// inside the run holds only a lower bound of its length.
    fn near(&self, before: &[u8], after: &[u8], length: usize) -> bool {
        let (left, right, exact_side) = self.compare(before, after);
        let shift = if before.is_empty() || after.is_empty() {
            length.saturating_sub(self.selected)
        } else {
            length.abs_diff(self.selected)
        };
        exact_side && left + right <= FLANK_MISMATCHES && shift <= MAX_LENGTH_SHIFT
    }

    /// Few mismatches over both longest observed contexts. Another haplotype
    /// usually also differs beyond the flank; strong read-through miscalls do not.
    fn near_context(&self, other: &Family) -> bool {
        let (l, r) = (
            other.before.len().min(self.before.len()),
            other.after.len().min(self.after.len()),
        );
        mismatches(
            &other.before[other.before.len() - l..],
            &self.before[self.before.len() - l..],
        ) + mismatches(&other.after[..r], &self.after[..r])
            <= FLANK_MISMATCHES
    }

    /// A weak flank-indel spelling must match the observed context on both
    /// sides of the run. One shared short flank alone can link distinct loci.
    fn near_indel_context(&self, other: &Family) -> bool {
        if self.selected.abs_diff(other.selected) > MAX_LENGTH_SHIFT {
            return false;
        }
        let (l, r) = (
            other.before.len().min(self.before.len()),
            other.after.len().min(self.after.len()),
        );
        (0..=FLANK_MISMATCHES).any(|left_edits| {
            super::similar_paths(
                &self.before[self.before.len() - l..],
                &other.before[other.before.len() - l..],
                left_edits,
            ) && super::similar_paths(
                &self.after[..r],
                &other.after[..r],
                FLANK_MISMATCHES - left_edits,
            )
        })
    }
}

#[derive(Default)]
pub(super) struct Families {
    families: Vec<Family>,
    removed_kmers: usize,
    enabled: bool,
}

/// Mean pair-record support over k-mers holding the whole run and both full
/// flanks. A missing k-mer counts as zero, so a broken length cannot win.
fn mean(observed: &Observed, length: usize, k: usize) -> f64 {
    observed.support as f64 / (k + 1 - length - 2 * FLANK) as f64
}

fn stronger(families: &[Family], a: usize, b: usize) -> Ordering {
    families[a]
        .strength
        .total_cmp(&families[b].strength)
        .then(b.cmp(&a))
}

/// Keep one spelling per tract: the consensus flanks at the consensus family's
/// selected length. Other lengths, flank variants and the run k-mers that only
/// longer lengths contain are removed from both strands.
pub(super) fn collapse(edges: &mut Vec<Edge>, k: usize, fraction: f64) -> Families {
    let mut observed = BTreeMap::<Key, BTreeMap<usize, Observed>>::new();
    let mut spelled = Vec::new();
    for (id, edge) in edges.iter().enumerate() {
        let sequence = graph_kmer::decode(edge.code, k);
        let found = runs(&sequence);
        if found.is_empty() {
            continue;
        }
        for &(base, start, length) in &found {
            if start < FLANK || start + length + FLANK > k {
                continue;
            }
            let key = (
                base,
                flank(&sequence[start - FLANK..start]),
                flank(&sequence[start + length..start + length + FLANK]),
            );
            let entry = observed.entry(key).or_default().entry(length).or_default();
            entry.support += u64::from(edge.count);
            entry.kmers += 1;
            if edge.count > entry.example.0 {
                entry.example = (edge.count, sequence.clone());
            }
            if start > entry.before.len() {
                entry.before = sequence[..start].to_vec();
            }
            if k - start - length > entry.after.len() {
                entry.after = sequence[start + length..].to_vec();
            }
        }
        spelled.push((id, sequence, found));
    }
    // Every flank pair is a candidate tract, so a haplotype seen at one run
    // length still protects its own k-mers below.
    let mut families = Vec::new();
    for (key, lengths) in observed {
        let (&selected, best) = lengths
            .iter()
            .max_by(|(a, x), (b, y)| {
                mean(x, **a, k)
                    .total_cmp(&mean(y, **b, k))
                    .then(x.support.cmp(&y.support))
                    .then(b.cmp(a))
            })
            .expect("one length");
        families.push(Family {
            key,
            strength: mean(best, selected, k),
            before: best.before.clone(),
            after: best.after.clone(),
            lengths: lengths
                .into_iter()
                .map(|(length, o)| Length {
                    length,
                    mean: mean(&o, length, k),
                    kmers: o.kmers,
                    example: o.example.1,
                })
                .collect(),
            selected,
            variant_of: None,
            supported: false,
        });
    }
    let mut left = HashMap::<Side, Vec<usize>>::new();
    let mut right = HashMap::<Side, Vec<usize>>::new();
    for (id, family) in families.iter().enumerate() {
        left.entry((family.key.0, family.key.1))
            .or_default()
            .push(id);
        right
            .entry((family.key.0, family.key.2))
            .or_default()
            .push(id);
    }
    // Each family joins the strongest stronger family it is near; following
    // these links gives its tract. A spelling above the error fraction must also
    // match the longer context. A weak spelling near two real tracts joins one
    // of them, so it cannot merge them.
    let mut parent: Vec<Option<usize>> = vec![None; families.len()];
    for id in 0..families.len() {
        let family = &families[id];
        let (base, before, after) = family.key;
        parent[id] = left[&(base, before)]
            .iter()
            .chain(&right[&(base, after)])
            .copied()
            .filter(|&other| {
                let f = &families[other];
                stronger(&families, other, id) == Ordering::Greater
                    && f.near(&family.key.1, &family.key.2, family.selected)
                    && (family.strength <= fraction * f.strength || f.near_context(family))
            })
            .max_by(|&a, &b| stronger(&families, a, b));
    }
    let mut tracts = BTreeMap::<usize, Vec<usize>>::new();
    for id in 0..families.len() {
        let mut top = id;
        while let Some(next) = parent[top] {
            top = next;
        }
        tracts.entry(top).or_default().push(id);
    }
    // Reads miscall only the side they reach after the run, so a single wrong
    // flank can outscore the correct pair. Choose each flank by the support
    // summed over all spellings of the other flank instead.
    let mut active = vec![false; families.len()];
    for ids in tracts.values() {
        let side = |pick: fn(&Key) -> [u8; FLANK]| {
            let mut support = BTreeMap::<[u8; FLANK], f64>::new();
            for &id in ids {
                *support.entry(pick(&families[id].key)).or_default() += families[id].strength;
            }
            support
                .into_iter()
                .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)))
                .expect("a tract has a family")
                .0
        };
        let (before, after) = (side(|key| key.1), side(|key| key.2));
        let strongest = ids
            .iter()
            .copied()
            .max_by(|&a, &b| stronger(&families, a, b))
            .expect("a tract has a family");
        let consensus = ids
            .iter()
            .copied()
            .find(|&id| families[id].key.1 == before && families[id].key.2 == after)
            .unwrap_or(strongest);
        for &id in ids {
            families[id].variant_of = (id != consensus).then_some(consensus);
        }
        // Only a long run whose consensus flanks show several run lengths
        // collapses; flank spellings alone, as in AT-rich arrays, do not.
        active[consensus] =
            families[consensus].selected >= MIN_RUN && families[consensus].lengths.len() > 1;
    }
    // A run belongs to the tract of every family it is near.
    let tracts_of = |base: u8, before: &[u8], after: &[u8], length: usize| {
        let mut candidates = Vec::new();
        if before.len() >= FLANK {
            candidates.extend(left.get(&(base, flank(&before[before.len() - FLANK..]))));
        }
        if after.len() >= FLANK {
            candidates.extend(right.get(&(base, flank(&after[..FLANK]))));
        }
        let mut tracts: Vec<_> = candidates
            .into_iter()
            .flatten()
            .copied()
            .filter(|&f| families[f].near(before, after, length))
            .map(|f| families[f].variant_of.unwrap_or(f))
            .collect();
        tracts.sort_unstable();
        tracts.dedup();
        tracts
    };
    // A k-mer is removed only when a collapsing tract claims it and no tract
    // keeps it as its own spelling: a collapsing tract, or a distinct tract above
    // the error fraction of the claim, such as another haplotype.
    let mut removed = graph_kmer::Set::default();
    for (id, sequence, found) in &spelled {
        for &(base, start, length) in found {
            let (before, after) = (&sequence[..start], &sequence[start + length..]);
            let tracts = tracts_of(base, before, after, length);
            let claim = tracts
                .iter()
                .filter(|&&t| active[t])
                .map(|&t| families[t].strength)
                .fold(0.0, f64::max);
            let kept = tracts.iter().any(|&t| {
                families[t].accepts(before, after, length)
                    && (active[t] || families[t].strength > fraction * claim)
            });
            if claim > 0.0 && !kept {
                removed.insert(edges[*id].code);
                removed.insert(graph_kmer::reverse_code(edges[*id].code, k));
            }
        }
    }
    edges.retain(|edge| !removed.contains(&edge.code));
    // A tract's own minor lengths and its flank variants are compared with its
    // selected support. A tract within the error fraction of a family of another
    // tract sharing one exact flank is an error branch only when both observed
    // contexts also match within a few edits (for example a flank indel). A
    // distinct locus sharing only the short flank retains its own evidence.
    for id in 0..families.len() {
        let top = families[id].variant_of.unwrap_or(id);
        if !active[top] {
            continue;
        }
        let family = &families[id];
        let support = if top == id {
            family.strongest_minor().map_or(0.0, |l| l.mean)
        } else {
            family.strength
        };
        let (base, before, after) = family.key;
        let dominated = left[&(base, before)]
            .iter()
            .chain(&right[&(base, after)])
            .any(|&other| {
                families[other].variant_of.unwrap_or(other) != id
                    && family.strength <= fraction * families[other].strength
                    && family.near_indel_context(&families[other])
            });
        if support > fraction * families[top].strength && !(top == id && dominated) {
            families[top].supported = true;
        }
    }
    // Report and mark only collapsing tracts, renumbered in key order.
    let kept: Vec<_> = (0..families.len())
        .filter(|&id| active[families[id].variant_of.unwrap_or(id)])
        .collect();
    let index: HashMap<_, _> = kept
        .iter()
        .enumerate()
        .map(|(new, &old)| (old, new))
        .collect();
    let mut slots: Vec<_> = families.into_iter().map(Some).collect();
    let families = kept
        .iter()
        .map(|&old| {
            let mut family = slots[old].take().expect("kept once");
            family.variant_of = family.variant_of.map(|v| index[&v]);
            family
        })
        .collect();
    Families {
        families,
        removed_kmers: removed.len() / 2,
        enabled: true,
    }
}

impl Families {
    pub fn removed_kmers(&self) -> usize {
        self.removed_kmers
    }

    fn tracts(&self) -> impl Iterator<Item = (usize, &Family)> {
        self.families
            .iter()
            .enumerate()
            .filter(|(_, f)| f.variant_of.is_none())
    }

    /// A tract whose minor length or flank variant exceeds the error fraction is
    /// retained evidence: like a permissive bubble, it vetoes final publication
    /// of its raw graph component.
    pub fn alternatives(&self) -> Vec<bubbles::Collapsed> {
        self.tracts()
            .filter(|(_, f)| f.supported)
            .map(|(id, f)| {
                let selected = f.selected();
                let variants = self
                    .families
                    .iter()
                    .filter(|v| v.variant_of == Some(id))
                    .map(Family::selected);
                let minor = f
                    .strongest_minor()
                    .into_iter()
                    .chain(variants)
                    .max_by(|a, b| a.mean.total_cmp(&b.mean))
                    .expect("a supported tract has an alternative");
                bubbles::Collapsed::homopolymer(
                    minor.example.clone(),
                    selected.example.clone(),
                    minor.mean,
                    selected.mean,
                )
            })
            .collect()
    }

    /// Tracts with supported alternatives (selected run and both flanks), located
    /// by their exact spelling in either orientation, including across a circle
    /// origin. This local interval replaces the component-wide bubble flag.
    pub fn regions(&self, candidate: &Candidate) -> Vec<Region> {
        let n = candidate.sequence.len();
        let mut sequence = candidate.sequence.clone();
        let mut regions = Vec::new();
        for (i, family) in self.tracts().filter(|(_, f)| f.supported) {
            let (base, left, right) = &family.key;
            let pattern = [&left[..], &vec![*base; family.selected], &right[..]].concat();
            if candidate.circular {
                sequence.truncate(n);
                sequence.extend(
                    candidate
                        .sequence
                        .iter()
                        .copied()
                        .cycle()
                        .take(pattern.len() - 1),
                );
            }
            let alternatives = family.lengths.len() - 1
                + self
                    .families
                    .iter()
                    .filter(|v| v.variant_of == Some(i))
                    .count();
            for spelling in [pattern.clone(), kmer::reverse_complement(&pattern)] {
                for (start, _) in sequence
                    .windows(spelling.len())
                    .enumerate()
                    .filter(|(_, w)| *w == spelling.as_slice())
                {
                    // The run and both flanks: flank variants were collapsed too.
                    let (from, width) = (start % n, spelling.len());
                    let segments = if from + width <= n {
                        vec![(from, from + width)]
                    } else {
                        vec![(0, from + width - n), (from, n)]
                    };
                    for (start, end) in segments {
                        regions.push(Region {
                            start,
                            end,
                            evidence: format!("homopolymer_family_{}", i + 1),
                            reason: "homopolymer_length_family",
                            alternatives,
                        });
                    }
                }
            }
        }
        regions
    }

    pub fn report(&self) -> Value {
        json!({"enabled":self.enabled,
            "method":"k-mers holding a run and both flanks, grouped by exact flank pair; each pair selects the run length with the highest mean pair-record support; a weaker pair with one flank exact and the other within the mismatch limit (above the error fraction, also over the longest observed context) is a flank variant of the strongest such pair; each tract takes the flanks with the most support summed over the other side and removes its other lengths and flank variants, keeping k-mers that spell a collapsing or distinct supported tract",
            "minimum_selected_run":MIN_RUN,"minimum_variant_run":MIN_VARIANT_RUN,"flank_bases":FLANK,
            "flank_variant_mismatches":FLANK_MISMATCHES,"flank_variant_length_shift":MAX_LENGTH_SHIFT,
            "supported_alternative_rule":"a minor length or flank variant above max_error_fraction of the selected support is retained like a permissive bubble: it vetoes final publication and marks primary sequence; a weak tract is excluded as a flank error only with an exact shared flank, compatible run lengths and matching observed contexts on both sides",
            "removed_kmers":self.removed_kmers,
            "families":self.families.iter().enumerate().map(|(i,f)| json!({
                "family":i+1,"base":(f.key.0 as char).to_string(),
                "left_flank":String::from_utf8_lossy(&f.key.1),"right_flank":String::from_utf8_lossy(&f.key.2),
                "selected_length":f.selected,"selected_mean_kmer_support":(f.strength*1000.0).round()/1000.0,
                "lengths":f.lengths.iter().map(|l| json!({"length":l.length,"mean_kmer_support":(l.mean*1000.0).round()/1000.0,"kmers":l.kmers})).collect::<Vec<_>>(),
                "max_minor_ratio":f.strongest_minor().map(|l| (l.mean/f.strength*1000.0).round()/1000.0),
                "flank_variant_of":f.variant_of.map(|v| v+1),
                "status":if f.variant_of.is_some() {"flank_variant_removed"} else if f.supported {"supported_alternatives"} else {"alternatives_within_error_fraction"},
            })).collect::<Vec<_>>(),
            "limitation":"modal spelling from short-read k-mer support; PCR and sequencing stutter shift run lengths and reads crossing a long run miscall the next bases, so this is a consensus, not a measured length heteroplasmy or a variant call; repeat copies with identical flanks share one tract"})
    }
}

#[cfg(test)]
mod tests {
    use super::super::{compact, graph};
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

    /// Random flanks that cannot extend an A or C run.
    fn flanks(seed: u64) -> (Vec<u8>, Vec<u8>) {
        let (mut left, mut right) = (dna(150, seed), dna(150, seed + 1));
        left[149] = b'G';
        right[0] = b'T';
        (left, right)
    }

    fn tract(left: &[u8], base: u8, length: usize, right: &[u8]) -> Vec<u8> {
        [left, &vec![base; length], right].concat()
    }

    fn codes(edges: &[Edge]) -> graph_kmer::Set {
        edges.iter().map(|e| e.code).collect()
    }

    fn present(codes: &graph_kmer::Set, sequence: &[u8], k: usize) -> bool {
        let mut all = true;
        graph_kmer::visit(sequence, None, k, 0, |_, code| all &= codes.contains(&code));
        all
    }

    fn counts(variants: &[(Vec<u8>, u32)], k: usize) -> graph_kmer::Map<u32> {
        let mut counts = graph_kmer::Map::default();
        for (sequence, n) in variants {
            graph_kmer::visit(sequence, None, k, 0, |_, code| {
                *counts.entry(graph_kmer::canonical(code, k)).or_default() += n;
            });
        }
        counts
    }

    fn spellings(edges: &[Edge], k: usize) -> Vec<Vec<u8>> {
        let (unitigs, _) = compact(edges, k);
        let mut spellings: Vec<_> = unitigs.into_iter().map(|u| u.sequence).collect();
        spellings.sort();
        spellings
    }

    #[test]
    fn interleaved_run_lengths_collapse_to_the_supported_length_on_both_strands() {
        let k = 63;
        let (left, right) = flanks(7);
        for reverse in [false, true] {
            let variants: Vec<_> = [(23, 50), (22, 20), (24, 15), (21, 5), (26, 4)]
                .into_iter()
                .map(|(length, n)| {
                    let sequence = tract(&left, b'A', length, &right);
                    let sequence = if reverse {
                        kmer::reverse_complement(&sequence)
                    } else {
                        sequence
                    };
                    (sequence, n)
                })
                .collect();
            let (mut edges, _, _) = graph(&counts(&variants, k), k, 3);
            let families = collapse(&mut edges, k, 0.2);
            let expected = tract(&left, b'A', 23, &right);
            let mut both = vec![expected.clone(), kmer::reverse_complement(&expected)];
            both.sort();
            assert_eq!(spellings(&edges, k), both, "reverse={reverse}");
            assert_eq!(families.families.len(), 1);
            let family = &families.families[0];
            assert_eq!(family.selected, 23);
            assert!(family.supported);
            assert_eq!(
                family.lengths.iter().map(|l| l.length).collect::<Vec<_>>(),
                vec![21, 22, 23, 24, 26]
            );
            let alternatives = families.alternatives();
            assert_eq!(alternatives.len(), 1);
            assert_eq!(families.report()["families"][0]["max_minor_ratio"], 0.4);
        }
    }

    #[test]
    fn stutter_below_the_error_fraction_collapses_without_retained_alternatives() {
        let k = 63;
        let (left, right) = flanks(3);
        let variants = vec![
            (tract(&left, b'C', 12, &right), 60),
            (tract(&left, b'C', 11, &right), 6),
            (tract(&left, b'C', 13, &right), 4),
        ];
        let (mut edges, _, _) = graph(&counts(&variants, k), k, 3);
        let families = collapse(&mut edges, k, 0.2);
        assert_eq!(families.families[0].selected, 12);
        assert!(!families.families[0].supported);
        assert!(families.alternatives().is_empty());
        assert_eq!(
            families.report()["families"][0]["status"],
            "alternatives_within_error_fraction"
        );
        let expected = tract(&left, b'C', 12, &right);
        let mut both = vec![expected.clone(), kmer::reverse_complement(&expected)];
        both.sort();
        assert_eq!(spellings(&edges, k), both);
    }

    #[test]
    fn read_through_flank_variants_join_the_dominant_tract() {
        let k = 63;
        let (left, right) = flanks(13);
        let substitute = |sequence: &[u8], positions: &[usize]| {
            let mut sequence = sequence.to_vec();
            for &p in positions {
                sequence[p] = if sequence[p] == b'G' { b'C' } else { b'G' };
            }
            sequence
        };
        // Substitutions after the run on each side, as reads crossing it make.
        let (left_variant, right_variant) =
            (substitute(&left, &[148, 145]), substitute(&right, &[3, 8]));
        let variants = vec![
            (tract(&left, b'A', 23, &right), 60),
            (tract(&left, b'A', 22, &right), 20),
            (tract(&left_variant, b'A', 23, &right), 25),
            (tract(&left_variant, b'A', 24, &right), 8),
            (tract(&left, b'A', 23, &right_variant), 18),
            (tract(&left, b'A', 22, &right_variant), 6),
            (tract(&left_variant, b'A', 23, &right_variant), 4),
            (tract(&left_variant, b'A', 22, &right_variant), 3),
        ];
        let (mut edges, _, _) = graph(&counts(&variants, k), k, 3);
        let families = collapse(&mut edges, k, 0.2);
        let roots: Vec<_> = families.tracts().map(|(_, f)| f.selected).collect();
        assert_eq!(roots, vec![23]);
        assert!(families.families.len() >= 3);
        let codes = codes(&edges);
        assert!(present(&codes, &tract(&left, b'A', 23, &right), k));
        // No remaining k-mer holds a run of the tract with a variant flank.
        for edge in &edges {
            let sequence = graph_kmer::decode(edge.code, k);
            for (base, start, length) in runs(&sequence) {
                assert_eq!(base, b'A');
                let consensus = tract(&left, b'A', 23, &right);
                let window = consensus
                    .windows(k)
                    .any(|w| w == sequence.as_slice() || w == kmer::reverse_complement(&sequence));
                assert!(
                    window,
                    "kept a variant k-mer with a {length}-base run at {start}"
                );
            }
        }
        let report = families.report();
        assert!(report["families"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f["status"] == "flank_variant_removed"));
        assert_eq!(families.alternatives().len(), 1);
    }

    #[test]
    fn each_flank_is_chosen_by_support_summed_over_the_other_flank() {
        let k = 63;
        let (left, right) = flanks(47);
        let mut wrong_left = left.clone();
        wrong_left[148] = if wrong_left[148] == b'G' { b'C' } else { b'G' };
        let (mut wrong_right, mut other_right) = (right.clone(), right.clone());
        wrong_right[3] = if wrong_right[3] == b'G' { b'C' } else { b'G' };
        other_right[8] = if other_right[8] == b'G' { b'C' } else { b'G' };
        // One miscalled left flank outscores the correct pair, but the correct
        // left flank also carries both right-flank miscalls.
        let mut variants = Vec::new();
        for (l, r, n) in [
            (&left, &right, 12),
            (&wrong_left, &right, 14),
            (&left, &wrong_right, 8),
            (&left, &other_right, 7),
        ] {
            variants.push((tract(l, b'A', 23, r), n));
            variants.push((tract(l, b'A', 22, r), 4));
        }
        let (mut edges, _, _) = graph(&counts(&variants, k), k, 3);
        let families = collapse(&mut edges, k, 0.2);
        let tracts: Vec<_> = families.tracts().map(|(_, f)| f.key).collect();
        assert_eq!(tracts.len(), 1);
        assert_eq!(
            (tracts[0].1.to_vec(), tracts[0].2.to_vec()),
            (left[138..].to_vec(), right[..12].to_vec())
        );
        assert!(present(&codes(&edges), &tract(&left, b'A', 23, &right), k));
    }

    #[test]
    fn a_second_haplotype_differing_beyond_the_flank_is_kept() {
        let k = 63;
        let (left, right) = flanks(53);
        // Three differences within 12 bases of the run and more beyond it.
        let mut other = left.clone();
        for p in [139, 142, 143, 130, 125, 120, 115] {
            other[p] = if other[p] == b'G' { b'C' } else { b'G' };
        }
        let variants = vec![
            (tract(&left, b'A', 10, &right), 150),
            (tract(&left, b'A', 9, &right), 4),
            (tract(&other, b'A', 9, &right), 50),
        ];
        let (mut edges, _, _) = graph(&counts(&variants, k), k, 3);
        let families = collapse(&mut edges, k, 0.2);
        let codes = codes(&edges);
        assert!(present(&codes, &tract(&left, b'A', 10, &right), k));
        assert!(present(&codes, &tract(&other, b'A', 9, &right), k));
        assert!(families.alternatives().is_empty());
    }

    #[test]
    fn weak_variants_and_dominated_tracts_supply_no_supported_alternatives() {
        let k = 63;
        let (left, right) = flanks(43);
        let mut variant = right.clone();
        variant[3] = if variant[3] == b'G' { b'C' } else { b'G' };
        // One extra base after the run: no longer a near flank, still one locus.
        let shifted = [b"C".as_slice(), &right].concat();
        let variants = vec![
            (tract(&left, b'A', 14, &right), 200),
            (tract(&left, b'A', 13, &right), 10),
            // A weak flank variant whose own two lengths are balanced.
            (tract(&left, b'A', 14, &variant), 6),
            (tract(&left, b'A', 15, &variant), 5),
            (tract(&left, b'A', 14, &shifted), 4),
            (tract(&left, b'A', 13, &shifted), 3),
        ];
        let (mut edges, _, _) = graph(&counts(&variants, k), k, 3);
        let families = collapse(&mut edges, k, 0.2);
        assert_eq!(families.families.len(), 3);
        assert_eq!(families.tracts().count(), 2);
        assert!(families.families.iter().all(|f| !f.supported));
        assert!(families.alternatives().is_empty());
    }

    #[test]
    fn distinct_tracts_sharing_one_flank_keep_their_local_minor_evidence() {
        let k = 63;
        for share_left in [false, true] {
            let (left, right) = flanks(43);
            let (mut other_left, mut other_right) = flanks(171);
            if share_left {
                other_left[138..].copy_from_slice(&left[138..]);
                assert!(mismatches(&right[..FLANK], &other_right[..FLANK]) > FLANK_MISMATCHES);
            } else {
                other_right[..FLANK].copy_from_slice(&right[..FLANK]);
                assert!(mismatches(&left[138..], &other_left[138..]) > FLANK_MISMATCHES);
            }
            let variants = vec![
                (tract(&left, b'A', 18, &right), 1000),
                (tract(&other_left, b'A', 14, &other_right), 40),
                (tract(&other_left, b'A', 13, &other_right), 30),
            ];
            let input = counts(&variants, k);
            let (mut edges, _, _) = graph(&input, k, 3);
            let families = collapse(&mut edges, k, 0.2);
            let active: Vec<_> = families.tracts().collect();
            assert_eq!(active.len(), 1);
            assert_eq!(active[0].1.selected, 14);
            assert!(!present(&codes(&edges), &variants[2].0, k));
            assert!(
                active[0].1.supported,
                "30/40 exceeds the local error fraction"
            );
            assert_eq!(families.report()["families"][0]["max_minor_ratio"], 0.75);
            assert_eq!(
                families.report()["families"][0]["status"],
                "supported_alternatives"
            );
            let alternatives = families.alternatives();
            assert_eq!(alternatives.len(), 1);
            let components = bubbles::Components::new(&input, &alternatives, k, 3);
            for sequence in [
                variants[1].0.clone(),
                kmer::reverse_complement(&variants[1].0),
            ] {
                let candidate = Candidate {
                    sequence,
                    circular: false,
                };
                assert_eq!(components.affected(&candidate, k), vec![1]);
                assert_eq!(families.regions(&candidate).len(), 1);
            }
        }
    }

    #[test]
    fn short_selected_runs_and_single_lengths_are_not_families() {
        let k = 63;
        let (left, right) = flanks(19);
        let (other_left, other_right) = flanks(29);
        let variants = vec![
            (tract(&left, b'A', 8, &right), 50),
            (tract(&left, b'A', 9, &right), 20),
            (tract(&other_left, b'C', 15, &other_right), 50),
        ];
        let (mut edges, _, _) = graph(&counts(&variants, k), k, 3);
        let before = edges.len();
        let families = collapse(&mut edges, k, 0.2);
        assert!(families.families.is_empty());
        assert_eq!(edges.len(), before);
    }

    #[test]
    fn regions_locate_the_selected_tract_in_either_orientation_and_across_the_origin() {
        let k = 63;
        let (left, right) = flanks(37);
        let variants = vec![
            (tract(&left, b'A', 15, &right), 40),
            (tract(&left, b'A', 14, &right), 20),
        ];
        let (mut edges, _, _) = graph(&counts(&variants, k), k, 3);
        let families = collapse(&mut edges, k, 0.2);
        let genome = tract(&left, b'A', 15, &right);
        let linear = Candidate {
            sequence: kmer::reverse_complement(&genome),
            circular: false,
        };
        let regions = families.regions(&linear);
        assert_eq!(regions.len(), 1);
        assert_eq!((regions[0].start, regions[0].end), (138, 177));
        assert_eq!(regions[0].reason, "homopolymer_length_family");
        assert_eq!(regions[0].alternatives, 1);
        let circle = Candidate {
            sequence: [&genome[157..], &genome[..157]].concat(),
            circular: true,
        };
        let mut wrapped: Vec<_> = families
            .regions(&circle)
            .iter()
            .map(|r| (r.start, r.end))
            .collect();
        wrapped.sort_unstable();
        let n = genome.len();
        assert_eq!(wrapped, vec![(0, 20), (n - 19, n)]);
    }
}
