//! Stratify already generated sample paths without using reference lengths or
//! annotation. Shared, exact sequence blocks identify traversal counts within
//! each k; unresolved homology between different k graphs is not assumed.
use super::{compact_at, Edge};
use crate::kmer::{self, graph as g};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

pub(crate) struct Candidate<'a> {
    pub sequence: &'a [u8],
    pub k: usize,
    pub circular: bool,
}

#[derive(Clone, Default, Eq, Ord, PartialEq, PartialOrd)]
struct Signature {
    k: usize,
    copies: Vec<(usize, usize)>,
    order: Vec<(usize, bool)>,
}

fn expanded(candidate: &Candidate<'_>) -> Vec<u8> {
    let mut sequence = candidate.sequence.to_vec();
    if candidate.circular {
        sequence.extend(
            candidate
                .sequence
                .iter()
                .copied()
                .cycle()
                .take(candidate.k - 1),
        );
    }
    sequence
}

fn signatures(candidates: &[Candidate<'_>]) -> Vec<Signature> {
    let mut result = vec![Signature::default(); candidates.len()];
    let ks: BTreeSet<_> = candidates.iter().map(|p| p.k).collect();
    for k in ks {
        let group: Vec<_> = candidates
            .iter()
            .enumerate()
            .filter(|(_, p)| p.k == k)
            .collect();
        let mut codes = BTreeSet::new();
        let mut boundaries = HashSet::new();
        for &(_, p) in &group {
            g::visit(&expanded(p), None, k, 0, |offset, c| {
                codes.insert(c);
                codes.insert(g::reverse_code(c, k));
                if offset == 0 {
                    boundaries.insert(c.prefix());
                    boundaries.insert(g::reverse_code(c.prefix(), k - 1));
                }
                if !p.circular && offset + k == p.sequence.len() {
                    boundaries.insert(c.suffix(k - 1));
                    boundaries.insert(g::reverse_code(c.suffix(k - 1), k - 1));
                }
            });
        }
        let edges: Vec<_> = codes
            .into_iter()
            .map(|code| Edge {
                code,
                from: code.prefix(),
                to: code.suffix(k - 1),
                count: 1,
            })
            .collect();
        let (units, _) = compact_at(&edges, k, &boundaries);
        let canonical: Vec<_> = units
            .iter()
            .map(|u| {
                let reverse = kmer::reverse_complement(&u.sequence);
                if reverse < u.sequence {
                    (reverse, true)
                } else {
                    (u.sequence.clone(), false)
                }
            })
            .collect();
        let ids: BTreeMap<_, _> = canonical
            .iter()
            .map(|(s, _)| s.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .enumerate()
            .map(|(i, s)| (s, i))
            .collect();
        let mut tokens = HashMap::new();
        for (u, (key, reverse)) in units.iter().zip(&canonical) {
            g::visit(&u.sequence, None, k, 0, |offset, c| {
                tokens.insert(c, (ids[key], *reverse, offset));
            });
        }
        let mut repeated = HashSet::new();
        for &(i, p) in &group {
            let signature = &mut result[i];
            signature.k = k;
            let mut visits = BTreeMap::<usize, usize>::new();
            // Boundaries split each candidate's start/end. Each first edge is
            // one complete block visit, including repeated self-loop visits.
            g::visit(&expanded(p), None, k, 0, |_, c| {
                if let Some(&(id, reverse, 0)) = tokens.get(&c) {
                    signature.order.push((id, reverse));
                    *visits.entry(id).or_default() += 1;
                }
            });
            repeated.extend(visits.iter().filter_map(|(&id, &n)| (n > 1).then_some(id)));
            signature.copies = visits.into_iter().collect();
        }
        for &(i, _) in &group {
            // Single-copy sequence bubbles are order/sequence alternatives,
            // not automatically new repeat-count families.
            result[i].copies.retain(|(id, _)| repeated.contains(id));
        }
    }
    result
}

pub(crate) struct Selection {
    pub indices: Vec<usize>,
    pub report: Value,
}

pub(crate) fn select(
    candidates: &[Candidate<'_>],
    coverage_order: &[usize],
    budget: usize,
) -> Selection {
    let signatures = signatures(candidates);
    let copy_key = |i: usize| (signatures[i].k, signatures[i].copies.clone());
    let mut copies = BTreeSet::new();
    let mut orders = BTreeSet::new();
    let mut selected = BTreeSet::new();
    let mut indices = Vec::new();
    // Preserve the coverage incumbent, then give each count family a turn,
    // then each distinct order; only then spend slots on further variants.
    for pass in 0..3 {
        for &i in coverage_order {
            if indices.len() == budget {
                break;
            }
            if selected.contains(&i)
                || (pass == 0 && copies.contains(&copy_key(i)))
                || (pass == 1 && orders.contains(&signatures[i]))
            {
                continue;
            }
            selected.insert(i);
            indices.push(i);
            copies.insert(copy_key(i));
            orders.insert(signatures[i].clone());
        }
    }
    let copy_count = (0..signatures.len())
        .map(copy_key)
        .collect::<BTreeSet<_>>()
        .len();
    let order_count = signatures.iter().collect::<BTreeSet<_>>().len();
    let report = json!({
        "method":"repeat_block_copies_then_order_then_sequence",
        "sample_path_budget":budget,"generated_paths":candidates.len(),"scored_paths":indices.len(),
        "copy_families":copy_count,"scored_copy_families":copies.len(),
        "order_families":order_count,"scored_order_families":orders.len(),
        "unscored_copy_families":copy_count.saturating_sub(copies.len()),
        "unscored_order_families":order_count.saturating_sub(orders.len()),
        "coverage_ranks":indices.iter().map(|i| coverage_order.iter().position(|x| x==i).unwrap()+1).collect::<Vec<_>>(),
        "baseline_is_additional":true,
        "scope":"exact sequence blocks in the generated path union within each k; graph traversal families, not proven biological copy counts"
    });
    Selection { indices, report }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random(n: usize) -> Vec<u8> {
        let mut state = 909177u64;
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn repeat_copies_are_not_starved_by_coverage_ranked_sequence_variants() {
        let source = random(300);
        let repeat = &source[80..140];
        let mut sequences = Vec::new();
        for variant in 0..9 {
            let mut s = [&source[..80], repeat, repeat, &source[140..]].concat();
            s[220 + variant] = if s[220 + variant] == b'A' { b'C' } else { b'A' };
            sequences.push(s);
        }
        sequences.push([&source[..80], repeat, repeat, repeat, &source[140..]].concat());
        let inputs: Vec<_> = sequences
            .iter()
            .map(|s| Candidate {
                sequence: s,
                k: 15,
                circular: false,
            })
            .collect();
        let selection = select(&inputs, &(0..10).collect::<Vec<_>>(), 8);
        assert_eq!(selection.indices[0], 0);
        assert!(selection.indices.contains(&9), "{}", selection.report);
        assert_eq!(selection.indices.len(), 8);
        assert!(!(0..8).contains(&9));
    }

    #[test]
    fn equal_edge_counts_keep_different_orders() {
        let source = random(300);
        let (s, r, a, b, t) = (
            &source[..50],
            &source[50..100],
            &source[100..160],
            &source[160..220],
            &source[220..],
        );
        let paths = [
            [s, r, a, r, b, r, t].concat(),
            [s, r, b, r, a, r, t].concat(),
        ];
        let inputs: Vec<_> = paths
            .iter()
            .map(|p| Candidate {
                sequence: p,
                k: 15,
                circular: false,
            })
            .collect();
        let sig = signatures(&inputs);
        assert_eq!(sig[0].copies, sig[1].copies);
        assert_ne!(sig[0].order, sig[1].order);
        let selected = select(&inputs, &[0, 1], 2);
        assert_eq!(selected.report["copy_families"], 1);
        assert_eq!(selected.report["order_families"], 2);
        assert_eq!(selected.indices, [0, 1]);
    }
}
