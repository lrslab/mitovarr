//! Remove weak disconnected sequence shadows using retained sample support.
//! Counts and connected alleles are untouched by this pass; no reference or
//! expected length participates. Low-support isolated homologs can be lost.
use super::{compact, graph_kmer, Edge, Unitig};
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};

pub(super) fn remove(
    edges: &mut Vec<Edge>,
    unitigs: Vec<Unitig>,
    k: usize,
    fraction: f64,
    retained_fragments: &HashSet<graph_kmer::Code>,
) -> (Vec<Unitig>, usize, Value) {
    if fraction <= 0.0 {
        let branches = compact(edges, k).1;
        return (unitigs, branches, Value::Null);
    }
    let n = unitigs.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let mut by_node: HashMap<graph_kmer::Code, usize> = HashMap::new();
    for (id, u) in unitigs.iter().enumerate() {
        let mut keys = vec![
            graph_kmer::canonical(u.from, k - 1),
            graph_kmer::canonical(u.to, k - 1),
        ];
        // Opposite circular unitigs can start at different vertices. Their
        // minimum canonical node still joins both strands into one component.
        if u.circular {
            let mut minimum = keys[0];
            graph_kmer::visit(&u.sequence, None, k - 1, 0, |_, node| {
                minimum = minimum.min(graph_kmer::canonical(node, k - 1));
            });
            keys.push(minimum);
        }
        for key in keys {
            if let Some(&other) = by_node.get(&key) {
                let (a, b) = (find(&mut parent, id), find(&mut parent, other));
                parent[a] = b;
            } else {
                by_node.insert(key, id);
            }
        }
    }
    let component: Vec<usize> = (0..n).map(|id| find(&mut parent, id)).collect();
    let mut maximum = HashMap::<usize, u32>::new();
    for (id, u) in unitigs.iter().enumerate() {
        let m = maximum.entry(component[id]).or_default();
        *m = (*m).max(u.max_count);
    }
    let mut order: Vec<usize> = maximum.keys().copied().collect();
    order.sort_by_key(|c| (std::cmp::Reverse(maximum[c]), *c));
    let mut members = HashMap::<usize, Vec<usize>>::new();
    for (id, &group) in component.iter().enumerate() {
        members.entry(group).or_default().push(id);
    }
    // Strongest retained support per canonical 31-mer, filled in decreasing
    // component strength so a shadow can only point at a retained competitor.
    let mut strongest = HashMap::<crate::kmer::Code, (f64, usize)>::new();
    let mut removed = HashSet::new();
    let (mut removed_kmers, mut removed_components) = (0usize, 0usize);
    let mut protected_weak_link_components = 0usize;
    for c in order {
        let ids = &members[&c];
        let mut competitor = 0.0f64;
        for &id in ids {
            crate::kmer::visit(&unitigs[id].sequence, None, 31, 0, |_, code| {
                if let Some(&(support, owner)) = strongest.get(&crate::kmer::canonical(code, 31)) {
                    if owner != c {
                        competitor = competitor.max(support);
                    }
                }
            });
        }
        let mut protected = false;
        for &id in ids {
            graph_kmer::visit(&unitigs[id].sequence, None, k, 0, |_, code| {
                protected |= retained_fragments.contains(&graph_kmer::canonical(code, k));
            });
        }
        let shadow = competitor > 0.0 && f64::from(maximum[&c]) <= fraction * competitor;
        if shadow && protected {
            protected_weak_link_components += 1;
        }
        if shadow && !protected {
            removed.insert(c);
            removed_components += 1;
            removed_kmers += ids
                .iter()
                .map(|&id| unitigs[id].sequence.len() + 1 - k)
                .sum::<usize>();
            continue;
        }
        for &id in ids {
            let support = unitigs[id].mean_count;
            crate::kmer::visit(&unitigs[id].sequence, None, 31, 0, |_, code| {
                let entry = strongest
                    .entry(crate::kmer::canonical(code, 31))
                    .or_insert((support, c));
                if support > entry.0 {
                    *entry = (support, c);
                }
            });
        }
    }
    let mut codes = HashSet::new();
    for (id, u) in unitigs.iter().enumerate() {
        if removed.contains(&component[id]) {
            graph_kmer::visit(&u.sequence, None, k, 0, |_, code| {
                codes.insert(code);
            });
        }
    }
    edges.retain(|edge| !codes.contains(&edge.code));
    let (unitigs, branches) = compact(edges, k);
    let report = json!({"method":"disconnected component max support <= max_error_fraction x mean support of a retained unitig sharing a canonical 31-mer",
        "removed_components":removed_components,"protected_weak_link_components":protected_weak_link_components,"protection_rule":"retain weak-link remnants with a stronger unitig between the same full flanks and a length difference above three bases","component_scope":"canonical (k-1) nodes; both strands together","max_error_fraction":fraction,"shared_seed_bases":31,"removed_directed_kmers":removed_kmers,
        "retained_canonical_kmers":edges.len()/2,"retained_directed_branch_nodes":branches});
    (unitigs, branches, report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{denovo::graph, kmer};

    fn remove(
        edges: &mut Vec<Edge>,
        unitigs: Vec<Unitig>,
        k: usize,
        fraction: f64,
    ) -> (Vec<Unitig>, usize, Value) {
        super::remove(edges, unitigs, k, fraction, &HashSet::new())
    }

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

    fn counts(sequences: &[(&[u8], u32)], k: usize) -> HashMap<graph_kmer::Code, u32> {
        let mut result = HashMap::new();
        for &(sequence, support) in sequences {
            graph_kmer::visit(sequence, None, k, 0, |_, code| {
                *result.entry(graph_kmer::canonical(code, k)).or_default() += support;
            });
        }
        result
    }

    fn codes(edges: &[Edge]) -> Vec<graph_kmer::Code> {
        edges.iter().map(|e| e.code).collect()
    }

    fn shadow(major: &[u8]) -> Vec<u8> {
        let mut sequence = major[100..300].to_vec();
        // Share 31-mers but no complete (k-1) node with the main graph.
        for i in (0..sequence.len()).step_by(50) {
            sequence[i] = if sequence[i] == b'A' { b'C' } else { b'A' };
        }
        sequence
    }

    #[test]
    fn removes_error_islands_symmetrically_without_changing_raw_counts() {
        let major = dna(500, 91427);
        let minor = shadow(&major);
        let unrelated = dna(180, 9216);
        for reverse in [false, true] {
            let sequences: Vec<_> = [&major, &minor, &unrelated]
                .into_iter()
                .map(|s| {
                    if reverse {
                        kmer::reverse_complement(s)
                    } else {
                        s.clone()
                    }
                })
                .collect();
            let input = counts(
                &[
                    (&sequences[0], 100),
                    (&sequences[1], 20),
                    (&sequences[2], 3),
                ],
                63,
            );
            let original = input.clone();
            let (mut edges, unitigs, _) = graph(&input, 63, 3);
            let (_, _, report) = remove(&mut edges, unitigs, 63, 0.2);
            assert_eq!(report["removed_components"], 1);
            assert_eq!(report["removed_directed_kmers"], 2 * (minor.len() + 1 - 63));
            let (expected, _, _) = graph(&counts(&[(&major, 100), (&unrelated, 3)], 63), 63, 3);
            assert_eq!(codes(&edges), codes(&expected));
            assert_eq!(input, original);
        }
    }

    #[test]
    fn retains_strong_and_thirty_percent_disconnected_homologs() {
        let major = dna(500, 91427);
        let minor = shadow(&major);
        for support in [30, 100] {
            let input = counts(&[(&major, 100), (&minor, support)], 63);
            let (mut edges, unitigs, _) = graph(&input, 63, 3);
            let before = codes(&edges);
            let (_, _, report) = remove(&mut edges, unitigs, 63, 0.2);
            assert_eq!(report["removed_components"], 0);
            assert_eq!(codes(&edges), before);
        }
    }

    #[test]
    fn connected_fifteen_percent_bubble_and_ten_percent_deletion_are_not_shadows() {
        let major = dna(700, 123491);
        let mut snp = major.clone();
        snp[300] = if snp[300] == b'A' { b'C' } else { b'A' };
        let deletion = [&major[..220], &major[400..]].concat();
        for (minor, support) in [(&snp, 15), (&deletion, 10)] {
            let input = counts(&[(&major, 100), (minor, support)], 63);
            let (mut edges, unitigs, branches) = graph(&input, 63, 3);
            assert!(branches > 0);
            let before = codes(&edges);
            let (_, after_branches, report) = remove(&mut edges, unitigs, 63, 0.2);
            assert_eq!(report["removed_components"], 0);
            assert_eq!(codes(&edges), before);
            assert_eq!(after_branches, branches);
        }
    }

    #[test]
    fn zero_fraction_does_not_remove_or_report_measured_shadows() {
        let major = dna(500, 91427);
        let minor = shadow(&major);
        let (mut edges, unitigs, _) = graph(&counts(&[(&major, 100), (&minor, 3)], 63), 63, 3);
        let before = codes(&edges);
        let (_, _, report) = remove(&mut edges, unitigs, 63, 0.0);
        assert!(report.is_null());
        assert_eq!(codes(&edges), before);
    }

    #[test]
    fn circular_strands_are_removed_as_one_component() {
        let major = dna(1000, 514279);
        let circle = major[200..250].to_vec();
        // k=63 needs more than one circumference of this small circle.
        let circular_sequence: Vec<_> = circle
            .iter()
            .cycle()
            .take(circle.len() + 62)
            .copied()
            .collect();
        let (mut edges, unitigs, _) = graph(
            &counts(&[(&major, 100), (&circular_sequence, 5)], 63),
            63,
            3,
        );
        assert!(unitigs.iter().filter(|u| u.circular).count() >= 2);
        let (_, _, report) = remove(&mut edges, unitigs, 63, 0.2);
        assert_eq!(report["removed_components"], 1);
        assert_eq!(report["removed_directed_kmers"], 100);
    }
    #[test]
    fn ten_percent_deletion_junction_survives_the_complete_cleanup_pipeline() {
        let major = dna(700, 123491);
        let deletion = [&major[..220], &major[400..]].concat();
        let input = counts(&[(&major, 100), (&deletion, 10)], 63);
        let mut junction = None;
        graph_kmer::visit(&deletion[190..253], None, 63, 0, |_, code| {
            junction = Some(code)
        });
        let junction = junction.unwrap();
        let (mut edges, _, _) = graph(&input, 63, 3);
        let mut retained = HashSet::new();
        let (unitigs, _, previous, _) =
            crate::denovo::simplify_retaining_links(&mut edges, 63, 0.2, 0.49, &mut retained);
        assert!(
            edges.iter().any(|e| e.code == junction),
            "baseline lost the junction"
        );
        let (_, _, shadows) = super::remove(&mut edges, unitigs, 63, 0.2, &retained);
        assert_eq!(shadows["protected_weak_link_components"], 1);
        assert!(
            edges.iter().any(|e| e.code == junction),
            "cleanup={previous}; shadows={shadows}"
        );
    }
}
