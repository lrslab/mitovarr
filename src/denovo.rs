//! Conservative sample-only de Bruijn graph. No reference or seed input exists here.
pub(crate) mod align;
pub(crate) mod anchor;
mod bubbles;
mod comparison;
pub(crate) mod evidence;
mod layout;
pub(crate) mod multik;
pub(crate) mod pairs;
mod paths;
mod shadows;
pub(crate) mod shortlist;
mod threading;
use crate::{
    cli::{GraphSettings, Limits, Reads},
    io::{fasta::write_fasta_records, fastq},
    kmer::{self, graph as graph_kmer},
};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};
use tracing::info;

#[derive(Debug)]
struct Edge {
    code: graph_kmer::Code,
    from: graph_kmer::Code,
    to: graph_kmer::Code,
    count: u32,
}

#[derive(Debug)]
struct Unitig {
    from: graph_kmer::Code,
    to: graph_kmer::Code,
    sequence: Vec<u8>,
    circular: bool,
    mean_count: f64,
    max_count: u32,
    min_count: u32,
    first_count: u32,
    last_count: u32,
}

#[derive(Clone, Debug)]
struct Candidate {
    sequence: Vec<u8>,
    circular: bool,
}

/// Booth's minimum rotation, O(n), so arbitrary FASTQ order never chooses the origin.
fn rotation_start(sequence: &[u8]) -> usize {
    let n = sequence.len();
    if n == 0 {
        return 0;
    }
    let (mut a, mut b, mut offset) = (0, 1, 0);
    while a < n && b < n && offset < n {
        let (x, y) = (sequence[(a + offset) % n], sequence[(b + offset) % n]);
        if x == y {
            offset += 1;
            continue;
        }
        if x > y {
            a += offset + 1;
            if a <= b {
                a = b + 1;
            }
        } else {
            b += offset + 1;
            if b <= a {
                b = a + 1;
            }
        }
        offset = 0;
    }
    a.min(b)
}

fn min_rotation(sequence: &[u8]) -> Vec<u8> {
    let start = rotation_start(sequence);
    [&sequence[start..], &sequence[..start]].concat()
}

fn normalize(sequence: &[u8], circular: bool) -> Vec<u8> {
    let reverse = kmer::reverse_complement(sequence);
    if circular {
        min_rotation(sequence).min(min_rotation(&reverse))
    } else {
        sequence.to_vec().min(reverse)
    }
}

fn walk(
    first: usize,
    edges: &[Edge],
    outgoing: &graph_kmer::Map<(usize, usize)>,
    incoming: &graph_kmer::Map<usize>,
    visited: &mut [bool],
    k: usize,
) -> Unitig {
    let mut sequence = graph_kmer::decode(edges[first].from, k - 1);
    let mut edge_id = first;
    let (mut total_count, mut max_count, mut length) = (0u64, 0u32, 0usize);
    let mut min_count = u32::MAX;
    loop {
        visited[edge_id] = true;
        let edge = &edges[edge_id];
        total_count += u64::from(edge.count);
        max_count = max_count.max(edge.count);
        min_count = min_count.min(edge.count);
        length += 1;
        sequence.push(edge.code.last_base());
        let next = outgoing.get(&edge.to);
        if next.is_none_or(|run| run.1 != 1) || incoming.get(&edge.to) != Some(&1) {
            return Unitig {
                from: edges[first].from,
                to: edge.to,
                sequence,
                circular: false,
                mean_count: total_count as f64 / length as f64,
                max_count,
                min_count,
                first_count: edges[first].count,
                last_count: edge.count,
            };
        }
        let next = next.expect("one outgoing edge").0;
        if visited[next] {
            return Unitig {
                from: edges[first].from,
                to: edge.to,
                sequence,
                circular: next == first,
                mean_count: total_count as f64 / length as f64,
                max_count,
                min_count,
                first_count: edges[first].count,
                last_count: edge.count,
            };
        }
        edge_id = next;
    }
}

fn graph(
    counts: &graph_kmer::Map<u32>,
    k: usize,
    min_count: u32,
) -> (Vec<Edge>, Vec<Unitig>, usize) {
    let mut edges = Vec::new();
    for (&code, &count) in counts.iter().filter(|(_, count)| **count >= min_count) {
        for code in [code, graph_kmer::reverse_code(code, k)] {
            edges.push(Edge {
                code,
                from: code.prefix(),
                to: code.suffix(k - 1),
                count,
            });
        }
    }
    edges.sort_by_key(|e| e.code);
    let (unitigs, branch_nodes) = compact(&edges, k);
    (edges, unitigs, branch_nodes)
}

fn compact(edges: &[Edge], k: usize) -> (Vec<Unitig>, usize) {
    compact_at(edges, k, &HashSet::new())
}

fn compact_at(
    edges: &[Edge],
    k: usize,
    boundaries: &HashSet<graph_kmer::Code>,
) -> (Vec<Unitig>, usize) {
    // Per node: (first outgoing edge id, outgoing edge count). Walks follow an
    // edge only when it is the node's single outgoing edge.
    let mut outgoing: graph_kmer::Map<(usize, usize)> =
        graph_kmer::Map::with_capacity_and_hasher(edges.len(), Default::default());
    let mut incoming: graph_kmer::Map<usize> =
        graph_kmer::Map::with_capacity_and_hasher(edges.len(), Default::default());
    for (id, edge) in edges.iter().enumerate() {
        outgoing
            .entry(edge.from)
            .and_modify(|run| run.1 += 1)
            .or_insert((id, 1));
        *incoming.entry(edge.to).or_default() += 1;
    }
    // Split unitigs at requested sample-derived corridor boundaries.
    for &node in boundaries {
        incoming.insert(node, 2);
    }
    let branch_nodes = outgoing
        .iter()
        .filter(|(node, run)| run.1 > 1 || incoming.get(node).is_some_and(|n| *n > 1))
        .count()
        + incoming
            .iter()
            .filter(|(node, n)| **n > 1 && !outgoing.contains_key(node))
            .count();
    let mut visited = vec![false; edges.len()];
    let mut unitigs = Vec::new();
    for (id, edge) in edges.iter().enumerate() {
        if !visited[id] && (outgoing[&edge.from].1 != 1 || incoming.get(&edge.from) != Some(&1)) {
            unitigs.push(walk(id, edges, &outgoing, &incoming, &mut visited, k));
        }
    }
    for id in 0..edges.len() {
        if !visited[id] {
            unitigs.push(walk(id, edges, &outgoing, &incoming, &mut visited, k));
        }
    }
    (unitigs, branch_nodes)
}

/// Small bounded edit distance for homologous local bubble paths, not repeat scaffolding.
/// Two reused rows; only the band and the cells it reads are ever relevant.
fn similar_paths(a: &[u8], b: &[u8], max_edits: usize) -> bool {
    if a.len().abs_diff(b.len()) > max_edits {
        return false;
    }
    let mut row: Vec<usize> = (0..=b.len()).collect();
    let mut next = vec![max_edits + 1; b.len() + 1];
    for (i, x) in a.iter().enumerate() {
        let start = i.saturating_sub(max_edits);
        next[0] = i + 1;
        // The band's left neighbour must read as out of range, not a stale value.
        if start > 0 {
            next[start] = max_edits + 1;
        }
        for j in start..b.len().min(i + max_edits + 1) {
            next[j + 1] = (row[j] + usize::from(*x != b[j]))
                .min(row[j + 1] + 1)
                .min(next[j] + 1);
        }
        std::mem::swap(&mut row, &mut next);
    }
    row[b.len()] <= max_edits
}

/// Compressed adjacency lists indexed by unitig id.
struct Lists {
    offsets: Vec<usize>,
    ids: Vec<usize>,
}

impl Lists {
    fn get(&self, id: usize) -> &[usize] {
        &self.ids[self.offsets[id]..self.offsets[id + 1]]
    }
}

/// Per-round bubble-search graph, built once for all candidate unitigs.
struct BubbleGraph<'a> {
    unitigs: &'a [Unitig],
    starts: &'a graph_kmer::Map<Vec<usize>>,
    /// Each unitig's (reverse of `to`, reverse of `from`): its twin's endpoints.
    reverse: Vec<(graph_kmer::Code, graph_kmer::Code)>,
    /// Unitigs starting at this unitig's end node (`starts[to]`, same order).
    successors: Lists,
    /// Unitigs whose reverse endpoints equal this unitig's endpoints.
    conjugates: Lists,
    k: usize,
}

/// State of the walk being expanded; marks are all zero between searches.
struct Walker {
    /// Path membership counts by unitig id.
    marks: Vec<u32>,
    path: Vec<usize>,
    pending: Vec<(usize, usize, f64)>,
    states: usize,
}

struct BubbleSearch<'a> {
    graph: BubbleGraph<'a>,
    walker: Walker,
}

/// The weak unitig and limits of one bubble search.
struct WeakPath<'w> {
    id: usize,
    unitig: &'w Unitig,
    reverse_from: graph_kmer::Code,
    reverse_to: graph_kmer::Code,
    max_length: usize,
    max_edits: usize,
    fraction: f64,
    removed: &'w [bool],
}

impl<'a> BubbleSearch<'a> {
    fn new(unitigs: &'a [Unitig], starts: &'a graph_kmer::Map<Vec<usize>>, k: usize) -> Self {
        let reverse: Vec<_> = unitigs
            .iter()
            .map(|u| {
                (
                    graph_kmer::reverse_code(u.to, k - 1),
                    graph_kmer::reverse_code(u.from, k - 1),
                )
            })
            .collect();
        let lists = |select: &dyn Fn(usize) -> Vec<usize>| {
            let mut offsets = vec![0];
            let mut ids = Vec::new();
            for id in 0..unitigs.len() {
                ids.extend(select(id));
                offsets.push(ids.len());
            }
            Lists { offsets, ids }
        };
        let starting = |node: &graph_kmer::Code| starts.get(node).map_or(&[][..], Vec::as_slice);
        let successors = lists(&|id| starting(&unitigs[id].to).to_vec());
        let conjugates = lists(&|id| {
            let (from, to) = reverse[id];
            starting(&from)
                .iter()
                .copied()
                .filter(|&p| unitigs[p].to == to)
                .collect()
        });
        Self {
            graph: BubbleGraph {
                unitigs,
                starts,
                reverse,
                successors,
                conjugates,
                k,
            },
            walker: Walker {
                marks: vec![0; unitigs.len()],
                path: Vec::new(),
                pending: Vec::new(),
                states: 0,
            },
        }
    }

    /// Compare complete compacted paths between the same junctions. Coverage is
    /// averaged over path edges, so a single high/low k-mer does not veto a bubble.
    /// Only short, similar spellings qualify; loops and conjugate traversals remain.
    fn stronger_path(
        &mut self,
        weak_id: usize,
        removed: &[bool],
        fraction: f64,
    ) -> Option<bubbles::Alternative> {
        let graph = &self.graph;
        let unitig = &graph.unitigs[weak_id];
        let (reverse_from, reverse_to) = graph.reverse[weak_id];
        let weak = WeakPath {
            id: weak_id,
            unitig,
            reverse_from,
            reverse_to,
            max_length: 3 * graph.k,
            max_edits: 3.max((unitig.sequence.len() - graph.k + 1) / 10),
            fraction,
            removed,
        };
        self.walker.states = 0;
        let first = graph
            .starts
            .get(&unitig.from)
            .map_or(&[][..], Vec::as_slice);
        // A bounded simple-path search on unitigs, rather than one visit per k-mer
        // node: different coverage histories reaching a node must remain distinct.
        graph
            .expand(&mut self.walker, &weak, first, 0, 0.0)
            .flatten()
    }
}

impl BubbleGraph<'_> {
    /// Generate (and count) every child of the current walk in order, then explore
    /// accepted children last-first: the order of the former LIFO state stack.
    /// `Some(result)` ends the search (an alternative, or `None` at the state cap).
    fn expand(
        &self,
        walker: &mut Walker,
        weak: &WeakPath,
        children: &[usize],
        length: usize,
        mass: f64,
    ) -> Option<Option<bubbles::Alternative>> {
        let k = self.k;
        let base = walker.pending.len();
        for &id in children {
            walker.states += 1;
            if walker.states > 4096 {
                walker.pending.truncate(base);
                return Some(None);
            }
            let u = &self.unitigs[id];
            let n = u.sequence.len() - k + 1;
            if id == weak.id
                || u.from == u.to
                || u.circular
                || u.to == weak.unitig.from
                || (u.from == weak.reverse_from && u.to == weak.reverse_to)
                || length + n > weak.max_length
                || self.on_path(walker, id)
            {
                continue;
            }
            if weak.removed[id] {
                continue;
            }
            let new_mass = mass + u.mean_count * n as f64;
            if u.to == weak.unitig.to {
                if weak.unitig.mean_count <= weak.fraction * new_mass / (length + n) as f64 {
                    let mut sequence = graph_kmer::decode(weak.unitig.from, k - 1);
                    for &p in walker.path.iter().chain([&id]) {
                        sequence.extend_from_slice(&self.unitigs[p].sequence[k - 1..]);
                    }
                    if similar_paths(&weak.unitig.sequence, &sequence, weak.max_edits) {
                        let minimum = walker
                            .path
                            .iter()
                            .chain([&id])
                            .map(|&p| self.unitigs[p].min_count)
                            .min()
                            .unwrap();
                        walker.pending.truncate(base);
                        return Some(Some(bubbles::Alternative {
                            sequence,
                            mean_count: new_mass / (length + n) as f64,
                            min_count: minimum,
                        }));
                    }
                }
            } else {
                walker.pending.push((id, length + n, new_mass));
            }
        }
        let end = walker.pending.len();
        for i in (base..end).rev() {
            let (id, length, mass) = walker.pending[i];
            walker.marks[id] += 1;
            walker.path.push(id);
            let done = self.expand(walker, weak, self.successors.get(id), length, mass);
            walker.path.pop();
            walker.marks[id] -= 1;
            if done.is_some() {
                walker.pending.truncate(base);
                return done;
            }
        }
        walker.pending.truncate(base);
        None
    }

    /// A unitig may not end at a node the walk already left (a path unitig starts
    /// there) or traverse the twin of a path unitig.
    fn on_path(&self, walker: &Walker, id: usize) -> bool {
        self.successors.get(id).iter().any(|&p| walker.marks[p] > 0)
            || self.conjugates.get(id).iter().any(|&p| walker.marks[p] > 0)
    }
}

/// A compacted unitig and its reverse-complement twin share a minimum
/// canonical edge, even when a circular twin starts at a different vertex.
/// Each directed edge belongs to exactly one unitig, so this key is unique
/// to that pair (or to one self-conjugate unitig).
fn unitig_twins(unitigs: &[Unitig], k: usize) -> Vec<usize> {
    let mut owners = graph_kmer::Map::default();
    let mut twins: Vec<_> = (0..unitigs.len()).collect();
    for (id, unitig) in unitigs.iter().enumerate() {
        let mut key = None;
        graph_kmer::visit(&unitig.sequence, None, k, 0, |_, code| {
            let code = graph_kmer::canonical(code, k);
            key = Some(key.map_or(code, |old: graph_kmer::Code| old.min(code)));
        });
        if let Some(other) = owners.insert(key.expect("unitig edge"), id) {
            twins[id] = other;
            twins[other] = id;
        }
    }
    twins
}

/// A longer/shorter alternative with the same full graph flanks is retained
/// structural evidence, unlike an arbitrary fragment detached during cleanup.
/// Three bases is the existing conservative sequencing-error edit allowance.
fn has_structural_alternative(
    weak: &Unitig,
    unitigs: &[Unitig],
    starts: &graph_kmer::Map<Vec<usize>>,
    fraction: f64,
) -> bool {
    starts.get(&weak.from).into_iter().flatten().any(|&id| {
        let strong = &unitigs[id];
        !strong.circular
            && strong.to == weak.to
            && strong.sequence.len().abs_diff(weak.sequence.len()) > 3
            && f64::from(weak.max_count) <= fraction * f64::from(strong.min_count)
    })
}

/// Topology-aware error removal uses local competing support, never a genome length,
/// reference, global coverage cutoff, or a choice of route through an unbridged repeat.
#[cfg(test)]
fn simplify(
    edges: &mut Vec<Edge>,
    k: usize,
    fraction: f64,
    bubble_fraction: f64,
) -> (Vec<Unitig>, usize, Value, Vec<bubbles::Collapsed>) {
    simplify_retaining_links(edges, k, fraction, bubble_fraction, &mut HashSet::new())
}

fn simplify_retaining_links(
    edges: &mut Vec<Edge>,
    k: usize,
    fraction: f64,
    bubble_fraction: f64,
    retained_fragments: &mut HashSet<graph_kmer::Code>,
) -> (Vec<Unitig>, usize, Value, Vec<bubbles::Collapsed>) {
    let (mut unitigs, mut branches) = compact(edges, k);
    let (raw_kmers, raw_branches) = (edges.len() / 2, branches);
    let (mut rounds, mut tip_kmers, mut bubble_kmers, mut weak_link_kmers) = (0, 0, 0, 0);
    let mut alternatives = std::collections::BTreeMap::new();
    if fraction > 0.0 {
        loop {
            let mut starts: graph_kmer::Map<Vec<usize>> = graph_kmer::Map::default();
            let mut ends: graph_kmer::Map<Vec<usize>> = graph_kmer::Map::default();
            for (id, unitig) in unitigs.iter().enumerate() {
                starts.entry(unitig.from).or_default().push(id);
                ends.entry(unitig.to).or_default().push(id);
            }
            let mut removed = graph_kmer::Set::default();
            let mut remove = |unitig: &Unitig| {
                graph_kmer::visit(&unitig.sequence, None, k, 0, |_, code| {
                    removed.insert(code);
                    removed.insert(graph_kmer::reverse_code(code, k));
                });
            };
            for unitig in &unitigs {
                if unitig.circular || unitig.sequence.len() - k + 1 > 2 * k {
                    continue;
                }
                let competitor = if !ends.contains_key(&unitig.from) {
                    ends.get(&unitig.to)
                        .into_iter()
                        .flatten()
                        .map(|id| unitigs[*id].last_count)
                        .max()
                } else if !starts.contains_key(&unitig.to) {
                    starts
                        .get(&unitig.from)
                        .into_iter()
                        .flatten()
                        .map(|id| unitigs[*id].first_count)
                        .max()
                } else {
                    None
                };
                if competitor
                    .is_some_and(|count| f64::from(unitig.max_count) <= f64::from(count) * fraction)
                {
                    remove(unitig);
                }
            }
            let tips = removed.len() / 2;
            // Recompact after tip clipping before interpreting bubbles.
            if removed.is_empty() {
                let twins = unitig_twins(&unitigs, k);
                let mut search = BubbleSearch::new(&unitigs, &starts, k);
                let mut unavailable = vec![false; unitigs.len()];
                let mut order: Vec<_> = (0..unitigs.len()).collect();
                order
                    .sort_by_key(|id| (unitigs[*id].max_count, unitigs[*id].from, unitigs[*id].to));
                for id in order {
                    let unitig = &unitigs[id];
                    if unitig.circular
                        || unitig.from == unitig.to
                        || unitig.sequence.len() - k + 1 > 3 * k
                        || starts[&unitig.from].len() < 2
                        || ends[&unitig.to].len() < 2
                    {
                        continue;
                    }
                    if let Some(alternative) =
                        search.stronger_path(id, &unavailable, bubble_fraction)
                    {
                        // The more permissive cleanup may improve the primary,
                        // but must not erase evidence needed for final status.
                        let conservative = unitig.sequence.len() - k < 2 * k
                            && alternative.sequence.len() - k < 2 * k
                            && f64::from(unitig.max_count)
                                <= fraction * f64::from(alternative.min_count)
                            && similar_paths(&unitig.sequence, &alternative.sequence, 3);
                        if !conservative {
                            let bubble = bubbles::Collapsed::new(unitig, alternative);
                            alternatives.entry(bubble.removed.clone()).or_insert(bubble);
                        }
                        unavailable[id] = true;
                        unavailable[twins[id]] = true;
                        graph_kmer::visit(&unitig.sequence, None, k, 0, |_, code| {
                            removed.insert(code);
                            removed.insert(graph_kmer::reverse_code(code, k));
                        });
                    }
                }
            }
            let mut weak_links = 0;
            if removed.is_empty() {
                // Intersecting error paths need not form a simple bubble. Disconnect
                // only their weak boundary kmer, leaving the rest as a separate
                // candidate. Relative support is local; tandem-repeat loops and
                // long competing paths are never pruned by this rule.
                for unitig in &unitigs {
                    if unitig.circular
                        || unitig.from == unitig.to
                        || unitig.sequence.len() - k + 1 > 2 * k
                    {
                        continue;
                    }
                    let outgoing_competitor = starts[&unitig.from]
                        .iter()
                        .map(|id| &unitigs[*id])
                        .filter(|u| u.from != u.to)
                        .map(|u| u.first_count)
                        .max()
                        .unwrap_or(0);
                    let incoming_competitor = ends[&unitig.to]
                        .iter()
                        .map(|id| &unitigs[*id])
                        .filter(|u| u.from != u.to)
                        .map(|u| u.last_count)
                        .max()
                        .unwrap_or(0);
                    for (window, competitor) in [
                        (&unitig.sequence[..k], outgoing_competitor),
                        (
                            &unitig.sequence[unitig.sequence.len() - k..],
                            incoming_competitor,
                        ),
                    ] {
                        if f64::from(unitig.max_count) <= f64::from(competitor) * fraction {
                            // Protect supported structural junction remnants,
                            // not every error branch detached by this pass.
                            if has_structural_alternative(unitig, &unitigs, &starts, fraction) {
                                graph_kmer::visit(&unitig.sequence, None, k, 0, |_, code| {
                                    retained_fragments.insert(graph_kmer::canonical(code, k));
                                });
                            }
                            graph_kmer::visit(window, None, k, 0, |_, code| {
                                removed.insert(code);
                                removed.insert(graph_kmer::reverse_code(code, k));
                            });
                        }
                    }
                }
                weak_links = removed.len() / 2;
            }
            if removed.is_empty() {
                break;
            }
            tip_kmers += tips;
            bubble_kmers += removed.len() / 2 - tips - weak_links;
            weak_link_kmers += weak_links;
            edges.retain(|edge| !removed.contains(&edge.code));
            (unitigs, branches) = compact(edges, k);
            rounds += 1;
        }
    }
    let alternatives: Vec<_> = alternatives.into_values().collect();
    let stats = json!({"retained_bubble_alternatives":bubbles::reports(&alternatives),"bubble_alternatives_file":if alternatives.is_empty(){None}else{Some("bubble_alternatives.fasta")},"method":"short tips, whole-path coverage bubbles, weak boundary disconnection; strand-symmetric local support", "max_error_fraction":fraction,"max_bubble_fraction":bubble_fraction, "max_path_kmers":3*k,"max_bubble_edits":3.max(3*k/10),"bubble_edit_rule":"max(3, weak_path_kmers / 10)","bubble_support":"length-weighted mean pair-record kmer support","max_bubble_search_states":4096,"rounds":rounds,"raw_canonical_kmers":raw_kmers,"raw_directed_branch_nodes":raw_branches,"tip_kmers_removed":tip_kmers,"bubble_kmers_removed":bubble_kmers,"weak_link_kmers_removed":weak_link_kmers,"cleaned_canonical_kmers":edges.len()/2,"cleaned_directed_branch_nodes":branches});
    (unitigs, branches, stats, alternatives)
}

const MIXED_MIN_MINOR: f64 = 0.2;
const MIXED_MIN_BUBBLES: usize = 3;

/// Minor-branch support fraction of each balanced allele bubble left after
/// cleanup: two non-circular unitigs between the same nodes, lengths within three
/// bases, counted once per strand pair. Cleanup keeps such pairs only when both
/// branches are well supported, so several of them signal a second haplotype.
fn balanced_bubbles(unitigs: &[Unitig], k: usize) -> Vec<f64> {
    let mut ends = std::collections::BTreeMap::<_, Vec<usize>>::new();
    for (id, u) in unitigs.iter().enumerate() {
        if !u.circular && u.from != u.to {
            ends.entry((u.from, u.to)).or_default().push(id);
        }
    }
    let mut minors = Vec::new();
    for ((from, to), mut ids) in ends {
        let twin = (
            graph_kmer::reverse_code(to, k - 1),
            graph_kmer::reverse_code(from, k - 1),
        );
        // One strand of each pair: the lexicographically smaller endpoint key.
        if ids.len() < 2 || twin < (from, to) {
            continue;
        }
        ids.sort_by(|&a, &b| {
            unitigs[b]
                .mean_count
                .total_cmp(&unitigs[a].mean_count)
                .then(unitigs[a].sequence.cmp(&unitigs[b].sequence))
        });
        let (major, minor) = (&unitigs[ids[0]], &unitigs[ids[1]]);
        if major.sequence.len().abs_diff(minor.sequence.len()) <= 3 {
            minors.push(minor.mean_count / (major.mean_count + minor.mean_count));
        }
    }
    minors.sort_by(f64::total_cmp);
    minors
}

/// A reference-free mixed-sample signal; not proof of heteroplasmy. Minor
/// sequences cleaned as errors or too divergent to share k-mers are not counted.
/// Supported structural alternatives (e.g. a deletion in part of the molecules)
/// come from path resolution.
fn mixed_sample(unitigs: &[Unitig], k: usize, structural: &[Value]) -> Value {
    let minors = balanced_bubbles(unitigs, k);
    let balanced = minors.iter().filter(|&&m| m >= MIXED_MIN_MINOR).count();
    let signal = balanced >= MIXED_MIN_BUBBLES || !structural.is_empty();
    json!({"method":"allele bubbles left after cleanup: parallel non-circular unitigs between the same nodes, lengths within 3 bases, one per strand pair",
        "balanced_allele_bubbles":balanced,"minimum_minor_fraction":MIXED_MIN_MINOR,"signal_threshold":MIXED_MIN_BUBBLES,
        "minor_fractions":minors.iter().map(|m| (m * 1000.0).round() / 1000.0).collect::<Vec<_>>(),
        "structural_alternatives":structural,
        "structural_method":"two routes from one path source, each with >= the junction fragment minimum and the minor >= 0.2 of both; one route's target returns to the other's target through at most six strict joins, including their intermediate segments; bypassed bases >= 50; one per bypassed length",
        "status":if signal {"mixed_sample_signal"} else {"not_detected"},
        "limitation":"signal, not proof of heteroplasmy or contamination; low-level minor alleles removed as errors and divergent minor sequences are not counted; a mixed sample can fragment the graph and shorten the primary"})
}

fn write_graph(out: &Path, unitigs: &[Unitig], k: usize) -> Result<()> {
    let mut writer = BufWriter::new(File::create(out.join("assembly_graph.gfa"))?);
    writeln!(writer, "H\tVN:Z:1.0")?;
    let mut starts: graph_kmer::Map<Vec<usize>> = graph_kmer::Map::default();
    for (id, unitig) in unitigs.iter().enumerate() {
        writeln!(
            writer,
            "S\tu{}\t{}\tdp:f:{:.3}",
            id + 1,
            String::from_utf8_lossy(&unitig.sequence),
            unitig.mean_count,
        )?;
        starts.entry(unitig.from).or_default().push(id);
    }
    for (id, unitig) in unitigs.iter().enumerate() {
        for next in starts.get(&unitig.to).into_iter().flatten() {
            writeln!(writer, "L\tu{}\t+\tu{}\t+\t{}M", id + 1, next + 1, k - 1)?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Full-query gapped placements in both orientations; count quality-filtered
/// matching bases once per fragment and preserve ambiguous placements.
fn validate(
    reads: &Reads,
    sequence: &[u8],
    supported_graph_closure: bool,
    out: &Path,
    settings: &GraphSettings,
    limits: &Limits,
) -> Result<(Value, bool)> {
    let n = sequence.len();
    let aligner = align::Index::new(sequence, true, sequence.len());
    let mut depth = vec![0u32; n];
    let mut names = HashSet::new();
    let (mut junction_fragments, mut exact_reads, mut repeated_names) = (0u64, 0u64, 0u64);
    let (mut aligned_reads, mut ambiguous_reads, mut alignment_limited_reads) =
        (0usize, 0usize, 0usize);
    let mut observed_depth = vec![0u32; n];
    let mut deletion_depth = vec![0u32; n];
    let anchor = settings.k.min(20);
    let scan = fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
        let name = a.name(1)?;
        if names.contains(name) {
            repeated_names += 1;
            return Ok(());
        }
        let mut covered = Vec::new();
        let mut observed = Vec::new();
        let mut deleted = Vec::new();
        let mut junction = false;
        let mut matches = 0;
        for read in [a, b] {
            let mut placements = aligner.align(read);
            if placements.exhausted {
                alignment_limited_reads += 1;
                continue;
            }
            placements
                .hits
                .sort_by(|a, b| b.logp.total_cmp(&a.logp).then(a.start.cmp(&b.start)));
            let Some(hit) = placements.hits.first() else {
                continue;
            };
            if placements
                .hits
                .get(1)
                .is_some_and(|other| hit.logp - other.logp < 5.0)
            {
                ambiguous_reads += 1;
                continue;
            }
            let trace = align::trace_alignment(read, hit, &aligner.sequence, limits.min_baseq);
            junction |= n - hit.start >= anchor
                && hit.end >= n + anchor
                && trace.matches_contiguous(n - anchor, n + anchor);
            observed.extend(trace.observed.into_iter().map(|p| p % n));
            covered.extend(trace.matching.into_iter().map(|p| p % n));
            deleted.extend(trace.deletions.into_iter().map(|p| p % n));
            matches += 1;
            exact_reads += u64::from(hit.edits == 0);
            aligned_reads += 1;
        }
        if matches > 0 {
            ensure!(
                names.len() < settings.max_validation_pairs,
                "validation exceeds --max-validation-pairs; no final assembly published"
            );
            names.insert(name.to_owned());
            covered.sort_unstable();
            covered.dedup();
            for pos in covered {
                depth[pos] = depth[pos].saturating_add(1);
            }
            observed.sort_unstable();
            observed.dedup();
            for p in observed {
                observed_depth[p] = observed_depth[p].saturating_add(1);
            }
            deleted.sort_unstable();
            deleted.dedup();
            for p in deleted {
                deletion_depth[p] = deletion_depth[p].saturating_add(1);
            }
            junction_fragments += u64::from(junction);
        }
        Ok(())
    })?;
    let minimum = *depth.iter().min().unwrap_or(&0);
    let low_depth = depth.iter().filter(|d| **d < settings.min_depth).count();
    let mut writer = BufWriter::new(File::create(out.join("coverage.tsv"))?);
    writeln!(writer, "pos\tdepth\tmethod")?;
    for (i, d) in depth.iter().enumerate() {
        writeln!(writer, "{}\t{}\tgapped_matching_fragment", i + 1, d)?;
    }
    writer.flush()?;
    let complete = low_depth == 0
        && (supported_graph_closure
            || junction_fragments >= u64::from(settings.min_junction_fragments));
    Ok((
        json!({"method":"unique full-query affine-gap placement; quality-filtered matching bases only; mates collapsed per position","aligned_reads":aligned_reads,"ambiguous_placement_reads":ambiguous_reads,"alignment_limited_reads":alignment_limited_reads,"minimum_observed_base_depth":observed_depth.iter().min(),"positions_with_deletion_support":deletion_depth.iter().filter(|&&n|n>0).count(), "input_r1":reads.r1.canonicalize()?,"input_r2":reads.r2.canonicalize()?,"scan":scan,"exact_reads":exact_reads,"unique_fragments":names.len(),"repeated_names_ignored":repeated_names,"minimum_depth":minimum,"low_depth_bases":low_depth,"literal_junction_fragments":junction_fragments,"read_supported_graph_closure":supported_graph_closure,"closure_basis":if supported_graph_closure{"read-supported graph joins before output rotation"}else{"literal output-origin read bridge"},"junction_anchor_bases":anchor,"passed":complete}),
        complete,
    ))
}

/// Count each canonical k-mer once per complete pair. `None` means the distinct
/// set exceeds `--max-kmers`. With `(seed, n)`, only pairs whose seeded name
/// priority lies in the lowest 1/2^n are counted; the returned scan describes them.
fn count_kmers(
    reads: &Reads,
    settings: &GraphSettings,
    limits: &Limits,
    sample: Option<(u64, u32)>,
) -> Result<Option<(graph_kmer::Map<u32>, fastq::Scan)>> {
    let mut counts = graph_kmer::Map::<u32>::default();
    let mut fragment = graph_kmer::Set::default();
    let mut retained = fastq::Scan { pairs: 0, bases: 0 };
    let mut exceeded = false;
    let scan = fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
        if let Some((seed, halvings)) = sample {
            if !retained_pair(a.name(1)?, seed, halvings) {
                return Ok(());
            }
        }
        retained.pairs += 1;
        retained.bases += (a.sequence.len() + b.sequence.len()) as u64;
        fragment.clear();
        for read in [a, b] {
            graph_kmer::visit(
                &read.sequence,
                Some(&read.quality),
                settings.k,
                limits.min_baseq,
                |_, code| {
                    fragment.insert(graph_kmer::canonical(code, settings.k));
                },
            );
        }
        for &code in &fragment {
            let count = if counts.len() < limits.max_kmers {
                counts.entry(code).or_default()
            } else if let Some(count) = counts.get_mut(&code) {
                count
            } else {
                exceeded = true;
                anyhow::bail!("assembly exceeds --max-kmers");
            };
            *count = count.saturating_add(1);
        }
        Ok(())
    });
    match scan {
        Ok(_) => Ok(Some((counts, retained))),
        Err(_) if exceeded => Ok(None),
        Err(error) => Err(error),
    }
}

const MAX_SUBSAMPLE_HALVINGS: u32 = 6;

fn retained_pair(name: &str, seed: u64, halvings: u32) -> bool {
    crate::subsample::priority(name.as_bytes(), seed).leading_zeros() >= halvings
}

/// Runs only after all complete pairs exceeded `--max-kmers`. Halve the retained
/// pairs until the graph fits, write them for every assembly stage and report the
/// trigger. Final read validation keeps using the caller's complete input.
fn subsample(
    reads: &Reads,
    out: &Path,
    settings: &GraphSettings,
    limits: &Limits,
    seed: u64,
) -> Result<(graph_kmer::Map<u32>, fastq::Scan, Option<Reads>, Value)> {
    let mut attempts = vec![json!({"fraction":1.0,"status":"exceeds_max_kmers"})];
    for halvings in 1..=MAX_SUBSAMPLE_HALVINGS {
        let fraction = 0.5f64.powi(halvings as i32);
        let Some((counts, scan)) = count_kmers(reads, settings, limits, Some((seed, halvings)))?
        else {
            attempts.push(json!({"fraction":fraction,"status":"exceeds_max_kmers"}));
            continue;
        };
        attempts.push(json!({"fraction":fraction,"status":"within_max_kmers",
            "retained_pairs":scan.pairs,"distinct_kmers":counts.len()}));
        let sampled = Reads {
            r1: out.join("subsampled.R1.fastq"),
            r2: out.join("subsampled.R2.fastq"),
        };
        let mut one = BufWriter::new(File::create(&sampled.r1)?);
        let mut two = BufWriter::new(File::create(&sampled.r2)?);
        let original = fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
            if retained_pair(a.name(1)?, seed, halvings) {
                a.write(&mut one)?;
                b.write(&mut two)?;
            }
            Ok(())
        })?;
        one.flush()?;
        two.flush()?;
        info!(
            original_pairs = original.pairs,
            retained_pairs = scan.pairs,
            fraction,
            "all pairs exceeded --max-kmers; assembling deterministic subsample"
        );
        let report = json!({"enabled":true,"applied":true,
            "trigger":{"limit":"--max-kmers","max_kmers":limits.max_kmers,
                "observed":"distinct k-mers from all complete pairs exceeded the limit"},
            "method":"complete pairs whose seeded fragment-name priority lies in the lowest 1/2^n; mates kept together; independent of read order",
            "seed":seed,"fraction":fraction,"attempts":attempts,
            "original_pairs":original.pairs,"original_bases":original.bases,
            "retained_pairs":scan.pairs,"retained_bases":scan.bases,
            "original_r1":reads.r1.canonicalize()?,"original_r2":reads.r2.canonicalize()?,
            "retained_r1":"subsampled.R1.fastq","retained_r2":"subsampled.R2.fastq",
            "evidence_scope":"graph construction, read threading, path scoring and optional layouts use retained pairs; final read validation uses all complete pairs",
            "limitation":"retained depth is lower everywhere; low-frequency alleles and weak junction or repeat support can be lost from graph assembly"});
        return Ok((counts, scan, Some(sampled), report));
    }
    anyhow::bail!(
        "assembly exceeds --max-kmers even after subsampling to 1/{} of complete pairs; no final assembly published",
        1u32 << MAX_SUBSAMPLE_HALVINGS
    )
}

#[allow(clippy::too_many_arguments)]
pub fn assemble(
    reads: &Reads,
    validation_reads: &Reads,
    out: &Path,
    settings: &GraphSettings,
    limits: &Limits,
    recruitment: Option<Value>,
    bam_baited: bool,
    subsample_seed: Option<u64>,
) -> Result<bool> {
    ensure!(
        (15..=127).contains(&settings.k) && settings.k % 2 == 1,
        "assembly k must be odd in 15..=127"
    );
    ensure!(
        settings.max_error_fraction.is_finite()
            && (0.0..0.5).contains(&settings.max_error_fraction),
        "--max-error-fraction must be at least 0 and less than 0.5"
    );
    ensure!(
        settings.max_bubble_fraction.is_finite()
            && (0.0..0.5).contains(&settings.max_bubble_fraction),
        "--max-bubble-fraction must be at least 0 and less than 0.5"
    );
    ensure!(
        settings.min_kmer_count > 0
            && settings.min_depth > 0
            && settings.min_junction_fragments > 0
            && limits.max_kmers > 0
            && settings.max_validation_pairs > 0,
        "assembly thresholds and limits must be positive"
    );
    let (counts, scan, sampled, subsampling) = match count_kmers(reads, settings, limits, None)? {
        Some((counts, scan)) => (
            counts,
            scan,
            None,
            json!({"enabled":subsample_seed.is_some(),"applied":false}),
        ),
        None => {
            let seed = subsample_seed.ok_or_else(|| {
                anyhow::anyhow!("assembly exceeds --max-kmers; no final assembly published")
            })?;
            subsample(reads, out, settings, limits, seed)?
        }
    };
    let reads = sampled.as_ref().unwrap_or(reads);
    info!(
        pairs = scan.pairs,
        distinct_kmers = counts.len(),
        "building sample-only graph"
    );
    let (mut edges, _, _) = graph(&counts, settings.k, settings.min_kmer_count);
    let mut retained_fragments = HashSet::new();
    let (unitigs, branches, mut simplification, collapsed_bubbles) = simplify_retaining_links(
        &mut edges,
        settings.k,
        settings.max_error_fraction,
        settings.max_bubble_fraction,
        &mut retained_fragments,
    );
    let (unitigs, branches) = if settings.max_error_fraction > 0.0 {
        let (unitigs, branches, report) = shadows::remove(
            &mut edges,
            unitigs,
            settings.k,
            settings.max_error_fraction,
            &retained_fragments,
        );
        simplification["shadow_components"] = report;
        (unitigs, branches)
    } else {
        simplification["shadow_components"] = Value::Null;
        (unitigs, branches)
    };
    simplification["cleaned_canonical_kmers"] = json!(edges.len() / 2);
    simplification["cleaned_directed_branch_nodes"] = json!(branches);
    info!(
        raw_branches = simplification["raw_directed_branch_nodes"].as_u64(),
        branches, "simplified sample-only graph"
    );
    bubbles::write(out, &collapsed_bubbles)?;
    let bubble_components = bubbles::Components::new(
        &counts,
        &collapsed_bubbles,
        settings.k,
        settings.min_kmer_count,
    );
    write_graph(out, &unitigs, settings.k)?;
    let threaded = threading::candidates(&unitigs, reads, settings, limits)?;
    let unresolved_branches = threaded.unresolved;
    let read_threading = threaded.report;
    let mut path_resolution = paths::resolve(&threaded.graph, reads, out, settings, limits)?;
    let mixed_sample = mixed_sample(
        &unitigs,
        settings.k,
        &path_resolution.structural_alternatives,
    );
    for primary in &mut path_resolution.primary_candidates {
        let affected = bubble_components.affected(&primary.candidate, settings.k);
        if !affected.is_empty() {
            primary.regions.push(paths::primary::Region {
                start: 0,
                end: primary.candidate.sequence.len(),
                evidence: format!(
                    "collapsed_bubbles_{}",
                    affected
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                reason: "collapsed_bubble_component_unresolved",
                alternatives: affected.len(),
            });
        }
    }
    let candidates = path_resolution.candidates;
    let records: Vec<_> = candidates
        .iter()
        .enumerate()
        .map(|(i, c)| {
            (
                format!(
                    "candidate_{} circular={} length={}",
                    i + 1,
                    c.circular,
                    c.sequence.len()
                ),
                c.sequence.clone(),
            )
        })
        .collect();
    write_fasta_records(&out.join("candidates.fasta"), &records)?;
    // The supported main circle is judged independently of residual components.
    // Canonical k-mer mass counts each graph observation once, not once per path
    // or repeat copy; length alone is never the selection criterion.
    let circles: Vec<_> = candidates
        .iter()
        .enumerate()
        .filter(|(_, c)| c.circular && c.sequence.len() >= settings.k)
        .map(|(i, _)| i)
        .collect();
    let total_mass: u64 = edges.iter().map(|edge| u64::from(edge.count)).sum::<u64>() / 2;
    let mut component_fraction = 0.0;
    let selected_id = if circles.len() == 1 {
        let id = circles[0];
        let c = &candidates[id];
        let extended = [c.sequence.as_slice(), &c.sequence[..settings.k - 1]].concat();
        let mut codes = HashSet::new();
        graph_kmer::visit(&extended, None, settings.k, 0, |_, code| {
            codes.insert(graph_kmer::canonical(code, settings.k));
        });
        let mass: u64 = codes
            .iter()
            .map(|code| u64::from(counts.get(code).copied().unwrap_or(0)))
            .sum();
        component_fraction = mass as f64 / total_mass.max(1) as f64;
        Some(id)
    } else {
        None
    };
    let selected = selected_id.map(|id| &candidates[id]);
    let selected_bubbles = selected
        .map(|c| bubble_components.affected(c, settings.k))
        .unwrap_or_default();
    let (mut validation, read_validated) = if let Some(candidate) = selected {
        validate(
            validation_reads,
            &candidate.sequence,
            path_resolution
                .supported_cycles
                .contains(&candidate.sequence),
            out,
            settings,
            limits,
        )?
    } else {
        (Value::Null, false)
    };
    let end_search_limited = path_resolution.report["end_assembly"]["status"] == "resource_limited";
    let complete = read_validated && selected_bubbles.is_empty() && !end_search_limited;
    if let Some(object) = validation.as_object_mut() {
        object.insert(
            "unresolved_collapsed_bubbles".into(),
            json!(selected_bubbles),
        );
    }
    let (primary, best_circular) = paths::primary::write_outputs(
        path_resolution.primary_candidates,
        &counts,
        settings.k,
        out,
        selected.filter(|_| complete),
    )?;
    let coverage_layout = if settings.coverage_layout && !primary.is_null() {
        let backbone = crate::io::fasta::load_reference(&out.join("primary.fasta"), None)?;
        layout::propose(
            &backbone.sequence,
            primary["circular"].as_bool().unwrap_or(false),
            reads,
            out,
            settings,
            limits,
        )?
    } else {
        json!({"enabled":settings.coverage_layout,"status":if settings.coverage_layout {"no_primary_backbone"} else {"disabled"}})
    };
    let reason = if complete {
        "supported_dominant_circle"
    } else if end_search_limited {
        "end_search_resource_limit"
    } else if !selected_bubbles.is_empty() {
        "supported_bubble_alternatives"
    } else if selected.is_some() {
        if bam_baited {
            "insufficient_extracted_read_or_junction_support"
        } else {
            "insufficient_original_read_or_junction_support"
        }
    } else if edges.is_empty() {
        "no_supported_kmers"
    } else if circles.len() > 1 {
        "competing_circular_candidates"
    } else if unresolved_branches {
        "unresolved_graph_branches"
    } else if candidates.len() != 1 {
        "multiple_components"
    } else if !candidates[0].circular {
        "linear_or_incomplete_graph"
    } else {
        "cycle_shorter_than_k"
    };
    let mut ambiguity = BufWriter::new(File::create(out.join("ambiguity.tsv"))?);
    writeln!(ambiguity, "candidate\tstart\tend\treason")?;
    if !complete {
        for (i, candidate) in candidates.iter().enumerate() {
            writeln!(
                ambiguity,
                "candidate_{}\t1\t{}\t{}",
                i + 1,
                candidate.sequence.len(),
                reason
            )?;
        }
    }
    ambiguity.flush()?;
    let report = json!({
        "schema":"mitovar-rust.assembly.v1", "backend":"rust-de-bruijn", "workflow":if bam_baited {"bam_baited_de_novo"} else if recruitment.is_some(){"seed_recruited_de_novo"}else{"reference_blind_de_novo"},
        "status":if complete{"complete"}else{"unresolved"},"reason":reason,"final_published":complete,
        "parameters":{"k":settings.k,"min_kmer_count":settings.min_kmer_count,"max_error_fraction":settings.max_error_fraction,"max_bubble_fraction":settings.max_bubble_fraction,"min_baseq":limits.min_baseq,"min_depth":settings.min_depth,"min_junction_fragments":settings.min_junction_fragments,"max_kmers":limits.max_kmers,"max_validation_pairs":settings.max_validation_pairs},
        "input_scope":if bam_baited {"mt-bam-extracted"} else {"original-paired-fastq"},"validation_scope":if bam_baited {"mt-bam-extracted"} else {"original-complete-paired-fastq"},
        "assembly_r1":reads.r1.canonicalize()?,"assembly_r2":reads.r2.canonicalize()?,"assembly_reads":scan,"subsampling":subsampling,
        "distinct_kmers":counts.len(),"retained_canonical_kmers":edges.len()/2,"minimum_retained_pair_count":edges.iter().map(|e|e.count).min(),"directed_branch_nodes":branches,
        "graph_simplification":simplification,"mixed_sample":mixed_sample,
        "read_threading":read_threading,
        "paired_path_resolution":path_resolution.report,
        "primary":primary,
        "best_circular":best_circular,
        "sequence_output":{"status":if complete {"complete"} else if !primary.is_null() {"partial_available"} else {"unavailable"},
            "fasta":if !primary.is_null() {Some("primary.fasta")} else {None},
            "gene_quality":"not_assessed", "ncr1_required_for_partial_output":false},
        "coverage_layout":coverage_layout,
        "candidate_count":candidates.len(),"candidates":candidates.iter().enumerate().map(|(i,c)|json!({"id":format!("candidate_{}",i+1),"length":c.sequence.len(),"circular":c.circular})).collect::<Vec<_>>(),
        "selected_candidate":selected_id.map(|id|format!("candidate_{}",id+1)),"selection_rule":"one circular candidate at least k bases, no unresolved collapsed alternatives in its raw graph component, with gapped read-supported coverage and closure; graph k-mer mass is diagnostic, other small components retained without global veto", "selected_component_support_fraction":component_fraction,
        "candidate_count_before_paired_paths":threaded.candidates.len(),
        "recruitment":recruitment,"validation":validation,"origin":"lexicographically minimal cyclic rotation across both strands",
        "limitations":["local error cleanup can remove low-frequency biological alleles; this is a primary-consensus assembler, not a heteroplasmy caller","optional coverage layouts estimate repeat copies and remain separate from final assembly; no annotation or reference repair","complete means this graph and gapped-read support criteria passed; mitochondrial identity requires independent annotation"],
    });
    crate::write_json(&out.join("report.json"), &report)?;
    if complete {
        write_fasta_records(
            &out.join("consensus.fasta"),
            &[(
                "mitochondrion".into(),
                selected.expect("validated selection").sequence.clone(),
            )],
        )?;
    }
    info!(
        reason,
        complete,
        candidates = candidates.len(),
        "assembly finished"
    );
    Ok(complete)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rotation_matches_exhaustive_search() {
        for sequence in [
            b"ACGACG".as_slice(),
            b"TTAGCCGA",
            b"AAAAAAAA",
            b"TGATG",
            b"A",
        ] {
            let expected = (0..sequence.len())
                .map(|i| [&sequence[i..], &sequence[..i]].concat())
                .min()
                .unwrap();
            assert_eq!(min_rotation(sequence), expected);
        }
    }

    /// The prior per-row allocating kernel, kept as a differential oracle.
    fn similar_paths_reference(a: &[u8], b: &[u8], max_edits: usize) -> bool {
        if a.len().abs_diff(b.len()) > max_edits {
            return false;
        }
        let mut row: Vec<usize> = (0..=b.len()).collect();
        for (i, x) in a.iter().enumerate() {
            let mut next = vec![max_edits + 1; b.len() + 1];
            next[0] = i + 1;
            for j in i.saturating_sub(max_edits)..b.len().min(i + max_edits + 1) {
                next[j + 1] = (row[j] + usize::from(*x != b[j]))
                    .min(row[j + 1] + 1)
                    .min(next[j] + 1);
            }
            row = next;
        }
        row[b.len()] <= max_edits
    }

    #[test]
    fn reused_row_edit_distance_matches_the_allocating_kernel() {
        let mut state = 0x9e3779b97f4a7c15u64;
        let mut next = |n: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % n
        };
        let (mut accepted, mut rejected) = (0, 0);
        for _ in 0..4000 {
            let length = 1 + next(160) as usize;
            let a: Vec<u8> = (0..length).map(|_| b"ACGT"[next(4) as usize]).collect();
            let mut b = a.clone();
            for _ in 0..next(10) {
                let at = next(b.len() as u64 + 1) as usize;
                match next(3) {
                    0 if at < b.len() => b[at] = b"ACGT"[next(4) as usize],
                    1 if at < b.len() => {
                        b.remove(at);
                    }
                    _ => b.insert(at, b"ACGT"[next(4) as usize]),
                }
            }
            for max_edits in 0..=8 {
                let expected = similar_paths_reference(&a, &b, max_edits);
                assert_eq!(similar_paths(&a, &b, max_edits), expected);
                assert_eq!(
                    similar_paths(&b, &a, max_edits),
                    similar_paths_reference(&b, &a, max_edits)
                );
                if expected {
                    accepted += 1;
                } else {
                    rejected += 1;
                }
            }
        }
        // Both outcomes are exercised near the edit bound.
        assert!(accepted > 1000 && rejected > 1000);
    }
}

#[cfg(test)]
mod cleanup_tests {
    use super::*;

    fn dna(n: usize) -> Vec<u8> {
        let mut state = 872398u64;
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 32) as usize & 3]
            })
            .collect()
    }

    fn counts(sequences: &[(&[u8], u32)], k: usize) -> graph_kmer::Map<u32> {
        let mut result = graph_kmer::Map::default();
        for &(s, n) in sequences {
            graph_kmer::visit(s, None, k, 0, |_, c| {
                *result.entry(graph_kmer::canonical(c, k)).or_default() += n;
            });
        }
        result
    }

    #[test]
    fn wider_bubbles_use_mean_support_and_preserve_balanced_alternatives() {
        let major = dna(650);
        let mut minor = major.clone();
        for i in [200, 220, 240, 260, 280] {
            minor[i] = if minor[i] == b'A' { b'C' } else { b'A' };
        }
        let k = 63;
        let expected = counts(&[(&major, 1)], k);
        let mut input = counts(&[(&major, 100), (&minor, 15)], k);
        // One locally inflated minor k-mer must not dominate whole-arm support.
        graph_kmer::visit(&minor[220..220 + k], None, k, 0, |_, c| {
            input.insert(graph_kmer::canonical(c, k), 35);
        });
        let (mut edges, _, _) = graph(&input, k, 3);
        let (_, branches, stats, _) = simplify(&mut edges, k, 0.2, 0.49);
        assert_eq!(branches, 0);
        assert!(stats["bubble_kmers_removed"].as_u64().unwrap() > 126);
        let actual: HashSet<_> = edges
            .iter()
            .map(|e| graph_kmer::canonical(e.code, k))
            .collect();
        assert_eq!(actual, expected.keys().copied().collect());
        let reversed = counts(
            &[
                (&kmer::reverse_complement(&major), 100),
                (&kmer::reverse_complement(&minor), 15),
            ],
            k,
        );
        let (mut reverse_edges, _, _) = graph(&reversed, k, 3);
        simplify(&mut reverse_edges, k, 0.2, 0.49);
        assert_eq!(
            actual,
            reverse_edges
                .iter()
                .map(|e| graph_kmer::canonical(e.code, k))
                .collect()
        );
        let balanced = counts(&[(&major, 100), (&minor, 100)], k);
        let (mut edges, _, _) = graph(&balanced, k, 3);
        let (_, branches, _, _) = simplify(&mut edges, k, 0.2, 0.49);
        assert!(branches > 0);
        assert!(edges
            .iter()
            .any(|e| !expected.contains_key(&graph_kmer::canonical(e.code, k))));
    }

    /// The prior path-cloning search, kept as a differential oracle.
    fn stronger_bubble_path_reference(
        weak_id: usize,
        unitigs: &[Unitig],
        starts: &graph_kmer::Map<Vec<usize>>,
        removed: &[bool],
        k: usize,
        fraction: f64,
    ) -> Option<bubbles::Alternative> {
        let weak = &unitigs[weak_id];
        let max_length = 3 * k;
        let max_edits = 3.max((weak.sequence.len() - k + 1) / 10);
        let reverse_from = graph_kmer::reverse_code(weak.to, k - 1);
        let reverse_to = graph_kmer::reverse_code(weak.from, k - 1);
        // A bounded simple-path search on unitigs, rather than one visit per k-mer
        // node: different coverage histories reaching a node must remain distinct.
        let mut stack = vec![(weak.from, Vec::<usize>::new(), 0usize, 0.0)];
        let mut states = 0;
        while let Some((node, path, length, mass)) = stack.pop() {
            for &id in starts.get(&node).into_iter().flatten() {
                states += 1;
                if states > 4096 {
                    return None;
                }
                let u = &unitigs[id];
                let n = u.sequence.len() - k + 1;
                if id == weak_id
                    || u.from == u.to
                    || u.circular
                    || u.to == weak.from
                    || (u.from == reverse_from && u.to == reverse_to)
                    || length + n > max_length
                    || path.iter().any(|&p| {
                        let old = &unitigs[p];
                        old.from == u.to
                            || (u.from == graph_kmer::reverse_code(old.to, k - 1)
                                && u.to == graph_kmer::reverse_code(old.from, k - 1))
                    })
                {
                    continue;
                }
                if removed[id] {
                    continue;
                }
                let new_mass = mass + u.mean_count * n as f64;
                let mut route = path.clone();
                route.push(id);
                if u.to == weak.to {
                    if weak.mean_count <= fraction * new_mass / (length + n) as f64 {
                        let mut sequence = graph_kmer::decode(weak.from, k - 1);
                        for &p in &route {
                            sequence.extend_from_slice(&unitigs[p].sequence[k - 1..]);
                        }
                        if similar_paths(&weak.sequence, &sequence, max_edits) {
                            let minimum =
                                route.iter().map(|&p| unitigs[p].min_count).min().unwrap();
                            return Some(bubbles::Alternative {
                                sequence,
                                mean_count: new_mass / (length + n) as f64,
                                min_count: minimum,
                            });
                        }
                    }
                } else {
                    stack.push((u.to, route, length + n, new_mass));
                }
            }
        }
        None
    }

    #[test]
    fn marked_bubble_search_matches_the_path_cloning_search() {
        let mut state = 0x2545f4914f6cdd1du64;
        let mut next = |n: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % n
        };
        let (mut found, mut capped, mut exhausted) = (0, 0, 0);
        for (round, k) in [11usize, 13, 15, 21, 11, 13, 11, 15, 21, 15]
            .into_iter()
            .enumerate()
        {
            let genome = dna(220 + 40 * round);
            let mut sequences = vec![(genome.clone(), 60u32)];
            // Random concatenations of a few short units: k-mers recur in many
            // contexts, so walks branch combinatorially and reach the state cap.
            if round >= 8 {
                let units: [&[u8]; 4] = [b"ACGTACG", b"GGTACCA", b"TTGACAT", b"CAGTTGC"];
                for _ in 0..80 {
                    let read: Vec<u8> = (0..20)
                        .flat_map(|_| units[next(4) as usize].iter().copied())
                        .collect();
                    sequences.push((read, 3 + next(20) as u32));
                }
            }
            // Many supported variants in one short window cross into a dense tangle.
            let dense = round >= 6;
            for _ in 0..(12 + 6 * round) {
                let mut copy = genome.clone();
                if dense {
                    for _ in 0..4 + next(8) {
                        copy[100 + next(50) as usize] = b"ACGT"[next(4) as usize];
                    }
                    sequences.push((copy, 4 + next(30) as u32));
                    continue;
                }
                // Error copies and short tandem duplications create bubbles, tips and tangles.
                for _ in 0..1 + next(3) {
                    let at = next(copy.len() as u64) as usize;
                    match next(4) {
                        0 => copy[at] = b"ACGT"[next(4) as usize],
                        1 => {
                            copy.remove(at);
                        }
                        2 => {
                            let end = (at + 5 + next(20) as usize).min(copy.len());
                            let unit = copy[at..end].to_vec();
                            copy.splice(at..at, unit);
                        }
                        _ => copy.insert(at, b"ACGT"[next(4) as usize]),
                    }
                }
                sequences.push((copy, 3 + next(40) as u32));
            }
            let pairs: Vec<(&[u8], u32)> =
                sequences.iter().map(|(s, n)| (s.as_slice(), *n)).collect();
            let (_, unitigs, _) = graph(&counts(&pairs, k), k, 3);
            let mut starts = graph_kmer::Map::<Vec<_>>::default();
            for (id, u) in unitigs.iter().enumerate() {
                starts.entry(u.from).or_default().push(id);
            }
            let mut search = BubbleSearch::new(&unitigs, &starts, k);
            for fraction in [0.2, 0.49, 0.9] {
                let removed: Vec<bool> = (0..unitigs.len()).map(|_| next(5) == 0).collect();
                for weak in 0..unitigs.len() {
                    if unitigs[weak].circular || unitigs[weak].from == unitigs[weak].to {
                        continue;
                    }
                    let expected = stronger_bubble_path_reference(
                        weak, &unitigs, &starts, &removed, k, fraction,
                    );
                    let actual = search.stronger_path(weak, &removed, fraction);
                    assert!(search.walker.marks.iter().all(|&m| m == 0));
                    match (&expected, &actual) {
                        (None, None) if search.walker.states > 4096 => capped += 1,
                        (None, None) => exhausted += 1,
                        (Some(a), Some(b)) => {
                            assert_eq!(a.sequence, b.sequence);
                            assert_eq!(a.mean_count.to_bits(), b.mean_count.to_bits());
                            assert_eq!(a.min_count, b.min_count);
                            found += 1;
                        }
                        _ => panic!("k={k} weak={weak} fraction={fraction}: results differ"),
                    }
                }
            }
        }
        assert!(
            found > 0 && capped > 0 && exhausted > 0,
            "{found} {capped} {exhausted}"
        );
    }

    #[test]
    fn balanced_haplotypes_leave_counted_allele_bubbles_but_errors_do_not() {
        let k = 21;
        let major = dna(600);
        let mut minor = major.clone();
        for at in [100, 250, 400] {
            minor[at] = if minor[at] == b'A' { b'C' } else { b'A' };
        }
        let cleaned = |minor_count: u32| {
            let (mut edges, _, _) = graph(&counts(&[(&major, 60), (&minor, minor_count)], k), k, 3);
            simplify(&mut edges, k, 0.2, 0.49).0
        };
        let unitigs = cleaned(40);
        let minors = balanced_bubbles(&unitigs, k);
        assert_eq!(minors.len(), 3, "{minors:?}");
        assert!(minors.iter().all(|&m| (m - 0.4).abs() < 0.01), "{minors:?}");
        assert_eq!(
            mixed_sample(&unitigs, k, &[])["status"],
            "mixed_sample_signal"
        );
        // A weak minor sequence is cleaned like sequencing error and not counted.
        let unitigs = cleaned(5);
        assert!(balanced_bubbles(&unitigs, k).is_empty());
        assert_eq!(mixed_sample(&unitigs, k, &[])["status"], "not_detected");
    }

    #[test]
    fn high_coverage_self_loop_cannot_inflate_a_bubble_alternative() {
        let k = 3;
        let make = |s: &[u8], mean_count: f64| {
            let mut codes = Vec::new();
            graph_kmer::visit(s, None, k, 0, |_, c| codes.push(c));
            Unitig {
                from: codes[0].prefix(),
                to: codes.last().unwrap().suffix(k - 1),
                sequence: s.to_vec(),
                circular: false,
                mean_count,
                max_count: mean_count as u32,
                min_count: mean_count as u32,
                first_count: mean_count as u32,
                last_count: mean_count as u32,
            }
        };
        let unitigs = vec![
            make(b"AAGCC", 10.0),
            make(b"AAT", 10.0),
            make(b"ATAT", 1000.0),
            make(b"ATCC", 10.0),
        ];
        let mut starts = graph_kmer::Map::<Vec<_>>::default();
        for (id, u) in unitigs.iter().enumerate() {
            starts.entry(u.from).or_default().push(id);
        }
        assert!(BubbleSearch::new(&unitigs, &starts, k)
            .stronger_path(0, &vec![false; unitigs.len()], 0.2)
            .is_none());
    }

    #[test]
    fn wide_bubble_cleanup_does_not_choose_large_deletions() {
        let major = dna(650);
        let short = [&major[..200], &major[300..]].concat();
        let input = counts(&[(&major, 100), (&short, 15)], 63);
        let (mut edges, _, _) = graph(&input, 63, 3);
        let (_, _, stats, _) = simplify(&mut edges, 63, 0.2, 0.49);
        assert_eq!(stats["bubble_kmers_removed"], 0);
    }
    #[test]
    fn twin_flags_match_deleted_edge_membership_for_linear_and_circular_unitigs() {
        for k in [31, 63, 95, 127] {
            let major = dna(700);
            let mut alternate = major.clone();
            alternate[310] = if alternate[310] == b'A' { b'C' } else { b'A' };
            let circle = &major[100..180];
            let circular: Vec<_> = circle
                .iter()
                .cycle()
                .take(circle.len() + k - 1)
                .copied()
                .collect();
            let (_, unitigs, _) = graph(
                &counts(&[(&major, 100), (&alternate, 15), (&circular, 5)], k),
                k,
                3,
            );
            let twins = unitig_twins(&unitigs, k);
            for (removed_id, unitig) in unitigs.iter().enumerate() {
                let mut removed = HashSet::new();
                graph_kmer::visit(&unitig.sequence, None, k, 0, |_, code| {
                    removed.insert(code);
                    removed.insert(graph_kmer::reverse_code(code, k));
                });
                for (id, other) in unitigs.iter().enumerate() {
                    let mut unavailable = false;
                    graph_kmer::visit(&other.sequence, None, k, 0, |_, code| {
                        unavailable |= removed.contains(&code)
                    });
                    assert_eq!(
                        unavailable,
                        id == removed_id || id == twins[removed_id],
                        "k={k}, removed={removed_id}, other={id}"
                    );
                }
            }
        }
    }
    #[test]
    fn structural_link_protection_requires_distinct_lengths_and_shared_full_flanks() {
        let major = dna(700);
        for (deleted, expected) in [(0, false), (2, false), (100, true), (180, true)] {
            let mut minor = [&major[..220], &major[220 + deleted..]].concat();
            if deleted == 0 {
                minor[220] = if minor[220] == b'A' { b'C' } else { b'A' };
            }
            let (_, unitigs, _) = graph(&counts(&[(&major, 100), (&minor, 10)], 63), 63, 3);
            let mut starts = graph_kmer::Map::<Vec<usize>>::default();
            for (id, unitig) in unitigs.iter().enumerate() {
                starts.entry(unitig.from).or_default().push(id);
            }
            let protected = unitigs
                .iter()
                .filter(|u| has_structural_alternative(u, &unitigs, &starts, 0.2))
                .count();
            assert_eq!(protected, if expected { 2 } else { 0 }, "deleted={deleted}");
        }
    }
}
