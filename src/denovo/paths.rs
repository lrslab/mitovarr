//! Joint fragment scoring of complete paths between sample-graph anchors.
use super::{
    align::{Index, Placements, Prepared},
    normalize,
    pairs::{self, Library},
    threading::PathGraph,
    Candidate,
};
use crate::{
    cli::{GraphSettings, Limits, Reads},
    io::fastq,
    kmer,
};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::{
    cmp::Reverse,
    collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet},
    path::Path,
};

mod circles;
mod ends;
pub(super) mod primary;
mod repeats;
mod seeds;
#[cfg(test)]
mod terminal_tests;

const MAX_PATHS: usize = 128;
const MAX_STATES: usize = 50_000;
const FLANK: usize = 20;

// Accumulated intermediate bases, complete node history, then next node.
// Keeping the history is essential: equal-cost arrivals at one node can spell
// different alleles or repeat copy counts. Ordering only schedules discovery;
// it supplies no length preference to the fragment scorer.
type PendingWalk = Reverse<(usize, Vec<usize>, usize)>;

#[derive(Default)]
struct Frontier {
    pending: BinaryHeap<PendingWalk>,
    peak: usize,
    limited: bool,
}

impl Frontier {
    fn push(&mut self, bases: usize, path: Vec<usize>, node: usize) -> bool {
        if self.pending.len() >= MAX_STATES {
            self.limited = true;
            return false;
        }
        self.pending.push(Reverse((bases, path, node)));
        self.peak = self.peak.max(self.pending.len());
        true
    }
}

struct Route {
    circle_length: Option<usize>,
    to: usize,
    nodes: Vec<usize>,
    sequence: Vec<u8>,
    left_end: usize,
    right_start: usize,
    // A pure SNP bubble contains shared k-1 overlaps on both sides. Those
    // shared bases are not the allele whose spanning support is being counted.
    snp_interval: Option<(usize, usize)>,
    // Additional allele support is for the partial primary only. Strict joins
    // retain their original graph-corridor / single-SNP spanning requirement.
    primary_interval: Option<(usize, usize)>,
}
impl Route {
    fn spanning_interval(&self) -> (usize, usize) {
        self.snp_interval
            .unwrap_or((self.left_end, self.right_start))
    }

    fn primary_spanning_interval(&self) -> (usize, usize) {
        self.primary_interval.unwrap_or(self.spanning_interval())
    }
}

/// Locate every difference within one shared source/target corridor, including
/// the full range of equivalent indel placements in a homopolymer. This never
/// changes the interval used to grant a strict join or complete publication.
fn identify_primary_allele(from: usize, routes: &mut [Route]) {
    let Some(first) = routes.first() else { return };
    if routes.len() < 2
        || first.to == from
        || routes
            .iter()
            .any(|r| r.to != first.to || r.circle_length.is_some())
        || routes.iter().all(|r| r.sequence == first.sequence)
    {
        return;
    }
    let prefix = routes
        .iter()
        .map(|r| {
            r.sequence
                .iter()
                .zip(&first.sequence)
                .take_while(|(a, b)| a == b)
                .count()
        })
        .min()
        .unwrap();
    let suffix = routes
        .iter()
        .map(|r| {
            r.sequence
                .iter()
                .rev()
                .zip(first.sequence.iter().rev())
                .take_while(|(a, b)| a == b)
                .count()
        })
        .min()
        .unwrap();
    let shortest = routes.iter().map(|r| r.sequence.len()).min().unwrap();
    let start = prefix.min(shortest - suffix);
    let ambiguous_overlap = prefix.saturating_sub(shortest - suffix);
    // Shared flanks must stay outside the differing block on every route.
    // Different exits and differences inside anchors are not local alleles.
    if routes.iter().any(|r| {
        let end = r.sequence.len() - suffix + ambiguous_overlap;
        start < r.left_end || end > r.right_start || start < FLANK || end + FLANK > r.sequence.len()
    }) {
        return;
    }
    for r in routes {
        r.primary_interval = Some((start, r.sequence.len() - suffix + ambiguous_overlap));
    }
}

fn identify_single_snp(from: usize, routes: &mut [Route]) {
    if routes.len() < 2 {
        return;
    }
    let first = &routes[0];
    if first.to == from
        || routes.iter().any(|r| {
            r.to != first.to
                || r.sequence.len() != first.sequence.len()
                || r.circle_length.is_some()
        })
    {
        return;
    }
    let Some(position) = (0..first.sequence.len()).find(|&p| {
        routes
            .iter()
            .skip(1)
            .any(|r| r.sequence[p] != first.sequence[p])
    }) else {
        return;
    };
    if routes.iter().any(|r| {
        r.sequence[position + 1..] != first.sequence[position + 1..]
            || position < r.left_end
            || position >= r.right_start
    }) {
        return;
    }
    for route in routes {
        route.snp_interval = Some((position, position + 1));
    }
}
struct Locus {
    from: usize,
    routes: Vec<Route>,
    indexes: Vec<Index>,
    states: usize,
    exhausted: bool,
    state_limited: bool,
    path_limited: bool,
    cycle_expansion_limited: bool,
    frontier_peak: usize,
    frontier_limited: bool,
    length_limited: bool,
    repeat_cycles: Vec<repeats::Cycle>,
    dead_ends: usize,
    fragments: usize,
    unplaced: usize,
    r1_only: usize,
    r2_only: usize,
    r2_anchored: usize,
    alignment_limited: usize,
    totals: Vec<f64>,
    offset: f64,
    delta_squares: Vec<Vec<f64>>,
    strong_votes: Vec<Vec<u32>>,
    spanning: Vec<u32>,
    primary_spanning: Vec<u32>,
    primary_strong_votes: Vec<Vec<u32>>,
    selected: Option<usize>,
}

fn enumerate(
    graph: &PathGraph,
    from: usize,
    anchors: &HashSet<usize>,
    overlap: usize,
    max_bases: usize,
    context: usize,
    max_read: usize,
) -> Locus {
    let left = &graph.segments[from].sequence;
    let left_start = left.len().saturating_sub(context);
    let left = &left[left_start..];
    let mut routes = Vec::new();
    let mut states = 0;
    let mut exhausted = false;
    let mut state_limited = false;
    let mut path_limited = false;
    let mut cycle_expansion_limited = false;
    let mut length_limited = false;
    let mut dead_ends = 0;
    let cycles = repeats::cycles(graph, anchors, overlap);
    let by_node: HashMap<_, _> = cycles
        .iter()
        .enumerate()
        .flat_map(|(i, c)| c.nodes.iter().map(move |&node| (node, i)))
        .collect();
    let mut used_cycles = BTreeSet::new();
    let mut frontier = Frontier::default();
    for &to in &graph.outgoing[from] {
        if !frontier.push(0, Vec::new(), to) {
            exhausted = true;
            break;
        }
    }
    let mut unique = BTreeSet::new();
    'search: while let Some(Reverse((bases, path, node))) = frontier.pending.pop() {
        states += 1;
        if states > MAX_STATES {
            exhausted = true;
            state_limited = true;
            break;
        }
        if anchors.contains(&node) {
            let mut sequence = left.to_vec();
            for &id in &path {
                sequence.extend_from_slice(&graph.segments[id].sequence[overlap..]);
            }
            let right_start = sequence.len();
            let target = &graph.segments[node].sequence;
            sequence.extend_from_slice(&target[overlap..target.len().min(context)]);
            if unique.insert((node, sequence.clone())) {
                routes.push(Route {
                    circle_length: None,
                    to: node,
                    nodes: path,
                    sequence,
                    left_end: left.len() - overlap,
                    right_start,
                    snp_interval: None,
                    primary_interval: None,
                });
                if routes.len() >= MAX_PATHS && !frontier.pending.is_empty() {
                    exhausted = true;
                    path_limited = true;
                    break;
                }
            }
            continue;
        }
        if let Some(&ci) = by_node.get(&node) {
            let cycle = &cycles[ci];
            used_cycles.insert(ci);
            // Every additional turn is an explicit integer parameter. Generate
            // exits around one turn for each admissible n instead of revisiting
            // cycle edges recursively. The higher-copy tail remains unknown.
            length_limited = true;
            let entry = cycle.nodes.iter().position(|&id| id == node).unwrap();
            let mut extension = Vec::new();
            let mut added = 0;
            let maximum_turns = max_bases.saturating_sub(bases) / cycle.period;
            let mut expanded = Vec::new();
            'turns: for _ in 0..=maximum_turns {
                for offset in 0..cycle.nodes.len() {
                    let id = cycle.nodes[(entry + offset) % cycle.nodes.len()];
                    added += graph.segments[id].sequence.len() - overlap;
                    if bases + added > max_bases {
                        break 'turns;
                    }
                    extension.push(id);
                    for &next in &graph.outgoing[id] {
                        if by_node.get(&next) == Some(&ci) {
                            continue;
                        }
                        let mut nodes = path.clone();
                        nodes.extend_from_slice(&extension);
                        expanded.push((next, nodes, bases + added));
                        if expanded.len() >= MAX_PATHS {
                            exhausted = true;
                            cycle_expansion_limited = true;
                            break 'turns;
                        }
                    }
                }
            }
            for (next, nodes, bases) in expanded {
                if !frontier.push(bases, nodes, next) {
                    exhausted = true;
                    break 'search;
                }
            }
            continue;
        }
        let added = graph.segments[node].sequence.len().saturating_sub(overlap);
        if bases + added > max_bases {
            length_limited = true;
            continue;
        }
        let mut path = path;
        path.push(node);
        if graph.outgoing[node].is_empty() {
            dead_ends += 1;
        }
        for &next in &graph.outgoing[node] {
            if !frontier.push(bases + added, path.clone(), next) {
                exhausted = true;
                break 'search;
            }
        }
    }
    identify_single_snp(from, &mut routes);
    identify_primary_allele(from, &mut routes);
    let n = routes.len();
    let indexes = routes
        .iter()
        .map(|r| Index::new(&r.sequence, false, max_read))
        .collect();
    Locus {
        from,
        routes,
        indexes,
        states,
        exhausted,
        state_limited,
        path_limited,
        cycle_expansion_limited,
        frontier_peak: frontier.peak,
        frontier_limited: frontier.limited,
        length_limited,
        repeat_cycles: used_cycles.into_iter().map(|i| cycles[i].clone()).collect(),
        dead_ends,
        fragments: 0,
        unplaced: 0,
        r1_only: 0,
        r2_only: 0,
        r2_anchored: 0,
        alignment_limited: 0,
        totals: vec![0.0; n],
        offset: 0.0,
        delta_squares: vec![vec![0.0; n]; n],
        strong_votes: vec![vec![0; n]; n],
        spanning: vec![0; n],
        primary_spanning: vec![0; n],
        primary_strong_votes: vec![vec![0; n]; n],
        selected: None,
    }
}

fn joint(
    a: &Placements,
    b: &Placements,
    route: &Route,
    model: Option<&Library>,
    baseline: f64,
    pivot: (isize, bool, usize),
) -> (f64, bool, bool) {
    if a.hits.is_empty() || b.hits.is_empty() {
        return (-12.0, false, false);
    }
    // A geometry state has a prior fixed by the library, not by how many hits
    // this candidate happens to provide. Duplicate alignments within a state
    // contribute only their best full-fragment sequence explanation.
    let mut states = BTreeMap::<usize, f64>::new();
    let mut strongest = 0.0f64;
    let mut bridge = 0.0f64;
    let mut primary_bridge = 0.0f64;
    for x in &a.hits {
        if x.reverse != pivot.1 || (x.start as isize).abs_diff(pivot.0) > pivot.2 {
            continue;
        }
        for y in &b.hits {
            if route
                .circle_length
                .is_some_and(|n| x.end.max(y.end) - x.start.min(y.start) > n)
            {
                continue;
            }
            let span = match model {
                Some(m) => match m.span(x, y) {
                    Some(d) => d,
                    None => continue,
                },
                None => 0,
            };
            let empirical = model.is_some_and(|m| m.geometry_usable);
            let weight = model
                .filter(|m| m.geometry_usable)
                .map_or(1.0, |m| m.probability(span));
            let sequence = (x.logp + y.logp - baseline).exp().min(1.0);
            states
                .entry(if empirical { span } else { 0 })
                .and_modify(|e| *e = e.max(sequence))
                .or_insert(sequence);
            let relative = sequence
                * model
                    .filter(|m| m.geometry_usable)
                    .map_or(1.0, |m| weight / m.peak_probability());
            strongest = strongest.max(relative);
            let crosses = |(left_end, right_start)| {
                let single = |h: &super::align::Hit| {
                    h.start + FLANK <= left_end && h.end >= right_start + FLANK
                };
                single(x)
                    || single(y)
                    || (model.is_some_and(|m| m.geometry_usable && span <= m.upper)
                        && x.start.min(y.start) + FLANK <= left_end
                        && x.end.max(y.end) >= right_start + FLANK)
            };
            if crosses(route.spanning_interval()) {
                bridge = bridge.max(relative);
            }
            if crosses(route.primary_spanning_interval()) {
                primary_bridge = primary_bridge.max(relative);
            }
        }
    }
    let likelihood: f64 = states
        .into_iter()
        .map(|(d, e)| {
            model
                .filter(|m| m.geometry_usable)
                .map_or(1.0, |m| m.probability(d))
                * e
        })
        .sum();
    (
        (likelihood.min(1.0) + (-12.0f64).exp()).ln(),
        bridge >= 0.1 && bridge >= strongest * 0.99,
        primary_bridge >= 0.1 && primary_bridge >= strongest * 0.99,
    )
}

fn single(
    a: &Placements,
    route: &Route,
    baseline: f64,
    pivot: (isize, bool, usize),
) -> (f64, bool, bool) {
    let hits: Vec<_> = a
        .hits
        .iter()
        .filter(|x| x.reverse == pivot.1 && (x.start as isize).abs_diff(pivot.0) <= pivot.2)
        .collect();
    let best = hits
        .iter()
        .map(|h| h.logp)
        .fold(f64::NEG_INFINITY, f64::max);
    let bridge = |(left_end, right_start)| {
        hits.iter().any(|h| {
            h.logp >= best - 2.3 && h.start + FLANK <= left_end && h.end >= right_start + FLANK
        })
    };
    (
        ((best - baseline).exp().min(1.0) + (-12.0f64).exp()).ln(),
        bridge(route.spanning_interval()),
        bridge(route.primary_spanning_interval()),
    )
}

fn score(
    a: &Placements,
    b: &Placements,
    route: &Route,
    model: Option<&Library>,
    baseline: f64,
    pivot: (isize, bool, usize),
    r1_only: bool,
) -> (f64, bool, bool) {
    if r1_only {
        single(a, route, baseline, pivot)
    } else {
        joint(a, b, route, model, baseline, pivot)
    }
}

fn unscored_contrast_bound(unknown: usize) -> f64 {
    // Every fragment score is in [-12, log(1 + exp(-12))]. Uniform
    // omission makes pools comparable but does not make its contrast zero.
    unknown as f64 * (12.0 + (-12.0f64).exp().ln_1p())
}

fn stable_margin(difference: f64, square_sum: f64, unknown: usize) -> bool {
    let bound = unscored_contrast_bound(unknown);
    let lower = difference - bound;
    // Triangle inequality also covers unknown fragments sharing an existing
    // start cluster. This is a support criterion, not a confidence interval.
    lower >= 10.0 && lower >= 3.0 * (square_sum.sqrt() + bound)
}

fn minor_supported(major: u32, minor: u32, unknown: usize, settings: &GraphSettings) -> bool {
    let upper = f64::from(minor) + unknown as f64;
    upper < f64::from(settings.min_junction_fragments)
        || upper <= f64::from(major) * settings.max_error_fraction
}

fn rescue_mates(
    indexes: &[Index],
    placements: &mut [(Placements, Placements)],
    prepared: &mut Prepared,
    pivots: &[(isize, bool, usize)],
    model: Option<&Library>,
    read_length: usize,
) {
    // If every candidate is equally unplaced, no candidate receives an
    // exclusive penalty and there is no decision-dependent rescue to perform.
    if !placements.iter().any(|(_, b)| !b.hits.is_empty()) {
        return;
    }
    let nominated: Vec<_> = placements
        .iter()
        .enumerate()
        .flat_map(|(i, (_, b))| b.hits.iter().map(move |h| (i, h.clone())))
        .collect();
    let longest = indexes.iter().map(|i| i.sequence.len()).max().unwrap();
    let shortest = indexes.iter().map(|i| i.sequence.len()).min().unwrap();
    let displacement = longest - shortest + (read_length / 10).max(3) + 12;
    for (i, (a, b)) in placements.iter_mut().enumerate() {
        if !a.hits.is_empty() && !b.exhausted {
            let pivot = pivots[i];
            let start = pivot.0.max(0) as usize;
            let allowed = if let Some(m) = model.filter(|m| m.geometry_usable) {
                (start.saturating_sub(m.upper), start + read_length + m.upper)
            } else {
                (0, indexes[i].sequence.len())
            };
            if b.hits.is_empty() {
                *b = indexes[i].rescue_near(prepared, allowed.0, allowed.1, !pivot.1);
                continue;
            }
            // An existing seeded decoy must not suppress a better seedless
            // explanation. Project better mate placements from the shared
            // candidate pool relative to its fixed anchor-mate pivots. Sequence length
            // differences and the edit allowance pad those nominated windows.
            let best = b
                .hits
                .iter()
                .map(|h| h.logp)
                .fold(f64::NEG_INFINITY, f64::max);
            let mut windows: Vec<_> = nominated
                .iter()
                .filter(|(_, h)| h.reverse != pivot.1 && h.logp > best + 1e-9)
                .filter_map(|(j, h)| {
                    let p = h.start as isize - pivots[*j].0 + pivot.0;
                    let left = p.saturating_sub(displacement as isize).max(0) as usize;
                    let right = (p + (h.end - h.start + displacement) as isize).max(0) as usize;
                    let bounds = (left.max(allowed.0), right.min(allowed.1));
                    (bounds.0 < bounds.1).then_some(bounds)
                })
                .collect();
            windows.sort_unstable();
            let mut merged = Vec::<(usize, usize)>::new();
            for (left, right) in windows {
                if let Some(last) = merged.last_mut() {
                    if left <= last.1 {
                        last.1 = last.1.max(right);
                        continue;
                    }
                }
                merged.push((left, right));
            }
            for (left, right) in merged {
                let rescued = indexes[i].rescue_near(prepared, left, right, !pivot.1);
                b.exhausted |= rescued.exhausted;
                b.hits.extend(rescued.hits);
            }
            b.hits.sort_by(|a, b| {
                (a.start, a.reverse)
                    .cmp(&(b.start, b.reverse))
                    .then(b.logp.total_cmp(&a.logp))
            });
            b.hits
                .dedup_by(|b, a| a.start == b.start && a.reverse == b.reverse);
            b.exhausted |= b.hits.len() > 512;
        }
    }
}

pub(super) struct Resolution {
    pub candidates: Vec<Candidate>,
    pub primary_candidates: Vec<primary::Primary>,
    pub supported_cycles: BTreeSet<Vec<u8>>,
    pub report: Value,
}

fn path_anchors(graph: &PathGraph, k: usize, seeds: &seeds::Seeds) -> HashSet<usize> {
    let mut has_incoming = vec![false; graph.segments.len()];
    for targets in &graph.outgoing {
        for &target in targets {
            has_incoming[target] = true;
        }
    }
    graph
        .segments
        .iter()
        .enumerate()
        .filter_map(|(id, s)| {
            // A real gene tail can end in a short sink behind a graph merge.
            // It needs unique sequence and room beyond the shared overlap, not
            // the long internal-anchor cutoff. The reverse source is eligible
            // too; actual joins still require full-fragment spanning evidence.
            let terminal = graph.outgoing[id].is_empty() || !has_incoming[id];
            let isolated = graph.outgoing[id].is_empty() && !has_incoming[id];
            let minimum_length = if terminal && !isolated {
                k - 1 + FLANK
            } else {
                4 * k
            };
            if s.circular || s.sequence.len() < minimum_length {
                return None;
            }
            seeds.anchor_length(&s.sequence, k).map(|_| id)
        })
        .collect()
}

pub(super) fn resolve(
    graph: &PathGraph,
    reads: &Reads,
    out: &Path,
    settings: &GraphSettings,
    limits: &Limits,
) -> Result<Resolution> {
    let training = pairs::estimate(graph, reads, settings, limits)?;
    let overlap = settings.k - 1;
    let context = training
        .model
        .as_ref()
        .map_or(4 * training.max_read, |m| m.upper + 2 * training.max_read)
        .max(4 * settings.k);
    let max_bases = training
        .model
        .as_ref()
        .map_or(4 * training.max_read, |m| 2 * m.upper)
        .max(4 * settings.k);
    let original_segments = graph.segments.len();
    let (graph, end_assembly) = ends::augment(graph, reads, out, settings, limits, max_bases)?;
    let graph = &graph;
    let seeds = seeds::Seeds::new(graph, settings.k.clamp(31, 63));
    let occurrence = seeds.occurrences(31);
    let mut anchors = path_anchors(graph, settings.k, &seeds);
    // Newly proposed connectors are alternatives between existing anchors,
    // never new independent anchors that could certify their own joins.
    anchors.retain(|&id| id < original_segments);
    let short_terminal_anchors = anchors
        .iter()
        .filter(|&&id| graph.segments[id].sequence.len() < 4 * settings.k)
        .count();
    let mut anchor_order: Vec<_> = anchors.iter().copied().collect();
    anchor_order.sort_unstable();
    let mut loci: Vec<_> = anchor_order
        .iter()
        .map(|&from| {
            enumerate(
                graph,
                from,
                &anchors,
                overlap,
                max_bases,
                context,
                training.max_read,
            )
        })
        .collect();
    tracing::info!(
        anchors = anchors.len(),
        local_paths = loci.iter().map(|l| l.routes.len()).sum::<usize>(),
        "scoring complete local paths against paired reads"
    );
    let records: Vec<_> = loci.iter().enumerate().flat_map(|(i,locus)| locus.routes.iter().enumerate().map(move |(j,route)|
        (format!("locus_{}_path_{} source_segment={} target_segment={} left_anchor_end={} right_anchor_start={}",i+1,j+1,locus.from+1,route.to+1,route.left_end,route.right_start),route.sequence.clone())
    )).collect();
    crate::io::fasta::write_fasta_records(&out.join("path_candidates.fasta"), &records)?;
    // Recruit a COMMON fragment pool from each source anchor's terminal region.
    // The pool is independent of which alternative is currently being scored.
    let mut lookup = BTreeMap::<usize, HashMap<_, Vec<(usize, isize, kmer::Code)>>>::new();
    let mut anchor_contexts = Vec::new();
    for (id, locus) in loci
        .iter()
        .enumerate()
        .filter(|(_, l)| !l.routes.is_empty() && !l.exhausted && !l.length_limited)
    {
        let sequence = &graph.segments[locus.from].sequence;
        let recruit = training
            .model
            .as_ref()
            .map_or(2 * training.max_read, |m| m.upper)
            .max(settings.k);
        let start = sequence.len().saturating_sub(recruit);
        let context_start = sequence.len().saturating_sub(context);
        let recruitment = seeds.recruitment(&sequence[start..], settings.k);
        anchor_contexts.push(json!({"source_segment":locus.from+1,
            "anchor_seed_bases":seeds.anchor_length(sequence,settings.k),
            "recruitment_seeds":recruitment.iter().map(|(k,p)|json!({"k":k,"unique_positions":p.len()})).collect::<Vec<_>>()}));
        for (k, positions) in recruitment {
            let index = lookup.entry(k).or_default();
            for (offset, code) in positions {
                index.entry(kmer::canonical(code, k)).or_default().push((
                    id,
                    (start + offset - context_start) as isize,
                    code,
                ));
            }
        }
    }
    let mut names = HashSet::new();
    let mut repeat_depth = repeats::Depth::new(&loci, graph);
    let scan = fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
        repeat_depth.observe_r1(a, limits.min_baseq);
        // Either mate may qualify through a shared unique anchor. Prefer R1
        // when both qualify, BEFORE any candidate alignment or geometry score.
        // The other mate need not map, so inclusion cannot favour a repeat length.
        let mut relevant = BTreeMap::<usize, [BTreeSet<(isize, bool)>; 2]>::new();
        for (mate, read) in [a, b].into_iter().enumerate() {
            for (&seed_k, index) in &lookup {
                kmer::visit(
                    &read.sequence,
                    Some(&read.quality),
                    seed_k,
                    limits.min_baseq.min(10),
                    |offset, code| {
                        if let Some(anchors) = index.get(&kmer::canonical(code, seed_k)) {
                            for &(id, position, forward) in anchors {
                                let reverse = code != forward;
                                let offset = if reverse {
                                    read.sequence.len() - seed_k - offset
                                } else {
                                    offset
                                };
                                relevant.entry(id).or_default()[mate]
                                    .insert((position - offset as isize, reverse));
                            }
                        }
                    },
                );
            }
        }
        if relevant.is_empty() {
            return Ok(());
        }
        let name = a.name(1)?;
        if names.contains(name) {
            return Ok(());
        }
        ensure!(
            names.len() < settings.max_validation_pairs,
            "path scoring exceeds --max-validation-pairs"
        );
        names.insert(name.to_owned());
        let mut prepared_a = Prepared::new(a);
        let mut prepared_b = Prepared::new(b);
        for (id, anchors) in relevant {
            let mate = usize::from(anchors[0].is_empty());
            let pivots = &anchors[mate];
            // Library::span uses genomic forward/reverse geometry, not mate
            // labels; swapping the working pair preserves FR/RF and distance.
            let (a, b, prepared_a, prepared_b) = if mate == 0 {
                (a, b, &mut prepared_a, &mut prepared_b)
            } else {
                (b, a, &mut prepared_b, &mut prepared_a)
            };
            let locus = &mut loci[id];
            locus.fragments += 1;
            locus.r2_anchored += mate;
            let &(lo, reverse) = pivots.first().unwrap();
            let &(hi, _) = pivots.last().unwrap();
            let tolerance = (a.sequence.len() / 10).max(3);
            if pivots.iter().any(|p| p.1 != reverse) || hi.abs_diff(lo) > 2 * tolerance {
                locus.alignment_limited += 1;
                continue;
            }
            let pivot = ((lo + hi) / 2, reverse, tolerance);
            let mut placements: Vec<_> = locus
                .indexes
                .iter()
                .map(|index| {
                    (
                        index.align_prepared(prepared_a),
                        index.align_prepared(prepared_b),
                    )
                })
                .collect();
            rescue_mates(
                &locus.indexes,
                &mut placements,
                prepared_b,
                &vec![pivot; locus.routes.len()],
                training.model.as_ref(),
                a.sequence.len(),
            );
            if placements.iter().any(|(a, b)| a.exhausted || b.exhausted) {
                // Omit the same molecule from EVERY hypothesis if its alignment
                // search is capped; it cannot cast an exclusive vote for any path.
                locus.alignment_limited += 1;
                continue;
            }
            let best_a = placements
                .iter()
                .flat_map(|(a, _)| a.hits.iter().map(|h| h.logp))
                .fold(f64::NEG_INFINITY, f64::max);
            let best_b = placements
                .iter()
                .flat_map(|(_, b)| b.hits.iter().map(|h| h.logp))
                .fold(f64::NEG_INFINITY, f64::max);
            let r1_only = !best_b.is_finite();
            locus.r1_only += usize::from(r1_only && mate == 0);
            locus.r2_only += usize::from(r1_only && mate == 1);
            if !best_a.is_finite() {
                locus.unplaced += 1;
            }
            let baseline = super::align::ideal_score(a)
                + if r1_only {
                    0.0
                } else {
                    super::align::ideal_score(b)
                };
            let scores: Vec<_> = placements
                .iter()
                .zip(&locus.routes)
                .map(|((a, b), route)| {
                    score(
                        a,
                        b,
                        route,
                        training.model.as_ref(),
                        baseline,
                        pivot,
                        r1_only,
                    )
                })
                .collect();
            let maximum = scores.iter().map(|s| s.0).fold(f64::NEG_INFINITY, f64::max);
            locus.offset += maximum;
            let relative: Vec<_> = scores.iter().map(|s| s.0 - maximum).collect();
            for (i, &score) in relative.iter().enumerate() {
                locus.totals[i] += score;
                locus.spanning[i] += u32::from(scores[i].1);
                locus.primary_spanning[i] += u32::from(scores[i].2);
                for (j, &other) in relative.iter().enumerate() {
                    let delta = score - other;
                    locus.delta_squares[i][j] += delta * delta;
                    locus.strong_votes[i][j] += u32::from(delta >= 2.3 && scores[i].1);
                    locus.primary_strong_votes[i][j] += u32::from(delta >= 2.3 && scores[i].2);
                }
            }
        }
        Ok(())
    })?;
    for locus in &mut loci {
        if locus.exhausted || locus.length_limited || locus.routes.is_empty() {
            continue;
        }
        let best = (0..locus.routes.len())
            .max_by(|&a, &b| locus.totals[a].total_cmp(&locus.totals[b]).then(b.cmp(&a)))
            .unwrap();
        let wins = (0..locus.routes.len()).filter(|&i| i != best).all(|i| {
            let difference = locus.totals[best] - locus.totals[i];
            stable_margin(
                difference,
                locus.delta_squares[best][i],
                locus.alignment_limited,
            ) && minor_supported(
                locus.strong_votes[best][i],
                locus.strong_votes[i][best],
                locus.alignment_limited,
                settings,
            )
        });
        if wins && locus.spanning[best] >= settings.min_junction_fragments {
            locus.selected = Some(best);
        }
    }
    let complete_cycles =
        circles::resolve(&mut loci, graph, reads, out, settings, limits, &training)?;
    let (mut primary_candidates, primary_choices) =
        primary::local_candidates(graph, &loci, &anchor_order, overlap, settings);
    primary_candidates.extend(complete_cycles.primaries);
    let mut connections = HashMap::new();
    let mut incoming = HashMap::<usize, usize>::new();
    for (i, locus) in loci.iter().enumerate() {
        if let Some(selected) = locus.selected {
            connections.insert(locus.from, (i, selected));
            *incoming.entry(locus.routes[selected].to).or_default() += 1;
        }
    }
    connections.retain(|_, route| incoming[&loci[route.0].routes[route.1].to] == 1);
    let mut consumed = HashSet::new();
    let mut candidates = BTreeSet::new();
    let mut supported_cycles = BTreeSet::new();
    let starts: Vec<_> = anchor_order
        .iter()
        .copied()
        .filter(|id| incoming.get(id) != Some(&1))
        .chain(anchor_order.iter().copied())
        .collect();
    for first in starts {
        if consumed.contains(&first) {
            continue;
        }
        let mut sequence = graph.segments[first].sequence.clone();
        let mut node = first;
        let mut circular = false;
        loop {
            consumed.insert(node);
            let Some(&(locus, route)) = connections.get(&node) else {
                break;
            };
            let route = &loci[locus].routes[route];
            if consumed.contains(&route.to) && route.to != first {
                break;
            }
            for &id in &route.nodes {
                sequence.extend_from_slice(&graph.segments[id].sequence[overlap..]);
                consumed.insert(id);
            }
            if route.to == first {
                circular = true;
                sequence.truncate(sequence.len() - overlap);
                break;
            }
            sequence.extend_from_slice(&graph.segments[route.to].sequence[overlap..]);
            node = route.to;
        }
        let sequence = normalize(&sequence, circular);
        if circular {
            supported_cycles.insert(sequence.clone());
        }
        candidates.insert((sequence, circular));
    }
    for (id, segment) in graph.segments.iter().enumerate() {
        if !consumed.contains(&id) {
            candidates.insert((
                normalize(&segment.sequence, segment.circular),
                segment.circular,
            ));
        }
    }
    let mut reports: Vec<_> = loci.iter().enumerate().map(|(li,l)|json!({"source_segment":l.from+1,"candidate_independent_log_offset":l.offset,"score_definition":"sum of per-fragment log scores minus that fragment maximum over all paths","source_length":graph.segments[l.from].sequence.len(),"search_states":l.states,"search_exhausted":l.exhausted,"search_cuts":{"visited_states":l.state_limited,"candidate_paths":l.path_limited,"cycle_pre_expansion":l.cycle_expansion_limited,"pending_states":l.frontier_limited,"intermediate_length":l.length_limited},"peak_pending_states":l.frontier_peak,"frontier_limited":l.frontier_limited,"length_limited":l.length_limited,"repeat_copy_model":repeat_depth.report(l,graph,occurrence),"dead_end_walks":l.dead_ends,"fragments":l.fragments,"unplaced_pairs":l.unplaced,"r1_only_fragments":l.r1_only,"r2_only_fragments":l.r2_only,"r2_anchored_fragments":l.r2_anchored,"alignment_limited_pairs":l.alignment_limited,"unscored_maximum_total_contrast":unscored_contrast_bound(l.alignment_limited),"selected_path":l.selected.map(|i|i+1),"primary_path":primary_choices[li].map(|i|i+1),"paths":l.routes.iter().enumerate().map(|(i,r)|json!({"path":i+1,"target_segment":r.to+1,"target_length":graph.segments[r.to].sequence.len(),"intermediate_segments":r.nodes.iter().map(|id|id+1).collect::<Vec<_>>(),"connecting_bases":r.right_start-r.left_end,"spanning_interval":{"start_0based":r.spanning_interval().0,"end_0based_exclusive":r.spanning_interval().1,"basis":if r.snp_interval.is_some(){"single_snp_shared_flanks"}else{"graph_corridor"}},"score":l.totals[i],"spanning_fragments":l.spanning[i],"strong_votes_against_alternatives":l.strong_votes[i],"primary_spanning_interval":{"start_0based":r.primary_spanning_interval().0,"end_0based_exclusive":r.primary_spanning_interval().1,"basis":if r.primary_interval.is_some(){"shared_flank_allele_partial_primary_only"}else if r.snp_interval.is_some(){"single_snp_shared_flanks"}else{"graph_corridor"}},"primary_spanning_fragments":l.primary_spanning[i],"primary_strong_votes_against_alternatives":l.primary_strong_votes[i]})).collect::<Vec<_>>()})).collect();
    for (locus, report) in loci.iter().zip(&mut reports) {
        let scored = !locus.exhausted && !locus.length_limited && !locus.routes.is_empty();
        report["read_scoring"] = json!(if scored {
            "computed"
        } else if locus.exhausted || locus.length_limited {
            "skipped_search_incomplete"
        } else {
            "skipped_no_candidates"
        });
        if !scored {
            for field in [
                "candidate_independent_log_offset",
                "fragments",
                "unplaced_pairs",
                "r1_only_fragments",
                "r2_only_fragments",
                "r2_anchored_fragments",
                "alignment_limited_pairs",
                "unscored_maximum_total_contrast",
            ] {
                report[field] = Value::Null;
            }
            for path in report["paths"].as_array_mut().expect("path reports") {
                for field in [
                    "score",
                    "spanning_fragments",
                    "strong_votes_against_alternatives",
                    "primary_spanning_fragments",
                    "primary_strong_votes_against_alternatives",
                ] {
                    path[field] = Value::Null;
                }
            }
            if let Some(cycles) = report["repeat_copy_model"]
                .get_mut("simple_cycles")
                .and_then(Value::as_array_mut)
            {
                for cycle in cycles {
                    for candidate in cycle["candidates"]
                        .as_array_mut()
                        .expect("repeat candidates")
                    {
                        candidate["fragment_score"] = Value::Null;
                        candidate["spanning_fragments"] = Value::Null;
                    }
                }
            }
        }
    }
    Ok(Resolution {
        supported_cycles,
        primary_candidates,
        candidates: candidates
            .into_iter()
            .map(|(sequence, circular)| Candidate { sequence, circular })
            .collect(),
        report: json!({"method":"complete graph paths between sample-only anchors; full-query affine-gap alignment and one empirical paired-fragment score per name","fragment_model":training.report,"end_assembly":end_assembly,"complete_cycle_comparison":complete_cycles.report,"primary_selection":{"local_rule":"highest fragment score with spanning support and observed strong-vote majority; reciprocal views prefer a competitive comparison over a sole route, then strict over relaxed; equal conflicts stay open","allele_spanning":"same-endpoint shared-flank allele interval for partial primary only; strict join intervals unchanged","minimum_spanning_fragments":settings.min_junction_fragments,"limited_alignments":"omitted from every candidate; reported as unknown, not counted as observed opposition for primary output","incomplete_path_search":"no relaxed join","cycle_rule":"one highest-scoring enumerated complete cycle with supported corridor boundaries; unresolved alleles and phase retained locally; strict final criteria unchanged"},"score_is_correctness_probability":false,"alignment":{"mode":"full-query affine-gap fitting; free target flanks, no query clipping","seed_length":17,"seed_start_margin":12,"maximum_seed_starts_per_strand":256,"gap_open_log_score":-8.0,"gap_extend_log_score":-2.0,"quality_error_floor":0.001,"outlier_log_floor":-12.0},"anchor_contexts":anchor_contexts,"anchor_seed_policy":"prefer 31; sample-unique 47/63 contexts rescue ambiguous short seeds without changing graph k or join evidence","fragment_pool":"identical names for all paths, recruited by unique source-anchor adaptive contexts on either mate; prefer R1 when both qualify, fixed anchor before alignment, no conditioning on the other mate placement; empirical weights are not renormalized over found hits","maximum_paths_per_anchor":MAX_PATHS,"maximum_search_states_per_anchor":MAX_STATES,"maximum_pending_states_per_anchor":MAX_STATES,"search_order":"increasing accumulated intermediate bases, then full node history and target ID; scheduling only, not a selection score","maximum_intermediate_bases":max_bases,"context_bases":context,"selection":{"minimum_score_margin":10.0,"minimum_margin_over_root_sum_squared_deltas":3.0,"strong_fragment_log_margin":2.3,"minimum_spanning_fragments":settings.min_junction_fragments,"max_strong_minor_major_ratio":settings.max_error_fraction},"scan":scan,"distinct_scored_fragment_names":names.len(),"anchor_count":anchors.len(),"short_terminal_anchor_count":short_terminal_anchors,"selected_connections":connections.len(),"loci":reports}),
    })
}

#[cfg(test)]
mod tests {
    use super::super::threading::Segment;
    use super::*;

    fn walk_graph(sequences: &[&[u8]], outgoing: Vec<Vec<usize>>) -> PathGraph {
        PathGraph {
            segments: sequences
                .iter()
                .map(|sequence| Segment {
                    sequence: sequence.to_vec(),
                    circular: false,
                })
                .collect(),
            outgoing,
        }
    }

    #[test]
    fn finite_walks_preserve_equal_cost_histories_and_ignore_edge_order() {
        let mut graph = walk_graph(
            &[b"GGGA", b"ACA", b"AGA", b"ATA", b"ACT", b"ATATA"],
            vec![vec![1, 2, 4, 5], vec![3], vec![3], vec![4], vec![], vec![3]],
        );
        let anchors = HashSet::from([0, 4]);
        let a = enumerate(&graph, 0, &anchors, 1, 100, 20, 20);
        assert!(!a.exhausted && !a.length_limited);
        let histories: BTreeSet<_> = a.routes.iter().map(|r| r.nodes.clone()).collect();
        assert_eq!(
            histories,
            BTreeSet::from([vec![], vec![1, 3], vec![2, 3], vec![5, 3]])
        );
        let sequences: Vec<_> = a.routes.iter().map(|r| r.sequence.clone()).collect();
        assert_eq!(
            sequences,
            [
                b"GGGACT".to_vec(),
                b"GGGACATACT".to_vec(),
                b"GGGAGATACT".to_vec(),
                b"GGGATATATACT".to_vec()
            ]
        );
        for next in &mut graph.outgoing {
            next.reverse();
        }
        let b = enumerate(&graph, 0, &anchors, 1, 100, 20, 20);
        assert_eq!(
            sequences,
            b.routes
                .iter()
                .map(|r| r.sequence.clone())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn complex_repeat_does_not_starve_short_exit_or_become_complete() {
        // The self edge and two-node return make this SCC non-simple. DFS
        // spends its path quota at the length boundary before this short exit.
        let graph = walk_graph(
            &[b"GGGA", b"AAA", b"ACA", b"ATG"],
            vec![vec![1], vec![1, 2, 3], vec![1], vec![]],
        );
        let anchors = HashSet::from([0, 3]);
        assert!(repeats::cycles(&graph, &anchors, 1).is_empty());
        let locus = enumerate(&graph, 0, &anchors, 1, 200, 20, 20);
        assert_eq!(locus.routes[0].nodes, [1]);
        assert_eq!(locus.routes[0].sequence, b"GGGAAATG");
        assert_eq!(locus.routes.len(), MAX_PATHS);
        assert!(locus.exhausted && locus.path_limited);
        assert!(!locus.cycle_expansion_limited && !locus.state_limited);
        assert!(locus.selected.is_none());
        assert!(locus.routes.iter().any(|r| r.nodes == [1, 2, 1]));
        assert!(locus
            .routes
            .windows(2)
            .all(|r| r[0].sequence.len() <= r[1].sequence.len()));
    }

    #[test]
    fn pending_walk_limit_is_explicitly_incomplete() {
        let graph = PathGraph {
            segments: (0..MAX_STATES + 2)
                .map(|_| Segment {
                    sequence: b"AA".to_vec(),
                    circular: false,
                })
                .collect(),
            outgoing: std::iter::once((1..MAX_STATES + 2).collect())
                .chain((0..MAX_STATES + 1).map(|_| Vec::new()))
                .collect(),
        };
        let locus = enumerate(&graph, 0, &HashSet::from([0, MAX_STATES + 1]), 1, 10, 2, 2);
        assert!(locus.exhausted && locus.frontier_limited);
        assert_eq!(locus.frontier_peak, MAX_STATES);
        assert!(locus.routes.is_empty() && locus.selected.is_none());
    }

    #[test]
    fn finite_candidate_set_matches_dag_oracle_and_reverse_views() {
        let graph = walk_graph(
            &[b"GGGA", b"ACA", b"AGA", b"ATA", b"ACT", b"ATATA", b"AAC"],
            vec![
                vec![1, 2, 4, 5],
                vec![3],
                vec![3],
                vec![4, 6],
                vec![],
                vec![3],
                vec![],
            ],
        );
        fn oracle(
            graph: &PathGraph,
            node: usize,
            seq: Vec<u8>,
            out: &mut BTreeSet<(usize, Vec<u8>)>,
        ) {
            if graph.outgoing[node].is_empty() {
                out.insert((node, seq));
            } else {
                for &next in &graph.outgoing[node] {
                    let next_seq = [seq.as_slice(), &graph.segments[next].sequence[1..]].concat();
                    oracle(graph, next, next_seq, out);
                }
            }
        }
        let mut expected = BTreeSet::new();
        oracle(&graph, 0, graph.segments[0].sequence.clone(), &mut expected);
        let anchors = HashSet::from([0, 4, 6]);
        let result = enumerate(&graph, 0, &anchors, 1, 100, 100, 20);
        assert!(!result.exhausted && !result.length_limited);
        assert_eq!(
            expected,
            result
                .routes
                .iter()
                .map(|r| (r.to, r.sequence.clone()))
                .collect()
        );
        let mut reverse = PathGraph {
            segments: graph
                .segments
                .iter()
                .map(|s| Segment {
                    sequence: kmer::reverse_complement(&s.sequence),
                    circular: false,
                })
                .collect(),
            outgoing: vec![Vec::new(); graph.segments.len()],
        };
        for (from, next) in graph.outgoing.iter().enumerate() {
            for &to in next {
                reverse.outgoing[to].push(from);
            }
        }
        for from in [4, 6] {
            let reversed = enumerate(&reverse, from, &anchors, 1, 100, 100, 20);
            let expected: BTreeSet<_> = expected
                .iter()
                .filter(|(to, _)| *to == from)
                .map(|(_, seq)| kmer::reverse_complement(seq))
                .collect();
            assert!(!reversed.exhausted && !reversed.length_limited);
            assert_eq!(
                expected,
                reversed.routes.iter().map(|r| r.sequence.clone()).collect()
            );
        }
    }

    #[test]
    fn visited_state_budget_cannot_claim_no_path_or_complete_search() {
        // A complex repeat with a dead exit has many histories but no target.
        // Its queued frontier stays below the independent pending-state cap.
        let graph = walk_graph(
            &[b"AA", b"AA", b"AA", b"AA", b"AA"],
            vec![vec![1], vec![1, 2, 3], vec![1], vec![], vec![]],
        );
        let locus = enumerate(&graph, 0, &HashSet::from([0, 4]), 1, 1_000_000, 2, 2);
        assert!(locus.exhausted && locus.state_limited);
        assert!(!locus.path_limited && !locus.frontier_limited && !locus.length_limited);
        assert!(locus.routes.is_empty());
    }

    #[test]
    fn incomplete_fragments_cannot_be_treated_as_zero_contrast() {
        assert!(stable_margin(83.0, 11.0 * (83.0f64 / 11.0).powi(2), 0));
        assert!(!stable_margin(83.0, 11.0 * (83.0f64 / 11.0).powi(2), 18));
        assert!(stable_margin(1000.0, 1000.0, 1));
    }

    #[test]
    fn existing_seeded_decoy_does_not_suppress_better_seedless_mate() {
        use crate::io::fastq::Record;
        let mut seed = 918731u64;
        let mut dna = |n| {
            (0..n)
                .map(|_| {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    b"ACGT"[(seed % 4) as usize]
                })
                .collect::<Vec<_>>()
        };
        let anchor = dna(100);
        let query = dna(80);
        let change = |s: &mut Vec<u8>, p: usize| s[p] = if s[p] == b'A' { b'C' } else { b'A' };
        let mut seedless = query.clone();
        let mut qualities = vec![b'I'; 80];
        for p in [15, 31, 47, 63] {
            change(&mut seedless, p);
            qualities[p] = b'+';
        }
        let mut decoy = query.clone();
        for p in [25, 26, 27] {
            change(&mut decoy, p);
        }
        let mut alternate = query.clone();
        for p in [25, 26] {
            change(&mut alternate, p);
        }
        let spacer = dna(100);
        let sequences = [seedless, alternate].map(|s| {
            [
                anchor.clone(),
                spacer.clone(),
                s,
                spacer.clone(),
                decoy.clone(),
            ]
            .concat()
        });
        let indexes: Vec<_> = sequences
            .iter()
            .map(|s| Index::new(s, false, 100))
            .collect();
        let record = |sequence, quality| Record {
            header: "@f".into(),
            plus: "+".into(),
            sequence,
            quality,
        };
        let a = record(anchor, vec![b'I'; 100]);
        qualities.reverse();
        let b = record(kmer::reverse_complement(&query), qualities);
        let mut placements: Vec<_> = indexes.iter().map(|i| (i.align(&a), i.align(&b))).collect();
        assert!(!placements[0].1.hits.is_empty());
        assert!(placements[0].1.hits.iter().all(|h| h.start != 200));
        rescue_mates(
            &indexes,
            &mut placements,
            &mut Prepared::new(&b),
            &[(0, false, 10); 2],
            None,
            100,
        );
        let best = placements[0]
            .1
            .hits
            .iter()
            .max_by(|a, b| a.logp.total_cmp(&b.logp))
            .unwrap();
        assert_eq!(best.start, 200);
        assert_eq!(best.edits, 4);
    }

    #[test]
    fn snp_spanning_interval_preserves_other_exits_and_structural_alternatives() {
        let route = |sequence| Route {
            circle_length: None,
            to: 1,
            nodes: vec![],
            sequence,
            left_end: 70,
            right_start: 131,
            snp_interval: None,
            primary_interval: None,
        };
        let a = vec![b'A'; 220];
        let mut c = a.clone();
        c[100] = b'C';
        let mut paths = vec![route(a.clone()), route(c.clone())];
        identify_single_snp(0, &mut paths);
        assert!(paths.iter().all(|p| p.spanning_interval() == (100, 101)));
        for case in 0..6 {
            let mut paths = vec![route(a.clone()), route(c.clone())];
            match case {
                0 => paths[1].sequence[120] = b'G',
                1 => paths[1].sequence.insert(101, b'T'),
                2 => paths[1].to = 2,
                3 => paths.iter_mut().for_each(|p| p.to = 0),
                4 => paths[1].sequence = a.clone(),
                _ => paths.iter_mut().for_each(|p| p.left_end = 110),
            }
            identify_single_snp(0, &mut paths);
            assert!(
                paths.iter().all(|p| p.snp_interval.is_none()),
                "case={case}"
            );
        }
    }

    #[test]
    fn duplicate_placements_cannot_increase_fragment_score() {
        let hit = super::super::align::Hit {
            start: 0,
            end: 100,
            reverse: false,
            logp: 0.0,
            edits: 0,
            traceback: None,
        };
        let route = Route {
            circle_length: None,
            to: 1,
            nodes: vec![],
            sequence: vec![b'A'; 120],
            left_end: 40,
            right_start: 60,
            snp_interval: None,
            primary_interval: None,
        };
        let one = Placements {
            hits: vec![hit.clone()],
            exhausted: false,
            complete: false,
        };
        let many = Placements {
            hits: vec![hit; 10],
            exhausted: false,
            complete: false,
        };
        assert_eq!(
            joint(&one, &one, &route, None, 0.0, (0, false, 3)),
            joint(&many, &many, &route, None, 0.0, (0, false, 3))
        );
        assert_eq!(
            joint(&one, &one, &route, None, 0.0, (100, false, 3)).0,
            -12.0
        );
    }
    #[test]
    fn repeat_revisits_are_candidates_and_length_bound_is_not_completion() {
        let segment = |n| Segment {
            sequence: vec![b'A'; n],
            circular: false,
        };
        let graph = PathGraph {
            segments: vec![segment(300), segment(82), segment(300)],
            outgoing: vec![vec![1], vec![1, 2], vec![]],
        };
        let locus = enumerate(&graph, 0, &HashSet::from([0, 2]), 62, 100, 300, 100);
        assert!(locus.length_limited);
        assert_eq!(locus.routes.len(), 5);
        assert!(locus.routes.iter().any(|r| r.nodes.len() == 5));
        assert!(locus.selected.is_none());
        let large = enumerate(&graph, 0, &HashSet::from([0, 2]), 62, 100_000, 300, 100);
        assert!(large.exhausted && large.length_limited && large.cycle_expansion_limited);
        assert_eq!(large.routes.len(), MAX_PATHS);
        assert!(large.states <= MAX_PATHS + 1);
    }
    fn allele_routes() -> Vec<Route> {
        [10, 13]
            .into_iter()
            .map(|n| Route {
                circle_length: None,
                to: 1,
                nodes: vec![],
                sequence: [vec![b'C'; 100], vec![b'T'; n], vec![b'G'; 150]].concat(),
                left_end: 20,
                right_start: 210 + n,
                snp_interval: None,
                primary_interval: None,
            })
            .collect()
    }

    #[test]
    fn primary_allele_spans_ambiguous_indel_positions_without_changing_strict_interval() {
        let mut routes = allele_routes();
        identify_primary_allele(0, &mut routes);
        assert_eq!(routes[0].primary_spanning_interval(), (100, 110));
        assert_eq!(routes[1].primary_spanning_interval(), (100, 113));
        assert_eq!(routes[0].spanning_interval(), (20, 220));
        assert_eq!(routes[1].spanning_interval(), (20, 223));
        let mut reverse = allele_routes();
        for r in &mut reverse {
            let n = r.sequence.len();
            r.sequence = kmer::reverse_complement(&r.sequence);
            (r.left_end, r.right_start) = (n - r.right_start, n - r.left_end);
        }
        identify_primary_allele(0, &mut reverse);
        for (f, r) in routes.iter().zip(&reverse) {
            let (a, b) = f.primary_spanning_interval();
            assert_eq!(
                r.primary_spanning_interval(),
                (f.sequence.len() - b, f.sequence.len() - a)
            );
        }
    }

    #[test]
    fn local_allele_evidence_does_not_grant_a_strict_join_or_cross_a_different_exit() {
        use super::super::align::Hit;
        let mut routes = allele_routes();
        identify_primary_allele(0, &mut routes);
        let hit = Hit {
            start: 80,
            end: 145,
            reverse: false,
            logp: 0.0,
            edits: 0,
            traceback: None,
        };
        let placed = Placements {
            hits: vec![hit],
            exhausted: false,
            complete: true,
        };
        for route in &routes {
            let single_score = single(&placed, route, 0.0, (80, false, 3));
            let joint_score = joint(&placed, &placed, route, None, 0.0, (80, false, 3));
            assert!(!single_score.1 && single_score.2);
            assert!(!joint_score.1 && joint_score.2);
        }
        for case in 0..5 {
            let mut routes = allele_routes();
            match case {
                0 => routes[1].to = 2,
                1 => routes.iter_mut().for_each(|r| r.to = 0),
                2 => routes[1].circle_length = Some(260),
                3 => routes[1].sequence = routes[0].sequence.clone(),
                _ => routes.iter_mut().for_each(|r| r.left_end = 105),
            }
            identify_primary_allele(0, &mut routes);
            assert!(
                routes.iter().all(|r| r.primary_interval.is_none()),
                "case={case}"
            );
        }
    }

    #[test]
    fn separated_alleles_require_one_fragment_across_all_differences() {
        let mut routes = allele_routes();
        routes[1].sequence = routes[0].sequence.clone();
        routes[1].sequence[100] = b'A';
        routes[1].sequence[200] = b'A';
        identify_primary_allele(0, &mut routes);
        assert_eq!(routes[0].primary_spanning_interval(), (100, 201));
        let placed = Placements {
            hits: vec![super::super::align::Hit {
                start: 80,
                end: 145,
                reverse: false,
                logp: 0.0,
                edits: 0,
                traceback: None,
            }],
            exhausted: false,
            complete: true,
        };
        assert!(!single(&placed, &routes[0], 0.0, (80, false, 3)).2);
    }
}
