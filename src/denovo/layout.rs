//! Optional, explicitly inferred repeat layouts. Coverage ranks sample-graph
//! walks; it never grants structural validation or final-assembly status.
use super::{compact_at, paths::primary, Edge, Unitig};
use crate::{
    cli::{GraphSettings, Limits, Reads},
    io::{fasta::write_fasta_records, fastq},
    kmer::{self, graph as g},
};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

const DEPTH_K: usize = 31;
const FALLBACK_DEPTH_K: usize = 15;
const MAX_STATES: usize = 50_000;
const MAX_PENDING: usize = 1024;
const MAX_COMPLETIONS: usize = 256;
const MAX_ADDED_BASES: usize = 10_000;
const MAX_CODE_COPIES: usize = 4096;
const ANCHOR_WINDOW: usize = 500;

type Occurrences = BTreeMap<kmer::Code, usize>;

fn occurrences(sequence: &[u8], circular: bool) -> Occurrences {
    occurrences_at(sequence, circular, DEPTH_K)
}

fn occurrences_at(sequence: &[u8], circular: bool, k: usize) -> Occurrences {
    let mut sequence = sequence.to_vec();
    if circular && !sequence.is_empty() {
        let prefix: Vec<_> = sequence.iter().copied().cycle().take(k - 1).collect();
        sequence.extend(prefix);
    }
    let mut counts = BTreeMap::new();
    kmer::visit(&sequence, None, k, 0, |_, code| {
        *counts.entry(kmer::canonical(code, k)).or_default() += 1;
    });
    counts
}

struct Depth {
    k: usize,
    counts: HashMap<kmer::Code, u64>,
    mean: f64,
}

impl Depth {
    fn copies(&self, code: kmer::Code) -> f64 {
        self.counts.get(&code).copied().unwrap_or(0) as f64 / self.mean
    }

    fn baseline(&self, sequence: &[u8], unique: &Occurrences) -> (usize, f64) {
        let codes: BTreeSet<_> = occurrences_at(sequence, false, self.k)
            .into_keys()
            .filter(|code| unique.get(code) == Some(&1))
            .collect();
        let sum: u64 = codes
            .iter()
            .map(|c| self.counts.get(c).copied().unwrap_or(0))
            .sum();
        (codes.len(), sum as f64 / codes.len().max(1) as f64)
    }
}

struct Proposal {
    sequence: Vec<u8>,
    junctions: [usize; 2],
    regions: Vec<primary::Region>,
    k: usize,
    walk: Vec<usize>,
    prefix_trim: usize,
    suffix_trim: usize,
    minimum_pair_support: u32,
    weakest_corridor_pair_count: u32,
    shortest_path_seed: bool,
    coverage_search: bool,
    source_graphs: BTreeSet<(usize, u32)>,
    generation_methods: BTreeSet<String>,
    graph_file: Option<String>,
    closing_overlap: usize,
}

struct Pending {
    cost: f64,
    node: g::Code,
    bases: usize,
    walk: Vec<usize>,
    used: Occurrences,
}

impl PartialEq for Pending {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Pending {}
impl PartialOrd for Pending {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Pending {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cost
            .total_cmp(&other.cost)
            .then(self.bases.cmp(&other.bases))
            .then(self.walk.cmp(&other.walk))
            .then(self.node.cmp(&other.node))
    }
}

fn code(sequence: &[u8]) -> g::Code {
    let mut result = g::Code::default();
    g::visit(sequence, None, sequence.len(), 0, |_, c| result = c);
    result
}

// The retained backbone and both flanks are sample sequences. Move a boundary
// inward only to escape a non-unique endpoint, never to match a reference.
fn unique_flanks(backbone: &[u8], k: usize) -> Option<(usize, usize)> {
    let overlap = k - 1;
    if backbone.len() < 2 * k {
        return None;
    }
    let mut counts = HashMap::<g::Code, usize>::new();
    g::visit(backbone, None, overlap, 0, |_, c| {
        *counts.entry(g::canonical(c, overlap)).or_default() += 1;
    });
    let bound = ANCHOR_WINDOW.min((backbone.len() - 2 * k) / 2);
    let left = (0..=bound).find(|&p| {
        counts.get(&g::canonical(code(&backbone[p..p + overlap]), overlap)) == Some(&1)
    })?;
    let right_trim = (0..=bound).find(|&p| {
        let end = backbone.len() - p;
        counts.get(&g::canonical(code(&backbone[end - overlap..end]), overlap)) == Some(&1)
    })?;
    Some((left, backbone.len() - right_trim))
}

fn corridor(
    counts: &HashMap<g::Code, u32>,
    backbone: &[u8],
    k: usize,
    minimum: u32,
) -> Vec<Unitig> {
    let overlap = k - 1;
    let start = code(&backbone[backbone.len() - overlap..]);
    let target = code(&backbone[..overlap]);
    // Cut one central unique backbone node to prevent a second whole-backbone
    // traversal. Blocking every unique interior word would also remove genuine
    // NCR copies of sequence seen once in the backbone, particularly at low k.
    let mut frequency = HashMap::<g::Code, usize>::new();
    g::visit(backbone, None, overlap, 0, |_, c| {
        *frequency.entry(g::canonical(c, overlap)).or_default() += 1;
    });
    let mut guard = None;
    g::visit(backbone, None, overlap, 0, |p, c| {
        if p >= ANCHOR_WINDOW.min(backbone.len() / 4)
            && p + overlap + ANCHOR_WINDOW.min(backbone.len() / 4) < backbone.len()
            && frequency[&g::canonical(c, overlap)] == 1
        {
            let rank = (p.abs_diff(backbone.len() / 2), p, c);
            if guard.is_none_or(|old| rank < old) {
                guard = Some(rank);
            }
        }
    });
    let blocked: HashSet<_> = guard
        .into_iter()
        .flat_map(|(_, _, c)| [c, g::reverse_code(c, overlap)])
        .collect();
    let mut edges = Vec::new();
    for (&canonical, &count) in counts.iter().filter(|(_, n)| **n >= minimum) {
        for c in [canonical, g::reverse_code(canonical, k)] {
            let (from, to) = (c.prefix(), c.suffix(overlap));
            if from != target && !blocked.contains(&from) && !blocked.contains(&to) {
                edges.push(Edge {
                    code: c,
                    from,
                    to,
                    count,
                });
            }
        }
    }
    edges.sort_by_key(|e| e.code);
    let (units, _) = compact_at(&edges, k, &HashSet::from([start, target]));
    let mut incoming = HashMap::<_, Vec<_>>::new();
    let mut outgoing = HashMap::<_, Vec<_>>::new();
    for (id, u) in units.iter().enumerate() {
        incoming.entry(u.to).or_default().push(id);
        outgoing.entry(u.from).or_default().push(id);
    }
    let mut reaches = HashSet::from([target]);
    let mut pending = vec![target];
    while let Some(node) = pending.pop() {
        for &id in incoming.get(&node).into_iter().flatten() {
            if reaches.insert(units[id].from) {
                pending.push(units[id].from);
            }
        }
    }
    let mut accessible = HashSet::from([start]);
    pending.push(start);
    while let Some(node) = pending.pop() {
        for &id in outgoing.get(&node).into_iter().flatten() {
            if reaches.contains(&units[id].to) && accessible.insert(units[id].to) {
                pending.push(units[id].to);
            }
        }
    }
    units
        .into_iter()
        .filter(|u| accessible.contains(&u.from) && reaches.contains(&u.to))
        .collect()
}

fn contributions(unitig: &Unitig, k: usize, depth_k: usize) -> Occurrences {
    // Only newly appended bases contribute. The k-1 overlap is not counted
    // again, even when the layout graph uses a k larger than the depth k.
    occurrences_at(&unitig.sequence[k - depth_k..], false, depth_k)
}

fn incremental_cost(before: usize, added: usize, estimated: f64) -> f64 {
    (before as f64 + added as f64 - estimated).powi(2) - (before as f64 - estimated).powi(2)
}

fn proposal(
    units: &[Unitig],
    backbone: &[u8],
    k: usize,
    trims: (usize, usize),
    walk: Vec<usize>,
    shortest_path_seed: bool,
) -> Proposal {
    let overlap = k - 1;
    let mut sequence = backbone.to_vec();
    for &id in &walk {
        sequence.extend_from_slice(&units[id].sequence[overlap..]);
    }
    sequence.truncate(sequence.len() - overlap);
    let (sequence, junctions, regions) = normalize_proposal(sequence, backbone.len(), overlap, k);
    Proposal {
        sequence,
        junctions,
        regions,
        k,
        prefix_trim: trims.0,
        suffix_trim: trims.1,
        minimum_pair_support: 0,
        weakest_corridor_pair_count: walk
            .iter()
            .map(|&id| units[id].min_count)
            .min()
            .unwrap_or(0),
        walk,
        shortest_path_seed,
        coverage_search: !shortest_path_seed,
        source_graphs: BTreeSet::new(),
        generation_methods: BTreeSet::from(["legacy_coverage_walk".to_owned()]),
        graph_file: None,
        closing_overlap: overlap,
    }
}

fn normalize_proposal(
    sequence: Vec<u8>,
    backbone_len: usize,
    overlap: usize,
    k: usize,
) -> (Vec<u8>, [usize; 2], Vec<primary::Region>) {
    let length = sequence.len();
    let regions = vec![
        primary::Region {
            start: backbone_len.saturating_sub(overlap).min(length - 1),
            end: length,
            evidence: format!("coverage_layout_k{k}"),
            reason: "coverage_inferred_repeat_layout",
            alternatives: 0,
        },
        primary::Region {
            start: 0,
            end: overlap.min(length),
            evidence: format!("coverage_layout_k{k}"),
            reason: "coverage_inferred_repeat_layout",
            alternatives: 0,
        },
    ];
    let reverse = kmer::reverse_complement(&sequence);
    let offset = super::rotation_start(&sequence);
    let rev_offset = super::rotation_start(&reverse);
    let forward = [&sequence[offset..], &sequence[..offset]].concat();
    let reversed = [&reverse[rev_offset..], &reverse[..rev_offset]].concat();
    let junctions = [backbone_len.min(length), 0].map(|p| {
        if reversed < forward {
            (2 * length - p - rev_offset) % length
        } else {
            (p + length - offset) % length
        }
    });
    let (sequence, regions) = primary::normalize_regions(sequence, true, regions);
    (sequence, junctions, regions)
}

fn merge_proposal(pool: &mut BTreeMap<Vec<u8>, Proposal>, candidate: Proposal) {
    if let Some(existing) = pool.get_mut(&candidate.sequence) {
        existing.shortest_path_seed |= candidate.shortest_path_seed;
        existing.coverage_search |= candidate.coverage_search;
        existing.source_graphs.extend(candidate.source_graphs);
        existing
            .generation_methods
            .extend(candidate.generation_methods);
    } else {
        pool.insert(candidate.sequence.clone(), candidate);
    }
}

fn shortest_walk(
    units: &[Unitig],
    outgoing: &HashMap<g::Code, Vec<usize>>,
    start: g::Code,
    target: g::Code,
    overlap: usize,
) -> Option<Vec<usize>> {
    let mut distances = HashMap::from([(start, 0usize)]);
    let mut previous = HashMap::<g::Code, usize>::new();
    let mut pending = BTreeSet::from([(0usize, start)]);
    while let Some((distance, node)) = pending.pop_first() {
        if distances.get(&node) != Some(&distance) {
            continue;
        }
        if node == target {
            let mut walk = Vec::new();
            let mut node = node;
            while node != start {
                let &id = previous.get(&node)?;
                walk.push(id);
                node = units[id].from;
            }
            walk.reverse();
            return Some(walk);
        }
        for &id in outgoing.get(&node).into_iter().flatten() {
            let next = units[id].to;
            let cost = distance + units[id].sequence.len() - overlap;
            if cost <= MAX_ADDED_BASES && distances.get(&next).is_none_or(|&old| cost < old) {
                distances.insert(next, cost);
                previous.insert(next, id);
                pending.insert((cost, next));
            }
        }
    }
    None
}

fn search(
    units: &[Unitig],
    backbone: &[u8],
    k: usize,
    depth: &Depth,
    trims: (usize, usize),
) -> (Vec<Proposal>, Value) {
    let overlap = k - 1;
    let start = code(&backbone[backbone.len() - overlap..]);
    let target = code(&backbone[..overlap]);
    let mut outgoing = HashMap::<_, Vec<usize>>::new();
    let additions: Vec<_> = units
        .iter()
        .enumerate()
        .map(|(id, u)| {
            outgoing.entry(u.from).or_default().push(id);
            contributions(u, k, depth.k)
        })
        .collect();
    let mut base = occurrences_at(backbone, false, depth.k);
    // k>31 appends some prefix 31-mers twice before the closing overlap is
    // removed. Subtract those occurrences from the starting state exactly once.
    kmer::visit(&backbone[..overlap], None, depth.k, 0, |_, c| {
        let count = base.get_mut(&kmer::canonical(c, depth.k)).unwrap();
        *count -= 1;
    });
    let mut pending = BTreeSet::from([Pending {
        cost: 0.0,
        node: start,
        bases: 0,
        walk: Vec::new(),
        used: BTreeMap::new(),
    }]);
    let mut states = 0;
    let mut peak = 1;
    let (mut frontier_cut, mut length_cut, mut copy_cut, mut completion_cut) =
        (false, false, false, false);
    let mut found = BTreeMap::new();
    let seed = shortest_walk(units, &outgoing, start, target, overlap);
    let shortest_seed_found = seed.is_some();
    if let Some(walk) = seed {
        let candidate = proposal(units, backbone, k, trims, walk, true);
        found.insert(candidate.sequence.clone(), candidate);
    }
    let mut completion_count = 0;
    while let Some(state) = pending.pop_first() {
        if states >= MAX_STATES {
            pending.insert(state);
            break;
        }
        states += 1;
        if state.node == target {
            let candidate = proposal(units, backbone, k, trims, state.walk, false);
            merge_proposal(&mut found, candidate);
            completion_count += 1;
            if completion_count >= MAX_COMPLETIONS {
                completion_cut = !pending.is_empty();
                break;
            }
            continue;
        }
        for &id in outgoing.get(&state.node).into_iter().flatten() {
            let u = &units[id];
            let bases = state.bases + u.sequence.len() - overlap;
            if bases > MAX_ADDED_BASES {
                length_cut = true;
                continue;
            }
            let mut cost = state.cost;
            let mut allowed = true;
            for (&c, &n) in &additions[id] {
                let before =
                    base.get(&c).copied().unwrap_or(0) + state.used.get(&c).copied().unwrap_or(0);
                let estimated = depth.copies(c);
                // A generous data-derived exploration range, explicitly not a
                // confidence interval or evidence against its unsearched tail.
                let bound = ((estimated * 1.5).ceil() as usize + 2)
                    .min(MAX_CODE_COPIES)
                    .max(base.get(&c).copied().unwrap_or(0));
                if before + n > bound {
                    allowed = false;
                    copy_cut = true;
                    break;
                }
                cost += incremental_cost(before, n, estimated);
            }
            if !allowed {
                continue;
            }
            let mut used = state.used.clone();
            for (&c, &n) in &additions[id] {
                *used.entry(c).or_default() += n;
            }
            let mut walk = state.walk.clone();
            walk.push(id);
            pending.insert(Pending {
                cost,
                node: u.to,
                bases,
                walk,
                used,
            });
            if pending.len() > MAX_PENDING {
                pending.pop_last();
                frontier_cut = true;
            }
            peak = peak.max(pending.len());
        }
    }
    let report = json!({"k":k,"corridor_segments":units.len(),"states":states,"peak_pending":peak,
        "candidates":found.len(),"shortest_seed_found":shortest_seed_found,"search_complete":false,
        "cuts":{"states":states>=MAX_STATES&&!pending.is_empty(),"pending":frontier_cut,
            "added_bases":length_cut,"copy_exploration_bound":copy_cut,"completions":completion_cut},
        "retained_candidates_survive_other_k_limits":true,
        "status":if !found.is_empty(){"inferred_candidates"}else if units.is_empty(){"no_graph_connection"}else{"no_layout_within_search_bounds"}});
    (found.into_values().collect(), report)
}

fn write_graph(
    out: &Path,
    units: &[Unitig],
    backbone: &[u8],
    k: usize,
    depth: &Depth,
    minimum: u32,
) -> Result<()> {
    let mut file = BufWriter::new(File::create(
        out.join(format!("coverage_layout_k{k}_min{minimum}.gfa")),
    )?);
    writeln!(file, "H\tVN:Z:1.0")?;
    writeln!(file, "S\tbackbone\t{}", String::from_utf8_lossy(backbone))?;
    let overlap = k - 1;
    let mut starts = HashMap::<_, Vec<_>>::new();
    let start = code(&backbone[backbone.len() - overlap..]);
    let target = code(&backbone[..overlap]);
    for (id, u) in units.iter().enumerate() {
        starts.entry(u.from).or_default().push(id);
        let c = contributions(u, k, depth.k);
        let windows: usize = c.values().sum();
        let mean = c
            .iter()
            .map(|(code, n)| *n as f64 * depth.counts.get(code).copied().unwrap_or(0) as f64)
            .sum::<f64>()
            / windows.max(1) as f64;
        writeln!(
            file,
            "S\tu{}\t{}\tdp:f:{mean:.6}\tFC:f:{:.6}",
            id + 1,
            String::from_utf8_lossy(&u.sequence),
            u.mean_count
        )?;
        if u.from == start {
            writeln!(file, "L\tbackbone\t+\tu{}\t+\t{overlap}M", id + 1)?;
        }
        if u.to == target {
            writeln!(file, "L\tu{}\t+\tbackbone\t+\t{overlap}M", id + 1)?;
        }
    }
    for (id, u) in units.iter().enumerate() {
        for &next in starts.get(&u.to).into_iter().flatten() {
            writeln!(file, "L\tu{}\t+\tu{}\t+\t{overlap}M", id + 1, next + 1)?;
        }
    }
    file.flush()?;
    Ok(())
}

fn fit(counts: &Occurrences, domain: &BTreeSet<kmer::Code>, depth: &Depth, scale: f64) -> f64 {
    domain
        .iter()
        .map(|&c| (counts.get(&c).copied().unwrap_or(0) as f64 - depth.copies(c) / scale).powi(2))
        .sum::<f64>()
        / domain.len().max(1) as f64
}

pub(super) fn propose(
    backbone: &[u8],
    circular: bool,
    reads: &Reads,
    out: &Path,
    settings: &GraphSettings,
    limits: &Limits,
) -> Result<Value> {
    if circular {
        return Ok(json!({"enabled":true,"status":"primary_already_circular"}));
    }
    let legacy_ks: BTreeSet<_> = [settings.k, 47, 31, 21, 15]
        .into_iter()
        .filter(|&k| k <= settings.k && k >= FALLBACK_DEPTH_K)
        .collect();
    let mut ks = legacy_ks.clone();
    if settings.multi_k_layout {
        ks.extend(super::multik::KS);
    }
    if ks.is_empty() {
        return Ok(json!({"enabled":true,"status":"requires_k_at_least_15"}));
    }
    let mut graphs: BTreeMap<_, HashMap<g::Code, u32>> =
        ks.iter().map(|&k| (k, HashMap::new())).collect();
    let mut depth = Depth {
        k: DEPTH_K,
        counts: HashMap::new(),
        mean: 0.0,
    };
    let mut fallback_depth = Depth {
        k: FALLBACK_DEPTH_K,
        counts: HashMap::new(),
        mean: 0.0,
    };
    let scan = fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
        kmer::visit(
            &a.sequence,
            Some(&a.quality),
            DEPTH_K,
            limits.min_baseq,
            |_, c| {
                let n = depth.counts.entry(kmer::canonical(c, DEPTH_K)).or_default();
                *n = n.saturating_add(1);
            },
        );
        ensure!(
            depth.counts.len() <= limits.max_kmers,
            "coverage layout depth exceeds --max-kmers"
        );
        kmer::visit(
            &a.sequence,
            Some(&a.quality),
            FALLBACK_DEPTH_K,
            limits.min_baseq,
            |_, c| {
                let n = fallback_depth
                    .counts
                    .entry(kmer::canonical(c, FALLBACK_DEPTH_K))
                    .or_default();
                *n = n.saturating_add(1);
            },
        );
        ensure!(
            fallback_depth.counts.len() <= limits.max_kmers,
            "coverage layout fallback depth exceeds --max-kmers"
        );
        for (&k, counts) in &mut graphs {
            let mut fragment = HashSet::new();
            for read in [a, b] {
                g::visit(
                    &read.sequence,
                    Some(&read.quality),
                    k,
                    limits.min_baseq,
                    |_, c| {
                        fragment.insert(g::canonical(c, k));
                    },
                );
            }
            for c in fragment {
                ensure!(
                    counts.contains_key(&c) || counts.len() < limits.max_kmers,
                    "coverage layout graph exceeds --max-kmers at k={k}"
                );
                let n = counts.entry(c).or_default();
                *n = n.saturating_add(1);
            }
        }
        Ok(())
    })?;
    let backbone_counts = occurrences(backbone, false);
    let (baseline_windows, mean) = depth.baseline(backbone, &backbone_counts);
    depth.mean = mean;
    if baseline_windows < DEPTH_K || mean <= 0.0 {
        return Ok(
            json!({"enabled":true,"status":"insufficient_single_copy_depth","baseline_unique_31mers":baseline_windows,"baseline_mean_r1_depth":mean}),
        );
    }
    let left_depth = depth.baseline(
        &backbone[..backbone.len().min(ANCHOR_WINDOW)],
        &backbone_counts,
    );
    let right_depth = depth.baseline(
        &backbone[backbone.len().saturating_sub(ANCHOR_WINDOW)..],
        &backbone_counts,
    );
    // Calibrate against nearby unique flanks where possible; keep the whole
    // backbone mean too so local sampling bias is visible in the report.
    let local_windows = left_depth.0 + right_depth.0;
    if backbone.len() > 2 * ANCHOR_WINDOW && local_windows >= DEPTH_K {
        depth.mean = (left_depth.0 as f64 * left_depth.1 + right_depth.0 as f64 * right_depth.1)
            / local_windows as f64;
    }
    if depth.mean <= 0.0 {
        depth.mean = mean;
    }
    let fallback_unique = occurrences_at(backbone, false, FALLBACK_DEPTH_K);
    let fl = fallback_depth.baseline(
        &backbone[..backbone.len().min(ANCHOR_WINDOW)],
        &fallback_unique,
    );
    let fr = fallback_depth.baseline(
        &backbone[backbone.len().saturating_sub(ANCHOR_WINDOW)..],
        &fallback_unique,
    );
    fallback_depth.mean = if backbone.len() > 2 * ANCHOR_WINDOW && fl.0 + fr.0 >= FALLBACK_DEPTH_K {
        (fl.0 as f64 * fl.1 + fr.0 as f64 * fr.1) / (fl.0 + fr.0) as f64
    } else {
        fallback_depth.baseline(backbone, &fallback_unique).1
    };
    let mut proposals = BTreeMap::<Vec<u8>, Proposal>::new();
    let mut searches = Vec::new();
    // The forced-layout fallback may use singleton-supported sample edges.
    // This is isolated from every strict graph and publication decision.
    for (pass, minimum) in [settings.min_kmer_count, 1].into_iter().enumerate() {
        if pass > 0 && (settings.min_kmer_count == 1 || !proposals.is_empty()) {
            break;
        }
        for (&k, counts) in graphs.iter().rev() {
            if !legacy_ks.contains(&k) {
                continue;
            }
            if k < DEPTH_K && (!proposals.is_empty() || fallback_depth.mean <= 0.0) {
                searches.push(json!({"k":k,"status":"short_k_fallback_not_needed_or_no_depth"}));
                continue;
            }
            let guiding_depth = if k >= DEPTH_K {
                &depth
            } else {
                &fallback_depth
            };
            let Some((left, right)) = unique_flanks(backbone, k) else {
                searches
                    .push(json!({"k":k,"status":"no_unique_flank_pair","search_complete":false}));
                continue;
            };
            let kept = &backbone[left..right];
            let units = corridor(counts, kept, k, minimum);
            tracing::info!(
                k,
                minimum,
                corridor_segments = units.len(),
                "fitting optional coverage layout"
            );
            write_graph(out, &units, kept, k, guiding_depth, minimum)?;
            let (found, mut report) = search(
                &units,
                kept,
                k,
                guiding_depth,
                (left, backbone.len() - right),
            );
            report["minimum_pair_support"] = json!(minimum);
            report["guidance_depth_k"] = json!(guiding_depth.k);
            report["guidance_mean_r1_depth"] = json!(guiding_depth.mean);
            report["prefix_trim_bases"] = json!(left);
            report["suffix_trim_bases"] = json!(backbone.len() - right);
            searches.push(report);
            for mut p in found {
                p.minimum_pair_support = minimum;
                p.source_graphs.insert((k, minimum));
                merge_proposal(&mut proposals, p);
            }
        }
    }
    if settings.multi_k_layout {
        for (pass, minimum) in [settings.min_kmer_count, 1].into_iter().enumerate() {
            if pass > 0
                && (settings.min_kmer_count == 1
                    || proposals
                        .values()
                        .any(|p| p.source_graphs.iter().any(|(_, m)| *m > 1)))
            {
                break;
            }
            // One common long-flank pair for the additional k searches. Actual
            // read windows determine whether a high k graph is usable.
            let anchor_k = (graphs
                .iter()
                .rev()
                .find(|(_, c)| c.values().any(|&n| n >= minimum))
                .map(|(&k, _)| k)
                .unwrap_or(31))
            .max(31);
            if let Some((left, right)) = unique_flanks(backbone, anchor_k) {
                let kept = &backbone[left..right];
                let flank_len = anchor_k - 1;
                let from = &kept[kept.len() - flank_len..];
                let to = &kept[..flank_len];
                let mut outside = occurrences(kept, false);
                for flank in [from, to] {
                    for (c, n) in occurrences(flank, false) {
                        *outside
                            .get_mut(&c)
                            .expect("flank belongs to retained backbone") -= n;
                    }
                }
                let guide = super::multik::Guide {
                    k: DEPTH_K,
                    observed: &depth.counts,
                    outside,
                    mean: depth.mean,
                };
                let guard = if kept.len() >= 4 * flank_len {
                    let start = kept.len() / 2;
                    let candidate = &kept[start..start + flank_len];
                    let key = g::canonical(code(candidate), flank_len);
                    let mut count = 0;
                    g::visit(kept, None, flank_len, 0, |_, c| {
                        count += usize::from(g::canonical(c, flank_len) == key)
                    });
                    (count == 1).then_some(candidate)
                } else {
                    None
                };
                for &k in ks.iter().rev() {
                    let graph_file = format!("coverage_multik_k{k}_min{minimum}.gfa");
                    let path = out.join(&graph_file);
                    let (found, mut report) = super::multik::search(super::multik::Request {
                        graphs: &graphs,
                        k,
                        minimum,
                        left: from,
                        right: to,
                        guard,
                        guide: &guide,
                        max_span: MAX_ADDED_BASES + 2 * flank_len,
                        max_paths: MAX_COMPLETIONS,
                        graph_file: Some(&path),
                    })?;
                    report["graph_file"] = if path.exists() {
                        json!(graph_file)
                    } else {
                        Value::Null
                    };
                    report["prefix_trim_bases"] = json!(left);
                    report["suffix_trim_bases"] = json!(backbone.len() - right);
                    searches.push(report);
                    for p in found {
                        let mut sequence = kept.to_vec();
                        sequence.extend_from_slice(&p.sequence[flank_len..]);
                        sequence.truncate(sequence.len() - flank_len);
                        let (sequence, junctions, regions) =
                            normalize_proposal(sequence, kept.len(), flank_len, k);
                        merge_proposal(
                            &mut proposals,
                            Proposal {
                                sequence,
                                junctions,
                                regions,
                                k,
                                walk: p.walk,
                                prefix_trim: left,
                                suffix_trim: backbone.len() - right,
                                minimum_pair_support: minimum,
                                weakest_corridor_pair_count: p.support,
                                shortest_path_seed: p.schedules.contains("nearby_exit_seed"),
                                coverage_search: p.schedules.contains("coverage_residual"),
                                source_graphs: BTreeSet::from([(k, minimum)]),
                                generation_methods: p
                                    .schedules
                                    .into_iter()
                                    .map(|s| format!("multi_k_{s}"))
                                    .collect(),
                                graph_file: Some(graph_file.clone()),
                                closing_overlap: flank_len,
                            },
                        );
                    }
                }
            } else {
                searches.push(json!({"method":"long_flanks_cross_k_read_context","status":"no_shared_long_flanks","search_complete":false}));
            }
        }
    }
    let mut proposals: Vec<_> = proposals.into_values().collect();
    let counts: Vec<_> = proposals
        .iter()
        .map(|p| occurrences(&p.sequence, true))
        .collect();
    // One fixed union, including unchosen graph spellings. A short candidate
    // cannot improve its fit by omitting inconvenient motif k-mers from its own
    // denominator. Constant backbone codes are excluded from the ranking.
    let mut domain = BTreeSet::new();
    for c in &counts {
        domain.extend(c.keys().copied());
    }
    domain.retain(|c| {
        counts
            .iter()
            .any(|n| n.get(c).copied().unwrap_or(0) != backbone_counts.get(c).copied().unwrap_or(0))
    });
    let mut order: Vec<_> = (0..proposals.len()).collect();
    order.sort_by(|&a, &b| {
        fit(&counts[a], &domain, &depth, 1.0)
            .total_cmp(&fit(&counts[b], &domain, &depth, 1.0))
            .then(proposals[a].sequence.cmp(&proposals[b].sequence))
    });
    let sensitivity:Vec<_>=[0.8,1.0,1.2].into_iter().filter_map(|scale| {
        let winner=(0..proposals.len()).min_by(|&a,&b|fit(&counts[a],&domain,&depth,scale).total_cmp(&fit(&counts[b],&domain,&depth,scale)).then(proposals[a].sequence.cmp(&proposals[b].sequence)))?;
        Some(json!({"baseline_multiplier":scale,"selected_rank":order.iter().position(|&i|i==winner).unwrap()+1,"length":proposals[winner].sequence.len()}))
    }).collect();
    let total = proposals.len();
    let shortlist = super::shortlist::select(
        &proposals
            .iter()
            .map(|p| super::shortlist::Candidate {
                sequence: &p.sequence,
                k: p.k,
                circular: true,
            })
            .collect::<Vec<_>>(),
        &order,
        8,
    );
    let compared = &shortlist.indices;
    let mut flank_seeds = HashSet::new();
    let kept_left = proposals.iter().map(|p| p.prefix_trim).max().unwrap_or(0);
    let kept_right = proposals
        .iter()
        .map(|p| backbone.len() - p.suffix_trim)
        .min()
        .unwrap_or(backbone.len());
    let kept = if kept_left < kept_right {
        &backbone[kept_left..kept_right]
    } else {
        &[]
    };
    for flank in [
        &kept[..kept.len().min(ANCHOR_WINDOW)],
        &kept[kept.len().saturating_sub(ANCHOR_WINDOW)..],
    ] {
        kmer::visit(flank, None, 21, 0, |_, c| {
            flank_seeds.insert(kmer::canonical(c, 21));
        });
    }
    // Freeze eligibility on the full generated union, independently of the
    // shortlist. Changing which eight paths are compared cannot change reads.
    for p in &proposals {
        let counts = occurrences_at(&p.sequence, true, 21);
        flank_seeds.retain(|c| counts.get(c) == Some(&1));
    }
    let mut hypotheses = vec![super::comparison::Hypothesis {
        sequence: backbone,
        circular: false,
        junctions: [0, backbone.len()],
    }];
    hypotheses.extend(compared.iter().map(|&i| super::comparison::Hypothesis {
        sequence: &proposals[i].sequence,
        circular: true,
        junctions: proposals[i].junctions,
    }));
    let mut paired_comparison = if total > 0 {
        super::comparison::compare(
            &hypotheses,
            &flank_seeds,
            backbone,
            reads,
            settings,
            limits,
            out,
        )?
    } else {
        Value::Null
    };
    if !paired_comparison.is_null() {
        paired_comparison["coverage_ranks"] = json!(std::iter::once(0)
            .chain(
                compared
                    .iter()
                    .map(|i| order.iter().position(|p| p == i).unwrap() + 1)
            )
            .collect::<Vec<_>>());
        paired_comparison["cohort_generated_paths"] = json!(total);
        paired_comparison["shortlist"] = shortlist.report;
        let best = (paired_comparison["status"] == "evaluated")
            .then(|| paired_comparison["best_candidate"].as_u64().unwrap() as usize);
        paired_comparison["read_preferred_coverage_rank"] = best.map_or(Value::Null, |i| {
            paired_comparison["coverage_ranks"][i].clone()
        });
        paired_comparison["read_preferred_is_linear_baseline"] = json!(best.map(|i| i == 0));
    }
    // This first comparison is diagnostic. Preserve the coverage representative
    // and open primary while exposing competing read evidence, including a
    // baseline that is not a circular assembly hypothesis.
    let rankings:Vec<_>=order.iter().enumerate().map(|(rank,&i)| {
        let p=&proposals[i];
        json!({"id":format!("coverage_layout_{}",rank+1),"length":p.sequence.len(),"k":p.k,
            "unobserved_candidate_31mers":counts[i].keys().filter(|c|depth.counts.get(c).copied().unwrap_or(0)==0).count(),
            "copy_rmse":fit(&counts[i],&domain,&depth,1.0).sqrt(),"shortest_path_seed":p.shortest_path_seed,"coverage_search":p.coverage_search,"source_graphs":p.source_graphs.iter().map(|&(k,m)| json!({"k":k,"minimum_pair_support":m})).collect::<Vec<_>>(),"minimum_pair_support":p.minimum_pair_support,"weakest_corridor_pair_count":p.weakest_corridor_pair_count,
            "graph_file":p.graph_file.clone().unwrap_or_else(||format!("coverage_layout_k{}_min{}.gfa",p.k,p.minimum_pair_support)),
            "generation_methods":p.generation_methods,
            "closing_overlap_bases":p.closing_overlap,
            "walk":p.graph_file.is_none().then(||"backbone+".to_owned()).into_iter().chain(p.walk.iter().map(|n|format!("u{}+",n+1))).collect::<Vec<_>>(),
            "walk_spelling":"append each unitig after its k-1 overlap to the retained input backbone, then remove closing_overlap_bases and normalize the circle",
            "prefix_trim_bases":p.prefix_trim,"suffix_trim_bases":p.suffix_trim})
    }).collect();
    let records:Vec<_>=order.iter().enumerate().map(|(rank,&i)|
        (format!("coverage_layout_{} circular=true length={} final_verified=false evidence=coverage_inferred",rank+1,proposals[i].sequence.len()),proposals[i].sequence.clone())).collect();
    write_fasta_records(&out.join("coverage_layout_candidates.fasta"), &records)?;
    let mut ambiguity = BufWriter::new(File::create(out.join("coverage_layout_ambiguity.tsv"))?);
    writeln!(
        ambiguity,
        "sequence\tstart\tend\tevidence\treason\talternatives"
    )?;
    if let Some(&i) = order.first() {
        let best = &mut proposals[i];
        write_fasta_records(&out.join("coverage_layout.fasta"),&[(format!("coverage_layout circular=true length={} final_verified=false evidence=coverage_inferred",best.sequence.len()),best.sequence.clone())])?;
        for r in &best.regions {
            writeln!(
                ambiguity,
                "coverage_layout\t{}\t{}\t{}\t{}\t{}",
                r.start + 1,
                r.end,
                r.evidence,
                r.reason,
                total.saturating_sub(1)
            )?;
        }
        let mut coverage = BufWriter::new(File::create(out.join("coverage_layout_depth.tsv"))?);
        writeln!(coverage,"position\tinstances_in_layout\tobserved_r1_31mer_windows\tdepth_per_layout_copy\tbaseline_r1_depth")?;
        let extended = [
            best.sequence.clone(),
            best.sequence
                .iter()
                .copied()
                .cycle()
                .take(DEPTH_K - 1)
                .collect(),
        ]
        .concat();
        let mut rows = Vec::new();
        kmer::visit(&extended, None, DEPTH_K, 0, |p, c| {
            let c = kmer::canonical(c, DEPTH_K);
            let instances = counts[i][&c];
            let observed = depth.counts.get(&c).copied().unwrap_or(0);
            rows.push((
                p + 1,
                instances,
                observed,
                observed as f64 / instances as f64,
            ));
        });
        for (p, n, observed, per_copy) in rows {
            writeln!(
                coverage,
                "{p}\t{n}\t{observed}\t{per_copy:.6}\t{:.6}",
                depth.mean
            )?;
        }
        coverage.flush()?;
    }
    ambiguity.flush()?;
    let read_validation = if let Some(&i) = order.first() {
        let directory = out.join("coverage_layout_validation");
        std::fs::create_dir_all(&directory)?;
        let mut validation = super::validate(
            reads,
            &proposals[i].sequence,
            false,
            &directory,
            settings,
            limits,
        )?
        .0;
        validation["candidate_id"] = json!("coverage_layout_1");
        validation["sequence_length"] = json!(proposals[i].sequence.len());
        validation["scope"]=json!("coverage representative only; matching depth and output origin, not every inferred join or structural uniqueness");
        validation
    } else {
        Value::Null
    };
    let report = json!({"schema":"mitovar-rust.coverage-layout.v1","enabled":true,"status":if total>0 {"coverage_inferred_layout"} else {"no_connected_layout"},
        "method":"sample-derived unique flanks; multi-k raw supported graph corridors; bounded coverage-guided walks; common 31-mer occurrence copy fit",
        "file":if total>0 {Some("coverage_layout.fasta")}else{None},"candidates_file":"coverage_layout_candidates.fasta",
        "depth_file":if total>0 {Some("coverage_layout_depth.tsv")}else{None},
        "ambiguity_file":"coverage_layout_ambiguity.tsv","final_verified":false,"structural_uniqueness_claim":false,
        "coverage_statistic":"R1 31-mer window occurrences, repeated offsets retained; R2 omitted to avoid overlapping-mate double counting",
        "baseline":"arithmetic mean R1 depth of unique 31-mers in both 500-base backbone flanks; whole-backbone fallback when flanks are unavailable",
        "whole_backbone_unique_31mers":baseline_windows,"whole_backbone_mean_r1_depth":mean,"baseline_mean_r1_depth":depth.mean,
        "left_flank_unique_31mers":left_depth.0,"left_flank_mean_r1_depth":left_depth.1,
        "right_flank_unique_31mers":right_depth.0,"right_flank_mean_r1_depth":right_depth.1,
        "copy_fit_domain_31mers":domain.len(),"candidate_count":total,"retained_candidates":rankings,
        "paired_comparison":paired_comparison,"read_validation":read_validation,"layout_selection":"coverage representative; paired comparison diagnostic only",
        "baseline_sensitivity":sensitivity,"searches":searches,"scan":scan,
        "multi_k_layout":{"enabled":settings.multi_k_layout,"requested_k":ks,"candidate_k_eligibility":"observed valid sample windows and sufficient flank length","additional_state_budget_per_k":if settings.multi_k_layout {50_000}else{0},"legacy_candidates_retained":true,"common_fragment_comparison_after_union":true,"main_backbone_k":settings.k},
        "limits":{"states_per_k":MAX_STATES,"pending_per_k":MAX_PENDING,"completions_per_k":MAX_COMPLETIONS,
            "additional_bases":MAX_ADDED_BASES,"maximum_code_copies":MAX_CODE_COPIES,"retained_fasta_records":total,"export_scope":"all generated unique spellings"},
        "search_complete":false,"shortest_path_seed":"one shortest connected graph walk per k and support threshold is retained independently of coverage exploration caps, then ranked by the common coverage fit","weak_support_policy":"only when the normal-support candidate union is empty, retry sample graphs with minimum pair-record support 1; all such layouts remain inferred","short_k_policy":"k21/15 used only if all larger k graphs found no layout; final cross-k ranking always uses the common 31-mer model","copy_bound":"ceil(1.5 * observed R1 depth / baseline) + 2; exploration heuristic, unsearched copy tails remain unresolved",
        "limitations":["coverage estimates copy usage, not repeat ordering or same-fragment phase","singleton-edge fallback can follow sequencing errors; minimum pair support and weakest corridor support are reported","AT-dependent sampling, PCR duplicates and other genomic copies can bias depth","no missing sequence or graph edge is invented; no reference, annotation or expected genome length enters selection","a bounded candidate union is ranked; the best candidate is not a certified optimum"]});
    crate::write_json(&out.join("coverage_layout.json"), &report)?;
    Ok(report)
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
                b"ACGT"[(seed % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn graph_walk_and_occurrence_depth_recover_a_tandem_layout_but_not_a_missing_bridge() {
        let backbone = dna(1100, 172839);
        // The NCR reuses a segment that occurs once inside the backbone; it
        // must remain traversable despite being locally unique there.
        let genome = [
            backbone.clone(),
            backbone[320..400].to_vec(),
            b"AT".repeat(90),
        ]
        .concat();
        let mut depth = Depth {
            k: DEPTH_K,
            counts: HashMap::new(),
            mean: 0.0,
        };
        let mut counts = HashMap::<g::Code, u32>::new();
        for start in 0..genome.len() {
            let read: Vec<_> = (0..100)
                .map(|p| genome[(start + p) % genome.len()])
                .collect();
            kmer::visit(&read, None, DEPTH_K, 0, |_, c| {
                *depth.counts.entry(kmer::canonical(c, DEPTH_K)).or_default() += 1;
            });
            let mut fragment = HashSet::new();
            g::visit(&read, None, 63, 0, |_, c| {
                fragment.insert(g::canonical(c, 63));
            });
            for c in fragment {
                *counts.entry(c).or_default() += 1;
            }
        }
        depth.mean = depth
            .baseline(&backbone[..200], &occurrences(&backbone, false))
            .1;
        assert_eq!(depth.mean, 70.0);
        let units = corridor(&counts, &backbone, 63, 3);
        let (found, report) = search(&units, &backbone, 63, &depth, (0, 0));
        let truth = super::super::normalize(&genome, true);
        assert!(found.iter().any(|p| p.sequence == truth));
        assert_eq!(report["search_complete"], false);
        let domain: BTreeSet<_> = found
            .iter()
            .flat_map(|p| occurrences(&p.sequence, true).into_keys())
            .collect();
        let best = found
            .iter()
            .min_by(|a, b| {
                fit(&occurrences(&a.sequence, true), &domain, &depth, 1.0).total_cmp(&fit(
                    &occurrences(&b.sequence, true),
                    &domain,
                    &depth,
                    1.0,
                ))
            })
            .unwrap();
        assert_eq!(best.sequence, truth);
        let low_depth = Depth {
            k: DEPTH_K,
            counts: HashMap::new(),
            mean: 1.0,
        };
        let (forced, limited) = search(&units, &backbone, 63, &low_depth, (0, 0));
        assert!(!forced.is_empty());
        assert!(forced.iter().any(|p| p.shortest_path_seed));
        assert_eq!(limited["shortest_seed_found"], true);
        assert_eq!(limited["search_complete"], false);
        // Every inferred base and join is spelled by the reported unitig walk.
        let mut spelled = backbone.clone();
        for &id in &best.walk {
            spelled.extend_from_slice(&units[id].sequence[62..]);
        }
        spelled.truncate(spelled.len() - 62);
        assert_eq!(super::super::normalize(&spelled, true), truth);
        // A genuine sequence gap remains a gap, regardless of the depth model.
        let mut weak = counts.clone();
        g::visit(&genome[1080..1150], None, 63, 0, |_, c| {
            weak.insert(g::canonical(c, 63), 1);
        });
        assert!(search(
            &corridor(&weak, &backbone, 63, 3),
            &backbone,
            63,
            &depth,
            (0, 0)
        )
        .0
        .is_empty());
        let supported_once = search(
            &corridor(&weak, &backbone, 63, 1),
            &backbone,
            63,
            &depth,
            (0, 0),
        )
        .0;
        let recovered = supported_once.iter().find(|p| p.sequence == truth).unwrap();
        assert_eq!(recovered.weakest_corridor_pair_count, 1);
        let mut missing = counts;
        g::visit(&genome[1080..1150], None, 63, 0, |_, c| {
            missing.remove(&g::canonical(c, 63));
        });
        let disconnected = corridor(&missing, &backbone, 63, 3);
        assert!(search(&disconnected, &backbone, 63, &depth, (0, 0))
            .0
            .is_empty());
    }

    #[test]
    fn mean_depth_counts_repeat_occurrences_not_presence_or_mates() {
        let counts = occurrences(&[b'A'; 80], false);
        assert_eq!(counts.values().sum::<usize>(), 50);
        let depth = Depth {
            k: DEPTH_K,
            counts: HashMap::from([(0, 250)]),
            mean: 50.0,
        };
        assert_eq!(depth.copies(0), 5.0);
        assert!(incremental_cost(4, 1, 5.0) < 0.0);
        assert!(incremental_cost(5, 1, 5.0) > 0.0);
    }

    #[test]
    fn common_fit_domain_penalizes_missing_sequence_and_distinguishes_copy_counts() {
        let depth = Depth {
            k: DEPTH_K,
            counts: HashMap::from([(0, 300), (1, 100)]),
            mean: 100.0,
        };
        let domain = BTreeSet::from([0, 1]);
        let correct = BTreeMap::from([(0, 3), (1, 1)]);
        assert_eq!(fit(&correct, &domain, &depth, 1.0), 0.0);
        for wrong in [
            BTreeMap::from([(0, 1), (1, 1)]),
            BTreeMap::from([(0, 4), (1, 1)]),
            BTreeMap::from([(0, 3)]),
        ] {
            assert!(fit(&wrong, &domain, &depth, 1.0) > 0.0);
        }
    }
}
