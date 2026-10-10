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
    collections::{HashMap, HashSet},
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
    outgoing: &HashMap<graph_kmer::Code, Vec<usize>>,
    incoming: &HashMap<graph_kmer::Code, usize>,
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
        if next.is_none_or(|ids| ids.len() != 1) || incoming.get(&edge.to) != Some(&1) {
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
        let next = next.expect("one outgoing edge")[0];
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
    counts: &HashMap<graph_kmer::Code, u32>,
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
    let mut outgoing: HashMap<graph_kmer::Code, Vec<usize>> = HashMap::new();
    let mut incoming: HashMap<graph_kmer::Code, usize> = HashMap::new();
    for (id, edge) in edges.iter().enumerate() {
        outgoing.entry(edge.from).or_default().push(id);
        *incoming.entry(edge.to).or_default() += 1;
    }
    // Split unitigs at requested sample-derived corridor boundaries.
    for &node in boundaries {
        incoming.insert(node, 2);
    }
    let branch_nodes = outgoing
        .keys()
        .chain(incoming.keys())
        .copied()
        .collect::<HashSet<_>>()
        .into_iter()
        .filter(|node| {
            outgoing.get(node).is_some_and(|ids| ids.len() > 1)
                || incoming.get(node).is_some_and(|n| *n > 1)
        })
        .count();
    let mut visited = vec![false; edges.len()];
    let mut unitigs = Vec::new();
    for (id, edge) in edges.iter().enumerate() {
        if !visited[id] && (outgoing[&edge.from].len() != 1 || incoming.get(&edge.from) != Some(&1))
        {
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
fn similar_paths(a: &[u8], b: &[u8], max_edits: usize) -> bool {
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

/// Compare complete compacted paths between the same junctions. Coverage is
/// averaged over path edges, so a single high/low k-mer does not veto a bubble.
/// Only short, similar spellings qualify; loops and conjugate traversals remain.
fn stronger_bubble_path(
    weak_id: usize,
    unitigs: &[Unitig],
    starts: &HashMap<graph_kmer::Code, Vec<usize>>,
    removed: &HashSet<graph_kmer::Code>,
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
            let mut unavailable = false;
            graph_kmer::visit(&u.sequence, None, k, 0, |_, c| {
                unavailable |= removed.contains(&c);
            });
            if unavailable {
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
                        let minimum = route.iter().map(|&p| unitigs[p].min_count).min().unwrap();
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

/// Topology-aware error removal uses local competing support, never a genome length,
/// reference, global coverage cutoff, or a choice of route through an unbridged repeat.
fn simplify(
    edges: &mut Vec<Edge>,
    k: usize,
    fraction: f64,
    bubble_fraction: f64,
) -> (Vec<Unitig>, usize, Value, Vec<bubbles::Collapsed>) {
    let (mut unitigs, mut branches) = compact(edges, k);
    let (raw_kmers, raw_branches) = (edges.len() / 2, branches);
    let (mut rounds, mut tip_kmers, mut bubble_kmers, mut weak_link_kmers) = (0, 0, 0, 0);
    let mut alternatives = std::collections::BTreeMap::new();
    if fraction > 0.0 {
        loop {
            let mut starts: HashMap<graph_kmer::Code, Vec<usize>> = HashMap::new();
            let mut ends: HashMap<graph_kmer::Code, Vec<usize>> = HashMap::new();
            for (id, unitig) in unitigs.iter().enumerate() {
                starts.entry(unitig.from).or_default().push(id);
                ends.entry(unitig.to).or_default().push(id);
            }
            let mut removed = HashSet::new();
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
                        stronger_bubble_path(id, &unitigs, &starts, &removed, k, bubble_fraction)
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

fn write_graph(out: &Path, unitigs: &[Unitig], k: usize) -> Result<()> {
    let mut writer = BufWriter::new(File::create(out.join("assembly_graph.gfa"))?);
    writeln!(writer, "H\tVN:Z:1.0")?;
    let mut starts: HashMap<graph_kmer::Code, Vec<usize>> = HashMap::new();
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

pub fn assemble(
    reads: &Reads,
    validation_reads: &Reads,
    out: &Path,
    settings: &GraphSettings,
    limits: &Limits,
    recruitment: Option<Value>,
    bam_baited: bool,
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
    let mut counts = HashMap::<graph_kmer::Code, u32>::new();
    let scan = fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
        let mut fragment = HashSet::new();
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
        for code in fragment {
            ensure!(
                counts.contains_key(&code) || counts.len() < limits.max_kmers,
                "assembly exceeds --max-kmers; no final assembly published"
            );
            let count = counts.entry(code).or_default();
            *count = count.saturating_add(1);
        }
        Ok(())
    })?;
    info!(
        pairs = scan.pairs,
        distinct_kmers = counts.len(),
        "building sample-only graph"
    );
    let (mut edges, _, _) = graph(&counts, settings.k, settings.min_kmer_count);
    let (unitigs, branches, simplification, collapsed_bubbles) = simplify(
        &mut edges,
        settings.k,
        settings.max_error_fraction,
        settings.max_bubble_fraction,
    );
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
    let complete = read_validated && selected_bubbles.is_empty();
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
        "assembly_r1":reads.r1.canonicalize()?,"assembly_r2":reads.r2.canonicalize()?,"assembly_reads":scan,
        "distinct_kmers":counts.len(),"retained_canonical_kmers":edges.len()/2,"minimum_retained_pair_count":edges.iter().map(|e|e.count).min(),"directed_branch_nodes":branches,
        "graph_simplification":simplification,
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

    fn counts(sequences: &[(&[u8], u32)], k: usize) -> HashMap<graph_kmer::Code, u32> {
        let mut result = HashMap::new();
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
        let mut starts = HashMap::<_, Vec<_>>::new();
        for (id, u) in unitigs.iter().enumerate() {
            starts.entry(u.from).or_default().push(id);
        }
        assert!(stronger_bubble_path(0, &unitigs, &starts, &HashSet::new(), k, 0.2).is_none());
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
}
