//! Small-k paths repair missing connections between existing graph ends.
//! This proposes sequence only; the shared fragment scorer decides each join.
use super::*;
use crate::denovo::threading::Segment;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

type Bridge = (usize, usize, Vec<u8>);
const MAX_ENDPOINT_PAIRS: usize = 10_000;

struct SearchGraph<'a> {
    outgoing: &'a HashMap<kmer::Code, Vec<(kmer::Code, u8)>>,
    reaches_target: HashSet<kmer::Code>,
}

impl<'a> SearchGraph<'a> {
    fn new(
        outgoing: &'a HashMap<kmer::Code, Vec<(kmer::Code, u8)>>,
        targets: &HashMap<kmer::Code, Vec<usize>>,
        k: usize,
    ) -> Self {
        // Reverse reachability is an over-approximation: the exact main-k
        // target flank is still checked when a walk arrives. A node with no
        // route to even the target's small-k suffix cannot spell a connector.
        // Derive its four possible predecessors instead of duplicating all
        // edges in another adjacency map.
        let mut reaches_target: HashSet<_> = targets.keys().copied().collect();
        let mut pending: Vec<_> = reaches_target.iter().copied().collect();
        while let Some(node) = pending.pop() {
            for base in 0..4u128 {
                let predecessor = (base << (2 * (k - 2))) | (node >> 2);
                if outgoing
                    .get(&predecessor)
                    .is_some_and(|edges| edges.iter().any(|&(next, _)| next == node))
                    && reaches_target.insert(predecessor)
                {
                    pending.push(predecessor);
                }
            }
        }
        Self {
            outgoing,
            reaches_target,
        }
    }
}

#[derive(Default)]
struct Search {
    bridges: BTreeSet<Bridge>,
    states: usize,
    limited: bool,
    unreachable_branches_pruned: usize,
    limit_reasons: BTreeSet<&'static str>,
}

fn code(sequence: &[u8]) -> kmer::Code {
    let mut result = 0;
    kmer::visit(sequence, None, sequence.len(), 0, |_, c| result = c);
    result
}

#[derive(Clone, Copy)]
struct SearchLimits {
    k: usize,
    overlap: usize,
    max_bases: usize,
    max_states: usize,
    max_paths: usize,
}

fn target_budget(reachable_targets: usize) -> (usize, usize) {
    if reachable_targets == 0 {
        (0, 0)
    } else {
        (
            MAX_STATES / reachable_targets,
            MAX_PATHS / reachable_targets,
        )
    }
}

impl Search {
    fn mark_limit(&mut self, reason: &'static str) {
        self.limited = true;
        self.limit_reasons.insert(reason);
    }
}

fn search(
    from: usize,
    graph: &PathGraph,
    targets: &HashMap<kmer::Code, Vec<usize>>,
    all_targets: &HashMap<kmer::Code, Vec<usize>>,
    small_graph: &SearchGraph<'_>,
    limits: SearchLimits,
) -> Search {
    let SearchLimits {
        k,
        overlap,
        max_bases,
        max_states,
        max_paths,
    } = limits;
    let source = &graph.segments[from].sequence;
    let initial = source[source.len() - overlap..].to_vec();
    let start = code(&source[source.len() - k + 1..]);
    let mut result = Search::default();
    if !small_graph.reaches_target.contains(&start) {
        result.unreachable_branches_pruned = 1;
        return result;
    }
    if max_states == 0 {
        result.mark_limit("target_state_budget");
    }
    if max_paths == 0 {
        result.mark_limit("target_candidate_budget");
    }
    if result.limited {
        return result;
    }
    let mut stack = vec![(start, initial)];
    while let Some((node, sequence)) = stack.pop() {
        if !small_graph.reaches_target.contains(&node) {
            result.unreachable_branches_pruned += 1;
            continue;
        }
        if result.states == max_states {
            result.mark_limit("target_state_budget");
            break;
        }
        result.states += 1;
        let arrived = sequence.len() > overlap
            && all_targets.get(&node).is_some_and(|ids| {
                ids.iter()
                    .any(|&to| sequence.ends_with(&graph.segments[to].sequence[..overlap]))
            });
        if arrived {
            // A corridor still stops at the FIRST full graph endpoint, even
            // when this search is collecting a different target. Splitting
            // budgets must not introduce walks through intervening anchors.
            if let Some(ids) = targets.get(&node) {
                for &to in ids {
                    if sequence.ends_with(&graph.segments[to].sequence[..overlap]) {
                        result.bridges.insert((from, to, sequence.clone()));
                    }
                }
            }
            if result.bridges.len() >= max_paths && !stack.is_empty() {
                result.mark_limit("target_candidate_budget");
                break;
            }
            continue;
        }
        if let Some(next) = small_graph.outgoing.get(&node) {
            let next: Vec<_> = next
                .iter()
                .filter(|(to, _)| small_graph.reaches_target.contains(to))
                .collect();
            result.unreachable_branches_pruned += small_graph.outgoing[&node].len() - next.len();
            if next.is_empty() {
                continue;
            }
            if sequence.len() - overlap >= max_bases {
                result.mark_limit("maximum_added_bases");
                continue;
            }
            for &&(to, base) in next.iter().rev() {
                let mut extended = sequence.clone();
                extended.push(base);
                stack.push((to, extended));
            }
        }
    }
    result
}

fn reachable_segments(graph: &PathGraph, start: usize) -> HashSet<usize> {
    let mut seen = HashSet::from([start]);
    let mut pending = vec![start];
    while let Some(id) = pending.pop() {
        for &next in &graph.outgoing[id] {
            if seen.insert(next) {
                pending.push(next);
            }
        }
    }
    seen
}

fn connection_role(
    from: usize,
    to: usize,
    reverse_segments: &[Vec<usize>],
    target_reaches_source: bool,
) -> &'static str {
    if reverse_segments[from].contains(&to) {
        "same_end_reverse_complement"
    } else if target_reaches_source {
        "backbone_closure"
    } else {
        "other_end_connection"
    }
}

struct Admission {
    eligible: bool,
    reason: &'static str,
}

fn admission(
    bridge: &Bridge,
    proposed: &BTreeMap<Bridge, BTreeSet<usize>>,
    limited_sources: &BTreeSet<usize>,
    reverse_segments: &[Vec<usize>],
) -> Admission {
    let (from, to, sequence) = bridge;
    if limited_sources.contains(from) {
        return Admission {
            eligible: false,
            reason: "source_competing_search_incomplete",
        };
    }
    let reverse = kmer::reverse_complement(sequence);
    let mut limited_reciprocal = false;
    for &a in &reverse_segments[*to] {
        // Absence of a spelling in a capped reciprocal search is unknown,
        // not evidence that the reciprocal spelling is missing.
        limited_reciprocal |= limited_sources.contains(&a);
        for &b in &reverse_segments[*from] {
            if !limited_sources.contains(&a) && proposed.contains_key(&(a, b, reverse.clone())) {
                return Admission {
                    eligible: true,
                    reason: "retained_for_fragment_scoring",
                };
            }
        }
    }
    Admission {
        eligible: false,
        reason: if limited_reciprocal {
            "reciprocal_competing_search_incomplete"
        } else {
            "reciprocal_spelling_missing"
        },
    }
}

pub(super) fn augment(
    graph: &PathGraph,
    reads: &Reads,
    out: &Path,
    settings: &GraphSettings,
    limits: &Limits,
    max_bases: usize,
) -> Result<(PathGraph, Value)> {
    let overlap = settings.k - 1;
    let mut incoming = vec![0; graph.segments.len()];
    for targets in &graph.outgoing {
        for &id in targets {
            incoming[id] += 1;
        }
    }
    let seeds = seeds::Seeds::new(graph, settings.k.clamp(31, 63));
    let anchors = path_anchors(graph, settings.k, &seeds);
    let mut sources: Vec<_> = anchors
        .iter()
        .copied()
        .filter(|&id| graph.outgoing[id].is_empty())
        .collect();
    let mut sinks: Vec<_> = anchors
        .iter()
        .copied()
        .filter(|&id| incoming[id] == 0)
        .collect();
    sources.sort_unstable();
    sinks.sort_unstable();
    let small_ks: Vec<_> = [47, 31].into_iter().filter(|&k| k < settings.k).collect();
    let endpoint_pairs = sources.len().saturating_mul(sinks.len());
    let endpoint_limited = !small_ks.is_empty() && endpoint_pairs > MAX_ENDPOINT_PAIRS;
    // The main graph has already consumed all input reads. An optional bridge
    // metadata limit must retain that graph, without allocating all pairs or
    // pretending that the omitted competitors have been searched.
    let mut by_sequence = HashMap::<_, Vec<_>>::new();
    for (id, segment) in graph.segments.iter().enumerate() {
        by_sequence
            .entry(segment.sequence.clone())
            .or_default()
            .push(id);
    }
    let reverse_segments: Vec<_> = graph
        .segments
        .iter()
        .map(|segment| {
            by_sequence
                .get(&kmer::reverse_complement(&segment.sequence))
                .cloned()
                .unwrap_or_default()
        })
        .collect();
    let mut roles = BTreeMap::new();
    if !small_ks.is_empty() && !endpoint_limited {
        for &to in &sinks {
            let backbone = reachable_segments(graph, to);
            for &from in &sources {
                roles.insert(
                    (from, to),
                    connection_role(from, to, &reverse_segments, backbone.contains(&from)),
                );
            }
        }
    }
    let mut proposed = BTreeMap::<Bridge, BTreeSet<usize>>::new();
    let mut limited_sources: BTreeSet<_> = sources
        .iter()
        .copied()
        .filter(|_| endpoint_limited)
        .collect();
    let mut limited_pairs = BTreeSet::new();
    let mut searches = Vec::new();
    let mut target_searches = Vec::new();
    for &small_k in small_ks.iter().filter(|_| !endpoint_limited) {
        if sources.is_empty() || sinks.is_empty() {
            break;
        }
        let mut counts = HashMap::<kmer::Code, u32>::new();
        let mut fragment = HashSet::new();
        fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
            fragment.clear();
            for read in [a, b] {
                kmer::visit(
                    &read.sequence,
                    Some(&read.quality),
                    small_k,
                    limits.min_baseq,
                    |_, c| {
                        fragment.insert(kmer::canonical(c, small_k));
                    },
                );
            }
            for &c in &fragment {
                let n = if counts.len() < limits.max_kmers {
                    counts.entry(c).or_default()
                } else {
                    counts.get_mut(&c).ok_or_else(|| {
                        anyhow::anyhow!("end assembly exceeds --max-kmers at k={small_k}")
                    })?
                };
                *n = n.saturating_add(1);
            }
            Ok(())
        })?;
        let mask = (1u128 << (2 * (small_k - 1))) - 1;
        let mut outgoing = HashMap::<_, Vec<_>>::new();
        for (&c, &count) in &counts {
            if count < settings.min_kmer_count {
                continue;
            }
            for c in [c, kmer::reverse_code(c, small_k)] {
                outgoing
                    .entry(c >> 2)
                    .or_default()
                    .push((c & mask, b"ACGT"[(c & 3) as usize]));
            }
        }
        for next in outgoing.values_mut() {
            next.sort_unstable();
            next.dedup();
        }
        let mut targets = HashMap::<_, Vec<_>>::new();
        for &to in &sinks {
            targets
                .entry(code(
                    &graph.segments[to].sequence[overlap - small_k + 1..overlap],
                ))
                .or_default()
                .push(to);
        }
        // One adjacency map and one target reachability set at a time. The
        // first pass records only which sources can reach each target suffix;
        // the second reuses each target's set across every source search.
        // This avoids retaining a graph-sized set for every endpoint.
        let mut plans = Vec::new();
        let mut reachable_counts = HashMap::<usize, usize>::new();
        for &to in &sinks {
            let target = HashMap::from([(
                code(&graph.segments[to].sequence[overlap - small_k + 1..overlap]),
                vec![to],
            )]);
            let small_graph = SearchGraph::new(&outgoing, &target, small_k);
            let reachable: HashSet<_> = sources
                .iter()
                .copied()
                .filter(|&from| {
                    let s = &graph.segments[from].sequence;
                    small_graph
                        .reaches_target
                        .contains(&code(&s[s.len() - small_k + 1..]))
                })
                .collect();
            for &from in &reachable {
                *reachable_counts.entry(from).or_default() += 1;
            }
            plans.push((to, target, reachable));
        }
        let mut totals: BTreeMap<_, Search> = sources
            .iter()
            .map(|&from| (from, Search::default()))
            .collect();
        for (to, target, reachable) in plans {
            let small_graph = if reachable.is_empty() {
                None
            } else {
                Some(SearchGraph::new(&outgoing, &target, small_k))
            };
            for &from in &sources {
                let count = reachable_counts.get(&from).copied().unwrap_or(0);
                let suffix_reachable = reachable.contains(&from);
                // Equal floor shares, without distributing the remainder by
                // segment ID. Total traversal-state and candidate quotas
                // never exceed the original per-source/per-k limits.
                let (max_states, max_paths) =
                    target_budget(if suffix_reachable { count } else { 0 });
                let result = if suffix_reachable {
                    search(
                        from,
                        graph,
                        &target,
                        &targets,
                        small_graph.as_ref().expect("reachable target"),
                        SearchLimits {
                            k: small_k,
                            overlap,
                            max_bases,
                            max_states,
                            max_paths,
                        },
                    )
                } else {
                    Search {
                        unreachable_branches_pruned: 1,
                        ..Search::default()
                    }
                };
                target_searches.push(json!({
                    "k":small_k,"source_segment":from+1,"target_segment":to+1,
                    "role":roles[&(from,to)],"target_suffix_reachable":suffix_reachable,
                    "reachable_targets_for_source":count,"state_budget":max_states,"candidate_budget":max_paths,
                    "states":result.states,"candidate_connections":result.bridges.len(),
                    "status":if !suffix_reachable {"unreachable"} else if result.limited {"limited"} else {"complete"},
                    "search_complete":!result.limited,"search_limited":result.limited,
                    "limit_reasons":result.limit_reasons,"unreachable_branches_pruned":result.unreachable_branches_pruned,
                }));
                let total = totals.get_mut(&from).expect("source summary");
                total.states += result.states;
                total.unreachable_branches_pruned += result.unreachable_branches_pruned;
                total.limited |= result.limited;
                total.limit_reasons.extend(&result.limit_reasons);
                total.bridges.extend(result.bridges.iter().cloned());
                if result.limited {
                    limited_sources.insert(from);
                    limited_pairs.insert((from, to));
                }
                for bridge in result.bridges {
                    proposed.entry(bridge).or_default().insert(small_k);
                }
            }
        }
        for (from, result) in totals {
            searches.push(json!({"k":small_k,"source_segment":from+1,"states":result.states,
                "candidate_connections":result.bridges.len(),"search_limited":result.limited,
                "limit_reasons":result.limit_reasons,"unreachable_branches_pruned":result.unreachable_branches_pruned,
                "reachable_targets":reachable_counts.get(&from).copied().unwrap_or(0),"target_search_count":sinks.len()}));
        }
    }
    // Target budgets isolate discovery, not structural confidence. An
    // incomplete competing target at ANY k still vetoes every join from the
    // source, as does an incomplete reciprocal source. Retained proposals
    // alone cannot hide unenumerated inversions or repeat layouts.
    let retained: Vec<_> = proposed
        .iter()
        .filter(|(bridge, _)| {
            admission(bridge, &proposed, &limited_sources, &reverse_segments).eligible
        })
        .collect();
    let mut proposal_records = Vec::new();
    let mut proposal_reports = Vec::new();
    for (index, (bridge, ks)) in proposed.iter().enumerate() {
        let (from, to, sequence) = bridge;
        let admitted = admission(bridge, &proposed, &limited_sources, &reverse_segments);
        let search_complete = !limited_pairs.contains(&(*from, *to));
        let id = format!("end_bridge_proposal_{}", index + 1);
        let role = roles[&(*from, *to)];
        proposal_records.push((format!(
            "{id} source_segment={} target_segment={} role={role} candidate_only=true search_complete={search_complete} join_eligible={}",
            from+1,to+1,admitted.eligible),sequence.clone()));
        proposal_reports.push(json!({"id":id,"source_segment":from+1,"target_segment":to+1,
            "role":role,"sequence_length":sequence.len(),"k_values":ks,"search_complete":search_complete,
            "join_eligible":admitted.eligible,"admission_reason":admitted.reason}));
    }
    crate::io::fasta::write_fasta_records(
        &out.join("end_bridge_proposals.fasta"),
        &proposal_records,
    )?;
    let mut augmented = graph.clone();
    let mut records = Vec::new();
    let mut connections = Vec::new();
    for ((from, to, sequence), ks) in retained {
        let id = augmented.segments.len();
        augmented.segments.push(Segment {
            sequence: sequence.clone(),
            circular: false,
        });
        augmented.outgoing[*from].push(id);
        augmented.outgoing.push(vec![*to]);
        records.push((
            format!(
                "end_bridge_{} source_segment={} target_segment={} candidate_only=true",
                id + 1,
                from + 1,
                to + 1
            ),
            sequence.clone(),
        ));
        connections.push(
            json!({"segment":id+1,"source_segment":from+1,"target_segment":to+1,
            "sequence_length":sequence.len(),"k_values":ks,"role":roles[&(*from,*to)]}),
        );
    }
    crate::io::fasta::write_fasta_records(&out.join("end_bridge_candidates.fasta"), &records)?;
    let mut gfa = BufWriter::new(File::create(out.join("path_graph.gfa"))?);
    writeln!(gfa, "H\tVN:Z:1.0")?;
    for (id, s) in augmented.segments.iter().enumerate() {
        writeln!(
            gfa,
            "S\t{}\t{}\tTP:Z:{}",
            id + 1,
            String::from_utf8_lossy(&s.sequence),
            if id < graph.segments.len() {
                "sample_graph"
            } else {
                "end_bridge_candidate"
            }
        )?;
    }
    for (from, targets) in augmented.outgoing.iter().enumerate() {
        for &to in targets {
            // Threading already removed the closing overlap from an isolated
            // circular segment. Its self-link therefore consumes zero bases.
            let shared = if from == to && augmented.segments[from].circular {
                0
            } else {
                overlap
            };
            writeln!(gfa, "L\t{}\t+\t{}\t+\t{}M", from + 1, to + 1, shared)?;
        }
    }
    gfa.flush()?;
    Ok((
        augmented,
        json!({"method":"sample-only small-k paths between unique graph ends; full main-k overlap at both flanks; reciprocal candidate spellings; joint fragment scoring decides joins",
        "k_values":small_ks,
        "status":if endpoint_limited{"resource_limited"}else{"searched"},
        "resource_limit":if endpoint_limited{Some(json!({"resource":"endpoint_pairs","observed":endpoint_pairs,"limit":MAX_ENDPOINT_PAIRS,"action":"retain main graph; defer end bridges and complete publication"}))}else{None},
        "source_ends":sources.len(),"target_ends":sinks.len(),"endpoint_pairs":endpoint_pairs,"maximum_endpoint_pairs":MAX_ENDPOINT_PAIRS,"maximum_added_bases":max_bases,
        "maximum_states_per_end_per_k":MAX_STATES,"maximum_candidates_per_end_per_k":MAX_PATHS,
        "reachability_pruning":"target-specific suffix reachability is an over-approximation; require full main-k flank at arrival and stop at the first complete endpoint",
        "target_budget_policy":"equal floor shares of the existing per-source per-k state and candidate budgets among suffix-reachable targets; no ID-based remainder allocation",
        "role_policy":"diagnostic only; reverse-complement endpoints remain eligible proposals, not automatically rejected inversions",
        "proposal_search_complete_scope":"all tested k values for this source-target pair; join eligibility additionally requires complete competing searches from both reciprocal sources",
        "limited_sources":limited_sources.iter().map(|id|id+1).collect::<Vec<_>>(),
        "incomplete_search_policy":"save discovered proposals but withhold graph joins from both views of any limited source; retain the open end",
        "candidate_file":"end_bridge_candidates.fasta","proposal_file":"end_bridge_proposals.fasta","proposals":proposal_reports,"graph_file":"path_graph.gfa",
        "directed_connections_added":connections.len(),"connections":connections,"searches":searches,"target_searches":target_searches}),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn limits(k: usize, overlap: usize, max_bases: usize) -> SearchLimits {
        SearchLimits {
            k,
            overlap,
            max_bases,
            max_states: MAX_STATES,
            max_paths: MAX_PATHS,
        }
    }

    fn small_edges(sequences: &[&[u8]], k: usize) -> HashMap<kmer::Code, Vec<(kmer::Code, u8)>> {
        let mask = (1 << (2 * (k - 1))) - 1;
        let mut edges = HashMap::<_, Vec<_>>::new();
        for sequence in sequences {
            kmer::visit(sequence, None, k, 0, |_, c| {
                edges
                    .entry(c >> 2)
                    .or_default()
                    .push((c & mask, b"ACGT"[(c & 3) as usize]));
            });
        }
        for next in edges.values_mut() {
            next.sort_unstable();
            next.dedup();
        }
        edges
    }

    fn segments(sequences: &[&[u8]]) -> PathGraph {
        PathGraph {
            segments: sequences
                .iter()
                .map(|s| Segment {
                    sequence: s.to_vec(),
                    circular: false,
                })
                .collect(),
            outgoing: vec![vec![]; sequences.len()],
        }
    }

    #[test]
    fn excessive_endpoint_metadata_preserves_graph_without_reading_fastq() {
        let mut seed = 783411u64;
        let mut sequences = Vec::new();
        for _ in 0..51 {
            let sequence: Vec<_> = (0..252)
                .map(|_| {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    b"ACGT"[(seed % 4) as usize]
                })
                .collect();
            sequences.push(kmer::reverse_complement(&sequence));
            sequences.push(sequence);
        }
        let graph = segments(&sequences.iter().map(Vec::as_slice).collect::<Vec<_>>());
        let out = std::env::temp_dir().join(format!(
            "mitovarr-endpoint-cap-missing-input-{}",
            std::process::id()
        ));
        let reads = Reads {
            r1: out.join("missing-r1.fq"),
            r2: out.join("missing-r2.fq"),
        };
        assert!(!reads.r1.exists() && !reads.r2.exists());
        let settings = GraphSettings {
            coverage_layout: false,
            multi_k_layout: false,
            k: 63,
            min_kmer_count: 3,
            max_error_fraction: 0.2,
            max_bubble_fraction: 0.49,
            min_depth: 10,
            min_junction_fragments: 3,
            max_validation_pairs: 10_000,
        };
        let resource_limits = Limits {
            max_kmers: 100_000,
            min_baseq: 0,
        };
        std::fs::create_dir(&out).unwrap();
        let (unchanged, report) =
            augment(&graph, &reads, &out, &settings, &resource_limits, 1000).unwrap();
        assert_eq!(unchanged.outgoing, graph.outgoing);
        assert!(unchanged
            .segments
            .iter()
            .zip(&graph.segments)
            .all(|(a, b)| a.sequence == b.sequence));
        assert_eq!(report["status"], "resource_limited");
        assert_eq!(report["resource_limit"]["observed"], 10_404);
        assert_eq!(report["resource_limit"]["limit"], 10_000);
        assert_eq!(report["limited_sources"].as_array().unwrap().len(), 102);
        assert_eq!(report["target_searches"], json!([]));
        assert_eq!(report["directed_connections_added"], 0);
        assert!(out.join("path_graph.gfa").exists());
        assert!(std::fs::read(out.join("end_bridge_candidates.fasta"))
            .unwrap()
            .is_empty());

        // At main k=31 no optional smaller-k search is requested.
        let shorter_k = GraphSettings { k: 31, ..settings };
        let (unchanged, report) =
            augment(&graph, &reads, &out, &shorter_k, &resource_limits, 1000).unwrap();
        assert_eq!(unchanged.segments.len(), graph.segments.len());
        assert_eq!(unchanged.outgoing, graph.outgoing);
        assert_eq!(report["endpoint_pairs"], 10_404);
        assert_eq!(report["maximum_endpoint_pairs"], 10_000);
        assert_eq!(report["k_values"], json!([]));
        assert_eq!(report["target_searches"], json!([]));
        assert_eq!(report["directed_connections_added"], 0);

        // Exercise the complete output path with full synthetic input: this
        // optional cap must no longer erase the already assembled main graph.
        let mut r1 = File::create(&reads.r1).unwrap();
        let mut r2 = File::create(&reads.r2).unwrap();
        for (i, sequence) in sequences.iter().enumerate() {
            for copy in 0..3 {
                for (mate, writer, sequence) in [
                    (1, &mut r1, sequence.clone()),
                    (2, &mut r2, kmer::reverse_complement(sequence)),
                ] {
                    writeln!(
                        writer,
                        "@end_{i}_{copy}/{mate}\n{}\n+\n{}",
                        String::from_utf8_lossy(&sequence),
                        "I".repeat(sequence.len())
                    )
                    .unwrap();
                }
            }
        }
        drop(r1);
        drop(r2);
        let assembled = out.join("assembled");
        std::fs::create_dir(&assembled).unwrap();
        assert!(!crate::denovo::assemble(
            &reads,
            &reads,
            &assembled,
            &settings,
            &resource_limits,
            None,
            false
        )
        .unwrap());
        let report: Value =
            serde_json::from_slice(&std::fs::read(assembled.join("report.json")).unwrap()).unwrap();
        assert_eq!(report["reason"], "end_search_resource_limit");
        assert_eq!(report["assembly_reads"]["pairs"], sequences.len() * 3);
        assert_eq!(report["final_published"], false);
        assert_eq!(report["sequence_output"]["status"], "partial_available");
        assert!(assembled.join("primary.fasta").exists());
        assert!(assembled.join("candidates.fasta").exists());
        assert!(!assembled.join("consensus.fasta").exists());
        std::fs::remove_dir_all(out).unwrap();
    }

    #[test]
    fn target_shares_keep_the_existing_source_resource_bounds() {
        for count in [0, 1, 2, 3, MAX_PATHS, MAX_PATHS + 1, MAX_STATES + 1] {
            let (states, paths) = target_budget(count);
            assert!(states * count <= MAX_STATES);
            assert!(paths * count <= MAX_PATHS);
            if count > 0 {
                assert_eq!(states, MAX_STATES / count);
                assert_eq!(paths, MAX_PATHS / count);
            }
        }
        let graph = segments(&[b"AAC", b"GTT"]);
        let edges = small_edges(&[b"AACGT"], 3);
        let target = HashMap::from([(code(b"GT"), vec![1])]);
        let small_graph = SearchGraph::new(&edges, &target, 3);
        let (max_states, max_paths) = target_budget(MAX_PATHS + 1);
        let result = search(
            0,
            &graph,
            &target,
            &target,
            &small_graph,
            SearchLimits {
                max_states,
                max_paths,
                ..limits(3, 2, 12)
            },
        );
        assert!(result.limited);
        assert_eq!(result.states, 0);
        assert!(result.bridges.is_empty());
        assert!(result.limit_reasons.contains("target_candidate_budget"));
    }

    #[test]
    fn competing_target_repeat_cannot_consume_a_completed_target_budget() {
        let graph = segments(&[b"AAC", b"GTG", b"TTT"]);
        let edges = small_edges(&[b"AACGT", b"ACCCCTT"], 3);
        let good_target = HashMap::from([(code(b"GT"), vec![1])]);
        let repeat_target = HashMap::from([(code(b"TT"), vec![2])]);
        let all_targets = HashMap::from([(code(b"GT"), vec![1]), (code(b"TT"), vec![2])]);
        let quota = SearchLimits {
            max_states: 25,
            max_paths: 2,
            ..limits(3, 2, 12)
        };
        let good_graph = SearchGraph::new(&edges, &good_target, 3);
        let repeat_graph = SearchGraph::new(&edges, &repeat_target, 3);
        let good = search(0, &graph, &good_target, &all_targets, &good_graph, quota);
        let repeated = search(
            0,
            &graph,
            &repeat_target,
            &all_targets,
            &repeat_graph,
            quota,
        );
        assert!(!good.limited);
        assert_eq!(good.bridges, BTreeSet::from([(0, 1, b"ACGT".to_vec())]));
        assert!(repeated.limited);
        assert!(repeated.limit_reasons.contains("target_candidate_budget"));
        assert!(good.states + repeated.states <= 50);
        assert!(good.bridges.len() + repeated.bridges.len() <= 4);
        assert!(repeated.bridges.iter().all(|(_, to, _)| *to == 2));

        // Isolation preserves the completed proposal, not an unsupported
        // uniqueness claim. A different target's limit still vetoes its join.
        let bridge = good.bridges.first().unwrap().clone();
        let reciprocal = (4, 3, kmer::reverse_complement(&bridge.2));
        let reverse_segments = vec![vec![3], vec![4], vec![5], vec![0], vec![1], vec![2]];
        let proposed = BTreeMap::from([
            (bridge.clone(), BTreeSet::from([47])),
            (reciprocal, BTreeSet::from([47])),
        ]);
        // Discovery at k=47 cannot erase a source limit at another k.
        let rejected = admission(&bridge, &proposed, &BTreeSet::from([0]), &reverse_segments);
        assert!(!rejected.eligible);
        assert_eq!(rejected.reason, "source_competing_search_incomplete");
        let rejected = admission(&bridge, &proposed, &BTreeSet::from([4]), &reverse_segments);
        assert!(!rejected.eligible);
        assert_eq!(rejected.reason, "reciprocal_competing_search_incomplete");
        assert!(admission(&bridge, &proposed, &BTreeSet::new(), &reverse_segments).eligible);
        let without_reciprocal = BTreeMap::from([(bridge.clone(), BTreeSet::from([47]))]);
        let unknown = admission(
            &bridge,
            &without_reciprocal,
            &BTreeSet::from([4]),
            &reverse_segments,
        );
        assert_eq!(unknown.reason, "reciprocal_competing_search_incomplete");
        let missing = admission(
            &bridge,
            &without_reciprocal,
            &BTreeSet::new(),
            &reverse_segments,
        );
        assert_eq!(missing.reason, "reciprocal_spelling_missing");
    }

    #[test]
    fn target_isolation_still_stops_at_the_first_full_endpoint() {
        let graph = segments(&[b"AAC", b"GTT", b"TTG"]);
        let edges = small_edges(&[b"AACGTT"], 3);
        let target = HashMap::from([(code(b"TT"), vec![2])]);
        let all_targets = HashMap::from([(code(b"GT"), vec![1]), (code(b"TT"), vec![2])]);
        let small_graph = SearchGraph::new(&edges, &target, 3);
        let result = search(
            0,
            &graph,
            &target,
            &all_targets,
            &small_graph,
            limits(3, 2, 12),
        );
        assert!(result.bridges.is_empty());
        assert!(!result.limited);
        // Without the intervening endpoint, exactly the same graph has a path.
        let without_other_endpoint =
            search(0, &graph, &target, &target, &small_graph, limits(3, 2, 12));
        assert_eq!(
            without_other_endpoint.bridges,
            BTreeSet::from([(0, 2, b"ACGTT".to_vec())])
        );
    }

    #[test]
    fn short_overlapping_bridge_uses_the_existing_source_flank_context() {
        let graph = segments(&[b"AAAC", b"AACG"]);
        let edges = small_edges(&[b"AACG"], 3);
        let target = HashMap::from([(code(b"CG"), vec![1])]);
        let small_graph = SearchGraph::new(&edges, &target, 3);
        let result = search(0, &graph, &target, &target, &small_graph, limits(3, 4, 12));
        assert!(!result.limited);
        assert_eq!(result.bridges, BTreeSet::from([(0, 1, b"AAACG".to_vec())]));
    }

    #[test]
    fn reverse_complement_counterpart_is_a_role_not_a_prohibition() {
        let graph = segments(&[b"AAC", b"GTT"]);
        let edges = small_edges(&[b"AACGT"], 3);
        let target = HashMap::from([(code(b"GT"), vec![1])]);
        let small_graph = SearchGraph::new(&edges, &target, 3);
        let result = search(0, &graph, &target, &target, &small_graph, limits(3, 2, 12));
        let bridge = result.bridges.first().unwrap();
        let reverse_segments = vec![vec![1], vec![0]];
        assert_eq!(
            connection_role(0, 1, &reverse_segments, false),
            "same_end_reverse_complement"
        );
        assert_eq!(bridge.2, kmer::reverse_complement(&bridge.2));
        let proposed = BTreeMap::from([(bridge.clone(), BTreeSet::from([31]))]);
        assert!(admission(bridge, &proposed, &BTreeSet::new(), &reverse_segments).eligible);
    }

    #[test]
    fn lower_k_connects_only_observed_sequence_and_marks_an_unsearched_tail() {
        let mut seed = 617113u64;
        let sequence: Vec<_> = (0..220)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed % 4) as usize]
            })
            .collect();
        let graph = PathGraph {
            segments: vec![
                Segment {
                    sequence: sequence[..90].to_vec(),
                    circular: false,
                },
                Segment {
                    sequence: sequence[130..].to_vec(),
                    circular: false,
                },
            ],
            outgoing: vec![vec![], vec![]],
        };
        let k = 31;
        let mask = (1u128 << (2 * (k - 1))) - 1;
        let mut outgoing = HashMap::<_, Vec<_>>::new();
        kmer::visit(&sequence, None, k, 0, |_, c| {
            outgoing
                .entry(c >> 2)
                .or_default()
                .push((c & mask, b"ACGT"[(c & 3) as usize]));
        });
        let targets = HashMap::from([(code(&sequence[162..192]), vec![1])]);
        let small_graph = SearchGraph::new(&outgoing, &targets, k);
        let joined = search(
            0,
            &graph,
            &targets,
            &targets,
            &small_graph,
            limits(k, 62, 200),
        );
        assert!(!joined.limited);
        assert_eq!(
            joined.bridges,
            BTreeSet::from([(0, 1, sequence[28..192].to_vec())])
        );
        let capped = search(
            0,
            &graph,
            &targets,
            &targets,
            &small_graph,
            limits(k, 62, 40),
        );
        assert!(capped.limited);
        assert!(capped.bridges.is_empty());
        outgoing.remove(&code(&sequence[100..130]));
        let small_graph = SearchGraph::new(&outgoing, &targets, k);
        let absent = search(
            0,
            &graph,
            &targets,
            &targets,
            &small_graph,
            limits(k, 62, 200),
        );
        assert!(!absent.limited);
        assert!(absent.bridges.is_empty());
        assert_eq!(absent.states, 0);
    }

    #[test]
    fn targetless_repeat_cannot_exhaust_a_valid_connector_search() {
        let graph = PathGraph {
            segments: vec![
                Segment {
                    sequence: b"AAC".to_vec(),
                    circular: false,
                },
                Segment {
                    sequence: b"GTT".to_vec(),
                    circular: false,
                },
            ],
            outgoing: vec![vec![], vec![]],
        };
        // AAC -> ACG -> CGT reaches the target; ACC -> CCC -> CCC
        // is a targetless repeat which used to invalidate the whole source.
        let mut edges = HashMap::<_, Vec<_>>::new();
        for sequence in [b"AACGT".as_slice(), b"ACCCC"] {
            kmer::visit(sequence, None, 3, 0, |_, c| {
                let edge = (c & 15, b"ACGT"[(c & 3) as usize]);
                let next = edges.entry(c >> 2).or_default();
                if !next.contains(&edge) {
                    next.push(edge);
                }
            });
        }
        let targets = HashMap::from([(code(b"GT"), vec![1])]);
        let small_graph = SearchGraph::new(&edges, &targets, 3);
        let result = search(
            0,
            &graph,
            &targets,
            &targets,
            &small_graph,
            limits(3, 2, 12),
        );
        assert!(!result.limited);
        assert_eq!(result.bridges, BTreeSet::from([(0, 1, b"ACGT".to_vec())]));
        assert_eq!(result.unreachable_branches_pruned, 1);

        // If the repeat gains an exit to the target, its unsearched copy tail
        // must still block completeness. Reachability must not select one copy.
        edges
            .entry(code(b"CC"))
            .or_default()
            .push((code(b"CG"), b'G'));
        let small_graph = SearchGraph::new(&edges, &targets, 3);
        let repeat = search(
            0,
            &graph,
            &targets,
            &targets,
            &small_graph,
            limits(3, 2, 12),
        );
        assert!(repeat.limited);
        assert!(repeat.bridges.len() > 1);
    }

    #[test]
    fn reachable_suffix_does_not_replace_full_flank_matching() {
        let graph = PathGraph {
            segments: vec![
                Segment {
                    sequence: b"AAAC".to_vec(),
                    circular: false,
                },
                Segment {
                    sequence: b"TTGT".to_vec(),
                    circular: false,
                },
            ],
            outgoing: vec![vec![], vec![]],
        };
        let mut edges = HashMap::<_, Vec<_>>::new();
        kmer::visit(b"AACGT", None, 3, 0, |_, c| {
            edges
                .entry(c >> 2)
                .or_default()
                .push((c & 15, b"ACGT"[(c & 3) as usize]));
        });
        let targets = HashMap::from([(code(b"GT"), vec![1])]);
        let small_graph = SearchGraph::new(&edges, &targets, 3);
        let result = search(
            0,
            &graph,
            &targets,
            &targets,
            &small_graph,
            limits(3, 4, 12),
        );
        assert!(result.bridges.is_empty());
        assert!(!result.limited);
    }
}
