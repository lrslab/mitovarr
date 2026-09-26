//! Reference-guided local assembly. No reference enters the de novo workflow.
mod discovery;
mod graph;
mod scoring;
pub(crate) use discovery::Signals;

use crate::{
    cli::{GraphSettings, Limits, PolishArgs, Reads},
    consensus::ConsensusResult,
    denovo::pairs,
    io::fasta::write_consensus_fasta,
    kmer,
    repeat::AmbiguityInterval,
};
use anyhow::{ensure, Context, Result};
use graph::Graph;
use serde_json::{json, Value};
use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufWriter, Write},
};

struct Alternative {
    sequence: Vec<u8>,
    support: u32,
}
struct Patch {
    start: usize,
    end: usize,
    sequence: Vec<u8>,
    locus: usize,
}

pub(crate) struct ResultAssembly {
    pub sequence: Vec<u8>,
    pub resolved: bool,
    pub unresolved_ambiguity: usize,
    pub report: Value,
}

fn rotate(sequence: &[u8], cut: usize) -> Vec<u8> {
    [&sequence[cut..], &sequence[..cut]].concat()
}

fn anchor_interval(
    sequence: &[u8],
    interval: (usize, usize),
    graph: &Graph,
    minimum: u32,
) -> Option<(usize, usize)> {
    let k = graph.k;
    let mut occurrences = HashMap::<kmer::Code, usize>::new();
    kmer::visit(sequence, None, k, 0, |_, c| {
        *occurrences.entry(kmer::canonical(c, k)).or_default() += 1
    });
    let valid = |p: usize| {
        let s = &sequence[p..p + k];
        graph::code(s).is_some_and(|c| occurrences.get(&kmer::canonical(c, k)) == Some(&1))
            && graph.support(s) >= minimum
    };
    let left_end = interval.0.checked_sub(k)?;
    let left = (left_end.saturating_sub(350)..=left_end)
        .rev()
        .find(|&p| valid(p))?;
    let right = (interval.1..=(interval.1 + 350).min(sequence.len().checked_sub(k)?))
        .find(|&p| valid(p))?;
    Some((left, right + k))
}

fn path_support(graph: &Graph, sequence: &[u8]) -> u32 {
    if sequence.contains(&b'N') || sequence.len() < graph.k {
        return 0;
    }
    sequence
        .windows(graph.k)
        .map(|s| graph.support(s))
        .min()
        .unwrap_or(0)
}

fn ref_position(called: &ConsensusResult, baseline_position: usize) -> usize {
    called
        .reference_offsets
        .partition_point(|&offset| offset <= baseline_position)
        .saturating_sub(1)
        + 1
}

/// Split a zero-based interval at an output/circular origin for JSON consumers.
fn segments(start: usize, len: usize, total: usize) -> Vec<Value> {
    if len == 0 || total == 0 {
        return vec![];
    }
    let start = start % total;
    let end = start + len;
    let mut parts = vec![json!({"start":start+1,"end":end.min(total)})];
    if end > total {
        parts.push(json!({"start":1,"end":end-total}));
    }
    parts
}

pub(crate) fn assemble(
    args: &PolishArgs,
    called: &ConsensusResult,
    ambiguity: &[AmbiguityInterval],
    signals: Signals,
    contig: &str,
) -> Result<ResultAssembly> {
    ensure!(
        (15..=31).contains(&args.local_k)
            && args.local_max_span > 2 * args.local_k
            && args.local_max_paths > 1,
        "local k must be 15..31, local span > 2*k and local path limit > 1"
    );
    ensure!(
        !called.sequence.is_empty(),
        "local assembly requires a nonempty pileup backbone"
    );
    let reads = Reads {
        r1: args
            .r1
            .clone()
            .context("local assembly requires paired R1 FASTQ")?,
        r2: args
            .r2
            .clone()
            .context("local assembly requires paired R2 FASTQ")?,
    };
    let directory = args.out.join("local_assembly");
    fs::create_dir_all(&directory)?;
    write_consensus_fasta(
        &args.out.join("reference_candidate.fasta"),
        contig,
        &called.sequence,
    )?;
    let (mask, reasons) = signals.windows(called, ambiguity);
    let (cut, intervals) = discovery::rotate_windows(&mask);
    let backbone = rotate(&called.sequence, cut);
    let graph = Graph::build(&reads, args.local_k, args.min_baseq, args.max_records)?;
    let settings = GraphSettings {
        k: args.local_k,
        min_kmer_count: 3,
        max_error_fraction: 0.2,
        max_bubble_fraction: 0.49,
        min_depth: args.min_depth,
        min_junction_fragments: 3,
        max_validation_pairs: args.max_records,
    };
    let limits = Limits {
        max_kmers: 2_000_000,
        min_baseq: args.min_baseq,
    };
    let mut training_sequence = backbone.clone();
    for &(start, end) in &intervals {
        training_sequence[start..end].fill(b'N');
    }
    let training_sequences = training_sequence
        .split(|b| *b == b'N')
        .filter(|s| s.len() >= 31)
        .map(|s| s.to_vec())
        .collect();
    let training = pairs::estimate_sequences(training_sequences, &reads, &settings, &limits)?;
    let mut gfa = BufWriter::new(File::create(args.out.join("local_assembly.gfa"))?);
    writeln!(gfa, "H\tVN:Z:1.0")?;
    let mut events = BufWriter::new(File::create(args.out.join("structural_variants.tsv"))?);
    writeln!(events,"locus\tref_window_start\tref_window_end\twraps_origin\tbaseline_length\tselected_length\tdelta\tstatus\tselected_candidate\treason")?;
    let mut reports = Vec::new();
    let mut fallback_graphs = HashMap::<usize, Graph>::new();
    let mut patches = Vec::<Patch>::new();
    let mut resolved_mask = vec![false; backbone.len()];
    let mut previous_end = 0;
    let mut evaluated_intervals = Vec::new();
    for (number, &interval) in intervals.iter().enumerate() {
        let locus = number + 1;
        evaluated_intervals.push(interval);
        tracing::info!(locus, regions = intervals.len(), "assemble mapping region");
        let Some((start, end)) = anchor_interval(&backbone, interval, &graph, args.min_depth)
        else {
            reports.push(json!({"locus":locus,"status":"unresolved","reason":"no_supported_unique_flanks","rotated_baseline_interval":[interval.0+1,interval.1]}));
            writeln!(
                events,
                "{locus}\t.\t.\t.\t{}\t.\t.\tunresolved\t.\tno_supported_unique_flanks",
                interval.1 - interval.0
            )?;
            continue;
        };
        if start < previous_end || end - start > args.local_max_span {
            reports.push(json!({"locus":locus,"status":"unresolved","reason":"overlapping_flanks_or_span_limit"}));
            writeln!(
                events,
                "{locus}\t.\t.\t.\t{}\t.\t.\tunresolved\t.\toverlapping_flanks_or_span_limit",
                end - start
            )?;
            continue;
        }
        previous_end = end;
        evaluated_intervals[number] = (start, end);
        let baseline = backbone[start..end].to_vec();
        let mut search = graph.paths(
            &baseline[..graph.k],
            &baseline[baseline.len() - graph.k..],
            args.local_max_span,
            args.local_max_paths,
            args.min_depth,
        );
        let mut proposals = vec![
            json!({"k":graph.k,"paths":search.paths.len(),"states":search.states,"limits":search.limits}),
        ];
        // Smaller k can bridge a coverage break. Keep the original full flanks,
        // pool all proposed spellings before scoring, and retain every limit.
        // No reference bases are supplied to a gap: sequences come from FASTQ.
        if search.paths.is_empty() && search.limits.is_empty() {
            for k in [21, 15].into_iter().filter(|&k| k < graph.k) {
                if let std::collections::hash_map::Entry::Vacant(entry) = fallback_graphs.entry(k) {
                    entry.insert(Graph::build(&reads, k, args.min_baseq, args.max_records)?);
                }
                let fallback = &fallback_graphs[&k];
                let mut extra = fallback.anchored_paths(
                    &baseline[..graph.k],
                    &baseline[baseline.len() - graph.k..],
                    args.local_max_span,
                    args.local_max_paths,
                    args.min_depth,
                );
                proposals.push(json!({"k":k,"paths":extra.paths.len(),"states":extra.states,"limits":extra.limits}));
                search.states += extra.states;
                search.low_support_edges_excluded += extra.low_support_edges_excluded;
                search.limits.extend(extra.limits);
                search.paths.append(&mut extra.paths);
            }
            search.paths.sort_by(|a, b| a.0.cmp(&b.0));
            search.paths.dedup_by(|a, b| a.0 == b.0);
            if search.paths.len() >= args.local_max_paths {
                search.limits.insert("pooled_path_limit");
                search.paths.truncate(args.local_max_paths);
            }
        }
        let mut alternatives = vec![Alternative {
            support: path_support(&graph, &baseline),
            sequence: baseline.clone(),
        }];
        for (sequence, support) in search.paths {
            if sequence != baseline {
                alternatives.push(Alternative { sequence, support });
            }
        }
        let baseline_start = (cut + start) % backbone.len();
        let baseline_end = (cut + end - 1) % backbone.len();
        let ref_start = ref_position(called, baseline_start);
        let ref_end = ref_position(called, baseline_end);
        let prefix = format!("locus_{locus}");
        let mut fasta = BufWriter::new(File::create(directory.join(format!("{prefix}.fasta")))?);
        writeln!(
            gfa,
            "S\t{prefix}_left\t{}",
            String::from_utf8_lossy(&baseline[..graph.k])
        )?;
        writeln!(
            gfa,
            "S\t{prefix}_right\t{}",
            String::from_utf8_lossy(&baseline[baseline.len() - graph.k..])
        )?;
        for (i, p) in alternatives.iter().enumerate() {
            writeln!(
                fasta,
                ">{prefix}_candidate_{i} baseline={}\n{}",
                i == 0,
                String::from_utf8_lossy(&p.sequence)
            )?;
            writeln!(
                gfa,
                "S\t{prefix}_candidate_{i}\t{}\tKC:i:{}",
                String::from_utf8_lossy(&p.sequence),
                p.support
            )?;
            writeln!(
                gfa,
                "L\t{prefix}_left\t+\t{prefix}_candidate_{i}\t+\t{}M",
                graph.k
            )?;
            writeln!(
                gfa,
                "L\t{prefix}_candidate_{i}\t+\t{prefix}_right\t+\t{}M",
                graph.k
            )?;
        }
        fasta.flush()?;
        let evaluation = if search.limits.is_empty() {
            scoring::evaluate(
                &backbone,
                start,
                end,
                &alternatives,
                &reads,
                training.model.as_ref(),
                graph.max_read,
                args,
            )?
        } else {
            // An incomplete search cannot select a replacement. Preserve its
            // paths, but do not spend minutes ranking an ineligible subset.
            // Null means unmeasured, not zero fragments or no read support.
            scoring::Evaluation {
                selected: None,
                report: json!({"status":"unresolved","reason":"graph_search_limit",
                    "scoring_status":"not_run_graph_search_limit",
                    "selected_candidate":null,"best_candidate":null,"common_fragments":null,
                    "shared_unique_anchor_kmers":null,"unscored_fragments":null,
                    "candidates":alternatives.iter().enumerate().map(|(i,p)|json!({
                        "candidate":i,"length":p.sequence.len(),"minimum_kmer_fragments":p.support,
                        "score":null,"junction_fragments":null,"strong_votes_against":null
                    })).collect::<Vec<_>>() }),
            }
        };
        if let Some(selected) = evaluation.selected {
            for p in start..end {
                resolved_mask[(cut + p) % backbone.len()] = true;
            }
            if selected != 0 {
                patches.push(Patch {
                    start,
                    end,
                    sequence: alternatives[selected].sequence.clone(),
                    locus,
                });
            }
        }
        let selected_length = evaluation.selected.map(|i| alternatives[i].sequence.len());
        let delta = selected_length.map(|n| n as isize - baseline.len() as isize);
        writeln!(
            events,
            "{locus}\t{ref_start}\t{ref_end}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            baseline_start > baseline_end,
            baseline.len(),
            selected_length.map_or(".".into(), |x| x.to_string()),
            delta.map_or(".".into(), |x| x.to_string()),
            evaluation.report["status"].as_str().unwrap(),
            evaluation.selected.map_or(".".into(), |x| x.to_string()),
            evaluation.report["reason"].as_str().unwrap()
        )?;
        let mut entry = evaluation.report;
        if let Some(selected) = evaluation.selected.filter(|&i| i != 0) {
            let alt = &alternatives[selected].sequence;
            let prefix = baseline.iter().zip(alt).take_while(|(a, b)| a == b).count();
            let suffix = baseline[prefix..]
                .iter()
                .rev()
                .zip(alt[prefix..].iter().rev())
                .take_while(|(a, b)| a == b)
                .count();
            let replaced = baseline.len() - prefix - suffix;
            let inserted = alt.len() - prefix - suffix;
            let left = (cut + start + prefix + backbone.len() - 1) % backbone.len();
            let right = (cut + start + baseline.len() - suffix) % backbone.len();
            entry["edit"] = json!({"relative_to":"reference_candidate.fasta (after pileup indels)",
                "kind":if replaced==0{"insertion"}else if inserted==0{"deletion"}else{"replacement"},
                "left_flank_reference_pos":ref_position(called,left),"right_flank_reference_pos":ref_position(called,right),
                "replaced_baseline_bases":replaced,"inserted_sample_bases":inserted,
                "baseline_segments":segments((left+1)%backbone.len(),replaced,backbone.len()),
                "replacement_sequence":String::from_utf8_lossy(&alt[prefix..alt.len()-suffix])});
        }
        entry["locus"] = json!(locus);
        entry["reference_interval"] =
            json!({"start":ref_start,"end":ref_end,"wraps_origin":baseline_start>baseline_end});
        entry["baseline_segments"] = json!(segments(baseline_start, end - start, backbone.len()));
        entry["candidates_fasta"] = json!(format!("local_assembly/{prefix}.fasta"));
        entry["graph_proposals"] = json!(proposals);
        entry["search"] = json!({"states":search.states,"limits":search.limits,"low_support_edges_excluded":search.low_support_edges_excluded});
        reports.push(entry);
    }
    gfa.flush()?;
    events.flush()?;
    let origin = (backbone.len() - cut) % backbone.len();
    let mut assembled = Vec::new();
    let mut blocks = Vec::new();
    let mut cursor = 0;
    let mut output_origin = None;
    let mut shifted_origin = false;
    for patch in &patches {
        if cursor <= origin && origin < patch.start {
            output_origin = Some(assembled.len() + origin - cursor);
        }
        assembled.extend_from_slice(&backbone[cursor..patch.start]);
        let output_start = assembled.len();
        if patch.start <= origin && origin < patch.end {
            // The old origin can lie in a deleted/replaced sequence. Use its
            // supported left anchor instead of inventing a coordinate inside it.
            output_origin = Some(output_start);
            shifted_origin = true;
        }
        assembled.extend_from_slice(&patch.sequence);
        blocks.push((patch.locus, output_start, patch.sequence.len()));
        cursor = patch.end;
    }
    if origin >= cursor {
        output_origin = Some(assembled.len() + origin - cursor);
    }
    assembled.extend_from_slice(&backbone[cursor..]);
    let output_origin = output_origin.unwrap_or(0);
    let sequence = rotate(&assembled, output_origin);
    for &(locus, start, length) in &blocks {
        reports[locus - 1]["output_segments"] = json!(segments(
            (start + sequence.len() - output_origin) % sequence.len(),
            length,
            sequence.len()
        ));
    }
    let mut uncertainty = BufWriter::new(File::create(args.out.join("structural_ambiguity.tsv"))?);
    writeln!(
        uncertainty,
        "sequence\tstart\tend\tevidence\treason\talternatives"
    )?;
    let project_boundary = |position: usize, right: bool| {
        let mut shift = 0isize;
        for patch in &patches {
            if patch.start < position && position < patch.end {
                return (patch.start as isize + shift) as usize
                    + if right { patch.sequence.len() } else { 0 };
            }
            if patch.end <= position {
                shift += patch.sequence.len() as isize - (patch.end - patch.start) as isize;
            }
        }
        (position as isize + shift) as usize
    };
    for (i, report) in reports.iter_mut().enumerate() {
        let (start, end) = evaluated_intervals[i];
        if report.get("output_segments").is_none() {
            report["output_segments"] = json!(segments(
                (project_boundary(start, false) + sequence.len() - output_origin) % sequence.len(),
                project_boundary(end, true) - project_boundary(start, false),
                sequence.len()
            ));
        }
        if report["status"] == "unresolved" {
            for segment in report["output_segments"].as_array().unwrap() {
                writeln!(
                    uncertainty,
                    "{contig}\t{}\t{}\tlocal_locus_{}\t{}\t{}",
                    segment["start"],
                    segment["end"],
                    i + 1,
                    report["reason"].as_str().unwrap(),
                    report["candidates"].as_array().map_or(0, Vec::len)
                )?;
            }
        }
    }
    uncertainty.flush()?;
    let unresolved_ambiguity = ambiguity
        .iter()
        .filter(|a| {
            let start = called.reference_offsets[a.start - 1];
            let end = called.reference_offsets[a.end];
            resolved_mask[start..end].iter().any(|resolved| !*resolved)
        })
        .count();
    let resolved = reports.iter().all(|r| r["status"] == "resolved");
    let report = json!({"schema":"mitovar-rust.mapping-local.v1","workflow":"reference_guided_local_assembly",
        "r1":reads.r1.canonicalize()?,"r2":reads.r2.canonicalize()?,"read_scope":"provided complete pairs; may be mt-bam-extracted, not assumed whole-sample",
        "discovery_signals":reasons,"named_pairs":graph.pairs,"sample_kmers":graph.counts.len(),"weak_graph_edges_removed":graph.weak_edges_removed,"fragment_model":training.report,
        "parameters":{"k":graph.k,"minimum_kmer_fragments":3,"minimum_path_kmer_fragments":args.min_depth,"minimum_junction_fragments":3,"breakpoint_cluster_width":7,"minimum_breakpoint_depth_fraction":0.05,"weak_edge_filter":"support < 0.02 of strongest outgoing edge, independent of absolute depth","max_span":args.local_max_span,"max_paths":args.local_max_paths,"max_states":200000,"max_visits_per_node":4,"max_kmers":2000000,"fallback_k":[21,15],"fallback_rule":"only when main graph has no path and no search limit; retain full main-k flanks; pool and jointly score proposals"},
        "windows":reports,"applied_replacements":patches.len(),"all_windows_resolved":resolved,"remaining_repeat_intervals":unresolved_ambiguity,
        "internal_rotation_baseline_offset":cut,"output_origin_shifted_to_supported_anchor":shifted_origin,
        "coordinate_system":"one-based inclusive; baseline_segments refer to reference_candidate.fasta; output_segments refer to candidate.fasta; pileup TSVs remain reference-based",
        "limitations":["local candidates use the reference/pileup backbone outside each window","candidate windows are scored independently; distant variants are not phased","no claim of whole-genome circularity or exhaustive SV discovery","input FASTQ may omit reads not recovered by the upstream BAM extraction"],
        "sequence_source":"sample FASTQ graph paths; no reference-filled bases"});
    crate::write_json(&directory.join("report.json"), &report)?;
    Ok(ResultAssembly {
        sequence,
        resolved,
        unresolved_ambiguity,
        report,
    })
}
