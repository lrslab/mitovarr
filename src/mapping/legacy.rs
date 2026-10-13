//! Best-effort gap filling. Preserve called sequence and emit one maximum
//! usable consensus; structural validation is reported separately.
use super::{graph::Graph, scoring, Alternative, Patch};
use crate::{
    cli::{MappingArgs, Reads},
    denovo::{pairs::Library, shortlist},
    io::fasta::write_fasta_records,
};
use anyhow::Result;
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    io::{BufWriter, Write},
};

pub(super) const KS: [usize; 4] = [21, 17, 13, 11];
const SCORED: usize = 2;

fn preserves_calls(sequence: &[u8], baseline: &[u8]) -> bool {
    sequence.len() == baseline.len()
        && sequence
            .iter()
            .zip(baseline)
            .all(|(&a, &b)| b == b'N' || a == b)
}

struct Proposal {
    sequence: Vec<u8>,
    support: u32,
    k: usize,
    sources: Vec<Value>,
}

pub(super) struct Inferred {
    pub length_preferred: Option<Patch>,
    pub read_preferred: Option<Patch>,
    pub report: Value,
}

/// Fast mode only fills N at the existing span. Reuse the established masked
/// walk; do not enumerate alternative lengths or run paired-fragment scoring.
#[allow(clippy::too_many_arguments)]
pub(super) fn fill(
    args: &MappingArgs,
    reads: &Reads,
    backbone: &[u8],
    start: usize,
    end: usize,
    locus: usize,
    graphs: &mut HashMap<usize, Graph>,
) -> Result<Inferred> {
    let baseline = &backbone[start..end];
    let mut report = json!({"method":"fixed_span_gap_fill","reference_guided":true,
        "target_length":baseline.len(),"target_source":"existing pileup backbone span; reference-guided length prior",
        "called_bases_and_span_preserved":true,"reference_bases_inserted":false,
        "search_complete":false,"final_verified":false,"strict_replacement":false,
        "read_preferred_rank":null,"paired_comparison":{"status":"not_evaluated",
            "scoring_status":"not_run_fast_consensus","best_candidate":null,
            "common_fragments":null,"unscored_fragments":null,"candidates":null}});
    if !baseline.contains(&b'N') {
        report["status"] = json!("not_needed_no_masked_bases");
        return Ok(Inferred {
            length_preferred: None,
            read_preferred: None,
            report,
        });
    }
    let mut pool = BTreeMap::<Vec<u8>, Proposal>::new();
    let mut searches = Vec::new();
    for k in KS.into_iter().filter(|&k| k <= args.local_k) {
        if let std::collections::hash_map::Entry::Vacant(entry) = graphs.entry(k) {
            entry.insert(Graph::build(
                reads,
                k,
                args.min_baseq,
                args.max_records,
                args.local_max_kmers,
            )?);
        }
        let (filled, search) = graphs[&k].fill_masked_span(baseline, 2);
        searches.push(search);
        if let Some((sequence, support)) = filled {
            pool.entry(sequence.clone())
                .or_insert(Proposal {
                    sequence,
                    support,
                    k,
                    sources: Vec::new(),
                })
                .sources
                .push(json!({"method":"maximum_support_masked_span","k":k,
                "minimum_pair_support":2,"path_minimum_pair_support":support}));
        }
    }
    report["searches"] = json!(searches);
    let mut pool: Vec<_> = pool.into_values().collect();
    pool.sort_by(|a, b| {
        (b.support >= args.min_depth)
            .cmp(&(a.support >= args.min_depth))
            .then(b.k.cmp(&a.k))
            .then(b.support.cmp(&a.support))
            .then(a.sequence.cmp(&b.sequence))
    });
    report["candidate_count"] = json!(pool.len());
    report["ranked_candidates"] = json!(pool
        .iter()
        .enumerate()
        .map(|(i, p)| json!({
            "rank":i+1,"length":p.sequence.len(),"k":p.k,"minimum_kmer_fragments":p.support,
            "called_bases_and_span_preserved":true,"sources":p.sources,"paired_candidate":null
        }))
        .collect::<Vec<_>>());
    if pool.is_empty() {
        report["status"] = json!("no_fill_candidate");
        return Ok(Inferred {
            length_preferred: None,
            read_preferred: None,
            report,
        });
    }
    let file = format!("local_assembly/locus_{locus}.gap_fill.fasta");
    write_fasta_records(
        &args.out.join(&file),
        &pool
            .iter()
            .enumerate()
            .map(|(i, p)| {
                (
                    format!(
                        "locus_{locus}_fill_{} k={} final_verified=false",
                        i + 1,
                        p.k
                    ),
                    p.sequence.clone(),
                )
            })
            .collect::<Vec<_>>(),
    )?;
    report["status"] = json!("inferred_fill");
    report["candidates_file"] = json!(file);
    report["maximum_consensus_rank"] = json!(1);
    report["maximum_consensus_selection"] = json!("preserve called bases and span; prefer supported longer k then support; structural scoring deferred");
    Ok(Inferred {
        length_preferred: Some(Patch {
            start,
            end,
            locus,
            sequence: pool.remove(0).sequence,
        }),
        read_preferred: None,
        report,
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn infer(
    args: &MappingArgs,
    reads: &Reads,
    backbone: &[u8],
    start: usize,
    end: usize,
    target_len: usize,
    locus: usize,
    main: &Graph,
    graphs: &mut HashMap<usize, Graph>,
    strict: &[Alternative],
    strict_origins: &HashMap<Vec<u8>, Vec<usize>>,
    library: Option<&Library>,
) -> Result<Inferred> {
    let left = &backbone[start..start + main.k];
    let right = &backbone[end - main.k..end];
    let mut pool = BTreeMap::<Vec<u8>, Proposal>::new();
    for p in strict.iter().skip(1) {
        let ks = &strict_origins[&p.sequence];
        pool.insert(
            p.sequence.clone(),
            Proposal {
                sequence: p.sequence.clone(),
                support: p.support,
                k: *ks.iter().max().unwrap(),
                sources: ks
                    .iter()
                    .map(|&k| json!({"method":"strict_search_candidate","k":k}))
                    .collect(),
            },
        );
    }
    let mut thresholds = vec![
        args.min_depth.max(2),
        args.min_depth.div_ceil(2).max(2),
        args.min_depth.div_ceil(4).max(2),
        3,
        2,
    ];
    thresholds.retain(|&n| n <= args.min_depth.max(2));
    thresholds.sort_unstable_by(|a, b| b.cmp(a));
    thresholds.dedup();
    let mut searches = Vec::new();
    for k in KS.into_iter().filter(|&k| k <= main.k) {
        if let std::collections::hash_map::Entry::Vacant(entry) = graphs.entry(k) {
            entry.insert(Graph::build(
                reads,
                k,
                args.min_baseq,
                args.max_records,
                args.local_max_kmers,
            )?);
        }
        let (filled, fill_report) = graphs[&k].fill_masked_span(&backbone[start..end], 2);
        searches.push(fill_report);
        if let Some((sequence, support)) = filled {
            let p = pool.entry(sequence.clone()).or_insert(Proposal {
                sequence,
                support,
                k,
                sources: vec![],
            });
            p.sources
                .push(json!({"method":"maximum_support_masked_span","k":k,
                "minimum_pair_support":2,"path_minimum_pair_support":support,
                "called_bases_preserved":true}));
        }
        for &minimum in &thresholds {
            let (paths, report) =
                graphs[&k].legacy_paths(left, right, target_len, args.local_max_span, minimum);
            let found = !paths.is_empty();
            searches.push(report);
            for (sequence, support) in paths {
                if sequence == backbone[start..end] {
                    continue;
                }
                let p = pool.entry(sequence.clone()).or_insert(Proposal {
                    sequence,
                    support,
                    k,
                    sources: vec![],
                });
                p.sources
                    .push(json!({"method":"legacy_reference_span_beam","k":k,
                    "minimum_pair_support":minimum,"path_minimum_pair_support":support}));
            }
            // The original search relaxed support only when no bridge existed.
            // Other k values still run, and their alternatives remain inspectable.
            if found {
                break;
            }
        }
    }
    if pool.is_empty() {
        return Ok(Inferred {
            length_preferred: None,
            read_preferred: None,
            report: json!({"status":"no_bridge_candidate","reference_guided":true,
                "target_length":target_len,"target_source":"reference interval between sample flanks",
                "search_complete":false,"final_verified":false,"searches":searches}),
        });
    }
    let pool: Vec<_> = pool.into_values().collect();
    let mut order: Vec<_> = (0..pool.len()).collect();
    order.sort_by(|&a, &b| {
        let a_preserves = preserves_calls(&pool[a].sequence, &backbone[start..end]);
        let b_preserves = preserves_calls(&pool[b].sequence, &backbone[start..end]);
        b_preserves
            .cmp(&a_preserves)
            .then_with(|| {
                if a_preserves && b_preserves {
                    // At the same complete span, prefer supported longer context.
                    (pool[b].support >= args.min_depth)
                        .cmp(&(pool[a].support >= args.min_depth))
                        .then(pool[b].k.cmp(&pool[a].k))
                        .then(pool[b].support.cmp(&pool[a].support))
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .then_with(|| {
                pool[a]
                    .sequence
                    .len()
                    .abs_diff(target_len)
                    .cmp(&pool[b].sequence.len().abs_diff(target_len))
            })
            .then(pool[b].support.cmp(&pool[a].support))
            .then(pool[b].k.cmp(&pool[a].k))
            .then(pool[a].sequence.cmp(&pool[b].sequence))
    });
    let file = format!("local_assembly/locus_{locus}.legacy_bridge.fasta");
    write_fasta_records(&args.out.join(&file), &order.iter().enumerate().map(|(rank,&i)|
        (format!("locus_{locus}_legacy_{} k={} target_length={} reference_guided=true final_verified=false",rank+1,pool[i].k,target_len),pool[i].sequence.clone())).collect::<Vec<_>>())?;
    let shortlist = shortlist::select(
        &pool
            .iter()
            .map(|p| shortlist::Candidate {
                sequence: &p.sequence,
                k: p.k,
                circular: false,
            })
            .collect::<Vec<_>>(),
        &order,
        SCORED,
    );
    let mut alternatives = vec![Alternative {
        sequence: backbone[start..end].to_vec(),
        support: super::path_support(main, &backbone[start..end]),
    }];
    for &i in &shortlist.indices {
        alternatives.push(Alternative {
            sequence: pool[i].sequence.clone(),
            support: pool[i].support,
        });
    }
    tracing::info!(
        locus,
        candidates = pool.len(),
        scored = alternatives.len(),
        "score legacy bridge hypotheses"
    );
    let mut evaluation = scoring::evaluate(
        backbone,
        start,
        end,
        &alternatives,
        reads,
        library,
        main.max_read,
        args,
        Some(&args.out.join(format!(
            "local_assembly/locus_{locus}.legacy_bridge.fragments.jsonl.gz"
        ))),
        Some(
            &pool
                .iter()
                .map(|p| p.sequence.as_slice())
                .collect::<Vec<_>>(),
        ),
    )?
    .report;
    let best = evaluation["best_candidate"].as_u64().unwrap_or(0) as usize;
    let unique_measured_best = evaluation["best_candidates"]
        .as_array()
        .is_some_and(|v| v.len() == 1)
        && evaluation["common_fragments"].as_u64().unwrap_or(0) > 0;
    let read_preferred = (unique_measured_best && best > 0).then(|| Patch {
        start,
        end,
        locus,
        sequence: alternatives[best].sequence.clone(),
    });
    // Never expose a bounded, reference-targeted comparison as a resolved locus.
    evaluation["within_shortlist_status"] = evaluation["status"].clone();
    evaluation["within_shortlist_selected_candidate"] = evaluation["selected_candidate"].clone();
    evaluation["status"] = json!("diagnostic_only");
    evaluation["selected_candidate"] = Value::Null;
    evaluation["final_verified"] = json!(false);
    let ranking:Vec<_>=order.iter().enumerate().map(|(rank,&i)|json!({"rank":rank+1,
        "length":pool[i].sequence.len(),"target_length_difference":pool[i].sequence.len().abs_diff(target_len),
        "called_bases_and_span_preserved":preserves_calls(&pool[i].sequence,&backbone[start..end]),
        "sources":pool[i].sources,"paired_candidate":shortlist.indices.iter().position(|&j|j==i).map(|j|j+1)})).collect();
    let preferred_rank = (unique_measured_best && best > 0).then(|| {
        order
            .iter()
            .position(|&i| i == shortlist.indices[best - 1])
            .unwrap()
            + 1
    });
    Ok(Inferred {
        length_preferred: Some(Patch {
            start,
            end,
            locus,
            sequence: pool[order[0]].sequence.clone(),
        }),
        read_preferred,
        report: json!({"status":"inferred_bridge","reference_guided":true,"target_length":target_len,
            "target_source":"reference interval between sample flanks; length is a proposal heuristic",
            "candidate_count":pool.len(),"candidates_file":file,"length_preferred_rank":1,
            "maximum_consensus_rank":if backbone[start..end].contains(&b'N') && preserves_calls(&pool[order[0]].sequence,&backbone[start..end]){Some(1)}else{None},
            "maximum_consensus_selection":"fill N only at the existing span and preserve every called base; otherwise retain the backbone; other hypotheses remain diagnostic",
            "read_preferred_rank":preferred_rank,"read_preference_scope":"unique measured best within shortlist; unknown fragments and unsearched paths remain unresolved",
            "shortlist":shortlist.report,"ranked_candidates":ranking,"paired_comparison":evaluation,
            "searches":searches,"search_complete":false,"final_verified":false,
            "strict_replacement":false,"reference_bases_inserted":false,
            "limitations":["legacy strategy port, not a replay of the July binary","reference target length can bias repeat copy number","incidental short-k matches can propose false joins across real read gaps","complete sample flanks retained at every k","minimum support can relax to two named pairs","bounded candidates and read preference do not establish structural uniqueness"]}),
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn write_hypotheses(
    args: &MappingArgs,
    backbone: &[u8],
    cut: usize,
    contig: &str,
    strict_validated: bool,
    strict: &[Patch],
    length: &[Patch],
    preferred: &[Patch],
    intervals: &[(usize, usize)],
    reports: &[Value],
) -> Result<(Value, Vec<u8>)> {
    let mut outputs = Vec::new();
    let mut maximum_sequence = Vec::new();
    let maximum: Vec<_> = length
        .iter()
        .filter(|p| {
            backbone[p.start..p.end].contains(&b'N')
                && preserves_calls(&p.sequence, &backbone[p.start..p.end])
        })
        .cloned()
        .collect();
    for (label, proposed) in [
        ("maximum_consensus", maximum.as_slice()),
        ("legacy_bridge", length),
        ("legacy_bridge.read_preferred", preferred),
    ] {
        if proposed.is_empty() && label != "maximum_consensus" {
            continue;
        }
        let mut patches: Vec<_> = strict.iter().chain(proposed).collect();
        patches.sort_by_key(|p| p.start);
        let mut sequence = Vec::new();
        let mut cursor = 0;
        for p in &patches {
            sequence.extend_from_slice(&backbone[cursor..p.start]);
            sequence.extend_from_slice(&p.sequence);
            cursor = p.end;
        }
        sequence.extend_from_slice(&backbone[cursor..]);
        let file = format!("{label}.fasta");
        let name = if label == "maximum_consensus" {
            contig
        } else {
            label
        };
        let verified = label == "maximum_consensus" && strict_validated;
        write_fasta_records(&args.out.join(&file),&[(format!("{name} consensus_mode=best_effort reference_guided=true final_verified={verified} circularity_verified=false topology=reference_assumed"),sequence.clone())])?;
        let project = |position: usize, right: bool| {
            let mut shift = 0isize;
            for p in &patches {
                if p.start < position && position < p.end {
                    return (p.start as isize + shift) as usize
                        + if right { p.sequence.len() } else { 0 };
                }
                if p.end <= position {
                    shift += p.sequence.len() as isize - (p.end - p.start) as isize;
                }
            }
            (position as isize + shift) as usize
        };
        let uncertainty_file = format!("{label}.ambiguity.tsv");
        let mut uncertainty = BufWriter::new(File::create(args.out.join(&uncertainty_file))?);
        writeln!(
            uncertainty,
            "sequence\tstart\tend\tevidence\treason\talternatives"
        )?;
        for (i, &(start, end)) in intervals.iter().enumerate() {
            if reports[i]["status"] == "unresolved" {
                writeln!(
                    uncertainty,
                    "{name}\t{}\t{}\tlocal_locus_{}\t{}\t{}",
                    project(start, false) + 1,
                    project(end, true),
                    i + 1,
                    if args.deep_validation {
                        "legacy_bridge_structure_unverified"
                    } else {
                        "structural_validation_deferred"
                    },
                    reports[i]["legacy_bridge"]["candidate_count"]
                        .as_u64()
                        .unwrap_or(0)
                )?;
            }
        }
        uncertainty.flush()?;
        outputs.push(json!({"file":file,"ambiguity_file":uncertainty_file,"length":sequence.len(),
            "ambiguous_bases":sequence.iter().filter(|&&b|b==b'N').count(),"inferred_replacements":proposed.len(),"final_verified":verified}));
        if label == "maximum_consensus" {
            maximum_sequence = sequence;
        }
    }
    let report = json!({"enabled":true,"status":if maximum.is_empty(){"no_inferred_windows"}else{"best_effort_consensus"},
        "k_values":KS,"outputs":outputs,"final_verified":false,"search_complete":false,
        "maximum_consensus":"maximum_consensus.fasta","consensus_mode":"best_effort",
        "maximum_consensus_rejected_patches":length.iter().filter(|p|backbone[p.start..p.end].contains(&b'N') && !preserves_calls(&p.sequence,&backbone[p.start..p.end])).map(|p|json!({"locus":p.locus,"reason":"would_change_called_bases_or_span","baseline_length":p.end-p.start,"proposed_length":p.sequence.len()})).collect::<Vec<_>>(),
        "internal_rotation_baseline_offset":cut,"coordinate_system":"FASTA and ambiguity TSV share the internal local-assembly rotation; not reference coordinates",
        "cross_window_phase_verified":false,"strict_sequence_changed_by_legacy_bridge":false});
    crate::write_json(&args.out.join("legacy_bridge.json"), &report)?;
    Ok((report, maximum_sequence))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Commands};
    use clap::Parser;

    #[test]
    fn maximum_output_rejects_called_base_changes_and_length_changes() {
        let root =
            std::env::temp_dir().join(format!("mitovarr-protected-maximum-{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let backbone = b"GGACNNTGCC";
        for (i, proposed) in [b"ATCCTG".as_slice(), b"ACCTG", b"ACCCCTG", b"ACCCTG"]
            .into_iter()
            .enumerate()
        {
            let out = root.join(i.to_string());
            std::fs::create_dir(&out).unwrap();
            let cli = Cli::parse_from([
                "mitovarr",
                "assemble-mapping",
                "--ref",
                "unused",
                "--bam",
                "unused",
                "--out",
                out.to_str().unwrap(),
                "--local-assembly",
                "--r1",
                "unused",
                "--r2",
                "unused",
                "--deep-validation",
                "--best-effort",
            ]);
            let Commands::AssembleMapping(args) = cli.command else {
                unreachable!()
            };
            let patch = Patch {
                start: 2,
                end: 8,
                locus: 1,
                sequence: proposed.to_vec(),
            };
            let (report, maximum) = write_hypotheses(
                &args,
                backbone,
                0,
                "MtDNA",
                false,
                &[],
                &[patch],
                &[],
                &[(2, 8)],
                &[json!({"status":"unresolved","legacy_bridge":{"candidate_count":1}})],
            )
            .unwrap();
            let expected = if i == 3 {
                b"GGACCCTGCC".as_slice()
            } else {
                backbone
            };
            assert_eq!(maximum, expected);
            assert_eq!(
                crate::io::fasta::load_reference(&out.join("maximum_consensus.fasta"), None)
                    .unwrap()
                    .sequence,
                expected
            );
            assert_eq!(
                report["maximum_consensus_rejected_patches"]
                    .as_array()
                    .unwrap()
                    .len(),
                usize::from(i != 3)
            );
            assert_eq!(report["outputs"][0]["final_verified"], false);
            // Alternative structures stay explicitly diagnostic and available.
            assert_eq!(
                crate::io::fasta::load_reference(&out.join("legacy_bridge.fasta"), None)
                    .unwrap()
                    .sequence,
                [b"GG".as_slice(), proposed, b"CC"].concat()
            );
        }
        std::fs::remove_dir_all(root).unwrap();
    }
}
