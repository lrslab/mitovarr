//! Optional mapping-layout hypotheses. Coverage shortlists graph paths; paired
//! scores compare the same fragments, but a bounded pool never grants finality.
use super::{graph::Graph, scoring, Alternative, Patch};
use crate::{
    cli::{MappingArgs, Reads},
    denovo::pairs::Library,
    io::fasta::write_fasta_records,
    kmer,
};
use anyhow::Result;
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fs::File,
    io::{BufWriter, Write},
};

const FLANK_WINDOW: usize = 500;
const MAX_SCORED_PATHS: usize = 8;
type Counts = BTreeMap<kmer::Code, usize>;

pub(super) fn occurrences(sequence: &[u8], circular: bool, k: usize) -> Counts {
    let mut extended = sequence.to_vec();
    if circular && !sequence.is_empty() {
        extended.extend(sequence.iter().copied().cycle().take(k - 1));
    }
    let mut counts = BTreeMap::new();
    kmer::visit(&extended, None, k, 0, |_, c| {
        *counts.entry(kmer::canonical(c, k)).or_default() += 1;
    });
    counts
}

pub(super) struct Depth<'a> {
    graph: &'a Graph,
    outside: Counts,
    pub mean: f64,
    report: Value,
}

impl<'a> Depth<'a> {
    fn new(graph: &'a Graph, backbone: &[u8], start: usize, end: usize) -> Option<Self> {
        let all = occurrences(backbone, true, graph.k);
        let baseline = |sequence: &[u8]| {
            let unique: Vec<_> = occurrences(sequence, false, graph.k)
                .into_keys()
                .filter(|c| all.get(c) == Some(&1))
                .collect();
            let sum: u64 = unique
                .iter()
                .map(|c| graph.r1_occurrences.get(c).copied().unwrap_or(0))
                .sum();
            (unique.len(), sum as f64 / unique.len().max(1) as f64)
        };
        // Flanks include the supported anchor, but never the inferred interior.
        let left = baseline(&backbone[start.saturating_sub(FLANK_WINDOW)..start + graph.k]);
        let right = baseline(&backbone[end - graph.k..(end + FLANK_WINDOW).min(backbone.len())]);
        let global = baseline(backbone);
        let n = left.0 + right.0;
        let local = n >= graph.k;
        let mean = if local {
            (left.0 as f64 * left.1 + right.0 as f64 * right.1) / n as f64
        } else {
            global.1
        };
        if !mean.is_finite() || mean <= 0.0 || global.0 < graph.k {
            return None;
        }
        let mut outside = all;
        // Full flanks remain unchanged, so boundary-crossing k-mers are fixed.
        // Subtract only complete local windows; outside + path is then the
        // exact circular whole-contig count, including other homologous loci.
        for (c, count) in occurrences(&backbone[start..end], false, graph.k) {
            *outside.get_mut(&c).unwrap() -= count;
        }
        let report = json!({"k":graph.k,"mean_r1_window_depth":mean,
            "baseline":if local{"unique_local_flanks"}else{"unique_backbone"},
            "left_unique_kmers":left.0,"left_mean":left.1,"right_unique_kmers":right.0,"right_mean":right.1,
            "backbone_unique_kmers":global.0,"backbone_mean":global.1,
            "counting":"R1 window occurrences including repeated offsets; one record per named pair; R2 excluded"});
        Some(Self {
            graph,
            outside,
            mean,
            report,
        })
    }

    fn expected(&self, c: kmer::Code) -> f64 {
        self.graph.r1_occurrences.get(&c).copied().unwrap_or(0) as f64 / self.mean
    }

    pub fn increment(&self, c: kmer::Code, used: usize, added: usize) -> f64 {
        let before = (self.outside.get(&c).copied().unwrap_or(0) + used) as f64;
        (before + added as f64 - self.expected(c)).powi(2) - (before - self.expected(c)).powi(2)
    }

    pub fn exploration_bound(&self, c: kmer::Code) -> usize {
        ((1.5 * self.expected(c)).ceil() as usize + 2)
            .saturating_sub(self.outside.get(&c).copied().unwrap_or(0))
            .clamp(2, 4096)
    }

    fn fit(&self, local: &Counts, domain: &BTreeSet<kmer::Code>, scale: f64) -> f64 {
        domain
            .iter()
            .map(|&c| {
                let count = self.outside.get(&c).copied().unwrap_or(0)
                    + local.get(&c).copied().unwrap_or(0);
                (count as f64 - self.expected(c) / scale).powi(2)
            })
            .sum::<f64>()
            / domain.len().max(1) as f64
    }
}

struct Proposal {
    sequence: Vec<u8>,
    graph_k: usize,
    support: u32,
    minimum: u32,
    sources: BTreeSet<(usize, u32, String)>,
}

pub(super) struct Inferred {
    pub choices: Vec<Choice>,
    pub mean: f64,
    pub report: Value,
}

pub(super) struct Choice {
    sequence: Vec<u8>,
    paired_index: usize,
    coverage_rank: usize,
}

pub(super) struct WindowChoices {
    pub start: usize,
    pub end: usize,
    pub locus: usize,
    pub choices: Vec<Choice>,
    pub mean: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn infer(
    args: &MappingArgs,
    reads: &Reads,
    backbone: &[u8],
    start: usize,
    end: usize,
    locus: usize,
    main: &Graph,
    smaller: &mut HashMap<usize, Graph>,
    library: Option<&Library>,
    cross_k_graphs: &mut Option<crate::denovo::multik::Graphs>,
) -> Result<Inferred> {
    let Some(depth) = Depth::new(main, backbone, start, end) else {
        return Ok(Inferred {
            choices: vec![],
            mean: 0.0,
            report: json!({"status":"insufficient_unique_flank_depth","final_verified":false}),
        });
    };
    let left = &backbone[start..start + main.k];
    let right = &backbone[end - main.k..end];
    let mut searches = Vec::new();
    let mut pool = BTreeMap::<Vec<u8>, Proposal>::new();
    // Weak edges can propose hypotheses, never modify the strict graph. Try
    // ordinary support first; record singleton fallback only if disconnected.
    for minimum in [3u32, 1] {
        if minimum == 1 && !pool.is_empty() {
            break;
        }
        for k in [main.k, 21, 15]
            .into_iter()
            .filter(|&k| k <= main.k)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .rev()
        {
            let graph = if k == main.k {
                main
            } else {
                if let std::collections::hash_map::Entry::Vacant(entry) = smaller.entry(k) {
                    entry.insert(Graph::build(
                        reads,
                        k,
                        args.min_baseq,
                        args.max_records,
                        args.local_max_kmers,
                    )?);
                }
                &smaller[&k]
            };
            let Some(guidance) = Depth::new(graph, backbone, start, end) else {
                continue;
            };
            let (paths, mut report) = graph.coverage_paths(
                left,
                right,
                args.local_max_span,
                args.local_max_paths,
                minimum,
                &guidance,
            );
            report["coverage_model"] = guidance.report;
            searches.push(report);
            for (sequence, support) in paths {
                let p = pool.entry(sequence.clone()).or_insert(Proposal {
                    sequence,
                    graph_k: k,
                    support,
                    minimum,
                    sources: BTreeSet::new(),
                });
                p.sources
                    .insert((k, minimum, "legacy_coverage_walk".to_owned()));
            }
        }
    }
    if args.multi_k_layout {
        for minimum in [3, 1] {
            if minimum == 1
                && pool
                    .values()
                    .any(|p| p.sources.iter().any(|(_, m, _)| *m > 1))
            {
                break;
            }
            if cross_k_graphs.is_none() {
                *cross_k_graphs = Some(crate::denovo::multik::build(
                    reads,
                    args.min_baseq,
                    args.max_records,
                )?);
            }
            let graphs = cross_k_graphs.as_ref().unwrap();
            let flank_len = graphs
                .iter()
                .rev()
                .find(|(_, c)| c.values().any(|&n| n >= minimum))
                .map(|(&k, _)| k - 1)
                .unwrap_or(main.k)
                .max(main.k);
            let wide_start = start.saturating_sub(flank_len.saturating_sub(main.k));
            let wide_end = (end + flank_len.saturating_sub(main.k)).min(backbone.len());
            let left = &backbone[wide_start..start + main.k];
            let right = &backbone[end - main.k..wide_end];
            let mut outside = occurrences(backbone, true, main.k);
            for (c, n) in occurrences(&backbone[wide_start..wide_end], false, main.k) {
                *outside.get_mut(&c).expect("window in backbone") -= n;
            }
            let guide = crate::denovo::multik::Guide {
                k: main.k,
                observed: &main.r1_occurrences,
                outside,
                mean: depth.mean,
            };
            for &k in graphs.keys().rev() {
                let graph_file =
                    format!("local_assembly/locus_{locus}.multik_k{k}_min{minimum}.gfa");
                let path = args.out.join(&graph_file);
                let (found, mut report) =
                    crate::denovo::multik::search(crate::denovo::multik::Request {
                        graphs,
                        k,
                        minimum,
                        left,
                        right,
                        guard: None,
                        guide: &guide,
                        max_span: args.local_max_span + (start - wide_start) + (wide_end - end),
                        max_paths: args.local_max_paths,
                        graph_file: Some(&path),
                    })?;
                report["graph_file"] = if path.exists() {
                    json!(graph_file)
                } else {
                    Value::Null
                };
                report["flank_padding_left"] = json!(start - wide_start);
                report["flank_padding_right"] = json!(wide_end - end);
                searches.push(report);
                for p in found {
                    // Padding belongs to unchanged sequence outside the patch.
                    // A path that overlaps through that padding cannot replace
                    // this interval while retaining both original flank words.
                    let padding = (start - wide_start) + (wide_end - end);
                    if p.sequence.len() < padding + main.k {
                        let report = searches.last_mut().unwrap();
                        let n = report["overlapping_padding_paths"].as_u64().unwrap_or(0);
                        report["overlapping_padding_paths"] = json!(n + 1);
                        continue;
                    }
                    let sequence = p.sequence
                        [start - wide_start..p.sequence.len() - (wide_end - end)]
                        .to_vec();
                    if !sequence.starts_with(&backbone[start..start + main.k])
                        || !sequence.ends_with(&backbone[end - main.k..end])
                    {
                        continue;
                    }
                    let entry = pool.entry(sequence.clone()).or_insert(Proposal {
                        sequence,
                        graph_k: k,
                        support: p.support,
                        minimum,
                        sources: BTreeSet::new(),
                    });
                    for schedule in p.schedules {
                        entry
                            .sources
                            .insert((k, minimum, format!("multi_k_{schedule}")));
                    }
                }
            }
        }
    }
    if pool.is_empty() {
        return Ok(Inferred {
            choices: vec![],
            mean: depth.mean,
            report: json!({"status":"no_connected_layout","searches":searches,"coverage_model":depth.report,"final_verified":false}),
        });
    }
    let pool: Vec<_> = pool.into_values().collect();
    let counts: Vec<_> = pool
        .iter()
        .map(|p| occurrences(&p.sequence, false, main.k))
        .collect();
    let baseline_counts = occurrences(&backbone[start..end], false, main.k);
    let mut domain: BTreeSet<_> = baseline_counts.keys().copied().collect();
    for counts in &counts {
        domain.extend(counts.keys().copied());
    }
    domain.retain(|c| {
        counts
            .iter()
            .any(|n| n.get(c).copied().unwrap_or(0) != baseline_counts.get(c).copied().unwrap_or(0))
    });
    let mut order: Vec<_> = (0..pool.len()).collect();
    order.sort_by(|&a, &b| {
        depth
            .fit(&counts[a], &domain, 1.0)
            .total_cmp(&depth.fit(&counts[b], &domain, 1.0))
            .then(pool[a].sequence.cmp(&pool[b].sequence))
    });
    let sensitivity: Vec<_> = [0.8, 1.0, 1.2].into_iter().map(|scale| {
        let &i = order.iter().min_by(|&&a, &&b| depth.fit(&counts[a], &domain, scale).total_cmp(&depth.fit(&counts[b], &domain, scale)).then(pool[a].sequence.cmp(&pool[b].sequence))).unwrap();
        json!({"baseline_multiplier":scale,"coverage_rank":order.iter().position(|&p|p==i).unwrap()+1,"length":pool[i].sequence.len()})
    }).collect();
    let records: Vec<_> = order.iter().enumerate().map(|(rank, &i)| {
        (format!("locus_{locus}_coverage_{} graph_k={} minimum_pair_support={} final_verified=false",rank+1,pool[i].graph_k,pool[i].minimum),pool[i].sequence.clone())
    }).collect();
    let candidate_file = format!("local_assembly/locus_{locus}.coverage.fasta");
    write_fasta_records(&args.out.join(&candidate_file), &records)?;
    let mut alternatives = vec![Alternative {
        sequence: backbone[start..end].to_vec(),
        support: super::path_support(main, &backbone[start..end]),
    }];
    let shortlist = crate::denovo::shortlist::select(
        &pool
            .iter()
            .map(|p| crate::denovo::shortlist::Candidate {
                sequence: &p.sequence,
                k: p.graph_k,
                circular: false,
            })
            .collect::<Vec<_>>(),
        &order,
        MAX_SCORED_PATHS,
    );
    let mut compared = Vec::new();
    for &i in &shortlist.indices {
        let index = if pool[i].sequence == alternatives[0].sequence {
            0
        } else {
            alternatives.push(Alternative {
                sequence: pool[i].sequence.clone(),
                support: pool[i].support,
            });
            alternatives.len() - 1
        };
        compared.push((index, i));
    }
    tracing::info!(
        locus,
        candidates = pool.len(),
        scored = alternatives.len(),
        "compare coverage-shortlisted mapping paths"
    );
    let evaluation = scoring::evaluate(
        backbone,
        start,
        end,
        &alternatives,
        reads,
        library,
        main.max_read,
        args,
        Some(&args.out.join(format!(
            "local_assembly/locus_{locus}.coverage.fragments.jsonl.gz"
        ))),
        Some(
            &pool
                .iter()
                .map(|p| p.sequence.as_slice())
                .collect::<Vec<_>>(),
        ),
    )?;
    let measured = evaluation.report["common_fragments"].as_u64().unwrap_or(0)
        > evaluation.report["unscored_fragments"]
            .as_u64()
            .unwrap_or(0);
    let candidate_score = |index: usize| {
        evaluation.report["candidates"][index]["score"]
            .as_f64()
            .unwrap_or(f64::NEG_INFINITY)
    };
    // Pair scores rank the frozen coverage shortlist. Coverage is only a
    // proposal/ranking heuristic, not another likelihood added to these reads.
    // If the baseline is not a sample-graph path, retain its competing score
    // explicitly instead of silently claiming it lost to the forced hypothesis.
    let &(selected, winner) = if measured {
        compared
            .iter()
            .max_by(|a, b| {
                candidate_score(a.0)
                    .total_cmp(&candidate_score(b.0))
                    .then_with(|| {
                        depth.fit(&counts[b.1], &domain, 1.0).total_cmp(&depth.fit(
                            &counts[a.1],
                            &domain,
                            1.0,
                        ))
                    })
                    .then(b.1.cmp(&a.1))
            })
            .unwrap()
    } else {
        &compared[0]
    };
    let unknown = evaluation.report["unscored_fragments"]
        .as_u64()
        .unwrap_or(0) as usize;
    let omitted_bound = crate::denovo::evidence::unknown_contrast_bound(unknown);
    let ranking: Vec<_> = order.iter().enumerate().map(|(rank,&i)| json!({
        "coverage_rank":rank+1,"length":pool[i].sequence.len(),"graph_k":pool[i].graph_k,
        "minimum_pair_support":pool[i].minimum,"weakest_edge_pair_count":pool[i].support,
        "source_graphs":pool[i].sources.iter().map(|(k,m,s)|json!({"k":k,"minimum_pair_support":m,"method":s})).collect::<Vec<_>>(),
        "copy_rmse":depth.fit(&counts[i],&domain,1.0).sqrt(),
        "unobserved_r1_kmers":counts[i].keys().filter(|c|!main.r1_occurrences.contains_key(c)).count(),
        "paired_candidate_index":compared.iter().find(|(_,p)|*p==i).map(|(index,_)|index),
    })).collect();
    let report = json!({"status":"coverage_inferred_layout","final_verified":false,"search_complete":false,
        "candidate_count":pool.len(),"scored_sample_paths":compared.len(),"shortlist":shortlist.report,"candidate_file":candidate_file,
        "selected_candidate":selected,"selected_coverage_rank":order.iter().position(|&p|p==winner).unwrap()+1,
        "selected_length":pool[winner].sequence.len(),"graph_k":pool[winner].graph_k,"minimum_pair_support":pool[winner].minimum,
        "selection":if measured{"paired_score_within_structural_shortlist"}else{"coverage_only_no_scored_fragments"},
        "baseline_preferred_by_pairs":measured&&candidate_score(0)>candidate_score(selected)+5.0+omitted_bound,
        "coverage_model":depth.report,"common_scoring_domain_kmers":domain.len(),"baseline_sensitivity":sensitivity,
        "ranked_candidates":ranking,"searches":searches,"paired_comparison":evaluation.report,
        "multi_k_layout":{"enabled":args.multi_k_layout,"requested_k":crate::denovo::multik::KS,"additional_state_budget_per_k":if args.multi_k_layout {50_000}else{0},"legacy_candidates_retained":true,"common_fragment_comparison_after_union":true},
        "limitations":["bounded graph and structural shortlist are not exhaustive","pair and coverage scores reuse the same reads and are not independent evidence","copy depth assumes comparable sampling and no unmodelled motif sources","local layouts do not phase distant windows or verify circularity"]});
    let best_graph_score = compared
        .iter()
        .map(|&(index, _)| candidate_score(index))
        .fold(f64::NEG_INFINITY, f64::max);
    let best_pair_score = best_graph_score.max(candidate_score(0));
    // Omitted fragments widen the pairwise contrast by a known finite bound.
    // They cannot erase arbitrarily strong evidence against a path.
    let retention_margin = 5.0 + omitted_bound;
    let mut choices: Vec<_> = compared
        .iter()
        .filter(|&&(index, _)| {
            !measured || candidate_score(index) >= best_pair_score - retention_margin
        })
        .map(|&(paired_index, i)| Choice {
            sequence: pool[i].sequence.clone(),
            paired_index,
            coverage_rank: order.iter().position(|&p| p == i).unwrap() + 1,
        })
        .collect();
    // The pileup baseline remains a real competitor even when coverage ranks
    // it below the sample-path shortlist or graph search never rediscovers it.
    // Otherwise a forced contraction could be chosen despite opposing reads.
    if !compared.iter().any(|&(index, _)| index == 0)
        && (!measured || candidate_score(0) >= best_pair_score - retention_margin)
    {
        choices.push(Choice {
            sequence: backbone[start..end].to_vec(),
            paired_index: 0,
            coverage_rank: 0,
        });
    }
    let mut report = report;
    report["choice_retention_log_margin"] = json!(retention_margin);
    report["baseline_pair_margin_to_best_graph"] = json!(candidate_score(0) - best_graph_score);
    Ok(Inferred {
        choices,
        mean: depth.mean,
        report,
    })
}

fn replace_counts(counts: &mut Counts, before: &[u8], after: &[u8], k: usize) {
    for (c, n) in occurrences(before, false, k) {
        *counts.get_mut(&c).unwrap() -= n;
    }
    for (c, n) in occurrences(after, false, k) {
        *counts.entry(c).or_default() += n;
    }
}

fn patched_sequence(backbone: &[u8], patches: &[Patch]) -> Vec<u8> {
    let mut sequence = Vec::new();
    let mut cursor = 0;
    for patch in patches {
        sequence.extend_from_slice(&backbone[cursor..patch.start]);
        sequence.extend_from_slice(&patch.sequence);
        cursor = patch.end;
    }
    sequence.extend_from_slice(&backbone[cursor..]);
    sequence
}

/// Allocate occurrence depth once across all inferred windows. Independently
/// fitting each locus would spend the same repeat-copy evidence multiple times.
/// Pair scores only screen choices *within* each window; they are not summed
/// across windows, whose fragment cohorts may overlap and do not phase loci.
pub(super) fn joint_select(
    args: &MappingArgs,
    backbone: &[u8],
    graph: &Graph,
    strict: &[Patch],
    options: &[WindowChoices],
    reports: &mut [Value],
) -> Result<(Vec<Patch>, Value)> {
    if options.is_empty() {
        return Ok((strict.to_vec(), json!({"status":"no_choices"})));
    }
    let mut fixed = occurrences(backbone, true, graph.k);
    for p in strict {
        replace_counts(&mut fixed, &backbone[p.start..p.end], &p.sequence, graph.k);
    }
    let mut varying = BTreeMap::<kmer::Code, (f64, usize)>::new();
    let mut differences = Vec::new();
    for window in options {
        let before = occurrences(&backbone[window.start..window.end], false, graph.k);
        let mut local_codes: BTreeSet<_> = before.keys().copied().collect();
        let after: Vec<_> = window
            .choices
            .iter()
            .map(|c| occurrences(&c.sequence, false, graph.k))
            .collect();
        for counts in &after {
            local_codes.extend(counts.keys().copied());
        }
        local_codes.retain(|c| {
            after
                .iter()
                .any(|n| n.get(c).copied().unwrap_or(0) != before.get(c).copied().unwrap_or(0))
        });
        for &c in &local_codes {
            let entry = varying.entry(c).or_default();
            entry.0 += window.mean;
            entry.1 += 1;
        }
        differences.push(
            after
                .iter()
                .map(|a| {
                    local_codes
                        .iter()
                        .map(|&c| {
                            (
                                c,
                                a.get(&c).copied().unwrap_or(0) as i32
                                    - before.get(&c).copied().unwrap_or(0) as i32,
                            )
                        })
                        .collect::<BTreeMap<_, _>>()
                })
                .collect::<Vec<_>>(),
        );
    }
    let codes: Vec<_> = varying.keys().copied().collect();
    let expected: Vec<_> = codes
        .iter()
        .map(|c| {
            graph.r1_occurrences.get(c).copied().unwrap_or(0) as f64
                / (varying[c].0 / varying[c].1 as f64)
        })
        .collect();
    let base: Vec<_> = codes
        .iter()
        .map(|c| fixed.get(c).copied().unwrap_or(0) as i32)
        .collect();
    let delta: Vec<Vec<Vec<i32>>> = differences
        .iter()
        .map(|window| {
            window
                .iter()
                .map(|choice| {
                    codes
                        .iter()
                        .map(|c| choice.get(c).copied().unwrap_or(0))
                        .collect()
                })
                .collect()
        })
        .collect();
    let mut remaining = vec![vec![(0i32, 0i32); codes.len()]; options.len() + 1];
    for i in (0..options.len()).rev() {
        for j in 0..codes.len() {
            remaining[i][j] = (
                remaining[i + 1][j].0 + delta[i].iter().map(|d| d[j]).min().unwrap(),
                remaining[i + 1][j].1 + delta[i].iter().map(|d| d[j]).max().unwrap(),
            );
        }
    }
    let lower_bound = |counts: &[i32], remaining: &[(i32, i32)]| {
        counts
            .iter()
            .zip(remaining)
            .zip(&expected)
            .map(|((&n, &(low, high)), &observed)| {
                let closest = observed.clamp((n + low) as f64, (n + high) as f64);
                (closest - observed).powi(2)
            })
            .sum::<f64>()
            / codes.len().max(1) as f64
    };
    // An optimistic count interval for unassigned windows keeps partial choices
    // from being penalized for copies that a later window may still explain.
    let mut beam = vec![(lower_bound(&base, &remaining[0]), base, Vec::<usize>::new())];
    let mut beam_truncated = false;
    for i in 0..options.len() {
        let mut expanded = Vec::new();
        for (_, counts, selected) in beam {
            for (j, change) in delta[i].iter().enumerate() {
                let counts: Vec<_> = counts.iter().zip(change).map(|(a, b)| a + b).collect();
                let mut selected = selected.clone();
                selected.push(j);
                expanded.push((lower_bound(&counts, &remaining[i + 1]), counts, selected));
            }
        }
        expanded.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.2.cmp(&b.2)));
        if expanded.len() > 512 {
            beam_truncated = true;
            expanded.truncate(512);
        }
        beam = expanded;
    }
    let best = &beam[0];
    let make_patches = |selected: &[usize]| {
        let mut patches = strict.to_vec();
        for (window, &j) in options.iter().zip(selected) {
            patches.push(Patch {
                start: window.start,
                end: window.end,
                locus: window.locus,
                sequence: window.choices[j].sequence.clone(),
            });
        }
        patches.sort_by_key(|p| p.start);
        patches
    };
    let configurations: Vec<_> = beam.iter().take(32).enumerate().map(|(rank,(fit,_,selected))| {
        json!({"rank":rank+1,"copy_rmse":fit.sqrt(),"choices":options.iter().zip(selected).map(|(w,&j)|json!({"locus":w.locus,"coverage_rank":w.choices[j].coverage_rank,"paired_candidate_index":w.choices[j].paired_index})).collect::<Vec<_>>()})
    }).collect();
    let records: Vec<_> = beam.iter().take(32).enumerate().map(|(rank,(_,_,selected))|(
        format!("mapping_layout_{} final_verified=false circularity_verified=false evidence=coverage_inferred",rank+1),patched_sequence(backbone,&make_patches(selected)))).collect();
    write_fasta_records(&args.out.join("coverage_layout_candidates.fasta"), &records)?;
    for (window, &j) in options.iter().zip(&best.2) {
        let choice = &window.choices[j];
        let report = &mut reports[window.locus - 1]["coverage_layout"];
        report["selected_candidate"] = json!(choice.paired_index);
        report["selected_coverage_rank"] = if choice.coverage_rank > 0 {
            json!(choice.coverage_rank)
        } else {
            Value::Null
        };
        report["selected_length"] = json!(choice.sequence.len());
        if choice.coverage_rank > 0 {
            let ranked = report["ranked_candidates"][choice.coverage_rank - 1].clone();
            report["graph_k"] = ranked["graph_k"].clone();
            report["minimum_pair_support"] = ranked["minimum_pair_support"].clone();
        } else {
            report["graph_k"] = Value::Null;
            report["minimum_pair_support"] = Value::Null;
        }
        report["selected_source"] = json!(if choice.coverage_rank > 0 {
            "sample_graph"
        } else {
            "pileup_baseline"
        });
        report["status"] = json!(if choice.sequence == backbone[window.start..window.end] {
            "baseline_retained"
        } else {
            "coverage_inferred_layout"
        });
        report["selection"] = json!("joint_copy_depth_among_retained_local_hypotheses");
        report["pair_compatibility_verified"] = json!(false);
        report["joint_eligible_choices"] = json!(window.choices.len());
        report["local_pair_score_tolerance"] = json!(5.0);
        let paired = &report["paired_comparison"];
        report["baseline_preferred_by_pairs"] = json!(
            paired["common_fragments"].as_u64().unwrap_or(0)
                > paired["unscored_fragments"].as_u64().unwrap_or(0)
                && paired["candidates"][0]["score"].as_f64().unwrap_or(0.0)
                    > paired["candidates"][choice.paired_index]["score"]
                        .as_f64()
                        .unwrap_or(0.0)
                        + 5.0
        );
    }
    let report = json!({"status":"joint_coverage_hypothesis","depth_k":graph.k,"common_domain_kmers":codes.len(),
        "copy_rmse":best.0.sqrt(),"depth_baseline":"arithmetic mean local baseline across windows varying each kmer",
        "retained_configurations":beam.len(),"exported_configurations":configurations.len(),"max_configurations":512,"beam_truncated":beam_truncated,
        "equally_scored_retained_configurations":beam.iter().filter(|b|(b.0-best.0).abs()<1e-9).count(),
        "configurations":configurations,"final_verified":false,"cross_window_phase_verified":false,
        "pair_scores_summed_across_windows":false});
    Ok((make_patches(&best.2), report))
}

/// Mapping hypotheses retain the local-search rotation. Every reported interval
/// is projected into this FASTA, separately from strict candidate coordinates.
pub(super) fn write_layout(
    args: &MappingArgs,
    backbone: &[u8],
    cut: usize,
    mut patches: Vec<Patch>,
    intervals: &[(usize, usize)],
    windows: &[Value],
    joint: Value,
) -> Result<Value> {
    let inferred = windows
        .iter()
        .filter(|w| w["coverage_layout"]["status"] == "coverage_inferred_layout")
        .count();
    if inferred == 0 {
        let report = json!({"enabled":true,"status":if windows.iter().any(|w|w["coverage_layout"]["status"]=="baseline_retained"){"baseline_retained"}else{"no_inferred_windows"},"final_verified":false,"joint_copy_model":joint});
        crate::write_json(&args.out.join("coverage_layout.json"), &report)?;
        return Ok(report);
    }
    patches.sort_by_key(|p| p.start);
    let mut sequence = Vec::new();
    let mut cursor = 0;
    for patch in &patches {
        sequence.extend_from_slice(&backbone[cursor..patch.start]);
        sequence.extend_from_slice(&patch.sequence);
        cursor = patch.end;
    }
    sequence.extend_from_slice(&backbone[cursor..]);
    let project = |position: usize, right: bool| {
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
    let mut ambiguity = BufWriter::new(File::create(
        args.out.join("coverage_layout_ambiguity.tsv"),
    )?);
    writeln!(
        ambiguity,
        "sequence\tstart\tend\tevidence\treason\talternatives"
    )?;
    let mut projected = Vec::new();
    for (i, &(start, end)) in intervals.iter().enumerate() {
        let (a, b) = (project(start, false), project(end, true));
        if windows[i]["status"] == "unresolved" {
            let count = windows[i]["coverage_layout"]["candidate_count"]
                .as_u64()
                .unwrap_or(0);
            writeln!(
                ambiguity,
                "coverage_layout\t{}\t{}\tmapping_local_{}\tbounded_coverage_layout_unverified\t{}",
                a + 1,
                b,
                i + 1,
                count
            )?;
        }
        projected.push(json!({"locus":i+1,"start":a+1,"end":b,"coverage_layout":windows[i].get("coverage_layout")}));
    }
    ambiguity.flush()?;
    write_fasta_records(&args.out.join("coverage_layout.fasta"),&[(
        "coverage_layout final_verified=false circularity_verified=false evidence=coverage_inferred topology=reference_assumed".to_owned(),sequence.clone())])?;
    let report = json!({"schema":"mitovar-rust.mapping-coverage-layout.v1","enabled":true,"status":"coverage_inferred_layout",
        "final_verified":false,"circularity_verified":false,"length":sequence.len(),"ambiguous_bases":sequence.iter().filter(|&&b|b==b'N').count(),
        "inferred_windows":inferred,"baseline_rotation_offset":cut,"windows":projected,"joint_copy_model":joint,
        "coordinate_system":"coverage_layout.fasta retains the internal local-assembly rotation; no reference coordinates",
        "strict_candidate_and_final_independent":true});
    crate::write_json(&args.out.join("coverage_layout.json"), &report)?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn repeated_r1_windows_recover_six_copies_without_borrowing_flank_bases() {
        use std::fs;
        let root =
            std::env::temp_dir().join(format!("mapping-coverage-copies-{}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        let mut state = 7142823u64;
        let backbone: Vec<_> = (0..1800)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 32) as usize & 3]
            })
            .collect();
        let repeat = backbone[150..230].repeat(6);
        let sample = [&backbone[..1000], repeat.as_slice(), &backbone[1000..]].concat();
        let mut baseline = sample.clone();
        baseline[1000..1480].fill(b'N');
        let reads = Reads {
            r1: root.join("r1.fq"),
            r2: root.join("r2.fq"),
        };
        let mut one = File::create(&reads.r1).unwrap();
        let mut two = File::create(&reads.r2).unwrap();
        for start in 0..sample.len() {
            let r: Vec<_> = (0..100)
                .map(|i| sample[(start + i) % sample.len()])
                .collect();
            for (mate, sequence, out) in [
                (1, r.clone(), &mut one),
                (2, kmer::reverse_complement(&r), &mut two),
            ] {
                writeln!(
                    out,
                    "@p{start}/{mate}\n{}\n+\n{}",
                    String::from_utf8_lossy(&sequence),
                    "I".repeat(100)
                )
                .unwrap();
            }
        }
        let (start, end) = (969, 1511);
        let (left, right) = (&baseline[start..start + 31], &baseline[end - 31..end]);
        for k in [31, 15] {
            let graph = Graph::build(&reads, k, 0, 10000, 2_000_000).unwrap();
            let depth = Depth::new(&graph, &baseline, start, end).unwrap();
            let (paths, report) = graph.coverage_paths(left, right, 1000, 32, 3, &depth);
            assert!(paths.iter().any(|(p, _)| p == &sample[start..end]));
            assert!(paths
                .iter()
                .all(|(p, n)| p.starts_with(left) && p.ends_with(right) && *n >= 3));
            assert_eq!(report["search_complete"], false);
            let counts: Vec<_> = paths
                .iter()
                .map(|(p, _)| occurrences(p, false, k))
                .collect();
            let domain = counts.iter().flat_map(|c| c.keys().copied()).collect();
            let best = (0..paths.len())
                .min_by(|&a, &b| {
                    depth
                        .fit(&counts[a], &domain, 1.0)
                        .total_cmp(&depth.fit(&counts[b], &domain, 1.0))
                })
                .unwrap();
            assert_eq!(paths[best].0, sample[start..end]);
            let mut wrong = left.to_vec();
            wrong[10] = if wrong[10] == b'A' { b'C' } else { b'A' };
            assert!(graph
                .coverage_paths(&wrong, right, 1000, 32, 3, &depth)
                .0
                .is_empty());
            if k == 15 {
                let c = kmer::canonical(super::super::graph::code(&backbone[150..165]).unwrap(), k);
                assert!(graph.r1_occurrences[&c] > u64::from(graph.counts[&c]));
            }
        }
        fs::remove_dir_all(root).unwrap();
    }
}
