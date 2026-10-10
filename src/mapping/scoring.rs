use super::Alternative;
use crate::{
    cli::{MappingArgs, Reads},
    denovo::{
        align::{self, Index, Prepared},
        anchor,
        evidence::{joint, observations, rescue_pool, unknown_reason, JointEvidence, RescueStats},
        pairs::Library,
    },
    io::fastq,
    kmer,
};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

pub(super) struct Evaluation {
    pub selected: Option<usize>,
    pub report: Value,
}

struct Fragment {
    reads: [fastq::Record; 2],
    pivot: anchor::Pivot,
    observed_anchors: HashSet<kmer::Code>,
}

struct FragmentEvidence {
    rescued: RescueStats,
    unknown: Option<&'static str>,
    results: Vec<JointEvidence>,
    contractions: Vec<bool>,
    diagnostic: Option<Value>,
}

/// Keep a bounded set of jobs in flight on persistent workers. Only the
/// calling thread reduces scores or writes diagnostics, in input order.
fn ordered_stream<T: Send, R: Send>(
    threads: usize,
    calculate: &(impl Fn(&T) -> Result<R> + Sync),
    mut consume: impl FnMut(R) -> Result<()>,
    produce: impl FnOnce(&mut dyn FnMut(T) -> Result<()>) -> Result<()>,
) -> Result<()> {
    if threads == 1 {
        return produce(&mut |item| consume(calculate(&item)?));
    }
    std::thread::scope(|scope| {
        use std::sync::{mpsc, Arc, Mutex};
        let capacity = threads * 4;
        let (input_tx, input_rx) = mpsc::sync_channel::<(usize, T)>(capacity);
        let input_rx = Arc::new(Mutex::new(input_rx));
        let (output_tx, output_rx) = mpsc::channel::<(usize, Result<R>)>();
        let mut jobs = Vec::new();
        for _ in 0..threads {
            let input = Arc::clone(&input_rx);
            let output = output_tx.clone();
            jobs.push(
                std::thread::Builder::new().spawn_scoped(scope, move || loop {
                    // Release the queue lock before doing any fragment work.
                    let task = { input.lock().unwrap().recv() };
                    let Ok((ordinal, item)) = task else { break };
                    let value =
                        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| calculate(&item)))
                            .unwrap_or_else(|_| {
                                Err(anyhow::anyhow!("fragment scoring worker panicked"))
                            });
                    if output.send((ordinal, value)).is_err() {
                        break;
                    }
                })?,
            );
        }
        drop(output_tx);
        let mut ready = BTreeMap::new();
        let mut next = 0;
        let mut consume_next = || -> Result<()> {
            while !ready.contains_key(&next) {
                let (ordinal, value) = output_rx.recv().map_err(|_| {
                    anyhow::anyhow!("fragment scoring workers ended before their results")
                })?;
                ready.insert(ordinal, value);
            }
            consume(ready.remove(&next).unwrap()?)?;
            next += 1;
            Ok(())
        };
        let (mut submitted, mut outstanding) = (0, 0);
        let produced = produce(&mut |item| {
            input_tx.send((submitted, item)).map_err(|_| {
                anyhow::anyhow!("fragment scoring workers ended before their inputs")
            })?;
            submitted += 1;
            outstanding += 1;
            if outstanding == capacity {
                consume_next()?;
                outstanding -= 1;
            }
            Ok(())
        });
        // Also close the input on producer or consumer error, before joining.
        drop(input_tx);
        produced?;
        for _ in 0..outstanding {
            consume_next()?;
        }
        for job in jobs {
            job.join()
                .map_err(|_| anyhow::anyhow!("fragment scoring worker panicked"))?;
        }
        Ok(())
    })
}

fn changed_interval(baseline: &[u8], alternative: &[u8], start: usize, k: usize) -> [usize; 2] {
    if baseline == alternative {
        return [start + k, start + alternative.len() - k];
    }
    let prefix = baseline
        .iter()
        .zip(alternative)
        .take_while(|(a, b)| a == b)
        .count();
    let suffix = baseline[prefix..]
        .iter()
        .rev()
        .zip(alternative[prefix..].iter().rev())
        .take_while(|(a, b)| a == b)
        .count();
    [start + prefix, start + alternative.len() - suffix]
}

fn contraction_flanks(
    sequence: &[u8],
    window: std::ops::Range<usize>,
    junctions: [usize; 2],
    unique: &HashSet<kmer::Code>,
) -> Option<[usize; 2]> {
    let is_unique = |p: usize| {
        super::graph::code(&sequence[p..p + 21])
            .is_some_and(|c| unique.contains(&kmer::canonical(c, 21)))
    };
    let left_end = junctions[0].checked_sub(21)?;
    let right_end = window.end.checked_sub(21)?;
    let left = (window.start..=left_end).rev().find(|&p| is_unique(p))?;
    let right = (junctions[1]..=right_end).find(|&p| is_unique(p))?;
    Some([left, right])
}

#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate(
    backbone: &[u8],
    start: usize,
    end: usize,
    alternatives: &[Alternative],
    reads: &Reads,
    model: Option<&Library>,
    max_read: usize,
    args: &MappingArgs,
    diagnostic_path: Option<&Path>,
    cohort_paths: Option<&[&[u8]]>,
) -> Result<Evaluation> {
    let mut pool_seeds = HashSet::new();
    let local = &backbone[start..end];
    let cohort: std::collections::BTreeSet<&[u8]> = alternatives
        .iter()
        .map(|p| p.sequence.as_slice())
        .chain(cohort_paths.into_iter().flatten().copied())
        .chain(std::iter::once(local))
        .collect();
    let prefix = cohort
        .iter()
        .map(|p| {
            local
                .iter()
                .zip(p.iter())
                .take_while(|(a, b)| a == b)
                .count()
        })
        .min()
        .unwrap_or(0);
    let suffix = cohort
        .iter()
        .map(|p| {
            local
                .iter()
                .rev()
                .zip(p.iter().rev())
                .take_while(|(a, b)| a == b)
                .count()
        })
        .min()
        .unwrap_or(0)
        .min(local.len() - prefix);
    for flank in [
        &backbone[start.saturating_sub(500)..start + prefix],
        &backbone[end - suffix..(end + 500).min(backbone.len())],
    ] {
        kmer::visit(flank, None, 21, 0, |_, c| {
            pool_seeds.insert(kmer::canonical(c, 21));
        });
    }
    let sequences: Vec<Vec<u8>> = alternatives
        .iter()
        .map(|p| [&backbone[..start], p.sequence.as_slice(), &backbone[end..]].concat())
        .collect();
    // Freeze one flank-anchored cohort before scoring. Internal repeat-only
    // fragments cannot locate a new copy and must not become copy-count votes.
    // An anchor must be single-copy in *every* whole-contig hypothesis, which
    // includes the original source locus and the potential inserted copy.
    for path in &cohort {
        let sequence = [&backbone[..start], *path, &backbone[end..]].concat();
        let mut occurrences = HashMap::<kmer::Code, usize>::new();
        let extended = [sequence.as_slice(), &sequence[..20.min(sequence.len())]].concat();
        kmer::visit(&extended, None, 21, 0, |_, c| {
            *occurrences.entry(kmer::canonical(c, 21)).or_default() += 1
        });
        pool_seeds.retain(|c| occurrences.get(c) == Some(&1));
    }
    let indexes: Vec<_> = sequences
        .iter()
        .map(|s| Index::new(s, true, max_read))
        .collect();
    let locations: Vec<_> = indexes
        .iter()
        .map(|i| anchor::locations(i, &pool_seeds))
        .collect();
    let mut log = diagnostic_path
        .map(|p| {
            File::create(p).map(|f| {
                flate2::write::GzEncoder::new(BufWriter::new(f), flate2::Compression::default())
            })
        })
        .transpose()?;
    let mut unknown_reasons = BTreeMap::<&str, usize>::new();
    let junctions: Vec<_> = alternatives
        .iter()
        .map(|p| changed_interval(&backbone[start..end], &p.sequence, start, args.local_k))
        .collect();
    // A contraction can turn two repeated junctions into the same local match.
    // Separate junction counts then establish neither their linkage nor copy
    // count, especially against an N-masked baseline. Require the same fragment
    // to observe unique unchanged sequence on both sides of a >=k contraction.
    let requires_contraction: Vec<_> = alternatives
        .iter()
        .map(|p| (end - start).saturating_sub(p.sequence.len()) >= args.local_k)
        .collect();
    let flanks: Vec<_> = sequences
        .iter()
        .enumerate()
        .map(|(i, s)| {
            requires_contraction[i]
                .then(|| {
                    contraction_flanks(
                        s,
                        start..start + alternatives[i].sequence.len(),
                        junctions[i],
                        &pool_seeds,
                    )
                })
                .flatten()
        })
        .collect();
    let n = alternatives.len();
    let mut scores = vec![0.0f64; n];
    let mut observed_bridges = vec![[0usize; 2]; n];
    let mut mate_brackets = vec![[0usize; 2]; n];
    let mut votes = vec![vec![0usize; n]; n];
    let mut bridges = vec![[0usize; 2]; n];
    let mut contraction_bridges = vec![0usize; n];
    let mut names = HashSet::new();
    let (mut fragments, mut unknown) = (0usize, 0usize);
    let mut rescued = RescueStats::default();
    let threads = usize::from(args.scoring_threads)
        .min(std::thread::available_parallelism().map_or(1, usize::from))
        .max(1);
    let diagnostics = log.is_some();
    let calculate = |fragment: &Fragment| -> Result<FragmentEvidence> {
        let [a, b] = [&fragment.reads[0], &fragment.reads[1]];
        let pivot = &fragment.pivot;
        let anchored = [a, b][pivot.mate];
        if indexes
            .iter()
            .zip(&locations)
            .any(|(i, map)| !pivot.consistent(map, anchored.sequence.len(), i.circle_length()))
        {
            return Ok(FragmentEvidence {
                rescued: RescueStats::default(),
                unknown: Some("incompatible_common_anchor_projections"),
                results: Vec::new(),
                contractions: Vec::new(),
                diagnostic: diagnostics.then(|| json!({"fragment":a.name(1).unwrap(),"pivot":pivot.report(),"unknown_reason":"incompatible_common_anchor_projections"})),
            });
        }
        let ideal = align::ideal_score(a) + align::ideal_score(b);
        let mut one = Prepared::new(a);
        let mut two = Prepared::new(b);
        let mut placements: Vec<_> = indexes
            .iter()
            .map(|_| (align::Placements::default(), align::Placements::default()))
            .collect();
        let mut removed = Vec::new();
        for (i, (x, y)) in placements.iter_mut().enumerate() {
            let (hits, prepared) = if pivot.mate == 0 {
                (x, &mut one)
            } else {
                (y, &mut two)
            };
            removed.push(pivot.constrain(
                hits,
                anchored,
                prepared,
                &indexes[i],
                &locations[i],
                args.min_baseq,
            ));
        }
        let rs = rescue_pool(
            &indexes,
            &mut placements,
            &mut one,
            &mut two,
            model,
            pivot.mate,
        );
        if let Some(reason) = unknown_reason(&placements, &rs.states, pivot.mate) {
            let diagnostic = diagnostics.then(|| json!({"fragment":a.name(1).unwrap(),"pivot":pivot.report(),"search_status":rs.states,"unknown_reason":reason}));
            return Ok(FragmentEvidence {
                rescued: rs,
                unknown: Some(reason),
                results: Vec::new(),
                contractions: Vec::new(),
                diagnostic,
            });
        }
        let mut results = Vec::with_capacity(n);
        for (i, ((x, y), index)) in placements.iter().zip(&indexes).enumerate() {
            results.push(joint(
                &observations(x, a, index, args.min_baseq),
                &observations(y, b, index, args.min_baseq),
                model,
                ideal,
                sequences[i].len(),
                junctions[i],
                flanks[i],
            ));
        }
        let contractions = results
            .iter()
            .enumerate()
            .map(|(i, result)| {
                result.contraction
                    && flanks[i].is_some_and(|anchors| {
                        anchors.iter().all(|&p| {
                            let c = super::graph::code(&sequences[i][p..p + 21]).unwrap();
                            fragment.observed_anchors.contains(&kmer::canonical(c, 21))
                        })
                    })
            })
            .collect();
        let diagnostic = diagnostics.then(|| {
            let candidates: Vec<_> = results.iter_mut().enumerate().map(|(i, result)| json!({"candidate":i,"pivot_projection":pivot.projection(&locations[i],anchored.sequence.len()),"anchor_hits_removed":removed[i],"r1_placements":placements[i].0.hits.len(),"r2_placements":placements[i].1.hits.len(),"search_status":rs.states[i],"joint":std::mem::take(&mut result.diagnostics),"observed_junctions":result.observed_junctions,"mate_brackets":result.mate_brackets})).collect();
            json!({"fragment":a.name(1).unwrap(),"pivot":pivot.report(),"unknown_reason":null,"candidates":candidates})
        });
        Ok(FragmentEvidence {
            rescued: rs,
            unknown: None,
            results,
            contractions,
            diagnostic,
        })
    };
    let consume = |measured: FragmentEvidence| -> Result<()> {
        rescued.searches += measured.rescued.searches;
        rescued.added_placements += measured.rescued.added_placements;
        rescued.seeded_decoy_searches += measured.rescued.seeded_decoy_searches;
        if let Some(reason) = measured.unknown {
            unknown += 1;
            *unknown_reasons.entry(reason).or_default() += 1;
        } else {
            let results = &measured.results;
            for i in 0..n {
                scores[i] += results[i].score;
                contraction_bridges[i] += usize::from(measured.contractions[i]);
                for j in 0..2 {
                    bridges[i][j] += usize::from(results[i].junctions[j]);
                    observed_bridges[i][j] += usize::from(results[i].observed_junctions[j]);
                    mate_brackets[i][j] += usize::from(results[i].mate_brackets[j]);
                }
                for j in 0..n {
                    votes[i][j] += usize::from(results[i].score - results[j].score >= 2.3);
                }
            }
        }
        if let (Some(log), Some(diagnostic)) = (log.as_mut(), measured.diagnostic) {
            serde_json::to_writer(&mut *log, &diagnostic)?;
            writeln!(log)?;
        }
        Ok(())
    };
    // Only eligible fragments enter the bounded queue. No reads are sampled,
    // and duplicate names and input limits are checked in their original order.
    ordered_stream(threads, &calculate, consume, |emit| {
        fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
            if !names.insert(a.name(1)?.to_owned()) {
                return Ok(());
            }
            ensure!(
                names.len() <= args.max_records,
                "local scoring exceeds --max-records named FASTQ pairs"
            );
            let Some(pivot) = anchor::Pivot::choose([a, b], &pool_seeds, args.min_baseq) else {
                return Ok(());
            };
            fragments += 1;
            let mut observed_anchors = HashSet::new();
            if requires_contraction.iter().any(|&required| required) {
                for read in [a, b] {
                    kmer::visit(
                        &read.sequence,
                        Some(&read.quality),
                        21,
                        args.min_baseq,
                        |_, c| {
                            let c = kmer::canonical(c, 21);
                            if pool_seeds.contains(&c) {
                                observed_anchors.insert(c);
                            }
                        },
                    );
                }
            }
            emit(Fragment {
                reads: [a.clone(), b.clone()],
                pivot,
                observed_anchors,
            })
        })?;
        Ok(())
    })?;
    if let Some(log) = log {
        log.finish()?.flush()?;
    }
    let best = (0..n)
        .max_by(|&a, &b| scores[a].total_cmp(&scores[b]).then(b.cmp(&a)))
        .unwrap();
    let fraction = args.min_alt_frac.unwrap_or(args.min_indel_frac);
    let discriminated = (0..n).filter(|&i| i != best).all(|i| {
        votes[best][i] >= 3
            && votes[best][i] as f64 >= fraction * (votes[best][i] + votes[i][best]) as f64
            && scores[best] - scores[i]
                > 5.0 + crate::denovo::evidence::unknown_contrast_bound(unknown)
    });
    let supported = bridges[best].iter().all(|&c| c >= 3)
        && alternatives[best].support >= args.min_depth
        && !alternatives[best].sequence.contains(&b'N');
    let contraction_supported = !requires_contraction[best] || contraction_bridges[best] >= 3;
    let selected =
        (discriminated && supported && contraction_supported && unknown == 0).then_some(best);
    Ok(Evaluation {
        selected,
        report: json!({"scoring_status":"evaluated","alignment_search_exhaustive":false,"alignment_domain":align::profile::report(),"cohort_generated_paths_including_baseline":cohort.len(),"unknown_score_contrast_bound":crate::denovo::evidence::unknown_contrast_bound(unknown),"best_candidates":(0..n).filter(|&i|(scores[i]-scores[best]).abs()<=1e-8).collect::<Vec<_>>(),"unknown_reasons":unknown_reasons,"fragment_diagnostics":diagnostic_path.and_then(|p|p.file_name()).map(|p|p.to_string_lossy()),"mate_rescue":{"searched_hypothesis_mates":rescued.searches,"seeded_decoy_searches":rescued.seeded_decoy_searches,"added_placements":rescued.added_placements,"scope":"fixed common-flank mate, query offset and projected target word; complete endpoint/edit profiles within declared empirical neighborhoods; work-limit omissions shared across candidates"},"common_fragments":fragments,"shared_unique_anchor_kmers":pool_seeds.len(),"cohort":"fixed observed 21-mer in actually retained common backbone/flanks, unique in every whole-contig candidate; anchored mate scoring constrained by query-target traceback","unscored_fragments":unknown,"best_candidate":best,
            "selected_candidate":selected,"status":if selected.is_some(){"resolved"}else{"unresolved"},
            "reason":if !discriminated{"competing_paths"}else if !supported{"insufficient_path_or_junction_support"}else if !contraction_supported{"contraction_unique_flanks_not_linked"}else if unknown>0{"alignment_search_limit"}else{"supported"},
            "candidates":alternatives.iter().enumerate().map(|(i,p)|json!({"candidate":i,"length":p.sequence.len(),"minimum_kmer_fragments":p.support,"score":scores[i],"strong_votes_against":votes[i],"junction_fragments":bridges[i],"observed_junction_fragments":observed_bridges[i],"mate_bracket_fragments":mate_brackets[i],"junction_positions_rotated_baseline_start":junctions[i],"requires_contraction_flanks":requires_contraction[i],"contraction_flank_starts_0based":flanks[i],"contraction_flank_fragments":contraction_bridges[i]})).collect::<Vec<_>>(),
            "model":"one named fragment once; full-query affine alignment; best explanation per empirical span; original copy retained in every whole-contig candidate",
            "junction_evidence":"one read continuously matches 20 bases on each side, or geometry-compatible mates bracket the junction with observed 20-base terminal anchors; mate gaps are not sequenced junctions",
            "contraction_evidence":"the same joint placement observes both unique unchanged 21-base flanks in matching traceback",
            "scores_are_not_correctness_probabilities":true}),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Commands};
    use clap::Parser;
    use std::{fs, io::Write};

    #[test]
    fn persistent_workers_bound_jobs_and_preserve_order_on_success_and_error() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let mut seen = Vec::new();
        let mut workers = HashSet::new();
        ordered_stream(
            3,
            &|&value| {
                if value % 7 == 0 {
                    std::thread::yield_now();
                }
                Ok((value, std::thread::current().id()))
            },
            |(value, thread)| {
                seen.push(value);
                workers.insert(thread);
                Ok(())
            },
            |emit| {
                for value in 0..100 {
                    emit(value)?;
                }
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(seen, (0..100).collect::<Vec<_>>());
        assert!(workers.len() <= 3);

        let calculated = AtomicUsize::new(0);
        let failed = ordered_stream(
            3,
            &|&value| {
                calculated.fetch_add(1, Ordering::Relaxed);
                Ok(value)
            },
            |value| {
                ensure!(value != 5, "consumer stopped");
                Ok(())
            },
            |emit| {
                for value in 0..1000 {
                    emit(value)?;
                }
                Ok(())
            },
        );
        assert!(failed.unwrap_err().to_string().contains("consumer stopped"));
        assert!(calculated.load(Ordering::Relaxed) <= 18);

        let failed = ordered_stream(
            3,
            &|&value| Ok(value),
            |_| Ok(()),
            |emit| {
                for value in 0..7 {
                    emit(value)?;
                }
                anyhow::bail!("producer stopped");
            },
        );
        assert!(failed.unwrap_err().to_string().contains("producer stopped"));
    }

    #[test]
    fn worker_failures_reach_the_caller_without_waiting_for_a_missing_result() {
        for panic in [false, true] {
            let mut seen = Vec::new();
            let failed = ordered_stream(
                3,
                &|&value| {
                    if value == 7 {
                        assert!(!panic, "test worker panic");
                        anyhow::bail!("test worker error");
                    }
                    Ok(value)
                },
                |value| {
                    seen.push(value);
                    Ok(())
                },
                |emit| {
                    for value in 0..40 {
                        emit(value)?;
                    }
                    Ok(())
                },
            );
            assert!(failed.is_err());
            assert_eq!(seen, (0..7).collect::<Vec<_>>());
        }
    }

    #[test]
    fn parallel_fragments_preserve_scores_unknowns_and_diagnostic_order() {
        let dir =
            std::env::temp_dir().join(format!("mapping-scoring-threads-{}", std::process::id()));
        for name in ["one", "four"] {
            fs::create_dir_all(dir.join(name)).unwrap();
        }
        let mut seed = 73u64;
        let backbone: Vec<_> = (0..500)
            .map(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(seed >> 32) as usize % 4]
            })
            .collect();
        let reads = Reads {
            r1: dir.join("r1.fq"),
            r2: dir.join("r2.fq"),
        };
        let mut one = File::create(&reads.r1).unwrap();
        let mut two = File::create(&reads.r2).unwrap();
        for i in 0..48 {
            // Repeated names deliberately have later, different observations.
            let name = if i >= 40 { i - 40 } else { i };
            let mut a = backbone[120 + i % 20..180 + i % 20].to_vec();
            let mut b = kmer::reverse_complement(&backbone[320 + i % 20..380 + i % 20]);
            if i % 5 == 0 {
                a.insert(39, b'A');
            }
            if i % 7 == 0 {
                b[30] = b'N';
            }
            if i % 11 == 0 {
                a.fill(b'N');
                b.fill(b'N');
            }
            for (mate, sequence, output) in [(1, a, &mut one), (2, b, &mut two)] {
                let quality = (0..sequence.len()).map(|p| b"I:F+"[(i + p) % 4]).collect();
                fastq::Record {
                    header: format!("@fragment-{name}/{mate}"),
                    sequence,
                    quality,
                    plus: "+".into(),
                }
                .write(output)
                .unwrap();
            }
        }
        let Commands::AssembleMapping(mut args) = Cli::try_parse_from([
            "mitovarr",
            "assemble-mapping",
            "--ref",
            "unused",
            "--bam",
            "unused",
            "--out",
            "unused",
        ])
        .unwrap()
        .command
        else {
            panic!()
        };
        let mut alternate = backbone[150..350].to_vec();
        alternate[85] = if alternate[85] == b'A' { b'C' } else { b'A' };
        let alternatives = [
            Alternative {
                sequence: backbone[150..350].to_vec(),
                support: 10,
            },
            Alternative {
                sequence: alternate,
                support: 10,
            },
        ];
        let mut measured = Vec::new();
        for (threads, name) in [(1, "one"), (4, "four")] {
            args.scoring_threads = threads;
            let diagnostic = dir.join(name).join("fragments.jsonl.gz");
            measured.push(
                evaluate(
                    &backbone,
                    150,
                    350,
                    &alternatives,
                    &reads,
                    None,
                    61,
                    &args,
                    Some(&diagnostic),
                    None,
                )
                .unwrap(),
            );
        }
        assert!(measured[0].report["common_fragments"].as_u64().unwrap() > 16);
        assert_eq!(measured[0].selected, measured[1].selected);
        assert_eq!(measured[0].report, measured[1].report);
        assert_eq!(
            fs::read(dir.join("one/fragments.jsonl.gz")).unwrap(),
            fs::read(dir.join("four/fragments.jsonl.gz")).unwrap()
        );
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn cohort_eligibility_is_frozen_before_the_shortlist() {
        let dir = std::env::temp_dir().join(format!("frozen-cohort-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut state = 1197183u64;
        let backbone: Vec<_> = (0..500)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state % 4) as usize]
            })
            .collect();
        let reads = Reads {
            r1: dir.join("one.fq"),
            r2: dir.join("two.fq"),
        };
        for (mate, path) in [(1, &reads.r1), (2, &reads.r2)] {
            fs::write(
                path,
                format!(
                    "@pair/{mate}\n{}\n+\n{}\n",
                    String::from_utf8_lossy(&backbone[..40]),
                    "I".repeat(40)
                ),
            )
            .unwrap();
        }
        let Commands::AssembleMapping(args) = Cli::try_parse_from([
            "mitovarr",
            "assemble-mapping",
            "--ref",
            "unused",
            "--bam",
            "unused",
            "--out",
            "unused",
        ])
        .unwrap()
        .command
        else {
            panic!()
        };
        let mut paths = vec![Alternative {
            sequence: backbone[150..350].to_vec(),
            support: 10,
        }];
        let duplicate = [&backbone[150..250], &backbone[..70], &backbone[250..350]].concat();
        let frozen = [paths[0].sequence.clone(), duplicate.clone()];
        let frozen: Vec<_> = frozen.iter().map(Vec::as_slice).collect();
        let original = evaluate(
            &backbone, 150, 350, &paths, &reads, None, 40, &args, None, None,
        )
        .unwrap();
        let subset = evaluate(
            &backbone,
            150,
            350,
            &paths,
            &reads,
            None,
            40,
            &args,
            None,
            Some(&frozen),
        )
        .unwrap();
        paths.push(Alternative {
            sequence: duplicate,
            support: 10,
        });
        let expanded = evaluate(
            &backbone,
            150,
            350,
            &paths,
            &reads,
            None,
            40,
            &args,
            None,
            Some(&frozen),
        )
        .unwrap();
        assert_eq!(original.report["common_fragments"], 1);
        assert_eq!(subset.report["common_fragments"], 0);
        assert_eq!(expanded.report["common_fragments"], 0);
        assert_eq!(
            subset.report["shared_unique_anchor_kmers"],
            expanded.report["shared_unique_anchor_kmers"]
        );
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn masked_repeat_contraction_needs_linked_unique_flanks() {
        let dir = std::env::temp_dir().join(format!("contraction-scoring-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut state = 498723u64;
        let mut dna = |n| {
            (0..n)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state >> 32) as usize & 3]
                })
                .collect::<Vec<_>>()
        };
        let left = dna(80);
        let repeat = dna(25).repeat(4);
        let right = dna(80);
        let short = [left, repeat, right].concat();
        let baseline = [&short[..100], &[b'N'; 170], &short[100..]].concat();
        let alternatives = vec![
            Alternative {
                sequence: baseline.clone(),
                support: 0,
            },
            Alternative {
                sequence: short.clone(),
                support: 20,
            },
        ];
        let reads = Reads {
            r1: dir.join("r1.fq"),
            r2: dir.join("r2.fq"),
        };
        let Commands::AssembleMapping(args) = Cli::try_parse_from([
            "mitovarr",
            "assemble-mapping",
            "--ref",
            "unused",
            "--bam",
            "unused",
            "--out",
            "unused",
        ])
        .unwrap()
        .command
        else {
            panic!()
        };
        for linked in [false, true] {
            let mut one = fs::File::create(&reads.r1).unwrap();
            let mut two = fs::File::create(&reads.r2).unwrap();
            for id in 0..12 {
                // Left-side fragments cross BOTH local edit coordinates (100),
                // but cannot establish the repeat's distant right exit.
                let read = if linked {
                    &short[55..215]
                } else if id < 6 {
                    &short[55..135]
                } else {
                    &short[135..215]
                };
                for (mate, out, sequence) in [
                    (1, &mut one, read.to_vec()),
                    (2, &mut two, kmer::reverse_complement(read)),
                ] {
                    writeln!(
                        out,
                        "@p{id}/{mate}\n{}\n+\n{}",
                        String::from_utf8_lossy(&sequence),
                        "I".repeat(sequence.len())
                    )
                    .unwrap();
                }
            }
            drop(one);
            drop(two);
            let evaluation = evaluate(
                &baseline,
                0,
                baseline.len(),
                &alternatives,
                &reads,
                None,
                160,
                &args,
                None,
                None,
            )
            .unwrap();
            assert_eq!(
                evaluation.report["best_candidate"], 1,
                "{}",
                evaluation.report
            );
            assert!(evaluation.report["candidates"][1]["junction_fragments"]
                .as_array()
                .unwrap()
                .iter()
                .all(|v| v.as_u64().unwrap() >= 3));
            if linked {
                assert_eq!(evaluation.selected, Some(1), "{}", evaluation.report);
            } else {
                assert_eq!(evaluation.selected, None, "{}", evaluation.report);
                assert_eq!(
                    evaluation.report["reason"],
                    "contraction_unique_flanks_not_linked"
                );
                assert_eq!(
                    evaluation.report["candidates"][1]["contraction_flank_fragments"],
                    0
                );
            }
        }
        fs::remove_dir_all(dir).unwrap();
    }
}
