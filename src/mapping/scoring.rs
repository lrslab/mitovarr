use super::Alternative;
use crate::{
    cli::{PolishArgs, Reads},
    denovo::{
        align::{self, Hit, Index, Placements, Prepared},
        pairs::Library,
    },
    io::fastq,
    kmer,
};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashMap, HashSet};

pub(super) struct Evaluation {
    pub selected: Option<usize>,
    pub report: Value,
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

fn crosses(hit: &Hit, junction: usize) -> bool {
    hit.start + 20 <= junction && hit.end >= junction + 20
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

/// A repeated hit never increases a fragment's mass: retain only the best
/// sequence explanation within each library span state, as in de novo scoring.
#[allow(clippy::too_many_arguments)]
fn joint(
    a: &Placements,
    b: &Placements,
    model: Option<&Library>,
    ideal: f64,
    length: usize,
    junctions: [usize; 2],
    contraction: Option<[usize; 2]>,
) -> (f64, [bool; 2], bool) {
    let mut states = BTreeMap::<usize, f64>::new();
    let mut best = 0.0f64;
    let mut bridge = [0.0f64; 2];
    let mut contraction_bridge = 0.0f64;
    for x in &a.hits {
        for original_y in &b.hits {
            // Include a pair crossing the arbitrary FASTA origin. Moving one
            // end by a whole circle is geometry only, never an extra observation.
            for shift in [0, length] {
                let mut x = x.clone();
                let mut y = original_y.clone();
                if shift > 0 {
                    let h = if x.start < y.start { &mut x } else { &mut y };
                    h.start += shift;
                    h.end += shift;
                }
                if x.end.max(y.end) - x.start.min(y.start) > length {
                    continue;
                }
                let span = match model {
                    Some(m) => match m.span(&x, &y) {
                        Some(d) => d,
                        None => continue,
                    },
                    None => 0,
                };
                let empirical = model.filter(|m| m.geometry_usable);
                let weight = empirical.map_or(1.0, |m| m.probability(span));
                let sequence = (x.logp + y.logp - ideal).exp().min(1.0);
                states
                    .entry(if empirical.is_some() { span } else { 0 })
                    .and_modify(|v| *v = v.max(sequence))
                    .or_insert(sequence);
                let relative = sequence * empirical.map_or(1.0, |m| weight / m.peak_probability());
                best = best.max(relative);
                for (j, &junction) in junctions.iter().enumerate() {
                    let paired = empirical.is_some_and(|m| span <= m.upper)
                        && x.start.min(y.start) + 20 <= junction
                        && x.end.max(y.end) >= junction + 20;
                    if crosses(&x, junction) || crosses(&y, junction) || paired {
                        bridge[j] = bridge[j].max(relative);
                    }
                }
                if let Some(anchors) = contraction {
                    let contains = |h: &Hit, p: usize| h.start <= p && h.end >= p + 21;
                    let one_read = [&x, &y]
                        .iter()
                        .any(|h| anchors.iter().all(|&p| contains(h, p)));
                    let paired = empirical.is_some_and(|m| span <= m.upper)
                        && anchors.iter().all(|&p| contains(&x, p) || contains(&y, p));
                    if one_read || paired {
                        contraction_bridge = contraction_bridge.max(relative);
                    }
                }
            }
        }
    }
    let mass: f64 = states
        .into_iter()
        .map(|(span, probability)| {
            probability
                * model
                    .filter(|m| m.geometry_usable)
                    .map_or(1.0, |m| m.probability(span))
        })
        .sum();
    (
        (mass.min(1.0) + (-12.0f64).exp()).ln(),
        bridge.map(|p| p >= 0.1 && p >= 0.99 * best),
        contraction_bridge >= 0.1 && contraction_bridge >= 0.99 * best,
    )
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
    args: &PolishArgs,
) -> Result<Evaluation> {
    let mut pool_seeds = HashSet::new();
    for candidate in alternatives {
        kmer::visit(&candidate.sequence, None, 21, 0, |_, c| {
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
    for sequence in &sequences {
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
    let mut votes = vec![vec![0usize; n]; n];
    let mut bridges = vec![[0usize; 2]; n];
    let mut contraction_bridges = vec![0usize; n];
    let mut names = HashSet::new();
    let (mut fragments, mut unknown) = (0usize, 0usize);
    fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
        if !names.insert(a.name(1)?.to_owned()) {
            return Ok(());
        }
        ensure!(
            names.len() <= args.max_records,
            "local scoring exceeds --max-records named FASTQ pairs"
        );
        let mut relevant = false;
        let mut observed_anchors = HashSet::new();
        for read in [a, b] {
            kmer::visit(
                &read.sequence,
                Some(&read.quality),
                21,
                args.min_baseq,
                |_, c| {
                    let c = kmer::canonical(c, 21);
                    if pool_seeds.contains(&c) {
                        relevant = true;
                        observed_anchors.insert(c);
                    }
                },
            );
        }
        if !relevant {
            return Ok(());
        }
        fragments += 1;
        let ideal = align::ideal_score(a) + align::ideal_score(b);
        let mut one = Prepared::new(a);
        let mut two = Prepared::new(b);
        let mut results = Vec::with_capacity(n);
        for (i, index) in indexes.iter().enumerate() {
            let mut x = index.align_prepared(&mut one);
            let mut y = index.align_prepared(&mut two);
            // Rescue a mate without a seed only near a uniquely anchored partner;
            // apply the same operation to baseline and every alternative.
            if let Some(m) = model.filter(|m| m.geometry_usable) {
                for swap in [false, true] {
                    let (pivot, missing, prepared) = if swap {
                        (&y, &mut x, &mut one)
                    } else {
                        (&x, &mut y, &mut two)
                    };
                    if pivot.hits.len() == 1 && missing.hits.is_empty() && !missing.exhausted {
                        let p = &pivot.hits[0];
                        *missing = index.rescue_near(
                            prepared,
                            p.start.saturating_sub(m.upper),
                            p.end.saturating_add(m.upper),
                            !p.reverse,
                        );
                    }
                }
            }
            if x.exhausted || y.exhausted {
                unknown += 1;
                return Ok(());
            }
            results.push(joint(
                &x,
                &y,
                model,
                ideal,
                sequences[i].len(),
                junctions[i],
                flanks[i],
            ));
        }
        for i in 0..n {
            scores[i] += results[i].0;
            if results[i].2
                && flanks[i].is_some_and(|anchors| {
                    anchors.iter().all(|&p| {
                        let c = super::graph::code(&sequences[i][p..p + 21]).unwrap();
                        observed_anchors.contains(&kmer::canonical(c, 21))
                    })
                })
            {
                contraction_bridges[i] += 1;
            }
            for j in 0..2 {
                bridges[i][j] += usize::from(results[i].1[j]);
            }
            for j in 0..n {
                votes[i][j] += usize::from(results[i].0 - results[j].0 >= 2.3);
            }
        }
        Ok(())
    })?;
    let best = (0..n)
        .max_by(|&a, &b| scores[a].total_cmp(&scores[b]).then(b.cmp(&a)))
        .unwrap();
    let fraction = args.min_alt_frac.unwrap_or(args.min_indel_frac);
    let discriminated = (0..n).filter(|&i| i != best).all(|i| {
        votes[best][i] >= 3
            && votes[best][i] as f64 >= fraction * (votes[best][i] + votes[i][best]) as f64
            && scores[best] - scores[i] > 5.0 + 12.0 * unknown as f64
    });
    let supported = bridges[best].iter().all(|&c| c >= 3)
        && alternatives[best].support >= args.min_depth
        && !alternatives[best].sequence.contains(&b'N');
    let contraction_supported = !requires_contraction[best] || contraction_bridges[best] >= 3;
    let selected =
        (discriminated && supported && contraction_supported && unknown == 0).then_some(best);
    Ok(Evaluation {
        selected,
        report: json!({"scoring_status":"evaluated","common_fragments":fragments,"shared_unique_anchor_kmers":pool_seeds.len(),"cohort":"at least one mate contains a quality-filtered local 21-mer unique in every whole-contig candidate; both mates jointly scored","unscored_fragments":unknown,"best_candidate":best,
            "selected_candidate":selected,"status":if selected.is_some(){"resolved"}else{"unresolved"},
            "reason":if !discriminated{"competing_paths"}else if !supported{"insufficient_path_or_junction_support"}else if !contraction_supported{"contraction_unique_flanks_not_linked"}else if unknown>0{"alignment_search_limit"}else{"supported"},
            "candidates":alternatives.iter().enumerate().map(|(i,p)|json!({"candidate":i,"length":p.sequence.len(),"minimum_kmer_fragments":p.support,"score":scores[i],"strong_votes_against":votes[i],"junction_fragments":bridges[i],"junction_positions_rotated_baseline_start":junctions[i],"requires_contraction_flanks":requires_contraction[i],"contraction_flank_starts_0based":flanks[i],"contraction_flank_fragments":contraction_bridges[i]})).collect::<Vec<_>>(),
            "model":"one named fragment once; full-query affine alignment; best explanation per empirical span; original copy retained in every whole-contig candidate",
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
        let Commands::Polish(args) = Cli::try_parse_from([
            "mitovarr", "polish", "--ref", "unused", "--bam", "unused", "--out", "unused",
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
