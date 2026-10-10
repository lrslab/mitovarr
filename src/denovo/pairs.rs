//! Estimate an empirical fragment model from sample-only, unambiguous sequence.
use super::{
    align::{Hit, Index},
    normalize,
    threading::PathGraph,
};
use crate::{
    cli::{GraphSettings, Limits, Reads},
    io::fastq::{self, Record},
    kmer,
};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::collections::{BTreeSet, HashMap, HashSet};

#[derive(Debug)]
pub(crate) struct Library {
    pub inward: bool,
    pub geometry_usable: bool,
    peak: f64,
    pub upper: usize,
    probabilities: Vec<f64>,
}

impl Library {
    #[cfg(test)]
    pub(crate) fn test_model(inward: bool, geometry_usable: bool, probabilities: Vec<f64>) -> Self {
        Self {
            inward,
            geometry_usable,
            peak: probabilities.iter().copied().fold(0.0f64, f64::max),
            upper: probabilities.len().saturating_sub(1),
            probabilities,
        }
    }

    pub fn peak_probability(&self) -> f64 {
        self.peak
    }
    pub fn probability(&self, span: usize) -> f64 {
        self.probabilities.get(span).copied().unwrap_or(0.0)
    }
    pub fn span(&self, a: &Hit, b: &Hit) -> Option<usize> {
        self.span_coordinates((a.start, a.end, a.reverse), (b.start, b.end, b.reverse))
    }
    pub fn span_coordinates(
        &self,
        a: (usize, usize, bool),
        b: (usize, usize, bool),
    ) -> Option<usize> {
        if a.2 == b.2 {
            return None;
        }
        let (forward, reverse) = if a.2 { (b, a) } else { (a, b) };
        let inward = forward.0 < reverse.1;
        (inward == self.inward).then_some(a.1.max(b.1) - a.0.min(b.0))
    }
}

pub(crate) struct Training {
    pub model: Option<Library>,
    pub max_read: usize,
    pub report: Value,
}

struct ExactIndex {
    sequences: Vec<Vec<u8>>,
    seeds: HashMap<kmer::Code, Option<(usize, usize)>>,
}
impl ExactIndex {
    fn new(sequences: Vec<Vec<u8>>) -> Self {
        let mut seeds = HashMap::new();
        for (id, sequence) in sequences.iter().enumerate() {
            kmer::visit(sequence, None, 31, 0, |pos, code| {
                seeds
                    .entry(code)
                    .and_modify(|hit| *hit = None)
                    .or_insert(Some((id, pos)));
            });
        }
        Self { sequences, seeds }
    }
    fn place(&self, read: &Record, minq: u8) -> Option<(usize, Hit)> {
        let mut found = BTreeSet::new();
        for reverse in [false, true] {
            let query = if reverse {
                kmer::reverse_complement(&read.sequence)
            } else {
                read.sequence.iter().map(u8::to_ascii_uppercase).collect()
            };
            let quality: Vec<_> = if reverse {
                read.quality.iter().rev().copied().collect()
            } else {
                read.quality.clone()
            };
            let mut proposed = BTreeSet::new();
            kmer::visit(&query, Some(&quality), 31, minq, |q, code| {
                if let Some(Some((id, p))) = self.seeds.get(&code) {
                    if *p >= q {
                        proposed.insert((*id, p - q));
                    }
                }
            });
            for (id, start) in proposed {
                if self.sequences[id].get(start..start + query.len()) == Some(query.as_slice()) {
                    found.insert((id, start, reverse));
                }
            }
        }
        if found.len() != 1 {
            return None;
        }
        let (id, start, reverse) = found.into_iter().next()?;
        Some((
            id,
            Hit {
                start,
                end: start + read.sequence.len(),
                reverse,
                logp: 0.0,
                edits: 0,
                traceback: None,
            },
        ))
    }
}

pub(super) fn estimate(
    graph: &PathGraph,
    reads: &Reads,
    settings: &GraphSettings,
    limits: &Limits,
) -> Result<Training> {
    let sequences: BTreeSet<_> = graph
        .segments
        .iter()
        .filter(|s| s.sequence.len() >= 31)
        .map(|s| normalize(&s.sequence, s.circular))
        .collect();
    estimate_sequences(sequences.into_iter().collect(), reads, settings, limits)
}

/// The caller supplies training sequences; mapping uses the unambiguous pileup
/// backbone while de novo continues to use sample-only graph segments.
pub(crate) fn estimate_sequences(
    sequences: Vec<Vec<u8>>,
    reads: &Reads,
    settings: &GraphSettings,
    limits: &Limits,
) -> Result<Training> {
    let index = ExactIndex::new(sequences);
    let mut names = HashSet::new();
    let mut observations = Vec::new();
    let mut orientation_counts = [0usize; 3];
    let mut max_read = 0;
    fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
        max_read = max_read.max(a.sequence.len()).max(b.sequence.len());
        let Some((x, p)) = index.place(a, limits.min_baseq) else {
            return Ok(());
        };
        let Some((y, q)) = index.place(b, limits.min_baseq) else {
            return Ok(());
        };
        if x != y {
            return Ok(());
        }
        let name = a.name(1)?;
        if names.contains(name) {
            return Ok(());
        }
        ensure!(
            names.len() < settings.max_validation_pairs,
            "fragment calibration exceeds --max-validation-pairs"
        );
        names.insert(name.to_owned());
        let (forward, reverse) = if p.reverse { (&q, &p) } else { (&p, &q) };
        let layout = if p.reverse == q.reverse {
            2
        } else {
            usize::from(forward.start >= reverse.end)
        };
        orientation_counts[layout] += 1;
        if layout < 2 {
            let span = p.end.max(q.end) - p.start.min(q.start);
            observations.push((span, layout, index.sequences[x].len()));
        }
        Ok(())
    })?;
    let layout = usize::from(orientation_counts[1] > orientation_counts[0]);
    let chosen = orientation_counts[layout];
    if chosen < 100 || chosen * 10 < names.len() * 9 {
        return Ok(Training {
            model: None,
            max_read,
            report: json!({"status":"not_estimable","unique_same_segment_fragments":names.len(),"orientation_counts_fr_rf_same":orientation_counts,"minimum_training_fragments":100,"reason":"insufficient unique pairs or no dominant orientation; distance evidence disabled"}),
        });
    }
    let mut spans: Vec<_> = observations
        .iter()
        .filter(|o| o.1 == layout)
        .map(|o| o.0)
        .collect();
    spans.sort_unstable();
    let quantile = |f: f64| spans[((spans.len() - 1) as f64 * f).round() as usize];
    let pilot_median = quantile(0.5);
    let upper = quantile(0.999).saturating_add(100).max(max_read);
    // Freeze an R1-pivot cohort with a common observable mate range. Inclusion
    // must not depend on whether R2 maps successfully or on its measured span.
    let aligned_indexes: Vec<_> = index
        .sequences
        .iter()
        .map(|s| Index::new(s, false, max_read))
        .collect();
    let mut cohort = HashSet::new();
    let mut measured = Vec::new();
    let mut windows: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    let mut unmeasured = 0;
    let mut missing_reasons = std::collections::BTreeMap::<&str, usize>::new();
    fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
        let Some((id, p)) = index.place(a, limits.min_baseq) else {
            return Ok(());
        };
        let length = index.sequences[id].len();
        let right = (layout == 0) != p.reverse;
        if p.start < 20
            || p.end + 20 > length
            || (right && p.start + upper + 20 > length)
            || (!right && p.end < upper + 20)
        {
            return Ok(());
        }
        let name = a.name(1)?;
        if cohort.contains(name) {
            return Ok(());
        }
        ensure!(
            cohort.len() < settings.max_validation_pairs,
            "fragment calibration exceeds --max-validation-pairs"
        );
        cohort.insert(name.to_owned());
        let q = match index.place(b, limits.min_baseq) {
            Some((other, q)) if other == id => Some(q),
            Some(_) => {
                *missing_reasons.entry("exact_other_segment").or_default() += 1;
                None
            }
            None => {
                let mut targets = BTreeSet::new();
                for query in [b.sequence.clone(), kmer::reverse_complement(&b.sequence)] {
                    kmer::visit(&query, None, 31, 0, |_, code| {
                        if let Some(Some((target, _))) = index.seeds.get(&code) {
                            targets.insert(*target);
                        }
                    });
                }
                if targets.is_empty() || targets.len() > 64 {
                    *missing_reasons
                        .entry(if targets.is_empty() {
                            "no_unique_r2_seed"
                        } else {
                            "too_many_r2_seed_targets"
                        })
                        .or_default() += 1;
                    None
                } else {
                    // A seed in another fragment is not a competing full-read
                    // placement. Compare accepted gapped placements before
                    // deciding that an R2 with several seed targets is ambiguous.
                    let mut hits = Vec::new();
                    let mut exhausted = false;
                    for target in targets {
                        let found = aligned_indexes[target].align(b);
                        exhausted |= found.exhausted;
                        hits.extend(found.hits.into_iter().map(|hit| (target, hit)));
                    }
                    hits.sort_by(|a, b| {
                        b.1.logp
                            .total_cmp(&a.1.logp)
                            .then(a.0.cmp(&b.0))
                            .then(a.1.start.cmp(&b.1.start))
                    });
                    if !exhausted
                        && !hits.is_empty()
                        && hits[0].0 == id
                        && (hits.len() == 1 || hits[0].1.logp - hits[1].1.logp >= 5.0)
                    {
                        Some(hits.remove(0).1)
                    } else {
                        *missing_reasons
                            .entry(if exhausted {
                                "alignment_limit"
                            } else if hits.is_empty() {
                                "no_accepted_gapped_alignment"
                            } else if hits[0].0 != id {
                                "gapped_other_segment"
                            } else {
                                "ambiguous_gapped_alignment"
                            })
                            .or_default() += 1;
                        None
                    }
                }
            }
        };
        let Some(q) = q else {
            unmeasured += 1;
            return Ok(());
        };
        let (forward, reverse) = if p.reverse { (&q, &p) } else { (&p, &q) };
        let orientation = if p.reverse == q.reverse {
            2
        } else {
            usize::from(forward.start >= reverse.end)
        };
        let span = p.end.max(q.end) - p.start.min(q.start);
        if orientation != layout || span > upper {
            *missing_reasons
                .entry(if orientation != layout {
                    "other_orientation"
                } else {
                    "span_outside_search_range"
                })
                .or_default() += 1;
            unmeasured += 1;
            return Ok(());
        }
        measured.push(span);
        windows
            .entry((id, p.start / upper.max(1)))
            .or_default()
            .push(span);
        Ok(())
    })?;
    if measured.len() < 100 {
        return Ok(Training {
            model: None,
            max_read,
            report: json!({"status":"not_estimable","reason":"too few measured fragments in the fixed R1 internal cohort; no default distance substituted","pilot_median_span":pilot_median,"r1_pivot_cohort":cohort.len(),"measured_fragments":measured.len(),"unmeasured_fragments":unmeasured,"search_span_limit":upper}),
        });
    }
    measured.sort_unstable();
    let quantile = |f: f64| measured[((measured.len() - 1) as f64 * f).round() as usize];
    let mut mass = vec![0.0; upper + 1];
    for &span in &measured {
        let lo = span.saturating_sub(10).max(1);
        let hi = (span + 10).min(upper);
        for value in &mut mass[lo..=hi] {
            *value += 1.0 / (hi - lo + 1) as f64 / cohort.len() as f64;
        }
    }
    let mut window_summary:Vec<_>=windows.into_iter().map(|((segment,window),mut spans)| {
        spans.sort_unstable();
        json!({"segment":segment+1,"window":window,"fragments":spans.len(),"median_span":spans[spans.len()/2]})
    }).collect();
    window_summary.sort_by_key(|w| (w["segment"].as_u64(), w["window"].as_u64()));
    let unknown_fraction = unmeasured as f64 / cohort.len() as f64;
    let mut window_medians: Vec<_> = window_summary
        .iter()
        .filter(|w| w["fragments"].as_u64().unwrap_or(0) >= 100)
        .map(|w| w["median_span"].as_u64().unwrap() as usize)
        .collect();
    window_medians.sort_unstable();
    let relative_window_spread = if window_medians.len() >= 3 {
        let lo = window_medians[(window_medians.len() - 1) / 10];
        let hi = window_medians[(window_medians.len() - 1) * 9 / 10];
        (hi - lo) as f64 / quantile(0.5).max(1) as f64
    } else {
        0.0
    };
    let geometry_usable = unknown_fraction <= 0.05 && relative_window_spread <= 0.25;
    let report = json!({"status":"estimated","geometry_use":if geometry_usable{"empirical_compatibility"}else{"diagnostic_only"},"geometry_eligibility":{"maximum_unmeasured_fraction":0.05,"maximum_relative_window_median_spread":0.25,"observed_relative_window_median_spread":relative_window_spread},"quantile_scope":"successfully measured subset, not the entire R1 cohort","mass_scope":"P(span and successfully measured | qualifying R1), not a missingness-corrected population PMF","method":"fixed, unique R1 pivot inside a common observable interval; include unsuccessful R2 in denominator; full-query R2 gapped rescue; 21-bp smoothed empirical mass","pilot_unique_pairs":names.len(),"pilot_median_span":pilot_median,"r1_pivot_cohort":cohort.len(),"training_fragments":measured.len(),"unmeasured_fragments":unmeasured,"unmeasured_fraction":unmeasured as f64/cohort.len() as f64,"unmeasured_reasons":missing_reasons,"orientation":if layout==0{"FR"}else{"RF"},"orientation_counts_fr_rf_same":orientation_counts,"median_span":quantile(0.5),"q01_span":quantile(0.01),"q99_span":quantile(0.99),"search_span_limit":upper,"max_read_length":max_read,"maximum_gapped_calibration_seed_targets":64,"internal_windows":window_summary,"limitations":"conditional on extracted FASTQ and uniquely anchored R1; pooled library without read-group identity; unmeasured mass is retained, not renormalized into successful pairs"});
    Ok(Training {
        model: Some(Library {
            inward: layout == 0,
            geometry_usable,
            peak: mass.iter().copied().fold(0.0, f64::max),
            upper,
            probabilities: mass,
        }),
        max_read,
        report,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn same_start_alternative_end_can_win_the_joint_score() {
        use super::super::{align, evidence};
        let read = Record {
            header: "@end/1".into(),
            sequence: b"ACGTA".to_vec(),
            quality: vec![b'I'; 5],
            plus: "+".into(),
        };
        let index = Index::new(b"ACGTAG", false, 5);
        let mut prepared = align::Prepared::new(&read);
        let mut placements = index.profile_all(&mut prepared);
        assert!(placements.complete);
        placements.hits.retain(|h| h.start == 0 && !h.reverse);
        let observations = evidence::observations(&placements, &read, &index, 0);
        let best = observations
            .iter()
            .max_by(|a, b| a.0.logp.total_cmp(&b.0.logp))
            .unwrap();
        assert_eq!(best.0.end, 5);
        assert!(observations
            .iter()
            .any(|h| h.0.end == 6 && h.0.logp < best.0.logp));
        let mate = [(
            Hit {
                start: 1,
                end: 4,
                reverse: true,
                logp: 0.0,
                edits: 0,
                traceback: None,
            },
            align::BaseTrace::default(),
        )];
        let model = Library {
            inward: true,
            geometry_usable: true,
            peak: 1.0,
            upper: 6,
            probabilities: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.000001, 0.999999],
        };
        let score = |hits: &[(Hit, align::BaseTrace)]| {
            evidence::joint(
                hits,
                &mate,
                Some(&model),
                align::ideal_score(&read),
                0,
                [0, 0],
                None,
            )
        };
        let complete = score(&observations);
        let sequence_only = score(std::slice::from_ref(best));
        assert!(complete.score > sequence_only.score + 2.0);
        assert_eq!(
            complete.diagnostics["span_states"]
                .as_array()
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn mate_gap_is_a_bracket_but_never_an_observed_junction() {
        use super::super::{align, evidence};
        let mut seed = 81829u64;
        let sequence: Vec<_> = (0..220)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed % 4) as usize]
            })
            .collect();
        let record = |sequence: Vec<u8>| Record {
            header: "@fragment".into(),
            quality: vec![b'I'; sequence.len()],
            sequence,
            plus: "+".into(),
        };
        let a = record(sequence[10..60].to_vec());
        let b = record(kmer::reverse_complement(&sequence[140..190]));
        let index = Index::new(&sequence, false, 50);
        let mut probabilities = vec![0.0; 181];
        probabilities[180] = 1.0;
        let model = Library {
            inward: true,
            geometry_usable: true,
            peak: 1.0,
            upper: 180,
            probabilities,
        };
        let result = evidence::joint(
            &evidence::observations(&index.align(&a), &a, &index, 0),
            &evidence::observations(&index.align(&b), &b, &index, 0),
            Some(&model),
            align::ideal_score(&a) + align::ideal_score(&b),
            0,
            [100, 100],
            Some([20, 150]),
        );
        assert!(result.score > -1.0);
        assert_eq!(result.observed_junctions, [false, false]);
        assert_eq!(result.mate_brackets, [true, true]);
        assert!(result.contraction); // The two flanks, unlike the gap, were sequenced.
    }
    #[test]
    fn whole_hypothesis_rescue_keeps_seeded_decoy_and_is_rotation_invariant() {
        check_seedless_rescue(false);
    }
    #[test]
    fn sequence_best_decoy_in_geometry_tail_does_not_suppress_better_joint_explanation() {
        check_seedless_rescue(true);
    }
    fn check_seedless_rescue(tail: bool) {
        use super::super::{
            align::{self, Prepared},
            evidence,
        };
        let mut seed = 129173u64;
        let mut sequence: Vec<_> = (0..1000)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed % 4) as usize]
            })
            .collect();
        let a = Record {
            header: "@test/1".into(),
            sequence: sequence[850..930].to_vec(),
            plus: "+".into(),
            quality: vec![b'I'; 80],
        };
        let mut query: Vec<_> = (970..1050).map(|p| sequence[p % 1000]).collect();
        let mut quality = vec![b'I'; 80];
        for p in (6..80).step_by(14) {
            query[p] = if query[p] == b'A' { b'C' } else { b'A' };
            quality[79 - p] = if tail { b'$' } else { b'!' };
        }
        sequence[400..480].copy_from_slice(&query);
        let b = Record {
            header: "@test/2".into(),
            sequence: kmer::reverse_complement(&query),
            plus: "+".into(),
            quality,
        };
        let upper = if tail { 650 } else { 220 };
        let mut probabilities = vec![0.0; upper + 1];
        probabilities[200] = if tail { 0.999999 } else { 1.0 };
        if tail {
            probabilities[630] = 0.000001;
        }
        let model = Library {
            inward: true,
            geometry_usable: true,
            peak: 1.0,
            upper,
            probabilities,
        };
        let mut scores = Vec::new();
        for reverse in [false, true] {
            let oriented = if reverse {
                kmer::reverse_complement(&sequence)
            } else {
                sequence.clone()
            };
            for offset in [0, 30, 400, 900] {
                let rotated = [&oriented[offset..], &oriented[..offset]].concat();
                let indexes = [Index::new(&rotated, true, 80)];
                let mut one = Prepared::new(&a);
                let mut two = Prepared::new(&b);
                let mut placements = vec![(
                    indexes[0].align_prepared(&mut one),
                    indexes[0].align_prepared(&mut two),
                )];
                assert!(!placements[0].1.hits.is_empty());
                let before = placements[0].1.hits.clone();
                let before_score = evidence::joint(
                    &evidence::observations(&placements[0].0, &a, &indexes[0], 0),
                    &evidence::observations(&placements[0].1, &b, &indexes[0], 0),
                    Some(&model),
                    align::ideal_score(&a) + align::ideal_score(&b),
                    1000,
                    [0, 0],
                    None,
                )
                .score;
                let seeds = {
                    let mut seeds = HashSet::new();
                    kmer::visit(&a.sequence, None, 21, 0, |_, c| {
                        seeds.insert(kmer::canonical(c, 21));
                    });
                    seeds
                };
                let pivot = super::super::anchor::Pivot::choose([&a, &b], &seeds, 0).unwrap();
                assert_eq!(pivot.mate, 0);
                let locations = super::super::anchor::locations(&indexes[0], &seeds);
                pivot.constrain(
                    &mut placements[0].0,
                    &a,
                    &mut one,
                    &indexes[0],
                    &locations,
                    0,
                );
                let rescued = evidence::rescue_pool(
                    &indexes,
                    &mut placements,
                    &mut one,
                    &mut two,
                    Some(&model),
                    0,
                );
                assert!(rescued.seeded_decoy_searches > 0);
                assert!(
                    rescued.added_placements > 0,
                    "rotation {offset}, reverse {reverse}"
                );
                assert!(
                    placements[0].0.complete && placements[0].1.complete,
                    "{:?}",
                    rescued.states
                );
                for hit in before.into_iter().filter(|_| tail) {
                    assert!(placements[0]
                        .1
                        .hits
                        .iter()
                        .any(|h| h.start == hit.start && h.reverse == hit.reverse));
                }
                let result = evidence::joint(
                    &evidence::observations(&placements[0].0, &a, &indexes[0], 0),
                    &evidence::observations(&placements[0].1, &b, &indexes[0], 0),
                    Some(&model),
                    align::ideal_score(&a) + align::ideal_score(&b),
                    1000,
                    [0, 0],
                    None,
                );
                if tail {
                    assert!(
                        result.score > before_score + 2.0,
                        "{} -> {}",
                        before_score,
                        result.score
                    );
                }
                scores.push(result.score);
            }
        }
        assert!(
            scores.iter().all(|s| (*s - scores[0]).abs() < 1e-8),
            "{scores:?}"
        );
        assert!(scores[0] > -12.0);
    }
    #[test]
    fn overlapping_mates_use_outer_span_and_orientation() {
        let library = Library {
            inward: true,
            geometry_usable: true,
            peak: 1.0,
            upper: 300,
            probabilities: vec![0.0; 301],
        };
        let a = Hit {
            start: 10,
            end: 110,
            reverse: false,
            logp: 0.0,
            edits: 0,
            traceback: None,
        };
        let b = Hit {
            start: 80,
            end: 180,
            reverse: true,
            logp: 0.0,
            edits: 0,
            traceback: None,
        };
        assert_eq!(library.span(&a, &b), Some(170));
        assert_eq!(library.span(&b, &a), Some(170));
        assert_eq!(
            library.span(
                &a,
                &Hit {
                    reverse: false,
                    ..b.clone()
                }
            ),
            None
        );
        assert_eq!(
            library.span(
                &a,
                &Hit {
                    start: 0,
                    end: 5,
                    ..b
                }
            ),
            None
        );
    }
    #[test]
    fn a_short_seed_decoy_does_not_make_a_full_read_ambiguous() {
        use super::super::threading::Segment;
        use std::io::Write;
        let mut state = 171823u64;
        let sequence: Vec<_> = (0..2200)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state % 4) as usize]
            })
            .collect();
        let mut decoy = sequence[900..970].to_vec();
        decoy[30] = if decoy[30] == b'A' { b'C' } else { b'A' };
        let graph = PathGraph {
            segments: vec![
                Segment {
                    sequence: sequence.clone(),
                    circular: false,
                },
                Segment {
                    sequence: decoy.clone(),
                    circular: false,
                },
            ],
            outgoing: vec![vec![], vec![]],
        };
        let out =
            std::env::temp_dir().join(format!("mitovar-calibration-decoy-{}", std::process::id()));
        std::fs::create_dir_all(&out).unwrap();
        let reads = Reads {
            r1: out.join("r1.fq"),
            r2: out.join("r2.fq"),
        };
        let mut r1 = std::fs::File::create(&reads.r1).unwrap();
        let mut r2 = std::fs::File::create(&reads.r2).unwrap();
        for start in 500..1000 {
            let a = &sequence[start..start + 100];
            let mut b = sequence[start + 200..start + 300].to_vec();
            if start == 700 {
                b[..70].copy_from_slice(&decoy);
            }
            writeln!(
                r1,
                "@r{start}/1\n{}\n+\n{}",
                String::from_utf8_lossy(a),
                "I".repeat(100)
            )
            .unwrap();
            writeln!(
                r2,
                "@r{start}/2\n{}\n+\n{}",
                String::from_utf8_lossy(&kmer::reverse_complement(&b)),
                "I".repeat(100)
            )
            .unwrap();
        }
        drop((r1, r2));
        let settings = GraphSettings {
            coverage_layout: false,
            multi_k_layout: false,
            k: 63,
            min_kmer_count: 3,
            max_error_fraction: 0.2,
            max_bubble_fraction: 0.49,
            min_depth: 10,
            min_junction_fragments: 3,
            max_validation_pairs: 10000,
        };
        let result = estimate(
            &graph,
            &reads,
            &settings,
            &Limits {
                max_kmers: 10000,
                min_baseq: 20,
            },
        )
        .unwrap();
        assert_eq!(
            result.report["training_fragments"], 500,
            "{}",
            result.report
        );
        assert_eq!(result.report["unmeasured_fragments"], 0);
        assert_eq!(result.report["median_span"], 300);
        std::fs::remove_dir_all(&out).unwrap();
    }
}
