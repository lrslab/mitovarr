//! Shared whole-hypothesis fragment evidence. An alignment span is not a
//! sequenced junction, and alternative placements never become extra molecules.
use super::{
    align::{self, BaseTrace, Hit, Index, Placements},
    pairs::Library,
};
use crate::io::fastq;
use std::collections::BTreeMap;

mod indexed;

#[derive(Default)]
pub(crate) struct RescueStats {
    pub searches: usize,
    pub added_placements: usize,
    pub seeded_decoy_searches: usize,
    pub states: Vec<&'static str>,
}

/// Every hypothesis searches the bounded mate neighborhood of the same fixed
/// common-flank pivot. Sequence-only scores cannot bound the joint score.
/// Complete states describe endpoint/edit enumeration inside the declared
/// window. Resource-limited searches are omitted from every hypothesis.
pub(crate) fn rescue_pool(
    indexes: &[Index],
    placements: &mut [(Placements, Placements)],
    one: &mut align::Prepared,
    two: &mut align::Prepared,
    model: Option<&Library>,
    pivot_mate: usize,
) -> RescueStats {
    let mut stats = RescueStats::default();
    for (i, (a, b)) in placements.iter_mut().enumerate() {
        let (anchor, other, prepared) = if pivot_mate == 0 {
            (&*a, b, &mut *two)
        } else {
            (&*b, a, &mut *one)
        };
        let status = if anchor.exhausted {
            "alignment_limit"
        } else if !anchor.complete {
            "anchor_constraint_search_incomplete"
        } else if anchor.hits.is_empty() {
            "complete_no_anchor_alignment"
        } else if let Some(model) = model.filter(|m| m.geometry_usable) {
            // All retained anchor alignments observe the same query word at
            // the same physical target word; cover their small indel spread.
            let mut pivot = anchor.hits[0].clone();
            let circle = indexes[i].circle_length().unwrap_or(0) as isize;
            let bounds: Vec<_> = anchor
                .hits
                .iter()
                .map(|h| {
                    let mut delta = h.start as isize - pivot.start as isize;
                    if circle > 0 {
                        if delta > circle / 2 {
                            delta -= circle;
                        }
                        if delta < -circle / 2 {
                            delta += circle;
                        }
                    }
                    let start = pivot.start as isize + delta;
                    (start, start + (h.end - h.start) as isize)
                })
                .collect();
            let left = bounds.iter().map(|h| h.0).min().unwrap();
            let right = bounds.iter().map(|h| h.1).max().unwrap();
            let shift = if left < 0 { circle } else { 0 };
            pivot.start = (left + shift) as usize;
            pivot.end = (right + shift) as usize;
            stats.searches += 1;
            stats.seeded_decoy_searches += usize::from(!other.hits.is_empty());
            let before = other.hits.len();
            *other = indexes[i].profile_around(prepared, &pivot, model.upper, !pivot.reverse);
            stats.added_placements += other.hits.len().saturating_sub(before);
            if other.exhausted {
                "rescue_limit"
            } else if other.hits.is_empty() {
                "complete_no_mate_alignment"
            } else {
                "complete_mate_search"
            }
        } else {
            stats.searches += 1;
            *other = indexes[i].profile_all(prepared);
            if other.complete {
                "complete_sequence_search"
            } else {
                "sequence_search_limit"
            }
        };
        stats.states.push(status);
    }
    stats
}

pub(crate) fn unknown_reason(
    placements: &[(Placements, Placements)],
    states: &[&str],
    pivot_mate: usize,
) -> Option<&'static str> {
    if placements.iter().any(|(a, b)| a.exhausted || b.exhausted) {
        return Some("alignment_or_rescue_limit");
    }
    if placements.iter().zip(states).any(|((a, b), state)| {
        let (pivot, mate) = if pivot_mate == 0 { (a, b) } else { (b, a) };
        !pivot.complete
            || (!pivot.hits.is_empty() && !mate.complete)
            || *state == "anchor_no_full_query_alignment"
    }) {
        return Some("incomplete_constrained_alignment_domain");
    }
    None
}

/// Maximum change in a pairwise score contrast from u omitted fragments.
pub(crate) fn unknown_contrast_bound(unknown: usize) -> f64 {
    unknown as f64 * (12.0 + (-12.0f64).exp().ln_1p())
}

fn observed(trace: &BaseTrace, start: usize, end: usize, shift: usize) -> bool {
    start
        .checked_sub(shift)
        .zip(end.checked_sub(shift))
        .is_some_and(|(a, b)| trace.matches_contiguous(a, b))
}

fn crosses(hit: &Hit, trace: &BaseTrace, junction: usize, shift: usize) -> bool {
    hit.start + 20 <= junction
        && hit.end >= junction + 20
        && observed(trace, junction - 20, junction + 20, shift)
}

pub(crate) fn observations(
    placements: &Placements,
    read: &fastq::Record,
    index: &Index,
    minq: u8,
) -> Vec<(Hit, BaseTrace)> {
    placements
        .hits
        .iter()
        .map(|hit| {
            (
                hit.clone(),
                align::trace_alignment(read, hit, &index.sequence, minq),
            )
        })
        .collect()
}

pub(crate) struct JointEvidence {
    pub score: f64,
    pub junctions: [bool; 2],
    pub observed_junctions: [bool; 2],
    pub mate_brackets: [bool; 2],
    pub contraction: bool,
    pub diagnostics: serde_json::Value,
}

/// A repeated hit never increases a fragment's mass: retain only the best
/// sequence explanation within each library span state, as in de novo scoring.
#[allow(clippy::too_many_arguments)]
pub(crate) fn joint(
    a: &[(Hit, BaseTrace)],
    b: &[(Hit, BaseTrace)],
    model: Option<&Library>,
    ideal: f64,
    length: usize,
    junctions: [usize; 2],
    contraction: Option<[usize; 2]>,
) -> JointEvidence {
    if a.len().saturating_mul(b.len()) <= 64
        || !ideal.is_finite()
        || a.iter().chain(b).any(|(h, _)| !h.logp.is_finite())
    {
        joint_exhaustive(a, b, model, ideal, length, junctions, contraction)
    } else {
        indexed::joint(a, b, model, ideal, length, junctions, contraction)
    }
}

/// Small-input fallback and independent oracle for the indexed implementation.
#[allow(clippy::too_many_arguments)]
fn joint_exhaustive(
    a: &[(Hit, BaseTrace)],
    b: &[(Hit, BaseTrace)],
    model: Option<&Library>,
    ideal: f64,
    length: usize,
    junctions: [usize; 2],
    contraction: Option<[usize; 2]>,
) -> JointEvidence {
    let mut states = BTreeMap::<usize, (f64, usize)>::new();
    let mut best = 0.0f64;
    let mut bridge = [0.0f64; 2];
    let mut observed_bridge = [0.0f64; 2];
    let mut bracket_bridge = [0.0f64; 2];
    let mut contraction_bridge = 0.0f64;
    let empirical = model.filter(|m| m.geometry_usable);
    for (x, trace_x) in a {
        for (original_y, trace_y) in b {
            // Include a pair crossing the arbitrary FASTA origin. Moving one
            // end by a whole circle is geometry only, never an extra observation.
            for shift in [0, length].into_iter().take(if length > 0 { 2 } else { 1 }) {
                let (mut shift_x, mut shift_y) = (0, 0);
                if shift > 0 {
                    if x.start < original_y.start {
                        shift_x = shift;
                    } else {
                        shift_y = shift;
                    }
                }
                let cx = (x.start + shift_x, x.end + shift_x, x.reverse);
                let cy = (
                    original_y.start + shift_y,
                    original_y.end + shift_y,
                    original_y.reverse,
                );
                if length > 0 && cx.1.max(cy.1) - cx.0.min(cy.0) > length {
                    continue;
                }
                let span = match model {
                    Some(m) => match m.span_coordinates(cx, cy) {
                        Some(d) => d,
                        None => continue,
                    },
                    None => 0,
                };
                let weight = empirical.map_or(1.0, |m| m.probability(span));
                if weight == 0.0 {
                    continue;
                }
                let log_sequence = (x.logp + original_y.logp - ideal).min(0.0);
                states
                    .entry(if empirical.is_some() { span } else { 0 })
                    .and_modify(|v| {
                        if log_sequence > v.0 {
                            *v = (
                                log_sequence,
                                trace_x.ambiguous_contacts + trace_y.ambiguous_contacts,
                            );
                        } else if log_sequence == v.0 {
                            v.1 =
                                v.1.max(trace_x.ambiguous_contacts + trace_y.ambiguous_contacts);
                        }
                    })
                    .or_insert((
                        log_sequence,
                        trace_x.ambiguous_contacts + trace_y.ambiguous_contacts,
                    ));
                // Geometry weight / peak <= 1. Such a placement can add span
                // mass, but cannot reach the junction support threshold.
                if log_sequence < 0.1f64.ln() {
                    continue;
                }
                let relative =
                    log_sequence.exp() * empirical.map_or(1.0, |m| weight / m.peak_probability());
                best = best.max(relative);
                if relative < 0.1 {
                    continue;
                }
                let x = Hit {
                    start: cx.0,
                    end: cx.1,
                    ..x.clone()
                };
                let y = Hit {
                    start: cy.0,
                    end: cy.1,
                    ..original_y.clone()
                };
                for (j, &junction) in junctions.iter().enumerate() {
                    let mut witnessed = false;
                    for p in [junction, junction + length] {
                        let direct =
                            crosses(&x, trace_x, p, shift_x) || crosses(&y, trace_y, p, shift_y);
                        // A mate bracket is geometry, not a sequenced junction.
                        // Both ends must actually observe their boundary anchors
                        // and lie on opposite sides. A deletion across p cannot
                        // masquerade as a bracket using its alignment span.
                        let bracket = |left: &Hit,
                                       lt: &BaseTrace,
                                       ls: usize,
                                       right: &Hit,
                                       rt: &BaseTrace,
                                       rs: usize| {
                            left.end <= p
                                && right.start >= p
                                && left.end >= left.start + 20
                                && right.end >= right.start + 20
                                && observed(lt, left.end - 20, left.end, ls)
                                && observed(rt, right.start, right.start + 20, rs)
                        };
                        let paired = empirical.is_some_and(|m| span <= m.upper)
                            && (bracket(&x, trace_x, shift_x, &y, trace_y, shift_y)
                                || bracket(&y, trace_y, shift_y, &x, trace_x, shift_x));
                        if direct {
                            observed_bridge[j] = observed_bridge[j].max(relative);
                        }
                        if paired {
                            bracket_bridge[j] = bracket_bridge[j].max(relative);
                        }
                        witnessed |= direct || paired;
                    }
                    if witnessed {
                        bridge[j] = bridge[j].max(relative);
                    }
                }
                if let Some(anchors) = contraction {
                    let linked = [0, length].into_iter().any(|offset| {
                        let contains = |trace: &BaseTrace, shift, p| {
                            observed(trace, p + offset, p + offset + 21, shift)
                        };
                        let one_read = [(trace_x, shift_x), (trace_y, shift_y)]
                            .iter()
                            .any(|&(t, s)| anchors.iter().all(|&p| contains(t, s, p)));
                        let paired = empirical.is_some_and(|m| span <= m.upper)
                            && anchors.iter().all(|&p| {
                                contains(trace_x, shift_x, p) || contains(trace_y, shift_y, p)
                            });
                        one_read || paired
                    });
                    if linked {
                        contraction_bridge = contraction_bridge.max(relative);
                    }
                }
            }
        }
    }
    let span_states:Vec<_>=states.into_iter().map(|(span,(log_sequence,ambiguous))| {
        let sequence = log_sequence.exp();
        let geometry=model.filter(|m|m.geometry_usable).map_or(1.0,|m|m.probability(span));
        serde_json::json!({"span":span,"sequence_weight":sequence,"geometry_weight":geometry,"contribution":sequence*geometry,"ambiguous_base_contacts":ambiguous})
    }).collect();
    let mass: f64 = span_states
        .iter()
        .map(|s| s["contribution"].as_f64().unwrap())
        .sum();
    let floor = (-12.0f64).exp();
    let supported = |p| p >= 0.1 && p >= 0.99 * best;
    JointEvidence {
        score: (mass.min(1.0) + (-12.0f64).exp()).ln(),
        junctions: bridge.map(supported),
        observed_junctions: observed_bridge.map(supported),
        mate_brackets: bracket_bridge.map(supported),
        contraction: supported(contraction_bridge),
        diagnostics: serde_json::json!({"span_states":span_states,"mass_before_floor":mass.min(1.0),"outlier_mass_fraction":floor/(mass.min(1.0)+floor),"score":(mass.min(1.0)+floor).ln()}),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn omitted_fragments_have_a_finite_exact_contrast_bound() {
        let worst = (-12.0f64).exp().ln();
        let best = (1.0 + (-12.0f64).exp()).ln();
        assert!((unknown_contrast_bound(3) - 3.0 * (best - worst)).abs() < 1e-12);
        assert_eq!(unknown_contrast_bound(0), 0.0);
        assert!(unknown_contrast_bound(1) < 12.00001);
    }

    #[test]
    fn missing_mate_without_rescue_is_unknown_and_duplicate_hits_do_not_add_mass() {
        let hit = Hit {
            start: 10,
            end: 70,
            reverse: false,
            logp: 0.0,
            edits: 0,
            traceback: None,
        };
        let mut pool = vec![(
            Placements {
                hits: vec![hit.clone()],
                exhausted: false,
                complete: false,
            },
            Placements::default(),
        )];
        assert_eq!(
            unknown_reason(&pool, &["geometry_unavailable_seeded_only"], 0),
            Some("incomplete_constrained_alignment_domain")
        );
        pool[0].0.complete = true;
        pool[0].1.complete = true;
        assert_eq!(
            unknown_reason(&pool, &["complete_no_mate_alignment"], 0),
            None
        );
        pool[0].1.exhausted = true;
        assert_eq!(
            unknown_reason(&pool, &["rescue_limit"], 0),
            Some("alignment_or_rescue_limit")
        );
        let a = vec![(hit.clone(), BaseTrace::default())];
        let b = vec![(
            Hit {
                start: 90,
                end: 150,
                reverse: true,
                ..hit.clone()
            },
            BaseTrace::default(),
        )];
        let single = joint(&a, &b, None, 0.0, 0, [80, 80], None);
        let duplicate = vec![
            (hit.clone(), BaseTrace::default()),
            (hit, BaseTrace::default()),
        ];
        let repeated = joint(&duplicate, &b, None, 0.0, 0, [80, 80], None);
        assert_eq!(single.score, repeated.score);
        assert_eq!(single.diagnostics, repeated.diagnostics);
    }
}
