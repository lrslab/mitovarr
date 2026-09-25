use std::cmp::Ordering;

use serde::Serialize;

use crate::pileup::{
    base_to_index, dominant_base, index_to_base, normalize_base, IndelStats, PileupPosition,
};

#[derive(Debug, Clone, Copy)]
pub struct ConsensusSettings {
    pub min_depth: u32,
    pub min_snp_frac: f64,
    pub min_indel_frac: f64,
    pub repeat_min_dominance: f64,
    pub repeat_max_entropy: f64,
}

const MIN_REPORTED_MINOR_FRACTION: f64 = 0.05;

#[derive(Debug, Clone, Serialize)]
pub struct AlleleVote {
    pub base: char,
    pub depth: u32,
    pub fraction: f64,
    pub weighted_fraction: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct BaseVoteRecord {
    pub pos: usize,
    pub reference: char,
    pub consensus: char,
    pub depth: u32,
    pub alleles: Vec<AlleleVote>,
    pub is_repeat_ambiguous: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct VariantRecord {
    pub pos: usize,
    pub r#ref: String,
    pub alt: String,
    pub depth: u32,
    pub alt_depth: u32,
    pub alt_frac: f64,
    pub baseq_mean: f64,
    pub mapq_mean: f64,
    pub flags: String,
    pub region_flag: String,
    pub is_repeat_ambiguous: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct CoverageRecord {
    pub pos: usize,
    pub depth: u32,
    pub mapq_median: f64,
    pub entropy: f64,
    pub is_repeat_ambiguous: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct UncertaintyRecord {
    pub pos: usize,
    pub r#ref: String,
    pub major_base: String,
    pub depth: u32,
    pub total_obs: u32,
    pub primary_obs: u32,
    pub sec_supp_obs: u32,
    pub major_frac: f64,
    pub mapq_median: f64,
    pub entropy: f64,
    pub is_repeat_ambiguous: bool,
    pub reasons: String,
}

#[derive(Debug, Clone)]
pub enum IndelCallKind {
    Insertion(String),
    Deletion(usize),
}

#[derive(Debug, Clone)]
pub struct ChosenIndel {
    pub kind: IndelCallKind,
    pub stats: IndelStats,
    pub frac: f64,
}

#[derive(Debug, Clone)]
pub struct ConsensusDraft {
    pub base_calls: Vec<u8>,
    pub indels: Vec<Option<ChosenIndel>>,
    pub variants: Vec<VariantRecord>,
    pub coverage: Vec<CoverageRecord>,
    pub uncertainty: Vec<UncertaintyRecord>,
    pub base_votes: Vec<BaseVoteRecord>,
}

#[derive(Debug, Clone)]
pub struct ConsensusResult {
    pub sequence: Vec<u8>,
    /// Zero-based output boundary before each reference base, plus the final boundary.
    pub reference_offsets: Vec<usize>,
    pub variants: Vec<VariantRecord>,
    pub coverage: Vec<CoverageRecord>,
    pub uncertainty: Vec<UncertaintyRecord>,
    pub base_votes: Vec<BaseVoteRecord>,
}

pub fn build_consensus_draft(
    reference: &[u8],
    pileup: &[PileupPosition],
    repeat_mask: &[bool],
    settings: ConsensusSettings,
) -> ConsensusDraft {
    let mut base_calls = Vec::with_capacity(reference.len());
    let mut coverage = Vec::with_capacity(reference.len());
    let mut uncertainty = Vec::new();
    let mut indels = vec![None; reference.len()];

    for (i, pos) in pileup.iter().enumerate() {
        let repeat = repeat_mask.get(i).copied().unwrap_or(false);

        let total_weight = pos.total_base_weight();
        let (major_idx, major_weight) = dominant_base(&pos.base_weights);
        let major_base = index_to_base(major_idx);
        let major_frac = if total_weight > 0.0 {
            major_weight / total_weight
        } else {
            0.0
        };

        let entropy = pos.entropy();

        let threshold = if repeat {
            settings.repeat_min_dominance.max(settings.min_snp_frac)
        } else {
            settings.min_snp_frac
        };

        let call = if repeat {
            if pos.depth < settings.min_depth || total_weight == 0.0 {
                b'N'
            } else if major_frac >= threshold && entropy <= settings.repeat_max_entropy {
                major_base
            } else {
                b'N'
            }
        } else if pos.depth < settings.min_depth || total_weight == 0.0 {
            b'N'
        } else if major_frac > 0.5 && major_frac >= threshold {
            major_base
        } else {
            b'N'
        };

        base_calls.push(call);

        coverage.push(CoverageRecord {
            pos: i + 1,
            depth: pos.depth,
            mapq_median: pos.mapq_median(),
            entropy,
            is_repeat_ambiguous: repeat,
        });

        if call == b'N' {
            let mut reasons = Vec::new();
            if pos.depth < settings.min_depth {
                reasons.push("low_depth");
            }
            if total_weight == 0.0 {
                reasons.push("no_weighted_observations");
            }
            if total_weight > 0.0 && (major_frac < threshold || (!repeat && major_frac <= 0.5)) {
                reasons.push("low_dominance");
            }
            if repeat && entropy > settings.repeat_max_entropy {
                reasons.push("high_entropy");
            }
            if major_base == b'N' {
                reasons.push("ambiguous_major_base");
            }
            if reasons.is_empty() {
                reasons.push("ambiguous_call");
            }

            uncertainty.push(UncertaintyRecord {
                pos: i + 1,
                r#ref: char::from(normalize_base(reference[i])).to_string(),
                major_base: char::from(major_base).to_string(),
                depth: pos.depth,
                total_obs: pos.total_obs,
                primary_obs: pos.primary_obs,
                sec_supp_obs: pos.sec_supp_obs,
                major_frac,
                mapq_median: pos.mapq_median(),
                entropy,
                is_repeat_ambiguous: repeat,
                reasons: reasons.join("|"),
            });
        }

        if let Some(chosen_indel) = choose_indel(pos, repeat, settings) {
            indels[i] = Some(chosen_indel);
        }
    }

    let variants = collect_variants(
        reference,
        pileup,
        repeat_mask,
        &base_calls,
        &indels,
        settings,
    );
    let base_votes = collect_base_votes(reference, pileup, repeat_mask, &base_calls, settings);

    ConsensusDraft {
        base_calls,
        indels,
        variants,
        coverage,
        uncertainty,
        base_votes,
    }
}

pub fn finalize_consensus(draft: ConsensusDraft) -> ConsensusResult {
    let mut sequence = Vec::with_capacity(draft.base_calls.len());
    let mut reference_offsets = vec![0; draft.base_calls.len() + 1];
    let mut retained_reference_positions = vec![false; draft.base_calls.len()];
    let mut skip_until = 0usize;

    for (i, retained) in retained_reference_positions.iter_mut().enumerate() {
        reference_offsets[i] = sequence.len();
        if i < skip_until {
            continue;
        }

        sequence.push(draft.base_calls[i]);
        *retained = true;

        if let Some(call) = draft.indels[i].as_ref() {
            match &call.kind {
                IndelCallKind::Insertion(seq) => {
                    sequence.extend(seq.as_bytes());
                }
                IndelCallKind::Deletion(len) => {
                    let end = i.saturating_add(*len + 1);
                    if end > skip_until {
                        skip_until = end;
                    }
                }
            }
        }
    }
    reference_offsets[draft.base_calls.len()] = sequence.len();

    let uncertainty = draft
        .uncertainty
        .into_iter()
        .filter(|record| {
            retained_reference_positions
                .get(record.pos.saturating_sub(1))
                .copied()
                .unwrap_or(false)
        })
        .collect();
    let base_votes = draft
        .base_votes
        .into_iter()
        .filter(|record| retained_reference_positions[record.pos - 1])
        .collect();

    ConsensusResult {
        sequence,
        reference_offsets,
        variants: draft.variants,
        coverage: draft.coverage,
        uncertainty,
        base_votes,
    }
}

fn collect_base_votes(
    reference: &[u8],
    pileup: &[PileupPosition],
    repeat_mask: &[bool],
    base_calls: &[u8],
    settings: ConsensusSettings,
) -> Vec<BaseVoteRecord> {
    reference
        .iter()
        .zip(pileup)
        .enumerate()
        .filter_map(|(i, (&reference, pos))| {
            let weight = pos.total_base_weight();
            if pos.depth < settings.min_depth || weight == 0.0 {
                return None;
            }
            // Report mixed A/C/G/T support independently of which allele matches
            // the reference or is written to the consensus; N remains in denominators.
            let supported = pos.base_weights[..4]
                .iter()
                .filter(|&&w| w / weight >= MIN_REPORTED_MINOR_FRACTION)
                .count();
            if supported < 2 {
                return None;
            }
            let alleles = (0..5)
                .filter(|&j| pos.base_counts[j] > 0)
                .map(|j| AlleleVote {
                    base: char::from(index_to_base(j)),
                    depth: pos.base_counts[j],
                    fraction: f64::from(pos.base_counts[j]) / f64::from(pos.depth),
                    weighted_fraction: pos.base_weights[j] / weight,
                })
                .collect();
            Some(BaseVoteRecord {
                pos: i + 1,
                reference: char::from(normalize_base(reference)),
                consensus: char::from(base_calls[i]),
                depth: pos.depth,
                alleles,
                is_repeat_ambiguous: repeat_mask.get(i).copied().unwrap_or(false),
            })
        })
        .collect()
}

fn collect_variants(
    reference: &[u8],
    pileup: &[PileupPosition],
    repeat_mask: &[bool],
    base_calls: &[u8],
    indels: &[Option<ChosenIndel>],
    settings: ConsensusSettings,
) -> Vec<VariantRecord> {
    let mut variants = Vec::new();

    for (i, (ref_base_raw, pos)) in reference.iter().copied().zip(pileup.iter()).enumerate() {
        let repeat = repeat_mask.get(i).copied().unwrap_or(false);
        let ref_base = normalize_base(ref_base_raw);
        let call = base_calls
            .get(i)
            .copied()
            .map(normalize_base)
            .unwrap_or(ref_base);
        let total_weight = pos.total_base_weight();
        let (major_idx, major_weight) = dominant_base(&pos.base_weights);
        let major_frac = if total_weight > 0.0 {
            major_weight / total_weight
        } else {
            0.0
        };

        let (minor_idx, minor_weight) = second_best_base(&pos.base_weights, major_idx);
        let minor_frac = if total_weight > 0.0 {
            minor_weight / total_weight
        } else {
            0.0
        };
        let minor_base = index_to_base(minor_idx);

        if pos.depth >= settings.min_depth {
            if call != ref_base {
                let alt_idx = if call == b'N' {
                    major_idx
                } else {
                    base_to_index(call)
                };
                let alt_frac = if call == b'N' {
                    major_frac
                } else if total_weight > 0.0 {
                    pos.base_weights[alt_idx] / total_weight
                } else {
                    0.0
                };

                variants.push(make_substitution_variant(
                    i,
                    ref_base,
                    call,
                    pos,
                    alt_idx,
                    alt_frac,
                    repeat,
                    "consensus_change",
                ));
            } else if minor_frac >= MIN_REPORTED_MINOR_FRACTION && minor_base != ref_base {
                variants.push(make_substitution_variant(
                    i,
                    ref_base,
                    minor_base,
                    pos,
                    minor_idx,
                    minor_frac,
                    repeat,
                    "heteroplasmy",
                ));
            }
        }

        if let Some(chosen_indel) = indels.get(i).and_then(Option::as_ref) {
            add_indel_variant(i, reference, pos.depth, repeat, chosen_indel, &mut variants);
        }
    }

    variants
}

fn choose_indel(
    pileup_pos: &PileupPosition,
    repeat: bool,
    settings: ConsensusSettings,
) -> Option<ChosenIndel> {
    let mut candidates = Vec::new();

    for (seq, stats) in &pileup_pos.insertions {
        candidates.push((
            stats.weight,
            IndelCallKind::Insertion(seq.clone()),
            stats.clone(),
            format!("ins:{}", seq),
        ));
    }

    for (len, stats) in &pileup_pos.deletions {
        candidates.push((
            stats.weight,
            IndelCallKind::Deletion(*len),
            stats.clone(),
            format!("del:{}", len),
        ));
    }

    if candidates.is_empty() {
        return None;
    }

    let entropy = pileup_pos.entropy();

    let threshold = if repeat {
        settings.repeat_min_dominance.max(settings.min_indel_frac)
    } else {
        settings.min_indel_frac
    };

    if pileup_pos.depth < settings.min_depth {
        return None;
    }

    let denom = f64::from(pileup_pos.depth).max(1.0);

    candidates.sort_by(|a, b| {
        let fa = f64::from(a.2.count) / denom;
        let fb = f64::from(b.2.count) / denom;

        fb.partial_cmp(&fa)
            .unwrap_or(Ordering::Equal)
            .then_with(|| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal))
            .then_with(|| a.3.cmp(&b.3))
    });

    let (_, best_kind, best_stats, _) = candidates[0].clone();
    let frac = f64::from(best_stats.count) / denom;

    if frac < threshold {
        return None;
    }

    if repeat && entropy > settings.repeat_max_entropy {
        return None;
    }

    Some(ChosenIndel {
        kind: best_kind,
        stats: best_stats,
        frac,
    })
}

fn second_best_base(base_weights: &[f64; 5], major_idx: usize) -> (usize, f64) {
    base_weights
        .iter()
        .copied()
        .enumerate()
        .filter(|(idx, _)| *idx != major_idx)
        .max_by(|(ia, wa), (ib, wb)| {
            if wa == wb {
                ia.cmp(ib)
            } else {
                wa.partial_cmp(wb).unwrap_or(Ordering::Equal)
            }
        })
        .unwrap_or((4, 0.0))
}

#[allow(clippy::too_many_arguments)]
fn make_substitution_variant(
    pos0: usize,
    ref_base: u8,
    alt_base: u8,
    pileup_pos: &PileupPosition,
    alt_idx: usize,
    alt_frac: f64,
    repeat: bool,
    call_flag: &str,
) -> VariantRecord {
    let alt_depth = pileup_pos.base_counts[alt_idx];
    let baseq_mean = if alt_depth > 0 {
        pileup_pos.baseq_sums[alt_idx] as f64 / f64::from(alt_depth)
    } else {
        0.0
    };
    let mapq_mean = if alt_depth > 0 {
        pileup_pos.mapq_sums[alt_idx] as f64 / f64::from(alt_depth)
    } else {
        0.0
    };

    let mut flags = vec![call_flag.to_string(), "snv".to_string()];
    if repeat {
        flags.push("repeat_ambiguous".to_string());
    }

    let (pos, ref_allele, alt_allele) = normalize_variant(
        pos0 + 1,
        String::from_utf8_lossy(&[ref_base]).to_string(),
        String::from_utf8_lossy(&[alt_base]).to_string(),
    );

    VariantRecord {
        pos,
        r#ref: ref_allele,
        alt: alt_allele,
        depth: pileup_pos.depth,
        alt_depth,
        alt_frac,
        baseq_mean,
        mapq_mean,
        flags: flags.join(";"),
        region_flag: if repeat { "repeat" } else { "normal" }.to_string(),
        is_repeat_ambiguous: repeat,
    }
}

fn add_indel_variant(
    anchor: usize,
    reference: &[u8],
    depth: u32,
    repeat: bool,
    chosen_indel: &ChosenIndel,
    variants: &mut Vec<VariantRecord>,
) {
    let Some(&anchor_base) = reference.get(anchor) else {
        return;
    };

    let (pos, ref_allele, alt_allele, flag) = match &chosen_indel.kind {
        IndelCallKind::Insertion(seq) => {
            let pos = anchor + 1;
            let ref_allele = String::from_utf8_lossy(&[anchor_base]).to_string();
            let alt_allele = format!("{}{}", ref_allele, seq);
            (pos, ref_allele, alt_allele, "insertion")
        }
        IndelCallKind::Deletion(len) => {
            if *len == 0 {
                return;
            }

            let end = anchor
                .saturating_add(*len)
                .min(reference.len().saturating_sub(1));
            let ref_bytes = &reference[anchor..=end];
            let ref_allele = String::from_utf8_lossy(ref_bytes).to_string();
            let alt_allele = String::from_utf8_lossy(&[anchor_base]).to_string();
            (anchor + 1, ref_allele, alt_allele, "deletion")
        }
    };

    let (pos, ref_allele, alt_allele) = normalize_variant(pos, ref_allele, alt_allele);

    let alt_depth = chosen_indel.stats.count;
    let baseq_mean = if alt_depth > 0 {
        chosen_indel.stats.baseq_sum as f64 / f64::from(alt_depth)
    } else {
        0.0
    };
    let mapq_mean = if alt_depth > 0 {
        chosen_indel.stats.mapq_sum as f64 / f64::from(alt_depth)
    } else {
        0.0
    };

    let mut flags = vec![flag.to_string()];
    if repeat {
        flags.push("repeat_ambiguous".to_string());
    }

    variants.push(VariantRecord {
        pos,
        r#ref: ref_allele,
        alt: alt_allele,
        depth,
        alt_depth,
        alt_frac: chosen_indel.frac,
        baseq_mean,
        mapq_mean,
        flags: flags.join(";"),
        region_flag: if repeat { "repeat" } else { "normal" }.to_string(),
        is_repeat_ambiguous: repeat,
    });
}

pub fn normalize_variant(
    pos: usize,
    ref_allele: String,
    alt_allele: String,
) -> (usize, String, String) {
    if ref_allele == alt_allele {
        return (pos, ref_allele, alt_allele);
    }

    let mut p = pos;
    let mut r = ref_allele.into_bytes();
    let mut a = alt_allele.into_bytes();

    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }

    while r.len() > 1 && a.len() > 1 && r.first() == a.first() {
        r.remove(0);
        a.remove(0);
        p = p.saturating_add(1);
    }

    (
        p,
        String::from_utf8_lossy(&r).to_string(),
        String::from_utf8_lossy(&a).to_string(),
    )
}

#[cfg(test)]
mod tests {
    use crate::pileup::init_pileup;

    use super::*;

    fn majority_settings() -> ConsensusSettings {
        ConsensusSettings {
            min_depth: 10,
            min_snp_frac: 0.5,
            min_indel_frac: 0.7,
            repeat_min_dominance: 0.85,
            repeat_max_entropy: 0.8,
        }
    }

    #[test]
    fn majority_snp_preserves_minor_votes_without_lowering_indel_support() {
        let reference = b"CC";
        let mut pileup = init_pileup(reference);
        for pos in &mut pileup {
            for i in 0..100 {
                pos.add_base(if i < 54 { b'C' } else { b'A' }, 40, 60);
            }
        }
        for _ in 0..60 {
            pileup[0].add_insertion("G".into(), 40, 60, 2400.0);
            pileup[1].add_deletion(1, 60, 2400.0);
        }
        let draft = build_consensus_draft(reference, &pileup, &[false; 2], majority_settings());
        assert_eq!(draft.base_calls, b"CC");
        assert!(draft.indels.iter().all(Option::is_none));
        assert!(draft.uncertainty.is_empty());
        for record in &draft.base_votes {
            assert_eq!(record.consensus, 'C');
            assert_eq!(record.depth, 100);
            assert_eq!(record.alleles.len(), 2);
            for allele in &record.alleles {
                let expected = if allele.base == 'C' { 54 } else { 46 };
                assert_eq!(allele.depth, expected);
                assert!((allele.fraction - f64::from(expected) / 100.0).abs() < 1e-12);
                assert!((allele.weighted_fraction - allele.fraction).abs() < 1e-12);
            }
        }
        let stricter_snp = ConsensusSettings {
            min_snp_frac: 0.7,
            ..majority_settings()
        };
        let draft = build_consensus_draft(reference, &pileup, &[false; 2], stricter_snp);
        assert_eq!(draft.base_calls, b"NN");
        assert!(draft.indels.iter().all(Option::is_none));
        let lower_indel = ConsensusSettings {
            min_indel_frac: 0.5,
            ..majority_settings()
        };
        let draft = build_consensus_draft(reference, &pileup, &[false; 2], lower_indel);
        assert_eq!(draft.base_calls, b"CC");
        assert!(draft.indels.iter().all(Option::is_some));
    }

    #[test]
    fn majority_uses_quality_and_keeps_ties_repeat_and_low_depth_unresolved() {
        let reference = b"CCCCC";
        let mut pileup = init_pileup(reference);
        for i in 0..10 {
            pileup[0].add_base(if i < 6 { b'C' } else { b'A' }, 40, 60);
            pileup[1].add_base(if i < 5 { b'C' } else { b'A' }, 40, 60);
            pileup[2].add_base(if i < 6 { b'C' } else { b'A' }, 40, 60);
            pileup[4].add_base(
                if i < 6 { b'C' } else { b'A' },
                if i < 6 { 20 } else { 40 },
                60,
            );
        }
        for _ in 0..9 {
            pileup[3].add_base(b'C', 40, 60);
        }
        let mut draft = build_consensus_draft(
            reference,
            &pileup,
            &[false, false, true, false, false],
            majority_settings(),
        );
        assert_eq!(draft.base_calls, b"CNNNA");
        assert_eq!(draft.base_votes.len(), 4);
        assert_eq!(draft.base_votes[3].consensus, 'A');
        // Votes from accepted deleted reference positions must not be described
        // as mixed sites in the emitted sequence.
        draft.indels[0] = Some(ChosenIndel {
            kind: IndelCallKind::Deletion(2),
            stats: IndelStats::default(),
            frac: 1.0,
        });
        let result = finalize_consensus(draft);
        assert_eq!(result.sequence, b"CNA");
        assert_eq!(
            result.base_votes.iter().map(|r| r.pos).collect::<Vec<_>>(),
            vec![1, 5]
        );
    }

    #[test]
    fn test_normalize_variant_minimal_representation() {
        let (p, r, a) = normalize_variant(10, "AAT".to_string(), "ACT".to_string());
        assert_eq!(p, 11);
        assert_eq!(r, "A");
        assert_eq!(a, "C");
    }

    #[test]
    fn test_repeat_dominance_gate_emits_n_when_weak() {
        let reference = b"AAAA";
        let mut pileup = init_pileup(reference);
        pileup[0].depth = 20;
        pileup[0].base_weights = [10.0, 9.0, 0.0, 0.0, 0.0];

        let repeat_mask = vec![true, false, false, false];

        let draft = build_consensus_draft(
            reference,
            &pileup,
            &repeat_mask,
            ConsensusSettings {
                min_depth: 10,
                min_snp_frac: 0.7,
                min_indel_frac: 0.7,
                repeat_min_dominance: 0.85,
                repeat_max_entropy: 1.2,
            },
        );

        assert_eq!(draft.base_calls[0], b'N');
        assert!(draft.uncertainty[0].reasons.contains("low_dominance"));
        assert!(draft.uncertainty[0].is_repeat_ambiguous);
    }

    #[test]
    fn test_non_repeat_weak_support_emits_n() {
        let reference = b"AAAA";
        let mut pileup = init_pileup(reference);
        pileup[0].depth = 20;
        pileup[0].base_weights = [10.0, 9.0, 0.0, 0.0, 0.0];

        let draft = build_consensus_draft(
            reference,
            &pileup,
            &[false, false, false, false],
            ConsensusSettings {
                min_depth: 10,
                min_snp_frac: 0.7,
                min_indel_frac: 0.7,
                repeat_min_dominance: 0.85,
                repeat_max_entropy: 1.2,
            },
        );

        assert_eq!(draft.base_calls[0], b'N');
        assert!(draft.uncertainty[0].reasons.contains("low_dominance"));
        assert!(!draft.uncertainty[0].is_repeat_ambiguous);
    }

    #[test]
    fn test_repeat_low_depth_emits_n() {
        let reference = b"AAAA";
        let mut pileup = init_pileup(reference);
        pileup[0].depth = 2;
        pileup[0].base_weights = [2.0, 0.0, 0.0, 0.0, 0.0];

        let draft = build_consensus_draft(
            reference,
            &pileup,
            &[true, false, false, false],
            ConsensusSettings {
                min_depth: 10,
                min_snp_frac: 0.7,
                min_indel_frac: 0.7,
                repeat_min_dominance: 0.85,
                repeat_max_entropy: 1.2,
            },
        );

        assert_eq!(draft.base_calls[0], b'N');
        assert!(draft.uncertainty[0].reasons.contains("low_depth"));
        assert_eq!(draft.uncertainty[0].depth, 2);
    }

    #[test]
    fn test_finalize_consensus_applies_insertion() {
        let draft = ConsensusDraft {
            base_calls: b"AC".to_vec(),
            indels: vec![
                Some(ChosenIndel {
                    kind: IndelCallKind::Insertion("GG".to_string()),
                    stats: IndelStats::default(),
                    frac: 1.0,
                }),
                None,
            ],
            variants: Vec::new(),
            coverage: Vec::new(),
            uncertainty: Vec::new(),
            base_votes: Vec::new(),
        };

        let result = finalize_consensus(draft);
        assert_eq!(result.sequence, b"AGGC");
    }

    #[test]
    fn test_finalize_consensus_applies_deletion() {
        let draft = ConsensusDraft {
            base_calls: b"ACGT".to_vec(),
            indels: vec![
                Some(ChosenIndel {
                    kind: IndelCallKind::Deletion(2),
                    stats: IndelStats::default(),
                    frac: 1.0,
                }),
                None,
                None,
                None,
            ],
            variants: Vec::new(),
            coverage: Vec::new(),
            uncertainty: Vec::new(),
            base_votes: Vec::new(),
        };

        let result = finalize_consensus(draft);
        assert_eq!(result.sequence, b"AT");
    }

    #[test]
    fn test_finalize_consensus_excludes_uncertainty_at_deleted_reference_positions() {
        let reference = b"ACGT";
        let pileup = init_pileup(reference);
        let mut draft = build_consensus_draft(
            reference,
            &pileup,
            &[false; 4],
            ConsensusSettings {
                min_depth: 10,
                min_snp_frac: 0.7,
                min_indel_frac: 0.7,
                repeat_min_dominance: 0.85,
                repeat_max_entropy: 1.2,
            },
        );
        draft.indels[0] = Some(ChosenIndel {
            kind: IndelCallKind::Deletion(2),
            stats: IndelStats::default(),
            frac: 1.0,
        });

        let result = finalize_consensus(draft);
        let positions = result
            .uncertainty
            .iter()
            .map(|record| record.pos)
            .collect::<Vec<_>>();

        assert_eq!(positions, vec![1, 4]);
    }

    #[test]
    fn test_indel_call_requires_support_fraction() {
        let reference = b"AC";
        let mut pileup = init_pileup(reference);

        for _ in 0..100 {
            pileup[0].add_base(b'A', 30, 60);
        }

        pileup[0].add_deletion(1, 60, 1200.0);

        let draft = build_consensus_draft(
            reference,
            &pileup,
            &[false, false],
            ConsensusSettings {
                min_depth: 10,
                min_snp_frac: 0.7,
                min_indel_frac: 0.7,
                repeat_min_dominance: 0.85,
                repeat_max_entropy: 1.2,
            },
        );

        assert!(draft.indels[0].is_none(), "expected no deletion call at 1%");

        for _ in 0..80 {
            pileup[0].add_deletion(1, 60, 1200.0);
        }

        let draft = build_consensus_draft(
            reference,
            &pileup,
            &[false, false],
            ConsensusSettings {
                min_depth: 10,
                min_snp_frac: 0.7,
                min_indel_frac: 0.7,
                repeat_min_dominance: 0.85,
                repeat_max_entropy: 1.2,
            },
        );

        let call = draft.indels[0].as_ref().expect("expected deletion call");
        match call.kind {
            IndelCallKind::Deletion(len) => assert_eq!(len, 1),
            _ => panic!("expected deletion call"),
        }
        assert!(call.frac >= 0.7);
    }
}
