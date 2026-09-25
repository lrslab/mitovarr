use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet},
};

use noodles_sam::alignment::{
    record::cigar::{op::Kind, Op},
    RecordBuf,
};

#[derive(Debug, Clone, Copy)]
pub struct PileupSettings {
    pub min_mapq: u8,
    pub min_baseq: u8,
}

#[derive(Debug, Clone, Default)]
pub struct PileupRunStats {
    pub total_records: u64,
    pub filtered_records: u64,
    pub primary_records: u64,
    pub secondary_or_supp_records: u64,
    /// Additional base/indel observations removed because two records from the
    /// same physical fragment covered the same assembly coordinate.
    pub overlapping_fragment_observations_collapsed: u64,
    /// Overlapping fragment observations that disagreed before the explicit
    /// base-quality, then mapping-quality conflict rule was applied.
    pub conflicting_fragment_observations: u64,
}

/// Controls how alignment coordinates are projected onto the pileup target.
///
/// Circular mode folds coordinates modulo the target length and can recover a
/// conservatively matched terminal soft clip across the literal origin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PileupCoordinateMode {
    Linear,
    Circular,
}

#[derive(Debug, Clone, Default)]
pub struct IndelStats {
    pub count: u32,
    pub weight: f64,
    pub baseq_sum: u64,
    pub mapq_sum: u64,
}

#[derive(Debug, Clone)]
pub struct PileupPosition {
    pub ref_base: u8,
    pub depth: u32,
    pub primary_obs: u32,
    pub total_obs: u32,
    pub sec_supp_obs: u32,
    pub mapq_hist: [u32; 61],
    pub base_counts: [u32; 5],
    pub base_weights: [f64; 5],
    pub baseq_sums: [u64; 5],
    pub mapq_sums: [u64; 5],
    pub insertions: BTreeMap<String, IndelStats>,
    pub deletions: BTreeMap<usize, IndelStats>,
}

impl PileupPosition {
    pub fn new(ref_base: u8) -> Self {
        Self {
            ref_base,
            depth: 0,
            primary_obs: 0,
            total_obs: 0,
            sec_supp_obs: 0,
            mapq_hist: [0; 61],
            base_counts: [0; 5],
            base_weights: [0.0; 5],
            baseq_sums: [0; 5],
            mapq_sums: [0; 5],
            insertions: BTreeMap::new(),
            deletions: BTreeMap::new(),
        }
    }

    pub fn add_observation(&mut self, mapq: u8, sec_supp: bool) {
        self.total_obs = self.total_obs.saturating_add(1);
        let mapq_idx = usize::from(mapq.min(60));
        self.mapq_hist[mapq_idx] = self.mapq_hist[mapq_idx].saturating_add(1);
        if sec_supp {
            self.sec_supp_obs = self.sec_supp_obs.saturating_add(1);
        } else {
            self.primary_obs = self.primary_obs.saturating_add(1);
        }
    }

    pub fn add_base(&mut self, base: u8, baseq: u8, mapq: u8) {
        let idx = base_to_index(base);
        self.depth = self.depth.saturating_add(1);
        self.base_counts[idx] = self.base_counts[idx].saturating_add(1);

        let weight = quality_weight(baseq, mapq);
        self.base_weights[idx] += weight;
        self.baseq_sums[idx] = self.baseq_sums[idx].saturating_add(u64::from(baseq));
        self.mapq_sums[idx] = self.mapq_sums[idx].saturating_add(u64::from(mapq));
    }

    pub fn add_insertion(&mut self, sequence: String, baseq_mean: u8, mapq: u8, weight: f64) {
        let stats = self.insertions.entry(sequence).or_default();
        stats.count = stats.count.saturating_add(1);
        stats.weight += weight;
        stats.baseq_sum = stats.baseq_sum.saturating_add(u64::from(baseq_mean));
        stats.mapq_sum = stats.mapq_sum.saturating_add(u64::from(mapq));
    }

    pub fn add_deletion(&mut self, len: usize, mapq: u8, weight: f64) {
        let stats = self.deletions.entry(len).or_default();
        stats.count = stats.count.saturating_add(1);
        stats.weight += weight;
        stats.baseq_sum = stats.baseq_sum.saturating_add(0);
        stats.mapq_sum = stats.mapq_sum.saturating_add(u64::from(mapq));
    }

    pub fn total_base_weight(&self) -> f64 {
        self.base_weights.iter().sum()
    }

    pub fn entropy(&self) -> f64 {
        shannon_entropy(&self.base_weights)
    }

    pub fn mapq_median(&self) -> f64 {
        histogram_median(&self.mapq_hist)
    }
}

pub fn init_pileup(reference: &[u8]) -> Vec<PileupPosition> {
    reference
        .iter()
        .map(|&base| PileupPosition::new(normalize_base(base)))
        .collect()
}

const MIN_CIRCULAR_SOFTCLIP_ANCHOR: usize = 8;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct QualityObservation {
    baseq: u8,
    mapq: u8,
}

impl QualityObservation {
    fn score(self) -> (u8, u8) {
        (self.baseq, self.mapq)
    }
}

#[derive(Debug, Clone, Default)]
struct BaseEvidence {
    by_base: [Option<QualityObservation>; 5],
}

impl BaseEvidence {
    fn add(&mut self, base: u8, baseq: u8, mapq: u8) {
        let idx = base_to_index(base);
        let candidate = QualityObservation { baseq, mapq };
        let replace = self.by_base[idx].is_none_or(|current| candidate.score() > current.score());
        if replace {
            self.by_base[idx] = Some(candidate);
        }
    }

    fn observed_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.by_base
            .iter()
            .enumerate()
            .filter_map(|(idx, observation)| observation.map(|_| idx))
    }

    fn disagrees_with(&self, other: &Self) -> bool {
        let lhs = self.observed_indices().collect::<BTreeSet<_>>();
        let rhs = other.observed_indices().collect::<BTreeSet<_>>();
        !lhs.is_empty() && !rhs.is_empty() && lhs.is_disjoint(&rhs)
    }

    fn merge(&mut self, other: Self) {
        for (idx, observation) in other.by_base.into_iter().enumerate() {
            let Some(observation) = observation else {
                continue;
            };
            let base = index_to_base(idx);
            self.add(base, observation.baseq, observation.mapq);
        }
    }

    fn resolve(&self) -> Option<(u8, u8, u8)> {
        let best_score = self
            .by_base
            .iter()
            .flatten()
            .map(|observation| observation.score())
            .max()?;
        let mut winners = self
            .by_base
            .iter()
            .enumerate()
            .filter_map(|(idx, observation)| {
                observation
                    .filter(|observation| observation.score() == best_score)
                    .map(|observation| (idx, observation))
            });
        let (winner_idx, winner) = winners.next()?;
        if winners.next().is_some() {
            // Equal baseQ and MAPQ but different alleles: retain one fragment
            // observation, explicitly marked ambiguous rather than double count.
            Some((b'N', winner.baseq, winner.mapq))
        } else {
            Some((index_to_base(winner_idx), winner.baseq, winner.mapq))
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct CoverageObservation {
    mapq: u8,
    sec_supp: bool,
}

#[derive(Debug, Clone, Default)]
struct FragmentContribution {
    observations: BTreeMap<usize, CoverageObservation>,
    bases: BTreeMap<usize, BaseEvidence>,
    insertions: BTreeMap<usize, BTreeMap<String, QualityObservation>>,
    deletions: BTreeMap<usize, BTreeMap<usize, QualityObservation>>,
}

impl FragmentContribution {
    fn observe(&mut self, pos: usize, mapq: u8, sec_supp: bool) {
        self.observations
            .entry(pos)
            .and_modify(|observation| {
                observation.mapq = observation.mapq.max(mapq);
                // A fragment is primary at a position when at least one of its
                // contributing alignments is primary there.
                observation.sec_supp &= sec_supp;
            })
            .or_insert(CoverageObservation { mapq, sec_supp });
    }

    fn add_base(&mut self, pos: usize, base: u8, baseq: u8, mapq: u8) {
        self.bases.entry(pos).or_default().add(base, baseq, mapq);
    }

    fn add_insertion(&mut self, anchor: usize, sequence: String, baseq: u8, mapq: u8) {
        add_indel_evidence(
            self.insertions.entry(anchor).or_default(),
            sequence,
            QualityObservation { baseq, mapq },
        );
    }

    fn add_deletion(&mut self, anchor: usize, len: usize, baseq: u8, mapq: u8) {
        add_indel_evidence(
            self.deletions.entry(anchor).or_default(),
            len,
            QualityObservation { baseq, mapq },
        );
    }

    fn merge(&mut self, other: Self, stats: &mut PileupRunStats) {
        for (pos, observation) in other.observations {
            if self.observations.contains_key(&pos) {
                stats.overlapping_fragment_observations_collapsed = stats
                    .overlapping_fragment_observations_collapsed
                    .saturating_add(1);
            }
            self.observe(pos, observation.mapq, observation.sec_supp);
        }

        for (pos, evidence) in other.bases {
            let target = self.bases.entry(pos).or_default();
            if target.disagrees_with(&evidence) {
                stats.conflicting_fragment_observations =
                    stats.conflicting_fragment_observations.saturating_add(1);
            }
            target.merge(evidence);
        }

        merge_indel_maps(
            &mut self.insertions,
            other.insertions,
            &mut stats.conflicting_fragment_observations,
        );
        merge_indel_maps(
            &mut self.deletions,
            other.deletions,
            &mut stats.conflicting_fragment_observations,
        );
    }

    fn apply(self, pileup: &mut [PileupPosition]) {
        for (pos, observation) in self.observations {
            if let Some(pileup_pos) = pileup.get_mut(pos) {
                pileup_pos.add_observation(observation.mapq, observation.sec_supp);
            }
        }

        for (pos, evidence) in self.bases {
            let Some((base, baseq, mapq)) = evidence.resolve() else {
                continue;
            };
            if let Some(pileup_pos) = pileup.get_mut(pos) {
                pileup_pos.add_base(base, baseq, mapq);
            }
        }

        for (anchor, evidence) in self.insertions {
            let Some((sequence, observation)) = resolve_indel_evidence(evidence) else {
                continue;
            };
            let weight = quality_weight(observation.baseq, observation.mapq);
            if let Some(pileup_pos) = pileup.get_mut(anchor) {
                pileup_pos.add_insertion(sequence, observation.baseq, observation.mapq, weight);
            }
        }

        for (anchor, evidence) in self.deletions {
            let Some((len, observation)) = resolve_indel_evidence(evidence) else {
                continue;
            };
            let weight = quality_weight(observation.baseq, observation.mapq);
            if let Some(pileup_pos) = pileup.get_mut(anchor) {
                pileup_pos.add_deletion(len, observation.mapq, weight);
            }
        }
    }
}

fn add_indel_evidence<K: Ord>(
    evidence: &mut BTreeMap<K, QualityObservation>,
    allele: K,
    candidate: QualityObservation,
) {
    let replace = evidence
        .get(&allele)
        .is_none_or(|current| candidate.score() > current.score());
    if replace {
        evidence.insert(allele, candidate);
    }
}

fn merge_indel_maps<K: Ord>(
    target: &mut BTreeMap<usize, BTreeMap<K, QualityObservation>>,
    source: BTreeMap<usize, BTreeMap<K, QualityObservation>>,
    conflicts: &mut u64,
) {
    for (anchor, source_evidence) in source {
        let target_evidence = target.entry(anchor).or_default();
        if !target_evidence.is_empty()
            && source_evidence
                .keys()
                .any(|allele| !target_evidence.contains_key(allele))
        {
            *conflicts = conflicts.saturating_add(1);
        }
        for (allele, observation) in source_evidence {
            add_indel_evidence(target_evidence, allele, observation);
        }
    }
}

fn resolve_indel_evidence<K: Ord>(
    evidence: BTreeMap<K, QualityObservation>,
) -> Option<(K, QualityObservation)> {
    let best_score = evidence
        .values()
        .map(|observation| observation.score())
        .max()?;
    let mut winners = evidence
        .into_iter()
        .filter(|(_, observation)| observation.score() == best_score);
    let winner = winners.next()?;
    // Equally strong but different indels are not converted into a structural
    // call. The ambiguity counter records the conflict.
    winners.next().is_none().then_some(winner)
}

/// Stateful pileup builder that emits at most one allele observation per
/// physical fragment and assembly coordinate.
///
/// All named records are retained until `finish`. This is intentionally not an
/// early-flush design: coordinate-sorted supplementary alignments can arrive
/// after both primary mates, and releasing the pair early would make results
/// depend on record order. Memory is therefore proportional to the number of
/// named fragments and their merged covered coordinates. A future bounded
/// implementation must use query-name grouping or an external spill/merge,
/// never infer completeness merely from seeing first and last primary mates.
pub struct FragmentPileupAccumulator {
    pileup: Vec<PileupPosition>,
    settings: PileupSettings,
    coordinate_mode: PileupCoordinateMode,
    stats: PileupRunStats,
    pending: BTreeMap<Vec<u8>, FragmentContribution>,
}

impl FragmentPileupAccumulator {
    pub fn new(
        reference: &[u8],
        settings: PileupSettings,
        coordinate_mode: PileupCoordinateMode,
    ) -> Self {
        Self {
            pileup: init_pileup(reference),
            settings,
            coordinate_mode,
            stats: PileupRunStats::default(),
            pending: BTreeMap::new(),
        }
    }

    pub fn push_record(&mut self, record: &RecordBuf) {
        self.stats.total_records = self.stats.total_records.saturating_add(1);
        let flags = record.flags();

        if flags.is_unmapped() || flags.is_duplicate() || flags.is_qc_fail() {
            self.stats.filtered_records = self.stats.filtered_records.saturating_add(1);
            return;
        }

        let is_sec_supp = flags.is_secondary() || flags.is_supplementary();
        if is_sec_supp {
            self.stats.secondary_or_supp_records =
                self.stats.secondary_or_supp_records.saturating_add(1);
        } else {
            self.stats.primary_records = self.stats.primary_records.saturating_add(1);
        }

        if record.alignment_start().is_none() {
            self.stats.filtered_records = self.stats.filtered_records.saturating_add(1);
            return;
        }

        let contribution =
            record_contribution(&self.pileup, record, self.settings, self.coordinate_mode);

        let name = record
            .name()
            .map(|name| name.iter().copied().collect::<Vec<_>>());
        let Some(name) = name else {
            contribution.apply(&mut self.pileup);
            return;
        };

        self.pending
            .entry(name)
            .or_default()
            .merge(contribution, &mut self.stats);
    }

    pub fn finish(mut self) -> (Vec<PileupPosition>, PileupRunStats) {
        for contribution in std::mem::take(&mut self.pending).into_values() {
            contribution.apply(&mut self.pileup);
        }
        (self.pileup, self.stats)
    }
}

fn record_contribution(
    pileup: &[PileupPosition],
    record: &RecordBuf,
    settings: PileupSettings,
    coordinate_mode: PileupCoordinateMode,
) -> FragmentContribution {
    let Some(alignment_start) = record.alignment_start() else {
        return FragmentContribution::default();
    };
    let flags = record.flags();
    let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
    let sequence = record.sequence().as_ref();
    let quality_scores = record.quality_scores().as_ref();
    let alignment_start = usize::from(alignment_start).saturating_sub(1);

    let should_revcomp = flags.is_reverse_complemented()
        && should_reverse_complement_with_mode(
            pileup,
            record.cigar().as_ref(),
            sequence,
            quality_scores,
            alignment_start,
            settings.min_baseq,
            coordinate_mode,
        );
    let revcomp_sequence: Vec<u8>;
    let reversed_qualities: Vec<u8>;
    let (sequence, quality_scores) = if should_revcomp {
        revcomp_sequence = reverse_complement_sequence(sequence);
        reversed_qualities = quality_scores.iter().copied().rev().collect();
        (revcomp_sequence.as_slice(), reversed_qualities.as_slice())
    } else {
        (sequence, quality_scores)
    };

    let mut contribution = FragmentContribution::default();
    let is_sec_supp = flags.is_secondary() || flags.is_supplementary();
    let allele_eligible = !flags.is_secondary() && mapq >= settings.min_mapq;
    let mut ref_pos = alignment_start;
    let mut read_pos = 0usize;

    for op in record.cigar().as_ref() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    let raw_pos = ref_pos.saturating_add(i);
                    let Some(pos) =
                        project_reference_position(raw_pos, pileup.len(), coordinate_mode)
                    else {
                        continue;
                    };
                    contribution.observe(pos, mapq, is_sec_supp);
                    if !allele_eligible {
                        continue;
                    }
                    let query_pos = read_pos.saturating_add(i);
                    let Some(&raw_base) = sequence.get(query_pos) else {
                        continue;
                    };
                    let baseq = quality_scores.get(query_pos).copied().unwrap_or(0);
                    if baseq >= settings.min_baseq {
                        contribution.add_base(pos, normalize_base(raw_base), baseq, mapq);
                    }
                }
                ref_pos = ref_pos.saturating_add(len);
                read_pos = read_pos.saturating_add(len);
            }
            Kind::Insertion => {
                let start = read_pos;
                let end = read_pos.saturating_add(len).min(sequence.len());
                let inserted = sequence[start..end]
                    .iter()
                    .map(|base| normalize_base(*base))
                    .collect::<Vec<_>>();
                if allele_eligible && !inserted.is_empty() {
                    let mean_baseq = mean_u8(quality_scores.get(start..end).unwrap_or(&[]));
                    if mean_baseq >= settings.min_baseq {
                        if let Some(anchor) =
                            project_anchor_before(ref_pos, pileup.len(), coordinate_mode)
                        {
                            contribution.add_insertion(
                                anchor,
                                String::from_utf8_lossy(&inserted).to_string(),
                                mean_baseq,
                                mapq,
                            );
                        }
                    }
                }
                read_pos = read_pos.saturating_add(len);
            }
            Kind::Deletion | Kind::Skip => {
                for i in 0..len {
                    let raw_pos = ref_pos.saturating_add(i);
                    if let Some(pos) =
                        project_reference_position(raw_pos, pileup.len(), coordinate_mode)
                    {
                        contribution.observe(pos, mapq, is_sec_supp);
                    }
                }
                if allele_eligible && len > 0 {
                    if let Some(anchor) =
                        project_anchor_before(ref_pos, pileup.len(), coordinate_mode)
                    {
                        contribution.add_deletion(anchor, len, settings.min_baseq, mapq);
                    }
                }
                ref_pos = ref_pos.saturating_add(len);
            }
            Kind::SoftClip => {
                read_pos = read_pos.saturating_add(len);
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    if coordinate_mode == PileupCoordinateMode::Circular && allele_eligible {
        recover_circular_terminal_softclips(
            &mut contribution,
            pileup,
            record.cigar().as_ref(),
            sequence,
            quality_scores,
            alignment_start,
            mapq,
            is_sec_supp,
            settings.min_baseq,
        );
    }

    contribution
}

fn project_reference_position(
    raw_pos: usize,
    reference_len: usize,
    coordinate_mode: PileupCoordinateMode,
) -> Option<usize> {
    if reference_len == 0 {
        return None;
    }
    match coordinate_mode {
        PileupCoordinateMode::Linear => (raw_pos < reference_len).then_some(raw_pos),
        PileupCoordinateMode::Circular => Some(raw_pos % reference_len),
    }
}

fn project_anchor_before(
    raw_ref_pos: usize,
    reference_len: usize,
    coordinate_mode: PileupCoordinateMode,
) -> Option<usize> {
    if reference_len == 0 {
        return None;
    }
    match coordinate_mode {
        PileupCoordinateMode::Linear => Some(raw_ref_pos.saturating_sub(1).min(reference_len - 1)),
        PileupCoordinateMode::Circular => {
            Some(raw_ref_pos.wrapping_add(reference_len - 1) % reference_len)
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn recover_circular_terminal_softclips(
    contribution: &mut FragmentContribution,
    pileup: &[PileupPosition],
    cigar: &[Op],
    sequence: &[u8],
    quality_scores: &[u8],
    alignment_start: usize,
    mapq: u8,
    is_sec_supp: bool,
    min_baseq: u8,
) {
    if pileup.is_empty() || sequence.is_empty() {
        return;
    }
    let (leading, trailing) = terminal_softclip_lengths(cigar);
    let reference_end = alignment_start.saturating_add(cigar_reference_span(cigar));

    if leading >= MIN_CIRCULAR_SOFTCLIP_ANCHOR && alignment_start == 0 {
        let clip_len = leading.min(sequence.len());
        let reference_start = pileup.len().wrapping_sub(clip_len % pileup.len()) % pileup.len();
        recover_supported_softclip(
            contribution,
            pileup,
            sequence,
            quality_scores,
            0,
            clip_len,
            reference_start,
            mapq,
            is_sec_supp,
            min_baseq,
        );
    }

    if trailing >= MIN_CIRCULAR_SOFTCLIP_ANCHOR && reference_end >= pileup.len() {
        let clip_len = trailing.min(sequence.len());
        let read_start = sequence.len().saturating_sub(clip_len);
        let reference_start = reference_end % pileup.len();
        recover_supported_softclip(
            contribution,
            pileup,
            sequence,
            quality_scores,
            read_start,
            clip_len,
            reference_start,
            mapq,
            is_sec_supp,
            min_baseq,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn recover_supported_softclip(
    contribution: &mut FragmentContribution,
    pileup: &[PileupPosition],
    sequence: &[u8],
    quality_scores: &[u8],
    read_start: usize,
    clip_len: usize,
    reference_start: usize,
    mapq: u8,
    is_sec_supp: bool,
    min_baseq: u8,
) {
    let mut compared = 0usize;
    let mut matches = 0usize;
    for offset in 0..clip_len {
        let query_pos = read_start.saturating_add(offset);
        let Some(&raw_base) = sequence.get(query_pos) else {
            continue;
        };
        let baseq = quality_scores.get(query_pos).copied().unwrap_or(0);
        if baseq < min_baseq {
            continue;
        }
        let pos = reference_start.saturating_add(offset) % pileup.len();
        compared = compared.saturating_add(1);
        if normalize_base(raw_base) == pileup[pos].ref_base {
            matches = matches.saturating_add(1);
        }
    }
    if compared < MIN_CIRCULAR_SOFTCLIP_ANCHOR
        || matches.saturating_mul(5) < compared.saturating_mul(4)
    {
        return;
    }

    for offset in 0..clip_len {
        let query_pos = read_start.saturating_add(offset);
        let Some(&raw_base) = sequence.get(query_pos) else {
            continue;
        };
        let baseq = quality_scores.get(query_pos).copied().unwrap_or(0);
        if baseq < min_baseq {
            continue;
        }
        let pos = reference_start.saturating_add(offset) % pileup.len();
        contribution.observe(pos, mapq, is_sec_supp);
        contribution.add_base(pos, normalize_base(raw_base), baseq, mapq);
    }
}

fn terminal_softclip_lengths(cigar: &[Op]) -> (usize, usize) {
    let leading = cigar
        .iter()
        .find(|op| !matches!(op.kind(), Kind::HardClip | Kind::Pad))
        .filter(|op| op.kind() == Kind::SoftClip)
        .map_or(0, |op| op.len());
    let trailing = cigar
        .iter()
        .rev()
        .find(|op| !matches!(op.kind(), Kind::HardClip | Kind::Pad))
        .filter(|op| op.kind() == Kind::SoftClip)
        .map_or(0, |op| op.len());
    (leading, trailing)
}

fn cigar_reference_span(cigar: &[Op]) -> usize {
    cigar
        .iter()
        .filter(|op| {
            matches!(
                op.kind(),
                Kind::Match
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
                    | Kind::Deletion
                    | Kind::Skip
            )
        })
        .fold(0usize, |span, op| span.saturating_add(op.len()))
}

#[allow(clippy::too_many_arguments)]
fn should_reverse_complement_with_mode(
    pileup: &[PileupPosition],
    cigar: &[Op],
    sequence: &[u8],
    quality_scores: &[u8],
    alignment_start: usize,
    min_baseq: u8,
    coordinate_mode: PileupCoordinateMode,
) -> bool {
    if sequence.is_empty() {
        return false;
    }

    let mut ref_pos = alignment_start;
    let mut read_pos = 0usize;
    let mut direct_matches = 0u32;
    let mut revcomp_matches = 0u32;
    let mut compared = 0u32;

    'ops: for op in cigar {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    if compared >= 128 {
                        break 'ops;
                    }

                    let raw_pos = ref_pos.saturating_add(i);
                    let q = read_pos.saturating_add(i);

                    let Some(p) =
                        project_reference_position(raw_pos, pileup.len(), coordinate_mode)
                    else {
                        continue;
                    };
                    if q >= sequence.len() {
                        continue;
                    }

                    let rq = sequence.len().saturating_sub(1).saturating_sub(q);
                    let direct_bq = quality_scores.get(q).copied().unwrap_or(0);
                    let rev_bq = quality_scores.get(rq).copied().unwrap_or(0);

                    if direct_bq < min_baseq && rev_bq < min_baseq {
                        continue;
                    }

                    let ref_base = pileup[p].ref_base;

                    if direct_bq >= min_baseq && normalize_base(sequence[q]) == ref_base {
                        direct_matches = direct_matches.saturating_add(1);
                    }

                    if rev_bq >= min_baseq && reverse_complement_base(sequence[rq]) == ref_base {
                        revcomp_matches = revcomp_matches.saturating_add(1);
                    }

                    compared = compared.saturating_add(1);
                }

                ref_pos = ref_pos.saturating_add(len);
                read_pos = read_pos.saturating_add(len);
            }
            Kind::Insertion | Kind::SoftClip => {
                read_pos = read_pos.saturating_add(len);
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos = ref_pos.saturating_add(len);
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    // Reverse-complementing a sequence that is already in reference orientation will
    // catastrophically inflate entropy (A<->T / C<->G mixtures) across the genome.
    //
    // We therefore only flip when there is *strong* evidence that the stored sequence is
    // in the original (uncomplemented) orientation despite the REVERSE_COMPLEMENTED flag.
    //
    // Using a strict threshold also avoids random flips when the provided reference FASTA
    // does not match the reference used to generate the BAM/CRAM alignments.
    if compared < 8 {
        return false;
    }

    // Require >=80% of the sampled bases to match after reverse-complementing.
    // (revcomp_matches / compared) >= 0.8, in integer arithmetic.
    let revcomp_strong = u64::from(revcomp_matches) * 5 >= u64::from(compared) * 4;
    revcomp_strong && revcomp_matches > direct_matches
}

pub(crate) fn reverse_complement_sequence(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .copied()
        .map(reverse_complement_base)
        .collect()
}

pub(crate) fn reverse_complement_base(base: u8) -> u8 {
    match normalize_base(base) {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

pub fn quality_weight(baseq: u8, mapq: u8) -> f64 {
    f64::from(baseq.min(40)) * f64::from(mapq.min(60))
}

pub fn shannon_entropy(weights: &[f64; 5]) -> f64 {
    let total = weights.iter().sum::<f64>();
    if total <= f64::EPSILON {
        return 0.0;
    }

    weights
        .iter()
        .copied()
        .filter(|w| *w > 0.0)
        .map(|w| {
            let p = w / total;
            -(p * p.log2())
        })
        .sum()
}

pub fn histogram_median(hist: &[u32; 61]) -> f64 {
    let total: u32 = hist.iter().sum();
    if total == 0 {
        return 0.0;
    }

    let target = total.div_ceil(2);
    let mut running = 0u32;
    for (i, c) in hist.iter().copied().enumerate() {
        running = running.saturating_add(c);
        if running >= target {
            return i as f64;
        }
    }

    60.0
}

pub fn base_to_index(base: u8) -> usize {
    match normalize_base(base) {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

pub fn index_to_base(index: usize) -> u8 {
    match index {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

pub fn normalize_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' | b'N' => base.to_ascii_uppercase(),
        b'U' => b'T',
        _ => b'N',
    }
}

pub fn dominant_base(base_weights: &[f64; 5]) -> (usize, f64) {
    base_weights
        .iter()
        .copied()
        .enumerate()
        .max_by(|(ia, wa), (ib, wb)| {
            if wa == wb {
                base_tie_break(*ia, *ib)
            } else {
                wa.partial_cmp(wb).unwrap_or(Ordering::Equal)
            }
        })
        .unwrap_or((4, 0.0))
}

fn base_tie_break(lhs_idx: usize, rhs_idx: usize) -> Ordering {
    // Higher rank wins in ties: A > C > G > T > N.
    let rank = |idx: usize| -> i32 {
        match idx {
            0 => 5,
            1 => 4,
            2 => 3,
            3 => 2,
            _ => 1,
        }
    };

    rank(lhs_idx).cmp(&rank(rhs_idx))
}

fn mean_u8(values: &[u8]) -> u8 {
    if values.is_empty() {
        return 0;
    }

    let sum: u64 = values.iter().map(|v| u64::from(*v)).sum();
    let mean = sum / u64::try_from(values.len()).unwrap_or(1);
    u8::try_from(mean).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;
    use noodles_sam::alignment::{
        record::{cigar::op::Kind, Flags, MappingQuality},
        record_buf::{Cigar, QualityScores, Sequence},
        RecordBuf,
    };

    use super::*;

    fn mapped_record(
        name: &str,
        flags: Flags,
        start: usize,
        cigar: Cigar,
        sequence: &[u8],
        qualities: Vec<u8>,
    ) -> RecordBuf {
        RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(start).expect("valid position"))
            .set_mapping_quality(MappingQuality::new(60).expect("valid MAPQ"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(sequence))
            .set_quality_scores(QualityScores::from(qualities))
            .build()
    }

    #[test]
    fn test_quality_weight() {
        assert_eq!(quality_weight(20, 30), 600.0);
        assert_eq!(quality_weight(60, 100), 2400.0);
    }

    #[test]
    fn test_shannon_entropy() {
        let low = shannon_entropy(&[10.0, 0.0, 0.0, 0.0, 0.0]);
        let high = shannon_entropy(&[5.0, 5.0, 0.0, 0.0, 0.0]);
        assert!(high > low);
    }

    #[test]
    fn test_dominant_base_tie_break() {
        let (idx, _) = dominant_base(&[10.0, 10.0, 10.0, 10.0, 0.0]);
        assert_eq!(index_to_base(idx), b'A');
    }

    #[test]
    fn fragment_accumulator_handles_match_insertion_and_deletion() {
        let reference = b"ACGTACGT";
        let cigar: Cigar = [
            Op::new(Kind::Match, 4),
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::Match, 2),
            Op::new(Kind::Deletion, 1),
        ]
        .into_iter()
        .collect();
        let record = mapped_record("indels", Flags::empty(), 1, cigar, b"ACGTGGTT", vec![30; 8]);
        let mut accumulator = FragmentPileupAccumulator::new(
            reference,
            PileupSettings {
                min_mapq: 20,
                min_baseq: 20,
            },
            PileupCoordinateMode::Linear,
        );
        accumulator.push_record(&record);
        let (pileup, stats) = accumulator.finish();
        assert_eq!(stats.total_records, 1);
        assert_eq!(pileup[3].insertions["GG"].count, 1);
        assert_eq!(pileup[5].deletions[&1].count, 1);
        assert_eq!(pileup[5].depth, 1);
        assert_eq!(pileup[6].total_obs, 1);
    }

    #[test]
    fn fragment_accumulator_accepts_both_reverse_read_sequence_conventions() {
        let reference = b"ACGAACGA";
        for sequence in [b"TCGTTCGT", b"ACGAACGA"] {
            let cigar: Cigar = [Op::new(Kind::Match, 8)].into_iter().collect();
            let record = mapped_record(
                "reverse",
                Flags::REVERSE_COMPLEMENTED,
                1,
                cigar,
                sequence,
                vec![35; 8],
            );
            let mut accumulator = FragmentPileupAccumulator::new(
                reference,
                PileupSettings {
                    min_mapq: 20,
                    min_baseq: 20,
                },
                PileupCoordinateMode::Linear,
            );
            accumulator.push_record(&record);
            let (pileup, _) = accumulator.finish();
            for (i, &base) in reference.iter().enumerate() {
                assert_eq!(
                    pileup[i].base_counts[base_to_index(base)],
                    1,
                    "position {}",
                    i + 1
                );
            }
        }
    }

    #[test]
    fn test_fragment_accumulator_collapses_overlapping_mates() {
        let reference = b"ACGTACGTACGT";
        let cigar: Cigar = [Op::new(Kind::Match, 8)].into_iter().collect();
        let first = mapped_record(
            "fragment-1",
            Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            1,
            cigar.clone(),
            &reference[..8],
            vec![35; 8],
        );
        let last = mapped_record(
            "fragment-1",
            Flags::SEGMENTED | Flags::LAST_SEGMENT,
            5,
            cigar,
            &reference[4..],
            vec![35; 8],
        );

        let mut accumulator = FragmentPileupAccumulator::new(
            reference,
            PileupSettings {
                min_mapq: 20,
                min_baseq: 20,
            },
            PileupCoordinateMode::Circular,
        );
        accumulator.push_record(&first);
        accumulator.push_record(&last);
        let (pileup, stats) = accumulator.finish();

        assert!(pileup.iter().all(|position| position.depth == 1));
        assert_eq!(stats.overlapping_fragment_observations_collapsed, 4);
        assert_eq!(stats.conflicting_fragment_observations, 0);
    }

    #[test]
    fn test_fragment_conflicts_prefer_baseq_then_mapq_and_tie_to_n() {
        let reference = b"AAAAAAAAAAAA";
        let cigar: Cigar = [Op::new(Kind::Match, 8)].into_iter().collect();
        let high_quality = mapped_record(
            "quality-wins",
            Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            1,
            cigar.clone(),
            b"AAAAAAAA",
            vec![35; 8],
        );
        let low_quality = mapped_record(
            "quality-wins",
            Flags::SEGMENTED | Flags::LAST_SEGMENT,
            5,
            cigar.clone(),
            b"CCCCCCCC",
            vec![25; 8],
        );
        let tied_a = mapped_record(
            "exact-tie",
            Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            1,
            cigar.clone(),
            b"AAAAAAAA",
            vec![30; 8],
        );
        let tied_c = mapped_record(
            "exact-tie",
            Flags::SEGMENTED | Flags::LAST_SEGMENT,
            1,
            cigar,
            b"CAAAAAAA",
            vec![30; 8],
        );

        let mut accumulator = FragmentPileupAccumulator::new(
            reference,
            PileupSettings {
                min_mapq: 20,
                min_baseq: 20,
            },
            PileupCoordinateMode::Circular,
        );
        for record in [&high_quality, &low_quality, &tied_a, &tied_c] {
            accumulator.push_record(record);
        }
        let (pileup, stats) = accumulator.finish();

        // At positions 5-8 the higher-baseQ A observation wins over C.
        for position in &pileup[4..8] {
            assert_eq!(position.base_counts[base_to_index(b'A')], 2);
            assert_eq!(position.base_counts[base_to_index(b'C')], 0);
        }
        // Equal baseQ and MAPQ with different alleles becomes one ambiguous
        // molecule, never two allele observations.
        assert_eq!(pileup[0].depth, 2);
        assert_eq!(pileup[0].base_counts[base_to_index(b'A')], 1);
        assert_eq!(pileup[0].base_counts[base_to_index(b'N')], 1);
        assert_eq!(stats.conflicting_fragment_observations, 5);
    }

    #[test]
    fn test_circular_accumulator_folds_alignment_past_origin() {
        let reference = b"ACGTTGCAAGTCCTAG";
        let sequence = [&reference[12..], &reference[..4]].concat();
        let cigar: Cigar = [Op::new(Kind::Match, sequence.len())].into_iter().collect();
        let record = mapped_record(
            "origin-spanning",
            Flags::empty(),
            13,
            cigar,
            &sequence,
            vec![35; sequence.len()],
        );

        let mut accumulator = FragmentPileupAccumulator::new(
            reference,
            PileupSettings {
                min_mapq: 20,
                min_baseq: 20,
            },
            PileupCoordinateMode::Circular,
        );
        accumulator.push_record(&record);
        let (pileup, _) = accumulator.finish();

        for pos in [12usize, 13, 14, 15, 0, 1, 2, 3] {
            assert_eq!(pileup[pos].depth, 1, "position {}", pos + 1);
            assert_eq!(
                pileup[pos].base_counts[base_to_index(reference[pos])],
                1,
                "position {}",
                pos + 1
            );
        }
        assert!(pileup[4..12].iter().all(|position| position.depth == 0));
    }

    #[test]
    fn test_circular_accumulator_recovers_supported_terminal_softclip() {
        let reference = b"ACGTTGCAAGTCCTAGGATCCGTA";
        let sequence = [&reference[12..], &reference[..8]].concat();
        let cigar: Cigar = [Op::new(Kind::Match, 12), Op::new(Kind::SoftClip, 8)]
            .into_iter()
            .collect();
        let record = mapped_record(
            "origin-softclip",
            Flags::empty(),
            13,
            cigar,
            &sequence,
            vec![35; sequence.len()],
        );

        let mut accumulator = FragmentPileupAccumulator::new(
            reference,
            PileupSettings {
                min_mapq: 20,
                min_baseq: 20,
            },
            PileupCoordinateMode::Circular,
        );
        accumulator.push_record(&record);
        let (pileup, _) = accumulator.finish();

        for pos in 0..8 {
            assert_eq!(pileup[pos].depth, 1, "recovered position {}", pos + 1);
            assert_eq!(pileup[pos].base_counts[base_to_index(reference[pos])], 1);
        }
        assert!(pileup[8..12].iter().all(|position| position.depth == 0));
        assert!(pileup[12..].iter().all(|position| position.depth == 1));
    }

    #[test]
    fn test_supplementary_piece_contributes_nonoverlapping_origin_segment_once() {
        let reference = b"ACGTTGCAAGTCCTAGGATCCGTA";
        let cigar: Cigar = [Op::new(Kind::Match, 8)].into_iter().collect();
        let supplementary = mapped_record(
            "split-fragment",
            Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::SUPPLEMENTARY,
            1,
            cigar.clone(),
            &reference[..8],
            vec![35; 8],
        );
        let last = mapped_record(
            "split-fragment",
            Flags::SEGMENTED | Flags::LAST_SEGMENT,
            9,
            cigar.clone(),
            &reference[8..16],
            vec![35; 8],
        );
        let first = mapped_record(
            "split-fragment",
            Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            17,
            cigar,
            &reference[16..],
            vec![35; 8],
        );

        let mut accumulator = FragmentPileupAccumulator::new(
            reference,
            PileupSettings {
                min_mapq: 20,
                min_baseq: 20,
            },
            PileupCoordinateMode::Circular,
        );
        for record in [&supplementary, &last, &first] {
            accumulator.push_record(record);
        }
        let (pileup, stats) = accumulator.finish();

        assert!(pileup.iter().all(|position| position.depth == 1));
        assert_eq!(stats.primary_records, 2);
        assert_eq!(stats.secondary_or_supp_records, 1);
        assert!(pileup[..8]
            .iter()
            .all(|position| position.sec_supp_obs == 1));
    }

    #[test]
    fn test_fragment_accumulator_collapses_identical_mate_indels() {
        let reference = b"ACGTACGT";
        let sequence = b"ACGTGGACGT";
        let cigar: Cigar = [
            Op::new(Kind::Match, 4),
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::Match, 4),
        ]
        .into_iter()
        .collect();
        let first = mapped_record(
            "indel-fragment",
            Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            1,
            cigar.clone(),
            sequence,
            vec![35; sequence.len()],
        );
        let last = mapped_record(
            "indel-fragment",
            Flags::SEGMENTED | Flags::LAST_SEGMENT,
            1,
            cigar,
            sequence,
            vec![35; sequence.len()],
        );

        let mut accumulator = FragmentPileupAccumulator::new(
            reference,
            PileupSettings {
                min_mapq: 20,
                min_baseq: 20,
            },
            PileupCoordinateMode::Circular,
        );
        accumulator.push_record(&first);
        accumulator.push_record(&last);
        let (pileup, _) = accumulator.finish();

        assert!(pileup.iter().all(|position| position.depth == 1));
        assert_eq!(
            pileup[3].insertions.get("GG").map(|stats| stats.count),
            Some(1)
        );
    }

    #[test]
    fn test_secondary_alignment_changes_mapping_diagnostics_not_alleles() {
        let reference = b"AAAAAAAA";
        let cigar: Cigar = [Op::new(Kind::Match, 8)].into_iter().collect();
        let primary = mapped_record(
            "multi-mapper",
            Flags::empty(),
            1,
            cigar.clone(),
            b"AAAAAAAA",
            vec![35; 8],
        );
        let secondary = mapped_record(
            "multi-mapper",
            Flags::SECONDARY,
            1,
            cigar,
            b"CCCCCCCC",
            vec![35; 8],
        );

        let mut accumulator = FragmentPileupAccumulator::new(
            reference,
            PileupSettings {
                min_mapq: 20,
                min_baseq: 20,
            },
            PileupCoordinateMode::Circular,
        );
        accumulator.push_record(&primary);
        accumulator.push_record(&secondary);
        let (pileup, stats) = accumulator.finish();

        for position in pileup {
            assert_eq!(position.depth, 1);
            assert_eq!(position.base_counts[base_to_index(b'A')], 1);
            assert_eq!(position.base_counts[base_to_index(b'C')], 0);
            assert_eq!(position.total_obs, 1);
            assert_eq!(position.primary_obs, 1);
            assert_eq!(position.sec_supp_obs, 0);
        }
        assert_eq!(stats.primary_records, 1);
        assert_eq!(stats.secondary_or_supp_records, 1);
        assert_eq!(stats.overlapping_fragment_observations_collapsed, 8);
    }

    #[test]
    fn test_unmatched_terminal_softclip_is_not_folded() {
        let reference = b"ACGTTGCAAGTCCTAGGATCCGTA";
        let sequence = [&reference[12..], b"NNNNNNNN"].concat();
        let cigar: Cigar = [Op::new(Kind::Match, 12), Op::new(Kind::SoftClip, 8)]
            .into_iter()
            .collect();
        let record = mapped_record(
            "junk-softclip",
            Flags::empty(),
            13,
            cigar,
            &sequence,
            vec![35; sequence.len()],
        );

        let mut accumulator = FragmentPileupAccumulator::new(
            reference,
            PileupSettings {
                min_mapq: 20,
                min_baseq: 20,
            },
            PileupCoordinateMode::Circular,
        );
        accumulator.push_record(&record);
        let (pileup, _) = accumulator.finish();

        assert!(pileup[..12].iter().all(|position| position.depth == 0));
        assert!(pileup[12..].iter().all(|position| position.depth == 1));
    }

    #[test]
    fn test_late_supplementary_is_order_invariant_after_primary_pair() {
        let reference = b"ACGTTGCAAGTCCTAG";
        let cigar: Cigar = [Op::new(Kind::Match, 8)].into_iter().collect();
        let first = mapped_record(
            "order-invariant",
            Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            1,
            cigar.clone(),
            &reference[..8],
            vec![35; 8],
        );
        let last = mapped_record(
            "order-invariant",
            Flags::SEGMENTED | Flags::LAST_SEGMENT,
            9,
            cigar.clone(),
            &reference[8..],
            vec![35; 8],
        );
        let supplementary = mapped_record(
            "order-invariant",
            Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::SUPPLEMENTARY,
            1,
            cigar,
            &reference[..8],
            vec![35; 8],
        );

        let run = |records: [&RecordBuf; 3]| {
            let mut accumulator = FragmentPileupAccumulator::new(
                reference,
                PileupSettings {
                    min_mapq: 20,
                    min_baseq: 20,
                },
                PileupCoordinateMode::Circular,
            );
            for record in records {
                accumulator.push_record(record);
            }
            accumulator.finish()
        };

        let (late_pileup, late_stats) = run([&first, &last, &supplementary]);
        let (early_pileup, early_stats) = run([&supplementary, &last, &first]);

        assert_pileups_equal(&late_pileup, &early_pileup);
        assert_run_stats_equal(&late_stats, &early_stats);
        assert!(late_pileup.iter().all(|position| position.depth == 1));
        assert_eq!(late_stats.total_records, 3);
        assert_eq!(late_stats.primary_records, 2);
        assert_eq!(late_stats.secondary_or_supp_records, 1);
        assert_eq!(late_stats.overlapping_fragment_observations_collapsed, 8);
    }

    fn assert_pileups_equal(left: &[PileupPosition], right: &[PileupPosition]) {
        assert_eq!(left.len(), right.len());
        for (pos, (left, right)) in left.iter().zip(right).enumerate() {
            assert_eq!(left.ref_base, right.ref_base, "ref base at {pos}");
            assert_eq!(left.depth, right.depth, "depth at {pos}");
            assert_eq!(left.primary_obs, right.primary_obs, "primary obs at {pos}");
            assert_eq!(left.total_obs, right.total_obs, "total obs at {pos}");
            assert_eq!(left.sec_supp_obs, right.sec_supp_obs, "sec/supp at {pos}");
            assert_eq!(left.mapq_hist, right.mapq_hist, "MAPQ hist at {pos}");
            assert_eq!(left.base_counts, right.base_counts, "base counts at {pos}");
            assert_eq!(
                left.base_weights, right.base_weights,
                "base weights at {pos}"
            );
            assert_eq!(left.baseq_sums, right.baseq_sums, "baseQ sums at {pos}");
            assert_eq!(left.mapq_sums, right.mapq_sums, "mapQ sums at {pos}");
            assert_eq!(
                indel_signature(&left.insertions),
                indel_signature(&right.insertions),
                "insertions at {pos}"
            );
            assert_eq!(
                deletion_signature(&left.deletions),
                deletion_signature(&right.deletions),
                "deletions at {pos}"
            );
        }
    }

    fn indel_signature(indels: &BTreeMap<String, IndelStats>) -> Vec<(&str, u32, f64, u64, u64)> {
        indels
            .iter()
            .map(|(allele, stats)| {
                (
                    allele.as_str(),
                    stats.count,
                    stats.weight,
                    stats.baseq_sum,
                    stats.mapq_sum,
                )
            })
            .collect()
    }

    fn deletion_signature(
        indels: &BTreeMap<usize, IndelStats>,
    ) -> Vec<(usize, u32, f64, u64, u64)> {
        indels
            .iter()
            .map(|(allele, stats)| {
                (
                    *allele,
                    stats.count,
                    stats.weight,
                    stats.baseq_sum,
                    stats.mapq_sum,
                )
            })
            .collect()
    }

    fn assert_run_stats_equal(left: &PileupRunStats, right: &PileupRunStats) {
        assert_eq!(left.total_records, right.total_records);
        assert_eq!(left.filtered_records, right.filtered_records);
        assert_eq!(left.primary_records, right.primary_records);
        assert_eq!(
            left.secondary_or_supp_records,
            right.secondary_or_supp_records
        );
        assert_eq!(
            left.overlapping_fragment_observations_collapsed,
            right.overlapping_fragment_observations_collapsed
        );
        assert_eq!(
            left.conflicting_fragment_observations,
            right.conflicting_fragment_observations
        );
    }
}
