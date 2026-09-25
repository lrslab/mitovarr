use crate::{consensus::ConsensusResult, repeat::AmbiguityInterval};
use noodles_sam::alignment::{record::cigar::op::Kind, RecordBuf};
use std::collections::BTreeSet;

#[derive(Default)]
pub(crate) struct Signals {
    entries: Vec<(usize, usize, Vec<u8>, &'static str)>,
}

impl Signals {
    pub fn observe(&mut self, record: &RecordBuf, length: usize, min_mapq: u8) {
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
            || record.mapping_quality().map(u8::from).unwrap_or(0) < min_mapq
        {
            return;
        }
        let (Some(start), Some(name)) = (record.alignment_start(), record.name()) else {
            return;
        };
        let name: Vec<u8> = name.iter().copied().collect();
        let mut pos = usize::from(start) - 1;
        for op in record.cigar().as_ref() {
            let size = op.len();
            let reason = match op.kind() {
                Kind::SoftClip if size >= 15 => Some("soft_clip"),
                Kind::Insertion if size >= 20 => Some("large_insertion"),
                Kind::Deletion | Kind::Skip if size >= 20 => Some("large_deletion"),
                _ => None,
            };
            if let Some(reason) = reason {
                let end = if op.kind().consumes_reference() {
                    pos + size
                } else {
                    pos + 1
                };
                self.entries
                    .push((pos.min(length - 1), end.min(length), name.clone(), reason));
            }
            if op.kind().consumes_reference() {
                pos += size;
            }
        }
        if let Some(mate) = record.mate_alignment_start() {
            let a = usize::from(start) - 1;
            let b = usize::from(mate) - 1;
            let gap = a.abs_diff(b).min(length.saturating_sub(a.abs_diff(b)));
            if record.mate_reference_sequence_id() != record.reference_sequence_id() || gap > 1000 {
                self.entries
                    .push((a, pos.min(length), name, "discordant_mate"));
            }
        }
    }

    pub fn windows(
        mut self,
        called: &ConsensusResult,
        ambiguity: &[AmbiguityInterval],
    ) -> (Vec<bool>, Vec<&'static str>) {
        let n = called.sequence.len();
        let mut mask = vec![false; n];
        let mut reasons = BTreeSet::new();
        if n == 0 {
            return (mask, vec![]);
        }
        let mut mark = |start: usize, end: usize| {
            // Circular padding also joins the two halves of an origin event.
            let lo = start as isize - 60;
            let hi = end.max(start + 1) as isize + 60;
            for p in lo..hi {
                mask[p.rem_euclid(n as isize) as usize] = true;
            }
        };
        self.entries.sort_by_key(|x| (x.0, x.1));
        let mut i = 0;
        while i < self.entries.len() {
            let start = self.entries[i].0;
            let mut end = self.entries[i].1;
            let mut names = BTreeSet::new();
            let mut cluster_reasons = BTreeSet::new();
            // First require a concentrated breakpoint. Chaining arbitrary read
            // ends at high coverage can otherwise merge the entire mt genome.
            while i < self.entries.len() && self.entries[i].0 <= start + 6 {
                let e = &self.entries[i];
                names.insert(&e.2);
                cluster_reasons.insert(e.3);
                end = end.max(e.1);
                i += 1;
            }
            let local_depth = called.coverage
                [start.saturating_sub(1)..(end + 1).min(called.coverage.len())]
                .iter()
                .map(|p| p.depth as usize)
                .max()
                .unwrap_or(0);
            if names.len() >= 3 && names.len().saturating_mul(20) >= local_depth {
                mark(
                    called.reference_offsets[start],
                    called.reference_offsets[end],
                );
                reasons.extend(cluster_reasons);
            }
        }
        for a in ambiguity {
            mark(
                called.reference_offsets[a.start - 1],
                called.reference_offsets[a.end],
            );
            reasons.insert("repeat_ambiguity");
        }
        let mut i = 0;
        while i < n {
            if called.sequence[i] != b'N' {
                i += 1;
                continue;
            }
            let start = i;
            while i < n && called.sequence[i] == b'N' {
                i += 1;
            }
            if i - start >= 15 {
                mark(start, i);
                reasons.insert("uncovered_sequence");
            }
        }
        (mask, reasons.into_iter().collect())
    }
}

pub(super) fn rotate_windows(mask: &[bool]) -> (usize, Vec<(usize, usize)>) {
    let n = mask.len();
    if n == 0 || !mask.iter().any(|x| *x) {
        return (0, vec![]);
    }
    let (mut run, mut best, mut end) = (0, 0, 0);
    for p in 0..2 * n {
        if mask[p % n] {
            run = 0;
        } else {
            run += 1;
        }
        if run > best {
            best = run;
            end = p + 1;
        }
    }
    let cut = (end.saturating_sub(best / 2)) % n;
    let mut intervals: Vec<(usize, usize)> = vec![];
    let mut p = 0;
    while p < n {
        if !mask[(cut + p) % n] {
            p += 1;
            continue;
        }
        let start = p;
        while p < n && mask[(cut + p) % n] {
            p += 1;
        }
        if let Some(last) = intervals.last_mut().filter(|last| start <= last.1 + 150) {
            last.1 = p;
        } else {
            intervals.push((start, p));
        }
    }
    (cut, intervals)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::CoverageRecord;

    fn backbone() -> ConsensusResult {
        ConsensusResult {
            sequence: vec![b'A'; 2000],
            reference_offsets: (0..=2000).collect(),
            variants: vec![],
            uncertainty: vec![],
            base_votes: vec![],
            coverage: (1..=2000)
                .map(|pos| CoverageRecord {
                    pos,
                    depth: 1000,
                    mapq_median: 60.0,
                    entropy: 0.0,
                    is_repeat_ambiguous: false,
                })
                .collect(),
        }
    }

    #[test]
    fn scattered_clips_do_not_chain_into_whole_genome_window() {
        let mut signals = Signals::default();
        for pos in 0..2000 {
            signals.entries.push((
                pos,
                pos + 1,
                format!("noise{pos}").into_bytes(),
                "soft_clip",
            ));
        }
        for id in 0..60 {
            signals
                .entries
                .push((700, 701, format!("junction{id}").into_bytes(), "soft_clip"));
        }
        let (mask, _) = signals.windows(&backbone(), &[]);
        assert!(mask[700]);
        assert!(mask.iter().filter(|&&x| x).count() < 140);
        assert!(!mask[1200]);
    }

    #[test]
    fn repeated_alignments_of_one_name_do_not_create_a_breakpoint() {
        let signals = Signals {
            entries: (0..100)
                .map(|_| (700, 701, b"same_fragment".to_vec(), "soft_clip"))
                .collect(),
        };
        assert!(!signals.windows(&backbone(), &[]).0.iter().any(|&x| x));
    }
}
