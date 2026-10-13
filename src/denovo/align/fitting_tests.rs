// Frozen six-buffer implementation provides an independent layout oracle.
use super::*;

fn original_fitting(query: &[u8], quality: &[u8], target: &[u8]) -> Vec<Hit> {
    let empty = Cell {
        score: NEG,
        start: 0,
        edits: 0,
    };
    let mut m: Vec<_> = (0..=target.len())
        .map(|start| Cell {
            score: 0.0,
            start,
            edits: 0,
        })
        .collect();
    let mut ins = vec![empty; m.len()];
    let mut del = vec![empty; m.len()];
    let mut next_m = vec![empty; m.len()];
    let mut next_i = vec![empty; m.len()];
    let mut next_d = vec![empty; m.len()];
    let max_edits = (query.len() / 10).max(3);
    if query.len() > target.len() + max_edits {
        return Vec::new();
    }
    for (i, &base) in query.iter().enumerate() {
        next_m.fill(empty);
        next_i.fill(empty);
        next_d.fill(empty);
        next_i[0] = best(m[0].advance(GAP_OPEN, 1), ins[0].advance(GAP_EXTEND, 1));
        let match_score = emission(base, base, quality[i]);
        let mismatch_base = if base == b'A' { b'C' } else { b'A' };
        let mismatch_score = emission(base, mismatch_base, quality[i]);
        let lo = (i + 1).saturating_sub(max_edits).max(1);
        let hi = (target.len() + max_edits)
            .saturating_sub(query.len() - i - 1)
            .min(target.len());
        for j in lo..=hi {
            next_m[j] = best(best(m[j - 1], ins[j - 1]), del[j - 1]).advance(
                if !b"ACGT".contains(&target[j - 1]) {
                    0.25f64.ln()
                } else if base == target[j - 1] {
                    match_score
                } else {
                    mismatch_score
                },
                usize::from(base != target[j - 1]),
            );
            next_i[j] = best(m[j].advance(GAP_OPEN, 1), ins[j].advance(GAP_EXTEND, 1));
            next_d[j] = best(
                next_m[j - 1].advance(GAP_OPEN, 1),
                next_d[j - 1].advance(GAP_EXTEND, 1),
            );
        }
        std::mem::swap(&mut m, &mut next_m);
        std::mem::swap(&mut ins, &mut next_i);
        std::mem::swap(&mut del, &mut next_d);
    }
    let mut hits: Vec<_> = (1..=target.len())
        .filter_map(|end| {
            let cell = best(best(m[end], ins[end]), del[end]);
            (cell.score > NEG / 2.0 && end > cell.start && cell.edits <= max_edits).then_some(Hit {
                start: cell.start,
                end,
                reverse: false,
                logp: cell.score,
                edits: cell.edits,
                traceback: None,
            })
        })
        .collect();
    hits.sort_by(|a, b| {
        b.logp
            .total_cmp(&a.logp)
            .then(a.edits.cmp(&b.edits))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    let maximum = hits.first().map_or(NEG, |h| h.logp);
    let mut starts = BTreeSet::new();
    hits.retain(|h| h.logp >= maximum - 16.0 && starts.insert(h.start));
    hits
}

fn signature(hits: &[Hit]) -> Vec<(usize, usize, bool, u64, usize)> {
    hits.iter()
        .map(|h| {
            assert!(h.traceback.is_none());
            (h.start, h.end, h.reverse, h.logp.to_bits(), h.edits)
        })
        .collect()
}

fn check(query: &[u8], quality: &[u8], target: &[u8]) {
    let read = Record {
        header: "@oracle/1".into(),
        sequence: query.to_vec(),
        plus: "+".into(),
        quality: quality.to_vec(),
    };
    let prepared = Prepared::new(&read);
    for oriented in &prepared.orientations {
        assert_eq!(
            signature(&fitting_candidates(oriented, target)),
            signature(&original_fitting(
                &oriented.sequence,
                &oriented.quality,
                target
            )),
            "query={:?}, target={:?}",
            oriented.sequence,
            target,
        );
        let ideal: f64 = oriented
            .sequence
            .iter()
            .zip(&oriented.quality)
            .map(|(&b, &q)| emission(b, b, q))
            .sum();
        assert_eq!(oriented.ideal.to_bits(), ideal.to_bits());
    }
}

#[test]
fn rolling_fitting_preserves_band_boundaries_ties_and_unknowns() {
    for length in [0usize, 1, 2, 3, 4, 16, 30, 31, 63, 100, 151] {
        for delta in [-20isize, -4, -3, -1, 0, 1, 3, 4, 20, 150] {
            let target_len = length.saturating_add_signed(delta);
            for alphabet in [b"A".as_slice(), b"AT", b"ACGTN"] {
                let query: Vec<_> = alphabet.iter().copied().cycle().take(length).collect();
                let target: Vec<_> = alphabet.iter().copied().cycle().take(target_len).collect();
                let quality: Vec<_> = (0..length)
                    .map(|i| [b'!', b'"', b'+', b'I', b'~'][i % 5])
                    .collect();
                check(&query, &quality, &target);
            }
        }
    }
}

#[test]
fn rolling_fitting_matches_original_for_mutated_reads_and_competing_targets() {
    let mut state = 0x7502939a813789u64;
    let mut random = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state as usize
    };
    for case in 0..400 {
        let length = 1 + random() % 180;
        let template: Vec<_> = (0..length).map(|_| b"ACGT"[random() % 4]).collect();
        let flank: Vec<_> = (0..(random() % 70))
            .map(|_| b"ACGTN"[random() % 5])
            .collect();
        let mut target = [flank.as_slice(), &template, &flank].concat();
        if case % 3 == 0 {
            target.extend_from_slice(&template);
            target.extend_from_slice(&flank);
        }
        let mut query = template.clone();
        for _ in 0..(random() % (edit_budget(length) + 3)) {
            let at = random() % query.len().max(1);
            match random() % 3 {
                0 if !query.is_empty() => {
                    query[at] = b"ACGTN"[random() % 5];
                }
                1 if query.len() > 1 => {
                    query.remove(at);
                }
                _ => query.insert(at, b"ACGTN"[random() % 5]),
            }
        }
        if case % 7 == 0 {
            query.extend_from_slice(&flank);
        }
        let quality: Vec<_> = (0..query.len())
            .map(|_| 33 + (random() % 61) as u8)
            .collect();
        check(&query, &quality, &target);
    }
}

#[test]
#[ignore = "manual comparison of the original and in-place affine fitting kernels"]
fn fitting_kernel_benchmark() {
    use std::{hint::black_box, time::Instant};
    let target: Vec<_> = (0usize..900)
        .map(|i| b"ACGT"[(i * i / 7 + i * 3) % 4])
        .collect();
    let mut sequence = target[201..351].to_vec();
    sequence.insert(55, b'A');
    sequence.remove(92);
    sequence[132] = b'N';
    let read = Record {
        header: "@bench".into(),
        quality: vec![b'I'; sequence.len()],
        sequence,
        plus: "+".into(),
    };
    let prepared = Prepared::new(&read);
    let query = &prepared.orientations[0];
    assert_eq!(
        signature(&original_fitting(&query.sequence, &query.quality, &target)),
        signature(&fitting_candidates(query, &target))
    );
    let start = Instant::now();
    for _ in 0..300 {
        black_box(original_fitting(
            black_box(&query.sequence),
            &query.quality,
            black_box(&target),
        ));
    }
    let original = start.elapsed();
    let start = Instant::now();
    for _ in 0..300 {
        black_box(fitting_candidates(black_box(query), black_box(&target)));
    }
    eprintln!("original={original:?} rolling={:?}", start.elapsed());
}
