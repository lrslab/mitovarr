use super::*;
use std::collections::BTreeMap;

fn same_profiles(expected: &Placements, actual: &Placements) {
    assert_eq!(actual.complete, expected.complete);
    assert_eq!(actual.exhausted, expected.exhausted);
    assert_eq!(actual.hits.len(), expected.hits.len());
    for (a, b) in actual.hits.iter().zip(&expected.hits) {
        assert_eq!(
            (a.start, a.end, a.reverse, a.edits),
            (b.start, b.end, b.reverse, b.edits)
        );
        assert_eq!(a.logp.to_bits(), b.logp.to_bits());
        assert_eq!(a.traceback, b.traceback);
    }
}

#[test]
fn certified_endpoints_and_fallbacks_match_long_mixed_quality_profiles() {
    let mut seed = 91753u64;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 32) as usize
    };
    for case in 0..72 {
        let length = [80, 101, 150][case % 3];
        let target: Vec<_> = (0..length + 40)
            .map(|i| {
                if case % 5 == 0 {
                    b"ACGTAC"[i % 6]
                } else {
                    b"ACGT"[next() % 4]
                }
            })
            .collect();
        let mut query = target[20..20 + length].to_vec();
        for _ in 0..case % 10 {
            let p = 25 + next() % (query.len() - 25);
            query[p] = b"ACGTN"[next() % 5];
        }
        if case % 4 == 0 {
            query.insert(40, b'G');
        }
        if case % 7 == 0 {
            query.remove(60);
        }
        let quality: Vec<_> = (0..query.len()).map(|_| b"!$+5:FI~"[next() % 8]).collect();
        let constraint = (case % 2 == 0).then_some(Constraint {
            query: 3,
            target: 23,
            length: 12,
        });
        let budget = if case % 6 == 0 { 140_000 } else { MAX_CELLS };
        let edits = super::super::edit_budget(query.len());
        same_profiles(
            &fit(&query, &quality, &target, constraint, edits, budget),
            &fit_cached(
                &query,
                &quality,
                &target,
                constraint,
                edits,
                budget,
                &mut Subproblems::default(),
            ),
        );
    }
}

#[test]
fn unrestricted_certificate_is_used_and_preserves_each_accepted_endpoint() {
    let mut seed = 119u64;
    let target: Vec<_> = (0..160)
        .map(|_| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(seed >> 32) as usize % 4]
        })
        .collect();
    let query = &target[20..120];
    let quality = vec![b'I'; query.len()];
    let expected = fit(query, &quality, &target, None, 10, MAX_CELLS);
    let emissions: Vec<_> = query
        .iter()
        .map(|&q| b"ACGTN".map(|t| emission(q, t, b'I')))
        .collect();
    let bases: Vec<_> = target
        .iter()
        .map(|&b| b"ACGT".iter().position(|&c| c == b).unwrap())
        .collect();
    let mut certified_count = 0;
    let mut fallback_count = 0;
    for start in 10..=30 {
        let length = query.len() + 10;
        let bounds = edit_bounds(query, &target[start..start + length], 10, None);
        let found = certified::fit_start(
            query, &target, &emissions, &bases, start, length, None, 10, &bounds, 1e-7,
        );
        fallback_count += found.needs_profile.iter().filter(|&&v| v).count();
        for hit in found.hits.into_iter().flatten() {
            let original = expected
                .hits
                .iter()
                .find(|h| (h.start, h.end) == (hit.start, hit.end))
                .unwrap();
            assert_eq!(hit.logp.to_bits(), original.logp.to_bits());
            assert_eq!(hit.edits, original.edits);
            assert_eq!(hit.traceback, original.traceback);
            certified_count += 1;
        }
    }
    assert!(certified_count > 100);
    assert_eq!(fallback_count, 0);
}

#[test]
fn shared_start_subproblems_preserve_tracebacks_boundaries_and_work_limits() {
    let query = b"ACGTAC";
    for quality in [b"IIIIII", b"!!I!I!", b"$I$$I$"] {
        for edits in [1, 2, 3] {
            let mut cache = Subproblems::default();
            for target in [
                b"TTACGTACGGACGTACAA".as_slice(),
                b"ATACGTACGGACGTACAT",
                b"ACGTAC",
                b"ACGTACGGACGT",
            ] {
                for constraint in [
                    None,
                    Some(Constraint {
                        query: 1,
                        target: 3,
                        length: 2,
                    }),
                ] {
                    let expected = fit(query, quality, target, constraint, edits, MAX_CELLS);
                    same_profiles(
                        &expected,
                        &fit_cached(
                            query, quality, target, constraint, edits, MAX_CELLS, &mut cache,
                        ),
                    );
                    for budget in (0..1200).step_by(29).chain([MAX_CELLS]) {
                        same_profiles(
                            &fit(query, quality, target, constraint, edits, budget),
                            &fit_cached(
                                query, quality, target, constraint, edits, budget, &mut cache,
                            ),
                        );
                    }
                }
            }
            assert!(!cache.entries.is_empty());
            assert!(cache.entries.values().all(|e| e.work > 0));
        }
    }
}

#[test]
fn subproblem_cache_capacity_never_changes_the_search() {
    let query = b"ACGTAC";
    let target = b"TTACGTACGGACGTACAA";
    let mut cache = Subproblems {
        hits: 8192,
        ..Subproblems::default()
    };
    same_profiles(
        &fit(query, b"IIIIII", target, None, 2, MAX_CELLS),
        &fit_cached(query, b"IIIIII", target, None, 2, MAX_CELLS, &mut cache),
    );
    assert!(cache.entries.values().all(|entry| entry.hits.is_empty()));
}

#[test]
fn warm_subproblem_cache_preserves_the_endpoint_profile_limit() {
    let target = b"ACGT".repeat(2000);
    let expected = fit(b"ACGT", b"IIII", &target, None, 1, MAX_CELLS);
    assert!(expected.exhausted && !expected.complete);
    assert_eq!(expected.hits.len(), MAX_PROFILES + 1);
    let mut cache = Subproblems::default();
    for _ in 0..2 {
        same_profiles(
            &expected,
            &fit_cached(b"ACGT", b"IIII", &target, None, 1, MAX_CELLS, &mut cache),
        );
    }
}

#[test]
fn cached_profiles_keep_the_query_target_constraint() {
    use crate::{
        denovo::align::{trace_alignment, Index, Prepared},
        io::fastq::Record,
    };
    let read = Record {
        header: "@cache".into(),
        sequence: b"AC".to_vec(),
        quality: vec![b'I'; 2],
        plus: "+".into(),
    };
    let mut prepared = Prepared::new(&read);
    let index = Index::new(b"ACAC", false, 2);
    for target in [0, 2, 0] {
        let found = index.profile_anchor(&mut prepared, target, false, 0, 1);
        assert!(found.complete && !found.hits.is_empty());
        for hit in found.hits {
            let trace = trace_alignment(&read, &hit, &index.sequence, 0);
            assert!(trace.matching_pairs.contains(&(0, target)));
        }
    }
}

// Enumerate edit scripts without DP merging, bands or seed filtering. This
// independently checks the complete endpoint set and constrained optimum.
struct Oracle<'a> {
    query: &'a [u8],
    quality: &'a [u8],
    target: &'a [u8],
    constraint: Option<Constraint>,
    budget: usize,
    best: BTreeMap<(usize, usize), f64>,
}

impl Oracle<'_> {
    #[allow(clippy::too_many_arguments)]
    fn visit(
        &mut self,
        start: usize,
        i: usize,
        j: usize,
        edits: usize,
        last: u8,
        score: f64,
        matches: &mut Vec<(usize, usize)>,
    ) {
        if i == self.query.len()
            && j > start
            && self.constraint.is_none_or(|c| {
                (0..c.length).all(|d| matches.contains(&(c.query + d, c.target + d)))
            })
        {
            self.best
                .entry((start, j))
                .and_modify(|s| *s = s.max(score))
                .or_insert(score);
        }
        if i < self.query.len() && j < self.target.len() {
            let equal = self.query[i] == self.target[j];
            let next = edits + usize::from(!equal);
            if next <= self.budget {
                if equal {
                    matches.push((i, j));
                }
                self.visit(
                    start,
                    i + 1,
                    j + 1,
                    next,
                    0,
                    score + emission(self.query[i], self.target[j], self.quality[i]),
                    matches,
                );
                if equal {
                    matches.pop();
                }
            }
        }
        if edits < self.budget {
            if i < self.query.len() && last != 2 {
                self.visit(
                    start,
                    i + 1,
                    j,
                    edits + 1,
                    1,
                    score + if last == 1 { GAP_EXTEND } else { GAP_OPEN },
                    matches,
                );
            }
            if j < self.target.len() && last != 1 {
                self.visit(
                    start,
                    i,
                    j + 1,
                    edits + 1,
                    2,
                    score + if last == 2 { GAP_EXTEND } else { GAP_OPEN },
                    matches,
                );
            }
        }
    }
}

fn check(
    query: &[u8],
    quality: &[u8],
    target: &[u8],
    edits: usize,
    constraint: Option<Constraint>,
) {
    let mut oracle = Oracle {
        query,
        quality,
        target,
        budget: edits,
        constraint,
        best: BTreeMap::new(),
    };
    for start in 0..target.len() {
        oracle.visit(start, 0, start, 0, 0, 0.0, &mut Vec::new());
    }
    let actual = fit_cached(
        query,
        quality,
        target,
        constraint,
        edits,
        MAX_CELLS,
        &mut Subproblems::default(),
    );
    assert!(actual.complete && !actual.exhausted);
    let found: BTreeMap<_, _> = actual
        .hits
        .iter()
        .map(|h| ((h.start, h.end), h.logp))
        .collect();
    assert_eq!(
        found.keys().collect::<Vec<_>>(),
        oracle.best.keys().collect::<Vec<_>>(),
        "query={query:?} target={target:?} constraint={constraint:?} E={edits}"
    );
    for (endpoint, expected) in oracle.best {
        assert!(
            (found[&endpoint] - expected).abs() < 1e-8,
            "{endpoint:?}: got {}, expected {expected}",
            found[&endpoint]
        );
    }
    for hit in &actual.hits {
        let mut i = 0;
        let mut j = hit.start;
        let mut score = 0.0;
        let mut last = 0;
        let mut edits = 0;
        let mut matching = Vec::new();
        for &step in &hit.traceback.as_ref().unwrap().steps {
            match step {
                0 => {
                    score += emission(query[i], target[j], quality[i]);
                    edits += usize::from(query[i] != target[j]);
                    if query[i] == target[j] {
                        matching.push((i, j));
                    }
                    i += 1;
                    j += 1;
                }
                1 => {
                    score += if last == 1 { GAP_EXTEND } else { GAP_OPEN };
                    i += 1;
                    edits += 1;
                }
                _ => {
                    score += if last == 2 { GAP_EXTEND } else { GAP_OPEN };
                    j += 1;
                    edits += 1;
                }
            }
            last = step;
        }
        assert_eq!((i, j, edits), (query.len(), hit.end, hit.edits));
        assert!((hit.logp - score).abs() < 1e-8);
        if let Some(c) = constraint {
            assert!((0..c.length).all(|d| matching.contains(&(c.query + d, c.target + d))));
        }
    }
}

#[test]
fn endpoint_edit_and_pivot_profiles_match_exhaustive_edit_scripts() {
    for q in 0..8 {
        let query: Vec<_> = (0..3).map(|i| b"AC"[(q >> i) & 1]).collect();
        for t in 0..16 {
            let target: Vec<_> = (0..4).map(|i| b"AC"[(t >> i) & 1]).collect();
            for edits in 0..=2 {
                check(&query, b"I!I", &target, edits, None);
                check(
                    &query,
                    b"I!I",
                    &target,
                    edits,
                    Some(Constraint {
                        query: 1,
                        target: 1,
                        length: 1,
                    }),
                );
            }
        }
    }
    check(b"ACNAC", b"I!III", b"AACNACG", 3, None);
}

#[test]
fn lower_sequence_score_under_the_pivot_is_retained() {
    let constraint = Some(Constraint {
        query: 0,
        target: 0,
        length: 1,
    });
    check(b"ACAC", b"IIII", b"AACAC", 1, constraint);
    let anchored = fit(b"ACAC", b"IIII", b"AACAC", constraint, 1, MAX_CELLS);
    assert!(anchored
        .hits
        .iter()
        .any(|h| h.start == 0 && h.end == 5 && h.edits == 1));
    let free = fit(b"ACAC", b"IIII", b"AACAC", None, 1, MAX_CELLS);
    assert!(free
        .hits
        .iter()
        .any(|h| h.start == 1 && h.end == 5 && h.edits == 0));
}

#[test]
fn distinct_endpoints_survive_and_tied_traces_cannot_certify_adjacency() {
    let found = fit(b"AAA", b"III", b"AAAA", None, 1, MAX_CELLS);
    assert!(found.hits.iter().filter(|h| h.start == 0).count() > 1);
    let hit = found
        .hits
        .iter()
        .find(|h| h.start == 0 && h.end == 4)
        .unwrap();
    assert!(hit.traceback.as_ref().unwrap().ambiguous);
    let read = crate::io::fastq::Record {
        header: "@oracle/1".into(),
        sequence: b"AAA".to_vec(),
        quality: b"III".to_vec(),
        plus: "+".into(),
    };
    assert!(!super::super::trace_alignment(&read, hit, b"AAAA", 0).matches_contiguous(0, 2));
}

#[test]
fn incomplete_work_is_never_reported_as_proven_absence() {
    let found = fit(b"ACAC", b"IIII", b"AACAC", None, 2, 1);
    assert!(found.exhausted && !found.complete);
    let absent = fit(b"AAA", b"III", b"CCC", None, 0, MAX_CELLS);
    assert!(absent.complete && !absent.exhausted && absent.hits.is_empty());
}

#[test]
fn dominated_states_match_dense_dp_with_realistic_lengths_and_limits() {
    let mut seed = 73u64;
    for length in [12usize, 31, 70, 150] {
        let query: Vec<_> = (0..length)
            .map(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(seed >> 32) as usize % 4]
            })
            .collect();
        let quality: Vec<_> = (0..length).map(|i| b"F:,I"[i % 4]).collect();
        let mut target = [
            b"TTACGATTCGGA".as_slice(),
            query.as_slice(),
            b"GACTGACCTAGC",
        ]
        .concat();
        target.insert(12 + length / 2, b'A');
        let edits = (length / 10).max(3);
        for constraint in [
            None,
            Some(Constraint {
                query: 1,
                target: 13,
                length: 4,
            }),
        ] {
            for budget in [100_000, 300_000, MAX_CELLS] {
                same_profiles(
                    &fit(&query, &quality, &target, constraint, edits, budget),
                    &fit_cached(
                        &query,
                        &quality,
                        &target,
                        constraint,
                        edits,
                        budget,
                        &mut Subproblems::default(),
                    ),
                );
            }
        }
    }
}

#[test]
#[ignore = "manual release benchmark of exact edit-state dominance"]
fn benchmark_edit_state_dominance() {
    let mut seed = 73u64;
    let target: Vec<_> = (0..350)
        .map(|_| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(seed >> 32) as usize % 4]
        })
        .collect();
    let query = &target[90..240];
    let quality = vec![b'F'; query.len()];
    for constraint in [
        None,
        Some(Constraint {
            query: 65,
            target: 155,
            length: 21,
        }),
    ] {
        same_profiles(
            &fit(query, &quality, &target, constraint, 15, MAX_CELLS),
            &fit_cached(
                query,
                &quality,
                &target,
                constraint,
                15,
                MAX_CELLS,
                &mut Subproblems::default(),
            ),
        );
        for fast in [false, true, false, true] {
            let started = std::time::Instant::now();
            for _ in 0..100 {
                let result = if fast {
                    fit_cached(
                        query,
                        &quality,
                        &target,
                        constraint,
                        15,
                        MAX_CELLS,
                        &mut Subproblems::default(),
                    )
                } else {
                    fit(query, &quality, &target, constraint, 15, MAX_CELLS)
                };
                std::hint::black_box(result);
            }
            eprintln!(
                "constraint={constraint:?}, fast={fast}: {:?}",
                started.elapsed()
            );
        }
    }
}

#[test]
fn dominated_profiles_preserve_mixed_qualities_repeats_and_tied_traces() {
    let mut state = 71419u64;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        state >> 32
    };
    for case in 0..240 {
        let length = 4 + next() as usize % 65;
        let offset = 5;
        let target: Vec<_> = (0..length + 20)
            .map(|i| match case % 4 {
                0 => b"ACGTN"[next() as usize % 5],
                1 => b"ACGTAC"[i % 6],
                2 => b'A',
                _ => b"ACGT"[next() as usize % 4],
            })
            .collect();
        let mut query = target[offset..offset + length].to_vec();
        for _ in 0..case % 5 {
            let p = 2 + next() as usize % (query.len() - 2);
            query[p] = b"ACGTN"[next() as usize % 5];
        }
        if case % 7 == 0 {
            query.insert(query.len() / 2, b'C');
        }
        if case % 11 == 0 {
            query.remove(query.len() - 1);
        }
        let quality: Vec<_> = (0..query.len())
            .map(|_| b"!$+5:FI~"[next() as usize % 8])
            .collect();
        let edits = case % 9;
        let constraint = (case % 3 != 0).then_some(Constraint {
            query: 0,
            target: offset,
            length: 2,
        });
        let budget = if case % 6 == 0 { 20000 } else { MAX_CELLS };
        let expected = fit(&query, &quality, &target, constraint, edits, budget);
        let mut cache = Subproblems::default();
        for _ in 0..2 {
            same_profiles(
                &expected,
                &fit_cached(
                    &query, &quality, &target, constraint, edits, budget, &mut cache,
                ),
            );
        }
    }
    // Domains too wide for a u64 mask retain the dense fallback.
    same_profiles(
        &fit(b"AGCTA", b"!I+I~", b"TAGCTAAGCTAG", None, 64, MAX_CELLS),
        &fit_cached(
            b"AGCTA",
            b"!I+I~",
            b"TAGCTAAGCTAG",
            None,
            64,
            MAX_CELLS,
            &mut Subproblems::default(),
        ),
    );
}

#[test]
fn shared_suffix_bounds_preserve_every_feasible_banded_state() {
    let mut random = 571u64;
    for case in 0..24 {
        let query: Vec<_> = (0..10 + case)
            .map(|i| {
                random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
                if case % 3 == 0 {
                    b"ACN"[i % 3]
                } else {
                    b"ACGT"[(random >> 32) as usize % 4]
                }
            })
            .collect();
        let mut target = [b"ACGTACGT".as_slice(), &query, b"AGCTAGCTACGT"].concat();
        target.insert(query.len() / 2 + 8, b'N');
        let edits = case % 8;
        let suffix = SuffixBounds::new(&query, &target, edits);
        for start in 0..target.len() {
            let part = &target[start..target.len().min(start + query.len() + edits)];
            let old = edit_bounds(&query, part, edits, None);
            let new = edit_bounds(&query, part, edits, Some((&suffix, start)));
            for (a, b) in old.into_iter().zip(new) {
                assert_eq!(a[0], b[0]);
                assert_eq!(a[0] + a[1] <= edits, b[0] + b[1] <= edits);
                if a[0] + a[1] <= edits {
                    assert_eq!(a, b);
                }
            }
        }
    }
}
