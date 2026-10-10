use super::super::joint_exhaustive;
use super::*;

fn check(
    a: &[(Hit, BaseTrace)],
    b: &[(Hit, BaseTrace)],
    model: Option<&Library>,
    ideal: f64,
    length: usize,
    junctions: [usize; 2],
    contraction: Option<[usize; 2]>,
) {
    let expected = joint_exhaustive(a, b, model, ideal, length, junctions, contraction);
    let actual = joint(a, b, model, ideal, length, junctions, contraction);
    assert_eq!(actual.score.to_bits(), expected.score.to_bits());
    assert_eq!(actual.diagnostics, expected.diagnostics);
    assert_eq!(actual.junctions, expected.junctions);
    assert_eq!(actual.observed_junctions, expected.observed_junctions);
    assert_eq!(actual.mate_brackets, expected.mate_brackets);
    assert_eq!(actual.contraction, expected.contraction);
}

fn placement(
    start: usize,
    end: usize,
    reverse: bool,
    logp: f64,
    ambiguous: usize,
) -> (Hit, BaseTrace) {
    (
        Hit {
            start,
            end,
            reverse,
            logp,
            edits: 0,
            traceback: None,
        },
        BaseTrace {
            matching: (start..end).collect(),
            matching_junctions: (start + 1..end).collect(),
            ambiguous_contacts: ambiguous,
            ..BaseTrace::default()
        },
    )
}

#[test]
fn indexed_joint_matches_oracle_for_overlaps_containment_origins_and_evidence() {
    let mut seed = 19u64;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        seed >> 32
    };
    let probabilities: Vec<_> = (0..501)
        .map(|i| {
            if i % 13 == 0 {
                0.0
            } else {
                0.0001 * (1 + i % 29) as f64
            }
        })
        .collect();
    let models: Vec<_> = [true, false]
        .into_iter()
        .flat_map(|inward| {
            [true, false].map(|usable| Library::test_model(inward, usable, probabilities.clone()))
        })
        .collect();
    for case in 0..240 {
        let length = if case % 3 == 0 { 0 } else { 240 };
        let junctions = [next() as usize % 241, next() as usize % 241];
        let contraction = (case % 4 != 0).then(|| [next() as usize % 220, next() as usize % 220]);
        let mut make = |count| {
            (0..count)
                .map(|_| {
                    let start = next() as usize % 240;
                    let end = start + 1 + next() as usize % 160;
                    let score = [0.0, -0.01, -0.1, -0.5, -2.0, 0.1f64.ln(), -6.0, -14.0]
                        [next() as usize % 8];
                    let mut p = placement(start, end, next() % 2 != 0, score, next() as usize % 9);
                    match next() % 5 {
                        0 => p.1.ambiguous_trace = true,
                        1 => p.1.matching.retain(|p| p % 31 != 0),
                        2 => p.1.matching_junctions.retain(|p| p % 37 != 0),
                        _ => (),
                    }
                    p
                })
                .collect::<Vec<_>>()
        };
        let mut a = make(case % 19);
        let mut b = make(case % 23);
        for model in std::iter::once(None).chain(models.iter().map(Some)) {
            check(&a, &b, model, -0.5, length, junctions, contraction);
        }
        a.reverse();
        b.reverse();
        check(
            &a,
            &b,
            Some(&models[0]),
            -0.5,
            length,
            junctions,
            contraction,
        );
    }
}

#[test]
fn range_max_preserves_ambiguous_contacts_in_rounded_and_clamped_ties() {
    let a = [placement(0, 100, false, -1e16, 2)];
    let b = [
        placement(120, 200, true, 0.0, 0),
        placement(125, 200, true, -0.1, 9),
    ];
    check(&a, &b, None, -1e16, 0, [90, 110], None);
    let a = [placement(0, 100, false, 0.0, 2)];
    check(&a, &b, None, -1.0, 0, [90, 110], None);
    let tree = RangeMax::new(vec![
        (
            120,
            Entry {
                start: 120,
                end: 200,
                score: 0.0,
                ambiguous: 0,
            },
        ),
        (
            125,
            Entry {
                start: 125,
                end: 200,
                score: -0.1,
                ambiguous: 9,
            },
        ),
    ]);
    let best = tree.query(Range { lo: 0, hi: 200 }, -1e16, -1e16).unwrap();
    assert_eq!(best.score, 0.0);
    assert_eq!(best.ambiguous, 9);
}

#[test]
fn suboptimal_sequence_explanation_can_supply_junction_and_contraction_support() {
    let mut a = vec![
        placement(20, 110, false, 0.0, 0),
        placement(20, 110, false, -0.005, 0),
    ];
    a[0].1.ambiguous_trace = true;
    let b = [placement(125, 210, true, 0.0, 0)];
    let model = Library::test_model(true, true, vec![0.01; 401]);
    check(&a, &b, Some(&model), 0.0, 240, [70, 120], Some([30, 160]));
    let result = joint(&a, &b, Some(&model), 0.0, 240, [70, 120], Some([30, 160]));
    assert!(result.observed_junctions[0]);
    assert!(result.mate_brackets[1]);
    assert!(result.contraction);
}

#[test]
fn endpoint_grid_matches_every_span_and_support_flag() {
    let a: Vec<_> = (0..16)
        .flat_map(|i| {
            (0..16).map(move |j| {
                placement(20 + i, 155 + j, false, -((i + j) as f64) * 0.2, (i + j) % 3)
            })
        })
        .collect();
    let b: Vec<_> = (0..16)
        .flat_map(|i| {
            (0..16).map(move |j| {
                placement(190 + i, 325 + j, true, -((i + j) as f64) * 0.2, (i + j) % 5)
            })
        })
        .collect();
    let model = Library::test_model(
        true,
        true,
        (0..501)
            .map(|i| if i > 250 && i < 400 { 0.002 } else { 0.0 })
            .collect(),
    );
    check(&a, &b, Some(&model), -0.1, 600, [140, 185], Some([40, 260]));
}

#[test]
#[ignore = "manual release benchmark; compares the same full endpoint grid"]
fn benchmark_dense_endpoint_grid() {
    let a: Vec<_> = (0..31)
        .flat_map(|i| {
            (0..16).map(move |j| placement(20 + i, 155 + j, false, -((i + j) as f64) * 0.2, 0))
        })
        .collect();
    let b: Vec<_> = (0..31)
        .flat_map(|i| {
            (0..16).map(move |j| placement(190 + i, 325 + j, true, -((i + j) as f64) * 0.2, 0))
        })
        .collect();
    let model = Library::test_model(true, true, vec![0.002; 601]);
    check(&a, &b, Some(&model), 0.0, 600, [140, 185], Some([40, 260]));
    for (name, implementation) in [
        (
            "enumerated",
            joint_exhaustive as fn(_, _, _, _, _, _, _) -> _,
        ),
        ("indexed", joint),
    ] {
        let started = std::time::Instant::now();
        for _ in 0..100 {
            std::hint::black_box(implementation(
                &a,
                &b,
                Some(&model),
                0.0,
                600,
                [140, 185],
                Some([40, 260]),
            ));
        }
        eprintln!("{name}: {:?}", started.elapsed());
    }
}
