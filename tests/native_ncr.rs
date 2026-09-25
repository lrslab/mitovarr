use mitovarr::{
    annotation::{ncr, predict, rotate_to_pro, Feature, Panel, Prediction, Segment},
    cli::Species,
    io::fasta::load_reference,
};
use std::path::PathBuf;
fn gene(name: &str, start: usize, end: usize, strand: char) -> Feature {
    Feature {
        id: format!("{name}_{start}"),
        name: name.into(),
        kind: if name.starts_with("tRNA") {
            "tRNA"
        } else {
            "CDS"
        }
        .into(),
        product: String::new(),
        strand,
        segments: vec![Segment { start, end }],
        length: end - start + 1,
        input_five_prime_1based: if strand == '+' { start } else { end },
        five_prime_1based: if strand == '+' { start } else { end },
        model_length: end - start + 1,
        query_start: 1,
        query_end: end - start + 1,
        query_coverage: 1.0,
        identity: 1.0,
        alignment_score: 100,
        gap_bases: 0,
        five_prime_boundary: "aligned".into(),
        three_prime_boundary: "aligned".into(),
        status: "complete".into(),
        full_length_supported: true,
        ambiguous_bases: 0,
        competing_models: Vec::new(),
        rna_structure: None,
        cds: None,
        sequence: Vec::new(),
    }
}
fn prediction(features: Vec<Feature>) -> Prediction {
    Prediction {
        features,
        missing_models: Vec::new(),
        models_without_full_length_support: Vec::new(),
        search_limited_models: Vec::new(),
    }
}

#[test]
fn three_species_reference_ncr_catalog_and_rotation_invariance() {
    for (sp, name, expected) in [
        (Species::Cel, "cel", [Some(466), None, Some(109), None]),
        (
            Species::Cbr,
            "cbr",
            [Some(548), Some(222), Some(109), Some(324)],
        ),
        (
            Species::Ctro,
            "ctro",
            [Some(410), Some(130), Some(107), None],
        ),
    ] {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(format!("tests/fixtures/annotation/{name}.fasta"));
        let sequence = load_reference(&path, None).unwrap().sequence;
        let panel = Panel::bundled(sp).unwrap();
        let mut p = predict(&sequence, &panel, true);
        let before = ncr::derive(&sequence, &p, true, Some(sp), 20);
        let (rotated, rotation) = rotate_to_pro(&sequence, &mut p, true, true);
        assert!(rotation.applied);
        let after = ncr::derive(&rotated, &p, true, Some(sp), 20);
        for (i, length) in expected.iter().enumerate() {
            let locus = format!("NCR{}", i + 1);
            let a: Vec<_> = before.regions.iter().filter(|r| r.name == locus).collect();
            let b: Vec<_> = after.regions.iter().filter(|r| r.name == locus).collect();
            assert_eq!(a.len(), usize::from(length.is_some()), "{name} {locus}");
            assert_eq!(a.len(), b.len());
            if let Some(length) = length {
                assert_eq!(a[0].length, *length);
                assert_eq!(b[0].length, *length);
                assert_eq!(b[0].boundary_status, "supported");
            }
        }
        assert_eq!(
            after
                .regions
                .iter()
                .find(|r| r.name == "NCR1")
                .unwrap()
                .segments[0]
                .end,
            sequence.len()
        );
        if name == "cel" {
            assert_eq!(after.loci[1].observed_gap_length, Some(7));
            assert_eq!(after.loci[1].status, "below_min_length");
            assert_eq!(after.loci[3].overlap_bases, 1);
        }
        if name == "ctro" {
            assert_eq!(after.loci[1].flank_genes, ["tRNA-Arg", "tRNA-Gln"]);
            assert_eq!(after.loci[3].observed_gap_length, Some(1));
        }
    }
}

#[test]
fn circular_head_tail_are_one_ncr_before_length_filter() {
    let p = prediction(vec![
        gene("tRNA-Pro", 13, 30, '+'),
        gene("COX2", 31, 69, '+'),
        gene("tRNA-Ala", 70, 88, '+'),
    ]);
    let seq = vec![b'A'; 100];
    let a = ncr::derive(&seq, &p, true, Some(Species::Cel), 20);
    let r = a.regions.iter().find(|r| r.name == "NCR1").unwrap();
    assert_eq!(r.length, 24);
    assert_eq!(r.segments.len(), 2);
    assert_eq!(r.segments[0].start, 89);
    assert_eq!(r.segments[1].end, 12);
    let linear = ncr::derive(&seq, &p, false, Some(Species::Cel), 20);
    assert!(linear.regions.is_empty());
    assert_eq!(linear.loci[0].status, "unresolved_contig_ends");
    assert_eq!(linear.loci[0].observed_gap_length, None);
}

#[test]
fn linear_terminal_ncr_fragments_are_not_joined() {
    let p = prediction(vec![
        gene("tRNA-Pro", 22, 50, '+'),
        gene("COX2", 51, 99, '+'),
        gene("tRNA-Ala", 100, 129, '+'),
    ]);
    let a = ncr::derive(&[b'A'; 150], &p, false, Some(Species::Cel), 20);
    let pieces: Vec<_> = a
        .regions
        .iter()
        .filter(|r| r.locus.as_deref() == Some("NCR1"))
        .collect();
    assert_eq!(pieces.len(), 2);
    assert!(pieces.iter().all(|r| r.name == "NCR1_candidate"
        && r.terminal_fragment
        && r.length == 21
        && r.segments.len() == 1));
    assert_eq!(a.loci[0].status, "unresolved_contig_ends");
    assert_eq!(a.loci[0].observed_gap_length, None);
}

#[test]
fn no_genes_does_not_make_the_whole_contig_ncr1() {
    let a = ncr::derive(
        &[b'A'; 100],
        &prediction(Vec::new()),
        true,
        Some(Species::Cbr),
        20,
    );
    assert!(a.regions.is_empty());
    assert!(a.loci.iter().all(|l| l.status == "missing_flank_genes"));
}

#[test]
fn truncated_nd4_does_not_become_a_false_ncr3_insertion() {
    let mut panel = Panel::bundled(Species::Cel).unwrap();
    panel.genes.retain(|g| g.name == "ND4" || g.name == "COX1");
    let nd4 = panel
        .genes
        .iter()
        .find(|g| g.name == "ND4")
        .unwrap()
        .sequence
        .as_bytes();
    let cox1 = panel
        .genes
        .iter()
        .find(|g| g.name == "COX1")
        .unwrap()
        .sequence
        .as_bytes();
    let seq = [&nd4[..nd4.len() - 100], &[b'N'; 100], &[b'C'; 109], cox1].concat();
    let p = predict(&seq, &panel, false);
    let a = ncr::derive(&seq, &p, false, Some(Species::Cel), 20);
    let r = a
        .regions
        .iter()
        .find(|r| r.locus.as_deref() == Some("NCR3"))
        .unwrap();
    assert_eq!(r.length, 209);
    assert_eq!(r.name, "NCR3_candidate");
    assert_eq!(r.boundary_status, "uncertain");
    assert_eq!(a.loci[2].status, "boundary_uncertain");
}

#[test]
fn local_boundary_evidence_survives_unrelated_gene_quality_flags() {
    let mut gln = gene("tRNA-Gln", 1, 55, '+');
    gln.status = "boundary_complete_length_variant".into();
    gln.full_length_supported = false;
    gln.gap_bases = 1;
    let p = prediction(vec![gln, gene("tRNA-Phe", 278, 334, '+')]);
    let mut sequence = vec![b'A'; 334];
    sequence[55..277].fill(b'N');
    let a = ncr::derive(&sequence, &p, false, Some(Species::Cbr), 20);
    let r = &a.regions[0];
    assert_eq!(r.name, "NCR2");
    assert_eq!(r.length, 222);
    assert_eq!(r.ambiguous_bases, 222);
    assert_eq!(r.at_fraction, None);
    assert_eq!(a.loci[1].status, "annotated_with_unknown_bases");
}

#[test]
fn short_nd5_homolog_does_not_consume_ncr4() {
    let mut fragment = gene("ND5", 360, 620, '+');
    fragment.query_coverage = 261.0 / 1584.0;
    fragment.model_length = 1584;
    fragment.status = "partial".into();
    fragment.full_length_supported = false;
    let id = fragment.id.clone();
    let p = prediction(vec![
        gene("ND3", 1, 336, '+'),
        gene("ND5", 661, 2244, '+'),
        fragment,
    ]);
    let a = ncr::derive(&[b'A'; 2244], &p, false, Some(Species::Cbr), 20);
    let r = a.regions.iter().find(|r| r.name == "NCR4").unwrap();
    assert_eq!(r.length, 324);
    assert_eq!(r.embedded_partial_homologs, vec![id]);
}

#[test]
fn duplicate_loci_are_retained_and_cel_insert_is_a_positional_candidate() {
    let p = prediction(vec![
        gene("ND4", 1, 50, '+'),
        gene("COX1", 101, 150, '+'),
        gene("ND4", 151, 200, '+'),
        gene("COX1", 251, 300, '+'),
    ]);
    let a = ncr::derive(&[b'A'; 300], &p, false, Some(Species::Cel), 20);
    assert_eq!(a.loci[2].status, "multiple_candidates");
    assert_eq!(a.loci[2].region_ids.len(), 2);
    let p = prediction(vec![
        gene("tRNA-Gln", 1, 55, '+'),
        gene("tRNA-Phe", 106, 162, '+'),
    ]);
    let a = ncr::derive(&[b'A'; 162], &p, false, Some(Species::Cel), 20);
    assert_eq!(a.regions[0].name, "NCR2_candidate");
    assert_eq!(a.loci[1].status, "positional_candidate");
}

#[test]
fn negative_strand_uses_the_boundary_that_faces_the_ncr() {
    let mut cox = gene("COX1", 1, 60, '-');
    cox.three_prime_boundary = "partial_or_uncertain".into();
    cox.full_length_supported = false;
    let p = prediction(vec![cox, gene("ND4", 170, 240, '-')]);
    let a = ncr::derive(&[b'A'; 240], &p, false, Some(Species::Cel), 20);
    assert_eq!(a.regions[0].name, "NCR3");
    let mut p = p;
    p.features[0].five_prime_boundary = "partial_or_uncertain".into();
    let a = ncr::derive(&[b'A'; 240], &p, false, Some(Species::Cel), 20);
    assert_eq!(a.regions[0].name, "NCR3_candidate");
}

#[test]
fn terminal_label_requires_the_correct_gene_end() {
    for strand in ['+', '-'] {
        let p = prediction(vec![gene("ND5", 31, 90, strand)]);
        let a = ncr::derive(&[b'A'; 120], &p, false, Some(Species::Cbr), 20);
        assert_eq!(a.regions.len(), 2);
        let index = usize::from(strand == '-');
        assert_eq!(a.regions[index].name, "NCR4_candidate");
        assert_eq!(a.regions[1 - index].name, "unclassified_gap");
        assert_eq!(a.loci[3].observed_gap_length, None);
    }
}

#[test]
fn opposite_arc_is_not_a_second_copy_of_the_same_ncr() {
    let p = prediction(vec![gene("ND3", 1, 30, '+'), gene("ND5", 61, 90, '+')]);
    let a = ncr::derive(&[b'A'; 120], &p, true, Some(Species::Cbr), 20);
    assert_eq!(a.regions.len(), 2);
    assert_eq!(a.regions[0].name, "NCR4");
    assert_eq!(a.regions[1].name, "unclassified_gap");
    assert_eq!(a.loci[3].status, "annotated");
}

#[test]
fn adjacent_internal_flanks_are_not_open_contig_ends() {
    let p = prediction(vec![gene("ND3", 1, 30, '+'), gene("ND5", 31, 90, '+')]);
    let a = ncr::derive(&[b'A'; 90], &p, false, Some(Species::Cbr), 20);
    assert_eq!(a.loci[3].status, "adjacent_without_gap");
    assert_eq!(a.loci[3].observed_gap_length, Some(0));
    let p = prediction(vec![gene("ND5", 1, 60, '+'), gene("ND3", 61, 90, '+')]);
    let a = ncr::derive(&[b'A'; 90], &p, false, Some(Species::Cbr), 20);
    assert_eq!(a.loci[3].status, "unresolved_contig_ends");
}
