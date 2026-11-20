use mitovarr::{
    annotation::{predict, rotate_to_pro, Panel},
    cli::Species,
    io::fasta::load_reference,
    kmer::reverse_complement,
};
use serde_json::{json, Value};
use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
fn panel(names: &[&str]) -> Panel {
    let mut p = Panel::bundled(Species::Cel).unwrap();
    p.genes.retain(|g| names.contains(&g.name.as_str()));
    p
}
fn dna(p: &Panel, name: &str) -> Vec<u8> {
    p.genes
        .iter()
        .find(|g| g.name == name)
        .unwrap()
        .sequence
        .as_bytes()
        .to_vec()
}
fn reference(species: &str) -> Vec<u8> {
    load_reference(
        &PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(format!("tests/fixtures/annotation/{species}.fasta")),
        None,
    )
    .unwrap()
    .sequence
}

#[test]
fn all_three_pinned_references_recover_full_genes_boundaries_and_proteins() {
    for (name, species) in [
        ("cel", Species::Cel),
        ("cbr", Species::Cbr),
        ("ctro", Species::Ctro),
    ] {
        let p = Panel::bundled(species).unwrap();
        let sequence = reference(name);
        let prediction = predict(&sequence, &p, true);
        let raw: Value = serde_json::from_str(
            &fs::read_to_string(
                PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                    .join(format!("src/annotation/models/{name}.json")),
            )
            .unwrap(),
        )
        .unwrap();
        assert!(
            prediction.search_limited_models.is_empty(),
            "{name}: {:?}",
            prediction.search_limited_models
        );
        for g in raw["genes"].as_array().unwrap() {
            let hits: Vec<_> = prediction
                .features
                .iter()
                .filter(|f| f.name == g["name"].as_str().unwrap() && f.full_length_supported)
                .collect();
            assert_eq!(
                hits.len(),
                1,
                "{name} {}: {:?}",
                g["name"],
                prediction
                    .features
                    .iter()
                    .filter(|f| f.name == g["name"].as_str().unwrap())
                    .collect::<Vec<_>>()
            );
            let f = hits[0];
            assert_eq!(f.segments.len(), 1);
            assert_eq!(
                f.segments[0].start,
                g["source_start"].as_u64().unwrap() as usize
            );
            assert_eq!(
                f.segments[0].end,
                g["source_end"].as_u64().unwrap() as usize
            );
            assert_eq!(f.strand.to_string(), g["source_strand"].as_str().unwrap());
            if let Some(cds) = &f.cds {
                assert_eq!(
                    cds.protein,
                    g["translation"].as_str().unwrap(),
                    "{name} {}",
                    f.name
                );
            } else {
                assert_eq!(f.rna_structure.as_deref(), Some("not_assessed"));
            }
        }
        assert!(
            prediction.models_without_full_length_support.is_empty(),
            "{name}: {:?}",
            prediction.models_without_full_length_support
        );
    }
}

#[test]
fn partial_cds_and_rnas_are_not_counted_as_complete() {
    for (name, removed) in [("ND4", 100), ("rrnS", 104), ("tRNA-Ala", 8)] {
        let p = panel(&[name]);
        let mut q = dna(&p, name);
        q.truncate(q.len() - removed);
        let seq = [vec![b'N'; 30], q, vec![b'N'; 30]].concat();
        let prediction = predict(&seq, &p, false);
        let f = prediction.features.iter().find(|f| f.name == name).unwrap();
        assert!(!f.full_length_supported, "{name}: {f:?}");
        assert_eq!(f.three_prime_boundary, "partial_or_uncertain");
        assert_eq!(prediction.models_without_full_length_support, vec![name]);
    }
}

#[test]
fn rna_terminal_substitutions_are_full_length_but_missing_ends_are_not() {
    // CaenDR C. tropicalis strains carry A where NC_025756.1 tRNA-Trp starts with G.
    let mut p = Panel::bundled(Species::Ctro).unwrap();
    p.genes.retain(|g| g.name == "tRNA-Trp");
    let q = dna(&p, "tRNA-Trp");
    let flank = b"CCCCCCCCCCCCCCCCCCCC".to_vec();
    let call = |gene: Vec<u8>, left: &[u8], right: &[u8]| {
        let seq = [left.to_vec(), gene, right.to_vec()].concat();
        let prediction = predict(&seq, &p, false);
        prediction
            .features
            .into_iter()
            .find(|f| f.name == "tRNA-Trp")
            .unwrap()
    };
    let swap = |b: u8| if b == b'A' { b'G' } else { b'A' };

    let mut first = q.clone();
    first[0] = swap(first[0]);
    let f = call(first, &flank, &flank);
    assert!(f.full_length_supported, "{f:?}");
    assert_eq!(f.five_prime_boundary, "terminal_mismatch");
    assert_eq!(f.three_prime_boundary, "aligned");
    assert_eq!((f.segments[0].start, f.length), (flank.len() + 1, q.len()));

    let mut last = q.clone();
    *last.last_mut().unwrap() = swap(*last.last().unwrap());
    let f = call(last, &flank, &flank);
    assert!(f.full_length_supported, "{f:?}");
    assert_eq!(f.three_prime_boundary, "terminal_mismatch");
    assert_eq!(f.length, q.len());

    // An ambiguous or absent terminal base never completes an RNA.
    let mut ambiguous = q.clone();
    ambiguous[0] = b'N';
    let f = call(ambiguous, &flank, &flank);
    assert!(!f.full_length_supported);
    assert_eq!(f.five_prime_boundary, "partial_or_uncertain");
    let f = call(q[1..].to_vec(), b"NNNNNNNNNN", b"NNNNNNNNNN");
    assert!(!f.full_length_supported);
    let f = call(q[..q.len() - 1].to_vec(), b"NNNNNNNNNN", b"NNNNNNNNNN");
    assert!(!f.full_length_supported);
    assert_eq!(f.three_prime_boundary, "partial_or_uncertain");

    // Three clipped bases exceed the substitution allowance.
    let mut three = q.clone();
    for b in &mut three[..3] {
        *b = swap(*b);
    }
    let f = call(three, &flank, &flank);
    assert!(!f.full_length_supported);
}

#[test]
fn clipped_cds_tail_is_not_a_polyadenylated_stop() {
    let p = panel(&["ND4"]);
    let q = dna(&p, "ND4");
    for tail in [b"T".as_slice(), b"TA".as_slice()] {
        let end = (q.len() / 2..q.len() - 30)
            .find(|&end| end % 3 == tail.len() && q[..end].ends_with(tail))
            .unwrap();
        let prediction = predict(&q[..end], &p, false);
        let f = prediction
            .features
            .iter()
            .find(|f| f.name == "ND4")
            .unwrap();
        assert_eq!(f.three_prime_boundary, "partial_or_uncertain");
        assert!(!f.full_length_supported);
        let cds = f.cds.as_ref().unwrap();
        assert_eq!(cds.stop, "partial_or_uncertain");
        // An incomplete terminal codon alone is not evidence of a frameshift.
        assert!(cds.frame_valid);
    }
}

#[test]
fn cds_indels_ambiguities_and_internal_stops_are_reported() {
    let p = panel(&["ND4"]);
    let q = dna(&p, "ND4");
    let mut deletion = q.clone();
    deletion.remove(500);
    let a = predict(&deletion, &p, false);
    let f = &a.features[0];
    assert!(!f.full_length_supported);
    assert!(f.cds.as_ref().unwrap().non_triplet_alignment_gaps > 0);
    let mut uncertain = q.clone();
    uncertain[498..501].copy_from_slice(b"NNN");
    let a = predict(&uncertain, &p, false);
    assert!(!a.features[0].full_length_supported);
    assert_eq!(a.features[0].ambiguous_bases, 3);
    let mut stopped = q.clone();
    stopped[498..501].copy_from_slice(b"TAA");
    let a = predict(&stopped, &p, false);
    assert!(!a.features[0].full_length_supported);
    assert!(!a.features[0]
        .cds
        .as_ref()
        .unwrap()
        .internal_stop_codons
        .is_empty());
    let mut in_frame = q;
    in_frame.drain(498..501);
    let a = predict(&in_frame, &p, false);
    assert!(a.features[0].full_length_supported, "{:?}", a.features[0]);
}

#[test]
fn retained_duplicate_and_inverted_genes_prevent_ambiguous_rotation() {
    let p = panel(&["tRNA-Pro"]);
    let q = dna(&p, "tRNA-Pro");
    for second in [q.clone(), reverse_complement(&q)] {
        let seq = [q.clone(), vec![b'N'; 40], second, vec![b'N'; 40]].concat();
        let mut a = predict(&seq, &p, true);
        assert_eq!(a.features.len(), 2, "{:?}", a.features);
        let (out, r) = rotate_to_pro(&seq, &mut a, true, true);
        assert!(!r.applied);
        assert_eq!(r.reason, "multiple_trna_pro_loci");
        assert_eq!(out, seq);
    }
}

#[test]
fn truncated_trna_pro_cannot_set_a_false_five_prime_origin() {
    let p = panel(&["tRNA-Pro"]);
    let q = dna(&p, "tRNA-Pro");
    let seq = [vec![b'N'; 100], q[5..].to_vec(), vec![b'N'; 80]].concat();
    let mut a = predict(&seq, &p, true);
    assert_eq!(a.features.len(), 1);
    assert_eq!(a.features[0].query_start, 6);
    let (out, r) = rotate_to_pro(&seq, &mut a, true, true);
    assert!(!r.applied);
    assert_eq!(r.reason, "trna_pro_boundary_not_supported");
    assert_eq!(out, seq);
}

#[test]
fn circular_reverse_rotation_preserves_bases_and_all_feature_sequences() {
    let p = panel(&["tRNA-Pro", "ND4"]);
    let pro = dna(&p, "tRNA-Pro");
    let expected = [pro.clone(), vec![b'N'; 37], dna(&p, "ND4"), vec![b'N'; 50]].concat();
    for reverse in [false, true] {
        let mut seq = expected.clone();
        seq.rotate_left(23);
        if reverse {
            seq = reverse_complement(&seq);
        }
        let mut a = predict(&seq, &p, true);
        let (out, r) = rotate_to_pro(&seq, &mut a, true, true);
        assert!(r.applied, "{r:?}");
        assert_eq!(out, expected);
        assert_eq!(r.reverse_complemented, reverse);
        let anchor = a.features.iter().find(|f| f.name == "tRNA-Pro").unwrap();
        assert_eq!(anchor.five_prime_1based, 1);
        assert_eq!(anchor.strand, '+');
        for f in &a.features {
            let mut extracted = Vec::new();
            for s in &f.segments {
                let piece = &out[s.start - 1..s.end];
                if f.strand == '-' {
                    extracted.extend(reverse_complement(piece));
                } else {
                    extracted.extend(piece);
                }
            }
            assert_eq!(extracted, f.sequence, "{}", f.name);
        }
    }
    let mut a = predict(&expected, &p, false);
    let (out, r) = rotate_to_pro(&expected, &mut a, false, true);
    assert_eq!(r.reason, "linear_input");
    assert_eq!(out, expected);
}

static NEXT: AtomicU64 = AtomicU64::new(0);
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "mitovar-annotation-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&p).unwrap();
        Self(p)
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn quality_run(temp: &Temp, sequence: &[u8], evidence: Option<&str>) -> Value {
    let input = temp.0.join("partial.fasta");
    fs::write(
        &input,
        format!(
            ">sample circular=false status=unresolved\n{}\n",
            String::from_utf8_lossy(sequence)
        ),
    )
    .unwrap();
    let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_mitovarr"));
    command
        .env("PATH", "")
        .args(["annotate", "--species", "cel", "--fasta"])
        .arg(&input)
        .arg("--out")
        .arg(temp.0.join("annotation"));
    if let Some(rows) = evidence {
        let path = temp.0.join("ambiguity.tsv");
        fs::write(
            &path,
            format!("sequence\tstart\tend\tevidence\treason\talternatives\n{rows}"),
        )
        .unwrap();
        command.arg("--assembly-ambiguity").arg(path);
    }
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        load_reference(&temp.0.join("annotation/annotated.fasta"), None)
            .unwrap()
            .sequence,
        sequence
    );
    assert!(!temp.0.join("annotation/consensus.fasta").exists());
    serde_json::from_slice(&fs::read(temp.0.join("annotation/report.json")).unwrap()).unwrap()
}

#[test]
fn assembly_length_evidence_labels_ncr_loci_without_changing_sequence() {
    // Mapping maximum consensus reports filled N runs and unmeasured tandem-repeat
    // lengths inside a deferred window; a fill takes precedence.
    let window = "sample\t13300\t13794\tlocal_locus_1\tstructural_validation_deferred\t1\n";
    let fill = "sample\t13400\t13550\tlocal_locus_1\tfixed_span_graph_fill\t0\n";
    let tandem = "sample\t13420\t13600\ttandem_repeat_1\ttandem_repeat_length_unmeasured\t0\n";
    for (rows, origin, status) in [
        (
            format!("{window}{tandem}"),
            "not_assessed",
            "tandem_repeat_length_unmeasured",
        ),
        (
            format!("{window}{fill}{tandem}"),
            "fixed_span_graph_fill",
            "fixed_span_graph_fill",
        ),
        (
            format!("{window}{fill}"),
            "fixed_span_graph_fill",
            "fixed_span_graph_fill",
        ),
        (
            window.to_owned(),
            "not_assessed",
            "retained_assembly_alternatives",
        ),
        (
            format!(
                "{window}sample\t13420\t13600\tread_validation_1\tlow_matching_read_depth\t0\n"
            ),
            "not_assessed",
            "low_read_depth",
        ),
    ] {
        let temp = Temp::new();
        let sequence = reference("cel");
        let input = temp.0.join("consensus.fasta");
        fs::write(
            &input,
            format!(">sample\n{}\n", String::from_utf8_lossy(&sequence)),
        )
        .unwrap();
        let evidence = temp.0.join("ambiguity.tsv");
        fs::write(
            &evidence,
            format!("sequence\tstart\tend\tevidence\treason\talternatives\n{rows}"),
        )
        .unwrap();
        let out = temp.0.join("annotation");
        let result = std::process::Command::new(env!("CARGO_BIN_EXE_mitovarr"))
            .args(["annotate", "--species", "cel", "--topology", "circular"])
            .arg("--fasta")
            .arg(&input)
            .arg("--assembly-ambiguity")
            .arg(&evidence)
            .arg("--out")
            .arg(&out)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            load_reference(&out.join("annotated.fasta"), None)
                .unwrap()
                .sequence,
            sequence
        );
        let table = fs::read_to_string(out.join("ncr.tsv")).unwrap();
        let origins: Vec<(String, String)> = table
            .lines()
            .skip(1)
            .map(|l| {
                let c: Vec<_> = l.split('\t').collect();
                (c[3].to_owned(), c[12].to_owned())
            })
            .collect();
        assert!(
            origins.contains(&("NCR1".into(), origin.into())),
            "{origins:?}"
        );
        assert!(origins
            .iter()
            .filter(|(locus, _)| locus != "NCR1")
            .all(|(_, o)| o == "not_assessed"));
        let gff = fs::read_to_string(out.join("ncr.gff3")).unwrap();
        assert!(
            gff.lines()
                .any(|l| l.contains("Name=NCR1;")
                    && l.ends_with(&format!("sequence_origin={origin}")))
        );
        let report: Value =
            serde_json::from_slice(&fs::read(out.join("report.json")).unwrap()).unwrap();
        let locus = report["records"][0]["assembly_evidence"]["ncr_loci"]
            .as_array()
            .unwrap()
            .iter()
            .find(|l| l["locus"] == "NCR1")
            .unwrap()
            .clone();
        assert_eq!(locus["evidence_status"], status);
    }
}

#[test]
fn partial_genome_without_ncr1_retains_independent_gene_quality() {
    let temp = Temp::new();
    // The pinned cel fixture's last gene, tRNA-Ala, ends at 13328;
    // remove the entire terminal NCR1 and keep the molecule explicitly linear.
    let sequence = reference("cel")[..13328].to_vec();
    let report = quality_run(&temp, &sequence, None);
    let r = &report["records"][0];
    assert_eq!(r["circular_input_assumed"], false);
    assert_eq!(r["ncr"]["loci"][0]["status"], "unresolved_contig_ends");
    let q = &r["gene_quality"];
    assert_eq!(q["status"], "no_detected_gene_issues");
    assert_eq!(q["expected_models"], 36);
    assert_eq!(q["annotation_clear_models"].as_array().unwrap().len(), 36);
    assert_eq!(q["assembly_ambiguity_supplied"], false);
    assert_eq!(q["read_validation"], "not_performed");
    let table = fs::read_to_string(temp.0.join("annotation/gene_quality.tsv")).unwrap();
    assert_eq!(table.lines().count(), 37);
    assert!(table
        .lines()
        .skip(1)
        .all(|r| r.split('\t').nth(3) == Some("annotation_clear")));
}

#[test]
fn gene_quality_retains_coding_missing_copy_and_assembly_problems() {
    let temp = Temp::new();
    let p = Panel::bundled(Species::Cel).unwrap();
    let mut sequence = reference("cel");
    let prediction = predict(&sequence, &p, false);
    let start = |name: &str| {
        prediction
            .features
            .iter()
            .find(|f| f.name == name)
            .unwrap()
            .segments[0]
            .start
            - 1
    };
    let nd4 = start("ND4");
    sequence[nd4 + 90..nd4 + 93].copy_from_slice(b"TAA");
    let rrns = start("rrnS");
    sequence[rrns + 100..rrns + 112].fill(b'N');
    let ala = prediction
        .features
        .iter()
        .find(|f| f.name == "tRNA-Ala")
        .unwrap();
    for s in &ala.segments {
        sequence[s.start - 1..s.end].fill(b'N');
    }
    sequence.extend([b'N'; 30]);
    sequence.extend(dna(&p, "tRNA-Pro"));
    let cox1 = start("COX1") + 100;
    let evidence = format!(
        "sample\t{cox1}\t{}\tlocal_1\tcompeting_paths\t2\n",
        cox1 + 10
    );
    let report = quality_run(&temp, &sequence, Some(&evidence));
    let q = &report["records"][0]["gene_quality"];
    assert_eq!(q["status"], "gene_review_required");
    for (name, issue) in [
        ("ND4", "cds_boundary_or_frame_issue"),
        ("rrnS", "ambiguous_sequence"),
        ("tRNA-Ala", "not_detected"),
        ("tRNA-Pro", "multiple_homologous_placements"),
        ("COX1", "overlaps_assembly_ambiguity"),
    ] {
        let g = q["genes"]
            .as_array()
            .unwrap()
            .iter()
            .find(|g| g["gene"] == name)
            .unwrap();
        assert_eq!(g["status"], "review_required", "{name}: {g}");
        assert!(
            g["issues"].as_array().unwrap().contains(&json!(issue)),
            "{name}: {g}"
        );
    }
}

#[test]
fn annotation_round_trip_preserves_topology_and_rotates_assembly_uncertainty() {
    let reference = reference("cel");
    let n = reference.len();
    for reverse in [false, true] {
        let temp = Temp::new();
        let mut sequence = reference.clone();
        sequence.rotate_left(13700);
        if reverse {
            sequence = reverse_complement(&sequence);
        }
        let input = temp.0.join("input.fasta");
        fs::write(&input,format!(">sample circular=false status=unresolved final_verified=false annotation_topology=linear\n{}\n",String::from_utf8_lossy(&sequence))).unwrap();
        // This one NCR1 block straddles the input origin and becomes one
        // contiguous interval at 13680..13720 after tRNA-Pro normalization.
        let mut positions: Vec<_> = (13679..13720)
            .map(|p| {
                let p = (p + n - 13700) % n;
                if reverse {
                    n - 1 - p
                } else {
                    p
                }
            })
            .collect();
        positions.sort_unstable();
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        for p in positions {
            if let Some(last) = ranges.last_mut().filter(|last| last.1 + 1 == p + 1) {
                last.1 = p + 1;
            } else {
                ranges.push((p + 1, p + 1));
            }
        }
        let evidence = temp.0.join("primary_ambiguity.tsv");
        let rows: String = ranges
            .iter()
            .map(|(a, b)| format!("sample\t{a}\t{b}\tblock_1\tcompeting_cycle_paths\t2\n"))
            .collect();
        fs::write(
            &evidence,
            format!("sequence\tstart\tend\tevidence\treason\talternatives\n{rows}"),
        )
        .unwrap();
        let annotate = |fasta: &std::path::Path,
                        uncertainty: &std::path::Path,
                        out: &str,
                        topology: Option<&str>| {
            let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_mitovarr"));
            command
                .env("PATH", "")
                .args(["annotate", "--fasta"])
                .arg(fasta)
                .args([
                    "--species",
                    "cel",
                    "--rotate-to-trna-pro",
                    "--assembly-ambiguity",
                ])
                .arg(uncertainty)
                .arg("--out")
                .arg(temp.0.join(out));
            if let Some(t) = topology {
                command.args(["--topology", t]);
            }
            let result = command.output().unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            serde_json::from_slice::<Value>(
                &fs::read(temp.0.join(out).join("report.json")).unwrap(),
            )
            .unwrap()
        };
        let first = annotate(&input, &evidence, "first", Some("circular"));
        let record = &first["records"][0];
        assert_eq!(record["rotation"]["reverse_complemented"], reverse);
        assert_eq!(record["ncr"]["loci"][0]["status"], "annotated");
        assert_eq!(
            record["assembly_evidence"]["regions"]
                .as_array()
                .unwrap()
                .len(),
            1
        );
        let region = &record["assembly_evidence"]["regions"][0];
        assert_eq!(region["segments"], json!([{"start":13680,"end":13720}]));
        assert_eq!(region["overlapping_ncrs"], json!(["NCR1"]));
        assert_eq!(
            record["assembly_evidence"]["ncr_loci"][0]["evidence_status"],
            "unresolved_path_alternatives"
        );
        assert_eq!(first["assembly_status"], "not_modified");
        // Unresolved NCR1 evidence alone does not fail the gene-level checks.
        assert_eq!(record["gene_quality"]["status"], "no_detected_gene_issues");
        assert_eq!(record["gene_quality"]["ncr1_required"], false);
        assert_eq!(record["gene_quality"]["assembly_ambiguity_supplied"], true);
        let annotated = temp.0.join("first/annotated.fasta");
        assert_eq!(
            load_reference(&annotated, None).unwrap().sequence,
            reference
        );
        let text = fs::read_to_string(&annotated).unwrap();
        assert!(text.contains("circular=false status=unresolved final_verified=false"));
        assert_eq!(text.matches("annotation_topology=circular").count(), 1);
        let gff = fs::read_to_string(temp.0.join("first/assembly_ambiguity.gff3")).unwrap();
        assert!(gff.contains("\tregion\t13680\t13720\t"));
        let projected = temp.0.join("first/assembly_ambiguity.tsv");
        let second = annotate(&annotated, &projected, "second", None);
        assert_eq!(second["records"][0]["circular_input_assumed"], true);
        assert_eq!(second["records"][0]["ncr"]["loci"], record["ncr"]["loci"]);
        for (a, b) in second["records"][0]["ncr"]["regions"]
            .as_array()
            .unwrap()
            .iter()
            .zip(record["ncr"]["regions"].as_array().unwrap())
        {
            for key in [
                "locus",
                "segments",
                "length",
                "boundary_status",
                "ambiguous_bases",
                "left_genes",
                "right_genes",
            ] {
                assert_eq!(a[key], b[key]);
            }
        }
        assert_eq!(
            second["records"][0]["assembly_evidence"]["regions"][0]["segments"],
            region["segments"]
        );
        let forced = annotate(&annotated, &projected, "linear", Some("linear"));
        assert_eq!(forced["records"][0]["circular_input_assumed"], false);
        assert!(fs::read_to_string(temp.0.join("linear/annotated.fasta"))
            .unwrap()
            .contains("annotation_topology=linear"));
    }
}

#[test]
fn absent_and_short_ncr_loci_keep_junction_uncertainty_after_rotation() {
    let mut reference = reference("cel");
    // The same zero-gap geometry reproduced by the real CB4856 NCR1 holdout.
    reference.truncate(reference.len() - 466);
    let prediction = predict(&reference, &Panel::bundled(Species::Cel).unwrap(), true);
    let derived =
        mitovarr::annotation::ncr::derive(&reference, &prediction, true, Some(Species::Cel), 1);
    let short = derived
        .regions
        .iter()
        .find(|r| r.locus.as_deref() == Some("NCR2"))
        .unwrap();
    let n = reference.len();
    let mut marked: Vec<_> = [(0, "junction"), (n - 1, "junction"), (25, "gene_interior")].into();
    marked.extend(
        short
            .segments
            .iter()
            .flat_map(|s| (s.start - 1..s.end).map(|p| (p, "short_gap"))),
    );
    let nd3 = prediction
        .features
        .iter()
        .find(|f| f.name == "ND3")
        .unwrap();
    let nd5 = prediction
        .features
        .iter()
        .find(|f| f.name == "ND5")
        .unwrap();
    for a in &nd3.segments {
        for b in &nd5.segments {
            for position in a.start.max(b.start)..=a.end.min(b.end) {
                marked.push((position - 1, "overlap"));
            }
        }
    }
    for reverse in [false, true] {
        let temp = Temp::new();
        let mut sequence = reference.clone();
        sequence.rotate_left(73);
        if reverse {
            sequence = reverse_complement(&sequence);
        }
        fs::write(
            temp.0.join("input.fasta"),
            format!(
                ">sample status=unresolved\n{}\n",
                String::from_utf8_lossy(&sequence)
            ),
        )
        .unwrap();
        let mut tsv = "sequence\tstart\tend\tevidence\treason\talternatives\n".to_owned();
        for &(p, reason) in &marked {
            let p = (p + n - 73) % n;
            let p = if reverse { n - 1 - p } else { p } + 1;
            tsv +=
                &format!("sample\t{p}\t{p}\t{reason}\tinsufficient_path_or_junction_support\t1\n");
        }
        fs::write(temp.0.join("evidence.tsv"), tsv).unwrap();
        let result = std::process::Command::new(env!("CARGO_BIN_EXE_mitovarr"))
            .env("PATH", "")
            .args([
                "annotate",
                "--species",
                "cel",
                "--topology",
                "circular",
                "--rotate-to-trna-pro",
                "--fasta",
            ])
            .arg(temp.0.join("input.fasta"))
            .arg("--assembly-ambiguity")
            .arg(temp.0.join("evidence.tsv"))
            .arg("--out")
            .arg(temp.0.join("annotation"))
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let report: Value =
            serde_json::from_slice(&fs::read(temp.0.join("annotation/report.json")).unwrap())
                .unwrap();
        let record = &report["records"][0];
        assert_eq!(record["ncr"]["loci"][0]["status"], "adjacent_without_gap");
        assert_eq!(record["ncr"]["loci"][1]["status"], "below_min_length");
        assert_eq!(record["ncr"]["loci"][3]["status"], "overlapping_flanks");
        for index in [0, 1, 3] {
            assert_eq!(
                record["assembly_evidence"]["ncr_loci"][index]["evidence_status"],
                "retained_assembly_alternatives"
            );
        }
        for region in record["assembly_evidence"]["regions"].as_array().unwrap() {
            let expected = match region["evidence"].as_str().unwrap() {
                "junction" => json!(["NCR1"]),
                "short_gap" => json!(["NCR2"]),
                "overlap" => json!(["NCR4"]),
                _ => json!([]),
            };
            assert_eq!(region["overlapping_ncrs"], expected);
        }
        assert_eq!(
            load_reference(&temp.0.join("annotation/annotated.fasta"), None)
                .unwrap()
                .sequence,
            reference
        );
    }
}

#[test]
fn annotation_cli_works_offline_without_path_and_preserves_candidate_status() {
    let temp = Temp::new();
    let p = panel(&["tRNA-Pro"]);
    let q = dna(&p, "tRNA-Pro");
    fs::write(temp.0.join("input.fa"),format!(">candidate_1 circular=false status=unresolved\n{}\n>candidate_2 circular=true status=unresolved\n{}\n",String::from_utf8_lossy(&q),String::from_utf8_lossy(&[vec![b'N';17],q.clone(),vec![b'N';33]].concat()))).unwrap();
    fs::write(temp.0.join("panel.json"),serde_json::to_vec(&json!({"schema":p.schema,"accession":p.accession,"genetic_code":5,"genes":[{"name":"tRNA-Pro","kind":"tRNA","sequence":String::from_utf8_lossy(&q)}]})).unwrap()).unwrap();
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_mitovarr"))
        .env("PATH", "")
        .args(["annotate", "--fasta"])
        .arg(temp.0.join("input.fa"))
        .arg("--gene-models")
        .arg(temp.0.join("panel.json"))
        .arg("--out")
        .arg(temp.0.join("result"))
        .arg("--rotate-to-trna-pro")
        .status()
        .unwrap();
    assert!(status.success());
    let r: Value =
        serde_json::from_slice(&fs::read(temp.0.join("result/report.json")).unwrap()).unwrap();
    assert_eq!(r["assembly_status"], "not_modified");
    assert_eq!(r["records"][0]["rotation"]["reason"], "linear_input");
    assert_eq!(r["records"][1]["rotation"]["applied"], true);
    assert!(fs::read_to_string(temp.0.join("result/annotated.fasta"))
        .unwrap()
        .contains("status=unresolved"));
    assert!(!temp.0.join("result/consensus.fasta").exists());
    for file in [
        "input.fasta",
        "annotated.fasta",
        "annotations.gff3",
        "features.tsv",
        "genes.fna",
        "proteins.faa",
        "ncr.gff3",
        "ncr.tsv",
        "ncr_loci.tsv",
        "ncr.fasta",
        "report.json",
    ] {
        assert!(temp.0.join("result").join(file).exists());
    }
    // A custom model panel has no implicit species-specific positional catalog.
    assert_eq!(r["records"][1]["ncr"]["catalog"], Value::Null);
    assert_eq!(r["records"][1]["ncr"]["loci"], json!([]));
    assert_eq!(
        r["records"][1]["ncr"]["regions"][0]["name"],
        "unclassified_gap"
    );
    assert_eq!(r["records"][1]["ncr"]["regions"][0]["length"], 50);
    assert!(fs::read_to_string(temp.0.join("result/ncr.gff3"))
        .unwrap()
        .contains("Name=unclassified_gap"));
    assert_eq!(
        load_reference(&temp.0.join("result/ncr.fasta"), None)
            .unwrap()
            .sequence,
        vec![b'N'; 50]
    );
}

#[test]
fn a_complete_gene_does_not_hide_an_adjacent_partial_copy() {
    let p = panel(&["ND4"]);
    let q = dna(&p, "ND4");
    let seq = [
        q.clone(),
        vec![b'N'; 20],
        q[q.len() - 300..].to_vec(),
        vec![b'N'; 40],
    ]
    .concat();
    let a = predict(&seq, &p, false);
    assert_eq!(
        a.features
            .iter()
            .filter(|f| f.name == "ND4" && f.full_length_supported)
            .count(),
        1
    );
    let copy = a
        .features
        .iter()
        .find(|f| f.five_prime_1based > q.len())
        .unwrap();
    assert!(!copy.full_length_supported);
    assert_eq!(copy.length, 300);
    assert_eq!(copy.query_start, q.len() - 299);
}

#[test]
fn an_internal_rna_length_variant_is_not_certified_by_homology_alone() {
    let p = panel(&["tRNA-Pro"]);
    let mut q = dna(&p, "tRNA-Pro");
    q.remove(25);
    let mut a = predict(&q, &p, true);
    assert_eq!(a.features[0].status, "boundary_complete_length_variant");
    assert_eq!(a.features[0].rna_structure.as_deref(), Some("not_assessed"));
    assert!(!rotate_to_pro(&q, &mut a, true, true).1.applied);
}

#[test]
fn rotation_preserves_iupac_ambiguities_outside_the_anchor() {
    let p = panel(&["tRNA-Pro"]);
    let q = dna(&p, "tRNA-Pro");
    let input = [b"RYSWKMBDHVN".to_vec(), reverse_complement(&q)].concat();
    let mut a = predict(&input, &p, true);
    let (out, r) = rotate_to_pro(&input, &mut a, true, true);
    assert!(r.applied);
    assert!(r.reverse_complemented);
    assert_eq!(out, [q, b"NBDHVKMWSRY".to_vec()].concat());
}

/// The pinned cel reference cut inside NCR1 (reference bases 13329..13794), so
/// a linear contig has NCR1 at both ends: 294 bases before tRNA-Pro and 172
/// after tRNA-Ala.
fn contig_cut_in_ncr1() -> Vec<u8> {
    let reference = reference("cel");
    [&reference[13500..], &reference[..13500]].concat()
}

fn assembly_report(length: usize, circular: bool, closure: &str) -> Value {
    json!({"schema":"mitovar-rust.assembly.v1","backend":"rust-de-bruijn",
        "status":"unresolved","final_published":false,
        "mixed_sample":{"status":"not_detected","structural_alternatives":[]},
        "primary":{"file":"primary.fasta","length":length,"circular":circular,"final_verified":false,
            "closure":{"state":closure},"read_validation":{"status":"performed"},
            "collapsed_bubbles":{"local_intervals":[1],"not_placed":[2,3]}}})
}

struct LevelRun {
    level: Value,
    table: Option<String>,
}

/// Annotate `sequence` as the reported primary and return its claim level.
fn level_run(sequence: &[u8], topology: &str, rows: &str, report: Option<&Value>) -> LevelRun {
    level_run_named("primary", sequence, topology, rows, report)
}

fn level_run_named(
    name: &str,
    sequence: &[u8],
    topology: &str,
    rows: &str,
    report: Option<&Value>,
) -> LevelRun {
    let temp = Temp::new();
    let input = temp.0.join("primary.fasta");
    fs::write(
        &input,
        format!(
            ">{name} circular={} length={} final_verified=false\n{}\n",
            topology == "circular",
            sequence.len(),
            String::from_utf8_lossy(sequence)
        ),
    )
    .unwrap();
    let evidence = temp.0.join("primary_ambiguity.tsv");
    let ends = if topology == "linear" {
        format!(
            "{name}\t1\t1\tgraph_end\topen_graph_end\t0\n{name}\t{n}\t{n}\tgraph_end\topen_graph_end\t0\n",
            n = sequence.len()
        )
    } else {
        String::new()
    };
    fs::write(
        &evidence,
        format!("sequence\tstart\tend\tevidence\treason\talternatives\n{ends}{rows}"),
    )
    .unwrap();
    let out = temp.0.join("annotation");
    let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_mitovarr"));
    command
        .env("PATH", "")
        .args(["annotate", "--species", "cel", "--rotate-to-trna-pro"])
        .args(["--topology", topology])
        .arg("--fasta")
        .arg(&input)
        .arg("--assembly-ambiguity")
        .arg(&evidence)
        .arg("--out")
        .arg(&out);
    if let Some(report) = report {
        // A run writes the same intervals to its table and its report; a test
        // that sets the report's own intervals checks a mismatch.
        let mut report = report.clone();
        if report["primary"].is_object() && report["primary"].get("uncertainty_intervals").is_none()
        {
            report["primary"]["uncertainty_intervals"] = format!("{ends}{rows}")
                .lines()
                .map(|line| {
                    let c: Vec<_> = line.split('\t').collect();
                    json!({"start":c[1].parse::<u64>().unwrap(),"end":c[2].parse::<u64>().unwrap(),
                        "evidence":c[3],"reason":c[4],"alternatives":c[5].parse::<u64>().unwrap()})
                })
                .collect();
        }
        let path = temp.0.join("report.json");
        fs::write(&path, serde_json::to_vec(&report).unwrap()).unwrap();
        command.arg("--assembly-report").arg(path);
    }
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let annotation: Value =
        serde_json::from_slice(&fs::read(out.join("report.json")).unwrap()).unwrap();
    // The level never changes what annotation says about the assembly.
    assert_eq!(annotation["assembly_status"], "not_modified");
    assert_eq!(annotation["closure_validation"], "not_performed");
    LevelRun {
        level: annotation["records"][0]["claim_level"].clone(),
        table: fs::read_to_string(out.join("claim_level.tsv")).ok(),
    }
}

/// One-based position 100 bases into a gene of `sequence`, read as linear.
fn inside(sequence: &[u8], gene: &str) -> usize {
    let p = Panel::bundled(Species::Cel).unwrap();
    predict(sequence, &p, false)
        .features
        .iter()
        .find(|f| f.name == gene && f.full_length_supported)
        .unwrap()
        .segments[0]
        .start
        + 100
}

#[test]
fn complete_except_ncr1_is_reached_by_one_contig_whatever_its_closure() {
    let contig = contig_cut_in_ncr1();
    let n = contig.len();
    let cox1 = inside(&contig, "COX1");
    let nd4 = inside(&contig, "ND4");
    // Single-site alternatives in two genes; inside NCR1 a tract, an open
    // structure and the thin depth of the two open ends.
    let rows = format!(
        "primary\t{cox1}\t{cox1}\tcollapsed_bubble_7\tcollapsed_bubble_alternative\t1\n\
         primary\t{nd4}\t{}\tlocus_3\tminor_path_retained\t1\n\
         primary\t100\t130\thomopolymer_family_2\thomopolymer_length_family\t1\n\
         primary\t150\t190\tlocus_9\tmajority_path_structure_unresolved\t2\n\
         primary\t1\t9\tread_validation_1\tlow_matching_read_depth\t0\n\
         primary\t{}\t{n}\tread_validation_2\tlow_matching_read_depth\t0\n",
        nd4 + 2,
        n - 8
    );
    for (closure, words) in [
        ("not_observed", "not shown to be circular"),
        (
            "joined_by_reads",
            "overlapping reads join the two NCR1 ends",
        ),
    ] {
        let run = level_run(
            &contig,
            "linear",
            &rows,
            Some(&assembly_report(n, false, closure)),
        );
        let level = &run.level;
        assert_eq!(level["level"], "complete_except_ncr1");
        assert_eq!(level["status"], "reached", "{level:#}");
        assert_eq!(level["closure"], closure);
        assert_eq!(level["strict_final_published"], false);
        assert!(level["not_reached"].as_array().unwrap().is_empty());
        assert!(level["not_assessed_reason"].is_null());
        assert!(
            level["statement"].as_str().unwrap().contains(words),
            "{level:#}"
        );
        // Never a circle or a complete genome for an unjoined contig.
        if closure == "not_observed" {
            let statement = level["statement"].as_str().unwrap();
            assert!(
                !statement.contains("verified circle") && !statement.contains("complete genome")
            );
        }
        assert_eq!(level["flagged_sites_outside_ncr1"], 2);
        assert_eq!(level["flagged_sites_inside_ncr1"], 1);
        assert_eq!(
            level["flagged_site_reasons_outside_ncr1"],
            json!({"collapsed_bubble_alternative":1,"minor_path_retained":1})
        );
        assert_eq!(level["collapsed_bubbles_not_placed"], 2);
        assert_eq!(level["ncr1_assembled_bases"], json!([294, 172]));
        assert_eq!(level["ncr1_length"], "not_determined");
        assert!(level["rna_length_variants"].as_array().unwrap().is_empty());
        let table = run.table.unwrap();
        let rows: Vec<Vec<&str>> = table.lines().map(|l| l.split('\t').collect()).collect();
        assert_eq!(rows.len(), 2);
        assert_eq!(
            rows[0][..5],
            [
                "contig",
                "level",
                "status",
                "closure",
                "strict_final_published"
            ]
        );
        assert_eq!(
            rows[1][..5],
            [
                "primary",
                "complete_except_ncr1",
                "reached",
                closure,
                "false"
            ]
        );
        assert_eq!(rows[0].len(), rows[1].len());
        assert_eq!(rows[1][8], "2");
        assert_eq!(rows[1][11], "294,172");
        assert_eq!(rows[1][12], "not_determined");
    }
}

#[test]
fn open_structure_low_depth_mixture_and_gene_problems_keep_the_level_out() {
    let contig = contig_cut_in_ncr1();
    let n = contig.len();
    let cox1 = inside(&contig, "COX1");
    let clean = assembly_report(n, false, "not_observed");
    let interval = |reason: &str| format!("primary\t{cox1}\t{}\tx_1\t{reason}\t1\n", cox1 + 30);
    for (rows, expected) in [
        (
            interval("majority_path_structure_unresolved"),
            vec!["open_structure_outside_ncr1:majority_path_structure_unresolved"],
        ),
        (
            interval("topology_inferred_repeat_pairing"),
            vec!["open_structure_outside_ncr1:topology_inferred_repeat_pairing"],
        ),
        (
            interval("competing_cycle_paths"),
            vec!["open_structure_outside_ncr1:competing_cycle_paths"],
        ),
        (interval("low_matching_read_depth"), vec!["low_read_depth_outside_ncr1"]),
        // A whole-primary flag of an older assembly cannot be placed.
        (
            format!("primary\t1\t{n}\tcollapsed_bubbles_1,2\tcollapsed_bubble_component_unresolved\t2\n"),
            vec!["open_structure_outside_ncr1:collapsed_bubble_component_unresolved"],
        ),
        // A reason this build does not know is open, not a site.
        (interval("future_reason"), vec!["open_structure_outside_ncr1:future_reason"]),
        // An interval that starts in NCR1 and runs into tRNA-Pro is outside.
        (
            "primary\t280\t300\tread_validation_1\tlow_matching_read_depth\t0\n".to_owned(),
            vec!["low_read_depth_outside_ncr1"],
        ),
    ] {
        let level = level_run(&contig, "linear", &rows, Some(&clean)).level;
        assert_eq!(level["status"], "not_reached", "{rows}");
        assert_eq!(level["not_reached"], json!(expected), "{rows}");
        assert_eq!(level["closure"], "not_observed");
        assert!(level["statement"].as_str().unwrap().starts_with("not complete outside NCR1"));
    }
    // The same intervals wholly inside NCR1 do not count.
    for reason in [
        "majority_path_structure_unresolved",
        "low_matching_read_depth",
        "future_reason",
    ] {
        let rows = format!("primary\t50\t90\tx_1\t{reason}\t1\n");
        let level = level_run(&contig, "linear", &rows, Some(&clean)).level;
        assert_eq!(level["status"], "reached", "{reason}: {level:#}");
    }
    // A mixed sample.
    let mut mixed = clean.clone();
    mixed["mixed_sample"] =
        json!({"status":"mixed_sample_signal","structural_alternatives":[{"bypassed_bases":400}]});
    let level = level_run(&contig, "linear", "", Some(&mixed)).level;
    assert_eq!(
        level["not_reached"],
        json!(["mixed_sample_signal", "structural_mixture"])
    );
    // A coding gene with a lost base is a gene problem; nothing else changes.
    let mut frameshift = contig.clone();
    frameshift.remove(cox1);
    let level = level_run(
        &frameshift,
        "linear",
        "",
        Some(&assembly_report(n - 1, false, "not_observed")),
    )
    .level;
    assert_eq!(level["status"], "not_reached");
    assert_eq!(level["not_reached"], json!(["gene:COX1"]));
    // A contig that stops inside the genome does not end in NCR1 on both sides.
    let short = contig[..n - 2000].to_vec();
    let level = level_run(
        &short,
        "linear",
        "",
        Some(&assembly_report(short.len(), false, "not_observed")),
    )
    .level;
    assert_eq!(level["status"], "not_reached");
    let reasons: Vec<_> = level["not_reached"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| r.as_str().unwrap().to_owned())
        .collect();
    assert!(
        reasons.contains(&"not_one_contig_ending_in_ncr1".to_owned()),
        "{reasons:?}"
    );
    assert!(
        reasons.iter().any(|r| r.starts_with("gene:")),
        "{reasons:?}"
    );
}

#[test]
fn an_rna_length_variant_is_listed_and_does_not_remove_the_level() {
    let mut contig = contig_cut_in_ncr1();
    let p = Panel::bundled(Species::Cel).unwrap();
    let lys = predict(&contig, &p, false)
        .features
        .iter()
        .find(|f| f.name == "tRNA-Lys")
        .unwrap()
        .segments[0]
        .start;
    // One base lost inside the tRNA: both boundaries still align.
    contig.remove(lys + 25);
    let report = assembly_report(contig.len(), false, "not_observed");
    let level = level_run(&contig, "linear", "", Some(&report)).level;
    assert_eq!(level["status"], "reached", "{level:#}");
    assert_eq!(level["rna_length_variants"], json!(["tRNA-Lys"]));
}

#[test]
fn a_circle_reaches_the_level_and_only_the_strict_gate_calls_it_verified() {
    let circle = reference("cel");
    let n = circle.len();
    let cox1 = inside(&circle, "COX1");
    let mut report = assembly_report(n, true, "circle");
    let level = level_run(&circle, "circular", "", Some(&report)).level;
    assert_eq!(level["status"], "reached", "{level:#}");
    assert_eq!(level["closure"], "circle");
    assert_eq!(level["strict_final_published"], false);
    assert!(level["statement"]
        .as_str()
        .unwrap()
        .contains("not a verified circle"));
    assert_eq!(level["ncr1_assembled_bases"], json!([466]));
    assert_eq!(level["ncr1_length"], "not_determined");
    // Without rotation a circle cut inside NCR1 still reports one NCR1 total.
    let split = contig_cut_in_ncr1();
    let temp = Temp::new();
    let input = temp.0.join("primary.fasta");
    fs::write(
        &input,
        format!(
            ">primary circular=true\n{}\n",
            String::from_utf8_lossy(&split)
        ),
    )
    .unwrap();
    let evidence = temp.0.join("primary_ambiguity.tsv");
    fs::write(
        &evidence,
        "sequence\tstart\tend\tevidence\treason\talternatives\n",
    )
    .unwrap();
    let assembly = temp.0.join("report.json");
    let mut split_report = report.clone();
    split_report["primary"]["uncertainty_intervals"] = json!([]);
    fs::write(&assembly, serde_json::to_vec(&split_report).unwrap()).unwrap();
    let out = temp.0.join("annotation");
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_mitovarr"))
        .args([
            "annotate",
            "--species",
            "cel",
            "--topology",
            "circular",
            "--fasta",
        ])
        .arg(&input)
        .arg("--assembly-ambiguity")
        .arg(&evidence)
        .arg("--assembly-report")
        .arg(&assembly)
        .arg("--out")
        .arg(&out)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let annotation: Value =
        serde_json::from_slice(&fs::read(out.join("report.json")).unwrap()).unwrap();
    let level = &annotation["records"][0]["claim_level"];
    assert_eq!(level["status"], "reached", "{level:#}");
    assert_eq!(level["ncr1_assembled_bases"], json!([466]));
    // Competing cycle paths in a gene are open structure; a minor cycle path is a site.
    let rows = format!(
        "primary\t{cox1}\t{}\tcircle_family_1_block_1\tcompeting_cycle_paths\t1\n",
        cox1 + 4
    );
    let level = level_run(&circle, "circular", &rows, Some(&report)).level;
    assert_eq!(
        level["not_reached"],
        json!(["open_structure_outside_ncr1:competing_cycle_paths"])
    );
    let rows = format!(
        "primary\t{cox1}\t{}\tcircle_family_1_block_1\tminor_cycle_path_retained\t1\n",
        cox1 + 4
    );
    let level = level_run(&circle, "circular", &rows, Some(&report)).level;
    assert_eq!(level["status"], "reached");
    assert_eq!(level["flagged_sites_outside_ncr1"], 1);
    // Only the assembly's own strict result makes it a verified circle.
    report["final_published"] = json!(true);
    report["primary"]["final_verified"] = json!(true);
    let level = level_run(&circle, "circular", "", Some(&report)).level;
    assert_eq!(level["status"], "reached");
    assert_eq!(level["strict_final_published"], true);
    assert!(level["statement"]
        .as_str()
        .unwrap()
        .starts_with("verified circle"));
}

#[test]
fn the_level_is_not_assessed_without_the_matching_assembly_evidence() {
    let contig = contig_cut_in_ncr1();
    let n = contig.len();
    let clean = assembly_report(n, false, "not_observed");
    // No report: the field says so and no table is written.
    let run = level_run(&contig, "linear", "", None);
    assert_eq!(run.level["status"], "not_assessed");
    assert_eq!(
        run.level["not_assessed_reason"],
        "assembly_report_not_supplied"
    );
    assert!(run.level["closure"].is_null());
    assert!(run.table.is_none());
    let edit = |change: &dyn Fn(&mut Value)| {
        let mut report = clean.clone();
        change(&mut report);
        report
    };
    for (report, topology, reason) in [
        (
            edit(&|r| r["schema"] = json!("mitovar-rust.native-mapping.v1")),
            "linear",
            "not_a_de_novo_assembly_report",
        ),
        (
            edit(&|r| r["primary"] = Value::Null),
            "linear",
            "assembly_has_no_primary",
        ),
        (
            edit(&|r| r["primary"]["length"] = json!(n + 1)),
            "linear",
            "input_is_not_the_reported_primary",
        ),
        (
            clean.clone(),
            "circular",
            "annotation_topology_differs_from_assembly",
        ),
        // A report written before this level existed.
        (
            edit(&|r| {
                r["primary"].as_object_mut().unwrap().remove("closure");
            }),
            "linear",
            "assembly_report_lacks_closure",
        ),
        (
            edit(&|r| {
                r["primary"]
                    .as_object_mut()
                    .unwrap()
                    .remove("read_validation");
            }),
            "linear",
            "assembly_report_lacks_primary_read_validation",
        ),
        (
            edit(&|r| r["primary"]["read_validation"]["status"] = json!("resource_limited")),
            "linear",
            "primary_read_validation_incomplete",
        ),
        (
            edit(&|r| r["primary"]["closure"] = json!({"state":"not_assessed"})),
            "linear",
            "closure_not_assessed",
        ),
    ] {
        let run = level_run(&contig, topology, "", Some(&report));
        assert_eq!(run.level["status"], "not_assessed", "{reason}");
        assert_eq!(run.level["not_assessed_reason"], reason);
        assert!(run.level["not_reached"].as_array().unwrap().is_empty());
        let table = run.table.unwrap();
        let row: Vec<_> = table.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(row[2], "not_assessed");
        assert_eq!(row[3], "");
        assert_eq!(row[6], reason);
    }
    // Another record of the same length is not the reported primary.
    let run = level_run_named("candidate_1", &contig, "linear", "", Some(&clean));
    assert_eq!(
        run.level["not_assessed_reason"],
        "input_is_not_the_reported_primary"
    );
    // The report alone is refused: the level needs the interval table.
    let temp = Temp::new();
    let input = temp.0.join("primary.fasta");
    fs::write(
        &input,
        format!(">primary\n{}\n", String::from_utf8_lossy(&contig)),
    )
    .unwrap();
    let report = temp.0.join("report.json");
    fs::write(&report, serde_json::to_vec(&clean).unwrap()).unwrap();
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_mitovarr"))
        .args(["annotate", "--species", "cel", "--fasta"])
        .arg(&input)
        .arg("--assembly-report")
        .arg(&report)
        .arg("--out")
        .arg(temp.0.join("annotation"))
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("--assembly-ambiguity"));
    assert!(!temp.0.join("annotation").exists());
}

/// Review finding 1 (2026-11-13): a tRNA-Pro or tRNA-Ala that is an accepted RNA
/// length variant at a linear end still ends the contig in NCR1 (panel strain
/// cel ECA1210). Test written by the reviewing session.
#[test]
fn an_end_trna_length_variant_still_ends_the_contig_in_ncr1() {
    let p = Panel::bundled(Species::Cel).unwrap();
    for gene in ["tRNA-Pro", "tRNA-Ala"] {
        let mut contig = contig_cut_in_ncr1();
        let start = predict(&contig, &p, false)
            .features
            .iter()
            .find(|f| f.name == gene)
            .unwrap()
            .segments[0]
            .start;
        // One base lost inside the tRNA: both boundaries still align.
        contig.remove(start + 25);
        let features = predict(&contig, &p, false).features;
        let hit = features.iter().find(|f| f.name == gene).unwrap();
        assert_eq!(hit.status, "boundary_complete_length_variant", "{gene}");
        let report = assembly_report(contig.len(), false, "not_observed");
        let level = level_run(&contig, "linear", "", Some(&report)).level;
        assert_eq!(level["rna_length_variants"], json!([gene]), "{level:#}");
        assert_eq!(level["status"], "reached", "{gene}: {level:#}");
    }
}

/// Review finding 4 (2026-11-13): the table must be the report's own. A report
/// whose read validation found low-depth bases inside COX1, with a table that
/// lacks that row (for example the v0.2.2 table of the same primary), is not
/// assessed. Test written by the reviewing session.
#[test]
fn a_table_that_does_not_match_the_report_is_not_assessed() {
    let contig = contig_cut_in_ncr1();
    let n = contig.len();
    let cox1 = inside(&contig, "COX1");
    let mut report = assembly_report(n, false, "not_observed");
    report["primary"]["read_validation"] = json!({"status":"performed","low_depth_bases":31,
        "low_depth_intervals":[{"start":cox1,"end":cox1 + 30}]});
    report["primary"]["uncertainty_intervals"] = json!([
        {"start":1,"end":1,"evidence":"graph_end","reason":"open_graph_end","alternatives":0},
        {"start":cox1,"end":cox1 + 30,"evidence":"read_validation_1","reason":"low_matching_read_depth","alternatives":0},
        {"start":n,"end":n,"evidence":"graph_end","reason":"open_graph_end","alternatives":0}]);
    let level = level_run(&contig, "linear", "", Some(&report)).level;
    assert_eq!(level["status"], "not_assessed", "{level:#}");
    assert_eq!(
        level["not_assessed_reason"],
        "assembly_ambiguity_differs_from_report"
    );
    // With the matching table the low-depth bases keep the level out.
    let rows = format!(
        "primary\t{cox1}\t{}\tread_validation_1\tlow_matching_read_depth\t0\n",
        cox1 + 30
    );
    let level = level_run(&contig, "linear", &rows, Some(&report)).level;
    assert_eq!(level["status"], "not_reached", "{level:#}");
    assert_eq!(level["not_reached"], json!(["low_read_depth_outside_ncr1"]));
    // A report without intervals cannot be matched either.
    report["primary"]
        .as_object_mut()
        .unwrap()
        .remove("uncertainty_intervals");
    report["primary"]["uncertainty_intervals"] = Value::Null;
    let level = level_run(&contig, "linear", &rows, Some(&report)).level;
    assert_eq!(
        level["not_assessed_reason"],
        "assembly_ambiguity_differs_from_report"
    );
}
