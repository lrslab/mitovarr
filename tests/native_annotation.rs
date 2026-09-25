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
