use clap::Parser;
use flate2::{write::GzEncoder, Compression};
use mitovarr::{
    cli::Cli,
    io::{fasta::load_reference, fastq},
    kmer::reverse_complement,
    run,
};
use serde_json::Value;
use std::{
    fs,
    io::Write,
    num::NonZeroUsize,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};

static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "mitovar-native-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn dna(length: usize, mut seed: u64) -> Vec<u8> {
    (0..length)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed >> 32) as usize & 3]
        })
        .collect()
}

fn circle_reads(sequence: &[u8]) -> Vec<(String, Vec<u8>, Vec<u8>)> {
    circle_reads_with_length(sequence, 100)
}

fn circle_reads_with_length(
    sequence: &[u8],
    read_length: usize,
) -> Vec<(String, Vec<u8>, Vec<u8>)> {
    (0..sequence.len())
        .map(|start| {
            let r1 = (0..read_length)
                .map(|i| sequence[(start + i) % sequence.len()])
                .collect();
            let r2: Vec<_> = (0..read_length)
                .map(|i| sequence[(start + i + 70) % sequence.len()])
                .collect();
            (format!("read{start}"), r1, reverse_complement(&r2))
        })
        .collect()
}

fn write_reads(dir: &Fixture, records: &[(String, Vec<u8>, Vec<u8>)], gzip: bool) {
    for mate in [1, 2] {
        let mut bytes = Vec::new();
        for (name, a, b) in records {
            let sequence = if mate == 1 { a } else { b };
            writeln!(
                bytes,
                "@{name}/{mate} original metadata\n{}\n+\n{}",
                String::from_utf8_lossy(sequence),
                "I".repeat(sequence.len())
            )
            .unwrap();
        }
        let path = dir.path(&format!("r{mate}.fq"));
        if gzip {
            // Concatenated gzip members exercise native multi-member decompression.
            let split = bytes.len() / 2;
            let mut file = fs::File::create(path).unwrap();
            for chunk in [&bytes[..split], &bytes[split..]] {
                let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
                encoder.write_all(chunk).unwrap();
                file.write_all(&encoder.finish().unwrap()).unwrap();
            }
        } else {
            fs::write(path, bytes).unwrap();
        }
    }
}

fn assembly(dir: &Fixture, command: &str, out: &str, extra: &[&str]) -> bool {
    let mut args = vec![
        "mitovarr".to_owned(),
        command.into(),
        "--r1".into(),
        dir.path("r1.fq").display().to_string(),
        "--r2".into(),
        dir.path("r2.fq").display().to_string(),
        "--out".into(),
        dir.path(out).display().to_string(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    run(Cli::try_parse_from(args).unwrap()).unwrap()
}

fn report(dir: &Fixture, out: &str) -> Value {
    serde_json::from_slice(&fs::read(dir.path(&format!("{out}/report.json"))).unwrap()).unwrap()
}

fn equivalent_circle(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    for sequence in [a.to_vec(), reverse_complement(a)] {
        if [sequence.as_slice(), sequence.as_slice()]
            .concat()
            .windows(b.len())
            .any(|w| w == b)
        {
            return true;
        }
    }
    false
}

#[test]
fn circle_gzip_order_independent_and_no_external_runtime() {
    let dir = Fixture::new();
    let genome = dna(600, 9123);
    let mut reads = circle_reads(&genome);
    write_reads(&dir, &reads, true);
    assert!(assembly(&dir, "assemble-denovo", "first", &[]));
    let consensus = load_reference(&dir.path("first/consensus.fasta"), None).unwrap();
    assert!(equivalent_circle(&genome, &consensus.sequence));
    let primary = load_reference(&dir.path("first/primary.fasta"), None).unwrap();
    assert_eq!(primary.sequence, consensus.sequence);
    assert_eq!(report(&dir, "first")["primary"]["final_verified"], true);
    let original = fs::read(dir.path("first/consensus.fasta")).unwrap();
    reads.reverse();
    write_reads(&dir, &reads, false);
    assert!(assembly(&dir, "assemble-denovo", "second", &[]));
    assert_eq!(
        original,
        fs::read(dir.path("second/consensus.fasta")).unwrap()
    );
    assert_eq!(report(&dir, "first")["status"], "complete");
    assert_eq!(
        fs::read(dir.path("first/primary.fasta")).unwrap(),
        fs::read(dir.path("second/primary.fasta")).unwrap()
    );
    assert_eq!(report(&dir, "first")["validation"]["scan"]["pairs"], 600);
    assert!(assembly(&dir, "assemble-baited", "bam_baited", &[]));
    assert_eq!(
        original,
        fs::read(dir.path("bam_baited/consensus.fasta")).unwrap()
    );
    let baited = report(&dir, "bam_baited");
    assert_eq!(baited["workflow"], "bam_baited_de_novo");
    assert_eq!(baited["input_scope"], "mt-bam-extracted");
    assert_eq!(baited["validation_scope"], "mt-bam-extracted");
    assert!(baited["recruitment"].is_null());
    // The only executable started by this test is our own binary, with no tools on PATH.
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_mitovarr"))
        .env("PATH", dir.path("no-programs"))
        .args(["assemble-denovo", "--r1"])
        .arg(dir.path("r1.fq"))
        .arg("--r2")
        .arg(dir.path("r2.fq"))
        .arg("--out")
        .arg(dir.path("empty_path"))
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(
        original,
        fs::read(dir.path("empty_path/consensus.fasta")).unwrap()
    );
}

#[test]
fn seeded_recruitment_converges_and_limit_blocks_publication() {
    let dir = Fixture::new();
    let genome = dna(600, 9123);
    let mut reads = circle_reads(&genome);
    reads.push(("unrelated".into(), dna(100, 6353), dna(100, 8791)));
    write_reads(&dir, &reads, true);
    fs::write(
        dir.path("seed.fa"),
        format!(">seed\n{}\n", String::from_utf8_lossy(&genome[210..260])),
    )
    .unwrap();
    let seed = dir.path("seed.fa").display().to_string();
    assert!(assembly(
        &dir,
        "assemble-seeded",
        "complete",
        &["--seed", &seed]
    ));
    let result = report(&dir, "complete");
    assert_eq!(result["recruitment"]["selected_pairs"], 600);
    assert_eq!(result["validation"]["scan"]["pairs"], 601);
    assert_eq!(result["recruitment"]["reason"], "converged");
    assert!(!assembly(
        &dir,
        "assemble-seeded",
        "limited",
        &["--seed", &seed, "--max-rounds", "1"]
    ));
    assert!(!dir.path("limited/consensus.fasta").exists());
    assert_eq!(
        report(&dir, "limited")["recruitment"]["reason"],
        "max_rounds"
    );
    assert!(!assembly(
        &dir,
        "assemble-seeded",
        "pairs_limited",
        &["--seed", &seed, "--max-pairs", "2"]
    ));
    assert_eq!(
        report(&dir, "pairs_limited")["recruitment"]["reason"],
        "max_pairs"
    );
    assert!(!dir.path("pairs_limited/consensus.fasta").exists());
}

#[test]
fn short_calibration_and_insufficient_coverage_do_not_publish() {
    let dir = Fixture::new();
    let repeat = dna(120, 12533);
    let genome = [
        dna(120, 5231),
        repeat.clone(),
        dna(150, 9712),
        repeat,
        dna(140, 1563),
    ]
    .concat();
    write_reads(&dir, &circle_reads(&genome), false);
    assert!(!assembly(&dir, "assemble-denovo", "repeat", &[]));
    assert_eq!(
        report(&dir, "repeat")["paired_path_resolution"]["fragment_model"]["status"],
        "not_estimable"
    );
    assert!(
        report(&dir, "repeat")["directed_branch_nodes"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert!(!dir.path("repeat/consensus.fasta").exists());
    assert!(dir.path("repeat/assembly_graph.gfa").exists());
    write_reads(&dir, &circle_reads(&dna(600, 9123)), false);
    assert!(!assembly(
        &dir,
        "assemble-denovo",
        "coverage",
        &["--min-depth", "1000"]
    ));
    assert_eq!(
        report(&dir, "coverage")["reason"],
        "insufficient_original_read_or_junction_support"
    );
    assert!(!dir.path("coverage/consensus.fasta").exists());
}

#[test]
fn read_spanning_repeat_bridges_reconstruct_copies_and_reject_conflicts() {
    let dir = Fixture::new();
    let a = dna(300, 91234);
    let b = dna(300, 72351);
    let c = dna(300, 19467);
    let repeat = dna(130, 15767);
    let genome = [
        a.clone(),
        repeat.clone(),
        b.clone(),
        repeat.clone(),
        c.clone(),
        repeat.clone(),
    ]
    .concat();
    let alternative = [a, repeat.clone(), c, repeat.clone(), b, repeat].concat();
    write_reads(&dir, &circle_reads(&genome), false);
    assert!(!assembly(&dir, "assemble-denovo", "unspanned", &[]));
    let mut reads = circle_reads_with_length(&genome, 250);
    write_reads(&dir, &reads, true);
    assert!(assembly(&dir, "assemble-denovo", "bridged", &[]));
    let assembled = load_reference(&dir.path("bridged/consensus.fasta"), None).unwrap();
    assert!(equivalent_circle(&genome, &assembled.sequence));
    assert_eq!(
        report(&dir, "bridged")["read_threading"]["resolved_directed_repeat_unitigs"],
        2
    );
    let mut rare_conflicts = reads.clone();
    rare_conflicts.extend(
        circle_reads_with_length(&alternative, 250)
            .into_iter()
            .step_by(50)
            .map(|(name, a, b)| (format!("rare_{name}"), a, b)),
    );
    rare_conflicts.reverse();
    write_reads(&dir, &rare_conflicts, false);
    assert!(assembly(&dir, "assemble-denovo", "rare_conflicts", &[]));
    assert_eq!(
        fs::read(dir.path("bridged/consensus.fasta")).unwrap(),
        fs::read(dir.path("rare_conflicts/consensus.fasta")).unwrap()
    );
    reads.extend(
        circle_reads_with_length(&alternative, 250)
            .into_iter()
            .map(|(name, a, b)| (format!("alternate_{name}"), a, b)),
    );
    write_reads(&dir, &reads, false);
    assert!(!assembly(&dir, "assemble-denovo", "conflicting", &[]));
    assert_eq!(
        report(&dir, "conflicting")["read_threading"]["resolved_directed_repeat_unitigs"],
        0
    );
    assert!(!dir.path("conflicting/consensus.fasta").exists());

    let tandem = [dna(10, 34978), dna(80, 98123).repeat(20)].concat();
    write_reads(&dir, &circle_reads_with_length(&tandem, 250), false);
    assert!(!assembly(&dir, "assemble-denovo", "tandem", &[]));
    assert!(!dir.path("tandem/consensus.fasta").exists());
}

#[test]
fn graph_kmers_through_127_preserve_repeat_paths_and_remove_short_repeat_branches() {
    let dir = Fixture::new();
    let repeat = dna(110, 791827);
    let genome = [dna(550, 672817), repeat.clone(), dna(570, 198317), repeat].concat();
    let mut reads = circle_reads_with_length(&genome, 180);
    write_reads(&dir, &reads, true);
    for k in [63, 65, 95, 127] {
        let name = format!("k{k}");
        assert!(assembly(
            &dir,
            "assemble-denovo",
            &name,
            &["--k", &k.to_string()]
        ));
        let sequence = load_reference(&dir.path(&format!("{name}/consensus.fasta")), None).unwrap();
        assert!(equivalent_circle(&genome, &sequence.sequence), "k={k}");
    }
    assert!(
        report(&dir, "k63")["directed_branch_nodes"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert_eq!(report(&dir, "k127")["directed_branch_nodes"], 0);
    assert!(
        report(&dir, "k95")["read_threading"]["resolved_directed_repeat_unitigs"]
            .as_u64()
            .unwrap()
            > 0
    );
    reads.reverse();
    write_reads(&dir, &reads, false);
    assert!(assembly(
        &dir,
        "assemble-denovo",
        "reordered",
        &["--k", "127"]
    ));
    assert_eq!(
        fs::read(dir.path("k127/consensus.fasta")).unwrap(),
        fs::read(dir.path("reordered/consensus.fasta")).unwrap()
    );
    let gfa = fs::read_to_string(dir.path("k127/assembly_graph.gfa")).unwrap();
    assert!(gfa
        .lines()
        .any(|line| line.starts_with("L\t") && line.ends_with("126M")));
    for k in ["128", "129", "64"] {
        assert!(Cli::try_parse_from([
            "mitovarr",
            "assemble-denovo",
            "--r1",
            "a.fq",
            "--r2",
            "b.fq",
            "--out",
            "out",
            "--k",
            k
        ])
        .is_err());
    }
}

#[test]
fn larger_k_resolves_read_spanned_repeat_without_a_reference() {
    let dir = Fixture::new();
    let repeat = dna(50, 8712);
    let genome = [
        dna(170, 3478),
        repeat.clone(),
        dna(160, 1942),
        repeat,
        dna(140, 9831),
    ]
    .concat();
    write_reads(&dir, &circle_reads(&genome), true);
    assert!(assembly(&dir, "assemble-denovo", "k31", &["--k", "31"]));
    assert!(
        report(&dir, "k31")["read_threading"]["resolved_directed_repeat_unitigs"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert!(assembly(&dir, "assemble-denovo", "k63", &["--k", "63"]));
    let assembled = load_reference(&dir.path("k63/consensus.fasta"), None).unwrap();
    assert!(equivalent_circle(&genome, &assembled.sequence));
}

#[test]
fn paired_distance_connects_repeat_that_neither_read_spans() {
    let dir = Fixture::new();
    let repeat = dna(130, 77123);
    let genome = [
        dna(500, 27151),
        repeat.clone(),
        dna(500, 37621),
        repeat.clone(),
        dna(500, 83921),
        repeat,
    ]
    .concat();
    let reads: Vec<_> = (0..genome.len())
        .map(|start| {
            let a = (0..100)
                .map(|i| genome[(start + i) % genome.len()])
                .collect();
            let b: Vec<_> = (0..100)
                .map(|i| genome[(start + 200 + i) % genome.len()])
                .collect();
            (format!("paired_{start}"), a, reverse_complement(&b))
        })
        .collect();
    write_reads(&dir, &reads, false);
    let complete = assembly(&dir, "assemble-denovo", "paired", &[]);
    assert!(
        complete,
        "{}",
        serde_json::to_string_pretty(&report(&dir, "paired")["paired_path_resolution"]["loci"][0])
            .unwrap()
    );
    let assembled = load_reference(&dir.path("paired/consensus.fasta"), None).unwrap();
    assert!(equivalent_circle(&genome, &assembled.sequence));
    let result = report(&dir, "paired");
    assert_eq!(
        result["read_threading"]["resolved_directed_repeat_unitigs"],
        0
    );
    assert_eq!(
        result["paired_path_resolution"]["fragment_model"]["median_span"],
        300
    );
    assert!(
        result["paired_path_resolution"]["selected_connections"]
            .as_u64()
            .unwrap()
            >= 3
    );
}

#[test]
fn unrelated_small_linear_component_does_not_veto_supported_circle() {
    let dir = Fixture::new();
    let genome = dna(600, 72341);
    let mut reads = circle_reads(&genome);
    let other = dna(150, 927323);
    for i in 0..10 {
        reads.push((
            format!("other_{i}"),
            other[..100].to_vec(),
            reverse_complement(&other[50..]),
        ));
    }
    write_reads(&dir, &reads, false);
    assert!(assembly(&dir, "assemble-denovo", "components", &[]));
    let result = report(&dir, "components");
    assert!(result["candidate_count"].as_u64().unwrap() > 1);
    let assembled = load_reference(&dir.path("components/consensus.fasta"), None).unwrap();
    assert!(equivalent_circle(&genome, &assembled.sequence));
}

#[test]
fn local_graph_cleanup_recovers_noisy_circle_but_retains_balanced_alleles() {
    let dir = Fixture::new();
    let genome = dna(1600, 9273);
    let mut alternate = genome.clone();
    alternate[400] = if alternate[400] == b'A' { b'C' } else { b'A' };
    alternate.remove(800);
    alternate.insert(1100, b'A');
    let mut reads = circle_reads(&genome);
    let minority: Vec<_> = circle_reads(&alternate)
        .into_iter()
        .step_by(20)
        .map(|(name, a, b)| (format!("error_{name}"), a, b))
        .collect();
    reads.extend(minority);
    // Recurrent high-quality terminal substitutions create supported short tips.
    let (name, mut a, b) = reads[20].clone();
    let last = a.len() - 1;
    a[last] = if a[last] == b'A' { b'C' } else { b'A' };
    for i in 0..3 {
        reads.push((format!("tip_{name}_{i}"), a.clone(), b.clone()));
    }
    write_reads(&dir, &reads, true);
    // With cleanup disabled and a zero accepted minor-route ratio, newly
    // observed SNP-minority fragments retain strict ambiguity. The primary
    // still reconstructs the majority sequence from the unchanged reads.
    assert!(!assembly(
        &dir,
        "assemble-denovo",
        "raw",
        &["--max-error-fraction", "0"]
    ));
    assert!(!dir.path("raw/consensus.fasta").exists());
    let unsimplified = load_reference(&dir.path("raw/primary.fasta"), None).unwrap();
    assert!(equivalent_circle(&genome, &unsimplified.sequence));
    assert!(report(&dir, "raw")["paired_path_resolution"]["loci"]
        .as_array()
        .unwrap()
        .iter()
        .any(|l| l["selected_path"].is_null()
            && l["paths"]
                .as_array()
                .unwrap()
                .iter()
                .any(
                    |p| p["spanning_interval"]["basis"] == "single_snp_shared_flanks"
                        && p["strong_votes_against_alternatives"]
                            .as_array()
                            .unwrap()
                            .iter()
                            .any(|v| v.as_u64().unwrap() >= 3)
                )));
    assert!(assembly(&dir, "assemble-denovo", "clean", &[]));
    let assembled = load_reference(&dir.path("clean/consensus.fasta"), None).unwrap();
    assert!(equivalent_circle(&genome, &assembled.sequence));
    let result = report(&dir, "clean");
    assert!(
        result["graph_simplification"]["tip_kmers_removed"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert!(
        result["graph_simplification"]["bubble_kmers_removed"]
            .as_u64()
            .unwrap()
            > 0
    );
    reads.reverse();
    write_reads(&dir, &reads, false);
    assert!(assembly(&dir, "assemble-denovo", "reordered", &[]));
    assert_eq!(
        fs::read(dir.path("clean/consensus.fasta")).unwrap(),
        fs::read(dir.path("reordered/consensus.fasta")).unwrap()
    );

    let mut mixed = circle_reads(&genome);
    mixed.extend(
        circle_reads(&alternate)
            .into_iter()
            .map(|(name, a, b)| (format!("allele_{name}"), a, b)),
    );
    write_reads(&dir, &mixed, false);
    assert!(!assembly(&dir, "assemble-denovo", "balanced", &[]));
    assert!(
        report(&dir, "balanced")["directed_branch_nodes"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert!(!dir.path("balanced/consensus.fasta").exists());
}

#[test]
fn parser_rejects_mismatches_truncation_and_wrong_mates() {
    let dir = Fixture::new();
    for (r1, r2) in [
        ("@x/1\nACGT\n+\nIIII\n", "@y/2\nACGT\n+\nIIII\n"),
        ("@x/1\nACGT\n+\nIIII\n", "@x/2\nACGT\n+\n"),
        ("@x/2\nACGT\n+\nIIII\n", "@x/1\nACGT\n+\nIIII\n"),
        ("@x/1\nACGT\n+\nIIII\n", ""),
    ] {
        fs::write(dir.path("a"), r1).unwrap();
        fs::write(dir.path("b"), r2).unwrap();
        assert!(fastq::scan_pairs(&dir.path("a"), &dir.path("b"), |_, _, _| Ok(())).is_err());
    }
    assert!(Cli::try_parse_from([
        "mitovarr",
        "assemble-denovo",
        "--r1",
        "a",
        "--r2",
        "b",
        "--out",
        "o",
        "--seed",
        "s"
    ])
    .is_err());
    assert!(Cli::try_parse_from([
        "mitovarr",
        "assemble-denovo",
        "--r1",
        "a",
        "--r2",
        "b",
        "--out",
        "o",
        "--k",
        "32"
    ])
    .is_err());
}

fn write_bam(path: &Path, sequence: &[u8], minor_at: Option<usize>) {
    use noodles_sam::{
        self as sam,
        alignment::{
            io::Write,
            record::{
                cigar::{op::Kind, Op},
                Flags, MappingQuality,
            },
            record_buf::{Cigar, QualityScores, Sequence},
            RecordBuf,
        },
        header::record::value::{map::ReferenceSequence, Map},
    };
    let header = sam::Header::builder()
        .set_header(
            Map::<sam::header::record::value::map::Header>::builder()
                .insert(
                    sam::header::record::value::map::header::tag::SORT_ORDER,
                    "coordinate",
                )
                .build()
                .unwrap(),
        )
        .add_reference_sequence(
            "MtDNA",
            Map::<ReferenceSequence>::new(NonZeroUsize::new(sequence.len()).unwrap()),
        )
        .build();
    let mut writer = noodles_bam::io::Writer::new(fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for i in 0..20 {
        let mut observed = sequence.to_vec();
        if let Some(pos) = minor_at.filter(|_| i < 9) {
            observed[pos] = b'A';
        }
        let cigar: Cigar = [Op::new(Kind::Match, sequence.len())].into_iter().collect();
        let mut record = RecordBuf::builder()
            .set_name(format!("fragment{i}"))
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::MIN)
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(cigar)
            .set_sequence(Sequence::from(observed))
            .set_quality_scores(QualityScores::from(vec![40; sequence.len()]))
            .build();
        writer.write_alignment_record(&header, &record).unwrap();
        if minor_at.is_some() && i < 9 {
            // A second overlapping observation must not turn nine A fragments
            // into eighteen votes and outvote the eleven C fragments.
            *record.flags_mut() = Flags::SEGMENTED | Flags::LAST_SEGMENT;
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }
    writer.try_finish().unwrap();
    let index = noodles_bam::fs::index(path).unwrap();
    noodles_bam::bai::fs::write(path.with_extension("bam.bai"), &index).unwrap();
}

#[test]
fn indexed_bam_consensus_and_refusal_to_overwrite() {
    let dir = Fixture::new();
    let genome = dna(600, 9123);
    fs::write(
        dir.path("reference.fa"),
        format!(">MtDNA\n{}\n", String::from_utf8_lossy(&genome)),
    )
    .unwrap();
    write_bam(&dir.path("input.bam"), &genome, None);
    let args = vec![
        "mitovarr".into(),
        "polish".into(),
        "--ref".into(),
        dir.path("reference.fa").display().to_string(),
        "--bam".into(),
        dir.path("input.bam").display().to_string(),
        "--out".into(),
        dir.path("polished").display().to_string(),
        "--emit-coverage".into(),
        "--emit-variants".into(),
    ];
    assert!(run(Cli::try_parse_from(&args).unwrap()).unwrap());
    assert_eq!(
        load_reference(&dir.path("polished/consensus.fasta"), None)
            .unwrap()
            .sequence,
        genome
    );
    let before = fs::read(dir.path("polished/consensus.fasta")).unwrap();
    assert!(run(Cli::try_parse_from(&args).unwrap()).is_err());
    assert_eq!(
        before,
        fs::read(dir.path("polished/consensus.fasta")).unwrap()
    );
    assert_eq!(report(&dir, "polished")["circularity_verified"], false);
}

#[test]
fn mapping_majority_reports_both_alleles_and_legacy_threshold_remains_available() {
    let dir = Fixture::new();
    let mut genome = dna(600, 41293);
    genome[300] = b'C';
    fs::write(
        dir.path("reference.fa"),
        format!(">MtDNA\n{}\n", String::from_utf8_lossy(&genome)),
    )
    .unwrap();
    write_bam(&dir.path("input.bam"), &genome, Some(300));
    for (out, extra, expected) in [
        ("majority", vec![], b'C'),
        ("legacy", vec!["--min-alt-frac", "0.7"], b'N'),
        (
            "separate",
            vec!["--min-snp-frac", "0.8", "--min-indel-frac", "0.6"],
            b'N',
        ),
    ] {
        let mut args = vec![
            "mitovarr".to_owned(),
            "polish".into(),
            "--ref".into(),
            dir.path("reference.fa").display().to_string(),
            "--bam".into(),
            dir.path("input.bam").display().to_string(),
            "--out".into(),
            dir.path(out).display().to_string(),
        ];
        args.extend(extra.into_iter().map(str::to_owned));
        assert_eq!(
            run(Cli::try_parse_from(args).unwrap()).unwrap(),
            expected != b'N'
        );
        let sequence = load_reference(&dir.path(&format!("{out}/candidate.fasta")), None)
            .unwrap()
            .sequence;
        assert_eq!(sequence[300], expected);
        let votes = fs::read_to_string(dir.path(&format!("{out}/base_votes.tsv"))).unwrap();
        assert!(votes.contains(&format!(
            "301\tC\t{}\t20\tA\t9\t0.450000\t0.450000\tfalse",
            char::from(expected)
        )));
        assert!(votes.contains(&format!(
            "301\tC\t{}\t20\tC\t11\t0.550000\t0.550000\tfalse",
            char::from(expected)
        )));
        assert_eq!(report(&dir, out)["base_voting"]["mixed_sites"], 1);
    }
    assert_eq!(
        report(&dir, "majority")["parameters"]["min_indel_frac"],
        0.7
    );
    assert_eq!(report(&dir, "legacy")["parameters"]["min_snp_frac"], 0.7);
    assert_eq!(report(&dir, "legacy")["parameters"]["min_indel_frac"], 0.7);
    assert_eq!(
        report(&dir, "separate")["parameters"]["min_indel_frac"],
        0.6
    );
}

#[test]
fn three_haplotype_mixture_cannot_publish_absent_majority_haplotype() {
    let dir = Fixture::new();
    let base = dna(1600, 8817923);
    let positions = [1100, 1250, 1400];
    let mut reads = Vec::new();
    for (allele, p) in positions.into_iter().enumerate() {
        let mut genome = base.clone();
        genome[p] = if genome[p] == b'A' { b'C' } else { b'A' };
        for start in 0..genome.len() {
            let a = (0..100)
                .map(|i| genome[(start + i) % genome.len()])
                .collect();
            let b: Vec<_> = (0..100)
                .map(|i| genome[(start + 150 + i) % genome.len()])
                .collect();
            reads.push((format!("hap{allele}_{start}"), a, reverse_complement(&b)));
        }
    }
    write_reads(&dir, &reads, false);
    assert!(!assembly(&dir, "assemble-denovo", "mixture", &[]));
    assert!(!dir.path("mixture/consensus.fasta").exists());
    let r = report(&dir, "mixture");
    let families = r["paired_path_resolution"]["complete_cycle_comparison"]["families"]
        .as_array()
        .unwrap();
    assert!(!families.is_empty());
    assert!(families.iter().all(|f| f["selected_path"].is_null()
        && !f["unique_biological_haplotype_claim"].as_bool().unwrap()));
    assert_eq!(r["primary"]["final_verified"], false);
    assert_eq!(r["primary"]["structural_uniqueness_claim"], false);
    let primary = load_reference(&dir.path("mixture/primary.fasta"), None).unwrap();
    assert!(primary.sequence.len() >= 1000);
    assert!(dir.path("mixture/primary_ambiguity.tsv").exists());
}
