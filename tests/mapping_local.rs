use clap::{CommandFactory, Parser};
use mitovarr::{
    cli::{Cli, Commands},
    io::fasta::load_reference,
    kmer::reverse_complement,
    run,
};
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        io::Write as _,
        record::{
            cigar::{op::Kind, Op},
            Flags, MappingQuality,
        },
        record_buf::{Cigar, QualityScores, Sequence},
        RecordBuf,
    },
    header::record::value::{map::ReferenceSequence, Map},
};
use serde_json::Value;
use std::{
    fs,
    io::Write,
    num::NonZeroUsize,
    path::{Path, PathBuf},
    sync::atomic::{AtomicUsize, Ordering},
};

static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "mapping-local-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&p).unwrap();
        Self(p)
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
fn dna(n: usize, mut state: u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state >> 32) as usize & 3]
        })
        .collect()
}

/// Make actual sample reads and an intentionally reference-limited BAM: each
/// alignment retains the longest contiguous reference piece and soft clips its
/// new junction. Insert-interior reads map to the old homologous source copy.
fn inputs(
    dir: &Fixture,
    reference: &[u8],
    sample: &[u8],
    coordinates: &[usize],
    noisy: bool,
    extra_source: bool,
) {
    inputs_with_gap(
        dir,
        reference,
        sample,
        coordinates,
        noisy,
        extra_source,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn inputs_with_gap(
    dir: &Fixture,
    reference: &[u8],
    sample: &[u8],
    coordinates: &[usize],
    noisy: bool,
    extra_source: bool,
    gap: Option<std::ops::Range<usize>>,
) {
    fs::write(
        dir.path("ref.fa"),
        format!(">MtDNA\n{}\n", String::from_utf8_lossy(reference)),
    )
    .unwrap();
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
            Map::<ReferenceSequence>::new(NonZeroUsize::new(reference.len()).unwrap()),
        )
        .build();
    let mut records = Vec::new();
    let mut one = fs::File::create(dir.path("r1.fq")).unwrap();
    let mut two = fs::File::create(dir.path("r2.fq")).unwrap();
    let starts: Vec<_> = (0..sample.len())
        .step_by(2)
        .chain(if extra_source { vec![320; 200] } else { vec![] })
        .collect();
    for (id, start) in starts.into_iter().enumerate() {
        if gap
            .as_ref()
            .is_some_and(|g| (0..180).any(|i| g.contains(&((start + i) % sample.len()))))
        {
            continue;
        }
        for (mate, out) in [(0, &mut one), (1, &mut two)] {
            let pos = start + mate * 80;
            let sequence: Vec<_> = (0..100).map(|i| sample[(pos + i) % sample.len()]).collect();
            let mapped: Vec<_> = (0..100)
                .map(|i| coordinates[(pos + i) % sample.len()])
                .collect();
            let (mut best_start, mut best_end, mut run_start) = (0, 0, 0);
            for i in 1..=100 {
                if i == 100 || mapped[i] != mapped[i - 1] + 1 {
                    if i - run_start > best_end - best_start {
                        best_start = run_start;
                        best_end = i;
                    }
                    run_start = i;
                }
            }
            let mut ops = Vec::new();
            if best_start > 0 {
                ops.push(Op::new(Kind::SoftClip, best_start));
            }
            ops.push(Op::new(Kind::Match, best_end - best_start));
            if best_end < 100 {
                ops.push(Op::new(Kind::SoftClip, 100 - best_end));
            }
            let record = RecordBuf::builder()
                .set_name(format!("fragment{id}"))
                .set_flags(if mate == 0 {
                    Flags::SEGMENTED | Flags::FIRST_SEGMENT
                } else {
                    Flags::SEGMENTED | Flags::LAST_SEGMENT | Flags::REVERSE_COMPLEMENTED
                })
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(mapped[best_start] + 1).unwrap())
                .set_mapping_quality(MappingQuality::new(60).unwrap())
                .set_cigar(ops.into_iter().collect::<Cigar>())
                .set_sequence(Sequence::from(sequence.clone()))
                .set_quality_scores(QualityScores::from(vec![40; 100]))
                .build();
            records.push(record);
            let mut read = if mate == 0 {
                sequence
            } else {
                reverse_complement(&sequence)
            };
            if noisy && id % 13 == 0 {
                read.insert(45, b'A');
            }
            writeln!(
                out,
                "@fragment{id}/{}\n{}\n+\n{}",
                mate + 1,
                String::from_utf8_lossy(&read),
                "I".repeat(read.len())
            )
            .unwrap();
        }
    }
    records.sort_by_key(|r| r.alignment_start());
    let path = dir.path("mt.bam");
    let mut writer = noodles_bam::io::Writer::new(fs::File::create(&path).unwrap());
    writer.write_header(&header).unwrap();
    for r in records {
        writer.write_alignment_record(&header, &r).unwrap();
    }
    writer.try_finish().unwrap();
    let index = noodles_bam::fs::index(&path).unwrap();
    noodles_bam::bai::fs::write(path.with_extension("bam.bai"), &index).unwrap();
}
fn assemble(dir: &Fixture, out: &str, extra: &[&str]) -> (bool, Vec<u8>, Value) {
    assemble_mode(dir, out, extra, true)
}

fn mapping_cli(dir: &Fixture, out: &str, extra: &[&str], deep: bool) -> Cli {
    let mut args = vec![
        "mitovarr".into(),
        "assemble-mapping".into(),
        "--ref".into(),
        dir.path("ref.fa").display().to_string(),
        "--bam".into(),
        dir.path("mt.bam").display().to_string(),
        "--out".into(),
        dir.path(out).display().to_string(),
        "--local-assembly".into(),
        "--r1".into(),
        dir.path("r1.fq").display().to_string(),
        "--r2".into(),
        dir.path("r2.fq").display().to_string(),
    ];
    if deep {
        args.push("--deep-validation".into());
    }
    args.extend(extra.iter().map(|s| s.to_string()));
    Cli::try_parse_from(args).unwrap()
}

fn assemble_mode(dir: &Fixture, out: &str, extra: &[&str], deep: bool) -> (bool, Vec<u8>, Value) {
    let status = run(mapping_cli(dir, out, extra, deep)).unwrap();
    let sequence = load_reference(&dir.path(&format!("{out}/candidate.fasta")), None)
        .unwrap()
        .sequence;
    let report = serde_json::from_slice(
        &fs::read(dir.path(&format!("{out}/local_assembly/report.json"))).unwrap(),
    )
    .unwrap();
    (status, sequence, report)
}
fn circular_eq(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && [a, a].concat().windows(b.len()).any(|s| s == b)
}

#[test]
fn local_graph_cap_retains_full_bam_consensus_without_claiming_fastq_validation() {
    let dir = Fixture::new();
    let sample = dna(2000, 124891);
    inputs(
        &dir,
        &sample,
        &sample,
        &(0..sample.len()).collect::<Vec<_>>(),
        false,
        false,
    );
    for (out, deep) in [("limited-fast", false), ("limited-deep", true)] {
        let (available, candidate, local) =
            assemble_mode(&dir, out, &["--local-max-kmers", "5"], deep);
        assert_eq!(available, !deep);
        assert_eq!(candidate, sample);
        assert_eq!(local["status"], "resource_limited");
        assert_eq!(local["resource_limit"]["observed"], 6);
        assert_eq!(local["resource_limit"]["limit"], 5);
        assert_eq!(local["named_pairs"], Value::Null);
        assert_eq!(local["sample_kmers"], Value::Null);
        assert_eq!(local["fastq_scan_complete"], false);
        assert_eq!(local["all_windows_resolved"], false);
        let root: Value =
            serde_json::from_slice(&fs::read(dir.path(&format!("{out}/report.json"))).unwrap())
                .unwrap();
        assert_eq!(root["strict_validation_passed"], false);
        assert_eq!(root["final_published"], false);
        assert_eq!(root["records"], sample.len());
        assert_eq!(root["subsampling"]["enabled"], false);
        assert!(dir.path(&format!("{out}/base_votes.tsv")).exists());
        assert_eq!(dir.path(&format!("{out}/consensus.fasta")).exists(), !deep);
        if !deep {
            assert_eq!(
                load_reference(&dir.path(&format!("{out}/consensus.fasta")), None)
                    .unwrap()
                    .sequence,
                sample
            );
        }
        let fasta = dir.path(&format!("{out}/candidate.fasta"));
        let uncertainty = dir.path(&format!("{out}/structural_ambiguity.tsv"));
        let annotated = dir.path(&format!("{out}/annotated"));
        assert!(run(Cli::parse_from([
            "mitovarr",
            "annotate",
            "--fasta",
            fasta.to_str().unwrap(),
            "--species",
            "cel",
            "--topology",
            "linear",
            "--assembly-ambiguity",
            uncertainty.to_str().unwrap(),
            "--out",
            annotated.to_str().unwrap()
        ]))
        .is_ok());
    }
}

#[test]
fn local_input_errors_are_not_misreported_as_a_resource_limited_consensus() {
    let dir = Fixture::new();
    let sample = dna(1000, 91887);
    inputs(
        &dir,
        &sample,
        &sample,
        &(0..sample.len()).collect::<Vec<_>>(),
        false,
        false,
    );
    fs::write(dir.path("r1.fq"), b"invalid FASTQ\n").unwrap();
    assert!(run(mapping_cli(&dir, "malformed", &[], false)).is_err());
    assert!(!dir.path("malformed/consensus.fasta").exists());
    assert!(!dir.path("malformed/candidate.fasta").exists());
    let root: Value =
        serde_json::from_slice(&fs::read(dir.path("malformed/report.json")).unwrap()).unwrap();
    assert_eq!(root["status"], "failed");
}
fn insertion(at: usize, noisy: bool) -> (Fixture, Vec<u8>, Vec<u8>) {
    let dir = Fixture::new();
    let reference = dna(3000, 821423);
    let mut insert = reference[300..540].to_vec();
    // A new-copy-specific allele close enough to its junction to phase.
    insert[25] = if insert[25] == b'A' { b'C' } else { b'A' };
    let sample = [&reference[..at], insert.as_slice(), &reference[at..]].concat();
    let coordinates: Vec<_> = (0..at).chain(300..540).chain(at..reference.len()).collect();
    inputs(&dir, &reference, &sample, &coordinates, noisy, false);
    (dir, reference, sample)
}

#[test]
fn deep_local_assembly_keeps_supported_q19_fastq_paths() {
    let (dir, reference, sample) = insertion(2100, false);
    for mate in [1, 2] {
        let path = dir.path(&format!("r{mate}.fq"));
        let text = fs::read_to_string(&path).unwrap();
        let changed = text
            .lines()
            .enumerate()
            .map(|(i, line)| {
                if i % 4 == 3 {
                    "4".repeat(line.len())
                } else {
                    line.to_owned()
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        fs::write(path, changed + "\n").unwrap();
    }
    let (complete, sequence, result) = assemble(&dir, "default", &[]);
    assert!(complete, "{result:#}");
    assert!(circular_eq(&sequence, &sample));
    assert_eq!(result["applied_replacements"], 1);
    let parameters: Value =
        serde_json::from_slice(&fs::read(dir.path("default/report.json")).unwrap()).unwrap();
    assert_eq!(parameters["parameters"]["min_baseq"], 0);
    let (_, filtered, result) = assemble(&dir, "q20", &["--min-baseq", "20"]);
    assert_eq!(filtered, reference);
    assert_eq!(result["applied_replacements"], 0);
}

#[test]
fn homologous_insertion_is_recovered_with_source_copy_and_gapped_reads() {
    let (dir, reference, sample) = insertion(2100, true);
    let (complete, sequence, report) = assemble(&dir, "local", &[]);
    assert!(complete, "{report:#}");
    assert!(
        circular_eq(&sequence, &sample),
        "{} vs {}\n{report:#}",
        sequence.len(),
        sample.len()
    );
    assert_eq!(report["applied_replacements"], 1);
    assert_eq!(
        load_reference(&dir.path("local/reference_candidate.fasta"), None)
            .unwrap()
            .sequence,
        reference
    );
    assert!(dir.path("local/local_assembly.gfa").is_file());
    assert!(dir.path("local/structural_variants.tsv").is_file());
    // Input order does not choose among the graph paths or count mates twice.
    for mate in [1, 2] {
        let path = dir.path(&format!("r{mate}.fq"));
        let text = fs::read_to_string(&path).unwrap();
        let lines: Vec<_> = text.lines().collect();
        let mut reversed = lines
            .chunks(4)
            .rev()
            .map(|r| r.join("\n"))
            .collect::<Vec<_>>()
            .join("\n");
        reversed.push('\n');
        fs::write(path, reversed).unwrap();
    }
    let (complete, reordered, _) = assemble(&dir, "reordered", &[]);
    assert!(complete);
    assert_eq!(sequence, reordered);
    // Artificially incomplete candidate searches cannot publish.
    let (complete, _, report) = assemble(&dir, "bounded", &["--local-max-span", "200"]);
    assert!(!complete, "{report:#}");
}

#[test]
fn origin_spanning_local_replacement_preserves_circular_sample() {
    let (dir, _, sample) = insertion(20, false);
    let (complete, sequence, report) = assemble(&dir, "local", &[]);
    assert!(complete, "{report:#}");
    assert!(circular_eq(&sequence, &sample), "{report:#}");
    assert_eq!(report["output_origin_shifted_to_supported_anchor"], true);
}

#[test]
fn long_deletion_is_assembled_from_junction_not_from_missing_coverage() {
    let dir = Fixture::new();
    let reference = dna(3000, 182314);
    let sample = [&reference[..1800], &reference[2100..]].concat();
    let coordinates: Vec<_> = (0..1800).chain(2100..3000).collect();
    inputs(&dir, &reference, &sample, &coordinates, false, false);
    let (complete, sequence, report) = assemble(&dir, "local", &[]);
    assert!(complete, "{report:#}");
    assert!(circular_eq(&sequence, &sample), "{report:#}");
}

#[test]
fn extra_source_coverage_does_not_create_an_insert() {
    let dir = Fixture::new();
    let reference = dna(3000, 891714);
    inputs(
        &dir,
        &reference,
        &reference,
        &(0..3000).collect::<Vec<_>>(),
        false,
        true,
    );
    let (complete, sequence, report) = assemble(&dir, "local", &[]);
    assert!(complete, "{report:#}");
    assert_eq!(sequence, reference);
    assert_eq!(report["applied_replacements"], 0);
}

#[test]
fn indistinguishable_long_repeat_copies_remain_unresolved() {
    let dir = Fixture::new();
    let mut reference = dna(3000, 429381);
    let unit = dna(220, 715391);
    reference[300..520].copy_from_slice(&unit);
    let insert = [unit.as_slice(), unit.as_slice()].concat();
    let sample = [&reference[..2100], insert.as_slice(), &reference[2100..]].concat();
    let coordinates: Vec<_> = (0..2100)
        .chain(300..520)
        .chain(300..520)
        .chain(2100..3000)
        .collect();
    inputs(&dir, &reference, &sample, &coordinates, false, false);
    let (complete, strict, report) = assemble(&dir, "local", &[]);
    assert!(!complete, "{report:#}");
    assert!(!Path::new(&dir.path("local/consensus.fasta")).exists());
    assert!(report["windows"]
        .as_array()
        .unwrap()
        .iter()
        .any(|w| w["status"] == "unresolved"));
    let limited = report["windows"]
        .as_array()
        .unwrap()
        .iter()
        .find(|w| w["reason"] == "graph_search_limit")
        .unwrap();
    assert_eq!(limited["scoring_status"], "not_run_graph_search_limit");
    assert!(limited["common_fragments"].is_null());
    assert!(limited["best_candidate"].is_null());
    assert!(limited["candidates"].as_array().unwrap().len() > 1);
    assert!(limited["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .all(|c| c["score"].is_null() && c["junction_fragments"].is_null()));

    let (complete, with_layout, guided) = assemble(&dir, "guided", &["--coverage-layout"]);
    assert!(!complete);
    assert_eq!(strict, with_layout);
    assert!(!dir.path("guided/consensus.fasta").exists());
    let layout = &guided["coverage_layout"];
    assert_eq!(layout["status"], "coverage_inferred_layout", "{guided:#}");
    assert_eq!(layout["final_verified"], false);
    let inferred = load_reference(&dir.path("guided/coverage_layout.fasta"), None).unwrap();
    assert_eq!(inferred.sequence.len(), sample.len());
    // Total repeat use is identifiable from these uniform reads, but their
    // placement across two loci is not. Retain alternatives and no final.
    assert!(
        layout["joint_copy_model"]["equally_scored_retained_configurations"]
            .as_u64()
            .unwrap()
            > 1
    );
    assert_eq!(
        layout["joint_copy_model"]["cross_window_phase_verified"],
        false
    );
    let mut alternatives = noodles_fasta::io::Reader::new(std::io::BufReader::new(
        fs::File::open(dir.path("guided/coverage_layout_candidates.fasta")).unwrap(),
    ));
    assert!(alternatives
        .records()
        .any(|r| circular_eq(r.unwrap().sequence().as_ref(), &sample)));
    let window = guided["windows"]
        .as_array()
        .unwrap()
        .iter()
        .find(|w| w["coverage_layout"]["status"] == "coverage_inferred_layout")
        .unwrap();
    assert_eq!(window["status"], "unresolved");
    assert_eq!(
        window["coverage_layout"]["paired_comparison"]["scoring_status"],
        "evaluated"
    );
    assert!(
        window["coverage_layout"]["paired_comparison"]["common_fragments"]
            .as_u64()
            .unwrap()
            > 0
    );
    for row in fs::read_to_string(dir.path("guided/coverage_layout_ambiguity.tsv"))
        .unwrap()
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = row.split('\t').collect();
        let (a, b): (usize, usize) = (fields[1].parse().unwrap(), fields[2].parse().unwrap());
        assert!(a >= 1 && a <= b && b <= inferred.sequence.len());
    }
    let (complete, multi_sequence, multi) =
        assemble(&dir, "multik", &["--coverage-layout", "--multi-k-layout"]);
    assert!(!complete);
    assert_eq!(strict, multi_sequence);
    assert!(!dir.path("multik/consensus.fasta").exists());
    let cross: Vec<_> = multi["windows"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|w| w["coverage_layout"]["searches"].as_array())
        .flatten()
        .filter(|s| s["method"] == "long_flanks_cross_k_read_context")
        .collect();
    assert!(cross.iter().any(|s| s["k"] == 95));
    assert!(cross.iter().all(|s| s["search_complete"] == false));
}

#[test]
fn absent_coverage_is_not_a_deletion_or_reference_filled_sequence() {
    let dir = Fixture::new();
    let reference = dna(3000, 892414);
    inputs_with_gap(
        &dir,
        &reference,
        &reference,
        &(0..3000).collect::<Vec<_>>(),
        false,
        false,
        Some(1750..2050),
    );
    let (complete, sequence, report) = assemble(&dir, "local", &[]);
    assert!(!complete, "{report:#}");
    assert_eq!(sequence.len(), reference.len());
    assert!(sequence[1750..2050].iter().all(|b| *b == b'N'));
    assert_eq!(report["applied_replacements"], 0);
    let (complete, forced, report) = assemble(&dir, "guided", &["--coverage-layout"]);
    assert!(!complete);
    assert_eq!(forced, sequence);
    assert_eq!(report["coverage_layout"]["status"], "no_inferred_windows");
    assert!(!dir.path("guided/coverage_layout.fasta").exists());
    let (complete, legacy, report) = assemble(&dir, "legacy_missing", &["--best-effort"]);
    assert!(complete); // Best-effort output availability, not structural validation.
    assert_eq!(legacy, sequence);
    assert!(dir.path("legacy_missing/consensus.fasta").exists());
    let root: Value =
        serde_json::from_slice(&fs::read(dir.path("legacy_missing/report.json")).unwrap()).unwrap();
    assert_eq!(root["status"], "best_effort");
    assert_eq!(root["strict_validation_passed"], false);
    // At small k, incidental matches can create a graph route despite a true
    // coverage gap. Such a route is an explicit counterexample to N-free = solved.
    assert_eq!(report["legacy_bridge"]["final_verified"], false);
    let false_bridge = load_reference(&dir.path("legacy_missing/legacy_bridge.fasta"), None)
        .unwrap()
        .sequence;
    let mut doubled_reference = reference.clone();
    doubled_reference.extend_from_slice(&reference);
    assert!(!doubled_reference
        .windows(false_bridge.len())
        .any(|s| s == false_bridge));
    for window in report["windows"].as_array().unwrap() {
        if window["legacy_bridge"]["status"] == "inferred_bridge" {
            let comparison = &window["legacy_bridge"]["paired_comparison"];
            assert_eq!(comparison["status"], "diagnostic_only");
            assert_eq!(comparison["selected_candidate"], Value::Null);
            assert_eq!(
                comparison["within_shortlist_status"], "unresolved",
                "{window:#}"
            );
        }
    }
}

#[test]
fn best_effort_writes_consensus_after_strict_search_limits() {
    let dir = Fixture::new();
    let mut reference = dna(1400, 729381);
    let unit = dna(140, 915391);
    reference[200..340].copy_from_slice(&unit);
    let insert = unit.repeat(2);
    let sample = [&reference[..1000], insert.as_slice(), &reference[1000..]].concat();
    let coordinates: Vec<_> = (0..1000)
        .chain(200..340)
        .chain(200..340)
        .chain(1000..1400)
        .collect();
    inputs(&dir, &reference, &sample, &coordinates, false, false);
    let (_, strict, _) = assemble(&dir, "before_legacy", &["--local-max-paths", "2"]);
    let (complete, with_legacy, report) =
        assemble(&dir, "legacy", &["--best-effort", "--local-max-paths", "2"]);
    assert!(complete, "{report:#}");
    assert_eq!(with_legacy, strict);
    assert!(dir.path("legacy/consensus.fasta").exists());
    let root: Value =
        serde_json::from_slice(&fs::read(dir.path("legacy/report.json")).unwrap()).unwrap();
    assert_eq!(root["status"], "best_effort");
    assert_eq!(root["final_published"], false);
    assert_eq!(root["consensus_written"], true);
    assert_eq!(report["legacy_bridge"]["final_verified"], false);
    assert!(dir.path("legacy/legacy_bridge.fasta").exists());
    let window = report["windows"]
        .as_array()
        .unwrap()
        .iter()
        .find(|w| {
            w["reason"] == "graph_search_limit" && w["legacy_bridge"]["status"] == "inferred_bridge"
        })
        .unwrap();
    let inferred = &window["legacy_bridge"];
    assert_eq!(inferred["paired_comparison"]["scoring_status"], "evaluated");
    assert!(
        inferred["paired_comparison"]["common_fragments"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert_eq!(
        inferred["paired_comparison"]["selected_candidate"],
        Value::Null
    );
    assert_eq!(inferred["search_complete"], false);
    assert!(inferred["searches"]
        .as_array()
        .unwrap()
        .iter()
        .any(|s| s["k"] == 11));
    let sequence = load_reference(&dir.path("legacy/legacy_bridge.fasta"), None)
        .unwrap()
        .sequence;
    for line in fs::read_to_string(dir.path("legacy/legacy_bridge.ambiguity.tsv"))
        .unwrap()
        .lines()
        .skip(1)
    {
        let c: Vec<_> = line.split('\t').collect();
        let start: usize = c[1].parse().unwrap();
        let end: usize = c[2].parse().unwrap();
        assert!(start > 0 && end >= start && end <= sequence.len());
    }
}

#[test]
fn maximum_consensus_fills_six_repeat_copies_at_existing_span_despite_incomplete_strict_search() {
    let dir = Fixture::new();
    let donor = Fixture::new();
    let mut truth = dna(2000, 881392);
    truth[900..1122].copy_from_slice(&dna(37, 711).repeat(6));
    let coordinates: Vec<_> = (0..truth.len()).collect();
    inputs(&donor, &truth, &truth, &coordinates, false, false);
    inputs_with_gap(
        &dir,
        &truth,
        &truth,
        &coordinates,
        false,
        false,
        Some(900..1122),
    );
    // The BAM cannot cover this region; the supplied sample FASTQ can.
    for mate in [1, 2] {
        fs::copy(
            donor.path(&format!("r{mate}.fq")),
            dir.path(&format!("r{mate}.fq")),
        )
        .unwrap();
    }
    let (available, strict, report) = assemble(
        &dir,
        "maximum",
        &["--best-effort", "--local-max-paths", "2"],
    );
    assert!(available);
    assert!(strict.contains(&b'N'), "{report}");
    let maximum = load_reference(&dir.path("maximum/consensus.fasta"), None)
        .unwrap()
        .sequence;
    assert!(circular_eq(&maximum, &truth));
    let metadata: Value =
        serde_json::from_slice(&fs::read(dir.path("maximum/report.json")).unwrap()).unwrap();
    assert_eq!(metadata["ambiguous_bases"], 0);
    assert_eq!(metadata["sequence_output"]["fasta"], "consensus.fasta");
    assert_eq!(metadata["consensus_written"], true);
    assert_eq!(metadata["strict_validation_passed"], false);
    assert!(
        fs::read_to_string(dir.path("maximum/maximum_consensus.ambiguity.tsv"))
            .unwrap()
            .lines()
            .count()
            > 1
    );
    let (available, fast_strict, fast_report) = assemble_mode(&dir, "fast_maximum", &[], false);
    assert!(available);
    assert_eq!(fast_strict, strict);
    let fast = load_reference(&dir.path("fast_maximum/consensus.fasta"), None)
        .unwrap()
        .sequence;
    assert!(circular_eq(&fast, &truth));
    assert_eq!(fast_report["validation_mode"], "fast_consensus");
    assert!(fast_report["windows"]
        .as_array()
        .unwrap()
        .iter()
        .any(|w| w["legacy_bridge"]["status"] == "inferred_fill"));
    for window in fast_report["windows"].as_array().unwrap() {
        assert_eq!(window["scoring_status"], "not_run_fast_consensus");
        assert_eq!(window["common_fragments"], Value::Null);
    }
}

#[test]
fn fast_local_consensus_keeps_all_pairs_and_defers_structural_replacement() {
    let (dir, reference, sample) = insertion(2100, false);
    // A complete FASTQ pair absent from the BAM must still enter the graph.
    for mate in [1, 2] {
        let mut file = fs::OpenOptions::new()
            .append(true)
            .open(dir.path(&format!("r{mate}.fq")))
            .unwrap();
        writeln!(
            file,
            "@extra/{mate}\n{}\n+\n{}",
            String::from_utf8_lossy(&reference[400..500]),
            "I".repeat(100)
        )
        .unwrap();
    }
    let (available, candidate, report) = assemble_mode(&dir, "fast", &[], false);
    assert!(available);
    assert_eq!(candidate, reference);
    assert_eq!(report["named_pairs"], sample.len() / 2 + 1);
    assert_eq!(report["applied_replacements"], 0);
    assert_eq!(report["validation_mode"], "fast_consensus");
    assert_eq!(report["fragment_model"]["status"], "not_run_fast_consensus");
    let windows = report["windows"].as_array().unwrap();
    assert!(!windows.is_empty());
    for window in windows {
        assert_eq!(window["status"], "unresolved");
        assert_eq!(window["scoring_status"], "not_run_fast_consensus");
        assert_eq!(window["common_fragments"], Value::Null);
        assert_eq!(window["unscored_fragments"], Value::Null);
        assert_eq!(window["best_candidate"], Value::Null);
        assert_eq!(window["search"]["states"], Value::Null);
        for c in window["candidates"].as_array().into_iter().flatten() {
            assert_eq!(c["score"], Value::Null);
        }
    }
    assert!(!fs::read_dir(dir.path("fast/local_assembly"))
        .unwrap()
        .any(|p| p.unwrap().path().extension().is_some_and(|s| s == "gz")));
    let sequence = load_reference(&dir.path("fast/consensus.fasta"), None)
        .unwrap()
        .sequence;
    assert!(circular_eq(&sequence, &reference));
    let uncertainty = fs::read_to_string(dir.path("fast/maximum_consensus.ambiguity.tsv")).unwrap();
    assert!(uncertainty.contains("structural_validation_deferred"));
    for row in uncertainty.lines().skip(1) {
        let fields: Vec<_> = row.split('\t').collect();
        assert_eq!(fields[0], "MtDNA");
        let start = fields[1].parse::<usize>().unwrap();
        let end = fields[2].parse::<usize>().unwrap();
        assert!(start > 0 && end >= start && end <= sequence.len());
    }
    let root: Value =
        serde_json::from_slice(&fs::read(dir.path("fast/report.json")).unwrap()).unwrap();
    assert_eq!(root["status"], "best_effort");
    assert_eq!(root["validation_mode"], "fast_consensus");
    assert_eq!(root["strict_validation_passed"], false);
    assert_eq!(root["circularity_verified"], false);
    assert_eq!(root["subsampling"]["enabled"], false);
}

#[test]
fn best_effort_retains_a_complete_consensus_when_no_fill_is_needed() {
    let dir = Fixture::new();
    let truth = dna(1200, 798341);
    inputs(
        &dir,
        &truth,
        &truth,
        &(0..truth.len()).collect::<Vec<_>>(),
        false,
        false,
    );
    let (available, _, _) = assemble(&dir, "complete_maximum", &["--best-effort"]);
    assert!(available);
    let text = fs::read_to_string(dir.path("complete_maximum/consensus.fasta")).unwrap();
    assert!(text.lines().next().unwrap().contains("final_verified=true"));
    let seq = load_reference(&dir.path("complete_maximum/consensus.fasta"), None)
        .unwrap()
        .sequence;
    assert!(circular_eq(&seq, &truth));
    let report: Value =
        serde_json::from_slice(&fs::read(dir.path("complete_maximum/report.json")).unwrap())
            .unwrap();
    assert_eq!(report["strict_validation_passed"], true);
    assert_eq!(report["consensus_written"], true);
    let (available, _, _) = assemble_mode(&dir, "fast_complete", &[], false);
    assert!(available);
    let text = fs::read_to_string(dir.path("fast_complete/consensus.fasta")).unwrap();
    assert!(text
        .lines()
        .next()
        .unwrap()
        .contains("final_verified=false"));
    let report: Value =
        serde_json::from_slice(&fs::read(dir.path("fast_complete/report.json")).unwrap()).unwrap();
    assert_eq!(report["strict_validation_passed"], false);
    assert_eq!(report["validation_mode"], "fast_consensus");
}

#[test]
fn mapping_command_and_best_effort_options_use_current_names() {
    let top_help = Cli::command().render_long_help().to_string();
    for command in ["assemble-baited", "assemble-seeded", "assemble-mapping"] {
        assert!(top_help.contains(command));
    }
    assert!(!top_help.contains("assemble-denovo"));
    assert!(!top_help.contains("polish"));
    for command in ["assemble-mapping", "polish"] {
        let base = [
            "mitovarr",
            command,
            "--ref",
            "ref.fa",
            "--bam",
            "x.bam",
            "--out",
            "out",
            "--local-assembly",
            "--r1",
            "a.fq",
            "--r2",
            "b.fq",
        ];
        let Commands::AssembleMapping(default) = Cli::try_parse_from(base).unwrap().command else {
            panic!("expected mapping command");
        };
        assert!(!default.best_effort);
        assert!(!default.deep_validation);
        assert!(Cli::try_parse_from(base.into_iter().chain(["--coverage-layout"])).is_err());
        let Commands::AssembleMapping(deep) = Cli::try_parse_from(
            base.into_iter()
                .chain(["--deep-validation", "--coverage-layout"]),
        )
        .unwrap()
        .command
        else {
            panic!("expected mapping command")
        };
        assert!(deep.deep_validation && deep.coverage_layout);
        for option in ["--best-effort", "--legacy-bridge"] {
            let Commands::AssembleMapping(args) =
                Cli::try_parse_from(base.into_iter().chain([option]))
                    .unwrap()
                    .command
            else {
                panic!("expected mapping command");
            };
            assert!(args.best_effort);
        }
    }
    let help = Cli::try_parse_from(["mitovarr", "assemble-mapping", "--help"])
        .unwrap_err()
        .to_string();
    assert!(help.contains("--best-effort"));
    assert!(help.contains("--deep-validation"));
    assert!(!help.contains("--legacy-bridge"));
    assert!(help.contains("Usage: mitovarr assemble-mapping"));
    assert!(!help.contains("polish"));
    assert!(Cli::try_parse_from([
        "mitovarr",
        "assemble-mapping",
        "--ref",
        "ref.fa",
        "--bam",
        "x.bam",
        "--out",
        "out",
        "--best-effort"
    ])
    .is_err());
    assert!(Cli::try_parse_from([
        "mitovarr",
        "assemble-baited",
        "--r1",
        "a.fq",
        "--r2",
        "b.fq",
        "--out",
        "out",
        "--best-effort"
    ])
    .is_err());
}

#[test]
fn repeat_ambiguity_at_origin_remains_unresolved_without_coverage() {
    let dir = Fixture::new();
    let reference = dna(1000, 892414);
    inputs_with_gap(
        &dir,
        &reference,
        &reference,
        &(0..reference.len()).collect::<Vec<_>>(),
        false,
        false,
        Some(0..300),
    );
    let (complete, sequence, report) = assemble(&dir, "local", &[]);
    assert!(!complete, "{report:#}");
    assert_eq!(sequence.len(), reference.len());
    assert!(sequence[..300].iter().all(|&b| b == b'N'));
    let root_report: Value =
        serde_json::from_slice(&fs::read(dir.path("local/report.json")).unwrap()).unwrap();
    assert_eq!(
        root_report["sequence_output"]["status"],
        "partial_available"
    );
    assert_eq!(root_report["final_published"], false);
    assert_eq!(report["applied_replacements"], 0);
    assert!(report["remaining_repeat_intervals"].as_u64().unwrap() > 0);
    assert!(!dir.path("local/consensus.fasta").exists());
    let ambiguity = fs::read_to_string(dir.path("local/ambiguity.tsv")).unwrap();
    assert!(ambiguity.lines().skip(1).any(|row| row.starts_with("1\t")));
    let (available, fast_candidate, _) = assemble_mode(&dir, "fast_missing", &[], false);
    assert!(available);
    assert_eq!(fast_candidate, sequence);
    let fast = load_reference(&dir.path("fast_missing/consensus.fasta"), None)
        .unwrap()
        .sequence;
    assert!(fast.contains(&b'N'));
    let root: Value =
        serde_json::from_slice(&fs::read(dir.path("fast_missing/report.json")).unwrap()).unwrap();
    assert_eq!(root["strict_validation_passed"], false);
}

#[test]
fn explicit_subsampling_caps_bam_and_keeps_matching_complete_pairs() {
    let dir = Fixture::new();
    let genome = dna(600, 71923);
    let coordinates: Vec<_> = (0..genome.len()).collect();
    inputs(&dir, &genome, &genome, &coordinates, false, false);
    let (_, first, _) = assemble(&dir, "sampled", &["--subsample", "--max-records", "300"]);
    let (_, second, _) = assemble(&dir, "repeat", &["--subsample", "--max-records", "300"]);
    assert_eq!(first, second);
    let report: Value =
        serde_json::from_slice(&fs::read(dir.path("sampled/report.json")).unwrap()).unwrap();
    let sample = &report["subsampling"];
    assert_eq!(sample["input_target_records"], 600);
    assert_eq!(sample["retained_target_records"], 300);
    assert_eq!(sample["retained_fastq_pairs"], 150);
    assert_eq!(report["records"], 300);
    let scan = mitovarr::io::fastq::scan_pairs(
        &dir.path("sampled/subsampled.R1.fastq"),
        &dir.path("sampled/subsampled.R2.fastq"),
        |_, _, _| Ok(()),
    )
    .unwrap();
    assert_eq!(scan.pairs, 150);
    let (complete, _, _) = assemble(&dir, "low", &["--subsample", "--max-records", "2"]);
    assert!(!complete);
    assert!(!dir.path("low/consensus.fasta").exists());
}
