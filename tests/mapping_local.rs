use clap::Parser;
use mitovarr::{cli::Cli, io::fasta::load_reference, kmer::reverse_complement, run};
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
    args.extend(extra.iter().map(|s| s.to_string()));
    let status = run(Cli::try_parse_from(args).unwrap()).unwrap();
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
    let (complete, _, report) = assemble(&dir, "local", &[]);
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
