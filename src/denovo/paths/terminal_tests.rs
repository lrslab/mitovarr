use super::*;
use crate::denovo::threading::Segment;
use std::{fs, io::Write, path::PathBuf};

struct Fixture(PathBuf);
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn dna(n: usize, mut seed: u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed >> 32) as usize & 3]
        })
        .collect()
}

#[test]
fn short_reads_resolve_a_single_snp_without_a_fragment_distance_model() {
    let dir =
        Fixture(std::env::temp_dir().join(format!("mitovar-short-snp-{}", std::process::id())));
    fs::create_dir_all(&dir.0).unwrap();
    let left = dna(600, 172311);
    let right = dna(600, 318177);
    let expected = [left.clone(), vec![b'C'], right.clone()].concat();
    for k in [31, 47, 63] {
        let overlap = k - 1;
        let branch = |base| {
            [
                left[600 - overlap..].to_vec(),
                vec![base],
                right[..overlap].to_vec(),
            ]
            .concat()
        };
        let forward = [left.clone(), branch(b'C'), branch(b'A'), right.clone()];
        let graph = PathGraph {
            segments: forward
                .iter()
                .cloned()
                .chain(forward.iter().map(|s| kmer::reverse_complement(s)))
                .map(|sequence| Segment {
                    sequence,
                    circular: false,
                })
                .collect(),
            outgoing: vec![
                vec![1, 2],
                vec![3],
                vec![3],
                vec![],
                vec![],
                vec![4],
                vec![4],
                vec![5, 6],
            ],
        };
        let settings = GraphSettings {
            coverage_layout: false,
            multi_k_layout: false,
            k,
            min_kmer_count: 3,
            max_error_fraction: 0.2,
            max_bubble_fraction: 0.49,
            min_depth: 10,
            min_junction_fragments: 3,
            max_validation_pairs: 1000,
        };
        let limits = Limits {
            max_kmers: 100_000,
            min_baseq: 20,
        };
        for covered in [true, false] {
            let out = dir.0.join(format!("k{k}-{covered}"));
            fs::create_dir(&out).unwrap();
            let reads = Reads {
                r1: out.join("r1.fq"),
                r2: out.join("r2.fq"),
            };
            let mut r1 = fs::File::create(&reads.r1).unwrap();
            let mut r2 = fs::File::create(&reads.r2).unwrap();
            for (allele, copies) in [(b'C', 3), (b'A', 1)] {
                let genome = [left.clone(), vec![allele], right.clone()].concat();
                for start in 560..570 {
                    for copy in 0..copies {
                        let a = if covered {
                            &genome[start..start + 100]
                        } else {
                            &left[500..]
                        };
                        let b = if covered { a } else { &right[..100] };
                        for (mate, writer, sequence) in [
                            (1, &mut r1, a.to_vec()),
                            (2, &mut r2, kmer::reverse_complement(b)),
                        ] {
                            writeln!(
                                writer,
                                "@snp_{allele}_{start}_{copy}/{mate}\n{}\n+\n{}",
                                String::from_utf8(sequence).unwrap(),
                                "I".repeat(100)
                            )
                            .unwrap();
                        }
                    }
                }
            }
            drop(r1);
            drop(r2);
            let result = resolve(&graph, &reads, &out, &settings, &limits).unwrap();
            assert_eq!(result.report["fragment_model"]["status"], "not_estimable");
            assert_eq!(
                result
                    .primary_candidates
                    .iter()
                    .any(|p| p.candidate.sequence == normalize(&expected, false)),
                covered,
                "k={k}, covered={covered}: {}",
                result.report
            );
            assert!(result.supported_cycles.is_empty());
        }
    }
}

#[test]
fn multi_k_end_candidates_join_sample_sequence_but_leave_a_real_read_gap_open() {
    let dir = Fixture(
        std::env::temp_dir().join(format!("mitovar-end-reassembly-{}", std::process::id())),
    );
    fs::create_dir_all(&dir.0).unwrap();
    let genome = dna(1240, 718373);
    let forward = [genome[..600].to_vec(), genome[640..].to_vec()];
    let graph = PathGraph {
        segments: forward
            .iter()
            .cloned()
            .chain(forward.iter().map(|s| kmer::reverse_complement(s)))
            .map(|sequence| Segment {
                sequence,
                circular: false,
            })
            .collect(),
        outgoing: vec![vec![]; 4],
    };
    let settings = GraphSettings {
        coverage_layout: false,
        multi_k_layout: false,
        k: 63,
        min_kmer_count: 3,
        max_error_fraction: 0.2,
        max_bubble_fraction: 0.49,
        min_depth: 10,
        min_junction_fragments: 3,
        max_validation_pairs: 10_000,
    };
    let limits = Limits {
        max_kmers: 100_000,
        min_baseq: 20,
    };
    for covered in [false, true] {
        let out = dir.0.join(if covered { "covered" } else { "missing" });
        fs::create_dir(&out).unwrap();
        let reads = Reads {
            r1: out.join("r1.fq"),
            r2: out.join("r2.fq"),
        };
        let mut r1 = fs::File::create(&reads.r1).unwrap();
        let mut r2 = fs::File::create(&reads.r2).unwrap();
        for start in 0..=genome.len() - 220 {
            if !covered && start < 640 && start + 220 > 600 {
                continue;
            }
            for (mate, writer, sequence) in [
                (1, &mut r1, genome[start..start + 140].to_vec()),
                (
                    2,
                    &mut r2,
                    kmer::reverse_complement(&genome[start + 80..start + 220]),
                ),
            ] {
                writeln!(
                    writer,
                    "@f{start}/{mate}\n{}\n+\n{}",
                    String::from_utf8(sequence).unwrap(),
                    "I".repeat(140)
                )
                .unwrap();
            }
        }
        drop(r1);
        drop(r2);
        let result = resolve(&graph, &reads, &out, &settings, &limits).unwrap();
        assert_eq!(
            result.report["end_assembly"]["directed_connections_added"],
            if covered { 2 } else { 0 },
            "{}",
            result.report
        );
        assert_eq!(
            result
                .primary_candidates
                .iter()
                .any(|p| p.candidate.sequence == normalize(&genome, false)),
            covered,
            "{}",
            result.report
        );
        assert!(result.candidates.iter().all(|c| !c.circular));
        assert!(out.join("path_graph.gfa").is_file());
    }
}

#[test]
fn limited_end_search_preserves_proposals_without_closing_the_graph() {
    let dir =
        Fixture(std::env::temp_dir().join(format!("mitovar-limited-end-{}", std::process::id())));
    fs::create_dir_all(&dir.0).unwrap();
    let backbone = dna(800, 417833);
    // Reads cover the gap, but its repeat permits an unsearched copy-number
    // tail. Discovered circle proposals must not become strict graph joins.
    let genome = [backbone.clone(), vec![b'A'; 70]].concat();
    let extended = [genome.clone(), genome[..220].to_vec()].concat();
    let graph = PathGraph {
        segments: [backbone.clone(), kmer::reverse_complement(&backbone)]
            .into_iter()
            .map(|sequence| Segment {
                sequence,
                circular: false,
            })
            .collect(),
        outgoing: vec![vec![]; 2],
    };
    let settings = GraphSettings {
        coverage_layout: false,
        multi_k_layout: false,
        k: 63,
        min_kmer_count: 3,
        max_error_fraction: 0.2,
        max_bubble_fraction: 0.49,
        min_depth: 10,
        min_junction_fragments: 3,
        max_validation_pairs: 10_000,
    };
    let limits = Limits {
        max_kmers: 100_000,
        min_baseq: 20,
    };
    let reads = Reads {
        r1: dir.0.join("r1.fq"),
        r2: dir.0.join("r2.fq"),
    };
    let mut r1 = fs::File::create(&reads.r1).unwrap();
    let mut r2 = fs::File::create(&reads.r2).unwrap();
    for start in 0..genome.len() {
        for (mate, writer, sequence) in [
            (1, &mut r1, extended[start..start + 140].to_vec()),
            (
                2,
                &mut r2,
                kmer::reverse_complement(&extended[start + 80..start + 220]),
            ),
        ] {
            writeln!(
                writer,
                "@f{start}/{mate}\n{}\n+\n{}",
                String::from_utf8(sequence).unwrap(),
                "I".repeat(140)
            )
            .unwrap();
        }
    }
    drop(r1);
    drop(r2);
    let result = resolve(&graph, &reads, &dir.0, &settings, &limits).unwrap();
    let ends = &result.report["end_assembly"];
    assert!(!ends["limited_sources"].as_array().unwrap().is_empty());
    let proposals = fs::read_to_string(dir.0.join("end_bridge_proposals.fasta")).unwrap();
    assert!(proposals.starts_with('>'));
    let proposed = ends["proposals"].as_array().unwrap();
    assert!(proposed.iter().any(|p| p["role"] == "backbone_closure"));
    assert!(proposed.iter().all(|p| p["join_eligible"] == false));
    assert_eq!(ends["directed_connections_added"], 0, "{ends}");
    assert!(
        fs::read_to_string(dir.0.join("end_bridge_candidates.fasta"))
            .unwrap()
            .is_empty()
    );
    assert!(!fs::read_to_string(dir.0.join("path_graph.gfa"))
        .unwrap()
        .lines()
        .any(|line| line.starts_with("L\t")));
    assert!(result.supported_cycles.is_empty());
    assert!(result.candidates.iter().all(|c| !c.circular));
    assert!(result
        .primary_candidates
        .iter()
        .all(|p| !p.candidate.circular));
}

#[test]
fn either_mate_can_anchor_a_fragment_without_counting_it_twice() {
    let dir =
        Fixture(std::env::temp_dir().join(format!("mitovar-either-anchor-{}", std::process::id())));
    fs::create_dir_all(&dir.0).unwrap();
    let genome = dna(770, 937117);
    let forward = [genome[..600].to_vec(), genome[538..].to_vec()];
    let graph = PathGraph {
        segments: forward
            .iter()
            .cloned()
            .chain(forward.iter().map(|s| kmer::reverse_complement(s)))
            .map(|sequence| Segment {
                sequence,
                circular: false,
            })
            .collect(),
        outgoing: vec![vec![1], vec![], vec![], vec![2]],
    };
    let settings = GraphSettings {
        coverage_layout: false,
        multi_k_layout: false,
        k: 63,
        min_kmer_count: 3,
        max_error_fraction: 0.2,
        max_bubble_fraction: 0.49,
        min_depth: 10,
        min_junction_fragments: 3,
        max_validation_pairs: 1000,
    };
    let limits = Limits {
        max_kmers: 100_000,
        min_baseq: 20,
    };
    for mode in ["r1", "r2", "both"] {
        let out = dir.0.join(mode);
        fs::create_dir(&out).unwrap();
        let reads = Reads {
            r1: out.join("r1.fq"),
            r2: out.join("r2.fq"),
        };
        let mut r1 = fs::File::create(&reads.r1).unwrap();
        let mut r2 = fs::File::create(&reads.r2).unwrap();
        for start in 500..506 {
            let sequence = &genome[start..start + 150];
            for (mate, writer, sequence) in [
                (
                    1,
                    &mut r1,
                    if mode == "r2" {
                        vec![b'N'; 150]
                    } else {
                        sequence.to_vec()
                    },
                ),
                (
                    2,
                    &mut r2,
                    if mode == "r1" {
                        vec![b'N'; 150]
                    } else {
                        kmer::reverse_complement(sequence)
                    },
                ),
            ] {
                writeln!(
                    writer,
                    "@fragment{start}/{mate}\n{}\n+\n{}",
                    String::from_utf8(sequence).unwrap(),
                    "I".repeat(150)
                )
                .unwrap();
            }
        }
        drop(r1);
        drop(r2);
        let result = resolve(&graph, &reads, &out, &settings, &limits).unwrap();
        assert!(
            result
                .candidates
                .iter()
                .any(|c| c.sequence == normalize(&genome, false)),
            "mode={mode}: {}",
            result.report
        );
        assert_eq!(result.report["distinct_scored_fragment_names"], 6);
        for locus in result.report["loci"].as_array().unwrap() {
            if locus["paths"].as_array().unwrap().is_empty() {
                continue;
            }
            assert_eq!(locus["fragments"], 6, "mode={mode}");
            assert_eq!(
                locus["r2_anchored_fragments"],
                if mode == "r2" { 6 } else { 0 }
            );
            assert_eq!(locus["r2_only_fragments"], if mode == "r2" { 6 } else { 0 });
        }
    }
}

#[test]
fn short_terminal_gene_tail_joins_only_with_spanning_fragments() {
    let dir = Fixture(
        std::env::temp_dir().join(format!("mitovar-short-terminal-{}", std::process::id())),
    );
    fs::create_dir_all(&dir.0).unwrap();
    let main = dna(600, 319771);
    let tail = [main[538..].to_vec(), dna(170, 276443)].concat();
    // A short competing incoming tip keeps the tail separate during compaction.
    // Neither a reference nor a gene model is passed to the path resolver.
    let tip = [dna(15, 617411), main[538..].to_vec()].concat();
    let genome = [main.clone(), tail[62..].to_vec()].concat();
    let graph = PathGraph {
        segments: [main.clone(), tail.clone(), tip.clone()]
            .into_iter()
            .chain(
                [main, tail, tip]
                    .iter()
                    .map(|s| kmer::reverse_complement(s)),
            )
            .map(|sequence| Segment {
                sequence,
                circular: false,
            })
            .collect(),
        outgoing: vec![vec![1], vec![], vec![1], vec![], vec![3, 5], vec![]],
    };
    let settings = GraphSettings {
        coverage_layout: false,
        multi_k_layout: false,
        k: 63,
        min_kmer_count: 3,
        max_error_fraction: 0.2,
        max_bubble_fraction: 0.49,
        min_depth: 10,
        min_junction_fragments: 3,
        max_validation_pairs: 10_000,
    };
    let limits = Limits {
        max_kmers: 100_000,
        min_baseq: 20,
    };
    for bridge in [false, true] {
        let out = dir.0.join(if bridge { "spanned" } else { "unspanned" });
        fs::create_dir(&out).unwrap();
        let reads = Reads {
            r1: out.join("r1.fq"),
            r2: out.join("r2.fq"),
        };
        let mut r1 = fs::File::create(&reads.r1).unwrap();
        let mut r2 = fs::File::create(&reads.r2).unwrap();
        for start in 0..=genome.len() - 160 {
            if !bridge && start + 160 > 538 && start < 600 {
                continue;
            }
            for (mate, writer, sequence) in [
                (1, &mut r1, genome[start..start + 90].to_vec()),
                (
                    2,
                    &mut r2,
                    kmer::reverse_complement(&genome[start + 70..start + 160]),
                ),
            ] {
                writeln!(
                    writer,
                    "@fragment{start}/{mate}\n{}\n+\n{}",
                    String::from_utf8(sequence).unwrap(),
                    "I".repeat(90)
                )
                .unwrap();
            }
        }
        drop(r1);
        drop(r2);
        let result = resolve(&graph, &reads, &out, &settings, &limits).unwrap();
        assert_eq!(result.report["short_terminal_anchor_count"], 2);
        assert_eq!(
            result
                .candidates
                .iter()
                .any(|c| c.sequence == normalize(&genome, false)),
            bridge,
            "{:#}",
            result.report
        );
        assert!(result.supported_cycles.is_empty());
        assert!(result.candidates.iter().all(|c| !c.circular));
        assert_eq!(
            result
                .primary_candidates
                .iter()
                .any(|p| p.candidate.sequence == normalize(&genome, false)),
            bridge
        );
    }
}

#[test]
fn spanning_majority_keeps_a_primary_contig_across_a_real_minor_indel() {
    let dir = Fixture(
        std::env::temp_dir().join(format!("mitovar-primary-bubble-{}", std::process::id())),
    );
    fs::create_dir_all(&dir.0).unwrap();
    let left = dna(600, 717233);
    let right = dna(600, 113833);
    let major = dna(40, 277311);
    let mut minor = major.clone();
    minor.splice(20..20, [b'A', b'C']);
    let inner = |s: &[u8]| [left[538..].to_vec(), s.to_vec(), right[..62].to_vec()].concat();
    let forward = [left.clone(), inner(&major), inner(&minor), right.clone()];
    let graph = PathGraph {
        segments: forward
            .iter()
            .cloned()
            .chain(forward.iter().map(|s| kmer::reverse_complement(s)))
            .map(|sequence| Segment {
                sequence,
                circular: false,
            })
            .collect(),
        outgoing: vec![
            vec![1, 2],
            vec![3],
            vec![3],
            vec![],
            vec![],
            vec![4],
            vec![4],
            vec![5, 6],
        ],
    };
    let settings = GraphSettings {
        coverage_layout: false,
        multi_k_layout: false,
        k: 63,
        min_kmer_count: 3,
        max_error_fraction: 0.2,
        max_bubble_fraction: 0.49,
        min_depth: 10,
        min_junction_fragments: 3,
        max_validation_pairs: 10_000,
    };
    let limits = Limits {
        max_kmers: 100_000,
        min_baseq: 20,
    };
    let reads = Reads {
        r1: dir.0.join("r1.fq"),
        r2: dir.0.join("r2.fq"),
    };
    let mut r1 = fs::File::create(&reads.r1).unwrap();
    let mut r2 = fs::File::create(&reads.r2).unwrap();
    let expected = [left.clone(), major.clone(), right.clone()].concat();
    for (allele, core, copies) in [(0, &major, 3), (1, &minor, 1)] {
        let genome = [left.clone(), core.clone(), right.clone()].concat();
        for start in (0..=genome.len() - 300).step_by(3) {
            for copy in 0..copies {
                for (mate, writer, sequence) in [
                    (1, &mut r1, genome[start..start + 125].to_vec()),
                    (
                        2,
                        &mut r2,
                        kmer::reverse_complement(&genome[start + 175..start + 300]),
                    ),
                ] {
                    writeln!(
                        writer,
                        "@f{allele}_{start}_{copy}/{mate}\n{}\n+\n{}",
                        String::from_utf8(sequence).unwrap(),
                        "I".repeat(125)
                    )
                    .unwrap();
                }
            }
        }
    }
    drop((r1, r2));
    let result = resolve(&graph, &reads, &dir.0, &settings, &limits).unwrap();
    let expected = normalize(&expected, false);
    assert!(
        !result.candidates.iter().any(|c| c.sequence == expected),
        "{:#}",
        result.report
    );
    let primary = result
        .primary_candidates
        .iter()
        .find(|p| p.candidate.sequence == expected)
        .unwrap_or_else(|| panic!("{:#}", result.report));
    assert!(!primary.candidate.circular);
    assert!(primary
        .regions
        .iter()
        .any(|r| r.reason == "majority_path_structure_unresolved" && r.alternatives == 1));
    assert!(result.supported_cycles.is_empty());
    let alternatives = fs::read_to_string(dir.0.join("path_candidates.fasta")).unwrap();
    assert_eq!(
        alternatives.lines().filter(|l| l.starts_with('>')).count(),
        4
    );
}
