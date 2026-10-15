use crate::{
    cli::{MappingArgs, Reads},
    consensus::{build_consensus_draft, finalize_consensus, ConsensusSettings},
    io::fasta::{load_reference, write_consensus_fasta},
    mapping,
    pileup::{FragmentPileupAccumulator, PileupCoordinateMode, PileupSettings},
    repeat::{self, RepeatSettings},
    report,
};
use anyhow::{ensure, Context, Result};
use noodles_bam as bam;
use noodles_core::Region;
use noodles_sam::alignment::RecordBuf;
use serde_json::json;
use std::{fs, io::Write, path::Path};

fn read_index(path: &Path) -> Result<bam::bai::Index> {
    let mut appended = path.as_os_str().to_owned();
    appended.push(".bai");
    let index_path = [appended.into(), path.with_extension("bai")]
        .into_iter()
        .find(|p: &std::path::PathBuf| p.is_file())
        .context("missing BAM index; provide .bam.bai or .bai")?;
    bam::bai::fs::read(&index_path)
        .with_context(|| format!("read BAM index {}", index_path.display()))
}

pub fn assemble(args: &MappingArgs) -> Result<bool> {
    // Fast local assembly always delivers the available sequence separately
    // from structural validation. Explicit deep mode retains the old gates.
    let mut effective = args.clone();
    effective.best_effort |= effective.local_assembly && !effective.deep_validation;
    let args = &effective;
    ensure!(
        args.min_depth > 0 && args.window > 0 && args.max_records > 0,
        "depth, window and record limit must be positive"
    );
    let min_snp_frac = args.min_alt_frac.unwrap_or(args.min_snp_frac);
    let min_indel_frac = args.min_alt_frac.unwrap_or(args.min_indel_frac);
    for (name, fraction) in [("SNP", min_snp_frac), ("indel", min_indel_frac)] {
        ensure!(
            fraction.is_finite() && (0.5..=1.0).contains(&fraction),
            "{name} fraction must be in 0.5..=1"
        );
    }
    let reference = load_reference(&args.reference, args.contig.as_deref())?;
    ensure!(
        !reference.sequence.is_empty() && reference.sequence.len() <= 1_000_000,
        "reference must contain 1..1,000,000 bases"
    );
    let mut reader = bam::io::indexed_reader::Builder::default()
        .set_index(read_index(&args.bam)?)
        .build_from_path(&args.bam)?;
    let header = reader.read_header()?;
    let (id, _, metadata) = header
        .reference_sequences()
        .get_full(reference.selected_name.as_bytes())
        .context("reference contig absent from BAM header")?;
    ensure!(
        usize::from(metadata.length()) == reference.sequence.len(),
        "reference and BAM contig lengths differ"
    );
    let mut accumulator = FragmentPileupAccumulator::new(
        &reference.sequence,
        PileupSettings {
            min_mapq: args.min_mapq,
            min_baseq: args.min_baseq,
        },
        PileupCoordinateMode::Circular,
    );
    let region = Region::new(reference.selected_name.as_bytes(), ..);
    let selection = if args.subsample {
        let mut selected = crate::subsample::Selection::new(args.max_records, args.subsample_seed);
        for record in reader.query(&header, &region)?.records() {
            let record = RecordBuf::try_from_alignment_record(&header, &record?)?;
            if record.reference_sequence_id() == Some(id) {
                selected.observe(
                    record
                        .name()
                        .context("subsampling requires BAM QNAME")?
                        .as_ref(),
                );
            }
        }
        ensure!(
            selected.retained > 0,
            "record cap retained no complete named fragments"
        );
        Some(selected)
    } else {
        None
    };
    let (sampled_args, sampling_report) = if let Some(selected) = &selection {
        let (sampled, report) = selected.paired_reads(args)?;
        (Some(sampled), report)
    } else {
        (None, json!({"enabled":false}))
    };
    let mut previous_start = 0;
    let mut records = 0;
    let mut signals = mapping::Signals::default();
    for record in reader.query(&header, &region)?.records() {
        let record = RecordBuf::try_from_alignment_record(&header, &record?)?;
        if record.reference_sequence_id() != Some(id) {
            continue;
        }
        if selection
            .as_ref()
            .is_some_and(|s| !s.contains(record.name().map_or(&[], |n| n.as_ref())))
        {
            continue;
        }
        records += 1;
        ensure!(
            records <= args.max_records,
            "BAM query exceeds --max-records; no final assembly published"
        );
        if let Some(start) = record.alignment_start() {
            let start = usize::from(start);
            ensure!(start >= previous_start, "BAM must be coordinate-sorted");
            previous_start = start;
        }
        if args.local_assembly {
            signals.observe(&record, reference.sequence.len(), args.min_mapq);
        }
        accumulator.push_record(&record);
    }
    let (pileup, stats) = accumulator.finish();
    let mut ambiguity = repeat::detect_ambiguity(
        &pileup,
        RepeatSettings {
            enabled: true,
            window: args.window,
            ambiguity_min_triggers: 2,
            repeat_min_mapq: args.min_mapq,
        },
    );
    ambiguity.sort_by_key(|a| (a.start, a.end));
    let mask = repeat::mask_from_intervals(reference.sequence.len(), &ambiguity);
    let called = finalize_consensus(build_consensus_draft(
        &reference.sequence,
        &pileup,
        &mask,
        ConsensusSettings {
            min_depth: args.min_depth,
            min_snp_frac,
            min_indel_frac,
            repeat_min_dominance: 0.85,
            repeat_max_entropy: 0.8,
        },
    ));
    let local = if args.local_assembly {
        Some(mapping::assemble(
            sampled_args.as_ref().unwrap_or(args),
            &called,
            &ambiguity,
            signals,
            &reference.selected_name,
        )?)
    } else {
        None
    };
    let sequence = local
        .as_ref()
        .map_or(called.sequence.as_slice(), |r| r.sequence.as_slice());
    let remaining_ambiguity = local
        .as_ref()
        .map_or(ambiguity.len(), |r| r.unresolved_ambiguity);
    let complete = remaining_ambiguity == 0
        && !sequence.contains(&b'N')
        && local.as_ref().is_none_or(|r| r.resolved);
    write_consensus_fasta(
        &args.out.join("candidate.fasta"),
        &reference.selected_name,
        sequence,
    )?;
    let strict_length = sequence.len();
    let sequence = local
        .as_ref()
        .and_then(|r| r.maximum_sequence.as_deref())
        .unwrap_or(sequence);
    let consensus_written = complete || args.best_effort;
    // Tandem-repeat lengths in the exported consensus are measured only by reads
    // carrying both flanks; label the rest in the consensus's own ambiguity TSV.
    let tandem_repeats = match (args.best_effort, local.is_some()) {
        (true, true) => {
            let local_args = sampled_args.as_ref().unwrap_or(args);
            let reads = Reads {
                r1: local_args.r1.clone().expect("local assembly reads"),
                r2: local_args.r2.clone().expect("local assembly reads"),
            };
            let (flagged, report) = mapping::tandem::assess(sequence, true, &reads)?;
            append_tandem_rows(
                &args.out.join("maximum_consensus.ambiguity.tsv"),
                &reference.selected_name,
                sequence.len(),
                &flagged,
            )?;
            report
        }
        _ => json!({"status":"not_assessed","reason":"requires best-effort local assembly reads"}),
    };
    report::write_ambiguity_tsv(&args.out.join("ambiguity.tsv"), &ambiguity)?;
    report::write_uncertainty_tsv(&args.out.join("uncertainty.tsv"), &called.uncertainty)?;
    report::write_base_votes_tsv(&args.out.join("base_votes.tsv"), &called.base_votes)?;
    if args.emit_variants {
        report::write_variants_tsv(&args.out.join("variants.tsv"), &called.variants)?;
    }
    if args.emit_coverage {
        report::write_coverage_tsv(&args.out.join("coverage.tsv"), &called.coverage)?;
    }
    crate::write_json(
        &args.out.join("report.json"),
        &json!({
            "schema":if local.is_some(){"mitovar-rust.polish.v2"}else{"mitovar-rust.polish.v1"},"backend":if local.is_some(){"rust-fragment-pileup-local-graph"}else{"rust-fragment-pileup"},"workflow":if local.is_some(){"reference_guided_local_assembly"}else{"reference_guided_bam_consensus"},
            "status":if complete{"complete"}else if args.best_effort{"best_effort"}else{"unresolved"},
            "validation_mode":if !args.local_assembly{"pileup_only"}else if args.deep_validation{"deep_validation"}else{"fast_consensus"},
            "final_published":complete,"consensus_written":consensus_written,"strict_validation_passed":complete,
            "reference":args.reference.canonicalize()?,"bam":args.bam.canonicalize()?,
            "contig":reference.selected_name,"reference_length":reference.sequence.len(),"consensus_length":sequence.len(),
            "sequence_output":{"status":if complete {"complete"} else if args.best_effort {"best_effort"} else if sequence.iter().any(|b| b"ACGT".contains(b)) {"partial_available"} else {"no_called_bases"},
                "fasta":if args.best_effort{"consensus.fasta"}else{"candidate.fasta"},
                "ambiguity_file":if args.best_effort{Some("maximum_consensus.ambiguity.tsv")}else if local.is_some(){Some("structural_ambiguity.tsv")}else{None},
                "gene_quality":"not_assessed", "ncr1_required_for_partial_output":false},
            "strict_candidate":{"fasta":"candidate.fasta","length":strict_length},
            "ambiguous_bases":sequence.iter().filter(|b| **b == b'N').count(),"ambiguity_intervals":remaining_ambiguity,
            "baseline_ambiguity_intervals":ambiguity.len(),"local_assembly":local.as_ref().map(|r|&r.report),
            "parameters":{"min_mapq":args.min_mapq,"min_baseq":args.min_baseq,"min_depth":args.min_depth,"min_alt_frac":min_indel_frac,"min_snp_frac":min_snp_frac,"min_indel_frac":min_indel_frac,"window":args.window,"max_records":args.max_records,"repeat_aware":true},
            "base_voting":{"method":"quality-weighted named fragments; overlapping mates counted once","weight":"min(baseQ,40) * min(MAPQ,60)","ordinary_snp_requires_strict_majority":true,"repeat_min_dominance":0.85,"repeat_max_entropy":0.8,"minimum_reported_minor_weighted_fraction":0.05,"mixed_sites":called.base_votes.len(),"called_mixed_sites":called.base_votes.iter().filter(|v| v.consensus != 'N').count(),"votes_file":"base_votes.tsv","mixed_sites_are_not_proof_of_heteroplasmy":true,"mixed_sample":mixed_sample(&called.base_votes)},
            "subsampling":sampling_report,"tandem_repeats":tandem_repeats,
            "records":stats.total_records,"filtered_records":stats.filtered_records,"overlapping_fragment_observations_collapsed":stats.overlapping_fragment_observations_collapsed,
            "coordinate_system":if args.best_effort{"pileup TSVs remain reference-based; structural_ambiguity.tsv matches strict candidate.fasta; maximum_consensus.ambiguity.tsv matches consensus.fasta and maximum_consensus.fasta"}else if local.is_some(){"pileup TSVs remain reference-based; structural_ambiguity.tsv is in candidate.fasta coordinates; see local_assembly/report.json for replacement coordinates"}else{"reference coordinates before accepted indels; TSV positions always reference-based"},
            "circularity_verified":false,"limitations":if local.is_some(){vec!["pileup reports describe the baseline; structural replacements are recorded separately","local reference-guided assembly does not certify whole-genome circularity","BAM and matching paired FASTQ required"]}else{vec!["uses existing BAM mappings; does not align FASTQ","no reference-filled N calls or structural repair","BAM accepted; CRAM and remote retrieval are outside this standalone version"]}
        }),
    )?;
    if args.best_effort {
        fs::copy(
            args.out.join("maximum_consensus.fasta"),
            args.out.join("consensus.fasta"),
        )?;
    } else if complete {
        fs::copy(
            args.out.join("candidate.fasta"),
            args.out.join("consensus.fasta"),
        )?;
    }
    Ok(consensus_written)
}

/// Append one row per flagged block, split at the origin of the circular consensus.
fn append_tandem_rows(
    path: &Path,
    contig: &str,
    length: usize,
    flagged: &[(usize, mapping::tandem::Block, &str)],
) -> Result<()> {
    let mut file = fs::OpenOptions::new().append(true).open(path)?;
    for (index, block, reason) in flagged {
        let segments = if block.end <= length {
            vec![(block.start + 1, block.end)]
        } else {
            vec![(block.start + 1, length), (1, block.end - length)]
        };
        for (start, end) in segments {
            writeln!(
                file,
                "{contig}\t{start}\t{end}\ttandem_repeat_{}\t{reason}\t0",
                index + 1
            )?;
        }
    }
    Ok(())
}

const MIXED_MIN_MINOR: f64 = 0.2;
const MIXED_MIN_SITES: usize = 10;

/// Minor-allele fractions of non-repeat mixed sites. Many sites with a substantial
/// second allele signal a mixed sample; this is not proof of heteroplasmy.
fn mixed_sample(votes: &[crate::consensus::BaseVoteRecord]) -> serde_json::Value {
    let mut bins = [0usize; 4];
    let mut repeat = 0;
    for vote in votes {
        if vote.is_repeat_ambiguous {
            repeat += 1;
            continue;
        }
        let mut called: Vec<f64> = vote
            .alleles
            .iter()
            .filter(|a| "ACGT".contains(a.base))
            .map(|a| a.weighted_fraction)
            .collect();
        called.sort_by(|a, b| b.total_cmp(a));
        let minor = called.get(1).copied().unwrap_or(0.0);
        let bin = match minor {
            m if m >= 0.3 => 3,
            m if m >= MIXED_MIN_MINOR => 2,
            m if m >= 0.1 => 1,
            m if m >= 0.05 => 0,
            _ => continue,
        };
        bins[bin] += 1;
    }
    let substantial = bins[2] + bins[3];
    json!({"method":"second-largest weighted A/C/G/T fraction at mixed sites, repeat-ambiguous sites excluded",
        "sites_by_minor_fraction":{"0.05-0.10":bins[0],"0.10-0.20":bins[1],"0.20-0.30":bins[2],"0.30+":bins[3]},
        "repeat_ambiguous_sites_excluded":repeat,"sites_with_minor_at_least_0.2":substantial,"signal_threshold":MIXED_MIN_SITES,
        "status":if substantial >= MIXED_MIN_SITES {"mixed_sample_signal"} else {"not_detected"},
        "limitation":"signal, not proof of heteroplasmy or contamination; mismapped divergent reads also create mixed sites"})
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::{AlleleVote, BaseVoteRecord};

    fn vote(minor: f64, repeat: bool) -> BaseVoteRecord {
        let allele = |base, fraction: f64| AlleleVote {
            base,
            depth: (fraction * 100.0) as u32,
            fraction,
            weighted_fraction: fraction,
        };
        BaseVoteRecord {
            pos: 1,
            reference: 'A',
            consensus: 'A',
            depth: 100,
            alleles: vec![
                allele('A', 1.0 - minor),
                allele('G', minor),
                allele('N', 0.0),
            ],
            is_repeat_ambiguous: repeat,
        }
    }

    #[test]
    fn many_substantial_minor_alleles_signal_a_mixed_sample() {
        let mut votes: Vec<_> = (0..9).map(|_| vote(0.35, false)).collect();
        votes.extend([vote(0.07, false), vote(0.15, false), vote(0.25, true)]);
        let report = mixed_sample(&votes);
        assert_eq!(report["sites_by_minor_fraction"]["0.30+"], 9);
        assert_eq!(report["sites_by_minor_fraction"]["0.05-0.10"], 1);
        assert_eq!(report["sites_by_minor_fraction"]["0.10-0.20"], 1);
        assert_eq!(report["repeat_ambiguous_sites_excluded"], 1);
        assert_eq!(report["status"], "not_detected");
        votes.push(vote(0.22, false));
        let report = mixed_sample(&votes);
        assert_eq!(report["sites_with_minor_at_least_0.2"], 10);
        assert_eq!(report["status"], "mixed_sample_signal");
    }
}
