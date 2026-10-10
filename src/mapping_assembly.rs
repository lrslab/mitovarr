use crate::{
    cli::MappingArgs,
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
use std::{fs, path::Path};

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
            "base_voting":{"method":"quality-weighted named fragments; overlapping mates counted once","weight":"min(baseQ,40) * min(MAPQ,60)","ordinary_snp_requires_strict_majority":true,"repeat_min_dominance":0.85,"repeat_max_entropy":0.8,"minimum_reported_minor_weighted_fraction":0.05,"mixed_sites":called.base_votes.len(),"called_mixed_sites":called.base_votes.iter().filter(|v| v.consensus != 'N').count(),"votes_file":"base_votes.tsv","mixed_sites_are_not_proof_of_heteroplasmy":true},
            "subsampling":sampling_report,
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
