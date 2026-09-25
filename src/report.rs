use crate::{
    consensus::{BaseVoteRecord, CoverageRecord, UncertaintyRecord, VariantRecord},
    repeat::AmbiguityInterval,
};
use anyhow::{Context, Result};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

pub fn write_base_votes_tsv(path: &Path, records: &[BaseVoteRecord]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "pos\tref\tconsensus\tdepth\tallele\tallele_depth\tallele_frac\tweighted_frac\tis_repeat_ambiguous")?;
    for record in records {
        for allele in &record.alleles {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}",
                record.pos,
                record.reference,
                record.consensus,
                record.depth,
                allele.base,
                allele.depth,
                allele.fraction,
                allele.weighted_fraction,
                record.is_repeat_ambiguous
            )?;
        }
    }
    writer.flush()?;
    Ok(())
}

pub fn write_variants_tsv(path: &Path, variants: &[VariantRecord]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("failed to create variants.tsv: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "pos\tref\talt\tdepth\talt_depth\talt_frac\tbaseq_mean\tmapq_mean\tflags\tregion_flag\tis_repeat_ambiguous"
    )?;

    for record in variants {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.3}\t{:.3}\t{}\t{}\t{}",
            record.pos,
            record.r#ref,
            record.alt,
            record.depth,
            record.alt_depth,
            record.alt_frac,
            record.baseq_mean,
            record.mapq_mean,
            record.flags,
            record.region_flag,
            record.is_repeat_ambiguous
        )?;
    }

    writer.flush()?;
    Ok(())
}

pub fn write_coverage_tsv(path: &Path, coverage: &[CoverageRecord]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("failed to create coverage.tsv: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "pos\tdepth\tmapq_median\tentropy\tis_repeat_ambiguous"
    )?;

    for record in coverage {
        writeln!(
            writer,
            "{}\t{}\t{:.3}\t{:.6}\t{}",
            record.pos,
            record.depth,
            record.mapq_median,
            record.entropy,
            record.is_repeat_ambiguous
        )?;
    }

    writer.flush()?;
    Ok(())
}

pub fn write_uncertainty_tsv(path: &Path, records: &[UncertaintyRecord]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("failed to create uncertainty.tsv: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "pos\tref\tmajor_base\tdepth\ttotal_obs\tprimary_obs\tsec_supp_obs\tmajor_frac\tmapq_median\tentropy\tis_repeat_ambiguous\treasons"
    )?;

    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.3}\t{:.6}\t{}\t{}",
            record.pos,
            record.r#ref,
            record.major_base,
            record.depth,
            record.total_obs,
            record.primary_obs,
            record.sec_supp_obs,
            record.major_frac,
            record.mapq_median,
            record.entropy,
            record.is_repeat_ambiguous,
            record.reasons,
        )?;
    }

    writer.flush()?;
    Ok(())
}

pub fn write_ambiguity_tsv(path: &Path, intervals: &[AmbiguityInterval]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("failed to create ambiguity.tsv: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "start\tend\tlen\ttriggers\tmedian_mapq\tsec_supp_frac\tdepth_ratio\tmean_entropy\taction"
    )?;

    for interval in intervals {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{:.3}\t{:.6}\t{:.6}\t{:.6}\t{}",
            interval.start + 1,
            interval.end + 1,
            interval.len,
            interval.triggers,
            interval.median_mapq,
            interval.sec_supp_frac,
            interval.depth_ratio,
            interval.mean_entropy,
            interval.action
        )?;
    }

    writer.flush()?;
    Ok(())
}
