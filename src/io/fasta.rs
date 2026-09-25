use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::{anyhow, Context, Result};
use noodles_fasta as fasta;

use crate::pileup::normalize_base;

#[derive(Debug, Clone)]
pub struct ReferenceGenome {
    pub selected_name: String,
    pub sequence: Vec<u8>,
}

pub fn load_reference(path: &Path, contig: Option<&str>) -> Result<ReferenceGenome> {
    let mut reader = fasta::io::Reader::new(super::fastq::open(path)?);

    let mut records = Vec::new();

    for result in reader.records() {
        let record = result.with_context(|| {
            format!(
                "failed to parse FASTA record from reference: {}",
                path.display()
            )
        })?;
        records.push(record);
    }

    if records.is_empty() {
        return Err(anyhow!("reference FASTA contains no records"));
    }

    let selected = if let Some(name) = contig {
        records
            .iter()
            .find(|record| record.name() == name.as_bytes())
            .cloned()
            .ok_or_else(|| anyhow!("contig '{}' not found in reference FASTA", name))?
    } else {
        if records.len() > 1 {
            return Err(anyhow!(
                "reference FASTA contains multiple records; pass --contig to select one"
            ));
        }
        records
            .first()
            .cloned()
            .ok_or_else(|| anyhow!("reference FASTA contains no records"))?
    };

    let selected_name = String::from_utf8_lossy(selected.name()).to_string();
    let sequence = selected
        .sequence()
        .as_ref()
        .iter()
        .map(|b| normalize_base(*b))
        .collect::<Vec<_>>();

    Ok(ReferenceGenome {
        selected_name,
        sequence,
    })
}

pub fn write_consensus_fasta(path: &Path, name: &str, sequence: &[u8]) -> Result<()> {
    write_fasta_records(path, &[(name.to_owned(), sequence.to_vec())])
}

pub fn write_fasta_records(path: &Path, records: &[(String, Vec<u8>)]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("failed to create FASTA output: {}", path.display()))?;
    let writer = BufWriter::new(file);
    let mut writer = fasta::io::Writer::new(writer);

    for (name, sequence) in records {
        let definition = fasta::record::Definition::new(name.as_str(), None);
        let sequence = fasta::record::Sequence::from(sequence.clone());
        writer
            .write_record(&fasta::Record::new(definition, sequence))
            .with_context(|| format!("failed to write FASTA output: {}", path.display()))?;
    }

    writer
        .into_inner()
        .flush()
        .with_context(|| format!("failed to flush FASTA output: {}", path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn unique_fasta_test_dir() -> std::path::PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!("mtasm-fasta-test-{}-{nanos}", std::process::id()))
    }

    #[test]
    fn test_normalize_base() {
        assert_eq!(normalize_base(b'a'), b'A');
        assert_eq!(normalize_base(b'u'), b'T');
        assert_eq!(normalize_base(b'X'), b'N');
    }

    #[test]
    fn test_load_reference_accepts_single_record_without_contig() {
        let dir = unique_fasta_test_dir();
        std::fs::create_dir_all(&dir).expect("create temp dir");
        let path = dir.join("single.fasta");
        std::fs::write(&path, ">MtDNA\nACGT\n").expect("write fasta");

        let reference = load_reference(&path, None).expect("load single-record fasta");
        assert_eq!(reference.selected_name, "MtDNA");
        assert_eq!(reference.sequence, b"ACGT");

        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn test_load_reference_requires_contig_for_multi_record_fasta() {
        let dir = unique_fasta_test_dir();
        std::fs::create_dir_all(&dir).expect("create temp dir");
        let path = dir.join("multi.fasta");
        std::fs::write(&path, ">chr1\nACGT\n>chr2\nTGCA\n").expect("write fasta");

        let err = load_reference(&path, None).expect_err("expected multi-record fasta to fail");
        assert!(
            err.to_string().contains("pass --contig to select one"),
            "unexpected error: {err:#}"
        );

        let reference = load_reference(&path, Some("chr2")).expect("load selected contig");
        assert_eq!(reference.selected_name, "chr2");
        assert_eq!(reference.sequence, b"TGCA");

        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_dir_all(dir);
    }
}
