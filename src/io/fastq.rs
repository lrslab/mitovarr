use anyhow::{bail, ensure, Context, Result};
use flate2::read::MultiGzDecoder;
use serde::Serialize;
use std::{
    fs::File,
    io::{BufRead, BufReader, Read, Write},
    path::Path,
};

const MAX_LINE: u64 = 1_048_576;

pub fn open(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut reader =
        BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
    if reader.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

#[derive(Clone, Debug)]
pub struct Record {
    pub header: String,
    pub sequence: Vec<u8>,
    pub plus: String,
    pub quality: Vec<u8>,
}

impl Record {
    pub fn write(&self, out: &mut impl Write) -> Result<()> {
        writeln!(out, "{}", self.header)?;
        out.write_all(&self.sequence)?;
        writeln!(out, "\n{}", self.plus)?;
        out.write_all(&self.quality)?;
        writeln!(out)?;
        Ok(())
    }

    pub fn name(&self, mate: u8) -> Result<&str> {
        let mut fields = self.header[1..].split_ascii_whitespace();
        let token = fields.next().context("empty FASTQ name")?;
        let id = if let Some(id) = token.strip_suffix("/1") {
            ensure!(mate == 1, "R2 contains /1 read: {token}");
            id
        } else if let Some(id) = token.strip_suffix("/2") {
            ensure!(mate == 2, "R1 contains /2 read: {token}");
            id
        } else {
            token
        };
        ensure!(!id.is_empty(), "empty FASTQ pair name");
        if let Some(label) = fields
            .next()
            .filter(|label| label.starts_with("1:") || label.starts_with("2:"))
        {
            ensure!(
                label.as_bytes()[0] == b'0' + mate,
                "wrong CASAVA mate: {}",
                self.header
            );
        }
        Ok(id)
    }
}

fn line(reader: &mut dyn BufRead) -> Result<Option<String>> {
    let mut text = String::new();
    let size = reader.take(MAX_LINE + 1).read_line(&mut text)?;
    ensure!(size as u64 <= MAX_LINE, "FASTQ line exceeds 1 MiB");
    if size == 0 {
        return Ok(None);
    }
    if text.ends_with('\n') {
        text.pop();
    }
    if text.ends_with('\r') {
        text.pop();
    }
    Ok(Some(text))
}

fn next(reader: &mut dyn BufRead) -> Result<Option<Record>> {
    let Some(header) = line(reader)? else {
        return Ok(None);
    };
    ensure!(
        header.starts_with('@') && !header.contains('\0'),
        "invalid FASTQ header"
    );
    let sequence = line(reader)?
        .context("truncated FASTQ sequence")?
        .into_bytes();
    let plus = line(reader)?.context("truncated FASTQ separator")?;
    let quality = line(reader)?
        .context("truncated FASTQ quality")?
        .into_bytes();
    ensure!(
        plus.starts_with('+'),
        "invalid FASTQ separator for {header}"
    );
    ensure!(
        plus.len() == 1 || plus[1..] == header[1..],
        "FASTQ + header disagrees for {header}"
    );
    ensure!(
        !sequence.is_empty() && sequence.len() == quality.len(),
        "FASTQ sequence/quality length mismatch for {header}"
    );
    ensure!(
        sequence
            .iter()
            .all(|b| b"ACGTNRYKMSWBDHVacgtnrykmswbdhv".contains(b)),
        "invalid DNA in {header}"
    );
    ensure!(
        quality.iter().all(|q| (33..=126).contains(q)),
        "expected Phred+33 quality for {header}"
    );
    Ok(Some(Record {
        header,
        sequence,
        plus,
        quality,
    }))
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Scan {
    pub pairs: u64,
    pub bases: u64,
}

/// Stream complete synchronized pairs.
pub fn scan_pairs(
    r1: &Path,
    r2: &Path,
    mut visit: impl FnMut(u64, &Record, &Record) -> Result<()>,
) -> Result<Scan> {
    ensure!(
        r1.canonicalize()? != r2.canonicalize()?,
        "R1 and R2 must be different files"
    );
    let (mut one, mut two) = (open(r1)?, open(r2)?);
    let (mut pairs, mut bases) = (0, 0);
    loop {
        match (next(&mut one)?, next(&mut two)?) {
            (None, None) => break,
            (Some(a), Some(b)) => {
                ensure!(
                    a.name(1)? == b.name(2)?,
                    "FASTQ mates disagree at pair {}: {} vs {}",
                    pairs + 1,
                    a.header,
                    b.header
                );
                visit(pairs, &a, &b)?;
                pairs += 1;
                bases += (a.sequence.len() + b.sequence.len()) as u64;
            }
            _ => bail!("FASTQ files have different pair counts"),
        }
    }
    ensure!(pairs > 0, "paired FASTQ is empty");
    Ok(Scan { pairs, bases })
}
