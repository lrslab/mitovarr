//! Indexed remote BAM -> available mitochondrial pairs, without external tools.
use std::{
    collections::BTreeMap,
    fs::File,
    io::{self, BufWriter, Read, Seek, SeekFrom, Write},
    path::Path,
    time::{Duration, Instant},
};

use anyhow::{ensure, Context, Result};
use flate2::{write::GzEncoder, Compression};
use noodles_bam as bam;
use noodles_core::{Position, Region};
use noodles_sam::{
    self as sam,
    alignment::{io::Write as _, RecordBuf},
};
use serde_json::json;
use tracing::{debug, info};

use crate::cli::FetchCaendrArgs;

const BUCKET: &str = "https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/bam";
const RANGE_BYTES: u64 = 1_048_576;
const MATE_WINDOW: usize = 16_384;
const FASTQS: [&str; 4] = [
    "reads.R1.fastq.gz",
    "reads.R2.fastq.gz",
    "reads.singletons.fastq.gz",
    "reads.unpaired.fastq.gz",
];
type Pair = [Option<RecordBuf>; 2];

/// The only adapter needed by noodles: seekable HTTP with one bounded buffer.
struct RemoteBam {
    agent: ureq::Agent,
    url: String,
    position: u64,
    length: Option<u64>,
    buffer_start: u64,
    buffer: Vec<u8>,
    range_bytes: u64,
    requests: u64,
    bytes: u64,
}

impl RemoteBam {
    fn new(agent: ureq::Agent, url: &str) -> Self {
        Self {
            agent,
            url: url.into(),
            position: 0,
            length: None,
            buffer_start: 0,
            buffer: Vec::new(),
            range_bytes: RANGE_BYTES,
            requests: 0,
            bytes: 0,
        }
    }

    fn refill(&mut self) -> Result<()> {
        let started = Instant::now();
        let end = self.position.saturating_add(self.range_bytes - 1);
        let mut response = self
            .agent
            .get(&self.url)
            .header("Range", format!("bytes={}-{}", self.position, end))
            .header("Accept-Encoding", "identity")
            .call()
            .with_context(|| format!("fetch BAM bytes at {}", self.position))?;
        // This check prevents a server ignoring Range from returning the full BAM.
        ensure!(
            response.status().as_u16() == 206,
            "server ignored HTTP Range; refusing a full BAM download"
        );
        let range = response
            .headers()
            .get("Content-Range")
            .context("missing Content-Range")?
            .to_str()?;
        let (span, total) = range
            .strip_prefix("bytes ")
            .and_then(|s| s.split_once('/'))
            .context("invalid Content-Range")?;
        let (start, last) = span.split_once('-').context("invalid byte range")?;
        let (start, last, total) = (
            start.parse::<u64>()?,
            last.parse::<u64>()?,
            total.parse::<u64>()?,
        );
        ensure!(
            start == self.position && last >= start && last <= end && last < total,
            "unexpected BAM byte range: {range}"
        );
        self.buffer.clear();
        response
            .body_mut()
            .as_reader()
            .take(self.range_bytes + 1)
            .read_to_end(&mut self.buffer)?;
        ensure!(
            self.buffer.len() as u64 == last - start + 1,
            "truncated BAM range response"
        );
        self.buffer_start = start;
        self.length = Some(total);
        self.requests += 1;
        self.bytes += self.buffer.len() as u64;
        debug!(
            start,
            last,
            request = self.requests,
            seconds = started.elapsed().as_secs_f64(),
            "BAM range fetched"
        );
        Ok(())
    }
}

impl Read for RemoteBam {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() || self.length.is_some_and(|n| self.position >= n) {
            return Ok(0);
        }
        if self.position < self.buffer_start
            || self.position >= self.buffer_start + self.buffer.len() as u64
        {
            self.refill().map_err(io::Error::other)?;
        }
        let offset = (self.position - self.buffer_start) as usize;
        let n = out.len().min(self.buffer.len() - offset);
        out[..n].copy_from_slice(&self.buffer[offset..offset + n]);
        self.position += n as u64;
        Ok(n)
    }
}

impl Seek for RemoteBam {
    fn seek(&mut self, from: SeekFrom) -> io::Result<u64> {
        let position = match from {
            SeekFrom::Start(n) => Some(n),
            SeekFrom::Current(n) => self.position.checked_add_signed(n),
            SeekFrom::End(n) => {
                if self.length.is_none() {
                    self.refill().map_err(io::Error::other)?;
                }
                self.length.and_then(|length| length.checked_add_signed(n))
            }
        }
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid BAM seek"))?;
        self.position = position;
        Ok(position)
    }
}

fn mate_slot(record: &RecordBuf) -> Option<usize> {
    let flags = record.flags();
    if !flags.is_segmented() {
        return None;
    }
    match (flags.is_first_segment(), flags.is_last_segment()) {
        (true, false) => Some(0),
        (false, true) => Some(1),
        _ => None,
    }
}

fn primary(record: &RecordBuf) -> bool {
    !record.flags().is_secondary() && !record.flags().is_supplementary()
}

fn name(record: &RecordBuf) -> Result<&[u8]> {
    let name = record.name().context("BAM read has no QNAME")?;
    Ok(name.as_ref())
}

fn select_contig(header: &sam::Header, explicit: Option<&str>) -> Result<(usize, String)> {
    let choices: Vec<_> = header
        .reference_sequences()
        .iter()
        .enumerate()
        .filter_map(|(id, (name, _))| {
            let name = String::from_utf8_lossy(name).into_owned();
            let matches = explicit.map_or_else(
                || {
                    matches!(
                        name.to_ascii_lowercase().as_str(),
                        "mtdna"
                            | "chrm"
                            | "mt"
                            | "mitochondrion_genome"
                            | "mitochondria"
                            | "mitochondrion"
                    )
                },
                |wanted| wanted == name,
            );
            matches.then_some((id, name))
        })
        .collect();
    ensure!(
        choices.len() == 1,
        "cannot select one mitochondrial contig; use --contig with a BAM contig name"
    );
    Ok(choices.into_iter().next().expect("one contig"))
}

/// BAM SEQ is reference-oriented; reverse both SEQ and QUAL for reverse reads.
fn write_read(writer: &mut impl Write, record: &RecordBuf) -> Result<bool> {
    let mut sequence = record.sequence().as_ref().to_vec();
    ensure!(
        !sequence.contains(&b'='),
        "BAM uses '=' in SEQ; original bases cannot be recovered without a reference"
    );
    let missing_quality = record.quality_scores().is_empty();
    let mut quality = if missing_quality {
        vec![0; sequence.len()]
    } else {
        record.quality_scores().as_ref().to_vec()
    };
    ensure!(
        quality.len() == sequence.len() && quality.iter().all(|q| *q <= 93),
        "invalid BAM quality scores"
    );
    if record.flags().is_reverse_complemented() {
        sequence.reverse();
        quality.reverse();
        for base in &mut sequence {
            *base = match *base {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                b'M' => b'K',
                b'K' => b'M',
                b'R' => b'Y',
                b'Y' => b'R',
                b'V' => b'B',
                b'B' => b'V',
                b'H' => b'D',
                b'D' => b'H',
                b => b,
            };
        }
    }
    for q in &mut quality {
        *q += 33;
    }
    writer.write_all(b"@")?;
    writer.write_all(name(record)?)?;
    writer.write_all(b"\n")?;
    writer.write_all(&sequence)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(&quality)?;
    writer.write_all(b"\n")?;
    Ok(missing_quality)
}

fn write_bam(
    out: &Path,
    header: &sam::Header,
    pairs: &BTreeMap<Vec<u8>, Pair>,
    unpaired: &[RecordBuf],
) -> Result<()> {
    let mut records: Vec<_> = pairs
        .values()
        .flat_map(|pair| pair.iter().flatten())
        .chain(unpaired)
        .collect();
    records.sort_by_key(|record| {
        (
            record.reference_sequence_id().unwrap_or(usize::MAX),
            record
                .alignment_start()
                .map(usize::from)
                .unwrap_or(usize::MAX),
        )
    });
    let mut header = header.clone();
    header
        .header_mut()
        .get_or_insert_with(Default::default)
        .other_fields_mut()
        .insert(
            sam::header::record::value::map::header::tag::SORT_ORDER,
            "coordinate".into(),
        );
    let path = out.join("mt.bam");
    let mut writer = bam::io::Writer::new(File::create(&path)?);
    writer.write_header(&header)?;
    for record in records {
        writer.write_alignment_record(&header, record)?;
    }
    writer.try_finish()?;
    drop(writer);
    bam::bai::fs::write(out.join("mt.bam.bai"), &bam::fs::index(path)?)?;
    Ok(())
}

pub fn fetch(args: &FetchCaendrArgs) -> Result<bool> {
    ensure!(
        !args.strain.is_empty()
            && args
                .strain
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"_-".contains(&b)),
        "strain must contain only letters, digits, '_' or '-'"
    );
    let bam_url = format!("{}/{}/{}.bam", BUCKET, args.species.name(), args.strain);
    fetch_from_urls(args, &bam_url, &format!("{bam_url}.bai"))
}

fn fetch_from_urls(args: &FetchCaendrArgs, bam_url: &str, bai_url: &str) -> Result<bool> {
    ensure!(
        args.timeout_seconds > 0 && args.max_templates > 0,
        "timeout and template limit must be positive"
    );
    let agent = ureq::Agent::config_builder()
        .timeout_global(Some(Duration::from_secs(args.timeout_seconds)))
        .build()
        .new_agent();
    let mut response = agent
        .get(bai_url)
        .call()
        .context("fetch CaenDR BAM index (check species and strain)")?;
    let mut index_data = Vec::new();
    response
        .body_mut()
        .as_reader()
        .take(64 * 1024 * 1024 + 1)
        .read_to_end(&mut index_data)?;
    ensure!(index_data.len() <= 64 * 1024 * 1024, "BAI exceeds 64 MiB");
    let index = bam::bai::io::Reader::new(index_data.as_slice()).read_index()?;
    let mut reader = bam::io::IndexedReader::new(RemoteBam::new(agent, bam_url), index.clone());
    let header = reader.read_header()?;
    let (target_id, contig) = select_contig(&header, args.contig.as_deref())?;
    info!(
        species = args.species.name(),
        strain = args.strain,
        contig,
        "fetching indexed mitochondrial reads"
    );
    let mut pairs = BTreeMap::<Vec<u8>, Pair>::new();
    let mut unpaired = Vec::new();
    let region = Region::new(contig.as_bytes(), ..);
    reader.get_mut().get_mut().range_bytes = 8 * RANGE_BYTES;
    for record in reader.query(&header, &region)?.records() {
        let record = RecordBuf::try_from_alignment_record(&header, &record?)?;
        if !primary(&record) {
            continue;
        }
        let name = name(&record)?.to_vec();
        if let Some(slot) = mate_slot(&record) {
            ensure!(
                pairs.contains_key(&name) || pairs.len() + unpaired.len() < args.max_templates,
                "extraction exceeds --max-templates"
            );
            let pair = pairs.entry(name).or_default();
            if pair[slot].as_ref().is_none_or(|r| r.sequence().is_empty()) {
                pair[slot] = Some(record);
            }
        } else {
            ensure!(
                pairs.len() + unpaired.len() < args.max_templates,
                "extraction exceeds --max-templates"
            );
            unpaired.push(record);
        }
    }
    ensure!(
        !pairs.is_empty() || !unpaired.is_empty(),
        "no primary mitochondrial reads found"
    );

    // One non-iterative mate pass. Adjacent mate positions share a BAI window.
    let mut windows = BTreeMap::<_, usize>::new();
    for pair in pairs.values() {
        for slot in 0..2 {
            if pair[1 - slot].is_none() {
                if let Some(record) = &pair[slot] {
                    if let (Some(id), Some(pos)) = (
                        record.mate_reference_sequence_id(),
                        record.mate_alignment_start(),
                    ) {
                        // All primary records on this contig were already read.
                        if id != target_id {
                            *windows
                                .entry((id, (usize::from(pos) - 1) / MATE_WINDOW))
                                .or_default() += 1;
                        }
                    }
                }
            }
        }
    }
    info!(
        templates = pairs.len(),
        mate_windows = windows.len(),
        "recovering available mates"
    );
    let mut recovered = 0u64;
    for ((id, window), &needed) in &windows {
        let (_, metadata) = header
            .reference_sequences()
            .get_index(*id)
            .context("mate reference ID absent from BAM header")?;
        let start = window * MATE_WINDOW + 1;
        let end = ((window + 1) * MATE_WINDOW).min(usize::from(metadata.length()));
        if start > end {
            continue;
        }
        let interval = (Position::try_from(start)?..=Position::try_from(end)?).into();
        let reference_index = &index.reference_sequences()[*id];
        // Exact mate starts usually live in the finest BAI bin. Search it
        // first; long alignments in parent bins are only needed as fallback.
        let chunk_passes = [
            reference_index
                .bins()
                .get(&(4681 + window))
                .into_iter()
                .flat_map(|bin| bin.chunks().iter().copied())
                .collect::<Vec<_>>(),
            reader.index().query(*id, interval)?,
        ];
        let offsets = reference_index.index();
        let min_offset = offsets.get(*window).copied().unwrap_or_default();
        // Fetch through the next window plus one maximum BGZF block. Sparse
        // mate lookups should not each download an entire MiB unnecessarily.
        reader.get_mut().get_mut().range_bytes =
            offsets.get(window + 1).map_or(RANGE_BYTES, |next| {
                next.compressed()
                    .saturating_sub(min_offset.compressed())
                    .saturating_add(65_536)
                    .clamp(65_536, RANGE_BYTES)
            });
        let mut record = bam::Record::default();
        let mut remaining = needed;
        // BAI parent bins can include distant chunks beyond this window. On a
        // coordinate-sorted BAM, stop at the first record past its end instead
        // of letting a filtered query fetch those chunks again for every mate.
        'window: for chunks in chunk_passes {
            'pass: for chunk in chunks {
                reader.get_mut().seek(chunk.start().max(min_offset))?;
                while reader.get_ref().virtual_position() < chunk.end() {
                    if reader.read_record(&mut record)? == 0 {
                        break 'pass;
                    }
                    let reference_id = record.reference_sequence_id().transpose()?;
                    let position = record.alignment_start().transpose()?;
                    if reference_id.is_none_or(|r| r > *id)
                        || (reference_id == Some(*id)
                            && position.is_some_and(|p| usize::from(p) > end))
                    {
                        break 'pass;
                    }
                    if reference_id != Some(*id)
                        || position.is_none_or(|p| usize::from(p) < start)
                        || record.flags().is_secondary()
                        || record.flags().is_supplementary()
                        || !record
                            .name()
                            .is_some_and(|n| pairs.contains_key::<[u8]>(n.as_ref()))
                    {
                        continue;
                    }
                    let record = RecordBuf::try_from_alignment_record(&header, &record)?;
                    let Some(slot) = mate_slot(&record) else {
                        continue;
                    };
                    if let Some(pair) = pairs.get_mut(name(&record)?) {
                        if pair[slot].is_none()
                            && pair[1 - slot].as_ref().is_some_and(|other| {
                                other.mate_reference_sequence_id() == record.reference_sequence_id()
                                    && other.mate_alignment_start() == record.alignment_start()
                            })
                        {
                            pair[slot] = Some(record);
                            recovered += 1;
                            remaining -= 1;
                            if remaining == 0 {
                                break 'window;
                            }
                        }
                    }
                }
            }
        }
    }
    let remote = reader.get_ref().get_ref();
    let (requests, downloaded) = (remote.requests, remote.bytes);
    let mut outputs = Vec::new();
    for name in FASTQS {
        outputs.push(GzEncoder::new(
            BufWriter::new(File::create(args.out.join(name))?),
            Compression::fast(),
        ));
    }
    let (mut paired, mut singletons, mut unpaired_count, mut missing_sequence, mut missing_quality) =
        (0u64, 0u64, 0u64, 0u64, 0u64);
    for pair in pairs.values() {
        missing_sequence += pair
            .iter()
            .flatten()
            .filter(|r| r.sequence().is_empty())
            .count() as u64;
        let first = pair[0].as_ref().filter(|r| !r.sequence().is_empty());
        let second = pair[1].as_ref().filter(|r| !r.sequence().is_empty());
        match (first, second) {
            (Some(a), Some(b)) => {
                missing_quality += u64::from(write_read(&mut outputs[0], a)?)
                    + u64::from(write_read(&mut outputs[1], b)?);
                paired += 1;
            }
            (Some(read), None) | (None, Some(read)) => {
                missing_quality += u64::from(write_read(&mut outputs[2], read)?);
                singletons += 1;
            }
            (None, None) => (),
        }
    }
    for record in &unpaired {
        if record.sequence().is_empty() {
            missing_sequence += 1;
            continue;
        }
        missing_quality += u64::from(write_read(&mut outputs[3], record)?);
        unpaired_count += 1;
    }
    for writer in outputs {
        writer.finish()?.flush()?;
    }
    write_bam(&args.out, &header, &pairs, &unpaired)?;
    let out = args.out.canonicalize()?;
    let report = json!({
        "schema":"mitovar-rust.caendr.v1", "status":if paired > 0 {"complete"} else {"no_pairs"},
        "backend":"rust-noodles-http-range", "read_source":"mt-bam-extracted",
        "species":args.species.name(), "strain":args.strain,
        "source_bam":{"url":bam_url,"index_url":bai_url,"contig":contig},
        "counts":{"paired_templates":paired,"singletons":singletons,"unpaired":unpaired_count,"mates_recovered":recovered,"missing_sequence_records":missing_sequence,"reads_without_quality":missing_quality},
        "outputs":{"r1":{"path":out.join(FASTQS[0])},"r2":{"path":out.join(FASTQS[1])},"singletons":{"path":out.join(FASTQS[2])},"unpaired":{"path":out.join(FASTQS[3])}},
        "extracted_bam":out.join("mt.bam"),"extracted_bam_index":out.join("mt.bam.bai"),
        "bam_range_requests":requests,"bam_bytes_received":downloaded,"mate_windows":windows.len(),
        "parameters":{"timeout_seconds":args.timeout_seconds,"max_templates":args.max_templates},
        "filters":"primary alignments; retain duplicate/QC-fail flags; no MAPQ cutoff",
        "quality_source":"stored BAM QUAL; absent quality becomes Phred 0",
        "limitations":["Only indexed mitochondrial reads and their available mates; not original whole-sample FASTQ","Unplaced/absent mates and hard-clipped sequence cannot be recovered; no second baiting pass"]
    });
    crate::write_json(&out.join("recovery.json"), &report)?;
    crate::write_json(&out.join("report.json"), &report)?;
    info!(
        paired,
        singletons,
        unpaired = unpaired_count,
        recovered,
        "CaenDR extraction complete"
    );
    Ok(paired > 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use sam::{
        alignment::{
            record::{
                cigar::{op::Kind, Op},
                Flags, MappingQuality,
            },
            record_buf::{Cigar, QualityScores, Sequence},
        },
        header::record::value::{map::ReferenceSequence, Map},
    };
    use std::{
        fs,
        io::{BufRead, BufReader},
        net::{TcpListener, TcpStream},
        num::NonZeroUsize,
        path::PathBuf,
        sync::{
            atomic::{AtomicBool, AtomicU64, Ordering},
            Arc,
        },
        thread,
    };

    static NEXT: AtomicU64 = AtomicU64::new(0);

    struct Fixture(PathBuf);
    impl Fixture {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "mtasm-caendr-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir_all(&path).unwrap();
            Self(path)
        }
        fn args(&self) -> FetchCaendrArgs {
            let out = self.0.join("reads");
            fs::create_dir(&out).unwrap();
            FetchCaendrArgs {
                species: crate::cli::Species::Cbr,
                strain: "AF16".into(),
                out,
                contig: None,
                timeout_seconds: 5,
                max_templates: 100,
            }
        }
    }
    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    struct Server {
        url: String,
        stop: Arc<AtomicBool>,
        thread: Option<thread::JoinHandle<()>>,
    }
    impl Server {
        fn new(bam: Vec<u8>, bai: Vec<u8>, ignore_range: bool) -> Self {
            let listener = TcpListener::bind("127.0.0.1:0").unwrap();
            let url = format!("http://{}", listener.local_addr().unwrap());
            let stop = Arc::new(AtomicBool::new(false));
            let done = stop.clone();
            let thread = thread::spawn(move || {
                for stream in listener.incoming() {
                    let mut stream = stream.unwrap();
                    if done.load(Ordering::Relaxed) {
                        break;
                    }
                    stream
                        .set_read_timeout(Some(Duration::from_secs(5)))
                        .unwrap();
                    let mut input = BufReader::new(stream.try_clone().unwrap());
                    let mut first = String::new();
                    input.read_line(&mut first).unwrap();
                    let is_index = first.split_whitespace().nth(1).unwrap().ends_with(".bai");
                    let mut range = None;
                    loop {
                        let mut line = String::new();
                        input.read_line(&mut line).unwrap();
                        if line.trim().is_empty() {
                            break;
                        }
                        if let Some(value) = line.to_ascii_lowercase().strip_prefix("range: bytes=")
                        {
                            let (a, b) = value.trim().split_once('-').unwrap();
                            range =
                                Some((a.parse::<usize>().unwrap(), b.parse::<usize>().unwrap()));
                        }
                    }
                    let data = if is_index { &bai } else { &bam };
                    if !is_index {
                        assert!(range.is_some(), "BAM must never use an unbounded GET");
                    }
                    let (status, body, content_range) =
                        if let Some((start, end)) = range.filter(|_| !ignore_range) {
                            let end = end.min(data.len() - 1);
                            (
                                "206 Partial Content",
                                &data[start..=end],
                                format!("Content-Range: bytes {start}-{end}/{}\r\n", data.len()),
                            )
                        } else {
                            ("200 OK", data.as_slice(), String::new())
                        };
                    write!(stream,"HTTP/1.1 {status}\r\nContent-Length: {}\r\n{content_range}Connection: close\r\n\r\n",body.len()).unwrap();
                    let _ = stream.write_all(body);
                }
            });
            Self {
                url,
                stop,
                thread: Some(thread),
            }
        }
    }
    impl Drop for Server {
        fn drop(&mut self) {
            self.stop.store(true, Ordering::Relaxed);
            let _ = TcpStream::connect(self.url.trim_start_matches("http://"));
            self.thread.take().unwrap().join().unwrap();
        }
    }

    fn record(
        name: &str,
        flags: Flags,
        reference: usize,
        start: usize,
        mate: Option<(usize, usize)>,
    ) -> RecordBuf {
        let mut record = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_reference_sequence_id(reference)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_mapping_quality(MappingQuality::new(0).unwrap())
            .set_cigar([Op::new(Kind::Match, 8)].into_iter().collect::<Cigar>())
            .set_sequence(Sequence::from(b"ACGTRMKB".to_vec()))
            .set_quality_scores(QualityScores::from((1..=8).collect::<Vec<_>>()))
            .build();
        if let Some((id, pos)) = mate {
            *record.mate_reference_sequence_id_mut() = Some(id);
            *record.mate_alignment_start_mut() = Position::new(pos);
        }
        record
    }

    fn bam_fixture(dir: &Path, distant_data: bool) -> (Vec<u8>, Vec<u8>) {
        let first = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        let second = Flags::SEGMENTED | Flags::LAST_SEGMENT | Flags::REVERSE_COMPLEMENTED;
        // Cover both a finest-bin mate and one crossing a 16 KiB bin boundary.
        let mate_pos = if distant_data { 1010 } else { 16380 };
        let mut records = vec![
            record("mtpair", first, 0, 10, Some((0, 60))),
            record("mtpair", second | Flags::SECONDARY, 0, 20, Some((0, 10))),
            record(
                "off",
                first | Flags::DUPLICATE | Flags::QC_FAIL,
                0,
                30,
                Some((1, mate_pos)),
            ),
            record(
                "orphan",
                first | Flags::MATE_UNMAPPED,
                0,
                40,
                Some((0, 100)),
            ),
            record("unpaired", Flags::empty(), 0, 50, None),
            record("mtpair", second, 0, 60, Some((0, 10))),
            record("noquality", Flags::empty(), 0, 70, None),
            record("unrelated", second, 1, 1010, Some((1, 3000))),
            record("off", second, 1, mate_pos, Some((0, 30))),
        ];
        *records[6].quality_scores_mut() = QualityScores::default();
        if distant_data {
            // A missing indexed mate must end its parent-bin fallback at the
            // window boundary, even after the other mate has been recovered.
            *records[3].mate_reference_sequence_id_mut() = Some(1);
            *records[3].mate_alignment_start_mut() = Position::new(1011);
            // A distant alignment in a BAI parent bin is returned as a possible
            // chunk for the mate window, even though its POS is past the end.
            let mut boundary = record("past_window", first, 1, 20_000, None);
            *boundary.cigar_mut() = [Op::new(Kind::Match, 8), Op::new(Kind::Deletion, 2_000_000)]
                .into_iter()
                .collect();
            records.push(boundary);
            let mut distant = record("distant", first, 1, 30_000, None);
            let length = 1_500_000;
            let mut state = 1u32;
            let quality: Vec<_> = (0..length)
                .map(|_| {
                    state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    ((state >> 16) % 94) as u8
                })
                .collect();
            *distant.sequence_mut() = Sequence::from(vec![b'A'; length]);
            *distant.quality_scores_mut() = QualityScores::from(quality);
            *distant.cigar_mut() = [Op::new(Kind::Match, length)].into_iter().collect();
            records.push(distant);
        }
        let mut header = sam::Header::builder()
            .add_reference_sequence(
                "MtDNA",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .add_reference_sequence(
                "I",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(4_000_000).unwrap()),
            )
            .build();
        header
            .header_mut()
            .get_or_insert_with(Default::default)
            .other_fields_mut()
            .insert(
                sam::header::record::value::map::header::tag::SORT_ORDER,
                "coordinate".into(),
            );
        let path = dir.join("input.bam");
        let mut writer = bam::io::Writer::new(File::create(&path).unwrap());
        writer.write_header(&header).unwrap();
        for record in records {
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish().unwrap();
        drop(writer);
        let index = bam::fs::index(&path).unwrap();
        bam::bai::fs::write(dir.join("input.bam.bai"), &index).unwrap();
        (
            fs::read(path).unwrap(),
            fs::read(dir.join("input.bam.bai")).unwrap(),
        )
    }

    #[test]
    fn http_extraction_recovers_mates_and_original_orientation() {
        let fixture = Fixture::new();
        let (bam, bai) = bam_fixture(&fixture.0, false);
        let server = Server::new(bam, bai, false);
        let args = fixture.args();
        assert!(fetch_from_urls(
            &args,
            &format!("{}/input.bam", server.url),
            &format!("{}/input.bam.bai", server.url)
        )
        .unwrap());
        let scan = crate::io::fastq::scan_pairs(
            &args.out.join(FASTQS[0]),
            &args.out.join(FASTQS[1]),
            |_, a, b| {
                assert_eq!(&a.sequence, b"ACGTRMKB");
                assert_eq!(&b.sequence, b"VMKYACGT");
                assert_eq!(b.quality, (1..=8).rev().map(|q| q + 33).collect::<Vec<_>>());
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(scan.pairs, 2);
        let report: serde_json::Value =
            serde_json::from_slice(&fs::read(args.out.join("recovery.json")).unwrap()).unwrap();
        assert_eq!(report["counts"]["mates_recovered"], 1);
        assert_eq!(report["mate_windows"], 1);
        assert_eq!(report["counts"]["singletons"], 1);
        assert_eq!(report["counts"]["unpaired"], 2);
        assert_eq!(report["counts"]["reads_without_quality"], 1);
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(args.out.join("mt.bam"))
            .unwrap();
        let header = reader.read_header().unwrap();
        let region = Region::new("I", ..);
        let off: Vec<_> = reader
            .query(&header, &region)
            .unwrap()
            .records()
            .collect::<std::result::Result<_, _>>()
            .unwrap();
        assert_eq!(off.len(), 1);
    }

    #[test]
    fn mate_window_does_not_fetch_distant_parent_bin_data() {
        let fixture = Fixture::new();
        let (bam, bai) = bam_fixture(&fixture.0, true);
        assert!(bam.len() > RANGE_BYTES as usize);
        let server = Server::new(bam, bai, false);
        let args = fixture.args();
        assert!(fetch_from_urls(
            &args,
            &format!("{}/input.bam", server.url),
            &format!("{}/input.bam.bai", server.url)
        )
        .unwrap());
        let report: serde_json::Value =
            serde_json::from_slice(&fs::read(args.out.join("recovery.json")).unwrap()).unwrap();
        assert_eq!(report["counts"]["mates_recovered"], 1);
        assert_eq!(report["bam_range_requests"], 1);
    }

    #[test]
    fn range_ignoring_server_cannot_trigger_full_download() {
        let fixture = Fixture::new();
        let (bam, bai) = bam_fixture(&fixture.0, false);
        let server = Server::new(bam, bai, true);
        let args = fixture.args();
        let error = fetch_from_urls(
            &args,
            &format!("{}/input.bam", server.url),
            &format!("{}/input.bam.bai", server.url),
        )
        .unwrap_err();
        assert!(format!("{error:#}").contains("refusing a full BAM download"));
        assert!(!args.out.join(FASTQS[0]).exists());
    }
}
