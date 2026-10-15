//! Tandem-repeat length evidence in the exported consensus. A block's length is
//! measured only when reads carry both of its unique flanking anchors; otherwise
//! the consensus length is a reference-guided span and is labelled unmeasured,
//! with the longest anchored read entry as a lower bound. Depth is never used to
//! infer copy number.
use crate::{cli::Reads, io::fastq, kmer::reverse_complement};
use anyhow::Result;
use serde_json::{json, Value};
use std::collections::BTreeMap;

const MIN_PERIOD: usize = 2;
const MAX_PERIOD: usize = 75;
const WINDOW: usize = 20;
/// At least 17 of 20 window positions match the base one period later.
const MIN_WINDOW_MATCHES: usize = 17;
/// Shorter blocks are spanned by ordinary short reads and called by pileup.
const MIN_BLOCK: usize = 60;
const ANCHOR: usize = 12;
/// Unique anchors may sit this far outside the block.
const ANCHOR_SEARCH: usize = 20;
/// Reads needed for a gap or a lower bound.
const MIN_READS: usize = 3;

const MASK: u32 = (1 << (2 * ANCHOR)) - 1;

fn two_bit(base: u8) -> Option<u32> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn code(kmer: &[u8]) -> u32 {
    kmer.iter().fold(0, |value, &base| {
        (value << 2) | two_bit(base).expect("called anchor base")
    })
}

pub(crate) const UNMEASURED: &str = "tandem_repeat_length_unmeasured";
pub(crate) const CONFLICT: &str = "tandem_repeat_length_conflict";

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Block {
    /// 0-based start and exclusive end on the consensus; `end` may exceed the
    /// sequence length for a block crossing the origin of a circular sequence.
    pub start: usize,
    pub end: usize,
    pub period: usize,
}

/// Maximal runs of period-p self-similarity, merged across periods. The period
/// reported is the one contributing the longest run.
pub(crate) fn blocks(sequence: &[u8], circular: bool) -> Vec<Block> {
    let n = sequence.len();
    let text: Vec<u8> = if circular {
        sequence.iter().chain(sequence).copied().collect()
    } else {
        sequence.to_vec()
    };
    let mut found = Vec::new();
    for period in MIN_PERIOD..=MAX_PERIOD {
        if text.len() < period + WINDOW {
            break;
        }
        let matches: Vec<bool> = (0..text.len() - period)
            .map(|i| text[i] == text[i + period] && b"ACGT".contains(&text[i]))
            .collect();
        let mut count = matches[..WINDOW].iter().filter(|&&m| m).count();
        let mut run: Option<(usize, usize)> = None;
        for start in 0..=matches.len() - WINDOW {
            if start > 0 {
                count += usize::from(matches[start + WINDOW - 1]);
                count -= usize::from(matches[start - 1]);
            }
            if count >= MIN_WINDOW_MATCHES {
                run = Some(run.map_or((start, start), |(a, _)| (a, start)));
            } else if let Some((a, b)) = run.take() {
                found.push(Block {
                    start: a,
                    end: b + WINDOW + period,
                    period,
                });
            }
        }
        if let Some((a, b)) = run {
            found.push(Block {
                start: a,
                end: (b + WINDOW + period).min(text.len()),
                period,
            });
        }
    }
    found.retain(|b| b.end - b.start >= MIN_BLOCK.max(2 * b.period) && (!circular || b.start < n));
    found.sort_by_key(|b| (b.start, std::cmp::Reverse(b.end)));
    let mut merged: Vec<(Block, usize)> = Vec::new();
    for b in found {
        let length = b.end - b.start;
        match merged.last_mut() {
            Some((m, best)) if b.start < m.end => {
                if length > *best {
                    m.period = b.period;
                    *best = length;
                }
                m.end = m.end.max(b.end);
            }
            _ => merged.push((b, length)),
        }
    }
    let mut blocks: Vec<Block> = merged.into_iter().map(|(b, _)| b).collect();
    if circular {
        // A block found again one turn later is the same block.
        blocks.retain(|b| b.end - b.start < n);
        if blocks.len() > 1 {
            let last = blocks.last().unwrap().clone();
            if last.end > n && last.end - n > blocks[0].start {
                blocks.remove(0);
            }
        }
    }
    blocks
}

struct Anchors {
    left: Vec<u8>,
    right: Vec<u8>,
    /// Bases between the anchors in the consensus.
    gap: usize,
    /// 0-based consensus position of the first left-anchor base.
    left_start: usize,
}

fn occurrences(text: &[u8], needle: &[u8]) -> usize {
    text.windows(needle.len()).filter(|w| *w == needle).count()
}

/// The nearest unique, fully called 12-mers outside a block, on both strands.
fn anchors(sequence: &[u8], circular: bool, block: &Block) -> Option<Anchors> {
    let n = sequence.len();
    let at = |i: usize| -> Option<u8> {
        if circular {
            Some(sequence[i % n])
        } else {
            sequence.get(i).copied()
        }
    };
    let kmer = |start: usize| -> Option<Vec<u8>> {
        let s: Option<Vec<u8>> = (start..start + ANCHOR).map(at).collect();
        s.filter(|s| s.iter().all(|b| b"ACGT".contains(b)))
    };
    let search: Vec<u8> = if circular {
        sequence
            .iter()
            .chain(&sequence[..ANCHOR - 1])
            .copied()
            .collect()
    } else {
        sequence.to_vec()
    };
    let reverse = reverse_complement(&search);
    let unique = |k: &[u8]| occurrences(&search, k) + occurrences(&reverse, k) == 1;
    let offset = if circular { n } else { 0 };
    let left_start = (0..=ANCHOR_SEARCH).find_map(|shift| {
        let start = (block.start + offset).checked_sub(ANCHOR + shift)?;
        kmer(start).filter(|k| unique(k)).map(|k| (start, k))
    })?;
    let right_start = (0..=ANCHOR_SEARCH).find_map(|shift| {
        let start = block.end + offset + shift;
        if !circular && start + ANCHOR > n {
            return None;
        }
        kmer(start).filter(|k| unique(k)).map(|k| (start, k))
    })?;
    Some(Anchors {
        gap: right_start.0 - (left_start.0 + ANCHOR),
        left_start: left_start.0 % n.max(1),
        left: left_start.1,
        right: right_start.1,
    })
}

#[derive(Default)]
struct Evidence {
    gaps: BTreeMap<usize, usize>,
    /// Reads entering from one anchor without the other, per extension length.
    left: Vec<usize>,
    right: Vec<usize>,
}

/// The `MIN_READS`-th largest extension: supported by at least that many reads.
fn supported(extensions: &[usize]) -> Option<usize> {
    let mut sorted = extensions.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    sorted.get(MIN_READS - 1).copied()
}

/// Locate blocks, read their anchors from all supplied pairs and classify each.
/// Returns (report index, block, reason) for blocks needing an ambiguity row, and the report.
#[allow(clippy::type_complexity)]
pub(crate) fn assess(
    sequence: &[u8],
    circular: bool,
    reads: &Reads,
) -> Result<(Vec<(usize, Block, &'static str)>, Value)> {
    let candidates = blocks(sequence, circular);
    let anchored: Vec<_> = candidates
        .iter()
        .map(|b| anchors(sequence, circular, b))
        .collect();
    let mut evidence: Vec<Evidence> = candidates.iter().map(|_| Evidence::default()).collect();
    let mut max_read = 0;
    if anchored.iter().any(Option::is_some) {
        let codes: Vec<Option<[u32; 4]>> = anchored
            .iter()
            .map(|a| {
                a.as_ref().map(|a| {
                    [
                        code(&a.left),
                        code(&a.right),
                        code(&reverse_complement(&a.left)),
                        code(&reverse_complement(&a.right)),
                    ]
                })
            })
            .collect();
        // Per block: first L, last rc(L), first R, last R, first rc(R), last rc(R).
        let mut hits = vec![[None::<usize>; 6]; codes.len()];
        fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
            for read in [a, b] {
                let t = &read.sequence;
                let n = t.len();
                max_read = max_read.max(n);
                hits.iter_mut().for_each(|h| *h = [None; 6]);
                let (mut value, mut valid) = (0u32, 0usize);
                for (i, &base) in t.iter().enumerate() {
                    let Some(bits) = two_bit(base) else {
                        valid = 0;
                        continue;
                    };
                    value = ((value << 2) | bits) & MASK;
                    valid += 1;
                    if valid < ANCHOR {
                        continue;
                    }
                    let at = i + 1 - ANCHOR;
                    for (h, c) in hits.iter_mut().zip(&codes) {
                        let Some(c) = c else { continue };
                        if value == c[0] {
                            h[0].get_or_insert(at);
                        }
                        if value == c[2] {
                            h[1] = Some(at);
                        }
                        if value == c[1] {
                            h[2].get_or_insert(at);
                            h[3] = Some(at);
                        }
                        if value == c[3] {
                            h[4].get_or_insert(at);
                            h[5] = Some(at);
                        }
                    }
                }
                for (h, evidence) in hits.iter().zip(&mut evidence) {
                    // Orient the read so it runs left to right through the block; in
                    // a reverse read, rc(X) at p is X at n - p - ANCHOR.
                    if let Some(l) = h[0] {
                        match h[3].filter(|&r| r >= l + ANCHOR) {
                            Some(r) => *evidence.gaps.entry(r - l - ANCHOR).or_default() += 1,
                            None => evidence.left.push(n - l - ANCHOR),
                        }
                    } else if let Some(l) = h[1] {
                        match h[4].filter(|&r| r + ANCHOR <= l) {
                            Some(r) => *evidence.gaps.entry(l - r - ANCHOR).or_default() += 1,
                            None => evidence.left.push(l),
                        }
                    } else if let Some(r) = h[2] {
                        evidence.right.push(r);
                    } else if let Some(r) = h[5] {
                        evidence.right.push(n - r - ANCHOR);
                    }
                }
            }
            Ok(())
        })?;
    }
    let mut flagged = Vec::new();
    let mut reports = Vec::new();
    for (index, ((block, anchors), evidence)) in
        candidates.iter().zip(&anchored).zip(&evidence).enumerate()
    {
        let n = sequence.len();
        let mut report = json!({
            "start":block.start % n + 1,"end":(block.end - 1) % n + 1,"length":block.end - block.start,
            "period":block.period,"crosses_origin":block.end > n,
        });
        let Some(anchors) = anchors else {
            report["status"] = json!("anchors_not_unique");
            report["length_measured"] = json!(false);
            reports.push(report);
            flagged.push((index, block.clone(), UNMEASURED));
            continue;
        };
        let spanning: usize = evidence.gaps.values().sum();
        let modal = evidence
            .gaps
            .iter()
            .max_by_key(|(gap, reads)| (**reads, std::cmp::Reverse(**gap)))
            .map(|(gap, reads)| (*gap, *reads));
        // Reads with room to span the modal gap but lacking the far anchor argue
        // against it; spanning must be the usual outcome, not a rare artefact.
        let measured = modal.filter(|&(gap, reads)| {
            let room = evidence.left.iter().filter(|&&e| e >= gap + ANCHOR).count();
            reads >= MIN_READS && 2 * spanning >= spanning + room
        });
        // A read that enters by `e` bases without the far anchor shows gap > e - ANCHOR.
        let lower_bound = supported(&evidence.left)
            .max(supported(&evidence.right))
            .map(|e| e.saturating_sub(ANCHOR) + 1);
        let status = match measured {
            Some((gap, _)) if gap == anchors.gap => "spanned",
            Some(_) => CONFLICT,
            None => UNMEASURED,
        };
        report["status"] = json!(status);
        report["length_measured"] = json!(measured.is_some());
        report["anchors"] = json!({"left":String::from_utf8_lossy(&anchors.left),"right":String::from_utf8_lossy(&anchors.right),
            "left_start":anchors.left_start + 1,"consensus_gap":anchors.gap});
        report["spanning_reads"] = json!(spanning);
        report["spanning_gaps"] = json!(evidence
            .gaps
            .iter()
            .map(|(g, r)| json!({"gap":g,"reads":r}))
            .collect::<Vec<_>>());
        report["measured_gap"] = json!(measured.map(|(gap, _)| gap));
        report["entering_reads"] = json!({"left":evidence.left.len(),"right":evidence.right.len()});
        report["gap_lower_bound"] = json!(if measured.is_some() {
            None
        } else {
            lower_bound
        });
        reports.push(report);
        if status != "spanned" {
            flagged.push((index, block.clone(), status));
        }
    }
    let report = json!({"status":"assessed","method":"period 2-75 self-similarity blocks of at least 60 bases in the exported consensus; unique 12-mer anchors; reads carrying both anchors measure the gap",
        "measurement_rule":"at least 3 reads share the modal gap and spanning reads are at least half of the reads with room to span it",
        "lower_bound_rule":"third-longest one-sided anchored entry minus the anchor length, plus one",
        "copy_number":"not inferred from depth","max_read_length":max_read,"blocks":reports});
    Ok((flagged, report))
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn finds_tandem_blocks_but_not_unique_sequence() {
        let unit = b"TGAATATAAATATATATCATATTTATTATATAAATATGATAT";
        let flank = dna(400, 7);
        let sequence = [&flank[..200], &unit.repeat(5), &flank[200..]].concat();
        let found = blocks(&sequence, false);
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].period, 42);
        assert!(found[0].start <= 200 && found[0].start + 5 >= 200);
        assert!(found[0].end >= 410 && found[0].end <= 430, "{found:?}");
        assert!(blocks(&dna(2000, 11), false).is_empty());
        // Circular sequences report a block crossing the origin once.
        let rotated = [&sequence[300..], &sequence[..300]].concat();
        let crossing = blocks(&rotated, true);
        assert_eq!(crossing.len(), 1, "{crossing:?}");
        assert!(crossing[0].end > rotated.len());
    }

    /// Tile reads over `sample` (circular), alternating strands; mates are tiles too.
    fn write_reads(dir: &std::path::Path, sample: &[u8], read_length: usize) -> Reads {
        let reads = Reads {
            r1: dir.join("r1.fq"),
            r2: dir.join("r2.fq"),
        };
        let tile = |start: usize| -> Vec<u8> {
            (0..read_length)
                .map(|i| sample[(start + i) % sample.len()])
                .collect()
        };
        let (mut one, mut two) = (String::new(), String::new());
        for (i, start) in (0..sample.len()).step_by(2).enumerate() {
            let mut a = tile(start);
            if i % 2 == 1 {
                a = reverse_complement(&a);
            }
            let b = reverse_complement(&tile(start + 37));
            for (out, seq, mate) in [(&mut one, a, 1), (&mut two, b, 2)] {
                out.push_str(&format!(
                    "@r{i}/{mate}\n{}\n+\n{}\n",
                    String::from_utf8_lossy(&seq),
                    "I".repeat(seq.len())
                ));
            }
        }
        std::fs::write(&reads.r1, one).unwrap();
        std::fs::write(&reads.r2, two).unwrap();
        reads
    }

    #[test]
    fn reads_measure_spanned_blocks_and_bound_unspanned_ones() {
        let dir = std::env::temp_dir().join(format!("mitovarr-tandem-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let unit = b"TGAATATAAATATATATCATATTTATTATATAAATATGATAT";
        let flank = dna(600, 23);
        let with = |copies: usize| [&flank[..300], &unit.repeat(copies), &flank[300..]].concat();
        let consensus = with(3);
        let block = blocks(&consensus, true);
        assert_eq!(block.len(), 1);

        // Reads longer than the block carry both anchors: measured, no row.
        let (flagged, report) =
            assess(&consensus, true, &write_reads(&dir, &consensus, 220)).unwrap();
        assert!(flagged.is_empty(), "{report}");
        let b = &report["blocks"][0];
        assert_eq!(b["status"], "spanned");
        assert_eq!(b["measured_gap"], b["anchors"]["consensus_gap"]);

        // Short reads only enter the block: unmeasured with a read-derived bound.
        let (flagged, report) =
            assess(&consensus, true, &write_reads(&dir, &consensus, 90)).unwrap();
        assert_eq!(flagged.len(), 1);
        assert_eq!(flagged[0].2, UNMEASURED);
        let b = &report["blocks"][0];
        assert_eq!(b["spanning_reads"], 0);
        assert_eq!(b["length_measured"], false);
        // Entries reach at most 90 - 12 = 78 bases (a bound of 67); the third-longest
        // supported entry gives a bound that never exceeds the true gap.
        let bound = b["gap_lower_bound"].as_u64().unwrap();
        assert!((60..=67).contains(&bound), "{bound}");
        assert!(bound < b["anchors"]["consensus_gap"].as_u64().unwrap());
        assert!(b["entering_reads"]["left"].as_u64().unwrap() > 0);
        assert!(b["entering_reads"]["right"].as_u64().unwrap() > 0);

        // A sample with one unit fewer is spanned at a different gap: conflict.
        let (flagged, report) =
            assess(&consensus, true, &write_reads(&dir, &with(2), 220)).unwrap();
        assert_eq!(flagged.len(), 1);
        assert_eq!(flagged[0].2, CONFLICT);
        let b = &report["blocks"][0];
        assert_eq!(
            b["measured_gap"].as_u64().unwrap() + 42,
            b["anchors"]["consensus_gap"].as_u64().unwrap()
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
