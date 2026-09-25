//! Synchronous iterative recruitment. Only sample reads extend the frontier.
use crate::{
    cli::{BaitSettings, Limits, Reads},
    io::fastq,
    kmer,
};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};
use tracing::info;

pub struct Recruitment {
    pub reads: Reads,
    pub complete: bool,
    pub report: Value,
}

pub fn recruit(
    reads: &Reads,
    seed: &Path,
    out: &Path,
    settings: &BaitSettings,
    limits: &Limits,
) -> Result<Recruitment> {
    ensure!(
        (15..=31).contains(&settings.bait_k),
        "bait k must be in 15..=31"
    );
    ensure!(
        settings.max_rounds > 0 && settings.max_pairs > 0 && limits.max_kmers > 0,
        "recruitment limits must be positive"
    );
    crate::prepare_output(out)?;
    let mut keys = HashSet::new();
    let mut reader = noodles_fasta::io::Reader::new(fastq::open(seed)?);
    let mut seed_records = 0;
    for record in reader.records() {
        let record = record?;
        ensure!(
            record.sequence().len() <= 1_000_000,
            "seed record exceeds 1 Mb; provide mitochondrial seed sequences"
        );
        let mut exceeded = false;
        kmer::visit(
            record.sequence().as_ref(),
            None,
            settings.bait_k,
            0,
            |_, code| {
                if keys.len() < limits.max_kmers {
                    keys.insert(kmer::canonical(code, settings.bait_k));
                } else if !keys.contains(&kmer::canonical(code, settings.bait_k)) {
                    exceeded = true;
                }
            },
        );
        ensure!(!exceeded, "seed exceeds --max-kmers");
        seed_records += 1;
    }
    ensure!(
        !keys.is_empty(),
        "seed contains no usable A/C/G/T bait kmers"
    );
    let mut selected = HashSet::new();
    let recruited = Reads {
        r1: out.join("reads.R1.fastq"),
        r2: out.join("reads.R2.fastq"),
    };
    let (mut one, mut two) = (
        BufWriter::new(File::create(&recruited.r1)?),
        BufWriter::new(File::create(&recruited.r2)?),
    );

    let mut originals = None;
    let mut rounds = Vec::new();
    let mut complete = false;
    let mut reason = "max_rounds";
    for round in 1..=settings.max_rounds {
        let mut additions = HashSet::new();
        let mut new_pairs = 0;
        let mut limited = false;
        let scan = fastq::scan_pairs(&reads.r1, &reads.r2, |ordinal, a, b| {
            if selected.contains(&ordinal) || limited {
                return Ok(());
            }
            let mut matches = false;
            for read in [a, b] {
                kmer::visit(
                    &read.sequence,
                    Some(&read.quality),
                    settings.bait_k,
                    limits.min_baseq,
                    |_, code| {
                        matches |= keys.contains(&kmer::canonical(code, settings.bait_k));
                    },
                );
            }
            if matches {
                if selected.len() >= settings.max_pairs {
                    limited = true;
                    reason = "max_pairs";
                    return Ok(());
                }
                selected.insert(ordinal);
                a.write(&mut one)?;
                b.write(&mut two)?;
                new_pairs += 1;
                for read in [a, b] {
                    kmer::visit(
                        &read.sequence,
                        Some(&read.quality),
                        settings.bait_k,
                        limits.min_baseq,
                        |_, code| {
                            let code = kmer::canonical(code, settings.bait_k);
                            if !keys.contains(&code) && !additions.contains(&code) {
                                if keys.len() + additions.len() >= limits.max_kmers {
                                    limited = true;
                                    reason = "max_kmers";
                                } else {
                                    additions.insert(code);
                                }
                            }
                        },
                    );
                }
            }
            Ok(())
        })?;
        if originals.is_none() {
            originals = Some(scan);
        }
        rounds.push(json!({"round":round,"new_pairs":new_pairs,"total_pairs":selected.len(),"new_kmers":additions.len()}));
        info!(
            round,
            new_pairs,
            selected_pairs = selected.len(),
            "recruitment pass"
        );
        if limited {
            break;
        }
        keys.extend(additions);
        if new_pairs == 0 {
            complete = true;
            reason = if selected.is_empty() {
                "no_matching_pairs"
            } else {
                "converged"
            };
            break;
        }
    }
    let originals = originals.expect("positive round limit");
    one.flush()?;
    two.flush()?;
    let usable = complete && !selected.is_empty();
    let report = json!({
        "schema":"mitovar-rust.recruitment.v1", "backend":"rust-canonical-kmer-recruitment",
        "status":if usable {"complete"} else {"unresolved"}, "complete":complete,
        "usable_for_assembly":usable, "reason":reason,
        "seed":seed.canonicalize()?, "seed_records":seed_records,
        "seed_usage":"recruitment only; no seed sequence is passed to assembly",
        "original_r1":reads.r1.canonicalize()?, "original_r2":reads.r2.canonicalize()?, "originals":originals,
        "parameters":{"bait_k":settings.bait_k,"min_baseq":limits.min_baseq,"max_rounds":settings.max_rounds,"max_pairs":settings.max_pairs,"max_kmers":limits.max_kmers},
        "selected_pairs":selected.len(), "rounds":rounds,
        "r1":recruited.r1,"r2":recruited.r2,
        "records_preserved":"original header, sequence, plus header and qualities; line endings normalized to LF",
    });
    crate::write_json(&out.join("recruitment.json"), &report)?;
    Ok(Recruitment {
        reads: recruited,
        complete: usable,
        report,
    })
}
