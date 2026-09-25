//! Carry existing assembly evidence through annotation's coordinate transform.
//! This module neither selects paths nor infers new structural variants.
use super::{ncr::NcrAnalysis, Prediction, Record, Rotation, Segment};
use anyhow::{ensure, Context, Result};
use serde::Serialize;
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

pub(super) struct Input {
    contig: String,
    segment: Segment,
    evidence: String,
    reason: String,
    alternatives: usize,
}

#[derive(Debug, Serialize)]
pub(super) struct Region {
    pub id: String,
    pub evidence: String,
    pub reason: String,
    pub alternatives: usize,
    pub input_segments: Vec<Segment>,
    pub segments: Vec<Segment>,
    pub overlapping_genes: Vec<String>,
    pub overlapping_ncrs: Vec<String>,
}

#[derive(Debug, Serialize)]
pub(super) struct LocusEvidence {
    pub locus: String,
    pub annotation_status: String,
    pub evidence_status: &'static str,
    pub region_ids: Vec<String>,
}

#[derive(Debug, Serialize)]
pub(super) struct Analysis {
    pub coordinate_system: &'static str,
    pub structural_resolution_from_reads: &'static str,
    pub regions: Vec<Region>,
    pub ncr_loci: Vec<LocusEvidence>,
}

pub(super) fn read(path: &Path, records: &[Record]) -> Result<Vec<Input>> {
    let mut lines = BufReader::new(File::open(path)?).lines();
    let header = lines.next().context("empty assembly ambiguity TSV")??;
    let expected = [
        "sequence",
        "start",
        "end",
        "evidence",
        "reason",
        "alternatives",
    ];
    ensure!(header.split('\t').take(6).eq(expected),
        "--assembly-ambiguity requires primary_ambiguity.tsv, structural_ambiguity.tsv or annotation assembly_ambiguity.tsv in input FASTA coordinates; mapping ambiguity.tsv uses reference coordinates");
    let mut result = Vec::new();
    for (line_no, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let columns: Vec<_> = line.split('\t').collect();
        ensure!(
            columns.len() >= 6,
            "assembly ambiguity row {} has fewer than six columns",
            line_no + 2
        );
        let record = records
            .iter()
            .find(|r| r.name == columns[0])
            .with_context(|| {
                format!(
                    "assembly ambiguity names unknown input sequence {}",
                    columns[0]
                )
            })?;
        let start: usize = columns[1].parse().context("invalid ambiguity start")?;
        let end: usize = columns[2].parse().context("invalid ambiguity end")?;
        ensure!(
            start > 0 && end >= start && end <= record.sequence.len(),
            "assembly ambiguity coordinates outside {} (length {})",
            record.name,
            record.sequence.len()
        );
        result.push(Input {
            contig: record.name.clone(),
            segment: Segment { start, end },
            evidence: columns[3].into(),
            reason: columns[4].into(),
            alternatives: columns[5].parse()?,
        });
    }
    Ok(result)
}

fn merge(mut segments: Vec<Segment>) -> Vec<Segment> {
    segments.sort_by_key(|s| (s.start, s.end));
    let mut result: Vec<Segment> = Vec::new();
    for s in segments {
        if let Some(last) = result.last_mut().filter(|last| s.start <= last.end + 1) {
            last.end = last.end.max(s.end);
        } else {
            result.push(s);
        }
    }
    result
}

fn transform(s: &Segment, length: usize, rotation: &Rotation) -> Vec<Segment> {
    let start = if rotation.reverse_complemented {
        length - s.end
    } else {
        s.start - 1
    };
    let start = (start + length - rotation.offset_0based) % length;
    let end = start + s.end - s.start + 1;
    if end <= length {
        vec![Segment {
            start: start + 1,
            end,
        }]
    } else {
        vec![
            Segment {
                start: start + 1,
                end: length,
            },
            Segment {
                start: 1,
                end: end - length,
            },
        ]
    }
}

fn overlaps(a: &[Segment], b: &[Segment]) -> bool {
    a.iter()
        .any(|x| b.iter().any(|y| x.start <= y.end && y.start <= x.end))
}

pub(super) fn project(
    inputs: &[Input],
    record: &Record,
    record_index: usize,
    rotation: &Rotation,
    prediction: &Prediction,
    ncr: &NcrAnalysis,
) -> Analysis {
    let mut groups = BTreeMap::<(String, String, usize), Vec<Segment>>::new();
    for i in inputs.iter().filter(|i| i.contig == record.name) {
        groups
            .entry((i.evidence.clone(), i.reason.clone(), i.alternatives))
            .or_default()
            .push(i.segment.clone());
    }
    let regions: Vec<_> = groups
        .into_iter()
        .enumerate()
        .map(|(i, ((evidence, reason, alternatives), input))| {
            let input_segments = merge(input);
            let segments = merge(
                input_segments
                    .iter()
                    .flat_map(|s| transform(s, record.sequence.len(), rotation))
                    .collect(),
            );
            let mut overlapping_genes: Vec<_> = prediction
                .features
                .iter()
                .filter(|f| overlaps(&segments, &f.segments))
                .map(|f| f.name.clone())
                .collect();
            overlapping_genes.sort();
            overlapping_genes.dedup();
            let mut overlapping_ncrs: Vec<_> = ncr
                .loci
                .iter()
                .filter(|l| overlaps(&segments, &l.evidence_segments))
                .map(|l| l.name.clone())
                .collect();
            overlapping_ncrs.sort();
            overlapping_ncrs.dedup();
            Region {
                id: format!("contig{}_assembly_region{}", record_index + 1, i + 1),
                evidence,
                reason,
                alternatives,
                input_segments,
                segments,
                overlapping_genes,
                overlapping_ncrs,
            }
        })
        .collect();
    let ncr_loci = ncr
        .loci
        .iter()
        .map(|l| {
            let matched: Vec<_> = regions
                .iter()
                .filter(|r| r.overlapping_ncrs.contains(&l.name))
                .collect();
            let unknown = ncr
                .regions
                .iter()
                .any(|r| l.region_ids.contains(&r.id) && r.ambiguous_bases > 0);
            let evidence_status = if l.status == "unresolved_contig_ends" {
                "open_contig_end"
            } else if unknown {
                "unknown_sequence"
            } else if matched.iter().any(|r| {
                matches!(
                    r.reason.as_str(),
                    "competing_cycle_paths" | "majority_path_structure_unresolved"
                )
            }) {
                "unresolved_path_alternatives"
            } else if !matched.is_empty() {
                "retained_assembly_alternatives"
            } else {
                "no_reported_path_ambiguity"
            };
            LocusEvidence {
                locus: l.name.clone(),
                annotation_status: l.status.clone(),
                evidence_status,
                region_ids: matched.iter().map(|r| r.id.clone()).collect(),
            }
        })
        .collect();
    Analysis {
        coordinate_system: "one-based inclusive in annotated.fasta",
        structural_resolution_from_reads: "not_certified_by_coordinate_projection",
        regions,
        ncr_loci,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn interval_transform_preserves_every_base_on_both_strands_across_origin() {
        for reverse in [false, true] {
            for offset in 0..20 {
                let rotation = Rotation {
                    requested: true,
                    applied: true,
                    reason: "test".into(),
                    reverse_complemented: reverse,
                    offset_0based: offset,
                    input_anchor_1based: None,
                };
                for (start, end) in [(1, 1), (1, 20), (3, 9), (19, 20)] {
                    let s = Segment { start, end };
                    let actual: std::collections::BTreeSet<_> = transform(&s, 20, &rotation)
                        .iter()
                        .flat_map(|s| s.start..=s.end)
                        .collect();
                    let expected = (start - 1..end)
                        .map(|p| {
                            let p = if reverse { 19 - p } else { p };
                            (p + 20 - offset) % 20 + 1
                        })
                        .collect();
                    assert_eq!(actual, expected);
                }
            }
        }
    }
}
