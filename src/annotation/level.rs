//! A reporting level for a frozen de novo primary: complete outside NCR1.
//! It reads the assembly report and the projected assembly intervals. It never
//! changes a sequence, an assembly status or the strict circle criteria, and
//! it never states an NCR1 length.
use super::{
    ncr::NcrAnalysis, quality::Summary, uncertainty::Analysis, Prediction, Record, Segment,
};
use serde::Serialize;
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};

pub(super) const LEVEL: &str = "complete_except_ncr1";
const NCR1: &str = "NCR1";

/// Alternatives kept on record: collapsed bubbles, homopolymer families,
/// single-site majority choices and minor routes the strict rules rejected.
/// They are counted and reported; they do not remove the level.
const SITE_REASONS: [&str; 5] = [
    "collapsed_bubble_alternative",
    "homopolymer_length_family",
    "majority_path_site_unresolved",
    "minor_cycle_path_retained",
    "minor_path_retained",
];
const OPEN_END: &str = "open_graph_end";
const LOW_DEPTH: &str = "low_matching_read_depth";

const RULE: &str = "one contig (a circle, or a linear primary whose outermost genes counted as present are tRNA-Ala and tRNA-Pro with NCR1 at its ends); every expected gene present once and full length, an RNA with both boundaries aligned and a length differing from its model counting as present; no mixed-sample signal or structural mixture; no interval outside NCR1 other than alternatives kept on record (collapsed bubbles, homopolymer families, single-site majority choices, minor routes the strict rules rejected); the primary read-validated with no base outside NCR1 below the depth minimum";
const SCOPE: &str = "a statement about the frozen primary outside NCR1 only; closure is reported beside it and does not decide it; NCR1 length and repeat copy number are not determined; the strict final status is unchanged";

#[derive(Debug, Serialize)]
pub(super) struct ClaimLevel {
    pub level: &'static str,
    /// `reached`, `not_reached` or `not_assessed`.
    pub status: &'static str,
    /// `circle`, `joined_by_reads` or `not_observed`; None when not assessed.
    pub closure: Option<String>,
    pub strict_final_published: Option<bool>,
    pub statement: String,
    pub not_reached: Vec<String>,
    pub not_assessed_reason: Option<&'static str>,
    pub rna_length_variants: Vec<String>,
    pub flagged_sites_outside_ncr1: usize,
    pub flagged_sites_inside_ncr1: usize,
    pub flagged_site_reasons_outside_ncr1: BTreeMap<String, usize>,
    pub collapsed_bubbles_not_placed: Option<usize>,
    /// Assembled NCR1 bases on the primary: two flanks of a linear contig, or
    /// the spelled region of a circle. Never a measured NCR1 length.
    pub ncr1_assembled_bases: Vec<usize>,
    pub ncr1_length: &'static str,
    pub rule: &'static str,
    pub scope: &'static str,
}

pub(super) struct Inputs<'a> {
    pub record: &'a Record,
    pub circular: bool,
    pub prediction: &'a Prediction,
    pub ncr: &'a NcrAnalysis,
    pub evidence: Option<&'a Analysis>,
    pub quality: &'a Summary,
    /// The de novo `report.json` that wrote the input primary, when supplied.
    pub assembly: Option<&'a Value>,
    /// This record's rows of the supplied ambiguity table, as read.
    pub table_rows: Option<Vec<IntervalRow>>,
}

/// One uncertainty interval: one-based inclusive start and end, evidence,
/// reason and alternative count, as in `primary_ambiguity.tsv`.
pub(super) type IntervalRow = (usize, usize, String, String, usize);

/// The report's own `primary.uncertainty_intervals`, in table form.
fn report_rows(primary: &Value) -> Option<Vec<IntervalRow>> {
    let mut rows = primary["uncertainty_intervals"]
        .as_array()?
        .iter()
        .map(|r| {
            Some((
                usize::try_from(r["start"].as_u64()?).ok()?,
                usize::try_from(r["end"].as_u64()?).ok()?,
                r["evidence"].as_str()?.to_owned(),
                r["reason"].as_str()?.to_owned(),
                usize::try_from(r["alternatives"].as_u64()?).ok()?,
            ))
        })
        .collect::<Option<Vec<_>>>()?;
    rows.sort_unstable();
    Some(rows)
}

fn not_assessed(reason: &'static str) -> ClaimLevel {
    ClaimLevel {
        level: LEVEL,
        status: "not_assessed",
        closure: None,
        strict_final_published: None,
        statement: format!("not assessed: {reason}"),
        not_reached: Vec::new(),
        not_assessed_reason: Some(reason),
        rna_length_variants: Vec::new(),
        flagged_sites_outside_ncr1: 0,
        flagged_sites_inside_ncr1: 0,
        flagged_site_reasons_outside_ncr1: BTreeMap::new(),
        collapsed_bubbles_not_placed: None,
        ncr1_assembled_bases: Vec::new(),
        ncr1_length: "not_determined",
        rule: RULE,
        scope: SCOPE,
    }
}

fn overlaps(a: &Segment, b: &Segment) -> bool {
    a.start <= b.end && b.start <= a.end
}

/// Spelled NCR1 bases: each open-end flank of a linear contig, or one total
/// for a circle even when its NCR1 is split at the sequence origin.
fn assembled_ncr1(
    segments: &[Segment],
    circular: bool,
    in_gene: impl Fn(&Segment) -> bool,
) -> Vec<usize> {
    let lengths = segments
        .iter()
        .filter(|s| !in_gene(s))
        .map(|s| s.end - s.start + 1);
    if circular {
        vec![lengths.sum()]
    } else {
        lengths.collect()
    }
}

pub(super) fn assess(input: Inputs) -> ClaimLevel {
    let Some(assembly) = input.assembly else {
        return not_assessed("assembly_report_not_supplied");
    };
    if assembly["schema"] != "mitovar-rust.assembly.v1" || assembly["backend"] != "rust-de-bruijn" {
        return not_assessed("not_a_de_novo_assembly_report");
    }
    let primary = &assembly["primary"];
    if !primary.is_object() {
        return not_assessed("assembly_has_no_primary");
    }
    if input.record.name != "primary"
        || primary["length"].as_u64() != Some(input.record.sequence.len() as u64)
    {
        return not_assessed("input_is_not_the_reported_primary");
    }
    if primary["circular"].as_bool() != Some(input.circular) {
        return not_assessed("annotation_topology_differs_from_assembly");
    }
    let Some(evidence) = input.evidence else {
        return not_assessed("assembly_ambiguity_not_supplied");
    };
    // The table must be this run's own: the level reads its intervals.
    let mut table = input.table_rows.clone().unwrap_or_default();
    table.sort_unstable();
    if report_rows(primary) != Some(table) {
        return not_assessed("assembly_ambiguity_differs_from_report");
    }
    let Some(ncr1) = input.ncr.loci.iter().find(|l| l.name == NCR1) else {
        return not_assessed("no_ncr_catalog");
    };
    let closure = if input.circular {
        "circle"
    } else {
        match primary["closure"]["state"].as_str() {
            Some("joined_by_reads") => "joined_by_reads",
            Some("not_observed") => "not_observed",
            Some("not_assessed") => return not_assessed("closure_not_assessed"),
            _ => return not_assessed("assembly_report_lacks_closure"),
        }
    };
    match primary["read_validation"]["status"].as_str() {
        Some("performed") => {}
        Some("resource_limited") => return not_assessed("primary_read_validation_incomplete"),
        _ => return not_assessed("assembly_report_lacks_primary_read_validation"),
    }

    let features = &input.prediction.features;
    let in_gene = |s: &Segment| {
        features
            .iter()
            .any(|f| f.segments.iter().any(|g| overlaps(g, s)))
    };
    let in_ncr1 =
        |s: &Segment| !in_gene(s) && ncr1.evidence_segments.iter().any(|g| overlaps(g, s));
    let mut not_reached = BTreeSet::new();

    // Genes (rule 2), first: rule 1 uses the RNA length variants it accepts.
    // Overlap with an assembly interval is judged below, by reason.
    let mut rna_length_variants = Vec::new();
    for gene in &input.quality.genes {
        let issues: Vec<_> = gene
            .issues
            .iter()
            .copied()
            .filter(|issue| *issue != "overlaps_assembly_ambiguity")
            .collect();
        if gene.placements == 1 && gene.full_length_placements == 1 && issues.is_empty() {
            continue;
        }
        let mut hits = features.iter().filter(|f| f.name == gene.gene);
        let length_variant = matches!(
            (hits.next(), hits.next()),
            (Some(hit), None) if hit.status == "boundary_complete_length_variant"
        ) && matches!(gene.kind.as_str(), "tRNA" | "rRNA")
            && issues == ["no_full_length_supported_placement"];
        if length_variant {
            rna_length_variants.push(gene.gene.clone());
        } else {
            not_reached.insert(format!("gene:{}", gene.gene));
        }
    }

    // 1. One contig. A circle qualifies as it is; a linear primary must end in
    // NCR1 on both sides: of the genes counted as present (full length, or an
    // accepted RNA length variant), the one starting first and the one ending
    // last are tRNA-Ala and tRNA-Pro.
    if !input.circular {
        let present: Vec<_> = features
            .iter()
            .filter(|f| {
                f.full_length_supported
                    || (f.status == "boundary_complete_length_variant"
                        && rna_length_variants.contains(&f.name))
            })
            .filter_map(|f| {
                let start = f.segments.iter().map(|s| s.start).min()?;
                let end = f.segments.iter().map(|s| s.end).max()?;
                Some((start, end, f.name.as_str()))
            })
            .collect();
        let first = present
            .iter()
            .min_by_key(|(start, end, name)| (*start, *end, *name));
        let last = present
            .iter()
            .max_by_key(|(start, end, name)| (*end, *start, std::cmp::Reverse(*name)));
        let ends: BTreeSet<_> = [first, last]
            .into_iter()
            .flatten()
            .map(|(_, _, name)| *name)
            .collect();
        if ncr1.status != "unresolved_contig_ends"
            || ends != BTreeSet::from(["tRNA-Ala", "tRNA-Pro"])
        {
            not_reached.insert("not_one_contig_ending_in_ncr1".to_owned());
        }
    }

    // 3. One sample.
    if assembly["mixed_sample"]["status"] == "mixed_sample_signal" {
        not_reached.insert("mixed_sample_signal".to_owned());
    }
    if assembly["mixed_sample"]["structural_alternatives"]
        .as_array()
        .is_some_and(|events| !events.is_empty())
    {
        not_reached.insert("structural_mixture".to_owned());
    }

    // 4. Assembly intervals, each segment placed on its own. Alternatives
    // kept on record are counted; any other interval outside NCR1 is open
    // structure or missing read support.
    let (mut sites_outside, mut sites_inside) = (0, 0);
    let mut site_reasons = BTreeMap::<String, usize>::new();
    for region in &evidence.regions {
        if region.reason == OPEN_END {
            continue;
        }
        let outside = region.segments.iter().any(|s| !in_ncr1(s));
        if SITE_REASONS.contains(&region.reason.as_str()) {
            if outside {
                sites_outside += 1;
                *site_reasons.entry(region.reason.clone()).or_default() += 1;
            } else {
                sites_inside += 1;
            }
        } else if outside {
            not_reached.insert(if region.reason == LOW_DEPTH {
                "low_read_depth_outside_ncr1".to_owned()
            } else {
                format!("open_structure_outside_ncr1:{}", region.reason)
            });
        }
    }

    let strict = assembly["final_published"] == true && primary["final_verified"] == true;
    let reached = not_reached.is_empty();
    let statement = if !reached {
        "not complete outside NCR1: see not_reached".to_owned()
    } else if strict {
        "verified circle: the strict final criteria passed; also complete outside NCR1".to_owned()
    } else {
        match closure {
            "circle" => "complete outside NCR1; circular graph path or cycle candidate, not a verified circle; NCR1 length not determined",
            "joined_by_reads" => "complete outside NCR1; overlapping reads join the two NCR1 ends; the join is not applied and NCR1 length is not determined; not a verified circle",
            _ => "complete outside NCR1: one contig with every expected gene and both ends inside NCR1; reads do not join the ends, so the sequence is not shown to be circular",
        }
        .to_owned()
    };
    ClaimLevel {
        level: LEVEL,
        status: if reached { "reached" } else { "not_reached" },
        closure: Some(closure.to_owned()),
        strict_final_published: Some(strict),
        statement,
        not_reached: not_reached.into_iter().collect(),
        not_assessed_reason: None,
        rna_length_variants,
        flagged_sites_outside_ncr1: sites_outside,
        flagged_sites_inside_ncr1: sites_inside,
        flagged_site_reasons_outside_ncr1: site_reasons,
        collapsed_bubbles_not_placed: primary["collapsed_bubbles"]["not_placed"]
            .as_array()
            .map(Vec::len),
        ncr1_assembled_bases: assembled_ncr1(&ncr1.evidence_segments, input.circular, in_gene),
        ncr1_length: "not_determined",
        rule: RULE,
        scope: SCOPE,
    }
}
