//! Gene-level triage of a frozen sequence. NCR completeness is not a gate here;
//! these checks do not certify read support, repeat structure or circularity.
use super::{uncertainty::Analysis, Panel, Prediction};
use serde::Serialize;

#[derive(Debug, Serialize)]
pub(super) struct GeneQuality {
    pub gene: String,
    pub kind: String,
    pub status: &'static str,
    pub placements: usize,
    pub full_length_placements: usize,
    pub feature_ids: Vec<String>,
    pub assembly_region_ids: Vec<String>,
    pub issues: Vec<&'static str>,
}

#[derive(Debug, Serialize)]
pub(super) struct Summary {
    pub status: &'static str,
    pub ncr1_required: bool,
    pub assembly_ambiguity_supplied: bool,
    pub read_validation: &'static str,
    pub expected_models: usize,
    pub annotation_clear_models: Vec<String>,
    pub review_required_models: Vec<String>,
    pub genes: Vec<GeneQuality>,
}

pub(super) fn assess(
    panel: &Panel,
    prediction: &Prediction,
    evidence: Option<&Analysis>,
) -> Summary {
    let genes: Vec<_> = panel
        .genes
        .iter()
        .map(|model| {
            let hits: Vec<_> = prediction
                .features
                .iter()
                .filter(|f| f.name == model.name)
                .collect();
            let full_length = hits.iter().filter(|f| f.full_length_supported).count();
            let assembly_region_ids: Vec<_> = evidence
                .into_iter()
                .flat_map(|e| &e.regions)
                .filter(|r| r.overlapping_genes.contains(&model.name))
                .map(|r| r.id.clone())
                .collect();
            let mut issues = Vec::new();
            if hits.is_empty() {
                issues.push("not_detected");
            } else if full_length == 0 {
                issues.push("no_full_length_supported_placement");
            }
            // An extra partial homolog is not automatically a second functional
            // gene, but still needs review before assigning reads to this locus.
            if hits.len() > 1 {
                issues.push("multiple_homologous_placements");
            }
            if prediction.search_limited_models.contains(&model.name) {
                issues.push("annotation_search_limited");
            }
            if hits.iter().any(|f| f.ambiguous_bases > 0) {
                issues.push("ambiguous_sequence");
            }
            if hits.iter().any(|f| !f.competing_models.is_empty()) {
                issues.push("competing_gene_identity");
            }
            if hits.iter().any(|f| {
                f.cds.as_ref().is_some_and(|c| {
                    !c.valid_start
                        || matches!(c.stop.as_str(), "missing" | "partial_or_uncertain")
                        || !c.frame_valid
                        || !c.internal_stop_codons.is_empty()
                })
            }) {
                issues.push("cds_boundary_or_frame_issue");
            }
            if !assembly_region_ids.is_empty() {
                issues.push("overlaps_assembly_ambiguity");
            }
            GeneQuality {
                gene: model.name.clone(),
                kind: model.kind.clone(),
                status: if issues.is_empty() {
                    "annotation_clear"
                } else {
                    "review_required"
                },
                placements: hits.len(),
                full_length_placements: full_length,
                feature_ids: hits.iter().map(|f| f.id.clone()).collect(),
                assembly_region_ids,
                issues,
            }
        })
        .collect();
    Summary {
        status: if genes.iter().all(|g| g.issues.is_empty()) {
            "no_detected_gene_issues"
        } else {
            "gene_review_required"
        },
        ncr1_required: false,
        assembly_ambiguity_supplied: evidence.is_some(),
        read_validation: "not_performed",
        expected_models: panel.genes.len(),
        annotation_clear_models: genes
            .iter()
            .filter(|g| g.issues.is_empty())
            .map(|g| g.gene.clone())
            .collect(),
        review_required_models: genes
            .iter()
            .filter(|g| !g.issues.is_empty())
            .map(|g| g.gene.clone())
            .collect(),
        genes,
    }
}
