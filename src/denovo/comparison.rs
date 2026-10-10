//! Bounded whole-sequence comparison for optional layout hypotheses. The open
//! backbone stays linear; neither a score nor a single candidate proves closure.
use super::{align, anchor, evidence, pairs};
use crate::{
    cli::{GraphSettings, Limits, Reads},
    io::fastq,
    kmer,
};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

pub(super) struct Hypothesis<'a> {
    pub sequence: &'a [u8],
    pub circular: bool,
    pub junctions: [usize; 2],
}

pub(super) fn compare(
    hypotheses: &[Hypothesis<'_>],
    flank_seeds: &HashSet<kmer::Code>,
    training_sequence: &[u8],
    reads: &Reads,
    settings: &GraphSettings,
    limits: &Limits,
    out: &Path,
) -> Result<Value> {
    let training =
        pairs::estimate_sequences(vec![training_sequence.to_vec()], reads, settings, limits)?;
    let mut seeds = flank_seeds.clone();
    for h in hypotheses {
        let mut counts = HashMap::new();
        let mut sequence = h.sequence.to_vec();
        if h.circular {
            sequence.extend(h.sequence.iter().copied().cycle().take(20));
        }
        kmer::visit(&sequence, None, 21, 0, |_, c| {
            *counts.entry(kmer::canonical(c, 21)).or_insert(0usize) += 1
        });
        seeds.retain(|c| counts.get(c) == Some(&1));
    }
    let indexes: Vec<_> = hypotheses
        .iter()
        .map(|h| align::Index::new(h.sequence, h.circular, training.max_read))
        .collect();
    let locations: Vec<_> = indexes
        .iter()
        .map(|i| anchor::locations(i, &seeds))
        .collect();
    let mut log = flate2::write::GzEncoder::new(
        BufWriter::new(File::create(
            out.join("coverage_layout_fragments.jsonl.gz"),
        )?),
        flate2::Compression::default(),
    );
    let mut unknown_reasons = BTreeMap::<&str, usize>::new();
    let n = hypotheses.len();
    let mut totals = vec![0.0; n];
    let mut junctions = vec![[0usize; 2]; n];
    let mut observed_junctions = vec![[0usize; 2]; n];
    let mut mate_brackets = vec![[0usize; 2]; n];
    let mut unplaced = vec![0usize; n];
    let mut names = HashSet::new();
    let (mut fragments, mut unknown) = (0usize, 0usize);
    let mut rescue = evidence::RescueStats::default();
    fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
        if !names.insert(a.name(1)?.to_owned()) {
            return Ok(());
        }
        ensure!(
            names.len() <= settings.max_validation_pairs,
            "layout comparison exceeds --max-validation-pairs"
        );
        let Some(pivot) = anchor::Pivot::choose([a, b], &seeds, limits.min_baseq) else {
            return Ok(());
        };
        fragments += 1;
        let anchored = [a, b][pivot.mate];
        if indexes
            .iter()
            .zip(&locations)
            .any(|(i, map)| !pivot.consistent(map, anchored.sequence.len(), i.circle_length()))
        {
            unknown += 1;
            *unknown_reasons
                .entry("incompatible_common_anchor_projections")
                .or_default() += 1;
            serde_json::to_writer(
                &mut log,
                &json!({"fragment":a.name(1)?,"pivot":pivot.report(),"unknown_reason":"incompatible_common_anchor_projections"}),
            )?;
            writeln!(log)?;
            return Ok(());
        }
        let mut one = align::Prepared::new(a);
        let mut two = align::Prepared::new(b);
        let mut placements: Vec<_> = indexes
            .iter()
            .map(|_| (align::Placements::default(), align::Placements::default()))
            .collect();
        let mut removed = Vec::new();
        for (i, (x, y)) in placements.iter_mut().enumerate() {
            let (hits, prepared) = if pivot.mate == 0 {
                (x, &mut one)
            } else {
                (y, &mut two)
            };
            removed.push(pivot.constrain(
                hits,
                anchored,
                prepared,
                &indexes[i],
                &locations[i],
                limits.min_baseq,
            ));
        }
        let rs = evidence::rescue_pool(
            &indexes,
            &mut placements,
            &mut one,
            &mut two,
            training.model.as_ref(),
            pivot.mate,
        );
        rescue.searches += rs.searches;
        rescue.seeded_decoy_searches += rs.seeded_decoy_searches;
        rescue.added_placements += rs.added_placements;
        if let Some(reason) = evidence::unknown_reason(&placements, &rs.states, pivot.mate) {
            unknown += 1;
            *unknown_reasons.entry(reason).or_default() += 1;
            serde_json::to_writer(
                &mut log,
                &json!({"fragment":a.name(1)?,"pivot":pivot.report(),"search_status":rs.states,"unknown_reason":reason}),
            )?;
            writeln!(log)?;
            return Ok(());
        }
        let mut diagnostics = Vec::new();
        let ideal = align::ideal_score(a) + align::ideal_score(b);
        for (i, ((x, y), index)) in placements.iter().zip(&indexes).enumerate() {
            unplaced[i] += usize::from(x.hits.is_empty() || y.hits.is_empty());
            let h = &hypotheses[i];
            let result = evidence::joint(
                &evidence::observations(x, a, index, limits.min_baseq),
                &evidence::observations(y, b, index, limits.min_baseq),
                training.model.as_ref(),
                ideal,
                if h.circular { h.sequence.len() } else { 0 },
                h.junctions,
                None,
            );
            diagnostics.push(json!({"candidate":i,"pivot_projection":pivot.projection(&locations[i],anchored.sequence.len()),"anchor_hits_removed":removed[i],"r1_placements":x.hits.len(),"r2_placements":y.hits.len(),"search_status":rs.states[i],"joint":result.diagnostics,"observed_junctions":result.observed_junctions,"mate_brackets":result.mate_brackets}));
            totals[i] += result.score;
            if h.circular {
                for j in 0..2 {
                    junctions[i][j] += usize::from(result.junctions[j]);
                    observed_junctions[i][j] += usize::from(result.observed_junctions[j]);
                    mate_brackets[i][j] += usize::from(result.mate_brackets[j]);
                }
            }
        }
        serde_json::to_writer(
            &mut log,
            &json!({"fragment":a.name(1)?,"pivot":pivot.report(),"unknown_reason":null,"candidates":diagnostics}),
        )?;
        writeln!(log)?;
        Ok(())
    })?;
    log.finish()?.flush()?;
    let best = (0..n)
        .max_by(|&a, &b| totals[a].total_cmp(&totals[b]).then(b.cmp(&a)))
        .unwrap();
    let best_candidates: Vec<_> = (0..n)
        .filter(|&i| (totals[i] - totals[best]).abs() <= 1e-8)
        .collect();
    let next_score = (0..n)
        .filter(|i| !best_candidates.contains(i))
        .map(|i| totals[i])
        .max_by(f64::total_cmp);
    let unknown_bound = evidence::unknown_contrast_bound(unknown);
    Ok(
        json!({"status":if fragments > unknown {"evaluated"} else {"no_scored_common_fragments"},
        "alignment_search_exhaustive":false,"alignment_domain":align::profile::report(),"unknown_score_contrast_bound":evidence::unknown_contrast_bound(unknown),"common_fragments":fragments,"unscored_fragments":unknown,"shared_unique_flank_21mers":seeds.len(),
        "best_candidate":best,"best_candidates":best_candidates,"candidates_within_5_log_units":(0..n).filter(|&i|totals[i]>=totals[best]-5.0).collect::<Vec<_>>(),"candidates_within_5_plus_unknown_bound":(0..n).filter(|&i|totals[i]>=totals[best]-5.0-unknown_bound).collect::<Vec<_>>(),"score_margin_to_next":next_score.map(|s|totals[best]-s),"measured_preference_robust_to_omitted_fragments":fragments>unknown && best_candidates.len()==1 && next_score.is_some_and(|s|totals[best]-s>unknown_bound),"unknown_reasons":unknown_reasons,"fragment_diagnostics":"coverage_layout_fragments.jsonl.gz","final_verified":false,"phase_verified":false,
        "cohort":"one named pair once; fixed mate and query offset of an observed retained-backbone flank 21-mer unique in every whole hypothesis; traceback-constrained scoring and rescue",
        "baseline_topology":"linear; no artificial closure or filled bases",
        "fragment_model":training.report,
        "mate_rescue":{"searched_hypothesis_mates":rescue.searches,"seeded_decoy_searches":rescue.seeded_decoy_searches,"added_placements":rescue.added_placements},
        "candidates":hypotheses.iter().enumerate().map(|(i,h)|json!({"candidate":i,"length":h.sequence.len(),"circular_hypothesis":h.circular,"score":totals[i],"junction_fragments":junctions[i],"observed_junction_fragments":observed_junctions[i],"mate_bracket_fragments":mate_brackets[i],"junction_positions_0based":h.junctions,"pairs_with_missing_mate_placement":unplaced[i]})).collect::<Vec<_>>(),
        "limitations":["bounded candidate subset; not an exhaustive structural comparison","mate brackets are geometry evidence, not observed bases in the gap","common-flank ascertainment; no calibrated correctness probability","scores do not establish whole-sequence depth or cross-locus phase"]}),
    )
}
