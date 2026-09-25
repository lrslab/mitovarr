//! Post-assembly, homology-assisted gene prediction. No assembly code calls this
//! module: candidates and selection are frozen before annotation is introduced.
mod align;
mod genetic_code;
pub mod ncr;
mod output;
mod uncertainty;

use crate::cli::{AnnotationArgs, Species, Topology};
use align::{Hit, Search};
use anyhow::{ensure, Context, Result};
use genetic_code::{check, is_start, is_stop, CdsCheck};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{collections::BTreeSet, fs, io::BufRead, path::Path};

#[derive(Clone, Debug, Deserialize)]
pub struct GeneModel {
    pub name: String,
    pub kind: String,
    pub sequence: String,
    #[serde(default)]
    pub product: String,
}
#[derive(Clone, Debug, Deserialize)]
pub struct Panel {
    pub schema: String,
    pub accession: String,
    pub genetic_code: u8,
    pub genes: Vec<GeneModel>,
}
impl Panel {
    pub fn bundled(species: Species) -> Result<Self> {
        Self::parse(match species {
            Species::Cel => include_str!("models/cel.json"),
            Species::Cbr => include_str!("models/cbr.json"),
            Species::Ctro => include_str!("models/ctro.json"),
        })
    }
    pub fn parse(text: &str) -> Result<Self> {
        let panel: Self = serde_json::from_str(text)?;
        ensure!(
            panel.schema == "mitovar-rust.gene-models.v1",
            "unsupported gene-model schema"
        );
        ensure!(
            panel.genetic_code == 5,
            "annotation currently implements genetic code 5 only"
        );
        ensure!(
            !panel.genes.is_empty() && panel.genes.len() <= 512,
            "panel must have 1..512 models"
        );
        let mut names = BTreeSet::new();
        for gene in &panel.genes {
            ensure!(
                names.insert(&gene.name),
                "duplicate model name: {}",
                gene.name
            );
            ensure!(
                matches!(gene.kind.as_str(), "CDS" | "rRNA" | "tRNA"),
                "unsupported model type: {}",
                gene.kind
            );
            ensure!(
                (15..=5000).contains(&gene.sequence.len())
                    && gene.sequence.bytes().all(|b| b"ACGT".contains(&b)),
                "model {} must contain 15..5000 unambiguous DNA bases",
                gene.name
            );
        }
        Ok(panel)
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Segment {
    pub start: usize,
    pub end: usize,
}
#[derive(Clone, Debug, Serialize)]
pub struct Feature {
    pub id: String,
    pub name: String,
    pub kind: String,
    pub product: String,
    pub strand: char,
    /// 1-based inclusive segments, in biological (5′ to 3′) order.
    pub segments: Vec<Segment>,
    pub length: usize,
    pub input_five_prime_1based: usize,
    pub five_prime_1based: usize,
    pub model_length: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub query_coverage: f64,
    pub identity: f64,
    pub alignment_score: i32,
    pub gap_bases: usize,
    pub five_prime_boundary: String,
    pub three_prime_boundary: String,
    pub status: String,
    pub full_length_supported: bool,
    pub ambiguous_bases: usize,
    pub competing_models: Vec<String>,
    pub rna_structure: Option<String>,
    pub cds: Option<CdsCheck>,
    #[serde(skip)]
    pub sequence: Vec<u8>,
}
impl Feature {
    fn update_segments(&mut self, genome_length: usize) {
        let five = self.five_prime_1based - 1;
        let low = if self.strand == '+' {
            five
        } else {
            (five + genome_length + 1 - self.length) % genome_length
        };
        self.segments = if low + self.length <= genome_length {
            vec![Segment {
                start: low + 1,
                end: low + self.length,
            }]
        } else {
            vec![
                Segment {
                    start: low + 1,
                    end: genome_length,
                },
                Segment {
                    start: 1,
                    end: (low + self.length) % genome_length,
                },
            ]
        };
        if self.strand == '-' {
            self.segments.reverse();
        }
    }
}
#[derive(Debug, Serialize)]
pub struct Prediction {
    pub features: Vec<Feature>,
    pub missing_models: Vec<String>,
    pub models_without_full_length_support: Vec<String>,
    pub search_limited_models: Vec<String>,
}
#[derive(Debug, Serialize)]
pub struct Rotation {
    pub requested: bool,
    pub applied: bool,
    pub reason: String,
    pub reverse_complemented: bool,
    /// Rotate left by this number of bases AFTER the optional reverse complement.
    pub offset_0based: usize,
    pub input_anchor_1based: Option<usize>,
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b,
    }
}
pub(super) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}
fn overlap(a: &Feature, b: &Feature) -> usize {
    a.segments
        .iter()
        .flat_map(|x| {
            b.segments
                .iter()
                .map(move |y| (x.end.min(y.end) + 1).saturating_sub(x.start.max(y.start)))
        })
        .sum()
}
fn matching_fraction(a: &[u8], b: &[u8]) -> f64 {
    a.iter()
        .zip(b)
        .filter(|(x, y)| x == y && b"ACGT".contains(x))
        .count() as f64
        / a.len().max(1) as f64
}

fn feature(model: &GeneModel, hit: Hit, target: &[u8], n: usize, reverse: bool) -> Option<Feature> {
    let query = model.sequence.as_bytes();
    let covered = hit.qend - hit.qstart;
    let fraction = covered as f64 / query.len() as f64;
    let (min_coverage, min_identity, min_length) = match model.kind.as_str() {
        "tRNA" => (0.80, 0.85, 30),
        "rRNA" => (0.20, 0.75, 80),
        _ => (0.20, 0.70, 80),
    };
    if fraction < min_coverage
        || hit.identity() < min_identity
        || covered < min_length.min(query.len())
    {
        return None;
    }
    let (mut start, mut end) = (hit.start, hit.end);
    let mut five = if hit.qstart == 0 && hit.terminal_supported(true) {
        "aligned"
    } else {
        "partial_or_uncertain"
    };
    let mut three = if hit.qend == query.len() && hit.terminal_supported(false) {
        "aligned"
    } else {
        "partial_or_uncertain"
    };
    // A clipped CDS terminus can be extended only onto sample bases with a start/
    // stop codon and nearby homology. Long missing ends cannot be manufactured.
    if model.kind == "CDS" {
        let window = query.len().min(12);
        if hit.qstart > 0 && hit.qstart <= 9 {
            if let Some(projected) = start.checked_sub(hit.qstart) {
                if projected + window <= target.len()
                    && is_start(&target[projected..projected + 3])
                    && matching_fraction(&query[..window], &target[projected..projected + window])
                        >= 0.75
                {
                    start = projected;
                    five = "orf_inferred";
                }
            }
        }
        let missing = query.len() - hit.qend;
        if missing > 0 && missing <= 9 {
            let projected = end + missing;
            if projected <= target.len()
                && projected >= window
                && is_stop(&target[projected - 3..projected])
                && matching_fraction(
                    &query[query.len() - window..],
                    &target[projected - window..projected],
                ) >= 0.75
            {
                end = projected;
                three = "orf_inferred";
            }
        }
    }
    if end <= start || end - start > n {
        return None;
    }
    let sequence = target[start..end].to_vec();
    let ambiguous_bases = sequence.iter().filter(|b| !b"ACGT".contains(b)).count();
    let complete_five = five != "partial_or_uncertain";
    let complete_three = three != "partial_or_uncertain";
    let phase = if complete_five {
        0
    } else {
        (3 - hit.qstart % 3) % 3
    };
    let cds = (model.kind == "CDS").then(|| {
        check(
            &sequence,
            phase,
            complete_five,
            complete_three,
            hit.non_triplet_gaps(),
        )
    });
    let mut status = if complete_five && complete_three {
        "complete_homology"
    } else {
        "partial"
    };
    if let Some(c) = &cds {
        if complete_five && complete_three {
            status = if !c.valid_start
                || c.stop == "missing"
                || !c.frame_valid
                || !c.internal_stop_codons.is_empty()
            {
                "disrupted"
            } else {
                "complete"
            };
        }
    } else if status == "complete_homology" && hit.gap_bases() > 0 {
        status = "boundary_complete_length_variant";
    }
    if ambiguous_bases > 0 && matches!(status, "complete" | "complete_homology") {
        status = "uncertain_bases";
    }
    if model.kind == "tRNA" && hit.identity() < 0.90 && status == "complete_homology" {
        status = "uncertain_homology";
    }
    let five_prime = if reverse {
        n - 1 - start % n
    } else {
        start % n
    };
    let mut feature = Feature {
        id: String::new(),
        name: model.name.clone(),
        kind: model.kind.clone(),
        product: model.product.clone(),
        strand: if reverse { '-' } else { '+' },
        segments: Vec::new(),
        length: end - start,
        input_five_prime_1based: five_prime + 1,
        five_prime_1based: five_prime + 1,
        model_length: query.len(),
        query_start: hit.qstart + 1,
        query_end: hit.qend,
        query_coverage: fraction,
        identity: hit.identity(),
        alignment_score: hit.score,
        gap_bases: hit.gap_bases(),
        five_prime_boundary: five.into(),
        three_prime_boundary: three.into(),
        status: status.into(),
        full_length_supported: matches!(status, "complete" | "complete_homology"),
        ambiguous_bases,
        competing_models: Vec::new(),
        rna_structure: (model.kind != "CDS").then(|| "not_assessed".into()),
        cds,
        sequence,
    };
    feature.update_segments(n);
    Some(feature)
}

pub fn predict(sequence: &[u8], panel: &Panel, circular: bool) -> Prediction {
    let n = sequence.len();
    let mut all = Vec::new();
    let mut limited = BTreeSet::new();
    for reverse in [false, true] {
        let mut target = if reverse {
            reverse_complement(sequence)
        } else {
            sequence.to_vec()
        };
        if circular {
            target.extend_from_within(..);
        }
        let short_index = Search::new(&target, 7);
        let long_index = Search::new(&target, 9);
        for model in &panel.genes {
            let index = if model.kind == "tRNA" {
                &short_index
            } else {
                &long_index
            };
            let (hits, search_limited) = index.hits(model.sequence.as_bytes(), &target);
            if search_limited {
                limited.insert(model.name.clone());
            }
            for hit in hits {
                if let Some(f) = feature(model, hit, &target, n, reverse) {
                    all.push(f);
                }
            }
        }
    }
    all.sort_by(|a, b| {
        b.alignment_score
            .cmp(&a.alignment_score)
            .then(b.length.cmp(&a.length))
            .then(a.five_prime_1based.cmp(&b.five_prime_1based))
            .then(a.name.cmp(&b.name))
    });
    let mut features: Vec<Feature> = Vec::new();
    for f in all {
        // Merge repeated seed-window observations of one physical locus, never
        // select a single genome-wide hit per gene or impose reference synteny.
        if features.iter().any(|old| {
            old.name == f.name
                && old.strand == f.strand
                && overlap(old, &f) * 5 >= old.length.min(f.length) * 4
        }) {
            continue;
        }
        features.push(f);
    }
    for i in 0..features.len() {
        for j in i + 1..features.len() {
            if features[i].name != features[j].name
                && features[i].kind == features[j].kind
                && overlap(&features[i], &features[j]) * 5
                    >= features[i].length.max(features[j].length) * 4
            {
                let ni = features[i].name.clone();
                let nj = features[j].name.clone();
                features[i].competing_models.push(nj);
                features[j].competing_models.push(ni);
            }
        }
    }
    for f in &mut features {
        if !f.competing_models.is_empty() {
            f.full_length_supported = false;
            f.status = "ambiguous_gene_identity".into();
        }
    }
    features.sort_by(|a, b| {
        a.five_prime_1based
            .cmp(&b.five_prime_1based)
            .then(a.name.cmp(&b.name))
            .then(a.strand.cmp(&b.strand))
    });
    let missing_models = panel
        .genes
        .iter()
        .filter(|g| !features.iter().any(|f| f.name == g.name))
        .map(|g| g.name.clone())
        .collect();
    let models_without_full_length_support = panel
        .genes
        .iter()
        .filter(|g| {
            !features
                .iter()
                .any(|f| f.name == g.name && f.full_length_supported)
        })
        .map(|g| g.name.clone())
        .collect();
    Prediction {
        features,
        missing_models,
        models_without_full_length_support,
        search_limited_models: limited.into_iter().collect(),
    }
}

pub fn rotate_to_pro(
    sequence: &[u8],
    prediction: &mut Prediction,
    circular: bool,
    requested: bool,
) -> (Vec<u8>, Rotation) {
    let mut rotation = Rotation {
        requested,
        applied: false,
        reason: "not_requested".into(),
        reverse_complemented: false,
        offset_0based: 0,
        input_anchor_1based: None,
    };
    let mut output = sequence.to_vec();
    if !requested {
        return (output, rotation);
    }
    let anchors: Vec<_> = prediction
        .features
        .iter()
        .filter(|f| f.name == "tRNA-Pro" && f.kind == "tRNA")
        .collect();
    let reason = if !circular {
        Some("linear_input")
    } else if prediction
        .search_limited_models
        .iter()
        .any(|x| x == "tRNA-Pro")
    {
        Some("trna_pro_search_incomplete")
    } else if anchors.is_empty() {
        Some("trna_pro_not_found")
    } else if anchors.len() != 1 {
        Some("multiple_trna_pro_loci")
    } else if !anchors[0].full_length_supported
        || anchors[0].five_prime_boundary != "aligned"
        || anchors[0].query_start != 1
    {
        Some("trna_pro_boundary_not_supported")
    } else {
        None
    };
    if let Some(reason) = reason {
        rotation.reason = reason.into();
        return (output, rotation);
    }
    let n = sequence.len();
    let anchor = anchors[0].five_prime_1based - 1;
    let reverse = anchors[0].strand == '-';
    let offset = if reverse { n - 1 - anchor } else { anchor };
    if reverse {
        output = reverse_complement(&output);
    }
    output.rotate_left(offset);
    rotation.applied = true;
    rotation.reason = "unique_full_length_trna_pro".into();
    rotation.reverse_complemented = reverse;
    rotation.offset_0based = offset;
    rotation.input_anchor_1based = Some(anchor + 1);
    for f in &mut prediction.features {
        let p = f.five_prime_1based - 1;
        let p = if reverse { n - 1 - p } else { p };
        f.five_prime_1based = (p + n - offset) % n + 1;
        if reverse {
            f.strand = if f.strand == '+' { '-' } else { '+' };
        }
        f.update_segments(n);
    }
    (output, rotation)
}

pub(super) struct Record {
    pub name: String,
    pub header: String,
    pub sequence: Vec<u8>,
}
fn read_records(path: &Path) -> Result<Vec<Record>> {
    let mut records: Vec<Record> = Vec::new();
    let mut names = BTreeSet::new();
    let mut total = 0;
    for line in crate::io::fastq::open(path)?.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            let name = header
                .split_whitespace()
                .next()
                .context("empty FASTA header")?
                .to_owned();
            ensure!(names.insert(name.clone()), "duplicate FASTA ID: {name}");
            records.push(Record {
                name,
                header: header.into(),
                sequence: Vec::new(),
            });
        } else if !line.trim().is_empty() {
            let record = records
                .last_mut()
                .context("sequence before first FASTA header")?;
            for b in line
                .bytes()
                .filter(|b| !b.is_ascii_whitespace())
                .map(|b| b.to_ascii_uppercase())
            {
                ensure!(
                    b"ACGTRYSWKMBDHVN".contains(&b),
                    "invalid DNA base in {}",
                    record.name
                );
                record.sequence.push(b);
                total += 1;
            }
            ensure!(
                record.sequence.len() <= 200_000 && total <= 5_000_000,
                "annotation supports mtDNA candidates up to 200 kb each / 5 Mb total"
            );
        }
    }
    ensure!(
        !records.is_empty() && records.iter().all(|r| !r.sequence.is_empty()),
        "empty FASTA or empty sequence"
    );
    Ok(records)
}

pub fn run(args: &AnnotationArgs) -> Result<bool> {
    let panel = if let Some(path) = &args.gene_models {
        Panel::parse(&fs::read_to_string(path)?)?
    } else {
        Panel::bundled(
            args.species
                .context("a species or gene-model panel is required")?,
        )?
    };
    let records = read_records(&args.fasta)?;
    let uncertainty = args
        .assembly_ambiguity
        .as_deref()
        .map(|path| uncertainty::read(path, &records))
        .transpose()?;
    let mut results = Vec::new();
    for (i, record) in records.iter().enumerate() {
        tracing::info!(contig = %record.name, length = record.sequence.len(), "predicting mitochondrial genes");
        let circular = match args.topology {
            Topology::Auto => record
                .header
                .split_whitespace()
                .skip(1)
                .find_map(|x| match x {
                    "annotation_topology=circular" => Some(true),
                    "annotation_topology=linear" => Some(false),
                    _ => None,
                })
                .unwrap_or_else(|| {
                    record
                        .header
                        .split_whitespace()
                        .skip(1)
                        .any(|x| x == "circular=true")
                }),
            Topology::Linear => false,
            Topology::Circular => true,
        };
        let mut prediction = predict(&record.sequence, &panel, circular);
        for (j, f) in prediction.features.iter_mut().enumerate() {
            f.id = format!("contig{}_gene{}", i + 1, j + 1);
        }
        let (sequence, rotation) = rotate_to_pro(
            &record.sequence,
            &mut prediction,
            circular,
            args.rotate_to_trna_pro,
        );
        let mut ncr = ncr::derive(
            &sequence,
            &prediction,
            circular,
            args.species,
            args.ncr_min_len as usize,
        );
        for region in &mut ncr.regions {
            region.id = format!("contig{}_{}", i + 1, region.id);
        }
        for locus in &mut ncr.loci {
            for id in &mut locus.region_ids {
                *id = format!("contig{}_{}", i + 1, id);
            }
        }
        let assembly_evidence = uncertainty
            .as_ref()
            .map(|input| uncertainty::project(input, record, i, &rotation, &prediction, &ncr));
        results.push(output::Annotated {
            sequence,
            prediction,
            rotation,
            circular,
            ncr,
            assembly_evidence,
        });
    }
    output::write(&args.out, &records, &results, &panel)?;
    let reports: Vec<_> = records
        .iter()
        .zip(&results)
        .map(|(r, a)| {
            json!({
                "contig": r.name, "input_header": r.header, "length": r.sequence.len(),
                "circular_input_assumed": a.circular, "rotation": a.rotation,
                "feature_counts": output::counts(&a.prediction.features, false),
                "full_length_counts": output::counts(&a.prediction.features, true),
                "prediction": a.prediction,
                "ncr": a.ncr,
                "assembly_evidence": a.assembly_evidence,
            })
        })
        .collect();
    crate::write_json(
        &args.out.join("report.json"),
        &json!({
        "schema": "mitovar-rust.annotation.v1", "status": "completed",
        "software_version": env!("CARGO_PKG_VERSION"),
            "method": "native_homology_assisted_gene_prediction", "model_accession": panel.accession,
            "genetic_code": 5, "input_fasta": args.fasta,
            "input_assembly_ambiguity": args.assembly_ambiguity,
            "assembly_status": "not_modified", "closure_validation": "not_performed",
            "rna_structure": "not_assessed", "records": reports,
        }),
    )?;
    Ok(true)
}
