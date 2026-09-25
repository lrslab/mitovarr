//! NCR labels describe conserved positions, not inferred sequence origins.
//! Boundaries come from observed sample gene placements after rotation.
use super::{Feature, Prediction, Segment};
use crate::cli::Species;
use serde::Serialize;
use std::collections::BTreeSet;

struct Spec {
    name: &'static str,
    a: &'static str,
    b: &'static str,
    expected: bool,
}
fn catalog(species: Option<Species>) -> Vec<Spec> {
    let Some(species) = species else {
        return Vec::new();
    };
    vec![
        Spec {
            name: "NCR1",
            a: "tRNA-Ala",
            b: "tRNA-Pro",
            expected: true,
        },
        Spec {
            name: "NCR2",
            a: if matches!(species, Species::Ctro) {
                "tRNA-Arg"
            } else {
                "tRNA-Gln"
            },
            b: if matches!(species, Species::Ctro) {
                "tRNA-Gln"
            } else {
                "tRNA-Phe"
            },
            expected: !matches!(species, Species::Cel),
        },
        Spec {
            name: "NCR3",
            a: "ND4",
            b: "COX1",
            expected: true,
        },
        Spec {
            name: "NCR4",
            a: "ND3",
            b: "ND5",
            expected: matches!(species, Species::Cbr),
        },
    ]
}
#[derive(Debug, Serialize)]
pub struct NcrRegion {
    pub id: String,
    pub name: String,
    pub locus: Option<String>,
    /// 1-based inclusive, in forward genomic order through a possible circular origin.
    pub segments: Vec<Segment>,
    pub length: usize,
    pub left_genes: Vec<String>,
    pub right_genes: Vec<String>,
    pub left_feature_ids: Vec<String>,
    pub right_feature_ids: Vec<String>,
    pub boundary_status: String,
    pub terminal_fragment: bool,
    pub ambiguous_bases: usize,
    /// A/T divided by unambiguous A/C/G/T; None when every base is ambiguous.
    pub at_fraction: Option<f64>,
    pub sequence_origin: &'static str,
    pub embedded_partial_homologs: Vec<String>,
    pub notes: Vec<String>,
    #[serde(skip)]
    pub sequence: Vec<u8>,
}
#[derive(Debug, Serialize)]
pub struct NcrLocus {
    pub name: String,
    pub expected_in_species_reference: bool,
    pub flank_genes: [String; 2],
    pub status: String,
    /// Gap span on the input sequence, not a read-validated SV length. None for open ends.
    pub observed_gap_length: Option<usize>,
    pub overlap_bases: usize,
    pub region_ids: Vec<String>,
    pub missing_flank_genes: Vec<String>,
    /// Actual gap or touching/overlapping flank bases, including loci without
    /// an emitted NCR region. Used only to associate existing assembly evidence.
    #[serde(skip)]
    pub(super) evidence_segments: Vec<Segment>,
}
#[derive(Debug, Serialize)]
pub struct NcrAnalysis {
    pub catalog: Option<String>,
    pub min_length: usize,
    pub regions: Vec<NcrRegion>,
    pub loci: Vec<NcrLocus>,
}
struct Gap {
    start: usize,
    length: usize,
    left: Vec<usize>,
    right: Vec<usize>,
    label: Option<usize>,
    terminal: bool,
}
fn parts(start: usize, length: usize, n: usize) -> Vec<Segment> {
    if start + length <= n {
        vec![Segment {
            start: start + 1,
            end: start + length,
        }]
    } else {
        vec![
            Segment {
                start: start + 1,
                end: n,
            },
            Segment {
                start: 1,
                end: start + length - n,
            },
        ]
    }
}
fn contains(f: &Feature, position: usize) -> bool {
    f.segments
        .iter()
        .any(|s| s.start <= position && position <= s.end)
}
fn names(indices: &[usize], p: &Prediction) -> Vec<String> {
    indices
        .iter()
        .map(|&i| p.features[i].name.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}
fn ids(indices: &[usize], p: &Prediction) -> Vec<String> {
    indices.iter().map(|&i| p.features[i].id.clone()).collect()
}
fn pair(spec: &Spec, left: &[usize], right: &[usize], p: &Prediction) -> bool {
    left.iter().any(|&l| {
        right.iter().any(|&r| {
            let (a, b) = (&p.features[l], &p.features[r]);
            let names_match =
                (a.name == spec.a && b.name == spec.b) || (a.name == spec.b && b.name == spec.a);
            // The catalog describes a -> b in transcriptional direction. The
            // opposite arc of a circle is not the same positional NCR. Mixed
            // strands retain a candidate, but cannot have supported boundaries.
            names_match && (a.strand != b.strand || (a.name == spec.a) == (a.strand == '+'))
        })
    })
}
fn supported(indices: &[usize], left: bool, p: &Prediction) -> bool {
    if indices.len() != 1 {
        return false;
    }
    let f = &p.features[indices[0]];
    let boundary = if left == (f.strand == '+') {
        &f.three_prime_boundary
    } else {
        &f.five_prime_boundary
    };
    f.query_coverage >= 0.8
        && f.competing_models.is_empty()
        && boundary != "partial_or_uncertain"
        && !p.search_limited_models.contains(&f.name)
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
fn touching(a: &Feature, b: &Feature, n: usize, circular: bool) -> bool {
    if a.strand != b.strand {
        return false;
    }
    let (Some(x), Some(y)) = (a.segments.last(), b.segments.first()) else {
        return false;
    };
    if a.strand == '+' {
        x.end + 1 == y.start || (circular && x.end == n && y.start == 1)
    } else {
        y.end + 1 == x.start || (circular && x.start == 1 && y.end == n)
    }
}

pub fn derive(
    sequence: &[u8],
    p: &Prediction,
    circular: bool,
    species: Option<Species>,
    min_length: usize,
) -> NcrAnalysis {
    let specs = catalog(species);
    let n = sequence.len();
    // A short secondary gene homolog (e.g. an ND5-derived repeat) must not consume
    // an NCR when a near-full copy of that gene is already represented elsewhere.
    // Keep all copies in the gene output, and keep partial sole copies in the mask.
    let near_full: BTreeSet<_> = p
        .features
        .iter()
        .filter(|f| f.query_coverage >= 0.8)
        .map(|f| f.name.as_str())
        .collect();
    let used: Vec<_> = p
        .features
        .iter()
        .enumerate()
        .filter(|(_, f)| f.query_coverage >= 0.8 || !near_full.contains(f.name.as_str()))
        .map(|(i, _)| i)
        .collect();
    let mut covered = vec![false; n];
    for &i in &used {
        for s in &p.features[i].segments {
            covered[s.start - 1..s.end].fill(true);
        }
    }
    let at = |position: usize| {
        used.iter()
            .copied()
            .filter(|&i| contains(&p.features[i], position))
            .collect::<Vec<_>>()
    };
    let mut runs = Vec::new();
    if covered.iter().any(|&b| b) {
        let mut i = 0;
        while i < n {
            if covered[i] {
                i += 1;
                continue;
            }
            let start = i;
            while i < n && !covered[i] {
                i += 1;
            }
            runs.push((start, i - start));
        }
        // Join circular head/tail BEFORE the size filter. Linear ends remain two
        // independent fragments and can never create a synthetic circular NCR.
        if circular
            && runs.len() >= 2
            && runs[0].0 == 0
            && runs.last().is_some_and(|&(s, l)| s + l == n)
        {
            let head = runs.remove(0);
            let tail = runs.pop().unwrap();
            runs.push((tail.0, tail.1 + head.1));
        }
    }
    let mut gaps = Vec::new();
    for (start, length) in runs {
        let end = (start + length) % n;
        let left = if start == 0 && !circular {
            Vec::new()
        } else {
            at((start + n - 1) % n + 1)
        };
        let right = if start + length == n && !circular {
            Vec::new()
        } else {
            at(end + 1)
        };
        let terminal = left.is_empty() || right.is_empty();
        let label = if terminal {
            let flank = if left.is_empty() { &right } else { &left };
            let hits: Vec<_> = specs
                .iter()
                .enumerate()
                .filter(|(_, s)| {
                    flank.iter().any(|&i| {
                        let f = &p.features[i];
                        let required = if left.is_empty() != (f.strand == '+') {
                            s.a
                        } else {
                            s.b
                        };
                        f.name == required
                    })
                })
                .map(|(i, _)| i)
                .collect();
            (hits.len() == 1).then(|| hits[0])
        } else {
            specs.iter().position(|s| pair(s, &left, &right, p))
        };
        gaps.push(Gap {
            start,
            length,
            left,
            right,
            label,
            terminal,
        });
    }
    let mut regions = Vec::new();
    let mut region_gap = Vec::new();
    for (g_index, g) in gaps
        .iter()
        .enumerate()
        .filter(|(_, g)| g.length >= min_length)
    {
        let mut notes = Vec::new();
        let same_strand = g.left.len() == 1
            && g.right.len() == 1
            && p.features[g.left[0]].strand == p.features[g.right[0]].strand;
        let boundary = if g.terminal {
            "open_contig_end"
        } else if supported(&g.left, true, p) && supported(&g.right, false, p) && same_strand {
            "supported"
        } else {
            "uncertain"
        };
        if !same_strand && !g.terminal {
            notes.push("flank_orientation_or_assignment_uncertain".into());
        }
        let is_canonical = g.label.is_some_and(|i| specs[i].expected) && boundary == "supported";
        let name = g
            .label
            .map(|i| {
                if is_canonical {
                    specs[i].name.to_owned()
                } else {
                    format!("{}_candidate", specs[i].name)
                }
            })
            .unwrap_or_else(|| "unclassified_gap".into());
        let segments = parts(g.start, g.length, n);
        let seq: Vec<_> = (0..g.length).map(|i| sequence[(g.start + i) % n]).collect();
        let known = seq.iter().filter(|&&b| b"ACGT".contains(&b)).count();
        let at_count = seq.iter().filter(|&&b| b == b'A' || b == b'T').count();
        let embedded = p
            .features
            .iter()
            .enumerate()
            .filter(|(i, f)| {
                !used.contains(i)
                    && f.segments.iter().any(|s| {
                        segments
                            .iter()
                            .any(|r| s.start <= r.end && r.start <= s.end)
                    })
            })
            .map(|(_, f)| f.id.clone())
            .collect();
        regions.push(NcrRegion {
            id: format!("ncr_region_{}", regions.len() + 1),
            name,
            locus: g.label.map(|i| specs[i].name.into()),
            segments,
            length: g.length,
            left_genes: names(&g.left, p),
            right_genes: names(&g.right, p),
            left_feature_ids: ids(&g.left, p),
            right_feature_ids: ids(&g.right, p),
            boundary_status: boundary.into(),
            terminal_fragment: g.terminal,
            ambiguous_bases: g.length - known,
            at_fraction: (known > 0).then(|| at_count as f64 / known as f64),
            sequence_origin: "not_assessed",
            embedded_partial_homologs: embedded,
            notes,
            sequence: seq,
        });
        region_gap.push(g_index);
    }
    let first = covered
        .iter()
        .position(|&x| x)
        .map(|x| at(x + 1))
        .unwrap_or_default();
    let last = covered
        .iter()
        .rposition(|&x| x)
        .map(|x| at(x + 1))
        .unwrap_or_default();
    let mut loci = Vec::new();
    for (s_index, spec) in specs.iter().enumerate() {
        let found: Vec<_> = regions
            .iter()
            .enumerate()
            .filter(|(i, _)| gaps[region_gap[*i]].label == Some(s_index))
            .map(|(_, r)| r)
            .collect();
        let bounded: Vec<_> = found.iter().filter(|r| !r.terminal_fragment).collect();
        let raw: Vec<_> = gaps
            .iter()
            .filter(|g| g.label == Some(s_index) && !g.terminal)
            .collect();
        let a: Vec<_> = used
            .iter()
            .copied()
            .filter(|&i| p.features[i].name == spec.a)
            .collect();
        let b: Vec<_> = used
            .iter()
            .copied()
            .filter(|&i| p.features[i].name == spec.b)
            .collect();
        let missing_flank_genes = [(spec.a, &a), (spec.b, &b)]
            .into_iter()
            .filter(|(_, v)| v.is_empty())
            .map(|(name, _)| name.into())
            .collect();
        let mut length = None;
        let mut overlap_bases = 0;
        let status = if bounded.len() > 1 {
            "multiple_candidates"
        } else if let Some(r) = bounded.first() {
            length = Some(r.length);
            if r.boundary_status != "supported" {
                "boundary_uncertain"
            } else if !spec.expected {
                "positional_candidate"
            } else if r.ambiguous_bases > 0 {
                "annotated_with_unknown_bases"
            } else {
                "annotated"
            }
        } else if !found.is_empty() || (!circular && pair(spec, &last, &first, p)) {
            "unresolved_contig_ends"
        } else if let [g] = raw.as_slice() {
            length = Some(g.length);
            if supported(&g.left, true, p)
                && supported(&g.right, false, p)
                && p.features[g.left[0]].strand == p.features[g.right[0]].strand
            {
                "below_min_length"
            } else {
                "boundary_uncertain"
            }
        } else if a.is_empty() || b.is_empty() {
            "missing_flank_genes"
        } else if a.len() != 1 || b.len() != 1 {
            "ambiguous_flank_copies"
        } else {
            let (fa, fb) = (&p.features[a[0]], &p.features[b[0]]);
            overlap_bases = overlap(fa, fb);
            if overlap_bases > 0 {
                length = Some(0);
                "overlapping_flanks"
            } else if touching(fa, fb, n, circular) {
                length = Some(0);
                "adjacent_without_gap"
            } else {
                "flanks_not_adjacent"
            }
        };
        let mut evidence_segments: Vec<_> = gaps
            .iter()
            .filter(|g| g.label == Some(s_index))
            .flat_map(|g| parts(g.start, g.length, n))
            .collect();
        if length == Some(0) && a.len() == 1 && b.len() == 1 {
            let (fa, fb) = (&p.features[a[0]], &p.features[b[0]]);
            if overlap_bases > 0 {
                for x in &fa.segments {
                    for y in &fb.segments {
                        let (start, end) = (x.start.max(y.start), x.end.min(y.end));
                        if start <= end {
                            evidence_segments.push(Segment { start, end });
                        }
                    }
                }
            } else if touching(fa, fb, n, circular) {
                // There are no intergenic bases at a zero-length locus. Both
                // sides of its literal junction still carry assembly evidence.
                let three = fa.segments.last().unwrap();
                let five = fb.segments.first().unwrap();
                for position in if fa.strand == '+' {
                    [three.end, five.start]
                } else {
                    [three.start, five.end]
                } {
                    evidence_segments.push(Segment {
                        start: position,
                        end: position,
                    });
                }
            }
        }
        loci.push(NcrLocus {
            name: spec.name.into(),
            expected_in_species_reference: spec.expected,
            flank_genes: [spec.a.into(), spec.b.into()],
            status: status.into(),
            observed_gap_length: length,
            overlap_bases,
            region_ids: found.iter().map(|r| r.id.clone()).collect(),
            missing_flank_genes,
            evidence_segments,
        });
    }
    NcrAnalysis {
        catalog: species.map(|s| s.name().into()),
        min_length,
        regions,
        loci,
    }
}
