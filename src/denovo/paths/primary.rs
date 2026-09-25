//! A usable primary graph path, with structural alternatives kept separately.
//! None of these choices grant the strict final-assembly status.
use super::*;
use crate::kmer::graph as graph_kmer;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

#[derive(Clone)]
pub(in crate::denovo) struct Region {
    pub start: usize,
    pub end: usize,
    pub evidence: String,
    pub reason: &'static str,
    pub alternatives: usize,
}

pub(in crate::denovo) struct Primary {
    pub candidate: Candidate,
    pub source: Value,
    pub regions: Vec<Region>,
}

/// Project an interval in a rotated, possibly reverse-complemented view back
/// onto the normalized FASTA. An empty allele marks its adjacent base; intervals
/// crossing the FASTA origin are split rather than spanning the entire genome.
pub(super) fn project_interval(
    start: usize,
    end: usize,
    length: usize,
    origin: usize,
    reverse: bool,
) -> Vec<(usize, usize)> {
    if length == 0 {
        return Vec::new();
    }
    let width = (end - start).max(1).min(length);
    let from = if reverse {
        (length - (origin + start + width) % length) % length
    } else {
        (origin + start) % length
    };
    if from + width <= length {
        vec![(from, from + width)]
    } else {
        vec![(0, from + width - length), (from, length)]
    }
}

fn normalize_regions(
    sequence: Vec<u8>,
    circular: bool,
    regions: Vec<Region>,
) -> (Vec<u8>, Vec<Region>) {
    let reverse = kmer::reverse_complement(&sequence);
    let offset = if circular {
        super::super::rotation_start(&sequence)
    } else {
        0
    };
    let rev_offset = if circular {
        super::super::rotation_start(&reverse)
    } else {
        0
    };
    let forward = [&sequence[offset..], &sequence[..offset]].concat();
    let reversed = [&reverse[rev_offset..], &reverse[..rev_offset]].concat();
    let is_reverse = reversed < forward;
    let normalized = if is_reverse { reversed } else { forward };
    let origin = if is_reverse {
        rev_offset
    } else {
        sequence.len() - offset
    };
    let mut mapped: Vec<_> = regions
        .into_iter()
        .flat_map(|r| {
            project_interval(r.start, r.end, sequence.len(), origin, is_reverse)
                .into_iter()
                .map(move |(start, end)| Region {
                    start,
                    end,
                    ..r.clone()
                })
        })
        .collect();
    mapped.sort_by_key(|r| (r.start, r.end, r.evidence.clone()));
    (normalized, mapped)
}

fn choose(locus: &Locus, settings: &GraphSettings) -> Option<usize> {
    if locus.selected.is_some() {
        return locus.selected;
    }
    if locus.exhausted || locus.length_limited || locus.routes.is_empty() {
        return None;
    }
    let best = (0..locus.routes.len())
        .max_by(|&a, &b| locus.totals[a].total_cmp(&locus.totals[b]).then(b.cmp(&a)))?;
    // Complete cycle families have their own shared-pool whole-path comparison.
    if locus.routes[best].to == locus.from || locus.spanning[best] < settings.min_junction_fragments
    {
        return None;
    }
    (0..locus.routes.len())
        .filter(|&i| i != best)
        .all(|i| {
            locus.totals[best] > locus.totals[i]
                && locus.strong_votes[best][i] >= settings.min_junction_fragments
                && locus.strong_votes[best][i] > locus.strong_votes[i][best]
        })
        .then_some(best)
}

fn route_key(graph: &PathGraph, locus: &Locus, route: &Route, overlap: usize) -> Vec<u8> {
    let mut sequence = graph.segments[locus.from].sequence.clone();
    for &id in route.nodes.iter().chain(std::iter::once(&route.to)) {
        sequence.extend_from_slice(&graph.segments[id].sequence[overlap..]);
    }
    normalize(&sequence, false)
}

/// A physical join has one spelling in both graph orientations. A qualified
/// strict choice takes precedence; conflicting relaxed views remain open.
fn reciprocal_choices(
    graph: &PathGraph,
    loci: &[Locus],
    overlap: usize,
    choices: &mut [Option<usize>],
) {
    let mut views = BTreeMap::<Vec<u8>, Vec<(usize, usize)>>::new();
    for (li, l) in loci.iter().enumerate() {
        for (ri, r) in l.routes.iter().enumerate().filter(|(_, r)| r.to != l.from) {
            views
                .entry(route_key(graph, l, r, overlap))
                .or_default()
                .push((li, ri));
        }
    }
    let keys: Vec<_> = loci
        .iter()
        .zip(choices.iter())
        .map(|(l, r)| r.map(|i| route_key(graph, l, &l.routes[i], overlap)))
        .collect();
    let original = choices.to_vec();
    for (li, key) in keys.iter().enumerate() {
        let Some(key) = key else {
            continue;
        };
        let Some(related) = views.get(key) else {
            continue;
        };
        if related.iter().any(|&(other, _)| {
            keys[other].as_ref().is_some_and(|k| k != key)
                && !(loci[li].selected == original[li] && loci[other].selected != original[other])
        }) {
            choices[li] = None;
        }
    }
    let retained = choices.to_vec();
    for (li, selected) in retained.iter().enumerate() {
        if selected.is_none() {
            continue;
        }
        if let Some(related) = views.get(keys[li].as_ref().unwrap()) {
            for &(other, ri) in related {
                if !loci[other].exhausted && !loci[other].length_limited {
                    choices[other] = Some(ri);
                }
            }
        }
    }
}

pub(super) fn local_candidates(
    graph: &PathGraph,
    loci: &[Locus],
    anchors: &[usize],
    overlap: usize,
    settings: &GraphSettings,
) -> (Vec<Primary>, Vec<Option<usize>>) {
    let mut choices: Vec<_> = loci.iter().map(|l| choose(l, settings)).collect();
    reciprocal_choices(graph, loci, overlap, &mut choices);
    let mut connections = HashMap::new();
    let mut incoming = HashMap::<usize, usize>::new();
    for (li, selected) in choices.iter().enumerate() {
        if let Some(ri) = *selected {
            connections.insert(loci[li].from, (li, ri));
            *incoming.entry(loci[li].routes[ri].to).or_default() += 1;
        }
    }
    connections.retain(|_, entry| incoming[&loci[entry.0].routes[entry.1].to] == 1);
    let starts: Vec<_> = anchors
        .iter()
        .copied()
        .filter(|id| incoming.get(id) != Some(&1))
        .chain(anchors.iter().copied())
        .collect();
    let mut consumed = HashSet::new();
    let mut assembled = BTreeMap::new();
    for first in starts {
        if consumed.contains(&first) {
            continue;
        }
        let mut sequence = graph.segments[first].sequence.clone();
        let mut node = first;
        let mut circular = false;
        let mut regions = Vec::new();
        let mut joins = Vec::new();
        loop {
            consumed.insert(node);
            let Some(&(li, ri)) = connections.get(&node) else {
                break;
            };
            let locus = &loci[li];
            let route = &locus.routes[ri];
            if consumed.contains(&route.to) && route.to != first {
                break;
            }
            if locus.routes.len() > 1 {
                let prefix = locus
                    .routes
                    .iter()
                    .map(|r| {
                        r.sequence
                            .iter()
                            .zip(&route.sequence)
                            .take_while(|(a, b)| a == b)
                            .count()
                    })
                    .min()
                    .unwrap();
                let suffix = locus
                    .routes
                    .iter()
                    .map(|r| {
                        r.sequence[prefix..]
                            .iter()
                            .rev()
                            .zip(route.sequence[prefix..].iter().rev())
                            .take_while(|(a, b)| a == b)
                            .count()
                    })
                    .min()
                    .unwrap();
                if prefix + suffix < route.sequence.len()
                    || locus.routes.iter().any(|r| r.sequence != route.sequence)
                {
                    let offset = sequence.len() - route.left_end - overlap;
                    regions.push(Region {
                        start: offset + prefix,
                        end: offset + route.sequence.len() - suffix,
                        evidence: format!("locus_{}", li + 1),
                        reason: if locus.selected == Some(ri) {
                            "minor_path_retained"
                        } else {
                            "majority_path_structure_unresolved"
                        },
                        alternatives: locus.routes.len() - 1,
                    });
                }
            }
            joins
                .push(json!({"locus":li+1,"path":ri+1,"strict_selected":locus.selected==Some(ri)}));
            for &id in &route.nodes {
                sequence.extend_from_slice(&graph.segments[id].sequence[overlap..]);
                consumed.insert(id);
            }
            if route.to == first {
                circular = true;
                sequence.truncate(sequence.len() - overlap);
                break;
            }
            sequence.extend_from_slice(&graph.segments[route.to].sequence[overlap..]);
            node = route.to;
        }
        if !circular {
            add_open_ends(&mut regions, sequence.len());
        }
        let (sequence, regions) = normalize_regions(sequence, circular, regions);
        assembled
            .entry((sequence.clone(), circular))
            .or_insert(Primary {
                candidate: Candidate { sequence, circular },
                source: json!({"type":"joined_graph_paths","joins":joins}),
                regions,
            });
    }
    for (id, segment) in graph.segments.iter().enumerate() {
        if consumed.contains(&id) {
            continue;
        }
        let sequence = normalize(&segment.sequence, segment.circular);
        let mut regions = Vec::new();
        if !segment.circular {
            add_open_ends(&mut regions, sequence.len());
        }
        assembled
            .entry((sequence.clone(), segment.circular))
            .or_insert(Primary {
                candidate: Candidate {
                    sequence,
                    circular: segment.circular,
                },
                source: json!({"type":"graph_segment","segment":id+1}),
                regions,
            });
    }
    (assembled.into_values().collect(), choices)
}

fn add_open_ends(regions: &mut Vec<Region>, length: usize) {
    if length == 0 {
        return;
    }
    for (start, end) in [(0, 1), (length - 1, length)] {
        regions.push(Region {
            start,
            end,
            evidence: "graph_end".into(),
            reason: "open_graph_end",
            alternatives: 0,
        });
    }
}

fn support_mass(candidate: &Candidate, counts: &HashMap<graph_kmer::Code, u32>, k: usize) -> u64 {
    let mut sequence = candidate.sequence.clone();
    if candidate.circular && !sequence.is_empty() {
        sequence.extend(candidate.sequence.iter().copied().cycle().take(k - 1));
    }
    let mut codes = HashSet::new();
    graph_kmer::visit(&sequence, None, k, 0, |_, code| {
        codes.insert(graph_kmer::canonical(code, k));
    });
    codes
        .iter()
        .map(|code| u64::from(counts.get(code).copied().unwrap_or(0)))
        .sum()
}

pub(in crate::denovo) fn write_primary(
    mut candidates: Vec<Primary>,
    counts: &HashMap<graph_kmer::Code, u32>,
    k: usize,
    out: &Path,
    final_candidate: Option<&Candidate>,
) -> Result<Value> {
    // Whole-cycle likelihood has already chosen ONE representative per family.
    // Mass chooses the main component, never the repeat count within a family.
    // Collapse duplicate spellings, preferring the more detailed cycle record.
    let mut unique = BTreeMap::new();
    for p in candidates.drain(..) {
        unique.insert((p.candidate.sequence.clone(), p.candidate.circular), p);
    }
    if let Some(c) = final_candidate {
        unique
            .entry((c.sequence.clone(), c.circular))
            .or_insert_with(|| Primary {
                candidate: c.clone(),
                source: json!({"type":"validated_candidate"}),
                regions: Vec::new(),
            });
    }
    let selected = unique
        .into_values()
        .filter(|p| p.candidate.sequence.len() >= k)
        .max_by_key(|p| {
            (
                final_candidate.is_some_and(|c| {
                    c.sequence == p.candidate.sequence && c.circular == p.candidate.circular
                }),
                support_mass(&p.candidate, counts, k),
                p.candidate.circular,
                std::cmp::Reverse(p.candidate.sequence.clone()),
            )
        });
    let mut ambiguity = BufWriter::new(File::create(out.join("primary_ambiguity.tsv"))?);
    writeln!(
        ambiguity,
        "sequence\tstart\tend\tevidence\treason\talternatives"
    )?;
    let Some(mut primary) = selected else {
        ambiguity.flush()?;
        return Ok(Value::Null);
    };
    primary
        .regions
        .sort_by_key(|r| (r.start, r.end, r.evidence.clone()));
    let c = &primary.candidate;
    let verified =
        final_candidate.is_some_and(|f| f.sequence == c.sequence && f.circular == c.circular);
    crate::io::fasta::write_fasta_records(
        &out.join("primary.fasta"),
        &[(
            format!(
                "primary circular={} length={} final_verified={verified}",
                c.circular,
                c.sequence.len()
            ),
            c.sequence.clone(),
        )],
    )?;
    for r in &primary.regions {
        writeln!(
            ambiguity,
            "primary\t{}\t{}\t{}\t{}\t{}",
            r.start + 1,
            r.end,
            r.evidence,
            r.reason,
            r.alternatives
        )?;
    }
    ambiguity.flush()?;
    Ok(json!({
        "file":"primary.fasta","ambiguity_file":"primary_ambiguity.tsv","length":c.sequence.len(),
        "circular":c.circular,"final_verified":verified,"source":primary.source,
        "selection":"whole-cycle fragment score within a family; distinct sample k-mer support mass between components; validated final candidate takes precedence",
        "canonical_kmer_support_mass":support_mass(c, counts, k),
        "uncertainty_intervals":primary.regions.iter().map(|r|json!({"start":r.start+1,"end":r.end,"evidence":r.evidence,"reason":r.reason,"alternatives":r.alternatives})).collect::<Vec<_>>(),
        "interval_coordinates":"one-based inclusive in primary.fasta, before annotation or tRNA-Pro rotation; empty alleles mark one adjacent base",
        "structural_uniqueness_claim":false,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uncertain_intervals_follow_rotation_reverse_complement_and_empty_alleles() {
        assert_eq!(
            project_interval(0, 4, 20, 18, false),
            vec![(0, 2), (18, 20)]
        );
        assert_eq!(project_interval(2, 6, 20, 3, true), vec![(11, 15)]);
        assert_eq!(project_interval(2, 2, 20, 3, true), vec![(14, 15)]);
        let sequence = b"TGTTACGGAACCTGT".to_vec();
        let region = Region {
            start: 3,
            end: 7,
            evidence: "test".into(),
            reason: "test",
            alternatives: 1,
        };
        for circular in [false, true] {
            let (normalized, intervals) =
                normalize_regions(sequence.clone(), circular, vec![region.clone()]);
            assert_eq!(normalized, normalize(&sequence, circular));
            let mut positions: Vec<_> = intervals.iter().flat_map(|r| r.start..r.end).collect();
            positions.sort_unstable();
            assert_eq!(positions.len(), 4);
            let bases: Vec<_> = positions.iter().map(|&p| normalized[p]).collect();
            let mut actual = bases.clone();
            actual.sort_unstable();
            let mut forward = sequence[3..7].to_vec();
            forward.sort_unstable();
            let mut reverse = kmer::reverse_complement(&sequence[3..7]);
            reverse.sort_unstable();
            assert!(actual == forward || actual == reverse);
        }
    }

    #[test]
    fn repeated_traversal_does_not_multiply_component_support() {
        let sequence = b"ACGTTAGGCATCTGAACCTGTAGTCCGATCAATGCG".to_vec();
        let k = 15;
        let extended = [sequence.clone(), sequence.clone()].concat();
        let mut counts = HashMap::new();
        graph_kmer::visit(&extended, None, k, 0, |_, code| {
            counts.insert(graph_kmer::canonical(code, k), 20);
        });
        let one = Candidate {
            sequence,
            circular: true,
        };
        let two = Candidate {
            sequence: extended,
            circular: true,
        };
        assert_eq!(
            support_mass(&one, &counts, k),
            support_mass(&two, &counts, k)
        );
    }

    #[test]
    fn validated_circle_remains_primary_when_relaxed_connections_change_candidates() {
        let out =
            std::env::temp_dir().join(format!("mitovar-primary-final-{}", std::process::id()));
        std::fs::create_dir_all(&out).unwrap();
        let final_candidate = Candidate {
            sequence: b"ACGTTAGGCATCTGAACCTGTAGTCCGATCAATGCG".to_vec(),
            circular: true,
        };
        let alternatives = vec![Primary {
            candidate: Candidate {
                sequence: [
                    final_candidate.sequence.clone(),
                    final_candidate.sequence.clone(),
                ]
                .concat(),
                circular: false,
            },
            source: json!({"type":"joined_graph_paths"}),
            regions: Vec::new(),
        }];
        let report = write_primary(
            alternatives,
            &HashMap::new(),
            15,
            &out,
            Some(&final_candidate),
        )
        .unwrap();
        assert_eq!(report["final_verified"], true);
        assert_eq!(report["circular"], true);
        assert_eq!(
            crate::io::fasta::load_reference(&out.join("primary.fasta"), None)
                .unwrap()
                .sequence,
            final_candidate.sequence
        );
        std::fs::remove_dir_all(&out).unwrap();
    }
}
