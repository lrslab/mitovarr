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

pub(in crate::denovo) fn normalize_regions(
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
    let (spanning, votes) = if locus.routes[best].primary_interval.is_some() {
        (&locus.primary_spanning, &locus.primary_strong_votes)
    } else {
        (&locus.spanning, &locus.strong_votes)
    };
    // Complete cycle families have their own shared-pool whole-path comparison.
    if locus.routes[best].to == locus.from || spanning[best] < settings.min_junction_fragments {
        return None;
    }
    (0..locus.routes.len())
        .filter(|&i| i != best)
        .all(|i| {
            locus.totals[best] > locus.totals[i]
                && votes[best][i] >= settings.min_junction_fragments
                && votes[best][i] > votes[i][best]
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

/// A physical join has one spelling in both graph orientations. A comparison
/// with alternatives outranks a single-route observation. Among competitive
/// comparisons a strict choice still takes precedence; equal conflicts stay
/// open. These priorities affect partial primary output only.
fn reciprocal_choices(
    graph: &PathGraph,
    loci: &[Locus],
    overlap: usize,
    choices: &mut [Option<usize>],
) -> HashSet<usize> {
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
    let mut blocked = HashSet::new();
    let priority = |li: usize| (loci[li].routes.len() > 1, loci[li].selected == original[li]);
    for (li, key) in keys.iter().enumerate() {
        let Some(key) = key else {
            continue;
        };
        let Some(related) = views.get(key) else {
            continue;
        };
        if related.iter().any(|&(other, _)| {
            keys[other].as_ref().is_some_and(|k| k != key) && priority(li) <= priority(other)
        }) {
            choices[li] = None;
            blocked.insert(li);
        }
    }
    let retained = choices.to_vec();
    let mut projected = vec![BTreeSet::new(); loci.len()];
    for (li, selected) in retained.iter().enumerate() {
        if selected.is_some() {
            if let Some(related) = views.get(keys[li].as_ref().unwrap()) {
                for &(other, ri) in related {
                    if !loci[other].exhausted && !loci[other].length_limited {
                        projected[other].insert(ri);
                    }
                }
            }
        }
    }
    // Two reverse-only observations can project different routes into an
    // undecided forward locus. Keep that conflict open in both orientations;
    // iteration order must never select whichever projection is visited last.
    let conflicted: BTreeSet<_> = retained
        .iter()
        .enumerate()
        .filter(|(_, selected)| selected.is_some())
        .filter_map(|(li, _)| {
            let key = keys[li].as_ref()?;
            views
                .get(key)?
                .iter()
                .any(|&(other, _)| projected[other].len() > 1)
                .then_some(key.clone())
        })
        .collect();
    for key in &conflicted {
        blocked.extend(views[key].iter().map(|&(li, _)| li));
    }
    for (li, selected) in retained.iter().enumerate() {
        if selected.is_none() {
            continue;
        }
        if conflicted.contains(keys[li].as_ref().unwrap()) {
            choices[li] = None;
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
    // A stronger reciprocal observation can resolve a rejected choice. Only
    // unresolved conflicts must remain distinct from a lack of read evidence.
    blocked.retain(|&li| choices[li].is_none());
    blocked
}

/// One strand view of a two-copy direct repeat: the entrance join into the loop
/// and the loop's return join to the exit, plus the other pairing's routes,
/// each as (locus, route) and sorted by locus.
struct Pairing {
    repeat: usize,
    segments: [usize; 4],
    selected: [(usize, usize); 2],
    other: [(usize, usize); 2],
}

/// Two loci A and Y whose two routes each cross the same non-anchor repeat R to
/// the same two exits. When one exit leads back to Y through joins already
/// chosen, R has two direct copies on one molecule: A → R → loop … Y → R → exit.
/// The other pairing would separate the loop and one copy of R as a subcircle.
/// Two loops through R likewise close one circle. Reads may already have chosen
/// either join; topology fills the other when nothing conflicts and exactly one
/// such pairing exists in a strand view. The joins then extend to the other
/// strand view by spelling, as reciprocal choices do. Strict joins are unchanged.
fn repeat_pairings(
    graph: &PathGraph,
    loci: &[Locus],
    overlap: usize,
    choices: &mut [Option<usize>],
    blocked: &HashSet<usize>,
    settings: &GraphSettings,
) -> (HashSet<(usize, usize)>, Vec<Value>) {
    let by_source: HashMap<_, _> = loci.iter().enumerate().map(|(i, l)| (l.from, i)).collect();
    let mut views = BTreeMap::<Vec<u8>, Vec<(usize, usize)>>::new();
    for (li, l) in loci.iter().enumerate() {
        for (ri, r) in l.routes.iter().enumerate() {
            views
                .entry(route_key(graph, l, r, overlap))
                .or_default()
                .push((li, ri));
        }
    }
    let mut crossing = BTreeMap::<usize, BTreeSet<usize>>::new();
    for (li, l) in loci.iter().enumerate() {
        for &node in l.routes.iter().flat_map(|r| &r.nodes) {
            crossing.entry(node).or_default().insert(li);
        }
    }
    let scored = |li: usize| !loci[li].exhausted && !loci[li].length_limited;
    // Both routes cross the repeat, to two different exits.
    let exits = |li: usize, repeat: usize| match &loci[li].routes[..] {
        [a, b] if a.to != b.to && a.nodes.contains(&repeat) && b.nodes.contains(&repeat) => {
            Some([a.to.min(b.to), a.to.max(b.to)])
        }
        _ => None,
    };
    let route_to = |li: usize, to: usize| {
        loci[li]
            .routes
            .iter()
            .position(|r| r.to == to)
            .expect("exit route")
    };
    // Keyed by the joins: the two loops of one circle find one pairing twice.
    let mut groups = BTreeMap::<Vec<u8>, Vec<BTreeMap<[(usize, usize); 2], Pairing>>>::new();
    for (&repeat, crossing) in &crossing {
        let entrances: Vec<_> = crossing
            .iter()
            .copied()
            .filter(|&li| exits(li, repeat).is_some())
            .collect();
        let reaches = |mut node: usize, target: usize| {
            for _ in 0..=loci.len() {
                if node == target {
                    return true;
                }
                let Some(route) = by_source
                    .get(&node)
                    .and_then(|&li| choices[li].map(|ri| &loci[li].routes[ri]))
                else {
                    return false;
                };
                if route.nodes.contains(&repeat) {
                    return false;
                }
                node = route.to;
            }
            false
        };
        let compatible = |li: usize, ri: usize| choices[li].map_or(scored(li), |c| c == ri);
        let mut found = BTreeMap::new();
        for &a in &entrances {
            for &y in &entrances {
                let Some(targets) =
                    exits(a, repeat).filter(|t| a != y && exits(y, repeat) == Some(*t))
                else {
                    continue;
                };
                for (inner, exit) in [(targets[0], targets[1]), (targets[1], targets[0])] {
                    let (ra, ry) = (route_to(a, inner), route_to(y, exit));
                    let claimed = loci.iter().enumerate().any(|(li, l)| {
                        li != a
                            && li != y
                            && choices[li].is_some_and(|c| [inner, exit].contains(&l.routes[c].to))
                    });
                    if compatible(a, ra)
                        && compatible(y, ry)
                        && !claimed
                        && reaches(inner, loci[y].from)
                    {
                        let mut selected = [(a, ra), (y, ry)];
                        let mut other = [(a, route_to(a, exit)), (y, route_to(y, inner))];
                        selected.sort_unstable();
                        other.sort_unstable();
                        found.entry(selected).or_insert(Pairing {
                            repeat,
                            segments: [loci[a].from, inner, loci[y].from, exit],
                            selected,
                            other,
                        });
                    }
                }
            }
        }
        if !found.is_empty() {
            groups
                .entry(normalize(&graph.segments[repeat].sequence, false))
                .or_default()
                .push(found);
        }
    }
    let key =
        |&(li, ri): &(usize, usize)| route_key(graph, &loci[li], &loci[li].routes[ri], overlap);
    let mut inferred = HashSet::new();
    let mut reports = Vec::new();
    for found in groups.into_values() {
        let ambiguous = found.iter().any(|f| f.len() > 1);
        let candidates: Vec<_> = found.iter().filter_map(|f| f.values().next()).collect();
        let first = candidates[0];
        let keys: BTreeSet<_> = first.selected.iter().map(key).collect();
        let agree = candidates
            .iter()
            .all(|p| p.selected.iter().map(key).collect::<BTreeSet<_>>() == keys);
        // Every strand view of both joins, as reciprocal choices would project them.
        let joins: Vec<_> = keys.iter().flat_map(|k| views[k].iter().copied()).collect();
        let conflict = joins
            .iter()
            .any(|&(li, ri)| blocked.contains(&li) || choices[li].map_or(!scored(li), |c| c != ri));
        let status = if ambiguous {
            "ambiguous"
        } else if !agree || conflict {
            "conflicting_choice"
        } else {
            "applied"
        };
        let spanning = |&(li, ri): &(usize, usize)| {
            let l = &loci[li];
            json!({"locus":li+1,"path":ri+1,"source_segment":l.from+1,"target_segment":l.routes[ri].to+1,
                "spanning_fragments":scored(li).then_some(l.spanning[ri]),
                "read_choice":choices[li]==Some(ri) && !inferred.contains(&(li, ri))})
        };
        // The same physical route can recruit different fragments from its
        // two source anchors. Check either view without adding duplicate counts.
        let observed = first
            .other
            .iter()
            .flat_map(|route| &views[&key(route)])
            .any(|&(li, ri)| {
                scored(li) && loci[li].spanning[ri] >= settings.min_junction_fragments
            });
        reports.push(json!({"repeat_segment":first.repeat+1,
            "repeat_length":graph.segments[first.repeat].sequence.len(),
            "entrance_segment":first.segments[0]+1,"loop_start_segment":first.segments[1]+1,
            "loop_end_segment":first.segments[2]+1,"exit_segment":first.segments[3]+1,
            "selected_routes":first.selected.iter().map(spanning).collect::<Vec<_>>(),
            "alternative_routes":first.other.iter().map(spanning).collect::<Vec<_>>(),
            "alternative_pairing_reads":if observed {"observed"} else {"not_observed"},
            "status":status}));
        if status == "applied" {
            for (li, ri) in joins {
                if choices[li].is_none() {
                    choices[li] = Some(ri);
                    inferred.insert((li, ri));
                }
            }
        }
    }
    (inferred, reports)
}

pub(super) fn local_candidates(
    graph: &PathGraph,
    loci: &[Locus],
    anchors: &[usize],
    overlap: usize,
    settings: &GraphSettings,
) -> (Vec<Primary>, Vec<Option<usize>>, Vec<Value>) {
    let mut choices: Vec<_> = loci.iter().map(|l| choose(l, settings)).collect();
    let blocked = reciprocal_choices(graph, loci, overlap, &mut choices);
    let (inferred, pairings) =
        repeat_pairings(graph, loci, overlap, &mut choices, &blocked, settings);
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
            let topology = inferred.contains(&(li, ri));
            if topology {
                // The whole repeat copy: which entrance pairs with which exit.
                let offset = sequence.len() - route.left_end - overlap;
                regions.push(Region {
                    start: offset + route.left_end,
                    end: offset + route.right_start,
                    evidence: format!("locus_{}", li + 1),
                    reason: "topology_inferred_repeat_pairing",
                    alternatives: locus.routes.len() - 1,
                });
            } else if locus.routes.len() > 1 {
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
                let (prefix, end) = route
                    .primary_interval
                    .unwrap_or((prefix, route.sequence.len() - suffix));
                if prefix < end || locus.routes.iter().any(|r| r.sequence != route.sequence) {
                    let offset = sequence.len() - route.left_end - overlap;
                    regions.push(Region {
                        start: offset + prefix,
                        end: offset + end,
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
            joins.push(
                json!({"locus":li+1,"path":ri+1,"strict_selected":locus.selected==Some(ri),
                "basis":if topology {"single_molecule_repeat_pairing"} else {"fragment_score"}}),
            );
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
    (assembled.into_values().collect(), choices, pairings)
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

fn support_mass(candidate: &Candidate, counts: &graph_kmer::Map<u32>, k: usize) -> u64 {
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

pub(in crate::denovo) fn write_outputs(
    candidates: Vec<Primary>,
    counts: &graph_kmer::Map<u32>,
    k: usize,
    out: &Path,
    final_candidate: Option<&Candidate>,
) -> Result<(Value, Value)> {
    // Whole-cycle likelihood has already chosen ONE representative per family.
    // Mass chooses the main component, never the repeat count within a family.
    // Collapse duplicate spellings, preferring the more detailed cycle record.
    let mut unique = BTreeMap::new();
    for p in candidates {
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
    let ranked: Vec<_> = unique
        .into_values()
        .filter(|p| p.candidate.sequence.len() >= k)
        .map(|mut p| {
            p.regions
                .sort_by_key(|r| (r.start, r.end, r.evidence.clone()));
            let mass = support_mass(&p.candidate, counts, k);
            let verified = final_candidate.is_some_and(|c| {
                c.sequence == p.candidate.sequence && c.circular == p.candidate.circular
            });
            (p, mass, verified)
        })
        .collect();
    let rank = |(p, mass, verified): &(Primary, u64, bool)| {
        (
            *verified,
            *mass,
            p.candidate.circular,
            std::cmp::Reverse(p.candidate.sequence.clone()),
        )
    };
    let selected = ranked.iter().max_by_key(|p| rank(p));
    let primary = if let Some((p, mass, verified)) = selected {
        write_selected(p, *mass, *verified, out, "primary")?
    } else {
        write_ambiguity(out, "primary", &[])?;
        Value::Null
    };
    // Select independently: a dominant linear primary must not hide a usable
    // circular candidate, and a small residual circle must not replace it.
    let circles: Vec<_> = ranked
        .iter()
        .filter(|(p, _, _)| p.candidate.circular)
        .collect();
    let selection = "whole-cycle fragment score within a family; distinct sample k-mer support mass between eligible circular components; validated final candidate takes precedence; exact support ties use normalized sequence lexicographic order";
    let best_circular =
        if let Some((p, mass, verified)) = circles.iter().copied().max_by_key(|p| rank(p)) {
            let mut report = write_selected(p, *mass, *verified, out, "best_circular")?;
            report["status"] = json!("selected");
            report["selection"] = json!(selection);
            report["eligible_representative_count"] = json!(circles.len());
            report
        } else {
            json!({
                "status":"unavailable","reason":"no_eligible_circular_candidate",
                "file":null,"ambiguity_file":null,"length":null,"circular":null,
                "final_verified":false,"source":null,"selection":selection,
                "eligible_representative_count":0,"canonical_kmer_support_mass":null,
                "uncertainty_intervals":[],"structural_uniqueness_claim":false,
            })
        };
    Ok((primary, best_circular))
}

fn write_ambiguity(out: &Path, id: &str, regions: &[Region]) -> Result<()> {
    let mut ambiguity = BufWriter::new(File::create(out.join(format!("{id}_ambiguity.tsv")))?);
    writeln!(
        ambiguity,
        "sequence\tstart\tend\tevidence\treason\talternatives"
    )?;
    for r in regions {
        writeln!(
            ambiguity,
            "{id}\t{}\t{}\t{}\t{}\t{}",
            r.start + 1,
            r.end,
            r.evidence,
            r.reason,
            r.alternatives
        )?;
    }
    ambiguity.flush()?;
    Ok(())
}

fn write_selected(
    primary: &Primary,
    mass: u64,
    verified: bool,
    out: &Path,
    id: &str,
) -> Result<Value> {
    let c = &primary.candidate;
    let file = format!("{id}.fasta");
    crate::io::fasta::write_fasta_records(
        &out.join(&file),
        &[(
            format!(
                "{id} circular={} length={} final_verified={verified}",
                c.circular,
                c.sequence.len()
            ),
            c.sequence.clone(),
        )],
    )?;
    write_ambiguity(out, id, &primary.regions)?;
    Ok(json!({
        "file":file,"ambiguity_file":format!("{id}_ambiguity.tsv"),"length":c.sequence.len(),
        "circular":c.circular,"final_verified":verified,"source":primary.source,
        "selection":"whole-cycle fragment score within a family; distinct sample k-mer support mass between components; validated final candidate takes precedence",
        "canonical_kmer_support_mass":mass,
        "uncertainty_intervals":primary.regions.iter().map(|r|json!({"start":r.start+1,"end":r.end,"evidence":r.evidence,"reason":r.reason,"alternatives":r.alternatives})).collect::<Vec<_>>(),
        "interval_coordinates":format!("one-based inclusive in {file}, before annotation or tRNA-Pro rotation; empty alleles mark one adjacent base"),
        "structural_uniqueness_claim":false,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn output_dir(name: &str) -> std::path::PathBuf {
        let out =
            std::env::temp_dir().join(format!("mitovar-primary-{name}-{}", std::process::id()));
        std::fs::create_dir_all(&out).unwrap();
        out
    }

    fn representative(sequence: &[u8], circular: bool, source: Value) -> Primary {
        Primary {
            candidate: Candidate {
                sequence: normalize(sequence, circular),
                circular,
            },
            source,
            regions: Vec::new(),
        }
    }

    fn add_counts(counts: &mut graph_kmer::Map<u32>, candidate: &Candidate, count: u32, k: usize) {
        let mut sequence = candidate.sequence.clone();
        if candidate.circular {
            sequence.extend(candidate.sequence.iter().copied().cycle().take(k - 1));
        }
        graph_kmer::visit(&sequence, None, k, 0, |_, code| {
            counts.insert(graph_kmer::canonical(code, k), count);
        });
    }

    #[test]
    fn scheduling_length_is_not_selection_and_incomplete_search_cannot_join() {
        use super::super::super::threading::Segment;
        let graph = PathGraph {
            segments: [b"GGGA".as_slice(), b"AAA", b"ACCCA", b"ATG"]
                .into_iter()
                .map(|s| Segment {
                    sequence: s.to_vec(),
                    circular: false,
                })
                .collect(),
            outgoing: vec![vec![1, 2], vec![3], vec![3], vec![]],
        };
        let mut locus = enumerate(&graph, 0, &HashSet::from([0, 3]), 1, 100, 20, 20);
        assert!(locus.routes[0].sequence.len() < locus.routes[1].sequence.len());
        locus.totals = vec![-100.0, 0.0];
        locus.spanning[1] = 20;
        locus.strong_votes[1][0] = 20;
        let settings = GraphSettings {
            coverage_layout: false,
            multi_k_layout: false,
            k: 15,
            min_kmer_count: 3,
            max_error_fraction: 0.2,
            max_bubble_fraction: 0.49,
            min_depth: 10,
            min_junction_fragments: 3,
            max_validation_pairs: 10000,
        };
        assert_eq!(choose(&locus, &settings), Some(1));
        locus.exhausted = true;
        assert_eq!(choose(&locus, &settings), None);
        locus.exhausted = false;
        locus.length_limited = true;
        assert_eq!(choose(&locus, &settings), None);
    }

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
        let mut counts = graph_kmer::Map::default();
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
    fn best_circle_is_selected_independently_of_a_heavier_linear_primary() {
        let out = output_dir("independent-circle");
        let k = 15;
        let linear = representative(
            b"GGCAGGACTCGGATCAGTTGGCCAGTAAGTCCAGCATAGCCAGT",
            false,
            json!({"type":"joined_graph_paths"}),
        );
        let mut best = representative(
            b"ACGTTAGGCATCTGAACCTGTAGTCCGATCAATGCG",
            true,
            json!({"type":"complete_cycle_family","family":1,"score":-10.0,
                "runner_up_score_margin":0.000005,"strict_selected":false}),
        );
        best.regions = vec![
            Region {
                start: 20,
                end: 23,
                evidence: "circle_family_1_block_2".into(),
                reason: "competing_cycle_paths",
                alternatives: 2,
            },
            Region {
                start: 4,
                end: 9,
                evidence: "circle_family_1_block_1".into(),
                reason: "competing_cycle_paths",
                alternatives: 1,
            },
        ];
        let weaker = representative(
            b"TCCGCTACTAACGTCCAGGTCTATGACCGTCATGCATA",
            true,
            json!({"type":"complete_cycle_family","family":2,"score":0.0}),
        );
        let expected_linear = linear.candidate.sequence.clone();
        let expected_circle = best.candidate.sequence.clone();
        let mut counts = graph_kmer::Map::default();
        add_counts(&mut counts, &weaker.candidate, 2, k);
        add_counts(&mut counts, &best.candidate, 8, k);
        add_counts(&mut counts, &linear.candidate, 100, k);
        assert!(
            support_mass(&linear.candidate, &counts, k) > support_mass(&best.candidate, &counts, k)
        );
        assert!(
            support_mass(&best.candidate, &counts, k) > support_mass(&weaker.candidate, &counts, k)
        );
        let (primary, circle) =
            write_outputs(vec![linear, weaker, best], &counts, k, &out, None).unwrap();
        assert_eq!(primary["circular"], false);
        assert_eq!(circle["status"], "selected");
        assert_eq!(circle["eligible_representative_count"], 2);
        assert_eq!(circle["source"]["family"], 1);
        assert_eq!(circle["source"]["runner_up_score_margin"], 0.000005);
        assert_eq!(circle["final_verified"], false);
        assert_eq!(circle["structural_uniqueness_claim"], false);
        assert_eq!(circle["uncertainty_intervals"][0]["start"], 5);
        assert_eq!(circle["uncertainty_intervals"][0]["end"], 9);
        assert_eq!(circle["uncertainty_intervals"][1]["start"], 21);
        assert_eq!(
            crate::io::fasta::load_reference(&out.join("primary.fasta"), None)
                .unwrap()
                .sequence,
            expected_linear
        );
        let fasta =
            crate::io::fasta::load_reference(&out.join("best_circular.fasta"), None).unwrap();
        assert_eq!(fasta.sequence, expected_circle);
        let header = std::fs::read_to_string(out.join("best_circular.fasta")).unwrap();
        assert!(header.starts_with(&format!(
            ">best_circular circular=true length={} final_verified=false\n",
            expected_circle.len()
        )));
        let ambiguity = std::fs::read_to_string(out.join("best_circular_ambiguity.tsv")).unwrap();
        assert!(ambiguity
            .contains("best_circular\t5\t9\tcircle_family_1_block_1\tcompeting_cycle_paths\t1\n"));
        assert!(ambiguity.contains(
            "best_circular\t21\t23\tcircle_family_1_block_2\tcompeting_cycle_paths\t2\n"
        ));
        assert!(!out.join("consensus.fasta").exists());
        std::fs::remove_dir_all(out).unwrap();
    }

    #[test]
    fn exact_circle_support_ties_use_normalized_sequence_independent_of_input_order() {
        let sequences = [
            b"ACGTTAGGCATCTGAACCTGTAGTCCGATCAATGCG".as_slice(),
            b"TCCGCTACTAACGTCCAGGTCTATGACCGTCATGCATA".as_slice(),
        ];
        let expected = sequences.iter().map(|s| normalize(s, true)).min().unwrap();
        let first = output_dir("tie-forward");
        let second = output_dir("tie-reversed");
        for (out, reversed) in [(&first, false), (&second, true)] {
            let mut candidates: Vec<_> = sequences
                .iter()
                .map(|s| representative(s, true, json!({"type":"graph_segment"})))
                .collect();
            if reversed {
                candidates.reverse();
            }
            let (_, circle) =
                write_outputs(candidates, &graph_kmer::Map::default(), 15, out, None).unwrap();
            assert_eq!(circle["canonical_kmer_support_mass"], 0);
            assert_eq!(
                crate::io::fasta::load_reference(&out.join("best_circular.fasta"), None)
                    .unwrap()
                    .sequence,
                expected
            );
        }
        assert_eq!(
            std::fs::read(first.join("best_circular.fasta")).unwrap(),
            std::fs::read(second.join("best_circular.fasta")).unwrap()
        );
        std::fs::remove_dir_all(first).unwrap();
        std::fs::remove_dir_all(second).unwrap();
    }

    #[test]
    fn no_eligible_circle_keeps_linear_primary_without_circular_files() {
        let out = output_dir("no-circle");
        let linear = representative(
            b"ACGTTAGGCATCTGAACCTGTAGTCCGATCAATGCG",
            false,
            json!({"type":"graph_segment"}),
        );
        let short = representative(b"ACGT", true, json!({"type":"graph_segment"}));
        let (primary, circle) = write_outputs(
            vec![linear, short],
            &graph_kmer::Map::default(),
            15,
            &out,
            None,
        )
        .unwrap();
        assert_eq!(primary["circular"], false);
        assert_eq!(circle["status"], "unavailable");
        assert_eq!(circle["reason"], "no_eligible_circular_candidate");
        assert_eq!(circle["eligible_representative_count"], 0);
        assert!(circle["file"].is_null());
        assert!(circle["ambiguity_file"].is_null());
        assert!(circle["source"].is_null());
        assert!(!out.join("best_circular.fasta").exists());
        assert!(!out.join("best_circular_ambiguity.tsv").exists());
        assert!(out.join("primary.fasta").exists());
        std::fs::remove_dir_all(out).unwrap();
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
        let mut alternatives = vec![Primary {
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
        let heavy_circle = representative(
            b"TCCGCTACTAACGTCCAGGTCTATGACCGTCATGCATA",
            true,
            json!({"type":"graph_segment"}),
        );
        let mut counts = graph_kmer::Map::default();
        add_counts(&mut counts, &heavy_circle.candidate, 100, 15);
        assert!(
            support_mass(&heavy_circle.candidate, &counts, 15)
                > support_mass(&final_candidate, &counts, 15)
        );
        alternatives.push(heavy_circle);
        let (report, circle) =
            write_outputs(alternatives, &counts, 15, &out, Some(&final_candidate)).unwrap();
        assert_eq!(report["final_verified"], true);
        assert_eq!(report["circular"], true);
        assert_eq!(circle["status"], "selected");
        assert_eq!(circle["final_verified"], true);
        assert_eq!(circle["source"]["type"], "validated_candidate");
        assert_eq!(
            crate::io::fasta::load_reference(&out.join("best_circular.fasta"), None)
                .unwrap()
                .sequence,
            final_candidate.sequence
        );
        assert_eq!(
            crate::io::fasta::load_reference(&out.join("primary.fasta"), None)
                .unwrap()
                .sequence,
            final_candidate.sequence
        );
        std::fs::remove_dir_all(&out).unwrap();
    }
    fn reciprocal_fixture() -> (PathGraph, Vec<Locus>) {
        use super::super::super::threading::Segment;
        let forward = [b"AGGCA".to_vec(), b"ATCGA".to_vec(), b"AACCT".to_vec()];
        let sequences = forward
            .iter()
            .cloned()
            .chain(forward.iter().map(|s| kmer::reverse_complement(s)));
        let graph = PathGraph {
            segments: sequences
                .map(|sequence| Segment {
                    sequence,
                    circular: false,
                })
                .collect(),
            outgoing: vec![vec![1, 2], vec![], vec![], vec![], vec![3], vec![3]],
        };
        let anchors = HashSet::from([0, 1, 2, 3, 4, 5]);
        let mut loci: Vec<_> = [0, 4, 5]
            .into_iter()
            .map(|from| enumerate(&graph, from, &anchors, 1, 100, 100, 100))
            .collect();
        loci[0].totals = vec![0.0, -100.0];
        loci[0].spanning = vec![100, 3];
        loci[0].strong_votes = vec![vec![0, 100], vec![3, 0]];
        loci[1].spanning[0] = 100;
        loci[1].selected = Some(0);
        loci[2].spanning[0] = 3;
        loci[2].selected = Some(0);
        (graph, loci)
    }

    #[test]
    fn reciprocal_competition_outranks_a_sole_route_without_changing_strict_choices() {
        for order in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            let (graph, mut loci) = reciprocal_fixture();
            let mut choices = vec![Some(0), Some(0), Some(0)];
            // Permute the data, not only the visit order, to detect incidental
            // ordering in reverse-view projection.
            let mut data: Vec<_> = loci.drain(..).zip(choices.drain(..)).map(Some).collect();
            let (loci, mut choices): (Vec<_>, Vec<_>) =
                order.into_iter().map(|i| data[i].take().unwrap()).unzip();
            reciprocal_choices(&graph, &loci, 1, &mut choices);
            for (l, choice) in loci.iter().zip(choices) {
                assert_eq!(choice, if l.from == 5 { None } else { Some(0) });
                if l.from == 5 {
                    assert_eq!(l.selected, Some(0));
                }
            }
        }
    }

    /// Both strand views of a forward segment graph (overlap 1); every segment
    /// except `repeats` is an anchor.
    fn repeat_fixture(
        segments: &[&[u8]],
        forward: &[Vec<usize>],
        repeats: &[usize],
    ) -> (PathGraph, Vec<Locus>) {
        use super::super::super::threading::Segment;
        let n = segments.len();
        let mut sequences: Vec<Vec<u8>> = segments.iter().map(|s| s.to_vec()).collect();
        sequences.extend(segments.iter().map(|s| kmer::reverse_complement(s)));
        let mut outgoing = vec![Vec::new(); 2 * n];
        for (from, targets) in forward.iter().enumerate() {
            for &to in targets {
                outgoing[from].push(to);
                outgoing[n + to].push(n + from);
            }
        }
        let graph = PathGraph {
            segments: sequences
                .into_iter()
                .map(|sequence| Segment {
                    sequence,
                    circular: false,
                })
                .collect(),
            outgoing,
        };
        let anchors: HashSet<_> = (0..2 * n)
            .filter(|&id| !repeats.contains(&(id % n)))
            .collect();
        let mut order: Vec<_> = anchors.iter().copied().collect();
        order.sort_unstable();
        let loci = order
            .iter()
            .map(|&from| enumerate(&graph, from, &anchors, 1, 100, 100, 100))
            .collect();
        (graph, loci)
    }

    fn settings() -> GraphSettings {
        GraphSettings {
            coverage_layout: false,
            multi_k_layout: false,
            k: 15,
            min_kmer_count: 3,
            max_error_fraction: 0.2,
            max_bubble_fraction: 0.49,
            min_depth: 10,
            min_junction_fragments: 3,
            max_validation_pairs: 10000,
        }
    }

    #[test]
    fn two_copy_direct_repeat_pairs_entrance_loop_and_exit_into_one_molecule() {
        let (graph, loci) = repeat_fixture(
            &[b"TTACC", b"CGTAG", b"GAAATC", b"GTCAA"],
            &[vec![1], vec![2, 3], vec![1], vec![]],
            &[1],
        );
        let anchors: Vec<_> = loci.iter().map(|l| l.from).collect();
        let (primaries, choices, reports) =
            local_candidates(&graph, &loci, &anchors, 1, &settings());
        let expected = normalize(b"TTACCGTAGAAATCGTAGTCAA", false);
        let joined: Vec<_> = primaries
            .iter()
            .filter(|p| p.candidate.sequence.len() > 6)
            .collect();
        assert_eq!(joined.len(), 1);
        assert_eq!(joined[0].candidate.sequence, expected);
        let inferred: Vec<_> = joined[0]
            .regions
            .iter()
            .filter(|r| r.reason == "topology_inferred_repeat_pairing")
            .map(|r| (r.start, r.end))
            .collect();
        assert_eq!(inferred, vec![(4, 9), (13, 18)]);
        assert!(joined[0].source["joins"]
            .as_array()
            .unwrap()
            .iter()
            .all(
                |j| j["basis"] == "single_molecule_repeat_pairing" && j["strict_selected"] == false
            ));
        assert_eq!(reports.len(), 1);
        assert_eq!(reports[0]["status"], "applied");
        assert_eq!(reports[0]["alternative_pairing_reads"], "not_observed");
        assert_eq!(choices.iter().filter(|c| c.is_some()).count(), 4);
        assert!(loci.iter().all(|l| l.selected.is_none()));
    }

    #[test]
    fn a_read_supported_route_or_one_qualifying_view_keeps_topology_out() {
        let shape = (
            [b"TTACC".as_slice(), b"CGTAG", b"GAAATC", b"GTCAA"],
            [vec![1], vec![2, 3], vec![1], vec![]],
        );
        let (graph, mut loci) = repeat_fixture(&shape.0, &shape.1, &[1]);
        // Reads favour the route that separates the loop as a subcircle.
        let entrance = loci.iter().position(|l| l.from == 0).unwrap();
        let skip = loci[entrance]
            .routes
            .iter()
            .position(|r| r.to == 3)
            .unwrap();
        loci[entrance].totals[skip] = 100.0;
        loci[entrance].spanning[skip] = 20;
        loci[entrance].strong_votes[skip][1 - skip] = 20;
        let anchors: Vec<_> = loci.iter().map(|l| l.from).collect();
        let (_, choices, reports) = local_candidates(&graph, &loci, &anchors, 1, &settings());
        assert_eq!(choices[entrance], Some(skip));
        assert!(reports.is_empty());
        // An exhausted search in the other strand view of a join blocks it.
        let (graph, mut loci) = repeat_fixture(&shape.0, &shape.1, &[1]);
        let reverse = loci.iter().position(|l| l.from == 6).unwrap();
        loci[reverse].exhausted = true;
        let (_, choices, reports) = local_candidates(&graph, &loci, &anchors, 1, &settings());
        assert!(choices.iter().all(Option::is_none));
        assert_eq!(reports[0]["status"], "conflicting_choice");
    }

    #[test]
    fn a_loop_closed_by_chosen_joins_pairs_despite_a_stub_into_its_end() {
        // a → R → l1 → l2 → w → R → b, with a stub s also entering w.
        let (graph, mut loci) = repeat_fixture(
            &[
                b"TTACC", b"CGTAG", b"GAAAT", b"TCCGA", b"TTGTA", b"AGGC", b"GTCAA",
            ],
            &[
                vec![1],
                vec![2, 6],
                vec![3],
                vec![5],
                vec![5],
                vec![1],
                vec![],
            ],
            &[1, 5],
        );
        let select = |loci: &mut Vec<Locus>, from: usize, to: usize| {
            let li = loci.iter().position(|l| l.from == from).unwrap();
            let ri = loci[li].routes.iter().position(|r| r.to == to).unwrap();
            loci[li].selected = Some(ri);
            loci[li].spanning[ri] = 50;
        };
        select(&mut loci, 0, 2);
        select(&mut loci, 2, 3);
        let anchors: Vec<_> = loci.iter().map(|l| l.from).collect();
        let (primaries, _, reports) = local_candidates(&graph, &loci, &anchors, 1, &settings());
        assert_eq!(reports.len(), 1);
        assert_eq!(reports[0]["status"], "applied");
        assert_eq!(reports[0]["loop_end_segment"], 4);
        let expected = normalize(b"TTACCGTAGAAATCCGAGGCGTAGTCAA", false);
        let joined = primaries
            .iter()
            .find(|p| p.candidate.sequence == expected)
            .expect("one molecule");
        let joins = joined.source["joins"].as_array().unwrap();
        assert_eq!(
            joins
                .iter()
                .filter(|j| j["basis"] == "single_molecule_repeat_pairing")
                .count(),
            1
        );
        assert_eq!(
            joined
                .regions
                .iter()
                .filter(|r| r.reason == "topology_inferred_repeat_pairing")
                .count(),
            1
        );
    }

    #[test]
    fn two_loops_through_one_repeat_close_a_single_circle() {
        let (graph, loci) = repeat_fixture(
            &[b"GAACC", b"CGTAG", b"GATTC"],
            &[vec![1], vec![0, 2], vec![1]],
            &[1],
        );
        let anchors: Vec<_> = loci.iter().map(|l| l.from).collect();
        let (primaries, _, reports) = local_candidates(&graph, &loci, &anchors, 1, &settings());
        assert_eq!(reports[0]["status"], "applied");
        let circles: Vec<_> = primaries.iter().filter(|p| p.candidate.circular).collect();
        assert_eq!(circles.len(), 1);
        assert_eq!(
            circles[0].candidate.sequence,
            normalize(b"GAACCGTAGATTCGTA", true)
        );
    }

    #[test]
    fn conflicting_reverse_only_choices_stay_open_instead_of_last_projection_winning() {
        let (graph, loci) = reciprocal_fixture();
        let mut choices = vec![None, Some(0), Some(0)];
        reciprocal_choices(&graph, &loci, 1, &mut choices);
        assert_eq!(choices, vec![None, None, None]);
    }

    #[test]
    fn reciprocal_read_conflicts_block_topology_in_both_locus_orders() {
        for strict in [false, true] {
            for reverse_order in [false, true] {
                let (graph, mut loci) = repeat_fixture(
                    &[b"TTACC", b"CGTAG", b"GAAATC", b"GTCAA"],
                    &[vec![1], vec![2, 3], vec![1], vec![]],
                    &[1],
                );
                for from in [6, 7] {
                    let l = loci.iter_mut().find(|l| l.from == from).unwrap();
                    let ri = l.routes.iter().position(|r| r.to == 4).unwrap();
                    l.totals[ri] = 100.0;
                    l.spanning[ri] = 20;
                    l.strong_votes[ri][1 - ri] = 20;
                    l.selected = strict.then_some(ri);
                }
                if reverse_order {
                    loci.reverse();
                }
                let mut choices: Vec<_> = loci.iter().map(|l| choose(l, &settings())).collect();
                assert_eq!(choices.iter().filter(|c| c.is_some()).count(), 2);
                let blocked = reciprocal_choices(&graph, &loci, 1, &mut choices);
                assert!(choices.iter().all(Option::is_none));
                assert!(loci
                    .iter()
                    .enumerate()
                    .filter(|(_, l)| [0, 6, 7].contains(&l.from))
                    .all(|(li, _)| blocked.contains(&li)));
                let anchors: Vec<_> = loci.iter().map(|l| l.from).collect();
                let (primaries, choices, reports) =
                    local_candidates(&graph, &loci, &anchors, 1, &settings());
                assert!(choices.iter().all(Option::is_none));
                assert_eq!(reports.len(), 1);
                assert_eq!(reports[0]["status"], "conflicting_choice");
                assert_eq!(reports[0]["alternative_pairing_reads"], "observed");
                assert!(primaries
                    .iter()
                    .flat_map(|p| &p.regions)
                    .all(|r| r.reason != "topology_inferred_repeat_pairing"));
            }
        }
    }

    #[test]
    fn a_rejected_forward_read_choice_stays_open_after_repeat_pairing() {
        let (graph, mut loci) = repeat_fixture(
            &[b"TTACC", b"CGTAG", b"GAAATC", b"GTCAA"],
            &[vec![1], vec![2, 3], vec![1], vec![]],
            &[1],
        );
        for (from, to) in [(0, 3), (7, 6)] {
            let l = loci.iter_mut().find(|l| l.from == from).unwrap();
            let ri = l.routes.iter().position(|r| r.to == to).unwrap();
            l.totals[ri] = 100.0;
            l.spanning[ri] = 20;
            l.strong_votes[ri][1 - ri] = 20;
        }
        let entrance = loci.iter().position(|l| l.from == 0).unwrap();
        let mut choices: Vec<_> = loci.iter().map(|l| choose(l, &settings())).collect();
        let blocked = reciprocal_choices(&graph, &loci, 1, &mut choices);
        assert!(blocked.contains(&entrance));
        assert!(choices[entrance].is_none());
        let anchors: Vec<_> = loci.iter().map(|l| l.from).collect();
        let (_, choices, reports) = local_candidates(&graph, &loci, &anchors, 1, &settings());
        assert!(choices[entrance].is_none());
        assert_eq!(reports[0]["status"], "conflicting_choice");
    }

    #[test]
    fn alternative_pairing_evidence_checks_both_views_without_adding_counts() {
        for (forward, reverse, observed) in
            [(0, 20, true), (20, 0, true), (2, 2, false), (0, 0, false)]
        {
            let (graph, mut loci) = repeat_fixture(
                &[b"TTACC", b"CGTAG", b"GAAATC", b"GTCAA"],
                &[vec![1], vec![2, 3], vec![1], vec![]],
                &[1],
            );
            for (from, to, spanning) in [(0, 3, forward), (7, 4, reverse)] {
                let l = loci.iter_mut().find(|l| l.from == from).unwrap();
                let ri = l.routes.iter().position(|r| r.to == to).unwrap();
                l.spanning[ri] = spanning;
            }
            // Spanning evidence without a winning score/vote margin is not a join.
            assert!(loci.iter().all(|l| choose(l, &settings()).is_none()));
            let anchors: Vec<_> = loci.iter().map(|l| l.from).collect();
            let (_, _, reports) = local_candidates(&graph, &loci, &anchors, 1, &settings());
            assert_eq!(reports.len(), 1);
            assert_eq!(reports[0]["status"], "applied");
            assert_eq!(
                reports[0]["alternative_pairing_reads"],
                if observed { "observed" } else { "not_observed" }
            );
        }
    }
}
