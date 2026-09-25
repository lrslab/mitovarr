//! Resolve a repeat only when reads span its complete unitig and both flanks.
use super::{normalize, Candidate, Unitig};
use crate::{
    cli::{GraphSettings, Limits, Reads},
    io::fastq,
    kmer::graph as graph_kmer,
};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::collections::{BTreeSet, HashMap, HashSet};

#[derive(Clone, Copy)]
struct Run {
    unitig: usize,
    first: usize,
    last: usize,
    read_first: usize,
    read_last: usize,
}

struct Node {
    unitig: usize,
    entry: Option<usize>,
    exit: Option<usize>,
}

#[derive(Clone)]
pub(super) struct Segment {
    pub sequence: Vec<u8>,
    pub circular: bool,
}

#[derive(Clone)]
pub(super) struct PathGraph {
    pub segments: Vec<Segment>,
    pub outgoing: Vec<Vec<usize>>,
}

pub(super) struct Threaded {
    pub candidates: Vec<Candidate>,
    pub unresolved: bool,
    pub report: Value,
    pub graph: PathGraph,
}

type Route = (usize, usize, u32);

fn dominant_routes(routes: &[Route], minimum: u32, fraction: f64) -> Vec<Route> {
    routes
        .iter()
        .copied()
        .filter(|&(a, b, n)| {
            let row_other: u64 = routes
                .iter()
                .filter(|r| r.0 == a && r.1 != b)
                .map(|r| u64::from(r.2))
                .sum();
            let col_other: u64 = routes
                .iter()
                .filter(|r| r.1 == b && r.0 != a)
                .map(|r| u64::from(r.2))
                .sum();
            n >= minimum
                && row_other as f64 <= f64::from(n) * fraction
                && col_other as f64 <= f64::from(n) * fraction
        })
        .collect()
}

pub(super) fn candidates(
    unitigs: &[Unitig],
    reads: &Reads,
    settings: &GraphSettings,
    limits: &Limits,
) -> Result<Threaded> {
    let k = settings.k;
    let anchor_kmers = 5;
    let mut starts: HashMap<graph_kmer::Code, Vec<usize>> = HashMap::new();
    let mut ends: HashMap<graph_kmer::Code, Vec<usize>> = HashMap::new();
    for (id, unitig) in unitigs.iter().enumerate() {
        starts.entry(unitig.from).or_default().push(id);
        ends.entry(unitig.to).or_default().push(id);
    }
    let repeats: HashSet<_> = unitigs
        .iter()
        .enumerate()
        .filter_map(|(id, u)| {
            let incoming = ends.get(&u.from).map_or(0, Vec::len);
            let outgoing = starts.get(&u.to).map_or(0, Vec::len);
            (!u.circular
                && u.from != u.to
                && incoming > 1
                && incoming == outgoing
                && starts[&u.from].len() == 1
                && ends[&u.to].len() == 1)
                .then_some(id)
        })
        .collect();
    let mut support: HashMap<(usize, usize, usize), u32> = HashMap::new();
    let mut names = HashSet::new();
    let mut scan = None;
    if !repeats.is_empty() {
        let mut index = HashMap::new();
        for (id, unitig) in unitigs.iter().enumerate() {
            graph_kmer::visit(&unitig.sequence, None, k, 0, |offset, code| {
                index.insert(code, (id, offset));
            });
        }
        let reverse_ids: Vec<_> = unitigs
            .iter()
            .map(|unitig| {
                let mut last_code = graph_kmer::Code::default();
                graph_kmer::visit(
                    &unitig.sequence[unitig.sequence.len() - k..],
                    None,
                    k,
                    0,
                    |_, code| last_code = code,
                );
                index[&graph_kmer::reverse_code(last_code, k)].0
            })
            .collect();
        scan = Some(fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
            let mut fragment_routes = HashSet::new();
            for read in [a, b] {
                let mut runs: Vec<Run> = Vec::new();
                graph_kmer::visit(
                    &read.sequence,
                    Some(&read.quality),
                    k,
                    limits.min_baseq,
                    |position, code| {
                        if let Some(&(id, offset)) = index.get(&code) {
                            if let Some(last) = runs.last_mut() {
                                if last.unitig == id
                                    && last.last + 1 == offset
                                    && last.read_last + 1 == position
                                {
                                    last.last = offset;
                                    last.read_last = position;
                                    return;
                                }
                            }
                            runs.push(Run {
                                unitig: id,
                                first: offset,
                                last: offset,
                                read_first: position,
                                read_last: position,
                            });
                        }
                    },
                );
                for triple in runs.windows(3) {
                    let (a, r, b) = (triple[0], triple[1], triple[2]);
                    if repeats.contains(&r.unitig)
                        && r.first == 0
                        && r.last + k == unitigs[r.unitig].sequence.len()
                        && a.last + k == unitigs[a.unitig].sequence.len()
                        && b.first == 0
                        && a.read_last + 1 == r.read_first
                        && r.read_last + 1 == b.read_first
                        && a.last - a.first + 1 >= anchor_kmers
                        && b.last - b.first + 1 >= anchor_kmers
                    {
                        fragment_routes.insert((a.unitig, r.unitig, b.unitig));
                        fragment_routes.insert((
                            reverse_ids[b.unitig],
                            reverse_ids[r.unitig],
                            reverse_ids[a.unitig],
                        ));
                    }
                }
            }
            if !fragment_routes.is_empty() {
                let name = a.name(1)?;
                if names.contains(name) {
                    return Ok(());
                }
                ensure!(
                    names.len() < settings.max_validation_pairs,
                    "repeat threading exceeds --max-validation-pairs"
                );
                names.insert(name.to_owned());
                for route in fragment_routes {
                    let count = support.entry(route).or_default();
                    *count = count.saturating_add(1);
                }
            }
            Ok(())
        })?);
    }

    // A bijection is required. A dominant route must win at both its entrance and
    // exit. Rare conflicting read paths can arise from errors or chimeric molecules;
    // retain their counts in the report without giving them an absolute veto.
    let mut resolved = HashMap::new();
    let mut repeat_order: Vec<_> = repeats.iter().copied().collect();
    repeat_order.sort_unstable();
    let mut reported = Vec::new();
    for id in repeat_order {
        let unitig = &unitigs[id];
        let incoming = &ends[&unitig.from];
        let outgoing = &starts[&unitig.to];
        let routes: Vec<_> = incoming
            .iter()
            .flat_map(|&a| {
                outgoing
                    .iter()
                    .filter_map(|&b| {
                        let count = support.get(&(a, id, b)).copied().unwrap_or(0);
                        (count > 0).then_some((a, b, count))
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        let chosen = dominant_routes(
            &routes,
            settings.min_junction_fragments,
            settings.max_error_fraction,
        );
        let inputs: HashSet<_> = chosen.iter().map(|r| r.0).collect();
        let outputs: HashSet<_> = chosen.iter().map(|r| r.1).collect();
        let complete = chosen.len() == incoming.len()
            && inputs.len() == incoming.len()
            && outputs.len() == outgoing.len();
        if complete {
            resolved.insert(id, chosen.clone());
        }
        reported.push(json!({"unitig":format!("u{}",id+1),"length":unitig.sequence.len(),"resolved":complete,
            "observed_routes":routes.iter().map(|(a,b,n)|json!({"from":format!("u{}",a+1),"to":format!("u{}",b+1),"fragments":n,"selected":complete && chosen.contains(&(*a,*b,*n))})).collect::<Vec<_>>()}));
    }

    let mut nodes = Vec::new();
    let mut copies = vec![Vec::new(); unitigs.len()];
    for (id, ids) in copies.iter_mut().enumerate() {
        if let Some(routes) = resolved.get(&id) {
            for &(entry, exit, _) in routes {
                ids.push(nodes.len());
                nodes.push(Node {
                    unitig: id,
                    entry: Some(entry),
                    exit: Some(exit),
                });
            }
        } else {
            ids.push(nodes.len());
            nodes.push(Node {
                unitig: id,
                entry: None,
                exit: None,
            });
        }
    }
    let mut outgoing = vec![Vec::new(); nodes.len()];
    let mut incoming = vec![Vec::new(); nodes.len()];
    for (id, node) in nodes.iter().enumerate() {
        for &next in starts.get(&unitigs[node.unitig].to).into_iter().flatten() {
            if node.exit.is_some_and(|exit| exit != next) {
                continue;
            }
            for &copy in &copies[next] {
                if nodes[copy].entry.is_some_and(|entry| entry != node.unitig) {
                    continue;
                }
                outgoing[id].push(copy);
                incoming[copy].push(id);
            }
        }
    }
    let branch_nodes = (0..nodes.len())
        .filter(|&id| incoming[id].len() > 1 || outgoing[id].len() > 1)
        .count();
    let mut visited = vec![false; nodes.len()];
    let mut unique = BTreeSet::new();
    let mut segments = Vec::new();
    let mut owner = vec![usize::MAX; nodes.len()];
    let mut terminal = Vec::new();
    // Start at boundaries, then consume remaining closed cycles.
    let starts: Vec<_> = (0..nodes.len())
        .filter(|&id| incoming[id].len() != 1 || outgoing[incoming[id][0]].len() != 1)
        .chain(0..nodes.len())
        .collect();
    for first in starts {
        if visited[first] {
            continue;
        }
        let mut sequence = unitigs[nodes[first].unitig].sequence.clone();
        let mut node = first;
        let mut circular = false;
        loop {
            visited[node] = true;
            owner[node] = segments.len();
            if outgoing[node].len() != 1 {
                break;
            }
            let next = outgoing[node][0];
            if incoming[next].len() != 1 {
                break;
            }
            if visited[next] {
                circular = next == first;
                break;
            }
            sequence.extend_from_slice(&unitigs[nodes[next].unitig].sequence[k - 1..]);
            node = next;
        }
        if circular {
            sequence.truncate(sequence.len() - k + 1);
        }
        unique.insert((normalize(&sequence, circular), circular));
        terminal.push(node);
        segments.push(Segment { sequence, circular });
    }
    let candidates = unique
        .into_iter()
        .map(|(sequence, circular)| Candidate { sequence, circular })
        .collect();
    let outgoing = terminal
        .iter()
        .map(|&last| outgoing[last].iter().map(|&n| owner[n]).collect())
        .collect();
    Ok(Threaded {
        candidates,
        unresolved: branch_nodes > 0,
        graph: PathGraph { segments, outgoing },
        report: json!({"method":"exact read paths spanning an entire repeat unitig and both flanks; dominant one-to-one entrance/exit routes", "anchor_kmers_per_flank":anchor_kmers,"minimum_bridge_fragments":settings.min_junction_fragments,"max_alternative_route_ratio":settings.max_error_fraction,"scan":scan,"bridging_fragments":names.len(),"directed_repeat_unitigs":repeats.len(),"resolved_directed_repeat_unitigs":resolved.len(),"unresolved_unitig_junctions":branch_nodes,"repeats":reported}),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weak_alternatives_still_count_towards_relative_route_support() {
        assert!(dominant_routes(&[(0, 2, 3), (0, 3, 2), (1, 2, 2), (1, 3, 3)], 3, 0.2).is_empty());
        assert_eq!(
            dominant_routes(&[(0, 2, 100), (0, 3, 2), (1, 2, 1), (1, 3, 80)], 3, 0.2),
            vec![(0, 2, 100), (1, 3, 80)]
        );
        // Individually small alternatives can collectively make the route ambiguous.
        assert!(dominant_routes(&[(0, 2, 10), (0, 3, 2), (0, 4, 2)], 3, 0.2).is_empty());
    }
}
