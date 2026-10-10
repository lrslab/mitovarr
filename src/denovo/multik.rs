//! Bounded cross-k hypotheses over the same sample reads. Long k-mer contexts
//! schedule additional walks; they are never independent votes or hard vetoes.
use super::{compact_at, Edge, Unitig};
use crate::{
    cli::Reads,
    io::fastq,
    kmer::{self, graph as g},
};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

pub(crate) const KS: [usize; 7] = [15, 21, 31, 47, 63, 95, 127];
pub(crate) type Graphs = BTreeMap<usize, HashMap<g::Code, u32>>;
pub(crate) type Counts = BTreeMap<kmer::Code, usize>;
const STATES: usize = 50_000;
const PENDING: usize = 1024;

pub(crate) fn occurrences(sequence: &[u8], k: usize) -> Counts {
    let mut result = BTreeMap::new();
    kmer::visit(sequence, None, k, 0, |_, c| {
        *result.entry(kmer::canonical(c, k)).or_default() += 1;
    });
    result
}

pub(crate) fn add_pair(
    graphs: &mut Graphs,
    a: &fastq::Record,
    b: &fastq::Record,
    minq: u8,
    maximum: usize,
) -> Result<()> {
    for (&k, counts) in graphs {
        let mut observed = HashSet::new();
        for read in [a, b] {
            g::visit(&read.sequence, Some(&read.quality), k, minq, |_, c| {
                observed.insert(g::canonical(c, k));
            });
        }
        for c in observed {
            ensure!(
                counts.contains_key(&c) || counts.len() < maximum,
                "multi-k graph exceeds {maximum} distinct kmers at k={k}"
            );
            let n = counts.entry(c).or_default();
            *n = n.saturating_add(1);
        }
    }
    Ok(())
}

pub(crate) fn build(reads: &Reads, minq: u8, max_pairs: usize) -> Result<Graphs> {
    let mut graphs = KS.into_iter().map(|k| (k, HashMap::new())).collect();
    let mut names = HashSet::new();
    fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
        if !names.insert(a.name(1)?.to_owned()) {
            return Ok(());
        }
        ensure!(
            names.len() <= max_pairs,
            "multi-k layout exceeds named FASTQ pair limit"
        );
        add_pair(&mut graphs, a, b, minq, 2_000_000)
    })?;
    Ok(graphs)
}

pub(crate) struct Guide<'a> {
    pub k: usize,
    pub observed: &'a HashMap<kmer::Code, u64>,
    pub outside: Counts,
    pub mean: f64,
}
impl Guide<'_> {
    fn expected(&self, c: kmer::Code) -> f64 {
        self.observed.get(&c).copied().unwrap_or(0) as f64 / self.mean
    }
    fn increment(&self, c: kmer::Code, used: usize, n: usize) -> f64 {
        let before = (self.outside.get(&c).copied().unwrap_or(0) + used) as f64;
        (before + n as f64 - self.expected(c)).powi(2) - (before - self.expected(c)).powi(2)
    }
    fn bound(&self, c: kmer::Code) -> usize {
        ((1.5 * self.expected(c)).ceil() as usize + 2)
            .min(4096)
            .saturating_sub(self.outside.get(&c).copied().unwrap_or(0))
            .max(2)
    }
}

fn code(sequence: &[u8]) -> Option<g::Code> {
    let mut result = None;
    if !sequence.is_empty() && sequence.len() <= 127 {
        g::visit(sequence, None, sequence.len(), 0, |_, c| result = Some(c));
    }
    result
}

// The target is the END of the complete long right flank. Its outgoing edges
// remain available when the same short node is encountered inside a repeat.
fn corridor(
    counts: &HashMap<g::Code, u32>,
    k: usize,
    minimum: u32,
    start: g::Code,
    target: g::Code,
    guards: &[Vec<u8>],
) -> Vec<Unitig> {
    let mut split = HashSet::from([start, target]);
    for guard in guards {
        if guard.len() >= k - 1 {
            if let Some(c) = code(&guard[guard.len() - (k - 1)..]) {
                split.insert(c);
            }
        }
    }
    let mut edges = Vec::new();
    for (&c, &n) in counts.iter().filter(|(_, n)| **n >= minimum) {
        for c in [c, g::reverse_code(c, k)] {
            edges.push(Edge {
                code: c,
                from: c.prefix(),
                to: c.suffix(k - 1),
                count: n,
            });
        }
    }
    edges.sort_by_key(|e| e.code);
    let (units, _) = compact_at(&edges, k, &split);
    let mut incoming = HashMap::<_, Vec<_>>::new();
    let mut outgoing = HashMap::<_, Vec<_>>::new();
    for (id, u) in units.iter().enumerate() {
        incoming.entry(u.to).or_default().push(id);
        outgoing.entry(u.from).or_default().push(id);
    }
    let mut reaches = HashSet::from([target]);
    let mut pending = vec![target];
    while let Some(n) = pending.pop() {
        for &i in incoming.get(&n).into_iter().flatten() {
            if reaches.insert(units[i].from) {
                pending.push(units[i].from);
            }
        }
    }
    let mut accessible = HashSet::from([start]);
    pending.push(start);
    while let Some(n) = pending.pop() {
        for &i in outgoing.get(&n).into_iter().flatten() {
            if reaches.contains(&units[i].to) && accessible.insert(units[i].to) {
                pending.push(units[i].to);
            }
        }
    }
    units
        .into_iter()
        .filter(|u| accessible.contains(&u.from) && reaches.contains(&u.to))
        .collect()
}

/// At a branch, use the longest available read window that distinguishes at
/// least one outgoing base. No observed extension means no preference, and
/// unobserved alternatives are still explored by the independent depth queue.
fn context_penalties(
    graphs: &Graphs,
    graph_k: usize,
    minimum: u32,
    sequence: &[u8],
    extensions: &[u8],
) -> (Vec<usize>, Option<usize>) {
    for (&k, counts) in graphs.iter().rev().filter(|(k, _)| **k > graph_k) {
        if sequence.len() + 1 < k {
            continue;
        }
        let mut word = sequence[sequence.len() + 1 - k..].to_vec();
        word.push(b'A');
        let supported: Vec<_> = extensions
            .iter()
            .map(|&base| {
                word[k - 1] = base;
                code(&word)
                    .and_then(|c| counts.get(&g::canonical(c, k)))
                    .copied()
                    .unwrap_or(0)
                    >= minimum
            })
            .collect();
        if supported.iter().any(|x| *x) && supported.iter().any(|x| !*x) {
            return (supported.iter().map(|x| usize::from(!x)).collect(), Some(k));
        }
    }
    (vec![0; extensions.len()], None)
}

#[derive(Clone)]
struct State {
    // The first priority is nonzero only for the context-guided schedule.
    penalty: usize,
    cost: f64,
    node: g::Code,
    sequence: Vec<u8>,
    used: Counts,
    support: u32,
    walk: Vec<usize>,
}
impl PartialEq for State {
    fn eq(&self, b: &Self) -> bool {
        self.cmp(b) == Ordering::Equal
    }
}
impl Eq for State {}
impl PartialOrd for State {
    fn partial_cmp(&self, b: &Self) -> Option<Ordering> {
        Some(self.cmp(b))
    }
}
impl Ord for State {
    fn cmp(&self, b: &Self) -> Ordering {
        self.penalty
            .cmp(&b.penalty)
            .then(self.cost.total_cmp(&b.cost))
            .then(self.sequence.len().cmp(&b.sequence.len()))
            .then(self.sequence.cmp(&b.sequence))
            .then(self.walk.cmp(&b.walk))
    }
}

pub(crate) struct Candidate {
    pub sequence: Vec<u8>,
    pub support: u32,
    pub walk: Vec<usize>,
    pub schedules: BTreeSet<String>,
}

fn write_graph(path: &Path, units: &[Unitig], k: usize) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "H\tVN:Z:1.0")?;
    let mut outgoing = HashMap::<_, Vec<_>>::new();
    for (i, u) in units.iter().enumerate() {
        outgoing.entry(u.from).or_default().push(i);
        writeln!(
            writer,
            "S\tu{}\t{}\tFC:f:{:.6}",
            i + 1,
            String::from_utf8_lossy(&u.sequence),
            u.mean_count
        )?;
    }
    for (i, u) in units.iter().enumerate() {
        for &j in outgoing.get(&u.to).into_iter().flatten() {
            writeln!(writer, "L\tu{}\t+\tu{}\t+\t{}M", i + 1, j + 1, k - 1)?;
        }
    }
    writer.flush()?;
    Ok(())
}

pub(crate) struct Request<'a> {
    pub graphs: &'a Graphs,
    pub k: usize,
    pub minimum: u32,
    pub left: &'a [u8],
    pub right: &'a [u8],
    pub guard: Option<&'a [u8]>,
    pub guide: &'a Guide<'a>,
    pub max_span: usize,
    pub max_paths: usize,
    pub graph_file: Option<&'a Path>,
}

pub(crate) fn search(request: Request<'_>) -> Result<(Vec<Candidate>, Value)> {
    let Request {
        graphs,
        k,
        minimum,
        left,
        right,
        guard,
        guide,
        max_span,
        max_paths,
        graph_file,
    } = request;
    let mut report = json!({"k":k,"minimum_pair_support":minimum,"search_complete":false,
        "method":"long_flanks_cross_k_read_context","left_flank_bases":left.len(),"right_flank_bases":right.len(),
        "max_states_per_k":STATES,"max_pending":PENDING,"max_span":max_span,"max_candidates":2*max_paths.div_ceil(2).max(1)+1,
        "depth_count_k":guide.k,"depth_counting":"complete spelled suffix, including cross-edge windows",
        "context_rule":"longest observed distinguishing next-base window; scheduling only; absence never removes an edge",
        "context_counts_are_independent_evidence":false,"retained_candidates_survive_other_k_limits":true,
        "limitations":["bounded schedules and coverage copy limits","exact read contexts stop at k127","full backbone guard excluded when present"]});
    if left.len() < k - 1 || right.len() < k - 1 || graphs[&k].is_empty() {
        report["status"] = json!("insufficient_read_or_flank_context");
        return Ok((vec![], report));
    }
    let Some(start) = code(&left[left.len() - (k - 1)..]) else {
        report["status"] = json!("ambiguous_left_flank");
        return Ok((vec![], report));
    };
    let Some(target) = code(&right[right.len() - (k - 1)..]) else {
        report["status"] = json!("ambiguous_right_flank");
        return Ok((vec![], report));
    };
    let support = left
        .windows(k)
        .map(|s| {
            code(s)
                .and_then(|c| graphs[&k].get(&g::canonical(c, k)))
                .copied()
                .unwrap_or(0)
        })
        .min()
        .unwrap_or(u32::MAX);
    if support < minimum {
        report["status"] = json!("unsupported_complete_left_flank");
        return Ok((vec![], report));
    }
    let guards: Vec<_> = guard
        .into_iter()
        .flat_map(|g| [g.to_vec(), kmer::reverse_complement(g)])
        .collect();
    let units = corridor(&graphs[&k], k, minimum, start, target, &guards);
    if let Some(path) = graph_file {
        write_graph(path, &units, k)?;
    }
    let mut outgoing = HashMap::<_, Vec<_>>::new();
    for (i, u) in units.iter().enumerate() {
        outgoing.entry(u.from).or_default().push(i);
    }
    let mut found = BTreeMap::<Vec<u8>, Candidate>::new();
    let mut schedules = Vec::new();
    for schedule in [
        "nearby_exit_seed",
        "coverage_residual",
        "cross_k_read_context",
    ] {
        let seed = schedule == "nearby_exit_seed";
        let context = schedule == "cross_k_read_context";
        let budget = if seed { STATES / 10 } else { STATES * 9 / 20 };
        let completions_max = if seed {
            1
        } else {
            max_paths.div_ceil(2).max(1)
        };
        let used = occurrences(left, guide.k);
        let cost = if seed {
            left.len() as f64
        } else {
            used.iter().map(|(&c, &n)| guide.increment(c, 0, n)).sum()
        };
        let mut pending = BTreeSet::from([State {
            penalty: 0,
            cost,
            node: start,
            sequence: left.to_vec(),
            used,
            support,
            walk: vec![],
        }]);
        let mut seen_seed = HashSet::new();
        let (mut states, mut completions, mut peak, mut contextual_branches) = (0, 0, 1, 0);
        let mut context_ks = BTreeSet::new();
        let mut cuts = BTreeSet::new();
        while let Some(state) = pending.pop_first() {
            if states >= budget {
                cuts.insert("states");
                break;
            }
            states += 1;
            if state.sequence.len() > left.len()
                && state.node == target
                && state.sequence.ends_with(right)
            {
                let entry = found
                    .entry(state.sequence.clone())
                    .or_insert_with(|| Candidate {
                        sequence: state.sequence.clone(),
                        support: state.support,
                        walk: state.walk.clone(),
                        schedules: BTreeSet::new(),
                    });
                entry.schedules.insert(schedule.to_owned());
                completions += 1;
                if completions >= completions_max {
                    if !pending.is_empty() {
                        cuts.insert("completions");
                    }
                    break;
                }
                continue;
            }
            if seed {
                let tail = state.sequence[state
                    .sequence
                    .len()
                    .saturating_sub(right.len().max(guard.map_or(0, |g| g.len())))..]
                    .to_vec();
                if !seen_seed.insert((state.node, tail)) {
                    continue;
                }
            }
            let next = outgoing.get(&state.node).map(Vec::as_slice).unwrap_or(&[]);
            let (penalties, context_k) = if context && next.len() > 1 {
                context_penalties(
                    graphs,
                    k,
                    minimum,
                    &state.sequence,
                    &next
                        .iter()
                        .map(|&i| units[i].sequence[k - 1])
                        .collect::<Vec<_>>(),
                )
            } else {
                (vec![0; next.len()], None)
            };
            if let Some(k) = context_k {
                contextual_branches += 1;
                context_ks.insert(k);
            }
            for (position, &i) in next.iter().enumerate() {
                let u = &units[i];
                let appended = &u.sequence[k - 1..];
                if state.sequence.len() + appended.len() > max_span {
                    cuts.insert("span");
                    continue;
                }
                let mut child = state.clone();
                child.sequence.extend_from_slice(appended);
                if guards.iter().any(|g| child.sequence.ends_with(g)) {
                    cuts.insert("full_backbone_guard");
                    continue;
                }
                child.node = u.to;
                child.walk.push(i);
                child.support = child.support.min(u.min_count);
                child.penalty += penalties[position];
                if seed {
                    child.cost = child.sequence.len() as f64;
                } else {
                    let from = state.sequence.len().saturating_sub(guide.k - 1);
                    let additions = occurrences(&child.sequence[from..], guide.k);
                    let mut allowed = true;
                    for (&c, &n) in &additions {
                        let before = child.used.get(&c).copied().unwrap_or(0);
                        if before + n > guide.bound(c) {
                            cuts.insert("copy_bound");
                            allowed = false;
                            break;
                        }
                        child.cost += guide.increment(c, before, n);
                        *child.used.entry(c).or_default() += n;
                    }
                    if !allowed {
                        continue;
                    }
                }
                pending.insert(child);
                if pending.len() > PENDING {
                    pending.pop_last();
                    cuts.insert("pending");
                }
                peak = peak.max(pending.len());
            }
        }
        schedules.push(
            json!({"schedule":schedule,"states":states,"state_budget":budget,
            "completions":completions,"peak_pending":peak,"cuts":cuts,
            "contextual_branch_visits":contextual_branches,"context_k":context_ks}),
        );
    }
    report["status"] = json!(if found.is_empty() {
        "no_layout_within_search_bounds"
    } else {
        "inferred_candidates"
    });
    report["corridor_segments"] = json!(units.len());
    report["candidates"] = json!(found.len());
    report["schedules"] = json!(schedules);
    Ok((found.into_values().collect(), report))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(n: usize, mut seed: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed % 4) as usize]
            })
            .collect()
    }
    fn graphs(sequences: &[Vec<u8>]) -> Graphs {
        KS.into_iter()
            .map(|k| {
                let mut counts = HashMap::new();
                for s in sequences {
                    g::visit(s, None, k, 0, |_, c| {
                        counts.insert(g::canonical(c, k), 10);
                    });
                }
                (k, counts)
            })
            .collect()
    }
    fn observed(sequence: &[u8], depth: u64) -> HashMap<kmer::Code, u64> {
        occurrences(sequence, 31)
            .into_iter()
            .map(|(c, n)| (c, n as u64 * depth))
            .collect()
    }
    fn run(
        graphs: &Graphs,
        sequence: &[u8],
        k: usize,
        guide: &Guide<'_>,
        paths: usize,
    ) -> (Vec<Candidate>, Value) {
        search(Request {
            graphs,
            k,
            minimum: 3,
            left: &sequence[..126],
            right: &sequence[sequence.len() - 126..],
            guard: None,
            guide,
            max_span: sequence.len() + 300,
            max_paths: paths,
            graph_file: None,
        })
        .unwrap()
    }

    #[test]
    fn long_exit_identity_survives_an_internal_short_target_and_high_k() {
        let left = dna(126, 9);
        let right = dna(126, 199);
        let sequence = [left, dna(41, 87), right[112..].to_vec(), dna(61, 45), right].concat();
        let graphs = graphs(std::slice::from_ref(&sequence));
        let counts = observed(&sequence, 10);
        let guide = Guide {
            k: 31,
            observed: &counts,
            outside: Counts::new(),
            mean: 10.0,
        };
        for k in [15, 31, 95, 127] {
            let (paths, report) = run(&graphs, &sequence, k, &guide, 32);
            assert!(
                paths.iter().any(|p| p.sequence == sequence),
                "k={k}: {report}"
            );
            assert!(paths
                .iter()
                .all(|p| p.sequence.ends_with(&sequence[sequence.len() - 126..])));
            assert_eq!(report["search_complete"], false);
            assert!(
                report["schedules"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|s| s["states"].as_u64().unwrap())
                    .sum::<u64>()
                    <= STATES as u64
            );
        }
    }

    #[test]
    fn continuous_high_k_context_proposes_the_phased_route_without_vetoing_other_routes() {
        let left = dna(126, 5);
        let other = dna(126, 67);
        let repeat = dna(30, 84);
        let right = dna(126, 19);
        let mut a = dna(70, 899);
        let mut b = dna(70, 912);
        a[0] = b'A';
        b[0] = b'T';
        let truth = [left.clone(), repeat.clone(), a, right.clone()].concat();
        let decoy_source = [other, repeat.clone(), b.clone(), right.clone()].concat();
        let false_join = [left, repeat, b, right].concat();
        let graphs = graphs(&[truth.clone(), decoy_source]);
        let mut counts = observed(&false_join, 30);
        for (c, n) in observed(&truth, 3) {
            *counts.entry(c).or_default() += n;
        }
        let guide = Guide {
            k: 31,
            observed: &counts,
            outside: Counts::new(),
            mean: 30.0,
        };
        let (paths, report) = run(&graphs, &truth, 15, &guide, 1);
        assert!(
            paths
                .iter()
                .any(|p| p.sequence == truth && p.schedules.contains("cross_k_read_context")),
            "{report}"
        );
        assert!(
            paths.iter().any(|p| p.sequence == false_join),
            "absence of high-k evidence must not remove the competing route"
        );
        let guided = &report["schedules"][2];
        assert!(guided["contextual_branch_visits"].as_u64().unwrap() > 0);
        let (none, k) = context_penalties(
            &Graphs::from([(63, HashMap::new())]),
            15,
            3,
            &truth[..150],
            b"AC",
        );
        assert_eq!(none, vec![0, 0]);
        assert_eq!(k, None);
    }

    #[test]
    fn small_k_walk_counts_cross_edge_31mers_and_recovers_three_repeat_copies() {
        let left = dna(126, 33);
        let repeat = dna(44, 771);
        let right = dna(126, 39);
        let truth = [left, repeat.repeat(3), right].concat();
        let graphs = graphs(std::slice::from_ref(&truth));
        let counts = observed(&truth, 20);
        let guide = Guide {
            k: 31,
            observed: &counts,
            outside: Counts::new(),
            mean: 20.0,
        };
        let (paths, report) = run(&graphs, &truth, 15, &guide, 64);
        let best = paths
            .iter()
            .min_by_key(|p| {
                let mut domain: BTreeSet<_> = counts.keys().copied().collect();
                let c = occurrences(&p.sequence, 31);
                domain.extend(c.keys().copied());
                domain
                    .iter()
                    .map(|code| {
                        let a = c.get(code).copied().unwrap_or(0) as i64;
                        let b = counts.get(code).copied().unwrap_or(0) as i64 / 20;
                        (a - b).pow(2)
                    })
                    .sum::<i64>()
            })
            .expect("candidate");
        assert_eq!(best.sequence, truth, "{report}");
        assert!(paths.iter().any(|p| p.sequence.len() < truth.len()));
        assert!(paths.iter().any(|p| p.sequence.len() > truth.len()));
    }

    #[test]
    fn a_true_gap_stays_disconnected_even_with_high_k_context_elsewhere() {
        let sequence = dna(460, 617);
        let mut graphs = graphs(std::slice::from_ref(&sequence));
        let edge = code(&sequence[240..261]).unwrap();
        graphs.get_mut(&21).unwrap().remove(&g::canonical(edge, 21));
        let counts = observed(&sequence, 10);
        let guide = Guide {
            k: 31,
            observed: &counts,
            outside: Counts::new(),
            mean: 10.0,
        };
        assert!(run(&graphs, &sequence, 21, &guide, 32).0.is_empty());
    }

    #[test]
    fn input_order_and_reverse_complement_preserve_candidate_spellings() {
        let truth = [dna(126, 171), dna(44, 188).repeat(2), dna(126, 161)].concat();
        let reverse = kmer::reverse_complement(&truth);
        let a = graphs(&[truth.clone(), reverse.clone()]);
        let b = graphs(&[reverse.clone(), truth.clone()]);
        let counts = observed(&truth, 10);
        let guide = Guide {
            k: 31,
            observed: &counts,
            outside: Counts::new(),
            mean: 10.0,
        };
        let first = run(&a, &truth, 31, &guide, 64)
            .0
            .into_iter()
            .map(|p| p.sequence)
            .collect::<BTreeSet<_>>();
        let second = run(&b, &truth, 31, &guide, 64)
            .0
            .into_iter()
            .map(|p| p.sequence)
            .collect::<BTreeSet<_>>();
        assert_eq!(first, second);
        let reversed = run(&b, &reverse, 31, &guide, 64)
            .0
            .into_iter()
            .map(|p| kmer::reverse_complement(&p.sequence))
            .collect::<BTreeSet<_>>();
        assert_eq!(first, reversed);
    }
}
