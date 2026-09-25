//! Small-k paths repair missing connections between existing graph ends.
//! This proposes sequence only; the shared fragment scorer decides each join.
use super::*;
use crate::denovo::threading::Segment;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

type Bridge = (usize, usize, Vec<u8>);

struct Search {
    bridges: BTreeSet<Bridge>,
    states: usize,
    limited: bool,
}

fn code(sequence: &[u8]) -> kmer::Code {
    let mut result = 0;
    kmer::visit(sequence, None, sequence.len(), 0, |_, c| result = c);
    result
}

fn search(
    from: usize,
    graph: &PathGraph,
    targets: &HashMap<kmer::Code, Vec<usize>>,
    outgoing: &HashMap<kmer::Code, Vec<(kmer::Code, u8)>>,
    k: usize,
    overlap: usize,
    max_bases: usize,
) -> Search {
    let source = &graph.segments[from].sequence;
    let initial = source[source.len() - overlap..].to_vec();
    let start = code(&source[source.len() - k + 1..]);
    let mut stack = vec![(start, initial)];
    let mut result = Search {
        bridges: BTreeSet::new(),
        states: 0,
        limited: false,
    };
    while let Some((node, sequence)) = stack.pop() {
        result.states += 1;
        if result.states > MAX_STATES {
            result.limited = true;
            break;
        }
        let mut arrived = false;
        if sequence.len() > overlap {
            if let Some(ids) = targets.get(&node) {
                for &to in ids {
                    if sequence.ends_with(&graph.segments[to].sequence[..overlap]) {
                        result.bridges.insert((from, to, sequence.clone()));
                        arrived = true;
                    }
                }
            }
        }
        if arrived {
            if result.bridges.len() >= MAX_PATHS && !stack.is_empty() {
                result.limited = true;
                break;
            }
            continue;
        }
        if let Some(next) = outgoing.get(&node) {
            if sequence.len() - overlap >= max_bases {
                result.limited = true;
                continue;
            }
            for &(to, base) in next.iter().rev() {
                let mut extended = sequence.clone();
                extended.push(base);
                stack.push((to, extended));
            }
        }
    }
    result
}

pub(super) fn augment(
    graph: &PathGraph,
    reads: &Reads,
    out: &Path,
    settings: &GraphSettings,
    limits: &Limits,
    max_bases: usize,
) -> Result<(PathGraph, Value)> {
    let overlap = settings.k - 1;
    let mut incoming = vec![0; graph.segments.len()];
    for targets in &graph.outgoing {
        for &id in targets {
            incoming[id] += 1;
        }
    }
    let mut occurrence = HashMap::new();
    for s in &graph.segments {
        kmer::visit(&s.sequence, None, 31, 0, |_, c| {
            *occurrence.entry(kmer::canonical(c, 31)).or_default() += 1;
        });
    }
    let anchors = path_anchors(graph, settings.k, &occurrence);
    let mut sources: Vec<_> = anchors
        .iter()
        .copied()
        .filter(|&id| graph.outgoing[id].is_empty())
        .collect();
    let mut sinks: Vec<_> = anchors
        .iter()
        .copied()
        .filter(|&id| incoming[id] == 0)
        .collect();
    sources.sort_unstable();
    sinks.sort_unstable();
    let mut proposed = BTreeMap::<Bridge, BTreeSet<usize>>::new();
    let mut limited_sources = BTreeSet::new();
    let mut searches = Vec::new();
    for small_k in [47, 31].into_iter().filter(|&k| k < settings.k) {
        if sources.is_empty() || sinks.is_empty() {
            break;
        }
        let mut counts = HashMap::<kmer::Code, u32>::new();
        fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
            let mut fragment = HashSet::new();
            for read in [a, b] {
                kmer::visit(
                    &read.sequence,
                    Some(&read.quality),
                    small_k,
                    limits.min_baseq,
                    |_, c| {
                        fragment.insert(kmer::canonical(c, small_k));
                    },
                );
            }
            for c in fragment {
                ensure!(
                    counts.contains_key(&c) || counts.len() < limits.max_kmers,
                    "end assembly exceeds --max-kmers at k={small_k}"
                );
                let n = counts.entry(c).or_default();
                *n = n.saturating_add(1);
            }
            Ok(())
        })?;
        let mask = (1u128 << (2 * (small_k - 1))) - 1;
        let mut outgoing = HashMap::<_, Vec<_>>::new();
        for (&c, &count) in &counts {
            if count < settings.min_kmer_count {
                continue;
            }
            for c in [c, kmer::reverse_code(c, small_k)] {
                outgoing
                    .entry(c >> 2)
                    .or_default()
                    .push((c & mask, b"ACGT"[(c & 3) as usize]));
            }
        }
        for next in outgoing.values_mut() {
            next.sort_unstable();
            next.dedup();
        }
        let mut targets = HashMap::<_, Vec<_>>::new();
        for &to in &sinks {
            targets
                .entry(code(
                    &graph.segments[to].sequence[overlap - small_k + 1..overlap],
                ))
                .or_default()
                .push(to);
        }
        for &from in &sources {
            let result = search(
                from, graph, &targets, &outgoing, small_k, overlap, max_bases,
            );
            searches.push(
                json!({"k":small_k,"source_segment":from+1,"states":result.states,
                "candidate_connections":result.bridges.len(),"search_limited":result.limited}),
            );
            if result.limited {
                limited_sources.insert(from);
            }
            for bridge in result.bridges {
                proposed.entry(bridge).or_default().insert(small_k);
            }
        }
    }
    // Require the reciprocal physical spelling too. A capped search in either
    // direction cannot silently leave the other direction as a sole candidate.
    let mut by_sequence = HashMap::<_, Vec<_>>::new();
    for (id, s) in graph.segments.iter().enumerate() {
        by_sequence.entry(s.sequence.clone()).or_default().push(id);
    }
    let retained: Vec<_> = proposed
        .iter()
        .filter(|((from, to, sequence), _)| {
            if limited_sources.contains(from) {
                return false;
            }
            let Some(reverse_from) =
                by_sequence.get(&kmer::reverse_complement(&graph.segments[*to].sequence))
            else {
                return false;
            };
            let Some(reverse_to) =
                by_sequence.get(&kmer::reverse_complement(&graph.segments[*from].sequence))
            else {
                return false;
            };
            let reverse = kmer::reverse_complement(sequence);
            reverse_from.iter().any(|&a| {
                !limited_sources.contains(&a)
                    && reverse_to
                        .iter()
                        .any(|&b| proposed.contains_key(&(a, b, reverse.clone())))
            })
        })
        .collect();
    let mut augmented = graph.clone();
    let mut records = Vec::new();
    let mut connections = Vec::new();
    for ((from, to, sequence), ks) in retained {
        let id = augmented.segments.len();
        augmented.segments.push(Segment {
            sequence: sequence.clone(),
            circular: false,
        });
        augmented.outgoing[*from].push(id);
        augmented.outgoing.push(vec![*to]);
        records.push((
            format!(
                "end_bridge_{} source_segment={} target_segment={} candidate_only=true",
                id + 1,
                from + 1,
                to + 1
            ),
            sequence.clone(),
        ));
        connections.push(
            json!({"segment":id+1,"source_segment":from+1,"target_segment":to+1,
            "sequence_length":sequence.len(),"k_values":ks}),
        );
    }
    crate::io::fasta::write_fasta_records(&out.join("end_bridge_candidates.fasta"), &records)?;
    let mut gfa = BufWriter::new(File::create(out.join("path_graph.gfa"))?);
    writeln!(gfa, "H\tVN:Z:1.0")?;
    for (id, s) in augmented.segments.iter().enumerate() {
        writeln!(
            gfa,
            "S\t{}\t{}\tTP:Z:{}",
            id + 1,
            String::from_utf8_lossy(&s.sequence),
            if id < graph.segments.len() {
                "sample_graph"
            } else {
                "end_bridge_candidate"
            }
        )?;
    }
    for (from, targets) in augmented.outgoing.iter().enumerate() {
        for &to in targets {
            // Threading already removed the closing overlap from an isolated
            // circular segment. Its self-link therefore consumes zero bases.
            let shared = if from == to && augmented.segments[from].circular {
                0
            } else {
                overlap
            };
            writeln!(gfa, "L\t{}\t+\t{}\t+\t{}M", from + 1, to + 1, shared)?;
        }
    }
    gfa.flush()?;
    Ok((
        augmented,
        json!({"method":"sample-only small-k paths between unique graph ends; full main-k overlap at both flanks; reciprocal candidate spellings; joint fragment scoring decides joins",
        "k_values":([47,31].into_iter().filter(|&k|k<settings.k).collect::<Vec<_>>()),
        "source_ends":sources.len(),"target_ends":sinks.len(),"maximum_added_bases":max_bases,
        "maximum_states_per_end_per_k":MAX_STATES,"maximum_candidates_per_end_per_k":MAX_PATHS,
        "limited_sources":limited_sources.iter().map(|id|id+1).collect::<Vec<_>>(),
        "incomplete_search_policy":"discard proposals from both views of any limited source; retain the open end",
        "candidate_file":"end_bridge_candidates.fasta","graph_file":"path_graph.gfa",
        "directed_connections_added":connections.len(),"connections":connections,"searches":searches}),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn lower_k_connects_only_observed_sequence_and_marks_an_unsearched_tail() {
        let mut seed = 617113u64;
        let sequence: Vec<_> = (0..220)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed % 4) as usize]
            })
            .collect();
        let graph = PathGraph {
            segments: vec![
                Segment {
                    sequence: sequence[..90].to_vec(),
                    circular: false,
                },
                Segment {
                    sequence: sequence[130..].to_vec(),
                    circular: false,
                },
            ],
            outgoing: vec![vec![], vec![]],
        };
        let k = 31;
        let mask = (1u128 << (2 * (k - 1))) - 1;
        let mut outgoing = HashMap::<_, Vec<_>>::new();
        kmer::visit(&sequence, None, k, 0, |_, c| {
            outgoing
                .entry(c >> 2)
                .or_default()
                .push((c & mask, b"ACGT"[(c & 3) as usize]));
        });
        let targets = HashMap::from([(code(&sequence[162..192]), vec![1])]);
        let joined = search(0, &graph, &targets, &outgoing, k, 62, 200);
        assert!(!joined.limited);
        assert_eq!(
            joined.bridges,
            BTreeSet::from([(0, 1, sequence[28..192].to_vec())])
        );
        let capped = search(0, &graph, &targets, &outgoing, k, 62, 40);
        assert!(capped.limited);
        assert!(capped.bridges.is_empty());
        outgoing.remove(&code(&sequence[100..130]));
        let absent = search(0, &graph, &targets, &outgoing, k, 62, 200);
        assert!(!absent.limited);
        assert!(absent.bridges.is_empty());
    }
}
