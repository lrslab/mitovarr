//! Compare a bounded family of complete graph cycles on one shared backbone.
//! This is a topology-specific optimization, not a species-specific workflow.
use super::*;

// Keep anchor-mate identity in start clusters: R1 and R2 coordinates are not
// interchangeable molecule starts even though the likelihood is symmetric.
type Start = (usize, isize, bool);
struct Hypothesis {
    local: Route,
    circle: Vec<u8>,
    views: Vec<(usize, usize)>,
}
struct Block {
    intervals: Vec<(usize, usize)>,
    alleles: Vec<usize>,
    allele_count: usize,
}
struct Family {
    hypotheses: Vec<Hypothesis>,
    indexes: Vec<Index>,
    origins: Vec<(usize, bool)>,
    seeds: Vec<Vec<(usize, kmer::Code)>>,
    blocks: Vec<Block>,
    ordered: bool,
    represented: bool,
    complete: bool,
    fragments: usize,
    omitted: usize,
    totals: Vec<f64>,
    offset: f64,
    r1_only: usize,
    r2_only: usize,
    r2_anchored: usize,
    clusters: BTreeMap<Start, Vec<f64>>,
    strong: Vec<Vec<u32>>,
    boundaries: Vec<[BTreeSet<Start>; 2]>,
    block_support: Vec<Vec<BTreeSet<Start>>>,
    phases: Vec<Vec<BTreeSet<Start>>>,
    selected: Option<usize>,
}

pub(super) struct Resolved {
    pub report: Value,
    pub primaries: Vec<primary::Primary>,
}

fn common_seeds(hypotheses: &[Hypothesis]) -> Vec<Vec<(usize, kmer::Code)>> {
    let indexes: Vec<_> = hypotheses
        .iter()
        .map(|h| {
            let mut whole_counts = HashMap::<kmer::Code, usize>::new();
            let circular: Vec<_> = h
                .circle
                .iter()
                .copied()
                .cycle()
                .take(h.circle.len() + 30)
                .collect();
            kmer::visit(&circular, None, 31, 0, |_, code| {
                *whole_counts.entry(kmer::canonical(code, 31)).or_default() += 1;
            });
            let mut index = HashMap::<_, Option<(usize, kmer::Code)>>::new();
            kmer::visit(&h.local.sequence, None, 31, 0, |p, c| {
                index
                    .entry(kmer::canonical(c, 31))
                    .and_modify(|v| *v = None)
                    .or_insert(Some((p, c)));
            });
            index.retain(|code, _| whole_counts.get(code) == Some(&1));
            index
        })
        .collect();
    let mut shared: Vec<_> = indexes[0]
        .iter()
        .filter_map(|(code, first)| {
            let _ = (*first)?;
            indexes
                .iter()
                .map(|m| m.get(code).copied().flatten())
                .collect::<Option<Vec<_>>>()
        })
        .collect();
    shared.sort_unstable_by_key(|positions| positions[0].0);
    shared
}

fn variable_blocks(
    hypotheses: &[Hypothesis],
    seeds: &[Vec<(usize, kmer::Code)>],
) -> (Vec<Block>, bool) {
    let mut regions: Vec<Vec<(usize, usize)>> = Vec::new();
    for adjacent in seeds.windows(2) {
        if (0..hypotheses.len()).any(|i| {
            adjacent[0][i].0 >= adjacent[1][i].0
                || adjacent[0][i].1 != adjacent[0][0].1
                || adjacent[1][i].1 != adjacent[1][0].1
        }) {
            return (Vec::new(), false);
        }
        let first_delta = adjacent[1][0].0 - adjacent[0][0].0;
        let delta_changes =
            (0..hypotheses.len()).any(|i| adjacent[1][i].0 - adjacent[0][i].0 != first_delta);
        let overlaps = (0..hypotheses.len()).any(|i| adjacent[0][i].0 + 31 > adjacent[1][i].0);
        // Two identical shared seeds can overlap by different amounts in a
        // homopolymer indel. Their union, not a nonexistent gap, is the event.
        if overlaps && !delta_changes {
            continue;
        }
        let bounds: Vec<_> = (0..hypotheses.len())
            .map(|i| {
                if overlaps {
                    (adjacent[0][i].0, adjacent[1][i].0 + 31)
                } else {
                    (adjacent[0][i].0 + 31, adjacent[1][i].0)
                }
            })
            .collect();
        let first = &hypotheses[0].local.sequence[bounds[0].0..bounds[0].1];
        if !delta_changes
            && (1..hypotheses.len())
                .all(|i| hypotheses[i].local.sequence[bounds[i].0..bounds[i].1] == *first)
        {
            continue;
        }
        if let Some(previous) = regions.last_mut() {
            if previous.iter().zip(&bounds).any(|(a, b)| a.1 >= b.0) {
                for (a, b) in previous.iter_mut().zip(&bounds) {
                    a.0 = a.0.min(b.0);
                    a.1 = a.1.max(b.1);
                }
                continue;
            }
        }
        regions.push(bounds);
    }
    let blocks = regions
        .into_iter()
        .map(|intervals| {
            let strings: Vec<_> = hypotheses
                .iter()
                .zip(&intervals)
                .map(|(h, &(a, b))| h.local.sequence[a..b].to_vec())
                .collect();
            let unique: BTreeSet<_> = strings.iter().cloned().collect();
            let ids: BTreeMap<_, _> = unique
                .into_iter()
                .enumerate()
                .map(|(i, s)| (s, i))
                .collect();
            Block {
                intervals,
                alleles: strings.iter().map(|s| ids[s]).collect(),
                allele_count: ids.len(),
            }
        })
        .collect();
    (blocks, true)
}

fn blocks_cover_differences(hypotheses: &[Hypothesis], blocks: &[Block]) -> bool {
    let mut cursors = vec![0; hypotheses.len()];
    for block in blocks {
        let first = &hypotheses[0].local.sequence[cursors[0]..block.intervals[0].0];
        if (1..hypotheses.len())
            .any(|i| hypotheses[i].local.sequence[cursors[i]..block.intervals[i].0] != *first)
        {
            return false;
        }
        for (i, cursor) in cursors.iter_mut().enumerate() {
            *cursor = block.intervals[i].1;
        }
    }
    (1..hypotheses.len()).all(|i| {
        hypotheses[i].local.sequence[cursors[i]..] == hypotheses[0].local.sequence[cursors[0]..]
    })
}

fn interval_bridge(
    a: &Placements,
    b: &Placements,
    interval: (usize, usize, Option<usize>),
    model: Option<&Library>,
    baseline: f64,
    pivot: (isize, bool, usize),
    r1_only: bool,
) -> bool {
    let interval = Route {
        circle_length: interval.2,
        to: 0,
        nodes: Vec::new(),
        sequence: Vec::new(),
        left_end: interval.0,
        right_start: interval.1,
        snp_interval: None,
        primary_interval: None,
    };
    score(a, b, &interval, model, baseline, pivot, r1_only).1
}

// Phase calls require one common joint placement. Alternative placements
// may contribute to compatibility scores, but their calls cannot be pooled.
fn phase_pair(
    a: &Placements,
    b: &Placements,
    model: Option<&Library>,
    pivot: (isize, bool, usize),
    r1_only: bool,
    circle_length: Option<usize>,
) -> Option<(super::super::align::Hit, Option<super::super::align::Hit>)> {
    let mut options = Vec::new();
    for x in &a.hits {
        if x.reverse != pivot.1 || (x.start as isize).abs_diff(pivot.0) > pivot.2 {
            continue;
        }
        if r1_only {
            options.push((x.logp, x, None));
            continue;
        }
        for y in &b.hits {
            if circle_length.is_some_and(|n| x.end.max(y.end) - x.start.min(y.start) > n) {
                continue;
            }
            let geometry = if let Some(m) = model {
                let Some(span) = m.span(x, y) else {
                    continue;
                };
                if m.geometry_usable {
                    let p = m.probability(span);
                    if p == 0.0 {
                        continue;
                    }
                    p.ln()
                } else {
                    0.0
                }
            } else {
                0.0
            };
            options.push((x.logp + y.logp + geometry, x, Some(y)));
        }
    }
    options.sort_by(|a, b| b.0.total_cmp(&a.0));
    let best = options.first()?;
    if options.iter().skip(1).any(|p| {
        best.0 - p.0 < 2.3
            && (p.1.start, p.1.end, p.2.map(|h| (h.start, h.end, h.reverse)))
                != (
                    best.1.start,
                    best.1.end,
                    best.2.map(|h| (h.start, h.end, h.reverse)),
                )
    }) {
        return None;
    }
    Some((best.1.clone(), best.2.cloned()))
}
fn sees(
    pair: &[(super::super::align::Hit, super::super::align::BaseTrace)],
    interval: (usize, usize),
) -> bool {
    pair.iter().any(|(h, trace)| {
        h.start + FLANK <= interval.0
            && h.end >= interval.1 + FLANK
            && if interval.0 == interval.1 {
                trace.matching_junctions.contains(&interval.0)
            } else {
                (interval.0..interval.1).all(|p| trace.matching.contains(&p))
            }
    })
}

fn whole_circle_index(h: &Hypothesis, max_read: usize) -> Result<(Index, usize, bool)> {
    let n = h.circle.len();
    for (reverse, sequence) in [h.circle.clone(), kmer::reverse_complement(&h.circle)]
        .into_iter()
        .enumerate()
    {
        let expanded: Vec<_> = sequence
            .iter()
            .copied()
            .cycle()
            .take(n + h.local.sequence.len())
            .collect();
        if let Some(origin) = expanded
            .windows(h.local.sequence.len())
            .position(|w| w == h.local.sequence)
        {
            let rotated: Vec<_> = (0..n).map(|i| sequence[(origin + i) % n]).collect();
            let copies = ((n + h.local.sequence.len()).div_ceil(n) + 1).max(3);
            return Ok((
                Index::new(&rotated.repeat(copies), false, max_read),
                origin % n,
                reverse == 1,
            ));
        }
    }
    anyhow::bail!("local graph walk is inconsistent with its complete cycle");
}

pub(super) fn resolve(
    loci: &mut [Locus],
    graph: &PathGraph,
    reads: &Reads,
    out: &Path,
    settings: &GraphSettings,
    limits: &Limits,
    training: &pairs::Training,
) -> Result<Resolved> {
    let overlap = settings.k - 1;
    let mut grouped = BTreeMap::<Vec<u8>, BTreeMap<Vec<u8>, Hypothesis>>::new();
    for (li, locus) in loci.iter().enumerate() {
        let anchor = &graph.segments[locus.from].sequence;
        let canonical = normalize(anchor, false);
        let reverse = *anchor != canonical;
        for (ri, route) in locus
            .routes
            .iter()
            .enumerate()
            .filter(|(_, r)| r.to == locus.from)
        {
            let mut circle = anchor.clone();
            for &node in &route.nodes {
                circle.extend_from_slice(&graph.segments[node].sequence[overlap..]);
            }
            circle.truncate(circle.len() - overlap);
            let circle = normalize(&circle, true);
            let local = if reverse {
                Route {
                    circle_length: None,
                    to: route.to,
                    nodes: Vec::new(),
                    sequence: kmer::reverse_complement(&route.sequence),
                    left_end: route.sequence.len() - route.right_start,
                    right_start: route.sequence.len() - route.left_end,
                    snp_interval: None,
                    primary_interval: None,
                }
            } else {
                Route {
                    circle_length: None,
                    to: route.to,
                    nodes: Vec::new(),
                    sequence: route.sequence.clone(),
                    left_end: route.left_end,
                    right_start: route.right_start,
                    snp_interval: None,
                    primary_interval: None,
                }
            };
            let family = grouped.entry(canonical.clone()).or_default();
            if let Some(h) = family.get_mut(&circle) {
                h.views.push((li, ri));
            } else {
                family.insert(
                    circle.clone(),
                    Hypothesis {
                        local,
                        circle,
                        views: vec![(li, ri)],
                    },
                );
            }
        }
    }
    let mut families = Vec::new();
    for (_, hypotheses) in grouped {
        let mut hypotheses: Vec<_> = hypotheses.into_values().collect();
        if hypotheses.len() < 2 {
            continue;
        }
        // Once physical cycle alternatives are compared together, a prior
        // one-sided local decision cannot bypass their full-circle evidence.
        for &(li, ri) in hypotheses.iter().flat_map(|h| &h.views) {
            if loci[li].selected == Some(ri) {
                loci[li].selected = None;
            }
        }
        let complete = hypotheses
            .iter()
            .flat_map(|h| &h.views)
            .all(|&(li, _)| !loci[li].exhausted && !loci[li].length_limited);
        let seeds = common_seeds(&hypotheses);
        let (mut blocks, ordered) = variable_blocks(&hypotheses, &seeds);
        let represented = ordered && blocks_cover_differences(&hypotheses, &blocks);
        let n = hypotheses.len();
        let b = blocks.len();
        let indexes = hypotheses
            .iter()
            .map(|h| whole_circle_index(h, training.max_read))
            .collect::<Result<Vec<_>>>()?;
        let (indexes, origins): (Vec<_>, Vec<_>) = indexes
            .into_iter()
            .map(|(index, origin, reverse)| (index, (origin, reverse)))
            .unzip();
        for (i, h) in hypotheses.iter_mut().enumerate() {
            let n = h.circle.len();
            h.local.circle_length = Some(n);
            h.local.left_end += n;
            h.local.right_start += n;
            for block in &mut blocks {
                block.intervals[i].0 += n;
                block.intervals[i].1 += n;
            }
        }
        families.push(Family {
            hypotheses,
            indexes,
            origins,
            seeds,
            blocks,
            ordered,
            represented,
            complete,
            fragments: 0,
            omitted: 0,
            totals: vec![0.0; n],
            offset: 0.0,
            r1_only: 0,
            r2_only: 0,
            r2_anchored: 0,
            clusters: BTreeMap::new(),
            strong: vec![vec![0; n]; n],
            boundaries: (0..n)
                .map(|_| std::array::from_fn(|_| BTreeSet::new()))
                .collect(),
            block_support: vec![vec![BTreeSet::new(); b]; n],
            phases: vec![vec![BTreeSet::new(); b.saturating_sub(1)]; n],
            selected: None,
        });
    }
    let records: Vec<_> = families
        .iter()
        .enumerate()
        .flat_map(|(fi, f)| {
            f.hypotheses.iter().enumerate().map(move |(hi, h)| {
                (
                    format!(
                        "circle_family_{}_path_{} circular=true candidate_only=true",
                        fi + 1,
                        hi + 1
                    ),
                    h.circle.clone(),
                )
            })
        })
        .collect();
    crate::io::fasta::write_fasta_records(&out.join("complete_path_candidates.fasta"), &records)?;
    let mut lookup = HashMap::<_, Vec<(usize, usize)>>::new();
    for (fi, f) in families.iter().enumerate() {
        for (si, s) in f.seeds.iter().enumerate() {
            lookup
                .entry(kmer::canonical(s[0].1, 31))
                .or_default()
                .push((fi, si));
        }
    }
    let mut names = HashSet::new();
    if !families.is_empty() {
        fastq::scan_pairs(&reads.r1, &reads.r2, |_, a, b| {
            let mut relevant = BTreeMap::<usize, [Vec<(usize, usize, kmer::Code)>; 2]>::new();
            for (mate, read) in [a, b].into_iter().enumerate() {
                kmer::visit(
                    &read.sequence,
                    Some(&read.quality),
                    31,
                    limits.min_baseq.min(10),
                    |q, c| {
                        if let Some(hits) = lookup.get(&kmer::canonical(c, 31)) {
                            for &(fi, si) in hits {
                                relevant.entry(fi).or_default()[mate].push((si, q, c));
                            }
                        }
                    },
                );
            }
            if relevant.is_empty() {
                return Ok(());
            }
            let name = a.name(1)?;
            if !names.insert(name.to_owned()) {
                return Ok(());
            }
            ensure!(
                names.len() <= settings.max_validation_pairs,
                "cycle path scoring exceeds --max-validation-pairs"
            );
            let mut pa = Prepared::new(a);
            let mut pb = Prepared::new(b);
            for (fi, anchors) in relevant {
                let mate = usize::from(anchors[0].is_empty());
                let seeds = &anchors[mate];
                let (a, b, pa, pb) = if mate == 0 {
                    (a, b, &mut pa, &mut pb)
                } else {
                    (b, a, &mut pb, &mut pa)
                };
                let ideal_a = super::super::align::ideal_score(a);
                let ideal_b = super::super::align::ideal_score(b);
                let f = &mut families[fi];
                f.fragments += 1;
                f.r2_anchored += mate;
                let tolerance = (a.sequence.len() / 10).max(3);
                let mut pivots = Vec::new();
                for hi in 0..f.hypotheses.len() {
                    let positions: BTreeSet<_> = seeds
                        .iter()
                        .map(|&(si, q, c)| {
                            let (p, forward) = f.seeds[si][hi];
                            let reverse = c != forward;
                            let q = if reverse {
                                a.sequence.len() - 31 - q
                            } else {
                                q
                            };
                            (
                                (p + f.hypotheses[hi].circle.len()) as isize - q as isize,
                                reverse,
                            )
                        })
                        .collect();
                    let &(lo, reverse) = positions.first().unwrap();
                    let &(hi, _) = positions.last().unwrap();
                    if positions.iter().any(|p| p.1 != reverse) || hi.abs_diff(lo) > 2 * tolerance {
                        break;
                    }
                    pivots.push(((lo + hi) / 2, reverse, tolerance));
                }
                if pivots.len() != f.hypotheses.len() {
                    f.omitted += 1;
                    continue;
                }
                let mut placements: Vec<_> = f
                    .indexes
                    .iter()
                    .map(|index| (index.align_prepared(pa), index.align_prepared(pb)))
                    .collect();
                rescue_mates(
                    &f.indexes,
                    &mut placements,
                    pb,
                    &pivots,
                    training.model.as_ref(),
                    a.sequence.len(),
                );
                if placements.iter().any(|(a, b)| a.exhausted || b.exhausted) {
                    f.omitted += 1;
                    continue;
                }
                let r1_only = placements.iter().all(|(_, b)| b.hits.is_empty());
                f.r1_only += usize::from(r1_only && mate == 0);
                f.r2_only += usize::from(r1_only && mate == 1);
                let baseline = ideal_a + if r1_only { 0.0 } else { ideal_b };
                let scores: Vec<_> = placements
                    .iter()
                    .enumerate()
                    .map(|(i, (a, b))| {
                        score(
                            a,
                            b,
                            &f.hypotheses[i].local,
                            training.model.as_ref(),
                            baseline,
                            pivots[i],
                            r1_only,
                        )
                        .0
                    })
                    .collect();
                let maximum = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                f.offset += maximum;
                let relative: Vec<_> = scores.iter().map(|s| s - maximum).collect();
                let start = (mate, pivots[0].0, pivots[0].1);
                let cluster = f
                    .clusters
                    .entry(start)
                    .or_insert_with(|| vec![0.0; scores.len()]);
                let mut observed = vec![vec![false; f.blocks.len()]; scores.len()];
                for i in 0..scores.len() {
                    f.totals[i] += relative[i];
                    cluster[i] += relative[i];
                    let (ref x, ref y) = placements[i];
                    let h = &f.hypotheses[i];
                    for (end, p) in [h.local.left_end + overlap, h.local.right_start]
                        .into_iter()
                        .enumerate()
                    {
                        if interval_bridge(
                            x,
                            y,
                            (p, p, Some(h.circle.len())),
                            training.model.as_ref(),
                            baseline,
                            pivots[i],
                            r1_only,
                        ) {
                            f.boundaries[i][end].insert(start);
                        }
                    }
                    for (bi, block) in f.blocks.iter().enumerate() {
                        let (left, right) = block.intervals[i];
                        observed[i][bi] = interval_bridge(
                            x,
                            y,
                            (left, right, Some(h.circle.len())),
                            training.model.as_ref(),
                            baseline,
                            pivots[i],
                            r1_only,
                        );
                        if observed[i][bi] {
                            f.block_support[i][bi].insert(start);
                        }
                    }
                    // Strong sequence/geometry contrasts are diagnostic; the
                    // separate block and phase conditions supply bridging support.
                    for j in 0..scores.len() {
                        if relative[i] - relative[j] >= 2.3 {
                            f.strong[i][j] += 1;
                        }
                    }
                }
                let phase_pairs: Vec<_> = placements
                    .iter()
                    .enumerate()
                    .map(|(i, (placed_a, placed_b))| {
                        let pair = phase_pair(
                            placed_a,
                            placed_b,
                            training.model.as_ref(),
                            pivots[i],
                            r1_only,
                            Some(f.hypotheses[i].circle.len()),
                        );
                        pair.into_iter()
                            .flat_map(|(x, y)| std::iter::once((x, a)).chain(y.map(|h| (h, b))))
                            .map(|(hit, read)| {
                                let trace = super::super::align::trace_alignment(
                                    read,
                                    &hit,
                                    &f.indexes[i].sequence,
                                    limits.min_baseq,
                                );
                                (hit, trace)
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect();
                let mut calls = vec![None; f.blocks.len()];
                for (bi, block) in f.blocks.iter().enumerate() {
                    let mut allele_scores = vec![f64::NEG_INFINITY; block.allele_count];
                    for (i, &score) in relative.iter().enumerate() {
                        allele_scores[block.alleles[i]] =
                            allele_scores[block.alleles[i]].max(score);
                    }
                    let best = (0..allele_scores.len())
                        .max_by(|&a, &b| {
                            allele_scores[a]
                                .total_cmp(&allele_scores[b])
                                .then(b.cmp(&a))
                        })
                        .unwrap();
                    if (0..allele_scores.len())
                        .filter(|&i| i != best)
                        .all(|i| allele_scores[best] - allele_scores[i] >= 2.3)
                    {
                        // Require an actual read across the variable bases.
                        // An empty mate gap is not an allele observation.
                        if (0..scores.len()).any(|i| {
                            block.alleles[i] == best && sees(&phase_pairs[i], block.intervals[i])
                        }) {
                            calls[bi] = Some(best);
                        }
                    }
                }
                for (i, phase_pair) in phase_pairs.iter().enumerate() {
                    for bi in 0..f.blocks.len().saturating_sub(1) {
                        if calls[bi] == Some(f.blocks[bi].alleles[i])
                            && calls[bi + 1] == Some(f.blocks[bi + 1].alleles[i])
                            && sees(phase_pair, f.blocks[bi].intervals[i])
                            && sees(phase_pair, f.blocks[bi + 1].intervals[i])
                        {
                            f.phases[i][bi].insert(start);
                        }
                    }
                }
            }
            Ok(())
        })?;
    }
    let mut reports = Vec::new();
    let mut primaries = Vec::new();
    for (fi, f) in families.iter_mut().enumerate() {
        let best = (0..f.totals.len())
            .max_by(|&a, &b| f.totals[a].total_cmp(&f.totals[b]).then(b.cmp(&a)))
            .unwrap();
        let margin = (0..f.totals.len()).filter(|&i| i != best).all(|i| {
            let difference = f.totals[best] - f.totals[i];
            let square: f64 = f.clusters.values().map(|s| (s[best] - s[i]).powi(2)).sum();
            stable_margin(difference, square, f.omitted)
                && minor_supported(f.strong[best][i], f.strong[i][best], f.omitted, settings)
        });
        let connected = f.boundaries[best]
            .iter()
            .chain(&f.block_support[best])
            .chain(&f.phases[best])
            .all(|s| s.len() >= settings.min_junction_fragments as usize);
        // Other outgoing topologies retain the original shared source-pivot
        // comparison. A well-scored internal bubble cannot erase a competing exit.
        let exits = f.hypotheses[best].views.iter().all(|&(li, ri)| {
            let l = &loci[li];
            l.routes
                .iter()
                .enumerate()
                .filter(|(_, r)| r.to != l.from)
                .all(|(other, _)| {
                    let d = l.totals[ri] - l.totals[other];
                    stable_margin(d, l.delta_squares[ri][other], l.alignment_limited)
                        && minor_supported(
                            l.strong_votes[ri][other],
                            l.strong_votes[other][ri],
                            l.alignment_limited,
                            settings,
                        )
                })
        });
        if f.complete
            && f.ordered
            && f.represented
            && margin
            && connected
            && exits
            && !f.blocks.is_empty()
        {
            f.selected = Some(best);
            for &(li, ri) in &f.hypotheses[best].views {
                loci[li].selected = Some(ri);
            }
        }
        // A primary representative uses the score of ONE existing complete
        // graph walk. Missing long-range phase or a close runner-up remains
        // localized uncertainty; it need not fragment the shared backbone.
        let primary_path = (f.complete
            && f.ordered
            && f.represented
            && exits
            && f.fragments > f.omitted
            && !f.blocks.is_empty()
            && f.boundaries[best]
                .iter()
                .all(|s| s.len() >= settings.min_junction_fragments as usize))
        .then_some(best);
        let runner_up = (0..f.totals.len())
            .filter(|&i| i != best)
            .max_by(|&a, &b| f.totals[a].total_cmp(&f.totals[b]).then(b.cmp(&a)))
            .unwrap();
        if primary_path.is_some() {
            let h = &f.hypotheses[best];
            let (origin, reverse) = f.origins[best];
            let reason = if f.selected.is_some() {
                "minor_cycle_path_retained"
            } else {
                "competing_cycle_paths"
            };
            let regions = f
                .blocks
                .iter()
                .enumerate()
                .flat_map(|(bi, b)| {
                    primary::project_interval(
                        b.intervals[best].0,
                        b.intervals[best].1,
                        h.circle.len(),
                        origin,
                        reverse,
                    )
                    .into_iter()
                    .map(move |(start, end)| primary::Region {
                        start,
                        end,
                        evidence: format!("circle_family_{}_block_{}", fi + 1, bi + 1),
                        reason,
                        alternatives: b.allele_count - 1,
                    })
                })
                .collect();
            primaries.push(primary::Primary {
                candidate: Candidate { sequence: h.circle.clone(), circular: true },
                source: json!({"type":"complete_cycle_family","family":fi+1,"path":best+1,
                    "candidate_id":format!("circle_family_{}_path_{}",fi+1,best+1),
                    "candidate_file":"complete_path_candidates.fasta","score":f.totals[best],
                    "runner_up_path":runner_up+1,"runner_up_score_margin":f.totals[best]-f.totals[runner_up],
                    "strict_selected":f.selected==Some(best)}),
                regions,
            });
        }
        reports.push(json!({"family":fi+1,"search_complete":f.complete,"common_anchor_order_consistent":f.ordered,"all_candidate_differences_represented_by_blocks":f.represented,"shared_unique_31mers":f.seeds.len(),"variable_blocks":f.blocks.len(),"fragments":f.fragments,"r1_only_fragments":f.r1_only,"r2_only_fragments":f.r2_only,"r2_anchored_fragments":f.r2_anchored,"candidate_independent_log_offset":f.offset,"score_definition":"sum of fragment log scores minus each fragment maximum over every candidate","omitted_alignment_or_pivot_limit":f.omitted,"unscored_maximum_total_contrast":unscored_contrast_bound(f.omitted),"distinct_r1_start_clusters":f.clusters.keys().filter(|s|s.0==0).count(),"distinct_anchor_start_clusters":f.clusters.len(),"selected_path":f.selected.map(|i|i+1),"best_path":best+1,"primary_path":primary_path.map(|i|i+1),"runner_up_path":runner_up+1,"runner_up_score_margin":f.totals[best]-f.totals[runner_up],"inference_status":if f.selected.is_some(){"dominant_consensus_with_observed_adjacent_phase"}else if !f.complete{"search_incomplete"}else{"mixed_or_unphased"},"unique_biological_haplotype_claim":false,"score_margin_passed":margin,"block_phase_and_boundary_support_passed":connected,"competing_exit_comparison_passed":exits,"paths":f.hypotheses.iter().enumerate().map(|(i,h)|json!({"path":i+1,"length":h.circle.len(),"alignment_origin_in_oriented_fasta":f.origins[i].0,"alignment_reverse_complement":f.origins[i].1,"score":f.totals[i],"source_views":h.views.iter().map(|&(li,ri)|json!({"locus":li+1,"path":ri+1})).collect::<Vec<_>>(),"boundary_start_clusters":f.boundaries[i].iter().map(BTreeSet::len).collect::<Vec<_>>(),"block_start_clusters":f.block_support[i].iter().map(BTreeSet::len).collect::<Vec<_>>(),"adjacent_block_phase_start_clusters":f.phases[i].iter().map(BTreeSet::len).collect::<Vec<_>>(),"block_alleles":f.blocks.iter().map(|b|b.alleles[i]).collect::<Vec<_>>(),"block_intervals":f.blocks.iter().map(|b|b.intervals[i]).collect::<Vec<_>>(),"strong_votes_against_alternatives":f.strong[i]})).collect::<Vec<_>>()}));
    }
    Ok(Resolved {
        primaries,
        report: json!({"method":"physical complete cycles on a shared sample backbone; reverse-complement views collapsed; common unique pivots from either mate throughout both flanks and internal corridor; one fragment once per family","candidate_file":"complete_path_candidates.fasta","selection_scope":"primary consensus conditional on enumerated graph paths; not a unique biological haplotype claim","phase_evidence":"adjacent variable blocks called jointly by one fragment, with reads spanning variable bases; at least the existing junction threshold in distinct anchor-mate start clusters","alignment_scope":"complete circular hypotheses, explicit unwrapped coordinates and at most one circumference per joint fragment","block_coordinate_system":"zero-based half-open coordinates in the unwrapped candidate; pivot origin shifted by one circumference","family_count":families.len(),"families":reports}),
    })
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
    #[test]
    fn phase_calls_require_observed_alleles_and_a_contiguous_empty_junction() {
        use super::super::super::align::trace_alignment;
        let mut target = dna(100, 991817);
        target[49..52].copy_from_slice(b"ACG");
        for change in [0, 1, 2] {
            let mut sequence = target.clone();
            if change == 1 {
                sequence.remove(50);
            }
            if change == 2 {
                sequence.insert(50, b'T');
            }
            let read = fastq::Record {
                header: "@phase/1".into(),
                quality: vec![b'I'; sequence.len()],
                sequence,
                plus: "+".into(),
            };
            let index = Index::new(&target, false, read.sequence.len());
            let placements = index.align(&read);
            let hit = placements
                .hits
                .iter()
                .find(|h| h.start == 0 && h.end == target.len() && !h.reverse)
                .unwrap();
            let trace = trace_alignment(&read, hit, &target, 20);
            let observations = vec![(hit.clone(), trace)];
            assert_eq!(
                sees(&observations, (50, 51)),
                change != 1,
                "deleted target base is not an observed allele"
            );
            assert_eq!(
                sees(&observations, (50, 50)),
                change == 0,
                "insertion or deletion cannot certify an empty junction"
            );
        }
    }
    #[test]
    fn common_anchors_separate_linked_variants_without_using_allele_reads_as_pivots() {
        let base = dna(600, 18713);
        let hypotheses: Vec<_> = (0..4)
            .map(|i| {
                let mut s = base.clone();
                if i & 1 != 0 {
                    s[200] = if s[200] == b'A' { b'C' } else { b'A' };
                }
                if i & 2 != 0 {
                    s.insert(350, b'G');
                }
                Hypothesis {
                    local: Route {
                        circle_length: None,
                        to: 0,
                        nodes: vec![],
                        sequence: s.clone(),
                        left_end: 100,
                        right_start: 500,
                        snp_interval: None,
                        primary_interval: None,
                    },
                    circle: s,
                    views: vec![],
                }
            })
            .collect();
        let shared = common_seeds(&hypotheses);
        let (blocks, ordered) = variable_blocks(&hypotheses, &shared);
        assert!(ordered);
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].intervals[0], (200, 201));
        assert_eq!(blocks[1].intervals[0], (350, 350));
        assert_eq!(blocks[1].intervals[2], (350, 351));
        assert_ne!(blocks[0].alleles[0], blocks[0].alleles[1]);
        assert_eq!(blocks[0].alleles[0], blocks[0].alleles[2]);
    }
    #[test]
    fn internal_paired_evidence_phases_a_corridor_but_unlinked_marginals_do_not() {
        use super::super::super::threading::Segment;
        use std::io::Write;
        let anchor = dna(1000, 821817);
        let core = dna(600, 121812);
        let true_circle = [anchor.clone(), core.clone()].concat();
        let mut segments = vec![Segment {
            sequence: anchor.clone(),
            circular: false,
        }];
        for bits in 0..4 {
            let mut inner = core.clone();
            for (mask, p) in [(1, 100), (2, 300)] {
                if bits & mask != 0 {
                    inner[p] = if inner[p] == b'A' { b'C' } else { b'A' };
                }
            }
            segments.push(Segment {
                sequence: [
                    anchor[anchor.len() - 62..].to_vec(),
                    inner,
                    anchor[..62].to_vec(),
                ]
                .concat(),
                circular: false,
            });
        }
        let mut outgoing = vec![vec![1, 2, 3, 4], vec![0], vec![0], vec![0], vec![0]];
        let reverse_segments: Vec<_> = segments
            .iter()
            .map(|s| Segment {
                sequence: kmer::reverse_complement(&s.sequence),
                circular: false,
            })
            .collect();
        segments.extend(reverse_segments);
        outgoing.extend([vec![6, 7, 8, 9], vec![5], vec![5], vec![5], vec![5]]);
        let graph = PathGraph { segments, outgoing };
        let settings = GraphSettings {
            coverage_layout: false,
            multi_k_layout: false,
            k: 63,
            min_kmer_count: 3,
            max_error_fraction: 0.2,
            max_bubble_fraction: 0.49,
            min_depth: 10,
            min_junction_fragments: 3,
            max_validation_pairs: 10000,
        };
        let limits = Limits {
            max_kmers: 10000,
            min_baseq: 20,
        };
        for phased in [true, false] {
            let out = std::env::temp_dir().join(format!(
                "mitovar-cycle-phasing-{}-{phased}",
                std::process::id()
            ));
            std::fs::create_dir_all(&out).unwrap();
            let reads = Reads {
                r1: out.join("r1.fq"),
                r2: out.join("r2.fq"),
            };
            let mut r1 = std::fs::File::create(&reads.r1).unwrap();
            let mut r2 = std::fs::File::create(&reads.r2).unwrap();
            let mut count = 0;
            for start in 0..true_circle.len() {
                // Remove only joint observations, retaining both marginal allele
                // coverages. Their high scores must not manufacture a haplotype.
                if !phased && (1001..=1100).contains(&start) {
                    continue;
                }
                let a: Vec<_> = (0..100)
                    .map(|i| true_circle[(start + i) % true_circle.len()])
                    .collect();
                let b: Vec<_> = (0..100)
                    .map(|i| true_circle[(start + 200 + i) % true_circle.len()])
                    .collect();
                writeln!(
                    r1,
                    "@f{start}/1\n{}\n+\n{}",
                    String::from_utf8_lossy(&a),
                    "I".repeat(100)
                )
                .unwrap();
                writeln!(
                    r2,
                    "@f{start}/2\n{}\n+\n{}",
                    String::from_utf8_lossy(&kmer::reverse_complement(&b)),
                    "I".repeat(100)
                )
                .unwrap();
                count += 1;
            }
            // These fragments have a shared backbone anchor exclusively on
            // R2. Missing R1 must not exclude them or duplicate their votes.
            for start in 0..5 {
                writeln!(
                    r1,
                    "@r2_anchor{start}/1\n{}\n+\n{}",
                    "N".repeat(100),
                    "I".repeat(100)
                )
                .unwrap();
                writeln!(
                    r2,
                    "@r2_anchor{start}/2\n{}\n+\n{}",
                    String::from_utf8_lossy(&kmer::reverse_complement(
                        &true_circle[200 + start..300 + start]
                    )),
                    "I".repeat(100)
                )
                .unwrap();
                count += 1;
            }
            drop((r1, r2));
            let training = pairs::estimate(&graph, &reads, &settings, &limits).unwrap();
            let anchors = HashSet::from([0, 5]);
            let mut loci: Vec<_> = [0, 5]
                .into_iter()
                .map(|from| enumerate(&graph, from, &anchors, 62, 1000, 300, 100))
                .collect();
            for locus in &mut loci {
                locus.selected = Some(0);
            }
            let result = resolve(
                &mut loci, &graph, &reads, &out, &settings, &limits, &training,
            )
            .unwrap();
            assert_eq!(result.report["family_count"], 1);
            assert_eq!(result.primaries.len(), 1);
            assert_eq!(
                result.primaries[0].candidate.sequence,
                normalize(&true_circle, true)
            );
            assert_eq!(result.primaries[0].regions.len(), 2);
            let family = &result.report["families"][0];
            assert!(family["primary_path"].is_number());
            assert_eq!(family["paths"].as_array().unwrap().len(), 4);
            assert_eq!(family["variable_blocks"], 2);
            assert!(family["fragments"].as_u64().unwrap() <= count);
            assert!(family["r2_anchored_fragments"].as_u64().unwrap() >= 5);
            assert!(family["r2_only_fragments"].as_u64().unwrap() >= 5);
            assert_eq!(family["selected_path"].is_number(), phased, "{family}");
            if phased {
                let selected = family["selected_path"].as_u64().unwrap() as usize - 1;
                assert_eq!(
                    family["paths"][selected]["source_views"]
                        .as_array()
                        .unwrap()
                        .len(),
                    2
                );
                assert!(loci.iter().all(|l| l.selected == Some(0)));
            } else {
                assert_eq!(family["block_phase_and_boundary_support_passed"], false);
                assert!(loci.iter().all(|l| l.selected.is_none()));
            }
            std::fs::remove_dir_all(&out).unwrap();
        }
    }
    #[test]
    fn distinct_joint_placements_cannot_be_combined_into_phase_evidence() {
        let hit = |start, reverse, logp| super::super::super::align::Hit {
            start,
            end: start + 100,
            reverse,
            logp,
            edits: 0,
            traceback: None,
        };
        let a = Placements {
            hits: vec![hit(100, false, 0.0)],
            exhausted: false,
            complete: false,
        };
        let mut b = Placements {
            hits: vec![hit(300, true, 0.0), hit(500, true, 0.0)],
            exhausted: false,
            complete: false,
        };
        assert!(phase_pair(&a, &b, None, (100, false, 3), false, None).is_none());
        b.hits[1].logp = -8.0;
        let pair = phase_pair(&a, &b, None, (100, false, 3), false, None).unwrap();
        assert_eq!(pair.1.unwrap().start, 300);
    }
    #[test]
    fn overlapping_shared_anchors_represent_homopolymer_length_changes() {
        let left = dna(120, 749171);
        let right = dna(120, 723181);
        let hypotheses: Vec<_> = [15, 16, 18]
            .into_iter()
            .map(|copies| {
                let sequence = [left.clone(), vec![b'A'; copies], right.clone()].concat();
                Hypothesis {
                    local: Route {
                        circle_length: None,
                        to: 0,
                        nodes: vec![],
                        left_end: 100,
                        right_start: 150,
                        sequence: sequence.clone(),
                        snp_interval: None,
                        primary_interval: None,
                    },
                    circle: sequence,
                    views: vec![],
                }
            })
            .collect();
        let seeds = common_seeds(&hypotheses);
        let (blocks, ordered) = variable_blocks(&hypotheses, &seeds);
        assert!(ordered);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].allele_count, 3);
        assert!(blocks_cover_differences(&hypotheses, &blocks));
        assert!(!blocks_cover_differences(&hypotheses, &[]));
        assert!(blocks[0].intervals.iter().all(|(a, b)| b - a > 31));
    }

    #[test]
    fn complete_circle_anchor_uniqueness_and_unwrapped_coordinates() {
        let local = dna(260, 712899);
        let repeated = local[90..130].to_vec();
        let physical = [local.clone(), dna(180, 928118), repeated].concat();
        let make = |circle| Hypothesis {
            local: Route {
                circle_length: None,
                to: 0,
                nodes: vec![],
                left_end: 50,
                right_start: 200,
                sequence: local.clone(),
                snp_interval: None,
                primary_interval: None,
            },
            circle,
            views: vec![],
        };
        let h = make(normalize(&physical, true));
        let shared = common_seeds(&[h]);
        assert!(shared.iter().all(|s| !(90..=99).contains(&s[0].0)));
        assert!(shared.iter().any(|s| s[0].0 == 50));
        for circle in [
            normalize(&physical, true),
            kmer::reverse_complement(&physical),
        ] {
            let h = make(circle);
            let (index, _, _) = whole_circle_index(&h, 100).unwrap();
            let n = physical.len();
            assert_eq!(&index.sequence[n..n + local.len()], &local);
            assert_eq!(&index.sequence[n - 40..n], &physical[n - 40..]);
            assert_eq!(
                &index.sequence[n + local.len()..2 * n],
                &physical[local.len()..]
            );
        }
    }
}
