//! Seeded, full-query affine-gap alignment. Reference flanks are free; query ends are not.
use crate::{io::fastq::Record, kmer};
use std::collections::{BTreeSet, HashMap};
use std::sync::{Arc, LazyLock};

mod blocks;
pub(crate) mod profile;

#[cfg(test)]
mod fitting_tests;

const SEED: usize = 17;
const BAND: usize = 12;
const MAX_STARTS: usize = 256;
const NEG: f64 = -1e100;
const GAP_OPEN: f64 = -8.0;
const GAP_EXTEND: f64 = -2.0;

pub(crate) fn edit_budget(length: usize) -> usize {
    (length / 10).max(3)
}

#[derive(Clone, Debug)]
pub(crate) struct Hit {
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
    pub logp: f64,
    pub edits: usize,
    pub traceback: Option<Arc<profile::Traceback>>,
}

#[derive(Clone, Default, Debug)]
pub(crate) struct Placements {
    pub hits: Vec<Hit>,
    pub exhausted: bool,
    /// Complete only in the declared endpoint/edit domain, not across all
    /// possible alignments or ungenerated assembly structures.
    pub complete: bool,
}

struct OrientedQuery {
    sequence: Vec<u8>,
    quality: Vec<u8>,
    seeds: Vec<(usize, kmer::Code)>,
    emissions: Vec<[f64; 2]>,
    blocks: Option<blocks::ExactBlocks>,
    ideal: f64,
}

/// A cache belongs to one read, including its qualities. It lives for one
/// fragment comparison and is reused across candidates, never across reads.
pub(crate) struct Prepared {
    orientations: [OrientedQuery; 2],
    windows: [HashMap<Vec<u8>, Vec<Hit>>; 2],
    rescues: [HashMap<Vec<u8>, Vec<Hit>>; 2],
    profiles: [HashMap<profile::Key, Placements>; 2],
    subproblems: [profile::Subproblems; 2],
    length: usize,
}
impl Prepared {
    pub fn new(read: &Record) -> Self {
        let orientations = [false, true].map(|reverse| {
            let sequence = if reverse {
                kmer::reverse_complement(&read.sequence)
            } else {
                read.sequence.iter().map(u8::to_ascii_uppercase).collect()
            };
            let quality = if reverse {
                read.quality.iter().rev().copied().collect()
            } else {
                read.quality.clone()
            };
            let mut seeds = Vec::new();
            // Seeds propose placements; quality is evidence in emission scores,
            // not a hidden hard cutoff that can erase every seed in a read.
            kmer::visit(&sequence, None, SEED, 0, |p, c| seeds.push((p, c)));
            let emissions: Vec<_> = sequence
                .iter()
                .zip(&quality)
                .map(|(&b, &q)| {
                    static SCORES: LazyLock<[[f64; 2]; 256]> = LazyLock::new(|| {
                        std::array::from_fn(|q| {
                            [emission(b'A', b'A', q as u8), emission(b'A', b'C', q as u8)]
                        })
                    });
                    if b"ACGT".contains(&b) {
                        SCORES[usize::from(q)]
                    } else {
                        [0.25f64.ln(); 2]
                    }
                })
                .collect();
            let ideal = emissions.iter().map(|scores| scores[0]).sum();
            let blocks = blocks::ExactBlocks::new(&sequence, edit_budget(sequence.len()));
            OrientedQuery {
                sequence,
                quality,
                seeds,
                emissions,
                blocks,
                ideal,
            }
        });
        Self {
            orientations,
            windows: std::array::from_fn(|_| HashMap::new()),
            rescues: std::array::from_fn(|_| HashMap::new()),
            profiles: std::array::from_fn(|_| HashMap::new()),
            subproblems: std::array::from_fn(|_| profile::Subproblems::default()),
            length: read.sequence.len(),
        }
    }
}

pub(crate) struct Index {
    pub sequence: Vec<u8>,
    seeds: HashMap<kmer::Code, Vec<usize>>,
    circle_length: Option<usize>,
}

fn emission(a: u8, b: u8, quality: u8) -> f64 {
    if !b"ACGT".contains(&a) || !b"ACGT".contains(&b) {
        return 0.25f64.ln();
    }
    // Cap the nominal accuracy: a base quality is not a guarantee of zero errors.
    let error = 10f64
        .powf(-f64::from(quality.saturating_sub(33)) / 10.0)
        .clamp(0.001, 0.75);
    if a == b {
        (1.0 - error).ln()
    } else {
        (error / 3.0).ln()
    }
}

pub(crate) fn ideal_score(read: &Record) -> f64 {
    read.sequence
        .iter()
        .zip(&read.quality)
        .map(|(&b, &q)| emission(b.to_ascii_uppercase(), b.to_ascii_uppercase(), q))
        .sum()
}

impl Index {
    pub(crate) fn circle_length(&self) -> Option<usize> {
        self.circle_length
    }

    pub fn new(sequence: &[u8], circular: bool, max_read: usize) -> Self {
        let mut expanded = sequence.to_vec();
        if circular && !sequence.is_empty() {
            expanded.extend(
                (0..max_read + BAND.max(edit_budget(max_read)))
                    .map(|i| sequence[i % sequence.len()]),
            );
        }
        let mut seeds: HashMap<_, Vec<_>> = HashMap::new();
        kmer::visit(&expanded, None, SEED, 0, |p, code| {
            seeds.entry(code).or_default().push(p);
        });
        Self {
            sequence: expanded,
            seeds,
            circle_length: circular.then_some(sequence.len()),
        }
    }

    pub fn align(&self, read: &Record) -> Placements {
        self.align_prepared(&mut Prepared::new(read))
    }

    /// Bounded seedless rescue in the region permitted by an independently
    /// anchored mate. Used symmetrically before an absent candidate hit can
    /// become negative evidence. It never searches an entire genome by default.
    pub fn rescue_near(
        &self,
        prepared: &mut Prepared,
        left: usize,
        right: usize,
        reverse: bool,
    ) -> Placements {
        let mut result = Placements::default();
        let left = left.min(self.sequence.len());
        let right = right.min(self.sequence.len());
        if left >= right {
            return result;
        }
        let orientation = usize::from(reverse);
        let query = &prepared.orientations[orientation];
        let target = &self.sequence[left..right];
        let cache = &mut prepared.rescues[orientation];
        let hits = cache
            .entry(target.to_vec())
            .or_insert_with(|| fitting_candidates(query, target));
        for hit in hits {
            if hit.edits <= (prepared.length / 10).max(3) {
                let mut hit = hit.clone();
                hit.start += left;
                hit.end += left;
                hit.reverse = reverse;
                result.hits.push(hit);
            }
        }
        if result.hits.len() > MAX_STARTS {
            result.hits.clear();
            result.exhausted = true;
        }
        result
    }

    pub fn align_prepared(&self, prepared: &mut Prepared) -> Placements {
        let mut result = Placements::default();
        for (orientation, query) in prepared.orientations.iter().enumerate() {
            let reverse = orientation == 1;
            let sequence = &query.sequence;
            let mut starts = BTreeSet::new();
            for &(q, code) in &query.seeds {
                if let Some(positions) = self.seeds.get(&code) {
                    for &p in positions {
                        let start = p as isize - q as isize;
                        if start >= -(BAND as isize) && starts.len() <= MAX_STARTS {
                            starts.insert(start);
                        }
                    }
                }
            }
            if starts.len() > MAX_STARTS {
                result.exhausted = true;
                continue;
            }
            let mut exact = Vec::new();
            for &start in &starts {
                if start >= 0 {
                    let start = start as usize;
                    if self.sequence.get(start..start + sequence.len()) == Some(sequence.as_slice())
                    {
                        exact.push(Hit {
                            start,
                            end: start + sequence.len(),
                            reverse,
                            logp: query.ideal,
                            edits: 0,
                            traceback: None,
                        });
                    }
                }
            }
            // An exact hit is already the optimal full-query score. Other genomic
            // locations are retained below if they have independent seed clusters.
            let exact_starts: Vec<_> = exact.iter().map(|h| h.start as isize).collect();
            result.hits.extend(exact);
            let mut windows = Vec::<(usize, usize)>::new();
            for start in starts {
                if exact_starts.iter().any(|p| p.abs_diff(start) <= BAND) {
                    continue;
                }
                let left = start.saturating_sub(BAND as isize).max(0) as usize;
                let right = ((start + sequence.len() as isize + BAND as isize).max(0) as usize)
                    .min(self.sequence.len());
                if let Some(previous) = windows.last_mut() {
                    if left <= previous.1.saturating_sub(sequence.len()) {
                        previous.1 = previous.1.max(right);
                        continue;
                    }
                }
                windows.push((left, right));
            }
            for (left, right) in windows {
                let target = &self.sequence[left..right];
                let cache = &mut prepared.windows[orientation];
                let hits = cache
                    .entry(target.to_vec())
                    .or_insert_with(|| fitting_candidates(query, target));
                for hit in hits {
                    let mut hit = hit.clone();
                    hit.start += left;
                    hit.end += left;
                    hit.reverse = reverse;
                    result.hits.push(hit);
                }
            }
        }
        if let Some(n) = self.circle_length {
            result.hits.retain(|h| h.start < n && h.end - h.start <= n);
        }
        result.hits.sort_by(|a, b| {
            a.start
                .cmp(&b.start)
                .then(a.reverse.cmp(&b.reverse))
                .then(b.logp.total_cmp(&a.logp))
        });
        result
            .hits
            .dedup_by(|b, a| a.start == b.start && a.reverse == b.reverse);
        let best = result.hits.iter().map(|h| h.logp).fold(NEG, f64::max);
        result
            .hits
            .retain(|h| h.edits <= (prepared.length / 10).max(3) && h.logp >= best - 16.0);
        result
    }
}

impl Index {
    fn profile_window(
        &self,
        prepared: &mut Prepared,
        left: usize,
        right: usize,
        reverse: bool,
        constraint: Option<profile::Constraint>,
    ) -> Placements {
        let left = left.min(self.sequence.len());
        let right = right.min(self.sequence.len());
        if left >= right {
            return Placements {
                complete: true,
                ..Placements::default()
            };
        }
        let orientation = usize::from(reverse);
        let query = &prepared.orientations[orientation];
        let key = profile::Key {
            target: self.sequence[left..right].to_vec(),
            constraint: constraint.map(|c| profile::Constraint {
                target: c.target - left,
                ..c
            }),
            edits: edit_budget(prepared.length),
        };
        let cache = &mut prepared.profiles[orientation];
        let mut found = if let Some(found) = cache.get(&key) {
            found.clone()
        } else {
            let found = profile::fit_cached(
                &query.sequence,
                &query.quality,
                &key.target,
                key.constraint,
                key.edits,
                profile::MAX_CELLS,
                &mut prepared.subproblems[orientation],
            );
            // Per-fragment caches remain bounded when many hypotheses differ.
            if cache.len() < 32 {
                cache.insert(key, found.clone());
            }
            found
        };
        for hit in &mut found.hits {
            hit.start += left;
            hit.end += left;
            hit.reverse = reverse;
        }
        found
    }

    fn normalize_profiles(&self, found: &mut Placements) {
        if let Some(n) = self.circle_length {
            found.hits.retain(|h| h.end - h.start <= n);
            for hit in &mut found.hits {
                let shift = hit.start / n * n;
                hit.start -= shift;
                hit.end -= shift;
            }
        }
        profile::dedup(found);
    }

    /// Enumerate all endpoint/edit interpretations of the fixed physical word.
    /// The chosen query offset and target identity are enforced inside the DP.
    pub(crate) fn profile_anchor(
        &self,
        prepared: &mut Prepared,
        target: usize,
        reverse: bool,
        query: usize,
        length: usize,
    ) -> Placements {
        let mut result = Placements {
            complete: true,
            ..Placements::default()
        };
        let edits = edit_budget(prepared.length);
        let n = self.circle_length.unwrap_or(0);
        let copies: BTreeSet<_> = [target, target + n].into_iter().collect();
        for target in copies {
            if target + length > self.sequence.len() {
                continue;
            }
            let left = (target as isize - query as isize - edits as isize).max(0) as usize;
            let right = target + prepared.length - query + edits;
            if n > 0 && left >= n {
                continue;
            }
            let found = self.profile_window(
                prepared,
                left,
                right,
                reverse,
                Some(profile::Constraint {
                    query,
                    target,
                    length,
                }),
            );
            result.complete &= found.complete;
            result.exhausted |= found.exhausted;
            result.hits.extend(found.hits);
        }
        self.normalize_profiles(&mut result);
        result
    }

    pub(crate) fn profile_around(
        &self,
        prepared: &mut Prepared,
        pivot: &Hit,
        radius: usize,
        reverse: bool,
    ) -> Placements {
        let mut result = Placements {
            complete: true,
            ..Placements::default()
        };
        let n = self.circle_length.unwrap_or(0) as isize;
        let mut windows = BTreeSet::new();
        for shift in [-n, 0, n] {
            let lo = (pivot.start as isize - radius as isize + shift).max(0) as usize;
            let hi = (pivot.end as isize + radius as isize + shift).max(0) as usize;
            if lo < self.sequence.len() && hi > lo {
                windows.insert((lo, hi.min(self.sequence.len())));
            }
        }
        for (lo, hi) in windows {
            let found = self.profile_window(prepared, lo, hi, reverse, None);
            result.complete &= found.complete;
            result.exhausted |= found.exhausted;
            result.hits.extend(found.hits);
        }
        self.normalize_profiles(&mut result);
        result
    }

    /// Without usable geometry, sequence comparison must either search both
    /// strands of the whole hypothesis in the edit domain or remain unknown.
    pub(crate) fn profile_all(&self, prepared: &mut Prepared) -> Placements {
        let mut result = Placements {
            complete: true,
            ..Placements::default()
        };
        for reverse in [false, true] {
            let found = self.profile_window(prepared, 0, self.sequence.len(), reverse, None);
            result.complete &= found.complete;
            result.exhausted |= found.exhausted;
            result.hits.extend(found.hits);
        }
        self.normalize_profiles(&mut result);
        result
    }
}

#[derive(Clone, Copy)]
struct Cell {
    score: f64,
    start: usize,
    edits: usize,
}
impl Cell {
    fn advance(self, score: f64, edit: usize) -> Self {
        Self {
            score: self.score + score,
            start: self.start,
            edits: self.edits + edit,
        }
    }
}
fn best(a: Cell, b: Cell) -> Cell {
    if a.score > b.score || (a.score == b.score && (a.edits, a.start) <= (b.edits, b.start)) {
        a
    } else {
        b
    }
}

/// M/I/D recurrences with full-query ends. The edit bound limits the active
/// diagonals: an accepted path cannot consume >max_edits extra query or target
/// bases. This excludes only paths that the caller would reject anyway.
#[cfg(test)]
fn fitting(query: &[u8], quality: &[u8], target: &[u8]) -> Option<Hit> {
    let read = Record {
        header: "@fitting".into(),
        sequence: query.to_vec(),
        plus: "+".into(),
        quality: quality.to_vec(),
    };
    fitting_candidates(&Prepared::new(&read).orientations[0], target)
        .into_iter()
        .next()
}

fn fitting_candidates(query: &OrientedQuery, target: &[u8]) -> Vec<Hit> {
    let max_edits = edit_budget(query.sequence.len());
    if query.sequence.is_empty()
        || target.is_empty()
        || query.sequence.len() > target.len() + max_edits
        || query.blocks.as_ref().is_some_and(|b| !b.possible(target))
    {
        return Vec::new();
    }
    let empty = Cell {
        score: NEG,
        start: 0,
        edits: 0,
    };
    // One in-place row retains M, I and best(M, I, D). The previous diagonal
    // and the current left M/D are scalars, so no second row or full-row
    // clearing is needed. The recurrence and its tie order are unchanged.
    #[derive(Clone, Copy)]
    struct Column {
        matched: Cell,
        inserted: Cell,
        best: Cell,
    }
    let empty_column = Column {
        matched: empty,
        inserted: empty,
        best: empty,
    };
    let mut row: Vec<_> = (0..=target.len())
        .map(|start| {
            let initial = Cell {
                score: 0.0,
                start,
                edits: 0,
            };
            Column {
                matched: initial,
                inserted: empty,
                best: best(best(initial, empty), empty),
            }
        })
        .collect();
    let mut previous_hi = target.len();
    for (i, (&base, &[match_score, mismatch_score])) in
        query.sequence.iter().zip(&query.emissions).enumerate()
    {
        let lo = (i + 1).saturating_sub(max_edits).max(1);
        let hi = (target.len() + max_edits)
            .saturating_sub(query.sequence.len() - i - 1)
            .min(target.len());
        // Both band edges move by at most one column per query base. An
        // entering upper column had no previous-row state; a leaving lower
        // column is cleared only after its diagonal has been saved.
        if hi > previous_hi {
            row[hi] = empty_column;
        }
        previous_hi = hi;
        let mut diagonal = row[lo - 1].best;
        if lo > 1 {
            row[lo - 1] = empty_column;
        }
        let inserted = best(
            row[0].matched.advance(GAP_OPEN, 1),
            row[0].inserted.advance(GAP_EXTEND, 1),
        );
        row[0] = Column {
            matched: empty,
            inserted,
            best: best(best(empty, inserted), empty),
        };
        let mut left_m = empty;
        let mut left_d = empty;
        for (column, &target_base) in row[lo..=hi].iter_mut().zip(&target[lo - 1..hi]) {
            let matched = diagonal.advance(
                if !b"ACGT".contains(&target_base) {
                    0.25f64.ln()
                } else if base == target_base {
                    match_score
                } else {
                    mismatch_score
                },
                usize::from(base != target_base),
            );
            let inserted = best(
                column.matched.advance(GAP_OPEN, 1),
                column.inserted.advance(GAP_EXTEND, 1),
            );
            let deleted = best(left_m.advance(GAP_OPEN, 1), left_d.advance(GAP_EXTEND, 1));
            diagonal = column.best;
            *column = Column {
                matched,
                inserted,
                best: best(best(matched, inserted), deleted),
            };
            left_m = matched;
            left_d = deleted;
        }
    }
    let mut hits: Vec<_> = row
        .iter()
        .enumerate()
        .skip(1)
        .filter_map(|(end, column)| {
            let cell = column.best;
            (cell.score > NEG / 2.0 && end > cell.start && cell.edits <= max_edits).then_some(Hit {
                start: cell.start,
                end,
                reverse: false,
                logp: cell.score,
                edits: cell.edits,
                traceback: None,
            })
        })
        .collect();
    hits.sort_by(|a, b| {
        b.logp
            .total_cmp(&a.logp)
            .then(a.edits.cmp(&b.edits))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    let maximum = hits.first().map_or(NEG, |h| h.logp);
    let mut starts = BTreeSet::new();
    hits.retain(|h| h.logp >= maximum - 16.0 && starts.insert(h.start));
    hits
}

/// Trace one accepted placement. Only sequenced target bases with sufficient
/// quality contribute matching depth; insertions, deletions and mate gaps do not.
#[cfg(test)]
pub(crate) fn trace_bases(
    read: &Record,
    hit: &Hit,
    target: &[u8],
    minq: u8,
) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let trace = trace_alignment(read, hit, target, minq);
    (trace.observed, trace.matching, trace.deletions)
}

#[derive(Default)]
pub(crate) struct BaseTrace {
    pub ambiguous_trace: bool,
    pub observed: Vec<usize>,
    pub matching: Vec<usize>,
    /// Oriented query offset and unfolded target position for known matching bases.
    pub matching_pairs: Vec<(usize, usize)>,
    pub ambiguous_contacts: usize,
    pub deletions: Vec<usize>,
    /// Reference positions whose preceding base is adjacent in the read too.
    pub matching_junctions: Vec<usize>,
}

impl BaseTrace {
    /// Every target base and every intervening junction must be sequenced in
    /// one contiguous matching query run. A deletion, insertion, mismatch or
    /// quality-masked base cannot certify a new sequence adjacency.
    pub(crate) fn matches_contiguous(&self, start: usize, end: usize) -> bool {
        !self.ambiguous_trace
            && end > start
            && self
                .matching
                .iter()
                .filter(|&&p| (start..end).contains(&p))
                .count()
                == end - start
            && self
                .matching_junctions
                .iter()
                .filter(|&&p| (start + 1..end).contains(&p))
                .count()
                == end - start - 1
    }
}

pub(crate) fn trace_alignment(read: &Record, hit: &Hit, target: &[u8], minq: u8) -> BaseTrace {
    let query = if hit.reverse {
        kmer::reverse_complement(&read.sequence)
    } else {
        read.sequence.iter().map(u8::to_ascii_uppercase).collect()
    };
    let quality: Vec<_> = if hit.reverse {
        read.quality.iter().rev().copied().collect()
    } else {
        read.quality.clone()
    };
    if let Some(trace) = &hit.traceback {
        return certified_trace(&query, &quality, hit, target, minq, trace);
    }
    if hit.edits == 0 {
        let observed: Vec<_> = quality
            .iter()
            .enumerate()
            .filter(|(_, q)| q.saturating_sub(33) >= minq)
            .map(|(i, _)| hit.start + i)
            .collect();
        let matching_pairs: Vec<_> = observed
            .iter()
            .filter_map(|&p| {
                let q = p - hit.start;
                (b"ACGT".contains(&query[q]) && target.get(p) == Some(&query[q])).then_some((q, p))
            })
            .collect();
        let matching: Vec<_> = matching_pairs.iter().map(|p| p.1).collect();
        let matching_junctions = matching
            .windows(2)
            .filter_map(|p| (p[1] == p[0] + 1).then_some(p[1]))
            .collect();
        let ambiguous_contacts = observed
            .iter()
            .filter(|&&p| !b"ACGT".contains(&query[p - hit.start]) || !b"ACGT".contains(&target[p]))
            .count();
        return BaseTrace {
            ambiguous_trace: false,
            observed,
            matching,
            matching_pairs,
            matching_junctions,
            ambiguous_contacts,
            deletions: Vec::new(),
        };
    }
    let reference = &target[hit.start..hit.end];
    let width = reference.len() + 1;
    let mut trace = vec![[0u8; 3]; (query.len() + 1) * width];
    let mut previous = [vec![NEG; width], vec![NEG; width], vec![NEG; width]];
    previous[0][0] = 0.0;
    for (j, cell) in trace.iter_mut().enumerate().take(width).skip(1) {
        previous[2][j] = GAP_OPEN + (j - 1) as f64 * GAP_EXTEND;
        cell[2] = if j == 1 { 0 } else { 2 };
    }
    let choose = |values: [f64; 3]| {
        (0..3)
            .max_by(|&a, &b| values[a].total_cmp(&values[b]).then(b.cmp(&a)))
            .unwrap()
    };
    let band = (query.len() / 10).max(3);
    for i in 1..=query.len() {
        let mut current = [vec![NEG; width], vec![NEG; width], vec![NEG; width]];
        current[1][0] = GAP_OPEN + (i - 1) as f64 * GAP_EXTEND;
        trace[i * width][1] = if i == 1 { 0 } else { 1 };
        for j in i.saturating_sub(band).max(1)..=(i + band).min(reference.len()) {
            let source = choose([previous[0][j - 1], previous[1][j - 1], previous[2][j - 1]]);
            current[0][j] =
                previous[source][j - 1] + emission(query[i - 1], reference[j - 1], quality[i - 1]);
            trace[i * width + j][0] = source as u8;
            let source = if previous[0][j] + GAP_OPEN >= previous[1][j] + GAP_EXTEND {
                0
            } else {
                1
            };
            current[1][j] = previous[source][j] + if source == 0 { GAP_OPEN } else { GAP_EXTEND };
            trace[i * width + j][1] = source as u8;
            let source = if current[0][j - 1] + GAP_OPEN >= current[2][j - 1] + GAP_EXTEND {
                0
            } else {
                2
            };
            current[2][j] =
                current[source][j - 1] + if source == 0 { GAP_OPEN } else { GAP_EXTEND };
            trace[i * width + j][2] = source as u8;
        }
        previous = current;
    }
    let (mut i, mut j) = (query.len(), reference.len());
    let mut state = choose([previous[0][j], previous[1][j], previous[2][j]]);
    let (mut observed, mut matching, mut deletions) = (Vec::new(), Vec::new(), Vec::new());
    let mut matching_junctions = Vec::new();
    let mut matching_pairs = Vec::new();
    let mut ambiguous_contacts = 0;
    let mut previous_match = None;
    while i > 0 || j > 0 {
        let from = trace[i * width + j][state] as usize;
        match state {
            0 => {
                if i == 0 || j == 0 {
                    break;
                }
                if quality[i - 1].saturating_sub(33) >= minq {
                    observed.push(hit.start + j - 1);
                    ambiguous_contacts += usize::from(
                        !b"ACGT".contains(&query[i - 1]) || !b"ACGT".contains(&reference[j - 1]),
                    );
                    if b"ACGT".contains(&query[i - 1]) && query[i - 1] == reference[j - 1] {
                        matching_pairs.push((i - 1, hit.start + j - 1));
                        matching.push(hit.start + j - 1);
                        if previous_match == Some((i, j)) {
                            matching_junctions.push(hit.start + j);
                        }
                        previous_match = Some((i - 1, j - 1));
                    }
                }
                i -= 1;
                j -= 1;
            }
            1 => {
                if i == 0 {
                    break;
                }
                i -= 1;
            }
            _ => {
                if j == 0 {
                    break;
                }
                deletions.push(hit.start + j - 1);
                j -= 1;
            }
        }
        state = from;
    }
    if i > 0 || j > 0 {
        return BaseTrace::default();
    }
    BaseTrace {
        ambiguous_trace: false,
        observed,
        matching,
        matching_pairs,
        ambiguous_contacts,
        deletions,
        matching_junctions,
    }
}

fn certified_trace(
    query: &[u8],
    quality: &[u8],
    hit: &Hit,
    target: &[u8],
    minq: u8,
    trace: &profile::Traceback,
) -> BaseTrace {
    let mut result = BaseTrace {
        ambiguous_trace: trace.ambiguous,
        ..BaseTrace::default()
    };
    let (mut q, mut p) = (0, hit.start);
    let mut previous_match = None;
    for &step in &trace.steps {
        match step {
            0 => {
                if q >= query.len() || p >= target.len() {
                    return BaseTrace::default();
                }
                if quality[q].saturating_sub(33) >= minq {
                    result.observed.push(p);
                    result.ambiguous_contacts +=
                        usize::from(!b"ACGT".contains(&query[q]) || !b"ACGT".contains(&target[p]));
                    if b"ACGT".contains(&query[q]) && query[q] == target[p] {
                        result.matching.push(p);
                        result.matching_pairs.push((q, p));
                        if previous_match == Some((q.wrapping_sub(1), p.wrapping_sub(1))) {
                            result.matching_junctions.push(p);
                        }
                        previous_match = Some((q, p));
                    }
                }
                q += 1;
                p += 1;
            }
            1 => q += 1,
            _ => {
                result.deletions.push(p);
                p += 1;
            }
        }
    }
    if q != query.len() || p != hit.end {
        return BaseTrace::default();
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    fn read(s: &[u8]) -> Record {
        Record {
            header: "@a/1".into(),
            sequence: s.to_vec(),
            plus: "+".into(),
            quality: vec![b'I'; s.len()],
        }
    }
    const SEQUENCE: &[u8] = b"GACTACGTGATAGCTGACTGCATGCTAGCTACGTTAGCGATCGTACCGATGCTACGATCGATCATGCATCGTAGCTAGCATGCTACGTACGATCG";
    #[test]
    fn unknown_bases_cannot_supply_exact_or_gapped_junction_observations() {
        let mut target = SEQUENCE.to_vec();
        target[37] = b'N';
        let index = Index::new(&target, false, 80);
        for insertion in [false, true] {
            let mut query = target[7..80].to_vec();
            if insertion {
                query.insert(10, b'C');
            }
            let query = read(&query);
            let placements = index.align(&query);
            let hit = placements
                .hits
                .iter()
                .max_by(|a, b| a.logp.total_cmp(&b.logp))
                .unwrap();
            assert_eq!(hit.edits, usize::from(insertion));
            let trace = trace_alignment(&query, hit, &index.sequence, 0);
            assert!(!trace.matching.contains(&37));
            assert!(!trace.matches_contiguous(30, 45));
            assert!(trace.ambiguous_contacts > 0);
            assert!(trace.matches_contiguous(45, 60));
        }
    }
    #[test]
    fn observed_junction_requires_adjacent_query_and_target_bases() {
        let index = Index::new(SEQUENCE, false, 80);
        for mutation in 0..4 {
            let mut sequence = SEQUENCE[7..80].to_vec();
            match mutation {
                1 => {
                    sequence.remove(30);
                }
                2 => sequence.insert(30, b'C'),
                3 => sequence = kmer::reverse_complement(&sequence),
                _ => {}
            }
            let query = read(&sequence);
            let placements = index.align(&query);
            let hit = placements
                .hits
                .iter()
                .max_by(|a, b| a.logp.total_cmp(&b.logp))
                .unwrap();
            assert_eq!((hit.start, hit.end), (7, 80));
            let trace = trace_alignment(&query, hit, &index.sequence, 0);
            assert_eq!(
                trace.matches_contiguous(17, 57),
                mutation == 0 || mutation == 3
            );
        }
    }
    #[test]
    fn low_quality_seeds_reach_alignment_and_still_receive_lower_scores() {
        let index = Index::new(SEQUENCE, false, 80);
        let high = read(&SEQUENCE[7..80]);
        let mut low = read(&SEQUENCE[7..80]);
        low.quality.fill(b'&'); // Q5: previously every 17-mer seed was discarded.
        let high_hits = index.align(&high);
        let low_hits = index.align(&low);
        let hit = |hits: &Placements| {
            hits.hits
                .iter()
                .find(|h| h.start == 7 && h.end == 80 && !h.reverse && h.edits == 0)
                .map(|h| h.logp)
                .unwrap()
        };
        assert!(hit(&low_hits) < hit(&high_hits));
    }
    #[test]
    fn full_query_mismatch_insertion_deletion_and_reverse() {
        let index = Index::new(SEQUENCE, false, 80);
        for mutation in 0..4 {
            let mut s = SEQUENCE[7..80].to_vec();
            match mutation {
                0 => s[30] = b'A',
                1 => s.insert(30, b'C'),
                2 => {
                    s.remove(30);
                }
                _ => s = kmer::reverse_complement(&s),
            }
            let hits = index.align(&read(&s));
            let hit = hits
                .hits
                .iter()
                .max_by(|a, b| a.logp.total_cmp(&b.logp))
                .unwrap();
            assert_eq!((hit.start, hit.end), (7, 80));
            assert_eq!(hit.reverse, mutation == 3);
            assert_eq!(hit.edits, usize::from(mutation != 3));
        }
    }
    #[test]
    fn cannot_hide_unmatched_query_ends_by_soft_clipping() {
        let target = &SEQUENCE[15..70];
        let aligned = fitting(SEQUENCE, &vec![b'I'; SEQUENCE.len()], target);
        assert!(aligned.is_none_or(|hit| hit.edits >= SEQUENCE.len() - target.len()));
        assert!(Index::new(target, false, SEQUENCE.len())
            .align(&read(SEQUENCE))
            .hits
            .is_empty());
    }
    #[test]
    fn circular_coordinates_and_identical_repeat_placements() {
        let index = Index::new(SEQUENCE, true, 60);
        let q = [&SEQUENCE[70..], &SEQUENCE[..35]].concat();
        let result = index.align(&read(&q));
        assert!(result
            .hits
            .iter()
            .any(|h| h.start == 70 && h.end > SEQUENCE.len() && h.edits == 0));
        let repeated = [SEQUENCE, b"TTTTTTTTTTTTTTTTTTTT", SEQUENCE].concat();
        assert_eq!(
            Index::new(&repeated, false, 60)
                .align(&read(&SEQUENCE[10..65]))
                .hits
                .iter()
                .filter(|h| h.edits == 0)
                .count(),
            2
        );
    }

    #[test]
    fn shared_window_cache_preserves_offsets_and_read_qualities() {
        let mut query = read(&SEQUENCE[20..80]);
        query.sequence[30] = if query.sequence[30] == b'A' {
            b'C'
        } else {
            b'A'
        };
        let first = Index::new(SEQUENCE, false, 100);
        let shifted = Index::new(&[b"GGCCG".as_slice(), SEQUENCE].concat(), false, 100);
        let mut prepared = Prepared::new(&query);
        let a = first.align_prepared(&mut prepared);
        let b = shifted.align_prepared(&mut prepared);
        let a = a
            .hits
            .iter()
            .max_by(|a, b| a.logp.total_cmp(&b.logp))
            .unwrap();
        let b = b
            .hits
            .iter()
            .max_by(|a, b| a.logp.total_cmp(&b.logp))
            .unwrap();
        assert_eq!((a.start, a.end), (20, 80));
        assert_eq!((b.start, b.end), (25, 85));
        assert_eq!(a.logp, b.logp);
        query.quality[30] = b'#';
        let low_quality = first.align(&query);
        assert!(low_quality
            .hits
            .iter()
            .any(|h| h.start == 20 && h.logp > a.logp + 3.0));
    }
    #[test]
    fn anchored_seedless_rescue_handles_distributed_errors() {
        let mut state = 813789u64;
        let target: Vec<_> = (0..180)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state % 4) as usize]
            })
            .collect();
        let mut query = read(&target[40..140]);
        for i in (7..100).step_by(11) {
            query.sequence[i] = if query.sequence[i] == b'A' {
                b'C'
            } else {
                b'A'
            };
            query.quality[i] = b'+';
        }
        let index = Index::new(&target, false, 100);
        assert!(index.align(&query).hits.is_empty());
        let result = index.rescue_near(&mut Prepared::new(&query), 15, 165, false);
        assert!(result
            .hits
            .iter()
            .any(|h| (h.start, h.end) == (40, 140) && h.edits == 9));
        let repeated = [target.clone(), vec![b'T'; 20], target.clone()].concat();
        let index = Index::new(&repeated, false, 100);
        let result = index.rescue_near(&mut Prepared::new(&query), 0, repeated.len(), false);
        assert!(result.hits.iter().any(|h| (h.start, h.end) == (40, 140)));
        assert!(result.hits.iter().any(|h| (h.start, h.end) == (240, 340)));
    }
    #[test]
    fn gapped_trace_counts_only_sequenced_matching_high_quality_bases() {
        let index = Index::new(SEQUENCE, false, 100);
        for mutation in 0..3 {
            for reverse in [false, true] {
                let mut q = read(&SEQUENCE[7..80]);
                match mutation {
                    0 => {
                        q.sequence[30] = if q.sequence[30] == b'A' { b'C' } else { b'A' };
                    }
                    1 => {
                        q.sequence.insert(30, b'C');
                        q.quality.insert(30, b'I');
                    }
                    _ => {
                        q.sequence.remove(30);
                        q.quality.remove(30);
                    }
                }
                q.quality[10] = b'!';
                if reverse {
                    q.sequence = kmer::reverse_complement(&q.sequence);
                    q.quality.reverse();
                }
                let hits = index.align(&q);
                let hit = hits
                    .hits
                    .iter()
                    .max_by(|a, b| a.logp.total_cmp(&b.logp))
                    .unwrap();
                let (observed, matching, deleted) = trace_bases(&q, hit, &index.sequence, 20);
                assert!(!observed.contains(&17));
                assert!(deleted
                    .iter()
                    .all(|p| !observed.contains(p) && !matching.contains(p)));
                assert_eq!(matching.len(), if mutation == 1 { 72 } else { 71 });
                assert_eq!(observed.len(), if mutation == 2 { 71 } else { 72 });
                assert_eq!(deleted.len(), usize::from(mutation == 2));
            }
        }
    }
}
