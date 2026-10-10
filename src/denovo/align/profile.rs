//! Endpoint profiles for fragment comparison. A start, end and edit budget are
//! part of the state; neither a better sequence score at another endpoint nor
//! an unconstrained traceback can discard a legal fragment explanation.
use super::{emission, Hit, Placements, GAP_EXTEND, GAP_OPEN, NEG};
use std::{
    collections::{BTreeSet, HashMap},
    sync::Arc,
};

pub(crate) const MAX_CELLS: usize = 8_000_000;
pub(crate) const MAX_PROFILES: usize = 4096;
const MAX_MATRIX_SLOTS: usize = 1_000_000;
const TIE: f64 = 1e-10;

mod certified;

pub(crate) fn report() -> serde_json::Value {
    serde_json::json!({
        "algorithm":"full-query affine endpoint/edit profiles with in-DP exact pivot",
        "edit_budget":"max(3, floor(read_length / 10))",
        "max_dp_states_per_window":MAX_CELLS,"work_accounting":"edit states plus prefix/suffix feasibility cells","max_matrix_slots":MAX_MATRIX_SLOTS,
        "max_endpoint_profiles":MAX_PROFILES,
        "search":"all feasible starts and endpoints within declared mate window, strand and edit budget",
        "incomplete":"omit the fragment from every compared hypothesis; never score as absence",
        "tied_tracebacks":"retain endpoint score; suppress junction certification",
        "global_exhaustiveness":false
    })
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Constraint {
    pub query: usize,
    pub target: usize,
    pub length: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Traceback {
    /// 0 consumes both sequences, 1 only query, 2 only target.
    pub steps: Vec<u8>,
    /// Tied optimal traces cannot certify a junction using an arbitrary trace.
    pub ambiguous: bool,
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub(super) struct Key {
    pub target: Vec<u8>,
    pub constraint: Option<Constraint>,
    pub edits: usize,
}

/// Owned by one prepared read/orientation, so query and qualities are fixed.
/// Cache complete start subproblems only, including their original work charge.
/// Cache capacity controls reuse, never the alignment domain or unknown status.
#[derive(Default)]
pub(super) struct Subproblems {
    entries: HashMap<Key, CachedStart>,
    hits: usize,
}

struct CachedStart {
    work: usize,
    hits: Vec<Hit>,
}

impl Subproblems {
    fn remember(&mut self, key: Key, work: usize, hits: &[Hit], start: usize) {
        if self.entries.len() >= 256 || self.hits + hits.len() > 8192 {
            return;
        }
        let hits: Vec<_> = hits
            .iter()
            .map(|h| Hit {
                start: h.start - start,
                end: h.end - start,
                ..h.clone()
            })
            .collect();
        self.hits += hits.len();
        self.entries.insert(key, CachedStart { work, hits });
    }
}

#[derive(Clone, Copy)]
struct Cell {
    score: f64,
    from: u8,
    ambiguous: bool,
}

impl Default for Cell {
    fn default() -> Self {
        Self {
            score: NEG,
            from: 0,
            ambiguous: false,
        }
    }
}

fn choose(candidates: impl IntoIterator<Item = (Cell, f64, u8)>) -> Cell {
    let mut best = Cell::default();
    for (source, cost, from) in candidates {
        if source.score <= NEG / 2.0 {
            continue;
        }
        let score = source.score + cost;
        let tied = (score - best.score).abs() <= TIE;
        if score > best.score {
            best = Cell {
                score,
                from,
                ambiguous: source.ambiguous || tied,
            };
        } else if tied {
            best.ambiguous = true;
        }
    }
    best
}

/// E edits can disrupt at most E of E+1 disjoint query blocks. Every accepted
/// full-query alignment therefore contains one exact block. Prefix indels can
/// move its implied start by at most E. This is a feasibility filter, not a
/// sequence-score cutoff (ambiguous bases use the same equality as the DP).
fn starts(query: &[u8], target: &[u8], edits: usize) -> BTreeSet<usize> {
    let minimum = query.len().saturating_sub(edits).max(1);
    if target.len() < minimum {
        return BTreeSet::new();
    }
    let last = target.len() - minimum;
    if query.len() <= edits {
        return (0..=last).collect();
    }
    let mut found = BTreeSet::new();
    for block in 0..=edits {
        let q = block * query.len() / (edits + 1);
        let end = (block + 1) * query.len() / (edits + 1);
        for (p, word) in target.windows(end - q).enumerate() {
            if word == &query[q..end] {
                let raw = p as isize - q as isize;
                let lo = (raw - edits as isize).max(0) as usize;
                let hi = (raw + edits as isize).min(last as isize);
                if hi >= lo as isize {
                    found.extend(lo..=hi as usize);
                }
            }
        }
    }
    found
}

struct Matrix {
    band: usize,
    edits: usize,
    cells: Vec<[Cell; 3]>,
    retained: Vec<[u64; 3]>,
    trace: Vec<[u8; 3]>,
}

/// Unit-cost prefix/suffix distances relax affine transitions and the pivot.
/// They are lower bounds only: a state with e edits can be skipped if even
/// this relaxed suffix needs more than E-e. No score-based pruning is used.
fn edit_bounds(
    query: &[u8],
    target: &[u8],
    edits: usize,
    suffix: Option<(&SuffixBounds, usize)>,
) -> Vec<[usize; 2]> {
    let band = 2 * edits + 1;
    let slot = |i: usize, j: usize| i * band + j + edits - i;
    let mut bounds = vec![[edits + 1; 2]; (query.len() + 1) * band];
    let get = |bounds: &[[usize; 2]], i: usize, j: usize, side: usize| {
        if i.abs_diff(j) <= edits {
            bounds[slot(i, j)][side]
        } else {
            edits + 1
        }
    };
    for i in 0..=query.len() {
        for j in i.saturating_sub(edits)..=(i + edits).min(target.len()) {
            let mut best = if i == 0 {
                j
            } else if j == 0 {
                i
            } else {
                edits + 1
            };
            if i > 0 && j > 0 {
                best = best.min(
                    get(&bounds, i - 1, j - 1, 0) + usize::from(query[i - 1] != target[j - 1]),
                );
            }
            if i > 0 {
                best = best.min(get(&bounds, i - 1, j, 0) + 1);
            }
            if j > 0 {
                best = best.min(get(&bounds, i, j - 1, 0) + 1);
            }
            bounds[slot(i, j)][0] = best.min(edits + 1);
        }
    }
    if let Some((suffix, start)) = suffix {
        for i in 0..=query.len() {
            for j in i.saturating_sub(edits)..=(i + edits).min(target.len()) {
                bounds[slot(i, j)][1] = suffix.cells[i * suffix.width + start + j];
            }
        }
        return bounds;
    }
    for i in (0..=query.len()).rev() {
        for j in (i.saturating_sub(edits)..=(i + edits).min(target.len())).rev() {
            let mut best = if i == query.len() { 0 } else { edits + 1 };
            if i < query.len() && j < target.len() {
                best = best.min(get(&bounds, i + 1, j + 1, 1) + usize::from(query[i] != target[j]));
            }
            if i < query.len() {
                best = best.min(get(&bounds, i + 1, j, 1) + 1);
            }
            if j < target.len() {
                best = best.min(get(&bounds, i, j + 1, 1) + 1);
            }
            bounds[slot(i, j)][1] = best.min(edits + 1);
        }
    }
    bounds
}

/// Restrict the fallback DP to paths that can reach an uncertified endpoint.
/// Every accepted path to those endpoints still passes this edit lower bound.
fn restrict_suffix(
    query: &[u8],
    target: &[u8],
    edits: usize,
    endpoints: &[bool],
    bounds: &mut [[usize; 2]],
) {
    let band = 2 * edits + 1;
    let slot = |i: usize, j: usize| i * band + j + edits - i;
    let get = |bounds: &[[usize; 2]], i: usize, j: usize| {
        if i.abs_diff(j) <= edits {
            bounds[slot(i, j)][1]
        } else {
            edits + 1
        }
    };
    for i in (0..=query.len()).rev() {
        for j in (i.saturating_sub(edits)..=(i + edits).min(target.len())).rev() {
            let mut best = if i == query.len() && endpoints[j] {
                0
            } else {
                edits + 1
            };
            if i < query.len() && j < target.len() {
                best = best.min(get(bounds, i + 1, j + 1) + usize::from(query[i] != target[j]));
            }
            if i < query.len() {
                best = best.min(get(bounds, i + 1, j) + 1);
            }
            if j < target.len() {
                best = best.min(get(bounds, i, j + 1) + 1);
            }
            bounds[slot(i, j)][1] = best.min(edits + 1);
        }
    }
}

/// Free-end suffix distances are shared by all starts in one target window.
/// A suffix that can pass prefix + suffix <= E cannot leave the original
/// start-relative band: its displacement is bounded by the same total edits.
/// Thus every feasible state's suffix bound is identical to the banded DP;
/// other states still fail the bound. Scores and logical work stay unchanged.
struct SuffixBounds {
    width: usize,
    cells: Vec<usize>,
}

impl SuffixBounds {
    fn new(query: &[u8], target: &[u8], edits: usize) -> Self {
        let width = target.len() + 1;
        let mut cells = vec![0; (query.len() + 1) * width];
        for i in (0..query.len()).rev() {
            cells[i * width + target.len()] = (query.len() - i).min(edits + 1);
            for j in (0..target.len()).rev() {
                cells[i * width + j] = (cells[(i + 1) * width + j + 1]
                    + usize::from(query[i] != target[j]))
                .min(cells[(i + 1) * width + j] + 1)
                .min(cells[i * width + j + 1] + 1)
                .min(edits + 1);
            }
        }
        Self { width, cells }
    }
}

impl Matrix {
    fn score_slot(&self, i: usize, j: usize, e: usize) -> usize {
        let row = if self.trace.is_empty() { i } else { i & 1 };
        (row * self.band + j + self.edits - i) * (self.edits + 1) + e
    }

    fn slot(&self, i: usize, j: usize, e: usize) -> Option<usize> {
        (i.abs_diff(j) <= self.edits && e <= self.edits)
            .then(|| (i * self.band + j + self.edits - i) * (self.edits + 1) + e)
    }
    fn get(&self, i: usize, j: usize, e: usize) -> [Cell; 3] {
        self.slot(i, j, e).map_or([Cell::default(); 3], |p| {
            let score_slot = self.score_slot(i, j, e);
            if self.retained.is_empty() {
                self.cells[score_slot]
            } else {
                let mask = self.retained[p / (self.edits + 1)];
                std::array::from_fn(|s| {
                    if mask[s] & (1 << e) != 0 {
                        self.cells[score_slot][s]
                    } else {
                        Cell::default()
                    }
                })
            }
        })
    }

    fn set(&mut self, i: usize, j: usize, e: usize, value: [Cell; 3]) {
        let p = self.slot(i, j, e).unwrap();
        let score_slot = self.score_slot(i, j, e);
        self.cells[score_slot] = value;
        if !self.trace.is_empty() {
            self.trace[p] = value.map(|c| c.from);
        }
    }

    fn from(&self, i: usize, j: usize, e: usize, state: usize) -> usize {
        let p = self.slot(i, j, e).unwrap();
        if self.trace.is_empty() {
            self.cells[p][state].from as usize
        } else {
            self.trace[p][state] as usize
        }
    }
}

#[cfg(test)]
pub(super) fn fit(
    query: &[u8],
    quality: &[u8],
    target: &[u8],
    constraint: Option<Constraint>,
    edits: usize,
    cell_limit: usize,
) -> Placements {
    fit_impl::<false>(query, quality, target, constraint, edits, cell_limit, None)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn fit_cached(
    query: &[u8],
    quality: &[u8],
    target: &[u8],
    constraint: Option<Constraint>,
    edits: usize,
    cell_limit: usize,
    cache: &mut Subproblems,
) -> Placements {
    fit_impl::<true>(
        query,
        quality,
        target,
        constraint,
        edits,
        cell_limit,
        Some(cache),
    )
}

#[allow(clippy::too_many_arguments)]
fn fit_impl<const DOMINANCE: bool>(
    query: &[u8],
    quality: &[u8],
    target: &[u8],
    constraint: Option<Constraint>,
    edits: usize,
    cell_limit: usize,
    mut cache: Option<&mut Subproblems>,
) -> Placements {
    let mut result = Placements {
        complete: true,
        ..Placements::default()
    };
    if query.is_empty() || target.is_empty() {
        return result;
    }
    if let Some(c) = constraint {
        if c.length == 0
            || query.get(c.query..c.query.saturating_add(c.length))
                != target.get(c.target..c.target.saturating_add(c.length))
            || c.query.saturating_add(c.length) > query.len()
            || c.target.saturating_add(c.length) > target.len()
        {
            return result;
        }
    }
    let band = 2 * edits + 1;
    let Some(size) = (query.len() + 1)
        .checked_mul(band)
        .and_then(|n| n.checked_mul(edits + 1))
    else {
        result.exhausted = true;
        result.complete = false;
        return result;
    };
    if size > cell_limit.min(MAX_MATRIX_SLOTS) {
        result.exhausted = true;
        result.complete = false;
        return result;
    }
    let candidates = if let Some(c) = constraint {
        let raw = c.target as isize - c.query as isize;
        let lo = (raw - edits as isize).max(0) as usize;
        let hi = (raw + edits as isize).min(c.target as isize);
        if hi < lo as isize {
            BTreeSet::new()
        } else {
            (lo..=hi as usize).collect()
        }
    } else {
        starts(query, target, edits)
    };
    let mut matrix = Matrix {
        band,
        edits,
        cells: Vec::new(),
        retained: Vec::new(),
        trace: Vec::new(),
    };
    let emissions: Vec<_> = query
        .iter()
        .zip(quality)
        .map(|(&q, &v)| b"ACGTN".map(|t| emission(q, t, v)))
        .collect();
    // A lower-edit state can replay every continuation of a higher-edit state
    // at the same coordinates and affine state. Discard the latter only when
    // its deficit cannot close to TIE under any remaining additions. All costs
    // are nonpositive and bounded; this deliberately loose roundoff allowance
    // covers both paths over their entire possible lengths.
    let steps = (2 * query.len() + edits) as f64;
    let maximum_cost = emissions
        .iter()
        .flatten()
        .map(|v| v.abs())
        .fold(GAP_OPEN.abs().max(GAP_EXTEND.abs()), f64::max);
    let dominance_margin = 4.0 * TIE + 16.0 * f64::EPSILON * steps * steps * maximum_cost;
    let target_bases: Vec<_> = target
        .iter()
        .map(|b| match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        })
        .collect();
    let mut work = 0usize;
    let sparse = DOMINANCE && edits < 64;
    let share_suffix = DOMINANCE
        && query.len() >= 50
        && candidates.len().saturating_mul(band) > target.len()
        && (query.len() + 1)
            .checked_mul(target.len() + 1)
            .is_some_and(|size| size <= MAX_MATRIX_SLOTS);
    let mut shared_suffix = None;
    'starts: for start in candidates {
        let length = (query.len() + edits).min(target.len() - start);
        let key = cache.as_ref().map(|_| Key {
            target: target[start..start + length].to_vec(),
            constraint: constraint.map(|c| Constraint {
                target: c.target - start,
                ..c
            }),
            edits,
        });
        if let Some(saved) = key
            .as_ref()
            .and_then(|key| cache.as_ref().and_then(|c| c.entries.get(key)))
        {
            if saved.work <= cell_limit - work {
                // Preserve the original logical budget even when no DP cells
                // are recomputed. Warm caches must not turn unknown into known.
                work += saved.work;
                for hit in &saved.hits {
                    result.hits.push(Hit {
                        start: hit.start + start,
                        end: hit.end + start,
                        ..hit.clone()
                    });
                    if result.hits.len() > MAX_PROFILES {
                        result.exhausted = true;
                        result.complete = false;
                        break 'starts;
                    }
                }
                continue;
            }
        }
        let work_before = work;
        let hits_before = result.hits.len();
        if share_suffix && shared_suffix.is_none() {
            shared_suffix = Some(SuffixBounds::new(query, target, edits));
        }
        let mut bounds = edit_bounds(
            query,
            &target[start..start + length],
            edits,
            shared_suffix.as_ref().map(|suffix| (suffix, start)),
        );
        work += 2 * bounds.len();
        if work > cell_limit {
            result.exhausted = true;
            result.complete = false;
            break;
        }
        if bounds[edits][1] > edits {
            if let (Some(cache), Some(key)) = (cache.as_mut(), key) {
                cache.remember(key, work - work_before, &[], start);
            }
            continue;
        }
        let mut certified = None;
        if DOMINANCE {
            let mut charge = 0;
            for i in 0..=query.len() {
                for j in i.saturating_sub(edits)..=(i + edits).min(length) {
                    if i == 0 && j == 0 {
                        continue;
                    }
                    let [prefix, suffix] = bounds[i * band + j + edits - i];
                    if prefix + suffix <= edits {
                        charge += edits - suffix - prefix + 1;
                    }
                }
            }
            if charge > cell_limit - work {
                result.exhausted = true;
                result.complete = false;
                break;
            }
            work += charge;
            let found = certified::fit_start(
                query,
                target,
                &emissions,
                &target_bases,
                start,
                length,
                constraint,
                edits,
                &bounds,
                dominance_margin,
            );
            if !found.needs_profile.iter().any(|&needed| needed) {
                for hit in found.hits.into_iter().flatten() {
                    result.hits.push(hit);
                    if result.hits.len() > MAX_PROFILES {
                        result.exhausted = true;
                        result.complete = false;
                        break 'starts;
                    }
                }
                if let (Some(cache), Some(key)) = (cache.as_mut(), key) {
                    cache.remember(key, work - work_before, &result.hits[hits_before..], start);
                }
                continue;
            }
            restrict_suffix(
                query,
                &target[start..start + length],
                edits,
                &found.needs_profile,
                &mut bounds,
            );
            certified = Some(found);
        }
        if matrix.cells.is_empty() {
            // Recurrences read only the current and preceding query rows.
            // Retain compact predecessor bytes for traceback at older rows;
            // their floating-point scores are never read again.
            matrix.cells.resize(
                if sparse { 2 * band * (edits + 1) } else { size },
                [Cell::default(); 3],
            );
            if sparse {
                matrix.trace.resize(size, [0; 3]);
            }
        } else if !sparse {
            matrix.cells.fill([Cell::default(); 3]);
        }
        let mut origin = [Cell::default(); 3];
        origin[0].score = 0.0;
        matrix.set(0, 0, 0, origin);
        if sparse {
            if matrix.retained.is_empty() {
                matrix.retained.resize((query.len() + 1) * band, [0; 3]);
            } else {
                matrix.retained.fill([0; 3]);
            }
            matrix.retained[edits][0] = 1;
        }
        for i in 0..=query.len() {
            for j in i.saturating_sub(edits)..=(i + edits).min(length) {
                if i == 0 && j == 0 {
                    continue;
                }
                let q_pivot = constraint.is_some_and(|c| i > c.query && i <= c.query + c.length);
                let t_pivot = constraint.is_some_and(|c| {
                    j > 0 && start + j > c.target && start + j <= c.target + c.length
                });
                let paired = constraint.is_none_or(|c| {
                    !(q_pivot || t_pivot)
                        || (q_pivot && t_pivot && i - c.query == start + j - c.target)
                });
                let mismatch = usize::from(i > 0 && j > 0 && query[i - 1] != target[start + j - 1]);
                let score = if i > 0 && j > 0 {
                    emissions[i - 1][target_bases[start + j - 1]]
                } else {
                    0.0
                };
                let [prefix_edits, suffix_edits] = bounds[i * band + j + edits - i];
                if prefix_edits + suffix_edits > edits {
                    continue;
                }
                let mut mask = [u64::MAX; 3];
                if sparse {
                    let get = |i: usize, j: usize| {
                        if i.abs_diff(j) > edits {
                            [0; 3]
                        } else {
                            matrix.retained[i * band + j + edits - i]
                        }
                    };
                    mask = [0; 3];
                    if i > 0 && j > 0 && paired {
                        let p = get(i - 1, j - 1);
                        mask[0] = (p[0] | p[1] | p[2]) << mismatch;
                    }
                    if i > 0 && !q_pivot {
                        let p = get(i - 1, j);
                        mask[1] = (p[0] | p[1]) << 1;
                    }
                    if j > 0 && !t_pivot {
                        let p = get(i, j - 1);
                        mask[2] = (p[0] | p[2]) << 1;
                    }
                    let domain =
                        (u64::MAX << prefix_edits) & (u64::MAX >> (63 - (edits - suffix_edits)));
                    mask = mask.map(|m| m & domain);
                    // Original logical work was charged before certification;
                    // this fallback may visit only uncertified endpoint paths.
                }
                let mut lower_edit_best = [NEG; 3];
                let mut remaining = mask[0] | mask[1] | mask[2];
                let mut dense = prefix_edits..=edits - suffix_edits;
                while let Some(e) = if sparse {
                    (remaining != 0).then(|| {
                        let e = remaining.trailing_zeros() as usize;
                        remaining &= remaining - 1;
                        e
                    })
                } else {
                    dense.next()
                } {
                    if !DOMINANCE {
                        work += 1;
                        if work > cell_limit {
                            result.exhausted = true;
                            result.complete = false;
                            break 'starts;
                        }
                    }
                    let mut cell = [Cell::default(); 3];
                    if i > 0
                        && j > 0
                        && e >= mismatch
                        && paired
                        && (!sparse || mask[0] & (1 << e) != 0)
                    {
                        let prev = matrix.get(i - 1, j - 1, e - mismatch);
                        cell[0] = choose((0..3).map(|s| (prev[s], score, s as u8)));
                    }
                    if i > 0 && e > 0 && !q_pivot && (!sparse || mask[1] & (1 << e) != 0) {
                        let prev = matrix.get(i - 1, j, e - 1);
                        cell[1] = choose([(prev[0], GAP_OPEN, 0), (prev[1], GAP_EXTEND, 1)]);
                    }
                    if j > 0 && e > 0 && !t_pivot && (!sparse || mask[2] & (1 << e) != 0) {
                        let prev = matrix.get(i, j - 1, e - 1);
                        cell[2] = choose([(prev[0], GAP_OPEN, 0), (prev[2], GAP_EXTEND, 2)]);
                    }
                    if DOMINANCE {
                        for s in 0..3 {
                            if cell[s].score < lower_edit_best[s] - dominance_margin {
                                cell[s] = Cell::default();
                                if sparse {
                                    mask[s] &= !(1 << e);
                                }
                            } else {
                                lower_edit_best[s] = lower_edit_best[s].max(cell[s].score);
                            }
                        }
                    }
                    matrix.set(i, j, e, cell);
                }
                if sparse {
                    matrix.retained[i * band + j + edits - i] = mask;
                }
            }
        }
        for end in query.len().saturating_sub(edits).max(1)..=length {
            if constraint.is_some_and(|c| start + end < c.target + c.length) {
                continue;
            }
            if let Some(found) = certified.as_mut() {
                if !found.needs_profile[end] {
                    if let Some(hit) = found.hits[end].take() {
                        result.hits.push(hit);
                        if result.hits.len() > MAX_PROFILES {
                            result.exhausted = true;
                            result.complete = false;
                            break 'starts;
                        }
                    }
                    continue;
                }
            }
            let mut best = Cell::default();
            let (mut edit, mut state) = (0, 0);
            for e in query.len().abs_diff(end)..=edits {
                for (s, cell) in matrix.get(query.len(), end, e).into_iter().enumerate() {
                    if cell.score <= NEG / 2.0 {
                        continue;
                    }
                    let tied = (cell.score - best.score).abs() <= TIE;
                    if cell.score > best.score {
                        best = Cell {
                            ambiguous: cell.ambiguous || tied,
                            ..cell
                        };
                        (edit, state) = (e, s);
                    } else if tied {
                        best.ambiguous = true;
                    }
                }
            }
            if best.score <= NEG / 2.0 {
                continue;
            }
            let (mut i, mut j, mut e) = (query.len(), end, edit);
            let mut steps = Vec::new();
            while i > 0 || j > 0 {
                let from = matrix.from(i, j, e, state);
                steps.push(state as u8);
                match state {
                    0 => {
                        e -= usize::from(query[i - 1] != target[start + j - 1]);
                        i -= 1;
                        j -= 1;
                    }
                    1 => {
                        e -= 1;
                        i -= 1;
                    }
                    _ => {
                        e -= 1;
                        j -= 1;
                    }
                }
                state = from;
            }
            steps.reverse();
            result.hits.push(Hit {
                start,
                end: start + end,
                reverse: false,
                logp: best.score,
                edits: edit,
                traceback: Some(Arc::new(Traceback {
                    steps,
                    ambiguous: best.ambiguous,
                })),
            });
            if result.hits.len() > MAX_PROFILES {
                result.exhausted = true;
                result.complete = false;
                break 'starts;
            }
        }
        if let (Some(cache), Some(key)) = (cache.as_mut(), key) {
            cache.remember(key, work - work_before, &result.hits[hits_before..], start);
        }
    }
    result
}

/// Merge only identical endpoint states. Alternative endpoints stay separate;
/// conflicting optimal traces suppress junction certification conservatively.
pub(super) fn dedup(result: &mut Placements) {
    result.hits.sort_by(|a, b| {
        (a.start, a.end, a.reverse)
            .cmp(&(b.start, b.end, b.reverse))
            .then(b.logp.total_cmp(&a.logp))
    });
    result.hits.dedup_by(|b, a| {
        if (a.start, a.end, a.reverse) != (b.start, b.end, b.reverse) {
            return false;
        }
        if (a.logp - b.logp).abs() <= TIE {
            if let (Some(x), Some(y)) = (&mut a.traceback, &b.traceback) {
                if x.ambiguous || y.ambiguous || x.steps != y.steps {
                    Arc::make_mut(x).ambiguous = true;
                }
            }
        }
        true
    });
    if result.hits.len() > MAX_PROFILES {
        result.exhausted = true;
        result.complete = false;
    }
}

#[cfg(test)]
mod tests;
