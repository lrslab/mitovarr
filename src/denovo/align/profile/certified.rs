//! Try an affine DP without an edit-count dimension. A legal unrestricted
//! optimum is also the optimum in the bounded domain. Cross-edit near ties
//! cannot be merged safely, so their descendants use the full profile DP.
use super::{Cell, Constraint, Hit, Traceback, GAP_EXTEND, GAP_OPEN, NEG, TIE};
use std::sync::Arc;

pub(super) struct Endpoints {
    pub hits: Vec<Option<Hit>>,
    pub needs_profile: Vec<bool>,
}

#[derive(Clone, Copy)]
struct State {
    cell: Cell,
    edits: usize,
    certified: bool,
}

impl Default for State {
    fn default() -> Self {
        Self {
            cell: Cell::default(),
            edits: 0,
            certified: true,
        }
    }
}

fn choose<const N: usize>(sources: [(State, f64, usize, u8); N], margin: f64) -> State {
    let candidates = sources.map(|(mut value, cost, edits, from)| {
        if value.cell.score > NEG / 2.0 {
            value.cell.score += cost;
            value.cell.from = from;
            value.edits += edits;
        }
        value
    });
    let mut best = State::default();
    for value in candidates {
        if value.cell.score > best.cell.score
            || (value.cell.score == best.cell.score && value.edits < best.edits)
        {
            best = value;
        }
    }
    if best.cell.score <= NEG / 2.0 {
        return best;
    }
    for value in candidates {
        if value.cell.score <= NEG / 2.0 || value.cell.from == best.cell.from {
            continue;
        }
        let difference = (value.cell.score - best.cell.score).abs();
        // A losing near-tied source can hide a different-edit competitor
        // even when its own retained edit count equals the winning count.
        if difference <= margin && (value.edits != best.edits || !value.certified) {
            best.certified = false;
        }
        if value.edits == best.edits && difference <= TIE {
            best.cell.ambiguous = true;
        }
    }
    best
}

#[allow(clippy::too_many_arguments)]
pub(super) fn fit_start(
    query: &[u8],
    target: &[u8],
    emissions: &[[f64; 5]],
    target_bases: &[usize],
    start: usize,
    length: usize,
    constraint: Option<Constraint>,
    edits: usize,
    bounds: &[[usize; 2]],
    margin: f64,
) -> Endpoints {
    let band = 2 * edits + 1;
    let slot = |i: usize, j: usize| i * band + j + edits - i;
    let mut cells = vec![[State::default(); 3]; (query.len() + 1) * band];
    cells[edits][0].cell.score = 0.0;
    let get = |cells: &[[State; 3]], i: usize, j: usize| {
        if i.abs_diff(j) > edits {
            [State::default(); 3]
        } else {
            cells[slot(i, j)]
        }
    };
    for i in 0..=query.len() {
        for j in i.saturating_sub(edits)..=(i + edits).min(length) {
            if i == 0 && j == 0 {
                continue;
            }
            let [prefix, suffix] = bounds[slot(i, j)];
            if prefix + suffix > edits {
                continue;
            }
            let q_pivot = constraint.is_some_and(|c| i > c.query && i <= c.query + c.length);
            let t_pivot = constraint
                .is_some_and(|c| j > 0 && start + j > c.target && start + j <= c.target + c.length);
            let paired = constraint.is_none_or(|c| {
                !(q_pivot || t_pivot) || (q_pivot && t_pivot && i - c.query == start + j - c.target)
            });
            let mut cell = [State::default(); 3];
            if i > 0 && j > 0 && paired {
                let previous = get(&cells, i - 1, j - 1);
                let mismatch = usize::from(query[i - 1] != target[start + j - 1]);
                let score = emissions[i - 1][target_bases[start + j - 1]];
                cell[0] = choose(
                    std::array::from_fn::<_, 3, _>(|s| (previous[s], score, mismatch, s as u8)),
                    margin,
                );
            }
            if i > 0 && !q_pivot {
                let previous = get(&cells, i - 1, j);
                cell[1] = choose(
                    [
                        (previous[0], GAP_OPEN, 1, 0),
                        (previous[1], GAP_EXTEND, 1, 1),
                    ],
                    margin,
                );
            }
            if j > 0 && !t_pivot {
                let previous = get(&cells, i, j - 1);
                cell[2] = choose(
                    [
                        (previous[0], GAP_OPEN, 1, 0),
                        (previous[2], GAP_EXTEND, 1, 2),
                    ],
                    margin,
                );
            }
            cells[slot(i, j)] = cell;
        }
    }
    let mut result = Endpoints {
        hits: vec![None; length + 1],
        needs_profile: vec![false; length + 1],
    };
    for end in query.len().saturating_sub(edits).max(1)..=length {
        if constraint.is_some_and(|c| start + end < c.target + c.length) {
            continue;
        }
        let previous = cells[slot(query.len(), end)];
        let best = choose(
            std::array::from_fn::<_, 3, _>(|s| (previous[s], 0.0, 0, s as u8)),
            margin,
        );
        if best.cell.score <= NEG / 2.0 {
            continue;
        }
        if best.edits > edits || !best.certified {
            result.needs_profile[end] = true;
            continue;
        }
        let (mut i, mut j, mut state) = (query.len(), end, best.cell.from as usize);
        let mut steps = Vec::with_capacity(query.len() + edits);
        while i > 0 || j > 0 {
            let from = cells[slot(i, j)][state].cell.from as usize;
            steps.push(state as u8);
            match state {
                0 => {
                    i -= 1;
                    j -= 1;
                }
                1 => i -= 1,
                _ => j -= 1,
            }
            state = from;
        }
        steps.reverse();
        result.hits[end] = Some(Hit {
            start,
            end: start + end,
            reverse: false,
            logp: best.cell.score,
            edits: best.edits,
            traceback: Some(Arc::new(Traceback {
                steps,
                ambiguous: best.cell.ambiguous,
            })),
        });
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cross_edit_ties_require_full_profiles_but_same_edit_ties_remain_ambiguous() {
        let state = |edits, score| State {
            edits,
            cell: Cell {
                score,
                ..Cell::default()
            },
            certified: true,
        };
        let margin = 1e-7;
        for delta in [0.0, TIE / 2.0, margin / 2.0] {
            let found = choose(
                [
                    (state(2, -10.0), 0.0, 0, 0),
                    (state(4, -10.0 + delta), 0.0, 0, 1),
                ],
                margin,
            );
            assert!(!found.certified);
        }
        let found = choose(
            [(state(2, -10.0), 0.0, 0, 0), (state(2, -10.0), 0.0, 0, 1)],
            margin,
        );
        assert!(found.certified && found.cell.ambiguous);
        assert_eq!(found.cell.from, 0);
        for delta in [0.0, TIE / 2.0, margin / 2.0] {
            let mut unresolved = state(2, -10.0 - delta);
            unresolved.certified = false;
            let found = choose(
                [(state(2, -10.0), 0.0, 0, 0), (unresolved, 0.0, 0, 1)],
                margin,
            );
            assert!(!found.certified);
        }
    }
}
