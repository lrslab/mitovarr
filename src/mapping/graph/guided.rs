//! Bounded proposal searches over raw sample edges, separate from strict paths.
use super::{code, Code, Graph};
use crate::{kmer, mapping::coverage::Depth};
use serde_json::{json, Value};
use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet, HashMap, VecDeque},
};

const MAX_STATES: usize = 25_000;
const MAX_PENDING: usize = 512;
const MAX_REACHABLE: usize = 200_000;

#[derive(Clone)]
struct State {
    priority: f64,
    cost: f64,
    node: Code,
    sequence: Vec<u8>,
    used: BTreeMap<Code, usize>,
    minimum: u32,
}
impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for State {}
impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        self.priority
            .total_cmp(&other.priority)
            .then(self.sequence.len().cmp(&other.sequence.len()))
            .then(self.sequence.cmp(&other.sequence))
    }
}

impl Graph {
    /// Keep full main-k flanks even while traversing a smaller-k graph. First
    /// retain a nearby exit independently of coverage, then explore copy depth.
    /// No exhausted/copy-bounded domain can establish a unique structure.
    pub fn coverage_paths(
        &self,
        left: &[u8],
        right: &[u8],
        max_span: usize,
        max_paths: usize,
        minimum: u32,
        depth: &Depth<'_>,
    ) -> (Vec<(Vec<u8>, u32)>, Value) {
        assert!(left.len() >= self.k && right.len() >= self.k);
        let mask = (1u128 << (2 * (self.k - 1))) - 1;
        let start = code(&left[left.len() - self.k + 1..]).unwrap();
        let target = code(&right[right.len() - self.k + 1..]).unwrap();
        let support = left
            .windows(self.k)
            .map(|s| self.support(s))
            .min()
            .unwrap_or(0);
        if support < minimum {
            return (
                vec![],
                json!({"status":"unsupported_complete_left_flank","search_complete":false}),
            );
        }
        let count = |edge| {
            self.counts
                .get(&kmer::canonical(edge, self.k))
                .copied()
                .unwrap_or(0)
        };
        let mut distance = HashMap::from([(target, 0usize)]);
        let mut queue = VecDeque::from([target]);
        let mut reachability_complete = true;
        while let Some(node) = queue.pop_front() {
            if distance.len() >= MAX_REACHABLE {
                reachability_complete = false;
                break;
            }
            for base in 0..4u128 {
                let previous = (base << (2 * (self.k - 2))) | (node >> 2);
                let edge = (previous << 2) | (node & 3);
                if count(edge) >= minimum && !distance.contains_key(&previous) {
                    distance.insert(previous, distance[&node] + 1);
                    queue.push_back(previous);
                }
            }
        }
        // A truncated reverse search is not evidence that omitted nodes cannot
        // reach the target. Keep the bounded domain explicit and do not prune it.
        let reachable = |n: Code| !reachability_complete || distance.contains_key(&n);
        let initial_used = super::super::coverage::occurrences(left, false, self.k);
        let initial_cost = initial_used
            .iter()
            .map(|(&c, &n)| depth.increment(c, 0, n))
            .sum();
        let mut found = BTreeMap::<Vec<u8>, u32>::new();
        let mut schedules = Vec::new();
        for guided in [false, true] {
            let initial = State {
                priority: if guided {
                    initial_cost
                } else {
                    left.len() as f64
                },
                cost: initial_cost,
                node: start,
                sequence: left.to_vec(),
                used: initial_used.clone(),
                minimum: support,
            };
            let mut pending = BTreeSet::from([initial]);
            let (mut states, mut completions, mut peak) = (0, 0, 1);
            let mut limits = BTreeSet::new();
            while let Some(state) = pending.pop_first() {
                if states >= MAX_STATES {
                    limits.insert("state_limit");
                    break;
                }
                states += 1;
                if state.sequence.len() > left.len() && state.sequence.ends_with(right) {
                    found.entry(state.sequence).or_insert(state.minimum);
                    completions += 1;
                    if !guided || completions >= max_paths {
                        if !pending.is_empty() {
                            limits.insert("completion_limit");
                        }
                        break;
                    }
                    continue;
                }
                for base in 0..4u128 {
                    let edge = (state.node << 2) | base;
                    let next = edge & mask;
                    let support = count(edge);
                    if support < minimum || !reachable(next) {
                        continue;
                    }
                    if state.sequence.len() + 1 + distance.get(&next).copied().unwrap_or(0)
                        > max_span
                    {
                        limits.insert("sequence_span_limit");
                        continue;
                    }
                    let c = kmer::canonical(edge, self.k);
                    let before = state.used.get(&c).copied().unwrap_or(0);
                    if guided && before + 1 > depth.exploration_bound(c) {
                        limits.insert("coverage_copy_bound");
                        continue;
                    }
                    let mut child = state.clone();
                    child.node = next;
                    child.sequence.push(b"ACGT"[base as usize]);
                    if child.sequence.ends_with(left) {
                        limits.insert("unique_left_flank_revisit");
                        continue;
                    }
                    child.cost += depth.increment(c, before, 1);
                    *child.used.entry(c).or_default() += 1;
                    child.minimum = child.minimum.min(support);
                    child.priority = if guided {
                        child.cost
                    } else {
                        child.sequence.len() as f64
                    };
                    pending.insert(child);
                    if pending.len() > MAX_PENDING {
                        pending.pop_last();
                        limits.insert("pending_limit");
                    }
                    peak = peak.max(pending.len());
                }
            }
            schedules.push(
                json!({"schedule":if guided{"coverage_residual"}else{"nearby_exit_seed"},
                "states":states,"completions":completions,"peak_pending":peak,"limits":limits}),
            );
        }
        let report = json!({"k":self.k,"minimum_pair_support":minimum,"candidates":found.len(),
            "status":if found.is_empty(){"no_layout_within_search_bounds"}else{"inferred_candidates"},
            "search_complete":false,"reverse_reachable_nodes":distance.len(),"reachability_complete":reachability_complete,
            "max_states_per_schedule":MAX_STATES,"max_pending":MAX_PENDING,"max_reachable":MAX_REACHABLE,
            "max_span":max_span,"max_guided_completions":max_paths,"schedules":schedules,
            "graph":"raw sample-pair-supported edges; strict graph cleanup is not used"});
        (found.into_iter().collect(), report)
    }
}
