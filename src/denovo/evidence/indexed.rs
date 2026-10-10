//! Exact pair aggregation. For a fixed first mate, partition second-mate
//! intervals by which mate supplies the outer start/end. Within each partition
//! the span is fixed by one indexed endpoint (or interval length). Range maxima
//! replace individual pair visits. Separate trace masks retain every structural
//! evidence predicate; no alignment, span or low-weight contribution is dropped.
use super::{crosses, observed, BaseTrace, Hit, JointEvidence, Library};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd)]
struct Mask {
    direct: u8,
    left: u8,
    right: u8,
    anchors: u8,
}

impl Mask {
    fn new(
        hit: &Hit,
        trace: &BaseTrace,
        shift: usize,
        length: usize,
        junctions: [usize; 2],
        contraction: Option<[usize; 2]>,
    ) -> Self {
        let shifted = Hit {
            start: hit.start + shift,
            end: hit.end + shift,
            ..hit.clone()
        };
        let mut mask = Self::default();
        let terminal = shifted.end >= shifted.start + 20;
        let left_observed = terminal && observed(trace, shifted.end - 20, shifted.end, shift);
        let right_observed = terminal && observed(trace, shifted.start, shifted.start + 20, shift);
        for (j, junction) in junctions.into_iter().enumerate() {
            for (copy, p) in [junction, junction + length].into_iter().enumerate() {
                mask.direct |= u8::from(crosses(&shifted, trace, p, shift)) << j;
                mask.left |= u8::from(left_observed && shifted.end <= p) << (2 * j + copy);
                mask.right |= u8::from(right_observed && shifted.start >= p) << (2 * j + copy);
            }
        }
        if let Some(anchors) = contraction {
            for (copy, offset) in [0, length].into_iter().enumerate() {
                for (i, p) in anchors.into_iter().enumerate() {
                    mask.anchors |= u8::from(observed(trace, p + offset, p + offset + 21, shift))
                        << (2 * copy + i);
                }
            }
        }
        mask
    }
}

#[derive(Clone, Copy)]
struct Entry {
    start: usize,
    end: usize,
    score: f64,
    ambiguous: usize,
}

#[derive(Clone, Copy)]
struct Node {
    low: f64,
    high: f64,
    ambiguous: usize,
}

impl Default for Node {
    fn default() -> Self {
        Self {
            low: f64::INFINITY,
            high: f64::NEG_INFINITY,
            ambiguous: 0,
        }
    }
}

#[derive(Clone, Copy)]
struct Best {
    score: f64,
    ambiguous: usize,
}

impl Best {
    fn retain(&mut self, score: f64, ambiguous: usize) {
        if score > self.score {
            *self = Self { score, ambiguous };
        } else if score == self.score {
            self.ambiguous = self.ambiguous.max(ambiguous);
        }
    }
}

struct RangeMax {
    coordinates: Vec<usize>,
    nodes: Vec<Node>,
    leaves: usize,
}

impl RangeMax {
    fn new(mut entries: Vec<(usize, Entry)>) -> Self {
        entries.sort_unstable_by_key(|e| e.0);
        let leaves = entries.len().next_power_of_two();
        let mut nodes = vec![Node::default(); 2 * leaves];
        for (i, (_, e)) in entries.iter().enumerate() {
            nodes[leaves + i] = Node {
                low: e.score,
                high: e.score,
                ambiguous: e.ambiguous,
            };
        }
        for i in (1..leaves).rev() {
            let (a, b) = (nodes[2 * i], nodes[2 * i + 1]);
            nodes[i] = Node {
                low: a.low.min(b.low),
                high: a.high.max(b.high),
                ambiguous: a.ambiguous.max(b.ambiguous),
            };
        }
        Self {
            coordinates: entries.into_iter().map(|e| e.0).collect(),
            nodes,
            leaves,
        }
    }

    // Apply the original floating-point expression before deciding ties. A
    // lower mate score can round to the same pair score, or hit the zero clamp;
    // its ambiguous contacts must then survive too.
    fn visit(&self, node: usize, x: f64, ideal: f64, best: &mut Best) {
        let n = self.nodes[node];
        if n.high == f64::NEG_INFINITY {
            return;
        }
        let high = (x + n.high - ideal).min(0.0);
        if high < best.score || (high == best.score && n.ambiguous <= best.ambiguous) {
            return;
        }
        if (x + n.low - ideal).min(0.0) == high {
            best.retain(high, n.ambiguous);
        } else {
            self.visit(2 * node, x, ideal, best);
            self.visit(2 * node + 1, x, ideal, best);
        }
    }

    fn query(&self, range: Range, x: f64, ideal: f64) -> Option<Best> {
        if range.lo > range.hi {
            return None;
        }
        let mut left = self.coordinates.partition_point(|&p| p < range.lo) + self.leaves;
        let mut right = self.coordinates.partition_point(|&p| p <= range.hi) + self.leaves;
        if left == right {
            return None;
        }
        let mut best = Best {
            score: f64::NEG_INFINITY,
            ambiguous: 0,
        };
        while left < right {
            if left % 2 == 1 {
                self.visit(left, x, ideal, &mut best);
                left += 1;
            }
            if right % 2 == 1 {
                right -= 1;
                self.visit(right, x, ideal, &mut best);
            }
            left /= 2;
            right /= 2;
        }
        Some(best)
    }
}

#[derive(Clone, Copy)]
struct Range {
    lo: usize,
    hi: usize,
}

impl Range {
    fn intersect(self, lo: usize, hi: usize) -> Self {
        Self {
            lo: self.lo.max(lo),
            hi: self.hi.min(hi),
        }
    }
}

struct Bucket {
    outer: usize,
    values: RangeMax,
}

fn buckets(
    entries: &[Entry],
    outer: fn(&Entry) -> usize,
    inner: fn(&Entry) -> usize,
) -> Vec<Bucket> {
    let mut groups = BTreeMap::<usize, Vec<(usize, Entry)>>::new();
    for e in entries {
        groups.entry(outer(e)).or_default().push((inner(e), *e));
    }
    groups
        .into_iter()
        .map(|(outer, entries)| Bucket {
            outer,
            values: RangeMax::new(entries),
        })
        .collect()
}

struct Group {
    reverse: bool,
    mask: Mask,
    starts: Vec<Bucket>,
    ends: Vec<Bucket>,
    lengths: Vec<Bucket>,
}

fn index(
    b: &[(Hit, BaseTrace)],
    shift: usize,
    length: usize,
    junctions: [usize; 2],
    contraction: Option<[usize; 2]>,
) -> Vec<Group> {
    let mut groups = BTreeMap::<(bool, Mask), Vec<Entry>>::new();
    for (h, t) in b {
        let mask = Mask::new(h, t, shift, length, junctions, contraction);
        groups.entry((h.reverse, mask)).or_default().push(Entry {
            start: h.start + shift,
            end: h.end + shift,
            score: h.logp,
            ambiguous: t.ambiguous_contacts,
        });
    }
    groups
        .into_iter()
        .map(|((reverse, mask), entries)| Group {
            reverse,
            mask,
            starts: buckets(&entries, |e| e.start, |e| e.end),
            ends: buckets(&entries, |e| e.end, |e| e.start),
            lengths: buckets(&entries, |e| e.end - e.start, |e| e.start),
        })
        .collect()
}

struct Summary<'a> {
    states: BTreeMap<usize, (f64, usize)>,
    model: Option<&'a Library>,
    length: usize,
    best: f64,
    direct: [f64; 2],
    brackets: [f64; 2],
    contraction: f64,
}

impl Summary<'_> {
    fn weight(&self, span: usize) -> Option<f64> {
        if self.length > 0 && span > self.length {
            return None;
        }
        let w = self
            .model
            .filter(|m| m.geometry_usable)
            .map_or(1.0, |m| m.probability(span));
        (w != 0.0).then_some(w)
    }

    fn accept(
        &mut self,
        span: usize,
        weight: f64,
        best: Best,
        x: Mask,
        y: Mask,
        ambiguous_x: usize,
    ) {
        let empirical = self.model.filter(|m| m.geometry_usable);
        let ambiguous = best.ambiguous + ambiguous_x;
        self.states
            .entry(if empirical.is_some() { span } else { 0 })
            .and_modify(|v| {
                if best.score > v.0 {
                    *v = (best.score, ambiguous);
                } else if best.score == v.0 {
                    v.1 = v.1.max(ambiguous);
                }
            })
            .or_insert((best.score, ambiguous));
        if best.score < 0.1f64.ln() {
            return;
        }
        let relative = best.score.exp() * empirical.map_or(1.0, |m| weight / m.peak_probability());
        self.best = self.best.max(relative);
        if relative < 0.1 {
            return;
        }
        let paired = empirical.is_some_and(|m| span <= m.upper);
        for j in 0..2 {
            if (x.direct | y.direct) & (1 << j) != 0 {
                self.direct[j] = self.direct[j].max(relative);
            }
            if paired && ((x.left & y.right) | (y.left & x.right)) & (3 << (2 * j)) != 0 {
                self.brackets[j] = self.brackets[j].max(relative);
            }
        }
        if (0..2).any(|copy| {
            let bits = 3 << (2 * copy);
            x.anchors & bits == bits
                || y.anchors & bits == bits
                || (paired && (x.anchors | y.anchors) & bits == bits)
        }) {
            self.contraction = self.contraction.max(relative);
        }
    }

    fn finish(self) -> JointEvidence {
        let empirical = self.model.filter(|m| m.geometry_usable);
        let span_states: Vec<_> = self.states.into_iter().map(|(span, (log_sequence, ambiguous))| {
            let sequence = log_sequence.exp();
            let geometry = empirical.map_or(1.0, |m| m.probability(span));
            serde_json::json!({"span":span,"sequence_weight":sequence,"geometry_weight":geometry,"contribution":sequence*geometry,"ambiguous_base_contacts":ambiguous})
        }).collect();
        let mass: f64 = span_states
            .iter()
            .map(|s| s["contribution"].as_f64().unwrap())
            .sum();
        let floor = (-12.0f64).exp();
        let supported = |p| p >= 0.1 && p >= 0.99 * self.best;
        JointEvidence {
            score: (mass.min(1.0) + floor).ln(),
            junctions: std::array::from_fn(|i| supported(self.direct[i].max(self.brackets[i]))),
            observed_junctions: self.direct.map(supported),
            mate_brackets: self.brackets.map(supported),
            contraction: supported(self.contraction),
            diagnostics: serde_json::json!({"span_states":span_states,"mass_before_floor":mass.min(1.0),"outlier_mass_fraction":floor/(mass.min(1.0)+floor),"score":(mass.min(1.0)+floor).ln()}),
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn joint(
    a: &[(Hit, BaseTrace)],
    b: &[(Hit, BaseTrace)],
    model: Option<&Library>,
    ideal: f64,
    length: usize,
    junctions: [usize; 2],
    contraction: Option<[usize; 2]>,
) -> JointEvidence {
    let mut summary = Summary {
        states: BTreeMap::new(),
        model,
        length,
        best: 0.0,
        direct: [0.0; 2],
        brackets: [0.0; 2],
        contraction: 0.0,
    };
    let original = index(b, 0, length, junctions, contraction);
    let shifted = if length > 0 {
        index(b, length, length, junctions, contraction)
    } else {
        vec![]
    };
    for (hit, trace) in a {
        for mode in 0..if length > 0 { 3 } else { 1 } {
            let shift_x = if mode == 1 { length } else { 0 };
            let (xs, xe) = (hit.start + shift_x, hit.end + shift_x);
            let mask = Mask::new(hit, trace, shift_x, length, junctions, contraction);
            for group in if mode == 2 { &shifted } else { &original } {
                let mut starts = Range {
                    lo: 0,
                    hi: usize::MAX,
                };
                let mut ends = starts;
                if mode == 1 {
                    starts.lo = hit.start + 1;
                }
                if mode == 2 {
                    starts.hi = hit.start + length;
                }
                if let Some(m) = model {
                    if hit.reverse == group.reverse {
                        continue;
                    }
                    if hit.reverse {
                        if m.inward {
                            let Some(hi) = xe.checked_sub(1) else {
                                continue;
                            };
                            starts.hi = starts.hi.min(hi);
                        } else {
                            starts.lo = starts.lo.max(xe);
                        }
                    } else if m.inward {
                        ends.lo = xs + 1;
                    } else {
                        ends.hi = xs;
                    }
                }
                if starts.lo > starts.hi || ends.lo > ends.hi {
                    continue;
                }
                let mut query = |bucket: &Bucket, range: Range, span: usize| {
                    if range.lo > range.hi {
                        return;
                    }
                    if let Some(w) = summary.weight(span) {
                        if let Some(best) = bucket.values.query(range, hit.logp, ideal) {
                            summary.accept(
                                span,
                                w,
                                best,
                                mask,
                                group.mask,
                                trace.ambiguous_contacts,
                            );
                        }
                    }
                };
                // x supplies the start, y supplies the end (including overlaps).
                for bucket in &group.ends {
                    let end = bucket.outer;
                    if end > xe && end >= ends.lo && end <= ends.hi {
                        query(bucket, starts.intersect(xs, usize::MAX), end - xs);
                    }
                }
                // y supplies the start, x supplies the end.
                for bucket in &group.starts {
                    let start = bucket.outer;
                    if start < xs && start >= starts.lo && start <= starts.hi {
                        query(bucket, ends.intersect(0, xe), xe - start);
                    }
                }
                for bucket in &group.lengths {
                    let span = bucket.outer;
                    let Some(end_hi) = ends.hi.checked_sub(span) else {
                        continue;
                    };
                    let range = starts.intersect(ends.lo.saturating_sub(span), end_hi);
                    // x contains y, including identical endpoint pairs.
                    if let Some(last) = xe.checked_sub(span) {
                        query(bucket, range.intersect(xs, last), xe - xs);
                    }
                    // y strictly contains x: the other three cases own ties.
                    if let Some(last) = xs.checked_sub(1) {
                        query(
                            bucket,
                            range.intersect((xe + 1).saturating_sub(span), last),
                            span,
                        );
                    }
                }
            }
        }
    }
    summary.finish()
}

#[cfg(test)]
mod tests;
