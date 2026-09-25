//! Seeded, full-query affine-gap alignment. Reference flanks are free; query ends are not.
use crate::{io::fastq::Record, kmer};
use std::collections::{BTreeSet, HashMap};

const SEED: usize = 17;
const BAND: usize = 12;
const MAX_STARTS: usize = 256;
const NEG: f64 = -1e100;
const GAP_OPEN: f64 = -8.0;
const GAP_EXTEND: f64 = -2.0;

#[derive(Clone, Debug)]
pub(crate) struct Hit {
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
    pub logp: f64,
    pub edits: usize,
}

#[derive(Default, Debug)]
pub(crate) struct Placements {
    pub hits: Vec<Hit>,
    pub exhausted: bool,
}

struct OrientedQuery {
    sequence: Vec<u8>,
    quality: Vec<u8>,
    seeds: Vec<(usize, kmer::Code)>,
    ideal: f64,
}

/// A cache belongs to one read, including its qualities. It lives for one
/// fragment comparison and is reused across candidates, never across reads.
pub(crate) struct Prepared {
    orientations: [OrientedQuery; 2],
    windows: [HashMap<Vec<u8>, Vec<Hit>>; 2],
    rescues: [HashMap<Vec<u8>, Vec<Hit>>; 2],
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
            kmer::visit(&sequence, Some(&quality), SEED, 10, |p, c| {
                seeds.push((p, c))
            });
            let ideal = sequence
                .iter()
                .zip(&quality)
                .map(|(&b, &q)| emission(b, b, q))
                .sum();
            OrientedQuery {
                sequence,
                quality,
                seeds,
                ideal,
            }
        });
        Self {
            orientations,
            windows: std::array::from_fn(|_| HashMap::new()),
            rescues: std::array::from_fn(|_| HashMap::new()),
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
    pub fn new(sequence: &[u8], circular: bool, max_read: usize) -> Self {
        let mut expanded = sequence.to_vec();
        if circular && !sequence.is_empty() {
            expanded.extend((0..max_read + BAND).map(|i| sequence[i % sequence.len()]));
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
            .or_insert_with(|| fitting_candidates(&query.sequence, &query.quality, target));
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
            let quality = &query.quality;
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
                    .or_insert_with(|| fitting_candidates(sequence, quality, target));
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
    fitting_candidates(query, quality, target)
        .into_iter()
        .next()
}

fn fitting_candidates(query: &[u8], quality: &[u8], target: &[u8]) -> Vec<Hit> {
    let empty = Cell {
        score: NEG,
        start: 0,
        edits: 0,
    };
    let mut m: Vec<_> = (0..=target.len())
        .map(|start| Cell {
            score: 0.0,
            start,
            edits: 0,
        })
        .collect();
    let mut ins = vec![empty; m.len()];
    let mut del = vec![empty; m.len()];
    let mut next_m = vec![empty; m.len()];
    let mut next_i = vec![empty; m.len()];
    let mut next_d = vec![empty; m.len()];
    let max_edits = (query.len() / 10).max(3);
    if query.len() > target.len() + max_edits {
        return Vec::new();
    }
    for (i, &base) in query.iter().enumerate() {
        next_m.fill(empty);
        next_i.fill(empty);
        next_d.fill(empty);
        next_i[0] = best(m[0].advance(GAP_OPEN, 1), ins[0].advance(GAP_EXTEND, 1));
        let match_score = emission(base, base, quality[i]);
        let mismatch_base = if base == b'A' { b'C' } else { b'A' };
        let mismatch_score = emission(base, mismatch_base, quality[i]);
        let lo = (i + 1).saturating_sub(max_edits).max(1);
        let hi = (target.len() + max_edits)
            .saturating_sub(query.len() - i - 1)
            .min(target.len());
        for j in lo..=hi {
            next_m[j] = best(best(m[j - 1], ins[j - 1]), del[j - 1]).advance(
                if !b"ACGT".contains(&target[j - 1]) {
                    0.25f64.ln()
                } else if base == target[j - 1] {
                    match_score
                } else {
                    mismatch_score
                },
                usize::from(base != target[j - 1]),
            );
            next_i[j] = best(m[j].advance(GAP_OPEN, 1), ins[j].advance(GAP_EXTEND, 1));
            next_d[j] = best(
                next_m[j - 1].advance(GAP_OPEN, 1),
                next_d[j - 1].advance(GAP_EXTEND, 1),
            );
        }
        std::mem::swap(&mut m, &mut next_m);
        std::mem::swap(&mut ins, &mut next_i);
        std::mem::swap(&mut del, &mut next_d);
    }
    let mut hits: Vec<_> = (1..=target.len())
        .filter_map(|end| {
            let cell = best(best(m[end], ins[end]), del[end]);
            (cell.score > NEG / 2.0 && end > cell.start && cell.edits <= max_edits).then_some(Hit {
                start: cell.start,
                end,
                reverse: false,
                logp: cell.score,
                edits: cell.edits,
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
pub(crate) fn trace_bases(
    read: &Record,
    hit: &Hit,
    target: &[u8],
    minq: u8,
) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
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
    if hit.edits == 0 {
        let positions: Vec<_> = quality
            .iter()
            .enumerate()
            .filter(|(_, q)| q.saturating_sub(33) >= minq)
            .map(|(i, _)| hit.start + i)
            .collect();
        return (positions.clone(), positions, Vec::new());
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
    while i > 0 || j > 0 {
        let from = trace[i * width + j][state] as usize;
        match state {
            0 => {
                if i == 0 || j == 0 {
                    break;
                }
                if quality[i - 1].saturating_sub(33) >= minq {
                    observed.push(hit.start + j - 1);
                    if query[i - 1] == reference[j - 1] {
                        matching.push(hit.start + j - 1);
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
        return (Vec::new(), Vec::new(), Vec::new());
    }
    (observed, matching, deletions)
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
