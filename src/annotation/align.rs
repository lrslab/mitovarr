//! Seed-window Smith–Waterman with affine gaps. Both query ends remain explicit:
//! an internal high-scoring fragment is never relabelled as a full gene.
use std::collections::{BTreeMap, HashMap};

#[derive(Clone, Debug)]
pub(super) struct Hit {
    pub qstart: usize,
    pub qend: usize,
    pub start: usize,
    pub end: usize,
    pub score: i32,
    pub matches: usize,
    pub ops: Vec<u8>,
}
impl Hit {
    pub fn identity(&self) -> f64 {
        self.matches as f64 / self.ops.len().max(1) as f64
    }
    pub fn gap_bases(&self) -> usize {
        self.ops.iter().filter(|&&x| x == b'I' || x == b'D').count()
    }
    pub fn non_triplet_gaps(&self) -> usize {
        let mut count = 0;
        let mut i = 0;
        while i < self.ops.len() {
            let op = self.ops[i];
            let mut end = i + 1;
            while end < self.ops.len() && self.ops[end] == op {
                end += 1;
            }
            if matches!(op, b'I' | b'D') && !(end - i).is_multiple_of(3) {
                count += 1;
            }
            i = end;
        }
        count
    }
    pub fn terminal_supported(&self, five_prime: bool) -> bool {
        let n = self.ops.len().min(8);
        let end = if five_prime {
            &self.ops[..n]
        } else {
            &self.ops[self.ops.len() - n..]
        };
        n >= 6
            && end.iter().filter(|&&x| x == b'=').count() >= n - 1
            && end.iter().all(|&x| x == b'=' || x == b'X')
    }
}
fn code(seq: &[u8]) -> Option<u32> {
    let mut x = 0;
    for &b in seq {
        x = (x << 2)
            | match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
    }
    Some(x)
}

pub(super) struct Search {
    index: HashMap<u32, Vec<usize>>,
    k: usize,
}
impl Search {
    pub fn new(target: &[u8], k: usize) -> Self {
        let mut index: HashMap<u32, Vec<usize>> = HashMap::new();
        for (i, word) in target.windows(k).enumerate() {
            if let Some(c) = code(word) {
                index.entry(c).or_default().push(i);
            }
        }
        Self { index, k }
    }
    pub fn hits(&self, query: &[u8], target: &[u8]) -> (Vec<Hit>, bool) {
        // Count distinct query seed positions and their span, not just overlapping
        // kmers from one short AT-rich match. Such a match is not a gene anchor.
        let mut diagonals: BTreeMap<isize, (usize, usize, usize, usize)> = BTreeMap::new();
        for (i, word) in query.windows(self.k).enumerate() {
            if let Some(positions) = code(word).and_then(|c| self.index.get(&c)) {
                // A ubiquitous low-complexity word cannot locate a gene. Other seeds
                // in the same query still nominate its window.
                if positions.len() > 64 {
                    continue;
                }
                for &p in positions {
                    let entry = diagonals
                        .entry(p as isize - i as isize)
                        .or_insert((0, i, i, 0));
                    entry.3 += if entry.0 == 0 {
                        self.k
                    } else {
                        self.k.min(i - entry.2)
                    };
                    entry.0 += 1;
                    entry.2 = i;
                }
            }
        }
        let min_span = query.len().min(24.max(query.len() / 20));
        let min_seed_bases = query.len().min((3 * self.k).max(min_span / 2));
        let mut ranked: Vec<_> = diagonals
            .into_iter()
            .filter(|(_, (count, first, last, bases))| {
                *count >= 3 && last - first + self.k >= min_span && *bases >= min_seed_bases
            })
            .map(|(d, (count, _, _, _))| (d, count))
            .collect();
        ranked.sort_by_key(|&(d, n)| (std::cmp::Reverse(n), d));
        let mut chosen: Vec<isize> = Vec::new();
        let mut limited = false;
        for (diagonal, _) in ranked {
            if chosen.iter().any(|&d| (d - diagonal).abs() <= 32) {
                continue;
            }
            if chosen.len() == 64 {
                limited = true;
                break;
            }
            chosen.push(diagonal);
        }
        let mut hits = Vec::new();
        let min_target = if self.k == 7 {
            query.len() * 4 / 5
        } else {
            (query.len() / 5).max(80).min(query.len())
        };
        for d in chosen {
            let begin = (d - 48).max(0) as usize;
            let end = (d + query.len() as isize + 48).max(0) as usize;
            let end = end.min(target.len());
            if begin >= end {
                continue;
            }
            let mut remaining = vec![(begin, end)];
            while let Some((left, right)) = remaining.pop() {
                if right - left < min_target {
                    continue;
                }
                if let Some(mut hit) = local(query, &target[left..right]) {
                    hit.start += left;
                    hit.end += left;
                    // A strong main gene must not hide an adjacent shorter copy
                    // in the same window. Search disjoint residual sequence;
                    // never align across a masked/artificial gap.
                    remaining.push((left, hit.start));
                    remaining.push((hit.end, right));
                    hits.push(hit);
                }
            }
        }
        (hits, limited)
    }
}

fn local(query: &[u8], target: &[u8]) -> Option<Hit> {
    if let Some(start) = target
        .windows(query.len())
        .position(|window| window == query)
    {
        return Some(Hit {
            qstart: 0,
            qend: query.len(),
            start,
            end: start + query.len(),
            score: query.len() as i32 * 3,
            matches: query.len(),
            ops: vec![b'='; query.len()],
        });
    }
    let width = target.len() + 1;
    let mut trace = vec![0u8; (query.len() + 1) * width];
    let mut previous = vec![0i32; width];
    let mut current = vec![0i32; width];
    let mut vertical = vec![i32::MIN / 4; width];
    let mut best = (0, 0, 0);
    for i in 1..=query.len() {
        let mut horizontal = i32::MIN / 4;
        current[0] = 0;
        for j in 1..width {
            let v_extend = vertical[j] - 1;
            let v_open = previous[j] - 8;
            vertical[j] = v_extend.max(v_open);
            let h_extend = horizontal - 1;
            let h_open = current[j - 1] - 8;
            horizontal = h_extend.max(h_open);
            let is_match = query[i - 1] == target[j - 1] && b"ACGT".contains(&query[i - 1]);
            let diagonal = previous[j - 1] + if is_match { 3 } else { -4 };
            let mut score = 0;
            let mut state = 0;
            for (value, s) in [(diagonal, 1), (vertical[j], 2), (horizontal, 3)] {
                if value > score {
                    score = value;
                    state = s;
                }
            }
            current[j] = score;
            trace[i * width + j] = state
                | if v_extend > v_open { 4 } else { 0 }
                | if h_extend > h_open { 8 } else { 0 };
            if score > best.0 {
                best = (score, i, j);
            }
        }
        std::mem::swap(&mut previous, &mut current);
    }
    if best.0 == 0 {
        return None;
    }
    let (mut i, mut j) = (best.1, best.2);
    let mut ops = Vec::new();
    let mut state = 0;
    loop {
        let t = trace[i * width + j];
        if state == 0 {
            state = t & 3;
        }
        match state {
            0 => break,
            1 => {
                ops.push(
                    if query[i - 1] == target[j - 1] && b"ACGT".contains(&query[i - 1]) {
                        b'='
                    } else {
                        b'X'
                    },
                );
                i -= 1;
                j -= 1;
                state = 0;
            }
            2 => {
                ops.push(b'D');
                i -= 1;
                state = if t & 4 != 0 { 2 } else { 0 };
            }
            3 => {
                ops.push(b'I');
                j -= 1;
                state = if t & 8 != 0 { 3 } else { 0 };
            }
            _ => unreachable!(),
        }
    }
    ops.reverse();
    Some(Hit {
        qstart: i,
        qend: best.1,
        start: j,
        end: best.2,
        score: best.0,
        matches: ops.iter().filter(|&&x| x == b'=').count(),
        ops,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn traceback_retains_indel_and_partial_query() {
        let q = b"ACGTTGCACGATCGATCCGTAGCTAGCTAACG";
        let mut t = q.to_vec();
        t.insert(15, b'A');
        let h = local(q, &t).unwrap();
        assert_eq!((h.qstart, h.qend, h.start, h.end), (0, q.len(), 0, t.len()));
        assert_eq!(h.non_triplet_gaps(), 1);
        let h = local(q, &q[5..]).unwrap();
        assert_eq!(h.qstart, 5);
    }
}
