//! A necessary condition for an alignment within the existing edit budget.
//!
//! An edit can disrupt at most one of e+1 disjoint query blocks. Consequently
//! at least one whole block must occur exactly in any target accepting <=e
//! edits, including indels and query-end insertions. This only rejects windows
//! with no acceptable alignment; it neither proposes nor limits placements.

pub(super) struct ExactBlocks {
    length: usize,
    mask: u32,
    words: Vec<u32>,
}

fn code(base: u8) -> Option<u32> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

impl ExactBlocks {
    pub(super) fn new(query: &[u8], edits: usize) -> Option<Self> {
        let count = edits.checked_add(1)?;
        let length = query.len() / count;
        // Unknown query symbols must not erase one of the required blocks.
        // Wider or empty words also fall back to the unchanged affine search.
        if !(1..=15).contains(&length) || query.iter().any(|&b| code(b).is_none()) {
            return None;
        }
        let mut words: Vec<_> = query
            .chunks_exact(length)
            .take(count)
            .map(|block| {
                block
                    .iter()
                    .fold(0, |word, &b| (word << 2) | code(b).unwrap())
            })
            .collect();
        words.sort_unstable();
        words.dedup();
        Some(Self {
            length,
            mask: (1u32 << (2 * length)) - 1,
            words,
        })
    }

    pub(super) fn possible(&self, target: &[u8]) -> bool {
        let (mut word, mut valid) = (0, 0);
        for &base in target {
            if let Some(base) = code(base) {
                word = ((word << 2) | base) & self.mask;
                valid += 1;
                if valid >= self.length && self.words.binary_search(&word).is_ok() {
                    return true;
                }
            } else {
                word = 0;
                valid = 0;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn minimum_edits(query: &[u8], target: &[u8]) -> usize {
        let mut row = vec![0usize; target.len() + 1];
        for (i, &a) in query.iter().enumerate() {
            let mut diagonal = row[0];
            row[0] = i + 1;
            for (j, &b) in target.iter().enumerate() {
                let above = row[j + 1];
                row[j + 1] = (diagonal + usize::from(a != b))
                    .min(above + 1)
                    .min(row[j] + 1);
                diagonal = above;
            }
        }
        *row.iter().min().unwrap()
    }

    #[test]
    fn block_rejection_implies_no_edit_bounded_full_query_alignment() {
        let mut state = 0x32157983u64;
        let mut random = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state as usize
        };
        let mut rejected = 0;
        for case in 0..800 {
            let length = 16 + random() % 160;
            let query: Vec<_> = (0..length).map(|_| b"ACGT"[random() % 4]).collect();
            let edits = super::super::edit_budget(length);
            let filter = ExactBlocks::new(&query, edits).unwrap();
            let mut target = if case % 3 == 0 {
                (0..(length + 50)).map(|_| b"ACGTN"[random() % 5]).collect()
            } else {
                query.clone()
            };
            // Include leading/trailing edits, indels inside and between blocks,
            // and exactly-at-budget explanations rather than only easy matches.
            for _ in 0..(edits + usize::from(case % 5 == 0)) {
                let at = random() % target.len();
                match random() % 3 {
                    0 => target[at] = b"ACGTN"[random() % 5],
                    1 => {
                        target.remove(at);
                    }
                    _ => target.insert(at, b"ACGTN"[random() % 5]),
                }
            }
            if case % 2 == 0 {
                target.splice(..0, [b'N', b'T', b'A']);
                target.extend_from_slice(b"GNT");
            }
            if !filter.possible(&target) {
                rejected += 1;
                assert!(minimum_edits(&query, &target) > edits);
            }
        }
        assert!(rejected > 100);
        assert!(ExactBlocks::new(b"ACGTNACGT", 3).is_none());
        assert!(ExactBlocks::new(b"ACG", 3).is_none());
        assert!(ExactBlocks::new(b"", 3).is_none());
        assert!(!ExactBlocks::new(b"AAAAAAAA", 1)
            .unwrap()
            .possible(b"AAANAAA"));
    }
}
