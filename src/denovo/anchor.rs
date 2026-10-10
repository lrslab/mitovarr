//! One observed common-flank word fixes mate, query offset and target identity
//! before alignment scores are inspected. All hypotheses use that same pivot.
use super::align::{self, Hit, Index, Placements, Prepared};
use crate::{io::fastq::Record, kmer};
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};

pub(crate) const K: usize = 21;
pub(crate) type Locations = HashMap<kmer::Code, (usize, kmer::Code)>;

pub(crate) fn locations(index: &Index, common: &HashSet<kmer::Code>) -> Locations {
    let length = index
        .circle_length()
        .map_or(index.sequence.len(), |n| n + K - 1);
    let mut result = HashMap::new();
    kmer::visit(
        &index.sequence[..length.min(index.sequence.len())],
        None,
        K,
        0,
        |p, c| {
            let key = kmer::canonical(c, K);
            if common.contains(&key) {
                result.insert(key, (p, c));
            }
        },
    );
    result
}

pub(crate) struct Pivot {
    pub mate: usize,
    pub code: kmer::Code,
    pub offset: usize,
    word: kmer::Code,
    words: Vec<(usize, kmer::Code)>,
}

impl Pivot {
    pub fn choose(reads: [&Record; 2], common: &HashSet<kmer::Code>, minq: u8) -> Option<Self> {
        for (mate, read) in reads.into_iter().enumerate() {
            let mut words = Vec::new();
            kmer::visit(&read.sequence, Some(&read.quality), K, minq, |p, c| {
                if common.contains(&kmer::canonical(c, K)) {
                    words.push((p, c));
                }
            });
            if let Some(&(offset, word)) = words.first() {
                // Only overlapping words describe this same anchor block.
                // A distant flank in a deletion-spanning read is independent
                // structural evidence; its changed distance is not pivot drift.
                words.retain(|&(p, _)| p < offset + K);
                return Some(Self {
                    mate,
                    code: kmer::canonical(word, K),
                    offset,
                    word,
                    words,
                });
            }
        }
        None
    }

    pub fn report(&self) -> Value {
        json!({"mate":self.mate+1,"canonical_21mer_code":self.code,"query_offset_0based":self.offset})
    }

    pub fn projection(&self, map: &Locations, read_length: usize) -> (usize, bool, usize) {
        let &(p, c) = &map[&self.code];
        let reverse = c != self.word;
        (
            p,
            reverse,
            if reverse {
                read_length - self.offset - K
            } else {
                self.offset
            },
        )
    }

    /// Reject incompatible projections from other observed common words in the
    /// chosen mate. Indels within the accepted full-query edit bound are allowed.
    pub fn consistent(&self, map: &Locations, read_length: usize, circle: Option<usize>) -> bool {
        let (p, reverse, q) = self.projection(map, read_length);
        let origin = p as isize - q as isize;
        let bound = 2 * (read_length / 10).max(3) as isize;
        self.words.iter().all(|&(offset, word)| {
            let &(p, c) = &map[&kmer::canonical(word, K)];
            if (c != word) != reverse {
                return false;
            }
            let q = if reverse {
                read_length - offset - K
            } else {
                offset
            };
            let mut delta = p as isize - q as isize - origin;
            if let Some(n) = circle {
                delta = delta.rem_euclid(n as isize);
                if delta > n as isize / 2 {
                    delta -= n as isize;
                }
            }
            delta.abs() <= bound
        })
    }

    fn accepts(&self, hit: &Hit, read: &Record, index: &Index, map: &Locations, minq: u8) -> bool {
        let (p, reverse, q) = self.projection(map, read.sequence.len());
        if hit.reverse != reverse {
            return false;
        }
        let trace = align::trace_alignment(read, hit, &index.sequence, minq);
        let n = index.circle_length().unwrap_or(0);
        [p, p + n].into_iter().any(|start| {
            trace
                .matching_pairs
                .iter()
                .filter(|&&(query, target)| {
                    (q..q + K).contains(&query) && target == start + query - q
                })
                .count()
                == K
        })
    }

    /// Replace heuristic hits with endpoint profiles optimized under the pivot.
    /// An empty, complete result proves absence only in the declared edit domain.
    pub fn constrain(
        &self,
        placements: &mut Placements,
        read: &Record,
        prepared: &mut Prepared,
        index: &Index,
        map: &Locations,
        minq: u8,
    ) -> usize {
        let removed = placements
            .hits
            .iter()
            .filter(|h| !self.accepts(h, read, index, map, minq))
            .count();
        let (p, reverse, q) = self.projection(map, read.sequence.len());
        *placements = index.profile_anchor(prepared, p, reverse, q, K);
        removed
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_shared_query_anchor_survives_indel_rotation_and_candidate_specific_decoys() {
        let mut seed = 773929u64;
        let genome: Vec<_> = (0..1000)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed % 4) as usize]
            })
            .collect();
        let mut query = genome[850..950].to_vec();
        query.insert(10, b'G');
        let read = Record {
            header: "@test/1".into(),
            quality: vec![b'I'; query.len()],
            sequence: query.clone(),
            plus: "+".into(),
        };
        let mut common = HashSet::new();
        kmer::visit(&genome[890..911], None, K, 0, |_, c| {
            common.insert(kmer::canonical(c, K));
        });
        let pivot = Pivot::choose([&read, &read], &common, 0).unwrap();
        assert_eq!(pivot.offset, 41);
        let mut decoy = query.clone();
        decoy[50] = if decoy[50] == b'A' { b'C' } else { b'A' };
        let mut changed = genome.clone();
        changed[400..501].copy_from_slice(&decoy);
        let mut scores = Vec::new();
        let mut removed = 0;
        for sequence in [genome, changed] {
            for reverse in [false, true] {
                let oriented = if reverse {
                    kmer::reverse_complement(&sequence)
                } else {
                    sequence.clone()
                };
                for offset in [0, 400, 900] {
                    let rotated = [&oriented[offset..], &oriented[..offset]].concat();
                    let index = Index::new(&rotated, true, read.sequence.len());
                    let map = locations(&index, &common);
                    assert!(pivot.consistent(&map, read.sequence.len(), index.circle_length()));
                    let mut prepared = Prepared::new(&read);
                    let mut hits = index.align_prepared(&mut prepared);
                    removed += pivot.constrain(&mut hits, &read, &mut prepared, &index, &map, 0);
                    assert!(!hits.exhausted && !hits.hits.is_empty());
                    assert!(hits
                        .hits
                        .iter()
                        .all(|h| pivot.accepts(h, &read, &index, &map, 0)));
                    scores.push(
                        hits.hits
                            .iter()
                            .map(|h| h.logp)
                            .max_by(f64::total_cmp)
                            .unwrap(),
                    );
                }
            }
        }
        assert!(
            removed > 0,
            "candidate-specific decoy must really be exercised"
        );
        assert!(
            scores.iter().all(|s| (s - scores[0]).abs() < 1e-8),
            "{scores:?}"
        );
    }
}
