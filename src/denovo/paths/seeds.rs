//! Sample-only anchor contexts. Short duplicated graph fragments must not make
//! an otherwise unique long segment disappear from the path decomposition.
use super::*;

pub(super) struct Seeds {
    counts: BTreeMap<usize, HashMap<kmer::Code, usize>>,
}

impl Seeds {
    pub fn new(graph: &PathGraph, maximum: usize) -> Self {
        let counts = [31, 47, 63]
            .into_iter()
            .filter(|&k| k <= maximum)
            .map(|k| {
                let mut counts = HashMap::new();
                for segment in &graph.segments {
                    kmer::visit(&segment.sequence, None, k, 0, |_, code| {
                        *counts.entry(kmer::canonical(code, k)).or_default() += 1;
                    });
                }
                (k, counts)
            })
            .collect();
        Self { counts }
    }

    pub fn occurrences(&self, k: usize) -> &HashMap<kmer::Code, usize> {
        &self.counts[&k]
    }

    fn unique(&self, sequence: &[u8], k: usize) -> Vec<(usize, kmer::Code)> {
        let mut positions = Vec::new();
        kmer::visit(sequence, None, k, 0, |offset, code| {
            // PathGraph contains both orientations, including repeated copies.
            if self.counts[&k].get(&kmer::canonical(code, k)) == Some(&2) {
                positions.push((offset, code));
            }
        });
        positions
    }

    pub fn anchor_length(&self, sequence: &[u8], minimum: usize) -> Option<usize> {
        self.counts
            .keys()
            .copied()
            .find(|&k| self.unique(sequence, k).len() >= minimum)
    }

    /// Keep usable 31-mers for short/noisy reads. Add one longer unique context
    /// only when the terminal region has too few unique 31-mer positions.
    /// Selection uses graph sequence alone, before any candidate read scoring.
    pub fn recruitment(
        &self,
        sequence: &[u8],
        minimum: usize,
    ) -> Vec<(usize, Vec<(usize, kmer::Code)>)> {
        let short = self.unique(sequence, 31);
        let mut result = Vec::new();
        let mut best_count = short.len();
        if !short.is_empty() {
            result.push((31, short));
        }
        if best_count < minimum {
            let mut fallback = None;
            for &k in self.counts.keys().filter(|&&k| k > 31) {
                let positions = self.unique(sequence, k);
                if positions.len() > best_count {
                    best_count = positions.len();
                    fallback = Some((k, positions));
                }
                if best_count >= minimum {
                    break;
                }
            }
            result.extend(fallback);
        }
        result
    }
}
