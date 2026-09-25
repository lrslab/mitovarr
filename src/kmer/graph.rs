//! Exact graph k-mers through 127 bases. The short seed indexes keep their
//! existing u128 representation; only assembly graph keys use two words.

#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Code {
    high: u128,
    low: u128,
}

impl Code {
    fn shift_right(self, bits: usize) -> Self {
        match bits {
            0 => self,
            1..=127 => Self {
                high: self.high >> bits,
                low: (self.low >> bits) | (self.high << (128 - bits)),
            },
            128..=255 => Self {
                high: 0,
                low: self.high >> (bits - 128),
            },
            _ => Self::default(),
        }
    }

    pub fn prefix(self) -> Self {
        self.shift_right(2)
    }

    pub fn suffix(self, bases: usize) -> Self {
        assert!(bases <= 127);
        if bases < 64 {
            Self {
                high: 0,
                low: self.low & ((1u128 << (2 * bases)) - 1),
            }
        } else {
            Self {
                high: self.high & ((1u128 << (2 * (bases - 64))) - 1),
                low: self.low,
            }
        }
    }

    pub fn last_base(self) -> u8 {
        b"ACGT"[(self.low & 3) as usize]
    }
}

pub fn reverse_code(code: Code, k: usize) -> Code {
    assert!((1..=127).contains(&k));
    if k <= 63 {
        return Code {
            high: 0,
            low: super::reverse_code(code.low, k),
        };
    }
    let reverse_word = |word: u128| {
        let reversed = word.reverse_bits();
        const EVEN_BITS: u128 = u128::MAX / 3;
        !(((reversed >> 1) & EVEN_BITS) | ((reversed & EVEN_BITS) << 1))
    };
    Code {
        high: reverse_word(code.low),
        low: reverse_word(code.high),
    }
    .shift_right(256 - 2 * k)
}

pub fn canonical(code: Code, k: usize) -> Code {
    code.min(reverse_code(code, k))
}

pub fn decode(mut code: Code, k: usize) -> Vec<u8> {
    assert!(k <= 127);
    let mut sequence = vec![b'A'; k];
    for base in sequence.iter_mut().rev() {
        *base = code.last_base();
        code = code.prefix();
    }
    sequence
}

pub fn visit(
    sequence: &[u8],
    quality: Option<&[u8]>,
    k: usize,
    min_baseq: u8,
    mut f: impl FnMut(usize, Code),
) {
    assert!((1..=127).contains(&k));
    if k <= 63 {
        super::visit(sequence, quality, k, min_baseq, |p, low| {
            f(p, Code { high: 0, low })
        });
        return;
    }
    let high_mask = (1u128 << (2 * (k - 64))) - 1;
    let mut code = Code::default();
    let mut valid = 0;
    for (i, base) in sequence.iter().enumerate() {
        let value = match base.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        };
        if let Some(value) =
            value.filter(|_| quality.is_none_or(|q| q[i].saturating_sub(33) >= min_baseq))
        {
            code.high = ((code.high << 2) | (code.low >> 126)) & high_mask;
            code.low = (code.low << 2) | value;
            valid += 1;
            if valid >= k {
                f(i + 1 - k, code);
            }
        } else {
            code = Code::default();
            valid = 0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(length: usize) -> Vec<u8> {
        let mut state = 912871u64;
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 32) as usize & 3]
            })
            .collect()
    }

    #[test]
    fn every_supported_width_preserves_all_bases_and_graph_nodes() {
        let sequence = dna(290);
        for k in 1..=127 {
            visit(&sequence, None, k, 0, |position, code| {
                let window = &sequence[position..position + k];
                assert_eq!(decode(code, k), window, "k={k}");
                let reverse = reverse_code(code, k);
                assert_eq!(decode(reverse, k), crate::kmer::reverse_complement(window));
                assert_eq!(reverse_code(reverse, k), code);
                assert_eq!(canonical(reverse, k), canonical(code, k));
                assert_eq!(decode(code.prefix(), k - 1), window[..k - 1]);
                assert_eq!(decode(code.suffix(k - 1), k - 1), window[1..]);
                if k <= 63 {
                    assert_eq!(code.high, 0);
                    assert_eq!(decode(code, k), crate::kmer::decode(code.low, k));
                }
            });
        }
    }

    #[test]
    fn high_word_changes_do_not_collide_and_bad_bases_reset_both_words() {
        for k in [65, 95, 127] {
            let a = dna(k);
            let mut b = a.clone();
            b[0] = if b[0] == b'A' { b'C' } else { b'A' };
            let mut codes = Vec::new();
            for s in [&a, &b] {
                visit(s, None, k, 0, |_, c| codes.push(c));
            }
            assert_eq!(codes[0].low, codes[1].low);
            assert_ne!(codes[0], codes[1]);
            let mut sequence = dna(3 * k);
            sequence[k] = b'N';
            let mut quality = vec![b'I'; sequence.len()];
            quality[2 * k] = b'!';
            let mut positions = Vec::new();
            visit(&sequence, Some(&quality), k, 20, |p, c| {
                positions.push(p);
                assert_eq!(decode(c, k), sequence[p..p + k]);
            });
            assert_eq!(positions, vec![0]);
        }
    }
}
