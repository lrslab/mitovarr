//! Two-bit kmers. Ambiguous or low-quality bases break the rolling window.

pub mod graph;

/// Internal numeric-key tables. Iteration order must never choose a graph path.
/// This deterministic, non-cryptographic hasher mixes every word of a k-mer;
/// callers still sort graph edges and all externally ordered candidates.
#[derive(Default)]
pub(crate) struct KmerHasher(u64);

impl std::hash::Hasher for KmerHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        let mut words = bytes.chunks_exact(8);
        for word in &mut words {
            self.write_u64(u64::from_ne_bytes(word.try_into().expect("eight bytes")));
        }
        let mut tail = 0u64;
        for (i, &byte) in words.remainder().iter().enumerate() {
            tail |= u64::from(byte) << (8 * i);
        }
        if !words.remainder().is_empty() {
            self.write_u64(tail);
        }
    }

    fn write_u64(&mut self, word: u64) {
        self.0 = (self.0.rotate_left(5) ^ word).wrapping_mul(0x517c_c1b7_2722_0a95);
    }

    fn write_u128(&mut self, word: u128) {
        self.write_u64(word as u64);
        self.write_u64((word >> 64) as u64);
    }
}

pub(crate) type Map<V> =
    std::collections::HashMap<Code, V, std::hash::BuildHasherDefault<KmerHasher>>;

pub type Code = u128;

pub fn reverse_code(code: Code, k: usize) -> Code {
    const EVEN_BITS: Code = 0x55555555555555555555555555555555;
    let reversed = code.reverse_bits();
    let pairs = ((reversed >> 1) & EVEN_BITS) | ((reversed & EVEN_BITS) << 1);
    !pairs >> (128 - 2 * k)
}

pub fn canonical(code: Code, k: usize) -> Code {
    code.min(reverse_code(code, k))
}

pub fn decode(mut code: Code, k: usize) -> Vec<u8> {
    let mut sequence = vec![b'A'; k];
    for base in sequence.iter_mut().rev() {
        *base = b"ACGT"[(code & 3) as usize];
        code >>= 2;
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
    assert!((1..=63).contains(&k));
    let mask = (1u128 << (2 * k)) - 1;
    let (mut code, mut valid) = (0, 0);
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
            code = ((code << 2) | value) & mask;
            valid += 1;
            if valid >= k {
                f(i + 1 - k, code);
            }
        } else {
            code = 0;
            valid = 0;
        }
    }
}

pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

pub(crate) type Set = std::collections::HashSet<Code, std::hash::BuildHasherDefault<KmerHasher>>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wide_kmers_preserve_all_bases_and_reverse_orientation() {
        let sequence = b"ACTGTTCCAGTACCACTGATCGTACTGACGATCCTGATCAGGCTTATCGGATTCGAGTTACGGA";
        for k in 1..=63 {
            visit(sequence, None, k, 0, |position, code| {
                let window = &sequence[position..position + k];
                assert_eq!(decode(code, k), window);
                assert_eq!(decode(reverse_code(code, k), k), reverse_complement(window));
                assert_eq!(canonical(code, k), canonical(reverse_code(code, k), k));
            });
        }
    }
}
