//! Two-bit kmers. Ambiguous or low-quality bases break the rolling window.

pub mod graph;

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
