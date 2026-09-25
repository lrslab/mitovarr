//! NCBI invertebrate mitochondrial translation table 5.
use serde::Serialize;
const AMINO: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
fn amino(codon: &[u8]) -> u8 {
    let mut index = 0;
    for &b in codon {
        index = index * 4
            + match b {
                b'T' => 0,
                b'C' => 1,
                b'A' => 2,
                b'G' => 3,
                _ => return b'X',
            };
    }
    if codon.len() == 3 {
        AMINO[index]
    } else {
        b'X'
    }
}
pub(super) fn is_start(codon: &[u8]) -> bool {
    matches!(codon, b"TTG" | b"ATT" | b"ATC" | b"ATA" | b"ATG" | b"GTG")
}
pub(super) fn is_stop(codon: &[u8]) -> bool {
    matches!(codon, b"TAA" | b"TAG")
}
#[derive(Clone, Debug, Serialize)]
pub struct CdsCheck {
    pub translation_table: u8,
    pub phase: usize,
    pub start_codon: String,
    pub valid_start: bool,
    pub stop: String,
    pub internal_stop_codons: Vec<usize>,
    pub ambiguous_codons: usize,
    pub non_triplet_alignment_gaps: usize,
    pub frame_valid: bool,
    pub protein: String,
}
pub(super) fn check(
    sequence: &[u8],
    phase: usize,
    complete_five: bool,
    complete_three: bool,
    non_triplet_gaps: usize,
) -> CdsCheck {
    let coding = &sequence[phase.min(sequence.len())..];
    let remainder = coding.len() % 3;
    let stop = if !complete_three {
        "partial_or_uncertain".into()
    } else if remainder == 0 && coding.len() >= 3 && is_stop(&coding[coding.len() - 3..]) {
        String::from_utf8_lossy(&coding[coding.len() - 3..]).into_owned()
    } else if remainder == 1 && coding.ends_with(b"T") {
        "T_requires_polyadenylation".into()
    } else if remainder == 2 && coding.ends_with(b"TA") {
        "TA_requires_polyadenylation".into()
    } else {
        "missing".into()
    };
    let mut protein: Vec<_> = coding.chunks_exact(3).map(amino).collect();
    if matches!(stop.as_str(), "TAA" | "TAG") {
        protein.pop();
    }
    let valid_start = complete_five && coding.len() >= 3 && is_start(&coding[..3]);
    if valid_start && !protein.is_empty() {
        protein[0] = b'M';
    }
    let internal_stop_codons = protein
        .iter()
        .enumerate()
        .filter_map(|(i, &b)| (b == b'*').then_some(i + 1))
        .collect();
    CdsCheck {
        translation_table: 5,
        phase,
        start_codon: String::from_utf8_lossy(&coding[..coding.len().min(3)]).into_owned(),
        valid_start,
        // A clipped 3′ end can stop within a codon without implying a frameshift.
        // This does not establish gene completeness or erase observed gap evidence.
        frame_valid: non_triplet_gaps == 0
            && (!complete_three || remainder == 0 || stop != "missing"),
        stop,
        internal_stop_codons,
        ambiguous_codons: protein.iter().filter(|&&x| x == b'X').count(),
        non_triplet_alignment_gaps: non_triplet_gaps,
        protein: String::from_utf8(protein).unwrap(),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn table_five_starts_stops_and_polyadenylation() {
        let c = check(b"TTGTGAAGAAGGTAA", 0, true, true, 0);
        assert_eq!(c.protein, "MWSS");
        assert!(c.valid_start);
        assert!(c.internal_stop_codons.is_empty());
        for tail in ["T", "TA"] {
            let c = check(format!("ATTGCT{tail}").as_bytes(), 0, true, true, 0);
            assert!(c.stop.ends_with("requires_polyadenylation"));
            assert!(c.frame_valid);
            assert_eq!(c.protein, "MA");
        }
        assert_eq!(
            check(b"ATGTAGGCTTAA", 0, true, true, 0).internal_stop_codons,
            vec![2]
        );
        assert_eq!(check(b"ATGNNNTAA", 0, true, true, 0).ambiguous_codons, 1);
    }

    #[test]
    fn uncertain_three_prime_end_keeps_stop_and_gap_evidence() {
        let c = check(b"ATGGCTTAA", 0, true, false, 0);
        assert_eq!(c.stop, "partial_or_uncertain");
        assert_eq!(c.protein, "MA*");
        assert_eq!(c.internal_stop_codons, vec![3]);
        assert!(!check(b"ATGGCTT", 0, true, false, 1).frame_valid);
    }
}
