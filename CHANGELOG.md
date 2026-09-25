# Changelog

## 0.1.0 - Unreleased

### Package and interfaces

- Exclude local AI-assistant files, TODOs, maintainer notes and exploratory-test directories from Git; limit the source package to runtime code, active regression tests and user documentation.
- Rename the package, library, executable and launcher to `mitovarr`. Preserve existing `mitovar-rust.*` JSON schema identifiers and gene panels; use `mitovarr` as the GFF3 source.
- Provide an independently buildable Rust package, English and Chinese usage guides, CLI documentation, method descriptions and editable workflow Figure 1. Keep experiment records and analysis outputs outside the repository and source package.
- Remove unused legacy pileup, variant-refresh and repeat-summary helpers and redundant tests; retain coverage through the production fragment accumulator and portable regression fixtures.

### Assembly and consensus

- Prevent local contractions from being published on separate repeat-internal junction evidence alone; require three named fragments linking unique unchanged flanks and report this evidence explicitly.

- Add native CaenDR extraction through indexed HTTP Range requests, local mate recovery and paired FASTQ/BAM output, with a portable launcher.
- Keep reference-blind FASTQ, seed-recruited FASTQ, BAM-baited de novo and reference-guided indexed-BAM workflows separate. Runtime processing uses in-process Rust implementations.
- Build sample de Bruijn graphs with odd k values from 15 through 127, local tip/bubble cleanup, weak-edge filtering, read-supported repeat traversal and smaller-k end connections.
- Compare complete candidate paths with full-query affine-gap alignment and empirical paired-fragment geometry. Retain alternatives, unsearched repeat tails and alignment limits; do not treat graph circularity as validated closure.
- Write a primary assembly and localized uncertainty alongside candidate FASTAs, GFA graphs and read-support reports. Publish a final assembly only when its evidence requirements pass.
- Use fragment-aware mapping consensus with separate SNP-majority and indel thresholds, base-vote reports and explicit repeat ambiguity. Offer deterministic complete-fragment subsampling with matching FASTQ selection for local assembly.
- Add optional local NCR assembly to mapping: discover suspect windows from BAM, generate sample-read paths and compare whole-contig alternatives with both-junction support. Keep incomplete searches unresolved and export candidate sequences, graphs, structural edits and coordinate mappings.
- Correct single-SNP spanning intervals, short-terminal anchor eligibility, either-mate anchoring and matching-base coverage without counting unsequenced mate gaps as depth.

### Annotation

- Require a supported CDS 3′ boundary before interpreting terminal stops or T/TA polyadenylation; retain stop/gap evidence for partial genes and avoid treating a clipped codon alone as a frameshift.

- Add homology-assisted CDS, tRNA and rRNA prediction on both strands with bundled three-species panels, mitochondrial genetic code 5 and explicit detection/completeness status.
- Rotate eligible circular input hypotheses to a unique, complete tRNA-Pro anchor, preserving sample bases and assembly status while transforming feature coordinates and strands.
- Emit positional NCR1–4 labels, generic gaps and locus states, retaining short gaps, unknown sequence, missing boundaries and separate linear ends.
- Preserve topology during repeated annotation and project assembly uncertainty into annotated coordinates, including short, zero-length and overlapping NCR boundaries.
- Export GFF3, feature tables, nucleotide/protein FASTA, NCR FASTA and JSON reports. RNA structure, general repeat copy-number inference and automated cross-sample NCR SV classification remain outside the implemented scope.
