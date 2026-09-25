# Reference-guided local NCR assembly

This optional stage follows ordinary
BAM pileup consensus. The reference locates suspect windows and flanks, while
paired FASTQ supplies every graph edge and every inserted base. This is a
reference-guided workflow, separate from reference-blind and BAM-baited de novo.

When graph enumeration reaches a resource limit, candidate FASTA/GFA and path
support are retained, but the ineligible subset is not scored. The window has
`scoring_status=not_run_graph_search_limit`; fragment counts, scores, junction
counts and best-candidate fields are JSON `null` (unmeasured). An exhaustive
bounded search uses `scoring_status=evaluated`. Skipping an ineligible subset
does not make an incomplete search resolved or imply that the sample lacks
supporting reads.

1. Collect named-fragment breakpoint evidence while reading BAM. Combine soft
   clips, large CIGAR indels, discordant mate endpoints and unresolved consensus
   intervals. Handle intervals crossing the mitochondrial origin.
2. Find supported flanks and enumerate bounded sample graph paths between them.
   Include the unmodified pileup candidate in each comparison. Keep graph and
   candidate sequence outputs for inspection; report resource limits.
3. Reuse the native full-query affine-gap aligner and empirical fragment model.
   Compare the same fragments against whole-contig alternatives, retaining the
   original homologous locus in every alternative. Shared repeat reads must not
   multiply their contribution simply because a candidate contains more copies.
   Require read-supported connections on both sides of a replacement. A contraction
   of at least `--local-k` bases additionally requires three named fragments
   that each observe both unique unchanged 21-bp flanks in quality-filtered
   sequence, with a compatible joint placement (one spanning read or usable
   paired-fragment geometry). The nearest qualifying flanks must be unique in
   every whole-contig candidate. Separate repeat-internal junction counts are
   insufficient to establish contraction length, particularly against an
   N-masked baseline. Missing or unlinked unique flanks leave the window
   unresolved as `contraction_unique_flanks_not_linked`. Scored candidates add
   `requires_contraction_flanks`, `contraction_flank_starts_0based` (in the rotated
   whole-contig candidate), and `contraction_flank_fragments`; unavailable flank
   positions are `null`. This conservative gate can leave a real contraction
   unresolved when reads do not connect its unique flanks.
4. Apply only resolved replacements. Keep baseline pileup reports in reference
   coordinates and write structural events and the mapping to final coordinates
   separately. Ambiguous paths remain inspectable and prevent final publication.
5. Freeze sequences before annotation; annotation never selects an assembly.

No external assembler/mapper, hashes, second baiting pass or reference-filled
sequence is introduced. A bounded short-read search cannot guarantee resolution
of indistinguishable long repeat copies.

## CLI

```bash
mitovarr assemble-mapping --ref reference.fasta --bam mt.bam \
  --local-assembly --r1 reads.R1.fastq.gz --r2 reads.R2.fastq.gz \
  --emit-variants --emit-coverage --out mapping

mitovarr annotate --fasta mapping/candidate.fasta --species cbr \
  --topology circular --rotate-to-trna-pro \
  --assembly-ambiguity mapping/structural_ambiguity.tsv --out annotation
```

`assemble-mapping` remains an alias of `polish`. Omitting `--local-assembly`
preserves ordinary indexed-BAM pileup. Local assembly requires matching paired
FASTQ; the complete extracted pairs from `fetch-caendr` are suitable. No extra
recruitment or network download occurs. `--topology circular` is an annotation
hypothesis and does not establish read-supported closure.

Defaults: `--local-k 31`, `--local-max-span 2000`, `--local-max-paths 128`.
`--max-records` also bounds unique FASTQ pair names in this mode. The sample graph
has a 2,000,000-kmer cap; each local search has a 200,000-state and four-visits-per-
node cap. A reached search cap leaves the region unresolved. Increasing a cap
does not add sequencing evidence.

### Optional high-depth subsampling

`--subsample` explicitly selects whole named target-BAM fragments before pileup;
`--max-records` (default 500,000) caps their total alignment records, not their
number of pairs. `--subsample-seed` defaults to 42 and requires `--subsample`.
Selection uses read-name priorities and is independent of BAM record order.

With `--local-assembly`, only complete supplied FASTQ pairs whose names were
selected are written to `subsampled.R1.fastq` and `subsampled.R2.fastq` and passed
to graph construction and paired scoring. A retained BAM fragment without a
complete FASTQ pair can support pileup but cannot contribute local pair
evidence. The report's `subsampling` object records input/retained counts, seed
and evidence scope. Coverage and variant fractions describe the retained data;
rare alleles or weak structural support can be lost. No support threshold is
relaxed. Without `--subsample`, exceeding the record cap remains an error.

## Discovery and sequence generation

Mapping signals include clips of at least 15 bp, explicit CIGAR insertions or
deletions of at least 20 bp, and mates on another contig or more than 1,000 bp
away after circular distance correction. Secondary, duplicate and QC-failed
records do not contribute; the existing minimum MAPQ applies. Supplementary
records can contribute their CIGAR breakpoints, with names deduplicated. The
stage does not reconstruct absent supplementary records from `SA` tags.

Breakpoint starts are clustered within seven bases. A cluster needs three named
fragments and support equal to at least 5% of the maximum adjacent pileup depth.
This avoids joining dispersed clipped read ends into a whole-genome interval at
high coverage. Existing repeat-ambiguity intervals and runs of at least 15 N
bases are also considered. Intervals receive 60-bp circular padding, and gaps up
to 150 bp are merged. The search extends at most 350 bp beyond each side to find
unique, sample-supported flanks. An internal rotation places the linear cut away
from suspect regions, so origin-spanning events can be searched in one piece.

Every graph kmer comes from quality-filtered FASTQ and is counted once per named
pair. Edges need at least three fragments. Outgoing branches below 2% of the
strongest outgoing edge are excluded from this dominant-consensus graph,
independent of absolute depth. More strongly supported alternatives remain;
this filter is not a detector of rare structural haplotypes.
Paths are generated
between the two sample-supported flank sequences. The `--min-depth` requirement
also applies during path search: paths with an edge below that depth cannot
produce a consensus and are excluded before consuming the candidate budget.
The number of encountered excluded edges is reported; this mode is not a
low-frequency haplotype detector. Excursions that would need to
revisit the unique left anchor are excluded. No reference sequence is injected
as a graph edge. The unmodified pileup sequence is always candidate 0.

## Whole-contig competition

For every local candidate, substitute that region into the complete pileup
backbone, retaining all other copies, including the original ND5 or NCR source.
Freeze one common fragment cohort using a quality-filtered 21-mer on either mate
that occurs once in every whole-contig candidate and occurs in the local region.
Repeat-only fragments still contribute sample graph sequence but do not become
independent votes for a copy's location. All candidates score the same cohort.

Use the existing full-query affine-gap aligner, including mismatch/base-quality
scores and bounded mate rescue. Estimate the existing empirical fragment model
on unambiguous pileup segments outside suspect windows. If it cannot be estimated
or is marked diagnostic-only, do not substitute a nominal insert size or use
distance as evidence. Each fragment contributes once: alternative placements
within the same library span state retain their best sequence explanation.
Duplicating a placement cannot multiply its evidence. Whole-contig placement
also includes pairs crossing the arbitrary FASTA origin.

A selected path needs at least three fragments supporting each changed boundary,
and all path kmers must meet `--min-depth`. Against every competing path it needs
at least three strong fragment preferences (log-score difference >= 2.3), the
configured `--min-indel-frac` majority among directional preferences (default
0.7), and a total score margin > 5. Exhausted alignment or graph searches remain
unresolved. These are selection heuristics, not calibrated correctness or
heteroplasmy probabilities. SNPs/indels inside the local sequence are selected as
part of the complete sample path; outside these windows the pileup is retained.

Windows are compared independently on the common pileup backbone. This version
does not jointly phase distant structural variants or globally re-polish copies
after several separate replacements. A new source-copy ambiguity outside the
discovered windows may therefore require a further local analysis. Detection is
conditional on the reads supplied; no claim of exhaustive SV discovery is made.

## Outputs and coordinates

| File | Meaning |
| --- | --- |
| `reference_candidate.fasta` | Frozen ordinary pileup consensus, including its accepted CIGAR indels and N calls. |
| `candidate.fasta` | Pileup backbone with supported sample-path replacements. |
| `consensus.fasta` | Published only when every discovered local window is resolved, remaining repeat ambiguity is absent, and the sequence has no N. Circularity remains unverified. |
| `local_assembly/locus_N.fasta` | Baseline candidate 0 and all generated local alternative sequences, with flanks. |
| `local_assembly.gfa` | Local candidate-path graph with shared flank nodes and k-bp overlaps; it is not the complete underlying kmer graph. |
| `structural_variants.tsv` | One row per evaluated window, selected candidate, length delta and unresolved reason. `ref_window_start/end` include flanks, not a minimal breakpoint. |
| `structural_ambiguity.tsv` | Unresolved intervals in final `candidate.fasta` coordinates; accepted by native `annotate --assembly-ambiguity`. |
| `local_assembly/report.json` | Library model, same-cohort pair scores, junction support, search limits, baseline/final coordinate segments and selected minimal edits. |

Existing `variants.tsv`, `coverage.tsv`, `ambiguity.tsv`, `uncertainty.tsv` and
`base_votes.tsv` continue to describe the **baseline reference-coordinate pileup**.
They are not silently relabelled as coverage or alleles on inserted sequence.
The opt-in main report uses `mitovar-rust.polish.v2`; ordinary polish retains v1.
Coordinates are one-based inclusive unless named as offsets. Selected `edit`
objects trim the common prefix/suffix, report both reference flank positions and
the actual replacement string, and describe edits relative to the baseline
after its accepted CIGAR indels. An insertion has zero replaced baseline bases.
Equivalent breakpoint placements within a repeat are not a claim of unique
biological breakpoint resolution.

When a replacement crosses or removes the original origin, the output begins at
its supported left flank, with `output_origin_shifted_to_supported_anchor=true`.
Each region's `baseline_segments` and `output_segments` describe the transform,
splitting intervals at either circular origin. Annotation/tRNA-Pro rotation is a
separate subsequent operation and carries structural ambiguity with it.

## Verification

The portable integration tests cover a 240-bp homologous insertion with a
copy-specific allele and indel-bearing reads, a 300-bp deletion represented only
by clipped junction reads, an origin-spanning insertion, read-order invariance,
source-only coverage inflation, a true coverage gap without a deletion junction,
an indistinguishable tandem-repeat copy count and an explicitly bounded search.
Discovery regressions cover dispersed high-depth clipping and duplicated names.
