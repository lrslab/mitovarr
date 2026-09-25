# Native de novo assembly

The assembler receives paired FASTQ only. For `assemble-baited`, these are the reads already extracted from mitochondrial BAM alignments and their recovered mates; input and validation scope remain `mt-bam-extracted`. No reference, target coordinates, trusted contigs or expected genome length enter candidate generation or selection. All computation described here runs in Rust.

See the [paired-path method](paired_paths.md) for alignment, fragment scoring and circular-path selection.

## Graph construction and cleanup

Canonical k-mers are counted once per pair record across both mates. Bases below baseQ 20 break the rolling window. The defaults are k = 63 and minimum pair-record support 3; odd k values from 15 to 127 are available. Two native 128-bit words store the complete graph k-mer exactly, including its prefix/suffix nodes and reverse complement. Small-k codes retain their previous numeric ordering; read-alignment and seed-recruitment indexes retain their existing single-word representation.

Larger k can distinguish longer repeat contexts, but each read supplies fewer windows and a single low-quality base interrupts more of them. For example, a fully usable 100-bp read has 38 windows at k=63 and only six at k=95. Read length alone does not establish adequate k-mer coverage. Graph keys occupy 32 bytes; `--max-kmers` is an explicit entry limit, not a fixed memory budget. A limit failure does not yield a partial assembly.

Each invocation still builds one main-k graph, with the smaller-k end proposal stage below. The command does not automatically choose among k=63/95/127 outputs. Raw k-mer support mass and fragment scores from different runs are not a shared cross-k score, and annotation never selects a graph.

Retained k-mers become two directed edges, one for each strand, between their (k−1)-base prefix and suffix nodes. The graph is compacted into nonbranching unitigs. With default weak/strong ratio 0.2, the following operations are repeated until no further edges qualify:

- Remove a dead-end tip of at most 2k edges when its maximum support is at most 0.2 times the competing boundary-edge support.
- Remove a weak short path between two junctions when a stronger path joins the same endpoints within 2k edges. All edges of the alternative must satisfy the support ratio, and the two spellings must differ by at most three edits. The bounded search leaves uncertain cases unchanged.
- Disconnect the boundary k-mer of a short weak path at a dominant competing branch. This handles intersecting error paths that do not form a simple bubble. The remaining sequence is retained. This rule excludes long paths and self-loop alternatives, rather than assigning a tandem-repeat copy number.

Deletions are symmetric across reverse complements. Local competing support is used instead of a global high-coverage cutoff. Setting `--max-error-fraction 0` disables these operations. These are consensus heuristics: low-frequency biological alleles or structural paths can also be removed, so the output is not a heteroplasmy estimate.

Tip/bubble simplification and local weak-edge disconnection are established assembly ideas; this package implements its own bounded versions. See [Zerbino & Birney, 2008](https://doi.org/10.1101/gr.074492.107) and [Nurk et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411777/).

## Read-supported repeat traversal

For a repeat unitig with multiple entrances and exits, the assembler scans the assembly reads for exact paths covering the entire unitig and at least five consecutive k-mers in each flanking unitig. Each named fragment contributes at most once to a route, even if both mates span it. This exact-threading stage does not use mate distance; a separate paired-path comparison follows it.

A selected entrance/exit route must have at least `--min-junction-fragments` supporting fragments (default 3). At each end, the **sum** of all competing observed routes must be at most `--max-error-fraction` times its support. Counts below three still contribute to that competing sum. Every entrance and exit must participate in a unique dominant pairing; missing evidence or a conflicting pairing leaves the repeat unresolved. Selected routes create separate graph traversals through the repeat and permit longer contigs. Unsupported tandem-copy counts are not filled in.

The method supports a primary sequence path. Rare alternate molecules can be suppressed by the dominance rule; all observed repeat-route counts are retained in `report.json`.

## Complete paired-path comparison

The retained graph then supplies bounded paths between reliable sample anchors. Native full-query affine-gap alignment allows mismatches and indels. A fixed R1 cohort estimates empirical paired orientation and span mass, including unmeasured mates in the denominator. Competing paths use the same named-fragment pool, one joint explanation per fragment. Read and mate geometry support are not local k-mer votes.

Pure single-SNP alternatives between the same two anchors use the actual differing base plus 20 shared bases on each side for spanning evidence. This prevents the k-dependent graph corridor from requiring a read longer than the input to span one SNP. The full fragment still scores against the complete path; indel, multi-site, competing-exit and repeat evidence intervals are unchanged.

For local paths and cycles, either mate may supply the shared unique anchor. R1 takes precedence when both qualify, before candidate alignment or mate-distance scoring. Each named fragment enters once per comparison; the learned orientation/span definition is symmetric under mate exchange. Library calibration retains its fixed R1 cohort. For cycles returning to a common backbone, reverse-complement views are merged. Conserved sequence separates variable blocks; adjacent blocks need joint allele observations on the same fragments before a path is qualified as phased for final publication. The primary output can retain the highest-scoring complete candidate while recording its unresolved blocks and phase. An unbounded repeat tail or unresolved competing exit remains unresolved. See [paired-path details](paired_paths.md) for the score, search limits, assumptions and output fields.

## End connections and simple repeat copies

Before path enumeration, native k=47 and k=31 graphs propose connections between eligible open ends when those k values are smaller than the main k. Every retained k-mer uses the same base-quality and pair-support thresholds. A proposed connector must match the full main-k overlap at both ends and have a reciprocal spelling; alternatives found at different k are pooled, not ranked by length. A search cap discards affected new connections rather than treating a partial search as unique. Connections then undergo the same complete-path fragment scoring. This searches the current FASTQ only, without another recruitment round.

Simple cyclic strongly connected components with one internal successor/predecessor per segment are represented as a turn count plus an exit position. Candidate counts are bounded by the existing path domain; the unsearched high-copy tail stays explicit. The repeat report compares observed high-quality R1 31-mer occurrence counts with the number of motif instances in each candidate at local unique-flank depth. Repeated offsets in R1 are retained; R2 is excluded from this statistic to avoid mate-overlap inflation. These predictions assume similar local sampling and no extra motif source. They are diagnostic, do not add an independent likelihood to the same reads, and do not establish an upper copy bound.

## Outputs and status

`primary.fasta` is the default downstream representative. Observed-majority local paths can extend it across supported short alternatives; a common-backbone cycle uses one whole candidate's paired-fragment score. `primary_ambiguity.tsv` marks the corresponding variable regions and open ends, and `report.json.primary` records its source and `final_verified`. Alternative path FASTAs and their read scores remain available. A primary sequence is a hypothesis when final status is unresolved; a closed graph walk alone does not establish NCR repeat copy number. See [primary path rules](paired_paths.md#primary-path-output).

`assembly_graph.gfa` contains the cleaned, compacted graph **before** read-supported repeat traversal. Its `dp` tag is the mean pair-record k-mer support along a unitig, not calibrated molecule depth. `read_threading.repeats` identifies these same GFA unitigs and records their observed and selected entrance/exit routes. `candidates.fasta` retains sequences from strict repeat/path traversal; one candidate may span multiple GFA unitigs. Reverse-complement duplicates are merged.

`path_graph.gfa` exports the segment graph used by path enumeration, including small-k candidate connectors. Its numeric IDs match `source_segment`, `target_segment` and `intermediate_segments`. Links have main k−1 overlap, except already-trimmed isolated circles, whose self-links have zero overlap. `TP` distinguishes original sample-graph segments from end-bridge candidates. `end_bridge_candidates.fasta` and `paired_path_resolution.end_assembly` record proposed connections and search limits. Per-locus `repeat_copy_model` supplies explicit periods, turn counts and coverage predictions.

The assembly report keeps its `mitovar-rust.assembly.v1` schema and adds `graph_simplification`, `read_threading` and `paired_path_resolution`. `directed_branch_nodes` describes the cleaned k-mer graph before threading; `read_threading.unresolved_unitig_junctions` counts branching unitig instances after threading. These are different graph representations and must not be compared as the same statistic.

A unique supported main circle is judged independently of unrelated residual small components. It must have quality-filtered matching-base coverage from gapped alignments and support for the actual graph joins or literal output-origin junction. Its fraction of retained k-mer support mass is reported but is not a publication threshold. Competing complete circles remain unresolved. `path_candidates.fasta` preserves all local hypotheses; `complete_path_candidates.fasta` preserves complete cycles on a common backbone, including unselected cycles. Original-FASTQ and seeded workflows use all original pairs for this final validation; BAM-baited mode uses the extracted pairs. Unresolved runs retain their candidates and return exit code 2. No candidate is padded or ordered using a reference.

This implementation does not yet use multi-k iteration, general gapped read-to-graph alignment, reference scaffolding or general tandem-repeat copy inference. Its gapped aligner compares finite path spellings, and its paired model never inserts unobserved bases. Multiple components remain explicit.
