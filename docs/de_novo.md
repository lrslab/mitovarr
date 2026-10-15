# Native de novo assembly

De novo assembly has two input routes: `assemble-baited` for BAM-extracted
pairs and `assemble-seeded` for seed recruitment from original paired FASTQ.
The seed is used only for recruitment.
For reference-guided processing of existing alignments, use the separate
[mapping assembly](mapping_local_assembly.md) workflow.

The assembler receives paired FASTQ only. For `assemble-baited`, these are the reads already extracted from mitochondrial BAM alignments and their recovered mates; input and validation scope remain `mt-bam-extracted`. No reference, target coordinates, trusted contigs or expected genome length enter candidate generation or selection. All computation described here runs in Rust.

See the [paired-path method](paired_paths.md) for alignment, fragment scoring and circular-path selection.

## Graph construction and cleanup

Canonical k-mers are counted once per pair record across both mates. Hard base-quality filtering is disabled by default (`--min-baseq 0`); A/C/G/T bases remain eligible regardless of quality. Ambiguous bases still break the rolling window. An explicit positive `--min-baseq` restores the optional cutoff for recruitment, graph construction, threading and matching-base validation. Alignment seeds have no independent quality cutoff, while full-query alignment scores retain quality-dependent emissions. The defaults are k = 63 and minimum pair-record support 3; odd k values from 15 to 127 are available. Two native 128-bit words store the complete graph k-mer exactly, including its prefix/suffix nodes and reverse complement. Small-k codes retain their previous numeric ordering; read-alignment and seed-recruitment indexes retain their existing single-word representation.

Larger k can distinguish longer repeat contexts, but each read supplies fewer windows; ambiguous bases, or low-quality bases when a cutoff is explicitly enabled, interrupt more of them. For example, a fully usable 100-bp read has 38 windows at k=63 and only six at k=95. Read length alone does not establish adequate k-mer coverage. Graph keys occupy 32 bytes; `--max-kmers` is an explicit entry limit, not a fixed memory budget. A limit failure does not yield a partial assembly.

Each invocation still builds one main-k graph, with the smaller-k end proposal stage below. The command does not automatically choose among k=63/95/127 outputs. Raw k-mer support mass and fragment scores from different runs are not a shared cross-k score, and annotation never selects a graph.

Main and end-graph counting reuse the per-pair deduplication buffer and update
each retained key with one table lookup. Every pair is still streamed in full,
each canonical k-mer still contributes at most once per pair, and a full key
table accepts further observations of existing keys while rejecting a new key.

Retained k-mers become two directed edges, one for each strand, between their (k−1)-base prefix and suffix nodes. The graph is compacted into nonbranching unitigs. The following operations repeat until no further edges qualify. Tip/weak-link cleanup and repeat routes retain `--max-error-fraction 0.2`; local bubbles use the independent `--max-bubble-fraction 0.49`.

- Remove a dead-end tip of at most 2k edges when its maximum support is at most 0.2 times the competing boundary-edge support.
- Remove a weak local path between two junctions when a stronger simple path joins the same endpoints within 3k edges. Compare length-weighted mean pair-record k-mer support across the complete paths, not the weakest competing k-mer against the strongest weak k-mer. Require weak/strong mean support at most `--max-bubble-fraction` and at most `max(3, weak_path_edges / 10)` sequence edits. Search compacted unitigs deterministically with a 4,096-state bound; retain the weak path if that bound is reached before finding a qualifying alternative. Exclude circular unitigs, self-loop traversals, repeated vertices and conjugate traversals. This is a bounded local consensus simplification, not unrestricted repeat resolution.
- Disconnect the boundary k-mer of a short weak path at a dominant competing branch. This handles intersecting error paths that do not form a simple bubble. The remaining sequence is retained. This rule excludes long paths and self-loop alternatives, rather than assigning a tandem-repeat copy number.

Deletions are symmetric across reverse complements. Local competing support is used instead of a global high-coverage cutoff. Setting `--max-error-fraction 0` disables these operations. Setting `--max-bubble-fraction 0` disables bubble collapse alone. Both ratios must be at least zero and below 0.5. The report records both thresholds and the bubble coverage/edit/search rules. These are consensus heuristics: low-frequency biological alleles or structural paths can also be removed, so the output is not a heteroplasmy estimate.

After local simplification, a separate pass removes disconnected sequencing-error
shadows. Components are joined through canonical (k−1)-base nodes so both
strands are decided together; circular strand views share a canonical node.
Process components by decreasing maximum k-mer support. Remove a component
when its maximum is at most `--max-error-fraction` (default 0.2) times the mean
support of an already retained unitig sharing a canonical 31-mer. Components
without a shared seed, and stronger isolated homologs, remain. Filter the
removed edges and compact again; original sample k-mer counts remain intact.
This pass is skipped when `--max-error-fraction 0` is used.

`graph_simplification.shadow_components` records the rule, strand-collapsed
component count and removed directed k-mer count; it is null when disabled.
The existing `cleaned_canonical_kmers` and `cleaned_directed_branch_nodes` now
refer to the graph after this pass. A connected bubble or deletion is not a
separate component for this rule. Earlier tip/bubble/weak-link rules still
apply: these checks do not guarantee retention of every low-frequency allele.
A weak-link remnant remains protected only when its original unitig had a
stronger unitig between the same full graph flanks, a length difference above
three bases, and sufficient support at every stronger-arm k-mer; `protected_weak_link_components` counts qualifying shadows kept
for that reason. This preserves minor deletion-junction candidates instead of
erasing them merely because the earlier cleanup disconnected their boundaries.
Isolated homologous biological sequence at or below the threshold, including
NUMT fragments, can also be removed. No reference or quality-only branch
classification is used, and complete-publication checks remain unchanged.

Permissive collapse does not establish structural uniqueness. If a removed bubble fails the conservative conditions (both arms at most 2k edges, at most three edits, and the weak arm's maximum support at most `--max-error-fraction` times every retained edge's support), both spellings are saved in `bubble_alternatives.fasta` and reported in `graph_simplification.retained_bubble_alternatives`. These are retained hypotheses, not independent variant calls. The primary may use the simplified graph, but a final candidate touching the original supported graph component of any such bubble remains unresolved with `supported_bubble_alternatives`. A disconnected component does not veto an otherwise supported main circle. Component membership is determined from the existing in-memory sample k-mers, before cleanup, so subsequent overlapping collapses cannot erase the uncertainty. No additional input scan is needed.

`primary_ambiguity.tsv` adds `collapsed_bubble_component_unresolved` across the affected primary sequence. This is deliberately a component-wide flag, not an inferred local variant interval. `validation.passed` continues to describe read coverage and junction checks; a final assembly also requires an empty `validation.unresolved_collapsed_bubbles` list. Balanced or phased supported alternatives still use the existing path-resolution rules. The conservative cleanup remains a consensus heuristic and is not proof that every low-frequency allele was an error.

Tip/bubble simplification and local weak-edge disconnection are established assembly ideas; this package implements its own bounded versions. Whole-path coverage comparison is also used by the [SPAdes bulge remover](https://github.com/ablab/spades/blob/v4.3.0/src/common/modules/simplification/bulge_remover.hpp); our implementation retains an explicit sequence-edit bound and does not import or invoke SPAdes code. See [Zerbino & Birney, 2008](https://doi.org/10.1101/gr.074492.107) and [Nurk et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411777/).

### Mixed-sample signal

`report.json.mixed_sample` counts allele bubbles left after cleanup: pairs of
non-circular unitigs between the same graph nodes whose lengths differ by at
most three bases, counted once per strand pair, with the minor branch's share
of mean k-mer support. Cleanup removes weak branches as errors, so pairs whose
minor share is at least 0.2 indicate a well-supported second haplotype;
`status` becomes `mixed_sample_signal` at three or more such bubbles and
`minor_fractions` lists every pair. Such a sample can fragment the graph and
leave a short primary (e.g. 1.9 kb) without any assembler fault. This is a
signal, not proof of heteroplasmy or contamination: minor alleles cleaned as
errors and minor sequences too divergent to share k-mers are not counted. The
reference-guided mapping report has a separate pileup-based summary
(`base_voting.mixed_sample`).

## Read-supported repeat traversal

For a repeat unitig with multiple entrances and exits, the assembler scans the assembly reads for exact paths covering the entire unitig and at least five consecutive k-mers in each flanking unitig. Each named fragment contributes at most once to a route, even if both mates span it. This exact-threading stage does not use mate distance; a separate paired-path comparison follows it.

A selected entrance/exit route must have at least `--min-junction-fragments` supporting fragments (default 3). At each end, the **sum** of all competing observed routes must be at most `--max-error-fraction` times its support. Counts below three still contribute to that competing sum. Every entrance and exit must participate in a unique dominant pairing; missing evidence or a conflicting pairing leaves the repeat unresolved. Selected routes create separate graph traversals through the repeat and permit longer contigs. Unsupported tandem-copy counts are not filled in.

The method supports a primary sequence path. Rare alternate molecules can be suppressed by the dominance rule; all observed repeat-route counts are retained in `report.json`.

## Complete paired-path comparison

The retained graph then supplies bounded paths between reliable sample anchors. Native full-query affine-gap alignment allows mismatches and indels. A fixed R1 cohort estimates empirical paired orientation and span mass, including unmeasured mates in the denominator. Competing paths use the same named-fragment pool, one joint explanation per fragment. Read and mate geometry support are not local k-mer votes.

Pure single-SNP alternatives between the same two anchors use the actual differing base plus 20 shared bases on each side for spanning evidence. This prevents the k-dependent graph corridor from requiring a read longer than the input to span one SNP. The full fragment still scores against the complete path; indel, multi-site, competing-exit and repeat evidence intervals are unchanged.

For local paths and cycles, either mate may supply the shared unique anchor. R1 takes precedence when both qualify, before candidate alignment or mate-distance scoring. Each named fragment enters once per comparison; the learned orientation/span definition is symmetric under mate exchange. Library calibration retains its fixed R1 cohort. For cycles returning to one common backbone, reverse-complement views are merged. Conserved sequence separates variable blocks; adjacent blocks need joint allele observations on the same fragments before a path in this family is qualified as phased for final publication. This pass does not provide general cross-locus phasing for multi-anchor circles. The primary output can retain the highest-scoring complete candidate while recording its unresolved blocks and phase. An unbounded repeat tail or unresolved competing exit remains unresolved. See [paired-path details](paired_paths.md) for the score, search limits, assumptions and output fields.

## End connections and simple repeat copies

Before path enumeration, native k=47 and k=31 graphs propose connections between eligible open ends when those k values are smaller than the main k. Every retained k-mer uses the same optional base-quality cutoff (disabled by default) and pair-support threshold. A proposed connector must match the full main-k overlap at both ends and have a reciprocal spelling; alternatives found at different k are pooled, not ranked by length. Discovered proposals are saved even when a search cap withholds them from the graph. Any incomplete competing target or k keeps the affected source ineligible for new joins, including its reciprocal view; partial enumeration is not treated as uniqueness. Admitted connections then undergo the same complete-path fragment scoring. This searches the current FASTQ only, without another recruitment round.

Before walking a smaller-k graph, reverse reachability identifies nodes that can reach each particular target suffix. The source/k budgets of 50,000 states and 128 candidate paths are divided equally between reachable targets. Unreachable targets consume no traversal budget, and a repetitive route to one target cannot spend another target's share. Same-end reverse-complement connections are classified explicitly but remain possible proposals. Reachability is only an over-approximation: every arriving sequence must still match the complete main-k target flank, including short bridges whose target prefix already overlaps the initial source context. Repeats with an exit retain their unresolved-copy behavior. Source summaries and target-specific searches, budgets and cutoffs are recorded in `paired_path_resolution.end_assembly`; a limit remains distinct from an unreachable target.

Budget shares use integer floor division; spare capacity is not allocated by
segment ID. Searches still stop on arrival at the first complete endpoint,
including another target, so target isolation does not create connectors that
skip intervening anchors. A proposal's `search_complete` means all tested k
searches for that source/target pair finished. `join_eligible` additionally
requires that no competing target/k search from the source or the matching
reciprocal source was limited. An eligible proposal still needs downstream
fragment scoring. Independent existing circular components retain their
existing validation rules.

Target-specific reachability adds preprocessing work per endpoint. To bound
the source/target metadata as well as the walk states, end assembly allows at
most 10,000 directed endpoint pairs (source count multiplied by target count)
when a smaller-k end search is enabled.
This bounds the two-k target report to 20,000 entries. Exceeding this limit
skips optional end search before allocating pair metadata or rescanning reads.
The complete main graph and available primary/candidate sequences are retained;
`end_assembly.status=resource_limited` reports the observed count and cap.
Complete publication is withheld with reason `end_search_resource_limit`.
This is a graph-search metadata limit, not a FASTQ input count, and it does not
trigger automatic subsampling.

### Subsampling after a `--max-kmers` overflow

Every de novo run first counts k-mers from all complete pairs. Without
`--subsample`, more distinct k-mers than `--max-kmers` is an error, as before.
With `--subsample`, that full-input overflow is the recorded trigger: the run
keeps complete pairs whose seeded fragment-name priority lies in the lowest 1/2,
then 1/4, and so on to 1/64, stopping at the first fraction whose graph fits.
Mates stay together and the selection does not depend on read order; change it
with `--subsample-seed` (default 42). Retained pairs are written as
`subsampled.R1.fastq` / `subsampled.R2.fastq` and used for graph construction,
threading, path scoring and optional layouts. Final read validation still uses
**all** complete pairs. `report.json.subsampling` records the trigger and cap,
each attempted fraction, original and retained pair/base counts, the seed and
this evidence scope; when the full input fits it reports `applied=false` and no
file is written. Input and validation scopes are unchanged. Typical triggers are
very deep extracted libraries whose sequencing-error k-mers far outnumber the
genome's; lower retained depth can also lose weak junction or repeat support.

Simple cyclic strongly connected components with one internal successor/predecessor per segment are represented as a turn count plus an exit position. Candidate counts are bounded by the existing path domain; the unsearched high-copy tail stays explicit. The repeat report compares observed R1 31-mer occurrence counts passing the optional base-quality cutoff with the number of motif instances in each candidate at local unique-flank depth. Repeated offsets in R1 are retained; R2 is excluded from this statistic to avoid mate-overlap inflation. These predictions assume similar local sampling and no extra motif source. They are diagnostic, do not add an independent likelihood to the same reads, and do not establish an upper copy bound.

## Outputs and status

NCR1 completion is optional for a partial deliverable. The additive
`sequence_output` report describes sequence availability independently of the
strict completion status. Evaluate genes on the frozen primary sequence using
[the partial-result workflow](partial_results.md); unresolved NCR1 alone does
not fail gene-level annotation checks. Partial availability does not relax path validation.

`primary.fasta` is the default downstream representative. Observed-majority local paths can extend it across supported short alternatives; a common-backbone cycle uses one whole candidate's paired-fragment score. `primary_ambiguity.tsv` marks the corresponding variable regions and open ends, and `report.json.primary` records its source and `final_verified`. Alternative path FASTAs and their read scores remain available. A primary sequence is a hypothesis when final status is unresolved; a closed graph walk alone does not establish NCR repeat copy number. See [primary path rules](paired_paths.md#primary-path-output).

`best_circular.fasta` separately selects the best available eligible circular
representative, even when the main primary component is linear. Its matching
`best_circular_ambiguity.tsv` and `report.json.best_circular` retain selection
provenance and unresolved evidence. Existing whole-cycle scores choose within
families; distinct k-mer support mass ranks different components, with a
validated final circle taking precedence. No eligible circle means no
best-circular FASTA. Selection does not change strict publication or infer
mitochondrial identity. See [best circular candidate](paired_paths.md#best-circular-candidate).

`assembly_graph.gfa` contains the cleaned, compacted graph **before** read-supported repeat traversal. Its `dp` tag is the mean pair-record k-mer support along a unitig, not calibrated molecule depth. `read_threading.repeats` identifies these same GFA unitigs and records their observed and selected entrance/exit routes. `candidates.fasta` retains sequences from strict repeat/path traversal; one candidate may span multiple GFA unitigs. Reverse-complement duplicates are merged.

`path_graph.gfa` exports the segment graph used by path enumeration, including small-k candidate connectors. Its numeric IDs match `source_segment`, `target_segment` and `intermediate_segments`. Links have main k−1 overlap, except already-trimmed isolated circles, whose self-links have zero overlap. `TP` distinguishes original sample-graph segments from end-bridge candidates. `end_bridge_candidates.fasta` records admitted graph connections. `end_bridge_proposals.fasta` also preserves withheld spellings with separate proposal IDs and explicit eligibility; these records do not imply graph edges or supported joins. `paired_path_resolution.end_assembly` records their status and search limits. Per-locus `repeat_copy_model` supplies explicit periods, turn counts and coverage predictions.

The assembly report keeps its `mitovar-rust.assembly.v1` schema and adds `graph_simplification`, `read_threading` and `paired_path_resolution`. `directed_branch_nodes` describes the cleaned k-mer graph before threading; `read_threading.unresolved_unitig_junctions` counts branching unitig instances after threading. These are different graph representations and must not be compared as the same statistic.

A unique supported main circle is judged independently of unrelated residual small components. It must have matching-base coverage from gapped alignments, subject to any explicitly enabled base-quality cutoff, plus support for the actual graph joins or literal output-origin junction. Its fraction of retained k-mer support mass is reported but is not a publication threshold. Competing complete circles remain unresolved. `path_candidates.fasta` preserves all local hypotheses; `complete_path_candidates.fasta` preserves complete cycles on a common backbone, including unselected cycles. Seeded assembly uses all original pairs for this final validation; BAM-baited mode uses the extracted pairs. Unresolved runs retain their candidates and return exit code 2. No candidate is padded or ordered using a reference.

This implementation does not yet use multi-k iteration, general gapped read-to-graph alignment, reference scaffolding or general tandem-repeat copy inference. Its gapped aligner compares finite path spellings, and its paired model never inserts unobserved bases. Multiple components remain explicit.

The optional `--multi-k-layout --coverage-layout` extension below pools additional local graph walks and uses longer exact read contexts. It does not rerun the entire main assembly at multiple k values or implement iterative graph replacement.


## Coverage-inferred layouts

`--coverage-layout` is an opt-in hypothesis output for unresolved linear primary
backbones, including NCR repeat gaps. Assembly still receives only its declared
sample FASTQ scope. No species, NCR coordinate, reference or expected genome
length enters this procedure. A primary that is already circular is reported
without an additional inferred layout.

The two boundaries are unique sample-backbone k−1 sequences. Up to 500 terminal
bases may be removed to reach unique boundaries, with both trims reported.
The missing corridor is sought in raw sample graphs at the main k and smaller
k=47/31, retaining the normal minimum pair-record k-mer support. If no layout
is found, k=21 then k=15 provide exploratory fallbacks. One central unique backbone node is excluded to prevent traversing the whole
backbone a second time. Other interior words remain available: an NCR can
reuse sequence present once elsewhere in the backbone. Every added base follows an observed graph edge; disconnected
flanks remain unresolved. If all normally supported graphs yield no layout,
the procedure retries with minimum pair-record support 1. These singleton edges
are allowed only for the explicitly inferred layout, and may be sequencing errors.
The selected graph threshold and weakest corridor support are always reported.

Coverage counts **R1 k-mer window occurrences**, including repeated offsets;
R2 is excluded from this statistic to avoid double-counting overlapping mates.
It is not the GFA pair-presence count. The arithmetic mean of unique canonical
31-mer depths in the two 500-base backbone flanks supplies the single-copy
baseline. Short backbones or unavailable flank depth use the whole-backbone
mean. Both local means and the whole-backbone mean remain in the report.
For each canonical 31-mer, observed depth divided by this baseline predicts
its total genomic multiplicity. Candidate circular spellings predict integer
multiplicities. Mean squared residual over one common candidate-union domain
ranks all candidates, including unchosen motifs in each comparison; invariant
backbone k-mers do not affect the rank. This score is a heuristic fit, not a
correctness probability or an independent read likelihood.

A shortest connected walk is retained as one seed candidate per graph before
coverage-guided exploration. It is found by Dijkstra search over the finite
retained graph, with the same appended-base bound. The final rank uses coverage,
not path length; seed provenance is reported.

Coverage guides bounded graph walks. Per k the guided-search limits are 50,000 visited states,
1,024 pending states, 256 completions in addition to the seed, and 10,000 appended bases. Per-code usage
is explored up to `ceil(1.5 * observed_depth / baseline) + 2`, capped at 4,096,
without removing already present backbone copies. These are search bounds,
not confidence intervals. Limits apply separately to each k and support threshold. k=21/15 searches use a 15-mer occurrence model for
traversal; final cross-k ranking always uses 31-mers. A limited search never
invalidates candidates discovered at another k. No search claims completeness.

### Additional cross-k proposals

`--multi-k-layout` retains the original candidate generator and adds searches at
k=15,21,31,47,63,95,127 (plus the requested main k). The same FASTQ pass collects
the graph counts. A k without usable observed windows or sufficient flank length
is reported as unavailable. An existing layout does not disable another k.
All ordinary-support searches run; additional singleton-support searches run only
when the union contains no ordinary-support candidate. Legacy fallback candidates
are preserved even when a new ordinary-support path is discovered later.

The additional graphs share a long sample-backbone flank pair, up to 126 bases.
An internal occurrence of the shorter target node does not terminate a walk:
the full right flank must match. When available, a unique long central backbone
word prevents a second complete backbone traversal. Its short constituent words
remain traversable. This guard is an explicit excluded domain, not evidence that
all homologous repeat copies have been considered.

Each k/support search has 50,000 additional state visits in total: 5,000 for a
nearby-exit seed, 22,500 for coverage, and 22,500 for read-context traversal.
Each queue holds at most 1,024 states. The two guided schedules each permit 128
completions; the seed permits one. At most 10,000 bases plus the two long flanks
are spelled. These budgets are additional to the legacy searches. Graph tables
retain the existing per-k `--max-kmers` limit; memory can grow with the number of
k values. Every exhausted domain remains explicit and no search claims
exhaustiveness.

At a shorter-k branch, the context schedule checks the longest observed higher-k
window distinguishing the possible next bases, requiring the same pair-record
support threshold. It prefers extensions with that continuous read context.
Missing context never deletes an alternative; the independent coverage schedule
and legacy candidates remain. This is bounded proposal scheduling, not a read
likelihood, a hard phasing rule or arbitrary gapped read-to-graph alignment.
It cannot distinguish identical repeats beyond the available read context.

Additional walks count 31-mers from the actual spelled suffix, including windows
crossing k=15/21 edges. Fixed backbone copies are subtracted once from available
depth. Final coverage ranking uses the common 31-mer domain for the entire union;
the existing eight-path structural shortlist and baseline are then compared on
one cohort frozen against that full union. Sequence duplicates merge k/support
and generation-method provenance; they do not add fragment votes. Scores from
separate runs with different cohorts must not be directly subtracted.

Additional graph files are `coverage_multik_k*_min*.gfa`. Their reported unitig
walks append sequence after each k-1 overlap to the retained input backbone;
remove the reported long closing overlap, then normalize the circle. The
`multi_k_layout` metadata records the requested schedule, added budget and
unchanged main-backbone k. Coverage representatives, read-preferred diagnostics,
primary output and final-publication rules retain their separate roles.

### Outputs

- `coverage_layout.fasta`: best inferred circle, when a connected candidate exists;
  its header explicitly says `final_verified=false` and `evidence=coverage_inferred`.
- `coverage_layout_candidates.fasta`: all generated, ranked and normalized alternatives.
- `coverage_layout.json` and `report.json.coverage_layout`: coverage baselines,
  copy-fit error, exact GFA walks, candidate 31-mers without R1 observations, search limits,
  candidate count and sensitivity to baseline multipliers 0.8, 1.0 and 1.2.
- `coverage_layout_k*_min*.gfa`: sample corridor and retained backbone. `dp` is mean
  R1 occurrence depth at the reported guidance k; `FC` is mean pair-record
  graph-k support. Walk IDs refer to these graphs before FASTA normalization.
- `coverage_layout_ambiguity.tsv`: inferred corridor and boundary intervals in
  normalized FASTA coordinates, suitable for `annotate --assembly-ambiguity`.
- `coverage_layout_depth.tsv`: each circular 31-mer's start position, multiplicity,
  observed R1 windows and depth divided by candidate multiplicity. This is a
  copy-depth diagnostic, not independent per-base mapped coverage.

Candidate export no longer drops records below coverage rank 32. Rankings also
retain both shortest-seed and coverage-search provenance when the same spelling
is found twice, plus its source k/support graphs. `paired_comparison` records
the structure-diverse shortlist, all-pool cohort size, read-preferred coverage
rank (0 denotes the open baseline), ties, and the unknown score contrast bound.
The coverage representative is unchanged by these diagnostic scores.

The separate inferred layout does not replace `primary.fasta`, modify strict
candidate selection, create `consensus.fasta`, or change unresolved exit code 2.
Coverage can suggest copy usage but cannot resolve repeat ordering or phase.
AT-dependent sampling, PCR duplicates and other genomic motif sources can bias
it. Use the recorded sensitivity and read evidence to assess the hypothesis;
annotate frozen layouts separately to identify NCR1 and its boundaries.

Coverage-layout reports include `paired_comparison` against the linear backbone and
`read_validation` of the displayed layout. The separate matching-base depth file is
`coverage_layout_validation/coverage.tsv`. These are diagnostics and do not make
a coverage hypothesis final; see [shared layout evidence](paired_paths.md#shared-whole-layout-read-evidence).

Whole-layout scores use a fixed query/target common-flank pivot and report ties,
unknown reasons and compressed per-fragment span contributions. The displayed
coverage representative remains unchanged by diagnostic scores.

The de novo graph-count, node-index and per-pair deduplication tables use an
internal deterministic numeric-key hasher without new dependencies. Full k-mer
words remain the keys, and graph/candidate ordering is explicitly sorted. Bubble
search marks unavailable unitigs and their reverse-complement twins once per
cleanup round, instead of revisiting every k-mer for each competitor lookup.
