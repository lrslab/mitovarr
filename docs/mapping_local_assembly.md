# Mapping assembly and local reconstruction

`assemble-mapping` performs mapping assembly from a matching reference FASTA
and an existing coordinate-sorted, indexed BAM. Its baseline is a
fragment-aware pileup consensus. `--local-assembly` defaults to fast consensus
from the pileup and all supplied matching paired FASTQ. It writes the available
sequence and uncertainty without full structural scoring. Add `--deep-validation`
to use the previous local search and scoring workflow. See the
[mapping assembly examples](README.md#mapping-assembly).

## Fast consensus (default)

Build the sample graph from every supplied complete pair, with the existing
quality filters and record limits. Preserve called bases, and attempt only
sample-supported N fills at the existing span using the masked-span method
below. Graphs at smaller k are built lazily when a supported anchored window
contains N. A window without N does not trigger bridge searches. Full path
enumeration, empirical fragment-model fitting, alternative-length bridge
searches and paired-fragment comparisons are deferred.

Fast local mode enables best-effort output automatically; the existing
`--best-effort` flag is still accepted. It writes `consensus.fasta` and
`maximum_consensus.fasta`, with `maximum_consensus.ambiguity.tsv`, and exits 0
when the output is written. `candidate.fasta` retains the pileup sequence and
`structural_ambiguity.tsv` describes its suspect regions.

Reports say `validation_mode=fast_consensus`; each suspect window has
`scoring_status=not_run_fast_consensus`. Scores, best-candidate fields, fragment
counts and unknown-fragment counts are null because they were not measured.
This is separate from a measured zero or a search-limit failure. All suspect
regions remain unresolved, including any filled region. Fast mode always has
`strict_validation_passed=false`, `final_published=false` and
`circularity_verified=false`, even for an N-free sequence.

`local_assembly/locus_N.gap_fill.fasta` retains the fixed-span fill proposals when
available. Fast mode does not write paired-fragment diagnostic gzip files or
read-preferred structural alternatives. `legacy_bridge` report fields retain
their historical names for output compatibility; their `fixed_span_gap_fill`
method and deferred comparison status identify fast fills.

## Deep structural validation (opt-in)

Add `--deep-validation` and choose a new output directory to enumerate local
structures and run the full paired-fragment comparison. Add `--best-effort` to
retain a primary even when this validation fails. Mapping `--coverage-layout`
requires `--deep-validation`; `--multi-k-layout` still requires coverage layouts.

An unresolved NCR1 does not prevent retaining `candidate.fasta` and evaluating
the rest of the sequence. `sequence_output` reports this availability separately
from strict completion. Use [independent gene checks](partial_results.md) with
the matching structural uncertainty. Optional proposal methods below retain
strict replacement gates and complete-output claims.

Deep validation follows ordinary
BAM pileup consensus. The reference locates suspect windows and flanks, while
paired FASTQ supplies every graph edge and every inserted base. This is a
reference-guided workflow, separate from BAM-baited and seeded de novo assembly.

When ordinary graph enumeration reaches a resource limit, candidate FASTA/GFA and path
support are retained, but the ineligible subset is not scored. The window has
`scoring_status=not_run_graph_search_limit`; fragment counts, scores, junction
counts and best-candidate fields are JSON `null` (unmeasured). An exhaustive
bounded search uses `scoring_status=evaluated`. Skipping an ineligible subset
does not make an incomplete search resolved or imply that the sample lacks
supporting reads.

The separate optional `--best-effort` and `--coverage-layout` comparisons
can score hypotheses from incomplete searches without selecting strict replacements.

1. Collect named-fragment breakpoint evidence while reading BAM. Combine soft
   clips, large CIGAR indels, discordant mate endpoints and unresolved consensus
   intervals. Handle intervals crossing the mitochondrial origin.
2. Find supported flanks and enumerate bounded sample graph paths between them.
   Include the unmodified pileup candidate in each comparison. Keep graph and
   candidate sequence outputs for inspection; report resource limits. When the
   main graph finds no connecting path without hitting a limit, lazily build
   sample-only k=21 and k=15 graphs (only values below `--local-k`). Pool their
   paths, retaining the complete original main-k flanks on both sides. The
   traversal follows sample edges through the entire left flank and stops only
   at the complete right flank. Internal copies of shorter flank sequences do
   not terminate the search or count as revisiting the complete left anchor.
   Apply
   the same minimum depth and path/span bounds. Any fallback or pooled-path
   limit keeps the window unresolved. Smaller k never supplies reference-filled
   gap sequence. Each window reports `graph_proposals` with k, paths, states and
   limits; accepted proposals still require the fragment comparison below.
3. Reuse the native full-query affine-gap aligner and empirical fragment model.
   Compare the same fragments against whole-contig alternatives, retaining the
   original homologous locus in every alternative. Shared repeat reads must not
   multiply their contribution simply because a candidate contains more copies.
   Require read-supported connections on both sides of a replacement. A contraction
   of at least `--local-k` bases additionally requires three named fragments
   that each observe both unique unchanged 21-bp flanks in sequence passing
   the optional base-quality cutoff, with a compatible joint placement (one spanning read or usable
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

## Maximum consensus and legacy bridges

Fast local mode and deep mode with `--best-effort` prioritize a usable primary
consensus. Unresolved structure, graph limits and a read score favoring the
N-containing baseline do not prevent writing the maximum consensus. This is
reference-guided mapping; de novo commands do not accept this flag or a length
prior. `candidate.fasta` and its strict comparison remain available.

For an anchored window containing N, first search sample FASTQ graphs at
k=21/17/13/11, up to `--local-k`, for a path at the current backbone length that
preserves every called A/C/G/T position. Each edge needs at least two named
pairs. A length-layered dynamic program maximizes summed `log(1 + support)`;
keep the highest-scoring spelling per node and position, with up to 1,024
nodes per layer. Length and the called positions constrain the walk. There is
no fixed repeat-visit ceiling, so six copies are not discarded just because a
strict traversal allows four. Merged histories and beam pruning are reported.
The length is a reference-guided prior, not a copy-number measurement.

In deep mode, pool these fills with the legacy sample-edge beam proposals and ordinary
local candidates. Prefer a full fill preserving called bases and span, then
paths meeting the ordinary depth threshold, longer graph k and support.
The maximum consensus applies a patch only if its length exactly matches the
backbone interval and every called A/C/G/T is preserved. If no patch qualifies,
retain the backbone, including its N calls. Nearest-reference-length proposals
remain available only in the separate diagnostic legacy FASTAs.
`maximum_consensus_rejected_patches` records rejected loci and span changes.
Each run of backbone N bases filled in the exported maximum consensus gets its
own `maximum_consensus.ambiguity.tsv` row with reason `fixed_span_graph_fill`
(alternatives 0), in addition to the window's `structural_validation_deferred`
row. These bases are sample-graph inferences whose length is the existing span,
not a measurement; annotation carries the label into NCR `sequence_origin`.
Do not change a window that already has no N just to add repeat copies.
Existing strict replacements outside these windows are retained.

The deep-mode legacy fallback uses complete main-k flanks and k=21/17/13/11 raw
pair-counted edges. At each k, support starts at `--min-depth` (at least two)
and may relax through its half/quarter values, three and two, stopping after
success. Its beam width is 256, state cap 50,000, reverse-reachability cap
100,000 and retained-path cap eight. Length is bounded by `--local-max-span`
and target plus `clamp(target,250,2000)/8 + 2*(k-1)`. The primary span fill is
independent of these repeat-use and legacy traversal limits.

In deep mode, compare up to two structural representatives, including the primary choice,
plus the unchanged baseline using the existing full-query aligner. Eligibility
is frozen against the entire generated pool. Scores, junction counts and
unknown fragments are confidence information, not output gates. All proposals
remain in local FASTA; scoring a subset does not certify uniqueness.

Outputs:

- `consensus.fasta` and `maximum_consensus.fasta`: the same best-effort primary
  sequence with the original contig name and explicit best-effort header.
  When no fill exists, retain the strict sequence, including remaining N.
- `maximum_consensus.ambiguity.tsv`: uncertainty for that primary sequence.
  The internal local-assembly rotation is retained; these are not reference
  coordinates. Supply this TSV when annotating `consensus.fasta`.
- `candidate.fasta` and `structural_ambiguity.tsv`: the separate strict result.
- `local_assembly/locus_N.legacy_bridge.fasta`: deep mode's pooled ranked alternatives.
- `local_assembly/locus_N.legacy_bridge.fragments.jsonl.gz`: deep-mode diagnostic scores;
  candidate indices are mapped in the window report.
- `legacy_bridge.fasta` and optional `legacy_bridge.read_preferred.fasta`:
  proposal and measured-read-preference comparisons, with matching ambiguity
  files. The read-preferred file requires a unique measured non-baseline best.
- `legacy_bridge.json` and local window reports: proposals, maximum selection,
  called-base preservation, traversal limits, scores and output files.

A successful best-effort run exits 0, including when repeat structure remains
unverified. Top-level `status=best_effort`, `consensus_written=true` and
`sequence_output.fasta=consensus.fasta` describe availability.
`strict_validation_passed` and the retained `final_published` field describe the
old strict validation independently. A complete strict result in deep mode may still have
`status=complete`. Pipeline failures remain errors and remove `consensus.fasta`.

Every bridge base follows a sample k-mer. Incidental short-k matches can still
join unrelated reads, and the reference-guided span can imply a wrong repeat
count. Uncertainty accompanies the consensus; annotation or absence of N does
not imply validated repeat structure or circularity.

## CLI

```bash
mitovarr assemble-mapping --ref reference.fasta --bam mt.bam \
  --local-assembly --r1 reads.R1.fastq.gz --r2 reads.R2.fastq.gz \
  --emit-variants --emit-coverage --out mapping

mitovarr annotate --fasta mapping/candidate.fasta --species cbr \
  --topology circular --rotate-to-trna-pro \
  --assembly-ambiguity mapping/structural_ambiguity.tsv --out annotation
```

In mapping assembly, omitting `--local-assembly`
preserves ordinary indexed-BAM pileup. Local assembly requires matching paired
FASTQ; the complete extracted pairs from `fetch-caendr` are suitable. No extra
recruitment or network download occurs. `--topology circular` is an annotation
hypothesis and does not establish read-supported closure.

Defaults: `--local-k 31`, `--local-max-span 2000`, `--local-max-paths 128`,
`--deep-validation` off, `--scoring-threads 4` (used only in deep mode, capped by available CPUs). Fragment scoring retains all
reads and reduces results in input order; `--scoring-threads 1` is serial.
`--max-records` also bounds unique FASTQ pair names in this mode. The sample graph
has a configurable `--local-max-kmers` cap (default 2,000,000); each local search
has a 200,000-state and four-visits-per-node cap. A reached search cap leaves the region unresolved. Increasing a cap
does not add sequencing evidence.

If a local graph exceeds its distinct-kmer or named-FASTQ-pair cap after the BAM
pileup is complete, discard the incomplete graph and retain that full-BAM
`candidate.fasta`, base votes and requested coverage/variant reports. Fast mode
or explicit `--best-effort` also writes that baseline as the maximum consensus
and `consensus.fasta`. The local report records `status=resource_limited`, the
observed value and limit, `fastq_scan_complete=false`, and null uncomputed graph
and fragment counts. Structural uncertainty covers the retained sequence and
strict validation remains false. This does not subsample or imply validation
against all FASTQ. Malformed input and I/O errors still fail; a BAM record cap
before pileup completion cannot use this partial-output fallback.

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
relaxed. Without `--subsample`, exceeding the BAM record cap remains an error.

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

Every graph kmer comes from sample FASTQ and is counted once per named
pair. Hard base-quality filtering is disabled by default (`--min-baseq 0`),
including local graphs, shared-cohort anchors and matching-base support.
A positive explicit cutoff remains available. Alignment seeds have no hidden
quality cutoff; quality-dependent alignment scores and BAM vote weights remain.
Edges need at least three fragments. Outgoing branches below 2% of the
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
Freeze one common fragment cohort using a 21-mer passing the optional base-quality cutoff on either mate
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
does not jointly phase distant structural variants or globally recompute the
consensus of homologous copies after several separate replacements. A new
source-copy ambiguity outside the discovered windows may therefore require a
further local analysis. Detection is conditional on the reads supplied; no claim
of exhaustive SV discovery is made.

## Outputs and coordinates

### Coverage-guided layout hypotheses

`assemble-mapping --local-assembly --deep-validation --coverage-layout` additionally explores
unresolved windows. This does not change ordinary pileup or strict replacement
selection. Reference coordinates locate windows; complete supported sample
flanks and paired FASTQ supply graph sequence. There is no reference-filled gap.

Adding `--multi-k-layout` retains those legacy candidates and adds the shared
cross-k generator described in [de novo layouts](de_novo.md#additional-cross-k-proposals).
One cached FASTQ pass counts named-pair support at k=15/21/31/47/63/95/127 for all
windows. The original local graph and strict replacement search remain unchanged.
Each additional graph uses the same extended left/right flank context (up to
126 bases, bounded by the current linearized backbone ends). Both original
local-k flanks remain in the replacement. Unchanged padding is removed only after
complete-flank traversal; a path overlapping away that padding is not a valid
patch for this window. Missing/unsupported flank context is reported, not filled.

The extra searches share one fixed occurrence-depth model at `--local-k`, count
cross-edge windows from the spelled suffix, and account for other backbone copies.
Per k/support they use 50,000 additional states split 5,000/22,500/22,500 across
nearby-exit, coverage and continuous-read-context schedules, with 1,024 pending
states per queue. Each guided schedule permits `ceil(local-max-paths/2)`
completions; the seed permits one. The span bound adds only the unchanged padding
to `--local-max-span`. Graphs retain a two-million-entry limit per k. Legacy search
costs remain in addition to these budgets.

All k pools merge by sequence before shared-flank recruitment and the existing
eight-path paired comparison. `source_graphs`, `multi_k_layout` and `searches`
record sources and limits; extra GFA files are
`local_assembly/locus_*.multik_k*_min*.gfa`. Agreement across k values is not
independent read support. Context is a proposal-order heuristic, and this option
does not verify phase across windows or change strict final publication.

For each window, raw pair-supported graphs at the local k and smaller k=21/15
contribute a union of complete paths. A nearby-exit search retains a short
connection before a second search schedules paths by incremental coverage-copy
residual. Every path must spell the entire original flanks, including in smaller
k graphs. The raw proposal graph bypasses strict relative edge cleanup; edges
need three named pairs, with an explicitly recorded one-pair fallback only when
the ordinary proposal pool is empty. This weak fallback can introduce incorrect
connections and never supplies evidence for final publication.

Depth counts R1 k-mer window occurrences, retaining repeated offsets and one
record per named pair. R2 is omitted to avoid overlapping-mate inflation.
The baseline is the arithmetic mean depth of unique k-mers within 500 bases of
the two flanks, falling back to unique backbone k-mers if too few are available.
Local and whole-backbone means are reported. Other homologous copies remain in
the whole-contig count model. Candidate counts are compared over one fixed union
of varying k-mers at the main local k, even for paths generated at smaller k.
Unobserved R1 k-mers and baseline sensitivity are diagnostics, not proof that
R2 lacks support or that a copy count is identifiable.

Each graph search has two schedules, each capped at 25,000 states and 512 pending
histories, and a 200,000-node reverse-reachability bound. Length and completion
bounds use `--local-max-span` and `--local-max-paths`. Coverage exploration caps
the local use of each canonical k-mer at a generous depth-derived count (up to
4096), with a minimum of two. This is an exploration heuristic, not a confidence
interval. A truncated reverse search does not classify unseen nodes as
unreachable. Limits in one graph never discard paths from another graph.

Up to eight sample paths plus the original pileup baseline are compared by the
whole-contig paired-fragment scorer. Preserve the coverage incumbent, then
represent repeat-block count families, distinct traversal orders, and additional
sequence variants. Report unscored families when the budget is too small.
Compared paths use a unique-anchor fragment cohort frozen on the full generated
pool and baseline before shortlisting; read scores are retained even
though graph search is incomplete. Coverage is not added as a second independent
likelihood from the same reads. The report keeps baseline preference and
unscored/exhausted alignments explicit. Its nested subset selection can pass
local scoring checks but does not resolve the full competitor domain.

A final bounded combination step allocates coverage across all inferred windows
at once: the same repeat observations cannot fund extra copies independently at
two loci. Choices within five log-score units of the best scored path in
each window remain eligible, including the pileup baseline even when graph
search has not regenerated it or coverage excludes it from the shortlist.
With `u` omitted fragments, expand this five-unit margin by
`u * (12 + ln(1 + exp(-12)))`, the maximum possible change in score contrast.
A path outside that margin stays excluded even when some fragments are unknown.
With no measured fragments, all shortlisted paths and the baseline remain eligible. The common occurrence model uses
the arithmetic mean local baseline across windows varying each k-mer. Strict
replacements and other backbone copies stay fixed in that model. An optimistic
remaining-copy interval schedules at most 512 configurations; up to 32 full
layouts are exported. Local paired scores are not summed across windows, and
this step does not establish phase between windows. Equal-scoring retained
configurations remain explicit.

Additional outputs:

| File | Meaning |
| --- | --- |
| `coverage_layout.fasta` | One inferred whole-contig hypothesis; `final_verified=false` and `circularity_verified=false`. |
| `coverage_layout_candidates.fasta` | Up to 32 retained combinations in joint coverage order. |
| `coverage_layout.json` | Selected windows, their searches and paired scores, global copy residual, configuration ties and limits. |
| `coverage_layout_ambiguity.tsv` | Unresolved windows projected into this hypothesis, accepted by `annotate --assembly-ambiguity`. |
| `local_assembly/locus_N.coverage.fasta` | All generated sample paths for that window, in coverage order. |

The layout keeps the internal local-search rotation; its coordinates differ
from `candidate.fasta` and the reference. `baseline_rotation_offset` records
that cut. Annotation may explicitly assume `--topology circular` and rotate
the hypothesis after selection, but neither operation verifies closure.
No selected layout FASTA is written if no unresolved window yields a changed
choice. `baseline_retained` means proposed paths were compared but the original
local sequence was kept; their window FASTAs and scores remain available.
Repeat sampling bias, low-k false connections, missing extracted reads and
unsearched configurations remain limitations even with a low residual or no N.

### Strict mapping outputs

| File | Meaning |
| --- | --- |
| `reference_candidate.fasta` | Frozen ordinary pileup consensus, including its accepted CIGAR indels and N calls. |
| `candidate.fasta` | Pileup backbone with supported sample-path replacements. |
| `consensus.fasta` | Published only when every discovered local window is resolved, remaining repeat ambiguity is absent, and the sequence has no N. Circularity remains unverified. |
| `local_assembly/locus_N.fasta` | Baseline candidate 0 and all generated local alternative sequences, with flanks. |
| `local_assembly.gfa` | Local candidate-path graph with shared flank nodes and k-bp overlaps; it is not the complete underlying kmer graph. |
| `structural_variants.tsv` | One row per evaluated window, selected candidate, length delta and unresolved reason. `ref_window_start/end` use original reference coordinates and include flanks, not a minimal breakpoint. |
| `structural_ambiguity.tsv` | Unresolved intervals in final `candidate.fasta` coordinates; accepted by native `annotate --assembly-ambiguity`. |
| `local_assembly/report.json` | Library model, same-cohort pair scores, junction support, search limits, baseline/final coordinate segments and selected minimal edits. |

Existing `variants.tsv`, `coverage.tsv`, `ambiguity.tsv`, `uncertainty.tsv` and
`base_votes.tsv` continue to describe the **baseline reference-coordinate pileup**.
They are not silently relabelled as coverage or alleles on inserted sequence.
Local mapping reports retain `mitovar-rust.polish.v2`; baseline mapping reports
retain `mitovar-rust.polish.v1`. These are historical format identifiers kept
for report compatibility; the command is `assemble-mapping`.
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

Whole-candidate evaluation reports seeded-decoy mate rescue, continuously
observed junction fragments and mate brackets separately. Contraction flanks
require actual matching traceback bases in the same placement. See
[shared layout evidence](paired_paths.md#shared-whole-layout-read-evidence).

Common scoring anchors come from actually retained flanks, including up to 500
bases of adjoining backbone, rather than the union of changed candidate interiors.
The same query offset and target word constrain rescue and joint scoring.
Compressed fragment diagnostics accompany optional coverage comparisons.

Joint coverage selection reports
`selection=joint_copy_depth_among_retained_local_hypotheses` and
`pair_compatibility_verified=false`. The reported baseline margin, unknown count
and `choice_retention_log_margin` distinguish measured opposition from differences
that omitted fragments could reverse. The bound does not validate structure or
unsearched candidates. Zero unknowns means completed comparison inside the
reported alignment domain, not exhaustive assembly search.
