# Changelog

## Unreleased

- Label tandem repeats in the best-effort mapping consensus whose length reads
  do not measure. Blocks with unique flanking anchors are measured only by reads
  carrying both anchors; others get a `tandem_repeat_length_unmeasured` (or
  `_conflict`) row in `maximum_consensus.ambiguity.tsv`, and
  `report.json.tandem_repeats` gives a read-derived lower bound. Annotation
  reports the locus `evidence_status`. No copy number is inferred from depth;
  sequences and other outputs are unchanged.

- Speed up de novo graph cleanup without changing results: bubble searches
  keep per-round successor/twin lists and mark the current walk instead of
  cloning and rescanning paths; the edit-distance check reuses two rows;
  compaction stores one outgoing run per node; removed k-mers use the k-mer
  hasher. The prior search and edit-distance kernels remain as differential
  test oracles. All 27 cel/ctro panel strains are byte-identical; ctro QG5638
  at an 8M cap drops from 800 to 207 s and cel NIC262 from 438 to 143 s.

- Add opt-in de novo `--subsample` / `--subsample-seed`. Only after all complete
  pairs exceed `--max-kmers`, keep seeded fragment-name-priority pairs at 1/2,
  1/4 … 1/64 until the graph fits, write them as `subsampled.R*.fastq` and
  report the trigger, attempts and counts in `report.json.subsampling`. Final
  read validation still uses all pairs; without the option an overflow remains
  an error.

- Label each mapping maximum-consensus N run filled at the existing span as a
  `fixed_span_graph_fill` row in `maximum_consensus.ambiguity.tsv`. Annotation
  reports it as NCR `sequence_origin` and locus `evidence_status`, so a filled
  NCR1 (e.g. the cel read-dropout core) is not read as a measured length.
  Sequences and other outputs are unchanged.

- Restore tRNA/rRNA termini clipped by one or two called sample bases next to a
  supported aligned end, labelled `terminal_mismatch` (e.g. the CaenDR
  *C. tropicalis* tRNA-Trp first base). Missing or ambiguous ends stay partial;
  the tRNA-Pro rotation origin still requires an `aligned` 5' boundary.

- Recover supported partial primary paths across local alleles and competing
  reciprocal views. Report partial-primary allele evidence separately, retain
  unresolved conflicts, and preserve strict join/publication criteria.

- Use internal numeric-key hashing, shared reverse target-label propagation and
  constant-time deleted-unitig/twin checks in de novo graph work, preserving
  candidate ordering and per-target budgets without adding dependencies.

- Skip read scoring for paths/cycle families that cannot be selected; retain
  candidate exports and report unmeasured evidence as null with `read_scoring`,
  including fragment scores nested in repeat-copy diagnostics.

- Remove weak disconnected graph shadows using strand-symmetric sample support;
  require 4k bases for isolated path anchors while preserving attached terminal
  anchors. Report shadow removal and post-cleanup graph counts without changing
  default base quality, input scope or complete-publication gates.

- Reduce seeded alignment and mate-rescue work with in-place affine score rows,
  cached per-quality emissions and a lossless disjoint-block rejection check.
  Reuse per-fragment k-mer sets and use one count lookup in main/end graph builds.
  Preserve full input, alignment scores and tie order, candidate searches,
  evidence limits and publication decisions; retain the prior fitting kernel
  as a differential test oracle.
- Recover sample-only path anchors and shared fragment recruitment using adaptive
  31/47/63-base contexts when short duplicated graph fragments obscure uniqueness.
  Preserve the main graph k, input scope and spanning/structural evidence gates.
- Enforce called-base and span preservation at the mapping maximum-consensus
  output boundary, including deep best-effort mode; retain unsafe alternatives
  only in diagnostic files and report rejected patches.
- Retain the complete main de novo graph when optional endpoint metadata exceeds
  its limit, withholding complete publication. Add `--local-max-kmers` (default
  2,000,000) and retain full-BAM partial consensus when local graph work reaches
  a resource cap. Report incomplete FASTQ validation explicitly; input errors
  and incomplete BAM pileups still fail.

## 0.1.0 - 2026-10-10

### Package and interfaces

- Add tag-driven GitHub binary releases for macOS Apple Silicon, static Linux
  x86_64 and Linux ARM64, gated by reusable native-platform CI and unpacked
  archive smoke checks. Package usage documentation with each executable and
  document installation and the maintainer release procedure. The CaenDR shell
  wrapper uses the packaged binary when available, without requiring Cargo.

- Make local mapping output a fast best-effort consensus by default, using all
  supplied pairs and preserving called bases. Only fill N at the existing span;
  mark structural scoring as deferred with null scores and fragment counts.
  Add `--deep-validation` for the prior full search/scoring workflow and require
  it for mapping coverage layouts. `--best-effort` remains available in deep
  mode and is implicit in fast local mode. Fast output never passes strict
  structural validation or claims verified circularity.

- Certify unrestricted affine endpoint optima against the existing edit and
  traceback domain, falling back to bounded profiles for uncertified endpoints.
  Keep two score rows and compact traceback predecessors. Add
  `assemble-mapping --scoring-threads` (default four, capped by available CPUs)
  with persistent workers, bounded in-flight fragments and deterministic
  score/diagnostic reduction.
  All input reads, candidate scores, logical limits and uncertainty rules remain.

- Aggregate whole-layout mate combinations with endpoint range indexes while
  preserving per-span scores, floating-point ties, junction/bracket support
  and contraction evidence. Reuse identical per-start affine subproblems
  across candidate windows and omit provably dominated internal edit states,
  charging the original logical work to retain
  alignment limits and unknown states. Input scope and output schemas are unchanged.

- Isolate small-k end searches by target within the existing total source/k
  budgets. Report target-specific reachability and limits, and preserve all
  discovered spellings in `end_bridge_proposals.fasta`. Keep admitted graph
  connectors in `end_bridge_candidates.fasta`; incomplete competing searches
  still prevent new joins from the affected source and its reciprocal view.
  Bound endpoint-pair metadata to 10,000 combinations and report an explicit
  error before allocation when a highly fragmented graph exceeds the limit.

- Present de novo assembly through `assemble-baited` and `assemble-seeded`
  only. Hide the unused direct `assemble-denovo` entry for compatibility with
  existing commands; remove it from current guides, help and workflow figures.
  Preserve the two active input/validation scopes and shared assembly engine.

- Standardize the two workflow names as de novo assembly and mapping assembly.
  Make `assemble-mapping` the canonical CLI command, with `polish` retained
  only as a hidden compatibility alias. Rename the Rust entry point to
  `mapping_assembly::assemble`, its CLI variant to `AssembleMapping` and its
  arguments to `MappingArgs`. Update help, documentation, examples and tests;
  keep historical report schema identifiers and assembly behavior unchanged.

- Expand the documentation entry point with workflow selection, the three
  mapping modes, runnable examples, matching FASTA/uncertainty outputs and
  independent annotation. Explain the current `--local-assembly` requirement
  for `--best-effort` and distinguish consensus availability from validation.

- Rename the mapping maximum-consensus option to `--best-effort` and describe
  it by its output behavior. Keep `--legacy-bridge` as a hidden compatibility
  alias for existing commands.

- Select an explicit `best_circular.fasta` and matching uncertainty table from
  eligible de novo circular representatives. Preserve within-family fragment
  scoring, rank components by distinct k-mer support, retain deterministic ties
  and report selection/provenance separately from verified completion. Keep the
  existing primary sequence, strict final gates and exit codes unchanged.

- Align the README and figure labels with bounded layout scoring, component-wide uncertainty and independent gene-quality checks. Clarify the common-backbone phasing scope, distinguish the overview from the detailed method schematic, and add the verified GitHub clone URL and package repository metadata. Assembly algorithms and output contracts are unchanged.

- Support partial results without mandatory NCR1 completion. Add `sequence_output` availability metadata and independent per-model `gene_quality.tsv` / `records[].gene_quality` checks for missing, incomplete, duplicated, uncertain or assembly-overlapping gene placements. Preserve assembly algorithms, sequence choices, exit codes and strict complete-output claims. Annotation quality does not certify read support.

- Add opt-in `--multi-k-layout` alongside `--coverage-layout` for de novo and mapping hypotheses. Retain legacy candidates, add long-flank searches through k127, schedule low-k branches with longer observed read contexts, count full spelled depth windows, merge provenance and compare the union on shared fragments. Report added budgets and unsearched domains without changing strict primary/final output.

- Align the README with implemented FASTQ contracts, mapping resource limits and mixed-site reporting; distinguish reference, candidate and annotation coordinates and explain topology handling for unmarked final FASTA headers.
- Embed a five-panel scientific Figure 1 in the README, with journal-sized editable PDF/SVG, a 600-dpi PNG, a manuscript legend and a reproducible drawing source. Explain input scopes, alternative repeat structures, shared-fragment evidence, optional multi-k/depth hypotheses, conditional outputs and independent annotation. Mark all examples as schematic and keep inferred layouts distinct from strict final consensus.
- Remove navigation links to the retired Chinese README; keep the Chinese annotation guide.
- Exclude local AI-assistant files, TODOs, maintainer notes and exploratory-test directories from Git; limit the source package to runtime code, active regression tests and user documentation.
- Rename the package, library, executable and launcher to `mitovarr`. Preserve existing `mitovar-rust.*` JSON schema identifiers and gene panels; use `mitovarr` as the GFF3 source.
- Provide an independently buildable Rust package, English and Chinese usage guides, CLI documentation, method descriptions and editable workflow Figure 1. Keep experiment records and analysis outputs outside the repository and source package.
- Remove unused legacy pileup, variant-refresh and repeat-summary helpers and redundant tests; retain coverage through the production fragment accumulator and portable regression fixtures.

### Assembly and consensus

- Make mapping `--best-effort` retain a maximum usable consensus. Add
  sample-only masked-span walks that preserve called bases without a fixed
  repeat-visit ceiling; prefer them to shortened legacy paths. Write
  `consensus.fasta` with a matching uncertainty TSV and successful exit even
  when strict validation is unresolved. Report consensus availability and
  strict validation separately; compare two representative paths plus baseline
  for confidence without using that comparison as an output veto.

- Add opt-in mapping `--best-effort`: reference-span-guided sample paths at
  k=21/17/13/11 with support relaxation, bounded beam search, full sample flanks
  and shared-fragment scoring even after strict graph limits. Export local
  alternatives, length-guided and measured-read-preferred whole hypotheses,
  uncertainty coordinates and explicit limits without promoting inferred
  repeat structures or changing default assembly outputs.

- Compare whole hypotheses with endpoint/edit-aware affine profiles and an in-DP common pivot; bind scores and junction observations to the same traceback, suppress tied-trace junction claims, and distinguish complete bounded absence from shared unknowns. Preserve library calibration and primary-validation behavior.
- Allocate the existing eight sample-path comparison slots across repeat-block count and order families before sequence variants, retaining the coverage incumbent and extra baseline. Freeze common-flank eligibility across all generated paths, export the entire de novo candidate pool with merged provenance, and report read-preferred ranks and unscored families.
- Bound the score contrast from omitted fragments when retaining mapping coverage choices, so unknowns no longer disable all rejection. Inferred layouts remain separate from final assemblies; graph-generation limits and publication requirements are unchanged.

- Correct joint coverage selection metadata from `joint_copy_depth_among_local_pair_compatible_choices` to `joint_copy_depth_among_retained_local_hypotheses`, with `pair_compatibility_verified=false`. Add `alignment_search_exhaustive=false` to evaluated whole-layout comparisons and document endpoint/traceback pruning. This reporting correction does not change scores, selected sequences or publication gates.

- Freeze shared flank identity and query offset throughout whole-candidate scoring and rescue; remove sequence-only rescue skipping, preserve complete bounded negative versus unknown search states, and stream span-level diagnostics. Exclude unknown bases from matching traceback evidence and bind validation to the actual coverage representative. Keep candidate generation, coverage-first eight and publication thresholds unchanged.

- Share bounded whole-hypothesis mate rescue and traceback-based junction/contraction evidence between mapping and optional de novo layouts. Rescue existing seeded decoys across circular origins, distinguish observed junctions from mate brackets, and validate sole coverage hypotheses without promoting them to final.

- Extend opt-in `--coverage-layout` to reference-guided local mapping. Retain nearby exits and bounded coverage-guided paths across k values, compare a fixed coverage shortlist with paired reads even when strict searches are capped, and allocate repeated k-mer observations once across multiple windows. Export hypotheses, alternative configurations, paired scores and explicit limits without changing strict candidate/final selection.
- Correct zero-based inclusive repeat-interval projection during local mapping discovery and remaining-ambiguity checks, preventing origin underflow and retaining the final base after pileup indels. Add boundary, indel and missing-origin-coverage regressions.
- Add opt-in `--coverage-layout`: fit repeat usage to arithmetic mean unique-flank R1 occurrence depth, retain multi-k sample-graph layouts and coverage sensitivity, report singleton-edge fallback when ordinary support filtering disconnects the corridor, and export inferred FASTA, GFA walks and ambiguity/depth reports independently of strict primary/final selection.
- Require actual matching traceback bases for complete-cycle phase calls. An empty allele requires a sequenced contiguous junction; deletion spans and inserted query bases cannot supply that observation.

- Schedule anchor-path walks by accumulated intermediate bases so complex repeats do not consume every candidate before short exits are visited. Preserve full histories, fragment-based selection and incomplete-search vetoes; bound and report pending search states.

- Prune smaller-k end-search branches that cannot reach any target suffix; preserve full-flank matching, reciprocal spelling, target-reaching repeat limits and final publication checks. Report excluded branches separately from capped searches.
- Disable hard base-quality filtering by default (`--min-baseq 0`) across recruitment, de novo and mapping, and remove the aligner's implicit Q10 seed cutoff. Preserve quality-dependent alignment scores and mapping vote weights, explicit cutoff overrides, graph support limits and structural publication checks.
- Preserve complete flank context during smaller-k mapping traversal so repeated short flank sequences do not discard supported paths; retain graph-edge support checks and search limits.
- Preserve permissively collapsed bubble alternatives and block final publication in their original graph component; prevent mixtures from becoming falsely unique through cleanup while retaining longer primary sequences.

- Compare bounded de novo bubbles using whole-path mean support and longer local alternatives, with independent `--max-bubble-fraction` (default 0.49); retain strict repeat-route thresholds.
- Add pooled k=21/15 sample-graph fallback proposals for mapping windows with no main-k path, preserving complete original flanks, search limits and paired-fragment replacement gates.

- Prevent local contractions from being published on separate repeat-internal junction evidence alone; require three named fragments linking unique unchanged flanks and report this evidence explicitly.

- Add native CaenDR extraction through indexed HTTP Range requests, local mate recovery and paired FASTQ/BAM output, with a portable launcher.
- Keep seed-recruited FASTQ, BAM-baited de novo and reference-guided indexed-BAM workflows separate. Runtime processing uses in-process Rust implementations.
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
