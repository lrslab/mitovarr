# Paired-path assembly

The assembler builds a sample-only main graph (default k=63; odd k=15–127), proposes end connections at smaller k values, and compares complete paths between anchor segments using paired FASTQ. The graph, library model, candidate generation and selection never receive a reference, species identity or expected genome length. `assemble-baited` uses the already extracted pairs; neither model fitting nor scoring recruits additional reads. Each run has one main k; paths are not automatically pooled across main-k runs.

## Alignment and fragment model

The native aligner consumes the entire read. Reference-window flanks are free, but query ends cannot be clipped for free. Its three affine-gap states allow substitutions and indels, with log-score penalties 8 for opening a gap and 2 for each additional base. Mismatch scores use Phred qualities, with the nominal error probability bounded to 0.001–0.75. These are scoring assumptions, not a calibrated assembly-error model.

Seventeen-base seeds locate candidate windows in either orientation. Reverse reads also reverse their quality arrays. The dynamic program retains only diagonals compatible with the accepted edit bound (10% of read length, at least three edits). Seed starts are bounded at 256 per strand. When a fragment reaches an alignment or pivot limit, it remains an unknown contribution shared by every competing path. Its maximum possible contrast is bounded by 12 + log(1 + exp(−12)), from the score range. Strict final-path selection uses the worst-case lower margin and possible opposing votes. The cluster norm also includes the bound by the triangle inequality; uniform omission alone is not treated as evidence of zero contrast. A per-fragment cache reuses identical target-window placement profiles across candidates; it includes the read and its qualities by construction and is not reused across different fragments.

An exploratory set of uniquely placed read pairs estimates orientation and a candidate span range. The measured fragment span is the outer extent of the two alignments, including overlap only once. FR/RF orientation uses the two 5′ ends, not a sort of alignment left coordinates.

Final calibration uses a fixed R1 cohort in internal sample sequence. R1 must have a unique full-read placement and enough continuous sequence in the expected mate direction to cover the exploratory range. R2 does not have to map for its fragment to enter the cohort. R2 gets full-query gapped rescue when unique sample seeds identify its region. The empirical span mass uses 21-bp smoothing and retains unmeasured R2 in its denominator. An insufficient calibration returns `not_estimable`; there is no default insert mean.

Calibration compares full-read placements on seed-nominated segments; an isolated seed on a short competing segment does not itself make the read ambiguous. Quantiles describe the measured subset, and the histogram represents measured span mass rather than a missingness-corrected population distribution. The same automatic policy makes distance diagnostic-only if more than 5% of qualifying R1 fragments are unmeasured or the central spread of internal-window medians exceeds 25% of the measured median (windows need 100 measurements, with at least three windows for this diagnostic). These are engineering eligibility conditions, not calibration guarantees.

The model is conditional on these FASTQ and anchor locations. It cannot undo BAM-extraction bias or recover identities of libraries pooled into one FASTQ pair. Its internal-window summaries and unmeasured fraction are reported. Named fragments are not claimed to be PCR-independent molecules.

## Complete path hypotheses

Read-threaded unitig instances retain their entrance/exit restrictions when converted to a segment graph. This avoids reintroducing connections that exact repeat threading already resolved. Internal anchors require at least 4k bases. A terminal source or sink can instead qualify with k−1 overlap bases plus 20 flanking bases, so a short terminal segment behind a merge is not excluded solely by the internal-anchor length rule. Both kinds require at least k unique canonical positions at one seed length: prefer 31, then try 47 and 63 up to the main k (with 31 retained for smaller main k). Short duplicated graph fragments therefore cannot by themselves disqualify a long unique context. Circular segments are excluded. Reverse-complement sources and sinks follow the same rule. These are candidate anchors only: a join still needs the full-fragment spanning evidence below. Uniqueness in this graph is a localization aid, not proof of biological single-copy status.

The search follows a complete path to the next anchor. Simple cyclic components are enumerated by explicit integer turn count and exit position; complex components retain bounded graph traversal. Sequence is concatenated with the same k−1 overlap used by the graph. A cycle removes exactly one terminal copy of that overlap. No bases are filled from mate distance.

The traversal schedules pending walks by increasing accumulated intermediate bases, with deterministic full-history/node tie breaks. This prevents depth-first repeat expansion from consuming the candidate budget before nearby exits are visited. Different histories arriving at the same node and length remain separate. This is a discovery order, not a biological length prior or a fragment-score term; a truncated set can still omit a better-supported longer path.

Each anchor search currently has a budget of 128 paths, 50,000 visited states and 50,000 pending states. `peak_pending_states` and `frontier_limited` expose the pending-state bound; hitting it also sets `search_exhausted` and prevents selection. `search_order` and `maximum_pending_states_per_anchor` document the traversal and its bound. Per-locus `search_cuts` separately records the visited-state, candidate-path, cycle pre-expansion, pending-state and intermediate-length limits. Simple cycles still cap their generated exits before enqueueing, so this traversal does not promise the globally shortest 128 complete walks. Every such cap retains the existing incomplete-search veto. Its length domain derives from the sample fragment range and read length. `search_exhausted` and `length_limited` remain explicit; an open repeat tail does not establish a unique repeat count. All generated local path spellings are saved in `path_candidates.fasta`, including alternatives that are not selected.

End search uses the same FASTQ at k=47/31 below the main k, with the same optional base-quality cutoff (disabled by default), support thresholds and full main-k overlaps at both flanks. Each reachable target gets a separate search within the total source/k budgets of 128 paths and 50,000 states. Same-end reverse-complement proposals retain their own target budget; they are not globally forbidden. Reciprocal spellings must agree, and any limited competing target keeps the affected source ineligible for new joins. `end_bridge_proposals.fasta` preserves discovered spellings, including withheld proposals. `end_bridge_candidates.fasta` contains only admitted graph connectors, and `path_graph.gfa` exports the actual segment graph. See [end connections and simple repeat copies](de_novo.md#end-connections-and-simple-repeat-copies).

Per-locus `repeat_copy_model` records cycle period, full turns, partial traversal via the complete connecting length, fragment scores, and occurrence-depth predictions. It counts repeated R1 windows passing the optional base-quality cutoff rather than pair-deduplicated graph k-mers; R2 is excluded to avoid overlap inflation. Depth predictions are diagnostics under local sampling assumptions, not another independent score or a hard bound on the unsearched tail.

## A fragment contributes one joint explanation

Candidate membership is frozen before scoring. Recruitment retains unique 31-mers in each source terminal region. When fewer than k positions are unique, it adds one longer context: the shortest 47/63 length with at least k unique positions, or the length recovering the most positions if neither reaches k. This choice uses graph sequence alone. `anchor_contexts` records each source’s anchor length and recruitment seed counts. Short reads can still use the retained 31-mers, and multiple lengths never count the same named fragment twice. A unique seed on either mate in the common source anchor admits a fragment independently of the other mate’s placement. R1 is used when both qualify, chosen before alignment and geometry scoring. The same fragment names are used for every alternative, once per pair. The chosen mate’s pivot is retained, so another copy of a flank cannot substitute for it. The FR/RF and outer-span model is symmetric under swapping the working mate order; calibration itself remains the fixed R1 cohort described above.

For each allowed orientation/span state, the scorer retains the best full-fragment sequence explanation. It then sums over the **fixed empirical geometry mass**. It neither sums duplicate placements nor renormalizes by the number of hits a candidate happened to produce. A fixed log background of −12 keeps unexplained fragments in the comparison. Reported scores subtract each fragment’s best score across every hypothesis. The candidate-independent sum of these offsets is reported, so the contrast can be distinguished from the unshifted score. They are compatibility scores, not assembly correctness probabilities. If geometry is diagnostic-only, the fragment retains a joint sequence explanation and known orientation without a distance preference. If the non-anchor mate is unplaceable across the entire candidate set, the same anchor-only mode is used for all hypotheses; this never becomes a candidate-specific exemption.

The strict decision used for final publication requires a score margin of at least 10, stability relative to fragment score differences, at least three spanning named fragments, and no substantial coherent competing path under the existing primary-consensus ratio. A read can span both sides directly, or its mate can supply the opposite flank under the empirical geometry model. Unresolved alternatives remain available. A high local k-mer count is never substituted for this joint evidence.

For linear alternatives joining the same anchor pair, with equal lengths and exactly one differing base, spanning support covers that SNP plus 20 bases on each side. The shared k−1 graph overlaps remain in the sequence and graph coordinates, but do not enlarge the SNP evidence interval. Fragment membership, full-query alignment, whole-fragment scores and support thresholds are unchanged. Indels, multiple differing sites, competing exits and cycle/repeat hypotheses keep the graph corridor rules. Per-path `spanning_interval` records its zero-based start, exclusive end and `single_snp_shared_flanks` or `graph_corridor` basis; `connecting_bases` and FASTA anchor coordinates retain their existing meaning.

## Common-backbone cycles and phase

When graph walks return to the same oriented backbone, the assembler spells complete cycles and merges their reverse-complement views. No expected genome length is supplied. Traversal length schedules candidate discovery only; whole-cycle selection uses fragment scores and the existing completeness gates. All cycle alternatives on that backbone are compared with a common pool admitted by seeds on either mate present exactly once in every complete circular hypothesis, including origin-spanning seeds. This includes both flanks and conserved internal sequence. Common seeds provide a candidate-specific coordinate for the same physical pivot; the other mate’s alignment success never controls inclusion. Uniqueness within the generated hypotheses does not establish biological single-copy status. Both reads are aligned against the complete circle, represented in unwrapped coordinates around the fixed anchor-mate pivot; accepted joint placements span at most one circumference. Other repeat copies outside the local corridor remain competing placements. Once this comparison applies, its decision supersedes earlier one-sided local cycle choices.

Ordered shared 31-mers delimit variable sequence blocks. Changes in the distance between overlapping shared anchors also define a block, so homopolymer indels cannot disappear into anchor overlap. Touching intervals are merged. The conserved sequence between all blocks must be identical across hypotheses before phase evidence can support selection; every candidate difference must be represented. The full fragment score is still accumulated once per cycle hypothesis. A block allele is called only when its best explanation exceeds every other allele by 2.3 log units and a read spans that block with flanking sequence. An ambiguous joint placement cannot pool allele observations from its different locations. Adjacent blocks must have compatible calls on the same joint placement in at least three distinct anchor-mate/start/strand clusters; high independent marginal allele scores alone cannot establish phase. Boundary and block bridging support also need three start clusters. These are operational consensus criteria, not a proof that an entire biological haplotype is unique.

Scores are compared on all candidates, with a margin of 10 and three times the root sum of squared start-cluster score differences. A substantial strong minority retains ambiguity. The root-sum-of-squares condition measures concentration of support; it is not a three-sigma test or a calibrated confidence interval. Alternative outgoing graph topologies remain in the original common-source comparison; internal bubble scoring cannot erase another exit. The new complete-cycle pass currently applies to cycles on one shared backbone. It is not an exhaustive multi-backbone assembly solver.

Identical target windows reuse one query's affine-gap result across candidate comparisons. Before a candidate's absent non-anchor-mate placement becomes an exclusive penalty, a bounded seedless rescue searches the independently anchored mate range if another hypothesis maps that mate. An existing lower-scoring seeded placement also triggers rescue in windows nominated by better placements on other candidates, projected from the common anchor-mate pivot and padded for candidate length differences and accepted edits. The existing hit is retained alongside rescued placements. When geometry is diagnostic-only, its pilot span does not exclude mate placements; enumeration limits remain explicitly reported search budgets. Queries still consume their full length. This bounded rescue does not make the general placement search exhaustive: all scores remain conditional on the seeded and nominated placement domain. A better seedless placement unrelated to any nominated window can still be missed. These regressions do not establish complete real-genome or unique-haplotype recovery.

The seeded and nominated-window affine fitter uses one in-place score row and
reuses the same quality emission values across candidates. A necessary-condition
check may reject a window before fitting: a full-query alignment with at most
`e` edits must preserve at least one of `e+1` disjoint exact query blocks.
Unknown query bases and unsuitable block lengths bypass this check. It does not
filter reads, reduce the placement domain, replace seedless rescue, or change
scores, ties, limits or unknown states. The previous six-buffer implementation
is retained as a test oracle; a separate edit-distance oracle checks rejection.

An isolated segment with neither incoming nor outgoing links must be at least
4k bases long to become a path anchor. A genuine terminal source or sink attached
to another segment retains the shorter `(k−1) + 20`-base threshold and the same
sample-unique seed requirement. End-connection search uses the same anchor rule.

## Primary path output

The default run now also writes `primary.fasta`, a single usable representative of the main sample component. This is separate from the strict `consensus.fasta` qualification. An unresolved minor allele or missing phase link no longer has to break the primary backbone.

For a local path, the primary choice is the highest complete-fragment score with at least the existing minimum spanning support (default three fragments). Against each alternative it must have a positive score contrast and at least three strong supporting fragments, with more observed strong votes than the alternative. The strict minority/majority ratio and worst-case opposition assigned to omitted alignments are not primary-output vetoes. Limited alignments still contribute to no candidate and remain reported as unknown. An exhausted or length-limited path search supplies no relaxed join. This particularly prevents a repeat's arbitrarily truncated copy-number range from becoming a chosen traversal.

A physical local join must have one spelling in both graph orientations. Reciprocal views inherit the same join; an existing strict selection takes precedence over a conflicting relaxed view. If two relaxed views disagree, that join remains open. They are never counted as independent extra fragment votes. `loci[].primary_path` records these primary choices; `selected_path` retains its previous strict meaning. A primary choice is still subject to graph connectivity and unambiguous incoming connections when contigs are spelled.

For a common-backbone complete-cycle family, the primary path is the **one existing complete candidate with the highest joint fragment score**. Candidate order breaks an exact score tie deterministically. Enumeration, common-block representation and corridor boundary support must be available, and competing exits must still pass the existing comparison. However, a close runner-up, coherent minor alleles or missing adjacent-block phase does not remove this representative. `families[].primary_path`, `runner_up_path` and `runner_up_score_margin` identify it. This can be an unphased combination present in the graph: it is not a demonstrated physical haplotype or a repeat-copy-number call. No independent local-majority pieces are substituted into a complete-cycle hypothesis.

One representative per cycle family enters component selection. Distinct sample canonical k-mer support mass chooses the main component; a repeated traversal cannot gain mass by counting the same k-mer again. Mass does not choose among the cycle family's repeat lengths. No expected length, mitochondrial reference or annotation enters this choice. A strictly validated final candidate takes precedence when available. Sample support is not a mitochondrial identity classifier.

`primary_ambiguity.tsv` localizes retained alternatives and open graph ends on `primary.fasta`. Columns are `sequence, start, end, evidence, reason, alternatives`; coordinates are **one-based inclusive before tRNA-Pro rotation**. Complete-cycle variable blocks are projected through the reported origin and reverse-complement transform. A block crossing FASTA origin is split; an empty allele marks one adjacent base. Local-path intervals enclose sequence differences across the compared alternatives; they are not claimed to be minimal SV breakpoints. A `graph_end` row marks an open end, not an estimated gap length. The `evidence` field points to a `locus_N` or `circle_family_N_block_N` in `paired_path_resolution`; those entries retain path identities, scores and support. Alternative spellings remain in `path_candidates.fasta` and `complete_path_candidates.fasta`.

Permissively collapsed bubbles also add a `collapsed_bubble_component_unresolved`
interval across the affected primary sequence. This component-wide flag is not
a localized breakpoint or a claim that every base is erroneous. It remains
uncertainty when passed to gene-quality checks.

To carry these intervals into tRNA-Pro coordinates, supply the matching FASTA and `annotate --assembly-ambiguity primary_ambiguity.tsv`. The Rust annotation command writes projected TSV/GFF3 intervals in `annotated.fasta` coordinates and links them to observed genes/NCRs. Source evidence and assembly status are retained. See [the annotation interface](annotation.md#carrying-assembly-uncertainty-into-annotation).

The new `report.json.primary` object records the chosen source, sequence length, graph circularity, uncertainty intervals and `final_verified`. A literal `circular=true` means a complete graph walk is present; `final_verified=false` retains uncertainty about its structure or support. `primary` is null and no primary FASTA is written if there is no candidate of at least k bases. The existing exit codes and final-publication criteria are unchanged: an unresolved run can return 2 and still provide a useful primary sequence. Annotation can then consume this frozen `primary.fasta`; its gene and NCR completeness must not be used to select a different path.

## Best circular candidate

Every de novo run independently selects one best available eligible circle and
writes `best_circular.fasta` with `best_circular_ambiguity.tsv`. Selection uses
the same eligible representatives described above: one highest-scoring complete
cycle per qualifying family, plus existing circular graph/joined-path candidates
of at least the main k length. It does not reconsider every raw ring by length
or k-mer mass, relax family eligibility, or generate new paths.

A validated final circle takes precedence. Otherwise, distinct canonical k-mer
support mass ranks representatives from different components, with normalized
sequence order breaking an exact mass tie deterministically. Whole-cycle
fragment scores remain comparable only within their shared-fragment family;
they are not compared across unrelated families. This preserves the existing
read-based choice of repeat length within each family. No reference, annotation
or expected genome length enters selection.

`report.json.best_circular.status` is `selected` when a ring is available and
`unavailable` otherwise. A selected entry records the FASTA and ambiguity file,
source candidate, available score and runner-up margin, support mass,
uncertainty intervals and `final_verified`. The source of a graph segment may
have no paired-family score; this is not replaced with an invented score. An
unavailable entry has null file paths, and neither best-circular file is written.
`eligible_representative_count` counts these deduplicated representatives after
within-family selection, not every raw cycle enumerated in a family.

Both files use record ID `best_circular`. Coordinates are one-based inclusive
in that FASTA, before annotation/rotation. The FASTA header states `circular=true`
and `final_verified=true|false`. Use the matching pair with
`annotate --assembly-ambiguity best_circular_ambiguity.tsv` so unresolved paths
and component-wide bubble uncertainty remain visible in gene checks.

This output selects among the enumerated eligible circles; it does not establish
unique structure, mitochondrial identity or an exhaustive optimum. It can be a
small residual circular component when the main component is linear.
`primary.fasta` retains its existing main-component ranking, so a residual ring
does not displace that backbone. Strict `consensus.fasta` publication, run status
and exit codes retain their existing evidence requirements. No circle is forced
when only linear or ineligible candidates are available.

## Output and closure

Local loci and cycle families add `r2_anchored_fragments` and `r2_only_fragments`; `r1_only_fragments` retains R1-only meaning. Cycle families add `distinct_anchor_start_clusters`. The legacy `distinct_r1_start_clusters` now explicitly counts only the R1-anchored subset, while boundary/block/phase cluster counts include either anchor mate. Both-end anchoring never contributes two votes for one fragment.

`paired_path_resolution.short_terminal_anchor_count` counts eligible oriented anchors below the internal 4k length cutoff; it does not count selected joins or independent molecules.

`report.json` retains the `mitovar-rust.assembly.v1` schema and adds `paired_path_resolution`, including the fragment model, search status, local path scores, spanning support and selected connections. `path_candidates.fasta` contains local hypotheses; `complete_path_candidates.fasta` contains complete common-backbone cycle hypotheses with selection status in `complete_cycle_comparison` (FASTA headers identify them as hypotheses); `candidates.fasta` retains the sequences assembled under strict joins for compatibility; `primary.fasta` is the default downstream representative described above. A graph-circle flag alone does not mean read-validated closure. Per-path origin/orientation fields relate the unwrapped zero-based, half-open block intervals to the corresponding normalized FASTA sequence: orient the FASTA, rotate at the reported origin, then repeat that circle; the scored pivot is shifted by one circumference.

A single main circle is considered independently of unrelated residual small components. Retained canonical k-mer support mass is a reported diagnostic, not a hard eligibility threshold. A single supported circle proceeds to read validation; this is not a mitochondrial identity classifier. Competing circles are retained. Final coverage uses unique full-query gapped placements with overlapping mates collapsed per base and the original depth and base-quality thresholds. Affine traceback distinguishes observed bases, matching bases and deletion events. Only matching bases passing the optional `--min-baseq` cutoff (default 0, disabled) enter `coverage.tsv` depth; its columns remain `pos, depth, method`, with method `gapped_matching_fragment`. Mismatches, insertions, deletions and unsequenced mate gaps are not matching-base depth. Ambiguous and alignment-limited reads are counted separately in the report. Supported graph joins are recorded before output rotation; moving the FASTA origin does not undo that evidence. When no such graph evidence exists, the literal output-origin bridge remains required. Distance between mates is never counted as sequenced-base depth.

A scorer cannot recover missing sequence when the graph contains no complete path. Repeat-entry, continuation and exit observations can all belong to one repeated traversal and must not be interpreted as mutually exclusive haplotypes. References are introduced only after candidate generation and status are frozen.

## Methodological context

The design uses ideas from [Pathset Graphs](https://pmc.ncbi.nlm.nih.gov/articles/PMC3619201/) and [ALE](https://academic.oup.com/bioinformatics/article/29/4/435/199222), without implementing their full algorithms or inheriting their assumptions. The distinction between coverage and bridging evidence is discussed in [Optimal assembly for high throughput shotgun sequencing](https://pmc.ncbi.nlm.nih.gov/articles/PMC3706340/). Short-contig and spanning-pair selection biases must be considered when estimating distances; see [Improved gap size estimation for scaffolding algorithms](https://doi.org/10.1093/bioinformatics/bts441).

### Shared whole-layout read evidence

Mapping alternatives and optional de novo coverage layouts share a bounded
fragment-evidence implementation. Freeze R1 (or R2 when only R2 qualifies), its first observed common-flank 21-mer,
and that word's query offset before inspecting alignment scores. The word must
belong to retained backbone/flank sequence and be unique in every hypothesis.
Enforce that exact word at its projected target position and orientation inside
the alignment DP. Keep separate start, end, edit-count and affine-gap states;
an unconstrained sequence optimum cannot discard a legal anchored explanation. Overlapping
words check the same anchor block; a distant flank across a genuine contraction
remains structural evidence rather than an anchor-conflict veto.

Apply bounded mate rescue around the frozen pivot in every hypothesis, even
when an existing seeded placement has the best sequence score: empirical span
weight can make another explanation better jointly. Circular rescue considers windows
on both sides of the arbitrary origin and deduplicates only identical physical
start/end/strand states. Full-query profiles enumerate all feasible endpoints
with at most `max(3, floor(read_length/10))` edits inside each declared search
window. A disjoint exact-block filter rules out infeasible starts; no sequence
`best-16` cutoff or exact-hit-neighborhood shortcut applies to these profiles.
The retained traceback supplies both its score and base observations. Tied
optimal traces keep their endpoint score but cannot certify junctions.
Library calibration, strict de novo path scoring, and primary read validation
retain their existing aligner. The profile implementation applies to mapping
alternatives and optional whole de novo layout comparisons. Their bounded
completeness claims do not extend to the legacy path scorer.

Each window is limited to 8,000,000 DP edit states and 1,000,000 matrix slots;
more than 4,096 endpoint profiles also marks the search incomplete. These limits
are reported in `alignment_domain`. Incomplete searches omit that fragment from
every compared hypothesis, never score it as absence.

Identical per-start alignment subproblems can be reused even when the enclosing
candidate windows differ. The cache belongs to one read and orientation, with
fixed query and qualities; its key includes the complete local target slice,
relative pivot and edit budget. It retains at most 256 subproblems and 8,192
hits per orientation. Only completed subproblems are cached. Reuse charges the
same logical DP/feasibility work as the original computation; an insufficient
remaining budget falls back to that computation. Warm caches therefore preserve
the original incomplete-search state, partial profiles and endpoint limits.

Within a start, a lower-edit path at the same query/target coordinates and
affine state can replay every continuation of a higher-edit path. The latter
is omitted only if its score deficit exceeds the traceback tie tolerance
plus a conservative bound on floating-point rounding over both complete paths.
Equal and near-tied states are retained. Compact masks skip unreachable or
dominated states and invalidate old matrix entries without clearing the full
matrix. All original logical edit states are still charged against the work
limit. No endpoint, lower-scoring endpoint explanation or structural observation
is removed by comparing it to an explanation at different coordinates.

For sufficiently reused windows, free-end suffix edit bounds are computed once
for the entire window. Whenever prefix plus suffix edits can fit the budget,
such a suffix cannot leave the original start-relative band. Feasible bounds
and logical state counts therefore match the per-start calculation exactly;
the shared bound table has the same one-million-slot allocation ceiling.

Before enumerating edit counts, a three-state affine DP finds unrestricted
endpoint optima within the same coordinate band and pivot constraint. An
optimum is accepted only when its path obeys the edit budget and no discarded
cross-edit competitor can affect its score, traceback choice or ambiguity
under the conservative rounding guard. Same-edit tied traces remain ambiguous.
Uncertified endpoints use the complete bounded profile DP. Its suffix edit
bound may exclude states that cannot reach any of those remaining endpoints.
Original logical work is charged before certification, including on a fast
success; budget exhaustion therefore preserves the original partial result.
The fallback keeps only two rows of scores and stores predecessor bytes for
older rows, without changing the original matrix-size acceptance limit.

Mapping fragment scoring is opt-in through `--deep-validation`; fast local
mode reports scores as unmeasured. Deep scoring uses `--scoring-threads` workers (default four,
capped by available CPUs). Each comparison reuses a fixed set of workers,
with at most four eligible fragments per worker in flight and a bounded buffer
for out-of-order results. Name deduplication and input-limit checks retain input order; workers
compute independent fragment evidence. Candidate sums, support counts and gzip
diagnostics are then written in the same input order as serial execution.
No parallel floating-point reduction or read subsampling is used. When running
many samples concurrently, `--scoring-threads 1` keeps each scorer serial.

Joint evidence partitions second-mate intervals into four disjoint cases for
each first-mate placement: left extension, right extension, containment and
being contained. Indexes by start, end and length provide range maxima within
each fixed-span case. Separate observation masks preserve direct junctions,
mate brackets and contraction linkage, including circular offsets. A range
maximum applies the original floating-point pair-score expression before
deciding ties, preserving ambiguous-base contacts even when a lower individual
score rounds or clamps to the same joint score. Low-weight span mass is kept.
Small placement sets retain direct enumeration. This reduces repeated pair
visits when many alignments share endpoints; it does not claim a universal
linear bound when endpoint and evidence classes are all distinct.

Junction reports separate continuously observed 20-base anchors on each side
(`observed_junction_fragments`) from geometry-compatible mate brackets
(`mate_bracket_fragments`). The latter do not sequence the intervening gap.
Contraction linkage requires matching traceback observations of both unique
flanks in the same joint placement. Insertions or deletions through a boundary
cannot establish a continuously observed junction merely through Hit start/end.

Optional layouts compare at most eight sample hypotheses plus their baseline.
They retain the coverage incumbent, then cover repeat-block traversal-count
families, distinct block orders, and further sequence variants in that order.
Blocks come from exact shared sequence in the generated path union, grouped
within each k; they are graph traversal classes, not biological copy-number
estimates. Different k graphs are not assumed to have identical block identities.
Reports count unscored count/order families and identify coverage ranks.
Eligibility of the common unique-flank cohort is frozen using every generated
path before shortlisting, including the unchanged baseline. De novo uses an
open linear baseline. Scores are diagnostic; the coverage representative remains an inferred
layout. Even a sole candidate receives paired comparison and full-sequence read
validation in `coverage_layout_validation/coverage.tsv`. `read_validation.passed`
is not structural uniqueness or phase verification. Cross-locus phase for
variable multi-anchor circles remains a separate, incomplete method-development
task; this comparison never upgrades them to final.

Whole-layout diagnostics are streamed as gzip JSON Lines:
`coverage_layout_fragments.jsonl.gz` for de novo and
`local_assembly/locus_<id>.coverage.fragments.jsonl.gz` for mapping coverage
comparisons. Each named fragment records the common pivot, per-hypothesis
projection and accepted placement counts, search status, best sequence weight
per span, empirical span weight, their product, mass before the outlier floor,
floor fraction and ambiguous-base contacts. The floor fraction is a numerical
score diagnostic, not a calibrated biological posterior. It is not valid to
add independently maximized sequence and distance scores.

`anchor_hits_removed` and `seeded_decoy_searches` remain compatibility counters.
The profile workflow starts from empty placement sets, so these are normally
zero; they do not estimate how many decoys exist in the sequence.

Exhausted alignment/rescue or inconsistent overlapping anchor projections omit
a fragment uniformly from every hypothesis and record the unknown reason.
Without usable geometry, both mate strands are searched across the whole
hypothesis in the same edit domain; a resource limit remains unknown and no
insert-length prior is invented. Only a completed constrained search with no
anchor or mate alignment may retain the outlier score. Completeness refers to
the declared pivot, edit, strand and target-window domain; it does not cover
unrestricted alignments or ungenerated structures. Unknown query or target
bases never supply matching depth or contiguous junction evidence, even for
an N==N exact string match. Their emission model is unchanged.

`best_candidates` lists numerical ties; de novo also lists candidates within
five log-score units, without a confidence interpretation. The unchanged
coverage representative's validation records its own candidate identity and
sequence length. That validation cannot be transferred to another scored path.

`alignment_search_exhaustive=false` continues to describe the unrestricted
search. Completed profile searches are exhaustive only inside `alignment_domain`.
Neither zero unknowns nor a score margin proves structural completeness or phase.

A fragment's score is `ln(min(mass,1)+exp(-12))`. For `u` uniformly omitted
fragments, the maximum change in a pairwise score contrast is
`u * (12 + ln(1 + exp(-12)))`. Mapping retains local choices within five log
units plus this bound of the best measured score. With no measured fragments,
all shortlisted choices remain hypotheses. This bound cannot account for
ungenerated structures, model misspecification or missing biological information.

Read scoring skips exhausted or length-limited loci, and cycle families lacking
complete enumeration, consistent order, complete block representation or any
variable blocks. Those candidates remain exported. `read_scoring` is
`skipped_search_incomplete`; unmeasured scores, fragment/support counts and
score-derived rankings are null, including fragment scores and spanning counts
inside `repeat_copy_model.simple_cycles[].candidates`. Measured R1 occurrence
depths remain available. A locus without candidates uses
`skipped_no_candidates`. Eligible comparisons use `computed` and retain the
same read pool and scoring. The cycle read scan is omitted when its index is
empty; repeat-depth diagnostics continue to scan all input pairs.

End searches propagate all target labels together on one reverse worklist. A
separate bit for each target preserves its reachable domain and the existing
per-target state/path budget. This does not replace target-specific reachability
with a union or change search completion, proposal admission or reported states.

### Partial-primary allele support and reciprocal choices

Partial primary paths additionally record `primary_spanning_interval`,
`primary_spanning_fragments` and `primary_strong_votes_against_alternatives`.
When every noncircular route joins the same source and target, the interval
covers all sequence differences within shared flanks. Overlapping common
prefix/suffix sequence expands the interval across every equivalent indel
position, including the entire ambiguous homopolymer. Multiple differences
must be covered by one fragment. Other exits, anchor differences and incomplete
searches do not gain a join from this rule. The original spanning interval,
strict scores, strict selected paths and complete-publication checks are unchanged.
Skipped scoring reports the additional evidence counts as null.

For partial primary output, a comparison of alternatives takes precedence over
an opposite view containing only one route. Among competitive comparisons,
strict selection still takes precedence over an observed majority. Equal
conflicts remain open. Conflicting reverse-only projections also remain open
instead of choosing the last visited route. All strict selections remain
unchanged, and an uncertain primary join is reported as structural uncertainty.
