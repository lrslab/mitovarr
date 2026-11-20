# Partial results and gene quality

NCR1 completeness is not required to retain and independently evaluate a partial
mitochondrial sequence. Assembly and NCR layout algorithms are held at the
current baseline; coverage layouts remain optional hypotheses. The primary
workflow now emphasizes gene sequence, uncertainty and variant evidence.

## Select the existing sequence

Use `primary.fasta` with `primary_ambiguity.tsv` for de novo work. For default local
mapping, use `consensus.fasta` with `maximum_consensus.ambiguity.tsv`. The
separate pileup/strict `candidate.fasta` uses `structural_ambiguity.tsv`. Keep any optional
coverage layout separate. Annotation evaluates frozen sequences independently
of assembly selection; keep any later choice based on gene content explicit.

New assembly reports include `sequence_output.status`: `complete`,
`partial_available`, `best_effort` (mapping best-effort mode), `unavailable` (de novo),
or `no_called_bases` (mapping).
`sequence_output.fasta` names the retained representative. `partial_available`
only means a sequence is present. Its gene quality has not yet been assessed.
Old reports lack this field; their retained sequences and uncertainty reports
remain valid inputs to annotation.

Assembly exit code 2 still means unresolved under the existing completion rules.
It can also mean unusable recruitment/extraction or a CLI usage error. Inspect
the report and actual outputs before continuing. Do not suppress every nonzero
exit or substitute a reference for missing sequence.

Fast local mapping enables best-effort output by default and defers structural
scoring. Add `--deep-validation` for full scoring; in that mode `--best-effort`
explicitly prioritizes a maximum usable consensus.
It writes `consensus.fasta` and `maximum_consensus.fasta` with
`maximum_consensus.ambiguity.tsv`, and exits 0 when written even if validation
remains incomplete. Check `status=best_effort`, `consensus_written`, and
`strict_validation_passed`. The retained `final_published` field still refers
to strict validation. Use the maximum consensus's own ambiguity TSV for
annotation; `structural_ambiguity.tsv` belongs to strict `candidate.fasta`.

## Evaluate genes independently

When a circular candidate is specifically needed, de novo runs provide
`best_circular.fasta` with `best_circular_ambiguity.tsv` if an eligible ring was
found. Check `report.json.best_circular.status` and `final_verified`; a selected
ring may remain unresolved. Its main-component identity must still be checked,
especially when `primary.fasta` is linear. Pass its own matching ambiguity table:

```bash
mitovarr annotate --fasta runs/sample_denovo/best_circular.fasta --species cel \
  --assembly-ambiguity runs/sample_denovo/best_circular_ambiguity.tsv \
  --out runs/sample_circle_genes
```

```bash
mitovarr annotate --fasta runs/sample_denovo/primary.fasta --species cel \
  --assembly-ambiguity runs/sample_denovo/primary_ambiguity.tsv \
  --assembly-report runs/sample_denovo/report.json \
  --out runs/sample_denovo_genes

mitovarr annotate --fasta runs/sample_mapping/candidate.fasta --species cel \
  --assembly-ambiguity runs/sample_mapping/structural_ambiguity.tsv \
  --out runs/sample_mapping_genes
```

Topology defaults to the existing FASTA annotation/header hints; an unmarked
mapping candidate is linear. Only use `--topology circular` to evaluate an
explicitly stated circular hypothesis. It does not certify a join. Rotation is
optional and preserves uncertainty coordinates.

`gene_quality.tsv` gives one row per expected model per input record. Its columns
are `contig`, `gene`, `type`, `status`, `placements`, `full_length_placements`,
`feature_ids`, `assembly_region_ids`, `issues`, and
`assembly_ambiguity_supplied`. Feature and region IDs refer to the existing
annotation outputs in `annotated.fasta` coordinates.

`annotation_clear` means one full-length supported placement and no detected
gene-level issue. `review_required` records missing or incomplete models,
multiple homologous placements, limited annotation searches, ambiguous sequence,
competing gene identities, CDS boundary/frame problems, or overlap with supplied
assembly uncertainty. Extra partial homologs trigger review; this is not a claim
that a functional gene is duplicated. Observed disruptions also need read review
before being called biological variants.

Uncertainty is not always localized to NCR1. De novo `primary_ambiguity.tsv`
uses these local reasons:

| reason | what the interval marks |
| --- | --- |
| `collapsed_bubble_alternative` | the bases where the primary differs from the other side of a permissively collapsed bubble; evidence `collapsed_bubble_<n>` names its record in `bubble_alternatives.fasta` and `graph_simplification.retained_bubble_alternatives` |
| `homopolymer_length_family` | a homopolymer tract and both flanks whose selected length is a consensus with supported alternatives |
| `minor_path_retained` | a strictly selected local route whose minor route is kept on record |
| `majority_path_site_unresolved` | a majority-chosen route whose alternatives differ at one small site: up to three substituted bases, or an indel of up to six bases whose shorter side has at most two |
| `majority_path_structure_unresolved` | a majority-chosen route with a longer or divergent alternative |
| `topology_inferred_repeat_pairing` | a repeat copy whose pairing was inferred from topology rather than spanning reads |
| `low_matching_read_depth` | bases of the primary below `--min-depth` matching fragments in its read validation |
| `open_graph_end` | the first or last base of a linear primary |

Genes overlapping any interval remain `review_required`; an interval does not
establish that the gene contains an error. Outputs of v0.2.2 and earlier mark
every permissive bubble of the component with one
`collapsed_bubble_component_unresolved` interval over the whole primary; that
flag cannot be localized after the fact. A bubble whose junction flanks are not
found on the primary marks nothing and is listed in
`primary.collapsed_bubbles.not_placed`. Final publication still counts every
bubble of the component.

The additive `records[].gene_quality` JSON object contains the same rows, expected
model count, names in each class and a summary of `no_detected_gene_issues` or
`gene_review_required`. `ncr1_required` is false. Missing or uncertain NCR1 alone
does not fail these checks. Uncertainty extending into a gene still does.

These are annotation and supplied-interval checks. Annotation reads no
sequencing reads, so `gene_quality.read_validation` remains `not_performed`.
The de novo report's `primary.read_validation` gives the per-base check of the
primary itself; its low-depth bases reach annotation only as
`low_matching_read_depth` intervals of a supplied `primary_ambiguity.tsv`.
Without `--assembly-ambiguity`, `assembly_ambiguity_supplied` is false, not
evidence that the assembly has no uncertain regions. RNA secondary structure is
not assessed. Inspect original read support before treating a region as
reliable for analysis.

## Claim level: complete except NCR1

The AT-rich core of NCR1 is often not covered by short reads, and its length can
vary within one sample. A de novo primary can still be complete everywhere
else. `annotate --assembly-report <de novo report.json>` states this as a
separate claim level, `complete_except_ncr1`, for `primary.fasta` with its
matching `primary_ambiguity.tsv`. The level is a report; it changes no
sequence, no assembly status and no strict criterion.

The level is reached when all of these hold on the frozen primary:

1. **One contig.** A circle, or a linear primary with NCR1 at both ends: of
   the genes counted as present (full length, or an RNA length variant accepted
   under rule 2), the one starting first and the one ending last are tRNA-Ala
   and tRNA-Pro.
2. **Genes.** Every expected gene is placed once and full length. An RNA whose
   two boundaries align but whose length differs from its model counts as
   present and is listed in `rna_length_variants`.
3. **One sample.** No mixed-sample signal and no structural mixture.
4. **Nothing open outside NCR1.** Outside NCR1 the primary carries no
   interval other than alternatives kept on record: collapsed bubbles
   (`collapsed_bubble_alternative`), homopolymer families
   (`homopolymer_length_family`), single-site majority choices
   (`majority_path_site_unresolved`) and minor routes the strict rules rejected
   (`minor_path_retained`, `minor_cycle_path_retained`). Any other reason,
   including `majority_path_structure_unresolved`, `competing_cycle_paths`,
   `topology_inferred_repeat_pairing`, `low_matching_read_depth` and reasons
   this version does not know, keeps the level out.
5. **Read validated.** `primary.read_validation.status` is `performed`.

Alternatives kept on record outside NCR1 are counted in
`flagged_sites_outside_ncr1` (by reason in `flagged_site_reasons_outside_ncr1`);
they do not remove the level. Each points to its evidence in the assembly
report, for example the support of a collapsed bubble's two sides. Intervals
that lie wholly inside NCR1 never count against it.

The closure state is reported with the level and never decides it:

| `closure` | meaning | statement |
| --- | --- | --- |
| `circle` | the primary is a circular path or cycle; only `strict_final_published=true` makes it a verified circle | complete outside NCR1; circular graph path or cycle candidate, not a verified circle; NCR1 length not determined. With the strict pass: verified circle; also complete outside NCR1 |
| `joined_by_reads` | the de novo path graph leads from the last segment of the linear primary back to its first within 2,000 added bases, crossing neither strand of the primary's own segments; overlapping sample k-mers only, the join is not applied | complete outside NCR1; overlapping reads join the two NCR1 ends; NCR1 length not determined; not a verified circle |
| `not_observed` | no such path | complete outside NCR1: one contig with every expected gene and both ends inside NCR1; reads do not join the ends, so the sequence is not shown to be circular |

Do not describe a `not_observed` primary as circular or as a complete genome.
`ncr1_assembled_bases` lists the assembled NCR1 bases on the primary (two
flanks of a linear contig, one region of a circle); `ncr1_length` is always
`not_determined`.

`status` is `reached`, `not_reached` (with the reasons in `not_reached`) or
`not_assessed` (with `not_assessed_reason`). The level is not assessed without
`--assembly-report`, for a record other than the reported `primary`, when the
annotation topology differs from the assembly, without `--assembly-ambiguity`,
when that table's rows differ from the report's `primary.uncertainty_intervals`
(it must be the same run's `primary_ambiguity.tsv`),
without an NCR catalog (`--gene-models` instead of `--species`), or when the
report predates these fields or its read validation hit `--max-validation-pairs`.

`claim_level.tsv` has one row per input record: `contig`, `level`, `status`,
`closure`, `strict_final_published`, `not_reached`, `not_assessed_reason`,
`rna_length_variants`, `flagged_sites_outside_ncr1`,
`flagged_sites_inside_ncr1`, `collapsed_bubbles_not_placed`,
`ncr1_assembled_bases`, `ncr1_length` and `statement`. `records[].claim_level`
holds the same fields plus `flagged_site_reasons_outside_ncr1`, the `rule` and
its `scope`.

## Preserve the distinction between usable regions and complete genomes

Strict completion requirements remain separate from best-effort consensus
availability. Annotation does not modify or reselect sequence and does not publish
a complete assembly. A complete set of gene annotations is not evidence for
NCR1 copy number, order or circularity.

Report partial sequence length separately from total mitochondrial genome
length. Unknown gap length stays unknown; an N run is not a measured gap length.
Candidate length ranges are not confidence intervals. Do not infer deletions
from an absent annotation model alone.

For gene variation and phylogeny, review genes and mask evidence-defined uncertain
positions, including affected coding regions. Mapping `variants.tsv`,
`coverage.tsv`, `ambiguity.tsv`, `uncertainty.tsv` and `base_votes.tsv` remain in
original reference coordinates, even after indels or local replacements.
`structural_ambiguity.tsv` uses candidate coordinates. Annotation does not
recompute coverage or move variant calls into rotated coordinates.
