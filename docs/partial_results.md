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

Uncertainty is not always localized to NCR1. In particular,
`collapsed_bubble_component_unresolved` can cover the entire de novo primary
sequence. Genes overlapping that conservative flag remain `review_required`;
it does not establish that every overlapping gene contains an error. Keep the
flag until more specific evidence supports a narrower interval.

The additive `records[].gene_quality` JSON object contains the same rows, expected
model count, names in each class and a summary of `no_detected_gene_issues` or
`gene_review_required`. `ncr1_required` is false. Missing or uncertain NCR1 alone
does not fail these checks. Uncertainty extending into a gene still does.

These are annotation and supplied-interval checks, not independent read
validation. `read_validation` remains `not_performed`; without
`--assembly-ambiguity`, `assembly_ambiguity_supplied` is false, not evidence that
the assembly has no uncertain regions. RNA secondary structure is not assessed.
Inspect original read support before treating a region as reliable for analysis.

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
