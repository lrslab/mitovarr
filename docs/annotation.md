# Native gene prediction and tRNA-Pro rotation

NCR1 completeness is not required for gene-level evaluation. Every annotation
now writes `gene_quality.tsv` and `records[].gene_quality`, distinguishing
annotation-clear models from missing, incomplete, duplicated or uncertain gene
placements. Supplied assembly uncertainty still flags affected genes. See
[partial results and the limits of these checks](partial_results.md).

[中文：tRNA/CDS 检出与完整度规则](annotation.zh-CN.md).

`mitovarr annotate` is an offline, in-process Rust **homology-assisted gene predictor** for already generated mitochondrial candidates. It annotates CDS, rRNA and tRNA on both strands, retains distinct copies, and optionally standardizes circular coordinates using tRNA-Pro. The command does not select assembly paths, repair sequence bases, validate circular closure or change an assembly's publication status.

```bash
mitovarr annotate \
  --fasta runs/AF16_denovo/primary.fasta --species cbr \
  --assembly-ambiguity runs/AF16_denovo/primary_ambiguity.tsv \
  --rotate-to-trna-pro --out runs/AF16_annotation
```

All FASTA records are processed; plain and gzip inputs are accepted. Each contig has its own completeness and rotation result. Output directories must be new or empty. Exit 0 means annotation finished; missing genes or a skipped rotation are biological results in the report, not a failed assembly or an execution error. Exit 1 indicates an input/runtime error.

## Models and scope

The executable contains fixed panels from these NCBI records. Each has 12 CDS, 2 rRNA and 22 tRNA models, extracted in transcriptional orientation. No network or external annotation program is used at runtime.

| Option | Organism | Model source |
| --- | --- | --- |
| `--species cel` | C. elegans | [NC_001328.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001328.1) |
| `--species cbr` | C. briggsae | [NC_009885.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_009885.1) |
| `--species ctro` | C. tropicalis | [NC_025756.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_025756.1) |

The model sequences, source coordinates, products and reference protein translations are in `src/annotation/models/`. The complete accession sequences in `tests/fixtures/annotation/` are positive controls. C. briggsae's source genes are on the minus strand. Source sequence topology is not used to assert topology of an assembly candidate. The C. briggsae source LOCUS is marked linear; its test fixture explicitly assumes circular coordinates to exercise reverse-orientation normalization, without using that test as closure evidence. The models were retrieved/assembled into this package on 2026-09-18 (Hong Kong).

Leu/Ser isoacceptors retain the source codon-family labels, such as `tRNA-Leu-CUN` and `tRNA-Ser-UCN`, instead of potentially inconsistent numbered labels. The small and large rRNAs are `rrnS` and `rrnL`.

This implementation is intended for these species and related strain variation. It is not a reference-free discovery model for novel genes, a covariance-model RNA predictor, or an assay of RNA processing/function. RNA structural assessment is explicitly `not_assessed`; a complete homolog is not proof of an intact functional RNA. Divergent or very short fragments can remain unannotated.

## Prediction and completeness

1. Exact 9-mer seeds (7-mer for tRNA) nominate local windows on both strands. A word occurring more than 64 times in the searched sequence is not used to nominate windows. Seeds on a nominated diagonal must span at least 24 nt or 5% of the model length, whichever is larger (capped at the model length), so overlapping seeds in one tiny shared motif do not nominate an entire gene. The union of seeded model bases must cover at least three seed lengths or half of that required span, whichever is larger (capped at the model length). An affine-gap Smith–Waterman alignment then allows mismatches and indels. The scores are +3 match, −4 mismatch, −8 gap open, −1 extension. They are alignment scores, not calibrated probabilities.
2. Distinct physical placements are kept, including inversions and partial homologs such as ND5-derived copies. Disjoint sequence remaining on either side of a window hit is searched again, so a complete gene cannot mask an adjacent partial copy. Overlapping seed-window observations of the same gene are merged; one hit per gene and reference gene order are **not** enforced. Substantial overlap between different gene models is reported as ambiguous gene identity.
3. Both model endpoints and terminal alignment evidence are checked independently. A partial alignment is reported as partial, even when its coverage exceeds 85%. For a CDS terminus clipped by at most 9 nt, the predictor can extend onto existing sample bases when the projected boundary has an allowed start/stop and at least 75% identity in the terminal 12 nt. This boundary is labelled `orf_inferred`. RNA termini are not extended this way.
4. CDS checks use [NCBI translation table 5](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/genetic-codes/#5-the-invertebrate-mitochondrial-code-transl_table5): TTG/ATT/ATC/ATA/ATG/GTG starts, TAA/TAG stops, TGA→Trp and AGA/AGG→Ser. The initial allowed start is translated as Met. Internal stops, ambiguous codons and non-triplet alignment gaps are reported. An in-frame deletion need not invalidate a CDS. A terminal T or TA in an otherwise defined coding boundary is labelled as requiring polyadenylation; no A bases are added and RNA maturation is not demonstrated. If the 3′ boundary is unsupported, `cds.stop` is `partial_or_uncertain`: a trailing T/TA is not evidence of polyadenylation, and a trailing full stop remains in the translated evidence. A terminal partial codon alone does not invalidate `frame_valid` for a clipped gene; non-triplet alignment gaps still do. This field does not establish completeness.

`full_length_supported` requires supported endpoints and an unambiguous sample sequence. CDS additionally needs an allowed start, a terminal stop interpretation, no internal stop, and no non-triplet alignment gap. RNA additionally requires an ungapped full homolog; tRNA identity must reach 90%. RNA indels receive `boundary_complete_length_variant`, retaining their boundaries and sequence without claiming structural completeness. Unknown bases, disrupted CDS, partial genes, and competing gene identities remain explicit.

Detection thresholds deliberately differ from completeness: CDS/rRNA fragments require ≥20% model coverage and ≥80 aligned model bases (or the model length if shorter), at ≥70%/75% identity respectively; tRNA requires ≥80% coverage, ≥30 bases and ≥85% identity. An endpoint counts as supported when its terminal 8 alignment columns have no gaps and at most one mismatch. Coverage is the model interval spanned; internal gaps are reported separately. These are heuristics, not biological certainty thresholds.

The seeded search examines at most 64 distinct start neighborhoods per gene/strand. If more pass, `search_limited_models` identifies the affected models; their absence/uniqueness cannot be concluded from the search. A limited tRNA-Pro search cannot authorize rotation. Annotation accepts mitochondrial candidates up to 200 kb per record and 5 Mb total, and model sequences up to 5 kb. Limits are reported as errors or explicit incomplete searches, never silently treated as exhaustive discovery.

### Detection is distinct from completeness

The search compares nucleotide models with candidate DNA; CDS translation is a subsequent check, not a protein-homology search. tRNA identity is assigned through sequence homology, without an independent anticodon or RNA-structure check. Annotation does not reload reads to establish support for a gene or variant.

| Detection criterion | tRNA | CDS | rRNA |
| --- | ---: | ---: | ---: |
| Model interval coverage | ≥80% | ≥20% | ≥20% |
| Alignment identity | ≥85% | ≥70% | ≥75% |
| Model interval length | ≥30 nt | ≥80 nt | ≥80 nt |

For shorter custom models, the interval-length threshold is capped at model length. Identity counts matching bases divided by all alignment columns, including gaps. Model coverage measures the spanned model interval, including positions deleted from the sample.

`missing_models` means that no retained feature passed the detection threshold. It means **not detected**, rather than demonstrated genomic absence. `models_without_full_length_support` includes both undetected models and detected models without a complete qualifying copy. A `partial`, `disrupted` or `boundary_complete_length_variant` feature is still detected. The 20% CDS detection threshold retains fragments such as partial ND5 homologs; it does not qualify them as complete CDS.

`feature_counts` and `full_length_counts` count physical features and can include multiple copies of one model. To assess model completeness, count distinct model names with at least one qualifying feature. Subtracting the complete count from 22 does not give a tRNA deletion count.

### Known RNA classification limitations

The strict completeness rule accepts internal SNPs within the identity and endpoint thresholds. An RNA with supported boundaries but an internal indel receives `boundary_complete_length_variant` and remains detected, while being excluded from the strict complete count.

A terminal SNP can be clipped by local alignment and cause a `partial` classification even when the corresponding sample base is present. A shorter annotated interval therefore does not establish that the assembly is missing sequence. RNA terminal-mismatch handling remains an implementation follow-up; RNA structure and function remain unassessed. The [Chinese method guide](annotation.zh-CN.md) explains the field meanings.

## Topology and rotation

Default `--topology auto` first reads `annotation_topology=circular|linear` from the FASTA header. Without that annotation tag, only a literal `circular=true` asserts a circular input hypothesis; other records are linear. Explicit `--topology linear` or `--topology circular` overrides either header convention. This option asserts an input convention; it does not establish the missing connection of a linear contig.

De novo `primary.fasta` records `circular` and `final_verified` in its header and uses the record ID expected by `primary_ambiguity.tsv`. The independent `best_circular.fasta` output carries the same metadata and pairs with `best_circular_ambiguity.tsv`. De novo `consensus.fasta` instead has the plain header `>mitochondrion`; mapping `candidate.fasta` and `consensus.fasta` likewise contain only the contig name. These unmarked records default to linear even when a de novo final assembly passed circular validation. Supply `--topology circular` explicitly to rotate an unmarked output when supported by its assembly evidence. An ambiguity table must match the selected FASTA's record IDs and coordinates.

Every `annotated.fasta` records the topology used in `annotation_topology`, so annotating that output again with auto preserves its interpretation. Existing assembly tokens such as `circular`, `status` and `final_verified` are retained unchanged. An explicit circular annotation assumption does not turn `circular=false` or an unresolved assembly into a validated circle.

For a mapping consensus in a known circular coordinate system, an explicit example is:

```bash
mitovarr annotate \
  --fasta runs/mapping/candidate.fasta --species cel \
  --topology circular --rotate-to-trna-pro --out runs/mapping_annotation
```

A circular hypothesis can be normalized while its assembly remains unresolved. Annotation never upgrades it to a validated circle. In particular, do not use `--topology circular` to join arbitrary de novo fragments.

Rotation requires exactly one detected tRNA-Pro locus, full-length support, and an aligned model 5′ boundary. A partial/ambiguous tRNA-Pro, multiple copies, an incomplete search, or a linear input leaves the sequence unrotated and gives a specific reason in the report. No high-scoring internal tRNA segment is used as the gene's 5′ end.

If tRNA-Pro is on the minus strand, the entire sequence is reverse-complemented. It is then rotated left to place that tRNA's biological 5′ base at position 1 on the plus strand. Every feature's coordinates and strand are transformed together. Input bases are preserved (including IUPAC ambiguity symbols), apart from case normalization and this reverse-complement/permutation. No reference bases replace sample bases.

The report records `input_anchor_1based`, `reverse_complemented` and `offset_0based`. For an input base at 0-based coordinate `x` and sequence length `L`, first set `y = L−1−x` if reversed, otherwise `y = x`; its output coordinate is `(y + L − offset_0based) % L`. This maps candidate coordinates. Existing BAM pileup/variant tables use the original reference coordinate system and are **not** rewritten by this command; accepted consensus indels require their separate coordinate mapping before applying this formula.

## Carrying assembly uncertainty into annotation

Optionally supply `--assembly-ambiguity primary_ambiguity.tsv` with the matching de novo `primary.fasta`, `best_circular_ambiguity.tsv` with `best_circular.fasta`, or local mapping `structural_ambiguity.tsv` with its `candidate.fasta`. A previously projected annotation table is also accepted with its matching FASTA. Rust transforms the existing intervals with the sequence, splits at the new origin and merges adjacent/overlapping pieces from the same evidence. Each projected region keeps its source evidence, reason, alternative count and input segments, and lists overlapping genes and NCR loci. This operation uses existing evidence; it does not rescore reads or infer structural variants.

The optional outputs `assembly_ambiguity.tsv` and `assembly_ambiguity.gff3` use **one-based inclusive coordinates in `annotated.fasta`**. The TSV retains the six input columns (`sequence, start, end, evidence, reason, alternatives`) and appends `id, input_segments_1based, genes, ncr_loci`; split regions have one row per output segment. The projected TSV and `annotated.fasta` can be supplied together to a subsequent annotation. Inputs must already use the same coordinates: mapping `ambiguity.tsv` describes the reference and is not accepted by this option.

Each JSON record adds `assembly_evidence` when the option is supplied, otherwise null. Its `ncr_loci` keep the biological `annotation_status` separately from `evidence_status`: `open_contig_end`, `unknown_sequence`, `unresolved_path_alternatives`, `retained_assembly_alternatives` or `no_reported_path_ambiguity`. All source regions remain available even when one status takes precedence. Absence of a reported path ambiguity does not establish resolved structure; `structural_resolution_from_reads` is explicitly `not_certified_by_coordinate_projection`.

Evidence is associated with positional loci even when no NCR interval passes the length threshold. For a zero-length gap, it overlaps either base beside the touching-gene junction; for overlapping flanks, their shared bases; for a short gap, the actual intergenic bases. Thus `adjacent_without_gap` can correctly coexist with `retained_assembly_alternatives`. Ambiguity inside a flanking gene away from that junction is not automatically NCR ambiguity. Mapping local assembly's `structural_ambiguity.tsv` can also be supplied with its matching `candidate.fasta`; the output column formats are unchanged.

## Files

| Output | Contents |
| --- | --- |
| `input.fasta` | All input records, with original headers; DNA normalized to uppercase |
| `annotated.fasta` | The same records, with successful rotations applied and `annotation_topology` recorded; existing assembly status tokens retained |
| `annotations.gff3` | Gene, NCR and optional projected assembly regions in `annotated.fasta` coordinates; circular-origin splits use the same feature ID and strand-aware CDS phase |
| `features.tsv` | One row per physical feature: segments, model interval, identity, boundaries, completeness and coding defects |
| `gene_quality.tsv` | One row per expected model: annotation status, placement counts, issues and overlap with supplied assembly uncertainty; NCR1 completeness is not a gate |
| `genes.fna` | All detected feature sequences in biological 5′→3′ orientation, including flagged partial/alternate copies |
| `proteins.faa` | Translations with partial/defective status in the header; internal stops and X codons retained |
| `ncr.gff3`, `ncr.tsv`, `ncr.fasta` | Positional NCR intervals and generic gaps, their boundaries and sample sequences |
| `ncr_loci.tsv` | Four NCR position states per species contig, including unresolved or short intervals |
| `assembly_ambiguity.tsv`, `assembly_ambiguity.gff3` | Only with `--assembly-ambiguity`: projected source intervals, including open ends and retained path alternatives |
| `report.json` | Per-contig genes, missing and non-full-length models, `gene_quality`, search limits, coordinate transformation, NCR regions/loci and optional `assembly_evidence` |

GFF/TSV/JSON segment positions are 1-based inclusive. Gene segments follow biological order; NCR segments follow forward genomic order through a possible circular origin. For partial CDS the phase is inferred from the model interval, and completeness remains partial; translation is never silently skipped because the GFF phase is missing. Gene/protein FASTA IDs match the GFF feature IDs. A `feature_counts` total includes partial and duplicate gene placements; `full_length_counts` includes supported genes, and `models_without_full_length_support` asks whether each expected model has at least one supported placement. Counts alone do not validate an assembly.

## Custom panels

Use `--gene-models panel.json` instead of `--species`. A panel contains a schema, accession, genetic code 5, and unique named models in transcriptional orientation:

The `mitovar-rust.gene-models.v1` schema identifier is deliberately retained
after the package rename to `mitovarr`; existing panels remain valid.

```json
{
  "schema": "mitovar-rust.gene-models.v1",
  "accession": "my-curated-panel-v1",
  "genetic_code": 5,
  "genes": [
    {"name": "tRNA-Pro", "kind": "tRNA", "sequence": "ACGTACGTACGTACGT", "product": "tRNA-Pro"}
  ]
}
```

The short illustrative sequence above is a format example, not a biological tRNA model. Supply complete curated CDS/rRNA/tRNA nucleotide sequences (15–5000 A/C/G/T bases); provenance fields may be kept alongside them. Runtime GenBank parsing and external database searches are not required.

## NCR and SV interpretation

`annotate` derives positional NCR1–4 labels from observed flanking genes, with a configurable `--ncr-min-len` (default 20 bp). The [NCR method](ncr.md) documents species-specific flanks, candidate labels, ignored short secondary homologs, unknown sequence and separate linear end fragments. A gap next to an uncertain/partial endpoint is not promoted to a supported NCR insertion. NCR sequence ancestry and structural variants from read evidence remain unassessed. Neither complete annotation nor tRNA-Pro rotation resolves repeat copy number or circular closure.

## Verification

Regression tests use bundled reference fixtures and small constructed examples. They cover gene detection and completeness, both strands, origin-spanning intervals, repeated annotation, explicit topology overrides, positional NCR labels and assembly-uncertainty projection.
