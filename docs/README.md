# mitovarr documentation

`mitovarr` is a native Rust package for mitochondrial **de novo assembly** and
**mapping assembly**, read extraction, gene annotation and non-coding region
(NCR) analysis. The Cargo package, library and executable are all named
**`mitovarr`**. These guides describe this package's current commands and outputs.

Start with [build and installation](../README.md#build-and-install), then choose
a workflow below. The [CLI reference](cli.md) lists the available options and
defaults. Each `--out` directory must be new or empty.

## Choose a workflow

There are two assembly workflows:

- **De novo assembly** builds sequence from paired sample reads. Its two
  input routes are BAM-extracted FASTQ and seed-based recruitment from original
  FASTQ.
- **Mapping assembly** uses a matching reference and existing BAM alignments.
  Its command is `assemble-mapping`, with optional reconstruction from matching
  paired FASTQ.

`fetch-caendr` and `bait` prepare reads; `annotate` evaluates the resulting
sequences. Choose the command that matches your inputs:

| Your starting point or task | Command | What it does |
| --- | --- | --- |
| A CaenDR species and strain | `fetch-caendr` | Extract mitochondrial BAM records and recover available mates into local BAM and FASTQ files. |
| Paired FASTQ extracted from mitochondrial BAM records | `assemble-baited` | Assemble those sample reads directly and validate with those pairs. No second recruitment step is needed. |
| Original paired FASTQ plus a mitochondrial seed FASTA | `assemble-seeded` | Recruit complete pairs using the seed, assemble sample reads, and validate with all original pairs. |
| Recruitment only | `bait` | Recruit complete pairs from original FASTQ using an explicit seed; report convergence. |
| Reference FASTA plus existing coordinate-sorted, indexed BAM | `assemble-mapping` | Mapping assembly: generate a consensus and optionally reconstruct uncertain regions using matching paired FASTQ. |
| An existing candidate or consensus FASTA | `annotate` | Assess gene placements and completeness, identify positional NCRs, and optionally rotate eligible circular sequences to tRNA-Pro. |

The two de novo input routes retain separate evidence scopes. A seed only
recruits reads in seeded mode; it does not supply assembly sequence. BAM
extraction has `mt-bam-extracted` input and validation scope, since the upstream
alignments determine which reads are available.

FASTA and FASTQ may be plain text or gzip. Paired FASTQ uses separate,
synchronized R1/R2 files with matching read names. Mapping requires a BAM index
(`.bam.bai` or `.bai`); this package does not accept CRAM or globally align FASTQ
to a reference. See the [input contracts](../README.md#workflows-and-input-scope).

## De novo assembly

For original whole-sample paired FASTQ, use a mitochondrial seed to recruit pairs:

```bash
mitovarr assemble-seeded \
  --r1 original_R1.fastq.gz --r2 original_R2.fastq.gz \
  --seed mitochondrial_seed.fasta --out runs/sample_seeded
```

For CaenDR data, first extract the mitochondrial reads, then use
`assemble-baited`. For example, extract *C. elegans* N2 and assemble all recovered
complete pairs:

```bash
mitovarr fetch-caendr --species cel --strain N2 --out runs/N2_reads

mitovarr assemble-baited \
  --r1 runs/N2_reads/reads.R1.fastq.gz \
  --r2 runs/N2_reads/reads.R2.fastq.gz \
  --out runs/N2_denovo
```

Inspect `runs/N2_denovo/report.json` and the retained FASTA files before
annotation. `primary.fasta` is the main representative; `best_circular.fasta`
is also written when an eligible circular candidate exists. Either can remain
unresolved. See [de novo outputs](de_novo.md#outputs-and-status) and
[CaenDR extraction](caendr.md).

## Mapping assembly

**`assemble-mapping` generates a sample consensus from existing BAM
alignments and a matching reference**. It collects read evidence at each
reference position, calls bases and supported insertions/deletions, and marks
weak or conflicting sequence with `N`. Overlapping mates contribute one
fragment observation. Missing sequence is not copied from the reference.

The reference must match the selected BAM contig's name, length and coordinate
system. Replace the example paths and `MtDNA` contig name with your sample's
values. The three modes below use the same command and add different work.

### 1. Consensus from BAM alignments

```bash
mitovarr assemble-mapping \
  --ref reference.fasta --bam mt.bam --contig MtDNA \
  --emit-variants --emit-coverage --out runs/sample_pileup
```

This requires no FASTQ input. It writes `candidate.fasta` and evidence reports;
`consensus.fasta` is written when the strict checks pass. Unresolved sequence
remains available in the candidate.

### 2. Fast consensus with all supplied paired reads (default local mode)

```bash
mitovarr assemble-mapping \
  --ref reference.fasta --bam mt.bam --contig MtDNA \
  --local-assembly --r1 reads.R1.fastq.gz --r2 reads.R2.fastq.gz \
  --emit-variants --emit-coverage --out runs/sample_consensus
```

Local mode builds the sample graph from all supplied complete pairs. It preserves
called bases, attempts fixed-span N fills, and outputs a best-effort consensus
with matching uncertainty. It defers structural enumeration and candidate scoring;
`validation_mode=fast_consensus` and `scoring_status=not_run_fast_consensus` make
that explicit. Unmeasured scores and fragment counts are null. Fast output never
passes strict structural validation, including when it contains no N.

Use matching pairs from the same sample, such as those produced by `fetch-caendr`.
The existing `--best-effort` flag is optional in fast local mode. Missing evidence
can leave N or an unresolved region; sequence availability is separate from
structural verification.

### 3. Deep structural validation (opt-in)

```bash
mitovarr assemble-mapping \
  --ref reference.fasta --bam mt.bam --contig MtDNA \
  --local-assembly --deep-validation --best-effort \
  --r1 reads.R1.fastq.gz --r2 reads.R2.fastq.gz \
  --emit-variants --emit-coverage --out runs/sample_deep
```

`--deep-validation` restores the previous full local search, empirical fragment
model and paired-candidate scoring. It can recover structural insertions and
deletions supported by both junctions. Add `--best-effort` to keep a maximum
usable primary when validation remains incomplete; otherwise strict output
rules apply. Length-guided legacy alternatives are available only in deep mode.
The reference-guided span remains a prior, not a measured repeat count.

`--best-effort` still requires `--local-assembly` and both FASTQ mates; it does
not enable local assembly by itself. The old `--legacy-bridge` spelling remains
a hidden alias. Mapping `--coverage-layout` additionally requires
`--deep-validation`. These switches do not apply to the de novo commands.

A successful run writes identical primary sequences to `consensus.fasta` and
`maximum_consensus.fasta`, with `maximum_consensus.ambiguity.tsv`. The strict
`candidate.fasta` remains available for comparison. If strict validation is
incomplete, the report records `status=best_effort`, `consensus_written=true`
and `strict_validation_passed=false`; successful consensus output exits 0.
A strict pass can still report `status=complete`. Best-effort output can retain
Ns or unresolved repeat structure.

See the [maximum-consensus method](mapping_local_assembly.md#maximum-consensus-and-legacy-bridges)
for gap-filling rules, alternatives and limits. Optional `--coverage-layout`
and `--multi-k-layout` produce separate layout hypotheses; they are not needed
for the best-effort command above.

## Choose the sequence and its uncertainty report

| Result | FASTA | Matching uncertainty file |
| --- | --- | --- |
| Main de novo representative | `primary.fasta` | `primary_ambiguity.tsv` |
| Best available eligible de novo circle, when present | `best_circular.fasta` | `best_circular_ambiguity.tsv` |
| Mapping candidate with local assembly | `candidate.fasta` | `structural_ambiguity.tsv` |
| Best-effort mapping primary | `consensus.fasta` or `maximum_consensus.fasta` | `maximum_consensus.ambiguity.tsv` |

De novo runs also write `consensus.fasta` when their strict final checks pass.
For annotation with uncertainty intervals, retain the representative and its
matching table listed above.

Inspect `report.json.sequence_output` for sequence availability. For mapping,
`strict_validation_passed` and `final_published` record strict validation
separately from `consensus_written`. For a selected de novo circle, inspect
`report.json.best_circular.status` and `final_verified`. A FASTA filename alone
does not establish complete-genome status.

For open-end diagnostics, `end_bridge_proposals.fasta` preserves discovered
connections and `report.json.paired_path_resolution.end_assembly` explains
their eligibility and search limits. Withheld proposals are not graph edges.
`end_bridge_candidates.fasta` contains the connectors actually admitted to the
path graph; neither file alone certifies a complete assembly.

Mapping `variants.tsv`, `coverage.tsv`, `ambiguity.tsv`, `uncertainty.tsv` and
`base_votes.tsv` describe the reference-coordinate pileup. They are not
recomputed on locally inserted sequence. In particular, pileup `ambiguity.tsv`
cannot be passed directly to `annotate --assembly-ambiguity`; use the matching
candidate/primary table above. FASTA record names and coordinates must agree.

## Annotate and assess retained sequence

NCR1 completion is not required to evaluate genes and other usable regions.
For a *C. elegans* best-effort mapping result:

```bash
mitovarr annotate \
  --fasta runs/sample_consensus/consensus.fasta --species cel \
  --assembly-ambiguity runs/sample_consensus/maximum_consensus.ambiguity.tsv \
  --out runs/sample_genes
```

Use `cbr` for *C. briggsae* and `ctro` for *C. tropicalis*. Review
`gene_quality.tsv`, `report.json`, the annotated features and NCR locus tables.
Annotation assesses the supplied sequence and uncertainty; it does not select
an assembly or independently validate reads, repeat copy number or closure.

Rotation is optional. `--rotate-to-trna-pro` requires an input interpreted as
circular and a suitable unique, complete tRNA-Pro anchor. Unmarked mapping
FASTA records are treated as linear by default; `--topology circular` explicitly
sets an annotation hypothesis. See [topology and rotation](annotation.md#topology-and-rotation)
and [partial results and gene quality](partial_results.md).

## Run status and resource limits

- Exit **0** means the command completed under its rules. This includes a
  written best-effort consensus or annotation with gene issues to review.
- Exit **2** can mean unresolved strict assembly, incomplete recruitment or
  extraction without complete pairs. CLI usage errors also return 2; inspect
  the message and report before using any retained output.
- Exit **1** indicates an input or runtime error.

Use all available reads in the selected workflow. Consider subsampling only
after a full-input run actually reaches a documented resource limit. Preserve
complete pairs and record the triggering limit, retained counts and fixed seed.
An unresolved repeat or graph-search cap alone is not a reason to discard reads.
See [mapping subsampling](mapping_local_assembly.md#optional-high-depth-subsampling)
and [de novo subsampling](de_novo.md#subsampling-after-a---max-kmers-overflow)
for the explicit sampling options and their evidence scope.

## Detailed guides

| Topic | Guide |
| --- | --- |
| Installation, input formats and defaults | [Package README](../README.md), [CLI reference](cli.md) |
| Binary platforms and maintainer release procedure | [Publishing releases](publishing.md) |
| Extraction and mate recovery | [CaenDR interface](caendr.md) |
| Graph construction, cleanup and de novo outputs | [Native de novo assembly](de_novo.md) |
| Fragment scoring, repeat paths and circular selection | [Paired-path method](paired_paths.md) |
| Pileup, local reconstruction and maximum consensus | [Mapping local assembly](mapping_local_assembly.md) |
| Retained partial sequences and gene quality | [Partial results](partial_results.md) |
| Gene models, CDS/RNA checks and rotation | [Annotation](annotation.md), [中文规则](annotation.zh-CN.md) |
| Positional NCR1–4 labels and boundaries | [NCR definitions](ncr.md) |
| Figure 1 overview | [Legend and source](figures/README.md), [PDF](figures/package_workflow.pdf), [SVG](figures/package_workflow.svg) |
| Detailed method schematic | [Workflow figure and legend](figure_1_workflow/README.md) |

Existing `mitovar-rust.*` JSON schema strings remain format identifiers for
compatibility with saved reports and gene-model panels. Current commands use
`mitovarr`, and generated GFF3 files use `mitovarr` in the source column.
