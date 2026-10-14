# mitovarr

Native Rust tools for mitochondrial **de novo assembly** and **mapping assembly**, read recruitment, gene annotation, and non-coding region (NCR) analysis. The Cargo package, Rust library and executable are all named `mitovarr`.

Read extraction, gzip I/O, assembly, consensus and annotation run in-process. No SPAdes, external mapper, samtools or external annotation program is required. Mapping uses existing BAM alignments; this package does not globally align FASTQ.

Partial mitochondrial sequences are supported deliverables: NCR1 completion is
not required to retain and evaluate the other regions. Use `primary.fasta`
(de novo) or `candidate.fasta` (mapping), their uncertainty reports, and the
independent `annotate` gene checks. Mapping `--best-effort` can also write a
best-effort `consensus.fasta`; its uncertainty and strict validation status are
reported separately. Complete-assembly claims retain their evidence requirements.
See [partial results and gene quality](docs/partial_results.md).

De novo runs also select the best available eligible circular candidate into
`best_circular.fasta`, with matching `best_circular_ambiguity.tsv` and
`report.json.best_circular`. This gives a single ring for downstream review even
when it remains unverified. Within a complete-cycle family, the existing paired
fragment score chooses the representative; distinct k-mer support mass ranks
different components. `primary.fasta` still represents the main component,
which can be linear. If no eligible circle exists, no best-circular FASTA is
written. See [circular candidate selection](docs/paired_paths.md#best-circular-candidate).

Open-end searches divide their bounded effort between reachable target ends.
`end_bridge_proposals.fasta` retains discovered connections for inspection,
including proposals withheld after an incomplete search. Only reciprocal,
eligible connections enter `end_bridge_candidates.fasta` and the path graph;
a saved proposal is not a validated join or a complete circular assembly.

Path anchors prefer unique 31-mer contexts and use longer sample-only contexts
when short graph fragments obscure uniqueness. Join evidence is unchanged.
Mapping maximum consensus fills only N at the existing span; hypotheses that
change called bases remain separate. Optional local graph/endpoint resource
limits retain usable partial outputs and explicitly defer structural validation.
See [mapping limits](docs/mapping_local_assembly.md) and
[paired-path anchors](docs/paired_paths.md).

[Workflow Figure 1](docs/figures/README.md) · [Documentation](docs/README.md) · [CLI reference](docs/cli.md)

## Package overview

[![Figure 1. Five schematic panels explain input scope, sample-derived repeat alternatives, shared-fragment comparison and optional copy-depth inference, three expected output states including retained partial sequences, and independent biological annotation. Depth-inferred layouts remain unverified; unresolved structures do not produce a final consensus.](docs/figures/package_workflow.png)](docs/figures/package_workflow.pdf)

**Figure 1 | Evidence-aware reconstruction of repeat-rich mitochondrial genomes with mitovarr.** **a,** Distinct input routes preserve their read-evidence scopes. **b,** Sample graphs and local reconstruction retain alternative repeat structures. **c,** Bounded candidate comparisons use the same fragments; optional multi-k layouts use copy-depth fit while retaining competing evidence. **d,** Expected outputs distinguish a consensus that passes strict checks, partial or unresolved candidates available for evaluation, and an unverified depth-inferred hypothesis. **e,** Frozen sequences undergo gene-quality checks, eligible tRNA-Pro normalization and positional NCR interpretation without changing assembly status. NCR1 completeness is not a gene-quality requirement. All sequences and depth values are schematic, not experimental results.

[Vector PDF](docs/figures/package_workflow.pdf) · [Editable SVG](docs/figures/package_workflow.svg) · [Full legend and drawing source](docs/figures/README.md) · [Detailed method schematic](docs/figure_1_workflow/README.md)

## Build and install

### Prebuilt binaries

Download the archive for your platform from
[GitHub Releases](https://github.com/lrslab/mitovarr/releases).
Each release tagged `v<version>` provides:

| Platform | Archive |
| --- | --- |
| macOS M series / Apple Silicon (ARM64, macOS 11+) | `mitovarr-v<version>-aarch64-apple-darwin.tar.gz` |
| Linux x86_64 (static musl) | `mitovarr-v<version>-x86_64-unknown-linux-musl.tar.gz` |
| Linux ARM64 / AArch64 (glibc 2.35+, e.g. Ubuntu 22.04+) | `mitovarr-v<version>-aarch64-unknown-linux-gnu.tar.gz` |

For example, after downloading a release for an Apple Silicon Mac, replace
`<version>` below with its version number:

```bash
tar -xzf mitovarr-v<version>-aarch64-apple-darwin.tar.gz
cd mitovarr-v<version>-aarch64-apple-darwin
./mitovarr --version
./mitovarr --help
mkdir -p "$HOME/.local/bin"
install -m 755 mitovarr "$HOME/.local/bin/mitovarr"
export PATH="$HOME/.local/bin:$PATH"
```

Use the corresponding archive name on Linux. Add the `PATH` line to your shell
configuration to keep it after restarting the terminal. The archives include
the usage guides and figures. Gene models are embedded in the executable;
running it requires no Rust, Cargo, Python or external bioinformatics programs.
The macOS binary is not Developer ID signed or notarized.

### Build from source

The package declares Rust 1.90 or newer, edition 2021. A C compiler toolchain is needed to build the TLS dependency `ring`.

```bash
git clone https://github.com/lrslab/mitovarr.git
cd mitovarr
cargo build --release --locked
./target/release/mitovarr --help

cargo install --path . --locked
mitovarr --version
```

Rebuilding the overview figure requires Python, ReportLab and Arial font files; PNG export additionally uses Poppler (`pdftoppm`). These are optional documentation tools, independent of the Rust build and runtime. See the [figure reproduction instructions](docs/figures/README.md#artwork-and-reproduction).

## Workflows and input scope

The two assembly workflows are:

- **De novo assembly** builds sequence from paired sample reads. Use
  `assemble-baited` for BAM-extracted pairs or `assemble-seeded` for seed-based
  recruitment from original FASTQ.
- **Mapping assembly** uses a matching reference and existing BAM alignments.
  `assemble-mapping` generates the consensus and can reconstruct uncertain
  regions with matching paired FASTQ.

Read extraction/recruitment prepares inputs; annotation assesses the resulting
sequences. The commands retain their distinct input and validation scopes:

| Command | Inputs | Scope |
| --- | --- | --- |
| `fetch-caendr` | Species and strain | Query the indexed mitochondrial BAM region and recover available mates; no complete-BAM download |
| `assemble-baited` | BAM-extracted paired FASTQ | BAM-baited de novo assembly; no second recruitment; validation uses the extracted pairs |
| `assemble-seeded` | Original paired FASTQ plus explicit seed FASTA | Seed recruits complete pairs only; assembly uses sample sequence and validation uses all original pairs |
| `bait` | Original paired FASTQ plus explicit seed FASTA | Recruitment alone, with convergence reported |
| `assemble-mapping` | Matching reference FASTA plus coordinate-sorted, indexed BAM | Mapping assembly: fragment-aware consensus with optional local reconstruction from matching FASTQ |
| `annotate` | Frozen FASTA plus species or gene models | Homology annotation, optional tRNA-Pro rotation and positional NCR labels |

`cel`, `cbr` and `ctro` identify *C. elegans*, *C. briggsae* and *C. tropicalis*. FASTA/FASTQ can be plain text or gzip. BAM requires a `.bam.bai` or `.bai` index; CRAM is not supported. Each output directory must be new or empty.

Paired FASTQ must use separate, synchronized R1/R2 files with equal, nonzero record counts, matching normalized read names, four lines per record and Phred+33 qualities. Interleaved or wrapped FASTQ is not supported.

Hard base-quality filtering is disabled by default (`--min-baseq 0`) in recruitment, de novo assembly and mapping. Low-quality bases no longer break otherwise supported k-mers, and alignment seeds have no implicit quality cutoff. Read alignment and mapping consensus still use base qualities as scoring weights. An explicit positive `--min-baseq` enables the optional cutoff; read-support, ambiguity, repeat and closure checks remain required.

## CaenDR to BAM-baited de novo assembly

```bash
mitovarr fetch-caendr --species cbr --strain AF16 --out runs/AF16_reads

mitovarr assemble-baited \
  --r1 runs/AF16_reads/reads.R1.fastq.gz \
  --r2 runs/AF16_reads/reads.R2.fastq.gz --out runs/AF16_denovo
```

The extraction writes paired/singleton/unpaired FASTQ, `mt.bam`, `mt.bam.bai` and reports. It uses HTTP Range requests and the remote index. These inputs have **`mt-bam-extracted`** scope: the upstream mappings affect which reads are available, even though the graph assembler receives only FASTQ. Unplaced or unavailable mates and hard-clipped sequence cannot be recovered.

The packaged [`scripts/fetch_caendr.sh`](scripts/fetch_caendr.sh) launches the included binary in a release archive, or builds and launches the native extractor through Cargo in a source checkout. It forwards arguments unchanged and can be called from another working directory. See [CaenDR extraction](docs/caendr.md).

For original whole-sample reads, recruit mitochondrial pairs using a seed:

```bash
mitovarr assemble-seeded --r1 original_R1.fastq.gz --r2 original_R2.fastq.gz \
  --seed mitochondrial_seed.fasta --out runs/seeded
```

Incomplete seed recruitment cannot publish a final assembly. The main graph defaults to `--k 63`, no hard base-quality cutoff (`--min-baseq 0`) and minimum pair-record k-mer support 3. Odd k values from 15 to 127 are accepted; larger k is not inherently better for short reads. Native k=47/31 end searches are used only when smaller than the main k. Anchor-path searches schedule nearby exits before longer repeated walks; candidate length is not a selection score. Path, visited-state and pending-state limits remain explicit, and incomplete searches cannot supply a resolved join. See [assembly](docs/de_novo.md) and [paired-path scoring](docs/paired_paths.md).

Add `--multi-k-layout` together with `--coverage-layout` in de novo mode or local mapping with `--deep-validation` to extend the hypothesis pool across k=15/21/31/47/63/95/127. Additional searches retain complete long flanks and use longer observed read windows to guide shorter-k branch traversal. The full generated pool defines common-fragment eligibility; at most eight sample hypotheses plus the baseline receive whole-layout fragment scores. Unscored families remain explicit, and agreement between k values does not multiply read support. This option increases work and retains unresolved alternatives. It does not change the strict primary/final assembly workflow.

For an explicit coverage-based hypothesis across an unresolved repeat, add
`--coverage-layout` to a de novo assembly command. This writes a separate
`coverage_layout.fasta` when the sample graph connects both backbone flanks,
plus alternative layouts, GFA walks, copy-depth fit and uncertainty reports.
A reported singleton-edge fallback can force a layout when ordinary graph
support filtering disconnects the corridor.
Coverage estimates repeat usage; it does not certify repeat order or final
assembly status. See [coverage layouts](docs/de_novo.md#coverage-inferred-layouts).

## Mapping assembly

The reference must match the BAM contig name, length and coordinate system. Select `--contig` explicitly for a multi-record reference. Pileup accumulation uses circular target coordinates, including supported alignment/soft-clip evidence across the origin; this coordinate convention does not validate circular closure.

```bash
mitovarr assemble-mapping \
  --ref reference.fasta --bam mt.bam --contig MtDNA \
  --emit-coverage --emit-variants --out runs/mapping
```

Ordinary SNP calls require a strict quality-weighted fragment majority (`--min-snp-frac 0.5`); ties remain N. Indels use `--min-indel-frac 0.7`. Default minimum MAPQ, baseQ and depth are 20, 0 and 10; baseQ 0 disables hard base-quality filtering. Overlapping mates count once, repeat-aware rules are always enabled, and unsupported bases remain N rather than being filled from the reference. The legacy `--min-alt-frac` sets both fractions and conflicts with explicitly setting either separate option.

`base_votes.tsv` reports retained mixed sites meeting `--min-depth` where at least two A/C/G/T alleles each contribute at least 5% of total base weight. Base weights are `min(baseQ,40) * min(MAPQ,60)`. Mixed-site evidence, including a `heteroplasmy` flag in `variants.tsv`, does not establish biological heteroplasmy.

For a fast consensus and uncertainty report, supply matching paired FASTQ:

```bash
mitovarr assemble-mapping \
  --ref reference.fasta --bam mt.bam --contig MtDNA \
  --local-assembly --r1 reads.R1.fastq.gz --r2 reads.R2.fastq.gz \
  --emit-coverage --emit-variants --out runs/mapping_local
```

Local mode now defaults to fast consensus. It uses every supplied complete pair
in the sample graph, preserves called bases and attempts N fills only at the
existing span. It writes `consensus.fasta`, `maximum_consensus.fasta` and matching
`maximum_consensus.ambiguity.tsv` without waiting for structural enumeration or
paired candidate scoring. Reports have `validation_mode=fast_consensus`,
`scoring_status=not_run_fast_consensus` and null unmeasured scores/counts.
The sequence is best effort; strict validation and circularity remain unverified.

Add `--deep-validation` to run the previous full local searches and comparisons
in a new output directory. Deep candidates can contain sample sequence absent
from the reference. Whole-contig comparisons retain homologous source copies
and use the same paired fragments across candidates; both replacement junctions
need support. Search limits and unresolved competition remain explicit. See
[local assembly and coordinate contracts](docs/mapping_local_assembly.md).

In deep mode, whole-layout scoring reuses identical affine alignment subproblems across
candidate windows, removes provably dominated internal edit states, and
aggregates mate combinations through endpoint indexes.
All supplied reads, scoring rules, structural evidence and logical search
limits are retained; cache capacity does not change an unknown result.
An affine fast path certifies legal endpoint optima before using the full
edit-count DP for the remaining endpoints. Mapping fragment scores use four
workers by default, capped by available CPUs; `--scoring-threads 1` makes
scoring serial. Scores and diagnostics are always accumulated in input order.

Fast local mode enables best-effort output automatically. In deep mode, add
`--best-effort` to `assemble-mapping --local-assembly --deep-validation`.
It fills N positions from sample-read
graphs at k=21/17/13/11, first preserving called bases and the existing span.
Repeat use is bounded by that span rather than a fixed copy ceiling. Legacy
length-guided paths remain available as fallback proposals in deep mode.

This mode writes `consensus.fasta` and `maximum_consensus.fasta` even when
structural validation is incomplete. A successful write exits 0; `status` is
`best_effort`, `consensus_written` is true, and `strict_validation_passed`
records validation separately. `final_published` retains its prior strict
meaning. In deep mode, read scores describe the proposals and do not veto the primary
output. `maximum_consensus.ambiguity.tsv` marks uncertainty in the primary
sequence's coordinates. `candidate.fasta` remains the strict comparison.
The reference-guided span is a prior and does not prove repeat copy number.
See the [maximum consensus method](docs/mapping_local_assembly.md#maximum-consensus-and-legacy-bridges).

Add `--deep-validation --coverage-layout` to a local mapping run to compare bounded repeat layouts
even when the strict search is incomplete. The shortlist retains the best
R1 copy-depth fit, then represents repeat-count and traversal-order families
before additional sequence variants. Up to eight sample paths plus the baseline
are compared with the same paired fragments; retained window choices share
a copy budget. `coverage_layout.fasta` and its
alternatives are explicitly unverified; `candidate.fasta` and strict final
publication remain independent. See [mapping coverage layouts](docs/mapping_local_assembly.md#coverage-guided-layout-hypotheses).

`--max-records` defaults to **500000** target-BAM records. Without `--subsample`, exceeding the cap is an error. In local mode the same option also caps distinct FASTQ pair names; it does not limit only the BAM scan.

Only after a full-input run hits a documented input/memory limit, explicitly add `--subsample --max-records 1000000` (default sampling seed 42, adjustable by `--subsample-seed`). Selection retains all target-BAM records for each selected QNAME under the record cap. With local assembly, the corresponding complete FASTQ pairs are saved as `subsampled.R1.fastq` / `subsampled.R2.fastq` and used for the graph and scoring. A retained BAM fragment lacking a complete supplied pair contributes only pileup evidence. Counts and scope are reported under `report.json.subsampling`; coverage/variant fractions describe the retained evidence.

De novo `assemble-baited` / `assemble-seeded` accept the same explicit `--subsample` (seed `--subsample-seed`, default 42). It acts only after all complete pairs exceed `--max-kmers`, then halves the retained pairs (down to 1/64) until the graph fits; final read validation still uses all pairs. See [de novo subsampling](docs/de_novo.md#subsampling-after-a---max-kmers-overflow).

## Gene annotation, rotation and NCRs

```bash
mitovarr annotate \
  --fasta runs/AF16_denovo/primary.fasta --species cbr \
  --assembly-ambiguity runs/AF16_denovo/primary_ambiguity.tsv \
  --rotate-to-trna-pro --out runs/AF16_annotation
```

Annotation uses nucleotide homology on both strands, followed by CDS coding checks with translation table 5 and RNA boundary checks. Each bundled species panel contains 12 CDS, 22 tRNA and 2 rRNA models; detected, partial and full-length models are reported separately. RNA structure is not assessed. CDS with an unsupported 3′ boundary report `cds.stop=partial_or_uncertain`, including fragments ending in T/TA.

Rotation needs an input interpreted as circular, a unique complete tRNA-Pro with an aligned model 5′ boundary, and a completed anchor search. Auto topology uses `annotation_topology` when present, otherwise a literal `circular=true` header; unmarked records are treated as linear. `--topology circular` explicitly sets an annotation hypothesis; it does not validate closure. Features and compatible assembly-ambiguity coordinates follow the rotation.

The example uses `primary.fasta` because its header records `circular` and `final_verified` and its record ID matches `primary_ambiguity.tsv`. Ordinary de novo final and mapping strict outputs have unmarked headers. Mapping best-effort consensus headers identify their mode and validation state. Explicitly supply `--topology circular` for an assumed circular mapping sequence when rotating it. Changing FASTA files also requires a matching ambiguity table; record IDs and coordinates must agree.

`--assembly-ambiguity` accepts de novo `primary_ambiguity.tsv` or `best_circular_ambiguity.tsv`, local mapping `structural_ambiguity.tsv`, or a previously projected annotation table **with the matching FASTA**. Ordinary pileup `ambiguity.tsv` uses reference coordinates and is not accepted directly.

NCR1-4 are positional labels derived from observed flanking genes. Short/zero-length gaps, missing flanks and open ends remain in `ncr_loci.tsv`; a species need not have four long NCRs. Annotation does not alter assembly selection or resolve repeats, and automatic cross-sample NCR SV classification is not implemented. See [annotation](docs/annotation.md), [中文规则](docs/annotation.zh-CN.md) and [NCR definitions](docs/ncr.md).

## Outputs and status

Files are workflow-dependent; an incomplete or failed stage may produce only a subset.

| Output | Interpretation |
| --- | --- |
| `primary.fasta`, `primary_ambiguity.tsv` | De novo representative and uncertainty intervals, which may cover the entire sequence; may remain unresolved. No primary FASTA is written when no candidate reaches the main k length |
| `best_circular.fasta`, `best_circular_ambiguity.tsv` | Best available eligible de novo ring and matching uncertainty; absent when no eligible ring exists. Selection does not grant verified-complete status |
| `candidates.fasta`, path candidate FASTAs | Alternative de novo sequences with selection evidence in JSON |
| `end_bridge_proposals.fasta` | Discovered end connectors, including withheld proposals; eligibility and limits are in JSON |
| `bubble_alternatives.fasta` | Retained removed/kept spellings when permissive de novo cleanup cannot justify final status |
| `assembly_graph.gfa`, `path_graph.gfa` | Cleaned unitigs and the path-enumeration graph, respectively |
| `candidate.fasta` | Mapping consensus candidate, including accepted local replacements when enabled |
| `reference_candidate.fasta`, `local_assembly.gfa` | Local mode's original pileup backbone and local candidate-path graph |
| `structural_variants.tsv` | Local mapping decisions; `ref_window_start/end` are original reference coordinates and include flanks |
| `structural_ambiguity.tsv` | Unresolved local intervals in the output `candidate.fasta` coordinates |
| `consensus.fasta` | Strict result, or a best-effort primary with mapping `--best-effort`; inspect `status` and `strict_validation_passed` |
| `maximum_consensus.fasta`, `maximum_consensus.ambiguity.tsv` | Mapping best-effort mode's primary consensus and matching uncertainty, including unresolved fills |
| `report.json`, evidence TSVs | Input scope, thresholds, coverage, votes, selection and uncertainty |
| `gene_quality.tsv` | Annotation checks for every expected gene; NCR1 completeness is not a gate |
| `annotated.fasta`, GFF3, gene/protein FASTAs, NCR tables/FASTA | Independent post-assembly features, sequences and locus states |

Mapping `variants.tsv`, `coverage.tsv`, `ambiguity.tsv`, `uncertainty.tsv` and `base_votes.tsv` describe the original reference-coordinate pileup, even after accepted indels or local replacements change the candidate. They are not recomputed on locally inserted sequence or transformed by annotation rotation. Local `baseline_segments` refer to `reference_candidate.fasta`; `output_segments` refer to `candidate.fasta`. TSV positions and JSON segments are one-based inclusive unless explicitly named as offsets.

De novo final validation requires a uniquely supported main circle, matching-base coverage, supported joins and no unresolved permissively collapsed bubble in its original graph component; small unrelated components do not automatically veto it. Passing these checks does not establish mitochondrial identity; annotation is an independent subsequent step. Mapping strict validation requires no N or remaining flagged ambiguity and, in local mode, resolved windows; best-effort consensus availability is separate. Local contractions of at least `--local-k` bases also require linked unique-flank evidence from at least three named fragments. Full gene annotation or an N-free sequence does not by itself establish NCR structure.

De novo bubble cleanup compares whole-path mean coverage (`--max-bubble-fraction 0.49`) within bounded local alternatives; tip and repeat-route thresholds remain separately controlled by `--max-error-fraction 0.2`. Mapping local assembly can pool smaller-k read-graph connections across a main-k break, keeping the original full flanks and all fragment-validation requirements. See the method guides for limits and biological-allele caveats.

Exit codes: **0** completed under the command's rules, including a written mapping best-effort consensus; **2** unresolved strict assembly/recruitment or extraction with no complete pairs; **1** input/runtime error. Annotation can exit 0 while reporting missing genes or a skipped rotation. CLI usage errors also use Clap's exit code 2; inspect the message/report to distinguish them.

For partial workflows, inspect the assembly report's additive `sequence_output`
field. `partial_available` identifies a retained sequence, not a quality pass;
`unavailable` (de novo) or `no_called_bases` (mapping) cannot be treated as a useful
partial result. Exit code 2 alone is not permission to continue. Existing
strict output contracts remain in force outside the explicit mapping best-effort
mode described above. Annotation
adds `records[].gene_quality` and does not promote assemblies to complete.

## Development

```bash
cargo fmt --check
cargo test --release --locked
cargo clippy --all-targets --all-features --locked -- -D warnings
```

GitHub Actions runs these checks on Linux x86_64, Linux ARM64 and macOS Apple
Silicon. Pushing a tag that matches the Cargo version runs the checks again,
builds and smoke-tests all three binary archives, then publishes a GitHub
Release with changelog notes. See [publishing releases](docs/publishing.md)
for versioning, local packaging checks and the tag procedure.

The former package/binary names were `mitovar-rust` / `mtasm-rust`; current commands and Rust imports use `mitovarr`. Existing `mitovar-rust.*` JSON schema identifiers and saved gene panels remain compatible. New GFF3 source fields use `mitovarr`. The checkout directory name is independent of the package name.

Optional coverage layouts include common-fragment whole-sequence comparison and
separate observed-junction versus mate-bracket evidence. De novo layouts also
receive full matching-base depth validation. These reports retain the inferred
status; read fit alone does not establish repeat ordering or copy number.

Whole-hypothesis comparisons freeze the same observed flank pivot for every
candidate and export compressed fragment-level score contributions. `N` bases
never count as observed matching junction sequence.

Coverage comparisons preserve endpoint and edit states, enforce the shared pivot
inside alignment, and report incomplete searches as unknown for every candidate.
The eight-path shortlist first covers repeat-count and traversal-order families,
with the baseline added separately; its cohort is frozen on all generated paths.
Mapping choices use the explicit worst-case unknown-fragment score bound when
filtering opposed paths. These remain inferred hypotheses, with unresolved
structure and search tails reported. See the
[shared-evidence limits](docs/paired_paths.md#shared-whole-layout-read-evidence).

Disconnected low-support graph shadows are removed using retained sample k-mer
support; isolated anchors require 4k bases. See [de novo cleanup](docs/de_novo.md)
for the threshold, report fields and low-frequency-homolog limitation.

Partial primary paths use shared-flank allele evidence and reciprocal candidate
comparisons while retaining strict complete-publication checks. The additional
report fields and uncertainty rules are documented in `docs/paired_paths.md`.
