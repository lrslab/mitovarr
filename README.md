# mitovarr

Native Rust tools for mitochondrial read recruitment, **BAM-baited de novo assembly**, reference-guided consensus and local reconstruction, gene annotation, and non-coding region (NCR) analysis. The Cargo package, Rust library and executable are all named `mitovarr`.

Read extraction, gzip I/O, assembly, consensus and annotation run in-process. No SPAdes, external mapper, samtools or external annotation program is required. Mapping uses existing BAM alignments; this package does not globally align FASTQ.

[Workflow Figure 1](docs/figure_1_workflow/README.md) · [Documentation](docs/README.md) · [CLI reference](docs/cli.md)

## Build and install

The package declares Rust 1.90 or newer, edition 2021. A C compiler toolchain is needed to build the TLS dependency `ring`.

```bash
cargo build --release --locked
./target/release/mitovarr --help

cargo install --path . --locked
mitovarr --version
```

Python and Poppler are optional tools for rebuilding documentation figures. They are not Rust build or runtime dependencies.

## Workflows and input scope

| Command | Inputs | Scope |
| --- | --- | --- |
| `fetch-caendr` | Species and strain | Query the indexed mitochondrial BAM region and recover available mates; no complete-BAM download |
| `assemble-baited` | BAM-extracted paired FASTQ | BAM-baited de novo assembly; no second recruitment; validation uses the extracted pairs |
| `assemble-denovo` | Original complete paired FASTQ | Reference-blind graph assembly and validation |
| `assemble-seeded` | Original paired FASTQ plus explicit seed FASTA | Seed recruits complete pairs only; assembly uses sample sequence and validation uses all original pairs |
| `bait` | Original paired FASTQ plus explicit seed FASTA | Recruitment alone, with convergence reported |
| `polish` / `assemble-mapping` | Matching reference FASTA plus coordinate-sorted, indexed BAM | Fragment-aware consensus; optional matching FASTQ enables local structural reconstruction |
| `annotate` | Frozen FASTA plus species or gene models | Homology annotation, optional tRNA-Pro rotation and positional NCR labels |

`assemble-mapping` is an alias of `polish`. `cel`, `cbr` and `ctro` identify *C. elegans*, *C. briggsae* and *C. tropicalis*. FASTA/FASTQ can be plain text or gzip. BAM requires a `.bam.bai` or `.bai` index; CRAM is not supported. Each output directory must be new or empty.

## CaenDR to BAM-baited de novo assembly

```bash
mitovarr fetch-caendr --species cbr --strain AF16 --out runs/AF16_reads

mitovarr assemble-baited \
  --r1 runs/AF16_reads/reads.R1.fastq.gz \
  --r2 runs/AF16_reads/reads.R2.fastq.gz --out runs/AF16_denovo
```

The extraction writes paired/singleton/unpaired FASTQ, `mt.bam`, `mt.bam.bai` and reports. It uses HTTP Range requests and the remote index. These inputs have **`mt-bam-extracted`** scope: the upstream mappings affect which reads are available, even though the graph assembler receives only FASTQ. Unplaced or unavailable mates and hard-clipped sequence cannot be recovered.

The packaged [`scripts/fetch_caendr.sh`](scripts/fetch_caendr.sh) builds and launches the same native extractor through Cargo. It forwards arguments unchanged and can be called from another working directory. See [CaenDR extraction](docs/caendr.md).

For original whole-sample reads:

```bash
mitovarr assemble-denovo --r1 original_R1.fastq.gz --r2 original_R2.fastq.gz --out runs/blind

mitovarr assemble-seeded --r1 original_R1.fastq.gz --r2 original_R2.fastq.gz \
  --seed mitochondrial_seed.fasta --out runs/seeded
```

Incomplete seed recruitment cannot publish a final assembly. The main graph defaults to `--k 63`, baseQ 20 and minimum pair-record k-mer support 3. Odd k values from 15 to 127 are accepted; larger k is not inherently better for short reads. Native k=47/31 end searches are used only when smaller than the main k. See [assembly](docs/de_novo.md) and [paired-path scoring](docs/paired_paths.md).

## Mapping and local reconstruction

The reference must match the BAM contig name, length and coordinate system. Select `--contig` explicitly for a multi-record reference.

```bash
mitovarr assemble-mapping \
  --ref reference.fasta --bam mt.bam --contig MtDNA \
  --emit-coverage --emit-variants --out runs/mapping
```

Ordinary SNP calls require a strict quality-weighted fragment majority (`--min-snp-frac 0.5`); ties remain N. Indels use `--min-indel-frac 0.7`. Overlapping mates count once, repeat-aware rules are always enabled, and unsupported bases remain N rather than being filled from the reference. `base_votes.tsv` retains mixed-allele evidence. The legacy `--min-alt-frac` sets both fractions and conflicts with explicitly setting either separate option.

Enable local sample-graph reconstruction with matching paired FASTQ:

```bash
mitovarr assemble-mapping \
  --ref reference.fasta --bam mt.bam --contig MtDNA \
  --local-assembly --r1 reads.R1.fastq.gz --r2 reads.R2.fastq.gz \
  --emit-coverage --emit-variants --out runs/mapping_local
```

Local candidates can contain sample sequence absent from the reference. Whole-contig comparisons retain homologous source copies and use the same paired fragments across candidates; both replacement junctions need support. Search limits or unresolved competition remain explicit. This does not guarantee reconstruction of every NCR or repeat. See [local assembly and coordinate contracts](docs/mapping_local_assembly.md).

For high-depth BAMs, explicitly add `--subsample --max-records 1000000` (default sampling seed 42, adjustable by `--subsample-seed`). Selection retains all target-BAM records for each selected QNAME under the record cap. With local assembly, the corresponding complete FASTQ pairs are saved as `subsampled.R1.fastq` / `subsampled.R2.fastq` and used for the graph and scoring. A retained BAM fragment lacking a complete supplied pair contributes only pileup evidence. Counts and scope are reported under `report.json.subsampling`; coverage/variant fractions describe the retained evidence. Without `--subsample`, exceeding `--max-records` is an error, not silent truncation.

## Gene annotation, rotation and NCRs

```bash
mitovarr annotate \
  --fasta runs/AF16_denovo/primary.fasta --species cbr \
  --assembly-ambiguity runs/AF16_denovo/primary_ambiguity.tsv \
  --rotate-to-trna-pro --out runs/AF16_annotation
```

Annotation uses nucleotide homology on both strands, followed by CDS coding checks with translation table 5 and RNA boundary checks. Each bundled species panel contains 12 CDS, 22 tRNA and 2 rRNA models; detected, partial and full-length models are reported separately. RNA structure is not assessed. CDS with an unsupported 3′ boundary report `cds.stop=partial_or_uncertain`, including fragments ending in T/TA.

Rotation needs an input interpreted as circular and a unique complete tRNA-Pro. Auto topology uses `annotation_topology` when present, otherwise a literal `circular=true` header. `--topology circular` explicitly sets an annotation hypothesis; it does not validate closure. Features and compatible assembly-ambiguity coordinates follow the rotation.

`--assembly-ambiguity` accepts de novo `primary_ambiguity.tsv`, local mapping `structural_ambiguity.tsv`, or a previously projected annotation table **with the matching FASTA**. Ordinary pileup `ambiguity.tsv` uses reference coordinates and is not accepted directly.

NCR1-4 are positional labels derived from observed flanking genes. Short/zero-length gaps, missing flanks and open ends remain in `ncr_loci.tsv`; a species need not have four long NCRs. Annotation does not alter assembly selection or resolve repeats, and automatic cross-sample NCR SV classification is not implemented. See [annotation](docs/annotation.md), [中文规则](docs/annotation.zh-CN.md) and [NCR definitions](docs/ncr.md).

## Outputs and status

| Output | Interpretation |
| --- | --- |
| `primary.fasta`, `primary_ambiguity.tsv` | De novo representative and localized uncertainty; may remain unresolved |
| `candidates.fasta`, path candidate FASTAs | Alternative de novo sequences with selection evidence in JSON |
| `bubble_alternatives.fasta` | Retained removed/kept spellings when permissive de novo cleanup cannot justify final status |
| `assembly_graph.gfa`, `path_graph.gfa` | Cleaned unitigs and the path-enumeration graph, respectively |
| `candidate.fasta` | Mapping consensus candidate, including accepted local replacements when enabled |
| `reference_candidate.fasta`, `local_assembly.gfa` | Local mode's original pileup backbone and local candidate-path graph |
| `structural_variants.tsv`, `structural_ambiguity.tsv` | Local mapping decisions and unresolved intervals in candidate coordinates |
| `consensus.fasta` | Only when the workflow's final-output rules pass; mapping completion does not verify circularity |
| `report.json`, evidence TSVs | Input scope, thresholds, coverage, votes, selection and uncertainty |
| `annotated.fasta`, GFF3, gene/protein FASTAs, NCR tables/FASTA | Independent post-assembly features, sequences and locus states |

De novo final validation requires a uniquely supported main circle, matching-base coverage, supported joins and no unresolved permissively collapsed bubble in its original graph component; small unrelated components do not automatically veto it. Mapping final output requires no N or remaining flagged ambiguity and, in local mode, resolved windows. Local contractions of at least `--local-k` bases also require linked unique-flank evidence from at least three named fragments. Full gene annotation or an N-free sequence does not by itself establish NCR structure.

De novo bubble cleanup compares whole-path mean coverage (`--max-bubble-fraction 0.49`) within bounded local alternatives; tip and repeat-route thresholds remain separately controlled by `--max-error-fraction 0.2`. Mapping local assembly can pool smaller-k read-graph connections across a main-k break, keeping the original full flanks and all fragment-validation requirements. See the method guides for limits and biological-allele caveats.

Exit codes: **0** completed under the command's rules; **2** unresolved assembly/recruitment or extraction with no complete pairs; **1** input/runtime error. Annotation can exit 0 while reporting missing genes or a skipped rotation. CLI usage errors also use Clap's exit code 2; inspect the message/report to distinguish them.

## Development

```bash
cargo fmt --check
cargo test --locked
cargo clippy --all-targets --all-features --locked -- -D warnings
```

The former package/binary names were `mitovar-rust` / `mtasm-rust`; current commands and Rust imports use `mitovarr`. Existing `mitovar-rust.*` JSON schema identifiers and saved gene panels remain compatible. New GFF3 source fields use `mitovarr`. The checkout directory name is independent of the package name.
