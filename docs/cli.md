# mitovarr CLI reference

Captured from `mitovarr 0.1.0` with unreleased fixes on 2026-10-12. Assembly commands are grouped into de novo assembly (`assemble-baited`, `assemble-seeded`) and mapping assembly (`assemble-mapping`). Current examples are in the [README](../README.md); method details are in the [documentation index](README.md). Regenerate these help sections after changing CLI options or defaults.

## Top-level commands

```text
Native Rust mtDNA de novo and mapping assembly, annotation and tRNA-Pro rotation

Usage: mitovarr <COMMAND>

Commands:
  fetch-caendr      Fetch mitochondrial reads and available mates from an indexed CaenDR BAM
  bait              Iteratively recruit original complete read pairs with seed kmers
  assemble-baited   BAM-baited de novo assembly from extracted pairs; no second recruitment
  assemble-seeded   Seeded de novo assembly; seed recruits reads only, validation uses all originals
  assemble-mapping  Mapping assembly from a reference and indexed BAM, with optional local reconstruction
  annotate          Predict and assess genes on frozen candidates independently of NCR1 completeness
  help              Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

## fetch-caendr

```text
Fetch mitochondrial reads and available mates from an indexed CaenDR BAM

Usage: mitovarr fetch-caendr [OPTIONS] --species <SPECIES> --strain <STRAIN> --out <OUT>

Options:
      --species <SPECIES>
          [possible values: cel, ctro, cbr]
      --strain <STRAIN>

      --out <OUT>

      --contig <CONTIG>
          Override automatic mitochondrial contig-name detection
      --timeout-seconds <TIMEOUT_SECONDS>
          Timeout for each HTTP request [default: 120]
      --max-templates <MAX_TEMPLATES>
          Bound the number of mitochondrial templates retained for mate recovery [default: 1000000]
  -h, --help
          Print help
```

## bait

```text
Iteratively recruit original complete read pairs with seed kmers

Usage: mitovarr bait [OPTIONS] --r1 <R1> --r2 <R2> --seed <SEED> --out <OUT>

Options:
      --r1 <R1>
      --r2 <R2>
      --seed <SEED>
      --out <OUT>
      --bait-k <BAIT_K>          [default: 21]
      --max-rounds <MAX_ROUNDS>  [default: 30]
      --max-pairs <MAX_PAIRS>    [default: 2000000]
      --max-kmers <MAX_KMERS>    Hard limit on distinct kmers in each recruitment or assembly set [default: 2000000]
      --min-baseq <MIN_BASEQ>    Optional hard base-quality cutoff; 0 disables it. Qualities still inform read scores [default: 0]
  -h, --help                     Print help
```

## assemble-baited

```text
BAM-baited de novo assembly from extracted pairs; no second recruitment

Usage: mitovarr assemble-baited [OPTIONS] --r1 <R1> --r2 <R2> --out <OUT>

Options:
      --r1 <R1>

      --r2 <R2>

      --out <OUT>

      --coverage-layout
          Also propose a coverage-fitted circular layout from sample graph paths; never certifies a final assembly
      --multi-k-layout
          Add cross-k layout proposals and read-context-guided walks; hypotheses remain unverified
      --k <K>
          Main odd k in 15..127; smaller k=47/31 graphs also propose open-end connections [default: 63]
      --min-kmer-count <MIN_KMER_COUNT>
          Minimum pair-record support for a retained kmer (mates counted once) [default: 3]
      --max-error-fraction <MAX_ERROR_FRACTION>
          Maximum weak/strong support ratio for tips, weak links and repeat routes; 0 disables all cleanup [default: 0.2]
      --max-bubble-fraction <MAX_BUBBLE_FRACTION>
          Maximum weak/strong mean coverage ratio for bounded local bubbles; does not relax repeat routes [default: 0.49]
      --min-depth <MIN_DEPTH>
          Minimum matching-base fragment depth after gapped validation [default: 10]
      --min-junction-fragments <MIN_JUNCTION_FRAGMENTS>
          Minimum named fragments supporting each repeat bridge and the output junction [default: 3]
      --max-validation-pairs <MAX_VALIDATION_PAIRS>
          Hard limit on unique names retained by final read validation [default: 2000000]
      --max-kmers <MAX_KMERS>
          Hard limit on distinct kmers in each recruitment or assembly set [default: 2000000]
      --min-baseq <MIN_BASEQ>
          Optional hard base-quality cutoff; 0 disables it. Qualities still inform read scores [default: 0]
  -h, --help
          Print help
```

## assemble-seeded

```text
Seeded de novo assembly; seed recruits reads only, validation uses all originals

Usage: mitovarr assemble-seeded [OPTIONS] --r1 <R1> --r2 <R2> --out <OUT> --seed <SEED>

Options:
      --r1 <R1>

      --r2 <R2>

      --out <OUT>

      --coverage-layout
          Also propose a coverage-fitted circular layout from sample graph paths; never certifies a final assembly
      --multi-k-layout
          Add cross-k layout proposals and read-context-guided walks; hypotheses remain unverified
      --k <K>
          Main odd k in 15..127; smaller k=47/31 graphs also propose open-end connections [default: 63]
      --min-kmer-count <MIN_KMER_COUNT>
          Minimum pair-record support for a retained kmer (mates counted once) [default: 3]
      --max-error-fraction <MAX_ERROR_FRACTION>
          Maximum weak/strong support ratio for tips, weak links and repeat routes; 0 disables all cleanup [default: 0.2]
      --max-bubble-fraction <MAX_BUBBLE_FRACTION>
          Maximum weak/strong mean coverage ratio for bounded local bubbles; does not relax repeat routes [default: 0.49]
      --min-depth <MIN_DEPTH>
          Minimum matching-base fragment depth after gapped validation [default: 10]
      --min-junction-fragments <MIN_JUNCTION_FRAGMENTS>
          Minimum named fragments supporting each repeat bridge and the output junction [default: 3]
      --max-validation-pairs <MAX_VALIDATION_PAIRS>
          Hard limit on unique names retained by final read validation [default: 2000000]
      --max-kmers <MAX_KMERS>
          Hard limit on distinct kmers in each recruitment or assembly set [default: 2000000]
      --min-baseq <MIN_BASEQ>
          Optional hard base-quality cutoff; 0 disables it. Qualities still inform read scores [default: 0]
      --seed <SEED>

      --bait-k <BAIT_K>
          [default: 21]
      --max-rounds <MAX_ROUNDS>
          [default: 30]
      --max-pairs <MAX_PAIRS>
          [default: 2000000]
  -h, --help
          Print help
```

## assemble-mapping

```text
Mapping assembly from a reference and indexed BAM, with optional local reconstruction

Usage: mitovarr assemble-mapping [OPTIONS] --ref <REFERENCE> --bam <BAM> --out <OUT>

Options:
      --local-assembly
          Build a sample graph from paired FASTQ and output a fast consensus with unresolved regions
      --deep-validation
          Enumerate local structures and compute complete paired-fragment scores (slower; opt-in)
      --best-effort
          Retain a usable consensus even if structure remains unresolved; enabled by default in fast local mode
      --coverage-layout
          Infer separate coverage-guided local layouts; bounded candidates never certify final structure
      --multi-k-layout
          Add cross-k local layout proposals with long flanks and shared read scoring
      --r1 <R1>
          Complete paired reads for local assembly (e.g. the matching fetch-caendr FASTQ)
      --r2 <R2>

      --local-k <LOCAL_K>
          Sample graph k for local assembly; independent of the de novo main k [default: 31]
      --local-max-span <LOCAL_MAX_SPAN>
          Maximum sequence length between local anchors, including the anchors [default: 2000]
      --local-max-paths <LOCAL_MAX_PATHS>
          Maximum local paths; reaching a search limit leaves the window unresolved [default: 128]
      --local-max-kmers <LOCAL_MAX_KMERS>
          Maximum distinct kmers per local graph; on overflow retain the full-BAM partial consensus [default: 2000000]
      --scoring-threads <SCORING_THREADS>
          Worker threads for deep validation; results keep input order. Capped by available CPUs [default: 4]
      --ref <REFERENCE>

      --bam <BAM>

      --contig <CONTIG>

      --out <OUT>

      --min-mapq <MIN_MAPQ>
          [default: 20]
      --min-baseq <MIN_BASEQ>
          Optional hard base-quality cutoff; 0 disables it. Quality-weighted voting and scoring remain [default: 0]
      --min-depth <MIN_DEPTH>
          [default: 10]
      --min-snp-frac <MIN_SNP_FRAC>
          Minimum weighted SNP fraction; a strict majority (> 0.5) is always required [default: 0.5]
      --min-indel-frac <MIN_INDEL_FRAC>
          Minimum fragment support fraction for an insertion or deletion [default: 0.7]
      --min-alt-frac <MIN_ALT_FRAC>
          Legacy override for both SNP and indel fractions; prefer the separate options
      --window <WINDOW>
          [default: 100]
      --max-records <MAX_RECORDS>
          Bound records retained for deterministic fragment-aware pileup [default: 500000]
      --subsample
          Deterministically subsample complete named fragments to the record cap
      --subsample-seed <SUBSAMPLE_SEED>
          Seed for fragment priority sampling (not a data integrity fingerprint) [default: 42]
      --emit-variants

      --emit-coverage

  -h, --help
          Print help
```

## annotate

```text
Predict and assess genes on frozen candidates independently of NCR1 completeness

Usage: mitovarr annotate [OPTIONS] --fasta <FASTA> --out <OUT>

Options:
      --fasta <FASTA>
          Candidate/consensus FASTA (all records, plain or gzip); never fed back to assembly

      --species <SPECIES>
          Bundled, pinned CDS/rRNA/tRNA homology models

          [possible values: cel, ctro, cbr]

      --gene-models <GENE_MODELS>
          Custom mitovar-rust.gene-models.v1 JSON panel (translation table 5)

      --out <OUT>


      --topology <TOPOLOGY>
          Possible values:
          - auto:     Use annotation_topology when present, otherwise a literal circular=true header
          - linear
          - circular: Treat each input as a circular hypothesis; this does not validate closure

          [default: auto]

      --rotate-to-trna-pro
          Orient a circular input so its unique, complete tRNA-Pro starts at base 1 on +

      --assembly-ambiguity <ASSEMBLY_AMBIGUITY>
          Project assembly ambiguity intervals from matching FASTA coordinates

      --ncr-min-len <NCR_MIN_LEN>
          Minimum intergenic length for NCR intervals; shorter named loci remain in ncr_loci.tsv

          [default: 20]

  -h, --help
          Print help (see a summary with '-h')
```
