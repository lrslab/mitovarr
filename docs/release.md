# Development and release

Package, library and binary: **`mitovarr`**. Initial version: **0.1.0 (unreleased)**.

## Source repository

This directory is the independent Git root. Track `Cargo.lock`, source, bundled gene models, regression tests, portable fixtures and documentation. The checkout directory name does not determine the Cargo package or binary name.

Experimental reports, sample result tables and research logs belong outside this repository and source package. Generated sequencing runs and analysis outputs also stay in the separate local development tree; build caches, binaries and release archives are ignored. Tests construct small BAM/FASTQ inputs locally and use bundled annotation fixtures; they do not require the real-read analysis directories or the CaenDR service. CaenDR HTTP behavior is tested against a local server.

[Current guides](README.md) describe the implementation. [Figure 1](figure_1_workflow/README.md) illustrates de novo, optional local mapping, annotation and NCR interpretation.

## Checks and source package

From this package root:

```bash
cargo fmt --check
cargo test --locked
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo package --locked
```

Before the initial commit, use `cargo package --locked --allow-dirty`. Packaging builds the extracted source to check that the archive is independently usable. The archive is `target/package/mitovarr-0.1.0.crate`; `CARGO_TARGET_DIR` relocates build outputs. The manifest includes source, models, launcher, tests/fixtures, documentation, lockfile and release metadata, while excluding local runs, binaries, build outputs and archives.

Removed pre-release helpers (`pileup::apply_record`, `consensus::refresh_variants`, `repeat::summarize` and `RepeatSummary`) had no production callers; mapping uses `FragmentPileupAccumulator`. No compatibility executable for `mtasm-rust` is installed. Update commands and Rust imports to `mitovarr`. Existing `mitovar-rust.*` JSON schemas and gene panels remain compatible; new GFF3 source fields use `mitovarr`.

`docs/cli.md` is captured from the current executable's `--help` output. Refresh it along with method guides when CLI options or defaults change.

## Publication metadata

The project owner still needs to choose a license and supply the real GitHub repository URL. Add the matching license text and `license` / `repository` manifest fields after those decisions; no license or URL is assumed here. Confirm the release version/date and update `CHANGELOG.md` before a tagged release. Source cleanup and local verification do not publish the repository.

Once the release state is committed, a crates.io release can be rehearsed with `cargo publish --dry-run --locked`. GitHub publication, registry publication and release asset uploads are separate operations.

## Platform scope

Preparation uses macOS ARM64 and Rust 1.90.0, matching the declared minimum; edition remains 2021. Validate other targets before advertising their binaries. Source builds need a C compiler for `ring`, but the executable does not launch external bioinformatics programs.

Figure generation is optional; dependencies and editable sources are documented with [Figure 1](figure_1_workflow/README.md).
