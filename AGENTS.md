# mitovarr

This directory implements the user's 2026-09-09 explicit request for a separate
pure-Rust fallback. That request supersedes the parent non-goal against adding a
second in-process assembler only within this directory.

- Keep this crate independently buildable. No imports or path dependencies on the parent project.
- Package, library and executable are named `mitovarr`. Existing `mitovar-rust.*` JSON schema identifiers remain unchanged for compatibility with saved reports and gene panels.
- This directory is the independent Git root. Track `Cargo.lock`, source, tests and portable documentation; keep generated runs, executables and release archives out of Git.
- Keep experiment records, analysis outputs and local binaries outside this repository in the sibling `mitovar_rust_local/` directory. Builds and tests must not depend on them.
- Keep Cargo package contents explicit in `Cargo.toml` and update `CHANGELOG.md` for release-facing changes. Do not invent a license or repository URL.
- Current usage and methods are indexed in `docs/README.md`. Do not include experiment logs, sample result tables or historical research figures in the repository or source package. Keep necessary regression tests, portable fixtures and runtime gene models.
- Production Rust code must never spawn a process or call an external assembler, mapper, recruiter, decompressor or shell.
- Optional documentation figure scripts are separate from the Rust build and runtime.
- Keep reference-blind FASTQ, seed recruitment and reference-guided indexed BAM inputs separate.
- `fetch-caendr` produces `mt-bam-extracted` FASTQ and indexed BAM. These reads are not original whole-sample FASTQ; no second baiting step is needed.
- `assemble-baited` assembles those paired FASTQ directly with the native graph engine and reports `mt-bam-extracted` input and validation scope. No seed or reference enters assembly.
- A seed may only select complete original pairs. Final validation uses all original paired FASTQ.
- Incomplete recruitment, ambiguous graph structure or insufficient read support must not publish a final FASTA.
- Simplify local sequencing-error paths and use sample-read bridges before judging graph continuity. A fragmented graph alone does not establish that the sample cannot be assembled.
- Keep resource limits explicit; never silently truncate reads or graph branches and publish.
- User explicitly rejected hashes, repeated input-integrity scans and additional audit designs. Keep the implementation minimal.
- Update README, tests and todo.md when changing a CLI or output contract.
- Run cargo fmt --check, cargo test and cargo clippy --all-targets --all-features.
