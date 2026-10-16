# Publishing binary releases

The independent repository is `lrslab/mitovarr`; its root contains `Cargo.toml`
and `.github/workflows/`. GitHub builds the `mitovarr` executable from this
repository. The parent development directory is not part of the build.

## Platforms and archive contents

| Target | Build and test runner | Runtime baseline |
| --- | --- | --- |
| `aarch64-apple-darwin` | `macos-14` (Apple Silicon) | macOS 11+, ARM64 |
| `x86_64-unknown-linux-musl` | `ubuntu-22.04` with `musl-tools` | Linux x86_64; static, no glibc dependency |
| `aarch64-unknown-linux-gnu` | `ubuntu-22.04-arm` | Linux ARM64 with glibc 2.35+ |

Each archive is named `mitovarr-v<version>-<target>.tar.gz` and contains one
directory with the same name. It includes the executable, README, changelog,
portable documentation and figures, `scripts/fetch_caendr.sh`, and
`scripts/smoke_release.sh`. The built-in gene panels are embedded in the binary.
Source-code links in the packaged guides point to the matching GitHub tag;
guides and figures remain local. Research outputs and development runs are
excluded. Any license files explicitly included by Cargo are copied as well.

Users need no Cargo or Python to run `mitovarr`. The optional smoke script
requires Bash and Python 3. The macOS binary has no Developer ID signing or
notarization step. Download and installation commands are in the
[package README](../README.md#build-and-install).

The workflow follows `trackcluster-rs`'s version-tag trigger, target matrix,
reusable CI gate and archive naming. Linux ARM64 uses a native GitHub runner,
so every target can execute its own archive smoke test. Its glibc baseline is
2.35, rather than the cross-built 2.31 baseline used by `trackcluster-rs`.
Runner labels follow the [GitHub runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).

## Prepare a version

1. Update `package.version` in `Cargo.toml` and refresh the package entry in
   `Cargo.lock` with `cargo check`. Keep dependency changes intentional.
2. Give the corresponding `CHANGELOG.md` section a dated heading, for example
   `## 0.1.0 - YYYY-MM-DD`, replacing the placeholder with the actual release
   date. `Unreleased`, missing sections and empty notes are rejected.
3. Update installation and CLI documentation when behavior changes. New
   portable documents and release scripts must appear in the explicit Cargo
   `include` list.
4. Commit the intended source and documentation changes. Do not include
   `target/`, `dist/`, local data or experiment records.

Rust is pinned to 1.90.0 in both workflows, matching the declared minimum
version. Actions are pinned to commits. When changing the Rust version,
update the manifest and both workflows together.

## Validate and rehearse packaging

From the repository root:

```bash
cargo fmt --check
cargo clippy --locked --all-targets --all-features -- -D warnings
RUST_BACKTRACE=1 cargo test --locked --release --all-features -- --test-threads=1 --nocapture
cargo package --locked
```

`.github/workflows/ci.yml` runs formatting, Clippy and the release-profile test
suite on all three native runners, and verifies the standalone Cargo package
on Linux x86_64. Test cases run sequentially with live diagnostics and panic
backtraces; scoring-worker tests retain their internal concurrency. This avoids
overlapping assembly fixtures and identifies the active test if a process aborts.
On macOS failures, CI also waits briefly for system crash reports, prints any
matching test reports in the job log and saves their full contents as the
`macos-test-crash-reports` artifact. Absence of a crash report does not establish
the abort's cause.
CI runs on branch pushes, pull requests and manual dispatch, and is called again
by the release workflow for the tagged commit.

For Linux x86_64, use `musl-gcc` to compile C dependencies and `gcc` as the
Rust linker (`CC_x86_64_unknown_linux_musl=musl-gcc` and
`CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER=gcc`). Rust 1.90's musl target
links its bundled static runtime by default. Ubuntu 22.04's `musl-gcc` wrapper
does not correctly handle the `-static-pie` link; do not use it as the final
linker. Keep the release's ELF linkage check enabled.

To rehearse the archive on an Apple Silicon Mac:

```bash
target=aarch64-apple-darwin
version=$(cargo metadata --locked --no-deps --format-version 1 | python3 -c 'import json,sys; print(json.load(sys.stdin)["packages"][0]["version"])')
MACOSX_DEPLOYMENT_TARGET=11.0 cargo build --locked --release --bin mitovarr --target "$target"
python3 scripts/package_release.py "$target"
unpack=$(mktemp -d)
tar -xzf "dist/mitovarr-v$version-$target.tar.gz" -C "$unpack"
bash "$unpack/mitovarr-v$version-$target/scripts/smoke_release.sh" \
  "$unpack/mitovarr-v$version-$target" "$version" "$(mktemp -d)"
```

Packaging requires Python 3 and uses Cargo's explicit source-package file
list. For a development rehearsal with uncommitted changes, pass
`--allow-dirty` to the packaging script and `cargo package`; release CI does
not use that exception. The script checks documentation links before writing
the archive. Smoke checks run the unpacked executable from a separate working
directory, exercise command help and the extraction wrapper, and annotate a
small synthetic gzipped FASTA with all three embedded species panels offline.
The CI test suite covers the assembly and mapping algorithms separately.

## Publish the tag

Push the release commit to `main` and wait for its CI run to pass. From that
same clean commit, with `version` set to the Cargo version as above:

```bash
test "$(git branch --show-current)" = main
test -z "$(git status --porcelain)"
test -n "$version"
git tag -a "v$version" -m "Release v$version"
git push origin "v$version"
```

`.github/workflows/release.yml` requires exactly `v<Cargo package version>`.
It validates the dated changelog, waits for reusable CI, then builds, packages
and smoke-tests all three targets. Linux x86_64 additionally checks that the
executable has neither a dynamic loader nor shared-library dependencies. Only
after all targets succeed does the publish job create the GitHub Release and
attach the three archives. Release
notes contain the version's changelog followed by GitHub-generated notes.
Versions such as `0.2.0-rc.1` are marked as prereleases.

The default `GITHUB_TOKEN` is sufficient: the publishing job alone requests
`contents: write`; other jobs have read-only repository access. No separate
release token is needed. A manual release dispatch must select an existing
version tag (for example `gh workflow run release.yml --ref v0.1.0`); selecting
a branch fails validation. Failed runs can also be rerun from GitHub Actions.

Confirm all three archives are attached before announcing the release. Never
move a published tag to different source; use a new version for corrections.
This workflow publishes GitHub binary releases only. crates.io publication is
a separate maintainer action.
