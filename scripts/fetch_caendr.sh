#!/bin/sh
# Run this package's native CaenDR extractor; forward its CLI arguments unchanged.
set -eu

package_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
if [ -x "$package_dir/mitovarr" ]; then
    exec "$package_dir/mitovarr" fetch-caendr "$@"
fi
exec cargo run --release --locked --manifest-path "$package_dir/Cargo.toml" \
    --bin mitovarr -- fetch-caendr "$@"
