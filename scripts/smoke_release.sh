#!/usr/bin/env bash
# Check an unpacked binary release offline, without Cargo or source files.
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 <release-dir> <version> <work-dir>" >&2
  exit 2
fi

release_dir=$(cd "$1" && pwd)
version=$2
mkdir -p "$3"
work_dir=$(cd "$3" && pwd)
binary="$release_dir/mitovarr"
cd "$work_dir"

test -x "$binary"
test "$("$binary" --version)" = "mitovarr $version"
"$binary" --help > help.txt
for command in fetch-caendr bait assemble-baited assemble-seeded assemble-mapping annotate; do
  "$binary" "$command" --help > "$command-help.txt"
done
sh "$release_dir/scripts/fetch_caendr.sh" --help > wrapper-help.txt
cmp fetch-caendr-help.txt wrapper-help.txt

# A short synthetic sequence exercises gzip input and all three embedded panels.
python3 - <<'PY'
import gzip
with gzip.open("input.fasta.gz", "wt") as output:
    output.write(">release_smoke\n" + "ACGT" * 30 + "\n")
PY
for species in cel cbr ctro; do
  "$binary" annotate --fasta input.fasta.gz --species "$species" --out "$species"
done
python3 - "$version" <<'PY'
import json
from pathlib import Path
import sys

for species in ("cel", "cbr", "ctro"):
    root = Path(species)
    report = json.loads((root / "report.json").read_text())
    assert report["software_version"] == sys.argv[1]
    assert report["status"] == "completed"
    assert report["assembly_status"] == "not_modified"
    assert report["closure_validation"] == "not_performed"
    assert report["records"][0]["length"] == 120
    assert (root / "annotated.fasta").stat().st_size > 0
    assert (root / "gene_quality.tsv").stat().st_size > 0
print("Release smoke checks passed: CLI, wrapper, gzip input and embedded gene panels.")
PY
