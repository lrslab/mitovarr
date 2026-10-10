#!/usr/bin/env python3
"""Package an already-built mitovarr binary with the portable Cargo documentation."""

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile
from urllib.parse import quote, unquote, urlsplit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=[
        "aarch64-apple-darwin", "x86_64-unknown-linux-musl", "aarch64-unknown-linux-gnu",
    ])
    parser.add_argument("--allow-dirty", action="store_true", help="include local changes for a packaging rehearsal")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    metadata = json.loads(subprocess.check_output(
        ["cargo", "metadata", "--locked", "--no-deps", "--format-version", "1"], cwd=root,
    ))
    package = next(p for p in metadata["packages"] if p["name"] == "mitovarr")
    tag = f"v{package['version']}"
    archive_name = f"mitovarr-{tag}-{args.target}"
    target_dir = Path(metadata["target_directory"])
    binary = target_dir / args.target / "release/mitovarr"
    if not binary.is_file():
        parser.error(f"build the release binary first: {binary}")
    command = ["cargo", "package", "--locked", "--list"]
    if args.allow_dirty:
        command.append("--allow-dirty")
    source_files = subprocess.check_output(command, cwd=root, text=True).splitlines()
    assets = [p for p in source_files if p.startswith(("docs/", "LICENSE")) or p in {
        "README.md", "CHANGELOG.md", "scripts/fetch_caendr.sh", "scripts/smoke_release.sh",
    }]
    output_dir = root / "dist"
    output_dir.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="release-stage-", dir=target_dir) as temporary:
        stage = Path(temporary) / archive_name
        stage.mkdir()
        shutil.copy2(binary, stage / "mitovarr")
        for relative in assets:
            destination = stage / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root / relative, destination)

        # Keep guides and figures local; link source-code references to this release tag.
        for document in stage.rglob("*.md"):
            def rewrite_link(match):
                target = match.group(2)
                url = urlsplit(target)
                if url.scheme or url.netloc or not url.path:
                    return match.group(0)
                resolved = (document.parent / unquote(url.path)).resolve()
                relative = resolved.relative_to(stage).as_posix()
                if resolved.exists():
                    return match.group(0)
                if relative in source_files and (root / relative).is_file():
                    source_url = f"{package['repository']}/blob/{quote(tag, safe='')}/{quote(relative)}"
                    if url.fragment:
                        source_url += f"#{url.fragment}"
                    return f"{match.group(1)}{source_url}{match.group(3)}"
                raise ValueError(f"{document.relative_to(stage)}: missing documentation target {target}")

            text = document.read_text(encoding="utf-8")
            text = re.sub(r"(!?\[[^\]]*\]\()([^)]+)(\))", rewrite_link, text)
            document.write_text(text, encoding="utf-8")

        for required in ("mitovarr", "README.md", "CHANGELOG.md", "docs/README.md",
                         "docs/publishing.md", "scripts/fetch_caendr.sh", "scripts/smoke_release.sh"):
            if not (stage / required).is_file() or (stage / required).stat().st_size == 0:
                raise ValueError(f"missing release asset: {required}")
        archive = output_dir / f"{archive_name}.tar.gz"
        with tarfile.open(archive, "w:gz") as output:
            output.add(stage, arcname=archive_name)
    print(archive)


if __name__ == "__main__":
    main()
