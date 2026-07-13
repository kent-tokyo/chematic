#!/usr/bin/env python3
"""Update version strings across docs when releasing a new chematic version.

Usage:
    python scripts/bump_version.py               # auto-detect old from git tag
    python scripts/bump_version.py --old 0.4.23  # explicit old version
"""

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent


def workspace_version() -> str:
    text = (ROOT / "Cargo.toml").read_text()
    m = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if not m:
        sys.exit("Could not find workspace version in Cargo.toml")
    return m.group(1)


def previous_tag_version() -> str:
    try:
        out = subprocess.check_output(
            ["git", "describe", "--tags", "--abbrev=0", "HEAD^"],
            cwd=ROOT, text=True, stderr=subprocess.DEVNULL,
        ).strip().lstrip("v")
        return out
    except subprocess.CalledProcessError:
        sys.exit(
            "Could not determine previous version from git tags.\n"
            "Pass --old VERSION explicitly."
        )


def bump_file(path: Path, old: str, new: str, patterns: list[re.Pattern]) -> bool:
    text = path.read_text()
    changed = False
    for pat in patterns:
        new_text = pat.sub(lambda m: m.group(0).replace(old, new), text)
        if new_text != text:
            text = new_text
            changed = True
    if changed:
        path.write_text(text)
    return changed


def main() -> None:
    args = sys.argv[1:]
    old = None
    if "--old" in args:
        idx = args.index("--old")
        old = args[idx + 1]
    dry_run = "--dry-run" in args

    new = workspace_version()
    if old is None:
        old = previous_tag_version()

    if old == new:
        print(f"Already at {new}, nothing to do.")
        return

    print(f"Bumping {old} → {new}")

    # (file, list of patterns that contain the version to replace)
    targets: list[tuple[Path, list[re.Pattern]]] = [
        (ROOT / "README.md", [
            re.compile(r"chematic v" + re.escape(old)),
            re.compile(r"v" + re.escape(old) + r" vs RDKit"),
            # Repository Structure tree comment: "workspace root (vX.Y.Z)"
            re.compile(r"workspace root \(v" + re.escape(old) + r"\)"),
            # BibTeX citation block: "version   = {X.Y.Z},"
            re.compile(r"version\s*=\s*\{" + re.escape(old) + r"\}"),
        ]),
        (ROOT / "README_ja.md", [
            re.compile(r"chematic v" + re.escape(old)),
            re.compile(r"v" + re.escape(old) + r" vs RDKit"),
        ]),
        (ROOT / "README_zh.md", [
            re.compile(r"chematic v" + re.escape(old)),
            re.compile(r"v" + re.escape(old) + r" vs RDKit"),
        ]),
        (ROOT / "docs" / "validation.md", [
            re.compile(r"chematic v" + re.escape(old)),
        ]),
        (ROOT / "docs" / "benchmark.md", [
            re.compile(r"chematic v" + re.escape(old)),
        ]),
        (ROOT / "CITATION.cff", [
            re.compile(r"^version: " + re.escape(old), re.MULTILINE),
        ]),
        (ROOT / "SECURITY.md", [
            # current release row
            re.compile(r"\| v" + re.escape(old) + r" \| Yes \|[^\n]*"),
            # "Latest release (vX.Y.Z)" prose line
            re.compile(r"Latest release \(v" + re.escape(old) + r"\)"),
        ]),
        (ROOT / "demo" / "pkg" / "package.json", [
            re.compile(r'"version":\s*"' + re.escape(old) + r'"'),
        ]),
    ]

    # Intra-workspace path-dependency version pins — key order varies between
    # crates ("path, version" in most, "version, path" in chematic-inchi), so
    # match both:
    #   chematic-core = { path = "../chematic-core", version = "0.4.19" }
    #   chematic-core = { version = "0.4.19", path = "../chematic-core" }
    for cargo_toml in sorted((ROOT / "crates").glob("*/Cargo.toml")):
        targets.append((cargo_toml, [
            re.compile(r'path\s*=\s*"[^"]+",\s*version\s*=\s*"' + re.escape(old) + r'"'),
            re.compile(r'version\s*=\s*"' + re.escape(old) + r'",\s*path\s*=\s*"[^"]+"'),
        ]))

    updated = []
    for path, patterns in targets:
        if not path.exists():
            continue
        if dry_run:
            text = path.read_text()
            if any(p.search(text) for p in patterns):
                print(f"  [dry-run] would update {path.relative_to(ROOT)}")
            continue
        if bump_file(path, old, new, patterns):
            updated.append(path.relative_to(ROOT))

    if dry_run:
        return

    if updated:
        print("Updated:")
        for p in updated:
            print(f"  {p}")
        print(f"\nNext: bash scripts/check.sh  →  git add -p  →  git commit")
    else:
        print("No occurrences of old version found in target files.")


if __name__ == "__main__":
    main()
