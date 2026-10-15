#!/usr/bin/env python3
"""Check canonical product spelling in public text and package metadata."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_SUFFIXES = {
    ".cff",
    ".html",
    ".htm",
    ".md",
    ".rst",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
}
PUBLIC_FILENAMES = {"package.json", "pyproject.toml"}
HISTORICAL_PREFIXES = ("benchmarks/", "docs/archive/", "validation/results/")
MISSPELLING = re.compile(r"\b" + "s" + r"chematic(?:-draw)?\b", re.IGNORECASE)


def tracked_public_files(root: Path = ROOT) -> list[Path]:
    completed = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=root,
        check=True,
        capture_output=True,
    )
    paths: list[Path] = []
    for item in completed.stdout.split(b"\0"):
        if not item:
            continue
        relative = Path(item.decode())
        relative_text = relative.as_posix()
        if relative_text.startswith(HISTORICAL_PREFIXES):
            continue
        if relative_text.startswith(".github/workflows/"):
            # Workflow command lines retain versioned internal CLI/schema names.
            continue
        if relative.suffix.lower() in PUBLIC_SUFFIXES or relative.name in PUBLIC_FILENAMES:
            paths.append(root / relative)
    return paths


def scan_paths(root: Path, paths: list[Path]) -> list[str]:
    errors: list[str] = []
    for path in paths:
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        for line_number, line in enumerate(text.splitlines(), start=1):
            if MISSPELLING.search(line):
                errors.append(
                    f"{path.relative_to(root)}:{line_number}: "
                    "unexpected leading-s product spelling"
                )
    return errors


def main() -> int:
    errors = scan_paths(ROOT, tracked_public_files())
    if errors:
        print("Product naming check failed:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("Product naming check OK: chematic and chematic-draw are canonical")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
