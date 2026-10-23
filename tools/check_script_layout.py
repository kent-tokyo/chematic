#!/usr/bin/env python3
"""Keep new maintenance scripts out of the legacy flat scripts directory."""

from __future__ import annotations

import subprocess
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "tools/script_root_baseline.list"
ALLOWED_DIRECTORIES = frozenset(
    {
        "benchmarks",
        "checks",
        "diagnostics",
        "fixtures",
        "generators",
        "lib",
        "release",
        "tests",
        "validation",
    }
)


def validate_layout(tracked: list[str], baseline: list[str]) -> list[str]:
    """Return layout errors without needing a Git checkout for negative controls."""
    errors: list[str] = []
    if baseline != sorted(set(baseline)):
        errors.append("tools/script_root_baseline.list must be sorted and unique")
    if any(PurePosixPath(path).parent != PurePosixPath("scripts") for path in baseline):
        errors.append("baseline entries must be direct children of scripts/")

    root_files = {path for path in tracked if PurePosixPath(path).parent == PurePosixPath("scripts")}
    new_root = root_files - set(baseline)
    stale = set(baseline) - root_files
    if new_root:
        errors.append("new files in scripts/ root: " + ", ".join(sorted(new_root)))
    if stale:
        errors.append("remove migrated files from baseline: " + ", ".join(sorted(stale)))

    unexpected = {
        path for path in tracked
        if len(PurePosixPath(path).parts) > 2
        and PurePosixPath(path).parts[1] not in ALLOWED_DIRECTORIES
    }
    if unexpected:
        errors.append("unclassified scripts/ subdirectories: " + ", ".join(sorted(unexpected)))
    return errors


def existing_files(root: Path, tracked: list[str]) -> list[str]:
    """Ignore paths deleted in the working tree but still present in Git's index."""
    return [path for path in tracked if (root / path).is_file()]


def main() -> int:
    result = subprocess.run(
        ["git", "ls-files", "-z", "--", "scripts"],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    tracked = existing_files(
        ROOT,
        [path for path in result.stdout.decode().split("\0") if path],
    )
    baseline = [
        line.strip() for line in BASELINE.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    errors = validate_layout(tracked, baseline)
    if errors:
        for error in errors:
            print(error)
        return 1
    print(f"script layout OK: {len(baseline)} legacy root files, {len(tracked)} tracked scripts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
