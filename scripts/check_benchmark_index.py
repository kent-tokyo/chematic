#!/usr/bin/env python3
"""Check that the benchmark record index links every checked-in record."""

from __future__ import annotations

import re
import sys
from pathlib import Path
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks"
INDEX = BENCHMARKS / "README.md"
LINK_PATTERN = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
RECORD_SUFFIXES = {".json", ".md"}


def benchmark_records() -> set[Path]:
    return {
        path.relative_to(BENCHMARKS)
        for path in BENCHMARKS.iterdir()
        if path.is_file()
        and path.name != INDEX.name
        and path.suffix in RECORD_SUFFIXES
    }


def main() -> int:
    errors: list[str] = []
    try:
        text = INDEX.read_text(encoding="utf-8")
    except OSError as exc:
        print(f"benchmark index check failed: {exc}", file=sys.stderr)
        return 1

    linked_records: set[Path] = set()
    for raw_target in LINK_PATTERN.findall(text):
        target = raw_target.strip().split(None, 1)[0]
        if target.startswith(("https://", "http://", "mailto:")):
            continue
        target_without_fragment = target.split("#", 1)[0]
        if not target_without_fragment:
            continue
        target_path = Path(unquote(target_without_fragment))
        resolved = (INDEX.parent / target_path).resolve()
        try:
            relative = resolved.relative_to(BENCHMARKS.resolve())
        except ValueError:
            relative = None
        if not resolved.is_file():
            errors.append(f"broken local link: {target}")
        elif relative is not None and relative.suffix in RECORD_SUFFIXES:
            linked_records.add(relative)

    missing = sorted(benchmark_records() - linked_records)
    errors.extend(f"unindexed benchmark record: {path}" for path in missing)

    if errors:
        print("Benchmark index consistency failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1

    print(f"Benchmark index consistency OK: {len(linked_records)} records indexed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
