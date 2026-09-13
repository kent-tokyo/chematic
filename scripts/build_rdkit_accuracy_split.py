#!/usr/bin/env python3
"""Build a reproducible, explicitly non-sealed descriptor evaluation split.

This is preparation infrastructure only. Existing repository corpora are
already exposed to development, so the generated holdout is labelled
``exposed_holdout`` and can never be mistaken for the planned sealed set.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUTS = [
    ROOT / "scripts/descriptor_census_corpus.smi",
    ROOT / "scripts/chembl_accuracy_corpus_4999.smi",
]
DEFAULT_OUT = ROOT / "validation/rdkit_accuracy_split_v2"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_unique(paths: list[Path]) -> tuple[list[str], list[dict[str, object]]]:
    rows: list[str] = []
    seen: set[str] = set()
    sources: list[dict[str, object]] = []
    for path in paths:
        source_rows = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
        before = len(rows)
        for smiles in source_rows:
            if smiles not in seen:
                seen.add(smiles)
                rows.append(smiles)
        sources.append({
            "path": str(path.relative_to(ROOT)),
            "sha256": digest(path),
            "input_rows": len(source_rows),
            "new_unique_rows": len(rows) - before,
        })
    return rows, sources


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--development-rows", type=int, default=2000)
    parser.add_argument("inputs", type=Path, nargs="*", default=DEFAULT_INPUTS)
    args = parser.parse_args()
    if args.development_rows <= 0:
        parser.error("--development-rows must be positive")
    paths = [path.resolve() for path in args.inputs]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        parser.error("missing input: " + ", ".join(missing))
    rows, sources = read_unique(paths)
    if args.development_rows >= len(rows):
        parser.error("development split must leave at least one holdout row")
    development = rows[: args.development_rows]
    exposed_holdout = rows[args.development_rows :]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    development_path = output / "development.smi"
    holdout_path = output / "exposed_holdout.smi"
    development_path.write_text("\n".join(development) + "\n", encoding="utf-8")
    holdout_path.write_text("\n".join(exposed_holdout) + "\n", encoding="utf-8")
    metadata = {
        "schema_version": 1,
        "protocol": "rdkit-accuracy-split-v2",
        "status": "exposed_holdout_not_sealed",
        "sealed": False,
        "reason": "Inputs are existing development/repository corpora; independent acquisition and freeze are still required.",
        "deduplication": "exact SMILES, first occurrence wins",
        "sources": sources,
        "total_unique_rows": len(rows),
        "planned_sealed_rows": 8000,
        "sealed_shortfall_rows": max(0, 8000 - len(exposed_holdout)),
        "splits": {
            "development": {"path": str(development_path.relative_to(ROOT)), "rows": len(development), "sha256": digest(development_path)},
            "exposed_holdout": {"path": str(holdout_path.relative_to(ROOT)), "rows": len(exposed_holdout), "sha256": digest(holdout_path)},
        },
    }
    (output / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": metadata["status"], "development_rows": len(development), "exposed_holdout_rows": len(exposed_holdout), "sealed_shortfall_rows": metadata["sealed_shortfall_rows"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
