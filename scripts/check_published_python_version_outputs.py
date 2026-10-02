#!/usr/bin/env python3
"""Fail-closed rowwise comparison of two published-wheel operation exports.

The export format is produced by bench_python_op_matrix_vs_rdkit.py
--outputs-only --output-values-jsonl. Only HBA and its named bundle may
change across v1.0.29 -> v1.0.30. Every other value *and error* is compared.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

EXPECTED_CHANGES = {"hba", "lipinski_bundle(mw,logp,hbd,hba)"}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compare(old_path: Path, new_path: Path) -> dict:
    counts: Counter[str] = Counter()
    errors: Counter[str] = Counter()
    changes: dict[str, list[int]] = {}
    opener_old = gzip.open if old_path.suffix == ".gz" else Path.open
    opener_new = gzip.open if new_path.suffix == ".gz" else Path.open
    with opener_old(old_path, "rt", encoding="utf-8") as old_file, opener_new(new_path, "rt", encoding="utf-8") as new_file:
        for line_no, (old_line, new_line) in enumerate(zip(old_file, new_file, strict=True), 1):
            old = json.loads(old_line)
            new = json.loads(new_line)
            if old["op"] != new["op"] or old["op"] in counts:
                raise ValueError(f"operation name mismatch or duplicate on line {line_no}")
            name = old["op"]
            if len(old["rows"]) != len(new["rows"]):
                raise ValueError(f"{name}: row count differs")
            changed = []
            for index, (before, after) in enumerate(zip(old["rows"], new["rows"], strict=True)):
                if before.get("input_index") != index or after.get("input_index") != index:
                    raise ValueError(f"{name}: noncontiguous input index at {index}")
                if ("error" in before) == ("value" in before) or ("error" in after) == ("value" in after):
                    raise ValueError(f"{name}: row {index} must contain exactly one result or error")
                if "error" in before or "error" in after:
                    errors[name] += 1
                if before != after:
                    changed.append(index)
            counts[name] = len(old["rows"])
            if changed:
                changes[name] = changed
    if len(counts) != 63:
        raise ValueError(f"expected 63 operations, found {len(counts)}")
    unexpected = sorted(set(changes) - EXPECTED_CHANGES)
    if unexpected:
        raise ValueError(f"non-HBA changes: {unexpected}")
    return {"schema": "published-python-version-output-diff/v1",
            "artifacts": {"v1.0.29_outputs_sha256": file_sha256(old_path),
                          "v1.0.30_outputs_sha256": file_sha256(new_path)},
            "operation_count": len(counts), "output_count": sum(counts.values()),
            "rows_by_operation": dict(counts), "error_rows_by_operation": dict(errors),
            "changed_rows_by_operation": {name: {"count": len(indices), "indices": indices}
                                          for name, indices in changes.items()},
            "unchanged_non_hba_operations": len(counts) - len(changes),
            "unexpected_changed_operations": unexpected}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old", type=Path, required=True)
    parser.add_argument("--new", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.old, args.new)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"checked {result['operation_count']} operations and {result['output_count']} rows; "
          f"changed: {[(op, info['count']) for op, info in result['changed_rows_by_operation'].items()]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
