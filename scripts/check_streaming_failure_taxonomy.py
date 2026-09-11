#!/usr/bin/env python3
"""Record parser-error variant counts for the bounded streaming corpus.

This is diagnostic coverage evidence.  It deliberately keeps malformed-input
acceptance and parser-error taxonomy separate from any performance claim.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from benchmark_version import workspace_version

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "validation" / "streaming_format_safety_cases.json"
DEFAULT_OUTPUT = ROOT / "validation" / "results" / f"streaming-failure-taxonomy-v{workspace_version(ROOT)}.json"
FORMATS = ("sdf", "mol", "xyz", "extxyz", "v3000", "mol2", "cml", "cdxml", "mmcif", "pdb")


def run(binary: list[str], fmt: str, path: Path) -> dict[str, object]:
    options = ("--max-line-bytes", "1") if fmt in {"cml", "cdxml", "pdb"} else ()
    completed = subprocess.run(
        [*binary, "--format", fmt, "--path", str(path), "--repeats", "1", *options],
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    return json.loads(completed.stdout)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", nargs="+", default=["target/debug/examples/streaming_benchmark"])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    corpus = json.loads(CORPUS.read_text(encoding="utf-8"))
    if corpus.get("schema_version") != 1 or set(corpus.get("cases", {})) != set(FORMATS):
        raise SystemExit("streaming safety corpus schema or format set is invalid")

    rows: dict[str, dict[str, object]] = {}
    errors: list[str] = []
    total = 0
    with tempfile.TemporaryDirectory(prefix="chematic-streaming-taxonomy-") as directory:
        temp = Path(directory)
        for fmt in FORMATS:
            counts: Counter[str] = Counter()
            cases = corpus["cases"][fmt]
            for index, content in enumerate(cases):
                path = temp / f"case-{index}.{fmt}"
                path.write_text(content, encoding="utf-8")
                result = run(args.binary, fmt, path)
                total += 1
                if result.get("records") != 0 or result.get("failures") != 1:
                    errors.append(f"{fmt}[{index}] aggregate failure mismatch: {result}")
                kinds = result.get("failure_kinds")
                if not isinstance(kinds, dict) or sum(kinds.values()) != 1:
                    errors.append(f"{fmt}[{index}] failure_kinds must contain exactly one failure: {result}")
                else:
                    for kind, count in kinds.items():
                        counts[str(kind)] += int(count)
            rows[fmt] = {"cases": len(cases), "failure_kinds": dict(sorted(counts.items()))}

    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "status": "local-verified" if not errors else "failed",
        "gate": "streaming_failure_taxonomy",
        "corpus": {
            "path": str(CORPUS.relative_to(ROOT)),
            "sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
            "cases_per_format": 12,
            "formats": len(FORMATS),
            "total_cases": total,
        },
        "rows": rows,
        "boundary": "aggregate rejection remains the safety contract; failure_kinds is diagnostic variant evidence only",
        "errors": errors,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if errors:
        print("Streaming failure taxonomy failures:")
        print("\n".join(errors))
        return 1
    print(f"Streaming failure taxonomy OK: {total} cases, {len(FORMATS)} formats, report={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
