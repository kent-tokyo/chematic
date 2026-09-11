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

# The checked-in twelve-case corpus is intentionally stable. These additional
# cases exercise Extended XYZ metadata branches that cannot be reached by the
# old generic malformed shapes, while keeping the base-corpus hash and its
# historical comparisons unchanged.
SUPPLEMENTAL_CASES = {
    "extxyz": [
        "1\n",
        "1\n=bad\nC 0 0 0\n",
        "1\nfoo=\"bad\nC 0 0 0\n",
        "1\nfoo=1 foo=2\nC 0 0 0\n",
        "1\nLattice=\"1 2\"\nC 0 0 0\n",
        "1\nProperties=species:S:1:pos:R\nC 0 0 0\n",
    ],
}


def run(binary: list[str], fmt: str, path: Path, line_limited: bool = True) -> dict[str, object]:
    options = ("--max-line-bytes", "1") if line_limited and fmt in {"cml", "cdxml", "pdb"} else ()
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
    mmcif = (ROOT / "benchmarks" / "fixtures" / "minimal.mmcif").read_text(encoding="utf-8")
    mol2 = (ROOT / "benchmarks" / "fixtures" / "ethanol.mol2").read_text(encoding="utf-8")
    cml = (ROOT / "benchmarks" / "fixtures" / "ethanol.cml").read_text(encoding="utf-8")
    cdxml = (ROOT / "benchmarks" / "fixtures" / "ethanol.cdxml").read_text(encoding="utf-8")
    supplemental_cases = {
        **SUPPLEMENTAL_CASES,
        "mmcif": [
            mmcif.replace("0.0 0.0 0.0", "nope 0.0 0.0", 1),
            mmcif.replace("ATOM 1 C", "ATOM nope C", 1),
            mmcif.replace("ATOM 1 C C1", "ATOM 1 Xx C1", 1),
            mmcif.replace("_atom_site.Cartn_z\n", "", 1),
        ],
        "mol2": [
            mol2.replace("     1     1     2    1", "     1     1", 1),
            mol2.replace("C.3", "Xx", 1),
        ],
        "cml": [
            cml.replace('elementType="C"', 'elementType="Xx"', 1),
            cml.replace('atomRefs2="a1 a2"', 'atomRefs2="a1 z9"', 1),
            cml.replace('atomRefs2="a1 a2"', 'atomRefs2="a1"', 1),
            cml.replace('order="1"', 'order="wat"', 1),
            cml.replace('x2="0.0"', 'x2="nan"', 1),
        ],
        "cdxml": [
            cdxml.replace('Element="6"', 'Element="999"', 1),
            cdxml.replace('E="3"', 'E="9"', 1),
            cdxml.replace(' B="1" E="2"', ' B="1"', 1),
        ],
    }

    rows: dict[str, dict[str, object]] = {}
    errors: list[str] = []
    total = 0
    with tempfile.TemporaryDirectory(prefix="chematic-streaming-taxonomy-") as directory:
        temp = Path(directory)
        for fmt in FORMATS:
            counts: Counter[str] = Counter()
            cases = [*corpus["cases"][fmt], *supplemental_cases.get(fmt, [])]
            for index, content in enumerate(cases):
                path = temp / f"case-{index}.{fmt}"
                path.write_text(content, encoding="utf-8")
                base_count = len(corpus["cases"][fmt])
                line_limited = not (fmt in {"cml", "cdxml", "pdb"} and index >= base_count)
                result = run(args.binary, fmt, path, line_limited=line_limited)
                total += 1
                if result.get("records") != 0 or result.get("failures") != 1:
                    errors.append(f"{fmt}[{index}] aggregate failure mismatch: {result}")
                kinds = result.get("failure_kinds")
                if not isinstance(kinds, dict) or sum(kinds.values()) != 1:
                    errors.append(f"{fmt}[{index}] failure_kinds must contain exactly one failure: {result}")
                else:
                    for kind, count in kinds.items():
                        counts[str(kind)] += int(count)
            rows[fmt] = {
                "cases": len(cases),
                "base_cases": len(corpus["cases"][fmt]),
                "supplemental_cases": len(supplemental_cases.get(fmt, [])),
                "failure_kinds": dict(sorted(counts.items())),
            }

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
