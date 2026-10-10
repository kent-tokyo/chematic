#!/usr/bin/env python3
"""Check CML/CDXML/mmCIF/PDB record accounting against Open Babel's CLI boundary."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

from benchmark_version import workspace_version

ROOT = Path(__file__).resolve().parents[1]


def executable_version(executable: str) -> str:
    completed = subprocess.run(
        [executable, "-V"],
        check=True,
        text=True,
        capture_output=True,
    )
    output = (completed.stdout or completed.stderr).strip()
    return output.splitlines()[0] if output else "unknown"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("cml", "cdxml", "mmcif", "pdb"), default="cml")
    parser.add_argument("--cml", type=Path, default=Path("benchmarks/fixtures/ethanol.cml"))
    parser.add_argument("--cdxml", type=Path, default=Path("benchmarks/fixtures/ethanol.cdxml"))
    parser.add_argument("--mmcif", type=Path, default=Path("benchmarks/fixtures/minimal.mmcif"))
    parser.add_argument("--pdb", type=Path, default=Path("benchmarks/fixtures/minimal.pdb"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--openbabel", default="obabel")
    parser.add_argument("--binary", nargs="+", default=["cargo", "run", "-p", "chematic-mol", "--example", "streaming_benchmark", "--offline", "--"])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = {"cml": args.cml, "cdxml": args.cdxml, "mmcif": args.mmcif, "pdb": args.pdb}[args.format]
    path = path if path.is_absolute() else ROOT / path
    payload = path.read_bytes()
    expected = args.repeats
    chematic = json.loads(subprocess.run([*args.binary, "--format", args.format, "--path", str(path), "--repeats", str(args.repeats)], cwd=ROOT, check=True, text=True, capture_output=True).stdout)
    started = subprocess.run
    converted = 0
    failures = 0
    pattern = re.compile(r"^(\d+) molecules? converted$")
    for _ in range(args.repeats):
        result = started([args.openbabel, f"-i{args.format}", str(path), "-osmi", "-O", "/dev/null"], check=False, text=True, capture_output=True)
        match = next((pattern.match(line.strip()) for line in result.stderr.splitlines()), None)
        if match:
            converted += int(match.group(1))
        if result.returncode != 0 or not match:
            failures += 1
    errors = []
    for name, row in (("chematic", chematic), ("openbabel", {"records": converted, "failures": failures, "input_bytes": len(payload) * args.repeats})):
        if row.get("records") != expected or row.get("failures") != 0:
            errors.append(f"{name} count mismatch: {row}")
        if row.get("input_bytes") != len(payload) * args.repeats:
            errors.append(f"{name} input byte mismatch: {row.get('input_bytes')}")
    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "status": "local-verified" if not errors else "failed",
        "gate": "openbabel_record_accounting",
        "format": args.format,
        "fixture": {
            "path": str(path.relative_to(ROOT)),
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        },
        "repeats": args.repeats,
        "expected_records": expected,
        "rows": {
            "chematic": chematic,
            "openbabel": {
                "records": converted,
                "failures": failures,
                "input_bytes": len(payload) * args.repeats,
                "comparison_boundary": "Open Babel CLI process per repetition; startup and conversion included",
            },
        },
        "comparison_boundary": {
            "chematic": f"Rust {args.format.upper()} materialized one-shot parser",
            "openbabel": f"Open Babel {args.format.upper()} CLI conversion per repetition",
        },
        "tool_versions": {"openbabel": executable_version(args.openbabel)},
        "claim_boundary": "record and failure accounting only; not semantic round-trip or speed parity",
    }
    if errors:
        print(f"streaming {args.format.upper()} contract failures:", *errors, sep="\n")
        return 1
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else ROOT / args.output
        target.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
