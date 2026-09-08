#!/usr/bin/env python3
"""Check gzip record accounting against Open Babel on the same payload.

Open Babel 3.x does not accept the gzip wrapper as an input format. The
contract therefore keeps the boundaries explicit: chematic reads the gzip
file, while Open Babel reads the byte-identical decompressed payload.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def version(executable: str) -> str:
    result = subprocess.run([executable, "-V"], check=True, text=True, capture_output=True)
    return (result.stdout or result.stderr).strip().splitlines()[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("sdf", "v3000", "mol2"), default="sdf")
    parser.add_argument("--sdf", type=Path, default=Path("benchmarks/fixtures/streaming.sdf"))
    parser.add_argument("--v3000", type=Path, default=Path("benchmarks/fixtures/ethanol.v3000"))
    parser.add_argument("--mol2", type=Path, default=Path("benchmarks/fixtures/ethanol.mol2"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--openbabel", default="obabel")
    parser.add_argument("--binary", nargs="+", default=[
        "cargo", "run", "-p", "chematic-mol", "--example", "streaming_benchmark", "--offline", "--",
    ])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")

    source = {"sdf": args.sdf, "v3000": args.v3000, "mol2": args.mol2}[args.format]
    source = source if source.is_absolute() else ROOT / source
    payload = source.read_bytes()
    compressed = gzip.compress(payload, mtime=0)
    openbabel_format = "sdf" if args.format == "sdf" else ("mol" if args.format == "v3000" else "mol2")
    expected_records = (2 if args.format == "sdf" else 1) * args.repeats

    with tempfile.TemporaryDirectory(prefix="chematic-gzip-openbabel-") as directory:
        temp = Path(directory)
        gzip_path = temp / f"fixture.{args.format}.gz"
        decompressed_path = temp / f"fixture.{args.format}"
        gzip_path.write_bytes(compressed)
        decompressed_path.write_bytes(payload)
        chematic = json.loads(subprocess.run(
            [*args.binary, "--format", args.format, "--gzip", "--path", str(gzip_path), "--repeats", str(args.repeats)],
            cwd=ROOT, check=True, text=True, capture_output=True,
        ).stdout)
        converted = 0
        failures = 0
        for _ in range(args.repeats):
            result = subprocess.run(
                [args.openbabel, f"-i{openbabel_format}", str(decompressed_path), "-osmi", "-O", "/dev/null"],
                check=False, text=True, capture_output=True,
            )
            line = next((line.strip() for line in result.stderr.splitlines() if line.strip().endswith(" molecule converted") or line.strip().endswith(" molecules converted")), "")
            try:
                count = int(line.split()[0])
            except (IndexError, ValueError):
                count = 0
            converted += count
            if result.returncode != 0 or not line:
                failures += 1

    openbabel = {
        "engine": "openbabel", "format": openbabel_format, "repeats": args.repeats,
        "records": converted, "failures": failures,
        "input_bytes": len(payload) * args.repeats,
        "comparison_boundary": "Open Babel CLI over the identical decompressed payload; gzip wrapper excluded",
    }
    report = {
        "schema_version": 1, "target_version": "1.0.9", "format": args.format,
        "compression": "gzip", "fixture": {
            "path": str(source.relative_to(ROOT)), "decompressed_bytes": len(payload),
            "compressed_bytes": len(compressed), "sha256": hashlib.sha256(payload).hexdigest(),
        }, "repeats": args.repeats, "expected_records": expected_records,
        "rows": {"chematic": chematic, "openbabel": openbabel},
        "comparison_boundary": {
            "chematic": f"Rust gzip-decompressing file-backed {args.format.upper()} reader",
            "openbabel": "Open Babel CLI over the identical decompressed payload; process startup included",
        }, "tool_versions": {"openbabel": version(args.openbabel)},
    }
    errors = []
    for name, row in (("chematic", chematic), ("openbabel", openbabel)):
        if row.get("records") != expected_records or row.get("failures") != 0:
            errors.append(f"{name} count mismatch: {row}")
    if errors:
        print("gzip Open Babel contract failures:", *errors, sep="\n")
        return 1
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else ROOT / args.output
        target.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
