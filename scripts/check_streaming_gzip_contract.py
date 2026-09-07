#!/usr/bin/env python3
"""Check gzip SDF record accounting across chematic and RDKit boundaries."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from bench_streaming_formats import measure, sdf_blocks  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, default=Path("benchmarks/fixtures/streaming.sdf"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--binary", nargs="+", default=["cargo", "run", "-p", "chematic-mol", "--example", "streaming_benchmark", "--offline", "--"])
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = args.sdf if args.sdf.is_absolute() else ROOT / args.sdf
    payload = path.read_bytes()
    compressed = gzip.compress(payload, mtime=0)
    gzip_path = Path("/tmp/chematic-streaming-contract.sdf.gz")
    gzip_path.write_bytes(compressed)
    chematic = json.loads(subprocess.run([*args.binary, "--format", "sdf", "--gzip", "--path", str(gzip_path), "--repeats", str(args.repeats)], cwd=ROOT, check=True, text=True, capture_output=True).stdout)
    rdkit = measure("sdf", sdf_blocks(payload.decode()), args.repeats, len(payload))
    expected = 2 * args.repeats
    if chematic.get("records") != expected or chematic.get("failures") != 0 or rdkit.get("records") != expected or rdkit.get("failures") != 0:
        print("gzip SDF contract failed:", json.dumps({"chematic": chematic, "rdkit": rdkit}))
        return 1
    report = {"schema_version": 1, "target_version": "1.0.9", "format": "sdf", "compression": "gzip", "fixture": {"path": str(path.relative_to(ROOT)), "decompressed_bytes": len(payload), "compressed_bytes": len(compressed), "sha256": hashlib.sha256(payload).hexdigest()}, "repeats": args.repeats, "expected_records": expected, "rows": {"chematic": chematic, "rdkit": rdkit}, "comparison_boundary": {"chematic": "Rust gzip-decompressing file-backed BufRead reader", "rdkit": "RDKit Python block parser over the decompressed fixture"}}
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
