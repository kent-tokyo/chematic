#!/usr/bin/env python3
"""Check gzip SDF/XYZ record accounting across chematic and RDKit boundaries."""

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
from bench_streaming_formats import measure, sdf_blocks, xyz_frames  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("sdf", "xyz"), default="sdf")
    parser.add_argument("--sdf", type=Path, default=Path("benchmarks/fixtures/streaming.sdf"))
    parser.add_argument("--xyz", type=Path, default=Path("benchmarks/fixtures/streaming.xyz"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--binary", nargs="+", default=["cargo", "run", "-p", "chematic-mol", "--example", "streaming_benchmark", "--offline", "--"])
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = args.sdf if args.format == "sdf" else args.xyz
    path = path if path.is_absolute() else ROOT / path
    payload = path.read_bytes()
    compressed = gzip.compress(payload, mtime=0)
    gzip_path = Path(f"/tmp/chematic-streaming-contract.{args.format}.gz")
    gzip_path.write_bytes(compressed)
    chematic = json.loads(subprocess.run([*args.binary, "--format", args.format, "--gzip", "--path", str(gzip_path), "--repeats", str(args.repeats)], cwd=ROOT, check=True, text=True, capture_output=True).stdout)
    frames = sdf_blocks(payload.decode()) if args.format == "sdf" else xyz_frames(payload.decode())
    rdkit = measure(args.format, frames, args.repeats, len(payload))
    expected = 2 * args.repeats
    if chematic.get("records") != expected or chematic.get("failures") != 0 or rdkit.get("records") != expected or rdkit.get("failures") != 0:
        print("gzip SDF contract failed:", json.dumps({"chematic": chematic, "rdkit": rdkit}))
        return 1
    report = {"schema_version": 1, "target_version": "1.0.9", "format": args.format, "compression": "gzip", "fixture": {"path": str(path.relative_to(ROOT)), "decompressed_bytes": len(payload), "compressed_bytes": len(compressed), "sha256": hashlib.sha256(payload).hexdigest()}, "repeats": args.repeats, "expected_records": expected, "rows": {"chematic": chematic, "rdkit": rdkit}, "comparison_boundary": {"chematic": f"Rust gzip-decompressing file-backed {args.format.upper()} reader", "rdkit": f"RDKit Python {args.format.upper()} block parser over the decompressed fixture"}}
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
