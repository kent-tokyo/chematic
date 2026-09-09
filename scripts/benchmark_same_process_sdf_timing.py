#!/usr/bin/env python3
"""Measure equivalent SDF parsing in one Python process.

The semantic contracts remain the correctness gate.  This companion records
paired timing context only: both implementations see the same bytes, are
warmed up, and are timed in alternating order.  It is not a claim about other
formats, machines, or the installed Open Babel CLI boundary.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import statistics
import tempfile
import time
from pathlib import Path

from benchmark_version import workspace_version

ROOT = Path(__file__).resolve().parents[1]


def chematic_parse(path: Path) -> int:
    import chematic

    stream = chematic.iter_sdf_batched(str(path), batch_size=100)
    records = sum(len(batch) for batch in stream)
    manifest = json.loads(stream.manifest_json())
    if records != 2 or int(manifest["rejected_records"]) != 0:
        raise RuntimeError("chematic SDF timing fixture did not parse as two records")
    return records


def rdkit_parse(payload: bytes) -> int:
    from rdkit import Chem

    records = 0
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(payload), sanitize=True, removeHs=False)
    for molecule in supplier:
        if molecule is None:
            raise RuntimeError("RDKit SDF timing fixture yielded a rejected record")
        records += 1
    if records != 2:
        raise RuntimeError("RDKit SDF timing fixture did not parse as two records")
    return records


def percentile(values: list[int], fraction: float) -> int:
    ordered = sorted(values)
    index = min(len(ordered) - 1, int((len(ordered) - 1) * fraction))
    return ordered[index]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, default=Path("benchmarks/fixtures/streaming.sdf"))
    parser.add_argument("--repeats", type=int, default=40)
    parser.add_argument("--inner", type=int, default=25)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 5 or args.inner <= 0:
        raise SystemExit("--repeats must be at least 5 and --inner must be positive")

    source = args.sdf if args.sdf.is_absolute() else ROOT / args.sdf
    payload = source.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    with tempfile.TemporaryDirectory(prefix="chematic-same-process-sdf-timing-") as directory:
        input_path = Path(directory) / "input.sdf"
        input_path.write_bytes(payload)

        # Import and warm both implementations before the measured rounds.
        chematic_parse(input_path)
        rdkit_parse(payload)
        chematic_ns: list[int] = []
        rdkit_ns: list[int] = []
        for round_index in range(args.repeats):
            if round_index % 2 == 0:
                started = time.perf_counter_ns()
                for _ in range(args.inner):
                    chematic_parse(input_path)
                chematic_ns.append(time.perf_counter_ns() - started)
                started = time.perf_counter_ns()
                for _ in range(args.inner):
                    rdkit_parse(payload)
                rdkit_ns.append(time.perf_counter_ns() - started)
            else:
                started = time.perf_counter_ns()
                for _ in range(args.inner):
                    rdkit_parse(payload)
                rdkit_ns.append(time.perf_counter_ns() - started)
                started = time.perf_counter_ns()
                for _ in range(args.inner):
                    chematic_parse(input_path)
                chematic_ns.append(time.perf_counter_ns() - started)

    chematic_per_op = [value / args.inner for value in chematic_ns]
    rdkit_per_op = [value / args.inner for value in rdkit_ns]
    chematic_p50 = statistics.median(chematic_per_op)
    rdkit_p50 = statistics.median(rdkit_per_op)
    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "status": "local-verified",
        "gate": "same_process_sdf_equivalent_timing",
        "fixture": {"path": str(source.relative_to(ROOT)), "bytes": len(payload), "sha256": digest},
        "repeats": args.repeats,
        "inner_iterations": args.inner,
        "timing": {
            "unit": "nanoseconds per parse operation",
            "chematic": {"p50": round(chematic_p50, 1), "p95": percentile([int(v) for v in chematic_per_op], 0.95)},
            "rdkit": {"p50": round(rdkit_p50, 1), "p95": percentile([int(v) for v in rdkit_per_op], 0.95)},
            "rdkit_over_schematic_p50": round(rdkit_p50 / chematic_p50, 4),
        },
        "semantic_contract": {"chematic_records": 2, "rdkit_records": 2, "failure_count": 0},
        "boundaries": [
            "same CPython process and same fixture bytes",
            "chematic uses the current source-built Python extension and file-backed iterator",
            "RDKit uses ForwardSDMolSupplier over the same bytes",
            "single SDF fixture; timing is context, not a general ranking",
        ],
        "tool_versions": {},
    }
    from rdkit import rdBase

    import chematic

    report["tool_versions"] = {"rdkit": rdBase.rdkitVersion, "chematic": chematic.__version__}
    target = args.output if args.output.is_absolute() else ROOT / args.output
    target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
