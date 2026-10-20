#!/usr/bin/env python3
"""Measure Open Babel CLI conversion lanes for the shared benchmark corpus.

These numbers include process startup, CLI parsing, and output conversion.
They are intentionally kept separate from in-process chematic/RDKit API
measurements and must not be read as parser-only or writer-only timings.
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import time
from pathlib import Path


def read_smiles(path: Path, limit: int) -> bytes:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            rows.append(line.split()[0])
    if not rows:
        raise SystemExit("SMILES corpus is empty")
    return ("\n".join(rows[:limit]) + "\n").encode("utf-8")


def default_sdf() -> Path:
    spec = importlib.util.find_spec("rdkit")
    if spec is None or not spec.submodule_search_locations:
        raise SystemExit("RDKit is required to locate the shared egfr.sdf fixture")
    root = Path(list(spec.submodule_search_locations)[0])
    path = root / "Contrib" / "PBF" / "testData" / "egfr.sdf"
    if not path.exists():
        raise SystemExit(f"SDF fixture not found: {path}")
    return path


def measure(command: list[str], *, input_bytes: bytes | None, repeats: int) -> float:
    # One warm-up establishes that command construction and executable lookup
    # are valid; reported values still include a fresh CLI process each time.
    subprocess.run(command, input=input_bytes, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, check=True)
    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        subprocess.run(command, input=input_bytes, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=True)
        samples.append(time.perf_counter() - started)
    samples.sort()
    return samples[len(samples) // 2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles", type=Path, default=Path("scripts/chembl_accuracy_corpus_4999.smi"))
    parser.add_argument("--sdf", type=Path)
    parser.add_argument("--limit", type=int, default=5000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--obabel", default="obabel")
    args = parser.parse_args()
    if args.limit < 1 or args.repeats < 1:
        parser.error("--limit and --repeats must be positive")

    smiles = read_smiles(args.smiles, args.limit)
    sdf = args.sdf or default_sdf()
    rows = []
    for label, command in (
        ("smiles_to_smiles", [args.obabel, "-ismi", "-osmi", "-O", "/dev/null"]),
        ("smiles_canonical", [args.obabel, "-ismi", "-osmi", "--canonical", "-O", "/dev/null"]),
    ):
        elapsed = measure(command, input_bytes=smiles, repeats=args.repeats)
        rows.append((label, elapsed, args.limit))
    for label, command in (
        ("sdf_to_sdf", [args.obabel, "-isdf", str(sdf), "-osdf", "-O", "/dev/null"]),
    ):
        elapsed = measure(command, input_bytes=None, repeats=args.repeats)
        rows.append((label, elapsed, 365))

    print(f"Open Babel CLI benchmark; executable={args.obabel}; repeats={args.repeats}")
    print(f"SMILES corpus={args.smiles} records={args.limit}")
    print(f"SDF corpus={sdf} records=365")
    for label, elapsed, records in rows:
        print(f"{label}: {elapsed * 1000:.2f} ms total; {elapsed / records * 1e6:.2f} us/record; {records / elapsed:.0f} records/s")


if __name__ == "__main__":
    main()
