#!/usr/bin/env python3
"""Compare declared Tinker data, without deriving a chemical graph or force field.

Only original, in-profile valid records are fed to Open Babel. The comparator's
custom-type writer and single-bond reader options are explicit. Exact numeric
equality is a serialization contract, not a chemical accuracy claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
AXES = ("title", "indices", "names", "atom_types", "coordinates", "connectivity", "cell")


def observe(text: str) -> dict:
    """Small independent table observer for outputs, not a chemistry oracle."""
    lines = text.splitlines()
    header = lines.pop(0).split(maxsplit=1)
    count = int(header[0])
    cell = None
    if len(lines[0].split()) == 6:
        try:
            values = list(map(float, lines[0].split()))
        except ValueError:
            pass
        else:
            cell = dict(lengths_angstrom=values[:3], angles_degrees=values[3:])
            lines.pop(0)
    rows = [line.split() for line in lines if line.strip()]
    if count < 1 or len(rows) != count:
        raise ValueError("incomplete or multiple records")
    atoms = []
    for row in rows:
        xyz = list(map(float, row[2:5]))
        if len(xyz) != 3 or not all(map(math.isfinite, xyz)):
            raise ValueError("non-finite or missing coordinates")
        atoms.append(dict(index=int(row[0]), name=row[1], position_angstrom=xyz,
                          atom_type=int(row[5]), neighbors=list(map(int, row[6:]))))
    return dict(title=header[1] if len(header) == 2 else "", cell=cell, atoms=atoms)


def axes(record: dict) -> dict:
    atoms = record["atoms"]
    return dict(title=record["title"], indices=[a["index"] for a in atoms],
                names=[a["name"] for a in atoms], atom_types=[a["atom_type"] for a in atoms],
                coordinates=[a["position_angstrom"] for a in atoms],
                connectivity=[sorted(a["neighbors"]) for a in atoms], cell=record["cell"])


def run(command: list[str], text: str) -> dict:
    try:
        result = subprocess.run(command, input=text, text=True, encoding="utf-8",
                                capture_output=True, timeout=15, check=False)
    except (OSError, subprocess.TimeoutExpired) as error:
        return dict(status="execution_failed", error=str(error))
    receipt = dict(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
    if result.returncode != 0:
        return dict(status="rejected", **receipt)
    try:
        record = observe(result.stdout)
    except (ValueError, IndexError) as error:
        return dict(status="unobservable", error=str(error), **receipt)
    return dict(status="observed", record=record, **receipt)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--obabel", type=Path, default=Path("/opt/homebrew/bin/obabel"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    packet_path = ROOT / "validation/tinker_xyz/cases.json"
    packet = json.loads(packet_path.read_text(encoding="utf-8"))
    commands = dict(chematic=[str(args.cli.resolve()), "tinker-xyz"],
                    openbabel=[str(args.obabel.resolve()), "-itxyz", "-otxyz", "-as", "-xc"])
    rows = []
    totals = {engine: {axis: 0 for axis in AXES} for engine in commands}
    for case in packet["cases"]:
        expected = axes(case["expected"])
        row = dict(id=case["id"], input_sha256=hashlib.sha256(case["input"].encode()).hexdigest())
        for engine, command in commands.items():
            result = run(command, case["input"])
            observed = axes(result["record"]) if result["status"] == "observed" else {}
            result["preserved"] = {axis: axis in observed and observed[axis] == expected[axis] for axis in AXES}
            for axis, preserved in result["preserved"].items():
                totals[engine][axis] += preserved
            row[engine] = result
        rows.append(row)
    sources = ["crates/chematic-mol/src/tinker.rs", "crates/chematic-mol/src/tinker/reader.rs",
               "crates/chematic-cli/src/main.rs", "crates/chematic-cli/src/dispatch.rs"]
    report = dict(schema_version=1, scope="single-frame Tinker XYZ; exact declared fields; original synthetic packet",
                  case_count=len(rows), commands=commands,
                  versions={"chematic": subprocess.check_output([str(args.cli.resolve()), "--version"], text=True).strip(),
                            "openbabel": subprocess.check_output([str(args.obabel.resolve()), "-V"], text=True).strip()},
                  packet_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest(),
                  binary_sha256={k: hashlib.sha256(Path(v[0]).read_bytes()).hexdigest() for k,v in commands.items()},
                  source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                  preserved_counts=totals, cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(dict(case_count=len(rows), preserved_counts=totals), indent=2))
    return int(any(value != len(rows) for value in totals["chematic"].values()))


if __name__ == "__main__":
    raise SystemExit(main())
