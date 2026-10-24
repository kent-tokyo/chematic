#!/usr/bin/env python3
"""Compare declared CONFIG data using an independent output-table observer.

Names, units and force fields are never inferred. Missing velocity/force data
is distinct from a numeric difference; absent input vectors are not evidence
that an engine can preserve present vectors. Only valid cases reach Open Babel.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
AXES = ("title", "levcfg", "imcon", "declared_atom_count", "names", "indices",
        "positions", "velocities", "forces", "cell_vectors", "header_extra")


def observe(text: str) -> dict:
    lines = iter(text.splitlines())
    try:
        title = next(lines)
        header = next(lines).split()
        if len(header) not in (2, 3):
            raise ValueError("unrecognized control header")
        level, mode = map(int, header[:2])
        if level not in (0, 1, 2) or mode not in (0, 1, 2, 3):
            raise ValueError("unrecognized level/boundary")
        count = int(header[2]) if len(header) == 3 else None

        def vector() -> list[float]:
            values = [float(v.replace("D", "e").replace("d", "e")) for v in next(lines).split()]
            if len(values) != 3 or not all(map(math.isfinite, values)):
                raise ValueError("incomplete/nonfinite vector")
            return values

        cell = [vector() for _ in range(3)] if mode else None
        atoms = []
        for line in lines:
            fields = line.split()
            if not fields:
                raise ValueError("blank atom header")
            atoms.append(dict(name=fields[0], index=int(fields[1]) if len(fields)>1 else None,
                              header_extra=fields[2:], position=vector(),
                              velocity=vector() if level else None,
                              force=vector() if level == 2 else None))
        if not atoms or (count is not None and count != len(atoms)):
            raise ValueError("missing/mismatched atom table")
    except StopIteration as error:
        raise ValueError("truncated record") from error
    return dict(title=title, levcfg=level, imcon=mode, declared_atom_count=count,
                cell_vectors=cell, atoms=atoms)


def axes(record: dict) -> dict:
    result = {k: record[k] for k in AXES if k in record}
    for axis, key in [("names","name"),("indices","index"),("positions","position"),
                      ("velocities","velocity"),("forces","force"),("header_extra","header_extra")]:
        result[axis] = [a[key] for a in record["atoms"]]
    return result


def run(command: list[str], text: str) -> dict:
    try:
        output = subprocess.run(command, input=text, text=True, encoding="utf-8",
                                capture_output=True, timeout=60, check=False)
    except (OSError, subprocess.TimeoutExpired) as error:
        return dict(status="execution_failed", error=str(error))
    receipt = dict(returncode=output.returncode, stdout=output.stdout, stderr=output.stderr)
    if output.returncode != 0:
        return dict(status="rejected", **receipt)
    try:
        record = observe(output.stdout)
    except (ValueError, IndexError) as error:
        return dict(status="unobservable", error=str(error), **receipt)
    return dict(status="observed", record=record, **receipt)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--obabel", type=Path, default=Path("/opt/homebrew/bin/obabel"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    packet_path = ROOT / "validation/dlpoly_config/cases.json"
    packet = json.loads(packet_path.read_text(encoding="utf-8"))
    commands = dict(chematic=[str(args.cli.resolve()),"dlpoly-config"],
                    openbabel=[str(args.obabel.resolve()),"-iCONFIG","-oCONFIG","-l1"])
    totals = {engine: {axis: 0 for axis in AXES} for engine in commands}
    rows = []
    for case in packet["cases"]:
        expected = axes(case["expected"])
        row = dict(id=case["id"],input_sha256=hashlib.sha256(case["input"].encode()).hexdigest())
        for engine,command in commands.items():
            result = run(command,case["input"])
            got = axes(result["record"]) if result["status"] == "observed" else {}
            result["preserved"] = {axis: axis in got and got[axis] == expected[axis] for axis in AXES}
            for axis,passed in result["preserved"].items():
                totals[engine][axis] += passed
            row[engine] = result
        rows.append(row)
    source_names = ["crates/chematic-mol/src/dlpoly.rs","crates/chematic-mol/src/dlpoly/reader.rs",
                    "crates/chematic-cli/src/main.rs","crates/chematic-cli/src/dispatch.rs"]
    report = dict(schema_version=1,scope="original single-CONFIG packet; exact source-unit data; not force-field validation",
                  case_count=len(rows),commands=commands,timeout_seconds=60,
                  versions={"chematic":subprocess.check_output([str(args.cli.resolve()),"--version"],text=True).strip(),
                            "openbabel":subprocess.check_output([str(args.obabel.resolve()),"-V"],text=True).strip()},
                  packet_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest(),
                  source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_names},
                  binary_sha256={k:hashlib.sha256(Path(v[0]).read_bytes()).hexdigest() for k,v in commands.items()},
                  observed_counts={k:sum(row[k]["status"] == "observed" for row in rows) for k in commands},
                  preserved_counts=totals,cases=rows)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:report[k] for k in ["case_count","observed_counts","preserved_counts"]},indent=2))
    return int(any(n != len(rows) for n in totals["chematic"].values()))


if __name__ == "__main__":
    raise SystemExit(main())
