#!/usr/bin/env python3
"""Observe original nonperiodic BGF fields and indexed bonds across round trips.

Extra generated remarks and reciprocal connections are allowed. ID renumbering
and explicit/implicit ORDER changes are reported separately from graph meaning.
This independent fixed-column observer does not use the chematic parser.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
ATOM_FIELDS = ("kind", "name", "residue_name", "chain", "residue_id", "position_angstrom",
               "force_field_type", "max_bonds", "lone_pairs", "partial_charge")
AXES = (*ATOM_FIELDS, "source_ids", "bond_graph", "order_encoding", "description", "force_field", "remarks")


def observe(text):
    lines = text.splitlines()
    if not lines or lines[0].split() != ["BIOGRF", "200"]:
        raise ValueError("missing BIOGRF 200 output")
    record = dict(version=200, atoms=[], connectivity=[], description=None, force_field=None, remarks=[])
    ended = False
    for line in lines[1:]:
        if not line.strip():
            continue
        if ended:
            raise ValueError("content after END")
        if line.strip() == "END":
            ended = True
        elif line.startswith(("ATOM  ", "HETATM")):
            if len(line) < 80:
                raise ValueError("truncated atom output")
            a = dict(kind=line[:6].strip(), id=int(line[7:12]), name=line[13:18].strip(),
                     residue_name=line[19:22].strip(), chain=line[23:24].strip(), residue_id=line[25:30].strip(),
                     position_angstrom=[float(line[i:i+10]) for i in (30,40,50)],
                     force_field_type=line[61:66].strip(), max_bonds=int(line[66:69]),
                     lone_pairs=int(line[69:71]), partial_charge=float(line[72:80]))
            if not all(map(math.isfinite, [*a["position_angstrom"], a["partial_charge"]])):
                raise ValueError("nonfinite atom output")
            record["atoms"].append(a)
        elif line.startswith("CONECT"):
            values = [int(line[i:i+6]) for i in range(6, len(line.rstrip()), 6)]
            record["connectivity"].append(dict(atom=values[0], neighbors=values[1:], orders=None))
        elif line.startswith("ORDER "):
            values = [int(line[i:i+6]) for i in range(6, len(line.rstrip()), 6)]
            c = record["connectivity"][-1]
            if c["atom"] != values[0] or c["orders"] is not None or len(values)-1 != len(c["neighbors"]):
                raise ValueError("unmatched ORDER output")
            c["orders"] = values[1:]
        elif line.startswith("DESCRP "):
            record["description"] = line[7:].rstrip()
        elif line.startswith("FORCEFIELD "):
            record["force_field"] = line[11:].rstrip()
        elif line.startswith("REMARK "):
            record["remarks"].append(line[7:].rstrip())
        elif not line.startswith("FORMAT "):
            raise ValueError("unknown output record")
    if not ended or not record["atoms"]:
        raise ValueError("incomplete output")
    # Validate references even when no graph comparison is requested.
    axes(record)
    return record


def axes(record):
    atoms = record["atoms"]
    ids = {a["id"]: i for i,a in enumerate(atoms)}
    if len(ids) != len(atoms):
        raise ValueError("duplicate atom IDs")
    graph = {}
    centers = set()
    for c in record["connectivity"]:
        if c["atom"] in centers or len(set(c["neighbors"])) != len(c["neighbors"]):
            raise ValueError("duplicate connectivity center/neighbor")
        centers.add(c["atom"])
        center = ids[c["atom"]]
        orders = c["orders"] if c["orders"] is not None else [1]*len(c["neighbors"])
        if len(orders) != len(c["neighbors"]):
            raise ValueError("ORDER length")
        for neighbor, order in zip(c["neighbors"], orders):
            edge = tuple(sorted((center, ids[neighbor])))
            if order not in (0,1,2,3):
                raise ValueError("unsupported bond order")
            order = max(order, 1)
            if edge[0] == edge[1] or order > 3 or (edge in graph and graph[edge] != order):
                raise ValueError("invalid indexed bond graph")
            graph[edge] = order
    return {**{f:[a[f] for a in atoms] for f in ATOM_FIELDS}, "source_ids":[a["id"] for a in atoms],
            "bond_graph":sorted((a,b,o) for (a,b),o in graph.items()),
            "order_encoding":record["connectivity"], "description":record["description"],
            "force_field":record["force_field"], "remarks":record["remarks"]}


def run(command, text):
    try:
        p = subprocess.run(command, input=text, text=True, encoding="utf-8", capture_output=True, timeout=30, check=False)
    except (OSError, subprocess.TimeoutExpired) as e:
        return dict(status="execution_failed", error=str(e))
    receipt = dict(returncode=p.returncode, stdout=p.stdout, stderr=p.stderr)
    if p.returncode:
        return dict(status="rejected", **receipt)
    try:
        record = observe(p.stdout)
    except (ValueError, KeyError, IndexError) as e:
        return dict(status="unobservable", error=str(e), **receipt)
    return dict(status="observed", record=record, **receipt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--obabel", type=Path, default=Path("/opt/homebrew/bin/obabel"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    packet_path = ROOT / "validation/bgf/cases.json"
    packet = json.loads(packet_path.read_text(encoding="utf-8"))
    commands = dict(chematic=[str(args.cli.resolve()), "bgf"], openbabel=[str(args.obabel.resolve()), "-ibgf", "-obgf", "-l1"])
    totals = {k:{a:0 for a in AXES} for k in commands}
    eligible = {a:0 for a in AXES}
    rows = []
    for case in packet["valid"]:
        expected = axes(case["expected"])
        row = dict(id=case["id"], input_sha256=hashlib.sha256(case["input"].encode()).hexdigest())
        for a in AXES:
            eligible[a] += expected[a] is not None and expected[a] != []
        for engine, command in commands.items():
            result = run(command, case["input"])
            got = axes(result["record"]) if result["status"] == "observed" else {}
            result["preserved"] = {a:None if expected[a] is None or expected[a] == [] else
                                   (all(v in got.get(a,[]) for v in expected[a]) if a == "remarks" else got.get(a) == expected[a]) for a in AXES}
            for a,v in result["preserved"].items():
                totals[engine][a] += v is True
            result["all_source_fields_preserved"] = result["status"] == "observed" and all(v is not False for v in result["preserved"].values())
            row[engine] = result
        rows.append(row)
    sources = ["crates/chematic-mol/src/bgf.rs", "crates/chematic-mol/src/bgf/reader.rs", "crates/chematic-mol/src/bgf/writer.rs", "crates/chematic-cli/src/main.rs", "crates/chematic-cli/src/dispatch.rs"]
    report = dict(schema_version=1, scope="Original nonperiodic BGF profile; exact source fields; indexed graph normalized for reciprocal entries and ORDER 0/1; no general superiority claim",
                  case_count=len(rows), timeout_seconds=30, commands=commands, eligible_counts=eligible, preserved_counts=totals,
                  observed_counts={k:sum(c[k]["status"] == "observed" for c in rows) for k in commands},
                  complete_source_preservation={k:sum(c[k]["all_source_fields_preserved"] for c in rows) for k in commands},
                  versions={"chematic":subprocess.check_output([commands["chematic"][0], "--version"], text=True).strip(), "openbabel":subprocess.check_output([commands["openbabel"][0], "-V"], text=True).strip()},
                  packet_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest(), source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                  binary_sha256={k:hashlib.sha256(Path(v[0]).read_bytes()).hexdigest() for k,v in commands.items()}, cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k:report[k] for k in ("case_count", "observed_counts", "eligible_counts", "preserved_counts", "complete_source_preservation")}, indent=2))
    return int(report["complete_source_preservation"]["chematic"] != len(rows))


if __name__ == "__main__":
    raise SystemExit(main())
