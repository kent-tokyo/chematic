#!/usr/bin/env python3
"""Authored CJSON → MOL2 charge development packet, never a holdout or ranking.

Run with python3 -m scripts.check_cjson_mol2_formal_charges --probe BINARY
--obabel /opt/homebrew/bin/obabel --output NEW.json. The probe is the
chematic-mol cjson_mol2_report example. Both commands receive identical bytes.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform
import struct

from scripts.check_interchange_preservation import bounded_run, identity, require, strict_json


def packet():
    """Original minimal inputs; no upstream or sealed corpus is imported."""
    specs = [
        ("neutral-control", [2], [0], [0, 0, 0], [], [], None),
        ("sodium-chloride", [11, 17], [1, -1],
         [1.123456789012345, -0.0, 0, -2.123456789012345, 0, 0], [], [],
         [0.987654321012345, -0.987654321012345]),
        ("calcium-chloride", [20, 17, 17], [2, -1, -1],
         [0, 0, 0, -3, 0, 0, 3, 0, 0], [], [], None),
        ("oxide", [8], [-2], [0, 0, 0], [], [], [-0.7]),
        ("zinc", [30], [2], [0, 0, 0], [], [], None),
        ("nitrate", [7, 8, 8, 8], [1, 0, -1, -1],
         [0, 0, 0, 1.2, 0, 0, -0.6, 1.1, 0, -0.6, -1.1, 0],
         [0, 1, 0, 2, 0, 3], [2, 1, 1], [0.7, -0.5, -0.6, -0.6]),
    ]
    cases = []
    for name, elements, charges, coords, bonds, orders, partial in specs:
        source = {"chemicalJson": 1,
                  "atoms": {"elements": {"number": elements},
                            "formalCharges": charges, "coords": {"3d": coords}},
                  "bonds": {"connections": {"index": bonds}, "order": orders}}
        if partial is not None:
            source["partialCharges"] = {"USER_CHARGES": partial}
        cases.append((name, source))
    return cases


ELEMENTS = {"He": 2, "N": 7, "O": 8, "Na": 11, "Cl": 17, "Ca": 20, "Zn": 30}


def mol2_snapshot(text):
    """Independent, deliberately small oracle for this packet's MOL2 profile.

    Recognizes UNITY charges anywhere around BOND, with N.4 fallback. Rejects
    ambiguity instead of silently substituting a default for malformed input.
    It does not use either toolkit's reader or claim general MOL2 perception.
    """
    sections = {}
    section = None
    for line in text.splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line[9:].strip()
            require(section not in sections, "duplicate MOL2 section")
            sections[section] = []
        elif line.strip() and not line.lstrip().startswith("#"):
            require(section is not None, "data before MOL2 header")
            sections[section].append(line)
    header = sections["MOLECULE"]
    count, bond_count = map(int, header[1].split()[:2])
    atoms = [line.split() for line in sections["ATOM"]]
    require(len(atoms) == count, "atom count mismatch")
    ids = [int(a[0]) for a in atoms]
    require(len(set(ids)) == count, "duplicate atom ID")
    indices = {a: i for i, a in enumerate(ids)}
    charges = [1 if a[5] == "N.4" else 0 for a in atoms]
    attributes = iter(sections.get("UNITY_ATOM_ATTR", []))
    seen = set()
    for line in attributes:
        atom, n = map(int, line.split())
        require(atom in indices and atom not in seen and n > 0, "invalid UNITY header")
        seen.add(atom)
        charge_seen = False
        for _ in range(n):
            key, value = next(attributes).split()
            require(key == "charge" and not charge_seen, "unexpected UNITY attribute")
            charges[indices[atom]] = int(value)
            charge_seen = True
    bonds = []
    for line in sections.get("BOND", []):
        _, a, b, order = line.split()[:4]
        x, y = sorted((indices[int(a)], indices[int(b)]))
        bonds.append([x, y, int(order)])
    require(len(bonds) == bond_count, "bond count mismatch")
    return {"elements": [ELEMENTS[a[5].split('.')[0]] for a in atoms],
            "formal_charges": charges,
            "coordinates": [float(v) for a in atoms for v in a[2:5]],
            "bonds": sorted(bonds),
            "partial_charges": [float(a[8]) if len(a) > 8 else None for a in atoms],
            "charge_method": header[3]}


def same_bits(a, b):
    return len(a) == len(b) and all(
        x is not None and y is not None and struct.pack("!d", x) == struct.pack("!d", y)
        for x, y in zip(a, b))


def sdf_graph(text):
    """Read just the explicit indexed graph/charges of Open Babel's V2000 output."""
    lines = text.splitlines()
    require("V2000" in lines[3], "expected V2000 readback")
    count, bond_count = int(lines[3][:3]), int(lines[3][3:6])
    atoms = lines[4:4 + count]
    charge_codes = {0: 0, 1: 3, 2: 2, 3: 1, 5: -1, 6: -2, 7: -3}
    charges = [charge_codes[int(a[36:39])] for a in atoms]
    seen = set()
    for line in lines[4 + count + bond_count:]:
        if line.startswith("M  CHG"):
            row = line.split()
            n = int(row[2])
            require(len(row) == 3 + 2 * n, "invalid M CHG count")
            for i in range(n):
                atom, charge = map(int, row[3 + i * 2:5 + i * 2])
                require(1 <= atom <= count and atom not in seen, "invalid M CHG atom")
                seen.add(atom)
                charges[atom - 1] = charge
    bonds = []
    for line in lines[4 + count:4 + count + bond_count]:
        a, b, order = int(line[:3]) - 1, int(line[3:6]) - 1, int(line[6:9])
        bonds.append([*sorted((a, b)), order])
    require("M  END" in lines and lines.count("$$$$") == 1, "incomplete or multiple records")
    return {"elements": [ELEMENTS[a[31:34].strip()] for a in atoms],
            "formal_charges": charges, "bonds": sorted(bonds)}


def evaluate(source, snapshot):
    atoms = source["atoms"]
    ends, orders = source["bonds"]["connections"]["index"], source["bonds"]["order"]
    bonds = sorted([*sorted(ends[i * 2:i * 2 + 2]), order]
                   for i, order in enumerate(orders))
    result = {"elements_and_indexed_bonds": snapshot["elements"] == atoms["elements"]["number"]
              and snapshot["bonds"] == bonds,
              "formal_charges": snapshot["formal_charges"] == atoms["formalCharges"],
              "coordinate_bits": same_bits(snapshot["coordinates"], atoms["coords"]["3d"])}
    if "partialCharges" in source:
        result["partial_charge_bits"] = same_bits(
            snapshot["partial_charges"], source["partialCharges"]["USER_CHARGES"])
        result["charge_method"] = snapshot["charge_method"] == "USER_CHARGES"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", required=True, type=Path)
    parser.add_argument("--obabel", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    require(not args.output.exists(), "output already exists")
    version = bounded_run([str(args.obabel), "-V"])
    require(version["exit_code"] == 0 and b"Open Babel 3.2.1" in version["stdout"],
            "expected pinned Open Babel 3.2.1")
    rows = []
    cross_reads = []
    cases = []
    for case, source in packet():
        data = json.dumps(source, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
        cases.append({"id": case, "source": source, "sha256": hashlib.sha256(data).hexdigest()})
        for engine, command in [
            ("chematic", [str(args.probe)]),
            ("openbabel", [str(args.obabel), "-icjson", "-omol2"]),
        ]:
            run = bounded_run(command, data, timeout=10, max_output=1048576)
            row = {"case": case, "engine": engine, "command": command,
                   "resource_status": run["resource_status"], "exit_code": run["exit_code"],
                   "stdout": run["stdout"].decode("utf8", errors="replace"),
                   "stderr": run["stderr"].decode("utf8", errors="replace")}
            row["terminal"] = "process_failure"
            if run["resource_status"] != "completed":
                row["terminal"] = run["resource_status"]
            elif run["exit_code"] == 0:
                try:
                    output = row["stdout"]
                    if engine == "chematic":
                        report = strict_json(output)
                        output = report["output"]
                        if output is None:
                            row["terminal"] = "refused"
                    if output:
                        row["snapshot"] = mol2_snapshot(output)
                        row["fields"] = evaluate(source, row["snapshot"])
                        row["terminal"] = ("packet_fields_preserved" if all(row["fields"].values())
                                           else "packet_field_loss")
                        if engine == "chematic":
                            command = [str(args.obabel), "-imol2", "-osdf"]
                            readback = bounded_run(command, output.encode(), timeout=10,
                                                   max_output=1048576)
                            cross = {"case": case, "command": command,
                                     "resource_status": readback["resource_status"],
                                     "exit_code": readback["exit_code"],
                                     "stdout": readback["stdout"].decode("utf8", errors="replace"),
                                     "stderr": readback["stderr"].decode("utf8", errors="replace")}
                            cross["terminal"] = "process_failure"
                            if readback["resource_status"] != "completed":
                                cross["terminal"] = readback["resource_status"]
                            elif readback["exit_code"] == 0:
                                try:
                                    graph = sdf_graph(cross["stdout"])
                                    # Graph/charge only: the SDF writer rounds coordinates.
                                    expected = {key: row["snapshot"][key] for key in graph}
                                    cross["snapshot"] = graph
                                    cross["terminal"] = ("emitted_graph_preserved" if graph == expected
                                                         else "emitted_graph_loss")
                                except (ValueError, KeyError, IndexError) as error:
                                    cross["terminal"] = "unresolved_output"
                                    cross["error"] = str(error)
                            cross_reads.append(cross)
                    elif row["terminal"] != "refused":
                        row["terminal"] = "parse_write_failure"
                except (ValueError, KeyError, IndexError, StopIteration) as error:
                    row["terminal"] = "unresolved_output"
                    row["error"] = str(error)
            rows.append(row)
    result = {"schema": "chematic.cjson-mol2-charge-development.v1",
              "scope": "Authored development packet; indexed atom order and explicit bonds, formal charges, coordinate bits, and declared partial-charge bits/method only. Not stereochemistry, implicit hydrogens, general typing, full format preservation, an independent holdout, or toolkit superiority.",
              "provenance": "Original synthetic inputs authored for this regression, 2026-10-21; recipe: packet() in this script; no third-party or sealed corpus input.",
              "license": "MIT OR Apache-2.0", "redistribution": "allowed: original minimal inputs",
              "minimization": "One neutral control, mono/divalent ions, and one connected charged graph. Coordinates are supplied, not optimized molecular geometries.",
              "oracle": {"name": "independent Python MOL2 packet parser", "version": 1,
                         "python_version": platform.python_version(),
                         "script": identity(__file__), "tolerance": "exact integers and binary64 bits; angstrom coordinates; elementary-charge values; no numeric tolerance"},
              "platform": platform.platform(),
              "engines": {"chematic": identity(args.probe), "openbabel": identity(args.obabel)},
              "openbabel_version": version["stdout"].decode().strip(),
              "expected": "Each case's source arrays define exact required fields; compare explicit indexed graph and charges; binary64 coordinates and declared partial charges. Names/typing/stereo and any unrequested charge assignment are outside scoring.",
              "cases": cases, "rows": rows, "openbabel_reads_chematic_mol2": cross_reads,
              "counts": {e: dict(Counter(r["terminal"] for r in rows if r["engine"] == e))
                         for e in ["chematic", "openbabel"]}}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
