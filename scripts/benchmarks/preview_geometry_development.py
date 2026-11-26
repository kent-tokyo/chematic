#!/usr/bin/env python3
"""Compare coarse preview diagnostics with a separate development screen."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import chematic
from chematic import chematic as extension
from rdkit import Chem, rdBase

from platinum_development_probe import graph_matches
from platinum_geometry import check_geometry, prepare_geometry_checks


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("coordinates", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inputs = json.loads(args.input.read_text(encoding="utf-8"))["results"]
    indices = {row["index"]: row for row in inputs}
    if len(indices) != len(inputs):
        parser.error("duplicate input indices")
    results = []
    for source in args.coordinates:
        source_hash = digest(source)
        data = json.loads(source.read_text(encoding="utf-8"))
        if data["development_sha256"] != digest(args.input):
            parser.error("coordinate source does not match input")
        if {r["index"] for r in data["results"]} != set(indices) or len(data["results"]) != len(indices):
            parser.error("incomplete or duplicate coordinate source")
        for row in data["results"]:
            selection = indices[row["index"]]
            result = {"index": row["index"], "name": row["name"], "category": selection["category"],
                      "source_sha256": source_hash}
            try:
                if row["name"] != selection["name"] or row.get("error"):
                    raise ValueError("coordinate source failure or identity mismatch")
                coords = row["coordinates"]
                mol = Chem.AddHs(Chem.MolFromSmiles(selection["smiles"]))
                other = chematic.from_smiles(selection["smiles"]).add_hydrogens()
                if not graph_matches(other, mol):
                    raise ValueError("atom or bond graph mismatch")
                checks = prepare_geometry_checks(mol)
                count, worst, rate, clashes = other._timed_preview_geometry(coords)
                result["diagnostics"] = {"bond_count": count, "worst_bond_length": worst,
                                         "bond_violation_rate_50pct": rate, "gross_clash_count": clashes}
                result["reference"] = check_geometry(coords, checks)
                result["coverage"] = {"bonds_supported": len(checks[0]), "bonds_total": mol.GetNumBonds(),
                                      "angles_supported": len(checks[1]),
                                      "angles_total": sum(a.GetDegree() * (a.GetDegree() - 1) // 2 for a in mol.GetAtoms())}
            except Exception as error:
                result["error"] = repr(error)
            results.append(result)
    output = {"input_sha256": digest(args.input), "coordinate_sha256": [digest(p) for p in args.coordinates],
              "extension_sha256": digest(Path(extension.__file__)), "script_sha256": digest(Path(__file__)),
              "python": platform.python_version(), "platform": platform.platform(), "rdkit": rdBase.rdkitVersion,
              "results": results}
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
