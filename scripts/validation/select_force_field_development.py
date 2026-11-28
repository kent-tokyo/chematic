#!/usr/bin/env python3
"""Select a fixed element-stratified COD development sample, without scoring."""
import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase


CAPS = {"boron": 20, "selenium": 20, "phosphorus_sulfur": 30, "other": 30}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdf", type=Path)
    parser.add_argument("split", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    split = json.loads(args.split.read_text(encoding="utf-8"))
    if hashlib.sha256(args.sdf.read_bytes()).hexdigest() != split["source_sha256"]:
        parser.error("SDF hash differs from the split")
    RDLogger.DisableLog("rdApp.*")
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    groups = {name: [] for name in CAPS}
    failures = []
    for row in split["rows"]:
        if row.get("split") != "development":
            continue
        molecule = supplier[row["index"]]
        if molecule is None:
            failures.append({"index": row["index"], "name": row["name"], "error": "parse failed"})
            continue
        if molecule.GetProp("_Name") != row["name"]:
            raise ValueError("reference identity differs from split")
        molecule = Chem.RemoveHs(molecule)
        if molecule.GetNumHeavyAtoms() > 30 or any(len(ring) >= 12 for ring in molecule.GetRingInfo().AtomRings()):
            continue
        atoms = {atom.GetAtomicNum() for atom in molecule.GetAtoms()}
        category = ("boron" if 5 in atoms else "selenium" if 34 in atoms else
                    "phosphorus_sulfur" if atoms & {15, 16} else "other")
        groups[category].append({"index": row["index"], "name": row["name"],
                                 "smiles": Chem.MolToSmiles(molecule, isomericSmiles=True),
                                 "category": category})
    selected = []
    for category, cap in CAPS.items():
        groups[category].sort(key=lambda row: (hashlib.sha256(row["name"].encode()).digest(), row["name"]))
        selected.extend(groups[category][:cap])
    result = {"source_sha256": split["source_sha256"],
              "split_sha256": hashlib.sha256(args.split.read_bytes()).hexdigest(),
              "rdkit_version": rdBase.rdkitVersion,
              "selection": "development partition only; at most 30 heavy atoms; no atom-ring basis member >=12; category precedence B, Se, P/S, other; lowest SHA256(name) within category, caps 20/20/30/30",
              "counts": {name: min(len(rows), CAPS[name]) for name, rows in groups.items()},
              "selection_parse_failures": failures, "results": selected}
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["counts"]))


if __name__ == "__main__":
    main()
