#!/usr/bin/env python3
"""Classify a scaffold holdout by graph features without scoring coordinates."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors


def classify(mol):
    if any(len(ring) >= 12 for ring in mol.GetRingInfo().AtomRings()):
        return "macrocycle"
    if rdMolDescriptors.CalcNumRotatableBonds(mol) >= 8:
        return "flexible"
    return "ordinary"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdf", type=Path)
    parser.add_argument("split", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    split_bytes = args.split.read_bytes()
    split = json.loads(split_bytes)
    if hashlib.sha256(args.sdf.read_bytes()).hexdigest() != split["source_sha256"]:
        parser.error("SDF hash differs from split")
    RDLogger.DisableLog("rdApp.*")
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    rows = []
    for record in split["rows"]:
        if record.get("split") != "holdout":
            continue
        result = {"index": record["index"], "name": record["name"],
                  "scaffold_sha256": record["scaffold_sha256"]}
        try:
            molecule = supplier[record["index"]]
            if molecule is None or molecule.GetProp("_Name") != record["name"]:
                raise ValueError("SDF identity or parsing failure")
            molecule = Chem.RemoveHs(molecule)
            result["group"] = classify(molecule)
        except Exception as error:
            result["group"] = "unclassified"
            result["error"] = repr(error)
        rows.append(result)
    result = {"source_sha256": split["source_sha256"],
              "split_sha256": hashlib.sha256(split_bytes).hexdigest(),
              "rdkit_version": rdBase.rdkitVersion,
              "rule": "macrocycle if an RDKit atom-ring basis member has at least 12 atoms; otherwise flexible if at least 8 RDKit strict rotatable bonds; otherwise ordinary",
              "counts": dict(Counter(row["group"] for row in rows)), "rows": rows}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
