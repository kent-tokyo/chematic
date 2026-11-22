#!/usr/bin/env python3
"""Create force-field input rows from a sealed holdout split and reference SDF."""
import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, RDLogger


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdf", type=Path)
    parser.add_argument("split", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    split_bytes = args.split.read_bytes()
    split = json.loads(split_bytes)
    source_sha256 = hashlib.sha256(args.sdf.read_bytes()).hexdigest()
    if source_sha256 != split["source_sha256"]:
        raise ValueError("reference SDF hash does not match split")
    selected = [r for r in split["rows"] if r.get("split") == "holdout"]
    indices = [r["index"] for r in selected]
    if not indices or len(set(indices)) != len(indices):
        raise ValueError("empty or duplicated holdout indices")
    RDLogger.DisableLog("rdApp.*")
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    results = []
    for row in selected:
        record = {"index": row["index"], "name": row["name"]}
        try:
            mol = supplier[row["index"]]
            if mol is None or mol.GetProp("_Name") != row["name"]:
                raise ValueError("reference parse or identity failure")
            record["smiles"] = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
        except Exception as exc:
            record["input_error"] = repr(exc)
        results.append(record)
    output = {
        "source_sha256": source_sha256,
        "split_sha256": hashlib.sha256(split_bytes).hexdigest(),
        "partition": "holdout",
        "results": results,
    }
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"rows": len(results), "input_errors": sum("input_error" in r for r in results)}))


if __name__ == "__main__":
    main()
