#!/usr/bin/env python3
"""Audit RDKit MMFF partial-charge sums against input formal charge."""
import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    RDLogger.DisableLog("rdApp.*")
    raw = args.input.read_bytes()
    baseline = json.loads(raw)
    results = []
    for row in baseline["results"]:
        result = {"index": row["index"], "name": row["name"]}
        try:
            mol = Chem.MolFromSmiles(row["smiles"])
            if mol is None:
                raise ValueError("SMILES parse failed")
            mol = Chem.AddHs(mol)
            props = AllChem.MMFFGetMoleculeProperties(mol)
            if props is None:
                raise ValueError("RDKit MMFF parameters unavailable")
            formal = Chem.GetFormalCharge(mol)
            partial = sum(props.GetMMFFPartialCharge(i) for i in range(mol.GetNumAtoms()))
            result.update(formal_charge=formal, partial_charge_sum=partial, mismatch=abs(partial-formal)>1e-6)
        except Exception as exc:
            result["error"] = repr(exc)
        results.append(result)
    output = {"baseline_sha256": hashlib.sha256(raw).hexdigest(), "rdkit_version": rdBase.rdkitVersion, "tolerance": 1e-6, "results": results}
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"records {len(results)}, errors {sum('error' in r for r in results)}, mismatches {sum(r.get('mismatch', False) for r in results)}")


if __name__ == "__main__":
    main()
