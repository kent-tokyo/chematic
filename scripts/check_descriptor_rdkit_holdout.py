#!/usr/bin/env python3
"""Run the descriptor compatibility profile against its structural holdout."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("holdout", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem import Descriptors
    import chematic

    rows = [json.loads(line) for line in args.holdout.read_text().splitlines() if line.strip()]
    checks = []
    for row in rows:
        smiles = row["smiles"]
        rd_mol = Chem.MolFromSmiles(smiles)
        ch_mol = chematic.from_smiles(smiles)
        rdkit_mw = Descriptors.MolWt(rd_mol)
        compat_mw = ch_mol.rdkit_mw
        isotope = any(atom.GetIsotope() for atom in rd_mol.GetAtoms())
        if isotope:
            passed = math.isnan(compat_mw)
            rule = "explicit isotope fails closed until isotope-table parity"
        else:
            passed = math.isfinite(compat_mw) and abs(compat_mw - rdkit_mw) <= 0.01
            rule = "unlabelled molecular weight within 0.01 Da"
        checks.append({
            "id": row["id"],
            "category": row["category"],
            "smiles": smiles,
            "rule": rule,
            "passed": passed,
            "chematic": compat_mw,
            "rdkit": rdkit_mw,
            "absolute_error": None if isotope else abs(compat_mw - rdkit_mw),
        })

    result = {
        "schema_version": 1,
        "profile": "rdkit_compat_unlabelled_v1",
        "corpus": str(args.holdout),
        "rows": len(checks),
        "passed": sum(check["passed"] for check in checks),
        "failed": sum(not check["passed"] for check in checks),
        "checks": checks,
    }
    args.json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if result["failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
