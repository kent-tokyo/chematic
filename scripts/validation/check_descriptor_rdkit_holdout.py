#!/usr/bin/env python3
"""Run the descriptor compatibility profile against its structural holdout."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


FIELD_RULES = {
    "molecular_weight": ("rdkit_mw", 0.01, 1e-6),
    "hba": ("rdkit_hba", 0.0, 0.0),
    "hbd": ("hbd", 0.0, 0.0),
    "tpsa": ("tpsa", 0.1, 1e-6),
    "logp": ("logp", 0.01, 1e-6),
    "molar_refractivity": ("molar_refractivity", 0.01, 1e-6),
    "fsp3": ("fsp3", 0.001, 1e-6),
    "aromatic_ring_count": ("rdkit_aromatic_ring_count", 0.0, 0.0),
}


def descriptor_values(rd_mol):
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    return {
        "molecular_weight": Descriptors.MolWt(rd_mol),
        "hba": Lipinski.NumHAcceptors(rd_mol),
        "hbd": Lipinski.NumHDonors(rd_mol),
        "tpsa": rdMolDescriptors.CalcTPSA(rd_mol, includeSandP=True),
        "logp": Crippen.MolLogP(rd_mol),
        "molar_refractivity": Crippen.MolMR(rd_mol),
        "fsp3": rdMolDescriptors.CalcFractionCSP3(rd_mol),
        "aromatic_ring_count": rdMolDescriptors.CalcNumAromaticRings(rd_mol),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("holdout", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--fields", nargs="+", choices=tuple(FIELD_RULES),
                        default=("molecular_weight",))
    args = parser.parse_args()

    from rdkit import Chem
    import chematic

    rows = [json.loads(line) for line in args.holdout.read_text().splitlines() if line.strip()]
    checks = {field: [] for field in args.fields}
    for row in rows:
        smiles = row["smiles"]
        rd_mol = Chem.MolFromSmiles(smiles)
        if rd_mol is None:
            raise SystemExit(f"RDKit rejected holdout row {row.get('id', '<unknown>')}")
        ch_mol = chematic.from_smiles(smiles)
        expected = descriptor_values(rd_mol)
        for field in args.fields:
            attr, tolerance, strict_tolerance = FIELD_RULES[field]
            expected_value = expected[field]
            actual_value = getattr(ch_mol, attr)
            actual_finite = (isinstance(actual_value, (int, float))
                             and not isinstance(actual_value, bool)
                             and math.isfinite(actual_value))
            error = abs(float(actual_value) - float(expected_value)) if actual_finite else None
            checks[field].append({
                "id": row["id"], "category": row["category"], "smiles": smiles,
                "tolerance": tolerance, "strict_tolerance": strict_tolerance,
                "passed": actual_finite and error <= tolerance,
                "strict_passed": actual_finite and error <= strict_tolerance,
                "chematic": actual_value, "rdkit": expected_value,
                "absolute_error": error,
            })

    field_results = {
        field: {
            "rows": len(items),
            "passed": sum(item["passed"] for item in items),
            "strict_passed": sum(item["strict_passed"] for item in items),
            "failed": sum(not item["strict_passed"] for item in items),
            "checks": items,
        }
        for field, items in checks.items()
    }
    result = {
        "schema_version": 2,
        "profile": "rdkit_compat_v2",
        "corpus": str(args.holdout),
        "rows": len(rows),
        "fields": field_results,
        "passed": sum(item["passed"] for items in checks.values() for item in items),
        "strict_passed": sum(item["strict_passed"] for items in checks.values() for item in items),
        "failed": sum(item["failed"] for item in field_results.values()),
    }
    args.json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if result["failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
