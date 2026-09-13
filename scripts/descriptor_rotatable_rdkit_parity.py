#!/usr/bin/env python3
"""Measure the declared rotatable-bond descriptor family against RDKit.

This is separate from the A1 core eight-field gate. Unsupported or failed
inputs stay in the denominator and are reported explicitly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_smiles(path: Path) -> list[str]:
    return [
        line.split()[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("corpus", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()

    from rdkit import Chem
    import rdkit
    from rdkit.Chem import Descriptors
    import chematic

    rows = load_smiles(args.corpus)
    matched = failed = unsupported = 0
    mismatches: list[dict] = []
    for smiles in rows:
        rd_mol = Chem.MolFromSmiles(smiles)
        try:
            ch_mol = chematic.from_smiles(smiles)
        except (TypeError, ValueError, RuntimeError):
            ch_mol = None
        if rd_mol is None or ch_mol is None:
            failed += 1
            continue
        value = getattr(ch_mol, "rotatable_bonds", None)
        if value is None:
            unsupported += 1
            continue
        expected = Descriptors.NumRotatableBonds(rd_mol)
        if value == expected:
            matched += 1
        else:
            mismatches.append({"smiles": smiles, "chematic": value, "rdkit": expected})

    result = {
        "schema_version": 1,
        "profile": "rdkit_rotatable_bonds_strict",
        "corpus": str(args.corpus),
        "rdkit_version": rdkit.__version__,
        "chematic_version": getattr(chematic, "__version__", "unknown"),
        "rows": len(rows),
        "matched": matched,
        "mismatched": len(rows) - failed - unsupported - matched,
        "failed": failed,
        "unsupported": unsupported,
        "coverage": matched / len(rows) if rows else 0.0,
        "mismatches": mismatches,
        "gate_passed": failed == 0 and unsupported == 0 and not mismatches,
        "boundary": "RDKit Strict NumRotatableBonds; this does not establish parity for other descriptor families",
    }
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
