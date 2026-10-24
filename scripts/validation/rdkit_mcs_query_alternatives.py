#!/usr/bin/env python3
"""Record RDKit FMCS query alternatives for fixed comparator development cases."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdFMCS

CASES = (
    (("CO", "NO", "SO"), "Any", "OrderExact", None),
    (("CO", "NO", "SO"), "AnyHeavyAtom", "OrderExact", None),
    (("C=N", "CN", "C#N"), "Elements", "Any", None),
    (("CO", "NO"), "Any", "OrderExact", "PO"),
    (("C=N", "CN"), "Elements", "Any", "C#N"),
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1, found {rdBase.rdkitVersion}")
    rows = []
    for smiles, atom_mode, bond_mode, unseen in CASES:
        params = rdFMCS.MCSParameters()
        params.AtomTyper = getattr(rdFMCS.AtomCompare, f"Compare{atom_mode}")
        params.BondTyper = getattr(rdFMCS.BondCompare, f"Compare{bond_mode}")
        mols = [Chem.MolFromSmiles(value) for value in smiles]
        if any(mol is None for mol in mols):
            raise ValueError(f"RDKit rejected {smiles}")
        result = rdFMCS.FindMCS(mols, params)
        query = Chem.MolFromSmarts(result.smartsString)
        if result.canceled or query is None or not all(mol.HasSubstructMatch(query) for mol in mols):
            raise ValueError(f"FMCS query failed to re-match an input: {smiles}")
        unseen_matches = None
        if unseen is not None:
            other = Chem.MolFromSmiles(unseen)
            if other is None:
                raise ValueError(f"RDKit rejected {unseen}")
            unseen_matches = other.HasSubstructMatch(query)
            if unseen_matches:
                raise ValueError(f"unseen alternative unexpectedly matched: {unseen}")
        rows.append({
            "smiles": smiles,
            "atom_compare": atom_mode,
            "bond_compare": bond_mode,
            "atoms": result.numAtoms,
            "bonds": result.numBonds,
            "smarts": result.smartsString,
            "all_inputs_match": True,
            "unseen_smiles": unseen,
            "unseen_matches": unseen_matches,
        })
    report = {
        "rdkit_version": rdBase.rdkitVersion,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cases": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
