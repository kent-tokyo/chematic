#!/usr/bin/env python3
"""Compare a fixed mixed-version SDF recognition fixture with installed engines."""

import hashlib
import json
from pathlib import Path

import cosmolkit
from rdkit import Chem


FIXTURE = Path(__file__).parent / "corpora/mixed-v2000-v3000-v1.sdf"
EXPECTED = [
    ["legacy", "old", None, 1, 0, "C"],
    ["M  END in title", "new", "M  END is data", 2, 1, "C[O-]"],
]


def main() -> None:
    if cosmolkit.__version__ != "0.5.0" or Chem.rdBase.rdkitVersion != "2025.09.3":
        raise SystemExit("requires COSMolKit 0.5.0 and RDKit 2025.09.3")

    cosmolkit_rows = []
    for index, record in enumerate(cosmolkit.SdfDataset.open(str(FIXTURE))):
        mol = record.molecule()
        cosmolkit_rows.append([
            record.title(),
            record.data_field("id"),
            record.data_field("note") if index else None,
            mol.num_atoms(),
            mol.num_bonds(),
            mol.to_smiles(),
        ])

    rdkit_rows = []
    for index, mol in enumerate(Chem.SDMolSupplier(str(FIXTURE), removeHs=False)):
        if mol is None:
            raise AssertionError(f"RDKit rejected source record {index}")
        rdkit_rows.append([
            mol.GetProp("_Name"),
            mol.GetProp("id"),
            mol.GetProp("note") if index else None,
            mol.GetNumAtoms(),
            mol.GetNumBonds(),
            Chem.MolToSmiles(mol),
        ])

    if cosmolkit_rows != EXPECTED or rdkit_rows != EXPECTED:
        raise AssertionError("installed comparator output differs from the fixed expectation")
    print(json.dumps({
        "fixture_sha256": hashlib.sha256(FIXTURE.read_bytes()).hexdigest(),
        "cosmolkit_version": cosmolkit.__version__,
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "rows": EXPECTED,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
