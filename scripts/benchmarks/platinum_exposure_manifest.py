#!/usr/bin/env python3
"""Generate a scaffold exclusion manifest from every record in an SDF."""
import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdf", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    molecules = []
    for index, mol in enumerate(Chem.SDMolSupplier(str(args.sdf), removeHs=True)):
        if mol is None:
            raise ValueError(f"unparseable SDF record {index}")
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        if Chem.MolFromSmiles(smiles) is None:
            raise ValueError(f"unparseable SMILES for record {index}")
        molecules.append({"smiles": smiles})
    result = {"source_sha256": hashlib.sha256(args.sdf.read_bytes()).hexdigest(), "rdkit_version": rdBase.rdkitVersion, "molecules": molecules}
    args.output.write_text(json.dumps(result, sort_keys=True) + "\n", encoding="utf-8")
    print(f"records: {len(molecules)}")


if __name__ == "__main__":
    main()
