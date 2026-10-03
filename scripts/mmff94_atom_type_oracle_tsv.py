#!/usr/bin/env python3
"""Emit per-atom MMFF94 types from a pinned RDKit for the source Rust gate."""

from __future__ import annotations

import argparse
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.03.6":
        raise SystemExit(f"expected RDKit 2026.03.6, got {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")
    with args.output.open("w", encoding="utf-8") as output:
        for i, line in enumerate(args.corpus.read_text(encoding="utf-8").splitlines()):
            smiles = line.split()[0]
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                output.write(f"{i}\tparse\n")
                continue
            mol = Chem.AddHs(mol)
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
            if props is None:
                output.write(f"{i}\tunsupported\n")
                continue
            types = ",".join(str(props.GetMMFFAtomType(j)) for j in range(mol.GetNumAtoms()))
            output.write(f"{i}\t{types}\n")


if __name__ == "__main__":
    main()
