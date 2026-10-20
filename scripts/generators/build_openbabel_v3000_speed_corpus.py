#!/usr/bin/env python3
"""Build fixed V3000 MOL blocks; $$$$ only frames the native API batch."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("smiles", "output", "manifest"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--limit", type=int, default=1000)
    args = parser.parse_args()
    if args.offset < 0 or args.limit < 1:
        parser.error("nonnegative offset and positive limit required")
    source = args.smiles.read_bytes()
    rows = [line.split()[0] for line in source.decode().splitlines() if line.strip()]
    rows = rows[args.offset:args.offset + args.limit]
    if len(rows) != args.limit:
        raise ValueError("source shorter than requested slice")
    blocks, atoms = [], []
    for index, smiles in enumerate(rows, args.offset):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"source row {index} failed; no filtering allowed")
        mol.SetProp("_Name", f"frontier-{index:06}")
        rdDepictor.Compute2DCoords(mol, canonOrient=True, clearConfs=True, forceRDKit=True)
        block = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=True)
        if "V3000" not in block.splitlines()[3]:
            raise ValueError(f"source row {index} is not V3000")
        blocks.append(block + "$$$$\n")
        atoms.append(mol.GetNumAtoms())
    payload = "".join(blocks).encode()
    args.output.write_bytes(payload)
    args.manifest.write_text(json.dumps({
        "schema": "openbabel-v3000-block-corpus/v1",
        "source_path": str(args.smiles), "source_sha256": hashlib.sha256(source).hexdigest(),
        "source_offset": args.offset, "source_rows": len(rows), "rows_filtered": 0,
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "coordinate_generation": "Compute2DCoords(canonOrient=True, clearConfs=True, forceRDKit=True)",
        "mol_generation": "MolToMolBlock(includeStereo=True, forceV3000=True)",
        "framing": "native MOL API batch; delimiters split before timing in both harnesses; no SD fields",
        "output_sha256": hashlib.sha256(payload).hexdigest(), "output_bytes": len(payload),
        "total_atoms": sum(atoms), "minimum_atoms": min(atoms), "maximum_atoms": max(atoms),
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
