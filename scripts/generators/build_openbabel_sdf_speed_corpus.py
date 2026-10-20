#!/usr/bin/env python3
"""Build the fixed, unfiltered SDF throughput corpus from a SMILES prefix."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=1000)
    parser.add_argument("--offset", type=int, default=0)
    args = parser.parse_args()
    if args.limit <= 0 or args.offset < 0:
        parser.error("positive limit and nonnegative offset required")
    data = args.smiles.read_bytes()
    rows = [line.split()[0] for line in data.decode().splitlines() if line.strip()][args.offset:args.offset + args.limit]
    if len(rows) != args.limit:
        raise ValueError("source shorter than requested limit")
    blocks = []
    atom_counts = []
    for index, smiles in enumerate(rows):
        source_index = args.offset + index
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"source row {index} did not parse; no row filtering allowed")
        molecule.SetProp("_Name", f"frontier-{source_index:06}")
        rdDepictor.Compute2DCoords(molecule, canonOrient=True, clearConfs=True, forceRDKit=True)
        block = Chem.MolToMolBlock(molecule, includeStereo=True, forceV3000=False)
        if "V2000" not in block.splitlines()[3]:
            raise ValueError(f"source row {index} does not fit V2000")
        blocks.append(block + f"> <source_index>\n{source_index}\n\n$$$$\n")
        atom_counts.append(molecule.GetNumAtoms())
    payload = "".join(blocks).encode()
    args.output.write_bytes(payload)
    args.manifest.write_text(json.dumps({
        "schema": "openbabel-sdf-speed-corpus/v1", "source_path": str(args.smiles),
        "source_sha256": hashlib.sha256(data).hexdigest(), "source_rows": len(rows),
        "source_offset": args.offset,
        "rows_filtered": 0, "rdkit_version": Chem.rdBase.rdkitVersion,
        "coordinate_generation": "Compute2DCoords(canonOrient=True, clearConfs=True, forceRDKit=True)",
        "output_sha256": hashlib.sha256(payload).hexdigest(), "output_bytes": len(payload),
        "minimum_atoms": min(atom_counts), "maximum_atoms": max(atom_counts),
        "total_atoms": sum(atom_counts),
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
