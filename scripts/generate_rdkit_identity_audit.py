#!/usr/bin/env python3
"""Generate the RDKit-pinned identity keys required by the sealed cohort gate.

The result is JSONL rather than an opaque cache so the cohort preparer can
verify one source row against one canonical/parent/scaffold audit row.  Invalid
SMILES are fatal: silently dropping them would invalidate the source hash and
row-accounting contract.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, help="SMILES-first .smi file")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jsonl-smiles-key", help="read SMILES from this key instead of a SMILES-first text file")
    args = parser.parse_args()
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem.Scaffolds import MurckoScaffold
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError as error:
        parser.error(f"RDKit is required for this oracle audit: {error}")

    lines = args.source.read_text(encoding="utf-8").splitlines()
    if args.jsonl_smiles_key:
        rows = []
        for line_no, line in enumerate(lines, 1):
            if not line.strip():
                continue
            try:
                value = json.loads(line).get(args.jsonl_smiles_key)
            except json.JSONDecodeError as error:
                parser.error(f"source row {line_no}: invalid JSON: {error.msg}")
            if not isinstance(value, str) or not value:
                parser.error(f"source row {line_no}: missing string {args.jsonl_smiles_key!r}")
            rows.append(value)
    else:
        rows = [line.split(maxsplit=1)[0] for line in lines if line.strip()]
    rendered: list[str] = []
    for index, smiles in enumerate(rows, 1):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            parser.error(f"source row {index}: RDKit cannot parse {smiles!r}")
        parent = rdMolStandardize.FragmentParent(molecule)
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=parent)
        rendered.append(json.dumps({
            "input_smiles": smiles,
            "canonical_smiles": Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True),
            "parent_smiles": Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True),
            "scaffold_smiles": scaffold,
            "oracle": {"engine": "RDKit", "version": rdBase.rdkitVersion},
        }, sort_keys=True))
    args.output.write_text("\n".join(rendered) + "\n", encoding="utf-8")
    print(json.dumps({"rows": len(rows), "output": str(args.output), "oracle": {"engine": "RDKit", "version": rdBase.rdkitVersion}}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
