#!/usr/bin/env python3
"""Check a bounded RDKit -> chematic -> RDKit V3000 semantic round trip.

This deliberately measures only the ordinary V3000 chemistry that chematic
models today: atom/bond topology, element, formal charge, isotope, and RDKit
isomeric canonical SMILES.  SGROUP, COLLECTION, ENDPTS, and ATTACH metadata
have a separate opaque-preservation contract and are not silently counted as
typed semantic interoperability here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = ROOT / "target" / "debug" / "chematic"

# Stable, non-query examples covering neutral, charged, isotopic, tetrahedral,
# and alkene stereo structures that RDKit can write as ordinary V3000 CTABs.
CASES = {
    "ethanol": "CCO",
    "benzene": "c1ccccc1",
    "pyridine": "n1ccccc1",
    "acetic_acid": "CC(=O)O",
    "acetate": "CC(=O)[O-]",
    "ammonium": "C[NH3+]",
    "trimethylammonium": "C[N+](C)(C)C",
    "nitrobenzene": "O=[N+]([O-])c1ccccc1",
    "sulfonamide": "CS(=O)(=O)N",
    "phosphonate": "COP(=O)(O)O",
    "chloroethane": "CCCl",
    "bromoethane": "CCBr",
    "iodomethane": "CI",
    "deuterated_ethanol": "[2H]CO",
    "carbon13_ethanol": "[13CH3]CO",
    "alanine": "N[C@@H](C)C(=O)O",
    "lactic_acid": "C[C@H](O)C(=O)O",
    "e_2_butene": "C/C=C/C",
    "z_2_butene": "C/C=C\\C",
    "cyclohexanol": "OC1CCCCC1",
}


def signature(mol: Chem.Mol) -> dict:
    atoms = [
        (atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetIsotope())
        for atom in mol.GetAtoms()
    ]
    bonds = sorted(
        (
            min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            str(bond.GetBondType()),
            str(bond.GetBondDir()),
        )
        for bond in mol.GetBonds()
    )
    return {
        "atom_count": mol.GetNumAtoms(),
        "bond_count": mol.GetNumBonds(),
        "atoms": atoms,
        "bonds": bonds,
        "canonical_isomeric_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
    }


def cli_provenance(cli: Path) -> dict[str, str]:
    """Record a relocatable identity for the tested executable."""
    with cli.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {
        "cli_name": cli.name,
        "cli_sha256": digest,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        parser.error(f"chematic CLI not found: {cli}")

    rows = []
    with tempfile.TemporaryDirectory(prefix="chematic-v3000-rdkit-") as raw_tmp:
        tmp = Path(raw_tmp)
        for case_id, smiles in CASES.items():
            original = Chem.MolFromSmiles(smiles)
            if original is None:
                raise RuntimeError(f"RDKit could not parse fixture {case_id}: {smiles}")
            source = Chem.MolToMolBlock(original, forceV3000=True)
            source_path = tmp / f"{case_id}.source.v3000"
            output_path = tmp / f"{case_id}.chematic.v3000"
            source_path.write_text(source, encoding="utf-8")
            run = subprocess.run(
                [
                    str(cli), "convert", "--input-format", "mol_v3000", "--output-format",
                    "mol_v3000", "--input", str(source_path), "--output", str(output_path),
                ],
                text=True,
                capture_output=True,
            )
            row = {"id": case_id, "smiles": smiles, "cli_exit_code": run.returncode}
            if run.returncode:
                row["error"] = run.stderr.strip()
                rows.append(row)
                continue
            written = output_path.read_text(encoding="utf-8")
            reopened = Chem.MolFromMolBlock(written, sanitize=True, removeHs=False)
            if reopened is None:
                row["error"] = "RDKit rejected chematic V3000 output"
                rows.append(row)
                continue
            before = signature(original)
            after = signature(reopened)
            row.update(
                {
                    "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
                    "output_sha256": hashlib.sha256(written.encode()).hexdigest(),
                    "before": before,
                    "after": after,
                    "semantic_equal": before == after,
                }
            )
            rows.append(row)

    failures = [row["id"] for row in rows if not row.get("semantic_equal", False)]
    result = {
        "schema_version": 1,
        "profile": "v3000_rdkit_ordinary_semantic_roundtrip_v1",
        "rdkit_version": Chem.rdBase.rdkitVersion,
        **cli_provenance(cli),
        "cases": len(rows),
        "semantic_equal": len(rows) - len(failures),
        "failures": failures,
        "scope": "RDKit-written ordinary V3000: topology, atom number/formal charge/isotope, bond type/direction, and isomeric identity",
        "not_claimed": [
            "SGROUP/COLLECTION typed semantic interoperability",
            "ENDPTS/ATTACH interpretation or haptic chemistry",
            "Indigo interoperability",
            "query or polymer expansion semantics",
        ],
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("cases", "semantic_equal", "failures")}))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
