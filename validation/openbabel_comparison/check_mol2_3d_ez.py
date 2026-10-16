#!/usr/bin/env python3
"""Compare E/Z preservation on the same pinned MOL2 bytes.

Requires RDKit, Open Babel, and a locally built CheMatic CLI. The expected
isomers are the source SMILES used to generate the committed conformers;
RDKit canonicalizes all three outputs to remove harmless SMILES spelling
differences. This measures only MOL2 -> SMILES E/Z preservation.
"""

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem


SOURCES = [
    r"F/C=C/F", r"F/C=C\F", r"C/C=C/C", r"C/C=C\C",
    r"Cl/C(Br)=C(/I)F", r"Cl/C(Br)=C(\I)F",
    r"F/C=C/C=C/F", r"F/C=C\C=C/F",
    r"O/N=C/c1ccccc1", r"O/N=C\c1ccccc1",
    r"C/C=C(/F)Cl", r"C/C=C(\F)Cl",
]
FIXTURES = Path(__file__).resolve().parents[2] / "crates/chematic-mol/tests/data"


def canonical(value):
    mol = Chem.MolFromSmiles(value.split()[0]) if value.strip() else None
    return Chem.MolToSmiles(mol) if mol is not None else None


def run(args):
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    return result.stdout if result.returncode == 0 else ""


def compare(cli, source, fixture, label):
    contents = fixture.read_bytes()
    expected = canonical(source)
    chematic = canonical(run([
        str(cli), "convert", "--input-format", "mol2",
        "--output-format", "smiles", "--input", str(fixture),
    ]))
    openbabel = canonical(run(["obabel", "-imol2", str(fixture), "-osmi"]))
    return {
        "source": source,
        "fixture": label,
        "sha256": hashlib.sha256(contents).hexdigest(),
        "expected": expected,
        "chematic": chematic,
        "openbabel": openbabel,
        "chematic_preserved": chematic == expected,
        "openbabel_preserved": openbabel == expected,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--holdout", action="store_true",
                        help="also generate five new conformers per isomer")
    args = parser.parse_args()
    if not args.cli.is_file():
        parser.error(f"CLI missing: {args.cli}")

    rows = []
    for index, source in enumerate(SOURCES):
        fixture = FIXTURES / f"mol2-3d-ez-{index:02d}.mol2"
        rows.append(compare(args.cli, source, fixture, fixture.name))
    if args.holdout:
        with tempfile.TemporaryDirectory() as directory:
            mol_path = Path(directory) / "source.mol"
            mol2_path = Path(directory) / "converted.mol2"
            for seed_index in range(5):
                for index, source in enumerate(SOURCES):
                    seed = 400000 + seed_index * 100 + index
                    molecule = Chem.AddHs(Chem.MolFromSmiles(source))
                    if AllChem.EmbedMolecule(molecule, randomSeed=seed) != 0:
                        raise RuntimeError(f"Embedding failed for {source}, seed {seed}")
                    mol_path.write_text(Chem.MolToMolBlock(molecule, forceV3000=False))
                    source_readback = Chem.MolFromMolFile(str(mol_path), removeHs=True)
                    if source_readback is None or Chem.MolToSmiles(source_readback) != canonical(source):
                        raise RuntimeError(f"Source identity changed for {source}, seed {seed}")
                    converted = subprocess.run(
                        ["obabel", "-imol", str(mol_path), "-omol2", "-O", str(mol2_path)],
                        capture_output=True, text=True, check=False,
                    )
                    if converted.returncode != 0:
                        raise RuntimeError(f"MOL2 generation failed: {converted.stderr}")
                    rows.append(compare(args.cli, source, mol2_path, f"seed-{seed}"))
    summary = {
        "rdkit": rdBase.rdkitVersion,
        "openbabel": run(["obabel", "-V"]).strip(),
        "cases": len(rows),
        "chematic_preserved": sum(row["chematic_preserved"] for row in rows),
        "openbabel_preserved": sum(row["openbabel_preserved"] for row in rows),
        "rows": rows,
    }
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    if summary["chematic_preserved"] != len(rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
