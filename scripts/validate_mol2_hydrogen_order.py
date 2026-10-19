#!/usr/bin/env python3
"""Independent MOL2 atom/coordinate/stereo checks for explicit-H ordering."""

import argparse
import hashlib
import json
from pathlib import Path
import random
import subprocess

import chematic as ch
import chematic._mol2_hydrogen_order as implementation
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

SOURCES = [
    "C",
    "CC",
    "CCO",
    "CC(C)O",
    "CC(=O)O",
    "NCCO",
    "C1CCCCC1",
    "c1ccccc1",
    "c1ccncc1",
    "C[C@H](O)F",
    "C[C@@H](O)F",
    "F/C=C/F",
    "F/C=C\\F",
    "C[C@H]1CCCCO1",
]


def atom_type(atom):
    symbol = atom.GetSymbol()
    if symbol in ("C", "N", "O"):
        if atom.GetIsAromatic():
            return symbol + ".ar"
        if atom.GetHybridization() == Chem.HybridizationType.SP2:
            return symbol + ".2"
        return symbol + ".3"
    return symbol


def block(mol, order):
    new = {old: i + 1 for i, old in enumerate(order)}
    conf = mol.GetConformer()
    lines = [
        "@<TRIPOS>MOLECULE",
        "authored_geometry",
        f"{len(order)} {mol.GetNumBonds()} 1 0 0",
        "SMALL",
        "USER_CHARGES",
        "",
        "@<TRIPOS>ATOM",
    ]
    for old in order:
        atom = mol.GetAtomWithIdx(old)
        point = conf.GetAtomPosition(old)
        lines.append(
            f"{new[old]} {atom.GetSymbol()}{old} {point.x:.12g} {point.y:.12g} {point.z:.12g} "
            f"{atom_type(atom)} 1 RES {((old % 7) - 3) / 100:.5f}"
        )
    lines.append("@<TRIPOS>BOND")
    for i, bond in enumerate(mol.GetBonds(), 1):
        kind = "ar" if bond.GetIsAromatic() else str(int(bond.GetBondTypeAsDouble()))
        lines.append(
            f"{i} {new[bond.GetBeginAtomIdx()]} {new[bond.GetEndAtomIdx()]} {kind}"
        )
    lines.extend(
        ["@<TRIPOS>SUBSTRUCTURE", f"1 RES {new[order[-1]]} GROUP 0 **** **** 0 ROOT"]
    )
    return "\n".join(lines) + "\n"


def read(text):
    mol = Chem.MolFromMol2Block(text, sanitize=True, removeHs=False)
    assert mol is not None, "independent RDKit MOL2 read failed"
    return mol


def edges(mol):
    return {
        (
            min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
        ): (b.GetIsAromatic(), b.GetBondTypeAsDouble())
        for b in mol.GetBonds()
    }


def check(before, after, order, tolerance):
    assert [after.GetAtomWithIdx(i).GetAtomicNum() for i in range(len(order))] == [
        before.GetAtomWithIdx(i).GetAtomicNum() for i in order
    ]
    inverse = {old: new for new, old in enumerate(order)}
    expected = {
        (min(inverse[i], inverse[j]), max(inverse[i], inverse[j])): value
        for (i, j), value in edges(before).items()
    }
    assert edges(after) == expected, "mapped graph mismatch"
    maximum = 0.0
    for new, old in enumerate(order):
        a = before.GetConformer().GetAtomPosition(old)
        b = after.GetConformer().GetAtomPosition(new)
        maximum = max(maximum, *(abs(a[k] - b[k]) for k in range(3)))
    assert maximum <= tolerance, "mapped coordinate mismatch"
    assert Chem.MolToInchi(before) == Chem.MolToInchi(
        after
    ), "independent stereo identity mismatch"
    return maximum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obabel", default="/opt/homebrew/bin/obabel")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rng = random.Random(2733)
    rows = []
    for index, smiles in enumerate(SOURCES):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert AllChem.EmbedMolecule(mol, randomSeed=2733 + index) == 0
        for variant in range(12):
            source_order = list(range(mol.GetNumAtoms()))
            if variant:
                rng.shuffle(source_order)
            source = block(mol, source_order)
            result = ch.order_mol2_hydrogens_checked(source)
            before = read(source)
            assert Chem.MolToInchi(before) == Chem.MolToInchi(mol)
            after = read(result.mol2)
            parent = None
            for atom in after.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    assert atom.GetDegree() == 1
                    assert (
                        atom.GetNeighbors()[0].GetIdx() == parent
                    ), "H not next to parent group"
                else:
                    parent = atom.GetIdx()
            rdkit_error = check(before, after, result.new_to_old, 1e-12)
            converted = subprocess.run(
                [args.obabel, "-imol2", "-omol2"],
                input=result.mol2,
                text=True,
                capture_output=True,
                check=True,
            )
            assert "1 molecule converted" in converted.stderr
            ob_error = check(
                after, read(converted.stdout), list(range(mol.GetNumAtoms())), 5.01e-5
            )
            assert ch.order_mol2_hydrogens_checked(result.mol2).mol2 == result.mol2
            rows.append(
                {
                    "smiles": smiles,
                    "variant": variant,
                    "atoms": mol.GetNumAtoms(),
                    "rdkit_coordinate_max_error": rdkit_error,
                    "openbabel_coordinate_max_error": ob_error,
                    "new_to_old": list(result.new_to_old),
                    "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
                    "output_sha256": hashlib.sha256(result.mol2.encode()).hexdigest(),
                }
            )
    probe = subprocess.run(
        [args.obabel, "-ismi", "-:CCO", "-h", "-omol2"],
        text=True,
        capture_output=True,
        check=True,
    )
    baseline = ch.order_mol2_hydrogens_checked(probe.stdout)
    assert baseline.new_to_old != tuple(range(len(baseline.new_to_old)))
    evidence = {
        "scope": "Python source editing over published native reader; no native rebuild",
        "rdkit_version": rdBase.rdkitVersion,
        "openbabel_version": subprocess.check_output(
            [args.obabel, "-V"], text=True
        ).strip(),
        "implementation_sha256": hashlib.sha256(
            Path(implementation.__file__).read_bytes()
        ).hexdigest(),
        "comparisons": len(rows),
        "rdkit_coordinate_max_error": max(
            row["rdkit_coordinate_max_error"] for row in rows
        ),
        "openbabel_coordinate_max_error": max(
            row["openbabel_coordinate_max_error"] for row in rows
        ),
        "baseline_obabel_ethanol": {
            "source": probe.stdout,
            "grouped": baseline.mol2,
            "new_to_old": list(baseline.new_to_old),
        },
        "cases": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evidence, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in evidence.items()
                if k not in ("cases", "baseline_obabel_ethanol")
            }
        )
    )


if __name__ == "__main__":
    main()
