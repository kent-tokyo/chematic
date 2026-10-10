"""Regenerate the coverage contract fixtures with RDKit 2026.03.1.

Requires that exact RDKit version; Rust tests consume the checked-in JSON and
need neither Python nor RDKit. Run from any directory in the checkout.
"""
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor

if rdBase.rdkitVersion != "2026.03.1":
    raise RuntimeError(f"expected RDKit 2026.03.1, got {rdBase.rdkitVersion}")

validation = Path(__file__).resolve().parents[1] / "validation"
smiles = [
    "C1CC2CCC1C2", "C12C3C4C1C5C2C3C45", "C1CCC2(CC1)CCCC2",
    "C1CCC2CCCCC2C1", "C1CC2CCC3CC(C1)C23", "C1CCC2(C1)CCCCC2",
    "C1CCCCCCCCCCC1", "C1CCCCCCC1", "C1CC2(C1)CCC1(CC2)CCC1",
    "c1ccc2c(c1)ccc1ccccc12", "C[C@H]1CCCC[C@@H]1O",
    "C1=CC=C2C(=C1)C=CC1=CC=CC=C12", "CC(C)(C)C(C)(C)C(C)(C)C",
    "C1CC2CC3CC1CC(C2)C3", "C1CCC2C3CCC4CCCCC4C3CCC12",
    "C1CCCCC1.C1CCCCC1", "F/C=C/C=C/C=C/F", "C1CCCC2CCCC(C1)C2",
]
for kind, count, prefix, suffix in [
    ("SP", 3, "F[Pt@", "](Cl)(Br)I"),
    ("TB", 20, "F[P@", "](Cl)(Br)(I)N"),
    ("OH", 30, "F[Co@", "](Cl)(Br)(I)(N)O"),
    ("SP", 3, "CC[Pt@", "](NCC)(N(C)C)Cl"),
    ("TB", 20, "CC[P@", "](Cl)(Br)(NC)N(C)C"),
    ("OH", 30, "CC[Co@", "](Cl)(Br)(I)(NC)N(C)C"),
    ("SP", 3, "N1CCN[Pt@", "]1(Cl)Br"),
    ("OH", 30, "N1CCN[Co@", "]1(Cl)(Br)(I)O"),
    ("SP", 3, "CC[Pt@", "](Cl)Br"),
    ("TB", 20, "CC[P@", "](Cl)(Br)N"),
    ("OH", 30, "CC[Co@", "](Cl)(Br)(I)N"),
]:
    smiles.extend(f"{prefix}{kind}{n}{suffix}" for n in range(1, count + 1))
depictions = []
for text in smiles:
    mol = Chem.MolFromSmiles(text)
    rdDepictor.Compute2DCoords(mol, forceRDKit=True)
    conf = mol.GetConformer()
    depictions.append({"smiles": text, "coords": [
        [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
        for i in range(mol.GetNumAtoms())
    ]})

pdb_rows = []
def record_pdb(label, block):
    reference = Chem.MolFromPDBBlock(block, proximityBonding=False)
    if reference is None:
        raise RuntimeError(f"RDKit could not read {label}")
    pdb_rows.append({"sequence": label, "pdb": block,
                     "smiles": Chem.MolToSmiles(reference),
                     "atoms": reference.GetNumAtoms(), "bonds": reference.GetNumBonds()})

for sequence in [*"ARNDCQEGHILKMFPSTWYV", "ACDEFGHIKLMNPQRSTVWY"]:
    record_pdb(sequence, Chem.MolToPDBBlock(Chem.MolFromSequence(sequence)))
for flavor in range(2, 10):
    for sequence in (["A", "C", "G", "U", "ACGU"] if flavor < 6 else
                     ["A", "C", "G", "T", "ACGT"]):
        record_pdb(f"{sequence},flavor={flavor}",
                   Chem.MolToPDBBlock(Chem.MolFromSequence(sequence, flavor=flavor)))
for name in [" C1 ", " N1 ", " O1 ", " S1 ", "CL  ", "BR  ", "ZN  ", "FE  ",
             "NA  ", "MG  ", "CA  ", "CO  ", "CU  ", "NI  ", "SE  ", " I1 ", "1C1 ", " Cl "]:
    charges = ["  ", "+ ", "++", "--", "1 ", "+0", "-0", " 1", " +", " -"] if name == "ZN  " else ["  "]
    for charge in charges:
        block = (f"HETATM    1 {name:4} UNL A   1       0.000   0.000   0.000"
                 f"  1.00  0.00            {charge}\nEND\n")
        record_pdb(f"atom name={name!r}, charge={charge!r}", block)
for filename, rows in [("depict", depictions), ("pdb-residue", pdb_rows)]:
    (validation / f"rdkit-2026.03.1-{filename}-contract.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8")
