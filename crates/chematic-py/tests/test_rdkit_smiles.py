"""``Mol.rdkit_smiles``: RDKit 2026.03.1's canonical SMILES.

Expected strings are ``Chem.MolToSmiles(Chem.MolFromSmiles(input))`` from
RDKit 2026.03.1; RDKit itself is not needed to run these tests.
"""

import pytest

import chematic

CASES = [
    ("OC(=O)[C@@H]1CCCN1", "O=C(O)[C@@H]1CCCN1"),
    ("[H][C@](F)(Cl)Br", "F[C@@H](Cl)Br"),
    ("C[C@H]1CC[C@@H](C)CC1", "C[C@H]1CC[C@@H](C)CC1"),
    ("F/C=C\\C=C/C=C\\F", "F\\C=C/C=C\\C=C/F"),
    ("C1=CC=CC=C1", "c1ccccc1"),
    ("N1C=CC=C1", "c1cc[nH]c1"),
    ("CN(=O)=O", "C[N+](=O)[O-]"),
    ("[2H]C([2H])([2H])O", "[2H]C([2H])([2H])O"),
    ("O.O.O.[Na+].[Na+].[O-]C(=O)C(=O)[O-]", "O.O.O.O=C([O-])C(=O)[O-].[Na+].[Na+]"),
    ("C%(100)CC%(100)", "C1CC1"),
]


@pytest.mark.parametrize("smiles,expected", CASES)
def test_rdkit_smiles_matches_rdkit(smiles, expected):
    assert chematic.from_smiles(smiles).rdkit_smiles == expected


def test_rdkit_smiles_leaves_native_smiles_unchanged():
    mol = chematic.from_smiles("OC(=O)[C@@H]1CCCN1")
    native = mol.smiles
    assert mol.rdkit_smiles == "O=C(O)[C@@H]1CCCN1"
    assert mol.smiles == native


@pytest.mark.parametrize(
    "smiles,kind",
    [
        ("CN(C)(C)(C)C", "sanitization failed"),
        ("F[Pt@SP1](Cl)(Br)I", "unsupported input"),
    ],
)
def test_rdkit_smiles_raises_typed_value_error(smiles, kind):
    with pytest.raises(ValueError, match=f"RDKit-compatible SMILES: {kind}"):
        chematic.from_smiles(smiles).rdkit_smiles
