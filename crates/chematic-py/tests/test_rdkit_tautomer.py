"""Mol.rdkit_canonical_tautomer / rdkit_tautomers / rdkit_tautomer_score,
pinned to RDKit 2026.03.1's rdMolStandardize.TautomerEnumerator() values
(cross-checked against RDKit when it is installed)."""

import pytest

import chematic

CASES = [
    # input, Canonicalize, ScoreTautomer(input), sorted Enumerate SMILES
    ("Oc1ccccn1", "O=c1cccc[nH]1", 100, ["O=C1CC=CC=N1", "O=c1cccc[nH]1", "Oc1ccccn1"]),
    (
        "C/N=C(\\NC)NCc1ccc(I)cc1",
        "CN=C(NC)NCc1ccc(I)cc1",
        257,
        ["CN=C(NC)NCc1ccc(I)cc1", "CNC(=NCc1ccc(I)cc1)NC"],
    ),
    (
        "O=C1CCCCC/C1=C\\c1ccccc1",
        "O=C1CCCCC/C1=C\\c1ccccc1",
        253,
        ["O=C1CCCCC/C1=C\\c1ccccc1", "OC1=CCCCC/C1=C\\c1ccccc1"],
    ),
    ("NCCCC(N)OP(O)O", "NCCCC(N)O[PH](=O)O", 0, ["NCCCC(N)OP(O)O", "NCCCC(N)O[PH](=O)O"]),
]


@pytest.mark.parametrize("smi,canon,score,tauts", CASES)
def test_pinned(smi, canon, score, tauts):
    m = chematic.from_smiles(smi)
    assert m.rdkit_canonical_tautomer() == canon
    assert m.rdkit_tautomer_score() == score
    assert m.rdkit_tautomers() == tauts
    assert m.rdkit_tautomer_status() == "Completed"


def test_max_transforms_status():
    m = chematic.from_smiles("CC1CC(=O)C2=C(O)c3c(O)ccc(O)c3CC2C1")
    assert m.rdkit_tautomer_status() == "MaxTransformsReached"
    assert len(m.rdkit_tautomers()) == 237
    assert m.rdkit_canonical_tautomer() == "CC1CC(=O)C2C(=O)c3c(O)ccc(O)c3CC2C1"


@pytest.mark.parametrize("smi", [c[0] for c in CASES])
def test_against_rdkit(smi):
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem.MolStandardize import rdMolStandardize

    rm = Chem.MolFromSmiles(smi)
    te = rdMolStandardize.TautomerEnumerator()
    m = chematic.from_smiles(smi)
    assert m.rdkit_canonical_tautomer() == Chem.MolToSmiles(te.Canonicalize(rm))
    assert m.rdkit_tautomers() == sorted(Chem.MolToSmiles(t) for t in te.Enumerate(rm))
    assert m.rdkit_tautomer_score() == rdMolStandardize.TautomerEnumerator.ScoreTautomer(rm)
