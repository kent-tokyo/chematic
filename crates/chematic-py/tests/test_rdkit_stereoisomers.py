"""Mol.rdkit_stereoisomers against RDKit 2026.03.1 EnumerateStereoisomers."""

import pytest

import chematic


def test_unique_unassigned_isomers_match_rdkit():
    mol = chematic.from_smiles("BrC=C[C@H]1OC(C2)(F)C2(Cl)C1")
    smiles = [m.rdkit_smiles for m in mol.rdkit_stereoisomers()]
    assert len(smiles) == 8 and smiles == sorted(smiles)
    assert smiles[0] == "F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2"


def test_existing_api_unchanged_and_options_checked():
    mol = chematic.from_smiles("CC(O)CC")
    assert len(mol.enumerate_stereoisomers()) >= 2  # native API, unchanged
    assert [m.rdkit_smiles for m in mol.rdkit_stereoisomers()] == ["CC[C@@H](C)O", "CC[C@H](C)O"]
    with pytest.raises(ValueError):
        mol.rdkit_stereoisomers(only_unassigned=False)


def test_random_sample_matches_rdkit_default_seed():
    # 11 unassigned centres = 2048 combinations > maxIsomers=1024: RDKit
    # samples with random.Random(hash(tuple(sorted((degree, Z))))).
    import hashlib

    smiles = sorted(
        m.rdkit_smiles
        for m in chematic.from_smiles("Br" + "C(Cl)" * 11 + "F").rdkit_stereoisomers()
    )
    assert len(smiles) == 1024
    digest = hashlib.sha256("\n".join(smiles).encode()).hexdigest()
    assert digest == "62988c34ffba9a112fab1409b27b88bdbdcba052cd0831c48dfada0a50f06f06"
