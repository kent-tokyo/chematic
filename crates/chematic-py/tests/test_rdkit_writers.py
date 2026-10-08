"""RDKit-compatible writers against RDKit 2026.03.1 (MolFromSmiles input)."""

import chematic


def test_smarts_cx_smarts_pdb():
    m = chematic.from_smiles("N->[Cu]")
    assert m.rdkit_smarts() == "[#7]->[Cu]"
    assert m.rdkit_cx_smarts() == "[#7]-[Cu] |C:0.0|"
    block = chematic.from_smiles("CC(=O)[O-]").rdkit_pdb_block()
    assert block.endswith("CONECT    1    2\nCONECT    2    3    3    4\nEND\n")
    assert "O1-\n" in block


def test_murcko_chiral_centers_counts():
    m = chematic.from_smiles("O=C1CC[C@H](C)N1CCc1ccccc1")
    assert m.rdkit_murcko_scaffold() == "O=C1CCCN1CCc1ccccc1"
    assert m.rdkit_chiral_centers() == [(4, "S")]
    assert chematic.from_smiles("Cn1cccc1CC1CCC1").rdkit_murcko_scaffold() == "c1c[nH]c(CC2CCC2)c1"
    assert chematic.from_smiles("CC(F)C(Cl)Br").rdkit_stereoisomer_count() == 4
    smiles = chematic.from_smiles("Cc1n[s+]([O-])nc1N").rdkit_stereoisomer_smiles()
    assert smiles == ["Cc1n[s@+]([O-])nc1N", "Cc1n[s@@+]([O-])nc1N"]


def test_rdkit_mol_hash():
    # rdMolHash.MolHash(Chem.MolFromSmiles(s), f[, useCXSmiles]), RDKit 2026.03.1.
    m = chematic.from_smiles("Cc1ccccc1CC(=O)O")
    assert m.rdkit_mol_hash("ExtendedMurcko") == "*c1ccccc1*"
    assert m.rdkit_mol_hash("extendedmurcko") == "*c1ccccc1*"
    assert m.rdkit_mol_hash("Regioisomer") == "*C.*CC(=O)O.c1ccccc1"
    assert chematic.from_smiles("Oc1ccccn1").rdkit_mol_hash("HetAtomTautomer") == "[O][C]1[CH][CH][CH][CH][N]1_1_0"
    assert chematic.from_smiles("[NH3+]CC([O-])=O").rdkit_mol_hash("MolFormula") == "C2H5NO2"
    assert chematic.from_smiles("N->[Cu]").rdkit_mol_hash("CanonicalSmiles", use_cx_smiles=True) == "[NH3][Cu] |C:0.0|"
    try:
        m.rdkit_mol_hash("NoSuchHash")
    except ValueError as e:
        assert "ExtendedMurcko" in str(e)
    else:
        raise AssertionError("unknown hash function accepted")
