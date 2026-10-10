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


def test_rdkit_cx_random_fragment_smiles():
    # RDKit 2026.03.1: MolToCXSmiles, MolToRandomSmilesVect, MolFragmentToSmiles.
    assert chematic.from_smiles("C[CH2]").rdkit_cx_smiles() == "[CH2]C |^1:0|"
    m = chematic.from_smiles("OC(=O)c1ccccc1")
    assert m.rdkit_fragment_smiles([0, 1, 2]) == "O=CO"
    assert m.rdkit_num_atoms() == 9
    assert chematic.from_smiles("CCO").rdkit_random_smiles(5, 42) == [
        "CCO", "OCC", "CCO", "CCO", "CCO"]


def test_rdkit_problems_matrices_fingerprints():
    # RDKit 2026.03.1: DetectChemistryProblems (sanitize=False), distance
    # matrices, sparse count fingerprints.
    assert chematic.rdkit_detect_chemistry_problems("c1cccc1") == [
        ("KekulizeException", [0, 1, 2, 3, 4])]
    assert chematic.rdkit_detect_chemistry_problems("F(C)C") == [("AtomValenceException", [0])]
    m = chematic.from_smiles("CC(=O)O")
    assert m.rdkit_chemistry_problems() == []
    assert m.rdkit_distance_matrix()[0] == [0.0, 1.0, 2.0, 2.0]
    assert m.rdkit_distance_matrix(use_bo=True)[1][2] == 0.5
    d = m.rdkit_distance_matrix_3d([[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 0, 1]])
    assert d[0][3] == 2 ** 0.5
    assert m.rdkit_morgan_sparse_counts(2)[2246728737] == 1
    assert len(m.rdkit_legacy_torsion_counts()) == 0
    assert sum(m.rdkit_atom_pair_sparse_counts().values()) == 6


def test_rdkit_reaction_smarts_writer_and_grouping():
    # RDKit 2026.03.1 ReactionToSmarts(ReactionFromSmarts(s)) / MolToSmarts.
    s = "[C@@H:1]([NH2:2])([#6:3])[C:4]=[O:5]>>[C@@H:1]([NH:2]C(C)=O)([#6:3])[C:4]=[O:5]"
    assert chematic.rdkit_reaction_to_smarts(s) == (
        "[N&H2:2][C@&H1:1]([#6:3])[C:4]=[O:5]>>[N&H1:2]([C@&H1:1]([#6:3])[C:4]=[O:5])C(C)=O")
    lactam = "([C:1](=[O:2])[OH].[NH2:3])>>[C:1](=[O:2])[N:3]"
    assert chematic.rdkit_reaction_to_smarts(lactam) == (
        "([C:1](=[O:2])[O&H1].[N&H2:3])>>[C:1](=[O:2])[N:3]")
    assert chematic.rdkit_smarts_to_smarts("[N;H2,H1;!$(NC=O)]") == "[N;H2,H1;!$(NC=O)]"
    r = chematic.run_smirks_checked(lactam, [chematic.from_smiles("NCCCC(=O)O")], rdkit_compat=True)
    assert r["status"] == "products"
    assert sorted({p.rdkit_smiles for ps in r["products"] for p in ps}) == ["O=C1CCCN1"]


def test_rdkit_reaction_profile_maps_and_rejected_products():
    # RDKit copies unmatched reactant atoms with their atom maps.
    r = chematic.run_smirks_checked("[C:1](=[O:2])[OH]>>[C:1](=[O:2])N",
                                    [chematic.from_smiles("[CH3:7]C(=O)O")], rdkit_compat=True)
    assert [p.smiles for ps in r["products"] for p in ps] == ["NC(=O)[CH3:7]"]
    # RunReactants returns product sets sanitization rejects; chematic
    # reports them as rejected_products.
    r = chematic.run_smirks_checked("[CH2:1][OH:2]>>[CH:1]=[O:2]",
                                    [chematic.from_smiles("[OH+3](C)CO")], rdkit_compat=True)
    assert r["reason"] == "product_valence" and len(r["rejected_products"]) == 1
    diagnostic = r["rejection_diagnostics"][0][0]
    assert diagnostic["reason"] == "valence"
    assert diagnostic["accepted"] is False
    assert all(
        {"atom_index", "atom_map", "element", "explicit_hydrogen_count",
         "observed_valence", "max_allowed_valence"} <= atom.keys()
        for atom in diagnostic["atoms"]
    )
    # A * atom has atomic number 0: [#6] does not match it.
    r = chematic.run_smirks_checked("[C:1](=[O:2])([#6:3])[#6:4]>>[C:1]([OH:2])([#6:3])[#6:4]",
                                    [chematic.from_smiles("*C(=O)*")], rdkit_compat=True)
    assert r["status"] == "no_match"
