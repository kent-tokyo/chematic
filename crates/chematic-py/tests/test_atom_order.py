"""Issue #650: SMILES atom output order and fragment source indices."""

import chematic


def _atoms(mol):
    # (symbol, atomic number, charge, aromatic) per atom
    return [row[:4] for row in mol.atom_table]


def test_smiles_with_atom_order_matches_smiles_and_reparse():
    for smi in ["OCC", "c1ccccc1O", "N[C@@H](C)C(=O)O", "[Na+].[Cl-]", "CC(=O)Nc1ccc(O)cc1", "F/C=C/Cl"]:
        mol = chematic.from_smiles(smi)
        out, order = mol.smiles_with_atom_order()
        assert out == mol.smiles
        atoms, again = _atoms(mol), _atoms(chematic.from_smiles(out))
        assert sorted(order) == list(range(len(atoms)))
        assert again == [atoms[source] for source in order]


def test_visit_order_example():
    assert chematic.from_smiles("OCC").smiles_with_atom_order() == ("C(C)O", [1, 2, 0])


def test_connected_components_with_atom_indices():
    mol = chematic.from_smiles("O.CC.N")
    parts = mol.connected_components_with_atom_indices()
    assert [source for _, source in parts] == [[0], [1, 2], [3]]
    assert [frag.smiles for frag, _ in parts] == [m.smiles for m in mol.connected_components()]
