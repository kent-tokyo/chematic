import chematic


def test_recursive_smarts_chirality_uses_mapped_neighbor_parity():
    trans = chematic.from_smiles("C1[C@H]2[C@@H](C2)C1")
    cis = chematic.from_smiles("C1[C@H]2[C@H](C2)C1")
    query = "[$(C1[C@H]2[C@@H](C2)C1)]"

    assert trans.has_substructure(query, use_chirality=True)
    assert not cis.has_substructure(query, use_chirality=True)
    assert trans.find_matches(query, use_chirality=True)
    assert not cis.find_matches(query, use_chirality=True)


def test_recursive_smarts_chirality_can_still_be_ignored_explicitly():
    cis = chematic.from_smiles("C1[C@H]2[C@H](C2)C1")
    query = "[$(C1[C@H]2[C@@H](C2)C1)]"

    assert cis.has_substructure(query)
    assert cis.find_matches(query)
