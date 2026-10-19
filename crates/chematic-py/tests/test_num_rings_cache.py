"""Ring-count memoization follows the immutable graph, never mutable properties."""
import pickle

import pytest
import chematic as ch


@pytest.mark.parametrize('smiles,count', [
    ('CCO', 0), ('C1CC1', 1), ('C1CC1.C1CCC1', 2),
    ('C12C3C4C1C5C2C3C45', 6), ('C' * 129 + 'C1CC1', 1),
])
def test_counts_survive_property_edits_and_pickle(smiles, count):
    mol = ch.from_smiles(smiles)
    assert mol.num_rings == count
    mol.SetProp('name', 'cached')
    mol.SetIntProp('row', 42)
    mol.ClearProp('name')
    assert mol.num_rings == count
    restored = pickle.loads(pickle.dumps(mol))
    assert restored.num_rings == count
    assert restored.num_rings == count


def test_editable_snapshots_keep_their_own_graph_counts():
    original = ch.from_smiles('C1CC1')
    assert original.num_rings == 1
    edit = ch.RWMol(original)
    edit.RemoveBond(0, 1)
    opened = edit.GetMol()
    assert opened.num_rings == 0
    edit.AddBond(0, 1, 'SINGLE')
    closed = edit.GetMol()
    assert closed.num_rings == 1
    assert opened.num_rings == 0
    assert original.num_rings == 1


def test_graph_transform_does_not_inherit_the_source_count():
    original = ch.from_smiles('C1CC1.CCCCC')
    assert original.num_rings == 1
    assert original.largest_fragment().num_rings == 0
    assert original.num_rings == 1
