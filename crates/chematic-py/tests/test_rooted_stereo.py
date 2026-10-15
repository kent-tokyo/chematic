"""Stereo-checked rooting is opt-in; atom correspondences use original indices."""
import pytest
import chematic


@pytest.mark.parametrize("source", [
    "[CH:1]1=[CH:2][C@H:3]2[O:4][C@@H:5]1[CH2:6][CH2:7]2",
    "[CH:1]1=[CH:2][C@@H:3]2[O:4][C@H:5]1[CH2:6][CH2:7]2",
    "[13CH3:1][C@H:2]([OH:3])[F:4]",
    "[F:1][C@H:2]([OH:3])[CH3:4].[Cl-:5]",
    "[F:1]/[CH:2]=[CH:3]/[Cl:4]",
    "[N:1]->[Cu:2]",
])
def test_all_roots_preserve_identity_and_return_original_indices(source):
    mol = chematic.from_smiles(source)
    baseline = mol.smiles
    n = mol.atom_count
    for root in range(n):
        text, order = mol.smiles_rooted_checked(root)
        assert order[0] == root
        assert sorted(order) == list(range(n))
        assert chematic.from_smiles(text).smiles == baseline
        assert mol.smiles == baseline


@pytest.mark.parametrize("source", ["F[Pt@SP1](Cl)(Br)I", "FC=[C@AL1]=C(Cl)Br", "C[S@](=O)CC", "C~C"])
def test_other_stereo_and_query_bonds_refuse(source):
    with pytest.raises(ValueError):
        chematic.from_smiles(source).smiles_rooted_checked(0)


@pytest.mark.parametrize("root", [2, 100, 2**40])
def test_out_of_range_roots_refuse(root):
    with pytest.raises(ValueError, match="root atom"):
        chematic.from_smiles("CC").smiles_rooted_checked(root)


def test_negative_roots_are_not_wrapped():
    with pytest.raises(OverflowError):
        chematic.from_smiles("CC").smiles_rooted_checked(-1)


def test_atom_cap_is_explicit():
    assert len(chematic.from_smiles("C" * 256).smiles_rooted_checked(0)[1]) == 256
    with pytest.raises(ValueError, match="atom limit"):
        chematic.from_smiles("C" * 257).smiles_rooted_checked(0)
