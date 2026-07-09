"""Tests for chematic.o3a_align (O3A-style 3D atom correspondence + alignment)."""
import pytest
import chematic


def test_self_rotation_recovers_identity_and_zero_rmsd():
    """Aligning a molecule onto a rotated+translated copy of itself should
    recover the identity correspondence with RMSD ~0."""
    m = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    coords1 = m.generate_3d()
    coords2 = [[-p[1] + 5.0, p[0] - 3.0, p[2] + 2.0] for p in coords1]

    pairs, aligned, rmsd, score = chematic.o3a_align(m, coords1, m, coords2)
    assert len(pairs) == m.heavy_atoms
    assert all(i == j for i, j in pairs)
    assert rmsd < 1e-6
    assert score > 0


def test_pairs_are_injective():
    """No mol2 atom index should appear twice in the correspondence."""
    m = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    coords1 = m.generate_3d()
    coords2 = [[-p[1] + 5.0, p[0] - 3.0, p[2] + 2.0] for p in coords1]

    pairs, _aligned, _rmsd, _score = chematic.o3a_align(m, coords1, m, coords2)
    js = [j for _i, j in pairs]
    assert len(set(js)) == len(js)


def test_aligned_coords_cover_full_mol2():
    """aligned_coords2 should have one entry per mol2 heavy atom, not just
    the matched subset."""
    benzene = chematic.from_smiles("c1ccccc1")
    toluene = chematic.from_smiles("Cc1ccccc1")
    bc = benzene.generate_3d()
    tc = toluene.generate_3d()

    _pairs, aligned, _rmsd, _score = chematic.o3a_align(benzene, bc, toluene, tc)
    assert len(aligned) == toluene.heavy_atoms


def test_shared_scaffold_scores_higher_than_unrelated():
    """Toluene shares benzene's aromatic ring; cyclohexane (saturated, a
    different MMFF94 atom type and non-planar geometry) is a genuine
    negative control."""
    benzene = chematic.from_smiles("c1ccccc1")
    toluene = chematic.from_smiles("Cc1ccccc1")
    cyclohexane = chematic.from_smiles("C1CCCCC1")

    bc = benzene.generate_3d()
    tc = toluene.generate_3d()
    cc = cyclohexane.generate_3d()

    _p1, _a1, _r1, related_score = chematic.o3a_align(benzene, bc, toluene, tc)
    _p2, _a2, _r2, unrelated_score = chematic.o3a_align(benzene, bc, cyclohexane, cc)
    assert related_score > unrelated_score


def test_coordinate_length_mismatch_raises():
    m = chematic.from_smiles("CO")
    with pytest.raises(ValueError):
        chematic.o3a_align(m, [[0.0, 0.0, 0.0]], m, [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]])
