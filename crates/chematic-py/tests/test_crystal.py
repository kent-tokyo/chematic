"""Tests for the chematic-crystal Python bindings: Lattice, PeriodicStructure,
Site, PeriodicNeighbor, CifSymmetryStatus.

Light-gate scope per project policy (binding-level change): focused cases
exercising the Python surface, not a corpus run -- the parser/expansion
grammar itself is covered by chematic-mol's Rust-side fixture suite
(crates/chematic-mol/src/cif.rs, cif_symmetry.rs).
"""
import numpy as np
import pytest

from chematic import Lattice, PeriodicStructure, Site


NACL_P1_CIF = """\
data_NaCl
_cell_length_a 5.6402
_cell_length_b 5.6402
_cell_length_c 5.6402
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na1 Na 0.0 0.0 0.0
Cl1 Cl 0.5 0.5 0.5
"""

# Realistic C2/c (space group No. 15) symop list -- with expand_symmetry
# (the default), Ti1's expanded images collapse pairwise onto 4 distinct
# sites (Ti1 sits on a special position); with expand_symmetry=False, the
# single listed atom is returned as-is, the pre-expansion behavior. Mirrors
# chematic-mol's own cif.rs test fixture for this exact scenario.
SYMMETRIC_C2C_CIF = """\
data_synthetic_c2c
_cell_length_a 10.0
_cell_length_b 8.0
_cell_length_c 12.0
_cell_angle_alpha 90
_cell_angle_beta 105.0
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'C 2/c'
loop_
_space_group_symop_operation_xyz
'x, y, z'
'-x, y, -z+1/2'
'-x, -y, -z'
'x, -y, z+1/2'
'x+1/2, y+1/2, z'
'-x+1/2, y+1/2, -z+1/2'
'-x+1/2, -y+1/2, -z'
'x+1/2, -y+1/2, z+1/2'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Ti1 Ti 0.0 0.25 0.25
"""

DISORDERED_CIF = """\
data_disordered
_cell_length_a 3.0
_cell_length_b 3.0
_cell_length_c 3.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
M1A Fe 0.0 0.0 0.0 0.6
M1B Ni 0.0 0.0 0.0 0.4
"""

NACL_POSCAR = """\
NaCl test structure
1.0
   5.6400000000000000    0.0000000000000000    0.0000000000000000
   0.0000000000000000    5.6400000000000000    0.0000000000000000
   0.0000000000000000    0.0000000000000000    5.6400000000000000
Na Cl
1 1
Direct
  0.0 0.0 0.0
  0.5 0.5 0.5
"""

NACL_POSCAR_SELECTIVE = """\
NaCl test structure
1.0
   5.6400000000000000    0.0000000000000000    0.0000000000000000
   0.0000000000000000    5.6400000000000000    0.0000000000000000
   0.0000000000000000    0.0000000000000000    5.6400000000000000
Na Cl
1 1
Selective dynamics
Direct
  0.0 0.0 0.0 T T T
  0.5 0.5 0.5 F F F
"""


# ---------------------------------------------------------------------------
# CIF round trip
# ---------------------------------------------------------------------------

def test_cif_roundtrip_preserves_lattice_sites_labels():
    s1 = PeriodicStructure.from_cif(NACL_P1_CIF)
    assert s1.symmetry_status.is_p1
    text = s1.to_cif()
    s2 = PeriodicStructure.from_cif(text)

    assert s1.lattice.lengths == pytest.approx(s2.lattice.lengths, abs=1e-3)
    assert s1.lattice.angles_degrees == pytest.approx(s2.lattice.angles_degrees, abs=1e-2)
    assert s1.site_count() == s2.site_count() == 2

    sites1 = sorted((site.species[0][0], site.species[0][1], site.fractional, site.label) for site in s1.sites)
    sites2 = sorted((site.species[0][0], site.species[0][1], site.fractional, site.label) for site in s2.sites)
    for (el1, occ1, f1, label1), (el2, occ2, f2, label2) in zip(sites1, sites2):
        assert el1 == el2
        assert occ1 == pytest.approx(occ2, abs=1e-3)
        assert f1 == pytest.approx(f2, abs=1e-4)
        assert label1 == label2


def test_cif_disordered_site_not_collapsed():
    s = PeriodicStructure.from_cif(DISORDERED_CIF)
    assert s.site_count() == 1
    species = dict(s.sites[0].species)
    assert species.keys() == {"Fe", "Ni"}
    assert species["Fe"] == pytest.approx(0.6)
    assert species["Ni"] == pytest.approx(0.4)


def test_cif_declared_symmetry_expands_by_default():
    s = PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF)
    status = s.symmetry_status
    assert status is not None
    assert status.is_p1 is False
    assert status.is_expanded is True
    assert status.is_complete_cell is True
    assert status.space_group_name == "C 2/c"
    assert status.operation_count == 8
    assert status.asymmetric_site_count == 1
    assert status.expanded_site_count == 4
    assert s.site_count() == 4
    # A genuinely complete (expanded) cell writes back out freely.
    assert s.to_cif()


def test_cif_expand_symmetry_false_keeps_asymmetric_unit_only():
    s = PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF, expand_symmetry=False)
    status = s.symmetry_status
    assert status is not None
    assert status.is_p1 is False
    assert status.is_expanded is False
    assert status.is_complete_cell is False
    assert status.space_group_name == "C 2/c"
    assert status.operation_count == 8
    assert status.asymmetric_site_count is None
    assert status.expanded_site_count is None
    assert s.site_count() == 1
    # Writing back must not silently re-declare this asymmetric unit as P1.
    with pytest.raises(ValueError):
        s.to_cif()


def test_cif_expand_symmetry_true_is_same_as_default():
    default = PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF)
    explicit = PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF, expand_symmetry=True)
    assert default.site_count() == explicit.site_count() == 4


def test_wrap_into_cell_and_make_supercell_preserve_complete_cell_status():
    # Regression spot-check (PR #333's wrap_into_cell()/make_supercell()
    # combined with this PR's symmetry expansion): both transforms must
    # keep reporting a genuinely complete cell as complete.
    s = PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF)
    assert s.symmetry_status.is_complete_cell is True

    wrapped = s.wrap_into_cell()
    assert wrapped.symmetry_status.is_complete_cell is True
    assert wrapped.to_cif()

    supercell = s.make_supercell((2, 1, 1))
    assert supercell.symmetry_status.is_complete_cell is True
    assert supercell.site_count() == 8
    assert supercell.to_cif()

    # And the unexpanded (opted-out) case must still refuse after both
    # transforms too -- the asymmetric-unit problem doesn't go away.
    unexpanded = PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF, expand_symmetry=False)
    assert unexpanded.wrap_into_cell().symmetry_status.is_complete_cell is False
    assert unexpanded.make_supercell((2, 1, 1)).symmetry_status.is_complete_cell is False
    with pytest.raises(ValueError):
        unexpanded.wrap_into_cell().to_cif()
    with pytest.raises(ValueError):
        unexpanded.make_supercell((2, 1, 1)).to_cif()


# ---------------------------------------------------------------------------
# POSCAR round trip
# ---------------------------------------------------------------------------

def test_poscar_roundtrip_is_semantically_equal():
    # write_poscar always emits scale 1.0 with pre-scaled vectors, so a
    # byte-identical round trip isn't the right claim to test -- what must
    # hold is that re-parsing the written text yields the same lattice and
    # sites, not the same text.
    s1 = PeriodicStructure.from_poscar(NACL_POSCAR)
    s2 = PeriodicStructure.from_poscar(s1.to_poscar())

    assert s1.lattice.lengths == pytest.approx(s2.lattice.lengths, abs=1e-6)
    assert s1.site_count() == s2.site_count() == 2
    for site1, site2 in zip(s1.sites, s2.sites):
        assert site1.species == site2.species
        assert site1.fractional == pytest.approx(site2.fractional, abs=1e-6)


def test_poscar_roundtrip_preserves_extras():
    s = PeriodicStructure.from_poscar(NACL_POSCAR_SELECTIVE)
    out = s.to_poscar()
    lines = out.splitlines()
    assert lines[0] == "NaCl test structure"
    assert "Selective dynamics" in out
    # Per-site selective-dynamics T/T/T and F/F/F flags, in the order the
    # sites were originally read (Na then Cl) -- the two lines right after
    # the "Direct" coordinate-mode line.
    direct_idx = lines.index("Direct")
    assert lines[direct_idx + 1].split()[-3:] == ["T", "T", "T"]
    assert lines[direct_idx + 2].split()[-3:] == ["F", "F", "F"]


def test_poscar_rejects_disordered_structure():
    lattice = Lattice.cubic(3.0)
    site = Site([("Fe", 0.6), ("Ni", 0.4)], (0.0, 0.0, 0.0))
    s = PeriodicStructure(lattice, [site])
    with pytest.raises(ValueError):
        s.to_poscar()


# ---------------------------------------------------------------------------
# Validation / fail-closed
# ---------------------------------------------------------------------------

def test_invalid_occupancy_sum_raises():
    with pytest.raises(ValueError):
        Site([("Fe", 0.7), ("Ni", 0.5)], (0.0, 0.0, 0.0))


def test_nan_infinity_rejected():
    with pytest.raises(ValueError):
        Lattice.cubic(float("nan"))
    with pytest.raises(ValueError):
        Site([("Na", 1.0)], (float("inf"), 0.0, 0.0))


def test_invalid_degenerate_lattice_raises():
    # Three linearly dependent rows -> zero-volume (singular) matrix.
    with pytest.raises(ValueError):
        Lattice.from_matrix([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]])


# ---------------------------------------------------------------------------
# NumPy interop
# ---------------------------------------------------------------------------

def test_lattice_matrix_numpy_shape_dtype():
    lat = Lattice.cubic(4.0)
    m = lat.matrix
    assert isinstance(m, np.ndarray)
    assert m.shape == (3, 3)
    assert m.dtype == np.float64
    assert lat.volume == pytest.approx(64.0)


def test_cartesian_positions_numpy_shape_dtype():
    s = PeriodicStructure.from_cif(NACL_P1_CIF)
    pos = s.cartesian_positions()
    assert isinstance(pos, np.ndarray)
    assert pos.shape == (2, 3)
    assert pos.dtype == np.float64
    assert pos[1] == pytest.approx([2.8201, 2.8201, 2.8201], abs=1e-3)


# ---------------------------------------------------------------------------
# neighbors / make_supercell / wrap_into_cell
# ---------------------------------------------------------------------------

def test_neighbors_and_supercell_simple_cubic():
    lattice = Lattice.cubic(3.0)
    site = Site([("Ar", 1.0)], (0.0, 0.0, 0.0))
    s = PeriodicStructure(lattice, [site])

    neighbors = s.neighbors(cutoff=3.0)
    assert len(neighbors) == 6
    assert all(n.distance == pytest.approx(3.0) for n in neighbors)

    supercell = s.make_supercell((2, 2, 2))
    assert supercell.site_count() == 8
    assert supercell.lattice.volume == pytest.approx(s.lattice.volume * 8.0)


def test_wrap_into_cell_reduces_out_of_range_fractional():
    lattice = Lattice.cubic(4.0)
    site = Site([("Fe", 1.0)], (1.25, -0.5, 3.0), label="Fe1")
    s = PeriodicStructure(lattice, [site])

    wrapped = s.wrap_into_cell()
    fx, fy, fz = wrapped.sites[0].fractional
    assert fx == pytest.approx(0.25, abs=1e-9)
    assert fy == pytest.approx(0.5, abs=1e-9)
    assert fz == pytest.approx(0.0, abs=1e-9)
    # original left unmodified
    assert s.sites[0].fractional == pytest.approx((1.25, -0.5, 3.0))


# ---------------------------------------------------------------------------
# Non-orthogonal (triclinic) lattice
# ---------------------------------------------------------------------------

def test_triclinic_lattice_frac_cart_roundtrip_and_neighbors():
    lat = Lattice.from_parameters(4.0, 5.0, 6.0, 80.0, 85.0, 75.0)
    alpha, beta, gamma = lat.angles_degrees
    assert alpha != pytest.approx(90.0)
    assert beta != pytest.approx(90.0)
    assert gamma != pytest.approx(90.0)

    point = (0.3, 0.4, 0.5)
    cart = lat.frac_to_cart(point)
    back = lat.cart_to_frac(cart)
    assert back == pytest.approx(point, abs=1e-9)

    site = Site([("C", 1.0)], (0.0, 0.0, 0.0))
    s = PeriodicStructure(lat, [site])
    neighbors = s.neighbors(cutoff=6.5)
    assert len(neighbors) > 0
    assert all(0.0 < n.distance <= 6.5 for n in neighbors)


# ---------------------------------------------------------------------------
# formula / repr
# ---------------------------------------------------------------------------

def test_formula_property_unreduced_and_occupancy_weighted():
    s = PeriodicStructure.from_cif(NACL_P1_CIF)
    assert s.formula == "ClNa"

    supercell = s.make_supercell((2, 2, 2))
    assert supercell.formula == "Cl8Na8"

    disordered = PeriodicStructure.from_cif(DISORDERED_CIF)
    assert disordered.formula == "Fe0.6Ni0.4"


def test_repr_smoke():
    s = PeriodicStructure.from_cif(NACL_P1_CIF)
    assert "PeriodicStructure" in repr(s)
    assert "Lattice" in repr(s.lattice)
    assert "Site" in repr(s.sites[0])
