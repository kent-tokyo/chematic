"""XSF geometry projection contracts over existing native crystal objects."""

import math

import numpy as np
import pytest
from chematic import Lattice, PeriodicStructure, Site, write_xsf_checked
from .test_crystal import NACL_P1_CIF, SYMMETRIC_C2C_CIF, NACL_POSCAR_SELECTIVE


def direct(sites=None):
    lattice = Lattice.from_matrix([[2.0, 0.1, 0.0], [0.0, 3.0, 0.2], [0.3, 0.0, 4.0]])
    return PeriodicStructure(
        lattice,
        (
            sites
            if sites is not None
            else [Site([("Cl", 1.0)], (1.25, -0.5, 2.0), "identity")]
        ),
    )


def test_basis_cartesian_order_no_wrap_and_source_unchanged():
    s = direct(
        [
            Site(
                [("Cl", 1.0)],
                (1.25, -0.5, 2.0),
                "label" + chr(10) + "not an XSF command",
            ),
            Site([("Na", 1.0)], (0.0, 0.0, 0.0), "Na_label"),
        ]
    )
    source = s.fractional_positions().copy()
    text = write_xsf_checked(s, allow_metadata_loss=True)
    lines = text.splitlines()
    assert lines[:2] == ["CRYSTAL", "PRIMVEC"]
    assert (
        np.asarray([[float(x) for x in row.split()] for row in lines[2:5]]).tobytes()
        == s.lattice.matrix.tobytes()
    )
    assert lines[5:7] == ["PRIMCOORD", "2 1"]
    assert [row.split()[0] for row in lines[7:]] == ["17", "11"]
    assert [float(x) for x in lines[7].split()[1:]] == pytest.approx([3.1, -1.375, 7.9])
    assert np.array_equal(s.fractional_positions(), source)
    assert "not an XSF" not in text and s.sites[0].label.startswith("label\n")
    assert write_xsf_checked(s, allow_metadata_loss=True) == text


def test_geometry_omission_requires_explicit_boolean():
    with pytest.raises(ValueError, match="geometry only"):
        write_xsf_checked(direct())
    for value in (1, "yes", None, np.bool_(True)):
        with pytest.raises(TypeError):
            write_xsf_checked(direct(), allow_metadata_loss=value)


@pytest.mark.parametrize("source", [None, "CIF text", {}, Lattice.cubic(1.0)])
def test_wrong_structure_type(source):
    with pytest.raises(TypeError):
        write_xsf_checked(source, allow_metadata_loss=True)


@pytest.mark.parametrize(
    "species",
    [
        [("C", 0.5)],
        [("C", 0.999999999)],
        [("C", 1.0), ("H", 0.0)],
        [("C", 0.5), ("N", 0.5)],
    ],
)
def test_partial_mixed_and_zero_occupancy_components_not_collapsed(species):
    s = direct([Site(species, (0.0, 0.0, 0.0))])
    with pytest.raises(ValueError, match="fully occupied"):
        write_xsf_checked(s, allow_metadata_loss=True)


def test_p1_and_expanded_cif_are_writable():
    p1 = PeriodicStructure.from_cif(NACL_P1_CIF)
    assert p1.symmetry_status.is_complete_cell
    assert "2 1\n" in write_xsf_checked(p1, allow_metadata_loss=True)
    expanded = PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF)
    assert expanded.symmetry_status.is_expanded
    assert f"{expanded.site_count()} 1\n" in write_xsf_checked(
        expanded, allow_metadata_loss=True
    )


@pytest.mark.parametrize(
    "transform",
    [lambda s: s, lambda s: s.wrap_into_cell(), lambda s: s.make_supercell((2, 1, 1))],
)
def test_unexpanded_symmetry_remains_refused_after_transforms(transform):
    s = transform(PeriodicStructure.from_cif(SYMMETRIC_C2C_CIF, expand_symmetry=False))
    assert not s.symmetry_status.is_complete_cell
    with pytest.raises(ValueError, match="unexpanded"):
        write_xsf_checked(s, allow_metadata_loss=True)


def test_poscar_constraints_source_kept_with_geometry_projection():
    s = PeriodicStructure.from_poscar(NACL_POSCAR_SELECTIVE)
    before = s.to_poscar()
    out = write_xsf_checked(s, allow_metadata_loss=True)
    assert out.startswith("CRYSTAL\n") and "Selective dynamics" not in out
    assert s.to_poscar() == before


def test_finite_native_input_can_overflow_cartesian_projection():
    s = PeriodicStructure(Lattice.cubic(1e100), [Site([("C", 1.0)], (1e250, 0.0, 0.0))])
    assert math.isfinite(s.sites[0].fractional[0])
    with pytest.raises(ValueError, match="finite"):
        write_xsf_checked(s, allow_metadata_loss=True)


def test_empty_and_site_cap():
    with pytest.raises(ValueError, match="1..50000"):
        write_xsf_checked(direct([]), allow_metadata_loss=True)
    site = Site([("C", 1.0)], (0.0, 0.0, 0.0))
    lattice = Lattice.cubic(3.0)
    with pytest.raises(ValueError, match="1..50000"):
        write_xsf_checked(
            PeriodicStructure(lattice, [site] * 50001), allow_metadata_loss=True
        )
    out = write_xsf_checked(
        PeriodicStructure(lattice, [site] * 50000), allow_metadata_loss=True
    )
    assert out.splitlines()[6] == "50000 1" and len(out.splitlines()) == 50007
