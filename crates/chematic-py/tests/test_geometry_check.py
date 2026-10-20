"""Original geometric controls for Open Babel #2651 (no imported corpus)."""

import copy
import math

import pytest

import chematic


def test_finite_collapsed_and_stretched_coordinates_do_not_silently_pass():
    mol = chematic.from_smiles("CC")
    for length in (0.0, 0.1, 10.0):
        result = mol.check_geometry([[0, 0, 0], [length, 0, 0]])
        assert result["bonds_checked"] == 1
        assert result["skipped_bonds"] == []
        (violation,) = result["bond_violations"]
        assert list(violation["atoms"]) == [0, 1]
        assert violation["bond"] == 0
        assert violation["distance"] == length
        assert violation["reference_length"] == pytest.approx(1.54)
        assert violation["relative_deviation"] == pytest.approx(
            abs(length - violation["reference_length"]) / violation["reference_length"]
        )
        assert result["nonbonded_clashes"] == []
        assert "valid" not in result
        assert "sound" not in result
        assert "is_fully_satisfied" not in result["stereo"]


def test_nonbonded_clashes_preserve_original_atom_indices_and_distances():
    rows = [[0, 0, 0], [1.54, 0, 0], [0, 0, 0], [0, 0.6, 0.8]]
    result = chematic.from_smiles("CCC.C").check_geometry(rows)
    clashes = result["nonbonded_clashes"]
    assert [list(c["atoms"]) for c in clashes] == [[0, 2], [0, 3], [2, 3]]
    assert [c["distance"] for c in clashes] == [0.0, 1.0, 1.0]
    assert result["bond_violations"] == []
    for c in clashes:
        a, b = c["atoms"]
        assert c["distance"] == math.dist(rows[a], rows[b])


def test_caller_thresholds_are_retained_and_use_strict_inequalities():
    mol = chematic.from_smiles("C.C")
    rows = [[0, 0, 0], [0, 0, 1.2]]
    assert mol.check_geometry(rows)["nonbonded_clashes"] == []
    r = mol.check_geometry(
        rows, nonbonded_min_distance=1.3, bond_relative_tolerance=0.2
    )
    assert r["nonbonded_min_distance"] == 1.3
    assert r["bond_relative_tolerance"] == 0.2
    assert len(r["nonbonded_clashes"]) == 1
    assert r["max_relative_bond_deviation"] is None


@pytest.mark.parametrize("smiles", ["*C", "N->[Cu+2]"])
def test_unsupported_reference_is_explicit_and_not_a_nonbonded_pair(smiles):
    r = chematic.from_smiles(smiles).check_geometry([[0, 0, 0], [0, 0, 0]])
    assert r["skipped_bonds"] == [0]
    assert r["bonds_checked"] == 0
    assert r["max_relative_bond_deviation"] is None
    assert r["nonbonded_clashes"] == []


TETRAHEDRON = [[1, 1, 1], [0, 0, 0], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]


def test_mirror_stereo_violation_is_distinct_from_unevaluable_planar_stereo():
    mol = chematic.from_smiles("F[C@](Cl)(Br)I")
    original = mol.check_geometry(TETRAHEDRON)["stereo"]
    mirrored = mol.check_geometry([[-x, y, z] for x, y, z in TETRAHEDRON])["stereo"]
    assert original["n_declared"] == mirrored["n_declared"] == 1
    assert original["n_violations"] + mirrored["n_violations"] == 1
    assert original["n_satisfied"] + mirrored["n_satisfied"] == 1
    planar = mol.check_geometry([[x, y, 0] for x, y, z in TETRAHEDRON])["stereo"]
    assert planar["n_violations"] == 0
    assert planar["n_unevaluable"] == 1
    (row,) = planar["tetrahedral"]
    assert row["atom"] == 1
    assert row["status"] == "unevaluable"
    assert row["rejection_reason"] == "degenerate_geometry"
    unspecified = chematic.from_smiles("FC(Cl)(Br)I").check_geometry(TETRAHEDRON)
    assert unspecified["stereo"]["n_declared"] == 0


def test_ez_reports_and_nontetrahedral_coverage():
    coords = [[-1, 1, 0], [0, 0, 0], [1.34, 0, 0], [2.34, -1, 0]]
    e = chematic.from_smiles("C/C=C/C").check_geometry(coords)["stereo"]
    z = chematic.from_smiles("C/C=C\\C").check_geometry(coords)["stereo"]
    assert e["n_declared"] == z["n_declared"] == 1
    assert e["n_violations"] + z["n_violations"] == 1
    r = chematic.from_smiles("Cl[Pt@SP1](Cl)(N)N").check_geometry([[0, 0, 0]] * 5)
    assert r["stereo"]["n_unevaluable"] == 1
    assert (
        r["stereo"]["tetrahedral"][0]["rejection_reason"] == "unsupported_coordination"
    )


def test_enhanced_group_semantics_are_not_silently_evaluated_as_absolute():
    mol = chematic.from_smiles("F[C@](Cl)(Br)I")
    grouped = mol.with_stereo_groups([("or", 1, [1], [])])
    before = grouped.stereo_groups
    result = grouped.check_geometry(TETRAHEDRON)
    assert result["stereo_groups_not_evaluated"] == 1
    assert grouped.stereo_groups == before


def test_read_only_and_rigid_transform_invariance():
    mol = chematic.from_smiles("F[C@](Cl)(Br)I")
    rows = copy.deepcopy(TETRAHEDRON)
    mol.SetProp("name", "original")
    before = mol.smiles
    a = mol.check_geometry(rows)
    moved = [[-y + 20, x - 10, z + 5] for x, y, z in rows]
    assert a == mol.check_geometry(moved)
    assert rows == TETRAHEDRON
    assert mol.smiles == before
    assert mol.GetProp("name") == "original"


def test_numpy_and_generator_inputs():
    np = pytest.importorskip("numpy")
    mol = chematic.from_smiles("CC")
    rows = np.array([[0, 0, 0], [1.54, 0, 0]], dtype=float)
    assert mol.check_geometry(rows) == mol.check_geometry(iter(rows.tolist()))


def test_excess_input_consumption_is_bounded():
    seen = []

    def infinite():
        while True:
            seen.append(1)
            yield [0, 0, 0]

    with pytest.raises(ValueError, match="one entry per atom"):
        chematic.from_smiles("CC").check_geometry(infinite())
    assert len(seen) == 3


def test_source_exception_propagates():
    def source():
        yield [0, 0, 0]
        raise RuntimeError("source failed")

    with pytest.raises(RuntimeError, match="source failed"):
        chematic.from_smiles("CC").check_geometry(source())


@pytest.mark.parametrize("rows", [[], [[0, 0, 0]], [[0, 0, 0]] * 3])
def test_bad_coordinate_counts(rows):
    with pytest.raises(ValueError):
        chematic.from_smiles("CC").check_geometry(rows)


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf, 1e6 + 1, -1e6 - 1])
def test_nonfinite_and_unbounded_coordinates(value):
    with pytest.raises(ValueError):
        chematic.from_smiles("CC").check_geometry([[0, 0, 0], [value, 0, 0]])


@pytest.mark.parametrize("tolerance", [-1, 1, math.nan, math.inf])
def test_invalid_bond_tolerance(tolerance):
    with pytest.raises(ValueError, match="bond_relative_tolerance"):
        chematic.from_smiles("C").check_geometry(
            [[0, 0, 0]], bond_relative_tolerance=tolerance
        )


@pytest.mark.parametrize("cutoff", [0, -1, math.nan, math.inf])
def test_invalid_clash_cutoff(cutoff):
    with pytest.raises(ValueError, match="nonbonded_min_distance"):
        chematic.from_smiles("C").check_geometry(
            [[0, 0, 0]], nonbonded_min_distance=cutoff
        )


@pytest.mark.parametrize("row", [[0, 0], [0, 0, 0, 0], None, ["x", 0, 0]])
def test_malformed_coordinate_rows(row):
    with pytest.raises((TypeError, ValueError)):
        chematic.from_smiles("C").check_geometry([row])


def test_atom_limit_is_checked_before_consuming_coordinates():
    def forbidden():
        raise AssertionError("must not iterate")
        yield

    mol = chematic.from_smiles("C" * 2049)
    with pytest.raises(ValueError, match="2048"):
        mol.check_geometry(forbidden())


def test_single_atom_has_no_claim_of_stereo_or_bond_validation():
    r = chematic.from_smiles("[He]").check_geometry([[0, 0, 0]])
    assert r["atom_count"] == 1
    assert r["bonds_checked"] == 0
    assert r["max_relative_bond_deviation"] is None
    assert r["stereo"]["n_declared"] == 0
    assert r["stereo_groups_not_evaluated"] == 0


def test_tiny_finite_distances_are_not_rounded_to_coincident_atoms():
    coords = [[0, 0, 0], [1e-200, 0, 0]]
    bonded = chematic.from_smiles("CC").check_geometry(coords)
    assert bonded["bond_violations"][0]["distance"] == 1e-200
    separate = chematic.from_smiles("C.C").check_geometry(
        coords, nonbonded_min_distance=1e-250
    )
    assert separate["nonbonded_clashes"] == []
