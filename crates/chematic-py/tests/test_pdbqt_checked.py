"""Output field integrity, explicit counts and compatibility for Open Babel #2101."""
import math

import pytest
import chematic


def atom_lines(text):
    return [line for line in text.splitlines() if line.startswith("ATOM")]


def test_accepted_output_preserves_legacy_bytes_and_explicit_hydrogens():
    mol = chematic.from_smiles("[H]OC(Cl)Br")
    coords = [(float(i), -2.345, 6.789) for i in range(mol.atom_count)]
    charges = [-0.123] * mol.atom_count
    text = chematic.write_pdbqt_checked(mol, coords, charges, "explicit H")
    assert text == mol.to_pdbqt(coords, charges, "explicit H")
    for i, line in enumerate(atom_lines(text)):
        assert len(line) == 79
        assert int(line[6:11]) == i + 1
        assert tuple(float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]) == coords[i]
        assert float(line[70:76]) == charges[i]
        assert line[77:79].strip() == ["HD", "OA", "C", "Cl", "Br"][i]


@pytest.mark.parametrize("axis", range(3))
@pytest.mark.parametrize("bad", [2261608.729, -1000.0, 10000.0, 9999.9999, -999.9999, math.inf, -math.inf, math.nan])
def test_coordinates_refuse_overflow_and_nonfinite(axis, bad):
    coords = [[0.0] * 3 for _ in range(2)]
    coords[1][axis] = bad
    with pytest.raises(ValueError, match=f"atom 1: {'xyz'[axis]}"):
        chematic.write_pdbqt_checked(chematic.from_smiles("CC"), coords, [0, 0])


@pytest.mark.parametrize("bad", [100, -10, 99.9999, -9.9999, math.inf, -math.inf, math.nan])
def test_partial_charge_refuses_overflow_and_nonfinite(bad):
    with pytest.raises(ValueError, match="atom 0: charge"):
        chematic.write_pdbqt_checked(chematic.from_smiles("C"), [(0, 0, 0)], [bad])


@pytest.mark.parametrize("count", [0, 1, 2, 4])
@pytest.mark.parametrize("field", ["coordinates", "charges"])
def test_counts_do_not_default_missing_or_discard_excess_values(count, field):
    mol = chematic.from_smiles("[H]O[H]")
    coords = [(0, 0, 0)] * (count if field == "coordinates" else 3)
    charges = [0] * (count if field == "charges" else 3)
    with pytest.raises(ValueError, match=field):
        chematic.write_pdbqt_checked(mol, coords, charges)


@pytest.mark.parametrize("name", ["a\nb", "a\rb", "\0", "\t", "\x7f", "分子", "a" * 1025])
def test_name_cannot_inject_records_or_exceed_limit(name):
    with pytest.raises(ValueError, match="name"):
        chematic.write_pdbqt_checked(chematic.from_smiles("C"), [(0, 0, 0)], [0], name)


def test_width_edges_and_valid_names():
    mol = chematic.from_smiles("CC")
    for name in ["", "LIG", "with spaces", "x" * 1024]:
        text = chematic.write_pdbqt_checked(mol, [(9999.9994, -999.9994, -0.0001), (0, 0, 0)], [99.9994, -9.9994], name)
        lines = atom_lines(text)
        assert len(lines) == 2
        assert all(len(line) == 79 for line in lines)
        assert lines[0][30:54] == "9999.999-999.999  -0.000"
        assert lines[0][70:76] == "99.999"
        assert lines[1][70:76] == "-9.999"


@pytest.mark.parametrize("smiles", ["*", "[B]", "[Si]"])
def test_unknown_autodock_type_is_rejected(smiles):
    with pytest.raises(ValueError, match="AutoDock type"):
        chematic.write_pdbqt_checked(chematic.from_smiles(smiles), [(0, 0, 0)], [0])


def test_python_boundary_rejects_wrong_shapes_and_types():
    mol = chematic.from_smiles("C")
    with pytest.raises((TypeError, ValueError)):
        chematic.write_pdbqt_checked(mol, [(0, 0)], [0])
    with pytest.raises(TypeError):
        chematic.write_pdbqt_checked(mol, [(0, 0, 0)], ["bad"])
    with pytest.raises(TypeError):
        chematic.write_pdbqt_checked("C", [(0, 0, 0)], [0])
