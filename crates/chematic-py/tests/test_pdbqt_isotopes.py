"""Explicit isotope loss for rigid PDBQT, Open Babel #3005."""
import json
import math
from pathlib import Path

import pytest
import chematic

CASES = json.loads((Path(__file__).resolve().parents[3] / "validation/openbabel_issues/pdbqt-isotopes/cases.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["name"])
def test_labels_require_projection_with_complete_serial_sidecar(case):
    mol = chematic.from_smiles(case["smiles"])
    before = mol.smiles
    coords = [[float(i), -2.345, 6.789] for i in range(mol.atom_count)]
    charges = [-0.123] * mol.atom_count
    if case["isotopes_omitted"]:
        first = case["isotopes_omitted"][0]
        with pytest.raises(ValueError, match=f"atom {first['source_atom_index']}: isotope {first['mass_number']}"):
            chematic.write_pdbqt_checked(mol, coords, charges)
    report = chematic.project_pdbqt_isotopes(mol, coords, charges)
    assert report["isotopes_omitted"] == case["isotopes_omitted"]
    assert report["pdbqt"] == mol.to_pdbqt([tuple(p) for p in coords], charges, "LIG")
    if not case["isotopes_omitted"]:
        assert report["pdbqt"] == chematic.write_pdbqt_checked(mol, coords, charges)
    lines = [line for line in report["pdbqt"].splitlines() if line.startswith("ATOM")]
    assert len(lines) == mol.atom_count
    for i, line in enumerate(lines):
        assert len(line) == 79
        assert int(line[6:11]) == i + 1
        assert [float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]] == coords[i]
        assert float(line[70:76]) == charges[i]
    assert mol.smiles == before


def test_donor_type_does_not_depend_on_hydrogen_isotope():
    args = ([[i, 0, 0] for i in range(4)], [-0.2, 0.3, -0.4, 0.1])
    labelled = chematic.from_smiles("[2H]O[13CH2][3H]")
    plain = chematic.from_smiles("[H]O[CH2][H]")
    report = chematic.project_pdbqt_isotopes(labelled, *args)
    assert report["pdbqt"] == chematic.write_pdbqt_checked(plain, *args)
    assert [line[77:79].strip() for line in report["pdbqt"].splitlines() if line.startswith("ATOM")] == ["HD", "OA", "C", "H"]


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf, 10000, -1000])
def test_projection_does_not_bypass_fixed_column_checks(bad):
    mol = chematic.from_smiles("[2H]")
    with pytest.raises(ValueError):
        chematic.project_pdbqt_isotopes(mol, [[bad, 0, 0]], [0])


def test_acetophenone_d8_does_not_create_donor_hydrogens():
    labelled = next(c["smiles"] for c in CASES if c["name"] == "acetophenone-d8")
    mol = chematic.from_smiles(labelled)
    ordinary = chematic.from_smiles(labelled.replace("[2H]", "[H]"))
    coords, charges = [[1, 2, 3]] * 17, [i * 0.01 for i in range(17)]
    report = chematic.project_pdbqt_isotopes(mol, coords, charges)
    assert len(report["isotopes_omitted"]) == 8
    assert report["pdbqt"] == chematic.write_pdbqt_checked(ordinary, coords, charges)
    types = [line[77:79].strip() for line in report["pdbqt"].splitlines() if line.startswith("ATOM")]
    assert len(types) == 17
    assert types.count("H") == 8
    assert "HD" not in types
    assert "TORSDOF" not in report["pdbqt"]


def test_projection_does_not_bypass_other_validation():
    mol = chematic.from_smiles("[2H]O")
    for coords, charges, name in [([], [0, 0], ""), ([[0, 0, 0]] * 2, [], ""),
                                  ([[0, 0, 0]] * 2, [0, 100], ""),
                                  ([[0, 0, 0]] * 2, [0, 0], "\nATOM")]:
        with pytest.raises(ValueError):
            chematic.project_pdbqt_isotopes(mol, coords, charges, name)
    with pytest.raises(ValueError, match="AutoDock type"):
        chematic.project_pdbqt_isotopes(chematic.from_smiles("[10B]"), [[0, 0, 0]], [0])
