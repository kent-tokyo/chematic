import json
import math
from pathlib import Path

import chematic
import pytest


def test_selenium_stereo_reports_satisfied_violated_and_unevaluable():
    path = (Path(__file__).resolve().parents[2] / "chematic-3d" / "src" /
            "stereo_constraints_se_neutral_oracle.json")
    row = json.loads(path.read_text(encoding="utf-8"))["rows"][0]
    mol = chematic.from_smiles(row["smiles"]).add_hydrogens()
    coords = row["coords"]
    for geometry, expected in ((coords, "satisfied"),
                               ([[-x, y, z] for x, y, z in coords], "violated"),
                               ([[0.0, 0.0, 0.0]] * len(coords), "unevaluable")):
        report = mol.stereo_3d_diagnostics(geometry)
        assert report["reported"] == 1
        assert report[expected] == 1
        assert report["atom_stereo"][0]["status"] == expected
        reason = report["atom_stereo"][0]["reason"]
        assert (reason is not None) == (expected == "unevaluable")
        assert report["double_bond_stereo"] == []


def test_double_bond_and_square_planar_are_reported_separately():
    mol = chematic.from_smiles("F/C=C/F")
    trans = [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, -1, 0]]
    cis = [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]]
    assert mol.stereo_3d_diagnostics(trans)["double_bond_stereo"][0]["status"] == "satisfied"
    assert mol.stereo_3d_diagnostics(cis)["double_bond_stereo"][0]["status"] == "violated"
    pt = chematic.from_smiles("[Pt@SP1](F)(Cl)(Br)I")
    coords = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    row = pt.stereo_3d_diagnostics(coords)["atom_stereo"][0]
    assert row == {"atom": 0, "status": "unevaluable", "reason": "UnsupportedCoordination"}


def test_stereo_diagnostics_empty_and_invalid_coordinates():
    mol = chematic.from_smiles("CCO")
    report = mol.stereo_3d_diagnostics([[0, 0, 0]] * 3)
    assert (report["reported"], report["violated"], report["unevaluable"]) == (0, 0, 0)
    with pytest.raises(ValueError, match="finite and match"):
        mol.stereo_3d_diagnostics([[0, 0, 0]])
    with pytest.raises(ValueError, match="finite and match"):
        mol.stereo_3d_diagnostics([[math.nan, 0, 0]] * 3)
