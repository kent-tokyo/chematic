"""Regression examples from RDKit's fixed-coordinate selenium stereo oracle."""

import json
from pathlib import Path

import chematic


def test_neutral_selenium_stereo_and_mirrors_on_fixed_coordinates():
    fixture = (Path(__file__).resolve().parents[2] / "chematic-3d" / "src" /
               "stereo_constraints_se_neutral_oracle.json")
    rows = json.loads(fixture.read_text(encoding="utf-8"))["rows"]
    assert len(rows) == 24
    for row in rows:
        mol = chematic.from_smiles(row["smiles"]).add_hydrogens()
        assert mol._timed_stereo_valid(row["coords"]), (row["smiles"], row["seed"])
        mirror = [[-x, y, z] for x, y, z in row["coords"]]
        assert not mol._timed_stereo_valid(mirror), (row["smiles"], row["seed"])
