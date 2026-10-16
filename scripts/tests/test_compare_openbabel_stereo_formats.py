"""Controls for stereo-safe explicit hydrogen normalization in the paired gate."""

from __future__ import annotations

import sys
from pathlib import Path

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_openbabel_stereo_formats import evaluate  # noqa: E402


def score(source: str, output: str) -> dict:
    expected = Chem.MolFromSmiles(source)
    assert expected is not None
    return evaluate(expected, Chem.MolToSmiles(expected), 0, output)


def test_nonstereogenic_directional_h_is_a_representation_difference():
    for output in ("C(N)(N)=N[H]", "N=C(N)N"):
        result = score("[H]/N=C(/N)N", output)
        assert not result["strict_equal"]
        assert result["hydrogen_representation_equivalent"]
        assert result["semantic_equal"]


def test_stereo_and_isotope_loss_do_not_pass_hydrogen_normalization():
    for source, output in (
        ("[H]/N=C(/F)Cl", "N=C(F)Cl"),
        ("F/C=C/F", "F/C=C\\F"),
        ("[2H]N=C(N)N", "N=C(N)N"),
        ("[H:7]/N=C(/N)N", "N=C(N)N"),
    ):
        result = score(source, output)
        assert not result["semantic_equal"]
        assert not result["hydrogen_representation_equivalent"]
