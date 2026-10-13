import math

import pytest

import chematic


def test_eem_charged_molecules_conserve_formal_charge():
    hydroxide = chematic.from_smiles("[O-][H]")
    charges = hydroxide.eem_charges([[0.0, 0.0, 0.0], [0.97, 0.0, 0.0]])
    assert sum(charges) == pytest.approx(-1.0, abs=1e-12)
    assert charges[0] < -1.0 and charges[1] > 0.0

    ammonium = chematic.from_smiles("[N+]([H])([H])([H])[H]")
    a = 1.01 / math.sqrt(3.0)
    charges = ammonium.eem_charges(
        [[0.0, 0.0, 0.0], [a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]]
    )
    assert sum(charges) == pytest.approx(1.0, abs=1e-12)


def test_eem_rejects_invalid_coordinate_count():
    with pytest.raises(ValueError, match="expected 2 coordinates"):
        chematic.from_smiles("[H][H]").eem_charges([[0.0, 0.0, 0.0]])
