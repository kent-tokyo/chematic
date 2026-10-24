"""Python binding for the fixed-coordinate RDKit 2026.09.1 USR profile."""

import pytest

import chematic


def test_rdkit_usr_matches_pinned_four_point_oracle():
    coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)]
    expected = [
        1.5305013515761465, 0.5315690104056001, 0.733435705162006,
        1.5, 1.118033988749895, 0.0,
        2.4419572339080924, 1.427180741099017, -1.0227286856189084,
        1.9604048132409448, 1.2871724702702962, -0.7105779995726041,
    ]
    assert chematic.rdkit_usr_descriptors(coords) == pytest.approx(expected, abs=1e-12)


def test_rdkit_usr_rejects_invalid_coordinates():
    with pytest.raises(ValueError, match="at least three atoms"):
        chematic.rdkit_usr_descriptors([(0.0, 0.0, 0.0)])
    with pytest.raises(ValueError, match="not finite"):
        chematic.rdkit_usr_descriptors([(0.0, 0.0, 0.0), (1.0, float("nan"), 0.0), (2.0, 0.0, 0.0)])
