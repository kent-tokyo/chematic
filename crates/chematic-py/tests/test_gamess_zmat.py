"""GAMESS Cartesian frame and explicit IZMAT definitions, Open Babel #2523."""

import pytest
import chematic

ELEMENTS = ["C", "O", "H", "N"]
POSITIONS = [[2.0, -1.0, 3.0], [3.2, -1.0, 3.0], [2.0, 0.0, 3.0], [2.2, 0.3, 4.0]]
REFS = [[], [0], [0, 1], [0, 1, 2]]


def write(**kwargs):
    return chematic.write_gamess_zmat(
        ELEMENTS, POSITIONS, REFS, charge=-1, multiplicity=2, **kwargs
    )


def test_geometry_and_internal_definitions():
    text = write()
    assert "COORD=UNIQUE UNITS=ANGS" in text
    assert "ICHARG=-1 MULT=2 NZVAR=6" in text
    assert "3,4,1,2,3\n $END" in text
    rows = text.split("C1\n")[1].splitlines()[:4]
    assert [list(map(float, row.split()[2:])) for row in rows] == POSITIONS
    assert all(len(line) <= 80 for line in text.splitlines())


@pytest.mark.parametrize(
    "kwargs",
    [
        {"significant_digits": 0},
        {"significant_digits": 18},
        {"significant_digits": True},
        {"title": "$END"},
        {"title": "a\nb"},
        {"title": "!comment"},
        {"title": "x" * 81},
    ],
)
def test_invalid_options(kwargs):
    with pytest.raises(ValueError):
        write(**kwargs)


@pytest.mark.parametrize(
    "references",
    [
        [[], [True], [0, 1], [0, 1, 2]],
        [[], [0], [0, 0], [0, 1, 2]],
        [[], [0], [0, 1], [0, 1, 4]],
        [[], [0], [0, 1], [0, 1, 2, 3]],
    ],
)
def test_invalid_references(references):
    with pytest.raises(ValueError):
        chematic.write_gamess_zmat(
            ELEMENTS, POSITIONS, references, charge=0, multiplicity=1
        )


@pytest.mark.parametrize(
    "charge,spin", [(True, 1), (0, True), (0, 0), (0, -1), (2**40, 1)]
)
def test_strict_electronic_state_types(charge, spin):
    with pytest.raises(ValueError):
        chematic.write_gamess_zmat(
            ELEMENTS, POSITIONS, REFS, charge=charge, multiplicity=spin
        )


def test_empty_missing_and_nonfinite_geometry():
    for elements, positions, refs in [
        ([], [], []),
        (ELEMENTS, POSITIONS[:3], REFS),
        (["C"], [[float("inf"), 0.0, 0.0]], [[]]),
        (["C"] * 1025, [[0.0, 0.0, 0.0]] * 1025, [[]] * 1025),
    ]:
        with pytest.raises(ValueError):
            chematic.write_gamess_zmat(
                elements, positions, refs, charge=0, multiplicity=1
            )


def test_monatomic_template_has_no_internal_coordinates():
    text = chematic.write_gamess_zmat(
        ["He"], [[0.0, 0.0, 0.0]], [[]], charge=0, multiplicity=1
    )
    assert "NZVAR=0" in text
    assert "$ZMAT" not in text


def test_xyz_conversion_retains_all_atoms_and_refuses_multiple_frames():
    xyz = "3\nwater\nO 0 0 0\nH .96 0 0\nH -.24 .93 0\n"
    text = chematic.xyz_to_gamess_zmat(xyz, [[], [0], [0, 1]], charge=0, multiplicity=1)
    assert sum(line.startswith("H 1.0 ") for line in text.splitlines()) == 2
    for bad in [
        xyz + xyz,
        xyz + "garbage\n",
        "",
        " " * (1 << 20) + xyz,
        xyz.replace("water", "Lattice=cell"),
        xyz.replace("O 0 0 0", "O 0 0 0 1.2"),
    ]:
        with pytest.raises(ValueError):
            chematic.xyz_to_gamess_zmat(
                bad, [[], [0], [0, 1]], charge=0, multiplicity=1
            )
