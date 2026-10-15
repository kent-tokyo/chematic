"""Independent PQR/PDBQT field contract; no upstream source fixture copied."""
from pathlib import Path
import pytest
import chematic

FIXTURE = Path(__file__).resolve().parents[3] / 'validation/openbabel_issues/pqr/charged-records.pqr'
TYPES = ['N', 'C', 'HD', 'Zn', 'Cl']


def test_charges_identity_and_sidecar():
    text = FIXTURE.read_text(encoding="utf-8")
    parsed = chematic.parse_pqr(text)
    report = chematic.pqr_to_pdbqt(text, TYPES, allow_radius_loss=True)
    rows = [s for s in report['pdbqt'].splitlines() if s.startswith(('ATOM', 'HETATM'))]
    assert len(rows) == 5
    assert [float(r[70:76]) for r in rows] == pytest.approx([-0.412, 0.123, 0.289, 2, -1])
    assert [r[12:16].strip() for r in rows] == ['N', 'CA', 'H', 'ZN', 'CL']
    assert [int(r[6:11]) for r in rows] == [7, 19, 20, 31, 45]
    assert [r[77:79].strip() for r in rows] == TYPES
    assert rows[0][17:27] == 'ALA A  52A'
    assert rows[4][21] == ' '
    assert report['source_charges'] == [a['charge'] for a in parsed['atoms']]
    assert report['source_radii'] == [a['radius'] for a in parsed['atoms']]
    assert report['source_coords'][0] == pytest.approx((-0.966123, 1.523, 1.412))
    assert report['radii_omitted'] is True
    assert report['charge_provenance'] == 'pqr_atom_records'
    assert report['atom_type_provenance'] == 'caller_supplied'
    assert report['max_charge_rounding_error'] == pytest.approx(0.00045)
    assert report['max_coordinate_rounding_error'] == pytest.approx(0.000123)
    assert chematic.from_pdbqt(report['pdbqt']).rdkit_num_atoms() == 5


def test_default_requires_radius_loss_consent():
    with pytest.raises(ValueError, match='allow_radius_loss'):
        chematic.pqr_to_pdbqt(FIXTURE.read_text(encoding="utf-8"), TYPES)


@pytest.mark.parametrize('types', [[], TYPES[:-1], TYPES + ['C'], ['C'] + TYPES[1:], ['X'] + TYPES[1:]])
def test_invalid_type_arrays(types):
    with pytest.raises(ValueError):
        chematic.pqr_to_pdbqt(FIXTURE.read_text(encoding="utf-8"), types, allow_radius_loss=True)


@pytest.mark.parametrize('text', [
    '', 'REMARK empty',
    'ATOM 1 N ALA A 1 0 0 0 NaN 1.5',
    'ATOM 1 N ALA A 1 0 0 0 -0.4 -1.5',
    'ATOM 1 N A1AJ1 A 1 0 0 0 -0.4 1.5',
    'ATOM 1 N ALA AA 1 0 0 0 -0.4 1.5',
    'ATOM 100000 N ALA A 1 0 0 0 -0.4 1.5',
    'ATOM 1 N ALA A 1 9999.9996 0 0 -0.4 1.5',
    'ATOM 1 N ALA A 1 0 0 0 99.9996 1.5',
])
def test_invalid_input_is_rejected(text):
    with pytest.raises(ValueError):
        chematic.pqr_to_pdbqt(text, ['N'], allow_radius_loss=True)


def test_preserves_nonpolar_hydrogens_and_their_source_charge():
    text = 'ATOM 1 C LIG 1 0 0 0 -0.123 1.7\nATOM 2 H LIG 1 1 0 0 0.123 1.2\n'
    report = chematic.pqr_to_pdbqt(text, ['C', 'H'], allow_radius_loss=True)
    rows = [s for s in report['pdbqt'].splitlines() if s.startswith('ATOM')]
    assert len(rows) == 2
    assert [r[77:79].strip() for r in rows] == ['C', 'H']
    assert [float(r[70:76]) for r in rows] == [-0.123, 0.123]
