"""Checked coordinate projection; independently decode PDB fixed columns."""
from pathlib import Path
import pytest
import chematic

ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / 'validation/openbabel_issues/mmcif/branched-residues.cif'
LONG = ROOT / 'validation/openbabel_issues/mmcif/long-component-name.cif'


def changed(column, value):
    text = FIXTURE.read_text(encoding="utf-8")
    rows = [s.split() for s in text.splitlines() if s.startswith('HETATM')]
    rows[0][column] = value
    return text.split('HETATM')[0] + '\n'.join(' '.join(r) for r in rows) + '\n#\n'


def convert(text, **kwargs):
    return chematic.mmcif_to_pdb_checked(text, allow_metadata_loss=True, **kwargs)


def test_auth_names_coordinates_and_source_are_preserved():
    text = FIXTURE.read_text(encoding="utf-8")
    result = convert(text)
    lines = result['pdb'].splitlines()
    assert lines[-1] == 'END'
    assert all(len(s) == 80 for s in lines[:-1])
    assert [s[17:20] for s in lines[:-1]] == ['NAG'] * 4 + ['BMA'] * 2
    assert [int(s[22:26]) for s in lines[:-1]] == [1, 1, 2, 2, 3, 3]
    assert [float(s[30:38]) for s in lines[:-1]] == [20, 22, 40, 42, 60, 62]
    assert result['source'] == chematic.parse_mmcif_residues(text)
    assert result['source_text'] == text
    assert result['namespace'] == 'auth'
    assert result['losses'] == ['identity_namespaces', 'entity_identifiers']


def test_loss_requires_explicit_argument():
    with pytest.raises(ValueError, match='allow_metadata_loss'):
        chematic.mmcif_to_pdb_checked(FIXTURE.read_text(encoding="utf-8"))


@pytest.mark.parametrize('allow', [False, True])
@pytest.mark.parametrize('namespace', ['auth', 'label'])
def test_five_character_names_never_truncate(allow, namespace):
    text = LONG.read_text(encoding="utf-8")
    # Make both namespaces independently representable except for the name.
    rows = [s.split() for s in text.splitlines() if s.startswith('HETATM')]
    for r in rows: r[8] = r[15]
    text = text.split('HETATM')[0] + '\n'.join(' '.join(r) for r in rows) + '\n#\n'
    with pytest.raises(ValueError, match='residue name'):
        chematic.mmcif_to_pdb_checked(text, namespace=namespace, allow_metadata_loss=allow)


@pytest.mark.parametrize('column,value,field', [
    (16, 'LONG', 'residue name'), (17, 'AB', 'chain identifier'),
    (15, '10000', 'residue number'), (15, '01', 'residue number'),
    (15, '+1', 'residue number'), (9, 'AA', 'insertion code'),
    (18, 'ABCDE', 'atom name'), (4, 'AB', 'alternate location'),
    (10, '9999.9996', 'x'), (13, '1.1', 'occupancy'),
    (14, '999.999', 'B factor'), (1, '100000', 'serial'),
])
def test_unrepresentable_fields_fail(column, value, field):
    with pytest.raises(ValueError, match=field):
        convert(changed(column, value))


def test_namespace_selection_has_no_fallback():
    with pytest.raises(ValueError, match='residue number'):
        convert(FIXTURE.read_text(encoding="utf-8"), namespace='label')
    with pytest.raises(ValueError, match='atom name'):
        convert(changed(18, '?'))
    with pytest.raises(ValueError, match='namespace'):
        convert(FIXTURE.read_text(encoding="utf-8"), namespace='auto')


def test_projection_collisions_are_errors():
    # Distinct label components sharing an author residue key cannot merge.
    with pytest.raises(ValueError, match='same PDB residue identifier'):
        convert(changed(5, 'MAN'))
    # Distinct label atom aliases sharing an author atom name cannot merge.
    with pytest.raises(ValueError, match='same PDB atom name'):
        convert(changed(18, 'O6'))
    with pytest.raises(ValueError, match='duplicate PDB serial'):
        convert(changed(1, '2'))


def test_numeric_precision_is_reported():
    result = convert(changed(10, '20.123456'))
    assert result['max_coordinate_rounding_error'] == pytest.approx(0.000456)
    assert 'numeric_rounding' in result['losses']
    assert result['source']['atoms'][0]['x'] == 20.123456
    assert float(result['pdb'].splitlines()[0][30:38]) == 20.123


def test_only_model_one_is_supported():
    with pytest.raises(ValueError, match='model 1 only'):
        convert(changed(19, '2'))


def test_label_namespace_uses_only_label_aliases():
    text = FIXTURE.read_text(encoding="utf-8")
    rows = [s.split() for s in text.splitlines() if s.startswith('HETATM')]
    for r in rows:
        r[8] = str(int(r[15]) + 10)
        r[3] = 'L' + r[3]
    text = text.split('HETATM')[0] + '\n'.join(' '.join(r) for r in rows) + '\n#\n'
    lines = convert(text, namespace='label')['pdb'].splitlines()[:-1]
    assert [int(s[22:26]) for s in lines] == [11, 11, 12, 12, 13, 13]
    assert all(s[21] == 'E' for s in lines)
    assert lines[0][12:16].strip() == 'LC1'
