"""Source atom order and numeric preservation; fixtures are original project data."""
import hashlib
import json
import struct
from pathlib import Path

import pytest
import chematic

PACKET = json.loads((Path(__file__).parents[3] / 'validation/pdbqt_atom_data/cases.json').read_text(encoding="utf-8"))
CASES = PACKET['cases']


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['id'])
def test_numeric_arrays_and_formal_charge_separation(case):
    text, expected = case['input'], case['expected']
    assert hashlib.sha256(text.replace('\r\n', '\n').encode()).hexdigest() == case['sha256']
    for source in [text, text.encode('utf-8')]:
        mol, coords, charges = chematic.from_pdbqt_with_data(source)
        assert coords == expected['coords']
        assert charges == expected['partial_charges']
        actual = [v for point in coords for v in point] + charges
        wanted = [v for point in expected['coords'] for v in point] + expected['partial_charges']
        assert [struct.pack('>d', v) for v in actual] == [struct.pack('>d', v) for v in wanted]
        assert mol.atom_count == len(coords) == len(charges)
        graph = json.loads(mol.to_moljson())
        assert graph['bonds'] == []
        assert all(atom['charge'] == 0 for atom in graph['atoms'])
        assert mol.formal_charge == 0
        # The existing element-only API remains compatible.
        assert mol.smiles == chematic.from_pdbqt(text).smiles


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['id'])
def test_bridge_reuses_source_coordinates_and_charges(case):
    source, expected = case['input'], case['expected']
    rewritten = chematic.convert_format(source, 'pdbqt', 'pdbqt')
    _, coords, charges = chematic.from_pdbqt_with_data(rewritten)
    assert coords == expected['coords']
    assert charges == pytest.approx(expected['partial_charges'], abs=0.000500000001, rel=0)
    for fmt in ['xyz', 'pdb']:
        output = chematic.convert_format(source, 'pdbqt', fmt)
        _, coords = getattr(chematic, f'from_{fmt}')(output)
        assert coords == expected['coords']


def test_explicit_overrides_have_priority_without_defaulting_other_array():
    source = CASES[2]['input']
    new_coords = [[3, 2, 1], [6, 5, 4]]
    output = chematic.convert_format(source, 'pdbqt', 'pdbqt', coords=new_coords)
    assert chematic.from_pdbqt_with_data(output)[1:] == (new_coords, [2, -1])
    output = chematic.convert_format(source, 'pdbqt', 'pdbqt', charges=[0.2, -0.2])
    assert chematic.from_pdbqt_with_data(output)[1:] == (CASES[2]['expected']['coords'], [0.2, -0.2])
    for coords, charges in [([], []), ([[0, 0, 0]], [0]), ([[0, 0, 0]]*3, [0]*3)]:
        with pytest.raises(ValueError, match='coordinates'):
            chematic.convert_format(source, 'pdbqt', 'pdbqt', coords=coords, charges=charges)
    with pytest.raises(ValueError, match='does not fit'):
        chematic.convert_format(source, 'pdbqt', 'pdbqt', charges=[100, 0])


@pytest.mark.parametrize('field,replacement', [('charge', 'NaN'), ('charge', 'inf'), ('charge', 'bad'), ('type', 'Unknown'), ('x', 'NaN'), ('x', 'bad')])
def test_malformed_atom_rejects_whole_parse(field, replacement):
    rows = CASES[2]['input'].splitlines()
    i = next(i for i, row in enumerate(rows) if row.startswith('HETATM'))
    row = rows[i]
    if field == 'x':
        rows[i] = row[:30] + f'{replacement:>8}' + row[38:]
    else:
        tail = row[66:].split()
        tail[0 if field == 'charge' else 1] = replacement
        rows[i] = row[:66] + '    ' + ' '.join(tail)
    with pytest.raises(ValueError):
        chematic.from_pdbqt_with_data('\n'.join(rows))


def test_invalid_utf8_and_empty_compatibility():
    with pytest.raises(UnicodeDecodeError):
        chematic.from_pdbqt_with_data(b'\xff')
    mol, coords, charges = chematic.from_pdbqt_with_data('REMARK no atoms\n')
    assert mol.atom_count == 0
    assert coords == charges == []
    with pytest.raises(ValueError, match='at least one atom'):
        chematic.convert_format('REMARK no atoms\n', 'pdbqt', 'pdbqt')
