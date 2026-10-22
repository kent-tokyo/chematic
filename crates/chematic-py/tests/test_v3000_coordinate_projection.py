"""Source V3000 xyz values survive the ordinary CJSON bridge, not just copying."""
import json
from pathlib import Path
import struct
import pytest
import chematic

PACKET = json.loads((Path(__file__).parents[3] / 'validation/v3000_coordinate_projection/cases.json').read_text(encoding="utf-8"))

def bits(rows):
    return [struct.pack('>d', float(v)) for row in rows for v in row]

@pytest.mark.parametrize('case', PACKET['cases'], ids=lambda c: c['id'])
@pytest.mark.parametrize('variant', PACKET['variants'])
def test_source_xyz_and_graph_in_cjson(case, variant):
    source = case['input']
    if variant == 'crlf': source = source.replace('\n', '\r\n')
    if variant == 'no-final-newline': source = source.rstrip('\n')
    expected = case['expected']
    mol, name, xyz = chematic.from_mol_v3000_with_coords_3d(source.encode())
    assert name == case['id']
    assert bits(xyz) == bits(expected['coordinates'])
    assert len(xyz) == mol.atom_count
    _, _, xy = chematic.from_mol_v3000_with_coords(source)
    assert bits(xy) == bits([row[:2] for row in xyz])
    value = json.loads(chematic.convert_format(source, '.V3000', 'cjson'))
    atoms = value['atoms']; bonds = value['bonds']
    flat = atoms['coords']['3d']
    assert bits([flat[i:i+3] for i in range(0, len(flat), 3)]) == bits(expected['coordinates'])
    assert atoms['elements']['number'] == expected['numbers']
    assert atoms.get('formalCharges', [0]*len(xyz)) == expected['formal_charges']
    assert atoms.get('isotopes', [0]*len(xyz)) == expected['isotopes']
    pairs = bonds['connections']['index']
    assert [[pairs[2*i], pairs[2*i+1], order] for i, order in enumerate(bonds['order'])] == expected['bonds']


def test_explicit_coordinates_keep_precedence_and_legacy_2d_shape():
    case = PACKET['cases'][0]
    coords = [[9.125, 8.25, 7.5], [6.0, 5.0, 4.0]]
    value = json.loads(chematic.convert_format(case['input'], 'v3000', 'cjson', coords=coords))
    assert value['atoms']['coords']['3d'] == sum(coords, [])
    empty = json.loads(chematic.convert_format(case['input'], 'v3000', 'cjson', coords=[]))
    assert 'coords' not in empty['atoms']
    _, _, xy = chematic.from_mol_v3000_with_coords(case['input'])
    assert all(len(row) == 2 for row in xy)

@pytest.mark.parametrize('token', ['NaN', 'inf', '-inf'])
def test_nonfinite_source_is_refused(token):
    with pytest.raises(ValueError):
        chematic.from_mol_v3000_with_coords_3d(PACKET['cases'][0]['input'].replace('0.123456789012345', token))

def test_size_and_utf8_boundaries():
    with pytest.raises(ValueError, match='input bytes'):
        chematic.from_mol_v3000_with_coords_3d(' '*(16*1024*1024+1))
    with pytest.raises(UnicodeDecodeError):
        chematic.from_mol_v3000_with_coords_3d(b'\xff')


def test_coordinate_outputs_receive_source_xyz_without_an_override():
    source = next(c['input'] for c in PACKET['cases'] if c['id'] == 'offset_charges')
    _, _, xyz = chematic.from_mol_v3000_with_coords_3d(source)
    assert xyz == [[0.0, 0.0, 2.125], [3.0, 0.0, 2.125]]
    for target in ['pdb', 'xyz', 'pdbqt']:
        output = chematic.convert_format(source, 'v3000', target)
        assert '2.125' in output
