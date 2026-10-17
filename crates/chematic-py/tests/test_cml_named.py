"""Explicit CML IDs and graph export (Open Babel #2120)."""
import itertools
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest
import chematic as c

NS = {'c': 'http://www.xml-cml.org/schema'}
FIXTURE = Path(__file__).resolve().parents[3] / 'validation/openbabel_issues/cml-named/water.mol2'


def test_mol2_names_and_xyz_export():
    text = FIXTURE.read_text()
    record = json.loads(c.parse_mol2_record_json(text))
    mol = c.from_mol2(text)
    names = [a['name'] for a in record['atoms']]
    output = c.write_cml_named_checked(mol, names, record['coordinates'])
    root = ET.fromstring(output)
    atoms = root.findall('c:atomArray/c:atom', NS)
    assert [a.get('id') for a in atoms] == names
    assert [[float(a.get(k)) for k in ['x3', 'y3', 'z3']] for a in atoms] == record['coordinates']
    assert [b.get('atomRefs2') for b in root.findall('c:bondArray/c:bond', NS)] == ['O_site H_one', 'O_site H-two']
    assert mol.formula == 'H2O'


@pytest.mark.parametrize('id', ['', '1C', '_C', 'C C', 'C:1', 'C.1', 'C"/>', 'C&x', '炭素', 'C'*129])
def test_invalid_id(id):
    with pytest.raises(ValueError, match='invalid CML atom ID'):
        c.write_cml_named_checked(c.from_smiles('C'), [id])


def test_duplicates_and_count():
    mol = c.from_smiles('CO')
    with pytest.raises(ValueError, match='duplicate'):
        c.write_cml_named_checked(mol, ['X', 'X'])
    for ids in [[], ['C1'], ['C1', 'O1', 'Extra']]:
        with pytest.raises(ValueError, match='one ID per atom'):
            c.write_cml_named_checked(mol, ids)


@pytest.mark.parametrize('ids', ['C', b'C', [1], [None], [b'C']])
def test_id_types(ids):
    with pytest.raises(TypeError):
        c.write_cml_named_checked(c.from_smiles('C'), ids)


@pytest.mark.parametrize('number', [float('nan'), float('inf'), -float('inf')])
def test_nonfinite_xyz(number):
    with pytest.raises(ValueError, match='non-finite coordinate'):
        c.write_cml_named_checked(c.from_smiles('C'), ['C1'], [[0, number, 0]])


@pytest.mark.parametrize('rows', [[], [[0,0,0],[1,1,1]]])
def test_coordinate_count(rows):
    with pytest.raises(ValueError, match='one 3D row per atom'):
        c.write_cml_named_checked(c.from_smiles('C'), ['C1'], rows)


@pytest.mark.parametrize('row', [[0,0], [0,0,0,0], [None,0,0], ['x',0,0]])
def test_coordinate_shape_and_type(row):
    with pytest.raises((ValueError, TypeError)):
        c.write_cml_named_checked(c.from_smiles('C'), ['C1'], [row])


@pytest.mark.parametrize('smiles', ['C[C@H](F)Cl', 'C/C=C/C', '[CH3:7]O', '*', 'N->[Cu+2]', 'C~C'])
def test_unsupported(smiles):
    mol = c.from_smiles(smiles)
    with pytest.raises(ValueError, match='unsupported'):
        c.write_cml_named_checked(mol, [f'A{i}' for i in range(mol.atom_count)])


def test_bounded_iterators_and_exceptions():
    consumed = []
    def ids():
        for i in itertools.count():
            consumed.append(i)
            yield f'A{i}'
    with pytest.raises(ValueError, match='one ID'):
        c.write_cml_named_checked(c.from_smiles('C'), ids())
    assert consumed == [0, 1]
    def broken():
        yield 'C1'
        raise RuntimeError('source failure')
    with pytest.raises(RuntimeError, match='source failure'):
        c.write_cml_named_checked(c.from_smiles('C'), broken())
    used = []
    def coords():
        for i in itertools.count():
            used.append(i)
            yield [i, 0, 0]
    with pytest.raises(ValueError, match='one 3D row'):
        c.write_cml_named_checked(c.from_smiles('C'), ['C1'], coords())
    assert used == [0, 1]


def test_no_coordinates_and_case_sensitive_ids():
    output = c.write_cml_named_checked(c.from_smiles('CO'), iter(['Name', 'name']))
    root = ET.fromstring(output)
    assert [a.get('id') for a in root.findall('c:atomArray/c:atom', NS)] == ['Name', 'name']
    assert 'x3=' not in output and 'x2=' not in output


def test_graph_limits():
    for mol in [c.RWMol().GetMol(), c.from_smiles('C' * 4097)]:
        with pytest.raises(ValueError, match='4096 atoms'):
            c.write_cml_named_checked(mol, [])
