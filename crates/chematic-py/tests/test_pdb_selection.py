"""Original synthetic MIT OR Apache-2.0 identity-preserving PDB contracts."""
import ast
import hashlib
import json
from pathlib import Path

import pytest
import chematic

FOLDER = Path(__file__).resolve().parents[3]/'validation/openbabel_issues/pdb-selection'
TEXT = (FOLDER/'rich.pdb').read_text(encoding="utf-8")


def atom_lines(text):
    return [s for s in text.splitlines() if s.startswith(('ATOM  ', 'HETATM'))]


def test_manifest_and_ligand_removal_preserve_original_lines():
    manifest = json.loads((FOLDER/'fixtures.json').read_text(encoding="utf-8"))
    for record in manifest['records']:
        text = (FOLDER/record['file']).read_text(encoding="utf-8")
        assert hashlib.sha256(text.encode()).hexdigest() == record['sha256']
        r = chematic.select_pdb_atoms(text, [0, 1, 2, 3])
        assert atom_lines(r['pdb']) == atom_lines(text)[:4]
        assert r['old_to_new'] == [0, 1, 2, 3, None, None]
        assert [a['residue_name'] for a in r['atoms']] == ['VAL', 'VAL', 'LEU', 'LEU']
        assert r['source_text'] == text and r['source_atom_count'] == 6
        assert r['retargeted_ter_lines'] == []


def test_source_order_anisou_ter_and_map():
    r = chematic.select_pdb_atoms(TEXT, [4, 2, 0])
    assert r['selected_source_indices'] == [0, 2, 4]
    assert r['old_to_new'] == [0, None, 1, None, 2, None]
    assert r['retargeted_ter_lines'] == [5, 10, 15]
    assert [s[6:11].strip() for s in r['pdb'].splitlines() if s.startswith('TER')] == ['2', '6', '10']
    assert len([s for s in r['pdb'].splitlines() if s.startswith('ANISOU')]) == 3
    assert 'CONECT' not in r['pdb']
    assert r['atoms'][0]['alt_loc'] == 'A' and r['atoms'][0]['atom_name'] == ' CA '
    assert r['atoms'][1]['residue_number'] == -2 and r['atoms'][1]['insertion_code'] == 'I'
    assert r['atoms'][1]['position_angstrom'] == [10, 0, 0]
    assert chematic.select_pdb_atoms(r['pdb'])['pdb'] == r['pdb']


def test_partial_connectivity_retains_duplicate_targets():
    text = TEXT.replace('END\n', 'CONECT    1    2    5    5    9\nEND\n')
    r = chematic.select_pdb_atoms(text, [0, 1, 2])
    assert 'CONECT    1    2    5    5\n' in r['pdb']


def test_metadata_loss_is_explicit_and_raw_source_stays_available():
    text = 'HEADER    original\nSEQRES   original\n'+TEXT
    with pytest.raises(ValueError, match='consent'):
        chematic.select_pdb_atoms(text, [0])
    r = chematic.select_pdb_atoms(text, [0], allow_metadata_loss=True)
    assert r['losses'] == [{'line': 1, 'record_name': 'HEADER'}, {'line': 2, 'record_name': 'SEQRES'}]
    assert r['source_text'] == text and 'HEADER' not in r['pdb']


@pytest.mark.parametrize('model', [1, 7])
@pytest.mark.parametrize('crlf', [False, True])
def test_single_model_and_empty_selection(model, crlf):
    rows = [s for s in TEXT.splitlines() if not s.startswith('CONECT') and s != 'END']
    text = f'MODEL     {model:4d}\n'+'\n'.join(rows)+'\nENDMDL\nEND\n'
    if crlf:
        text = text.replace('\n', '\r\n')
    r = chematic.select_pdb_atoms(text, [0])
    assert r['model_num'] == model and r['pdb'].endswith('ENDMDL\nEND\n')
    assert r['source_text'] == text
    assert chematic.select_pdb_atoms(text, [])['pdb'] == 'END\n'


@pytest.mark.parametrize('indices', [[0, 0], [6], [-1]])
def test_invalid_indices(indices):
    with pytest.raises((ValueError, OverflowError)):
        chematic.select_pdb_atoms(TEXT, indices)


@pytest.mark.parametrize('kwargs', [dict(max_atoms=5), dict(max_input_bytes=2), dict(max_line_len=20), dict(max_records=2)])
def test_limits(kwargs):
    with pytest.raises(ValueError):
        chematic.select_pdb_atoms(TEXT, **kwargs)


@pytest.mark.parametrize('text', [
    'ENDMDL\n'+TEXT, 'MODEL        1\n'+TEXT, TEXT+'ATOM',
    TEXT.replace('CONECT    1    2    2', 'CONECT    1   99'),
    TEXT.replace('CONECT    1    2    2', 'CONECT    1    1'),
    TEXT.replace('TER       3', 'TER       4'),
    TEXT.replace('    100    200', '    bad    200'),
    TEXT.replace(' 23.45', '   inf'),
    'REMARK λ\n'+TEXT,
])
def test_malformed_records_fail(text):
    with pytest.raises(ValueError):
        chematic.select_pdb_atoms(text)


def test_empty_source_none_vs_empty_and_independent_results():
    assert chematic.select_pdb_atoms('')['pdb'] == 'END\n'
    assert chematic.select_pdb_atoms(TEXT)['pdb'] == TEXT
    assert chematic.select_pdb_atoms(TEXT, [])['pdb'] == 'END\n'
    r = chematic.select_pdb_atoms(TEXT, [0])
    r['atoms'][0]['residue_name'] = 'changed'
    r['old_to_new'][0] = 999
    fresh = chematic.select_pdb_atoms(TEXT, [0])
    assert fresh['atoms'][0]['residue_name'] == 'VAL' and fresh['old_to_new'][0] == 0


def test_stub_and_hard_limit():
    p = Path(__file__).resolve().parents[1]/'python/chematic/__init__.pyi'
    names = {n.name for n in ast.parse(p.read_text(encoding="utf-8")).body if isinstance(n, ast.FunctionDef)}
    assert 'select_pdb_atoms' in names
    with pytest.raises(ValueError):
        chematic.select_pdb_atoms('REMARK'+'a'*1025+'\nEND\n', max_line_len=2**32, allow_metadata_loss=True)
