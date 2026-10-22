"""The ordinary same-format bridge must use the immutable source record path."""
import hashlib
import json
from pathlib import Path

import pytest
import chematic

PACKET = json.loads(
    (Path(__file__).parents[3] / 'validation/cjson_source_conversion/cases.json')
    .read_text(encoding='utf-8')
)
CASES = PACKET['cases']


def variant(case, name):
    if name == 'crlf':
        return case['input'].replace('\n', '\r\n')
    if name == 'minified':
        return json.dumps(json.loads(case['input']), ensure_ascii=True, separators=(',', ':'))
    return case['input']


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['id'])
@pytest.mark.parametrize('layout', PACKET['variants'])
def test_exact_source_including_uninterpreted_fields(case, layout):
    assert hashlib.sha256(case['input'].encode()).hexdigest() == case['sha256']
    text = variant(case, layout)
    assert chematic.convert_format(text, 'cjson', 'cjson') == text


@pytest.mark.parametrize('source,target', [('CJSON', '.cjson'), ('.cjson', 'CJSON')])
def test_normalized_aliases_and_explicit_defaults(source, target):
    text = CASES[0]['input']
    assert chematic.convert_format(text, source, target, coords=None, charges=None,
                                   name='LIG', comment='') == text


@pytest.mark.parametrize('override', [dict(coords=[]), dict(coords=[[0, 0, 0]] * 2),
                                     dict(charges=[]), dict(charges=[0, 0]),
                                     dict(name='edited'), dict(comment='edited')])
def test_overrides_are_refused_instead_of_ignored(override):
    with pytest.raises(ValueError, match='CJSON source-preserving.*overrides'):
        chematic.convert_format(CASES[0]['input'], 'cjson', 'cjson', **override)


@pytest.mark.parametrize('mutate', [
    lambda v: v.update(chemicalJson=0),
    lambda v: v['atoms'].update(formalCharges=[0]),
    lambda v: v['atoms']['coords'].update({'3d': [0]}),
    lambda v: v.update(partialCharges={'method': [0]}),
    lambda v: v['bonds'].update(order=[2.2]),
    lambda v: v['bonds']['connections'].update(index=[0, 0]),
    lambda v: v['bonds']['connections'].update(index=[0, 2]),
])
def test_local_profile_errors_propagate(mutate):
    value = json.loads(CASES[0]['input'])
    mutate(value)
    with pytest.raises(ValueError):
        chematic.convert_format(json.dumps(value), 'cjson', 'cjson')


@pytest.mark.parametrize('suffix', [', "chemicalJson":1', ', "opaque":{"id":1,"id":2}',
                                   ', "opaque":{"id":1,"\\u0069d":2}'])
def test_duplicate_keys_are_rejected_even_in_opaque_fields(suffix):
    text = CASES[0]['input'].rstrip()[:-1] + suffix + '}'
    with pytest.raises(ValueError, match='duplicate object key'):
        chematic.convert_format(text, 'cjson', 'cjson')


def test_input_limit_still_precedes_retention():
    with pytest.raises(ValueError, match='maximum input size'):
        chematic.convert_format(' ' * (16 * 1024 * 1024 + 1), 'cjson', 'cjson')


def test_other_routes_still_convert_graphs():
    mol = chematic.from_smiles(chematic.convert_format(CASES[0]['input'], 'cjson', 'smiles'))
    assert mol.atom_count == 2
