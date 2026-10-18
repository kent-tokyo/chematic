"""Original synthetic MIT OR Apache-2.0 MOPAC Cartesian contracts."""
import ast
import hashlib
import json
import math
from pathlib import Path

import pytest
import chematic

ROOT = Path(__file__).resolve().parents[3]
FIELDS = ('keywords', 'title', 'comment', 'labels', 'optimization_flags', 'translation_vectors_angstrom', 'translation_flags')


def emit(g):
    return chematic.write_mopac_cartesian(g['elements'], g['positions_angstrom'], **{k: g[k] for k in FIELDS})


def test_manifest_and_metadata_roundtrip():
    fixture = json.loads((ROOT/'validation/openbabel_issues/mopac/contract-fixtures.json').read_text(encoding="utf-8"))['records'][0]
    text = fixture['payload']
    assert hashlib.sha256(text.encode()).hexdigest() == fixture['sha256']
    g = chematic.parse_mopac_cartesian(text)
    for key, value in fixture['expected'].items():
        if key not in ('units', 'tolerance'):
            assert g[key] == value
    assert g['source_text'] == text
    back = chematic.parse_mopac_cartesian(emit(g))
    for key in g.keys()-{'source_text'}:
        assert back[key] == g[key]
    g['positions_angstrom'][0] = [0.12345678901234566, 2, 3]
    g['optimization_flags'][0] = None
    g['labels'][0] = None
    g['keywords'] = 'PM6 XYZ'
    back = chematic.parse_mopac_cartesian(emit(g))
    for key in g.keys()-{'source_text'}:
        assert back[key] == g[key]


@pytest.mark.parametrize('dimension', range(4))
def test_partial_periodicity_and_flags_absence(dimension):
    vectors = [[4, 1, 0], [-1, 5, 2], [1, 0, -6]][:dimension]
    g = chematic.parse_mopac_cartesian(chematic.write_mopac_cartesian(['C', 'H'], [[0, 0, 0], [1, 2, 3]], optimization_flags=[None, ['+0', '+1', '+2']], translation_vectors_angstrom=vectors))
    assert g['optimization_flags'] == [None, ['0', '1', '2']]
    assert g['translation_flags'] == [None]*dimension
    assert g['translation_vectors_angstrom'] == vectors
    assert g['periodicity'] == dimension


@pytest.mark.parametrize('body', [
    '', 'C 0 1 0 1 0 1 0 0 0', 'C 0 3 0 1 0 1', 'C 0 t 0 1 0 1',
    'C14.0 0 0 0', 'X 0 0 0', 'XX 0 0 0', '0 0 0 0', '99 0 0 0',
    'At 0 0 0', 'Fr 0 0 0', 'Es 0 0 0', '+ 0 0 0', 'C NaN 0 0',
    'C(label)0 0 0', 'C(unclosed 0 0 0', 'C(nested(a)) 0 0 0', 'C(λ) 0 0 0',
    'C('+('a'*39)+') 0 0 0', 'C 0 0 0\n\nH 1 2 3',
    'Tv 1 0 0\nC 0 0 0', 'C 0 0 0\nTv(label) 1 0 0', 'C 0 0 0\nTv 0 0 0',
    'C 0 0 0\nTv 1 0 0\nTv 2 0 0', 'C 0 0 0\nTv 1 0 0\nTv 0 1 0\nTv 1 1 1e-12',
])
def test_invalid_reader_input(body):
    with pytest.raises(ValueError):
        chematic.parse_mopac_cartesian('PM7 XYZ\noriginal\n\n'+body)


@pytest.mark.parametrize('keyword', ['A0', 'PM7 a0=1', '++', 'PM7 &'])
def test_unsupported_units_and_continuations(keyword):
    with pytest.raises(ValueError, match='line 1'):
        chematic.parse_mopac_cartesian(keyword+'\n\n\nC 0 0 0')


@pytest.mark.parametrize('kwargs', [{'max_atoms': 0}, {'max_input_bytes': 2}, {'max_line_len': 4}])
def test_reader_limits(kwargs):
    with pytest.raises(ValueError):
        chematic.parse_mopac_cartesian('\n\n\nC 0 0 0', **kwargs)


@pytest.mark.parametrize('elements,positions,kwargs', [
    ([], [], {}), (['C'], [], {}), (['c'], [[0, 0, 0]], {}),
    (['Es'], [[0, 0, 0]], {}), (['C'], [[float('inf'), 0, 0]], {}),
    (['C'], [[0, 0, 0]], {'labels': []}), (['C'], [[0, 0, 0]], {'labels': ['bad)']}),
    (['C'], [[0, 0, 0]], {'optimization_flags': []}),
    (['C'], [[0, 0, 0]], {'optimization_flags': [['1', '3', '0']]}),
    (['C'], [[0, 0, 0]], {'translation_flags': [['0', '0', '0']]}),
    (['C'], [[0, 0, 0]], {'title': 'bad\nheader'}),
])
def test_invalid_writer_input(elements, positions, kwargs):
    with pytest.raises(ValueError):
        chematic.write_mopac_cartesian(elements, positions, **kwargs)


def test_crlf_fortran_precision_and_stub():
    text = 'PM7 XYZ\r\ntitle \r\ncomment\r\nH -0.0 4.9406564584124654D-324 1.7976931348623157e308\r\n'
    g = chematic.parse_mopac_cartesian(text)
    back = chematic.parse_mopac_cartesian(emit(g))
    assert back['positions_angstrom'] == g['positions_angstrom']
    assert math.copysign(1, back['positions_angstrom'][0][0]) == -1
    assert g['title'] == 'title ' and g['source_text'] == text
    p = Path(__file__).resolve().parents[1]/'python/chematic/__init__.pyi'
    names = {n.name for n in ast.parse(p.read_text(encoding="utf-8")).body if isinstance(n, ast.FunctionDef)}
    assert {'parse_mopac_cartesian', 'write_mopac_cartesian'} <= names


@pytest.mark.parametrize("keywords", ["PM7 XYZ +\nCHARGE=+1", "PM7 + XYZ\nCHARGE=+1 +\nPRECISE"])
def test_keyword_continuation_keeps_descriptions_and_roundtrips(keywords):
    text = keywords + "\ntitle\ncomment\nC 0 0 0\nO 1.4 0.2 -0.5\n"
    parsed = chematic.parse_mopac_cartesian(text)
    assert parsed["keywords"] == keywords
    assert parsed["title"] == "title"
    assert parsed["comment"] == "comment"
    assert parsed["source_text"] == text
    assert parsed["elements"] == ["C", "O"]
    output = chematic.write_mopac_cartesian(
        parsed["elements"], parsed["positions_angstrom"],
        keywords=parsed["keywords"], title=parsed["title"], comment=parsed["comment"],
    )
    reread = chematic.parse_mopac_cartesian(output)
    assert reread["keywords"] == keywords
    assert reread["positions_angstrom"] == parsed["positions_angstrom"]
