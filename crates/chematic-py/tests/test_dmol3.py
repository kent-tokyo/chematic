"""Original synthetic MIT OR Apache-2.0 DMol3 .incoor contracts."""
import ast
import hashlib
import json
import math
from pathlib import Path

import pytest
import chematic

ROOT = Path(__file__).resolve().parents[3]


def emit(g):
    return chematic.write_dmol3(g['elements'], g['positions_bohr'], cell_bohr=g['cell_bohr'])


def test_manifest_projection_and_roundtrip():
    fixture = json.loads((ROOT/'validation/openbabel_issues/dmol3/contract-fixtures.json').read_text(encoding="utf-8"))['records'][0]
    text = fixture['payload']
    assert hashlib.sha256(text.encode()).hexdigest() == fixture['sha256']
    g = chematic.parse_dmol3(text)
    for key, value in fixture['expected'].items():
        if key not in ('units', 'tolerance'):
            assert g[key] == value
    assert g['source_text'] == text
    back = chematic.parse_dmol3(emit(g))
    for key in g.keys()-{'source_text'}:
        assert back[key] == g[key]
    g['positions_bohr'][0][0] = 0.12345678901234566
    g['cell_bohr'] = None
    back = chematic.parse_dmol3(emit(g))
    assert back['positions_bohr'] == g['positions_bohr'] and back['periodicity'] == 0
    assert back['cell_bohr'] is None and back['cell_angstrom'] is None


@pytest.mark.parametrize('text', [
    '', '$coordinates\n$end', '$coordinates\nC 0 0 0',
    'preamble\n$coordinates\nC 0 0 0\n$end', '$coordinates angstrom\nC 0 0 0\n$end',
    '$coordinates\nC 0 0 0 0.5\n$end', '$coordinates\nC1 0 0 0\n$end',
    '$coordinates\n6 0 0 0\n$end', '$coordinates\nc 0 0 0\n$end',
    '$coordinates\nC NaN 0 0\n$end', '$coordinates\nC 1e309 0 0\n$end',
    '$coordinates\nC 0 0 0\n$end\n$end', '$coordinates\nC 0 0 0\n$cell vectors\n$end',
    '$cell vectors\n1 0 0\n0 1 0\n$coordinates\nC 0 0 0\n$end',
    '$cell vectors\n1 0 0\n0 1 0\n1 1 1e-12\n$coordinates\nC 0 0 0\n$end',
])
def test_invalid_reader_input(text):
    with pytest.raises(ValueError):
        chematic.parse_dmol3(text)


@pytest.mark.parametrize('kwargs', [dict(max_atoms=0), dict(max_line_len=5), dict(max_input_bytes=2)])
def test_reader_limits(kwargs):
    with pytest.raises(ValueError):
        chematic.parse_dmol3('$coordinates\nC 0 0 0\n$end', **kwargs)


@pytest.mark.parametrize('elements,positions,kwargs', [
    ([], [], {}), (['C'], [], {}), (['c'], [[0, 0, 0]], {}),
    (['C'], [[float('inf'), 0, 0]], {}), (['C'], [[0, 0, 0]], {'cell_bohr': [[1, 0, 0]]}),
    (['C'], [[0, 0, 0]], {'cell_bohr': [[1, 0, 0], [2, 0, 0], [0, 0, 1]]}),
])
def test_invalid_writer_input(elements, positions, kwargs):
    with pytest.raises(ValueError):
        chematic.write_dmol3(elements, positions, **kwargs)


def test_precision_d_exponents_crlf_and_stub():
    text = '\r\n$coordinates\r\n\r\nH -0.0 4.9406564584124654D-324 1.7976931348623157e308\r\n$end\r\n'
    g = chematic.parse_dmol3(text)
    back = chematic.parse_dmol3(emit(g))
    assert back['positions_bohr'] == g['positions_bohr']
    assert math.copysign(1, back['positions_bohr'][0][0]) == -1
    assert g['source_text'] == text
    assert g['bohr_angstrom'] == 0.529177210544
    p = Path(__file__).resolve().parents[1]/'python/chematic/__init__.pyi'
    names = {n.name for n in ast.parse(p.read_text(encoding="utf-8")).body if isinstance(n, ast.FunctionDef)}
    assert {'parse_dmol3', 'write_dmol3'} <= names
