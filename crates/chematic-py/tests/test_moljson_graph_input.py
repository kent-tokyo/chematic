"""Known local MolJSON graph fields must retain meaning or fail explicitly."""
import importlib.util
import json
from pathlib import Path

import chematic
import pytest

ROOT = Path(__file__).resolve().parents[3]
PACKET = json.loads((ROOT/'validation/moljson_graph_input/cases.json').read_text(encoding='utf-8'))
SPEC = importlib.util.spec_from_file_location('moljson_observer', ROOT/'scripts/validation/check_moljson_graph_input.py')
OBSERVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OBSERVER)


@pytest.mark.parametrize('case', PACKET['valid'], ids=lambda c: c['id'])
def test_supported_fields_survive_direct_and_common_routes(case):
    text = case['input']
    for _ in range(3):
        mol = chematic.from_moljson(text)
        direct = mol.to_moljson()
        common = chematic.convert_format(text, 'moljson', 'moljson')
        assert OBSERVER.snapshot(direct) == case['expected']
        assert OBSERVER.snapshot(common) == case['expected']
        text = direct


@pytest.mark.parametrize('case', PACKET['invalid'], ids=lambda c: c['id'])
def test_malformed_or_ambiguous_graphs_are_refused(case):
    with pytest.raises(ValueError):
        chematic.from_moljson(case['input'])
    with pytest.raises(ValueError):
        chematic.convert_format(case['input'], 'moljson', 'moljson')
