"""Checked projection refuses payload before the graph bridge discards it."""
import json
from pathlib import Path
import pytest
import chematic

P=json.loads((Path(__file__).parents[3]/'validation/cml_cjson_projection/cases.json').read_text())
def variants(s):return [s,s.replace('><','>\n<'),s.replace('><','>\r\n<')]
@pytest.mark.parametrize('c',P['valid'],ids=lambda c:c['id'])
def test_direct_and_common_outputs_match(c):
    for s in variants(c['input']):
        out=chematic.project_cml_to_cjson(s)
        assert out==chematic.project_cml_to_cjson(s.encode())
        assert out==chematic.convert_format(s,'.CML','cjson')
        assert json.loads(out)['atoms']['elements']['number']
@pytest.mark.parametrize('c',P['invalid'],ids=lambda c:c['id'])
def test_unrepresented_fields_do_not_return_output(c):
    for s in variants(c['input']):
        with pytest.raises(ValueError):chematic.project_cml_to_cjson(s)
        with pytest.raises(ValueError):chematic.convert_format(s,'cml','cjson')
def test_input_and_override_boundaries():
    s=P['valid'][0]['input']
    for opts in [dict(coords=[]),dict(charges=[]),dict(name='override'),dict(comment='x')]:
        with pytest.raises(ValueError,match='overrides'):chematic.convert_format(s,'cml','cjson',**opts)
    with pytest.raises(UnicodeDecodeError):chematic.project_cml_to_cjson(b'\xff')
    with pytest.raises(ValueError,match='input bytes'):chematic.project_cml_to_cjson(' '*(16*1024*1024+1))
