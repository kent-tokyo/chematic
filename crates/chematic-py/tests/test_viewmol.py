"""ViewMol geometry: source values, fixed axes and explicit unit conversion."""
import copy,json
from pathlib import Path
import pytest
import chematic
P=json.loads((Path(__file__).resolve().parents[3]/'validation/viewmol_preservation/cases.json').read_text(encoding='utf-8'))
def bits(rows):return [[v.hex() for v in row] for row in rows]
@pytest.mark.parametrize('c',P['valid'],ids=lambda c:c['id'])
def test_roundtrip_and_projection(c):
 r=chematic.parse_viewmol(c['text'].encode());assert r==c['expected']
 for _ in range(3):
  assert bits(chematic.viewmol_positions_angstrom(r))==bits(c['expected_angstrom'])
  r=chematic.parse_viewmol(chematic.write_viewmol(r));assert r==c['expected']
@pytest.mark.parametrize('c',P['refused'],ids=lambda c:c['id'])
def test_refusals(c):
 with pytest.raises(ValueError,match='^ViewMol line '):chematic.parse_viewmol(c['text'])
@pytest.mark.parametrize('c',P['projection_refused'],ids=lambda c:c['id'])
def test_raw_roundtrip_survives_unsafe_projection(c):
 r=chematic.parse_viewmol(c['text']);assert r==c['expected'];assert chematic.parse_viewmol(chematic.write_viewmol(r))==r
 with pytest.raises(ValueError,match='projection'):chematic.viewmol_positions_angstrom(r)
def test_edits_limits_and_strict_fields():
 r=copy.deepcopy(P['valid'][0]['expected']);r['title']='edited';r['factor_to_angstrom']='.5';r['atoms'][0]['position'][0]='-.12345678901234567';r['atoms'][0]['fixed_axes']='zy'
 assert chematic.parse_viewmol(chematic.write_viewmol(r))==r
 for key in ['max_atoms','max_input_bytes','max_line_len']:
  with pytest.raises(ValueError):chematic.write_viewmol(r,**{key:0})
  with pytest.raises(ValueError):chematic.parse_viewmol(P['valid'][0]['text'],**{key:0})
 for key,value in [('charge',1),('isotope',13),('stereo',None)]:
  bad=copy.deepcopy(r);bad['atoms'][0][key]=value
  with pytest.raises(ValueError,match='unknown field'):chematic.write_viewmol(bad)
 for value in [0.125,None,True,'NaN','1e-999']:
  bad=copy.deepcopy(r);bad['atoms'][0]['position'][0]=value
  with pytest.raises(ValueError):chematic.write_viewmol(bad)
 with pytest.raises(TypeError):chematic.parse_viewmol(object())
 with pytest.raises(UnicodeDecodeError):chematic.parse_viewmol(b'\xff')
