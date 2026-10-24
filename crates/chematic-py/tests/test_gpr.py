"""Checked GPR records on the actual native extension."""
import copy
import json
from pathlib import Path
import pytest
import chematic
PACKET=json.loads((Path(__file__).resolve().parents[3]/'validation/gpr_preservation/cases.json').read_text(encoding="utf-8"))
@pytest.mark.parametrize('case',PACKET['valid'],ids=lambda c:c['id'])
def test_roundtrip(case):
    record=chematic.parse_gpr(case['text'].encode())
    assert record==case['expected']
    for _ in range(3):
        record=chematic.parse_gpr(chematic.write_gpr(record))
        assert record==case['expected']
@pytest.mark.parametrize('case',PACKET['refused'],ids=lambda c:c['id'])
def test_refusal(case):
    with pytest.raises(ValueError,match='^GPR line '):chematic.parse_gpr(case['text'])
def test_edits_absence_limits_and_unrepresented_meaning():
    r=copy.deepcopy(PACKET['valid'][2]['expected'])
    r['coordinates'][0]['position_nm'][0]='.12345678901234567'
    r['partial_charges'][0]['value']='-.12345678901234567'
    r['bonds'][0]['order']='D'
    assert chematic.parse_gpr(chematic.write_gpr(r))==r
    for key in ['max_atoms','max_bonds','max_input_bytes','max_line_len']:
        with pytest.raises(ValueError):chematic.write_gpr(r,**{key:0})
        with pytest.raises(ValueError):chematic.parse_gpr(PACKET['valid'][2]['text'],**{key:0})
    for key,value in [('formal_charge',1),('isotope',13),('stereo',None)]:
        bad=copy.deepcopy(r);bad['atoms'][0][key]=value
        with pytest.raises(ValueError,match='unknown field'):chematic.write_gpr(bad)
    for value in [0.125,True,None,'NaN','1e-999']:
        bad=copy.deepcopy(r);bad['coordinates'][0]['position_nm'][0]=value
        with pytest.raises(ValueError):chematic.write_gpr(bad)
    r['partial_charges']=None
    assert '!Charges' not in chematic.write_gpr(r)
    with pytest.raises(TypeError):chematic.parse_gpr(object())
    with pytest.raises(UnicodeDecodeError):chematic.parse_gpr(b'\xff')
