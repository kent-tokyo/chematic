"""Immutable V3000 source retention includes 3D values below perception tolerance."""
import hashlib
import json
from pathlib import Path
import pytest
import chematic

PACKET=json.loads((Path(__file__).parents[3]/'validation/v3000_source_roundtrip/cases.json').read_text(encoding="utf-8"))


def variant(text,layout):
    return text.replace('\n','\r\n') if layout=='crlf' else text.rstrip('\n') if layout=='no-final-newline' else text


@pytest.mark.parametrize('case',PACKET['cases'],ids=lambda c:c['id'])
@pytest.mark.parametrize('layout',PACKET['variants'])
def test_exact_source_in_both_public_routes(case,layout):
    assert hashlib.sha256(case['input'].encode()).hexdigest()==case['sha256']
    text=variant(case['input'],layout)
    assert chematic.roundtrip_mol_v3000_source(text)==text
    assert chematic.roundtrip_mol_v3000_source(text.encode())==text
    assert chematic.convert_format(text,'.v3000','mol-v3000')==text


@pytest.mark.parametrize('kwargs',[dict(coords=[]),dict(coords=[[0,0]]*2),dict(charges=[]),dict(charges=[0,0]),dict(name='new'),dict(comment='new')])
def test_same_format_overrides_are_refused(kwargs):
    with pytest.raises(ValueError,match='V3000 source-preserving.*overrides'):
        chematic.convert_format(PACKET['cases'][0]['input'],'v3000','v3000',**kwargs)


@pytest.mark.parametrize('change',[
    lambda s:s.replace('M  END',''),lambda s:s+'$$$$\n',lambda s:s+s,
    lambda s:s.replace('END CTAB','END CTAB-'),
    lambda s:s.replace('M  V30 2 O','M  V30 1 O'),lambda s:s.replace('M  V30 1 C','M  V30 0 C'),
    lambda s:s.replace('M  V30 1 1 1 2','M  V30 0 1 1 2'),
    lambda s:s.replace('M  V30 1 1 1 2','M  V30 1 1 1 99'),
    lambda s:s.replace('0.123456789012345','NaN'),
    lambda s:s.replace('M  V30 END CTAB','M  V30 END CTAB\nM  V30 BEGIN CTAB\nM  V30 END CTAB'),
    lambda s:s.replace('M  V30 END CTAB','M  V30 END CTAB\nM  V30 ignored'),
])
def test_invalid_source_never_returns_partial_text(change):
    text=change(PACKET['cases'][0]['input'])
    for call in [chematic.roundtrip_mol_v3000_source,lambda s:chematic.convert_format(s,'v3000','v3000')]:
        with pytest.raises(ValueError):call(text)


def test_utf8_and_size_boundaries():
    with pytest.raises(UnicodeDecodeError): chematic.roundtrip_mol_v3000_source(b'\xff')
    with pytest.raises(ValueError,match='input bytes'): chematic.roundtrip_mol_v3000_source(' '*(16*1024*1024+1))
