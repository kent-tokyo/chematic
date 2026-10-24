"""GROMOS96 labels/units/blocks are not inferred from successful coordinate parsing."""
import copy
import json
from pathlib import Path
import pytest
import chematic

P=json.loads((Path(__file__).resolve().parents[3]/'validation/gromos96_preservation/cases.json').read_text(encoding='utf-8'))

@pytest.mark.parametrize('c',P['valid'],ids=lambda c:c['id'])
def test_cycles(c):
    unit=c['record']['length_unit'];r=chematic.parse_gromos96(c['text'].encode('utf-8'),length_unit=unit)
    for _ in range(3):
        assert r==c['record'];assert isinstance(chematic.gromos96_geometry_angstrom(r),list)
        r=chematic.parse_gromos96(chematic.write_gromos96(r),length_unit=unit)
    assert r==c['record']

@pytest.mark.parametrize('c',P['refused'],ids=lambda c:c['id'])
def test_refused(c):
    with pytest.raises(ValueError,match='GROMOS96 line '):chematic.parse_gromos96(c['text'])

@pytest.mark.parametrize('c',P['projection_refused'],ids=lambda c:c['id'])
def test_unsafe_projection_retains_raw_blocks(c):
    r=chematic.parse_gromos96(c['text'],length_unit=c['record']['length_unit'])
    assert chematic.parse_gromos96(chematic.write_gromos96(r),length_unit=r['length_unit'])==c['record']
    with pytest.raises(ValueError,match='GROMOS96 line '):chematic.gromos96_geometry_angstrom(r)

def test_edits_context_types_and_whole_document_limits():
    r=copy.deepcopy(P['valid'][0]['record']);line=r['blocks'][1]['lines'][0];r['blocks'][1]['lines'][0]=line[:24]+'.5'.rjust(15)+line[39:]
    g=chematic.gromos96_geometry_angstrom(chematic.parse_gromos96(chematic.write_gromos96(r),length_unit='nm'))
    assert g[0]['positions_angstrom'][0]==[5.,2.5,5.]
    assert g[0]['box_vectors_angstrom']==[[20.,0.,0.],[0.,40.,0.],[0.,0.,80.]]
    assert g[0]['identities'][1]['atom_name']=='Cl'
    unknown=chematic.parse_gromos96(P['valid'][0]['text']);assert unknown['length_unit']=='unspecified'
    with pytest.raises(ValueError):chematic.gromos96_geometry_angstrom(unknown)
    for kwargs in ({'max_atoms':1},{'max_frames':0},{'max_line_len':10},{'max_input_bytes':10}):
        with pytest.raises(ValueError):chematic.parse_gromos96(P['valid'][0]['text'],**kwargs)
        with pytest.raises(ValueError):chematic.write_gromos96(r,**kwargs)
    for bad in ({**r,'unknown':1},{**r,'length_unit':'bohr'},{**r,'blocks':[{'kind':'TITLE','leading':[],'lines':[1]}]}):
        with pytest.raises(ValueError):chematic.write_gromos96(bad)
    with pytest.raises(ValueError):chematic.parse_gromos96(P['valid'][0]['text'],length_unit='bohr')
    with pytest.raises(UnicodeDecodeError):chematic.parse_gromos96(b'\xff')
    with pytest.raises(TypeError):chematic.parse_gromos96(1)
