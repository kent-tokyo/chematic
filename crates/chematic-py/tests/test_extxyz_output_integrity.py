"""Preflight geometry/properties and keep metadata literal across readers."""
import json
import math
from pathlib import Path
import struct
import pytest
import chematic
P=json.loads((Path(__file__).parents[3]/'validation/extxyz_output_integrity/cases.json').read_text())
def bits(values):return [struct.pack('>d',float(v)) for row in values for v in row]
def write(c):return chematic.to_extxyz(chematic.from_smiles('CO'),c['coords'],lattice=c['lattice'],properties=c['properties'],info=c['info'])
@pytest.mark.parametrize('c',P['valid'],ids=lambda c:c['id'])
def test_valid_fields_are_preserved(c):
    out=write(c);r=chematic.from_extxyz(out)
    assert bits(r['coords'])==bits(c['coords'])
    assert r['info']==c['info']
    if c['lattice'] is None:assert r['lattice'] is None
    else:assert bits([r['lattice']])==bits([c['lattice']])
    assert r['properties'].keys()==c['properties'].keys()
    for name,values in c['properties'].items():assert bits(r['properties'][name])==bits(values)
@pytest.mark.parametrize('c',P['invalid'],ids=lambda c:c['id'])
def test_unsafe_names_are_errors(c):
    with pytest.raises(ValueError):write(c)
@pytest.mark.parametrize('value',[math.nan,math.inf,-math.inf])
@pytest.mark.parametrize('field',['coords','lattice','properties'])
def test_nonfinite_values_are_errors(value,field):
    kwargs=dict(coords=[[0.,0.,0.],[1.,0.,0.]])
    if field=='coords':kwargs['coords'][1][2]=value
    elif field=='lattice':kwargs['lattice']=[value]*9
    else:kwargs['properties']={'data':[[1.],[value]]}
    with pytest.raises(ValueError):chematic.to_extxyz(chematic.from_smiles('CO'),**kwargs)
