"""Structured settings and exact Cartesian numbers survive checked output."""
import copy
import json
import math
import struct
from pathlib import Path

import chematic
import pytest

CASES=json.loads((Path(__file__).parents[3]/'validation/orca-checked-2026-10-24/cases.json').read_text())['cases']


@pytest.mark.parametrize('case',[c for c in CASES if c['lane']!='transport'],ids=lambda c:c['id'])
def test_structured_packet(case):
    doc=json.loads(case['input_json'])
    if case['expected']=='refuse':
        with pytest.raises((ValueError,TypeError)):chematic.write_orca_input_checked(doc)
        return
    text=chematic.write_orca_input_checked(doc)
    got=chematic.parse_orca_input(text)
    assert got==case['fields']
    assert [[struct.pack('>d',a[k]).hex() for k in ['x','y','z']] for a in got['coords']['atoms']]==case['coordinate_bits']
    assert chematic.write_orca_input_checked(got)==text


@pytest.mark.parametrize('value',[math.nan,math.inf,-math.inf])
def test_nonfinite_values(value):
    doc=json.loads(CASES[0]['input_json']);doc['coords']['atoms'][0]['x']=value
    with pytest.raises(ValueError):chematic.write_orca_input_checked(doc)


@pytest.mark.parametrize('limit',[True,False,0,-1,67_108_865,1.5])
def test_invalid_output_limit(limit):
    with pytest.raises((ValueError,TypeError,OverflowError)):
        chematic.write_orca_input_checked(json.loads(CASES[0]['input_json']),max_output_bytes=limit)


def test_exact_output_limit_and_coordinate_tag_alias():
    doc=json.loads(CASES[0]['input_json']);text=chematic.write_orca_input_checked(doc)
    assert chematic.write_orca_input_checked(doc,max_output_bytes=len(text.encode()))==text
    with pytest.raises(ValueError):chematic.write_orca_input_checked(doc,max_output_bytes=len(text.encode())-1)
    alias=copy.deepcopy(doc);alias['coords']['type']=alias['coords'].pop('kind')
    assert chematic.write_orca_input_checked(alias)==text
    alias['coords']['kind']='xyz'
    with pytest.raises(ValueError):chematic.write_orca_input_checked(alias)
