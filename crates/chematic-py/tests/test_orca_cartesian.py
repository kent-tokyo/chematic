"""Explicit numeric units and refusal of implicit Cartesian job selection."""
import json
import struct
from pathlib import Path
import pytest
import chematic

CASES=json.loads((Path(__file__).parents[3]/'validation/orca-cartesian-2026-10-22/cases.json').read_text())['cases']


@pytest.mark.parametrize('case',CASES,ids=lambda c:c['id'])
def test_checked_units_and_geometry(case):
    if case['expected']=='refuse':
        with pytest.raises(ValueError):chematic.parse_orca_cartesian(case['source'])
        if case['core_refuse']:
            with pytest.raises(ValueError):chematic.parse_orca_input(case['source'])
        return
    actual=chematic.parse_orca_cartesian(case['source']);expected=case['fields']
    assert actual['schema']=='chematic.orca-cartesian.v1'
    for key in ['coordinate_units','units_explicit','elements','charge','multiplicity','frozen','atom_extras']:
        assert actual[key]==expected[key]
    assert [[struct.pack('>d',v).hex() for v in row] for row in actual['coords']]==expected['coordinate_bits']


def test_bohr_values_are_tagged_without_rescaling():
    source='! Bohrs\n* xyz 0 1\nHe 2.0 0 0\n*\n'
    report=chematic.parse_orca_cartesian(source)
    assert report['coordinate_units']=='bohr' and report['coords'][0][0]==2.0
    assert chematic.parse_orca_input(source)['keywords']==['Bohrs']


@pytest.mark.parametrize('options',[dict(max_input_bytes=8),dict(max_atoms=0)])
def test_limits(options):
    with pytest.raises(ValueError):chematic.parse_orca_cartesian('* xyz 0 1\nHe 1 2 3\n*\n',**options)


def test_atom_budget_boundary():
    source='! Angs\n* xyz 0 1\n'+'He 0 0 0\n'*100_000+'*\n'
    assert len(chematic.parse_orca_cartesian(source)['coords'])==100_000
    with pytest.raises(ValueError):chematic.parse_orca_cartesian(source[:-2]+'He 0 0 0\n*\n')
