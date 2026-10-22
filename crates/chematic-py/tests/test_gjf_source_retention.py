"""Exact source retention is distinct from the Cartesian geometry projection."""
import json
from pathlib import Path
import pytest
import chematic

_ROOT=Path(__file__).parents[3]
_CASES=json.loads((_ROOT/'validation/gjf-cartesian-2026-10-22/source-cases.json').read_text(encoding="utf-8"))['cases']
_SOURCE=(_ROOT/'validation/gjf-cartesian-2026-10-22/configured-source.gjf').read_text(encoding="utf-8")


@pytest.mark.parametrize('case',_CASES,ids=lambda case:case['id'])
def test_original_source_packet(case):
    source=case['source']
    functions=[lambda:chematic.roundtrip_gjf_record(source),
               lambda:chematic.convert_format(source,'.GJF','com'),
               lambda:chematic.convert_format(source,'.COM','gjf')]
    for function in functions:
        if case['expected']=='retain': assert function().encode('utf8')==source.encode('utf8')
        else:
            with pytest.raises(ValueError):function()


@pytest.mark.parametrize('override',[dict(coords=[]),dict(coords=[[0,0,0]]),dict(charges=[]),dict(name='new title'),dict(comment='changed'),dict(comment=None),dict(comment=False)])
def test_source_conversion_cannot_ignore_requested_edits(override):
    with pytest.raises(ValueError,match='overrides'):
        chematic.convert_format(_SOURCE,'gjf','com',**override)


def test_other_gjf_output_requires_explicit_geometry():
    with pytest.raises(ValueError,match='source-preserving'):
        chematic.convert_format('C','smiles','gjf')


def test_limits_include_original_bom_and_are_checked_before_source_copy():
    source='\ufeff'+_SOURCE
    size=len(source.encode('utf8'))
    assert chematic.roundtrip_gjf_record(source,max_input_bytes=size,max_atoms=2)==source
    for kwargs in [dict(max_input_bytes=size-1),dict(max_atoms=1)]:
        with pytest.raises(ValueError,match='exceeding limit'):
            chematic.roundtrip_gjf_record(source,**kwargs)


def test_atom_boundary_and_both_public_paths():
    prefix='# HF/STO-3G\n\nresource control\n\n0 1\n'
    source=prefix+'H 0 0 0\n'*100_000+'\n'
    assert chematic.roundtrip_gjf_record(source)==source
    assert chematic.convert_format(source,'com','gjf')==source
    too_large=prefix+'H 0 0 0\n'*100_001+'\n'
    for function in [lambda:chematic.roundtrip_gjf_record(too_large),lambda:chematic.convert_format(too_large,'gjf','gjf')]:
        with pytest.raises(ValueError):function()
