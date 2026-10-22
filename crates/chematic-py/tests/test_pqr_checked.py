"""Checked PQR atom fields, separate from original document/source retention."""
import json
import struct
from pathlib import Path
import pytest
import chematic

PACKET=json.loads((Path(__file__).parents[3]/'validation/pqr-checked-2026-10-22/cases.json').read_text(encoding="utf-8"))
CASES=[c for c in PACKET['cases'] if c.get('surfaces',['python','wasm'])!=['wasm']]
ERRORS=(ValueError,TypeError,OverflowError)


@pytest.mark.parametrize('case',CASES,ids=lambda c:c['id'])
def test_record_packet(case):
    atoms=json.loads(case['records_json'])
    if case['expected']=='refuse':
        with pytest.raises(ERRORS):chematic.write_pqr_checked(atoms)
        return
    text=chematic.write_pqr_checked(iter(atoms))
    actual=chematic.parse_pqr(text)['atoms']
    assert len(actual)==len(case['fields'])
    for a,expected in zip(actual,case['fields']):
        for key in ['group_pdb','atom_name','res_name','chain_id','icode','element']:
            assert a[key]==expected[key]
        for key in ['serial','res_seq']:assert str(a[key])==expected[key]
        for key,bits in expected['bits'].items():assert struct.pack('>d',a[key]).hex()==bits


@pytest.mark.parametrize('value',[float('nan'),float('inf'),-float('inf')])
@pytest.mark.parametrize('field',['x','y','z','charge','radius'])
def test_nonfinite_fields(value,field):
    a=json.loads(CASES[0]['records_json'])[0];a[field]=value
    with pytest.raises(ValueError):chematic.write_pqr_checked([a])


@pytest.mark.parametrize('source',['records',b'records',bytearray(b'records'),[1],[None],[[]]])
def test_records_are_a_nontext_iterable_of_dictionaries(source):
    with pytest.raises(ERRORS):chematic.write_pqr_checked(source)


@pytest.mark.parametrize('limit',[True,False,-1,0,1.5,'1024',64*1024*1024+1])
def test_invalid_output_options_do_not_consume_geometry(limit):
    def forbidden():
        raise AssertionError('invalid option consumed atoms')
        yield
    with pytest.raises(ERRORS):chematic.write_pqr_checked(forbidden(),max_output_bytes=limit)


def test_output_byte_boundary():
    atoms=json.loads(CASES[0]['records_json'])
    output=chematic.write_pqr_checked(atoms)
    size=len(output.encode())
    assert chematic.write_pqr_checked(atoms,max_output_bytes=size)==output
    with pytest.raises(ValueError):chematic.write_pqr_checked(atoms,max_output_bytes=size-1)


def test_bounded_generator_and_atom_count_boundary():
    atom=json.loads(CASES[0]['records_json'])[0]
    assert len(chematic.parse_pqr(chematic.write_pqr_checked(atom for _ in range(100_000)))['atoms'])==100_000
    seen=[]
    def unlimited():
        while True:
            seen.append(1)
            if len(seen)>100_001:raise AssertionError('unbounded consumption')
            yield atom
    with pytest.raises(ValueError,match='100000'):chematic.write_pqr_checked(unlimited())
    assert len(seen)==100_001


def test_legacy_defaults_and_bytes_stay_available():
    atom=dict(atom_name='CA',x=1.2345678901234567,y=0.,z=0.)
    assert chematic.write_pqr([atom])=='REMARK   PQR written by chematic\nATOM 1 CA UNL 1 1.235 0.000 0.000 0.0000 0.0000\n'
    with pytest.raises(ValueError):chematic.write_pqr_checked([atom])
