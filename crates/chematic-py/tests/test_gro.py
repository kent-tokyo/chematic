"""GRO records preserve names and explicit units; elements require caller data."""
import ast
from pathlib import Path
import pytest
import chematic

ROOT=Path(__file__).resolve().parents[3]
FIXTURE=ROOT/'validation/openbabel_issues/gro/two-letter-elements.gro'


def atom(name='CA',precision=3,velocity=False):
    width=precision+5
    out=f'{99999:5d}{"ALA":<5}{name:>5}{0:5d}'
    out+=''.join(f'{v:{width}.{precision}f}' for v in [-1.25,0.5,4.75])
    if velocity:out+=''.join(f'{v:{width}.{precision+1}f}' for v in [-.125,.25,.5])
    return out


def frame(row=None,box='1 2 3'):
    if row is None:row=atom()
    return f'synthetic t= 2.5\n1\n{row}\n{box}\n'


def test_two_letter_names_do_not_infer_truncated_elements():
    text=FIXTURE.read_text(encoding="utf-8")
    d=chematic.parse_gro(text)
    assert [a['atom_name'] for a in d['atoms']]==['BE','CE','FE']
    assert [a['element'] for a in d['atoms']]==[None]*3
    assert [a['position_nm'][0] for a in d['atoms']]==[0,2,4]
    assert d['element_policy']=='unresolved' and d['source_text']==text


def test_explicit_elements_are_canonical_and_follow_source_order():
    d=chematic.parse_gro(FIXTURE.read_text(encoding="utf-8"),elements=['Be','Ce','Fe'])
    assert [a['element'] for a in d['atoms']]==['Be','Ce','Fe']
    assert d['element_policy']=='caller_supplied'


@pytest.mark.parametrize('element',['C','Ca'])
def test_ca_name_can_be_resolved_explicitly_without_guessing(element):
    d=chematic.parse_gro(frame(),elements=[element])
    assert d['atoms'][0]['atom_name']=='CA' and d['atoms'][0]['element']==element


@pytest.mark.parametrize('elements',[[],['Be'],['Be','Ce','Fe','C']])
def test_element_counts_must_match(elements):
    with pytest.raises(ValueError,match='elements were supplied'):
        chematic.parse_gro(FIXTURE.read_text(encoding="utf-8"),elements=elements)


@pytest.mark.parametrize('symbol',['BE','be','Xx','', 'Fe2+', ' C'])
def test_bad_element_symbols_do_not_guess_or_normalize(symbol):
    with pytest.raises(ValueError,match='canonical element'):
        chematic.parse_gro(frame(),elements=[symbol])


@pytest.mark.parametrize('precision',[2,3,5,8,12])
@pytest.mark.parametrize('velocity',[False,True])
def test_decimal_precision_and_optional_velocities(precision,velocity):
    d=chematic.parse_gro(frame(atom(precision=precision,velocity=velocity)))
    a=d['atoms'][0]
    assert a['position_nm']==[-1.25,.5,4.75]
    assert a['velocity_nm_per_ps']==([-.125,.25,.5] if velocity else None)
    assert a['residue_number']==99999 and a['atom_number']==0
    assert a['residue_name']=='ALA'


def test_three_and_nine_component_box_order_is_explicit():
    d=chematic.parse_gro(frame(box='1 2 3 0 0 .4 0 .5 .6'))
    assert d['box_vectors_nm']==[[1,0,0],[.4,2,0],[.5,.6,3]]
    assert len(d['box_values_nm'])==9
    d=chematic.parse_gro(frame())
    assert d['box_vectors_nm']==[[1,0,0],[0,2,0],[0,0,3]]
    assert d['box_values_nm']==[1,2,3]


def test_empty_frame_and_zero_box_are_not_replaced():
    d=chematic.parse_gro('\n0\n0 0 0\n',elements=[])
    assert d['title']=='' and d['atoms']==[] and d['box_vectors_nm']==[[0,0,0]]*3


def test_full_input_crlf_title_and_precision_are_preserved():
    text=frame().replace('synthetic','構造').replace('\n','\r\n')+'\r\n'
    d=chematic.parse_gro(text)
    assert d['title']=='構造 t= 2.5' and d['source_text']==text


@pytest.mark.parametrize('tail',['another frame\n',frame()])
def test_additional_frames_and_trailing_data_fail(tail):
    with pytest.raises(ValueError,match='additional frame or trailing data'):
        chematic.parse_gro(frame()+tail)


@pytest.mark.parametrize('box',['1 2','1 2 3 4','NaN 1 1','1 inf 1','1e999 1 1'])
def test_malformed_or_nonfinite_box_fails(box):
    with pytest.raises(ValueError):chematic.parse_gro(frame(box=box))


@pytest.mark.parametrize('source',['','title\n','title\n-1\n','title\n2000000\n0 0 0\n'])
def test_incomplete_or_bad_counts_fail(source):
    with pytest.raises(ValueError):chematic.parse_gro(source)


@pytest.mark.parametrize('kwargs',[{'max_input_bytes':8},{'max_atoms':0},{'max_line_len':8}])
def test_resource_limits(kwargs):
    with pytest.raises(ValueError):chematic.parse_gro(frame(),**kwargs)


def test_partial_velocity_row_and_nonascii_atom_field_fail():
    with pytest.raises(ValueError):chematic.parse_gro(frame(atom(velocity=True)[:-8]))
    with pytest.raises(ValueError):chematic.parse_gro(frame(atom().replace('ALA','原子')))


def test_stub_syntax_and_function_declaration():
    tree=ast.parse((ROOT/'crates/chematic-py/python/chematic/__init__.pyi').read_text(encoding="utf-8"))
    assert any(isinstance(n,ast.FunctionDef) and n.name=='parse_gro' for n in tree.body)


def test_one_decimal_positions_and_two_decimal_velocities():
    d=chematic.parse_gro(frame(atom(precision=1,velocity=True)))
    assert d['atoms'][0]['position_nm']==[-1.2,.5,4.8]
    assert d['atoms'][0]['velocity_nm_per_ps']==[-.12,.25,.5]
