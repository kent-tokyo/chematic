"""Original analytic geometry; MIT OR Apache-2.0. No proprietary input."""
import math
import pytest
import chematic

E = ['C','N','O','F']
P = [[0.,0.,0.],[1.234567891234567,0.,0.],[1.1,1.33333333333333,0.],[0.1,0.3,1.444444444444444]]
R = [[],[0],[1,0],[2,1,0]]


def write(**kwargs):
    return chematic.write_gzmat(E,P,R,charge=0,multiplicity=1,**kwargs)


def test_rows_precision_and_nonmutation():
    before=repr((E,P,R));rows=chematic.cartesian_to_zmatrix(E,P,R)
    assert rows[1]=={'element':'N','references':[0],'values':[P[1][0]]}
    text=write();values=[float(line.split('=')[1]) for line in text.splitlines() if '=' in line]
    assert values==[v for row in rows for v in row['values']]
    assert 'r2= 1.2345679e0' in write(significant_digits=8)
    assert repr((E,P,R))==before
    rows[1]['references'][0]=99;assert chematic.cartesian_to_zmatrix(E,P,R)[1]['references']==[0]


@pytest.mark.parametrize('value',[True,False,1.0,'8',None,-1,0,18,256])
def test_invalid_precision(value):
    with pytest.raises((ValueError,TypeError,OverflowError)):write(significant_digits=value)


@pytest.mark.parametrize('field',['charge','multiplicity'])
@pytest.mark.parametrize('value',[True,1.0,'1',None,2**100])
def test_state_requires_integers(field,value):
    kw={'charge':0,'multiplicity':1};kw[field]=value
    with pytest.raises((ValueError,TypeError,OverflowError)):chematic.write_gzmat(E,P,R,**kw)


@pytest.mark.parametrize('value',[True,0.0,'0',None,-1,4,2**100])
def test_bad_references(value):
    refs=[[],[value],[1,0],[2,1,0]]
    with pytest.raises((ValueError,TypeError,OverflowError)):chematic.cartesian_to_zmatrix(E,P,refs)


@pytest.mark.parametrize('refs',[[[],[0],[1,0],[2,1,1]],[[],[0],[1,0]],R+[[]],[[0],[0],[1,0],[2,1,0]]])
def test_wrong_reference_shapes(refs):
    with pytest.raises(ValueError):chematic.cartesian_to_zmatrix(E,P,refs)


@pytest.mark.parametrize('symbol',['c','*','X','13C','D',''])
def test_exact_element_symbols(symbol):
    with pytest.raises(ValueError):chematic.cartesian_to_zmatrix([symbol],[[0,0,0]],[[]])


@pytest.mark.parametrize('v',[math.nan,math.inf,-math.inf])
def test_nonfinite(v):
    with pytest.raises(ValueError):chematic.cartesian_to_zmatrix(['He'],[[v,0,0]],[[]])


def test_singular_frame_and_output_headers():
    with pytest.raises(ValueError):chematic.cartesian_to_zmatrix(['C']*3,[[0,0,0],[1,0,0],[2,0,0]],[[],[0],[1,0]])
    with pytest.raises(ValueError):write(title='name\n0 1')
    with pytest.raises(ValueError):write(route='HF')
    assert '# HF/STO-3G' in write(route='# HF/STO-3G')
    assert '\n-1 2\n' in chematic.write_gzmat(E,P,R,charge=-1,multiplicity=2)


def test_reference_iterators_are_bounded():
    import itertools
    with pytest.raises(ValueError):chematic.cartesian_to_zmatrix(E,P,itertools.repeat([]))
    with pytest.raises(ValueError):chematic.cartesian_to_zmatrix(E,P,[[],itertools.repeat(0),[1,0],[2,1,0]])


def test_single_atom_and_tiny_bond():
    text=chematic.write_gzmat(['He'],[[4.,5.,6.]],[[]],charge=0,multiplicity=1)
    assert 'He\n' in text and 'Variables:' not in text
    text=chematic.write_gzmat(['H','H'],[[0.,0.,0.],[0.,0.,1e-200]],[[],[0]],charge=0,multiplicity=1)
    assert float(text.split('r2= ')[1].strip()) == 1e-200
