"""Checked Cartesian GJF generation; geometry projection, not complete jobs."""
import itertools
import math
import struct
import pytest
import chematic


def write(elements=("C",), coords=((1.2345678901234567,-2.5e-20,3.75),), **kwargs):
    options = dict(charge=-1, multiplicity=2, units="angstrom")
    options.update(kwargs)
    return chematic.write_gjf_cartesian(elements, coords, **options)


def test_read_rewrite_retains_selected_fields_and_input_is_not_mutated():
    source = "# HF/STO-3G\n\noriginal\n\n-1 2\nC -1 1.2345678901234567 -2.5e-20 3.75\nO 0 -4.5 5.25 -6.75\nH 7.25 -8.5 9.75\n\n"
    g = chematic.parse_gjf_geometry(source)
    assert g["elements"] == ["C", "O", "H"]
    snapshot = repr((g["elements"],g["coords"],g["freeze_codes"]))
    for unit, keyword in [("angstrom","Ang"),("bohr","AU")]:
        text = write(g["elements"],g["coords"],freeze_codes=g["freeze_codes"],units=unit,title="日本語 geometry",route="#p HF/STO-3G")
        assert text.startswith(f"#p HF/STO-3G Units={keyword}\n\n日本語 geometry\n\n-1 2\n")
        back = chematic.parse_gjf_geometry(text)
        for key in ("elements","coords","freeze_codes","charge","multiplicity"):
            assert back[key] == g[key]
    assert repr((g["elements"],g["coords"],g["freeze_codes"])) == snapshot


@pytest.mark.parametrize("value", [0.,-0.,float.fromhex('0x0.0000000000001p-1022'),-float.fromhex('0x0.0000000000001p-1022'),float.fromhex('0x1.fffffffffffffp+1023'),1e-200,1.2345678901234567])
def test_numeric_values_are_preserved_bitwise(value):
    text = write(coords=[[value,-value,0.]])
    actual = chematic.parse_gjf_geometry(text)["coords"][0]
    assert [struct.pack('>d',v) for v in actual] == [struct.pack('>d',v) for v in [value,-value,0.]]


@pytest.mark.parametrize("coords,flags", [([],None),([[0,0,0]]*2,None),([[0,0,0]],[]),([[0,0,0]],[None,None]),([[0,0]],None),([[0,0,0,1]],None)])
def test_shape_mismatch_is_refused(coords,flags):
    with pytest.raises(ValueError): write(coords=coords,freeze_codes=flags)


@pytest.mark.parametrize("value", [math.nan,math.inf,-math.inf])
@pytest.mark.parametrize("axis", [0,1,2])
def test_nonfinite_is_refused(value,axis):
    coords = [[0.,0.,0.]]; coords[0][axis] = value
    with pytest.raises(ValueError): write(coords=coords)


@pytest.mark.parametrize("field", ["charge","multiplicity"])
@pytest.mark.parametrize("value", [True,False,1.,"1",None,2**100])
def test_state_is_strict_integer(field,value):
    with pytest.raises((ValueError,TypeError,OverflowError)): write(**{field:value})


@pytest.mark.parametrize("value", [True,False,1.,"1",2**32,-2**32])
def test_freeze_codes_are_strict_signed_integers(value):
    with pytest.raises((ValueError,TypeError,OverflowError)): write(freeze_codes=[value])


@pytest.mark.parametrize("symbol", ["c","*","X","D","13C","C1","", "a"*10000])
def test_only_canonical_element_symbols(symbol):
    with pytest.raises(ValueError): write(elements=[symbol])


@pytest.mark.parametrize("field,value", [
    ("title",""),("title","a\nb"),("title","a\rb"),("title","a\0b"),("title","a\u2028b"),("title","a\u2029b"),("title","x"*1025),("route","#"+"x"*1024),
    ("route","HF"),("route","## HF"),("route","# HF ! comment"),("route","# Units=AU"),("route","# UN(AU)"),("route","# U=AU"),
    ("units","AU"),("units",""),("multiplicity",0),
])
def test_headers_units_and_state_validation(field,value):
    with pytest.raises(ValueError): write(**{field:value})


def test_boolean_coordinates_are_refused():
    with pytest.raises(ValueError): write(coords=[[True,0,0]])


def test_iterators_are_bounded_before_materialization():
    with pytest.raises(ValueError): write(elements=itertools.repeat("C"))
    with pytest.raises(ValueError): write(coords=itertools.repeat([0,0,0]))
    with pytest.raises(ValueError): write(coords=[itertools.repeat(0.)])
    with pytest.raises(ValueError): write(freeze_codes=itertools.repeat(None))
    with pytest.raises(ValueError): write(elements=["C"]*100001)
    text = write(elements=(x for x in ["C"]),coords=(x for x in [[0,0,0]]),freeze_codes=iter([None]))
    assert chematic.parse_gjf_geometry(text)["freeze_codes"] == [None]


def test_signed_integer_boundaries_and_absent_codes_are_distinct():
    g = chematic.parse_gjf_geometry(write(elements=["C"]*5,coords=[[0,0,0]]*5,freeze_codes=[None,0,-1,-2**31,2**31-1],charge=-2**31,multiplicity=2**32-1))
    assert g["freeze_codes"] == [None,0,-1,-2**31,2**31-1]
    assert (g["charge"],g["multiplicity"]) == (-2**31,2**32-1)
