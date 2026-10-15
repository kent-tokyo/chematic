"""Original MIT OR Apache-2.0 correlated atom/bond group enumeration tests."""
import ast
from pathlib import Path
import pytest
import chematic


def source():
    return chematic.from_smiles("C/C=C/[C@H](F)C[C@@H](Cl)/C=C/Br").with_stereo_groups([("and",1,[3],[1]),("or",2,[6],[8])])


def test_mask_order_uses_stored_relative_groups():
    m=source(); groups,variants=m.enumerate_stereo_group_alternatives()
    assert groups==[0,1] and len(variants)==4
    assert variants[0].smiles==m.smiles
    assert variants[1].smiles==m.invert_stereo_group(0).smiles
    assert variants[2].smiles==m.invert_stereo_group(1).smiles
    assert variants[3].smiles==m.invert_stereo_group(0).invert_stereo_group(1).smiles
    assert all(v.stereo_groups==m.stereo_groups for v in variants)


def test_typed_properties_are_independent_and_source_is_unchanged():
    m=source();m.SetName("correlated");m.SetIntProp("trial",17)
    before=m.smiles;_,variants=m.enumerate_stereo_group_alternatives()
    assert all(v.GetName()=="correlated" and v.GetIntProp("trial")==17 for v in variants)
    variants[0].SetIntProp("trial",99)
    assert variants[1].GetIntProp("trial")==17 and m.GetIntProp("trial")==17
    assert m.smiles==before


@pytest.mark.parametrize("limit",[0,65,True,False,-1,2.0,"4",2**40,None])
def test_invalid_limits(limit):
    with pytest.raises(ValueError):source().enumerate_stereo_group_alternatives(limit)


@pytest.mark.parametrize("limit",[1,2,3])
def test_budget_refuses_partial_results(limit):
    with pytest.raises(ValueError,match="4 group alternatives"):
        source().enumerate_stereo_group_alternatives(limit)


def test_absolute_groups_do_not_consume_bits():
    m=chematic.from_smiles("F[C@](Cl)(Br)I.N[C@H](F)[C@@H](Cl)Br").with_stereo_groups([("abs",0,[1],[]),("or",2,[6],[]),("and",3,[8],[])])
    groups,variants=m.enumerate_stereo_group_alternatives()
    assert groups==[1,2] and len(variants)==4
    assert all(v.stereo_groups==m.stereo_groups for v in variants)


@pytest.mark.parametrize("s",["CC(F)Cl","F[C@](Cl)(Br)I"])
def test_no_relative_group_is_single_unchanged_assignment(s):
    m=chematic.from_smiles(s);groups,variants=m.enumerate_stereo_group_alternatives(1)
    assert groups==[] and len(variants)==1 and variants[0].smiles==m.smiles


@pytest.mark.parametrize("n",[1,2,3,4,5,6])
def test_every_group_mask_is_returned(n):
    m=chematic.from_smiles(".".join(["[C@H](F)(Cl)Br"]*n)).with_stereo_groups([("and",i+1,[4*i],[]) for i in range(n)])
    groups,variants=m.enumerate_stereo_group_alternatives()
    assert groups==list(range(n)) and len(variants)==2**n
    assert all(v.atom_count==m.atom_count for v in variants)


def test_seven_groups_are_rejected():
    m=chematic.from_smiles(".".join(["[C@H](F)(Cl)Br"]*7)).with_stereo_groups([("and",i+1,[4*i],[]) for i in range(7)])
    with pytest.raises(ValueError,match="7 relative groups"):
        m.enumerate_stereo_group_alternatives()


@pytest.mark.parametrize("s,atoms,bonds",[("FC(Cl)(Br)I",[1],[]),("CC=CC",[],[1]),("Cl[Pt@SP1](Cl)(N)N",[1],[]),("C1/C=C/CCCC1",[],[1]),("CC",[],[0])])
def test_unsupported_relative_members_fail(s,atoms,bonds):
    m=chematic.from_smiles(s).with_stereo_groups([("and",1,atoms,bonds)])
    before=m.smiles
    with pytest.raises(ValueError,match="stereo group 0"):
        m.enumerate_stereo_group_alternatives()
    assert m.smiles==before


def test_symmetry_equivalent_masks_are_not_collapsed():
    m=chematic.from_smiles("C[C@H](O)[C@H](O)C").with_stereo_groups([("and",1,[1,3],[])])
    _,variants=m.enumerate_stereo_group_alternatives()
    assert len(variants)==2


def test_unspecified_centers_are_not_assigned():
    m=chematic.from_smiles("FC(Cl)(Br)I.N[C@H](F)Cl").with_stereo_groups([("or",1,[6],[])])
    _,variants=m.enumerate_stereo_group_alternatives()
    assert len(variants)==2
    assert all(v.atom_count==m.atom_count for v in variants)


def test_type_stub_declares_method_on_mol():
    p=Path(__file__).resolve().parents[1]/"python/chematic/__init__.pyi"
    tree=ast.parse(p.read_text(encoding="utf-8"));mol=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=="Mol")
    assert any(isinstance(n,ast.FunctionDef) and n.name=="enumerate_stereo_group_alternatives" for n in mol.body)
