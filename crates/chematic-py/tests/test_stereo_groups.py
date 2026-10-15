"""Checked atom/bond enhanced groups, including correlated V3000 export."""
import ast
import pickle
from pathlib import Path
import pytest
import chematic


@pytest.mark.parametrize('kind,number,token', [('abs',0,'ABS'),('and',7,'RAC7'),('or',7,'REL7')])
@pytest.mark.parametrize('smiles,atoms,bonds', [('C/C=C/C=C/C',[],[1,3]),('C/C=C/[C@H](O)CC',[3],[1])])
def test_groups_round_trip_without_dropping_bonds(kind,number,token,smiles,atoms,bonds):
    source=chematic.from_smiles(smiles)
    coords=source.rdkit_2d_coords()
    expected=[(kind,number,atoms,bonds)]
    grouped=source.with_stereo_groups(expected)
    assert source.stereo_groups==[]
    assert grouped.stereo_groups==expected
    assert grouped.atom_table==source.atom_table and grouped.bond_table==source.bond_table
    block=grouped.to_mol_v3000(coords,strict=True)
    assert f'MDLV30/STEB{token} BONDS=' in block
    if atoms:assert f'MDLV30/STE{token} ATOMS=' in block
    reread=chematic.from_mol_v3000(block)
    assert reread.stereo_groups==expected
    assert reread.smiles==source.smiles
    assert pickle.loads(pickle.dumps(grouped)).stereo_groups==expected


def test_inspection_is_a_detached_copy_and_empty_list_clears():
    source=chematic.from_smiles('C/C=C/[C@H](O)CC').with_stereo_groups([('and',1,[3],[1])])
    groups=source.stereo_groups
    groups[0][2].append(0);groups[0][3].clear();groups.clear()
    assert source.stereo_groups==[('and',1,[3],[1])]
    assert source.with_stereo_groups([]).stereo_groups==[]
    assert source.stereo_groups==[('and',1,[3],[1])]


def test_typed_properties_are_copied_independently():
    source=chematic.from_smiles('C/C=C/[C@H](O)CC')
    source.SetName('source');source.SetIntProp('int',-7)
    source.SetUnsignedProp('unsigned',9);source.SetDoubleProp('double',1.25)
    source.SetBoolProp('bool',True);source.SetProp('text','original')
    out=source.with_stereo_groups([('and',1,[3],[1])])
    assert out.GetPropsAsDict()==source.GetPropsAsDict()
    assert out.GetIntProp('int')==-7 and out.GetUnsignedProp('unsigned')==9
    assert out.GetDoubleProp('double')==1.25 and out.GetBoolProp('bool') is True
    out.SetName('copy');out.SetIntProp('int',3)
    assert source.GetName()=='source' and source.GetIntProp('int')==-7


@pytest.mark.parametrize('groups,error', [
    ([('and',0,[],[1])],'positive'),([('or',0,[3],[])],'positive'),
    ([('abs',1,[3],[])],'must be 0'),([('AND',1,[3],[])],'kind'),
    ([('and',1,[],[])],'no atom or bond'),
    ([('and',1,[7],[])],'out-of-range atom'),([('and',1,[],[6])],'out-of-range bond'),
    ([('and',1,[3],[]),('and',1,[],[1])],'repeat the same'),
    ([('and',1,[3],[]),('or',1,[3],[])],'belongs to stereo groups'),
    ([('and',1,[],[1]),('abs',0,[],[1])],'belongs to stereo groups'),
])
def test_bad_grouping_fails_atomically(groups,error):
    source=chematic.from_smiles('C/C=C/[C@H](O)CC').with_stereo_groups([('or',2,[3],[1])])
    with pytest.raises(ValueError,match=error):source.with_stereo_groups(groups)
    assert source.stereo_groups==[('or',2,[3],[1])]


@pytest.mark.parametrize('value',[-1,2**32,True,False,1.5,'1',None])
@pytest.mark.parametrize('position',['number','atom','bond'])
def test_indices_and_numbers_require_unsigned_integers(value,position):
    kind,number,atoms,bonds='and',1,[3],[1]
    if position=='number':number=value
    elif position=='atom':atoms=[value]
    else:bonds=[value]
    with pytest.raises(ValueError):
        chematic.from_smiles('C/C=C/[C@H](O)CC').with_stereo_groups([(kind,number,atoms,bonds)])


def test_duplicates_within_a_group_deduplicate_in_first_seen_order():
    source=chematic.from_smiles('C/C=C/C=C/C')
    assert source.with_stereo_groups([('and',3,[],[3,1,3])]).stereo_groups==[('and',3,[],[3,1])]


def test_memberships_follow_native_atom_permutation_bond_indices_stay_fixed():
    source=chematic.from_smiles('C/C=C/[C@H](O)CC').with_stereo_groups([('and',1,[3],[1])])
    order=[3,0,1,2,4,5,6]
    out=source.renumber_atoms(order)
    assert out.stereo_groups==[('and',1,[0],[1])]
    assert source.stereo_groups==[('and',1,[3],[1])]


def test_and_or_numbers_have_distinct_namespaces():
    source=chematic.from_smiles('C/C=C/C=C/C')
    expected=[('and',1,[],[1]),('or',1,[],[3])]
    result=source.with_stereo_groups(expected)
    assert result.stereo_groups==expected
    assert chematic.from_mol_v3000(result.to_mol_v3000(source.rdkit_2d_coords(),strict=True)).stereo_groups==expected


def test_old_cx_compatibility_api_still_refuses_enhanced_groups():
    source=chematic.from_smiles('C/C=C/C=C/C').with_stereo_groups([('and',1,[],[1,3])])
    with pytest.raises(ValueError):source.rdkit_cx_smiles()


def test_shipped_stub_parses_and_declares_the_new_surface():
    path=Path(__file__).resolve().parents[1]/'python/chematic/__init__.pyi'
    tree=ast.parse(path.read_text(encoding="utf-8"))
    mol=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='Mol')
    assert {'stereo_groups','with_stereo_groups'} <= {n.name for n in mol.body if isinstance(n,ast.FunctionDef)}
