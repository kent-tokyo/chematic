"""Atomic full-group alternatives and shared-marker geometry preservation."""
import pytest
import chematic


@pytest.mark.parametrize('smiles,atoms,bonds,expected',[
    ('C/C=C/[C@H](O)CC',[3],[1],r'C/C=C\[C@@H](O)CC'),
    ('C/C=C/C=C/C',[],[1],r'C/C=C\C=C\C'),
    ('C/C=C/C=C/C',[],[3],r'C/C=C/C=C\C'),
    ('C/C=C/C=C/C',[],[1,3],r'C/C=C\C=C/C'),
    ('N[C@@H](C)C(=O)O',[1],[],'N[C@H](C)C(=O)O'),
    ('[H][C@](F)(Cl)Br',[1],[],'[H][C@@](F)(Cl)Br'),
])
@pytest.mark.parametrize('kind',['and','or'])
def test_full_alternative_and_double_inversion(smiles,atoms,bonds,expected,kind):
    source=chematic.from_smiles(smiles).with_stereo_groups([(kind,1,atoms,bonds)])
    before=source.smiles
    out=source.invert_stereo_group(0)
    assert out.smiles==chematic.from_smiles(expected).smiles
    assert out.stereo_groups==source.stereo_groups
    assert out.invert_stereo_group(0).smiles==before
    assert source.smiles==before


def test_other_group_center_and_bond_are_unchanged():
    smiles='C/C=C/[C@H](O)CC.F[C@](Cl)(Br)I'
    source=chematic.from_smiles(smiles).with_stereo_groups([('and',1,[3],[1]),('abs',0,[8],[])])
    expected=chematic.from_smiles(r'C/C=C\[C@@H](O)CC.F[C@](Cl)(Br)I')
    out=source.invert_stereo_group(0)
    assert out.smiles==expected.smiles
    assert out.stereo_groups==source.stereo_groups


def test_typed_properties_are_copied_independently():
    source=chematic.from_smiles('C/C=C/[C@H](O)CC').with_stereo_groups([('and',1,[3],[1])])
    source.SetName('source');source.SetIntProp('int',-7);source.SetDoubleProp('float',1.25)
    source.SetBoolProp('bool',True);source.SetProp('text','original')
    out=source.invert_stereo_group(0)
    assert out.GetName()=='source' and out.GetIntProp('int')==-7
    assert out.GetDoubleProp('float')==1.25 and out.GetBoolProp('bool') is True
    assert out.GetProp('text')=='original'
    out.SetName('copy');out.SetIntProp('int',3)
    assert source.GetName()=='source' and source.GetIntProp('int')==-7


@pytest.mark.parametrize('index',[-1,2**32,True,False,0.5,'0',None,1])
def test_invalid_indices_fail_atomically(index):
    source=chematic.from_smiles('F[C@](Cl)(Br)I').with_stereo_groups([('and',1,[1],[])])
    before=source.smiles
    with pytest.raises(ValueError):source.invert_stereo_group(index)
    assert source.smiles==before


@pytest.mark.parametrize('smiles,kind,atoms,bonds,error',[
    ('F[C@](Cl)(Br)I','abs',[1],[],'absolute'),
    ('FC(Cl)(Br)I','and',[1],[],'specified tetrahedral'),
    ('Cl[Pt@SP1](Cl)(N)N','and',[1],[],'specified tetrahedral'),
    ('CC','and',[],[0],'axes are unsupported'),
    ('CC=CC','and',[],[1],'fully specified'),
    ('C1/C=C/CCCC1','and',[],[1],'ring smaller'),
])
def test_unsupported_configurations_fail_explicitly(smiles,kind,atoms,bonds,error):
    source=chematic.from_smiles(smiles).with_stereo_groups([(kind,0 if kind=='abs' else 1,atoms,bonds)])
    before=source.smiles
    with pytest.raises(ValueError,match=error):source.invert_stereo_group(0)
    assert source.smiles==before


def test_atom_permutations_keep_the_selected_group_meaning():
    source=chematic.from_smiles('C/C=C/[C@H](O)CC').with_stereo_groups([('and',1,[3],[1])])
    expected=source.invert_stereo_group(0).smiles
    for order in [list(range(7)),list(reversed(range(7))),[3,0,1,2,4,5,6]]:
        out=source.renumber_atoms(order).invert_stereo_group(0)
        assert out.smiles==expected
        assert out.stereo_groups==[('and',1,[order.index(3)],[1])]


@pytest.mark.parametrize('kind',['and','or'])
def test_eight_membered_ring_alternative_is_reversible(kind):
    source=chematic.from_smiles('C1/C=C/CCCCC1').with_stereo_groups([(kind,1,[],[1])])
    out=source.invert_stereo_group(0)
    assert out.invert_stereo_group(0).smiles==source.smiles
