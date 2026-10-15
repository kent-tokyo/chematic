"""Original MIT OR Apache-2.0 proof-backed relative-group reduction tests."""
import pytest
import chematic


def source(s='C[C@H](O)[C@H](O)C',atoms=(1,3),bonds=(),kind='and'):
    return chematic.from_smiles(s).with_stereo_groups([(kind,1,list(atoms),list(bonds))])


@pytest.mark.parametrize('kind',['and','or'])
def test_atom_only_and_mixed_meso_reduction(kind):
    for m in [source(kind=kind),source(r'C/C=C/[C@H](C)/C=C\C',(3,),(1,5),kind)]:
        r=m.reduce_meso_stereo_group(0)
        assert r['reduced'] is True and r['checked_pairs']==1
        assert r['molecule'].stereo_groups==[('abs',0,m.stereo_groups[0][2],m.stereo_groups[0][3])]
        assert len(r['witnesses'])==1 and sorted(r['witnesses'][0][1])==list(range(m.atom_count))
        assert r['molecule'].smiles==m.smiles


@pytest.mark.parametrize('s,atoms,bonds',[
    ('C[C@H](O)[C@@H](O)C',(1,3),()),
    (r'C/C=C/[C@H](C)/C=C\C',(3,),()),
    ('C/C=C/C=C/C',(),(1,3)),
    ('[13CH3][C@H](O)[C@H](O)C',(1,3),()),
    ('[CH3:11][C@H](O)[C@H](O)[CH3:12]',(1,3),()),
])
def test_conservative_retention(s,atoms,bonds):
    m=source(s,atoms,bonds);r=m.reduce_meso_stereo_group(0)
    assert not r['reduced'] and not r['witnesses']
    assert r['molecule'].stereo_groups==m.stereo_groups and r['molecule'].smiles==m.smiles


def test_abs_merge_and_independent_relative_contexts():
    m=chematic.from_smiles('C[C@H](O)[C@H](O)C.N[C@H](F)Cl.[C@H](F)(Cl)Br').with_stereo_groups([('abs',0,[7],[]),('and',8,[1,3],[]),('or',9,[10],[])])
    r=m.reduce_meso_stereo_group(1)
    assert r['reduced'] and r['checked_pairs']==2
    assert r['molecule'].stereo_groups==[('abs',0,[7,1,3],[]),('or',9,[10],[])]
    assert [mask for mask,_ in r['witnesses']]==[0,2]


def test_properties_source_and_returned_lists_are_independent():
    m=source();m.SetName('meso trial');m.SetIntProp('run',7)
    before=m.stereo_groups;r=m.reduce_meso_stereo_group(0);out=r['molecule']
    assert out.GetName()=='meso trial' and out.GetIntProp('run')==7
    out.SetIntProp('run',99);r['witnesses'][0][1][0]=999
    assert m.GetIntProp('run')==7 and m.stereo_groups==before
    assert source().reduce_meso_stereo_group(0)['witnesses'][0][1][0]!=999


@pytest.mark.parametrize('index',[True,False,-1,1,2.0,'0',2**40,None])
def test_invalid_group_indices(index):
    with pytest.raises(ValueError):source().reduce_meso_stereo_group(index)


@pytest.mark.parametrize('limit',[0,10000001,-1,True,False,2.0,'100',2**40,None])
def test_invalid_visit_limits(limit):
    with pytest.raises(ValueError):source().reduce_meso_stereo_group(0,max_visits=limit)


def test_visit_budget_is_not_a_negative_result():
    m=source();before=m.stereo_groups
    with pytest.raises(ValueError,match='exhausted'):m.reduce_meso_stereo_group(0,max_visits=1)
    assert m.stereo_groups==before


@pytest.mark.parametrize('s,atoms,bonds',[
    ('C[C@H](O)[C@H](O)C.FC(Cl)(Br)I',(1,3),()),
    ('C[C@H](O)[C@H](O)C.Cl[Pt@SP1](Cl)(N)N',(1,3),()),
    ('C[C@H](O)[C@H](O)C.C~C',(1,3),()),
    ('FC(Cl)(Br)I',(1,),()),
    ('C1/C=C/CCCC1',(),(1,)),
])
def test_unsupported_inputs_fail_closed(s,atoms,bonds):
    with pytest.raises(ValueError):source(s,atoms,bonds).reduce_meso_stereo_group(0)


def test_v3000_output_keeps_abs_atom_and_bond_members_and_drawn_geometry():
    base=chematic.from_smiles(r'C/C=C/[C@H](C)/C=C\C');coords=base.rdkit_2d_coords()
    m=base.with_stereo_groups([('and',1,[3],[1,5])]);out=m.reduce_meso_stereo_group(0)['molecule']
    block=out.to_mol_v3000(coords,strict=True)
    back=chematic.from_mol_v3000(block)
    assert back.stereo_groups==[('abs',0,[3],[1,5])]
    assert back.smiles==out.smiles


def test_explicit_hydrogen_frames_can_supply_a_real_witness():
    m=source('[H][C@](C)(O)[C@@]([H])(O)C',(1,4))
    r=m.reduce_meso_stereo_group(0)
    assert r['reduced'] and len(r['witnesses'][0][1])==8


def test_substituent_reference_choices_do_not_confuse_e_and_z():
    m=source('F/C(Cl)=C(/F)Cl',(),(2,))
    assert not m.reduce_meso_stereo_group(0)['reduced']
