"""Original minimal inputs; pinned RDKit oracle checked separately."""
import pytest
import chematic


def test_existing_default_and_detached_seed_state():
    m=chematic.from_smiles('CC(C)O');m.SetIntProp('tag',7)
    seeds=[10,20,30,40,99];roots=[0,0]
    out=m.rdkit_layered_config(from_atoms=roots,atom_counts=seeds)
    assert out['bits']==m.rdkit_layered_config(from_atoms=[0])['bits']
    assert out['atom_counts']==[18,28,34,44,99]
    assert seeds==[10,20,30,40,99] and roots==[0,0] and m.GetIntProp('tag')==7
    assert m.rdkit_layered_config()['bits']==bytes(m.rdkit_layered_fp())
    assert isinstance(out['bits'],bytes)


@pytest.mark.parametrize('field',['layer_flags','min_path','max_path','nbits','max_visits','max_paths'])
@pytest.mark.parametrize('value',[True,1.0,'1',-1,None,2**100])
def test_strict_integer_arguments(field,value):
    with pytest.raises((ValueError,TypeError,OverflowError)):chematic.from_smiles('CCO').rdkit_layered_config(**{field:value})


@pytest.mark.parametrize('value',[0,1,'true',None])
def test_branched_flag_must_be_bool(value):
    with pytest.raises((ValueError,TypeError)):chematic.from_smiles('CCO').rdkit_layered_config(branched_paths=value)


@pytest.mark.parametrize('field',['from_atoms','atom_counts'])
@pytest.mark.parametrize('value',[[True],[1.0],['0'],[-1],[2**40]])
def test_sequence_entries_checked(field,value):
    with pytest.raises((ValueError,TypeError,OverflowError)):chematic.from_smiles('CCO').rdkit_layered_config(**{field:value})


@pytest.mark.parametrize('kwargs,code',[({'nbits':0},'fp_size'),({'max_path':13},'path_range'),({'min_path':8,'max_path':7},'path_range'),({'from_atoms':[3]},'root_index'),({'atom_counts':[0]},'atom_counts'),({'set_only_bits':b''},'mask_size'),({'nbits':9,'set_only_bits':bytes([0,2])},'mask_padding'),({'branched_paths':False},'unrooted_linear_unsupported'),({'max_visits':1},'visit_limit'),({'max_paths':1},'path_limit'),({'branched_paths':False,'from_atoms':[1],'max_paths':1},'frontier_limit')])
def test_typed_failures(kwargs,code):
    with pytest.raises(chematic.ChematicInputError) as error: chematic.from_smiles('CCO').rdkit_layered_config(**kwargs)
    assert error.value.code==code and error.value.format=='rdkit_layered_config'


def test_masking_counts_padding_and_counter_wrap():
    m=chematic.from_smiles('CCO')
    assert m.rdkit_layered_config(nbits=1)['bits']==b'\x01'
    assert m.rdkit_layered_config(nbits=1,set_only_bits=b'\0')['atom_counts']==[0,0,0]
    out=m.rdkit_layered_config(atom_counts=[2**32-1]*3)
    assert out['atom_counts']==[1,2,1]
    assert len(m.rdkit_layered_config(nbits=9)['bits'])==2
    with pytest.raises(ValueError):m.rdkit_layered_config(set_only_bits=[0]*256)


def test_rooted_linear_and_bounded_iterators():
    m=chematic.from_smiles('CC(C)O')
    assert m.rdkit_layered_config(branched_paths=False,from_atoms=[1])['atom_counts']==[1,3,1,1]
    import itertools
    with pytest.raises(ValueError):m.rdkit_layered_config(from_atoms=itertools.repeat(0))
    with pytest.raises(ValueError):m.rdkit_layered_config(atom_counts=itertools.repeat(0))
