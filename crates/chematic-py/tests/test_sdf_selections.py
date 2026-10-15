"""Original minimal fixtures; see validation/cosmolkit_comparison/sdf_selections."""
from pathlib import Path
import gc
import pytest
import chematic

FIXTURE=Path(__file__).resolve().parents[3]/'validation/cosmolkit_comparison/sdf_selections/minimal.sdf'

@pytest.fixture
def dataset(tmp_path):
    path=tmp_path/'selected.sdf';path.write_bytes(FIXTURE.read_bytes())
    return chematic.SdfDataset.open(str(path)),path


def test_order_duplicates_fields_and_detached_results(dataset):
    d,_=dataset;b=d[[4,0,4,2]]
    assert isinstance(b,chematic.SdfSelectionBatch)
    assert b.source_indices()==[4,0,4,2] and b.valid_mask()==[True]*4
    assert b.to_smiles_list()==['F','C','F','O']
    assert b[-1].get('ID')=='oxygen-id' and b[-1].get('note')=='line1\nline2'
    b[0].mol.SetProp('test','changed');assert not b[2].mol.HasProp('test')
    indices=b.source_indices();indices[0]=99;assert b.source_indices()[0]==4
    records=b.records();records.clear();assert len(b)==4


@pytest.mark.parametrize('key,indices',[(slice(None,None,-1),[4,3,2,1,0]),(slice(4,0,-2),[4,2]),(slice(-100,100,2),[0,2,4]),(slice(3,1),[]),([-1,-5,-1],[4,0,4]),((2,0),[2,0]),([True,False,True,False,True],[0,2,4]),([],[])])
def test_python_selection_semantics(dataset,key,indices):
    d,_=dataset;b=d.select(key)
    assert b.source_indices()==indices
    assert b.valid_mask()==[i!=1 for i in indices]


def test_keep_raise_policies_and_errors(dataset):
    d,_=dataset;b=d.select([2,1,0,1])
    assert b.source_indices()==[2,1,0,1] and b.valid_mask()==[True,False,True,False]
    assert b.to_smiles_list()==['O',None,'C',None] and b[1] is None
    errors=b.errors();assert [(e['position'],e['source_index'],e['code']) for e in errors]==[(1,1,'malformed'),(3,1,'malformed')]
    assert all(e['message'] for e in errors)
    for operation in [lambda:d.select([2,1],errors='raise'),lambda:d[[2,1]],lambda:d[1]]:
        with pytest.raises(ValueError):operation()
    with pytest.raises(ValueError):d.select([],errors='skip')


@pytest.mark.parametrize('key',[True,False,0.0,'0',{0},iter([0]),[0,False],[True,0],[True],[-6],[5],[2**100],slice(None,None,0)])
def test_invalid_selectors(dataset,key):
    d,_=dataset
    with pytest.raises((ValueError,TypeError,IndexError,OverflowError)):d[key]


def test_indexing_keeps_scalar_api_and_boundaries(dataset):
    d,_=dataset;assert d[0].name=='carbon' and d[-1].name=='fluorine'
    assert d.get(-1).name=='fluorine';assert d.metadata(-1).title=='fluorine'
    assert 'V3000' in d.raw_record(-1)
    with pytest.raises(IndexError):d[999999999999999999999]
    assert [r.name for r in d.read_range(2)]==['oxygen','nitrogen','fluorine']


def test_batches_survive_wrapper_and_keep_failures(dataset):
    d,_=dataset;stream=d.batches(size=2);del d;gc.collect()
    assert len(stream)==3
    a=next(stream);assert a.source_indices()==[0,1] and a.valid_mask()==[True,False]
    assert len(stream)==2
    rest=list(stream);assert [b.source_indices() for b in rest]==[[2,3],[4]]
    assert len(stream)==0 and list(stream)==[]


def test_failure_and_cancellation_fuse_iterator(dataset):
    d,path=dataset;stream=d.batches(size=2,errors='raise')
    with pytest.raises(ValueError):next(stream)
    assert len(stream)==0 and list(stream)==[]
    stream=d.batches();stream.cancel();assert list(stream)==[]
    stream=d.batches();path.write_bytes(b'')
    with pytest.raises(OSError):next(stream)
    assert list(stream)==[]


@pytest.mark.parametrize('size',[True,0,-1,1.2,'2',10001,2**100])
def test_batch_size_limits(dataset,size):
    d,_=dataset
    with pytest.raises((ValueError,TypeError,OverflowError)):d.batches(size=size)


@pytest.mark.parametrize('limit',[True,0,-1,'5',1.2,2**100])
def test_byte_limits_validate_arguments(dataset,limit):
    d,_=dataset
    with pytest.raises((ValueError,TypeError,OverflowError)):d.select([],max_bytes=limit)
    with pytest.raises((ValueError,TypeError,OverflowError)):d.batches(max_bytes=limit)


def test_limits_checked_before_io(dataset):
    d,path=dataset
    byte_count=d.metadata(0).byte_len
    assert len(d.select([0,0],max_bytes=2*byte_count))==2
    with pytest.raises(ValueError):d.select([0,0],max_bytes=2*byte_count-1)
    path.write_bytes(b'')
    with pytest.raises(IndexError):d.select([0,5])
    with pytest.raises(ValueError):d.select([0],max_bytes=1)
    with pytest.raises(ValueError):d.select([0]*10001)
    with pytest.raises(OSError):d.select([0])


def test_empty_dataset_and_all_invalid(tmp_path):
    path=tmp_path/'empty.sdf';path.write_text('')
    d=chematic.SdfDataset(str(path));assert len(d[:])==0 and list(d.batches())==[]
    path.write_text('broken\n$$$$\n');d=chematic.SdfDataset(str(path));b=d.select([0,0])
    assert b.valid_mask()==[False,False] and b.to_smiles_list()==[None,None]
