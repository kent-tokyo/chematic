"""General residue grouping, independent synthetic Open Babel #3028 fixture."""
from pathlib import Path
import itertools
import pytest
import chematic

FIXTURE=Path(__file__).resolve().parents[3]/'validation/openbabel_issues/mmcif/branched-residues.cif'


def source_rows():
    return [s.split() for s in FIXTURE.read_text(encoding="utf-8").splitlines() if s.startswith('HETATM')]


def render(rows):
    return FIXTURE.read_text(encoding="utf-8").split('HETATM')[0]+'\n'.join(' '.join(r) for r in rows)+'\n#\n'


def test_three_sugars_are_distinct_despite_missing_label_sequence():
    text=FIXTURE.read_text(encoding="utf-8")
    d=chematic.parse_mmcif_residues(text)
    assert [r['label_comp_id'] for r in d['residues']]==['NAG','NAG','BMA']
    assert [r['auth_seq_id'] for r in d['residues']]==['1','2','3']
    assert [r['label_seq_id'] for r in d['residues']]==[None,None,None]
    assert d['atom_to_residue']==[0,0,1,1,2,2]
    assert [r['atom_indices'] for r in d['residues']]==[[0,1],[2,3],[4,5]]
    assert [a['res_name'] for a in d['atoms']]==['NAG']*4+['BMA']*2
    # The legacy API keeps its existing per-atom schema and behavior.
    assert 'residues' not in chematic.parse_mmcif(text)


def test_all_720_input_orders_preserve_partition():
    rows=source_rows()
    for order in itertools.permutations(range(6)):
        d=chematic.parse_mmcif_residues(render([rows[i] for i in order]))
        assert len(d['residues'])==3
        for i,old in enumerate(order):
            r=d['residues'][d['atom_to_residue'][i]]
            assert r['auth_seq_id']==str(old//2+1)
            assert i in r['atom_indices']
            assert len(r['atom_indices'])==2


@pytest.mark.parametrize('column,value',[(6,'F'),(7,'2'),(8,'99'),(9,'AB'),(19,'2'),(5,'MAN'),(16,'AUTHOR_NAME')])
def test_complete_namespaces_prevent_collisions(column,value):
    rows=source_rows()[:2]
    rows[1][3]=rows[1][18]='C1'
    rows[1][column]=value
    d=chematic.parse_mmcif_residues(render(rows))
    assert len(d['residues'])==2


def test_alternate_locations_and_full_insertion_code():
    rows=source_rows()[:2]
    rows[1][3]=rows[1][18]='C1'
    rows[0][4]='AA';rows[1][4]='BB'
    rows[0][9]=rows[1][9]='AB'
    d=chematic.parse_mmcif_residues(render(rows))
    assert len(d['residues'])==1
    assert [a['alt_loc'] for a in d['atoms']]==['AA','BB']
    assert [a['icode'] for a in d['atoms']]==['AB','AB']
    assert d['residues'][0]['insertion_code']=='AB'


def test_long_component_and_separate_atom_aliases_are_preserved():
    rows=source_rows()[:1]
    rows[0][5]='A1AJ1';rows[0][16]='AUTHOR_LONG';rows[0][18]='NAME_LONG'
    d=chematic.parse_mmcif_residues(render(rows))
    assert d['residues'][0]['label_comp_id']=='A1AJ1'
    assert d['residues'][0]['auth_comp_id']=='AUTHOR_LONG'
    assert d['atoms'][0]['label_atom_id']=='C1'
    assert d['atoms'][0]['auth_atom_id']=='NAME_LONG'


@pytest.mark.parametrize('columns',[(6,17),(8,15),(5,16),(3,18)])
def test_missing_identity_fails_instead_of_fabricating_zero(columns):
    rows=source_rows()[:1]
    for c in columns:rows[0][c]='?'
    with pytest.raises(ValueError,match='cannot be inferred'):
        chematic.parse_mmcif_residues(render(rows))


def test_duplicate_atom_is_an_error():
    rows=source_rows()[:2];rows[1][3]=rows[1][18]='C1'
    with pytest.raises(ValueError,match='repeats'):
        chematic.parse_mmcif_residues(render(rows))


@pytest.mark.parametrize('kwargs',[{'max_input_bytes':8},{'max_atoms':2},{'max_line_len':8}])
def test_resource_limits(kwargs):
    with pytest.raises(ValueError):chematic.parse_mmcif_residues(FIXTURE.read_text(encoding="utf-8"),**kwargs)


@pytest.mark.parametrize('column,value',[(19,'4294967297'),(19,'0'),(15,'12A'),(10,'NaN'),(2,'Unknown')])
def test_invalid_or_unsupported_values_fail_explicitly(column,value):
    rows=source_rows()[:1];rows[0][column]=value
    with pytest.raises(ValueError):chematic.parse_mmcif_residues(render(rows))


def test_incomplete_row_is_an_error():
    rows=source_rows()[:1];rows[0].pop()
    with pytest.raises(ValueError,match='incomplete'):chematic.parse_mmcif_residues(render(rows))
