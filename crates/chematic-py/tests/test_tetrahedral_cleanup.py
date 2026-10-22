import json
import random
from pathlib import Path
import pytest
import chematic

CASES=json.loads((Path(__file__).resolve().parents[3]/'validation/rdkit_issues/tetrahedral-cleanup/cases.json').read_text(encoding="utf-8"))['cases']

def mirror(s):
    return s.replace('@@','TEMP').replace('@','@@').replace('TEMP','@')

@pytest.mark.parametrize('case',CASES,ids=lambda c:c['id'])
@pytest.mark.parametrize('mirrored',[False,True])
def test_proof_decisions_and_atom_permutations(case,mirrored):
    text=mirror(case['smiles']) if mirrored else case['smiles']
    source=chematic.from_smiles(text)
    for seed in [None,0,11,53]:
        perm=list(range(source.atom_count))
        if seed is not None:random.Random(seed).shuffle(perm)
        m=source.renumber_atoms(perm);a=perm.index(case['atom_index']);before=m.smiles
        report=m.remove_redundant_tetrahedral_stereo(a)
        assert report['removed']==case['removed']
        assert report['atom_index']==a
        assert report['search_visits']>0
        result=report['molecule']
        assert result.atom_count==m.atom_count
        assert result.bond_table == m.bond_table
        assert result.formula==m.formula
        assert m.smiles==before
        if report['removed']:
            assert report['witness'][a]==a
            assert sorted(report['witness'])==list(range(m.atom_count))
        else:
            assert report['witness'] is None
            assert result.smiles==before

@pytest.mark.parametrize('value',[True,False,-1,1.5,'1',None,2**32])
def test_target_indices_are_strict(value):
    with pytest.raises(ValueError):chematic.from_smiles('F[C@H](F)Cl').remove_redundant_tetrahedral_stereo(value)

@pytest.mark.parametrize('value',[True,False,-1,0,1,1.5,'1000',None,10000001,2**32])
def test_bad_or_exhausted_budgets_have_no_partial_result(value):
    m=chematic.from_smiles('F[C@H](F)Cl');before=m.smiles
    with pytest.raises(ValueError):m.remove_redundant_tetrahedral_stereo(1,max_visits=value)
    assert m.smiles==before

@pytest.mark.parametrize('smiles,idx',[('F[C@H](F)C(Cl)(Br)I',1),('F[C@H](F)*',1),('FC=[C@AL1]=C(Cl)Br',2),('F[C@H](F)N->[Cu]',1),('C/[C@H](C)F',1)])
def test_unsupported_global_stereo_is_not_silently_ignored(smiles,idx):
    m=chematic.from_smiles(smiles)
    with pytest.raises(ValueError):m.remove_redundant_tetrahedral_stereo(idx)


def test_source_properties_are_copied_independently():
    m=chematic.from_smiles('C[C@H](C)C');m.SetProp('name','source');m.SetIntProp('count',7)
    out=m.remove_redundant_tetrahedral_stereo(1)['molecule']
    assert out.GetProp('name')=='source';assert out.GetIntProp('count')==7
    out.SetProp('name','result');assert m.GetProp('name')=='source'


def test_enhanced_groups_and_large_graphs_are_rejected():
    m=chematic.from_smiles('F[C@H](F)Cl').with_stereo_groups([('and',1,[1],[])])
    with pytest.raises(ValueError,match='enhanced'):m.remove_redundant_tetrahedral_stereo(1)
    with pytest.raises(ValueError,match='64 atoms'):chematic.from_smiles('C'*65).remove_redundant_tetrahedral_stereo(1)
    with pytest.raises(ValueError,match='outside'):chematic.from_smiles('C').remove_redundant_tetrahedral_stereo(99)
