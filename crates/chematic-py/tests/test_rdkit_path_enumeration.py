"""Original synthetic MIT OR Apache-2.0 contracts; RDKit 2026.09.1 oracle."""
import ast
import hashlib
import json
from pathlib import Path

import pytest
import chematic

MANIFEST = Path(__file__).resolve().parents[3]/'validation/rdkit_issues/path-enumeration/fixtures.json'


def test_golden_paths_and_fixture_integrity():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert manifest['sha256'] == hashlib.sha256(json.dumps(manifest['cases'], sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    for case in manifest['cases']:
        assert case['sha256'] == hashlib.sha256(case['smiles'].encode()).hexdigest()
    for case in manifest['golden']:
        mol = chematic.from_smiles(case['smiles'])
        assert mol.rdkit_paths_of_length(case['length'], **case['params']) == case['expected']['paths']


@pytest.mark.parametrize('use_bonds', [False, True])
@pytest.mark.parametrize('length', [0, 2**32-1])
def test_zero_and_impossible_lengths_need_no_search_budget(use_bonds, length):
    assert chematic.from_smiles('CCO').rdkit_paths_of_length(length, use_bonds=use_bonds, max_visits=0, max_frontier=0, max_results=0, max_path_indices=0) == []


@pytest.mark.parametrize('root', [3, 999])
def test_outside_roots_return_no_paths(root):
    assert chematic.from_smiles('CCO').rdkit_paths_of_length(1, rooted_at_atom=root, max_visits=0, max_frontier=0) == []


@pytest.mark.parametrize('kwargs', [dict(max_visits=0), dict(max_frontier=2), dict(max_results=1), dict(max_path_indices=2)])
def test_each_budget_raises_instead_of_truncating(kwargs):
    m = chematic.from_smiles('CCO')
    before = m.smiles
    with pytest.raises(ValueError, match='resource limit'):
        m.rdkit_paths_of_length(1, **kwargs)
    assert m.smiles == before
    assert m.rdkit_paths_of_length(1) == [[0], [1]]


def test_exact_budget_boundary_and_no_visits_for_singletons():
    m = chematic.from_smiles('CCO')
    assert m.rdkit_paths_of_length(1, rooted_at_atom=0, max_visits=1, max_frontier=1, max_results=1, max_path_indices=2) == [[0]]
    with pytest.raises(ValueError, match='frontier indices'):
        m.rdkit_paths_of_length(1, rooted_at_atom=0, max_path_indices=1)
    assert m.rdkit_paths_of_length(1, use_bonds=False, max_visits=0, max_results=3, max_path_indices=3) == [[0], [1], [2]]


def test_returned_lists_are_independent_and_source_properties_survive():
    m = chematic.from_smiles('C1CC1')
    m.SetIntProp('path_sentinel', 9710)
    first = m.rdkit_paths_of_length(2)
    first[0].clear()
    first.append([999])
    assert m.rdkit_paths_of_length(2) == [[0, 1], [2, 1], [0, 2]]
    assert m.GetIntProp('path_sentinel') == 9710


def test_hydrogens_and_sanitized_indexing():
    m = chematic.from_smiles('[H]C([H])([H])[H]')
    assert m.atom_count == 5 and m.rdkit_num_atoms() == 1
    assert m.rdkit_paths_of_length(1, use_bonds=False) == [[0]]
    m = m.rdkit_add_hydrogens()
    assert len(m.rdkit_paths_of_length(1, use_hs=True)) == 4
    assert m.rdkit_paths_of_length(1) == []
    # RDKit builds shortest-mode adjacency from its unfiltered distance matrix.
    assert m.rdkit_paths_of_length(1, only_shortest_paths=True) == m.rdkit_paths_of_length(1, use_hs=True)


def test_empty_model_input_cap_and_type_boundaries():
    assert chematic.RWMol().GetMol().rdkit_paths_of_length(1) == []
    with pytest.raises(ValueError, match='1026 atoms'):
        chematic.from_smiles('.'.join(['[H]C']*513)).rdkit_paths_of_length(0)
    for length in [-1, 2**32]:
        with pytest.raises((ValueError, OverflowError)):
            chematic.from_smiles('C').rdkit_paths_of_length(length)


def test_stub_declares_all_options():
    p = Path(__file__).resolve().parents[1]/'python/chematic/__init__.pyi'
    cls = next(n for n in ast.parse(p.read_text(encoding="utf-8")).body if isinstance(n, ast.ClassDef) and n.name == 'Mol')
    method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == 'rdkit_paths_of_length')
    assert {a.arg for a in method.args.kwonlyargs} == {'use_bonds', 'use_hs', 'rooted_at_atom', 'only_shortest_paths', 'max_visits', 'max_frontier', 'max_results', 'max_path_indices'}
