"""Reaction-origin comparison must preserve symmetric mappings and null origins."""

import pytest

pytest.importorskip("rdkit")

from rdkit import Chem

from scripts.reaction_atom_provenance_gate import labelled_smiles, rdkit_sets, rust_sets


def test_labelled_product_distinguishes_source_atom_from_created_atom() -> None:
    mol = Chem.MolFromSmiles("CO")
    retained = labelled_smiles(mol, [(0, 0), (0, 1)], [0, 2])
    created = labelled_smiles(mol, [(0, 0), None], [0, 2])
    assert retained != created
    with pytest.raises(ValueError, match="outside the reactant"):
        labelled_smiles(mol, [(0, 2), None], [0, 2])


def test_symmetric_embedding_can_match_graph_but_lose_provenance() -> None:
    case = {
        "smirks": "[C:1][O:2][C:3]>>[C:1][O:2].[C:3]",
        "reactants": ["COC"],
    }
    oracle_graph, oracle_origins, raw = rdkit_sets(case)
    assert len(oracle_graph) == 1
    assert len(oracle_origins) == 2
    assert raw == 2
    candidate_graph, candidate_origins = rust_sets(case, {
        "sets_with_atom_sources": [[
            {"smiles": "CO", "atom_sources": [[0, 0], [0, 1]]},
            {"smiles": "C", "atom_sources": [[0, 2]]},
        ]],
    })
    assert candidate_graph == oracle_graph
    assert len(candidate_origins) == 1
    assert candidate_origins != oracle_origins
