"""Tests for SimilarityIndex (LSH-based approximate nearest neighbor search)."""
import pytest
import chematic


SMILES_LIST = [
    "c1ccccc1",            # benzene
    "Cc1ccccc1",           # toluene
    "Clc1ccccc1",          # chlorobenzene
    "Oc1ccccc1",           # phenol
    "CCO",                  # ethanol
    "CCCO",                 # propanol
    "CC(=O)O",              # acetic acid
]


@pytest.fixture(scope="module")
def index():
    idx = chematic.SimilarityIndex()
    for smi in SMILES_LIST:
        idx.add(smi)
    return idx


def test_from_smiles_constructor():
    idx = chematic.SimilarityIndex.from_smiles(SMILES_LIST)
    assert idx is not None


def test_search_returns_results(index):
    query = chematic.from_smiles("c1ccccc1")
    results = index.search(query, threshold=0.3)
    assert isinstance(results, list)
    assert len(results) >= 1


def test_search_self_similar(index):
    # Benzene should find itself as most similar
    query = chematic.from_smiles("c1ccccc1")
    results = index.search(query, threshold=0.5)
    assert len(results) >= 1
    idx_found, sim = results[0]
    assert isinstance(idx_found, int)
    assert 0.0 <= sim <= 1.0
    assert sim >= 0.5


def test_search_with_k(index):
    query = chematic.from_smiles("c1ccccc1")
    results = index.search(query, threshold=0.0, k=3)
    assert len(results) <= 3


def test_get_smiles(index):
    smiles = index.get_smiles(0)
    assert isinstance(smiles, str)
    assert len(smiles) > 0


def test_get_smiles_out_of_bounds(index):
    result = index.get_smiles(9999)
    assert result is None


def test_add_and_search():
    idx = chematic.SimilarityIndex()
    idx.add("c1ccccc1")
    query = chematic.from_smiles("c1ccccc1")
    results = idx.search(query, threshold=0.0)
    assert len(results) == 1
