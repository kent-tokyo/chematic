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
    results = index.search("c1ccccc1", threshold=0.3)
    assert isinstance(results, list)
    assert len(results) >= 1


def test_search_self_similar(index):
    # Benzene should find itself as most similar
    results = index.search("c1ccccc1", threshold=0.5)
    assert len(results) >= 1
    idx_found, sim = results[0]
    assert isinstance(idx_found, int)
    assert 0.0 <= sim <= 1.0
    assert sim >= 0.5


def test_search_with_k(index):
    results = index.search("c1ccccc1", threshold=0.0, k=3)
    assert len(results) <= 3


def test_get_smiles(index):
    smiles = index.get_smiles(0)
    assert isinstance(smiles, str)
    assert len(smiles) > 0


def test_get_smiles_out_of_bounds(index):
    # Consistent with add()/search()/from_smiles(), which all raise ValueError
    # on invalid input rather than returning a sentinel.
    with pytest.raises(ValueError):
        index.get_smiles(9999)


def test_add_and_search():
    idx = chematic.SimilarityIndex()
    idx.add("c1ccccc1")
    results = idx.search("c1ccccc1", threshold=0.0)
    assert len(results) == 1
