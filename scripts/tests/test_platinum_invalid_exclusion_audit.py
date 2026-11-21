import importlib.util
from pathlib import Path

from rdkit import Chem


spec = importlib.util.spec_from_file_location(
    "platinum_invalid_exclusion_audit",
    Path(__file__).resolve().parents[1] / "benchmarks" / "platinum_invalid_exclusion_audit.py",
)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_relaxed_query_ignores_bond_orders_but_preserves_elements():
    query = audit.relaxed_graph_query(Chem.MolFromSmiles("C=C"))
    for smiles in ("CC", "C=C", "C#C", "c1ccccc1"):
        assert Chem.MolFromSmiles(smiles).HasSubstructMatch(query)
    assert not Chem.MolFromSmiles("NN").HasSubstructMatch(query)


def test_relaxed_ring_query_matches_unsanitized_invalid_graph():
    query = audit.relaxed_graph_query(Chem.MolFromSmiles("c1ccccc1"))
    target = Chem.MolFromSmiles("O[Hg]C1=CC=CC=C1", sanitize=False)
    assert target is not None
    target.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(target)
    assert target.HasSubstructMatch(query)
    assert not Chem.MolFromSmiles("C1CCCCC1O").HasSubstructMatch(
        audit.relaxed_graph_query(Chem.MolFromSmiles("n1ccccc1"))
    )
