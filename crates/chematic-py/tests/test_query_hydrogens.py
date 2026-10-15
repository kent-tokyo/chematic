"""RDKit #9675: folding query H preserves tetrahedral neighbour parity."""
import pytest
import chematic

PATTERNS = [
    "[H][C@@](C)(N)C(=O)O", "[C@@]([H])(C)(N)C(=O)O",
    "[C@](C)([H])(N)C(=O)O", "[C@@](C)(N)([H])C(=O)O",
    "[C@](C)(N)(C(=O)O)[H]", "C[C@]([H])(N)C(=O)O",
    "C[C@@](N)([H])C(=O)O",
]


@pytest.mark.parametrize("pattern", PATTERNS)
def test_chirality_survives_h_position_and_query_compaction(pattern):
    original = chematic.SmartsQuery(pattern)
    folded = original.merge_hydrogens()
    assert folded.atom_count == original.atom_count - 1
    assert folded.source == pattern
    assert repr(folded).endswith(".merge_hydrogens()")
    same = chematic.from_smiles("C[C@H](N)C(=O)O")
    other = chematic.from_smiles("C[C@@H](N)C(=O)O")
    assert folded.matches(same, use_chirality=True)
    assert not folded.matches(other, use_chirality=True)
    assert folded.matches(other)
    assert bool(folded.find_matches(same, use_chirality=True))
    assert folded.find_matches(other, use_chirality=True) == []
    assert folded.matches_many([same, other], use_chirality=True) == [True, False]
    batch = chematic.MoleculeBatch([same, other])
    assert folded.matches_batch(batch, use_chirality=True) == [True, False]


@pytest.mark.parametrize("pattern", ["[2H]C", "[H:7]C", "[H+]C", "[#1;D1]C", "[H]~C", "[H]=C", "[H][H]", "[H]"])
def test_h_with_extra_meaning_is_retained(pattern):
    query = chematic.compile_smarts(pattern)
    assert query.merge_hydrogens().atom_count == query.atom_count


@pytest.mark.parametrize("target", ["[H][C@@](C)(N)C(=O)O", "C[C@]([H])(N)C(=O)O", "[2H][C@@](C)(N)C(=O)O"])
def test_explicit_target_h_matches_correct_configuration(target):
    query = chematic.SmartsQuery(PATTERNS[1]).merge_hydrogens()
    assert query.matches(chematic.from_smiles(target), use_chirality=True)


def test_minimum_h_count_and_idempotence():
    query = chematic.SmartsQuery("C([H])[H]").merge_hydrogens()
    assert query.atom_count == query.merge_hydrogens().atom_count == 1
    assert query.matches_many([chematic.from_smiles(s) for s in ["C", "[CH2]", "[CH]", "[C]"]]) == [True, True, False, False]


@pytest.mark.parametrize("pattern", ["[C@]([H])(N)C", "[C@]([H])([H])(N)C", "[$(C[H])]"])
def test_unsupported_folding_fails_explicitly(pattern):
    query = chematic.SmartsQuery(pattern)
    with pytest.raises(ValueError):
        query.merge_hydrogens()
    assert query.source == pattern
