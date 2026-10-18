"""Owned typed score callbacks; all validation precedes user scoring."""

import gc
import itertools
import numpy as np
import pytest
import chematic as C
from chematic import _scored_tautomer as module


def mol(smiles, tag=""):
    result = C.from_smiles(smiles)
    result.SetName(tag)
    return result


def test_custom_score_can_choose_an_enol_and_does_not_enumerate():
    selected = C.rdkit_pick_tautomer_scored(
        (mol(s) for s in ["CC(C)=O", "C=C(C)O"]),
        lambda candidate: -candidate.rdkit_tautomer_score(),
    )
    assert selected.rdkit_smiles == "C=C(C)O"
    assert (
        C.rdkit_pick_tautomer_scored([mol("C=C(C)O")], lambda _: -10).rdkit_smiles
        == "C=C(C)O"
    )


def test_signed_integer_scores_numpy_and_index_limit():
    candidates = [mol("CC(C)=O", "first"), mol("C=C(C)O", "second")]
    selected = C.rdkit_pick_tautomer_scored(
        candidates,
        lambda m: np.int64(2**31 - 1 if m.GetName() == "second" else -(2**31)),
        max_candidates=np.int64(2),
    )
    assert selected.GetName() == "second"


def test_lexical_ties_and_first_identical_occurrence():
    a, b = mol("CC(C)=O", "first"), mol("C=C(C)O", "second")
    assert C.rdkit_pick_tautomer_scored([a, b], lambda _: 0).rdkit_smiles == min(
        a.rdkit_smiles, b.rdkit_smiles
    )
    b = mol("CC(C)=O", "second")
    assert C.rdkit_pick_tautomer_scored([a, b], lambda _: 0).GetName() == "first"


def test_callback_views_cannot_change_returned_metadata_or_originals():
    sources = [mol("CC(C)=O", "keto"), mol("C=C(C)O", "enol")]
    sources[1].SetIntProp("typed", 7)
    sources[1].SetBoolProp("accepted", True)
    views = []

    def score(view):
        original_name = view.GetName()
        views.append(view)
        view.SetName("callback edit")
        view.SetIntProp("typed", -99)
        return int(original_name == "enol")

    winner = C.rdkit_pick_tautomer_scored(sources, score)
    assert [s.GetName() for s in sources] == ["keto", "enol"]
    assert winner.GetName() == "enol"
    assert winner.GetIntProp("typed") == 7
    assert winner.GetBoolProp("accepted") is True
    assert winner is not sources[1] and all(winner is not view for view in views)
    winner.SetName("return edit")
    assert sources[1].GetName() == "enol"
    assert views[1].GetName() == "callback edit"


def test_generator_reuse_is_snapshotted_before_advancing():
    source = mol("CC(C)=O")

    def candidates():
        for index in range(3):
            source.SetIntProp("index", index)
            yield source
            gc.collect()

    seen = []

    def score(view):
        seen.append(view.GetIntProp("index"))
        return -view.GetIntProp("index")

    selected = C.rdkit_pick_tautomer_scored(candidates(), score)
    assert seen == [0, 1, 2]
    assert selected.GetIntProp("index") == 0
    assert source.GetIntProp("index") == 2


def test_mutating_other_sources_during_callback_does_not_change_snapshots():
    sources = [mol("CCO", "a"), mol("CCN", "b")]
    seen = []

    def score(view):
        sources[1].SetName("changed externally")
        seen.append(view.GetName())
        return 1 if view.GetName() == "b" else 0

    chosen = C.rdkit_pick_tautomer_scored(sources, score)
    assert seen == ["a", "b"] and chosen.GetName() == "b"


@pytest.mark.parametrize("candidate", [None, "CCO", 7])
def test_invalid_candidate_is_rejected_before_any_callback(candidate):
    calls = []
    with pytest.raises(TypeError, match="candidate 1"):
        C.rdkit_pick_tautomer_scored(
            [mol("CCO"), candidate], lambda m: calls.append(m) or 1
        )
    assert not calls


@pytest.mark.parametrize(
    "value,error",
    [
        (0, ValueError),
        (-1, ValueError),
        (4097, ValueError),
        (True, TypeError),
        (np.bool_(True), TypeError),
        (1.0, TypeError),
        (None, TypeError),
    ],
)
def test_limits_validate_before_input_is_touched(value, error):
    def broken():
        raise AssertionError("input touched")
        yield

    with pytest.raises(error):
        C.rdkit_pick_tautomer_scored(broken(), lambda _: 1, max_candidates=value)


@pytest.mark.parametrize("value", [None, False, 1, "not callable"])
def test_callback_type_is_checked_before_iterating(value):
    with pytest.raises(TypeError, match="callable"):
        C.rdkit_pick_tautomer_scored(None, value)


@pytest.mark.parametrize(
    "score,error",
    [
        (None, TypeError),
        (True, TypeError),
        (np.bool_(True), TypeError),
        (1.5, TypeError),
        ("1", TypeError),
        (2**31, ValueError),
        (-(2**31) - 1, ValueError),
    ],
)
def test_bad_score_never_returns_a_partial_winner(score, error):
    source = mol("CCO", "original")
    with pytest.raises(error):
        C.rdkit_pick_tautomer_scored([source], lambda _: score)
    assert source.GetName() == "original"


def test_empty_and_noniterable_sources():
    with pytest.raises(ValueError, match="one candidate"):
        C.rdkit_pick_tautomer_scored([], lambda _: 1)
    with pytest.raises(TypeError):
        C.rdkit_pick_tautomer_scored(None, lambda _: 1)


def test_iterator_exception_propagates_and_callbacks_have_not_started():
    failure = RuntimeError("iterator failure")
    calls = []

    def source():
        yield mol("CCO")
        raise failure

    with pytest.raises(RuntimeError) as error:
        C.rdkit_pick_tautomer_scored(source(), lambda m: calls.append(m) or 1)
    assert error.value is failure and not calls


def test_callback_exception_propagates_unchanged():
    failure = RuntimeError("callback failure")

    def score(_):
        raise failure

    with pytest.raises(RuntimeError) as error:
        C.rdkit_pick_tautomer_scored([mol("CCO")], score)
    assert error.value is failure


def test_count_cap_consumes_only_one_additional_item_and_no_callbacks():
    seen, calls = [], []

    def source():
        for index in itertools.count():
            seen.append(index)
            yield mol("C")

    with pytest.raises(ValueError, match="max_candidates"):
        C.rdkit_pick_tautomer_scored(
            source(), lambda m: calls.append(m) or 1, max_candidates=2
        )
    assert seen == [0, 1, 2] and not calls


def test_atom_cap_precedes_canonicalization_and_callback():
    source = mol(".".join(["C"] * 1025))
    calls = []
    with pytest.raises(ValueError, match="atom/bond"):
        C.rdkit_pick_tautomer_scored([source], lambda m: calls.append(m) or 1)
    assert not calls


def test_bond_and_total_atom_caps_precede_callback(monkeypatch):
    monkeypatch.setattr(module, "_MAX_BONDS", 1)
    with pytest.raises(ValueError, match="atom/bond"):
        C.rdkit_pick_tautomer_scored([mol("CCO")], lambda _: pytest.fail("callback"))
    monkeypatch.setattr(module, "_MAX_BONDS", 4096)
    monkeypatch.setattr(module, "_MAX_TOTAL_ATOMS", 5)
    with pytest.raises(ValueError, match="combined"):
        C.rdkit_pick_tautomer_scored(
            [mol("CCO"), mol("CCN")], lambda _: pytest.fail("callback")
        )


def test_snapshot_byte_caps_precede_callback(monkeypatch):
    source = mol("CCO")
    source.SetProp("large", "x" * (1 << 20))
    with pytest.raises(ValueError, match="snapshot byte"):
        C.rdkit_pick_tautomer_scored([source], lambda _: pytest.fail("callback"))
    source = mol("CCO")
    size = len(source.__reduce__()[1][0])
    monkeypatch.setattr(module, "_MAX_TOTAL_PAYLOAD_BYTES", size)
    with pytest.raises(ValueError, match="combined"):
        C.rdkit_pick_tautomer_scored(
            [source, source], lambda _: pytest.fail("callback")
        )


def test_empty_molecule_is_a_candidate():
    selected = C.rdkit_pick_tautomer_scored([C.RWMol().GetMol()], lambda _: 0)
    assert selected.atom_count == 0


@pytest.mark.parametrize(
    "smiles",
    [
        "N[C@@H](C)C(=O)O",
        "F/C=C\\Cl",
        "[13CH3:7][NH2+]C(=O)[O-]",
        "[Pt@SP2](Cl)(Cl)(N)N",
        "CC=[C@AL1]=CCC",
    ],
)
def test_stereo_isotopes_and_maps_survive_owned_selection(smiles):
    original = mol(smiles)
    selected = C.rdkit_pick_tautomer_scored([original], lambda _: 1)
    assert selected.rdkit_smiles == original.rdkit_smiles
    assert selected.smiles == original.smiles
    assert selected.formula == original.formula


def test_unsupported_canonicalization_prevents_all_callbacks():
    calls = []
    with pytest.raises(ValueError, match="candidate 1"):
        C.rdkit_pick_tautomer_scored(
            [mol("CCO"), mol("[CH5]")], lambda m: calls.append(m) or 1
        )
    assert not calls
