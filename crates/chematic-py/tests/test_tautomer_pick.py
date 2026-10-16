"""Checked supplied-candidate selection motivated by RDKit #9601."""

import gc
import itertools

import pytest

import chematic as c


@pytest.mark.parametrize("factory", [list, tuple, iter])
def test_supplied_candidates_only(factory):
    enol, keto = [c.from_smiles(s) for s in ["C=C(C)O", "CC(C)=O"]]
    assert c.rdkit_pick_tautomer(factory([enol, keto])).rdkit_smiles == "CC(C)=O"
    assert c.rdkit_pick_tautomer(factory([enol])).rdkit_smiles == "C=C(C)O"


@pytest.mark.parametrize("items", [[], (), iter(())])
def test_empty(items):
    with pytest.raises(ValueError, match="at least one"):
        c.rdkit_pick_tautomer(items)


@pytest.mark.parametrize("limit", [-1, 0, 4097])
def test_invalid_limits_do_not_consume(limit):
    def forbidden():
        raise AssertionError("must not consume")
        yield
    with pytest.raises(ValueError, match="max_candidates"):
        c.rdkit_pick_tautomer(forbidden(), max_candidates=limit)


@pytest.mark.parametrize("value", [None, 1, object()])
def test_non_iterable(value):
    with pytest.raises(TypeError):
        c.rdkit_pick_tautomer(value)


@pytest.mark.parametrize("value", [None, "CCO", b"CCO", 1, object()])
def test_wrong_item(value):
    with pytest.raises(TypeError, match="candidate 1 must be a Mol"):
        c.rdkit_pick_tautomer([c.from_smiles("CCO"), value])


def test_iterator_exception_propagates():
    def broken():
        yield c.from_smiles("CCO")
        raise RuntimeError("source failure")
    with pytest.raises(RuntimeError, match="source failure"):
        c.rdkit_pick_tautomer(broken())


def test_bounded_infinite_source():
    consumed = []
    def infinite():
        for i in itertools.count():
            consumed.append(i)
            yield c.from_smiles("CCO")
    with pytest.raises(ValueError, match="max_candidates"):
        c.rdkit_pick_tautomer(infinite(), max_candidates=2)
    assert consumed == [0, 1, 2]


def test_exact_limit_and_duplicate_tie_copy():
    a, b = c.from_smiles("CCO"), c.from_smiles("OCC")
    a.SetProp("source", "first")
    b.SetProp("source", "second")
    selected = c.rdkit_pick_tautomer([a, b], max_candidates=2)
    assert selected is not a
    assert selected.GetProp("source") == "first"
    selected.SetProp("source", "changed")
    assert a.GetProp("source") == "first"


def test_generator_temporaries_survive_collection():
    def generated():
        for i, smiles in enumerate(["C=C(C)O", "CC(C)=O", "CCO"]):
            gc.collect()
            m = c.from_smiles(smiles)
            m.SetProp("row", str(i))
            yield m
    selected = c.rdkit_pick_tautomer(generated())
    gc.collect()
    assert selected.rdkit_smiles == "CC(C)=O"
    assert selected.GetProp("row") == "1"


@pytest.mark.parametrize("order", [["CCC", "CC(C)O"], ["CC(C)O", "CCC"]])
def test_lexical_tie(order):
    assert c.rdkit_pick_tautomer(map(c.from_smiles, order)).rdkit_smiles == "CC(C)O"


def test_unsupported_candidate_not_ignored():
    with pytest.raises(ValueError, match="candidate 1"):
        c.rdkit_pick_tautomer([c.from_smiles("c1ccccc1"), c.from_smiles("[CH5]")])


def test_empty_molecule_is_valid():
    assert c.rdkit_pick_tautomer([c.RWMol().GetMol()]).rdkit_smiles == ""


def test_atom_limit():
    with pytest.raises(ValueError, match="1024 atoms"):
        c.rdkit_pick_tautomer([c.from_smiles("C" * 1025)])


def test_typed_properties_and_source_stereo_are_preserved():
    source = c.from_smiles("[13CH3:9][C@H:7](F)Cl")
    source.SetIntProp("number", 42)
    source.SetBoolProp("flag", True)
    source.SetDoubleProp("measurement", 1.25)
    source.SetProp("_Name", "selected source")
    before = source.rdkit_smiles
    chosen = c.rdkit_pick_tautomer([source])
    assert chosen.rdkit_smiles == before == source.rdkit_smiles
    assert chosen.GetIntProp("number") == 42
    assert chosen.GetBoolProp("flag") is True
    assert chosen.GetDoubleProp("measurement") == 1.25
    assert chosen.GetProp("_Name") == "selected source"
    chosen.ClearProp("number")
    assert source.GetIntProp("number") == 42


def test_lazy_sequence_temporary_ownership():
    class Lazy:
        def __getitem__(self, i):
            if i == 3:
                raise IndexError
            m = c.from_smiles(["C=C(C)O", "CC(C)=O", "CCO"][i])
            m.SetIntProp("row", i)
            return m
    chosen = c.rdkit_pick_tautomer(Lazy())
    assert chosen.GetIntProp("row") == 1


def test_default_budget_and_hard_maximum():
    mol = c.from_smiles("C")
    with pytest.raises(ValueError, match="max_candidates"):
        c.rdkit_pick_tautomer(itertools.repeat(mol))
    assert c.rdkit_pick_tautomer([mol] * 4096, max_candidates=4096).rdkit_smiles == "C"


def test_source_graph_is_not_cleaned_or_transformed():
    # RDKit's writer drops this non-stereogenic source tag. Selection itself
    # must preserve it, rather than silently replacing the original graph.
    source = c.from_smiles("C[C@H](C)C")
    before = source.smiles
    assert c.rdkit_pick_tautomer([source]).smiles == before


@pytest.mark.parametrize("limit", ["2", 2.5, None])
def test_limit_type(limit):
    with pytest.raises(TypeError):
        c.rdkit_pick_tautomer([c.from_smiles("C")], max_candidates=limit)


def test_limit_integer_overflow():
    with pytest.raises(OverflowError):
        c.rdkit_pick_tautomer([c.from_smiles("C")], max_candidates=1 << 100)
