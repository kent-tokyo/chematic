"""Detect legacy spiro oscillation without inventing a stereo correction."""

import pytest
import chematic as ch
import chematic._rdkit_smiles_checked as implementation


@pytest.mark.parametrize(
    "source",
    [
        "C[N@+]1(F)CCC2(CC[C@H](Cl)CC2)CC1",
        "C[N@+]1(F)CC[C@]2(CC[C@H](Cl)CC2)CC1",
        "C[C@]1(F)CCC2(CC[C@H](Cl)CC2)CC1",
    ],
)
def test_9391_rejects_oscillating_spiro_output(source):
    mol = ch.from_smiles(source)
    before = mol.to_moljson()
    first = mol.rdkit_smiles
    second = ch.from_smiles(first).rdkit_smiles
    assert first != second
    assert ch.from_smiles(second).rdkit_smiles == first
    with pytest.raises(ch.RdkitSmilesStabilityError) as caught:
        ch.rdkit_smiles_checked(mol)
    error = caught.value
    assert isinstance(error, ValueError)
    assert error.code == "non_idempotent"
    assert error.first_smiles == first
    assert error.roundtrip_smiles == second
    assert mol.to_moljson() == before


@pytest.mark.parametrize(
    "source",
    [
        "CCO",
        "C1=CC=CC=C1",
        "[13CH3]C(=O)[O-].[Na+]",
        "N[C@@H](C)C(=O)O",
        "N[C@H](C)C(=O)O",
        "F/C=C/F",
        "F/C=C\\F",
        "C[C@H]1CC[C@@H](C)CC1",
        "C[C@H]1CC[C@H](C)CC1",
        "C1CCC2(CC1)CCCCC2",
        "C1C2CC3CC1CC(C2)C3",
        "C[C@]1(F)CCC2(CC[C@](I)(Cl)CC2)CC1",
    ],
)
def test_accepts_stable_outputs_without_modifying_input(source):
    mol = ch.from_smiles(source)
    before = mol.to_moljson()
    expected = mol.rdkit_smiles
    assert ch.rdkit_smiles_checked(mol) == expected
    assert ch.rdkit_smiles_checked(ch.from_smiles(expected)) == expected
    assert mol.to_moljson() == before


def test_bounded_orders_are_deterministic_unique_and_broad():
    first = list(implementation._bounded_atom_orders(7, 32))
    second = list(implementation._bounded_atom_orders(7, 32))
    assert first == second
    assert len(first) == len(set(first)) == 32
    assert first[0] == tuple(range(7))
    assert tuple(reversed(range(7))) in first
    assert tuple(range(1, 7)) + (0,) in first
    assert (0, 2, 4, 6, 1, 3, 5) in first


def test_rejects_real_square_planar_permutation_dependence():
    mol = ch.from_smiles("F[Pt@SP1](Cl)(Br)I")
    with pytest.raises(ch.RdkitSmilesStabilityError) as caught:
        ch.rdkit_smiles_checked(mol)
    assert caught.value.code == "permutation_dependent"
    assert caught.value.first_smiles == "[F][Pt@SP1]([Cl])([Br])[I]"
    assert caught.value.roundtrip_smiles == "[F][Pt@SP3]([Cl])([Br])[I]"
    assert caught.value.atom_order == (0, 2, 4, 1, 3)


def test_small_graph_stops_after_all_permutations():
    assert set(implementation._bounded_atom_orders(2, 32)) == {(0, 1), (1, 0)}


def test_rejects_permutation_dependent_output_with_witness(monkeypatch):
    mol = ch.from_smiles("CCO")
    native_write = implementation._write_permuted

    def write(candidate, order):
        if order == (2, 1, 0):
            return "OCC"
        return native_write(candidate, order)

    monkeypatch.setattr(implementation, "_write_permuted", write)
    with pytest.raises(ch.RdkitSmilesStabilityError) as caught:
        ch.rdkit_smiles_checked(mol, max_permutations=2)
    assert caught.value.code == "permutation_dependent"
    assert caught.value.first_smiles == "CCO"
    assert caught.value.roundtrip_smiles == "OCC"
    assert caught.value.atom_order == (2, 1, 0)


def test_permutation_writer_failure_keeps_order_and_cause(monkeypatch):
    failure = RuntimeError("cannot write permutation")

    def fail(candidate, order):
        raise failure

    monkeypatch.setattr(implementation, "_write_permuted", fail)
    with pytest.raises(ch.RdkitSmilesStabilityError) as caught:
        ch.rdkit_smiles_checked(ch.from_smiles("CCO"), max_permutations=2)
    assert caught.value.code == "permutation_write_failed"
    assert caught.value.roundtrip_smiles is None
    assert caught.value.atom_order == (2, 1, 0)
    assert caught.value.__cause__ is failure


@pytest.mark.parametrize("value", [True, 1.5, "32", None])
def test_permutation_limit_requires_integer(value):
    with pytest.raises(TypeError, match="max_permutations"):
        ch.rdkit_smiles_checked(ch.from_smiles("CCO"), max_permutations=value)


@pytest.mark.parametrize("value", [-1, 0, 1, 257])
def test_permutation_limit_is_bounded(value):
    with pytest.raises(ValueError, match="between 2 and 256"):
        ch.rdkit_smiles_checked(ch.from_smiles("CCO"), max_permutations=value)


def test_preserves_distinct_tetrahedral_and_double_bond_isomers():
    pairs = [("N[C@H](C)C(=O)O", "N[C@@H](C)C(=O)O"), ("F/C=C/F", "F/C=C\\F")]
    for left, right in pairs:
        assert ch.rdkit_smiles_checked(ch.from_smiles(left)) != ch.rdkit_smiles_checked(
            ch.from_smiles(right)
        )


@pytest.mark.parametrize("value", [None, "CCO", 42, object()])
def test_requires_native_molecule(value):
    with pytest.raises(TypeError, match="chematic.Mol"):
        ch.rdkit_smiles_checked(value)


def test_initial_writer_failure_propagates():
    with pytest.raises(ValueError, match="sanitization failed") as caught:
        ch.rdkit_smiles_checked(ch.from_smiles("CN(C)(C)(C)C"))
    assert not isinstance(caught.value, ch.RdkitSmilesStabilityError)


@pytest.mark.parametrize(
    "backend_error", [ValueError("parse"), RuntimeError("backend")]
)
def test_unreadable_output_keeps_candidate_and_cause(monkeypatch, backend_error):
    def fail(text):
        raise backend_error

    monkeypatch.setattr(implementation, "from_smiles", fail)
    with pytest.raises(ch.RdkitSmilesStabilityError) as caught:
        ch.rdkit_smiles_checked(ch.from_smiles("CCO"))
    assert caught.value.code == "unreadable_output"
    assert caught.value.first_smiles == "CCO"
    assert caught.value.roundtrip_smiles is None
    assert caught.value.__cause__ is backend_error


def test_rewrite_failure_keeps_candidate_and_cause(monkeypatch):
    failure = ValueError("rewrite failed")

    class Unwritable:
        @property
        def rdkit_smiles(self):
            raise failure

    monkeypatch.setattr(implementation, "from_smiles", lambda text: Unwritable())
    with pytest.raises(ch.RdkitSmilesStabilityError) as caught:
        ch.rdkit_smiles_checked(ch.from_smiles("CCO"))
    assert caught.value.code == "rewrite_failed"
    assert caught.value.first_smiles == "CCO"
    assert caught.value.roundtrip_smiles is None
    assert caught.value.__cause__ is failure


def test_checks_once_even_if_later_iterations_would_converge(monkeypatch):
    calls = []

    class Changed:
        rdkit_smiles = "OCC"

    def read(text):
        calls.append(text)
        return Changed()

    monkeypatch.setattr(implementation, "from_smiles", read)
    with pytest.raises(ch.RdkitSmilesStabilityError, match="non_idempotent"):
        ch.rdkit_smiles_checked(ch.from_smiles("CCO"))
    assert calls == ["CCO"]


def test_unexpected_internal_error_is_not_reclassified(monkeypatch):
    def fail(text):
        raise AssertionError("internal failure")

    monkeypatch.setattr(implementation, "from_smiles", fail)
    with pytest.raises(AssertionError, match="internal failure"):
        ch.rdkit_smiles_checked(ch.from_smiles("CCO"))
