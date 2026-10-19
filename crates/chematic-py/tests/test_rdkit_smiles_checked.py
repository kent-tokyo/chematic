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
        "F[Pt@SP1](Cl)(Br)I",
    ],
)
def test_accepts_stable_outputs_without_modifying_input(source):
    mol = ch.from_smiles(source)
    before = mol.to_moljson()
    expected = mol.rdkit_smiles
    assert ch.rdkit_smiles_checked(mol) == expected
    assert ch.rdkit_smiles_checked(ch.from_smiles(expected)) == expected
    assert mol.to_moljson() == before


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
