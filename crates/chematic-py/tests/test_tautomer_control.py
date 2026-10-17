"""Real native enumeration checkpoints, cancellation and exception propagation."""

from dataclasses import FrozenInstanceError

import pytest
import chematic as ch


@pytest.mark.parametrize(
    "smiles",
    [
        "",
        "CCO",
        "CC=O",
        "CC(=O)C",
        "CC(=O)CC=O",
        "OC=C",
        "O=c1cccc[nH]1",
        "Oc1ccccn1",
        "NC(=N)N",
        "NC(=N)NC(=N)N",
        "CC(=O)/C=C/c1ccc(O)cc1",
        "C[C@@H](C=O)C(=O)O",
        "O=P(O)(O)O",
        "[NH4+].[Cl-]",
    ],
)
def test_continue_preserves_existing_enumeration(smiles):
    mol = ch.from_smiles(smiles) if smiles else ch.RWMol().GetMol()
    result = ch.enumerate_tautomers_checked(mol)
    assert sorted(result.smiles) == mol.rdkit_tautomers()
    assert result.status == mol.rdkit_tautomer_status()
    assert result.keys == tuple(sorted(set(result.keys)))
    assert len(result.keys) == len(result.smiles)


def test_cancel_first_checkpoint_returns_only_input():
    mol = ch.from_smiles("CC=O")
    trace = []

    def cancel(progress):
        trace.append(progress)
        return False

    result = ch.enumerate_tautomers_checked(mol, callback=cancel)
    assert result.status == "Canceled"
    assert result.smiles == result.keys == ("CC=O",)
    assert trace == [ch.TautomerProgress(1, ("CC=O",))]
    assert mol.smiles == "CC=O"


@pytest.mark.parametrize("stop", [2, 3, 5])
def test_cancel_later_and_restart(stop):
    mol = ch.from_smiles("CC(=O)CC=O")
    trace = []

    def cancel(progress):
        trace.append(progress)
        return len(trace) < stop

    result = ch.enumerate_tautomers_checked(mol, callback=cancel)
    assert len(trace) == stop
    assert result.status == "Canceled"
    assert len(result.smiles) > 1
    assert ch.enumerate_tautomers_checked(mol).status == "Completed"
    assert [p.transforms for p in trace] == sorted(p.transforms for p in trace)
    assert all(p.keys == tuple(sorted(set(p.keys))) for p in trace)


def test_snapshots_and_results_are_immutable():
    snapshots = []

    def proceed(progress):
        snapshots.append(progress)
        return True

    result = ch.enumerate_tautomers_checked(ch.from_smiles("CC=O"), callback=proceed)
    assert snapshots[0].keys == ("CC=O",)
    assert len(result.keys) == 2
    with pytest.raises(FrozenInstanceError):
        snapshots[0].transforms = 20
    with pytest.raises(FrozenInstanceError):
        result.status = "Canceled"
    with pytest.raises(TypeError):
        snapshots[0].keys[0] = "CCO"


@pytest.mark.parametrize("smiles", ["", "CCO", "[Na+].[Cl-]"])
def test_no_match_does_not_invoke_callback(smiles):
    def fail(_):
        pytest.fail("no transform can match")

    assert (
        ch.enumerate_tautomers_checked(
            ch.from_smiles(smiles) if smiles else ch.RWMol().GetMol(), callback=fail
        ).status
        == "Completed"
    )


@pytest.mark.parametrize(
    "kwargs,status",
    [
        ({"max_tautomers": 1}, "MaxTautomersReached"),
        ({"max_transforms": 1}, "MaxTransformsReached"),
        ({"max_tautomers": 1, "max_transforms": 1}, "MaxTransformsReached"),
    ],
)
def test_limits_precede_callback(kwargs, status):
    def fail(_):
        pytest.fail("limits must take precedence")

    result = ch.enumerate_tautomers_checked(
        ch.from_smiles("CC=O"), callback=fail, **kwargs
    )
    assert result.status == status
    assert result.smiles == ("CC=O",)


@pytest.mark.parametrize("value", [True, False, 1.0, "1", None])
@pytest.mark.parametrize("name", ["max_tautomers", "max_transforms"])
def test_limit_types(value, name):
    with pytest.raises(TypeError):
        ch.enumerate_tautomers_checked(ch.from_smiles("CC=O"), **{name: value})


@pytest.mark.parametrize("value", [-1, 0, 10001, 10**100])
@pytest.mark.parametrize("name", ["max_tautomers", "max_transforms"])
def test_limit_ranges(value, name):
    with pytest.raises(ValueError):
        ch.enumerate_tautomers_checked(ch.from_smiles("CC=O"), **{name: value})


@pytest.mark.parametrize("mol", [None, "CC=O", object()])
def test_wrong_molecule(mol):
    with pytest.raises(TypeError):
        ch.enumerate_tautomers_checked(mol)


@pytest.mark.parametrize("callback", [False, 42, [], object()])
def test_wrong_callback(callback):
    with pytest.raises(TypeError):
        ch.enumerate_tautomers_checked(ch.from_smiles("CC=O"), callback=callback)


@pytest.mark.parametrize("answer", [0, 1, None, [], "yes", object()])
def test_wrong_callback_result(answer):
    with pytest.raises(TypeError, match="callback must return bool"):
        ch.enumerate_tautomers_checked(
            ch.from_smiles("CC=O"), callback=lambda _: answer
        )


@pytest.mark.parametrize(
    "error",
    [
        ValueError("same object"),
        RuntimeError("stop"),
        KeyboardInterrupt(),
        SystemExit(17),
    ],
)
def test_callback_exception_propagates_unchanged(error):
    calls = []

    def fail(progress):
        calls.append(progress)
        raise error

    with pytest.raises(type(error)) as caught:
        ch.enumerate_tautomers_checked(ch.from_smiles("CC=O"), callback=fail)
    assert caught.value is error
    assert len(calls) == 1


def test_callable_instance_and_reentrancy():
    class Cancel:
        def __call__(self, progress):
            assert (
                ch.enumerate_tautomers_checked(ch.from_smiles("CCO")).status
                == "Completed"
            )
            return False

    assert (
        ch.enumerate_tautomers_checked(ch.from_smiles("CC=O"), callback=Cancel()).status
        == "Canceled"
    )


def test_atom_cap_before_native_conversion():
    mol = ch.from_smiles(".".join(["C"] * 1025))
    with pytest.raises(ValueError, match="1024"):
        ch.enumerate_tautomers_checked(mol)
