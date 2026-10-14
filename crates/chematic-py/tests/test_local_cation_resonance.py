"""Local resonance API and explicit limits, inspired by RDKit #9458."""
import chematic
import pytest


@pytest.mark.parametrize("pair", [
    ("COC(C)=[N+](C)C", "C[O+]=C(C)N(C)C"),
    ("CSC(C)=[N+](C)C", "C[S+]=C(C)N(C)C"),
])
def test_both_starting_forms_reach_the_same_pair(pair):
    expected = {chematic.from_smiles(s).smiles for s in pair}
    for smiles in pair:
        report = chematic.from_smiles(smiles).local_cation_resonance()
        assert report["status"] == "complete_within_scope"
        assert report["scope"] == "local_nos_cation_shift_v1"
        assert report["transforms_attempted"] == 2
        assert len(report["forms"]) == 2
        assert {m.smiles for m in report["forms"]} == expected


def test_budget_statuses_are_explicit():
    mol = chematic.from_smiles("CSC(C)=[N+](C)C")
    assert mol.local_cation_resonance(max_forms=2)["status"] == "complete_within_scope"
    limited = mol.local_cation_resonance(max_forms=1)
    assert limited["status"] == "form_limit_reached"
    assert len(limited["forms"]) == 1
    limited = mol.local_cation_resonance(max_transforms=1)
    assert limited["status"] == "transform_limit_reached"
    assert limited["transforms_attempted"] == 1


@pytest.mark.parametrize("kwargs", [
    {"max_forms": 0}, {"max_forms": 1025}, {"max_transforms": 0},
    {"max_transforms": 100001}, {"max_atoms": 0}, {"max_atoms": 2},
])
def test_invalid_limits_raise_valueerror(kwargs):
    with pytest.raises(ValueError):
        chematic.from_smiles("CCO").local_cation_resonance(**kwargs)


@pytest.mark.parametrize("smiles", [
    "C/C(SC)=[N+](C)/C", "C[C@H](F)CSC(C)=[N+](C)C", "*C",
])
def test_stereo_and_queries_are_not_silently_projected(smiles):
    with pytest.raises(ValueError):
        chematic.from_smiles(smiles).local_cation_resonance()


def test_indexed_symmetric_placements_are_retained():
    report = chematic.from_smiles("NC(N)=[NH2+]").local_cation_resonance()
    assert report["status"] == "complete_within_scope"
    assert len(report["forms"]) == 3
    assert len({m.smiles for m in report["forms"]}) == 1
