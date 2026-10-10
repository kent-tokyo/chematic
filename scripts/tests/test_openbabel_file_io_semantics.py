"""Unit tests for the Open Babel semantic and paired-speed gates."""

from __future__ import annotations

from scripts.bench_openbabel_file_io_paired import ci95_lower
from scripts.check_openbabel_file_io_semantics import compare, coordinate_comparison


def record(*, distances: list[float], metadata: dict | None = None) -> dict:
    return {
        "graph_key": "CCO",
        "atom_count": 3,
        "bond_count": 2,
        "total_formal_charge": 0,
        "isotopes": [],
        "pairwise_distances": distances,
        "metadata": metadata or {},
    }


def test_cdxml_uniform_scale_is_coordinate_equivalent():
    delta, scale = coordinate_comparison(
        [1.5, 3.0, 1.5], [3.0, 6.0, 3.0], allow_uniform_scale=True
    )
    assert delta == 0.0
    assert scale == 2.0


def test_non_cdxml_scale_change_is_not_coordinate_equivalent():
    result = compare(
        "cml",
        {"records": [record(distances=[1.5, 3.0, 1.5])]},
        {"records": [record(distances=[3.0, 6.0, 3.0])]},
    )
    assert result["status"] == "semantic_mismatch"
    assert any("coordinate distance" in error for error in result["errors"])


def test_mol2_partial_charge_change_is_a_semantic_mismatch():
    baseline = {"atoms": [{"partial_charge": -0.394}], "bonds": []}
    changed = {"atoms": [{"partial_charge": 0.0}], "bonds": []}
    result = compare(
        "mol2",
        {"records": [record(distances=[], metadata=baseline)]},
        {"records": [record(distances=[], metadata=changed)]},
    )
    assert result["status"] == "semantic_mismatch"
    assert result["errors"] == ["record 0: format metadata differs"]


def test_paired_speed_lower_bound_is_below_mean_for_variable_samples():
    values = [2.0 + index / 100 for index in range(21)]
    assert 1.0 < ci95_lower(values) < sum(values) / len(values)
