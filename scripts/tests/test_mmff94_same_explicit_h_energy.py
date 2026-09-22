from scripts.mmff94_same_explicit_h_energy import (
    gradient_diagnostic,
    percentile_nearest_rank,
    summarize,
    terminal,
)


class QuadraticEnergy:
    def mmff94_bounded_analytic_gradient(self, coords):
        return [[2.0 * value for value in point] for point in coords]

    def mmff94_energy_breakdown(self, coords):
        return {"total": sum(value * value for point in coords for value in point)}


def test_nearest_rank_percentile_is_deterministic():
    assert percentile_nearest_rank([], 9, 10) is None
    assert percentile_nearest_rank([3.0, 1.0, 2.0, 4.0], 9, 10) == 4.0


def test_summary_preserves_terminal_row_accounting():
    rows = [
        terminal(
            {"primary_category": "small"},
            0,
            "ok",
            abs_delta_kcal_mol=0.5,
        ),
        terminal({}, 1, "embed_failure"),
    ]
    result = summarize(rows)
    assert result["row_accounting"] == {
        "input_count": 2,
        "terminal_count": 2,
        "status_counts": {"ok": 1, "embed_failure": 1},
    }
    assert result["comparable_rows"] == 1
    assert result["within_1_kcal_mol"] == 1
    assert result["within_5_kcal_mol"] == 1
    assert result["gradient_diagnostic"] == {
        "rows": 0,
        "input_indices": [],
        "max_abs_error_kcal_mol_angstrom": None,
        "max_scaled_error": None,
    }


def test_gradient_diagnostic_matches_quadratic_energy():
    result = gradient_diagnostic(QuadraticEnergy(), [[1.0, -2.0, 0.5]], 1e-5)
    assert result["components"] == 3
    assert result["max_abs_error_kcal_mol_angstrom"] < 1e-9
    assert result["max_scaled_error"] < 1e-9
