from scripts.mmff94_same_explicit_h_energy import (
    percentile_nearest_rank,
    summarize,
    terminal,
)


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
