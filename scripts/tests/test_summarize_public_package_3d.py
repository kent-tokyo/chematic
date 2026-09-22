import pytest

from scripts.summarize_public_package_3d import paired_speedup, percentile


def row(name, elapsed):
    return {
        "tier": "A",
        "name": name,
        "status": "success",
        "elapsed_ms": elapsed,
    }


def test_percentile_interpolates():
    assert percentile([1.0, 3.0], 0.5) == pytest.approx(2.0)


def test_paired_speedup_uses_common_successes():
    result = paired_speedup(
        [row("a", 1.0), row("b", 2.0)],
        [row("a", 2.0), row("b", 4.0)],
    )

    assert result["common_successes"] == 2
    assert result["geometric_mean_rdkit_over_chematic"] == pytest.approx(2.0)
    assert result["chematic_faster_gate"] is True
