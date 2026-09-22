import pytest

from scripts.check_parse_morgan_rdkit_speed_gate import paired_summary


def measurement(chematic_values, rdkit_values):
    def arm(values):
        return [
            {
                "operations": {
                    "parse_fp": {"mean_ms": value},
                    "prepared_fp": {"mean_ms": value / 2},
                }
            }
            for value in values
        ]

    return {
        "schema_version": 2,
        "raw_runs": {
            "chematic": arm(chematic_values),
            "rdkit": arm(rdkit_values),
        },
    }


def test_paired_summary_selects_requested_operation():
    result = paired_summary(measurement([1.0, 1.0], [2.0, 2.0]), "prepared_fp")

    assert result["operation"] == "prepared_fp"
    assert result["geometric_mean_speedup"] == pytest.approx(2.0)


def test_paired_summary_rejects_missing_operation():
    with pytest.raises(ValueError, match="prepared_fp"):
        paired_summary(
            {
                "schema_version": 2,
                "raw_runs": {
                    "chematic": [
                        {"operations": {"parse_fp": {"mean_ms": 1.0}}},
                        {"operations": {"parse_fp": {"mean_ms": 1.0}}},
                    ],
                    "rdkit": [
                        {"operations": {"parse_fp": {"mean_ms": 2.0}}},
                        {"operations": {"parse_fp": {"mean_ms": 2.0}}},
                    ],
                },
            },
            "prepared_fp",
        )
