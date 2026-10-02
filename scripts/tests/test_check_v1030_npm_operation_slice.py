"""Published npm operation evidence cannot hide a chemical value mismatch."""

from scripts.check_v1030_artifact_packet import formula_composition, npm_python_outcome


def test_formula_order_is_spelling_only_when_counts_are_equal() -> None:
    assert formula_composition("C23H20N4O5SCl2") == formula_composition("C23H20Cl2N4O5S")
    assert npm_python_outcome(
        "formula",
        {"input_index": 2, "value": "C23H20N4O5SCl2"},
        {"input_index": 2, "value": "C23H20Cl2N4O5S"},
    ) == "formula_spelling_only"
    assert npm_python_outcome(
        "formula",
        {"input_index": 2, "value": "C23H20N4O5SCl"},
        {"input_index": 2, "value": "C23H20Cl2N4O5S"},
    ) == "different"


def test_roundoff_is_only_allowed_on_named_numeric_operations() -> None:
    actual = {"input_index": 5, "value": 0.777586729984065}
    expected = {"input_index": 5, "value": 0.7775867299840651}
    assert npm_python_outcome("qed", actual, expected) == "numeric_roundoff"
    assert npm_python_outcome("mw", actual, expected) == "different"
    assert npm_python_outcome("qed", {**actual, "input_index": 6}, expected) == "different"
