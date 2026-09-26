"""The #635 evidence checker must reject dropped or misclassified SMARTS cells."""

from scripts.check_rdkit_rebaseline_evidence import validate_smarts_cell_correspondence


def sample():
    rows = [
        {
            "input_index": 0,
            "status": "completed",
            "smarts": {
                "query_count": 31,
                "difference_count": 1,
                "differences": [
                    {
                        "query": "[R2]",
                        "chematic": [[1]],
                        "rdkit": [[1], [2]],
                        "chematic_error": None,
                    }
                ],
            },
        }
    ]
    classification = {
        "rows": [
            {
                "input_index": 0,
                "ring_facts": {"has_metal": False},
                "cells": [{"query": "[R2]", "family": "ring_semantics:symmetrized_rings"}],
            }
        ]
    }
    return rows, classification


def test_complete_cell_correspondence_passes():
    rows, classification = sample()
    errors = []
    validate_smarts_cell_correspondence(rows, classification, errors)
    assert errors == []


def test_missing_classified_query_is_rejected():
    rows, classification = sample()
    classification["rows"][0]["cells"] = []
    errors = []
    validate_smarts_cell_correspondence(rows, classification, errors)
    assert any("classified queries differ" in error for error in errors)


def test_refusal_cannot_be_silently_counted_as_ring_difference():
    rows, classification = sample()
    rows[0]["smarts"]["differences"][0]["chematic_error"] = "budget_exceeded"
    errors = []
    validate_smarts_cell_correspondence(rows, classification, errors)
    assert any("unaccounted refusal" in error for error in errors)


def test_organometallic_family_requires_a_metal():
    rows, classification = sample()
    classification["rows"][0]["cells"][0]["family"] = "ring_semantics:organometallic"
    errors = []
    validate_smarts_cell_correspondence(rows, classification, errors)
    assert any("unexpected organometallic" in error for error in errors)
