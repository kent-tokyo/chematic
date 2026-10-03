"""A6 evidence must distinguish usable geometry from convergence/readiness."""

import importlib.util
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "a6_quality_matrix.py"
sys.path.insert(0, str(SCRIPT.parent))
SPEC = importlib.util.spec_from_file_location("a6_quality_matrix", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_historical_quality_does_not_clear_current_a6():
    matrix = MODULE.build_matrix("1.0.31")
    assert matrix["ready_for_target_version"] is False
    assert matrix["speed_excluded"] is True
    assert matrix["dimensions"]["stereo"]["independent_scorer_clean"] == 265
    assert (
        matrix["dimensions"]["clash_and_geometry"]["independent_scorer_sound_and_clash_free"]
        == 265
    )
    assert matrix["dimensions"]["convergence"]["converged"] == 100
    assert matrix["dimensions"]["convergence"]["not_converged"] == 165
    assert (
        matrix["dimensions"]["typing_and_parameters"]["separate_census_heavy_type_mismatches"]
        == 2908
    )
    assert matrix["dimensions"]["same_coordinate_energy"]["over_5_kcal_mol_indices"] == [166, 231]
    assert matrix["dimensions"]["independent_conformer_quality"]["status"] == "not_measured"


def test_absent_convergence_result_is_not_treated_as_success():
    rows = [
        {
            "status": "success",
            "row_index": index,
            "force_field": {
                "coverage": {
                    f"{term}_missing": []
                    for term in ("bonds", "angles", "torsions", "oop", "stretch_bend")
                },
            },
            "final_validation": {"sound": True, "stereo_ok": True, "gross_clash_count": 0},
        }
        for index in range(265)
    ]
    with pytest.raises(ValueError, match="convergence"):
        MODULE.published_dimensions(rows)
