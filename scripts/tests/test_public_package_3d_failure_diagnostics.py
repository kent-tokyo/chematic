"""A failed A6 row retains typed diagnostics without becoming a scored conformer."""

from __future__ import annotations

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "public_package_3d_chematic.py"
SPEC = importlib.util.spec_from_file_location("public_package_3d_chematic", SCRIPT)
assert SPEC and SPEC.loader
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class FailingMol:
    def embed_pipeline_v2(self, config):
        error = ValueError("final stereo violation")
        error.diagnostics = {
            "stage": "final_stereo_verify",
            "cause": {"kind": "final_stereo_violation"},
            "last_known_coords": [[0.0, 1.0, 2.0]],
            "coords_are_diagnostic_only": True,
        }
        raise error


def test_typed_failure_keeps_only_diagnostic_coordinates():
    result = module.run_pipeline(FailingMol(), object())
    assert result["status"] == "typed_failure"
    assert result["failure_diagnostics"]["cause"] == {
        "kind": "final_stereo_violation"
    }
    assert result["failure_diagnostics"]["coords_are_diagnostic_only"] is True
    assert "coords" not in result
    assert "force_field" not in result


def test_untyped_failure_does_not_invent_diagnostics():
    class UnexpectedMol:
        def embed_pipeline_v2(self, config):
            raise RuntimeError("unexpected")

    result = module.run_pipeline(UnexpectedMol(), object())
    assert result["status"] == "typed_failure"
    assert "failure_diagnostics" not in result


def test_non_json_diagnostic_does_not_erase_failed_row():
    class MalformedMol:
        def embed_pipeline_v2(self, config):
            error = ValueError("failed")
            error.diagnostics = {"unsupported": object()}
            raise error

    result = module.run_pipeline(MalformedMol(), object())
    assert result["status"] == "typed_failure"
    assert "failure_diagnostics" not in result
    assert result["diagnostic_serialization_error"].startswith("TypeError:")
