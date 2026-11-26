import math

import pytest

import chematic


def test_timed_conformer_generation_reports_mmff_unavailable():
    result = chematic.generate_conformers_timed("CB(C)C", 2, timeout_seconds=15)
    assert result["status"] == "completed"
    assert result["embedded"] > 0
    assert result["force_field_error"]
    assert result["conformers"]
    assert all(row[0] == -1 and row[1] is None for row in result["conformers"])


def test_timed_conformer_generation_terminates_and_reports_timeout():
    result = chematic.generate_conformers_timed("C" * 80, 20, timeout_seconds=0.1)
    assert result["status"] == "timeout"
    assert result["elapsed_seconds"] < 10


def test_timed_conformer_generation_reports_input_failure():
    result = chematic.generate_conformers_timed("not smiles", 1, timeout_seconds=15)
    assert result["status"] == "error"
    assert result["error_type"]


def test_timed_conformer_generation_rejects_a_different_worker_package(monkeypatch):
    monkeypatch.setattr(chematic, "__file__", "/a/different/chematic/__init__.py")
    result = chematic.generate_conformers_timed("CCO", 1, timeout_seconds=15)
    assert result["status"] == "error"
    assert "different chematic package" in result["error"]


def test_timed_conformer_generation_reports_unevaluable_stereo():
    result = chematic.generate_conformers_timed(
        "O=C(O[I@SP3]1OB(O)c2c(F)cccc21)C(F)(F)F", 2, timeout_seconds=15
    )
    assert result["status"] == "completed"
    assert result["embedded"] > 0
    assert result["conformers"] == []
    assert result["stereo_unevaluable"] == result["embedded"]


@pytest.mark.parametrize("timeout", [0, -1, math.inf, math.nan])
def test_timed_conformer_generation_rejects_invalid_deadlines(timeout):
    with pytest.raises(ValueError):
        chematic.generate_conformers_timed("CCO", 1, timeout_seconds=timeout)


def test_partial_mode_returns_individually_checked_attempts():
    result = chematic.generate_conformers_timed("CCO", 3, timeout_seconds=15, return_partial=True)
    assert result["status"] == "completed"
    assert result["attempts_completed"] == 3
    mol = chematic.from_smiles("CCO").add_hydrogens()
    expected = []
    for seed in (42, 84, 126):
        expected.extend(mol.etkdg_embed_multiple_confs_mixed_checked_diagnostics(1, random_seed=seed)[0])
    # JSON converts the nested Rust tuples into lists.
    import json
    assert result["conformers"] == json.loads(json.dumps(expected))


def test_partial_timeout_keeps_completed_checkpoint(monkeypatch):
    import json
    import subprocess
    from chematic import _timed_conformers as timed

    class FakeProcess:
        returncode = None
        killed = False

        def communicate(self, request=None, timeout=None):
            if request is None:
                return "", ""
            path = json.loads(request)["checkpoint"]
            with open(path, "w", encoding="utf-8") as stream:
                json.dump({"conformers": [[1, 0.0, [[0, 0, 0]], None]],
                           "attempts_completed": 1, "embedded": 1}, stream)
            raise subprocess.TimeoutExpired("worker", timeout)

        def kill(self):
            self.killed = True

        def poll(self):
            return -9 if self.killed else None

    worker = FakeProcess()
    monkeypatch.setattr(timed.subprocess, "Popen", lambda *args, **kwargs: worker)
    result = timed.generate_conformers_timed("C", 10, timeout_seconds=1, return_partial=True)
    assert worker.killed
    assert result["status"] == "timeout"
    assert result["attempts_completed"] == 1
    assert len(result["conformers"]) == 1


def test_partial_timeout_before_first_attempt_is_empty():
    result = chematic.generate_conformers_timed("C" * 80, 20, timeout_seconds=0.001, return_partial=True)
    assert result["status"] == "timeout"
    assert result["conformers"] == []
    assert result["attempts_completed"] == 0


def test_partial_diagnostics_do_not_report_incomplete_rejection_totals(tmp_path):
    from chematic._timed_conformers import _run_partial

    class MixedDiagnostics:
        seeds = []

        def etkdg_embed_multiple_confs_mixed_checked_diagnostics(self, count, **kwargs):
            self.seeds.append(kwargs["random_seed"])
            if len(self.seeds) == 1:
                return [], 1, "missing MMFF", 1, 0, 0
            return [], 1, None, None, None, None

    molecule = MixedDiagnostics()
    request = {"num_confs": 2, "random_seed": 0, "max_iterations": 0,
               "max_mmff_iters": 200, "repair_mmff_iters": 0,
               "checkpoint": str(tmp_path / "completed.json")}
    result = _run_partial(molecule, request)
    assert molecule.seeds == [0, 2]
    assert result["attempts_completed"] == 2
    assert result["stereo_violated"] is None
    assert result["stereo_unevaluable"] is None
    assert result["invalid_coordinates"] is None
    import json
    assert json.loads((tmp_path / "completed.json").read_text()) == result


def test_partial_checkpoint_survives_real_worker_termination(tmp_path):
    import subprocess
    import sys
    import time
    from chematic._timed_conformers import _read_checkpoint

    checkpoint = tmp_path / "completed.json"
    script = '''
import sys, time
from chematic._timed_conformers import _run_partial
class SlowMol:
    def etkdg_embed_multiple_confs_mixed_checked_diagnostics(self, count, **kwargs):
        time.sleep(0.05)
        return [(1, 0.0, [[0.0, 0.0, 0.0]], None)], 1, None, None, None, None
_run_partial(SlowMol(), {"num_confs": 1000, "random_seed": 42, "max_iterations": 0,
    "max_mmff_iters": 200, "repair_mmff_iters": 0, "checkpoint": sys.argv[1]})
'''
    worker = subprocess.Popen([sys.executable, "-c", script, str(checkpoint)])
    try:
        deadline = time.monotonic() + 10
        while not checkpoint.exists() and time.monotonic() < deadline and worker.poll() is None:
            time.sleep(0.01)
        assert checkpoint.exists(), "worker did not publish its first checkpoint"
    finally:
        worker.kill()
        worker.wait()
    result = _read_checkpoint(str(checkpoint))
    assert result["attempts_completed"] >= 1
    assert len(result["conformers"]) == result["attempts_completed"]
