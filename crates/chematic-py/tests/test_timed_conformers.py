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
                json.dump({"attempt": 1, "rows": [[1, 0.0, [[0, 0, 0]], None]], "embedded": 1,
                           "error": None, "violated": None, "unevaluable": None, "invalid": None}, stream)
                stream.write("\n")
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
    from chematic._timed_conformers import _read_checkpoint
    assert _read_checkpoint(str(tmp_path / "completed.json")) == result


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


@pytest.mark.parametrize("worker_output,exit_code,expected_status", [
    ('{"status":"error","error_type":"ValueError","error":"calculation failed"}', 0, "error"),
    ("", 1, "worker_error"),
    ("invalid json", 0, "worker_error"),
])
def test_partial_worker_failures_keep_completed_results(monkeypatch, worker_output, exit_code, expected_status):
    import json
    from chematic import _timed_conformers as timed

    class FakeProcess:
        returncode = exit_code

        def communicate(self, request=None, timeout=None):
            if request is not None:
                checkpoint = json.loads(request)["checkpoint"]
                with open(checkpoint, "w", encoding="utf-8") as stream:
                    json.dump({"attempt": 1, "rows": [[1, 0.0, [[0, 0, 0]], None]], "embedded": 1,
                               "error": None, "violated": None, "unevaluable": None, "invalid": None}, stream)
                    stream.write("\n")
            return worker_output, "worker crashed"

        def poll(self):
            return self.returncode

    monkeypatch.setattr(timed.subprocess, "Popen", lambda *args, **kwargs: FakeProcess())
    result = timed.generate_conformers_timed("C", 10, timeout_seconds=1, return_partial=True)
    assert result["status"] == expected_status
    assert result["error"]
    assert result["attempts_completed"] == 1
    assert len(result["conformers"]) == 1


def test_invalid_checkpoint_preserves_timeout_status(monkeypatch):
    import json
    import subprocess
    from chematic import _timed_conformers as timed

    class FakeProcess:
        def communicate(self, request=None, timeout=None):
            if request is None:
                return "", ""
            checkpoint = json.loads(request)["checkpoint"]
            with open(checkpoint, "w", encoding="utf-8") as stream:
                stream.write("{\n")
            raise subprocess.TimeoutExpired("worker", timeout)

        def kill(self):
            pass

        def poll(self):
            return -9

    monkeypatch.setattr(timed.subprocess, "Popen", lambda *args, **kwargs: FakeProcess())
    result = timed.generate_conformers_timed("C", 10, timeout_seconds=1, return_partial=True)
    assert result["status"] == "timeout"
    assert "JSONDecodeError" in result["checkpoint_error"]
    assert "conformers" not in result


def test_partial_real_worker_error_keeps_prior_completed_attempt():
    result = chematic.generate_conformers_timed(
        "CCO", 2, random_seed=2**31 - 1, timeout_seconds=15, return_partial=True
    )
    assert result["status"] == "error"
    assert result["error_type"] == "ValueError"
    assert "seed" in result["error"]
    assert result["attempts_completed"] == 1
    assert result["embedded"] == 1
    assert len(result["conformers"]) == 1


def test_partial_checkpoint_ignores_interrupted_last_write(tmp_path):
    import json
    from chematic._timed_conformers import _read_checkpoint

    path = tmp_path / "journal.jsonl"
    entry = {"attempt": 1, "rows": [[1, 0.0, [[0, 0, 0]], None]], "embedded": 1,
             "error": None, "violated": None, "unevaluable": None, "invalid": None}
    with open(path, "wb") as stream:
        stream.write((json.dumps(entry) + "\n").encode("utf-8"))
        stream.write(b'{"attempt":2,"rows":[[1')
    result = _read_checkpoint(str(path))
    assert result["attempts_completed"] == 1
    assert len(result["conformers"]) == 1


def test_partial_checkpoint_rejects_out_of_order_records(tmp_path):
    import json
    from chematic._timed_conformers import _read_checkpoint

    path = tmp_path / "journal.jsonl"
    entry = {"attempt": 2, "rows": [], "embedded": 0, "error": None,
             "violated": None, "unevaluable": None, "invalid": None}
    path.write_text(json.dumps(entry) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="invalid checkpoint entry"):
        _read_checkpoint(str(path))


def test_preview_requires_partial_mode():
    with pytest.raises(ValueError, match="requires return_partial"):
        chematic.generate_conformers_timed("CCO", 1, timeout_seconds=1, return_preview=True)


@pytest.mark.parametrize("smiles", ["CCO", "CB(C)C", "F[C@](Cl)(Br)I"])
def test_preview_completed_attempt_matches_existing_partial_mode(smiles):
    options = {"timeout_seconds": 15, "return_partial": True, "num_confs": 2}
    baseline = chematic.generate_conformers_timed(smiles, **options)
    preview = chematic.generate_conformers_timed(smiles, return_preview=True, **options)
    assert baseline["status"] == preview["status"] == "completed"
    assert "in_progress" not in preview
    for key in ("conformers", "attempts_completed", "embedded", "force_field_error",
                "stereo_violated", "stereo_unevaluable", "invalid_coordinates"):
        assert preview[key] == baseline[key]


def test_preview_checkpoint_survives_termination_during_optimization(tmp_path):
    import subprocess
    import sys
    import time
    from chematic._timed_conformers import _read_checkpoint

    checkpoint = tmp_path / "completed.jsonl"
    script = '''
import sys, time
from chematic._timed_conformers import _run_partial
class SlowMol:
    def _timed_embed_corrected(self, **kwargs): return [[0.0, 0.0, 0.0]]
    def _timed_stereo_valid(self, coords): return True
    def _timed_optimize_embedded(self, coords, **kwargs):
        time.sleep(30)
        return [(1, 0.0, coords, None)], 1, None, None, None, None
_run_partial(SlowMol(), {"num_confs": 1, "random_seed": 42, "max_iterations": 0,
    "max_mmff_iters": 200, "repair_mmff_iters": 0, "return_preview": True, "checkpoint": sys.argv[1]})
'''
    worker = subprocess.Popen([sys.executable, "-c", script, str(checkpoint)])
    try:
        deadline = time.monotonic() + 10
        while not checkpoint.exists() and time.monotonic() < deadline and worker.poll() is None:
            time.sleep(0.01)
        assert checkpoint.exists(), "worker did not publish its preview"
    finally:
        worker.kill()
        worker.wait()
    result = _read_checkpoint(str(checkpoint))
    assert result["attempts_completed"] == 0
    assert result["conformers"] == []
    assert result["in_progress"] == {"attempt": 1, "coordinates": [[0.0, 0.0, 0.0]]}


def test_preview_is_withheld_for_unevaluable_stereo():
    mol = chematic.from_smiles("O=C(O[I@SP3]1OB(O)c2c(F)cccc21)C(F)(F)F").add_hydrogens()
    coords = mol._timed_embed_corrected(random_seed=42)
    assert coords is not None
    assert not mol._timed_stereo_valid(coords)
    with pytest.raises(ValueError):
        mol._timed_stereo_valid([[0.0, 0.0, 0.0]])


def test_preview_never_publishes_failed_stereo_check(tmp_path):
    from chematic._timed_conformers import _read_checkpoint, _run_partial

    class InvalidStereoMol:
        def _timed_embed_corrected(self, **kwargs):
            return [[0.0, 0.0, 0.0]]

        def _timed_stereo_valid(self, coords):
            return False

        def _timed_optimize_embedded(self, coords, **kwargs):
            raise RuntimeError("interrupted before optimization completed")

    path = tmp_path / "journal.jsonl"
    request = {"num_confs": 1, "random_seed": 42, "max_iterations": 0,
               "max_mmff_iters": 200, "repair_mmff_iters": 0,
               "return_preview": True, "checkpoint": str(path)}
    with pytest.raises(RuntimeError, match="interrupted"):
        _run_partial(InvalidStereoMol(), request)
    result = _read_checkpoint(str(path))
    assert result["attempts_completed"] == 0
    assert "in_progress" not in result


def test_completed_attempt_supersedes_preview_in_journal(tmp_path):
    import json
    from chematic._timed_conformers import _read_checkpoint

    path = tmp_path / "journal.jsonl"
    preview = {"preview_attempt": 1, "coordinates": [[0.0, 0.0, 0.0]]}
    completed = {"attempt": 1, "rows": [[0, 1.0, [[1.0, 0.0, 0.0]], None]],
                 "embedded": 1, "error": None, "violated": None, "unevaluable": None, "invalid": None}
    with open(path, "w", encoding="utf-8") as stream:
        stream.write(json.dumps(preview) + "\n" + json.dumps(completed) + "\n")
    result = _read_checkpoint(str(path))
    assert result["attempts_completed"] == 1
    assert "in_progress" not in result
    assert result["conformers"] == completed["rows"]


def test_preview_stereo_check_rejects_a_mirrored_tetrahedral_center():
    mol = chematic.from_smiles("F[C@](Cl)(Br)I").add_hydrogens()
    coords = mol._timed_embed_corrected(random_seed=42)
    assert coords is not None
    assert mol._timed_stereo_valid(coords)
    mirrored = [[-x, y, z] for x, y, z in coords]
    assert not mol._timed_stereo_valid(mirrored)
