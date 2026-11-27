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


def test_checked_mmff_reports_discarded_stereo_violation():
    mol = chematic.from_smiles("C[C@H](O)F").add_hydrogens()
    coords = mol._timed_prepare_corrected().embed(42)
    assert coords is not None and mol._timed_stereo_valid(coords)
    mirror = [[-x, y, z] for x, y, z in coords]
    assert not mol._timed_stereo_valid(mirror)
    rows, embedded, error, violated, unevaluable, invalid = mol._timed_optimize_embedded(
        mirror, max_mmff_iters=0
    )
    assert rows == [] and embedded == 1 and error is None
    assert (violated, unevaluable, invalid) == (1, 0, 0)


def test_timed_mmff_reports_zero_rejection_counts_when_attempts_succeed():
    for preview in (False, True):
        result = chematic.generate_conformers_timed(
            "CCO", 2, timeout_seconds=15, return_partial=True, return_preview=preview
        )
        assert result["status"] == "completed"
        assert result["embedded"] == 2 and result["force_field_error"] is None
        assert (result["stereo_violated"], result["stereo_unevaluable"],
                result["invalid_coordinates"]) == (0, 0, 0)


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
    assert [item["attempt"] for item in result["attempt_diagnostics"]] == [1, 2, 3]
    assert all(item["embedded"] == item["returned"] == 1 for item in result["attempt_diagnostics"])
    assert sum(sum(item["conformer_status_counts"].values())
               for item in result["attempt_diagnostics"]) == len(expected)


def test_attempt_diagnostics_distinguish_statuses_and_reject_bad_rows():
    from chematic._timed_conformers import _apply_checkpoint_entry, _empty_checkpoint

    entry = {"attempt": 1, "rows": [[0, 1.0, [], None], [1, 2.0, [], None],
             [-1, 3.0, [], None], [-1, None, [], "unsupported"], [42, 4.0, [], None]],
             "embedded": 9, "error": None, "violated": 2, "unevaluable": 1, "invalid": 0}
    result = _empty_checkpoint()
    known = {key: True for key in ("stereo_violated", "stereo_unevaluable", "invalid_coordinates",
                                   "unclassified_rejections")}
    _apply_checkpoint_entry(result, known, entry)
    assert result["attempt_diagnostics"] == [{
        "attempt": 1, "embedded": 9, "returned": 5, "force_field_error": None,
        "stereo_violated": 2, "stereo_unevaluable": 1, "invalid_coordinates": 0,
        "unclassified_rejections": 1,
        "conformer_status_counts": {"converged": 1, "unfinished": 1, "minimizer_failed": 1,
                                    "force_field_unavailable": 1, "other": 1},
    }]
    invalid_entry = {**entry, "attempt": 2, "rows": [[1]]}
    with pytest.raises(ValueError, match="invalid checkpoint conformer"):
        _apply_checkpoint_entry(result, known, invalid_entry)
    assert result["attempts_completed"] == len(result["attempt_diagnostics"]) == 1
    assert result["conformers"] == entry["rows"]
    assert result["unclassified_rejections"] == 1


@pytest.mark.parametrize("embedded,rows,counts,expected", [
    (3, [], (1, 1, 0), 1),
    (1, [], (None, 0, 0), None),
    (1, [[1, 1.0, [], None]], (None, None, None), 0),
    (0, [], (None, None, None), 0),
])
def test_unclassified_rejections_are_unknown_only_when_needed(embedded, rows, counts, expected):
    from chematic._timed_conformers import _apply_checkpoint_entry, _empty_checkpoint

    entry = {"attempt": 1, "rows": rows, "embedded": embedded, "error": None,
             "violated": counts[0], "unevaluable": counts[1], "invalid": counts[2]}
    result = _empty_checkpoint()
    known = {key: True for key in ("stereo_violated", "stereo_unevaluable", "invalid_coordinates",
                                   "unclassified_rejections")}
    _apply_checkpoint_entry(result, known, entry)
    assert result["attempt_diagnostics"][0]["unclassified_rejections"] == expected
    assert result["unclassified_rejections"] == expected


@pytest.mark.parametrize("embedded,rows,counts", [
    (0, [[1, 1.0, [], None]], (0, 0, 0)),
    (1, [], (1, 1, None)),
])
def test_inconsistent_candidate_counts_are_rejected_atomically(embedded, rows, counts):
    from chematic._timed_conformers import _apply_checkpoint_entry, _empty_checkpoint

    entry = {"attempt": 1, "rows": rows, "embedded": embedded, "error": None,
             "violated": counts[0], "unevaluable": counts[1], "invalid": counts[2]}
    result = _empty_checkpoint()
    known = {key: True for key in ("stereo_violated", "stereo_unevaluable", "invalid_coordinates",
                                   "unclassified_rejections")}
    with pytest.raises(ValueError, match="inconsistent checkpoint candidate counts"):
        _apply_checkpoint_entry(result, known, entry)
    assert result == _empty_checkpoint()
    assert all(known.values())


def test_unclassified_aggregate_remains_unknown_after_an_unknown_attempt(tmp_path):
    import json
    from chematic._timed_conformers import _read_checkpoint

    entries = [
        {"attempt": 1, "rows": [], "embedded": 1, "error": None,
         "violated": None, "unevaluable": 0, "invalid": 0},
        {"attempt": 2, "rows": [], "embedded": 2, "error": None,
         "violated": 1, "unevaluable": 0, "invalid": 0},
    ]
    path = tmp_path / "attempts.jsonl"
    path.write_text("".join(json.dumps(entry) + "\n" for entry in entries), encoding="utf-8")
    result = _read_checkpoint(str(path))
    assert [entry["unclassified_rejections"] for entry in result["attempt_diagnostics"]] == [None, 1]
    assert result["unclassified_rejections"] is None


def test_inconsistent_checkpoint_counts_preserve_previous_attempt(tmp_path):
    import json
    from chematic._timed_conformers import _attach_checkpoint, _read_checkpoint

    first = {"attempt": 1, "rows": [], "embedded": 2, "error": None,
             "violated": 1, "unevaluable": 0, "invalid": 0}
    bad = {**first, "attempt": 2, "embedded": 0}
    path = tmp_path / "attempts.jsonl"
    path.write_text(json.dumps(first) + "\n" + json.dumps(bad) + "\n", encoding="utf-8")
    response = {"status": "timeout"}
    _attach_checkpoint(response, str(path))
    assert response["status"] == "timeout"
    assert response["attempts_completed"] == 1
    assert response["unclassified_rejections"] == 1
    assert "line 2: ValueError: inconsistent checkpoint candidate counts" in response["checkpoint_error"]
    with pytest.raises(ValueError, match="inconsistent checkpoint candidate counts"):
        _read_checkpoint(str(path))


def test_attempt_diagnostics_record_an_empty_embedding(tmp_path):
    from chematic._timed_conformers import _read_checkpoint, _run_partial

    class NoEmbedding:
        def etkdg_embed_multiple_confs_mixed_checked_diagnostics(self, count, **kwargs):
            return [], 0, None, None, None, None

    checkpoint = str(tmp_path / "attempts.jsonl")
    request = {"num_confs": 1, "random_seed": 42, "max_iterations": 0,
               "max_mmff_iters": 200, "repair_mmff_iters": 0, "checkpoint": checkpoint}
    result = _run_partial(NoEmbedding(), request)
    attempt = result["attempt_diagnostics"][0]
    assert attempt["embedded"] == attempt["returned"] == 0
    assert attempt["stereo_violated"] is None
    assert sum(attempt["conformer_status_counts"].values()) == 0
    assert _read_checkpoint(checkpoint) == result


def test_partial_attempt_diagnostics_report_unsupported_mmff():
    result = chematic.generate_conformers_timed("CB(C)C", 1, timeout_seconds=15, return_partial=True)
    assert result["status"] == "completed"
    attempt = result["attempt_diagnostics"][0]
    assert attempt["embedded"] == attempt["returned"] == 1
    assert attempt["force_field_error"]
    assert attempt["conformer_status_counts"]["force_field_unavailable"] == 1
    assert attempt["conformer_status_counts"]["minimizer_failed"] == 0


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
    assert result["attempt_diagnostics"][0]["conformer_status_counts"]["unfinished"] == 1


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


@pytest.mark.parametrize("failure", ["timeout", "worker_error"])
@pytest.mark.parametrize("bad_record", [b"{\n", b'{"attempt":4}\n', b'\xff\n'])
def test_corrupt_checkpoint_keeps_verified_prefix(monkeypatch, failure, bad_record):
    import json
    import subprocess
    from chematic import _timed_conformers as timed

    completed = {"attempt": 1, "rows": [[1, 0.0, [[0, 0, 0]], None]], "embedded": 1,
                 "error": None, "violated": None, "unevaluable": None, "invalid": None}
    preview = {"preview_attempt": 2, "coordinates": [[1, 2, 3]]}

    class FakeProcess:
        returncode = 1

        def communicate(self, request=None, timeout=None):
            if request is None:
                return "", ""
            checkpoint = json.loads(request)["checkpoint"]
            with open(checkpoint, "wb") as stream:
                for record in (completed, preview):
                    stream.write((json.dumps(record) + "\n").encode("utf-8"))
                stream.write(bad_record)
                stream.write((json.dumps({**completed, "attempt": 2}) + "\n").encode("utf-8"))
            if failure == "timeout":
                raise subprocess.TimeoutExpired("worker", timeout)
            return "", "worker crashed"

        def kill(self):
            pass

        def poll(self):
            return self.returncode

    monkeypatch.setattr(timed.subprocess, "Popen", lambda *args, **kwargs: FakeProcess())
    result = timed.generate_conformers_timed("C", 10, timeout_seconds=1,
                                             return_partial=True, return_preview=True)
    assert result["status"] == failure
    assert result["attempts_completed"] == 1
    assert result["conformers"] == completed["rows"]
    assert result["in_progress"] == {"attempt": 2, "coordinates": [[1, 2, 3]]}
    assert result["checkpoint_error"].startswith("line 3:")
    if failure == "worker_error":
        assert result["error"] == "worker crashed"


def test_checkpoint_read_error_keeps_verified_prefix(monkeypatch):
    import json
    from chematic import _timed_conformers as timed

    entry = {"attempt": 1, "rows": [], "embedded": 0, "error": None,
             "violated": None, "unevaluable": None, "invalid": None}

    class FailingStream:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def __iter__(self):
            yield (json.dumps(entry) + "\n").encode("utf-8")
            raise OSError("disk read failed")

    monkeypatch.setattr(timed, "open", lambda *args: FailingStream(), raising=False)
    response = {"status": "worker_error", "error": "worker crashed"}
    timed._attach_checkpoint(response, "ignored")
    assert response["status"] == "worker_error"
    assert response["error"] == "worker crashed"
    assert response["attempts_completed"] == 1
    assert response["checkpoint_error"] == "OSError: disk read failed"
    with pytest.raises(OSError, match="disk read failed"):
        timed._read_checkpoint("ignored")


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
    def _timed_prepare_corrected(self, **kwargs): return self
    def embed(self, seed): return [[0.0, 0.0, 0.0]]
    def _timed_stereo_valid(self, coords): return True
    def _timed_preview_geometry(self, coords): return (0, None, None, 0)
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
    assert result["in_progress"] == {"attempt": 1, "coordinates": [[0.0, 0.0, 0.0]],
        "geometry": {"bond_count": 0, "worst_bond_length": None,
                     "bond_violation_rate_50pct": None, "gross_clash_count": 0},
        "geometry_error": None}


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
        def _timed_prepare_corrected(self, **kwargs):
            return self

        def embed(self, seed):
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


def test_preview_is_retained_if_geometry_diagnostics_fail(tmp_path):
    from chematic._timed_conformers import _read_checkpoint, _run_partial

    class GeometryFailure:
        def _timed_prepare_corrected(self, **kwargs):
            return self

        def embed(self, seed):
            return [[0.0, 0.0, 0.0]]

        def _timed_stereo_valid(self, coords):
            return True

        def _timed_preview_geometry(self, coords):
            raise ValueError("diagnostic unavailable")

        def _timed_optimize_embedded(self, coords, **kwargs):
            raise RuntimeError("optimizer interrupted")

    path = tmp_path / "preview.jsonl"
    request = {"num_confs": 1, "random_seed": 42, "max_iterations": 0,
               "max_mmff_iters": 200, "repair_mmff_iters": 0,
               "return_preview": True, "checkpoint": str(path)}
    with pytest.raises(RuntimeError, match="optimizer interrupted"):
        _run_partial(GeometryFailure(), request)
    preview = _read_checkpoint(str(path))["in_progress"]
    assert preview["coordinates"] == [[0.0, 0.0, 0.0]]
    assert preview["geometry"] is None
    assert preview["geometry_error"] == "ValueError: diagnostic unavailable"


def test_preview_geometry_reports_clashes_and_stretched_bonds():
    molecule = chematic.from_smiles("CCC")
    count, worst, rate, clashes = molecule._timed_preview_geometry(
        [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.2, 0.0, 0.0]]
    )
    assert (count, worst, rate, clashes) == (2, 1.5, 0.0, 1)
    count, worst, rate, clashes = molecule._timed_preview_geometry(
        [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [20.0, 0.0, 0.0]]
    )
    assert (count, worst, rate, clashes) == (2, 18.5, 0.5, 0)
    with pytest.raises(ValueError):
        molecule._timed_preview_geometry([[0.0, 0.0, 0.0]])
    with pytest.raises(ValueError):
        molecule._timed_preview_geometry([[float("nan"), 0.0, 0.0]] * 3)


def test_corrupt_preview_diagnostics_are_rejected(tmp_path):
    import json
    from chematic._timed_conformers import _read_checkpoint

    path = tmp_path / "journal.jsonl"
    preview = {"preview_attempt": 1, "coordinates": [[0.0, 0.0, 0.0]],
               "geometry": {"bond_count": -1, "worst_bond_length": None,
                            "bond_violation_rate_50pct": None, "gross_clash_count": 0}}
    path.write_text(json.dumps(preview) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="invalid preview geometry diagnostics"):
        _read_checkpoint(str(path))


@pytest.mark.parametrize("smiles", ["CCO", "CB(C)C", "F[C@](Cl)(Br)I", "CC.O"])
def test_prepared_corrected_embedding_matches_single_attempts(smiles):
    molecule = chematic.from_smiles(smiles).add_hydrogens()
    prepared = molecule._timed_prepare_corrected(max_iterations=0)
    for seed in (0, 2, 42, 84):
        assert prepared.embed(seed) == molecule._timed_embed_corrected(random_seed=seed)
    with pytest.raises(ValueError):
        prepared.embed(-1)


def test_prepared_preview_zero_seed_matches_existing_partial_mode():
    options = {"timeout_seconds": 15, "return_partial": True, "num_confs": 3, "random_seed": 0}
    baseline = chematic.generate_conformers_timed("CCO", **options)
    preview = chematic.generate_conformers_timed("CCO", return_preview=True, **options)
    assert baseline["status"] == preview["status"] == "completed"
    assert baseline["conformers"] == preview["conformers"]
    assert baseline["attempts_completed"] == preview["attempts_completed"] == 3


def test_timed_preview_regular_geometry_screen_rejects_distorted_stereo():
    import json
    from pathlib import Path

    path = Path(__file__).with_name("stereo_distorted_prototype_719.json")
    fixture = json.loads(path.read_text(encoding="utf-8"))
    molecule = chematic.from_smiles(fixture["smiles"])
    assert molecule._timed_stereo_valid(fixture["coords"])
    assert not molecule._timed_stereo_valid(fixture["coords"], require_regular_geometry=True)


@pytest.mark.parametrize("screen_passes", [False, True])
def test_strict_preview_checkpoint_contains_only_screened_coordinates(tmp_path, screen_passes):
    from chematic._timed_conformers import _read_checkpoint, _run_partial

    class FakeMol:
        def _timed_prepare_corrected(self, **kwargs):
            return self

        def embed(self, seed):
            return [[0.0, 0.0, 0.0]]

        def _timed_stereo_valid(self, coords, *, require_regular_geometry=False):
            assert require_regular_geometry is True
            return screen_passes

        def _timed_preview_geometry(self, coords):
            return (0, None, None, 0)

        def _timed_optimize_embedded(self, coords, **kwargs):
            raise RuntimeError("interrupted")

    path = tmp_path / "journal.jsonl"
    request = {"num_confs": 1, "random_seed": 42, "max_iterations": 0,
               "max_mmff_iters": 200, "repair_mmff_iters": 0,
               "return_preview": True, "preview_require_regular_geometry": True,
               "checkpoint": str(path)}
    with pytest.raises(RuntimeError, match="interrupted"):
        _run_partial(FakeMol(), request)
    result = _read_checkpoint(str(path))
    if screen_passes:
        assert result["in_progress"]["coordinates"] == [[0.0, 0.0, 0.0]]
        assert result["in_progress"]["stereo_screen"] == "regular_geometry"
    else:
        assert "in_progress" not in result


def test_invalid_preview_screen_checkpoint_is_rejected(tmp_path):
    import json
    from chematic._timed_conformers import _read_checkpoint

    path = tmp_path / "journal.jsonl"
    entry = {"preview_attempt": 1, "coordinates": [[0., 0., 0.]], "stereo_screen": "invalid"}
    path.write_text(json.dumps(entry) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="invalid preview stereo screen"):
        _read_checkpoint(str(path))


def test_timed_preview_screen_options_and_completed_ensemble():
    with pytest.raises(ValueError, match="requires return_preview"):
        chematic.generate_conformers_timed("CCO", 1, timeout_seconds=1,
                                          preview_require_regular_geometry=True)
    with pytest.raises(TypeError, match="must be a bool"):
        chematic.generate_conformers_timed("CCO", 1, timeout_seconds=1,
                                          preview_require_regular_geometry=1)
    options = {"num_confs": 2, "timeout_seconds": 15, "return_partial": True,
               "return_preview": True, "random_seed": 42}
    baseline = chematic.generate_conformers_timed("CCO", **options)
    checked = chematic.generate_conformers_timed("CCO", preview_require_regular_geometry=True, **options)
    assert baseline["status"] == checked["status"] == "completed"
    assert baseline["conformers"] == checked["conformers"]
    assert baseline["attempts_completed"] == checked["attempts_completed"] == 2


def test_strict_preview_checkpoint_survives_worker_termination(tmp_path):
    import subprocess
    import sys
    import time
    from chematic._timed_conformers import _read_checkpoint

    path = tmp_path / "strict-journal.jsonl"
    script = '''
import sys, time
from chematic._timed_conformers import _run_partial
class SlowMol:
    def _timed_prepare_corrected(self, **kwargs): return self
    def embed(self, seed): return [[0.0, 0.0, 0.0]]
    def _timed_stereo_valid(self, coords, *, require_regular_geometry=False): return require_regular_geometry
    def _timed_preview_geometry(self, coords): return (0, None, None, 0)
    def _timed_optimize_embedded(self, coords, **kwargs):
        time.sleep(30)
        return [], 1, None, None, None, None
_run_partial(SlowMol(), {"num_confs": 1, "random_seed": 42, "max_iterations": 0,
    "max_mmff_iters": 200, "repair_mmff_iters": 0, "return_preview": True,
    "preview_require_regular_geometry": True, "checkpoint": sys.argv[1]})
'''
    worker = subprocess.Popen([sys.executable, "-c", script, str(path)])
    try:
        deadline = time.monotonic() + 10
        while not path.exists() and time.monotonic() < deadline and worker.poll() is None:
            time.sleep(0.01)
        assert path.exists(), "worker did not publish a screened preview"
    finally:
        worker.kill()
        worker.wait()
    result = _read_checkpoint(str(path))
    assert result["attempts_completed"] == 0 and result["conformers"] == []
    assert result["in_progress"]["stereo_screen"] == "regular_geometry"
    assert result["in_progress"]["coordinates"] == [[0., 0., 0.]]
