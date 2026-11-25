import multiprocessing

import pytest

from chematic import timed_conformer_batch
from chematic._timed_conformers import _stereo_checked


def test_timeout_replaces_worker_and_preserves_following_results():
    before = {child.pid for child in multiprocessing.active_children()}
    results = timed_conformer_batch(["C" * 300, "CCO", "invalid", "CCO"], seconds=0.2)
    assert [result["status"] for result in results] == ["timeout", "ok", "error", "ok"]
    assert results[0]["conformers"] == []
    assert len(results[1]["conformers"]) == len(results[3]["conformers"]) == 1
    assert {child.pid for child in multiprocessing.active_children()} == before


def test_timeout_keeps_completed_candidates():
    result = timed_conformer_batch(["CCO"], seconds=1.0, num_confs=10000)[0]
    assert result["status"] == "timeout"
    assert 0 < len(result["conformers"]) < 10000
    assert all(len(candidate) == 3 for candidate in result["conformers"])


@pytest.mark.parametrize("seconds", [0, -1, float("inf"), float("nan")])
def test_invalid_deadline(seconds):
    with pytest.raises(ValueError):
        timed_conformer_batch(["CCO"], seconds=seconds)


def test_invalid_input_cleans_up_worker():
    before = {child.pid for child in multiprocessing.active_children()}
    with pytest.raises(TypeError):
        timed_conformer_batch(["CCO", None], seconds=2)
    assert {child.pid for child in multiprocessing.active_children()} == before


def test_unoptimized_fallback_is_explicit_and_retains_provenance():
    smiles = "COc1ccc(C2=Cc3cc(Cl)nc[n+]3[B-](c3ccccc3)(c3ccccc3)O2)cc1"
    unsupported = timed_conformer_batch([smiles], seconds=5)[0]
    assert unsupported["status"] == "error"
    assert "type assignment failed" in unsupported["error"]
    result = timed_conformer_batch([smiles], seconds=5, allow_unoptimized=True)[0]
    assert result["status"] == "ok"
    assert "type assignment failed" in result["fallback_reason"]
    assert len(result["conformers"]) == 1
    status, energy, coords = result["conformers"][0]
    assert status == -1 and energy is None and coords


def test_rejects_center_outside_its_ligand_tetrahedron():
    class FakeMolecule:
        bond_table = [(0, i, "SINGLE", False) for i in range(1, 5)]

        def check_geometry(self, coords):
            return {"stereo_groups_not_evaluated": False,
                    "stereo": {"n_violations": 0, "n_unevaluable": 0,
                               "tetrahedral": [{"atom": 0}]}}

    ligands = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert _stereo_checked(FakeMolecule(), [[0.2, 0.2, 0.2], *ligands])
    assert not _stereo_checked(FakeMolecule(), [[0.4, 0.4, 0.4], *ligands])


def _deadline_test_worker(connection):
    """A deterministic worker that would write after the deadline if left alive."""
    import time
    from pathlib import Path

    connection.send(None)
    while True:
        request = connection.recv()
        if request is None:
            return
        text = request[0]
        if text.startswith("slow:"):
            connection.send(("candidate", [(0, 0.0, [[0.0, 0.0, 0.0]])]))
            time.sleep(0.4)
            Path(text.removeprefix("slow:")).write_text("worker continued")
            connection.send(("candidate", [(0, 1.0, [[1.0, 0.0, 0.0]])]))
        else:
            connection.send(("candidate", [(0, 2.0, [[2.0, 0.0, 0.0]])]))
        connection.send(("ok", None))


def test_deadline_discards_late_result_and_stops_worker(monkeypatch, tmp_path):
    import time
    import chematic._timed_conformers as timed

    marker = tmp_path / "continued"
    before = {child.pid for child in multiprocessing.active_children()}
    monkeypatch.setattr(timed, "_worker", _deadline_test_worker)
    results = timed.timed_conformer_batch([f"slow:{marker}", "next"], seconds=0.1)
    assert [result["status"] for result in results] == ["timeout", "ok"]
    assert [candidate[1] for candidate in results[0]["conformers"]] == [0.0]
    assert [candidate[1] for candidate in results[1]["conformers"]] == [2.0]
    time.sleep(0.45)
    assert not marker.exists()
    assert {child.pid for child in multiprocessing.active_children()} == before
