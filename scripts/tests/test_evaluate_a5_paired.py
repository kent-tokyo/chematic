import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts/evaluate_a5_paired.py"


def run(packet: Path, *args: str):
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(packet), *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )


def write_rows(path: Path, rows: list[dict]) -> None:
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")


def row(case_id: str, cluster: str, schematic: bool, rdkit: bool, adjudicated: bool = True):
    return {
        "id": case_id,
        "cluster": cluster,
        "adjudicated": adjudicated,
        "chematic_correct": schematic,
        "rdkit_correct": rdkit,
        "chematic_status": "ok" if schematic else "refused",
        "rdkit_status": "ok" if rdkit else "refused",
    }


def test_paired_report_retains_unresolved_and_bootstraps(tmp_path):
    packet = tmp_path / "packet.jsonl"
    report = tmp_path / "report.json"
    write_rows(packet, [row("a", "s1", True, True), row("b", "s2", True, False), row("c", "s3", False, False, False)])
    result = run(packet, "--output", str(report), "--repetitions", "200")
    assert result.returncode == 0, result.stderr
    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["category_counts"] == {"both_correct": 1, "chematic_only_correct": 1, "unresolved": 1}
    assert data["all_rows"]["rows"] == 3
    assert data["adjudicated_valid"]["adjudicated_rows"] == 2
    assert data["paired_bootstrap"]["clusters"] == 2
    assert data["verdict"] in {"equivalent", "schematic_superior", "not_equivalent_or_inconclusive"}


def test_rejects_duplicate_ids_and_non_ok_correctness(tmp_path):
    packet = tmp_path / "packet.jsonl"
    write_rows(packet, [row("a", "s1", True, True), row("a", "s2", False, False)])
    result = run(packet)
    assert result.returncode != 0
    assert "duplicate id" in result.stderr

    write_rows(packet, [{**row("b", "s1", True, True), "chematic_status": "timeout"}])
    result = run(packet)
    assert result.returncode != 0
    assert "non-ok chematic status" in result.stderr


def test_empty_or_single_cluster_is_insufficient(tmp_path):
    packet = tmp_path / "packet.jsonl"
    write_rows(packet, [row("a", "s1", True, True)])
    data = json.loads(run(packet, "--repetitions", "200").stdout)
    assert data["verdict"] == "insufficient_evidence"


def test_reports_all_adjudication_categories(tmp_path):
    packet = tmp_path / "packet.jsonl"
    write_rows(
        packet,
        [
            row("both", "s1", True, True),
            row("chematic", "s2", True, False),
            row("rdkit", "s3", False, True),
            row("neither", "s4", False, False),
        ],
    )
    data = json.loads(run(packet, "--repetitions", "200").stdout)
    assert data["category_counts"] == {
        "both_correct": 1,
        "both_incorrect": 1,
        "rdkit_only_correct": 1,
        "chematic_only_correct": 1,
    }
