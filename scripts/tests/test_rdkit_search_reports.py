import json
from pathlib import Path

from scripts.check_rdkit_search_reports import validate


ROOT = Path(__file__).parents[2]
REPORTS = (
    ROOT / "validation/results/rdkit-search-oracle-current-rerun.json",
    ROOT / "validation/results/rdkit-search-threshold-gate-v1.0.13.json",
    ROOT / "validation/results/rdkit-search-cross-binding-parity-v1.0.13.json",
)


def copied_reports(tmp_path: Path) -> list[Path]:
    paths = []
    for index, source in enumerate(REPORTS):
        target = tmp_path / f"report-{index}.json"
        target.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
        paths.append(target)
    return paths


def test_search_report_packet_accepts_current_reports(tmp_path):
    paths = copied_reports(tmp_path)
    assert validate(tuple(paths)) == []


def test_search_report_packet_rejects_corpus_mismatch(tmp_path):
    paths = copied_reports(tmp_path)
    data = json.loads(paths[1].read_text(encoding="utf-8"))
    data["corpus"]["sha256"] = "0" * 64
    paths[1].write_text(json.dumps(data), encoding="utf-8")
    errors = validate(tuple(paths))
    assert any("same corpus identity" in error for error in errors)


def test_search_report_packet_rejects_provenance_mismatch(tmp_path):
    paths = copied_reports(tmp_path)
    data = json.loads(paths[2].read_text(encoding="utf-8"))
    data["provenance"]["tracked_diff_sha256"] = "f" * 64
    paths[2].write_text(json.dumps(data), encoding="utf-8")
    errors = validate(tuple(paths))
    assert any("same source/diff provenance" in error for error in errors)


def test_search_report_packet_rejects_nested_mismatch(tmp_path):
    paths = copied_reports(tmp_path)
    data = json.loads(paths[0].read_text(encoding="utf-8"))
    data["comparison"]["1"]["mismatches"] = 1
    paths[0].write_text(json.dumps(data), encoding="utf-8")
    errors = validate(tuple(paths))
    assert any("mismatch count" in error for error in errors)
