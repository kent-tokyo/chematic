import json
import subprocess
import sys

from pathlib import Path


ROOT = Path(__file__).parents[2]
VALIDATOR = ROOT / "scripts/check_descriptor_raw_accounting.py"
SOURCE = ROOT / "validation/results/descriptor-rdkit-diagnostics-v1.0.13-candidate-raw.json"


def test_raw_accounting_accepts_current_candidate():
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(SOURCE)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_raw_accounting_rejects_stale_summary(tmp_path):
    document = json.loads(SOURCE.read_text(encoding="utf-8"))
    document["fields"]["hba"]["strict_matches"] -= 1
    path = tmp_path / "stale.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(VALIDATOR), str(path)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "hba.strict_matches" in result.stderr
