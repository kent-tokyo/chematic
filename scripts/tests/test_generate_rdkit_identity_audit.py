import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts/generate_rdkit_identity_audit.py"


def test_fragment_parent_ring_cache_is_initialized_before_scaffold_generation(tmp_path):
    source = tmp_path / "cu-complex.smi"
    source.write_text("CCC1=[O+][Cu]2([O+]=C(CC)C1)[O+]=C(CC)CC(=[O+]2)CC\n", encoding="utf-8")
    output = tmp_path / "identity.jsonl"
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), str(source), "--output", str(output)],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    row = json.loads(output.read_text(encoding="utf-8"))
    assert row["input_smiles"] == source.read_text(encoding="utf-8").strip()
    assert isinstance(row["scaffold_smiles"], str)


def test_allow_parse_errors_retains_an_explicit_empty_key_row(tmp_path):
    source = tmp_path / "mixed.smi"
    source.write_text("CCO\nC1(CC\n", encoding="utf-8")
    output = tmp_path / "identity.jsonl"
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), str(source), "--output", str(output), "--allow-parse-errors"],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    rows = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]
    assert len(rows) == 2
    assert rows[1]["input_smiles"] == "C1(CC"
    assert rows[1]["parse_error"] == "RDKit cannot parse input"
    assert {key: rows[1][key] for key in ("canonical_smiles", "parent_smiles", "scaffold_smiles")} == {
        "canonical_smiles": "",
        "parent_smiles": "",
        "scaffold_smiles": "",
    }
