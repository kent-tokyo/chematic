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
