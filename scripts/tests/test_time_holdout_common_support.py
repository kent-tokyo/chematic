import hashlib
import json
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks" / "summarize_platinum_time_holdout.py"


def test_common_mmff_convergence_uses_only_audited_records(tmp_path):
    split = {"source_sha256": "source", "rows": [
        {"index": i, "name": str(i), "scaffold_sha256": str(i), "split": "holdout"}
        for i in range(3)
    ]}
    split_path = tmp_path / "split.json"
    split_path.write_text(json.dumps(split))
    hit = {"rmsd": [0.2], "stereo": [{"atoms": [], "bonds": []}], "statuses": [0]}
    miss = {"rmsd": [], "stereo": [], "statuses": []}
    data = {"source_sha256": "source", "seconds": 1.0, "geometry_checked": False,
            "results": [{"index": i, "name": str(i), "chematic": hit,
                         "rdkit": hit if i == 1 else miss} for i in range(3)]}
    data_path = tmp_path / "result.json"
    data_path.write_text(json.dumps(data))
    support = {"source_sha256": "source", "partition": "holdout",
               "split_sha256": hashlib.sha256(split_path.read_bytes()).hexdigest(),
               "results": [{"index": i, "name": str(i), "rdkit_mmff": i < 2,
                            "chematic_mmff": True, "graph_matches": True} for i in range(3)]}
    support_path = tmp_path / "support.json"
    support_path.write_text(json.dumps(support))
    groups = {"source_sha256": "source", "split_sha256": support["split_sha256"],
              "rows": [{"index": i, "name": str(i), "scaffold_sha256": str(i),
                        "group": "ordinary" if i < 2 else "macrocycle"} for i in range(3)]}
    groups_path = tmp_path / "groups.json"
    groups_path.write_text(json.dumps(groups))
    output_path = tmp_path / "summary.json"
    subprocess.run([sys.executable, str(SCRIPT), str(data_path), str(split_path),
                    "--support", str(support_path), "--groups", str(groups_path), "--output", str(output_path),
                    "--bootstrap", "100"], check=True, capture_output=True)
    result = json.loads(output_path.read_text())
    assert result["common_mmff_records"] == 2
    assert result["stages"]["common_mmff_converged_recovery"]["difference"] == 0.5
    assert result["stages"]["ordinary_recovery"]["records"] == 2
    assert result["stages"]["ordinary_common_mmff_converged_recovery"]["difference"] == 0.5
    assert result["stages"]["macrocycle_recovery"]["records"] == 1
    strict = subprocess.run([sys.executable, str(SCRIPT), str(data_path), str(split_path),
                             "--require-strict", "--output", str(output_path)], capture_output=True)
    assert strict.returncode != 0
    support["results"].pop()
    support_path.write_text(json.dumps(support))
    failed = subprocess.run([sys.executable, str(SCRIPT), str(data_path), str(split_path),
                             "--support", str(support_path), "--output", str(output_path)],
                            capture_output=True)
    assert failed.returncode != 0
