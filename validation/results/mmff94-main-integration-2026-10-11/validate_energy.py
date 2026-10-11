"""Measure shared-coordinate energies using the transplanted source wheel."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
EXT = next((ROOT / "target/integration-python/chematic").glob("*.so"))
RDKIT = Path("/Users/k_tanabe/Documents/Documents/oss_rust/chematic/.venv-v1gate/lib/python3.13/site-packages/rdkit/rdBase.so")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    env = dict(os.environ, PYTHONPATH=str(ROOT / "target/integration-python") + os.pathsep + str(ROOT))
    before = {str(path): sha(path) for path in (EXT, RDKIT)}
    command = [
        sys.executable, "-B", "scripts/mmff94_same_explicit_h_energy.py",
        "--output", str(OUT / "energy.jsonl"),
        "--summary", str(OUT / "energy-summary.json"),
        "--expected-schematic", "1.1.0", "--expected-rdkit", "2026.03.6",
        "--schematic-artifact", str(EXT), "--rdkit-artifact", str(RDKIT),
        "--schematic-build-command", "CARGO_PROFILE_DEV_DEBUG=0 CARGO_INCREMENTAL=0 maturin build --profile dev --offline -j 2 -m crates/chematic-py/Cargo.toml -i /private/tmp/chematic-mmff-development-20261011/bin/python --out target/integration-wheels -- -C strip=none -C debuginfo=0",
    ]
    with (OUT / "energy.stdout.json").open("w") as out, (OUT / "energy.stderr.log").open("w") as err:
        run = subprocess.run(command, cwd=ROOT, env=env, stdout=out, stderr=err)
    assert {str(path): sha(path) for path in (EXT, RDKIT)} == before
    summary = json.loads((OUT / "energy-summary.json").read_bytes())
    assert run.returncode == 2
    assert {key: value for key, value in summary["gate"].items() if not value} == {"source_tree_clean": False}
    (OUT / "energy-validation.json").write_text(json.dumps({"command":command, "exit_code":run.returncode, "artifact_hashes":before, "expected_nonzero_reason":"uncommitted transplant candidate"},indent=2) + "\n")
    print(json.dumps(summary["measurements"]), flush=True)


if __name__ == "__main__":
    main()
