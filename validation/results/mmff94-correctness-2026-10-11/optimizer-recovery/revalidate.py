"""Recheck the final source extension and score both candidates with exact JSON decoding."""
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
EXTENSION = Path("crates/chematic-py/python/chematic/chematic.cpython-313-darwin.so")
BINARY = Path("target/debug/examples/pipeline_v2_vs_rdkit_common_scorer")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    source = json.loads((OUT / "build-source.json").read_bytes())
    for record in [source, *source["additional_sources"]]:
        assert digest(record["path"]) == record["sha256"]
    extension_hash = digest(EXTENSION)
    scorer = json.loads((OUT / "scorer-build.json").read_bytes())["binary"]
    assert digest(BINARY) == scorer["sha256"]
    for filename in ["rust-tests.log", "clippy.log", "build.log", "scorer-build.log"]:
        shutil.copyfile(Path("/private/tmp/chematic-mmff-development-20261011") / ("recovery-" + filename), OUT / filename)
    commands = [
        ("python-tests.log", [sys.executable, "-B", "-m", "pytest",
            "scripts/tests/test_mmff94_same_explicit_h_energy.py",
            "crates/chematic-py/tests/test_pipeline_v2.py",
            "crates/chematic-py/tests/test_conformer_ensemble.py",
            "crates/chematic-py/tests/test_conformer_ensemble_v2.py", "-p", "no:cacheprovider", "-q"]),
        ("holdout.log", [sys.executable, "-B", str(OUT / "measure_holdout.py"), "candidate"]),
        ("energy-run.log", [sys.executable, "-B", str(OUT / "run_frozen_energy.py")]),
        ("quality-run.log", [sys.executable, "-B", str(OUT.parent / "run_quality.py"),
                             "--output", str(OUT / "quality.jsonl")]),
        ("gradient-run.log", [sys.executable, "-B", str(OUT / "compare_gradients.py")]),
        ("scorer.log", [str(BINARY), "--chematic-rows", str(OUT / "quality.jsonl"),
            "--rdkit-rows", "validation/results/public-package-3d-rdkit-2026.3.6-2026-09-22.jsonl",
            "--output", str(OUT / "quality-scored.jsonl")]),
    ]
    completed = []
    env = dict(os.environ, PYTHONPATH=".")
    for logfile, command in commands:
        with (OUT / logfile).open("w") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env)
        completed.append({"log": logfile, "command": command, "exit_code": result.returncode})
        print(logfile, result.returncode, flush=True)
        assert result.returncode == 0, logfile
    assert digest(EXTENSION) == extension_hash
    assert digest(BINARY) == scorer["sha256"]
    for record in [source, *source["additional_sources"]]:
        assert digest(record["path"]) == record["sha256"]
    (OUT / "validation-run.json").write_text(json.dumps({
        "source": source, "extension_sha256_before_and_after": extension_hash,
        "scorer_binary": scorer, "commands": completed[:4],
    }, indent=2) + "\n")
    (OUT / "gradient-scorer-run.json").write_text(json.dumps(completed[4:], indent=2) + "\n")


if __name__ == "__main__":
    main()
