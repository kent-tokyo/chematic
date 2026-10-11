"""Validate the transplanted main candidate in its own unoptimized source wheel."""
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
BUILD = Path("/private/tmp/chematic-main-integration-20261011")


def artifact(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    return {"path": str(path), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def main():
    wheels = list((ROOT / "target/integration-wheels").glob("*.whl"))
    assert len(wheels) == 1
    wheel = wheels[0]
    raw = wheel.read_bytes()
    wheel_hash = hashlib.sha256(raw).hexdigest()
    extract = ROOT / "target/integration-python"
    extract.mkdir(exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        for name in archive.namelist():
            member = PurePosixPath(name)
            assert not member.is_absolute() and ".." not in member.parts
        archive.extractall(extract)
    env = dict(os.environ, PYTHONPATH=str(extract) + os.pathsep + str(ROOT))
    version_check = """import importlib.metadata, json, pathlib, chematic
from rdkit import rdBase
assert importlib.metadata.version("chematic") == "1.1.0"
assert pathlib.Path(chematic.__file__).is_relative_to(pathlib.Path("target/integration-python").resolve())
assert rdBase.rdkitVersion == "2026.03.6"
print(json.dumps({"schematic_version":importlib.metadata.version("chematic"),"schematic_module":chematic.__file__,"rdkit_version":rdBase.rdkitVersion,"rdkit_binary":rdBase.__file__}))
"""
    commands = [
        ("imports.json", [sys.executable, "-B", "-c", version_check]),
        ("python-tests.log", [sys.executable, "-B", "-m", "pytest",
            "scripts/tests/test_mmff94_same_explicit_h_energy.py",
            "crates/chematic-py/tests/test_pipeline_v2.py",
            "crates/chematic-py/tests/test_conformer_ensemble.py",
            "crates/chematic-py/tests/test_conformer_ensemble_v2.py",
            "-p", "no:cacheprovider", "-q"]),
        ("quality-run.log", [sys.executable, "-B",
            "validation/results/mmff94-correctness-2026-10-11/run_quality.py",
            "--output", str(OUT / "quality.jsonl")]),
        ("scorer.log", [str(ROOT / "target/debug/examples/pipeline_v2_vs_rdkit_common_scorer"),
            "--chematic-rows", str(OUT / "quality.jsonl"),
            "--rdkit-rows", "validation/results/public-package-3d-rdkit-2026.3.6-2026-09-22.jsonl",
            "--output", str(OUT / "quality-scored.jsonl")]),
    ]
    baseline = json.loads((ROOT / "validation/results/mmff94-correctness-2026-10-11/optimizer-recovery/candidate-summary.json").read_bytes())
    for manifest in baseline["manifests"]:
        current = ROOT / "validation/manifests" / Path(manifest["path"]).name
        assert artifact(current)["sha256"] == manifest["sha256"]
    completed = []
    binary = ROOT / "target/debug/examples/pipeline_v2_vs_rdkit_common_scorer"
    binary_before = artifact(binary)
    for logfile, command in commands:
        with (OUT / logfile).open("w") as log:
            result = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        completed.append({"command": command, "exit_code": result.returncode, "log": logfile})
        print(logfile, result.returncode, flush=True)
        assert result.returncode == 0, logfile
    assert artifact(binary) == binary_before
    assert hashlib.sha256(wheel.read_bytes()).hexdigest() == wheel_hash
    imports = json.loads((OUT / "imports.json").read_bytes())
    assert artifact(imports["rdkit_binary"])["sha256"] == baseline["artifacts"]["rdkit"]["sha256"]
    quality = [json.loads(line) for line in (OUT / "quality.jsonl").read_bytes().splitlines()]
    scored = [json.loads(line) for line in (OUT / "quality-scored.jsonl").read_bytes().splitlines()]
    scored = [row for row in scored if row["engine"] == "chematic"]
    assert len(quality) == len(scored) == 265
    counts = {
        "inputs": len(quality),
        "success": sum(row["status"] == "success" for row in quality),
        "independently_sound": sum(row.get("independently_sound", False) for row in scored),
        "stereo_violated": sum(row.get("stereo", {}).get("violated", 0) for row in scored),
        "stereo_unevaluable": sum(row.get("stereo", {}).get("unevaluable", 0) for row in scored),
        "reported_converged": sum(row.get("force_field", {}).get("converged", False) for row in quality),
    }
    (OUT / "python-validation.json").write_text(json.dumps({
        "candidate_kind": "transplant_on_main_v1.1.0",
        "base_revision": subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "transplanted_original_commit": "0bae2e0bdec5ef299d1c4dd6e16e1c06cb164c61",
        "build_profile": "dev, unoptimized, debug=0, incremental=0; no runtime comparison",
        "imports": imports,
        "wheel": {"path":str(wheel),"bytes":len(raw),"sha256":wheel_hash},
        "scorer": binary_before,
        "rdkit": artifact(imports["rdkit_binary"]),
        "commands": completed,
        "quality": counts,
    }, indent=2) + "\n")
    print(json.dumps(counts), flush=True)
    assert counts["success"] == counts["independently_sound"] == 265
    assert counts["stereo_violated"] == counts["stereo_unevaluable"] == 0


if __name__ == "__main__":
    main()
