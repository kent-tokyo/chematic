import hashlib
import json
import sys
from pathlib import Path

import pytest
from rdkit import Chem, rdBase

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
from reserved_timed_shards import merge, prepare


def test_reserved_shards_cover_every_record_and_reject_mismatches(tmp_path):
    sdf = tmp_path / "source.sdf"
    writer = Chem.SDWriter(str(sdf))
    for i, smiles in enumerate(("CCO", "CCN", "c1ccccc1")):
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", str(i))
        writer.write(mol)
    writer.close()
    split = {"source_sha256": hashlib.sha256(sdf.read_bytes()).hexdigest(),
             "rows": [{"index": i, "name": str(i), "split": "holdout",
                       "scaffold_sha256": hashlib.sha256(str(i).encode()).hexdigest()}
                      for i in range(3)]}
    split_path = tmp_path / "split.json"
    split_path.write_text(json.dumps(split))
    directory = tmp_path / "plan"
    manifest = prepare(sdf, split_path, directory, 2)
    assert sum(row["count"] for row in manifest["shards"]) == 3
    bench = Path(__file__).resolve().parents[1] / "benchmarks"
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    paths = []
    for shard in manifest["shards"]:
        result = {"source_sha256": split["source_sha256"], "wheel_sha256": "wheel",
                  "source_commit": "commit", "rdkit": rdBase.rdkitVersion,
                  "profile": "mixed-fallback", "stereo_policy": "strict", "hard_timeout": True,
                  "geometry_checked": True, "seconds": 0.2, "max_candidates": 20,
                  "repair_iters": 400, "script_sha256": digest(bench / "platinum_time_budget.py"),
                  "stereo_verifier_sha256": digest(bench / "platinum_stereo_strict.py"),
                  "platform": "test", "chematic": "candidate", "results": shard["records"]}
        path = tmp_path / f"result-{shard['index']}.json"
        path.write_text(json.dumps(result))
        paths.append(path)
    output = tmp_path / "merged.json"
    merge(directory / "manifest.json", paths, output, "wheel", "commit")
    assert [row["index"] for row in json.loads(output.read_text())["results"]] == [0, 1, 2]
    bad = json.loads(paths[0].read_text())
    bad["wheel_sha256"] = "different"
    paths[0].write_text(json.dumps(bad))
    with pytest.raises(ValueError, match="settings differ"):
        merge(directory / "manifest.json", paths, output, "wheel", "commit")


def test_mixed_platforms_require_opt_in_and_preserve_provenance(tmp_path):
    bench = Path(__file__).resolve().parents[1] / "benchmarks"
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {"source_sha256": "source", "split_sha256": "split", "rdkit_version": "version",
                "shard_count": 2, "shards": [{"count": 1, "records": [{"index": i, "name": str(i)}]} for i in range(2)]}
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    paths = []
    for i, platform in enumerate(("Linux-7.0-azure-x86_64-with-glibc2.43", "Linux-7.0-aws-x86_64-with-glibc2.43")):
        data = {"source_sha256": "source", "wheel_sha256": "wheel", "source_commit": "commit",
                "rdkit": "version", "profile": "mixed-fallback", "stereo_policy": "strict",
                "hard_timeout": True, "geometry_checked": True, "seconds": 0.2,
                "max_candidates": 20, "repair_iters": 400,
                "script_sha256": digest(bench / "platinum_time_budget.py"),
                "stereo_verifier_sha256": digest(bench / "platinum_stereo_strict.py"),
                "platform": platform, "chematic": "candidate",
                "results": [{"index": i, "name": str(i)}]}
        path = tmp_path / f"result-{i}.json"
        path.write_text(json.dumps(data))
        paths.append(path)
    out = tmp_path / "merged.json"
    with pytest.raises(ValueError, match="scorer mismatch"):
        merge(manifest_path, paths, out, "wheel", "commit")
    merge(manifest_path, paths, out, "wheel", "commit", allow_mixed_platforms=True)
    merged = json.loads(out.read_text())
    assert len(merged["platforms"]) == 2
    assert [r["evaluation_platform"] for r in merged["results"]] == [json.loads(p.read_text())["platform"] for p in paths]
    changed = json.loads(paths[1].read_text())
    changed["platform"] = "Linux-7.0-aws-aarch64-with-glibc2.43"
    paths[1].write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="mixed operating systems or architectures"):
        merge(manifest_path, paths, out, "wheel", "commit", allow_mixed_platforms=True)
