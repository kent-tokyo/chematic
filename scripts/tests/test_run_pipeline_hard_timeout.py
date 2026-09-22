import json

from scripts import run_pipeline_hard_timeout as runner


def test_run_one_preserves_original_index(tmp_path, monkeypatch):
    manifest = tmp_path / "tier-a.json"
    manifest.write_text(
        json.dumps(
            {
                "molecules": [
                    {
                        "name": "ethane",
                        "smiles": "CC",
                        "primary_category": "small",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    binary = tmp_path / "fake-pipeline"
    binary.write_text(
        "#!/bin/sh\nprintf '%s\\n' '{\"status\":\"success\",\"name\":\"ethane\"}'\n",
        encoding="utf-8",
    )
    binary.chmod(0o755)
    monkeypatch.setattr(runner, "ROOT", tmp_path)
    monkeypatch.setattr(runner, "MANIFESTS", {"A": manifest})

    row = runner.run_one(binary, "A", 0, "test-arm", 1.0)

    assert row["status"] == "success"
    assert row["original_index"] == 0
    assert row["external_timeout_s"] == 1.0


def test_metadata_records_binary_manifest_and_runtime(tmp_path, monkeypatch):
    manifest = tmp_path / "tier-a.json"
    manifest.write_text('{"molecules": [{"name": "ethane", "smiles": "CC"}]}\n')
    binary = tmp_path / "pipeline"
    binary.write_bytes(b"candidate-binary")
    monkeypatch.setattr(runner, "ROOT", tmp_path)
    monkeypatch.setattr(runner, "MANIFESTS", {"A": manifest})
    monkeypatch.setenv("SCHEMATIC_MMFF94_MAX_ITERATIONS", "300")

    metadata = runner.build_metadata(binary, "test-arm", 25.0, ["A"])

    assert metadata["schema_version"] == 1
    assert metadata["binary"]["sha256"] == runner.sha256_file(binary)
    assert metadata["manifests"]["A"]["row_count"] == 1
    assert metadata["manifests"]["A"]["sha256"] == runner.sha256_file(manifest)
    assert metadata["environment"]["SCHEMATIC_MMFF94_MAX_ITERATIONS"] == "300"
    assert metadata["source_revision"] is None
