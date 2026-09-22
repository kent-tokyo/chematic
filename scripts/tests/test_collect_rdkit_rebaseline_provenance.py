import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

from scripts.collect_rdkit_rebaseline_provenance import (
    collect_npm_lane,
    collect_python_lane,
)


class FakeDistribution:
    version = "2026.3.6"
    metadata = {"Name": "rdkit"}

    def read_text(self, name: str) -> str | None:
        values = {
            "METADATA": "Name: rdkit\nVersion: 2026.3.6\n",
            "WHEEL": "Tag: cp313-cp313-macosx_11_0_arm64\n",
            "RECORD": "rdkit/__init__.py,sha256=fake,1\n",
            "direct_url.json": json.dumps(
                {"archive_info": {"hashes": {"sha256": "a" * 64}}}
            ),
        }
        return values.get(name)


def test_python_lane_keeps_runtime_and_backend_evidence_separate():
    def importer(name: str) -> object:
        if name == "rdkit":
            return SimpleNamespace(__file__="/tmp/rdkit/__init__.py")
        if name == "rdkit.rdBase":
            return SimpleNamespace(rdkitVersion="2026.03.6")
        raise ImportError(name)

    lane = collect_python_lane(
        "nanobind",
        "distributed extension symbol audit",
        distribution_lookup=lambda _name: FakeDistribution(),
        module_import=importer,
    )

    assert lane["availability"] == "measured"
    assert lane["package"]["artifact_archive_sha256"] == "a" * 64
    assert lane["runtime"]["version"] == "2026.03.6"
    assert lane["wrapper_backend"] == {
        "value": "nanobind",
        "evidence": "distributed extension symbol audit",
        "status": "declared",
    }
    assert lane["missing_dimensions"] == []


def test_python_lane_does_not_promote_installed_record_to_archive_hash():
    distribution = FakeDistribution()
    distribution.read_text = lambda name: (
        None if name == "direct_url.json" else FakeDistribution().read_text(name)
    )
    lane = collect_python_lane(
        "unknown",
        None,
        distribution_lookup=lambda _name: distribution,
        module_import=lambda _name: SimpleNamespace(
            __file__=None, rdkitVersion="2026.03.6"
        ),
    )

    assert lane["package"]["artifact_archive_sha256"] is None
    assert "artifact_archive_sha256" in lane["missing_dimensions"]
    assert "wrapper_backend" in lane["missing_dimensions"]


def test_npm_lane_records_lock_integrity_assets_and_explicit_runtime(tmp_path: Path):
    node_modules = tmp_path / "node_modules"
    package_dir = node_modules / "@rdkit" / "rdkit"
    dist = package_dir / "dist"
    dist.mkdir(parents=True)
    (package_dir / "package.json").write_text(
        json.dumps({"name": "@rdkit/rdkit", "version": "2026.3.6"}),
        encoding="utf-8",
    )
    wasm = dist / "RDKit_minimal.wasm"
    wasm.write_bytes(b"wasm")
    lock = {
        "packages": {
            "node_modules/@rdkit/rdkit": {
                "resolved": "https://registry.example/rdkit.tgz",
                "integrity": "sha512-example",
            }
        }
    }
    lock_path = node_modules / ".package-lock.json"
    lock_path.write_text(json.dumps(lock), encoding="utf-8")

    lane = collect_npm_lane(package_dir, "2026.03.6")

    assert lane["availability"] == "measured"
    assert lane["runtime"] == {"version": "2026.03.6", "status": "measured"}
    assert lane["package"]["registry"]["integrity"] == "sha512-example"
    assert lane["package"]["assets"] == [
        {
            "path": "dist/RDKit_minimal.wasm",
            "bytes": 4,
            "sha256": hashlib.sha256(b"wasm").hexdigest(),
        }
    ]
    assert lane["missing_dimensions"] == []


def test_npm_lane_never_infers_runtime_from_package_version(tmp_path: Path):
    package_dir = tmp_path / "node_modules" / "@rdkit" / "rdkit"
    package_dir.mkdir(parents=True)
    (package_dir / "package.json").write_text(
        json.dumps({"name": "@rdkit/rdkit", "version": "2026.3.6"}),
        encoding="utf-8",
    )

    lane = collect_npm_lane(package_dir, None)

    assert lane["availability"] == "installed_not_executed"
    assert lane["runtime"] == {"version": None, "status": "not_measured"}
    assert "runtime" in lane["missing_dimensions"]
