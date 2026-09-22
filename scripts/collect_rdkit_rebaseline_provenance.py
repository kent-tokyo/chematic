#!/usr/bin/env python3
"""Collect offline provenance for one RDKit rebaseline environment.

The collector deliberately separates installed-package identity from the
original wheel or npm tarball identity.  Missing archive hashes, runtime
versions, and wrapper-backend evidence remain explicit instead of being
inferred from a package version or an upstream release announcement.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.validation_provenance import collect as collect_source_provenance


PYTHON_BACKENDS = ("unknown", "boost_python", "nanobind")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def file_record(path: Path, root: Path) -> dict[str, object]:
    resolved = path.resolve()
    try:
        label = str(resolved.relative_to(root.resolve()))
    except ValueError:
        label = str(resolved)
    body = resolved.read_bytes()
    return {"path": label, "bytes": len(body), "sha256": sha256_bytes(body)}


def command_version(command: str) -> str | None:
    try:
        result = subprocess.run(
            [command, "--version"],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=10,
        )
    except (FileNotFoundError, subprocess.SubprocessError):
        return None
    return result.stdout.strip() or None


def distribution_text_record(
    distribution: object, name: str
) -> dict[str, object] | None:
    text = distribution.read_text(name)  # type: ignore[attr-defined]
    if text is None:
        return None
    body = text.encode("utf-8")
    return {"name": name, "bytes": len(body), "sha256": sha256_bytes(body)}


def direct_url_archive_hash(distribution: object) -> str | None:
    text = distribution.read_text("direct_url.json")  # type: ignore[attr-defined]
    if text is None:
        return None
    try:
        document = json.loads(text)
    except json.JSONDecodeError:
        return None
    archive = document.get("archive_info")
    if not isinstance(archive, dict):
        return None
    hashes = archive.get("hashes")
    if isinstance(hashes, dict) and isinstance(hashes.get("sha256"), str):
        return hashes["sha256"]
    legacy = archive.get("hash")
    if isinstance(legacy, str) and legacy.startswith("sha256="):
        return legacy.removeprefix("sha256=")
    return None


def collect_python_lane(
    wrapper_backend: str,
    backend_evidence: str | None,
    artifact_path: Path | None = None,
    *,
    distribution_lookup: Callable[[str], object] = importlib.metadata.distribution,
    module_import: Callable[[str], object] = importlib.import_module,
) -> dict[str, object]:
    try:
        distribution = distribution_lookup("rdkit")
    except importlib.metadata.PackageNotFoundError:
        return {
            "channel": "python_native",
            "availability": "unavailable",
            "reason": "rdkit distribution is not installed",
            "missing_dimensions": [
                "package",
                "artifact_archive_sha256",
                "runtime",
                "wrapper_backend",
            ],
        }

    metadata = getattr(distribution, "metadata", {})
    package_name = metadata.get("Name", "rdkit")
    package_version = getattr(distribution, "version", None)
    metadata_records = [
        record
        for name in ("METADATA", "WHEEL", "RECORD", "direct_url.json")
        if (record := distribution_text_record(distribution, name)) is not None
    ]
    direct_url_sha256 = direct_url_archive_hash(distribution)
    artifact = file_record(artifact_path, ROOT) if artifact_path is not None else None
    archive_sha256 = (
        str(artifact["sha256"]) if artifact is not None else direct_url_sha256
    )

    runtime_version: str | None = None
    runtime_error: str | None = None
    try:
        rdkit_module = module_import("rdkit")
        rd_base = module_import("rdkit.rdBase")
        runtime_version = getattr(rd_base, "rdkitVersion", None)
        module_path = getattr(rdkit_module, "__file__", None)
    except (ImportError, OSError) as exc:
        runtime_error = f"{type(exc).__name__}: {exc}"
        module_path = None

    missing_dimensions: list[str] = []
    if archive_sha256 is None:
        missing_dimensions.append("artifact_archive_sha256")
    if runtime_version is None:
        missing_dimensions.append("runtime")
    if wrapper_backend == "unknown":
        missing_dimensions.append("wrapper_backend")

    return {
        "channel": "python_native",
        "availability": (
            "measured" if runtime_version is not None else "installed_not_executed"
        ),
        "package": {
            "name": package_name,
            "version": package_version,
            "installed_metadata": metadata_records,
            "artifact": artifact,
            "artifact_archive_sha256": archive_sha256,
            "direct_url_archive_sha256": direct_url_sha256,
        },
        "runtime": {
            "version": runtime_version,
            "module_path": module_path,
            "error": runtime_error,
        },
        "wrapper_backend": {
            "value": None if wrapper_backend == "unknown" else wrapper_backend,
            "evidence": backend_evidence,
            "status": "unverified" if wrapper_backend == "unknown" else "declared",
        },
        "missing_dimensions": missing_dimensions,
    }


def find_node_modules(path: Path) -> Path | None:
    return next(
        (parent for parent in (path, *path.parents) if parent.name == "node_modules"),
        None,
    )


def collect_npm_lane(
    package_dir: Path,
    runtime_version: str | None,
    artifact_path: Path | None = None,
) -> dict[str, object]:
    package_dir = package_dir.resolve()
    package_json = package_dir / "package.json"
    if not package_json.is_file():
        return {
            "channel": "npm_wasm",
            "availability": "unavailable",
            "reason": f"missing package.json at {package_json}",
            "missing_dimensions": ["package", "registry_integrity", "runtime"],
        }

    package = json.loads(package_json.read_text(encoding="utf-8"))
    package_name = package.get("name")
    package_version = package.get("version")
    node_modules = find_node_modules(package_dir)
    lock_record: dict[str, object] = {"status": "not_found"}
    if node_modules is not None:
        lock_path = node_modules / ".package-lock.json"
        if lock_path.is_file():
            try:
                lock = json.loads(lock_path.read_text(encoding="utf-8"))
                entry = lock.get("packages", {}).get(f"node_modules/{package_name}")
            except (json.JSONDecodeError, OSError):
                entry = None
                lock_record = {"status": "unreadable"}
            if isinstance(entry, dict):
                lock_record = {
                    "status": "recorded",
                    "path": str(lock_path),
                    "sha256": sha256_bytes(lock_path.read_bytes()),
                    "resolved": entry.get("resolved"),
                    "integrity": entry.get("integrity"),
                }
            elif lock_record["status"] == "not_found":
                lock_record = {"status": "entry_not_found"}

    assets = [
        file_record(path, package_dir)
        for path in sorted(package_dir.rglob("*"))
        if path.is_file() and path.suffix in {".wasm", ".js", ".mjs", ".ts"}
    ]
    artifact = file_record(artifact_path, ROOT) if artifact_path is not None else None
    missing_dimensions: list[str] = []
    if artifact is None:
        missing_dimensions.append("artifact_archive_sha256")
    if lock_record.get("integrity") is None:
        missing_dimensions.append("registry_integrity")
    if runtime_version is None:
        missing_dimensions.append("runtime")

    return {
        "channel": "npm_wasm",
        "availability": (
            "measured" if runtime_version is not None else "installed_not_executed"
        ),
        "package": {
            "name": package_name,
            "version": package_version,
            "package_json": file_record(package_json, package_dir),
            "artifact": artifact,
            "artifact_archive_sha256": (
                artifact["sha256"] if artifact is not None else None
            ),
            "registry": lock_record,
            "assets": assets,
        },
        "runtime": {
            "version": runtime_version,
            "status": "measured" if runtime_version is not None else "not_measured",
        },
        "wrapper_backend": {
            "value": "official_rdkit_js" if package_name == "@rdkit/rdkit" else None,
            "status": (
                "package_identity" if package_name == "@rdkit/rdkit" else "unverified"
            ),
        },
        "missing_dimensions": missing_dimensions,
    }


def collect_file_group(paths: list[Path], root: Path) -> list[dict[str, object]]:
    return [file_record(path, root) for path in paths]


def build_record(args: argparse.Namespace) -> dict[str, object]:
    lanes: list[dict[str, object]] = []
    if args.python:
        lanes.append(
            collect_python_lane(
                args.python_wrapper_backend,
                args.python_backend_evidence,
                args.python_artifact,
            )
        )
    if args.npm_package is not None:
        lanes.append(
            collect_npm_lane(
                args.npm_package,
                args.npm_runtime_version,
                args.npm_artifact,
            )
        )
    return {
        "schema_version": 1,
        "profile": "rdkit_rebaseline_provenance_v1",
        "collected_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": collect_source_provenance(ROOT),
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version,
            "python_executable": sys.executable,
            "node": command_version("node"),
            "npm": command_version("npm"),
            "rustc": command_version("rustc"),
            "cargo": command_version("cargo"),
        },
        "inputs": {
            "corpora": collect_file_group(args.corpus, ROOT),
            "operation_configurations": collect_file_group(args.operation_config, ROOT),
            "sealed_accuracy_cohort_reused": False,
        },
        "lanes": lanes,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--python",
        action="store_true",
        help="inspect the installed rdkit Python distribution",
    )
    parser.add_argument(
        "--python-wrapper-backend", choices=PYTHON_BACKENDS, default="unknown"
    )
    parser.add_argument("--python-backend-evidence")
    parser.add_argument(
        "--python-artifact",
        type=Path,
        help="exact installed wheel/archive; hashed directly rather than inferred",
    )
    parser.add_argument(
        "--npm-package", type=Path, help="installed @rdkit/rdkit package directory"
    )
    parser.add_argument(
        "--npm-runtime-version",
        help="version returned by executing rdkit.version(); never inferred",
    )
    parser.add_argument(
        "--npm-artifact",
        type=Path,
        help="exact installed npm tarball; hashed directly rather than inferred",
    )
    parser.add_argument("--corpus", type=Path, action="append", default=[])
    parser.add_argument("--operation-config", type=Path, action="append", default=[])
    parser.add_argument(
        "--output", type=Path, help="write JSON to this path; stdout when omitted"
    )
    args = parser.parse_args()
    if not args.python and args.npm_package is None:
        parser.error("at least one of --python or --npm-package is required")
    if args.python_wrapper_backend != "unknown" and not args.python_backend_evidence:
        parser.error(
            "--python-backend-evidence is required for a declared wrapper backend"
        )
    if args.python_artifact is not None and not args.python_artifact.is_file():
        parser.error(f"Python artifact not found: {args.python_artifact}")
    if args.npm_artifact is not None and not args.npm_artifact.is_file():
        parser.error(f"npm artifact not found: {args.npm_artifact}")
    if args.npm_artifact is not None and args.npm_package is None:
        parser.error("--npm-artifact requires --npm-package")
    return args


def main() -> int:
    args = parse_args()
    try:
        record = build_record(args)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"RDKit provenance collection failed: {exc}", file=sys.stderr)
        return 1
    body = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(body, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(body, encoding="utf-8")
        print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
