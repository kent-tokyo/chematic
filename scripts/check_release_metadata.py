#!/usr/bin/env python3
"""Validate checked-in release metadata without requiring third-party packages."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "release-metadata" / "v1.0.6.json"
SCHEMA = ROOT / "docs" / "release-metadata-schema.json"


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def main() -> int:
    errors: list[str] = []
    try:
        document = json.loads(METADATA.read_text(encoding="utf-8"))
        schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"release metadata read failure: {exc}", file=sys.stderr)
        return 1

    require(schema.get("$schema") == "https://json-schema.org/draft/2020-12/schema", "schema is not draft 2020-12", errors)
    require(document.get("schema_version") == 1, "schema_version must be 1", errors)
    require(document.get("product") == "chematic", "product must be chematic", errors)
    release = document.get("release", {})
    version = release.get("version")
    require(version == "1.0.6", "checked-in metadata must describe v1.0.6", errors)
    require(release.get("tag") == f"v{version}", "release tag/version mismatch", errors)
    commit = release.get("commit")
    require(
        commit is None or bool(re.fullmatch(r"[0-9a-f]{40}", commit)),
        "release commit must be a full lowercase SHA-1 or null before tag publication",
        errors,
    )
    if commit is None:
        require(bool(release.get("commit_source")), "null release commit must document its source", errors)
    for key in ("rust", "python", "npm"):
        package = document.get("packages", {}).get(key, {})
        require(package.get("version") == version, f"{key} package version mismatch", errors)
        require(bool(package.get("registry_url")), f"{key} registry URL missing", errors)
    mcp = document.get("mcp", {})
    require(mcp.get("tool_count") == 20, "MCP tool count drifted from the documented registry", errors)
    require(mcp.get("transport") == ["stdio"], "MCP transport must remain stdio", errors)
    require(mcp.get("network_enabled_tools") == ["pubchem_lookup"], "MCP network tool declaration drifted", errors)
    for entry in document.get("benchmarks", {}).get("historical", []):
        require(entry.get("status") == "historical", f"benchmark is not marked historical: {entry.get('path')}", errors)
        require(bool(entry.get("path")), "historical benchmark path missing", errors)
    if errors:
        print("release metadata validation failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"release metadata OK: {METADATA.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
