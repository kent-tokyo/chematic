"""Read the current workspace version for generated benchmark metadata."""

from __future__ import annotations

import re
import os
from pathlib import Path


def workspace_version(root: Path) -> str:
    evidence_version = os.environ.get("SCHEMATIC_BENCHMARK_VERSION")
    if evidence_version:
        if not re.fullmatch(r"\d+\.\d+\.\d+", evidence_version):
            raise RuntimeError("SCHEMATIC_BENCHMARK_VERSION must be semantic version")
        return evidence_version
    cargo = (root / "Cargo.toml").read_text(encoding="utf-8")
    match = re.search(r"(?ms)^\[workspace\.package\]\s*$.*?^version\s*=\s*\"([^\"]+)\"\s*$", cargo)
    if not match:
        raise RuntimeError("Cargo.toml has no [workspace.package] version")
    return match.group(1)
