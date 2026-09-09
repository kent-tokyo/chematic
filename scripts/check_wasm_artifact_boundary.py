#!/usr/bin/env python3
"""Validate the intentional source-vs-checked-in Node/WASM boundary."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "crates" / "chematic-wasm" / "src"
ARTIFACT = ROOT / "crates" / "chematic-wasm" / "pkg-node" / "chematic_wasm.js"


def main() -> int:
    declaration = re.compile(
        r"#\[wasm_bindgen\](?:\s*#\[[^\]]+\])*\s*pub fn\s+(\w+)", re.MULTILINE
    )
    exported = re.compile(r"^exports\.([a-z_]\w*)\s*=", re.MULTILINE)
    source = "\n".join(p.read_text(encoding="utf-8") for p in sorted(SOURCE.rglob("*.rs")))
    artifact = ARTIFACT.read_text(encoding="utf-8")
    source_names = sorted(set(declaration.findall(source)) - {"start"})
    artifact_names = sorted(set(exported.findall(artifact)) - {"start"})
    missing = sorted(set(source_names) - set(artifact_names))
    stale = sorted(set(artifact_names) - set(source_names))
    report_path = ROOT / "validation" / "results" / f"wasm-node-artifact-rebuild-v{workspace_version(ROOT)}.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    recorded = report.get("existing_node_artifact_surface", {})
    errors: list[str] = []
    if recorded.get("source_export_count") != len(source_names):
        errors.append("recorded source export count is stale")
    if recorded.get("artifact_export_count") != len(artifact_names):
        errors.append("recorded artifact export count is stale")
    if recorded.get("missing_from_artifact") != missing:
        errors.append("recorded missing export list is stale")
    if recorded.get("stale_in_artifact") != stale:
        errors.append("recorded stale export list is stale")
    if not missing and not stale:
        errors.append("boundary unexpectedly claims no artifact drift")
    if errors:
        print("WASM artifact boundary failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(
        "WASM artifact boundary OK: "
        f"source={len(source_names)}, artifact={len(artifact_names)}, missing={len(missing)}, stale={len(stale)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
