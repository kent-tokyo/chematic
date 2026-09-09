#!/usr/bin/env python3
"""Compare declared WASM exports with the checked-in Node glue surface."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "crates" / "chematic-wasm" / "src"
ARTIFACT = ROOT / "crates" / "chematic-wasm" / "pkg-node" / "chematic_wasm.js"
DECLARED = re.compile(
    r"#\[wasm_bindgen\](?:\s*#\[[^\]]+\])*\s*pub fn\s+(\w+)",
    re.MULTILINE,
)
EXPORTED = re.compile(r"^exports\.([a-z_]\w*)\s*=", re.MULTILINE)


def main() -> int:
    source = "\n".join(p.read_text(encoding="utf-8") for p in sorted(SOURCE.rglob("*.rs")))
    artifact = ARTIFACT.read_text(encoding="utf-8")
    declared = sorted(set(DECLARED.findall(source)) - {"start"})
    exported = sorted(set(EXPORTED.findall(artifact)) - {"start"})
    missing = sorted(set(declared) - set(exported))
    stale = sorted(set(exported) - set(declared))
    result = {
        "schema_version": 1,
        "source_export_count": len(declared),
        "artifact_export_count": len(exported),
        "missing_from_artifact": missing,
        "stale_in_artifact": stale,
        "status": "match" if not missing and not stale else "drift",
    }
    print(json.dumps(result, indent=2))
    return 0 if not missing and not stale else 1


if __name__ == "__main__":
    raise SystemExit(main())
