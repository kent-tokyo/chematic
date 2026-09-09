#!/usr/bin/env python3
"""Audit the Rust/Python/WASM binding surface without building or networking.

The result is intentionally an inventory, not a parity claim.  It gives the
roadmap's shared-manifest work a deterministic starting point by listing
exported function names and the Python functions registered on the module.
"""

from __future__ import annotations

import json
import hashlib
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "validation" / "results" / "binding_surface_inventory-v1.0.10.json"
WASM = ROOT / "crates" / "chematic-wasm" / "src"
PY = ROOT / "crates" / "chematic-py" / "src"
DECL_RE = re.compile(
    r"#\[wasm_bindgen\](?:\s*#\[[^\]]+\])*\s*pub fn\s+(\w+)",
    re.MULTILINE,
)
PY_DECL_RE = re.compile(r"#\[pyfunction\](?:\s*#\[[^\]]+\])*\s*fn\s+(\w+)", re.MULTILINE)
PY_REGISTER_RE = re.compile(r"wrap_pyfunction!\((\w+)", re.MULTILINE)


def source_text(root: Path) -> str:
    return "\n".join(path.read_text(encoding="utf-8") for path in sorted(root.rglob("*.rs")))


def main() -> int:
    wasm_names = sorted(set(DECL_RE.findall(source_text(WASM))) - {"start"})
    python_declared = sorted(set(PY_DECL_RE.findall(source_text(PY))))
    python_registered = sorted(set(PY_REGISTER_RE.findall(source_text(PY))))
    missing_registration = sorted(set(python_declared) - set(python_registered))
    surface_payload = json.dumps(
        {
            "python_pyfunctions": python_declared,
            "python_registered": python_registered,
            "wasm_exports": wasm_names,
        },
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    surface_sha256 = hashlib.sha256(surface_payload).hexdigest()
    result = {
        "schema_version": 1,
        "scope": "static_binding_surface_inventory",
        "wasm_export_count": len(wasm_names),
        "wasm_exports": wasm_names,
        "python_pyfunction_count": len(python_declared),
        "python_pyfunctions": python_declared,
        "python_registered_count": len(python_registered),
        "python_registered": python_registered,
        "python_declared_without_registration": missing_registration,
        "surface_sha256": surface_sha256,
    }
    try:
        recorded = json.loads(REPORT.read_text(encoding="utf-8"))["results"]
    except (OSError, json.JSONDecodeError, KeyError) as exc:
        print(f"binding surface report read failure: {exc}", file=sys.stderr)
        return 1
    for key in ("wasm_export_count", "python_pyfunction_count", "python_registered_count"):
        if recorded.get(key) != result[key]:
            print(f"binding surface report drift for {key}: recorded={recorded.get(key)} current={result[key]}", file=sys.stderr)
            return 1
    if recorded.get("python_declared_without_registration") != missing_registration:
        print("binding surface report drift for unregistered Python functions", file=sys.stderr)
        return 1
    if recorded.get("surface_sha256") != surface_sha256:
        print(
            "binding surface report drift for exported-name set: "
            f"recorded={recorded.get('surface_sha256')} current={surface_sha256}",
            file=sys.stderr,
        )
        return 1
    if missing_registration:
        print(json.dumps(result, indent=2), file=sys.stderr)
        print("Python #[pyfunction] declarations missing module registration", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
