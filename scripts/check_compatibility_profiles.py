#!/usr/bin/env python3
"""Validate the machine-readable API/profile Compatibility Contract."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROFILES = ROOT / "validation" / "compatibility_profiles.json"
REQUIRED_IDS = {
    "smiles_parse_write", "canonical_identity", "aromaticity", "cip_labels",
    "smarts_substructure", "morgan_ecfp", "mol_sdf_v2000", "mol_sdf_v3000",
}
REQUIRED_KEYS = {"id", "api", "bindings", "support_status", "profile", "comparison", "oracle", "evidence", "coverage"}


def main() -> int:
    try:
        document = json.loads(PROFILES.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"Compatibility profiles invalid: {exc}", file=sys.stderr)
        return 1
    errors: list[str] = []
    if document.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    operations = document.get("operations")
    if not isinstance(operations, list):
        errors.append("operations must be an array")
        operations = []
    ids = set()
    for operation in operations:
        if not isinstance(operation, dict) or set(operation) != REQUIRED_KEYS:
            errors.append("every operation must contain exactly the contract keys")
            continue
        operation_id = operation["id"]
        if not isinstance(operation_id, str) or operation_id in ids:
            errors.append("operation ids must be unique non-empty strings")
        ids.add(operation_id)
        if not isinstance(operation["bindings"], list) or set(operation["bindings"]) != {"rust", "python", "node", "wasm"}:
            errors.append(f"{operation_id}: all four bindings are required")
        if operation["support_status"] not in {"stable", "experimental", "unsupported"}:
            errors.append(f"{operation_id}: invalid support_status")
        if operation["profile"] not in {"native", "rdkit_compatibility", "binding_consistency"}:
            errors.append(f"{operation_id}: invalid profile")
        if operation["coverage"] not in {"not_measured", "partial"} and not str(operation["coverage"]).startswith(tuple(str(n) for n in range(10))):
            errors.append(f"{operation_id}: coverage must be explicit")
        if not isinstance(operation["evidence"], list) or not operation["evidence"]:
            errors.append(f"{operation_id}: evidence is required")
        else:
            for path in operation["evidence"]:
                if not isinstance(path, str) or not (ROOT / path).is_file():
                    errors.append(f"{operation_id}: evidence path is missing: {path}")
    missing = REQUIRED_IDS - ids
    if missing:
        errors.append("required operations missing: " + ", ".join(sorted(missing)))
    if errors:
        print("Compatibility profiles invalid:", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 1
    print(f"Compatibility profiles OK: {len(operations)} API/profile operations")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
