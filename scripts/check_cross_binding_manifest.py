#!/usr/bin/env python3
"""Validate the shared cross-binding operation inventory.

This gate checks the manifest's bookkeeping only.  It does not claim that the
listed operations cover every public API, nor does it execute a binding test.
Each listed operation must point at an existing contract section and at least
one checked-in local test/source anchor; the four binding names are explicit
so a future expansion cannot silently reduce the declared surface.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "validation" / "cross_binding_contract.json"
REQUIRED_BINDINGS = {"rust", "python", "node", "wasm"}
SOURCE_ONLY_BINDINGS = {"rust", "python", "wasm"}
CONTRACT_SECTIONS = {
    "fixtures",
    "mol_v2000_contract",
    "mol_v3000_contract",
    "xyz_contract",
    "batch_canonicalization_contract",
    "smiles_validity_contract",
    "inchi_contract",
    "fragment_parent_contract",
    "charge_parent_contract",
    "isotope_parent_contract",
    "stereo_parent_contract",
    "canonical_tautomer_contract",
    "neutralize_charges_contract",
    "largest_fragment_contract",
    "tautomer_parent_contract",
    "super_parent_contract",
    "hydrogen_contract",
    "super_parent_report_contract",
    "extxyz_contract",
    "pdb_contract",
    "qcschema_contract",
    "qcschema_atomic_input_contract",
    "qcschema_atomic_result_contract",
    "orca_input_contract",
    "orca_output_contract",
    "mol2_contract",
    "cml_contract",
    "cdxml_contract",
    "mmcif_contract",
    "moljson_contract",
    "xyz_batch_contract",
    "rxn_document_contract",
    "reaction_application_contract",
    "reaction_balance_contract",
    "reaction_center_contract",
    "semantic_expansion_contract",
    "descriptor_contract",
    "standardization_contract",
    "fingerprint_contract",
    "fingerprint_detail_contract",
    "adversarial",
}


def main() -> int:
    try:
        document = json.loads(MANIFEST.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"cross-binding manifest read failure: {exc}", file=sys.stderr)
        return 1

    inventory = document.get("operation_manifest")
    operations = inventory.get("operations") if isinstance(inventory, dict) else None
    errors: list[str] = []
    if not isinstance(inventory, dict) or inventory.get("schema_version") != 1:
        errors.append("operation_manifest schema_version must be 1")
    if inventory.get("scope") != "currently_shared_stable_surface":
        errors.append("operation_manifest scope must identify the current shared surface")
    if not isinstance(operations, list) or not operations:
        errors.append("operation_manifest.operations must be a non-empty list")
        operations = []

    ids: set[str] = set()
    referenced_sections: set[str] = set()
    for index, operation in enumerate(operations):
        prefix = f"operation {index}"
        if not isinstance(operation, dict):
            errors.append(f"{prefix} must be an object")
            continue
        identifier = operation.get("id")
        if not isinstance(identifier, str) or not identifier:
            errors.append(f"{prefix} has no non-empty id")
        elif identifier in ids:
            errors.append(f"duplicate operation id: {identifier}")
        else:
            ids.add(identifier)
        section = operation.get("contract")
        if section not in CONTRACT_SECTIONS or section not in document:
            errors.append(f"{prefix} references unknown contract section: {section!r}")
        else:
            referenced_sections.add(section)
        bindings = operation.get("bindings")
        if set(bindings or []) != REQUIRED_BINDINGS:
            errors.append(f"{prefix} must declare exactly rust/python/node/wasm bindings")
        anchors = operation.get("test_anchors")
        if not isinstance(anchors, list) or not anchors:
            errors.append(f"{prefix} must have at least one test anchor")
        else:
            for anchor in anchors:
                if not isinstance(anchor, str) or not (ROOT / anchor).is_file():
                    errors.append(f"{prefix} has missing test anchor: {anchor!r}")

    missing_sections = sorted(CONTRACT_SECTIONS - referenced_sections)
    errors.extend(f"contract section is not in operation manifest: {section}" for section in missing_sections)
    source_only = document.get("source_only_operations")
    if not isinstance(source_only, list):
        errors.append("source_only_operations must be a list")
        source_only = []
    source_ids: set[str] = set()
    for index, operation in enumerate(source_only):
        prefix = f"source-only operation {index}"
        if not isinstance(operation, dict):
            errors.append(f"{prefix} must be an object")
            continue
        identifier = operation.get("id")
        if not isinstance(identifier, str) or not identifier or identifier in source_ids:
            errors.append(f"{prefix} has a missing or duplicate id")
        else:
            source_ids.add(identifier)
        section = operation.get("contract")
        if section not in document or section in CONTRACT_SECTIONS:
            errors.append(f"{prefix} must reference a source-only contract section")
        bindings = set(operation.get("bindings") or [])
        if bindings != SOURCE_ONLY_BINDINGS:
            errors.append(f"{prefix} must declare exactly rust/python/wasm bindings")
        if operation.get("missing_binding") not in {
            "node-generated-artifact",
            "node-generated-artifact-and-export",
        }:
            errors.append(f"{prefix} must record the missing Node generated artifact")
        for anchor in operation.get("test_anchors") or []:
            if not isinstance(anchor, str) or not (ROOT / anchor).is_file():
                errors.append(f"{prefix} has missing test anchor: {anchor!r}")
    if errors:
        print("Cross-binding manifest failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(
        "Cross-binding operation manifest OK: "
        f"{len(operations)} four-binding operations, "
        f"{len(source_only)} source-only operations, 4 bindings declared"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
