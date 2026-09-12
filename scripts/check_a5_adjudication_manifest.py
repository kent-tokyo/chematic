#!/usr/bin/env python3
"""Fail-closed structural validator for the A5 external-review packet.

This checks that a gold/blind packet is reviewable and that blind cases do not
leak an expected answer. It deliberately does not score chematic or RDKit and
cannot replace the required non-maintainer adjudication.
"""

from __future__ import annotations

import json
import hashlib
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "validation/manifests/rdkit_accuracy_adjudication_v1.json"


def fail(message: str) -> None:
    raise SystemExit(f"A5 adjudication manifest invalid: {message}")


def referenced_smiles(input_ref: str) -> str | None:
    path_text, _, case_id = input_ref.partition("#")
    if not path_text or not case_id:
        return None
    path = ROOT / path_text
    if not path.is_file():
        return None
    for line in path.read_text(encoding="utf-8").splitlines():
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if value.get("id") == case_id and isinstance(value.get("smiles"), str):
            return value["smiles"]
    return None


def main() -> int:
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1:
        fail("schema_version must be 1")
    if data.get("status") != "prepared_external_review_required":
        fail("manifest must remain explicitly external-review pending")
    if data.get("decision_rules", {}).get("reviewer_must_be_non_maintainer") is not True:
        fail("reviewer independence rule is missing")
    gold = data.get("gold_cases")
    blind = data.get("blind_cases")
    if not isinstance(gold, list) or not gold:
        fail("gold_cases must be non-empty")
    if not isinstance(blind, list) or not blind:
        fail("blind_cases must be non-empty")
    ids = [case.get("id") for case in [*gold, *blind]]
    if any(not isinstance(case_id, str) or not case_id for case_id in ids):
        fail("every case needs a non-empty id")
    if len(set(ids)) != len(ids):
        fail("case ids must be unique")
    for case in gold:
        for key in ("input_ref", "atom_mapping", "expected", "rationale", "unresolved_reason"):
            if key not in case:
                fail(f"gold case {case.get('id')!r} lacks {key}")
    for case in blind:
        digest = case.get("input_digest")
        if (
            not isinstance(digest, str)
            or re.fullmatch(r"sha256:[0-9a-f]{64}", digest) is None
        ):
            fail(f"blind case {case.get('id')!r} lacks an input digest")
        input_ref = case.get("input_ref")
        if not isinstance(input_ref, str):
            fail(f"blind case {case.get('id')!r} lacks an input reference")
        smiles = referenced_smiles(input_ref)
        if smiles is None:
            fail(f"blind case {case.get('id')!r} references no resolvable SMILES")
        actual_digest = "sha256:" + hashlib.sha256(smiles.encode("utf-8")).hexdigest()
        if digest != actual_digest:
            fail(f"blind case {case.get('id')!r} input digest does not match its reference")
        if case.get("expected") is not None or case.get("rationale") is not None:
            fail(f"blind case {case.get('id')!r} leaks adjudication labels")
    print(f"A5 adjudication packet OK: {len(gold)} gold + {len(blind)} blind; external review required")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
