#!/usr/bin/env python3
"""Bridge independently rerun published v1.0.30 chemistry rows to adjudicated residuals.

The source-candidate classification is reusable only because the *complete*
decompressed 10k row stream, including every RDKit and chematic value, is
byte-identical to the published-wheel rerun. This is not a new oracle or a
claim that the two artifacts are interchangeable on other inputs.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

if __package__:
    from .check_rdkit_rebaseline_evidence import validate_smarts_cell_correspondence
else:
    from check_rdkit_rebaseline_evidence import validate_smarts_cell_correspondence

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation/results"
PUBLISHED = RESULTS / "v1.0.30-published-python-chemistry-rows.jsonl.gz"
PUBLISHED_SUMMARY = RESULTS / "v1.0.30-published-python-chemistry-summary.json"
ADJUDICATED = RESULTS / "rdkit-rebaseline-python-chemistry-rows-v1.0.27-issue634-vs-2026.03.6-2026-09-28.jsonl.gz"
CLASSIFICATION = RESULTS / "rdkit-rebaseline-residual-classification-v1.0.27-issue634-vs-2026.03.6-2026-09-28.json"
PUBLISHED_COMPRESSED_SHA256 = "4c429609828c193d84bde7200e0e8b36eaa91e7c01a8aaf4f6ebe73839412452"
CLASSIFICATION_SHA256 = "7f6b32d614eefb5604df17261e23d8f746fa17f2afb2307686d1d97b329d39ef"
EXPECTED_FAMILIES = {
    "ring_semantics:organometallic": 6,
    "ring_semantics:symmetrized_rings": 194,
}
EXPECTED_CIP = {"lone_pair_center": 1, "phosphorus_oracle_unstable": 4}


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def check() -> dict[str, object]:
    published_compressed = PUBLISHED.read_bytes()
    if digest(published_compressed) != PUBLISHED_COMPRESSED_SHA256:
        raise ValueError("published-wheel raw-row archive changed")
    published_rows = gzip.decompress(published_compressed)
    adjudicated_compressed = ADJUDICATED.read_bytes()
    if published_rows != gzip.decompress(adjudicated_compressed):
        raise ValueError("published-wheel rows are not byte-identical to adjudicated rows")

    summary = json.loads(PUBLISHED_SUMMARY.read_text(encoding="utf-8"))
    classification_bytes = CLASSIFICATION.read_bytes()
    if digest(classification_bytes) != CLASSIFICATION_SHA256:
        raise ValueError("source residual classification changed")
    classification = json.loads(classification_bytes)
    if (summary["chematic_version"] != "1.0.30"
            or summary["rdkit_version"] != "2026.03.6"
            or summary["rows"]["sha256"] != digest(published_rows)
            or classification["rdkit_version"] != summary["rdkit_version"]
            or classification["rows_file"]["sha256"] != digest(adjudicated_compressed)
            or classification["rows"] != 10_000
            or classification["smarts_cells"] != 310_000):
        raise ValueError("published summary, oracle, and adjudication inputs disagree")

    rows = [json.loads(line) for line in published_rows.splitlines()]
    if len(rows) != 10_000 or [row["input_index"] for row in rows] != list(range(10_000)):
        raise ValueError("published-row accounting or order changed")
    unresolved = Counter(
        difference["operation"]
        for row in rows
        for difference in row["differences"]
        if difference["classification"] == "unresolved"
    )
    if (unresolved != {"smarts": 64, "cip": 5}
            or sum(unresolved.values()) != summary["counts"]["difference_class_unresolved"]):
        raise ValueError("unresolved row families changed")

    smarts = classification["smarts"]
    errors: list[str] = []
    validate_smarts_cell_correspondence(rows, smarts, errors)
    if (errors or smarts["differing_rows"] != 64 or smarts["differing_cells"] != 200
            or smarts["unclassified_cells"] != 0 or smarts["cells_by_family"] != EXPECTED_FAMILIES
            or summary["counts"]["smarts_differences"] != 200):
        raise ValueError(f"published SMARTS classification mismatch: {errors[:3]}")

    cip = classification["cip"]
    cip_indices = {row["input_index"] for row in rows if row["cip"]["exact"] is not True}
    if (cip_indices != {item["input_index"] for item in cip["rows"]}
            or cip["differing_rows"] != 5 or cip["unclassified_labels"] != 0
            or cip["differing_labels_by_family"] != EXPECTED_CIP
            or summary["counts"]["cip_difference"] != 5):
        raise ValueError("published CIP classification mismatch")
    for item in cip["rows"]:
        row = rows[item["input_index"]]
        for label in item["labels"]:
            expected_reason = {"lone_pair_center": "lone_pair_center",
                               "phosphorus_oracle_unstable": "oracle_unstable"}[label["family"]]
            if (label["kind"] != "atom" or label["chematic"] is not None
                    or label["chematic_unresolved"] != expected_reason
                    or row["cip"]["chematic_unresolved"].get(str(label["atom"])) != expected_reason
                    or row["cip"]["rdkit_atoms"].get(str(label["atom"])) != label["rdkit"]):
                raise ValueError(f"published CIP row {item['input_index']} is not the classified typed abstention")

    return {"published_rows": len(rows), "unresolved_reclassified": sum(unresolved.values()),
            "smarts_rows": 64, "smarts_cells": smarts["differing_cells"],
            "smarts_families": smarts["cells_by_family"], "cip_typed_abstentions": EXPECTED_CIP,
            "raw_rows_sha256": digest(published_rows)}


if __name__ == "__main__":
    print(json.dumps(check(), sort_keys=True))
