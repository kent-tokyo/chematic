#!/usr/bin/env python3
"""Validate the bounded independent RDKit MMFF94 availability oracle."""

from __future__ import annotations

import json
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "validation" / "results" / f"mmff94-rdkit-availability-oracle-v{workspace_version(ROOT)}.json"


def fail(message: str) -> None:
    raise SystemExit(f"MMFF94 RDKit availability oracle invalid: {message}")


def main() -> int:
    try:
        report = json.loads(REPORT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"cannot read report: {error}")
    if report.get("schema_version") != 1:
        fail("schema_version must be 1")
    if report.get("target_version") != workspace_version(ROOT):
        fail("target version is stale")
    if report.get("result") != "pass":
        fail("result must be pass")
    corpus = report.get("corpus")
    if not isinstance(corpus, dict) or corpus.get("molecules") != 265:
        fail("corpus must contain 265 molecules")
    if corpus.get("tiers") != {"A": 65, "B": 200}:
        fail("unexpected tier counts")
    for key in ("manifest_a", "manifest_b"):
        path = ROOT / corpus.get(key, "")
        if not path.is_file():
            fail(f"missing corpus manifest: {path}")
    results = report.get("results")
    expected = {
        "parse_ok": 265,
        "embedding_ok": 265,
        "mmff_properties_available": 265,
        "force_field_constructed": 265,
        "finite_energy": 265,
        "minimize_return_code_0": 142,
        "minimize_return_code_1": 123,
    }
    if results != expected:
        fail(f"result counts changed: {results!r}")
    if not isinstance(report.get("boundary"), str) or not report["boundary"].strip():
        fail("boundary disclosure is missing")
    print("MMFF94 RDKit availability oracle OK: 265/265 finite-energy constructions; minimize split 142/123")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
