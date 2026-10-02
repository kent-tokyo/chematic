#!/usr/bin/env python3
"""Separate distinct-product mismatches from raw reaction embedding counts.

Consumes the pinned published-v1.0.30 83-case, three-binding RDKit comparison.
Raw embedding counts are diagnostics: the same unique product can arise from
multiple symmetry-related RDKit matches.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "validation/results/v1.0.30-published-reaction-83-all-bindings.json"


def classify(oracle: dict, candidate: dict) -> dict:
    if oracle["status"] != "products":
        return {"outcome": "oracle_invalid_or_unavailable", "candidate_status": candidate["status"]}
    if candidate["status"] != "products":
        return {"outcome": "candidate_refusal_or_invalid", "candidate_status": candidate["status"],
                "rdkit_unique_product_sets": len(oracle["sets"]),
                "rdkit_raw_product_sets": oracle["raw_product_sets"]}
    oracle_unique = oracle["sets"]
    candidate_unique = candidate["sets"]
    oracle_raw = oracle["raw_product_sets"]
    candidate_raw = candidate["raw_product_sets"]
    if not all(isinstance(value, int) and value >= 0 for value in (oracle_raw, candidate_raw)):
        raise ValueError("raw product counts must be nonnegative integers")
    if candidate_unique != oracle_unique:
        outcome = "distinct_product_sets_differ"
    elif candidate_raw == oracle_raw:
        outcome = "distinct_and_raw_counts_match"
    elif candidate_raw < oracle_raw:
        outcome = "distinct_match_raw_count_lower"
    else:
        outcome = "distinct_match_raw_count_higher"
    return {"outcome": outcome, "rdkit_unique_product_sets": len(oracle_unique),
            "candidate_unique_product_sets": len(candidate_unique),
            "rdkit_raw_product_sets": oracle_raw,
            "candidate_raw_product_sets": candidate_raw}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.source.read_bytes()
    source = json.loads(raw)
    if (source["schema"] != "stratified-reaction-compatibility/v2"
            or source["profile"] != "published-legacy"
            or source["artifact"]["version"] != "1.0.30"
            or source["rdkit_version"] != "2026.03.6"
            or len(source["rows"]) != 83):
        raise ValueError("expected the pinned published-v1.0.30 83-case report")
    rows = []
    for row in source["rows"]:
        rust = {**row["rust"]["raw"], "sets": row["rust"].get("sets")}
        npm = {**row["npm"]["raw"], "sets": row["npm"].get("sets")}
        results = {
            "python": classify(row["rdkit"], row["chematic"]),
            "rust": classify(row["rdkit"], rust),
            "npm": classify(row["rdkit"], npm),
        }
        # The existing graph oracle adjudication must agree with this gate.
        for binding, prior in (("python", row["outcome"]),
                               ("rust", row["rust"]["outcome"]),
                               ("npm", row["npm"]["outcome"])):
            outcome = results[binding]["outcome"]
            if prior == "semantic_match" and outcome not in {
                    "distinct_and_raw_counts_match", "distinct_match_raw_count_lower",
                    "distinct_match_raw_count_higher"}:
                raise ValueError(f"{row['id']}: {binding} contradicts the existing semantic match")
            if prior == "wrong_confident" and outcome != "distinct_product_sets_differ":
                raise ValueError(f"{row['id']}: {binding} contradicts the existing confident mismatch")
        rows.append({"id": row["id"], "strata": row["strata"], "results": results})
    counts = {binding: dict(sorted(Counter(row["results"][binding]["outcome"]
                                   for row in rows).items()))
              for binding in ("python", "rust", "npm")}
    report = {"schema": "published-reaction-distinct-vs-raw-multiplicity/v1",
              "version": "1.0.30", "rdkit_version": "2026.03.6",
              "source_sha256": hashlib.sha256(raw).hexdigest(),
              "fixtures": source["fixtures"],
              "accounting": {"input": len(rows), "by_binding": counts},
              "limits": ["Distinct-product equality is required for semantic compatibility",
                         "Raw embedding counts are diagnostic; unequal counts alone do not imply a missing distinct product",
                         "No reaction yield or selectivity claim"],
              "rows": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"reaction multiplicity: {len(rows)} cases, {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
