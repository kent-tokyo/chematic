#!/usr/bin/env python3
"""Compare new product-template atom maps from a pinned Rust build with RDKit.

This supplement preserves the original 83-case diagnostic byte-for-byte.
Product atom maps are diagnostic metadata, not a claim about reaction yields.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from rdkit import rdBase

if __package__:
    from .reaction_atom_provenance_gate import rdkit_sets, rust_sets
    from .reaction_template_map_gate import rdkit_map_sets, rust_map_sets
else:
    from reaction_atom_provenance_gate import rdkit_sets, rust_sets
    from reaction_template_map_gate import rdkit_map_sets, rust_map_sets

ROOT = Path(__file__).resolve().parents[1]


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "validation/reaction_product_parity_cases.json")
    parser.add_argument("--strata", type=Path, default=ROOT / "validation/reaction_product_parity_strata_v2.json")
    parser.add_argument("--supplement", type=Path, default=ROOT / "validation/reaction_product_new_maps_v3.json")
    parser.add_argument("--baseline-rows", type=Path, default=ROOT / "validation/results/v1.0.30-published-rust-reaction-83-template-map-rows.json")
    parser.add_argument("--rust-rows", type=Path, required=True)
    parser.add_argument("--rust-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rdkit", default="2026.03.6")
    parser.add_argument("--expected-crate", default="chematic 1.0.30 from crates.io")
    parser.add_argument("--expected-version", default="1.0.30")
    parser.add_argument("--expected-profile", default=None)
    args = parser.parse_args()
    if rdBase.rdkitVersion != args.expected_rdkit:
        parser.error(f"RDKit {rdBase.rdkitVersion} != {args.expected_rdkit}")
    supplement_bytes = args.supplement.read_bytes()
    supplement = json.loads(supplement_bytes)
    base_hash = sha256(args.base.read_bytes())
    strata_hash = sha256(args.strata.read_bytes())
    if supplement["schema_version"] != 3 or (supplement["base_fixture_sha256"], supplement["strata_fixture_sha256"]) != (base_hash, strata_hash):
        raise ValueError("supplement fixture does not pin the historical input set")
    cases = supplement["cases"]
    if not cases or len({case["id"] for case in cases}) != len(cases):
        raise ValueError("supplement case IDs must be nonempty and unique")
    raw = args.rust_rows.read_bytes()
    rust_rows = json.loads(raw)
    baseline = json.loads(args.baseline_rows.read_bytes())
    summary = json.loads(args.rust_summary.read_bytes())
    if (len(baseline) != 83 or rust_rows[:83] != baseline
            or len(rust_rows) != 83 + len(cases)
            or [row["id"] for row in rust_rows[83:]] != [case["id"] for case in cases]
            or summary["schema"] != "published-rust-reaction-template-maps/v1"
            or summary["crate"] != args.expected_crate
            or summary.get("compatibility_profile") != args.expected_profile
            or summary["input_count"] != len(rust_rows)
            or summary["rows_sha256"] != sha256(raw)
            or summary["base_cases_sha256"] != base_hash
            or summary["strata_sha256"] != strata_hash
            or summary["supplement_sha256"] != sha256(supplement_bytes)):
        raise ValueError("Rust rows or fixture provenance disagree")

    results = []
    for case, row in zip(cases, rust_rows[83:], strict=True):
        if row["status"] != "products":
            outcome = "typed_refusal" if row["status"] in {"typed_refusal", "diagnosed_valence_refusal"} else "unexpected_status"
            results.append({"id": case["id"], "strata": case["strata"], "outcome": outcome, "rust_status": row["status"]})
            continue
        born = [(source, label) for product_set in row["sets"] for product in product_set
                for source, label in zip(product["atom_sources"], product["template_map_numbers"], strict=True)
                if source is None]
        if not born or ("mapped" in case["strata"]) != any(label is not None for _, label in born):
            raise ValueError(f"{case['id']}: born atom map diagnostic does not match fixture stratum")
        oracle_graph, oracle_origin, raw_count = rdkit_sets(case)
        oracle_maps, map_raw_count = rdkit_map_sets(case)
        if raw_count != map_raw_count:
            raise ValueError(f"{case['id']}: RDKit oracle counts disagree")
        candidate_graph, candidate_origin = rust_sets(case, {"sets_with_atom_sources": row["sets"]})
        candidate_maps = rust_map_sets(row)
        if candidate_graph != oracle_graph:
            outcome = "graph_mismatch"
        elif candidate_origin != oracle_origin:
            outcome = "origin_mismatch"
        elif candidate_maps != oracle_maps:
            outcome = "map_label_mismatch"
        else:
            outcome = "graph_origin_map_match"
        results.append({"id": case["id"], "strata": case["strata"], "outcome": outcome,
                        "rdkit_graph": oracle_graph, "rust_graph": candidate_graph,
                        "rdkit_origins": oracle_origin, "rust_origins": candidate_origin,
                        "rdkit_maps": oracle_maps, "rust_maps": candidate_maps,
                        "rdkit_raw_product_sets": raw_count, "rust_raw_product_sets": len(row["sets"])})
    report = {"schema": "published-rust-reaction-new-product-maps/v1",
              "crate_version": args.expected_version, "rdkit_version": rdBase.rdkitVersion,
              "fixtures": {"base_sha256": base_hash, "strata_sha256": strata_hash,
                           "supplement_sha256": sha256(supplement_bytes)},
              "rust_rows_sha256": sha256(raw),
              "rust_summary_sha256": sha256(args.rust_summary.read_bytes()),
              "baseline_rows_sha256": sha256(args.baseline_rows.read_bytes()),
              "accounting": {"input": len(cases), "outcomes": dict(sorted(Counter(item["outcome"] for item in results).items()))},
              "limits": ["Exposed supplement, not broad SMIRKS parity",
                         "No atom maps added to ordinary product molecules",
                         "No reaction yield or selectivity claim"],
              "rows": results}
    if args.expected_profile is not None:
        report["crate"] = args.expected_crate
        report["compatibility_profile"] = args.expected_profile
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"new product maps: {len(results)} inputs, {report['accounting']['outcomes']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
