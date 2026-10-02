#!/usr/bin/env python3
"""Stratified, fail-closed product compatibility on published Python wheels.

All 57 historical fixtures remain in the denominator. Product graphs and
multiplicity are compared; map/provenance evidence is explicitly unavailable
in the published Python run_smirks result and is never implied by a match.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

from reaction_product_parity_gate import rdkit_products

ROOT = Path(__file__).resolve().parents[1]
RDLogger.DisableLog("rdApp.*")


def load_cases(base_path: Path, strata_path: Path) -> tuple[list[dict], dict]:
    base = json.loads(base_path.read_text(encoding="utf-8"))
    strata = json.loads(strata_path.read_text(encoding="utf-8"))
    if strata.get("schema_version") != 2 or strata.get("base_fixture_sha256") != hashlib.sha256(base_path.read_bytes()).hexdigest():
        raise ValueError("strata manifest does not pin the historical base fixture")
    cases = [{**case, "strata": ["legacy"]} for case in base] + strata["cases"]
    ids = [case["id"] for case in cases]
    if len(base) != 57 or len(ids) != len(set(ids)):
        raise ValueError("legacy fixture count or case IDs changed")
    if any(not case.get("strata") or not isinstance(case.get("reactants"), list) for case in cases):
        raise ValueError("every case needs strata and reactants")
    return cases, {"base_sha256": hashlib.sha256(base_path.read_bytes()).hexdigest(),
                   "strata_sha256": hashlib.sha256(strata_path.read_bytes()).hexdigest()}


def candidate_products(chematic, case: dict, *, checked_rdkit_compat: bool = False) -> dict:
    try:
        reactants = [chematic.from_smiles(value) for value in case["reactants"]]
    except ValueError as exc:
        return {"status": "typed_refusal" if checked_rdkit_compat else "untyped_refusal",
                "reason": "reactant_parse" if checked_rdkit_compat else None,
                "stage": "reactant_parse", "error_type": type(exc).__name__, "detail": str(exc)}
    try:
        if checked_rdkit_compat:
            checked = chematic.run_smirks_checked(case["smirks"], reactants, rdkit_compat=True)
            if checked["status"] in {"typed_refusal", "typed_unsupported", "partial_products"}:
                return {key: value for key, value in checked.items() if key != "products"}
            raw = checked["products"]
        else:
            raw = chematic.run_smirks(case["smirks"], reactants)
    except ValueError as exc:
        # The published Python binding maps all Rust TransformError variants
        # to ValueError. Message inspection cannot promote this to typed refusal.
        return {"status": "untyped_refusal", "stage": "transform", "error_type": type(exc).__name__, "detail": str(exc)}
    values = set()
    for product_set in raw:
        normalized = []
        for product in product_set:
            mol = Chem.MolFromSmiles(product.smiles)
            if mol is None:
                return {"status": "invalid_product", "detail": product.smiles}
            normalized.append(Chem.MolToSmiles(mol, canonical=True))
        values.add(tuple(sorted(normalized)))
    result = {"status": "products", "sets": [list(items) for items in sorted(values)], "raw_product_sets": len(raw)}
    if checked_rdkit_compat:
        result["diagnostics"] = {key: value for key, value in checked.items() if key != "products"}
    return result


def classify_case(chematic, case: dict, *, checked_rdkit_compat: bool = False) -> dict:
    candidate = candidate_products(chematic, case, checked_rdkit_compat=checked_rdkit_compat)
    try:
        oracle_sets, oracle_raw_count = rdkit_products(case["smirks"], case["reactants"], case["id"])
        oracle = {"status": "products", "sets": oracle_sets, "raw_product_sets": oracle_raw_count}
    except ValueError as exc:
        oracle = {"status": "invalid_input_or_product", "detail": str(exc)}
    if oracle["status"] != "products":
        outcome = "joint_invalid_input" if candidate["status"] in {"untyped_refusal", "typed_refusal"} else "oracle_invalid"
    elif candidate["status"] == "products":
        outcome = "semantic_match" if candidate["sets"] == oracle["sets"] else "wrong_confident"
    else:
        outcome = candidate["status"]
    return {**case, "outcome": outcome, "rdkit": oracle, "chematic": candidate,
            "map_provenance": "not_exposed_by_published_python_run_smirks"}


def classify_rust_row(rust: dict, oracle: dict) -> dict:
    status = rust["status"]
    if status == "products":
        normalized = set()
        for product_set in rust["sets"]:
            products = []
            for spelling in product_set:
                mol = Chem.MolFromSmiles(spelling)
                if mol is None:
                    return {"outcome": "invalid_product", "raw": rust}
                products.append(Chem.MolToSmiles(mol, canonical=True))
            normalized.add(tuple(sorted(products)))
        sets = [list(items) for items in sorted(normalized)]
        outcome = ("semantic_match" if sets == oracle["sets"] else "wrong_confident") if oracle["status"] == "products" else "oracle_invalid"
        return {"outcome": outcome, "sets": sets, "raw": rust}
    if status in {"typed_refusal", "diagnosed_valence_refusal"}:
        return {"outcome": status if oracle["status"] == "products" else "joint_invalid_input", "raw": rust}
    raise ValueError(f"unknown Rust status: {status}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "validation/reaction_product_parity_cases.json")
    parser.add_argument("--strata", type=Path, default=ROOT / "validation/reaction_product_parity_strata_v2.json")
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rdkit", required=True)
    parser.add_argument("--rust-rows", type=Path, help="optional exact-version crates.io reaction rows")
    parser.add_argument("--checked-rdkit-compat", action="store_true",
                        help="source-candidate typed RDKit-compat profile; not published v1.0.30 evidence")
    args = parser.parse_args()
    if rdBase.rdkitVersion != args.expected_rdkit:
        parser.error(f"RDKit {rdBase.rdkitVersion} != {args.expected_rdkit}")
    import chematic

    cases, hashes = load_cases(args.base, args.strata)
    rows = [classify_case(chematic, case, checked_rdkit_compat=args.checked_rdkit_compat) for case in cases]
    rust_counts = None
    if args.rust_rows:
        rust_rows = json.loads(args.rust_rows.read_text(encoding="utf-8"))
        if [row["id"] for row in rust_rows] != [row["id"] for row in rows]:
            raise ValueError("Rust rows and pinned fixtures differ in count or order")
        for row, rust_row in zip(rows, rust_rows, strict=True):
            row["rust"] = classify_rust_row(rust_row, row["rdkit"])
        rust_counts = dict(sorted(Counter(row["rust"]["outcome"] for row in rows).items()))
    counts = Counter(row["outcome"] for row in rows)
    by_stratum: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        for name in row["strata"]:
            by_stratum[name][row["outcome"]] += 1
    report = {
        "schema": "stratified-reaction-compatibility/v2",
        "profile": "checked-rdkit-compat-source" if args.checked_rdkit_compat else "published-legacy",
        "artifact": {"version": importlib.metadata.version("chematic"),
                     "filename": args.artifact.name,
                     "sha256": hashlib.sha256(args.artifact.read_bytes()).hexdigest()},
        "rdkit_version": rdBase.rdkitVersion,
        "fixtures": hashes,
        "accounting": {"input": len(cases), "outcomes": dict(sorted(counts.items())),
                       "by_stratum": {key: dict(sorted(value.items())) for key, value in sorted(by_stratum.items())},
                       "rust_outcomes": rust_counts},
        "rows": rows,
        "limits": ["No reaction yield or selectivity claim", "Product atom-map provenance was not compared with RDKit",
                   "Python result does not expose product atom provenance",
                   "A checked source-candidate profile is not published v1.0.30 evidence" if args.checked_rdkit_compat
                   else "ValueError from published binding is not a typed refusal category"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"reaction v2: {len(rows)} inputs, {dict(counts)}")
    return 1 if counts["wrong_confident"] or counts["invalid_product"] or counts["oracle_invalid"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
