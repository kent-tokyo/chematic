#!/usr/bin/env python3
"""Compare published Rust reaction template-map identity with pinned RDKit.

RDKit's product `old_mapno` property is compared with map numbers recovered
from each accepted CheMatic ReactionMatch and TracedProduct atom source. The
ordinary CheMatic product molecule intentionally clears atom-map annotations.
The exposed 83 cases cannot prove broad SMIRKS parity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

if __package__:
    from .reaction_atom_provenance_gate import rdkit_sets, rust_sets
    from .run_reaction_compatibility_v2 import load_cases
else:
    from reaction_atom_provenance_gate import rdkit_sets, rust_sets
    from run_reaction_compatibility_v2 import load_cases

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]


def map_labelled_smiles(mol: Chem.Mol, labels: list[int | None]) -> str:
    if mol.GetNumAtoms() != len(labels):
        raise ValueError("product atom count and template-map count differ")
    labelled = Chem.Mol(mol)
    for atom, label in zip(labelled.GetAtoms(), labels, strict=True):
        if label is not None and not 0 < label <= 65535:
            raise ValueError("template-map number is outside the Rust u16 domain")
        atom.SetAtomMapNum(label or 0)
    return Chem.MolToSmiles(labelled, canonical=True, isomericSmiles=True)


def rdkit_map_sets(case: dict) -> tuple[list[list[str]], int]:
    reaction = AllChem.ReactionFromSmarts(case["smirks"])
    if reaction is None:
        raise ValueError("RDKit cannot parse SMIRKS")
    reactants = [Chem.MolFromSmiles(text) for text in case["reactants"]]
    if any(mol is None for mol in reactants):
        raise ValueError("RDKit cannot parse a reactant")
    raw = reaction.RunReactants(tuple(reactants))
    map_sets = set()
    for product_set in raw:
        labelled = []
        for product in product_set:
            checked = Chem.Mol(product)
            Chem.SanitizeMol(checked)
            labels = []
            for atom in checked.GetAtoms():
                if atom.HasProp("react_idx") != atom.HasProp("react_atom_idx"):
                    raise ValueError("RDKit product atom has incomplete reactant origin")
                labels.append(atom.GetIntProp("old_mapno") if atom.HasProp("old_mapno") else None)
            labelled.append(map_labelled_smiles(checked, labels))
        map_sets.add(tuple(sorted(labelled)))
    return [list(values) for values in sorted(map_sets)], len(raw)


def rust_map_sets(row: dict) -> list[list[str]]:
    result = set()
    for product_set in row["sets"]:
        labelled = []
        for product in product_set:
            mol = Chem.MolFromSmiles(product["smiles"])
            if mol is None:
                raise ValueError("published Rust product is not RDKit-readable")
            if len(product["atom_sources"]) != mol.GetNumAtoms():
                raise ValueError("published Rust atom origin count changed")
            labelled.append(map_labelled_smiles(mol, product["template_map_numbers"]))
        result.add(tuple(sorted(labelled)))
    return [list(values) for values in sorted(result)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "validation/reaction_product_parity_cases.json")
    parser.add_argument("--strata", type=Path, default=ROOT / "validation/reaction_product_parity_strata_v2.json")
    parser.add_argument("--rust-rows", type=Path, required=True)
    parser.add_argument("--rust-summary", type=Path, required=True)
    parser.add_argument("--baseline-rows", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rdkit", required=True)
    parser.add_argument("--expected-crate", default="chematic 1.0.30 from crates.io")
    parser.add_argument("--expected-version", default="1.0.30")
    args = parser.parse_args()
    if rdBase.rdkitVersion != args.expected_rdkit:
        parser.error(f"RDKit {rdBase.rdkitVersion} != {args.expected_rdkit}")
    cases, fixture_hashes = load_cases(args.base, args.strata)
    rust_bytes = args.rust_rows.read_bytes()
    rust_summary_bytes = args.rust_summary.read_bytes()
    summary = json.loads(rust_summary_bytes)
    baseline_bytes = args.baseline_rows.read_bytes()
    rust_rows = json.loads(rust_bytes)
    baseline = json.loads(baseline_bytes)
    if (summary["schema"] != "published-rust-reaction-template-maps/v1"
            or summary["crate"] != args.expected_crate
            or summary["input_count"] != len(cases) == 83
            or summary["base_cases_sha256"] != fixture_hashes["base_sha256"]
            or summary["strata_sha256"] != fixture_hashes["strata_sha256"]
            or summary["rows_sha256"] != hashlib.sha256(rust_bytes).hexdigest()
            or [row["id"] for row in rust_rows] != [case["id"] for case in cases]
            or [row["id"] for row in baseline] != [case["id"] for case in cases]):
        raise ValueError("published Rust map rows, baseline, and fixtures disagree")

    rows = []
    for case, rust, plain in zip(cases, rust_rows, baseline, strict=True):
        if rust["status"] != plain["status"]:
            raise ValueError(f"{case['id']}: per-match diagnostic changed published status")
        try:
            oracle_graph, oracle_origins, oracle_raw = rdkit_sets(case)
            oracle_maps, map_raw = rdkit_map_sets(case)
            if oracle_raw != map_raw:
                raise ValueError("RDKit raw product count changed between oracle passes")
        except ValueError as exc:
            outcome = "joint_invalid_input" if rust["status"] == "typed_refusal" else "oracle_invalid"
            rows.append({"id": case["id"], "outcome": outcome, "detail": str(exc),
                         "rust_status": rust["status"]})
            continue
        if rust["status"] != "products":
            outcome = ("typed_unsupported" if rust["status"] == "typed_unsupported"
                       else "typed_or_diagnosed_refusal" if rust["status"] in {"typed_refusal", "diagnosed_valence_refusal"}
                       else "unexpected_status")
            rows.append({"id": case["id"], "outcome": outcome,
                         "rdkit_raw_product_sets": oracle_raw, "rust_status": rust["status"],
                         "reason": rust.get("reason")})
            continue
        candidate_graph, candidate_origins = rust_sets(case, {"sets_with_atom_sources": rust["sets"]})
        baseline_graph, baseline_origins = rust_sets(case, plain)
        if candidate_graph != baseline_graph or candidate_origins != baseline_origins:
            raise ValueError(f"{case['id']}: per-match map runner changed published graph or provenance")
        candidate_maps = rust_map_sets(rust)
        if candidate_graph != oracle_graph:
            outcome = "graph_mismatch"
        elif candidate_origins != oracle_origins:
            outcome = "provenance_mismatch"
        elif candidate_maps != oracle_maps:
            outcome = "map_label_mismatch"
        else:
            outcome = "graph_origin_map_match"
        rows.append({"id": case["id"], "outcome": outcome,
                     "rdkit_graph": oracle_graph, "rust_graph": candidate_graph,
                     "rdkit_origins": oracle_origins, "rust_origins": candidate_origins,
                     "rdkit_maps": oracle_maps, "rust_maps": candidate_maps,
                     "rdkit_raw_product_sets": oracle_raw,
                     "rust_accepted_matches": rust["accepted_matches"],
                     "rust_valence_rejected_matches": rust["valence_rejected_matches"],
                     "rust_raw_product_sets": len(rust["sets"])})
    counts = dict(sorted(Counter(row["outcome"] for row in rows).items()))
    report = {"schema": "published-rust-reaction-template-map-parity/v1",
              "rdkit_version": rdBase.rdkitVersion, "crate_version": args.expected_version,
              "crate_label": args.expected_crate,
              "rust_rows_sha256": hashlib.sha256(rust_bytes).hexdigest(),
              "rust_summary_sha256": hashlib.sha256(rust_summary_bytes).hexdigest(),
              "baseline_rows_sha256": hashlib.sha256(baseline_bytes).hexdigest(),
              "fixtures": fixture_hashes, "accounting": {"input": len(cases), "outcomes": counts},
              "limits": ["Exposed 83-case diagnostic, not broad SMIRKS parity",
                         "Compares RDKit old_mapno with map numbers reconstructed from the accepted CheMatic match and traced product origins",
                         "Does not add atom-map numbers to ordinary CheMatic product molecules",
                         "Raw embedding counts are recorded but not required to match",
                         "Python/npm published bindings do not expose per-product map/origin metadata",
                         "No reaction yield or selectivity claim"],
              "rows": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"reaction map labels: {len(rows)} inputs, {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
