#!/usr/bin/env python3
"""Compare published Rust reaction-product atom origins with pinned RDKit.

This is an exposed diagnostic, not a broad reaction-accuracy claim. The Rust
runner tags every input atom and serializes the tags in product SMILES order.
RDKit's `react_idx`/`react_atom_idx` properties provide the independent oracle.
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
    from .run_reaction_compatibility_v2 import load_cases
else:
    from run_reaction_compatibility_v2 import load_cases

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]


def labelled_smiles(mol: Chem.Mol, origins: list[tuple[int, int] | None], offsets: list[int]) -> str:
    if mol.GetNumAtoms() != len(origins):
        raise ValueError("product atom count and origin count differ")
    labelled = Chem.Mol(mol)
    for atom, origin in zip(labelled.GetAtoms(), origins, strict=True):
        if origin is None:
            atom.SetAtomMapNum(0)
        else:
            reactant, index = origin
            if reactant < 0 or reactant + 1 >= len(offsets) or not 0 <= index < offsets[reactant + 1] - offsets[reactant]:
                raise ValueError("source atom is outside the reactant fixture")
            atom.SetAtomMapNum(offsets[reactant] + index + 1)
    return Chem.MolToSmiles(labelled, canonical=True, isomericSmiles=True)


def rdkit_sets(case: dict) -> tuple[list[list[str]], list[list[str]], int]:
    reaction = AllChem.ReactionFromSmarts(case["smirks"])
    if reaction is None:
        raise ValueError("RDKit cannot parse SMIRKS")
    reactants = [Chem.MolFromSmiles(text) for text in case["reactants"]]
    if any(mol is None for mol in reactants):
        raise ValueError("RDKit cannot parse a reactant")
    offsets = [0]
    for mol in reactants:
        offsets.append(offsets[-1] + mol.GetNumAtoms())
    raw = reaction.RunReactants(tuple(reactants))
    graphs, origins = set(), set()
    for product_set in raw:
        graph_tuple, origin_tuple = [], []
        for product in product_set:
            checked = Chem.Mol(product)
            Chem.SanitizeMol(checked)
            graph_tuple.append(Chem.MolToSmiles(checked, canonical=True, isomericSmiles=True))
            source_values = []
            for atom in checked.GetAtoms():
                has_index = atom.HasProp("react_atom_idx")
                has_reactant = atom.HasProp("react_idx")
                if has_index != has_reactant:
                    raise ValueError("RDKit product has incomplete origin properties")
                source_values.append((atom.GetIntProp("react_idx"), atom.GetIntProp("react_atom_idx")) if has_index else None)
            origin_tuple.append(labelled_smiles(checked, source_values, offsets))
        graphs.add(tuple(sorted(graph_tuple)))
        origins.add(tuple(sorted(origin_tuple)))
    return [list(x) for x in sorted(graphs)], [list(x) for x in sorted(origins)], len(raw)


def rust_sets(case: dict, row: dict) -> tuple[list[list[str]], list[list[str]]]:
    offsets = [0]
    for text in case["reactants"]:
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            raise ValueError("RDKit cannot parse a reactant")
        offsets.append(offsets[-1] + mol.GetNumAtoms())
    graphs, origins = set(), set()
    for product_set in row["sets_with_atom_sources"]:
        graph_tuple, origin_tuple = [], []
        for product in product_set:
            mol = Chem.MolFromSmiles(product["smiles"])
            if mol is None:
                raise ValueError("published Rust product is not RDKit-readable")
            graph_tuple.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
            source_values = [tuple(x) if x is not None else None for x in product["atom_sources"]]
            origin_tuple.append(labelled_smiles(mol, source_values, offsets))
        graphs.add(tuple(sorted(graph_tuple)))
        origins.add(tuple(sorted(origin_tuple)))
    return [list(x) for x in sorted(graphs)], [list(x) for x in sorted(origins)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "validation/reaction_product_parity_cases.json")
    parser.add_argument("--strata", type=Path, default=ROOT / "validation/reaction_product_parity_strata_v2.json")
    parser.add_argument("--rust-rows", type=Path, required=True)
    parser.add_argument("--rust-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rdkit", required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != args.expected_rdkit:
        parser.error(f"RDKit {rdBase.rdkitVersion} != {args.expected_rdkit}")
    cases, fixture_hashes = load_cases(args.base, args.strata)
    rust_bytes = args.rust_rows.read_bytes()
    rust_summary = json.loads(args.rust_summary.read_text(encoding="utf-8"))
    if (rust_summary.get("schema") != "published-rust-reactions/v2"
            or rust_summary.get("crate") != "chematic 1.0.30 from crates.io"
            or rust_summary.get("base_cases_sha256") != fixture_hashes["base_sha256"]
            or rust_summary.get("strata_sha256") != fixture_hashes["strata_sha256"]
            or rust_summary.get("rows_sha256") != hashlib.sha256(rust_bytes).hexdigest()):
        raise ValueError("published-crate reaction summary does not pin these rows and fixtures")
    rust_rows = json.loads(rust_bytes)
    if [row["id"] for row in rust_rows] != [case["id"] for case in cases]:
        raise ValueError("Rust rows do not cover the frozen fixtures in order")
    rows = []
    for case, rust in zip(cases, rust_rows, strict=True):
        try:
            oracle_graph, oracle_sources, oracle_raw = rdkit_sets(case)
        except ValueError as exc:
            outcome = "joint_invalid_input" if rust["status"] in {"typed_refusal", "diagnosed_valence_refusal"} else "oracle_invalid"
            rows.append({"id": case["id"], "outcome": outcome, "detail": str(exc), "rust_status": rust["status"]})
            continue
        if rust["status"] != "products":
            outcome = "typed_or_diagnosed_refusal" if rust["status"] in {"typed_refusal", "diagnosed_valence_refusal"} else "unexpected_status"
            rows.append({"id": case["id"], "outcome": outcome, "rust_status": rust["status"], "rdkit_raw_product_sets": oracle_raw})
            continue
        candidate_graph, candidate_sources = rust_sets(case, rust)
        recorded_graph = sorted({tuple(sorted(Chem.MolToSmiles(Chem.MolFromSmiles(text), canonical=True,
                                                               isomericSmiles=True) for text in product_set))
                                 for product_set in rust["sets"]})
        if sorted(tuple(product_set) for product_set in candidate_graph) != recorded_graph:
            raise ValueError(f"{case['id']}: provenance rows changed the ordinary product graph")
        if candidate_graph != oracle_graph:
            outcome = "graph_mismatch"
        elif candidate_sources != oracle_sources:
            outcome = "provenance_mismatch"
        else:
            outcome = "graph_and_provenance_match"
        rows.append({"id": case["id"], "outcome": outcome, "strata": case["strata"],
                     "rdkit_graph": oracle_graph, "rust_graph": candidate_graph,
                     "rdkit_origins": oracle_sources, "rust_origins": candidate_sources,
                     "rdkit_raw_product_sets": oracle_raw,
                     "rust_raw_product_sets": rust["raw_product_sets"]})
    counts = dict(sorted(Counter(row["outcome"] for row in rows).items()))
    report = {"schema": "published-rust-reaction-atom-provenance/v1",
              "rdkit_version": rdBase.rdkitVersion, "crate_version": "1.0.30",
              "rust_rows_sha256": hashlib.sha256(rust_bytes).hexdigest(),
              "rust_summary_sha256": hashlib.sha256(args.rust_summary.read_bytes()).hexdigest(),
              "fixtures": fixture_hashes, "accounting": {"input": len(cases), "outcomes": counts},
              "limits": ["Exposed 83-case diagnostic, not broad SMIRKS parity",
                         "Does not compare original template map labels directly",
                         "Raw embedding count is recorded but not required to match",
                         "Python/npm published bindings do not expose per-product origins",
                         "No reaction yield or selectivity claim"],
              "rows": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"reaction provenance: {len(rows)} inputs, {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
