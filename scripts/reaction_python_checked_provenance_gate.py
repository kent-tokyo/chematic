#!/usr/bin/env python3
"""Gate checked Python reaction graphs, atom origins and template maps on 83 fixtures.

This is a source-candidate gate. It does not promote a local extension or wheel
to published-package evidence, and raw reaction embedding counts are diagnostic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from zipfile import ZipFile

from rdkit import Chem, RDLogger, rdBase

if __package__:
    from .reaction_atom_provenance_gate import rdkit_sets, rust_sets
    from .reaction_template_map_gate import rdkit_map_sets, rust_map_sets
    from .run_reaction_compatibility_v2 import load_cases
else:
    from reaction_atom_provenance_gate import rdkit_sets, rust_sets
    from reaction_template_map_gate import rdkit_map_sets, rust_map_sets
    from run_reaction_compatibility_v2 import load_cases

ROOT = Path(__file__).resolve().parents[1]
RDLogger.DisableLog("rdApp.*")
EXPECTED = {
    # v2_isotope_methanol_split matches RDKit since #734 (was a refusal).
    "graph_origin_map_match": 77,
    "typed_unsupported": 3,
    "joint_invalid_input": 3,
}


def source_artifact(chematic, wheel: Path | None) -> dict:
    if wheel is None:
        return {"kind": "editable_source_extension", "published": False}
    if wheel.suffix != ".whl" or not wheel.is_file():
        raise ValueError("--wheel must identify one existing wheel file")
    package_dir = Path(chematic.__file__).resolve().parent
    extensions = list(package_dir.glob("chematic*.so")) + list(package_dir.glob("chematic*.pyd"))
    if len(extensions) != 1:
        raise ValueError("installed chematic extension is missing or ambiguous")
    installed_digest = hashlib.sha256(extensions[0].read_bytes()).hexdigest()
    with ZipFile(wheel) as archive:
        entries = [name for name in archive.namelist()
                   if name.endswith("/" + extensions[0].name)]
        if len(entries) != 1:
            raise ValueError("wheel extension is missing or ambiguous")
        wheel_extension_digest = hashlib.sha256(archive.read(entries[0])).hexdigest()
    if installed_digest != wheel_extension_digest:
        raise ValueError("installed extension does not match the supplied wheel")
    return {"kind": "release_profile_source_wheel", "published": False,
            "wheel": wheel.name, "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
            "extension_sha256": installed_digest}


def python_row(chematic, case: dict) -> dict:
    try:
        reactants = [chematic.from_smiles(value) for value in case["reactants"]]
    except ValueError as exc:
        return {"status": "typed_refusal", "reason": "reactant_parse", "detail": str(exc)}
    checked = chematic.run_smirks_checked(case["smirks"], reactants, rdkit_compat=True)
    status = checked["status"]
    if status not in {"products", "no_match"}:
        return {"status": status, "reason": checked["reason"],
                "detail": checked["detail"], "diagnostics": {key: checked[key] for key in
                ("accepted_matches", "applied_products", "valence_rejected_matches", "truncated_matches")}}
    if status == "no_match" and checked["products"]:
        raise ValueError(f"{case['id']}: no_match returned products")
    if not (len(checked["products"]) == len(checked["product_atom_sources"])
            == len(checked["product_template_maps"])):
        raise ValueError(f"{case['id']}: product-set metadata length differs")
    sets = []
    for product_set, source_set, map_set in zip(
            checked["products"], checked["product_atom_sources"],
            checked["product_template_maps"], strict=True):
        if not len(product_set) == len(source_set) == len(map_set):
            raise ValueError(f"{case['id']}: product metadata length differs")
        items = []
        for product, sources, maps in zip(product_set, source_set, map_set, strict=True):
            smiles, order = product.smiles_with_atom_order()
            if (sorted(order) != list(range(len(order)))
                    or len(sources) != len(order) or len(maps) != len(order)):
                raise ValueError(f"{case['id']}: canonical atom order is not a metadata permutation")
            parsed = Chem.MolFromSmiles(smiles)
            if parsed is None or parsed.GetNumAtoms() != len(order):
                raise ValueError(f"{case['id']}: canonical product is not RDKit-readable")
            items.append({"smiles": smiles,
                          "atom_sources": [sources[index] for index in order],
                          "template_map_numbers": [maps[index] for index in order]})
        sets.append(items)
    return {"status": status, "sets": sets,
            "diagnostics": {key: checked[key] for key in
                            ("accepted_matches", "applied_products", "valence_rejected_matches", "truncated_matches")}}


def classify(case: dict, candidate: dict) -> dict:
    status = candidate["status"]
    try:
        oracle_graph, oracle_origins, oracle_count = rdkit_sets(case)
        oracle_maps, map_count = rdkit_map_sets(case)
        if oracle_count != map_count:
            raise ValueError("RDKit reaction count changed between oracle passes")
    except ValueError as exc:
        outcome = "joint_invalid_input" if status == "typed_refusal" else "oracle_invalid"
        return {"id": case["id"], "outcome": outcome, "status": status, "detail": str(exc)}
    if status == "typed_unsupported":
        outcome = "typed_unsupported"
    elif status in {"typed_refusal", "partial_products"}:
        outcome = "typed_or_diagnosed_refusal" if candidate.get("reason") == "product_valence" else "unexpected_refusal"
    elif status in {"products", "no_match"}:
        graph, origins = rust_sets(case, {"sets_with_atom_sources": candidate["sets"]})
        maps = rust_map_sets({"sets": candidate["sets"]})
        outcome = ("graph_mismatch" if graph != oracle_graph else
                   "provenance_mismatch" if origins != oracle_origins else
                   "map_label_mismatch" if maps != oracle_maps else
                   "graph_origin_map_match")
    else:
        outcome = "unexpected_status"
    return {"id": case["id"], "outcome": outcome, "status": status,
            "reason": candidate.get("reason"), "diagnostics": candidate.get("diagnostics"),
            "rdkit_raw_product_sets": oracle_count}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "validation/reaction_product_parity_cases.json")
    parser.add_argument("--strata", type=Path, default=ROOT / "validation/reaction_product_parity_strata_v2.json")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--wheel", type=Path,
                        help="release-profile source wheel installed in this interpreter")
    parser.add_argument("--expected-rdkit", default="2026.03.6")
    args = parser.parse_args()
    if rdBase.rdkitVersion != args.expected_rdkit:
        parser.error(f"RDKit {rdBase.rdkitVersion} != {args.expected_rdkit}")
    import chematic

    artifact = source_artifact(chematic, args.wheel)

    cases, hashes = load_cases(args.base, args.strata)
    if len(cases) != 83:
        raise ValueError(f"reaction fixture count changed: {len(cases)}")
    rows = [classify(case, python_row(chematic, case)) for case in cases]
    counts = dict(sorted(Counter(row["outcome"] for row in rows).items()))
    report = {"schema": "source-python-checked-reaction-provenance/v1",
              "rdkit_version": rdBase.rdkitVersion, "fixtures": hashes,
              "artifact": artifact,
              "accounting": {"input": len(cases), "outcomes": counts}, "rows": rows,
              "limits": ["83 pinned cases only", "Raw embeddings are not required to match",
                         "No yield, selectivity or broad SMIRKS parity claim"]}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"{args.output}: sha256={hashlib.sha256(rendered.encode()).hexdigest()}")
    else:
        print(rendered)
    if counts != EXPECTED:
        raise ValueError(f"checked Python outcome accounting differs: {counts} != {EXPECTED}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
