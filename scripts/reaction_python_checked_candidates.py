#!/usr/bin/env python3
"""Emit checked Python reaction rows for the 83 fixtures without importing RDKit.

Some published chematic wheels exist only for interpreters RDKit 2026.03.6
does not ship for (the v1.0.34 Linux x86-64 wheel is CPython 3.9 only). This
script runs ``run_smirks_checked(..., rdkit_compat=True)`` in that interpreter
and writes the candidate rows; ``reaction_python_checked_provenance_gate.py
--candidates`` classifies them against RDKit in another interpreter.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import sys
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
DIAGNOSTICS = ("accepted_matches", "applied_products", "valence_rejected_matches", "truncated_matches")


def load_cases(base_path: Path, strata_path: Path):
    """Same fixture pinning as ``run_reaction_compatibility_v2.load_cases``."""
    base = json.loads(base_path.read_text(encoding="utf-8"))
    strata = json.loads(strata_path.read_text(encoding="utf-8"))
    base_sha = hashlib.sha256(base_path.read_bytes()).hexdigest()
    if strata.get("schema_version") != 2 or strata.get("base_fixture_sha256") != base_sha:
        raise ValueError("strata manifest does not pin the historical base fixture")
    cases = [{**case, "strata": ["legacy"]} for case in base] + strata["cases"]
    ids = [case["id"] for case in cases]
    if len(base) != 57 or len(ids) != len(set(ids)):
        raise ValueError("legacy fixture count or case IDs changed")
    return cases, {"base_sha256": base_sha,
                   "strata_sha256": hashlib.sha256(strata_path.read_bytes()).hexdigest()}


def python_row(chematic, case: dict, validate=None) -> dict:
    """One checked row; ``validate(case_id, smiles, atom_count)`` may reject a product."""
    try:
        reactants = [chematic.from_smiles(value) for value in case["reactants"]]
    except ValueError as exc:
        return {"status": "typed_refusal", "reason": "reactant_parse", "detail": str(exc)}
    checked = chematic.run_smirks_checked(case["smirks"], reactants, rdkit_compat=True)
    status = checked["status"]
    if status not in {"products", "no_match"}:
        return {"status": status, "reason": checked["reason"], "detail": checked["detail"],
                "diagnostics": {key: checked[key] for key in DIAGNOSTICS}}
    if status == "no_match" and checked["products"]:
        raise ValueError(f"{case['id']}: no_match returned products")
    if not (len(checked["products"]) == len(checked["product_atom_sources"])
            == len(checked["product_template_maps"])):
        raise ValueError(f"{case['id']}: product-set metadata length differs")
    sets = []
    for product_set, source_set, map_set in zip(
            checked["products"], checked["product_atom_sources"], checked["product_template_maps"]):
        if not len(product_set) == len(source_set) == len(map_set):
            raise ValueError(f"{case['id']}: product metadata length differs")
        items = []
        for product, sources, maps in zip(product_set, source_set, map_set):
            smiles, order = product.smiles_with_atom_order()
            if (sorted(order) != list(range(len(order)))
                    or len(sources) != len(order) or len(maps) != len(order)):
                raise ValueError(f"{case['id']}: canonical atom order is not a metadata permutation")
            if validate is not None:
                validate(case["id"], smiles, len(order))
            items.append({"smiles": smiles,
                          "atom_sources": [list(sources[i]) if sources[i] is not None else None
                                           for i in order],
                          "template_map_numbers": [maps[i] for i in order]})
        sets.append(items)
    return {"status": status, "sets": sets,
            "diagnostics": {key: checked[key] for key in DIAGNOSTICS}}


def wheel_artifact(chematic, wheel: Path) -> dict:
    """Prove the imported extension is the one inside ``wheel``."""
    package_dir = Path(chematic.__file__).resolve().parent
    extensions = list(package_dir.glob("chematic*.so")) + list(package_dir.glob("chematic*.pyd"))
    if len(extensions) != 1:
        raise ValueError("installed chematic extension is missing or ambiguous")
    installed = hashlib.sha256(extensions[0].read_bytes()).hexdigest()
    with ZipFile(wheel) as archive:
        entries = [n for n in archive.namelist() if n.endswith("/" + extensions[0].name)]
        if len(entries) != 1:
            raise ValueError("wheel extension is missing or ambiguous")
        if hashlib.sha256(archive.read(entries[0])).hexdigest() != installed:
            raise ValueError("installed extension does not match the supplied wheel")
    return {"wheel": wheel.name, "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
            "extension_sha256": installed}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT / "validation/reaction_product_parity_cases.json")
    parser.add_argument("--strata", type=Path, default=ROOT / "validation/reaction_product_parity_strata_v2.json")
    parser.add_argument("--wheel", type=Path, required=True,
                        help="the wheel installed in this interpreter")
    parser.add_argument("--published-from", help="where the wheel was downloaded (e.g. PyPI URL)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    import chematic

    cases, hashes = load_cases(args.base, args.strata)
    artifact = {"kind": "published_wheel" if args.published_from else "wheel",
                "published": bool(args.published_from), "published_from": args.published_from,
                "chematic_version": importlib.metadata.version("chematic"),
                "python": platform.python_version(), "platform": sys.platform,
                **wheel_artifact(chematic, args.wheel)}
    rows = {case["id"]: python_row(chematic, case) for case in cases}
    report = {"schema": "python-checked-reaction-candidates/v1", "fixtures": hashes,
              "artifact": artifact, "rows": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{args.output}: {len(rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
