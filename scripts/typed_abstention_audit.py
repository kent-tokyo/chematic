#!/usr/bin/env python3
"""Probe how the Python binding reports each kind of abstention.

One probe per (operation family, category) with an input built to trigger
it: ``malformed`` (input that is not the format), ``unsupported`` (valid
input outside the implemented scope), ``ambiguous`` (input with more than one
reading) and ``resource_limit`` (input past a declared bound). Each probe
records whether the call returned, returned a typed envelope (``code`` /
``status`` / ``reason`` field), or raised, and with which class and fields.
``scripts/typed_abstention_audit.mjs`` runs the same probes on the WASM
package. The result is the matrix in ``docs/error-and-limits.md``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import chematic

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

NUCLEIC = json.loads((ROOT / "validation/nucleic_acid_document_contract.json").read_text(encoding="utf-8"))


def nucleic_case(code):
    return next(c["document"] for c in NUCLEIC["invalid_cases"] if c["expected_code"] == code)


def over_limit_nucleic():
    doc = {"schema": "chematic.nucleic-acid.v1", "atom_ids": [], "strands": [], "linkages": [], "annotations": {}}
    for i in range(65):
        doc["atom_ids"].append(f"a{i}")
        doc["strands"].append({"id": f"s{i}", "kind": "dna", "annotations": {}, "residues": [
            {"id": f"r{i}", "base": "A", "sugar": "deoxyribose", "atom_refs": [f"a{i}"], "annotations": {}}]})
    return doc


def pipeline(smiles):
    from public_package_3d_chematic import pipeline_config

    return chematic.from_smiles(smiles).embed_pipeline_v2(pipeline_config(chematic, "mmff94_bond_angle_strict",
                                                                          stereo_safe=True))


PROBES = [
    ("smiles", "malformed", lambda: chematic.from_smiles("C1CC(")),
    ("smarts", "malformed", lambda: chematic.from_smiles("CC").has_substructure("[C;")),
    ("mol_block", "malformed", lambda: chematic.from_mol_block("not a mol block")),
    ("inchi", "malformed", lambda: chematic.from_inchi("InChI=1S/garbage")),
    ("smirks", "malformed", lambda: chematic.run_smirks_checked("[C:1>>", [chematic.from_smiles("CC")], rdkit_compat=True)),
    ("smirks", "unsupported", lambda: chematic.run_smirks_checked(
        "[C:1].[O:2]>>[C:1][O:2]", [chematic.from_smiles("C"), chematic.from_smiles("O"), chematic.from_smiles("N")],
        rdkit_compat=True)),
    ("3d", "unsupported", lambda: pipeline("Cl[Pt](Cl)([NH3])[NH3]")),
    ("nucleic_acid", "ambiguous", lambda: json.loads(chematic.nucleic_acid_validate_json(
        json.dumps(nucleic_case("ambiguous_atom_mapping"))))),
    ("nucleic_acid", "unsupported", lambda: json.loads(chematic.nucleic_acid_validate_json(
        json.dumps(nucleic_case("unsupported_linkage_topology"))))),
    ("nucleic_acid", "resource_limit", lambda: json.loads(chematic.nucleic_acid_validate_json(
        json.dumps(over_limit_nucleic())))),
    ("mmcif", "resource_limit", lambda: chematic.parse_mmcif("data_x\n" * 100, max_input_bytes=10)),
]


def describe(value):
    if isinstance(value, dict):
        typed = {k: value[k] for k in ("ok", "status", "reason", "code") if k in value}
        if "error" in value and isinstance(value["error"], dict):
            typed["error.code"] = value["error"].get("code")
        return {"outcome": "envelope" if typed else "returned", "fields": typed}
    return {"outcome": "returned", "fields": {}}


def run():
    rows = []
    for family, category, call in PROBES:
        row = {"family": family, "category": category}
        try:
            row.update(describe(call()))
        except TypeError as error:
            row.update({"outcome": "probe_error", "message": str(error)[:160]})
        except Exception as error:  # the representation is what is audited
            fields = {k: str(getattr(error, k))[:80] for k in ("code", "category", "kind", "reason", "path")
                      if hasattr(error, k)}
            diagnostics = getattr(error, "diagnostics", None)
            if isinstance(diagnostics, dict):
                cause = diagnostics.get("cause")
                fields["diagnostics.cause"] = (cause.get("kind") if isinstance(cause, dict) else str(cause))[:80]
            row.update({"outcome": "raised", "class": type(error).__name__, "fields": fields,
                        "message": str(error)[:160]})
        rows.append(row)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    rows = run()
    args.output.write_text(json.dumps({"schema": "typed-abstention-audit/v1", "binding": "python",
                                       "chematic_version": chematic.__version__, "probes": rows},
                                      indent=1) + "\n", encoding="utf-8")
    for r in rows:
        print(f"{r['family']:13s} {r['category']:15s} {r['outcome']:9s} {r.get('class', '')} {r.get('fields')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
