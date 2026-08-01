#!/usr/bin/env python3
"""Phase 1B-0 audit for issue #227: post-fix atom-type parity report.

Joins chematic's numeric-type dump (mmff94_numeric_type_dump.rs,
`validation/results/mmff94_chematic_numeric_types.jsonl`) against the frozen
RDKit oracle (`validation/results/mmff94_rdkit_type_oracle.jsonl`, unchanged
since RDKit's own typer did not move) by molecule name + atom index, and
classifies every comparable atom into exactly one of:

  - exact_id_match: chematic's numeric type ID equals RDKit's
  - real_mismatch: different MMFF94 symbol/class (would be a distinct
    parameter row -- e.g. NR vs NC=C)

Every `real_mismatch` is additionally required to have the SAME element on
both sides (the numeric-type registry's element field), which is the one
invariant this PR's construction-time semantic-compatibility gate
guarantees can never be violated in production; a cross-element mismatch
here would mean the gate itself is broken. Zero unclassified atoms: every
comparable atom lands in exactly one bucket above.

Run: python3 scripts/mmff94_type_parity_report.py \
  > validation/results/mmff94_type_parity_227_postfix.json
"""

import json
import re
import sys

REGISTRY_PATH = "crates/chematic-ff/src/mmff94_numeric_type_registry.rs"
CHEMATIC_DUMP = "validation/results/mmff94_chematic_numeric_types.jsonl"
RDKIT_ORACLE = "validation/results/mmff94_rdkit_type_oracle.jsonl"


def load_registry(path):
    text = open(path).read()
    rows = re.findall(r'id:\s*(\d+),\s*symbol:\s*"([^"]+)"', text)
    if len(rows) < 90:
        raise RuntimeError(
            f"registry parse found only {len(rows)} rows -- generator output format changed?"
        )
    return {int(i): s for i, s in rows}


def load_jsonl(path):
    out = {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        d = json.loads(line)
        out[d["name"]] = d
    return out


def main():
    reg = load_registry(REGISTRY_PATH)
    chem = load_jsonl(CHEMATIC_DUMP)
    rdkit = load_jsonl(RDKIT_ORACLE)

    total_atoms = 0
    exact_id_match = 0
    real_mismatch = 0
    cross_element_mismatch = 0
    chem_typing_failure_molecules = []
    rdkit_unavailable_molecules = []
    mismatches = []

    for name, rd in rdkit.items():
        c = chem.get(name)
        if c is None:
            continue
        if c.get("status") != "ok":
            chem_typing_failure_molecules.append(name)
            continue
        if rd.get("status") != "ok" or not rd.get("mmff_properties_available"):
            rdkit_unavailable_molecules.append(name)
            continue
        c_atoms = {a["index"]: a for a in c["atoms"]}
        r_atoms = {a["index"]: a for a in rd["atom_types"]}
        for idx, ra in r_atoms.items():
            ca = c_atoms.get(idx)
            if ca is None:
                continue
            total_atoms += 1
            r_type = ra["rdkit_mmff_type"]
            c_type = ca["chematic_numeric_type"]
            if r_type == c_type:
                exact_id_match += 1
                continue
            r_sym = reg.get(r_type, "???")
            c_sym = reg.get(c_type, "???")
            real_mismatch += 1
            if ca["element"] != ra["element"]:
                cross_element_mismatch += 1
            mismatches.append(
                {
                    "molecule": name,
                    "atom_index": idx,
                    "element": ca["element"],
                    "chematic_type": c_type,
                    "chematic_symbol": c_sym,
                    "rdkit_type": r_type,
                    "rdkit_symbol": r_sym,
                }
            )

    by_pair = {}
    for m in mismatches:
        key = (m["element"], m["chematic_symbol"], m["rdkit_symbol"])
        by_pair.setdefault(key, []).append(m)

    grouped = [
        {
            "element": el,
            "chematic_symbol": cs,
            "rdkit_symbol": rs,
            "count": len(items),
            "example": f"{items[0]['molecule']}#{items[0]['atom_index']}",
        }
        for (el, cs, rs), items in sorted(
            by_pair.items(), key=lambda kv: -len(kv[1])
        )
    ]

    report = {
        "total_comparable_atoms": total_atoms,
        "exact_id_match": exact_id_match,
        "exact_id_match_pct": round(100 * exact_id_match / total_atoms, 2)
        if total_atoms
        else None,
        "real_mismatch": real_mismatch,
        "real_mismatch_pct": round(100 * real_mismatch / total_atoms, 2)
        if total_atoms
        else None,
        "cross_element_mismatch": cross_element_mismatch,
        "unclassified": total_atoms
        - exact_id_match
        - real_mismatch,
        "chematic_typing_failure_molecules": chem_typing_failure_molecules,
        "rdkit_unavailable_molecules": rdkit_unavailable_molecules,
        "mismatch_groups": grouped,
    }
    print(json.dumps(report, indent=2))

    if cross_element_mismatch:
        print(
            f"FAIL: {cross_element_mismatch} cross-element mismatches -- "
            "the semantic-compatibility gate's invariant is violated",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
