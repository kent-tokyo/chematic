#!/usr/bin/env python3
"""Issue #227 Phase 0.4: corpus-wide MMFF aromatic atom + bond parity report,
joining chematic's `compute_mmff94_aromatic_view` dump against the RDKit
oracle (both on the 265-molecule Wave 1 corpus). Separate from
`mmff94_type_parity_report.py` (numeric types) -- this measures the
upstream aromaticity re-perception step itself.

Run: .venv/bin/python scripts/mmff94_aromaticity_corpus_parity_report_227.py
"""

import json

RDKIT_PATH = "validation/results/mmff94_aromaticity_corpus_parity_227_rdkit.jsonl"
CHEMATIC_PATH = "validation/results/mmff94_aromaticity_corpus_parity_227_chematic.jsonl"


def load(path):
    rows = {}
    with open(path) as f:
        for line in f:
            row = json.loads(line)
            rows[row["name"]] = row
    return rows


def main():
    rdkit_rows = load(RDKIT_PATH)
    chematic_rows = load(CHEMATIC_PATH)
    names = sorted(set(rdkit_rows) | set(chematic_rows))

    atom_total = atom_match = 0
    bond_total = bond_match = 0
    unavailable_molecules = []
    atom_mismatches = []
    bond_mismatches = []

    for name in names:
        r = rdkit_rows.get(name)
        c = chematic_rows.get(name)
        if r is None or r.get("status") != "ok" or c is None or c.get("status") != "ok":
            unavailable_molecules.append(
                {"name": name, "rdkit_status": r.get("status") if r else "missing",
                 "chematic_status": c.get("status") if c else "missing"}
            )
            continue

        for idx, rdkit_val in r["atom_aromatic"].items():
            atom_total += 1
            chematic_val = c["atom_aromatic"].get(idx)
            if chematic_val == rdkit_val:
                atom_match += 1
            else:
                if len(atom_mismatches) < 20:
                    atom_mismatches.append(f"{name}#{idx}: rdkit={rdkit_val} chematic={chematic_val}")

        for key, rdkit_val in r["bond_aromatic"].items():
            bond_total += 1
            chematic_val = c["bond_aromatic"].get(key)
            if chematic_val == rdkit_val:
                bond_match += 1
            else:
                if len(bond_mismatches) < 20:
                    bond_mismatches.append(f"{name}#{key}: rdkit={rdkit_val} chematic={chematic_val}")

    result = {
        "molecules_compared": len(names) - len(unavailable_molecules),
        "molecules_oracle_unavailable": len(unavailable_molecules),
        "atom_aromatic_parity": {
            "total": atom_total,
            "match": atom_match,
            "pct": round(100 * atom_match / atom_total, 4) if atom_total else None,
        },
        "bond_aromatic_parity": {
            "total": bond_total,
            "match": bond_match,
            "pct": round(100 * bond_match / bond_total, 4) if bond_total else None,
        },
    }
    print(json.dumps(result, indent=2))
    print()
    if unavailable_molecules:
        print("-- oracle-unavailable molecules --")
        for u in unavailable_molecules:
            print(f"  {u}")
    if atom_mismatches:
        print("-- atom aromaticity mismatches (up to 20) --")
        for m in atom_mismatches:
            print(f"  {m}")
    if bond_mismatches:
        print("-- bond aromaticity mismatches (up to 20) --")
        for m in bond_mismatches:
            print(f"  {m}")


if __name__ == "__main__":
    main()
