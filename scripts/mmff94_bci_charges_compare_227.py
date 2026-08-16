#!/usr/bin/env python3
"""Issue #227 Phase 2 (BCI investigation): compares chematic's dumped MMFF94
partial charges (`mmff94_bci_charges_dump_227.rs` output) against the live
RDKit oracle (`mmff94_bci_charges_oracle_227.py` output), per heavy atom, for
all 264 typing-succeeded molecules (not a sample -- same "check everything
cheap" discipline Phase 1's torsion investigation used, since this is
topology-only).

Usage:
    .venv/bin/python scripts/mmff94_bci_charges_compare_227.py \\
        <chematic_dump.jsonl> <rdkit_oracle.jsonl> [label]
"""

import json
import sys

TOL = 1e-6  # exact-match tolerance (both sides are the same float arithmetic
            # given the same inputs; float noise is not the phenomenon under test)


def load(path):
    out = {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        row = json.loads(line)
        out[(row["tier"], row["name"])] = row
    return out


def main():
    ch_path, rd_path = sys.argv[1], sys.argv[2]
    label = sys.argv[3] if len(sys.argv) > 3 else "chematic"
    ch = load(ch_path)
    rd = load(rd_path)

    keys = sorted(set(ch) & set(rd))
    n_molecules_compared = 0
    n_atoms_compared = 0
    n_atoms_exact_match = 0
    abs_diffs = []
    per_molecule_max_diff = {}
    mismatched_molecules = []

    for key in keys:
        crow, rrow = ch[key], rd[key]
        if crow["status"] != "ok" or rrow["status"] != "ok":
            continue
        if crow["n_heavy"] != rrow["n_heavy"]:
            print(f"WARN atom-count mismatch {key}: chematic={crow['n_heavy']} rdkit={rrow['n_heavy']}", file=sys.stderr)
            continue
        n_molecules_compared += 1
        max_diff_this_mol = 0.0
        for ci, ri in zip(crow["charges"], rrow["charges"]):
            assert ci["element"] == ri["element"], f"{key}: element mismatch at idx {ci['index']}"
            n_atoms_compared += 1
            d = abs(ci["chematic_partial_charge"] - ri["rdkit_partial_charge"])
            abs_diffs.append(d)
            max_diff_this_mol = max(max_diff_this_mol, d)
            if d < TOL:
                n_atoms_exact_match += 1
        per_molecule_max_diff[f"{key[0]}:{key[1]}"] = max_diff_this_mol
        if max_diff_this_mol > TOL:
            mismatched_molecules.append((key, max_diff_this_mol))

    abs_diffs.sort()
    n = len(abs_diffs)

    def pct(p):
        if n == 0:
            return None
        idx = min(n - 1, int(p * n))
        return abs_diffs[idx]

    mismatched_molecules.sort(key=lambda x: -x[1])

    summary = {
        "label": label,
        "molecules_compared": n_molecules_compared,
        "atoms_compared": n_atoms_compared,
        "atoms_exact_match": n_atoms_exact_match,
        "atoms_mismatched": n_atoms_compared - n_atoms_exact_match,
        "molecules_with_any_mismatch": len(mismatched_molecules),
        "mean_abs_diff": sum(abs_diffs) / n if n else None,
        "median_abs_diff": pct(0.5),
        "p90_abs_diff": pct(0.9),
        "p99_abs_diff": pct(0.99),
        "max_abs_diff": abs_diffs[-1] if n else None,
        "top_20_mismatched_molecules": [
            {"molecule": f"{k[0]}:{k[1]}", "max_abs_diff": d} for k, d in mismatched_molecules[:20]
        ],
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
