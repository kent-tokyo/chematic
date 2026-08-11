#!/usr/bin/env python3
"""Row-level chematic-vs-RDKit comparison for the platinum benchmark.

Reads chematic's own harness output (`platinum_after_fix_chematic.jsonl`,
produced by `cargo run -p chematic-mol --example platinum_benchmark`) and
the independent RDKit oracle output (`platinum_rdkit_oracle.jsonl`,
produced by `scripts/platinum_rdkit_oracle.py`), and reports, per row:

- whether chematic's `formula_matches_expected` and RDKit's own (computed
  with a completely independent regex-based formula parser, not chematic's
  `chematic_chem::formula::parse_formula`) agree on the SAME
  `formula_expected` string -- if both independently confirm it, the
  corpus's expected value is doubly cross-checked, not merely
  self-consistent with chematic's own parser.
- for the 6 corpus rows with no independently-sourced
  `exact_mass_expected` (see `pt_corpus.jsonl`'s `source.note` on those
  rows), chematic's computed mass vs RDKit's computed mass on the exact
  same input SMILES, since no external expected value exists for them.

Usage: python scripts/platinum_compare_chematic_rdkit.py
"""

import json

CHEMATIC_PATH = "validation/results/platinum_after_fix_chematic.jsonl"
RDKIT_PATH = "validation/results/platinum_rdkit_oracle.jsonl"
CORPUS_PATH = "validation/platinum/pt_corpus.jsonl"


def load_jsonl(path):
    return {json.loads(line)["id"]: json.loads(line) for line in open(path) if line.strip()}


def main():
    chematic = load_jsonl(CHEMATIC_PATH)
    rdkit = load_jsonl(RDKIT_PATH)
    corpus = load_jsonl(CORPUS_PATH)

    print(f"{'id':40}{'chem_formula_ok':>17}{'rdkit_formula_ok':>18}{'agree':>8}")
    both_confirm = 0
    disagreements = []
    for cid in corpus:
        c = chematic.get(cid, {})
        r = rdkit.get(cid, {})
        c_ok = c.get("formula_matches_expected")
        r_ok = r.get("formula_matches_expected")
        agree = c_ok == r_ok
        if c_ok and r_ok:
            both_confirm += 1
        if not agree:
            disagreements.append(cid)
        print(f"{cid:40}{str(c_ok):>17}{str(r_ok):>18}{str(agree):>8}")

    print()
    print(f"formula_expected independently confirmed by BOTH chematic and RDKit: {both_confirm}/{len(corpus)}")
    if disagreements:
        print(f"DISAGREEMENTS (chematic's formula_matches_expected != RDKit's): {disagreements}")
    else:
        print("No disagreements: every row where chematic says the expected formula is right, RDKit's own independent formula parser says so too, and vice versa.")

    print()
    print("=== chematic vs RDKit computed mass, for the 6 rows with no independently-sourced exact_mass_expected ===")
    print(f"{'id':40}{'chematic_mass':>16}{'rdkit_mass':>14}{'abs_diff':>10}")
    null_mass_rows = [cid for cid, row in corpus.items() if row.get("exact_mass_expected") is None]
    max_diff = 0.0
    for cid in null_mass_rows:
        c_mass = chematic.get(cid, {}).get("exact_mass")
        r_mass = rdkit.get(cid, {}).get("exact_mass")
        if c_mass is None or r_mass is None:
            print(f"{cid:40}{'MISSING':>16}{'MISSING':>14}")
            continue
        diff = abs(c_mass - r_mass)
        max_diff = max(max_diff, diff)
        print(f"{cid:40}{c_mass:>16.4f}{r_mass:>14.4f}{diff:>10.4f}")
    print(f"\nmax |chematic - rdkit| across these {len(null_mass_rows)} rows: {max_diff:.4f} Da")


if __name__ == "__main__":
    main()
