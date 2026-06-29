#!/usr/bin/env python3
"""Canonical SMILES differential validation: chematic vs RDKit.

Canonical SMILES strings are algorithm-specific, so chematic and RDKit will
NOT produce identical strings. The meaningful test is *semantic round-trip
equivalence*: when chematic's canonical SMILES is re-parsed by RDKit, does it
canonicalize to the same molecule as RDKit's native canonicalization of the
original input?

For each input SMILES s:
  cm   = chematic.from_smiles(s).smiles          # chematic canonical
  rd*  = RDKit.MolToSmiles(RDKit.MolFromSmiles(s))   # RDKit canonical of input
  rd∘cm = RDKit.MolToSmiles(RDKit.MolFromSmiles(cm)) # RDKit canonical of cm
  round-trip OK  ⇔  rd* == rd∘cm

Also checks chematic idempotency: canonical(canonical(s)) == canonical(s).

Usage:
    python scripts/canonical_diff.py [SMILES.csv] [--limit N]

Writes validation/results/canonical_diff.jsonl (one row per divergence).
"""
import json
import os
import sys

import chematic

try:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    sys.exit("RDKit is required for canonical_diff.py (pip install rdkit)")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "validation", "results", "canonical_diff.jsonl")


def has_ez(*strings):
    return any("/" in s or "\\" in s for s in strings)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") \
        else os.path.expanduser("~/Downloads/SMILES.csv")
    limit = None
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])

    smis = [l.strip() for l in open(path) if l.strip()]
    if limit:
        smis = smis[:limit]

    n = 0
    unparseable = 0
    rt_match = 0
    rt_mismatch = 0
    rt_mismatch_ez = 0
    idem_ok = 0
    idem_fail = 0
    idem_fail_ez = 0
    rows = []

    for s in smis:
        rm = Chem.MolFromSmiles(s)
        if rm is None:
            continue  # RDKit can't parse the reference; skip
        try:
            cm = chematic.from_smiles(s).smiles
        except Exception:
            continue
        n += 1

        # idempotency (chematic-internal)
        try:
            cm2 = chematic.from_smiles(cm).smiles
        except Exception:
            cm2 = None
        if cm2 == cm:
            idem_ok += 1
        else:
            idem_fail += 1
            ez = has_ez(s, cm)
            idem_fail_ez += ez
            rows.append({"smiles": s, "metric": "idempotency",
                         "pass1": cm, "pass2": cm2, "ez": ez, "ok": False})

        # round-trip semantic equivalence vs RDKit
        rd_of_cm = Chem.MolFromSmiles(cm)
        if rd_of_cm is None:
            unparseable += 1
            rows.append({"smiles": s, "metric": "rdkit_parse",
                         "chematic": cm, "ez": has_ez(cm), "ok": False})
            continue
        rd_native = Chem.MolToSmiles(rm)
        rd_of_cm_canon = Chem.MolToSmiles(rd_of_cm)
        if rd_native == rd_of_cm_canon:
            rt_match += 1
        else:
            rt_mismatch += 1
            ez = has_ez(s, cm)
            rt_mismatch_ez += ez
            rows.append({"smiles": s, "metric": "roundtrip",
                         "chematic": cm, "rdkit_native": rd_native,
                         "rdkit_of_chematic": rd_of_cm_canon, "ez": ez,
                         "ok": False})

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")

    print(f"corpus (parsed by both):       {n}")
    print(f"chematic SMILES RDKit-parseable: {n - unparseable}/{n} "
          f"({100 * (n - unparseable) / max(n, 1):.2f}%)")
    print(f"round-trip MATCH:              {rt_match}/{n} "
          f"({100 * rt_match / max(n, 1):.2f}%)")
    print(f"round-trip MISMATCH:           {rt_mismatch} "
          f"(of which E/Z-related: {rt_mismatch_ez})")
    print(f"idempotent:                    {idem_ok}/{n} "
          f"({100 * idem_ok / max(n, 1):.2f}%)")
    print(f"idempotency fails:             {idem_fail} "
          f"(of which E/Z-related: {idem_fail_ez})")
    print(f"\nwrote {len(rows)} divergence rows to {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
