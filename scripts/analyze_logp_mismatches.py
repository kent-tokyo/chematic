#!/usr/bin/env python3
"""Analyze LogP (Crippen) mismatches between chematic and RDKit on a SMILES corpus.

Usage:
    python scripts/analyze_logp_mismatches.py ~/Downloads/SMILES.csv [--tolerance 0.01] [--limit 5000]

Output:
    - Summary table to stdout
    - scripts/logp_mismatches.tsv  (tab-separated: smiles, rd_logp, ch_logp, delta)
"""
import argparse
import csv
import sys
from collections import Counter

try:
    from rdkit import Chem
    from rdkit.Chem import Crippen
except ImportError:
    sys.exit("rdkit-pypi required: pip install rdkit-pypi")

try:
    import chematic
except ImportError:
    sys.exit("chematic required: pip install chematic")

OUT_TSV = "scripts/logp_mismatches.tsv"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("corpus", help="CSV/TSV with a SMILES column")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="Absolute delta threshold (default 0.01)")
    parser.add_argument("--limit", type=int, default=5000,
                        help="Max molecules to process (default 5000)")
    args = parser.parse_args()

    mismatches = []
    total = 0
    skip = 0

    with open(args.corpus) as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i >= args.limit:
                break
            smi = row.get("smiles") or row.get("SMILES") or list(row.values())[0]
            smi = smi.strip()
            if not smi:
                continue

            rd_mol = Chem.MolFromSmiles(smi)
            if rd_mol is None:
                skip += 1
                continue

            try:
                ch_mol = chematic.from_smiles(smi)
                ch_logp = ch_mol.logp
            except Exception:
                skip += 1
                continue

            rd_logp = Crippen.MolLogP(rd_mol)
            delta = ch_logp - rd_logp
            total += 1

            if abs(delta) > args.tolerance:
                mismatches.append({
                    "smiles": smi,
                    "rd_logp": round(rd_logp, 4),
                    "ch_logp": round(ch_logp, 4),
                    "delta": round(delta, 4),
                })

    pct = len(mismatches) / total * 100 if total else 0
    print(f"Processed: {total}  Skipped: {skip}  Mismatches (|Δ|>{args.tolerance}): {len(mismatches)} ({pct:.1f}%)")

    if not mismatches:
        print("No mismatches — LogP is 100% within tolerance.")
        return

    # Sort by |delta| descending
    mismatches.sort(key=lambda r: abs(r["delta"]), reverse=True)

    # Bucket summary
    buckets = Counter()
    for r in mismatches:
        d = abs(r["delta"])
        if d > 1.0:
            buckets[">1.0"] += 1
        elif d > 0.5:
            buckets["0.5–1.0"] += 1
        elif d > 0.1:
            buckets["0.1–0.5"] += 1
        else:
            buckets[f"0.01–0.1"] += 1

    print("\nDelta magnitude buckets:")
    for label, count in sorted(buckets.items(), reverse=True):
        print(f"  {label:12s}: {count}")

    print(f"\nTop 10 mismatches (|Δ| largest first):")
    for r in mismatches[:10]:
        sign = "+" if r["delta"] > 0 else ""
        print(f"  Δ={sign}{r['delta']:+.3f}  rd={r['rd_logp']:.3f}  ch={r['ch_logp']:.3f}  {r['smiles'][:80]}")

    # Write TSV
    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["smiles", "rd_logp", "ch_logp", "delta"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(mismatches)
    print(f"\nSaved {len(mismatches)} mismatches → {OUT_TSV}")


if __name__ == "__main__":
    main()
