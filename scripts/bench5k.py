#!/usr/bin/env python3
"""
Benchmark chematic hba_count + aromatic_ring_count + [nH] SMARTS
against the 5,000-molecule SMILES corpus, using RDKit as ground truth.

Usage:
    python3 scripts/bench5k.py ~/Downloads/SMILES.csv

Requires:  pip install rdkit
"""

import sys
import csv
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("smiles_csv", help="CSV with a 'SMILES' column (or first column)")
    parser.add_argument("--detail", action="store_true",
                        help="Print every mismatching molecule to stderr")
    parser.add_argument("--limit", type=int, default=None,
                        help="Only show first N mismatches per category in --detail mode")
    args = parser.parse_args()

    # --- load libraries ---
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        from rdkit.Chem import AllChem
    except ImportError:
        sys.exit("rdkit not installed. pip install rdkit")

    try:
        import chematic
    except ImportError:
        sys.exit("chematic not installed.")

    # --- read SMILES ---
    smiles_list = []
    with open(args.smiles_csv) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        col = "SMILES" if "SMILES" in fieldnames else fieldnames[0]
        for row in reader:
            smiles_list.append(row[col].strip())

    print(f"Loaded {len(smiles_list)} SMILES from {args.smiles_csv}", flush=True)

    # --- counters ---
    total = 0
    parse_fail_ch = 0
    parse_fail_rd = 0

    # hba
    hba_match = 0
    hba_over = 0   # chematic > rdkit
    hba_under = 0  # chematic < rdkit
    hba_detail_count = 0

    # aromatic ring count
    arc_match = 0
    arc_detail_count = 0

    # [nH] SMARTS: whether the molecule has ANY [nH] atom
    nh_tp = 0  # both say yes
    nh_tn = 0  # both say no
    nh_fp = 0  # chematic yes, rdkit no
    nh_fn = 0  # chematic no, rdkit yes

    # [nH] SMARTS query objects
    rd_nh_query = Chem.MolFromSmarts("[nH]")

    for smi in smiles_list:
        if not smi:
            continue

        # --- RDKit ---
        rd_mol = Chem.MolFromSmiles(smi)
        if rd_mol is None:
            parse_fail_rd += 1
            continue

        rd_hba = rdMolDescriptors.CalcNumHBA(rd_mol)
        rd_arc = sum(
            1 for ring in rd_mol.GetRingInfo().AtomRings()
            if all(rd_mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
        )
        rd_has_nh = rd_mol.HasSubstructMatch(rd_nh_query)

        # --- chematic ---
        try:
            ch_mol = chematic.from_smiles(smi)
        except Exception:
            parse_fail_ch += 1
            continue

        ch_hba = ch_mol.hba
        ch_arc = ch_mol.aromatic_ring_count
        ch_has_nh = chematic.smarts_match("[nH]", ch_mol)

        total += 1

        if rd_hba == ch_hba:
            hba_match += 1
        else:
            delta = ch_hba - rd_hba
            if delta > 0:
                hba_over += 1
            else:
                hba_under += 1
            if args.detail and (args.limit is None or hba_detail_count < args.limit):
                print(f"HBA {delta:+d}  rd={rd_hba} ch={ch_hba}  {smi}", file=sys.stderr)
                hba_detail_count += 1

        if rd_arc == ch_arc:
            arc_match += 1
        else:
            if args.detail and (args.limit is None or arc_detail_count < args.limit):
                print(f"ARC  rd={rd_arc} ch={ch_arc}  {smi}", file=sys.stderr)
                arc_detail_count += 1

        if rd_has_nh and ch_has_nh:
            nh_tp += 1
        elif not rd_has_nh and not ch_has_nh:
            nh_tn += 1
        elif ch_has_nh and not rd_has_nh:
            nh_fp += 1
        else:
            nh_fn += 1

    # --- report ---
    print(f"\n{'='*55}")
    print(f"  Molecules evaluated:   {total:>6}")
    print(f"  RDKit parse failures:  {parse_fail_rd:>6}")
    print(f"  chematic parse fails:  {parse_fail_ch:>6}")
    print(f"{'='*55}")
    hba_miss = total - hba_match
    print(f"  HBA agreement:         {hba_match/total*100:6.1f}%  ({hba_match}/{total})")
    print(f"    over-count (ch>rd):  {hba_over:>6}  ({hba_over/total*100:.1f}%)")
    print(f"    under-count(ch<rd):  {hba_under:>6}  ({hba_under/total*100:.1f}%)")
    print(f"  Aromatic ring count:   {arc_match/total*100:6.1f}%  ({arc_match}/{total})")
    nh_denom = nh_tp + nh_tn + nh_fp + nh_fn
    nh_agree = (nh_tp + nh_tn) / nh_denom * 100 if nh_denom else 0
    nh_prec = nh_tp / (nh_tp + nh_fp) * 100 if (nh_tp + nh_fp) else 0
    nh_rec  = nh_tp / (nh_tp + nh_fn) * 100 if (nh_tp + nh_fn) else 0
    print(f"  [nH] SMARTS overall:   {nh_agree:6.1f}%")
    print(f"    precision (no false positives): {nh_prec:.1f}%")
    print(f"    recall    (no false negatives): {nh_rec:.1f}%")
    print(f"    TP={nh_tp}  TN={nh_tn}  FP={nh_fp}  FN={nh_fn}")
    print(f"{'='*55}")

if __name__ == "__main__":
    main()
