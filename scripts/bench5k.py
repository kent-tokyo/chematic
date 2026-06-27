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
    parser.add_argument("--json", metavar="PATH",
                        help="Write results as JSON to PATH (for validation dashboard)")
    args = parser.parse_args()

    # --- load libraries ---
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors, Crippen
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

    # hbd
    hbd_match = 0
    hbd_over = 0
    hbd_under = 0
    hbd_detail_count = 0

    # aromatic ring count
    arc_match = 0
    arc_detail_count = 0

    # [nH] SMARTS: whether the molecule has ANY [nH] atom
    nh_tp = 0  # both say yes
    nh_tn = 0  # both say no
    nh_fp = 0  # chematic yes, rdkit no
    nh_fn = 0  # chematic no, rdkit yes

    # tpsa (tolerance ±0.1 Å²)
    tpsa_match = 0
    tpsa_over = 0   # chematic > rdkit
    tpsa_under = 0  # chematic < rdkit
    tpsa_detail_count = 0

    # logp (tolerance ±0.01)
    logp_match = 0
    logp_over = 0
    logp_under = 0
    logp_detail_count = 0

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
        rd_hbd = rdMolDescriptors.CalcNumHBD(rd_mol)
        rd_arc = sum(
            1 for ring in rd_mol.GetRingInfo().AtomRings()
            if all(rd_mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
        )
        rd_has_nh = rd_mol.HasSubstructMatch(rd_nh_query)
        rd_tpsa = rdMolDescriptors.CalcTPSA(rd_mol, includeSandP=True)
        rd_logp = Crippen.MolLogP(rd_mol)

        # --- chematic ---
        try:
            ch_mol = chematic.from_smiles(smi)
        except Exception:
            parse_fail_ch += 1
            continue

        ch_hba = ch_mol.hba
        ch_hbd = ch_mol.hbd
        ch_arc = ch_mol.aromatic_ring_count
        ch_has_nh = chematic.smarts_match("[nH]", ch_mol)
        ch_tpsa = ch_mol.tpsa
        ch_logp = ch_mol.logp

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

        if rd_hbd == ch_hbd:
            hbd_match += 1
        else:
            delta = ch_hbd - rd_hbd
            if delta > 0:
                hbd_over += 1
            else:
                hbd_under += 1
            if args.detail and (args.limit is None or hbd_detail_count < args.limit):
                print(f"HBD {delta:+d}  rd={rd_hbd} ch={ch_hbd}  {smi}", file=sys.stderr)
                hbd_detail_count += 1

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

        tpsa_delta = ch_tpsa - rd_tpsa
        if abs(tpsa_delta) <= 0.1:
            tpsa_match += 1
        elif tpsa_delta > 0:
            tpsa_over += 1
            if args.detail and (args.limit is None or tpsa_detail_count < args.limit):
                print(f"TPSA +{tpsa_delta:.2f}  rd={rd_tpsa:.2f} ch={ch_tpsa:.2f}  {smi}", file=sys.stderr)
                tpsa_detail_count += 1
        else:
            tpsa_under += 1
            if args.detail and (args.limit is None or tpsa_detail_count < args.limit):
                print(f"TPSA {tpsa_delta:.2f}  rd={rd_tpsa:.2f} ch={ch_tpsa:.2f}  {smi}", file=sys.stderr)
                tpsa_detail_count += 1

        logp_delta = ch_logp - rd_logp
        if abs(logp_delta) <= 0.01:
            logp_match += 1
        elif logp_delta > 0:
            logp_over += 1
            if args.detail and (args.limit is None or logp_detail_count < args.limit):
                print(f"LogP +{logp_delta:.4f}  rd={rd_logp:.4f} ch={ch_logp:.4f}  {smi}", file=sys.stderr)
                logp_detail_count += 1
        else:
            logp_under += 1
            if args.detail and (args.limit is None or logp_detail_count < args.limit):
                print(f"LogP {logp_delta:.4f}  rd={rd_logp:.4f} ch={ch_logp:.4f}  {smi}", file=sys.stderr)
                logp_detail_count += 1

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
    print(f"  HBD agreement:         {hbd_match/total*100:6.1f}%  ({hbd_match}/{total})")
    print(f"    over-count (ch>rd):  {hbd_over:>6}  ({hbd_over/total*100:.1f}%)")
    print(f"    under-count(ch<rd):  {hbd_under:>6}  ({hbd_under/total*100:.1f}%)")
    print(f"  Aromatic ring count:   {arc_match/total*100:6.1f}%  ({arc_match}/{total})")
    nh_denom = nh_tp + nh_tn + nh_fp + nh_fn
    nh_agree = (nh_tp + nh_tn) / nh_denom * 100 if nh_denom else 0
    nh_prec = nh_tp / (nh_tp + nh_fp) * 100 if (nh_tp + nh_fp) else 0
    nh_rec  = nh_tp / (nh_tp + nh_fn) * 100 if (nh_tp + nh_fn) else 0
    print(f"  [nH] SMARTS overall:   {nh_agree:6.1f}%")
    print(f"    precision (no false positives): {nh_prec:.1f}%")
    print(f"    recall    (no false negatives): {nh_rec:.1f}%")
    print(f"    TP={nh_tp}  TN={nh_tn}  FP={nh_fp}  FN={nh_fn}")
    tpsa_miss = total - tpsa_match
    print(f"  TPSA (±0.1 Å²):        {tpsa_match/total*100:6.1f}%  ({tpsa_match}/{total})")
    print(f"    over  (ch>rd):       {tpsa_over:>6}  ({tpsa_over/total*100:.1f}%)")
    print(f"    under (ch<rd):       {tpsa_under:>6}  ({tpsa_under/total*100:.1f}%)")
    logp_miss = total - logp_match
    print(f"  LogP (±0.01):          {logp_match/total*100:6.1f}%  ({logp_match}/{total})")
    print(f"    over  (ch>rd):       {logp_over:>6}  ({logp_over/total*100:.1f}%)")
    print(f"    under (ch<rd):       {logp_under:>6}  ({logp_under/total*100:.1f}%)")
    print(f"{'='*55}")

    if args.json:
        import json, datetime, subprocess
        try:
            ver = subprocess.check_output(
                ["python3", "-c", "import chematic; print(chematic.__version__)"],
                text=True
            ).strip()
        except Exception:
            ver = "unknown"
        results = {
            "generated_at": datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            "chematic_version": ver,
            "corpus": {"total": total, "rdkit_parse_failures": parse_fail_rd,
                       "chematic_parse_failures": parse_fail_ch},
            "metrics": {
                "hba":  {"agreement_pct": round(hba_match/total*100, 2),  "match": hba_match,  "over": hba_over,  "under": hba_under,  "tolerance": "exact"},
                "hbd":  {"agreement_pct": round(hbd_match/total*100, 2),  "match": hbd_match,  "over": hbd_over,  "under": hbd_under,  "tolerance": "exact"},
                "arc":  {"agreement_pct": round(arc_match/total*100, 2),  "match": arc_match,  "tolerance": "exact"},
                "nh_smarts": {"agreement_pct": round(nh_agree, 2), "precision_pct": round(nh_prec, 2),
                              "recall_pct": round(nh_rec, 2), "tp": nh_tp, "tn": nh_tn, "fp": nh_fp, "fn": nh_fn},
                "tpsa": {"agreement_pct": round(tpsa_match/total*100, 2), "match": tpsa_match, "over": tpsa_over, "under": tpsa_under, "tolerance": "±0.1 Å²"},
                "logp": {"agreement_pct": round(logp_match/total*100, 2), "match": logp_match, "over": logp_over, "under": logp_under, "tolerance": "±0.01"},
            },
        }
        with open(args.json, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nJSON results written to {args.json}")

if __name__ == "__main__":
    main()
