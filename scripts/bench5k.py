#!/usr/bin/env python3
"""
Benchmark chematic descriptors against the 5,000-molecule SMILES corpus,
using RDKit as ground truth.

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
        from rdkit.Chem import rdMolDescriptors, Crippen, Lipinski
        from rdkit.Chem import AllChem
    except ImportError:
        sys.exit("rdkit not installed. pip install rdkit")

    # FindPotentialStereo (new CIP perception) — available RDKit 2022+.
    _find_stereo = getattr(Chem, "FindPotentialStereo", None)

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

    def make_int_counter():
        return {"match": 0, "over": 0, "under": 0, "detail_count": 0}

    def make_float_counter():
        return {"match": 0, "over": 0, "under": 0, "detail_count": 0}

    hba   = make_int_counter()
    hbd   = make_int_counter()
    arc   = make_int_counter()
    tpsa  = make_float_counter()
    logp  = make_float_counter()
    mr    = make_float_counter()
    fsp3  = make_float_counter()
    rb    = make_int_counter()   # rotatable bonds
    hac   = make_int_counter()   # heavy atom count
    nhet  = make_int_counter()   # num heteroatoms
    nsc     = make_int_counter()   # num stereocenters (legacy CalcNumAtomStereoCenters)
    nsc_new = make_int_counter()   # num stereocenters (new RDKit CIP / FindPotentialStereo)
    nahe  = make_int_counter()   # num aromatic heterocycles
    nalhe = make_int_counter()   # num aliphatic heterocycles
    nsat  = make_int_counter()   # num saturated heterocycles
    nalr  = make_int_counter()   # num aliphatic rings
    nsatr = make_int_counter()   # num saturated rings
    nspi  = make_int_counter()   # num spiro atoms
    nbr   = make_int_counter()   # num bridgehead atoms
    namide= make_int_counter()   # num amide bonds

    # [nH] SMARTS
    nh_tp = 0; nh_tn = 0; nh_fp = 0; nh_fn = 0

    # [nH] SMARTS query
    rd_nh_query = Chem.MolFromSmarts("[nH]")

    for smi in smiles_list:
        if not smi:
            continue

        # --- RDKit ---
        rd_mol = Chem.MolFromSmiles(smi)
        if rd_mol is None:
            parse_fail_rd += 1
            continue

        rd_hba  = rdMolDescriptors.CalcNumHBA(rd_mol)
        rd_hbd  = rdMolDescriptors.CalcNumHBD(rd_mol)
        rd_arc  = sum(
            1 for ring in rd_mol.GetRingInfo().AtomRings()
            if all(rd_mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
        )
        rd_has_nh = rd_mol.HasSubstructMatch(rd_nh_query)
        rd_tpsa = rdMolDescriptors.CalcTPSA(rd_mol, includeSandP=True)
        rd_logp = Crippen.MolLogP(rd_mol)
        rd_mr   = Crippen.MolMR(rd_mol)
        rd_fsp3 = rdMolDescriptors.CalcFractionCSP3(rd_mol)
        rd_rb   = rdMolDescriptors.CalcNumRotatableBonds(rd_mol)
        rd_hac  = rd_mol.GetNumHeavyAtoms()
        rd_nhet = Lipinski.NumHeteroatoms(rd_mol)
        rd_nsc  = rdMolDescriptors.CalcNumAtomStereoCenters(rd_mol)
        if _find_stereo is not None:
            rd_nsc_new = sum(
                1 for s in _find_stereo(rd_mol)
                if str(s.type).endswith("Atom_Tetrahedral")
            )
        else:
            rd_nsc_new = rd_nsc  # fallback for older RDKit
        rd_nahe = rdMolDescriptors.CalcNumAromaticHeterocycles(rd_mol)
        rd_nalhe= rdMolDescriptors.CalcNumAliphaticHeterocycles(rd_mol)
        rd_nsat = rdMolDescriptors.CalcNumSaturatedHeterocycles(rd_mol)
        rd_nalr = rdMolDescriptors.CalcNumAliphaticRings(rd_mol)
        rd_nsatr= rdMolDescriptors.CalcNumSaturatedRings(rd_mol)
        rd_nspi = rdMolDescriptors.CalcNumSpiroAtoms(rd_mol)
        rd_nbr  = rdMolDescriptors.CalcNumBridgeheadAtoms(rd_mol)
        rd_namide=rdMolDescriptors.CalcNumAmideBonds(rd_mol)

        # --- chematic ---
        try:
            ch_mol = chematic.from_smiles(smi)
        except Exception:
            parse_fail_ch += 1
            continue

        ch_hba  = ch_mol.hba
        ch_hbd  = ch_mol.hbd
        ch_arc  = ch_mol.aromatic_ring_count
        ch_has_nh = chematic.smarts_match("[nH]", ch_mol)
        ch_tpsa = ch_mol.tpsa
        ch_logp = ch_mol.logp
        ch_mr   = ch_mol.molar_refractivity
        ch_fsp3 = ch_mol.fsp3
        ch_rb   = ch_mol.rotatable_bonds
        ch_hac  = ch_mol.heavy_atoms
        ch_nhet = ch_mol.num_heteroatoms
        ch_nsc  = ch_mol.num_stereocenters
        ch_nahe = ch_mol.num_aromatic_heterocycles
        ch_nalhe= ch_mol.num_aliphatic_heterocycles
        ch_nsat = ch_mol.num_saturated_heterocycles
        ch_nalr = ch_mol.num_aliphatic_rings
        ch_nsatr= ch_mol.num_saturated_rings
        ch_nspi = ch_mol.num_spiro_atoms
        ch_nbr  = ch_mol.num_bridgehead_atoms
        ch_namide=ch_mol.num_amide_bonds

        total += 1

        # --- integer metric helper ---
        def check_int(c, rd_val, ch_val, label, smi):
            if rd_val == ch_val:
                c["match"] += 1
            else:
                d = ch_val - rd_val
                if d > 0: c["over"] += 1
                else:     c["under"] += 1
                if args.detail and (args.limit is None or c["detail_count"] < args.limit):
                    print(f"{label} {d:+d}  rd={rd_val} ch={ch_val}  {smi}", file=sys.stderr)
                    c["detail_count"] += 1

        def check_float(c, rd_val, ch_val, tol, label, smi, fmt=".2f"):
            delta = ch_val - rd_val
            if abs(delta) <= tol:
                c["match"] += 1
            elif delta > 0:
                c["over"] += 1
                if args.detail and (args.limit is None or c["detail_count"] < args.limit):
                    print(f"{label} +{delta:{fmt}}  rd={rd_val:{fmt}} ch={ch_val:{fmt}}  {smi}", file=sys.stderr)
                    c["detail_count"] += 1
            else:
                c["under"] += 1
                if args.detail and (args.limit is None or c["detail_count"] < args.limit):
                    print(f"{label} {delta:{fmt}}  rd={rd_val:{fmt}} ch={ch_val:{fmt}}  {smi}", file=sys.stderr)
                    c["detail_count"] += 1

        check_int(hba,    rd_hba,   ch_hba,   "HBA",    smi)
        check_int(hbd,    rd_hbd,   ch_hbd,   "HBD",    smi)
        check_int(arc,    rd_arc,   ch_arc,   "ARC",    smi)
        check_float(tpsa, rd_tpsa,  ch_tpsa,  0.1,  "TPSA",  smi, ".2f")
        check_float(logp, rd_logp,  ch_logp,  0.01, "LogP",  smi, ".4f")
        check_float(mr,   rd_mr,    ch_mr,    0.01, "MR",    smi, ".2f")
        check_float(fsp3, rd_fsp3,  ch_fsp3,  0.001,"Fsp3",  smi, ".4f")
        check_int(rb,     rd_rb,    ch_rb,    "RotB",  smi)
        check_int(hac,    rd_hac,   ch_hac,   "HAC",   smi)
        check_int(nhet,   rd_nhet,  ch_nhet,  "NHet",  smi)
        check_int(nsc,     rd_nsc,     ch_nsc, "NSC",     smi)
        check_int(nsc_new, rd_nsc_new, ch_nsc, "NSC_new", smi)
        check_int(nahe,   rd_nahe,  ch_nahe,  "NAHet", smi)
        check_int(nalhe,  rd_nalhe, ch_nalhe, "NALHet",smi)
        check_int(nsat,   rd_nsat,  ch_nsat,  "NSatHet",smi)
        check_int(nalr,   rd_nalr,  ch_nalr,  "NALRing",smi)
        check_int(nsatr,  rd_nsatr, ch_nsatr, "NSatRing",smi)
        check_int(nspi,   rd_nspi,  ch_nspi,  "NSpiro",smi)
        check_int(nbr,    rd_nbr,   ch_nbr,   "NBridge",smi)
        check_int(namide, rd_namide,ch_namide,"NAmide",smi)

        if rd_has_nh and ch_has_nh:       nh_tp += 1
        elif not rd_has_nh and not ch_has_nh: nh_tn += 1
        elif ch_has_nh and not rd_has_nh: nh_fp += 1
        else:                             nh_fn += 1

    # --- report ---
    def pct(c): return c["match"] / total * 100
    def fmt_int(label, c, width=28):
        s = f"  {label:{width}} {pct(c):6.1f}%  ({c['match']}/{total})"
        if c["over"] or c["under"]:
            s += f"\n    over  (ch>rd):  {c['over']:>6}  ({c['over']/total*100:.1f}%)"
            s += f"\n    under (ch<rd):  {c['under']:>6}  ({c['under']/total*100:.1f}%)"
        return s

    sep = "=" * 57
    print(f"\n{sep}")
    print(f"  Molecules evaluated:     {total:>6}")
    print(f"  RDKit parse failures:    {parse_fail_rd:>6}")
    print(f"  chematic parse fails:    {parse_fail_ch:>6}")
    print(sep)

    print(fmt_int("HBA agreement:", hba))
    print(fmt_int("HBD agreement:", hbd))
    print(f"  {'Aromatic ring count:':<28} {pct(arc):6.1f}%  ({arc['match']}/{total})")

    nh_denom = nh_tp + nh_tn + nh_fp + nh_fn
    nh_agree = (nh_tp + nh_tn) / nh_denom * 100 if nh_denom else 0
    nh_prec  = nh_tp / (nh_tp + nh_fp) * 100 if (nh_tp + nh_fp) else 0
    nh_rec   = nh_tp / (nh_tp + nh_fn) * 100 if (nh_tp + nh_fn) else 0
    print(f"  {'[nH] SMARTS overall:':<28} {nh_agree:6.1f}%")
    print(f"    precision (no FP): {nh_prec:.1f}%   recall (no FN): {nh_rec:.1f}%")
    print(f"    TP={nh_tp}  TN={nh_tn}  FP={nh_fp}  FN={nh_fn}")

    print(fmt_int("TPSA (±0.1 Å²):", tpsa))
    print(fmt_int("LogP (±0.01):", logp))
    print(fmt_int("MR   (±0.01):", mr))
    print(fmt_int("Fsp3 (±0.001):", fsp3))
    print(fmt_int("Rotatable bonds:", rb))
    print(fmt_int("Heavy atom count:", hac))
    print(fmt_int("Num heteroatoms:", nhet))
    print(fmt_int("Num stereocenters (legacy):", nsc))
    print(fmt_int("Num stereocenters (new CIP):", nsc_new))
    print(fmt_int("Num arom heterocycles:", nahe))
    print(fmt_int("Num aliph heterocycles:", nalhe))
    print(fmt_int("Num sat  heterocycles:", nsat))
    print(fmt_int("Num aliphatic rings:", nalr))
    print(fmt_int("Num saturated rings:", nsatr))
    print(fmt_int("Num spiro atoms:", nspi))
    print(fmt_int("Num bridgehead atoms:", nbr))
    print(fmt_int("Num amide bonds:", namide))

    print(sep)

    if args.json:
        import json, datetime, subprocess
        try:
            ver = subprocess.check_output(
                ["python3", "-c", "import chematic; print(chematic.__version__)"],
                text=True
            ).strip()
        except Exception:
            ver = "unknown"
        def metric_dict(c, tol=None):
            d = {"agreement_pct": round(pct(c), 2), "match": c["match"],
                 "over": c["over"], "under": c["under"]}
            if tol is not None:
                d["tolerance"] = tol
            return d
        results = {
            "generated_at": datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            "chematic_version": ver,
            "corpus": {"total": total, "rdkit_parse_failures": parse_fail_rd,
                       "chematic_parse_failures": parse_fail_ch},
            "metrics": {
                "hba":   metric_dict(hba,  "exact"),
                "hbd":   metric_dict(hbd,  "exact"),
                "arc":   metric_dict(arc,  "exact"),
                "nh_smarts": {"agreement_pct": round(nh_agree, 2),
                              "precision_pct": round(nh_prec, 2),
                              "recall_pct": round(nh_rec, 2),
                              "tp": nh_tp, "tn": nh_tn, "fp": nh_fp, "fn": nh_fn},
                "tpsa":  metric_dict(tpsa, "±0.1 Å²"),
                "logp":  metric_dict(logp, "±0.01"),
                "mr":    metric_dict(mr,   "±0.01"),
                "fsp3":  metric_dict(fsp3, "±0.001"),
                "rotatable_bonds":          metric_dict(rb,     "exact"),
                "heavy_atom_count":         metric_dict(hac,    "exact"),
                "num_heteroatoms":          metric_dict(nhet,   "exact"),
                "num_stereocenters":         metric_dict(nsc,     "exact"),
                "num_stereocenters_new_cip": metric_dict(nsc_new, "exact"),
                "num_aromatic_heterocycles":metric_dict(nahe,   "exact"),
                "num_aliphatic_heterocycles":metric_dict(nalhe, "exact"),
                "num_saturated_heterocycles":metric_dict(nsat,  "exact"),
                "num_aliphatic_rings":      metric_dict(nalr,   "exact"),
                "num_saturated_rings":      metric_dict(nsatr,  "exact"),
                "num_spiro_atoms":          metric_dict(nspi,   "exact"),
                "num_bridgehead_atoms":     metric_dict(nbr,    "exact"),
                "num_amide_bonds":          metric_dict(namide, "exact"),
            },
        }
        with open(args.json, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nJSON results written to {args.json}")

if __name__ == "__main__":
    main()
