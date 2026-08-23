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
        from rdkit.Chem import rdMolDescriptors, Crippen, Lipinski, Descriptors, EState
        from rdkit.Chem import AllChem, FilterCatalog
    except ImportError:
        sys.exit("rdkit not installed. pip install rdkit")

    # RDKit's SA_Score contrib script isn't pip-installable; load it from the
    # rdkit package's own Contrib/ directory (present in any full RDKit install).
    import os
    _sa_score_dir = os.path.join(os.path.dirname(__import__("rdkit").__file__), "Contrib", "SA_Score")
    sys.path.append(_sa_score_dir)
    try:
        import sascorer
    except ImportError:
        sascorer = None

    _pains_params = FilterCatalog.FilterCatalogParams()
    _pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    _pains_catalog = FilterCatalog.FilterCatalog(_pains_params)

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

    # per-molecule signed deltas (ch - rd), keyed by SMILES, for drift
    # tracking across runs (see scripts/veridict_accuracy_report.py) —
    # separate from the match/over/under aggregate counters above.
    tpsa_deltas = {}
    logp_deltas = {}

    def make_int_counter():
        return {"match": 0, "over": 0, "under": 0, "detail_count": 0}

    def make_float_counter():
        return {"match": 0, "over": 0, "under": 0, "detail_count": 0}

    mw    = make_float_counter()
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
    nsc_consensus = make_int_counter()  # where legacy == new CIP == chematic
    logp_max_delta = 0.0
    nsc_oracle_disagree = 0        # molecules where legacy != new CIP
    nsc_oracle_legacy_under = 0    # legacy < new CIP
    nsc_oracle_new_over = 0        # new CIP > legacy (same as above, broken out)
    nahe  = make_int_counter()   # num aromatic heterocycles
    nalhe = make_int_counter()   # num aliphatic heterocycles
    nsat  = make_int_counter()   # num saturated heterocycles
    nalr  = make_int_counter()   # num aliphatic rings
    nsatr = make_int_counter()   # num saturated rings
    nspi  = make_int_counter()   # num spiro atoms
    nbr   = make_int_counter()   # num bridgehead atoms
    namide= make_int_counter()   # num amide bonds

    # --- easy wins: spot-checked exact matches, no bulk script existed yet ---
    formal_charge = make_int_counter()
    num_valence_e = make_int_counter()
    exact_mass    = make_float_counter()
    hba_lipinski  = make_int_counter()
    nusc          = make_int_counter()   # num unspecified stereocenters
    ring_count    = make_int_counter()
    chi0  = make_float_counter(); chi1  = make_float_counter()
    chi0v = make_float_counter(); chi1v = make_float_counter()
    chi2v = make_float_counter(); chi3v = make_float_counter(); chi4v = make_float_counter()

    # --- named 1:1 RDKit equivalents flagged as divergent by spot-check;
    #     bulk run here quantifies the real corpus-wide agreement rate ---
    kappa1 = make_float_counter(); kappa2 = make_float_counter(); kappa3 = make_float_counter()
    hall_kier_alpha = make_float_counter()
    bertz_ct  = make_float_counter()
    balaban_j = make_float_counter()
    labute_asa = make_float_counter()
    sa_score  = make_float_counter()
    qed       = make_float_counter()
    max_estate = make_float_counter(); min_estate = make_float_counter(); sum_estate = make_float_counter()

    # --- array families (BCUT2D 8, MQN 42, VSA 47) ---
    BCUT2D_NAMES = ["bcut2d_mwhi", "bcut2d_mwlo", "bcut2d_chghi", "bcut2d_chglo",
                    "bcut2d_logphi", "bcut2d_logplo", "bcut2d_mrhi", "bcut2d_mrlo"]
    MQN_NAMES = [f"MQN{i}" for i in range(1, 43)]
    SLOGP_VSA_NAMES = [f"SlogP_VSA{i}" for i in range(1, 13)]
    SMR_VSA_NAMES   = [f"SMR_VSA{i}"   for i in range(1, 11)]
    PEOE_VSA_NAMES  = [f"PEOE_VSA{i}"  for i in range(1, 15)]
    ESTATE_VSA_NAMES = [f"EState_VSA{i}" for i in range(1, 12)]

    def make_family_counters(names):
        return {n: make_float_counter() for n in names}

    bcut2d = make_family_counters(BCUT2D_NAMES)
    mqn    = make_family_counters(MQN_NAMES)
    slogp_vsa  = make_family_counters(SLOGP_VSA_NAMES)
    smr_vsa    = make_family_counters(SMR_VSA_NAMES)
    peoe_vsa   = make_family_counters(PEOE_VSA_NAMES)
    estate_vsa = make_family_counters(ESTATE_VSA_NAMES)

    # PAINS boolean agreement (corpus-wide, vs the 5-molecule spot-check done in review)
    pains_tp = pains_tn = pains_fp = pains_fn = 0
    extended_skipped = 0  # molecules where RDKit couldn't compute the extended metrics (e.g. rare elements)

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

        rd_mw   = Descriptors.MolWt(rd_mol)
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

        # BCUT2D/Kappa/LabuteASA etc. rely on Gasteiger charges + valence tables
        # that RDKit doesn't have parameters for on some rare elements (e.g. Te) —
        # skip the extended-metric block for just those molecules rather than
        # aborting the whole run.
        rd_extended_ok = True
        try:
            rd_formal_charge = Chem.GetFormalCharge(rd_mol)
            rd_num_valence_e = Descriptors.NumValenceElectrons(rd_mol)
            rd_exact_mass    = Descriptors.ExactMolWt(rd_mol)
            rd_hba_lipinski  = rdMolDescriptors.CalcNumLipinskiHBA(rd_mol)
            rd_nusc          = rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(rd_mol)
            rd_ring_count    = rdMolDescriptors.CalcNumRings(rd_mol)
            rd_chi0, rd_chi1 = Descriptors.Chi0(rd_mol), Descriptors.Chi1(rd_mol)
            rd_chi0v = rdMolDescriptors.CalcChi0v(rd_mol)
            rd_chi1v = rdMolDescriptors.CalcChi1v(rd_mol)
            rd_chi2v = rdMolDescriptors.CalcChi2v(rd_mol)
            rd_chi3v = rdMolDescriptors.CalcChi3v(rd_mol)
            rd_chi4v = rdMolDescriptors.CalcChi4v(rd_mol)

            rd_kappa1 = rdMolDescriptors.CalcKappa1(rd_mol)
            rd_kappa2 = rdMolDescriptors.CalcKappa2(rd_mol)
            rd_kappa3 = rdMolDescriptors.CalcKappa3(rd_mol)
            rd_hall_kier_alpha = rdMolDescriptors.CalcHallKierAlpha(rd_mol)
            rd_bertz_ct  = Descriptors.BertzCT(rd_mol)
            rd_balaban_j = Descriptors.BalabanJ(rd_mol)
            rd_labute_asa = rdMolDescriptors.CalcLabuteASA(rd_mol)
            rd_sa_score = sascorer.calculateScore(rd_mol) if sascorer else None
            rd_qed = Descriptors.qed(rd_mol)
            rd_max_estate = Descriptors.MaxEStateIndex(rd_mol)
            rd_min_estate = Descriptors.MinEStateIndex(rd_mol)
            rd_sum_estate = sum(EState.EStateIndices(rd_mol))

            rd_bcut2d_vals = rdMolDescriptors.BCUT2D(rd_mol)
            rd_mqn_vals    = rdMolDescriptors.MQNs_(rd_mol)
            rd_slogp_vsa = [getattr(Descriptors, f"SlogP_VSA{i}")(rd_mol) for i in range(1, 13)]
            rd_smr_vsa   = [getattr(Descriptors, f"SMR_VSA{i}")(rd_mol)   for i in range(1, 11)]
            rd_peoe_vsa  = [getattr(Descriptors, f"PEOE_VSA{i}")(rd_mol)  for i in range(1, 15)]
            rd_estate_vsa= [getattr(Descriptors, f"EState_VSA{i}")(rd_mol) for i in range(1, 12)]

            rd_pains_hit = _pains_catalog.HasMatch(rd_mol)
        except Exception:
            rd_extended_ok = False

        # --- chematic ---
        try:
            ch_mol = chematic.from_smiles(smi)
        except Exception:
            parse_fail_ch += 1
            continue

        ch_mw   = ch_mol.mw
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

        ch_desc = ch_mol.descriptors()

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

        check_float(mw,   rd_mw,    ch_mw,    0.01, "MW",     smi, ".4f")
        check_int(hba,    rd_hba,   ch_hba,   "HBA",    smi)
        check_int(hbd,    rd_hbd,   ch_hbd,   "HBD",    smi)
        check_int(arc,    rd_arc,   ch_arc,   "ARC",    smi)
        check_float(tpsa, rd_tpsa,  ch_tpsa,  0.1,  "TPSA",  smi, ".2f")
        check_float(logp, rd_logp,  ch_logp,  0.01, "LogP",  smi, ".4f")
        logp_max_delta = max(logp_max_delta, abs(ch_logp - rd_logp))
        tpsa_deltas[smi] = ch_tpsa - rd_tpsa
        logp_deltas[smi] = ch_logp - rd_logp
        check_float(mr,   rd_mr,    ch_mr,    0.01, "MR",    smi, ".2f")
        check_float(fsp3, rd_fsp3,  ch_fsp3,  0.001,"Fsp3",  smi, ".4f")
        check_int(rb,     rd_rb,    ch_rb,    "RotB",  smi)
        check_int(hac,    rd_hac,   ch_hac,   "HAC",   smi)
        check_int(nhet,   rd_nhet,  ch_nhet,  "NHet",  smi)
        check_int(nsc,     rd_nsc,     ch_nsc, "NSC",     smi)
        check_int(nsc_new, rd_nsc_new, ch_nsc, "NSC_new", smi)
        if rd_nsc == rd_nsc_new == ch_nsc:
            nsc_consensus["match"] += 1
        elif ch_nsc > max(rd_nsc, rd_nsc_new):
            nsc_consensus["over"] += 1
        else:
            nsc_consensus["under"] += 1
        if rd_nsc != rd_nsc_new:
            nsc_oracle_disagree += 1
            if rd_nsc < rd_nsc_new:
                nsc_oracle_legacy_under += 1
            else:
                nsc_oracle_new_over += 1
        check_int(nahe,   rd_nahe,  ch_nahe,  "NAHet", smi)
        check_int(nalhe,  rd_nalhe, ch_nalhe, "NALHet",smi)
        check_int(nsat,   rd_nsat,  ch_nsat,  "NSatHet",smi)
        check_int(nalr,   rd_nalr,  ch_nalr,  "NALRing",smi)
        check_int(nsatr,  rd_nsatr, ch_nsatr, "NSatRing",smi)
        check_int(nspi,   rd_nspi,  ch_nspi,  "NSpiro",smi)
        check_int(nbr,    rd_nbr,   ch_nbr,   "NBridge",smi)
        check_int(namide, rd_namide,ch_namide,"NAmide",smi)

        if rd_extended_ok:
            check_int(formal_charge, rd_formal_charge, ch_desc["formal_charge"], "FormalChg", smi)
            check_int(num_valence_e, rd_num_valence_e, ch_desc["num_valence_electrons"], "NValE", smi)
            check_float(exact_mass, rd_exact_mass, ch_desc["exact_mass"], 0.01, "ExactMass", smi, ".4f")
            check_int(hba_lipinski, rd_hba_lipinski, ch_desc["hba_lipinski"], "HBALip", smi)
            check_int(nusc, rd_nusc, ch_desc["num_unspecified_stereocenters"], "NUnspecSC", smi)
            check_int(ring_count, rd_ring_count, ch_desc["ring_count"], "RingCount", smi)
            check_float(chi0,  rd_chi0,  ch_desc["chi0"],  1e-6, "Chi0",  smi, ".6f")
            check_float(chi1,  rd_chi1,  ch_desc["chi1"],  1e-6, "Chi1",  smi, ".6f")
            check_float(chi0v, rd_chi0v, ch_desc["chi0v"], 1e-6, "Chi0v", smi, ".6f")
            check_float(chi1v, rd_chi1v, ch_desc["chi1v"], 1e-6, "Chi1v", smi, ".6f")
            check_float(chi2v, rd_chi2v, ch_desc["chi2v"], 1e-6, "Chi2v", smi, ".6f")
            check_float(chi3v, rd_chi3v, ch_desc["chi3v"], 1e-6, "Chi3v", smi, ".6f")
            check_float(chi4v, rd_chi4v, ch_desc["chi4v"], 1e-6, "Chi4v", smi, ".6f")

            check_float(kappa1, rd_kappa1, ch_desc["kappa1"], 0.01, "Kappa1", smi, ".3f")
            check_float(kappa2, rd_kappa2, ch_desc["kappa2"], 0.01, "Kappa2", smi, ".3f")
            check_float(kappa3, rd_kappa3, ch_desc["kappa3"], 0.01, "Kappa3", smi, ".3f")
            check_float(hall_kier_alpha, rd_hall_kier_alpha, ch_desc["hall_kier_alpha"], 0.01, "HallKierAlpha", smi, ".3f")
            check_float(bertz_ct,  rd_bertz_ct,  ch_desc["bertz_ct"],  0.01, "BertzCT",  smi, ".2f")
            check_float(balaban_j, rd_balaban_j, ch_desc["balaban_j"], 0.01, "BalabanJ", smi, ".3f")
            check_float(labute_asa, rd_labute_asa, ch_desc["labute_asa"], 0.01, "LabuteASA", smi, ".2f")
            if rd_sa_score is not None:
                check_float(sa_score, rd_sa_score, ch_desc["sa_score"], 0.01, "SAScore", smi, ".3f")
            check_float(qed, rd_qed, ch_desc["qed"], 0.01, "QED", smi, ".4f")
            check_float(max_estate, rd_max_estate, ch_desc["max_estate"], 0.01, "MaxEState", smi, ".3f")
            check_float(min_estate, rd_min_estate, ch_desc["min_estate"], 0.01, "MinEState", smi, ".3f")
            check_float(sum_estate, rd_sum_estate, ch_desc["sum_estate"], 0.01, "SumEState", smi, ".3f")

            for name, rv in zip(BCUT2D_NAMES, rd_bcut2d_vals):
                check_float(bcut2d[name], rv, ch_desc[name], 0.01, name, smi, ".3f")
            for name, rv in zip(MQN_NAMES, rd_mqn_vals):
                check_float(mqn[name], float(rv), float(ch_desc[name]), 0.5, name, smi, ".0f")
            for name, rv in zip(SLOGP_VSA_NAMES, rd_slogp_vsa):
                check_float(slogp_vsa[name], rv, ch_desc[name], 0.01, name, smi, ".3f")
            for name, rv in zip(SMR_VSA_NAMES, rd_smr_vsa):
                check_float(smr_vsa[name], rv, ch_desc[name], 0.01, name, smi, ".3f")
            for name, rv in zip(PEOE_VSA_NAMES, rd_peoe_vsa):
                check_float(peoe_vsa[name], rv, ch_desc[name], 0.01, name, smi, ".3f")
            for name, rv in zip(ESTATE_VSA_NAMES, rd_estate_vsa):
                check_float(estate_vsa[name], rv, ch_desc[name], 0.01, name, smi, ".3f")

            ch_pains_hit = not ch_desc["pains_passes"]
            if rd_pains_hit and ch_pains_hit:           pains_tp += 1
            elif not rd_pains_hit and not ch_pains_hit: pains_tn += 1
            elif ch_pains_hit and not rd_pains_hit:     pains_fp += 1
            else:                                       pains_fn += 1
        else:
            extended_skipped += 1

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
    print(f"  extended metrics skipped:{extended_skipped:>6}  (rare elements RDKit lacks Gasteiger/valence params for)")
    print(sep)

    print(fmt_int("Molecular weight (±0.01):", mw))
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
    print("  Easy wins (untested, spot-checked exact):")
    print(fmt_int("Formal charge:", formal_charge))
    print(fmt_int("Num valence electrons:", num_valence_e))
    print(fmt_int("Exact mass (±0.01):", exact_mass))
    print(fmt_int("HBA (Lipinski):", hba_lipinski))
    print(fmt_int("Num unspecified stereocenters:", nusc))
    print(fmt_int("Ring count:", ring_count))
    print(fmt_int("Chi0 (±1e-6):", chi0))
    print(fmt_int("Chi1 (±1e-6):", chi1))
    print(fmt_int("Chi0v-Chi4v (±1e-6):", chi0v))

    print(sep)
    print("  Named RDKit equivalents flagged divergent by spot-check (±0.01 unless noted):")
    print(fmt_int("Kappa1:", kappa1))
    print(fmt_int("Kappa2:", kappa2))
    print(fmt_int("Kappa3:", kappa3))
    print(fmt_int("HallKierAlpha:", hall_kier_alpha))
    print(fmt_int("BertzCT:", bertz_ct))
    print(fmt_int("BalabanJ:", balaban_j))
    print(fmt_int("LabuteASA:", labute_asa))
    if sascorer:
        print(fmt_int("SA Score:", sa_score))
    else:
        print("  SA Score:                    skipped (rdkit Contrib/SA_Score not found)")
    print(fmt_int("QED:", qed))
    print(fmt_int("Max EState:", max_estate))
    print(fmt_int("Min EState:", min_estate))
    print(fmt_int("Sum EState:", sum_estate))

    print(sep)
    print("  PAINS (boolean agreement, corpus-wide):")
    pains_denom = pains_tp + pains_tn + pains_fp + pains_fn
    pains_agree = (pains_tp + pains_tn) / pains_denom * 100 if pains_denom else 0
    print(f"    agreement: {pains_agree:6.1f}%   TP={pains_tp} TN={pains_tn} FP={pains_fp} FN={pains_fn}")

    def fmt_family(label, counters):
        fam_match = sum(c["match"] for c in counters.values())
        fam_total = sum(c["match"] + c["over"] + c["under"] for c in counters.values())
        pct_fam = fam_match / fam_total * 100 if fam_total else 0
        worst = min(counters.items(), key=lambda kv: kv[1]["match"])
        worst_pct = worst[1]["match"] / total * 100 if total else 0
        print(f"  {label:<28} {pct_fam:6.1f}%  ({fam_match}/{fam_total} molecule×sub-descriptor pairs)"
              f"   worst sub-descriptor: {worst[0]} ({worst_pct:.1f}%)")

    print(sep)
    print("  Array families (aggregated across all sub-descriptors):")
    fmt_family("BCUT2D (8 values):", bcut2d)
    fmt_family("MQN (42 values):", mqn)
    fmt_family("SlogP_VSA (12 values):", slogp_vsa)
    fmt_family("SMR_VSA (10 values):", smr_vsa)
    fmt_family("PEOE_VSA (14 values):", peoe_vsa)
    fmt_family("EState_VSA (11 values):", estate_vsa)

    print(sep)

    if args.json:
        import json, datetime
        try:
            ver = chematic.__version__
        except AttributeError:
            ver = "unknown"
        def metric_dict(c, tol=None):
            d = {"agreement_pct": round(pct(c), 2), "match": c["match"],
                 "over": c["over"], "under": c["under"]}
            if tol is not None:
                d["tolerance"] = tol
            return d
        def family_dict(counters):
            per_sub = {name: metric_dict(c, "±0.01") for name, c in counters.items()}
            fam_match = sum(c["match"] for c in counters.values())
            fam_total = sum(c["match"] + c["over"] + c["under"] for c in counters.values())
            return {
                "agreement_pct": round(fam_match / fam_total * 100, 2) if fam_total else 0,
                "per_sub_descriptor": per_sub,
            }
        import rdkit as _rdkit
        results = {
            "generated_at": datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            "chematic_version": ver,
            "rdkit_version": _rdkit.__version__,
            "corpus": {"total": total, "rdkit_parse_failures": parse_fail_rd,
                       "chematic_parse_failures": parse_fail_ch,
                       "extended_metrics_skipped": extended_skipped},
            "stereocenters": {
                "oracle_disagreements": nsc_oracle_disagree,
                "oracle_disagree_legacy_under": nsc_oracle_legacy_under,
                "oracle_disagree_new_cip_over": nsc_oracle_new_over,
            },
            "metrics": {
                "mw":    metric_dict(mw,   "±0.01"),
                "hba":   metric_dict(hba,  "exact"),
                "hbd":   metric_dict(hbd,  "exact"),
                "arc":   metric_dict(arc,  "exact"),
                "nh_smarts": {"agreement_pct": round(nh_agree, 2),
                              "precision_pct": round(nh_prec, 2),
                              "recall_pct": round(nh_rec, 2),
                              "tp": nh_tp, "tn": nh_tn, "fp": nh_fp, "fn": nh_fn},
                "tpsa":  metric_dict(tpsa, "±0.1 Å²"),
                "logp":  {**metric_dict(logp, "±0.01"), "max_delta": logp_max_delta},
                "mr":    metric_dict(mr,   "±0.01"),
                "fsp3":  metric_dict(fsp3, "±0.001"),
                "rotatable_bonds":           metric_dict(rb,            "exact"),
                "heavy_atom_count":          metric_dict(hac,           "exact"),
                "num_heteroatoms":           metric_dict(nhet,          "exact"),
                "num_stereocenters":         metric_dict(nsc,           "exact"),
                "num_stereocenters_new_cip": metric_dict(nsc_new,       "exact"),
                "num_stereocenters_consensus": metric_dict(nsc_consensus, "exact"),
                "num_aromatic_heterocycles":metric_dict(nahe,   "exact"),
                "num_aliphatic_heterocycles":metric_dict(nalhe, "exact"),
                "num_saturated_heterocycles":metric_dict(nsat,  "exact"),
                "num_aliphatic_rings":      metric_dict(nalr,   "exact"),
                "num_saturated_rings":      metric_dict(nsatr,  "exact"),
                "num_spiro_atoms":          metric_dict(nspi,   "exact"),
                "num_bridgehead_atoms":     metric_dict(nbr,    "exact"),
                "num_amide_bonds":          metric_dict(namide, "exact"),
                "formal_charge":            metric_dict(formal_charge, "exact"),
                "num_valence_electrons":    metric_dict(num_valence_e, "exact"),
                "exact_mass":               metric_dict(exact_mass,    "±0.01"),
                "hba_lipinski":             metric_dict(hba_lipinski,  "exact"),
                "num_unspecified_stereocenters": metric_dict(nusc,     "exact"),
                "ring_count":               metric_dict(ring_count,    "exact"),
                "chi0":  metric_dict(chi0,  "±1e-6"),
                "chi1":  metric_dict(chi1,  "±1e-6"),
                "chi0v": metric_dict(chi0v, "±1e-6"),
                "chi1v": metric_dict(chi1v, "±1e-6"),
                "chi2v": metric_dict(chi2v, "±1e-6"),
                "chi3v": metric_dict(chi3v, "±1e-6"),
                "chi4v": metric_dict(chi4v, "±1e-6"),
                "kappa1": metric_dict(kappa1, "±0.01"),
                "kappa2": metric_dict(kappa2, "±0.01"),
                "kappa3": metric_dict(kappa3, "±0.01"),
                "hall_kier_alpha": metric_dict(hall_kier_alpha, "±0.01"),
                "bertz_ct":  metric_dict(bertz_ct,  "±0.01"),
                "balaban_j": metric_dict(balaban_j, "±0.01"),
                "labute_asa": metric_dict(labute_asa, "±0.01"),
                "sa_score": metric_dict(sa_score, "±0.01") if sascorer else None,
                "qed": metric_dict(qed, "±0.01"),
                "max_estate": metric_dict(max_estate, "±0.01"),
                "min_estate": metric_dict(min_estate, "±0.01"),
                "sum_estate": metric_dict(sum_estate, "±0.01"),
                "pains_bool": {
                    "agreement_pct": round(pains_agree, 2),
                    "tp": pains_tp, "tn": pains_tn, "fp": pains_fp, "fn": pains_fn,
                },
                "bcut2d":     family_dict(bcut2d),
                "mqn":        family_dict(mqn),
                "slogp_vsa":  family_dict(slogp_vsa),
                "smr_vsa":    family_dict(smr_vsa),
                "peoe_vsa":   family_dict(peoe_vsa),
                "estate_vsa": family_dict(estate_vsa),
            },
            "deltas": {
                "tpsa": tpsa_deltas,
                "logp": logp_deltas,
            },
        }
        with open(args.json, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nJSON results written to {args.json}")

if __name__ == "__main__":
    main()
