#!/usr/bin/env python3
"""
Per-atom aromaticity parity: chematic vs RDKit.

rdkit_compat_diff.py already checks aromatic atom/bond *counts* (order-invariant).
This checks the aromatic flag *per atom*, in two scenarios:

1. **as-is**: the corpus SMILES fed straight to chematic.from_smiles() (no
   apply_aromaticity() call). For aromatic-form (lowercase) input this exercises
   the "trust the parser" path — atom.aromatic is set directly during parsing.
   Does not exercise the Hückel algorithm at all.
2. **kekulized (worst-of-N)**: for each molecule, N independently-traversed
   valid Kekulized SMILES (RDKit `doRandom=True`) are each fed through
   chematic.from_smiles(...).apply_aromaticity() — this exercises chematic's
   Huckel perception algorithm, and specifically its sensitivity to which
   valid representation of a molecule it's handed (ring perception and
   Pass 1/Pass 2 propagation can both be order-sensitive; measuring only
   RDKit's single canonical Kekule form hides this). A molecule counts as
   passing only if **every** one of the N variants agrees with RDKit —
   single-representation testing (N=1) is what let known order-dependent
   molecules go undetected in earlier rounds.

Both scenarios compare against RDKit's own (always re-perceived) aromatic
flags. Alignment for the kekulized scenario uses a permanent atom-map-number
identity set before generating variants (no substructure matching — a prior
version of this script's order-invariance test used substructure matching
and had a backwards-indexing bug that produced false mismatches).

Usage:
    python scripts/aromaticity_atom_parity.py [SMILES.csv] [--detail] [--limit N] [-n N_VARIANTS] [--seed S]
"""
import sys
import re
import argparse


def per_atom_arom(mol):
    return [row[3] for row in mol.atom_table]


def map_nums_in_order(smi):
    return [int(x) for x in re.findall(r":(\d+)\]", smi)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("smiles_csv", nargs="?",
                        default="~/Downloads/SMILES.csv",
                        help="One SMILES per line, or CSV with a SMILES column")
    parser.add_argument("--detail", action="store_true",
                        help="Print every mismatching molecule to stderr")
    parser.add_argument("--limit", type=int, default=None,
                        help="Only show first N mismatches per scenario in --detail mode")
    parser.add_argument("-n", "--variants", type=int, default=10,
                        help="Number of random Kekulized traversals per molecule for the "
                             "worst-of-N kekulized scenario (default 10)")
    args = parser.parse_args()

    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        sys.exit("rdkit not installed. pip install rdkit")

    import chematic

    import os
    path = os.path.expanduser(args.smiles_csv)
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    # Tolerate a "SMILES" header / CSV with extra columns, same as bench5k.py.
    smis = []
    for line in lines:
        first = line.split(",")[0].strip()
        if first.lower() == "smiles":
            continue
        smis.append(first)

    as_is = {"mol_ok": 0, "mol_bad": 0, "atom_ok": 0, "atom_bad": 0, "shown": 0}
    kek = {"mol_ok": 0, "mol_bad": 0, "variant_ok": 0, "variant_bad": 0, "shown": 0}
    n_mol = 0

    for smi in smis:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        n_mol += 1

        # --- as_is: original SMILES, no apply_aromaticity(), trust-the-parser path ---
        rd_arom = [a.GetIsAromatic() for a in rd.GetAtoms()]
        try:
            cm = chematic.from_smiles(smi)
            cm_arom = per_atom_arom(cm)
        except Exception:
            cm_arom = None
        if cm_arom is not None and len(cm_arom) == len(rd_arom):
            matches = sum(1 for a, b in zip(cm_arom, rd_arom) if a == b)
            as_is["atom_ok"] += matches
            as_is["atom_bad"] += len(rd_arom) - matches
            if matches == len(rd_arom):
                as_is["mol_ok"] += 1
            else:
                as_is["mol_bad"] += 1
                if args.detail and (args.limit is None or as_is["shown"] < args.limit):
                    as_is["shown"] += 1
                    print(f"[as_is] MISMATCH smiles={smi!r}", file=sys.stderr)
                    print(f"  rdkit    ={rd_arom}", file=sys.stderr)
                    print(f"  chematic ={cm_arom}", file=sys.stderr)
        else:
            as_is["mol_bad"] += 1

        # --- kekulized worst-of-N: N random Kekulized traversals, atom-map aligned ---
        rd_mapped = Chem.MolFromSmiles(smi)
        for i, a in enumerate(rd_mapped.GetAtoms()):
            a.SetAtomMapNum(i + 1)
        Chem.Kekulize(rd_mapped, clearAromaticFlags=True)
        ground_truth = {i + 1: a.GetIsAromatic() for i, a in enumerate(rd.GetAtoms())}

        molecule_ok = True
        first_bad_variant = None
        for _ in range(args.variants):
            variant = Chem.MolToSmiles(rd_mapped, doRandom=True, kekuleSmiles=True, canonical=False)
            order = map_nums_in_order(variant)
            variant_clean = re.sub(r":\d+", "", variant)
            try:
                cmv = chematic.from_smiles(variant_clean).apply_aromaticity()
                cmv_arom = per_atom_arom(cmv)
            except Exception:
                cmv_arom = None
            ok = (cmv_arom is not None and len(cmv_arom) == len(order)
                  and all(dict(zip(order, cmv_arom)).get(k) == v for k, v in ground_truth.items()))
            if ok:
                kek["variant_ok"] += 1
            else:
                kek["variant_bad"] += 1
                molecule_ok = False
                if first_bad_variant is None:
                    first_bad_variant = variant

        if molecule_ok:
            kek["mol_ok"] += 1
        else:
            kek["mol_bad"] += 1
            if args.detail and (args.limit is None or kek["shown"] < args.limit):
                kek["shown"] += 1
                print(f"[kekulized worst-of-{args.variants}] MISMATCH smiles={smi!r} "
                      f"first_bad_variant={first_bad_variant!r}", file=sys.stderr)

    print(f"corpus (parsed by RDKit): {n_mol}")
    print(f"\n[as_is]")
    as_is_total = as_is["mol_ok"] + as_is["mol_bad"]
    as_is_atom_total = as_is["atom_ok"] + as_is["atom_bad"]
    print(f"  molecule-level (all atoms match): {as_is['mol_ok']}/{as_is_total} "
          f"({100 * as_is['mol_ok'] / max(as_is_total, 1):.2f}%)")
    print(f"  atom-level agreement:             {as_is['atom_ok']}/{as_is_atom_total} "
          f"({100 * as_is['atom_ok'] / max(as_is_atom_total, 1):.2f}%)")

    print(f"\n[kekulized, worst-of-{args.variants}]")
    kek_total = kek["mol_ok"] + kek["mol_bad"]
    kek_variant_total = kek["variant_ok"] + kek["variant_bad"]
    print(f"  molecule-level (all {args.variants} variants match): {kek['mol_ok']}/{kek_total} "
          f"({100 * kek['mol_ok'] / max(kek_total, 1):.2f}%)")
    print(f"  variant-level agreement:                  {kek['variant_ok']}/{kek_variant_total} "
          f"({100 * kek['variant_ok'] / max(kek_variant_total, 1):.2f}%)")


if __name__ == "__main__":
    main()
