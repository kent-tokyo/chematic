#!/usr/bin/env python3
"""
How far does find_sssr's non-determinism/non-minimality (see project plan,
and ringinfo_parity.py) leak into descriptors currently advertised as "100%
RDKit agreement"?

Measures worst-of-N self-stability (does the SAME molecule give the SAME
value across N independently-traversed valid Kekulized SMILES) for a set of
metrics chosen to span "should be ring-independent" through "should be
maximally exposed to ring perception":

  MolWt, TPSA, HBA, HBD          -- predicted stable (no ring dependency)
  RingCount                      -- predicted stable (mu is a topological invariant)
  NumAromaticRings                -- predicted unstable (~3%, matches the
                                     aromaticity-flag order-dependence already found)
  NumSaturatedRings, NumAliphaticRings -- predicted unstable
  LogP, MolarRefractivity        -- predicted unstable (Crippen atom-typing
                                     uses ring membership)
  Murcko scaffold (canonical SMILES) -- predicted substantially unstable
  [r5] / [r6] SMARTS match count -- predicted unstable

Self-stability only (oracle-independent: any variation for the identical
molecule is a bug regardless of what RDKit says). RDKit agreement is a
separate, already-well-covered question (bench5k.py).

Usage:
    python scripts/ring_collateral_damage.py [SMILES.csv] [--limit N] [-n N_VARIANTS]
"""
import sys
import re
import argparse


def map_nums_in_order(smi):
    return [int(x) for x in re.findall(r":(\d+)\]", smi)]


METRICS = [
    "mol_wt", "tpsa", "hba", "hbd",
    "ring_count", "num_aromatic_rings", "num_saturated_rings", "num_aliphatic_rings",
    "logp", "mr", "scaffold", "r5_matches", "r6_matches",
]


def compute(mol):
    return {
        "mol_wt": round(mol.mw, 4),
        "tpsa": round(mol.tpsa, 4),
        "hba": mol.hba,
        "hbd": mol.hbd,
        "ring_count": mol.ring_count,
        "num_aromatic_rings": mol.aromatic_ring_count,
        "num_saturated_rings": mol.num_saturated_rings,
        "num_aliphatic_rings": mol.num_aliphatic_rings,
        "logp": round(mol.logp, 4),
        "mr": round(mol.molar_refractivity, 4),
        "scaffold": mol.scaffold(),
        "r5_matches": len(mol.find_matches("[r5]")),
        "r6_matches": len(mol.find_matches("[r6]")),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("smiles_csv", nargs="?", default="~/Downloads/SMILES.csv")
    parser.add_argument("--limit", type=int, default=None, help="Cap corpus size for a quick run")
    parser.add_argument("-n", "--variants", type=int, default=10)
    parser.add_argument("--detail", action="store_true")
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
    smis = [l.split(",")[0].strip() for l in lines if l.split(",")[0].strip().lower() != "smiles"]
    if args.limit:
        smis = smis[: args.limit]

    n_mol = 0
    stable = {m: 0 for m in METRICS}
    unstable = {m: 0 for m in METRICS}
    examples = {m: [] for m in METRICS}

    for smi in smis:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        n_mol += 1

        rd_mapped = Chem.MolFromSmiles(smi)
        for i, a in enumerate(rd_mapped.GetAtoms()):
            a.SetAtomMapNum(i + 1)
        Chem.Kekulize(rd_mapped, clearAromaticFlags=True)

        values = {m: set() for m in METRICS}
        for _ in range(args.variants):
            variant = Chem.MolToSmiles(rd_mapped, doRandom=True, kekuleSmiles=True, canonical=False)
            variant_clean = re.sub(r":\d+", "", variant)
            try:
                cm = chematic.from_smiles(variant_clean)
                result = compute(cm)
                for m in METRICS:
                    values[m].add(result[m])
            except Exception:
                for m in METRICS:
                    values[m].add(None)

        for m in METRICS:
            if len(values[m]) <= 1:
                stable[m] += 1
            else:
                unstable[m] += 1
                if args.detail and len(examples[m]) < 5:
                    examples[m].append((smi, values[m]))

    print(f"corpus: {n_mol}, variants per molecule: {args.variants}")
    print(f"\n{'metric':22s} {'stable':>8s} {'unstable':>8s} {'unstable%':>10s}")
    for m in METRICS:
        total = stable[m] + unstable[m]
        pct = 100 * unstable[m] / max(total, 1)
        print(f"{m:22s} {stable[m]:8d} {unstable[m]:8d} {pct:9.2f}%")

    if args.detail:
        print()
        for m in METRICS:
            if examples[m]:
                print(f"=== {m} examples ===")
                for smi, vals in examples[m]:
                    print(f"  {smi[:70]!r}: {vals}")


if __name__ == "__main__":
    main()
