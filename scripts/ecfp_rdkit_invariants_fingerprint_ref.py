#!/usr/bin/env python3
"""
Reference (NON-GATING) fingerprint-level measurements for
`ecfp4_rdkit_invariants` vs RDKit's `GetMorganFingerprintAsBitVect(radius=2)`.

These numbers are recorded for context, not as an acceptance gate: chematic
still hashes with FNV-1a (RDKit uses its own hash), so raw bit-vector
identity was never a goal of this milestone -- only the atom invariant
PARTITION (scripts/ecfp_rdkit_invariant_parity.py) is the gate. Same
methodology as the existing tier3/tier0 measurements in
scripts/ecfp4_agreement.py (Pearson correlation of pairwise Tanimoto
similarity, same default sample size), so the "before" (chematic's default,
aromaticity-inclusive invariant, r=0.94) and "after" (RdkitMorgan atom
invariant) numbers are directly comparable.

Consumes the TSV from `cargo run -p chematic-fp --release --example
ecfp4_rdkit_invariants_bits -- <SMILES.csv> <out.tsv>` (`smiles\\tbit,bit,...`).

Usage:
    .venv/bin/python scripts/ecfp_rdkit_invariants_fingerprint_ref.py <bits.tsv> [SMILES.csv] [--pairs-sample N]
"""

import argparse
import random
import statistics
import time

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs

RDLogger.DisableLog("rdApp.*")


def load_bits(path):
    out = {}
    with open(path) as f:
        for line in f:
            smi, bits = line.rstrip("\n").split("\t", 1)
            out[smi] = frozenset(int(b) for b in bits.split(",")) if bits else frozenset()
    return out


def tanimoto(a, b):
    if not a and not b:
        return 1.0
    inter = len(a & b)
    union = len(a) + len(b) - inter
    return inter / union if union else 0.0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bits_path")
    p.add_argument("csv_path", nargs="?", default="~/Downloads/SMILES.csv")
    p.add_argument("--pairs-sample", type=int, default=300)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    chem_bits = load_bits(args.bits_path)

    import os

    with open(os.path.expanduser(args.csv_path)) as f:
        smis = [line.strip() for line in f if line.strip() and line.strip() in chem_bits]

    rng = random.Random(args.seed)
    sample = smis if len(smis) <= args.pairs_sample else rng.sample(smis, args.pairs_sample)

    chem_fps = []
    rd_fps = []
    densities = []
    start = time.perf_counter()
    for smi in sample:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        chem_fps.append(chem_bits[smi])
        densities.append(len(chem_bits[smi]))
        rd_fps.append(AllChem.GetMorganFingerprintAsBitVect(rd, 2, 2048))
    elapsed = time.perf_counter() - start

    chem_sims, rd_sims = [], []
    n = len(chem_fps)
    for i in range(n):
        for j in range(i + 1, n):
            chem_sims.append(tanimoto(chem_fps[i], chem_fps[j]))
            rd_sims.append(DataStructs.TanimotoSimilarity(rd_fps[i], rd_fps[j]))

    corr = statistics.correlation(chem_sims, rd_sims) if len(chem_sims) > 1 else None
    mean_abs_diff = (
        sum(abs(a - b) for a, b in zip(chem_sims, rd_sims)) / len(chem_sims) if chem_sims else None
    )

    # Collision rate: distinct set-bit-count vs distinct-molecule count, as a
    # coarse proxy (folding 2048 raw hash values into a 2048-bit vector will
    # always show some collision; this just records the observed density).
    all_bit_counts = [len(v) for v in chem_bits.values()]

    print(f"n_molecules (sample): {n}")
    print(f"n_pairs: {len(chem_sims)}")
    print(f"pearson_correlation vs RDKit GetMorganFingerprintAsBitVect(radius=2): {corr:.4f}" if corr is not None else "pearson_correlation: n/a")
    print(f"  reference: chematic's default (Chematic mode) tier3 correlation is r=0.94 (scripts/ecfp4_agreement.py)")
    print(f"mean_abs_tanimoto_diff: {mean_abs_diff:.4f}" if mean_abs_diff is not None else "mean_abs_tanimoto_diff: n/a")
    print()
    print(f"bit density (full corpus, {len(all_bit_counts)} molecules):")
    print(f"  mean set bits / 2048: {statistics.mean(all_bit_counts):.1f}")
    print(f"  median: {statistics.median(all_bit_counts):.0f}")
    print(f"  min/max: {min(all_bit_counts)}/{max(all_bit_counts)}")
    print()
    print(f"runtime (RDKit fingerprint generation, {n} molecules): {elapsed:.3f}s ({1000*elapsed/n:.3f}ms/mol)")
    print("  (chematic-side runtime not separately timed here -- Rust snapshot generation for")
    print("   the full 5,000-molecule corpus completed in well under 1s; see the example's own output)")


if __name__ == "__main__":
    main()
