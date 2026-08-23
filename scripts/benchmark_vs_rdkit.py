#!/usr/bin/env python3
"""
Benchmark: chematic vs RDKit — batch ECFP4 fingerprint generation.

Usage:
    # chematic only (no RDKit needed):
    python scripts/benchmark_vs_rdkit.py

    # With RDKit comparison:
    pip install rdkit
    python scripts/benchmark_vs_rdkit.py --rdkit
"""

import sys
import time
import argparse
import numpy as np

# ---------------------------------------------------------------------------
# Sample molecules (drug-like, diverse)
# ---------------------------------------------------------------------------

SMILES_SAMPLES = [
    "c1ccccc1", "CCO", "CC(=O)O", "c1ccc(cc1)O", "CC(=O)Oc1ccccc1C(=O)O",
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", "CN1CCC[C@H]1c2cccnc2",
    "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", "OC(=O)c1ccccc1O",
    "c1ccc2ccccc2c1", "C1CCCCC1", "c1ccncc1", "c1ccc(Cl)cc1",
    "CC(N)C(=O)O", "OCC(O)C(O)C(O)C(O)CO", "c1ccc(cc1)C(=O)O",
    "NC(=O)c1ccccc1", "c1ccc2[nH]ccc2c1", "CCc1ccccc1",
    "COc1ccc(cc1OC)C2CC(=O)c3ccccc3O2",
]

def generate_smiles(n: int) -> list[str]:
    """Generate n SMILES by cycling through samples."""
    return [SMILES_SAMPLES[i % len(SMILES_SAMPLES)] for i in range(n)]


def benchmark_chematic(smiles: list[str], repeats: int = 3) -> float:
    import chematic
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fps = chematic.bulk.ecfp4(smiles)
        t1 = time.perf_counter()
        times.append(t1 - t0)
    n = fps.shape[0]
    med = sorted(times)[len(times) // 2]
    return med, n


def benchmark_rdkit(smiles: list[str], repeats: int = 3) -> float:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    import numpy as np
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fps = []
        for s in smiles:
            mol = Chem.MolFromSmiles(s)
            if mol is not None:
                fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, 2048)
                arr = np.zeros((2048,), dtype=np.uint8)
                from rdkit.DataStructs import ConvertToNumpyArray
                ConvertToNumpyArray(fp, arr)
                fps.append(arr)
        t1 = time.perf_counter()
        times.append(t1 - t0)
    n = len(fps)
    med = sorted(times)[len(times) // 2]
    return med, n


def fmt(seconds: float, n: int) -> str:
    µs_per_mol = seconds * 1e6 / n
    return f"{seconds:.3f}s  ({µs_per_mol:.2f} µs/mol)"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rdkit", action="store_true", help="Include RDKit benchmark")
    parser.add_argument("--n", type=int, default=10_000, help="Number of molecules")
    parser.add_argument("--repeats", type=int, default=3, help="Timing repetitions")
    parser.add_argument("--corpus", metavar="FILE",
                        help="Also time ECFP4 over the real, diverse molecules in this "
                             "SMILES-per-line file (not cycled), reported separately from "
                             "the repeated-fixture sweep above -- never blend the two.")
    args = parser.parse_args()

    if args.corpus:
        with open(args.corpus) as f:
            corpus_smiles = [line.strip() for line in f if line.strip()]
        print(f"\nBenchmark: ECFP4 fingerprint generation -- diverse corpus ({args.corpus}, "
              f"{len(corpus_smiles)} molecules, {args.repeats} runs, median)")
        ch_time, ch_n = benchmark_chematic(corpus_smiles, args.repeats)
        row = f"chematic: {fmt(ch_time, ch_n)}"
        if args.rdkit:
            try:
                rd_time, rd_n = benchmark_rdkit(corpus_smiles, args.repeats)
                speedup = rd_time / ch_time
                row += f"   RDKit: {fmt(rd_time, rd_n)}   speedup: {speedup:.1f}x"
            except ImportError:
                row += "  [RDKit not available]"
        print(row)
        return

    sizes = [100, 1_000, args.n]
    print(f"\nBenchmark: ECFP4 fingerprint generation ({args.repeats} runs, median)")
    print(f"Platform: {sys.platform}\n")
    print(f"{'N':>8}  {'chematic (Rayon parallel)':>28}", end="")
    if args.rdkit:
        print(f"  {'RDKit (Python loop)':>24}  {'speedup':>8}", end="")
    print()
    print("-" * (8 + 30 + (36 if args.rdkit else 0)))

    for n in sizes:
        smiles = generate_smiles(n)
        ch_time, ch_n = benchmark_chematic(smiles, args.repeats)
        row = f"{n:>8}  {fmt(ch_time, ch_n):>28}"

        if args.rdkit:
            try:
                rd_time, rd_n = benchmark_rdkit(smiles, args.repeats)
                speedup = rd_time / ch_time
                row += f"  {fmt(rd_time, rd_n):>24}  {speedup:>7.1f}×"
            except ImportError:
                row += "  [RDKit not available]"

        print(row)

    # Also benchmark bulk.descriptors
    print()
    print("Benchmark: bulk.descriptors() — 55+ descriptors per molecule")
    print("-" * 50)
    import chematic
    for n in [100, 1_000]:
        smiles = generate_smiles(n)
        times = []
        for _ in range(args.repeats):
            t0 = time.perf_counter()
            descs = chematic.bulk.descriptors(smiles)
            t1 = time.perf_counter()
            times.append(t1 - t0)
        med = sorted(times)[len(times) // 2]
        n_valid = len(descs)
        print(f"{n:>8}  {fmt(med, n_valid):>28}")


if __name__ == "__main__":
    main()
