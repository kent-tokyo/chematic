#!/usr/bin/env python3
"""
RDKit benchmark for comparison with chematic Rust benchmarks.

Usage:
    python3 scripts/rdkit_benchmark.py

Measures median time (µs) per molecule for the same 10-molecule set
used in chematic's criterion benchmarks.

Requirements: rdkit (pip install rdkit or conda install -c conda-forge rdkit)
"""

import sys
import timeit
import statistics

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
except ImportError:
    print("ERROR: RDKit not installed. Run: conda install -c conda-forge rdkit")
    sys.exit(1)

# Same 10-molecule set as chematic benchmarks
BENCH_SMILES = [
    "c1ccccc1",                              # benzene
    "Cc1ccccc1",                             # toluene
    "CC(=O)Oc1ccccc1C(=O)O",                # aspirin
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",           # caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",           # ibuprofen
    "c1ccncc1",                              # pyridine
    "C1CCNCC1",                              # piperidine
    "CC(=O)Nc1ccc(O)cc1",                   # paracetamol
    "NCC(=O)O",                              # glycine
    "CN(C)C(=N)NC(=N)N",                    # metformin
]

N = 1000  # number of repetitions for timeit


def measure(label: str, fn, n: int = N) -> float:
    """Run fn n times and return median µs per call."""
    times = timeit.repeat(fn, number=1, repeat=n)
    median_s = statistics.median(times)
    median_us = median_s * 1e6
    print(f"  {label:<35} {median_us:8.2f} µs/call  ({n} reps)")
    return median_us


def bench_parse() -> float:
    """Parse 10 SMILES → 10 Mol objects."""
    def _fn():
        for smi in BENCH_SMILES:
            _ = Chem.MolFromSmiles(smi)
    return measure("parse_smiles_10mol", _fn)


def bench_ecfp4(mols) -> float:
    """Generate ECFP4 (Morgan r=2, 2048 bits) for 10 molecules."""
    def _fn():
        for mol in mols:
            _ = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    return measure("ecfp4_10mol", _fn)


def bench_tanimoto(mols) -> float:
    """Pairwise Tanimoto (55 pairs from 10 mols)."""
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]
    from rdkit.DataStructs import TanimotoSimilarity

    def _fn():
        for i in range(len(fps)):
            for j in range(i, len(fps)):
                _ = TanimotoSimilarity(fps[i], fps[j])
    return measure("tanimoto_ecfp4_pairs (55 pairs)", _fn)


def bench_descriptors(mols) -> float:
    """Compute 5 descriptors (MW, LogP, TPSA, HBD, HBA) for 10 molecules."""
    def _fn():
        for mol in mols:
            _ = Descriptors.MolWt(mol)
            _ = Descriptors.MolLogP(mol)
            _ = Descriptors.TPSA(mol)
            _ = rdMolDescriptors.CalcNumHBD(mol)
            _ = rdMolDescriptors.CalcNumHBA(mol)
    return measure("descriptors_5x10mol", _fn)


def bench_qed(mols) -> float:
    """QED drug-likeness score for 10 molecules."""
    from rdkit.Chem import QED

    def _fn():
        for mol in mols:
            _ = QED.qed(mol)
    return measure("qed_10mol", _fn)


def bench_smarts_match(mols) -> float:
    """SMARTS match (benzene ring) against 10 molecules."""
    query = Chem.MolFromSmarts("c1ccccc1")

    def _fn():
        for mol in mols:
            _ = mol.HasSubstructMatch(query)
    return measure("smarts_match_10mol", _fn)


def main():
    print(f"RDKit Benchmark — {len(BENCH_SMILES)} molecules, {N} repetitions each")
    print("=" * 70)

    # Pre-parse molecules
    mols = [Chem.MolFromSmiles(s) for s in BENCH_SMILES]
    assert all(m is not None for m in mols), "Some SMILES failed to parse!"

    print("\n--- SMILES Parsing ---")
    t_parse = bench_parse()

    print("\n--- Fingerprints ---")
    t_ecfp4 = bench_ecfp4(mols)
    t_tan = bench_tanimoto(mols)

    print("\n--- Descriptors ---")
    t_desc = bench_descriptors(mols)
    t_qed = bench_qed(mols)

    print("\n--- SMARTS ---")
    t_smarts = bench_smarts_match(mols)

    print("\n" + "=" * 70)
    print("Summary (per molecule, median µs):")
    print(f"  SMILES parse:      {t_parse / len(BENCH_SMILES):8.2f} µs/mol")
    print(f"  ECFP4:             {t_ecfp4 / len(BENCH_SMILES):8.2f} µs/mol")
    print(f"  5 descriptors:     {t_desc / len(BENCH_SMILES):8.2f} µs/mol")
    print(f"  QED:               {t_qed / len(BENCH_SMILES):8.2f} µs/mol")
    print(f"  SMARTS match:      {t_smarts / len(BENCH_SMILES):8.2f} µs/mol")
    print()
    print("Update docs/benchmark_results.md with these values.")


if __name__ == "__main__":
    main()
