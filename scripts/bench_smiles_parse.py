#!/usr/bin/env python3
"""
Benchmark SMILES parsing throughput: chematic vs RDKit.

Uses a built-in 20-molecule corpus (simple → drug-like → stereo → large)
or a user-supplied CSV file.  Repeats the corpus to reach --n total parses.

Usage:
    python scripts/bench_smiles_parse.py
    python scripts/bench_smiles_parse.py --rdkit
    python scripts/bench_smiles_parse.py --n 10000 --rdkit
    python scripts/bench_smiles_parse.py ~/Downloads/SMILES.csv --n 5000 --rdkit --json
"""

import argparse
import csv
import json
import sys
import time


# Diverse corpus: simple alkanes → drug-like → stereocentres → polycyclic
BUILTIN_SMILES: list[str] = [
    # Simple
    "C", "CC", "CCC", "c1ccccc1", "c1ccncc1",
    # Drug-like
    "CC(=O)Oc1ccccc1C(=O)O",            # aspirin
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",      # caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",      # ibuprofen
    "CC(=O)Nc1ccc(O)cc1",               # paracetamol
    "c1ccc2cc3ccccc3cc2c1",             # pyrene
    # Stereocentres
    "N[C@@H](Cc1ccccc1)C(=O)O",        # L-phenylalanine
    "OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",  # glucose
    # Rings / fused
    "C1CCC2CCCCC2C1",                   # decalin
    "C1CC2CCCCC2CC1",                   # bicyclo
    "C1=CC2=CC3=CC=CC=C3C=C2C=C1",     # coronene-like
    # SMARTS-challenging
    "Clc1ccccc1Cl",
    "F[B-](F)(F)F.[Na+]",              # salt
    "CN(C)C(=N)NC(=N)N",               # metformin (multiple N)
    "[NH4+].[O-]C(=O)c1ccccc1",       # charged pair
    "OC(=O)[C@H](N)Cc1c[nH]cn1",      # histidine
]


def build_corpus(smiles_list: list[str], n: int) -> list[str]:
    reps = (n // len(smiles_list)) + 1
    return (smiles_list * reps)[:n]


def run_chematic(corpus: list[str]) -> tuple[float, int]:
    import chematic
    # Warm-up (fills intern caches)
    for s in corpus[:min(20, len(corpus))]:
        chematic.from_smiles(s)
    t0 = time.perf_counter()
    ok = sum(1 for s in corpus if chematic.from_smiles(s) is not None)
    return time.perf_counter() - t0, ok


def run_rdkit(corpus: list[str]) -> tuple[float, int]:
    from rdkit import Chem
    for s in corpus[:min(20, len(corpus))]:
        Chem.MolFromSmiles(s)
    t0 = time.perf_counter()
    ok = sum(1 for s in corpus if Chem.MolFromSmiles(s) is not None)
    return time.perf_counter() - t0, ok


def fmt_row(name: str, elapsed: float, n: int, ok: int) -> str:
    us   = elapsed / n * 1e6
    rate = n / elapsed
    return (f"  {name:<12}  {n:>6} SMILES  "
            f"{elapsed*1000:>7.1f} ms  "
            f"{us:>6.2f} µs/mol  "
            f"{rate:>9,.0f} mol/s"
            + (f"  ({ok}/{n} parsed)" if ok < n else ""))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("smiles_file", nargs="?",
                    help="CSV with SMILES column (or plain text, one SMILES per line)")
    ap.add_argument("--rdkit",  action="store_true", help="Also benchmark RDKit")
    ap.add_argument("--n",      type=int, default=5000, help="Total parses to time")
    ap.add_argument("--json",   action="store_true",   help="Machine-readable output")
    args = ap.parse_args()

    # --- load SMILES ---
    smiles_list = BUILTIN_SMILES
    if args.smiles_file:
        loaded: list[str] = []
        with open(args.smiles_file) as f:
            first = f.read(1024); f.seek(0)
            if "," in first or "\t" in first:
                reader = csv.DictReader(f)
                fields = reader.fieldnames or []
                col = "SMILES" if "SMILES" in fields else (fields[0] if fields else None)
                loaded = [r[col].strip() for r in reader if col and r.get(col,"").strip()]
            else:
                loaded = [l.strip() for l in f if l.strip() and not l.startswith("#")]
        smiles_list = loaded
        if not args.json:
            print(f"Loaded {len(smiles_list)} SMILES from {args.smiles_file}\n")

    corpus = build_corpus(smiles_list, args.n)
    results: dict[str, object] = {"n": args.n, "corpus_size": len(smiles_list)}

    if not args.json:
        print(f"SMILES parse benchmark  (n={args.n})\n")

    # --- chematic ---
    try:
        elapsed, ok = run_chematic(corpus)
        results["chematic"] = {
            "total_ms":   round(elapsed * 1000, 1),
            "us_per_mol": round(elapsed / args.n * 1e6, 2),
            "mol_per_sec": int(args.n / elapsed),
            "parsed_ok":  ok,
        }
        if not args.json:
            print(fmt_row("chematic", elapsed, args.n, ok))
    except ImportError:
        if not args.json:
            print("  chematic not installed")

    # --- rdkit ---
    if args.rdkit:
        try:
            elapsed, ok = run_rdkit(corpus)
            results["rdkit"] = {
                "total_ms":   round(elapsed * 1000, 1),
                "us_per_mol": round(elapsed / args.n * 1e6, 2),
                "mol_per_sec": int(args.n / elapsed),
                "parsed_ok":  ok,
            }
            if not args.json:
                print(fmt_row("rdkit", elapsed, args.n, ok))
        except ImportError:
            if not args.json:
                print("  rdkit not installed — skipping")

    # --- speedup summary ---
    if "chematic" in results and "rdkit" in results:
        ch = results["chematic"]["total_ms"]
        rd = results["rdkit"]["total_ms"]
        speedup = rd / ch
        results["speedup_x"] = round(speedup, 1)
        if not args.json:
            print(f"\n  chematic is {speedup:.1f}× faster to parse than RDKit")

    if args.json:
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
