#!/usr/bin/env python3
"""
Benchmark SMARTS substructure matching throughput: chematic vs RDKit.

Pairs with scripts/rdkit_benchmark.py's bench_smarts_match (benzene-ring
query) and mirrors the query tiers already covered by chematic-smarts'
Criterion benches (crates/chematic-smarts/benches/smarts_bench.rs):
a simple ring query and a recursive SMARTS query. Repeats a 10-molecule
corpus to reach --n total matches per query.

Usage:
    python scripts/bench_smarts.py
    python scripts/bench_smarts.py --rdkit
    python scripts/bench_smarts.py --n 5000 --rdkit --json
"""

import argparse
import json
import time


# Same 10-molecule set as crates/chematic-smarts/benches/smarts_bench.rs
BENCH_SMILES: list[str] = [
    "c1ccccc1",
    "Cc1ccccc1",
    "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "c1ccncc1",
    "C1CCNCC1",
    "CC(=O)Nc1ccc(O)cc1",
    "NCC(=O)O",
    "CN(C)C(=N)NC(=N)N",
]

# name -> SMARTS pattern. "ring" is a plain aromatic-ring query (cheap);
# "recursive" is a recursive SMARTS (amide N-H), the expensive tier.
QUERIES: dict[str, str] = {
    "ring": "c1ccccc1",
    "recursive": "[NH;$(NC=O)]",
}


def build_corpus(n: int) -> list[str]:
    reps = (n // len(BENCH_SMILES)) + 1
    return (BENCH_SMILES * reps)[:n]


def run_chematic(corpus: list[str], pattern: str) -> tuple[float, int]:
    import chematic
    mols = [chematic.from_smiles(s) for s in corpus]
    # Warm-up — fills the SmartsCache LRU for this pattern
    for m in mols[:min(20, len(mols))]:
        chematic.smarts_match(pattern, m)
    t0 = time.perf_counter()
    hits = sum(1 for m in mols if chematic.smarts_match(pattern, m))
    return time.perf_counter() - t0, hits


def run_rdkit(corpus: list[str], pattern: str) -> tuple[float, int]:
    from rdkit import Chem
    mols = [Chem.MolFromSmiles(s) for s in corpus]
    query = Chem.MolFromSmarts(pattern)  # compiled once, like a cached chematic query
    for m in mols[:min(20, len(mols))]:
        m.HasSubstructMatch(query)
    t0 = time.perf_counter()
    hits = sum(1 for m in mols if m.HasSubstructMatch(query))
    return time.perf_counter() - t0, hits


def fmt_row(name: str, elapsed: float, n: int, hits: int) -> str:
    us   = elapsed / n * 1e6
    rate = n / elapsed
    return (f"    {name:<12}  {n:>6} mols  "
            f"{elapsed*1000:>7.1f} ms  "
            f"{us:>6.2f} µs/mol  "
            f"{rate:>9,.0f} mol/s  "
            f"({hits}/{n} matched)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rdkit", action="store_true", help="Also benchmark RDKit")
    ap.add_argument("--n",     type=int, default=5000, help="Total matches to time, per query")
    ap.add_argument("--json",  action="store_true",    help="Machine-readable output")
    args = ap.parse_args()

    corpus = build_corpus(args.n)
    results: dict[str, object] = {"n": args.n, "corpus_size": len(BENCH_SMILES), "queries": {}}

    if not args.json:
        print(f"SMARTS match benchmark  (n={args.n} per query)\n")

    for label, pattern in QUERIES.items():
        query_results: dict[str, object] = {"pattern": pattern}
        if not args.json:
            print(f"  [{label}] {pattern}")

        try:
            elapsed, hits = run_chematic(corpus, pattern)
            query_results["chematic"] = {
                "total_ms":    round(elapsed * 1000, 1),
                "us_per_mol":  round(elapsed / args.n * 1e6, 2),
                "mol_per_sec": int(args.n / elapsed),
                "matched":     hits,
            }
            if not args.json:
                print(fmt_row("chematic", elapsed, args.n, hits))
        except ImportError:
            if not args.json:
                print("    chematic not installed")

        if args.rdkit:
            try:
                elapsed, hits = run_rdkit(corpus, pattern)
                query_results["rdkit"] = {
                    "total_ms":    round(elapsed * 1000, 1),
                    "us_per_mol":  round(elapsed / args.n * 1e6, 2),
                    "mol_per_sec": int(args.n / elapsed),
                    "matched":     hits,
                }
                if not args.json:
                    print(fmt_row("rdkit", elapsed, args.n, hits))
            except ImportError:
                if not args.json:
                    print("    rdkit not installed — skipping")

        if "chematic" in query_results and "rdkit" in query_results:
            ch = query_results["chematic"]["total_ms"]
            rd = query_results["rdkit"]["total_ms"]
            speedup = rd / ch
            query_results["speedup_x"] = round(speedup, 1)
            if not args.json:
                direction = "faster" if speedup >= 1 else "slower"
                factor = speedup if speedup >= 1 else 1 / speedup
                print(f"    → chematic is {factor:.1f}× {direction} than RDKit")

        results["queries"][label] = query_results
        if not args.json:
            print()

    if args.json:
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
