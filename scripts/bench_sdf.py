#!/usr/bin/env python3
"""
Benchmark SDF read/write throughput: chematic vs RDKit.

Uses RDKit's own bundled sample SDF (Contrib/PBF/testData/egfr.sdf, ~365
records) as a reproducible corpus with no extra download step. Falls back
to a user-supplied SDF path.

Read:  chematic.iter_sdf(path)       vs RDKit Chem.SDMolSupplier(path)
Write: chematic.SDWriter(path)       vs RDKit Chem.SDWriter(path)

Usage:
    python scripts/bench_sdf.py
    python scripts/bench_sdf.py --rdkit
    python scripts/bench_sdf.py path/to/other.sdf --rdkit --json
"""

import argparse
import importlib.util
import json
import os
import tempfile
import time


def default_sdf_path() -> str | None:
    """Locate RDKit's bundled egfr.sdf sample via the installed rdkit package."""
    spec = importlib.util.find_spec("rdkit")
    if spec is None or not spec.submodule_search_locations:
        return None
    rdkit_dir = list(spec.submodule_search_locations)[0]
    candidate = os.path.join(rdkit_dir, "Contrib", "PBF", "testData", "egfr.sdf")
    return candidate if os.path.exists(candidate) else None


def run_chematic_read(path: str) -> tuple[float, int]:
    import chematic
    # Warm-up (import/JIT-free, but keeps timing consistent with other scripts)
    ok = 0
    t0 = time.perf_counter()
    for mol in chematic.SDMolSupplier(path, strictParsing=False):
        # Match RDKit's SDMolSupplier contract: parse the molecule and fields,
        # without requesting chematic's optional stereo/3D diagnostics or
        # canonical SMILES generation.
        if mol is not None:
            ok += 1
    return time.perf_counter() - t0, ok


def run_rdkit_read(path: str) -> tuple[float, int]:
    from rdkit import Chem
    ok = 0
    t0 = time.perf_counter()
    for mol in Chem.SDMolSupplier(path):
        if mol is not None:
            ok += 1
    return time.perf_counter() - t0, ok


def run_chematic_write(path: str, out_path: str) -> tuple[float, int]:
    import chematic
    mols = [rec.mol for rec in chematic.iter_sdf(path) if rec.mol is not None]
    t0 = time.perf_counter()
    # Serialization-only mode is compared separately from depiction. RDKit's
    # writer does not run a new 2D layout for each molecule either.
    with chematic.SDWriter(out_path, compute2d=False) as w:
        for m in mols:
            w.write(m)
    return time.perf_counter() - t0, len(mols)


def run_rdkit_write(path: str, out_path: str) -> tuple[float, int]:
    from rdkit import Chem
    mols = [m for m in Chem.SDMolSupplier(path) if m is not None]
    for mol in mols:
        mol.RemoveAllConformers()
    t0 = time.perf_counter()
    w = Chem.SDWriter(out_path)
    for m in mols:
        w.write(m)
    w.close()
    return time.perf_counter() - t0, len(mols)


def fmt_row(name: str, elapsed: float, n: int) -> str:
    if n == 0:
        return f"    {name:<12}  no records"
    us   = elapsed / n * 1e6
    rate = n / elapsed
    return (f"    {name:<12}  {n:>6} mols  "
            f"{elapsed*1000:>7.1f} ms  "
            f"{us:>6.2f} µs/mol  "
            f"{rate:>9,.0f} mol/s")


def report_pair(label: str, results: dict, chematic_r, rdkit_r, args) -> None:
    section: dict[str, object] = {}
    if not args.json:
        print(f"  [{label}]")

    if chematic_r is not None:
        elapsed, n = chematic_r
        section["chematic"] = {
            "total_ms":    round(elapsed * 1000, 1),
            "us_per_mol":  round(elapsed / n * 1e6, 2) if n else None,
            "mol_per_sec": int(n / elapsed) if elapsed else None,
            "n":           n,
        }
        if not args.json:
            print(fmt_row("chematic", elapsed, n))

    if rdkit_r is not None:
        elapsed, n = rdkit_r
        section["rdkit"] = {
            "total_ms":    round(elapsed * 1000, 1),
            "us_per_mol":  round(elapsed / n * 1e6, 2) if n else None,
            "mol_per_sec": int(n / elapsed) if elapsed else None,
            "n":           n,
        }
        if not args.json:
            print(fmt_row("rdkit", elapsed, n))

    if "chematic" in section and "rdkit" in section:
        ch = section["chematic"]["total_ms"]
        rd = section["rdkit"]["total_ms"]
        speedup = rd / ch if ch else float("inf")
        section["speedup_x"] = round(speedup, 1)
        if not args.json:
            direction = "faster" if speedup >= 1 else "slower"
            factor = speedup if speedup >= 1 else 1 / speedup
            print(f"    → chematic is {factor:.1f}× {direction} than RDKit")

    results[label] = section
    if not args.json:
        print()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sdf_file", nargs="?", help="SDF corpus (default: RDKit's bundled egfr.sdf)")
    ap.add_argument("--rdkit", action="store_true", help="Also benchmark RDKit")
    ap.add_argument("--json",  action="store_true", help="Machine-readable output")
    args = ap.parse_args()

    path = args.sdf_file or default_sdf_path()
    if path is None or not os.path.exists(path):
        raise SystemExit(
            "No SDF corpus found. Pass one explicitly, e.g.\n"
            "  python scripts/bench_sdf.py path/to/some.sdf"
        )

    results: dict[str, object] = {"source": path}
    if not args.json:
        print(f"SDF read/write benchmark  (corpus: {path})\n")

    # --- read ---
    ch_read = None
    try:
        ch_read = run_chematic_read(path)
    except ImportError:
        if not args.json:
            print("  chematic not installed")
    rd_read = run_rdkit_read(path) if args.rdkit else None
    report_pair("read", results, ch_read, rd_read, args)

    # --- write ---
    with tempfile.TemporaryDirectory() as tmp:
        ch_write = None
        try:
            ch_write = run_chematic_write(path, os.path.join(tmp, "chematic_out.sdf"))
        except ImportError:
            if not args.json:
                print("  chematic not installed")
        rd_write = run_rdkit_write(path, os.path.join(tmp, "rdkit_out.sdf")) if args.rdkit else None
        report_pair("write", results, ch_write, rd_write, args)

    if args.json:
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
