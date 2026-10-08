#!/usr/bin/env python3
"""Time the same operations in chematic, COSMolKit and RDKit (Python bindings).

Each engine runs in its own interpreter. A block runs every engine once, in
an order that rotates between blocks; per operation the median over blocks is
reported. Molecules are re-parsed (untimed) before each operation so that a
result cached on a molecule object cannot be reused. Calls return each
library's native result object (no normalisation), and SMARTS queries are
compiled once where the library offers a compiled query.

Timings depend on the host; the report records it. Speed is only meaningful
next to the correctness summary from ``compare_corpus.py``: an operation whose
outputs differ is a different amount of work.

    python bench_corpus.py --corpus CORPUS.smi --output bench.json --blocks 5
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
OPS = ["parse", "canonical_smiles", "formula", "mol_wt", "tpsa", "logp", "qed",
       "morgan2_2048", "maccs", "cip", "smarts_31", "molblock", "inchi"]


def table(engine: str):
    from corpus_engines import SMARTS_QUERIES
    if engine == "rdkit":
        from rdkit import Chem, RDLogger
        from rdkit.Chem import QED, Crippen, Descriptors, MACCSkeys, rdCIPLabeler, rdFingerprintGenerator
        from rdkit.Chem import rdMolDescriptors as D
        RDLogger.DisableLog("rdApp.*")
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        qs = [Chem.MolFromSmarts(q) for q in SMARTS_QUERIES]
        return Chem.MolFromSmiles, {
            "canonical_smiles": Chem.MolToSmiles, "formula": D.CalcMolFormula,
            "mol_wt": Descriptors.MolWt, "tpsa": D.CalcTPSA, "logp": Crippen.MolLogP,
            "qed": QED.qed, "morgan2_2048": gen.GetFingerprint, "maccs": MACCSkeys.GenMACCSKeys,
            "cip": rdCIPLabeler.AssignCIPLabels,
            "smarts_31": lambda m: [m.GetSubstructMatches(q, maxMatches=100000) for q in qs],
            "molblock": Chem.MolToMolBlock, "inchi": Chem.MolToInchi,
        }
    if engine == "chematic":
        import chematic as c
        return c.from_smiles, {
            "canonical_smiles": lambda m: m.smiles, "formula": lambda m: m.formula,
            "mol_wt": lambda m: m.rdkit_mw, "tpsa": lambda m: m.rdkit_tpsa,
            "logp": lambda m: m.logp, "qed": lambda m: m.qed,
            "morgan2_2048": lambda m: m.rdkit_ecfp_config(2, 2048), "maccs": lambda m: m.maccs(),
            "cip": lambda m: m.cip_stereo(mode="accurate"),
            "smarts_31": lambda m: [m.find_matches_rdkit_parity(q) for q in SMARTS_QUERIES],
            "molblock": lambda m: m.to_mol_block(),
        }
    import cosmolkit as ck
    qs = [ck.parse_smarts(q) for q in SMARTS_QUERIES]
    return ck.Molecule.from_smiles, {
        "canonical_smiles": lambda m: m.to_smiles(), "formula": ck.calc_mol_formula,
        "mol_wt": ck.calc_mol_wt, "tpsa": ck.calc_tpsa,
        "logp": lambda m: ck.calc_crippen_descriptors(m)[0], "qed": ck.calc_qed,
        "morgan2_2048": lambda m: m.fingerprint_morgan(radius=2, n_bits=2048),
        "maccs": lambda m: m.maccs_fingerprint(), "cip": lambda m: m.with_cip_labels(),
        "smarts_31": lambda m: [ck.get_substruct_matches(m, q, max_matches=100000) for q in qs],
        "molblock": lambda m: m.to_2d_sdf_string(), "inchi": lambda m: m.to_inchi(),
    }


def measure(engine: str, corpus: Path, ops: list[str]) -> dict:
    warnings.simplefilter("ignore")
    parse, fns = table(engine)
    rows = [l.split()[0] for l in corpus.read_text().splitlines() if l.strip()]

    def parsed():
        out = []
        for s in rows:
            try:
                m = parse(s)
            except Exception:  # noqa: BLE001
                m = None
            if m is not None:
                out.append(m)
        return out

    res = {}
    t = time.perf_counter()
    mols = parsed()
    res["parse"] = time.perf_counter() - t
    for op in ops:
        if op == "parse":
            continue
        fn = fns.get(op)
        if fn is None:
            res[op] = None
            continue
        mols = parsed()
        errors = 0
        t = time.perf_counter()
        for m in mols:
            try:
                fn(m)
            except Exception:  # noqa: BLE001 - counted, still timed
                errors += 1
        res[op] = time.perf_counter() - t
        if errors:
            res[op + "__errors"] = errors
    res["parsed"] = len(mols)
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--corpus", type=Path, required=True)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--ops", default=",".join(OPS))
    ap.add_argument("--python", action="append", metavar="ENGINE=PYTHON", default=[],
                    help="interpreter per engine (default: this one)")
    ap.add_argument("--_child")
    args = ap.parse_args()
    ops = args.ops.split(",")
    if args._child:
        print(json.dumps(measure(args._child, args.corpus, ops)))
        return 0
    pythons = {"rdkit": sys.executable, "chematic": sys.executable, "cosmolkit": sys.executable}
    pythons.update(dict(p.split("=", 1) for p in args.python))
    engines = list(pythons)
    samples: dict[str, list[dict]] = {e: [] for e in engines}
    versions = {}
    for block in range(args.blocks):
        order = engines[block % len(engines):] + engines[: block % len(engines)]
        for e in order:
            out = subprocess.run([pythons[e], __file__, "--corpus", str(args.corpus), "--ops", args.ops,
                                  "--_child", e], capture_output=True, text=True, check=True)
            samples[e].append(json.loads(out.stdout))
        print(f"block {block + 1}/{args.blocks} done", file=sys.stderr)
    for e in engines:
        code = {"rdkit": "from rdkit import rdBase;print(rdBase.rdkitVersion)",
                "chematic": "import chematic;print(chematic.__version__)",
                "cosmolkit": "import cosmolkit;print(cosmolkit.__version__)"}[e]
        versions[e] = subprocess.run([pythons[e], "-c", code], capture_output=True, text=True).stdout.strip()
    report = {"schema": "cosmolkit-corpus-bench/v1", "corpus": args.corpus.name,
              "blocks": args.blocks, "host": {"platform": platform.platform(),
                                               "cpus": os.cpu_count(), "python": platform.python_version()},
              "versions": versions, "median_seconds": {}, "parsed": {}, "errors": {}}
    for e in engines:
        report["parsed"][e] = samples[e][0]["parsed"]
        report["median_seconds"][e] = {
            op: (None if samples[e][0].get(op) is None
                 else statistics.median(s[op] for s in samples[e])) for op in ops}
        report["errors"][e] = {k: v for k, v in samples[e][0].items() if k.endswith("__errors")}
    ck, ch, rd = (report["median_seconds"][e] for e in ("cosmolkit", "chematic", "rdkit"))
    report["chematic_speedup_vs_cosmolkit"] = {
        op: (round(ck[op] / ch[op], 2) if ck.get(op) and ch.get(op) else None) for op in ops}
    report["chematic_speedup_vs_rdkit"] = {
        op: (round(rd[op] / ch[op], 2) if rd.get(op) and ch.get(op) else None) for op in ops}
    report["cosmolkit_speedup_vs_rdkit"] = {
        op: (round(rd[op] / ck[op], 2) if rd.get(op) and ck.get(op) else None) for op in ops}
    text = json.dumps(report, indent=1, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(f"{'op':18}{'chematic ms':>13}{'COSMolKit ms':>14}{'RDKit ms':>11}{'ch/ck':>8}")
    for op in ops:
        f = lambda v: f"{v * 1000:.0f}" if v else "-"  # noqa: E731
        print(f"{op:18}{f(ch.get(op)):>13}{f(ck.get(op)):>14}{f(rd.get(op)):>11}"
              f"{str(report['chematic_speedup_vs_cosmolkit'][op] or '-') + 'x':>8}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
