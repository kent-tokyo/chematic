#!/usr/bin/env python3
"""Checks of COSMolKit 0.5 surfaces that do not run per corpus molecule.

* ``DetectChemistryProblems`` on ``corpus_engines.BROKEN_SMILES`` read with
  sanitize=False: problem types and atom indices against RDKit.
* Reaction SMARTS round trip (``ReactionToSmarts(ReactionFromSmarts(s))``)
  on ``corpus_engines.REACTIONS``.

    python check_new_surfaces.py [--engine chematic --engine cosmolkit] [--output out.json]
"""

from __future__ import annotations

import argparse
import json

import corpus_engines as ce


def run(engine: str, fn, items):
    out = []
    for item in items:
        try:
            out.append(["ok", fn(engine, item)])
        except ce.UnsupportedError as exc:
            out.append(["unsupported", str(exc)])
        except Exception as exc:  # noqa: BLE001 - engine failure is data
            out.append(["error", f"{type(exc).__name__}: {exc}"[:200]])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", action="append", default=None)
    ap.add_argument("--output")
    args = ap.parse_args()
    engines = args.engine or ["chematic", "cosmolkit"]
    checks = {
        "chemistry_problems_unsanitized": (ce.chemistry_problems_unsanitized, ce.BROKEN_SMILES),
        "reaction_smarts_roundtrip": (ce.reaction_smarts_roundtrip, [r[1] for r in ce.REACTIONS]),
    }
    report = {}
    for name, (fn, items) in checks.items():
        ref = run("rdkit", fn, items)
        report[name] = {"items": items, "rdkit": ref}
        for eng in engines:
            got = run(eng, fn, items)
            report[name][eng] = got
            counts = {}
            for r, g in zip(ref, got):
                k = g[0] if g[0] != "ok" else ("match" if r == g else "mismatch")
                counts[k] = counts.get(k, 0) + 1
            print(f"{name:32s} {eng:10s} {counts}")
            for item, r, g in zip(items, ref, got):
                if g[0] == "ok" and r != g:
                    print(f"    {item}\n      rdkit: {r[1]}\n      {eng}: {g[1]}")
    if args.output:
        with open(args.output, "w") as f:
            json.dump(report, f, indent=1)


if __name__ == "__main__":
    main()
