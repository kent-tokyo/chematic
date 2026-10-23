#!/usr/bin/env python3
"""Record chematic's products for the BioTransformer rule corpus, without RDKit.

Runs every rule of the three pinned BioTransformer tables on every reactant
request of ``validation/biotransformer-requests-400.tsv`` (the implicit-H and
explicit-H SMILES RDKit 2026.03.6 writes for the corpus's 400 reactants) and
writes one gzip JSON line per pair: ``[rule, request, status, detail,
products]`` after a header line with the chematic version. RDKit then
compares them on Linux with ``biotransformer_rule_corpus.py --replay``, so a
published wheel for a platform RDKit has no wheel for (macOS x86-64) is still
measured, and the comparison runs against one pinned RDKit build.

Only the standard library, json5 and chematic are imported (Python 3.9 syntax).
"""

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import chematic  # noqa: E402
from biotransformer_rules import load_requests, load_rules  # noqa: E402
from chematic_reaction_worker import run  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--rules", type=Path, nargs="+", required=True)
    ap.add_argument("--requests", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True, help="JSONL, gzip-compressed")
    ap.add_argument("--rules-slice", default=None,
                    help="START:END, record only rules[START:END] (indices stay those of the full list)")
    args = ap.parse_args()

    rules, sources = load_rules(args.rules)
    start, end = 0, len(rules)
    if args.rules_slice:
        a, _, b = args.rules_slice.partition(":")
        start, end = int(a or 0), int(b) if b else len(rules)
    smiles = [s for pair in load_requests(args.requests) for s in pair]
    started = time.time()
    with gzip.open(args.output, "wt", encoding="utf-8") as out:
        out.write(json.dumps({"version": chematic.__version__, "file": chematic.__file__,
                              "rules_sources": sources, "rules": len(rules),
                              "requests": len(smiles)}) + "\n")
        for i in range(start, min(end, len(rules))):
            rule = rules[i]
            for k, smi in enumerate(smiles):
                r = run({"smirks": rule["smirks"], "smiles": smi})
                out.write(json.dumps([i, k, r["status"], r["detail"], r["products"]],
                                     separators=(",", ":")) + "\n")
    print(json.dumps({"rules": len(rules), "requests": len(smiles),
                      "elapsed_seconds": round(time.time() - started, 1)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
