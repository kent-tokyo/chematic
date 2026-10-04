#!/usr/bin/env python3
"""Find SMARTS ring-query cells whose RDKit answer depends on input atom order.

For every corpus row with at least three rings, the molecule is renumbered
with ``--seeds`` random permutations, written and re-parsed, and the matched
atoms of each ring query are mapped back to the original numbering. A cell
whose matched-atom set differs between permutations has no order-independent
RDKit answer. Used to adjudicate the opt-in SMARTS profile's
``ring_model_ambiguous`` refusals.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

RDLogger.DisableLog("rdApp.*")
QUERIES = ["[R1]", "[R2]", "[R3]", "[r5]", "[r6]", "[k5]", "[k6]", "[R]", "[R0]", "[x2]", "[x3]"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--seeds", type=int, default=6)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    lines = [line.split()[0] for line in args.corpus.read_text(encoding="utf-8").splitlines() if line.strip()]
    patterns = [Chem.MolFromSmarts(q) for q in QUERIES]
    checked, dependent = 0, {}
    for row, smiles in enumerate(lines):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetRingInfo().NumRings() < 3:
            continue
        checked += 1
        variants = [set() for _ in patterns]
        for seed in range(args.seeds):
            order = list(range(mol.GetNumAtoms()))
            random.Random(seed).shuffle(order)
            renumbered = Chem.RenumberAtoms(mol, order)
            reparsed = Chem.MolFromSmiles(Chem.MolToSmiles(renumbered, canonical=False))
            written = json.loads(renumbered.GetProp("_smilesAtomOutputOrder").replace(",]", "]"))
            for q, pattern in enumerate(patterns):
                variants[q].add(frozenset(order[written[m[0]]] for m in reparsed.GetSubstructMatches(pattern)))
        cells = [QUERIES[q] for q, v in enumerate(variants) if len(v) > 1]
        if cells:
            dependent[row] = cells
    report = {"schema": "rdkit-ring-order-dependence/v1", "rdkit": rdBase.rdkitVersion,
              "corpus_sha256": hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
              "seeds": args.seeds, "queries": QUERIES, "rows": len(lines),
              "rows_with_three_or_more_rings": checked, "order_dependent_cells": dependent}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({"checked": checked, "order_dependent_rows": len(dependent)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
