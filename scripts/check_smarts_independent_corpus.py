#!/usr/bin/env python3
"""SMARTS match sets on a corpus independent of the exposed 10k lane.

Rows of ``--corpus`` that also occur in the exposed 10k corpus are dropped, so
no row used to develop the opt-in profile is counted. Three steps, so the
chematic side can run in an interpreter RDKit does not support (the PyPI
v1.0.34 Linux wheel is CPython 3.9 only):

  emit     chematic only: per cell, native ``smarts_find`` and opt-in
           ``find_matches_rdkit_parity`` results (JSON lines)
  oracle   RDKit only: ``GetSubstructMatches`` (uniquify, no chirality)
  classify compare and write row-level accounting

Each cell is ``exact``, ``typed_refusal``/``typed_unsupported`` (opt-in only)
or ``differs``; a parse error on either side is recorded as such. Match sets
are compared as sets of sorted atom tuples, as in the 310k gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXPOSED = ROOT / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
QUERIES = ROOT / "validation/rdkit_rebaseline_smarts_queries.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(corpus: Path, limit: int) -> list[str]:
    exposed = {line.split()[0] for line in EXPOSED.read_text(encoding="utf-8").splitlines() if line.strip()}
    out = []
    for line in corpus.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        smiles = line.split()[0]
        if smiles not in exposed:
            out.append(smiles)
    return out[:limit] if limit else out


def sets(matches) -> list[list[int]]:
    return sorted({tuple(sorted(m)) for m in matches})


def emit(args) -> None:
    import chematic

    queries = json.loads(QUERIES.read_text(encoding="utf-8"))["queries"]
    corpus = rows(args.corpus, args.limit)
    with args.output.open("w", encoding="utf-8") as out:
        out.write(json.dumps({"chematic": chematic.__version__, "python": platform.python_version(),
                              "corpus_sha256": sha256(args.corpus), "rows": len(corpus),
                              "queries_sha256": sha256(QUERIES)}) + "\n")
        for index, smiles in enumerate(corpus):
            try:
                mol = chematic.from_smiles(smiles)
            except ValueError as exc:
                out.write(json.dumps({"i": index, "smiles": smiles, "parse_error": str(exc)}) + "\n")
                continue
            native, optin = [], []
            for query in queries:
                try:
                    native.append(sets(chematic.smarts_find(query, mol)))
                except ValueError as exc:
                    native.append({"error": str(exc)})
                result = mol.find_matches_rdkit_parity(query)
                optin.append(sets(result["matches"]) if result["status"] == "ok"
                             else {"status": result["status"], "reason": result["reason"]})
            out.write(json.dumps({"i": index, "smiles": smiles, "native": native, "optin": optin}) + "\n")
    print(f"{args.output}: {len(corpus)} rows")


def oracle(args) -> None:
    from rdkit import Chem, RDLogger, rdBase

    RDLogger.DisableLog("rdApp.*")
    queries = json.loads(QUERIES.read_text(encoding="utf-8"))["queries"]
    patterns = [Chem.MolFromSmarts(q) for q in queries]
    corpus = rows(args.corpus, args.limit)
    with args.output.open("w", encoding="utf-8") as out:
        out.write(json.dumps({"rdkit": rdBase.rdkitVersion, "corpus_sha256": sha256(args.corpus),
                              "rows": len(corpus), "queries_sha256": sha256(QUERIES)}) + "\n")
        for index, smiles in enumerate(corpus):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                out.write(json.dumps({"i": index, "smiles": smiles, "parse_error": True}) + "\n")
                continue
            out.write(json.dumps({"i": index, "smiles": smiles, "matches": [
                sets(mol.GetSubstructMatches(p, uniquify=True, maxMatches=100000)) for p in patterns
            ]}) + "\n")
    print(f"{args.output}: {len(corpus)} rows")


def classify(args) -> int:
    queries = json.loads(QUERIES.read_text(encoding="utf-8"))["queries"]
    with args.chematic.open(encoding="utf-8") as c, args.rdkit.open(encoding="utf-8") as r:
        c_head, r_head = json.loads(next(c)), json.loads(next(r))
        if (c_head["corpus_sha256"], c_head["rows"], c_head["queries_sha256"]) != (
                r_head["corpus_sha256"], r_head["rows"], r_head["queries_sha256"]):
            raise ValueError("chematic and RDKit runs used different inputs")
        counts = {"native": Counter(), "optin": Counter()}
        cells = {"native": [], "optin": []}
        for c_line, r_line in zip(c, r):
            crow, rrow = json.loads(c_line), json.loads(r_line)
            if crow["smiles"] != rrow["smiles"]:
                raise ValueError("rows are not aligned")
            if "parse_error" in crow or "parse_error" in rrow:
                for profile in counts:
                    side = ("both" if "parse_error" in crow and "parse_error" in rrow else
                            "chematic" if "parse_error" in crow else "rdkit")
                    counts[profile][f"parse_error_{side}"] += len(queries)
                continue
            for q, query in enumerate(queries):
                want = [list(m) for m in rrow["matches"][q]]
                for profile in ("native", "optin"):
                    got = crow[profile][q]
                    if isinstance(got, dict):
                        outcome = got.get("status", "query_error")
                    else:
                        outcome = "exact" if got == want else "differs"
                    counts[profile][outcome] += 1
                    if outcome != "exact":
                        cells[profile].append({"i": crow["i"], "smiles": crow["smiles"], "query": query,
                                               "outcome": outcome, "chematic": got, "rdkit": want})
    report = {"schema": "smarts-independent-corpus/v1", "chematic_run": c_head, "rdkit_run": r_head,
              "corpus": str(args.corpus_label), "queries": len(queries),
              "counts": {p: dict(sorted(v.items())) for p, v in counts.items()},
              "non_exact_cells": cells}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["counts"]))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="step", required=True)
    for name in ("emit", "oracle"):
        p = sub.add_parser(name)
        p.add_argument("--corpus", type=Path, required=True)
        p.add_argument("--limit", type=int, default=0)
        p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("classify")
    p.add_argument("--chematic", type=Path, required=True)
    p.add_argument("--rdkit", type=Path, required=True)
    p.add_argument("--corpus-label", default="")
    p.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.step == "emit":
        emit(args)
    elif args.step == "oracle":
        oracle(args)
    else:
        return classify(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
