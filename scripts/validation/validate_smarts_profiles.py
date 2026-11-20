#!/usr/bin/env python3
"""Compare native and versioned SMARTS atom sets by input row and query."""

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase

from rdkit_cli_search_differential import QUERIES


def atom_sets(matches):
    return sorted({tuple(sorted(mapping)) for mapping in matches})


def audit(path, profile, output):
    rd_queries = [Chem.MolFromSmarts(query) for query in QUERIES]
    native_queries = [chematic.compile_smarts(query) for query in QUERIES]
    counts = {"parse": Counter(), **{query: {"native": Counter(), "parity": Counter()}
                                      for query in QUERIES}}
    with gzip.open(output, "wt", encoding="utf-8") as rows:
        for row, line in enumerate(path.read_text().splitlines(), 1):
            if not line.strip() or line.startswith("#"):
                continue
            smiles = line.split()[0]
            ref = Chem.MolFromSmiles(smiles)
            try:
                mol = chematic.from_smiles(smiles)
                error = None
            except ValueError as exc:
                mol, error = None, str(exc)
            if ref is None or mol is None:
                status = f"chematic_{mol is not None}/rdkit_{ref is not None}"
                counts["parse"][status] += 1
                rows.write(json.dumps({"row": row, "smiles": smiles, "parse": status,
                                       "error": error}) + "\n")
                continue
            counts["parse"]["jointly_parsed"] += 1
            faults = []
            for query, rd_query, native_query in zip(QUERIES, rd_queries, native_queries, strict=True):
                rd_matches = ref.GetSubstructMatches(rd_query, maxMatches=100_000)
                if len(rd_matches) == 100_000:
                    raise RuntimeError(f"oracle mapping limit at row {row} query {query}")
                expected = atom_sets(rd_matches)
                native = atom_sets(native_query.find_matches(mol))
                state = "match" if native == expected else "mismatch"
                counts[query]["native"][state] += 1
                if state != "match":
                    faults.append({"query": query, "lane": "native", "expected": expected,
                                   "actual": native})
                result = mol.find_matches_rdkit_parity(query, profile=profile)
                if result["status"] == "ok":
                    actual = atom_sets(result["matches"])
                    state = "match" if actual == expected else "mismatch"
                else:
                    actual = None
                    state = result["status"]
                counts[query]["parity"][state] += 1
                if state != "match":
                    faults.append({"query": query, "lane": "parity", "expected": expected,
                                   "actual": actual, "status": state, "reason": result["reason"]})
            rows.write(json.dumps({"row": row, "smiles": smiles, "differences": faults}) + "\n")
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("corpus", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", choices=("2026.03.6", "2026.09.1"), default="2026.09.1")
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")
    output = args.output.with_suffix(".jsonl.gz")
    counts = audit(args.corpus, args.profile, output)
    totals = {key: sum((counts[q][key] for q in QUERIES), Counter())
              for key in ("native", "parity")}
    args.output.write_text(json.dumps({
        "corpus": str(args.corpus), "corpus_sha256": hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
        "rdkit_version": rdBase.rdkitVersion,
        "chematic_import": chematic.__file__, "parity_profile": args.profile,
        "output_contract": "unordered unique atom index sets, chirality off",
        "queries": QUERIES, "counts": counts, "totals": totals,
        "row_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }, indent=2) + "\n")
    print(totals)


if __name__ == "__main__":
    main()
