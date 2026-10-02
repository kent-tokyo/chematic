#!/usr/bin/env python3
"""Emit all exposed SMARTS matches for a pinned RDKit version.

Unlike the chemistry-difference packet, this stores every query/row cell so
another binding can be compared without treating absence from a residual list
as evidence of equality.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase

from run_rdkit_python_chemistry_lane import load_queries, load_smiles, normalized_match_sets


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rdkit", required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != args.expected_rdkit:
        parser.error(f"RDKit {rdBase.rdkitVersion} != {args.expected_rdkit}")
    smiles = load_smiles(args.corpus, None)
    queries = load_queries(args.queries)
    prepared = [Chem.MolFromSmarts(value) for value in queries]
    if any(query is None for query in prepared):
        parser.error("a query is not supported by the pinned RDKit")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps({
            "_manifest": True,
            "schema": "rdkit-smarts-all-cells/v1",
            "rdkit_version": rdBase.rdkitVersion,
            "corpus_sha256": hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
            "queries_sha256": hashlib.sha256(args.queries.read_bytes()).hexdigest(),
            "input_count": len(smiles),
            "query_count": len(queries),
        }, sort_keys=True) + "\n")
        for index, value in enumerate(smiles):
            mol = Chem.MolFromSmiles(value)
            matches = (
                [normalized_match_sets(mol.GetSubstructMatches(query, uniquify=True)) for query in prepared]
                if mol is not None else None
            )
            handle.write(json.dumps({
                "input_index": index, "smiles": value, "matches": matches,
            }, separators=(",", ":")) + "\n")
    print(f"wrote {len(smiles)} rows x {len(queries)} queries to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
