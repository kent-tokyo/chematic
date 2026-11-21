#!/usr/bin/env python3
"""Assign a pinned PDB macrocycle index to non-overlapping scaffold partitions."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold


INDEX_SHA256 = "cf408bbcbeed78fe8862f8f1ccf48a2714f18dcaad5b9ac1fa6c482f4d7155dc"
SOURCE_COMMIT = "6afbb010770f8adaadb9db7f83f7a86538757a62"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("index", type=Path)
    parser.add_argument("--exclude-split", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error("the existing scaffold hashes were produced with RDKit 2026.09.1")
    index_bytes = args.index.read_bytes()
    if digest(index_bytes) != INDEX_SHA256:
        parser.error("index content does not match the pinned source commit")
    RDLogger.DisableLog("rdApp.*")
    known = set()
    invalid = {}
    exclusion_hashes = {}
    for path in args.exclude_split:
        raw = path.read_bytes()
        exclusion_hashes[str(path)] = digest(raw)
        data = json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)
        known.update(r["scaffold_sha256"] for r in data["rows"] if "scaffold_sha256" in r)
        for record in data.get("invalid_exclusions", []):
            invalid[(record["path"], record["smiles"])] = record
    with args.index.open(encoding="utf-8", newline="") as source:
        index = list(csv.DictReader(source))
    if len(index) != 1861 or len({r["case_id"] for r in index}) != len(index):
        raise ValueError("unexpected count or duplicate case identifiers")
    rows, canonical = [], defaultdict(set)
    for i, record in enumerate(index):
        row = {"index": i, "name": record["case_id"], "pdb_id": record["pdb_id"],
               "primary_ring_size": int(record["primary_ring_size"]),
               "input_smiles_sha256": digest(record["unprotonated_isomeric_smiles"].encode())}
        try:
            mol = Chem.MolFromSmiles(record["unprotonated_isomeric_smiles"])
            if mol is None or len(Chem.GetMolFrags(mol)) != 1:
                raise ValueError("SMILES parse failed or is not a single component")
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            graph = Chem.MolFromSmiles(smiles)
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=graph, includeChirality=False)
            key = digest(scaffold.encode())
            partition = "excluded_known_scaffold" if key in known else (
                "development" if int(key[:8], 16) % 5 == 0 else "holdout"
            )
            row.update(scaffold_sha256=key, split=partition)
            canonical[partition].add(smiles)
        except Exception as exc:
            row.update(status="graph_or_scaffold_failure", error=repr(exc))
        rows.append(row)
    output = {
        "source_url": "https://github.com/dptech-corp/MacroCycle-Benchmark",
        "source_commit": SOURCE_COMMIT,
        "index_sha256": INDEX_SHA256,
        "rdkit_version": rdBase.rdkitVersion,
        "excluded_splits_sha256": exclusion_hashes,
        "known_scaffold_count": len(known),
        "invalid_exclusions": list(invalid.values()),
        "rule": "Exclude every scaffold in supplied splits, then assign SHA-256 of non-chiral canonical Murcko scaffold by first 32 bits modulo 5; zero is development and others holdout. No reference coordinates or conformer outcomes are used.",
        "counts": dict(Counter(r.get("split", r.get("status")) for r in rows)),
        "scaffold_counts": {p: len({r["scaffold_sha256"] for r in rows if r.get("split") == p}) for p in canonical},
        "unique_isomeric_smiles_counts": {p: len(values) for p, values in canonical.items()},
        "rows": rows,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"counts": output["counts"], "scaffold_counts": output["scaffold_counts"], "invalid_exclusions": len(invalid)}, sort_keys=True))


if __name__ == "__main__":
    main()
