#!/usr/bin/env python3
"""Check every randomized SMILES against a pinned RDKit stereo identity."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from rdkit import Chem, RDLogger


def read_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", required=True, type=Path)
    parser.add_argument("--dump", required=True, type=Path)
    parser.add_argument("--seeds", type=int, default=16)
    args = parser.parse_args()
    if Chem.rdBase.rdkitVersion != "2025.09.3":
        raise SystemExit("requires RDKit 2025.09.3")
    if args.seeds < 1:
        raise SystemExit("--seeds must be positive")
    RDLogger.DisableLog("rdApp.*")

    source = [line.split()[0] for line in args.corpus.read_text().splitlines()]
    molecules = [Chem.MolFromSmiles(smiles) for smiles in source]
    if any(mol is None for mol in molecules):
        raise SystemExit("corpus contains a molecule RDKit cannot parse")
    expected = [Chem.MolToSmiles(mol) for mol in molecules]

    seen: set[tuple[int, int]] = set()
    mismatches = []
    with read_text(args.dump) as rows:
        for line in rows:
            index_text, seed_text, smiles = line.rstrip("\n").split("\t", 2)
            index, seed = int(index_text), int(seed_text)
            key = (index, seed)
            if not (0 <= index < len(source) and 1 <= seed <= args.seeds) or key in seen:
                raise SystemExit(f"unexpected or duplicate row: {key}")
            seen.add(key)
            mol = Chem.MolFromSmiles(smiles)
            actual = Chem.MolToSmiles(mol) if mol is not None else None
            if actual != expected[index]:
                mismatches.append({
                    "index": index,
                    "seed": seed,
                    "expected": expected[index],
                    "actual": actual,
                    "output": smiles,
                })
    required = len(source) * args.seeds
    if len(seen) != required:
        raise SystemExit(f"missing rows: saw {len(seen)} of {required}")
    print(json.dumps({
        "schema": "random-smiles-stereo-comparison/v1",
        "corpus_sha256": hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
        "dump_sha256": hashlib.sha256(args.dump.read_bytes()).hexdigest(),
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "rows": len(source),
        "seeds_per_row": args.seeds,
        "tested": required,
        "exact": required - len(mismatches),
        "mismatches": mismatches[:20],
    }, sort_keys=True))
    if mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
