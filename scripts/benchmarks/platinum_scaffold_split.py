#!/usr/bin/env python3
"""Assign Platinum SDF records to stable scaffold-disjoint development and holdout sets."""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sdf", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--exclude-manifest", type=Path, action="append", default=[])
    parser.add_argument("--exclude-smiles", type=Path, action="append", default=[])
    parser.add_argument("--exclude-jsonl", type=Path, action="append", default=[])
    parser.add_argument("--allow-invalid-exclusions", action="store_true")
    args = parser.parse_args()
    source_hash = hashlib.sha256(args.sdf.read_bytes()).hexdigest()
    excluded = set()
    invalid_exclusions = []
    for path in args.exclude_manifest:
        for row in json.loads(path.read_text())["molecules"]:
            mol = Chem.MolFromSmiles(row["smiles"])
            if mol is None:
                raise ValueError(f"invalid SMILES in {path}: {row['smiles']}")
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(
                mol=mol, includeChirality=False
            )
            excluded.add(hashlib.sha256(scaffold.encode()).hexdigest())
    for path in args.exclude_smiles:
        for line in path.read_text().splitlines():
            smiles = line.strip()
            if not smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                if not args.allow_invalid_exclusions:
                    raise ValueError(f"invalid SMILES in {path}: {smiles}")
                invalid_exclusions.append({"path": str(path), "smiles": smiles})
                continue
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(
                mol=mol, includeChirality=False
            )
            excluded.add(hashlib.sha256(scaffold.encode()).hexdigest())
    for path in args.exclude_jsonl:
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            smiles = row["smiles"]
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                if not args.allow_invalid_exclusions:
                    raise ValueError(f"invalid SMILES in {path}: {smiles}")
                invalid_exclusions.append({"path": str(path), "smiles": smiles})
                continue
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
            excluded.add(hashlib.sha256(scaffold.encode()).hexdigest())
    rows = []
    counts = Counter()
    for index, molecule in enumerate(Chem.SDMolSupplier(str(args.sdf), removeHs=False)):
        if molecule is None:
            rows.append({"index": index, "status": "parse_failure"})
            counts["parse_failure"] += 1
            continue
        heavy = Chem.RemoveHs(molecule)
        # Reparse the graph to discard SDF atom properties before decomposition.
        graph = Chem.MolFromSmiles(Chem.MolToSmiles(heavy, isomericSmiles=True))
        if graph is None:
            rows.append({"index": index, "status": "graph_parse_failure"})
            counts["graph_parse_failure"] += 1
            continue
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(
            mol=graph, includeChirality=False
        )
        scaffold_hash = hashlib.sha256(scaffold.encode()).hexdigest()
        if scaffold_hash in excluded:
            split = "excluded_known_scaffold"
        else:
            split = "development" if int(scaffold_hash[:8], 16) % 5 == 0 else "holdout"
        name = molecule.GetProp("_Name") if molecule.HasProp("_Name") else ""
        rows.append(
            {
                "index": index,
                "name": name,
                "scaffold_sha256": scaffold_hash,
                "split": split,
            }
        )
        counts[split] += 1
    exclusions = {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in args.exclude_manifest + args.exclude_smiles + args.exclude_jsonl
    }
    args.output.write_text(
        json.dumps(
            {
                "source_sha256": source_hash,
                "rdkit_version": rdBase.rdkitVersion,
                "excluded_manifests_sha256": exclusions,
                "invalid_exclusions": invalid_exclusions,
                "rule": "Exclude listed scaffolds, then SHA-256 of non-chiral canonical Murcko scaffold, first 32 bits modulo 5; zero is development, others holdout; empty scaffolds remain together",
                "counts": dict(counts),
                "rows": rows,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(dict(counts), sort_keys=True))


if __name__ == "__main__":
    main()
