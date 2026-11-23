#!/usr/bin/env python3
"""Conservatively audit invalid prior exclusions against macrocycle index scaffolds."""
import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("index", type=Path)
    parser.add_argument("split", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    split_bytes, index_bytes = args.split.read_bytes(), args.index.read_bytes()
    split = json.loads(split_bytes)
    if sha(index_bytes) != split["index_sha256"] or rdBase.rdkitVersion != split["rdkit_version"]:
        raise ValueError("index hash or RDKit version differs from the split")
    RDLogger.DisableLog("rdApp.*")
    helper_path = Path(__file__).resolve().parents[1] / "benchmarks" / "platinum_invalid_exclusion_audit.py"
    spec = importlib.util.spec_from_file_location("platinum_invalid_exclusion_audit", helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    targets = []
    for entry in split["invalid_exclusions"]:
        target = Chem.MolFromSmiles(entry["smiles"], sanitize=False)
        if target is None:
            raise ValueError("an invalid exclusion cannot be parsed without sanitization")
        target.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(target)
        targets.append(target)
    with args.index.open(encoding="utf-8", newline="") as source:
        index = list(csv.DictReader(source))
    if len(index) != len(split["rows"]):
        raise ValueError("index and split counts differ")
    cache, overlaps = {}, []
    for row in split["rows"]:
        if row.get("split") not in ("development", "holdout"):
            continue
        record = index[row["index"]]
        if record["case_id"] != row["name"] or sha(record["unprotonated_isomeric_smiles"].encode()) != row["input_smiles_sha256"]:
            raise ValueError("index identity differs from split")
        key = row["scaffold_sha256"]
        if key not in cache:
            mol = Chem.MolFromSmiles(record["unprotonated_isomeric_smiles"])
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))
            text = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
            if sha(text.encode()) != key:
                raise ValueError("index scaffold differs from split")
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            if scaffold.GetNumAtoms() == 0:
                matches = list(range(len(targets)))
            else:
                query = helper.relaxed_graph_query(scaffold)
                matches = [i for i, target in enumerate(targets) if target.HasSubstructMatch(query)]
            cache[key] = matches
        if cache[key]:
            overlaps.append({"index": row["index"], "split": row["split"], "invalid_indices": cache[key]})
    output = {
        "index_sha256": sha(index_bytes), "split_sha256": sha(split_bytes),
        "rdkit_version": rdBase.rdkitVersion, "invalid_exclusion_count": len(targets),
        "eligible_scaffolds": len(cache), "overlap_counts": dict(Counter(x["split"] for x in overlaps)),
        "overlaps": overlaps,
        "rule": "Match each eligible non-chiral scaffold as an atom-number graph substructure of each unsanitized prior invalid exclusion, ignoring bond order. Match requires excluding its complete scaffold group before scoring.",
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"eligible_scaffolds": len(cache), "overlap_counts": output["overlap_counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
