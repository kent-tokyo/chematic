#!/usr/bin/env python3
"""Conservatively audit invalid exclusion graphs against eligible scaffolds.

A matching scaffold must be a subgraph of its input molecule. Ignore bond
orders and aromaticity when matching so invalid valences do not conceal a
possible overlap; matches are exclusions, not proof of chemical equivalence.
"""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold


def relaxed_graph_query(mol):
    query = Chem.RWMol()
    for atom in mol.GetAtoms():
        query.AddAtom(Chem.AtomFromSmarts(f"[#{atom.GetAtomicNum()}]"))
    for bond in mol.GetBonds():
        query.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), Chem.BondType.SINGLE)
        query.ReplaceBond(query.GetNumBonds() - 1, Chem.BondFromSmarts("~"))
    return Chem.MolFromSmarts(Chem.MolToSmarts(query.GetMol()))


def audit(split, sdf):
    targets = []
    for entry in split["invalid_exclusions"]:
        molecule = Chem.MolFromSmiles(entry["smiles"], sanitize=False)
        if molecule is None:
            raise ValueError("Cannot parse an invalid exclusion even without sanitization")
        molecule.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(molecule)
        targets.append(molecule)
    supplier = Chem.SDMolSupplier(str(sdf), removeHs=False)
    if len(supplier) != len(split["rows"]):
        raise ValueError("SDF and split record counts differ")
    cache = {}
    overlaps = []
    for row, mol in zip(split["rows"], supplier):
        if row.get("split") not in ("development", "holdout"):
            continue
        if mol is None:
            raise ValueError(f"Cannot parse eligible SDF record {row['index']}")
        key = row["scaffold_sha256"]
        if key not in cache:
            heavy = Chem.RemoveHs(mol)
            graph = Chem.MolFromSmiles(Chem.MolToSmiles(heavy, isomericSmiles=True))
            if graph is None:
                raise ValueError(f"Cannot reparse eligible graph {row['index']}")
            scaffold = MurckoScaffold.GetScaffoldForMol(graph)
            if scaffold.GetNumAtoms() == 0:
                matches = list(range(len(targets)))
            else:
                query = relaxed_graph_query(scaffold)
                matches = [i for i, target in enumerate(targets) if target.HasSubstructMatch(query)]
            cache[key] = matches
        if cache[key]:
            overlaps.append({"index": row["index"], "split": row["split"], "invalid_indices": cache[key]})
    return {"eligible_scaffolds": len(cache), "overlap_counts": dict(Counter(x["split"] for x in overlaps)), "overlaps": overlaps}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdf", type=Path)
    parser.add_argument("split", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")
    split_bytes = args.split.read_bytes()
    split = json.loads(split_bytes)
    source_hash = hashlib.sha256(args.sdf.read_bytes()).hexdigest()
    if source_hash != split["source_sha256"]:
        raise ValueError("SDF does not match split source")
    result = audit(split, args.sdf)
    result.update({"source_sha256": source_hash, "split_sha256": hashlib.sha256(split_bytes).hexdigest(), "rdkit_version": rdBase.rdkitVersion, "invalid_exclusion_count": len(split["invalid_exclusions"]), "rule": "Match eligible scaffold as an atom-number graph substructure of each unsanitized exclusion, ignoring bond orders; any match requires exclusion before scoring."})
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["overlap_counts"], sort_keys=True))


if __name__ == "__main__":
    main()
