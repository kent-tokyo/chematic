#!/usr/bin/env python3
"""Audit element, formal charge and aromaticity identity for selected records."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sdf", type=Path)
    p.add_argument("split", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--partition", choices=("development", "holdout"), default="holdout")
    args = p.parse_args()
    RDLogger.DisableLog("rdApp.*")
    split_bytes = args.split.read_bytes()
    split = json.loads(split_bytes)
    if hashlib.sha256(args.sdf.read_bytes()).hexdigest() != split["source_sha256"]:
        raise ValueError("SDF hash mismatch")
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    results = []
    for row in split["rows"]:
        if row.get("split") != args.partition:
            continue
        result = {"index": row["index"], "name": row["name"]}
        try:
            source = supplier[row["index"]]
            if source is None or source.GetProp("_Name") != row["name"]:
                raise ValueError("reference parse or identity failure")
            smiles = Chem.MolToSmiles(Chem.RemoveHs(source), isomericSmiles=True)
            reference = Chem.MolFromSmiles(smiles)
            if reference is None:
                raise ValueError("reference graph parse failed")
            reference = Chem.AddHs(reference)
            other = chematic.from_smiles(smiles).add_hydrogens()
            observed = [(a[1], a[2], bool(a[3])) for a in other.atom_table]
            expected = [(a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsAromatic()) for a in reference.GetAtoms()]
            result["atom_identity_match"] = observed == expected
            result["isotopic_atoms"] = sum(a.GetIsotope() != 0 for a in reference.GetAtoms())
            if result["isotopic_atoms"]:
                result["isotopic_formula_match"] = other.formula_with_isotopes == rdMolDescriptors.CalcMolFormula(reference, separateIsotopes=True, abbreviateHIsotopes=False)
            if observed != expected:
                result["mismatches"] = [i for i, (a, b) in enumerate(zip(observed, expected)) if a != b]
                result["atom_count_match"] = len(observed) == len(expected)
        except Exception as exc:
            result["error"] = repr(exc)
        results.append(result)
    summary = dict(Counter("error" if "error" in r else "match" if r["atom_identity_match"] else "mismatch" for r in results))
    out = {"source_sha256": split["source_sha256"], "split_sha256": hashlib.sha256(split_bytes).hexdigest(), "rdkit_version": rdBase.rdkitVersion, "partition": args.partition, "summary": summary, "records_with_isotopic_atoms": sum(r.get("isotopic_atoms", 0) > 0 for r in results), "results": results}
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": summary, "records_with_isotopic_atoms": out["records_with_isotopic_atoms"]}, sort_keys=True))


if __name__ == "__main__":
    main()
