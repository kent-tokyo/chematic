#!/usr/bin/env python3
"""Audit per-engine MMFF/UFF construction support for selected COD records."""
import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from platinum_development_probe import graph_matches


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sdf", type=Path)
    p.add_argument("split", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--partition", choices=("development", "holdout"), default="development")
    args = p.parse_args()
    RDLogger.DisableLog("rdApp.*")
    split = json.loads(args.split.read_text(encoding="utf-8"))
    source_hash = hashlib.sha256(args.sdf.read_bytes()).hexdigest()
    if source_hash != split["source_sha256"]:
        raise ValueError("source hash mismatch")
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    rows = []
    for record in split["rows"]:
        if record.get("split") != args.partition:
            continue
        result = {"index": record["index"], "name": record["name"]}
        try:
            source = supplier[record["index"]]
            if source is None or source.GetProp("_Name") != record["name"]:
                raise ValueError("reference parse or identity failure")
            graph = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.RemoveHs(source), isomericSmiles=True))
            if graph is None or len(Chem.GetMolFrags(graph)) != 1:
                raise ValueError("invalid or multifragment graph")
            mol = Chem.AddHs(graph)
            result["rdkit_mmff"] = bool(AllChem.MMFFHasAllMoleculeParams(mol))
            result["rdkit_uff"] = bool(AllChem.UFFHasAllMoleculeParams(mol))
            try:
                other = chematic.from_smiles(Chem.MolToSmiles(graph, isomericSmiles=True)).add_hydrogens()
                result["graph_matches"] = graph_matches(other, mol)
                # Construction tests parameterization, not energy at meaningful coordinates.
                coords = [[1.5 * i, 0.0, 0.0] for i in range(mol.GetNumAtoms())]
                for variant in ("mmff", "uff"):
                    try:
                        (other.mmff_force_field(coords) if variant == "mmff" else other.uff_force_field(coords))
                        result["chematic_" + variant] = True
                    except Exception as exc:
                        result["chematic_" + variant] = False
                        result["chematic_" + variant + "_error"] = repr(exc)
            except Exception as exc:
                result["chematic_graph_error"] = repr(exc)
        except Exception as exc:
            result["input_error"] = repr(exc)
        rows.append(result)
    output = {"source_sha256": source_hash, "split_sha256": hashlib.sha256(args.split.read_bytes()).hexdigest(), "rdkit_version": rdBase.rdkitVersion, "partition": args.partition, "counts": dict(Counter((r.get("rdkit_mmff"), r.get("chematic_mmff")) for r in rows)), "results": rows}
    # JSON object keys must be strings; retain the paired support table as rows.
    output["counts"] = [{"rdkit_mmff": k[0], "chematic_mmff": k[1], "records": v} for k, v in Counter((r.get("rdkit_mmff"), r.get("chematic_mmff")) for r in rows).items()]
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output["counts"]))


if __name__ == "__main__":
    main()
