#!/usr/bin/env python3
"""Dump first/repeated descriptor results; preserve failures and exact float bits."""
import argparse
import hashlib
import json
from pathlib import Path
import struct


def evaluate(function, molecule):
    try:
        value = function(molecule)
        return ({"float_bits": struct.pack(">d", value).hex()} if isinstance(value, float)
                else {"integer": value})
    except Exception as error:
        return {"error": type(error).__name__ + ": " + str(error)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", choices=("chematic", "rdkit"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.engine == "chematic":
        import chematic
        parse = chematic.from_smiles
        version = chematic.__version__
        functions = {name: (lambda m, key=name: getattr(m, key))
                     for name in ("tpsa", "rdkit_tpsa", "labute_asa", "num_rings")}
    else:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdMolDescriptors as D
        RDLogger.DisableLog("rdApp.*")
        parse = Chem.MolFromSmiles
        version = Chem.rdBase.rdkitVersion
        functions = {"tpsa": lambda m: D.CalcTPSA(m, includeSandP=True),
                     "rdkit_tpsa": D.CalcTPSA, "labute_asa": D.CalcLabuteASA, "num_rings": D.CalcNumRings}
    root = Path(__file__).resolve().parents[1]
    corpora = [root / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi",
               root / "validation/cosmolkit_comparison/corpora/stress-v1.smi"]
    with args.output.open("w") as output:
        for path in corpora:
            data = path.read_bytes()
            rows = [line.split()[0] for line in data.decode().splitlines() if line.strip()]
            output.write(json.dumps({"corpus": path.name, "sha256": hashlib.sha256(data).hexdigest(),
                                     "rows": len(rows), "engine": args.engine, "version": version}) + "\n")
            for index, smiles in enumerate(rows):
                row = {"corpus": path.name, "index": index}
                try:
                    molecule = parse(smiles)
                    if molecule is None:
                        raise ValueError("parse failed")
                    row["first"] = {key: evaluate(fn, molecule) for key, fn in functions.items()}
                    row["repeated"] = {key: evaluate(fn, molecule) for key, fn in functions.items()}
                except Exception as error:
                    row["parse_error"] = type(error).__name__ + ": " + str(error)
                output.write(json.dumps(row, sort_keys=True) + "\n")
            print(path.name, len(rows), "done", flush=True)


if __name__ == "__main__":
    main()
