#!/usr/bin/env python3
"""Compare one attempt with a bounded square-planar retry on exposed inputs."""

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from pathlib import Path

import chematic

DEVELOPMENT = [
    *(f"[Pt@SP{i}](F)(Cl)(Br)I" for i in (1, 2, 3)),
    *(f"F[Pt@SP{i}](Cl)Br" for i in (1, 2, 3)),
    "O=C(O[I@SP3]1OB(O)c2c(F)cccc21)C(F)(F)F",
    "Brc1ccc2c3c(cccc13)[Se@@](c1ccccc1)=[Se@SP3]2Br",
]
SYNTHETIC = [pattern.format(i) for pattern in (
    "[Pd@SP{}](F)(Cl)(Br)N", "N[Pt@SP{}](F)Cl",
    "[Ni@SP{}](F)(Cl)(Br)S", "O[Pd@SP{}](Cl)Br",
) for i in (1, 2, 3)]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_one(mol, seed, retry):
    try:
        result = (mol.rdkit_embed_square_planar_with_retry(random_seed=seed) if retry
                  else {"coordinates": mol.rdkit_embed_square_planar(random_seed=seed),
                        "attempts": 1, "used_tightened_bounds": False})
        angles = mol.square_planar_geometry_diagnostics(result["coordinates"])
        return {"status": "returned", **result, "angles": angles}
    except (RuntimeError, ValueError) as error:
        return {"status": "error", "error_type": type(error).__name__, "error": str(error)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    records = []
    for group, inputs, seeds in (("development", DEVELOPMENT, (1, 42)),
                                 ("synthetic", SYNTHETIC, (7, 23))):
        for smiles in inputs:
            mol = chematic.from_smiles(smiles).add_hydrogens()
            for seed in seeds:
                records.append({"group": group, "smiles": smiles, "seed": seed,
                                "single": run_one(mol, seed, False),
                                "retry": run_one(mol, seed, True)})
    extension = importlib.util.find_spec("chematic.chematic").origin
    report = {"selection": "exposed development and authored synthetic molecules; neither is a holdout",
              "platform": platform.platform(), "python": sys.version,
              "wheel_sha256": digest(args.wheel), "extension_sha256": digest(extension),
              "script_sha256": digest(__file__), "records": records}
    Path(args.output).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
