#!/usr/bin/env python3
"""Compare the existing and opt-in square-planar embedding on exposed cases."""

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from pathlib import Path

import chematic

CASES = [
    *(f"[Pt@SP{i}](F)(Cl)(Br)I" for i in (1, 2, 3)),
    *(f"F[Pt@SP{i}](Cl)Br" for i in (1, 2, 3)),
    "O=C(O[I@SP3]1OB(O)c2c(F)cccc21)C(F)(F)F",
    "Brc1ccc2c3c(cccc13)[Se@@](c1ccccc1)=[Se@SP3]2Br",
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_one(mol, seed, constrained):
    try:
        coords = (mol.rdkit_embed_square_planar(random_seed=seed) if constrained
                  else mol.rdkit_embed(random_seed=seed))
        angles = mol.square_planar_geometry_diagnostics(coords)
        in_policy = all(report["unavailable"] is None
                        and all(angle >= 150.0 for _, _, angle in report["trans"])
                        and all(60.0 <= angle <= 120.0 for _, _, angle in report["cis"])
                        for report in angles)
        return {"status": "returned", "coordinates": coords, "angles": angles,
                "passes_angle_policy": in_policy}
    except (ValueError, RuntimeError) as error:
        return {"status": "error", "error_type": type(error).__name__, "error": str(error)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", required=True, help="local wheel used for the installed package")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    extension = importlib.util.find_spec("chematic.chematic").origin
    records = []
    for smiles in CASES:
        mol = chematic.from_smiles(smiles).add_hydrogens()
        for seed in (1, 42):
            records.append({"smiles": smiles, "seed": seed,
                            "baseline": run_one(mol, seed, False),
                            "constrained": run_one(mol, seed, True)})
    report = {"selection": "6 synthetic Pt examples and 2 previously exposed development cases; not a holdout",
              "platform": platform.platform(), "python": sys.version,
              "wheel_sha256": digest(args.wheel), "extension_sha256": digest(extension),
              "script_sha256": digest(__file__), "records": records}
    Path(args.output).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
