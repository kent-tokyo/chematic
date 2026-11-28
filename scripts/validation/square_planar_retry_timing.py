#!/usr/bin/env python3
"""Interleaved local timing of square-planar single and retry calls."""
import argparse
import hashlib
import importlib.util
import json
import platform
import sys
import time
from pathlib import Path

import chematic
from square_planar_retry_development import DEVELOPMENT, SYNTHETIC


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    cases = []
    for group, inputs, seeds in (("development", DEVELOPMENT, (1, 42)),
                                 ("synthetic", SYNTHETIC, (7, 23))):
        for smiles in inputs:
            mol = chematic.from_smiles(smiles).add_hydrogens()
            cases.extend((group, smiles, seed, mol) for seed in seeds)
    rows = []
    for block, method in enumerate(("single", "retry", "retry", "single"), 1):
        for group, smiles, seed, mol in cases:
            start = time.perf_counter_ns()
            try:
                result = (mol.rdkit_embed_square_planar_with_retry(random_seed=seed) if method == "retry"
                          else mol.rdkit_embed_square_planar(random_seed=seed))
                status = "returned"
                attempts = result["attempts"] if method == "retry" else 1
            except (ValueError, RuntimeError) as error:
                status = "error"
                attempts = 2 if method == "retry" and "tightened-bounds retry:" in str(error) else 1
            elapsed_ns = time.perf_counter_ns() - start
            rows.append({"block": block, "method": method, "group": group, "smiles": smiles,
                         "seed": seed, "status": status, "attempts": attempts, "elapsed_ns": elapsed_ns})
    extension = importlib.util.find_spec("chematic.chematic").origin
    result = {"selection": "exposed development and authored synthetic molecules; neither is a holdout",
              "platform": platform.platform(), "python": sys.version,
              "wheel_sha256": digest(args.wheel), "extension_sha256": digest(extension),
              "script_sha256": digest(__file__),
              "case_script_sha256": digest(Path(__file__).with_name("square_planar_retry_development.py")),
              "records": rows}
    Path(args.output).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
