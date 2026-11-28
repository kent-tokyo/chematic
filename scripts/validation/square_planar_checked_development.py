#!/usr/bin/env python3
"""Screen exposed square-planar attempts and independently remeasure planarity."""

import argparse
import hashlib
import importlib.util
import itertools
import json
import platform
import sys
from pathlib import Path

import chematic
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def shape_oracle(report, coordinates):
    origin = np.asarray(coordinates[report["center"]], dtype=float)
    vectors = [np.asarray(coordinates[index], dtype=float) - origin
               for index in report["neighbors"] if index is not None]
    units = [vector / np.linalg.norm(vector) for vector in vectors]
    return max(abs(float(np.linalg.det(np.stack(triple))))
               for triple in itertools.combinations(units, 3))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--wheel", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    source = json.loads(Path(args.input).read_text(encoding="utf-8"))
    records = []
    for old in source["records"]:
        mol = chematic.from_smiles(old["smiles"]).add_hydrogens()
        row = {"group": old["group"], "smiles": old["smiles"], "seed": old["seed"]}
        try:
            retry = mol.rdkit_embed_square_planar_with_retry(random_seed=old["seed"])
            row["retry"] = {"status": "returned", **retry}
            prior = old["retry"]
            row["matches_previous"] = (prior["status"] == "returned" and
                                       retry["coordinates"] == prior["coordinates"] and
                                       retry["attempts"] == prior["attempts"])
            reports = mol.square_planar_geometry_diagnostics(retry["coordinates"])
            row["planarity"] = [{"center": report["center"],
                                  "native": report["max_unit_triple_product"],
                                  "numpy": shape_oracle(report, retry["coordinates"])}
                                 for report in reports]
        except (RuntimeError, ValueError) as error:
            row["retry"] = {"status": "error", "error": str(error)}
            row["matches_previous"] = (old["retry"]["status"] == "error" and
                                       str(error) == old["retry"]["error"])
            row["planarity"] = []
        try:
            checked = mol.rdkit_embed_square_planar_checked(
                random_seed=old["seed"], max_unit_triple_product=0.3)
            row["checked"] = {"status": "returned", **checked}
            row["checked_coordinates_match"] = (row["retry"]["status"] == "returned" and
                                                 checked["coordinates"] == row["retry"]["coordinates"])
        except (RuntimeError, ValueError) as error:
            message = str(error)
            category = ("planarity" if "unit triple product" in message else
                        "other_stereo" if "other declared stereo rejected" in message else "embedding")
            row["checked"] = {"status": "error", "category": category, "error": message}
        records.append(row)
    extension = importlib.util.find_spec("chematic.chematic").origin
    report = {"selection": "previously exposed development and authored synthetic inputs; not a holdout",
              "planarity_limit": 0.3, "platform": platform.platform(), "python": sys.version,
              "numpy_version": np.__version__, "input_sha256": digest(args.input),
              "wheel_sha256": digest(args.wheel), "extension_sha256": digest(extension),
              "script_sha256": digest(__file__), "records": records}
    Path(args.output).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
