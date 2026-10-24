#!/usr/bin/env python3
"""Score frozen explicit-H coordinates with RDKit's USR implementation."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from rdkit import rdBase
from rdkit.Chem import rdMolDescriptors
from rdkit.Geometry import Point3D


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fixture",
        type=Path,
        default=Path("validation/rdkit_3d_g1_fixtures.json"),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = args.fixture.read_bytes()
    records = json.loads(source)["fixtures"]
    rows = []
    for record in records:
        points = [
            Point3D(atom["x"], atom["y"], atom["z"])
            for atom in record["atoms"]
        ]
        moments = rdMolDescriptors.GetUSRFromDistributions(
            rdMolDescriptors.GetUSRDistributions(points)
        )
        rows.append(
            {
                "name": record["name"],
                "conf_id": record["conf_id"],
                "expected": list(moments),
            }
        )

    result = {
        "rdkit_version": rdBase.rdkitVersion,
        "source_fixture": str(args.fixture),
        "source_sha256": hashlib.sha256(source).hexdigest(),
        "method": "GetUSRFromDistributions(GetUSRDistributions(frozen explicit-H points))",
        "rows": rows,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
