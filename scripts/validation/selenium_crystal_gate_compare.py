#!/usr/bin/env python3
"""Evaluate a fixed development coordinate packet with one installed wheel."""

import argparse
import hashlib
import importlib.util
import json
import platform
from pathlib import Path

import chematic


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--wheel", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    data = json.loads(Path(args.input).read_text(encoding="utf-8"))
    records = []
    for row in data["results"]:
        if row["status"] != "ok":
            records.append({"index": row["index"], "status": "source_error"})
            continue
        try:
            mol = chematic.from_smiles(row["smiles"])
            coords = row["coordinates"]
            for mirror in (False, True):
                xyz = [[-x, y, z] if mirror else [x, y, z] for x, y, z in coords]
                records.append({"index": row["index"], "mirror": mirror, "status": "ok",
                                "valid": mol._timed_stereo_valid(xyz)})
        except (ValueError, RuntimeError) as error:
            records.append({"index": row["index"], "status": "error", "error": str(error)})
    result = {"input_sha256": digest(args.input), "wheel_sha256": digest(args.wheel),
              "extension_sha256": digest(importlib.util.find_spec("chematic.chematic").origin),
              "script_sha256": digest(__file__), "platform": platform.platform(), "records": records}
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
