#!/usr/bin/env python3
"""Write `Mol.formula` for the first N corpus rows (reference for
`check_published_npm_formula_ez_json.mjs`). Unparseable rows are `null`."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=5000)
    parser.add_argument("--wheel-sha256", help="digest of the installed wheel, recorded as given")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    import chematic

    lines = args.corpus.read_text(encoding="utf-8").rstrip("\n").split("\n")
    smiles = [line.strip().split()[0] for line in lines][: args.limit]
    formulas = []
    for text in smiles:
        try:
            formulas.append(chematic.from_smiles(text).formula)
        except ValueError:
            formulas.append(None)
    report = {"artifact": {"chematic_version": importlib.metadata.version("chematic"),
                           "python": platform.python_version(), "wheel_sha256": args.wheel_sha256},
              "smiles_sha256": hashlib.sha256("\n".join(smiles).encode()).hexdigest(),
              "formulas": formulas}
    args.output.write_text(json.dumps(report) + "\n", encoding="utf-8")
    print(f"{args.output}: {len(formulas)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
