#!/usr/bin/env python3
"""Write a deterministic Markdown comparison report for two result files."""

import argparse
import json
from pathlib import Path

from score import COMPARE, load
from validate import validate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    errors = validate(args.left) + validate(args.right)
    if errors:
        raise SystemExit("invalid comparison input: " + "; ".join(errors))
    left, right = load(args.left), load(args.right)
    rows = []
    for ident in sorted(set(left) | set(right)):
        cells = []
        for operation in COMPARE:
            a = left.get(ident, {}).get("operations", {}).get(operation, {"status": "unsupported"})
            b = right.get(ident, {}).get("operations", {}).get(operation, {"status": "unsupported"})
            if a["status"] == "unsupported" or b["status"] == "unsupported":
                result = "unsupported"
            elif a["status"] != "ok" or b["status"] != "ok":
                result = "failure"
            elif a.get("value") == b.get("value"):
                result = "match"
            else:
                result = "mismatch"
            cells.append(result)
        rows.append((ident, cells))
    lines = ["# Direct comparison report", "", f"- left: `{args.left}`", f"- right: `{args.right}`", "",
             "| Corpus id | canonical_smiles | formula | molecular_weight | rdkit_morgan_bits |",
             "|---|---:|---:|---:|---:|"]
    lines.extend("| " + ident + " | " + " | ".join(cells) + " |" for ident, cells in rows)
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
