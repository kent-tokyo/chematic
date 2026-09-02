#!/usr/bin/env python3
"""Classify common-schema differences without counting unsupported cases as mismatches."""

import argparse
import json
from collections import Counter
from validate import validate

COMPARE = ("canonical_smiles", "formula", "molecular_weight", "rdkit_morgan_bits")


def load(path):
    return {row["id"]: row for row in (json.loads(line) for line in path.read_text().splitlines() if line.strip())}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("left", type=__import__("pathlib").Path)
    parser.add_argument("right", type=__import__("pathlib").Path)
    parser.add_argument("--fail-on-mismatch", action="store_true")
    args = parser.parse_args()
    validation_errors = validate(args.left) + validate(args.right)
    if validation_errors:
        print(json.dumps({"valid": False, "errors": validation_errors}, sort_keys=True, indent=2))
        raise SystemExit(2)
    left, right = load(args.left), load(args.right)
    counts = Counter()
    mismatches = []
    for ident in sorted(set(left) | set(right)):
        if ident not in left or ident not in right:
            counts["missing_record"] += 1
            continue
        for operation in COMPARE:
            a = left[ident]["operations"].get(operation, {"status": "unsupported"})
            b = right[ident]["operations"].get(operation, {"status": "unsupported"})
            if a["status"] == "unsupported" or b["status"] == "unsupported":
                counts["unsupported"] += 1
            elif a["status"] != "ok" or b["status"] != "ok":
                counts["failure"] += 1
            elif a.get("value") == b.get("value"):
                counts["match"] += 1
            else:
                counts["mismatch"] += 1
                mismatches.append({"id": ident, "operation": operation,
                                   "left": a.get("value"), "right": b.get("value")})
    report = {"valid": True, "counts": counts, "mismatches": mismatches}
    print(json.dumps(report, sort_keys=True, indent=2))
    if args.fail_on_mismatch and (counts["mismatch"] or counts["failure"] or counts["missing_record"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
