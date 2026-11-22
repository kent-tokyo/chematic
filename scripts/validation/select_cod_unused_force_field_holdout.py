#!/usr/bin/env python3
"""Reserve all untouched holdout scaffold groups for a later field candidate."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("split", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.split.read_bytes()
    if hashlib.sha256(raw).hexdigest() != "c9d185e075f83dacc48620d1b843727cbaf6729406aee0568fada9e30bfbff9d":
        raise ValueError("unexpected prior split")
    data = json.loads(raw)
    groups = {}
    for row in data["rows"]:
        if row.get("scaffold_sha256"):
            groups.setdefault(row.get("split"), set()).add(row["scaffold_sha256"])
    candidate = groups["unused_holdout"]
    for previous in ("holdout", "development", "unused_development", "excluded_known_scaffold", "excluded_invalid_exclusion_overlap"):
        if candidate & groups[previous]:
            raise ValueError(f"scaffold overlap with {previous}")
    for row in data["rows"]:
        if row.get("split") in ("development", "holdout"):
            row["split"] = "excluded_previously_evaluated"
        elif row.get("split") == "unused_holdout":
            row["split"] = "holdout"
    data["counts"] = dict(Counter(r.get("split", r.get("status", "unknown")) for r in data["rows"]))
    data["rule"] += "; reserve all prior unused_holdout scaffold groups, excluding previous development and holdout"
    data["field_followup_selection"] = {
        "input_split_sha256": hashlib.sha256(raw).hexdigest(),
        "holdout_scaffolds": len(candidate),
        "selection": "all previously unused_holdout scaffold groups; no outcome-based selection",
    }
    args.output.write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps({"counts": data["counts"], "holdout_scaffolds": len(candidate)}, sort_keys=True))


if __name__ == "__main__":
    main()
