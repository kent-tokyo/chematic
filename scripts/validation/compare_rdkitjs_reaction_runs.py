#!/usr/bin/env python3
"""Compare two ``rdkitjs_reaction_corpus.mjs`` runs (e.g. RDKit.js 2026.03.6
and 2026.09.1 on the same rules and reactants).

A row is (rule, reactant index, hydrogen mode) with at least one raw
product in either run. Reports how many rows give the same sanitized
product sets, which differ, and the rules that ran in one build only.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def load(path: Path):
    lines = path.read_text().splitlines()
    header = json.loads(lines[0])
    rows, errors, parse_errors = {}, Counter(), set()
    for line in lines[1:]:
        r = json.loads(line)
        if "seconds" in r or r.get("start"):
            continue
        if r.get("error") == "parse":
            parse_errors.add(r["rule"])
            continue
        key = (r["rule"], r["r"], r["mode"])
        if "error" in r:
            errors[r["rule"]] += 1
            rows[key] = ("error",)
            continue
        rows[key] = (r["raw"], tuple(r["sets"]))
    return header, rows, errors, parse_errors


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", type=Path, required=True)
    ap.add_argument("--new", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    old_h, old, old_err, old_parse = load(args.old)
    new_h, new, new_err, new_parse = load(args.new)
    counts = Counter()
    differing = []
    for key in sorted(set(old) | set(new)):
        a, b = old.get(key, (0, ())), new.get(key, (0, ()))
        if a == b:
            counts["same"] += 1
        elif a[1:] == b[1:] and a != ("error",) and b != ("error",):
            counts["same_sets_raw_count_differs"] += 1
        else:
            counts["differ"] += 1
            differing.append({"rule": key[0], "reactant": key[1], "mode": key[2],
                              "old": list(a), "new": list(b)})
    report = {
        "schema": "rdkitjs-reaction-differential/v1",
        "old": old_h,
        "new": new_h,
        "rows_with_raw_products": dict(counts),
        "differing_by_rule": Counter(d["rule"] for d in differing).most_common(),
        "parse_errors": {"old_only": sorted(old_parse - new_parse), "new_only": sorted(new_parse - old_parse),
                         "both": len(old_parse & new_parse)},
        "run_errors": {"old": sum(old_err.values()), "new": sum(new_err.values())},
        "differing": differing,
    }
    args.output.write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps({k: report[k] for k in ("rows_with_raw_products", "parse_errors", "run_errors")}))
    print(report["differing_by_rule"][:20])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
