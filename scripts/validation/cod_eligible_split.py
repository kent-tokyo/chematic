#!/usr/bin/env python3
"""Restrict a COD holdout to the independently audited common MMFF domain."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("split", type=Path)
    p.add_argument("audit", type=Path)
    p.add_argument("support", type=Path)
    p.add_argument("--output-split", type=Path, required=True)
    p.add_argument("--output-audit", type=Path, required=True)
    args = p.parse_args()
    split_bytes, audit_bytes, support_bytes = args.split.read_bytes(), args.audit.read_bytes(), args.support.read_bytes()
    split, audit, support = json.loads(split_bytes), json.loads(audit_bytes), json.loads(support_bytes)
    split_hash = sha(split_bytes)
    if audit["split_sha256"] != split_hash or audit["source_sha256"] != split["source_sha256"] or audit["overlaps"]:
        raise ValueError("original exclusion audit is not clean or does not match")
    if support["split_sha256"] != split_hash or support["source_sha256"] != split["source_sha256"] or support["partition"] != "holdout":
        raise ValueError("support audit does not match the holdout")
    expected = {row["index"]: row for row in split["rows"] if row.get("split") == "holdout"}
    observed = {row["index"]: row for row in support["results"]}
    if set(expected) != set(observed) or len(observed) != len(support["results"]):
        raise ValueError("support audit is incomplete or duplicated")
    for index, original in expected.items():
        row = observed[index]
        if row["name"] != original["name"]:
            raise ValueError("support audit record identity mismatch")
        eligible = row.get("rdkit_mmff") is True and row.get("chematic_mmff") is True and row.get("graph_matches") is True
        if not eligible:
            original["split"] = "excluded_noncommon_mmff_or_graph"
    split["counts"] = dict(Counter(r.get("split", r.get("status", "unknown")) for r in split["rows"]))
    split["eligibility"] = {"rule": "both engines construct MMFF and graphs match in the independent support audit", "support_sha256": sha(support_bytes), "full_holdout_records": len(expected), "eligible_holdout_records": sum(row.get("split") == "holdout" for row in expected.values())}
    output_bytes = (json.dumps(split, indent=2) + "\n").encode()
    args.output_split.write_bytes(output_bytes)
    # No new record can be introduced by filtering an already clean audit.
    audit["split_sha256"] = sha(output_bytes)
    audit["eligible_scaffolds"] = len({r["scaffold_sha256"] for r in split["rows"] if r.get("split") in ("development", "holdout")})
    audit["derived_from_audit_sha256"] = sha(audit_bytes)
    audit["rule"] += " Filtered to a subset using the independently recorded support audit."
    args.output_audit.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(split["eligibility"], sort_keys=True))


if __name__ == "__main__":
    main()
