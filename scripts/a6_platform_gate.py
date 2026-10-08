#!/usr/bin/env python3
"""Gate one platform's A6 rows on the independent geometry/stereo/clash scorer (#739).

Reads the stereo-safe MMFF94 rows written by `public_package_3d_chematic.py`
and the rows `pipeline_v2_vs_rdkit_common_scorer` scored from them, and
requires, for the chosen arm:

- every input row a success (no `FinalStereoViolation` or other typed failure);
- every success independently sound, stereo-clean (no violated or
  unevaluable declared centre) and free of gross clashes.

It also writes a digest of the coordinates. With `--reference-rows` (another
platform's rows) it reports whether the coordinates are bit-identical row by
row. That comparison is reported, not gated: the gate is the quality above.

Exit status 1 when a quality check fails, unless `--report-only`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ARM = "chematic_pipeline_v2_mmff94_strict_stereo_safe"


def load_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as stream:
        return [json.loads(line) for line in stream if line.strip()]


def arm_rows(rows: list[dict], arm: str) -> list[dict]:
    return [row for row in rows if row.get("arm") == arm]


def coordinates_digest(rows: list[dict]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        payload = [row.get("tier"), row.get("name"), row.get("coords")]
        digest.update(json.dumps(payload, separators=(",", ":")).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def gate(rows: list[dict], scored: list[dict], expected_rows: int) -> dict:
    failures = [
        {"row_index": row.get("row_index"), "name": row.get("name"), "status": row.get("status"),
         "failure_cause": row.get("failure_cause")}
        for row in rows
        if row.get("status") != "success"
    ]
    by_key = {
        (row.get("tier"), row.get("name")): row
        for row in scored
        if row.get("engine") == "chematic"
    }
    unsound, stereo_bad, clashing, unscored = [], [], [], []
    for row in rows:
        if row.get("status") != "success":
            continue
        key = (row.get("tier"), row.get("name"))
        result = by_key.get(key)
        if result is None or result.get("status") != "scored":
            unscored.append(row.get("name"))
            continue
        if not result.get("independently_sound"):
            unsound.append(row.get("name"))
        stereo = result.get("stereo") or {}
        if stereo.get("violated", 0) or stereo.get("unevaluable", 0):
            stereo_bad.append(row.get("name"))
        if result.get("gross_clash_count", 0):
            clashing.append(row.get("name"))
    successes = len(rows) - len(failures)
    passed = (
        len(rows) == expected_rows
        and not failures
        and not unscored
        and not unsound
        and not stereo_bad
        and not clashing
    )
    return {
        "rows": len(rows),
        "expected_rows": expected_rows,
        "successes": successes,
        "failures": failures,
        "independently_sound": successes - len(unsound) - len(unscored),
        "stereo_clean": successes - len(stereo_bad) - len(unscored),
        "with_gross_clash": len(clashing),
        "unscored": unscored,
        "unsound": unsound,
        "stereo_violated_or_unevaluable": stereo_bad,
        "clashing": clashing,
        "passed": passed,
    }


def compare(rows: list[dict], reference: list[dict]) -> dict:
    ref = {(row.get("tier"), row.get("name")): row for row in reference}
    identical, different, missing = 0, [], []
    for row in rows:
        other = ref.get((row.get("tier"), row.get("name")))
        if other is None:
            missing.append(row.get("name"))
        elif other.get("status") == row.get("status") and other.get("coords") == row.get("coords"):
            identical += 1
        else:
            different.append(row.get("name"))
    return {
        "rows_identical": identical,
        "rows_different": different,
        "rows_missing_in_reference": missing,
        "reference_coordinates_sha256": coordinates_digest(reference),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--rows", required=True, type=Path, help="public_package_3d_chematic.py output")
    parser.add_argument("--scored", required=True, type=Path, help="pipeline_v2_vs_rdkit_common_scorer output")
    parser.add_argument("--arm", default=ARM)
    parser.add_argument("--expected-rows", type=int, default=265)
    parser.add_argument("--reference-rows", type=Path, help="another platform's rows to compare coordinates with")
    parser.add_argument("--reference-sha256", help="another platform's coordinates_sha256 to compare with")
    parser.add_argument("--scope", default="")
    parser.add_argument("--report-only", action="store_true", help="exit 0 even when a quality check fails")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows = arm_rows(load_jsonl(args.rows), args.arm)
    scored = arm_rows(load_jsonl(args.scored), args.arm)
    report = {
        "schema": "chematic.a6_platform_gate.v1",
        "issue": "https://github.com/kent-tokyo/chematic/issues/739",
        "scope": args.scope,
        "arm": args.arm,
        **gate(rows, scored, args.expected_rows),
        "coordinates_sha256": coordinates_digest(rows),
    }
    if args.reference_sha256:
        report["coordinates_equal_reference_sha256"] = report["coordinates_sha256"] == args.reference_sha256
    if args.reference_rows is not None:
        report["reference"] = compare(rows, arm_rows(load_jsonl(args.reference_rows), args.arm))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("scope", "rows", "successes", "independently_sound",
                                             "stereo_clean", "with_gross_clash", "passed")}))
    if "coordinates_equal_reference_sha256" in report:
        print(json.dumps({"coordinates_equal_reference_sha256": report["coordinates_equal_reference_sha256"]}))
    if "reference" in report:
        print(json.dumps({"rows_identical": report["reference"]["rows_identical"],
                          "rows_different": len(report["reference"]["rows_different"])}))
    return 0 if report["passed"] or args.report_only else 1


if __name__ == "__main__":
    raise SystemExit(main())
