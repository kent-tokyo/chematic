#!/usr/bin/env python3
"""Build an operation-level scorecard from common-schema JSONL results."""

import argparse
import json
from collections import Counter
from pathlib import Path

from validate import validate


def load(path: Path) -> dict[str, dict]:
    return {
        row["id"]: row
        for row in (json.loads(line) for line in path.read_text().splitlines() if line.strip())
    }


def parse_result(value: str) -> tuple[str, Path]:
    engine, separator, filename = value.partition("=")
    if not separator or not engine or not filename:
        raise argparse.ArgumentTypeError("result must use ENGINE=PATH")
    return engine, Path(filename)


def status_for(row: dict | None, operation: str) -> dict:
    if row is None:
        return {"status": "missing"}
    return row["operations"].get(operation, {"status": "unsupported"})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--result", action="append", required=True, type=parse_result,
        metavar="ENGINE=PATH", help="common-schema JSONL result (repeatable)",
    )
    parser.add_argument("--reference", default="rdkit", help="engine used as comparison reference")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fail-on-mismatch", action="store_true")
    args = parser.parse_args()

    paths = dict(args.result)
    if len(paths) != len(args.result):
        parser.error("each engine may be supplied only once")
    if args.reference not in paths:
        parser.error(f"reference engine {args.reference!r} was not supplied")

    errors = []
    results = {}
    for engine, path in paths.items():
        file_errors = validate(path)
        errors.extend(f"{engine}: {error}" for error in file_errors)
        if not file_errors:
            results[engine] = load(path)
    if errors:
        print(json.dumps({"valid": False, "errors": errors}, indent=2, sort_keys=True))
        raise SystemExit(2)

    rows = list(results.values())
    corpus_hashes = {row["corpus_sha256"] for records in rows for row in records.values()}
    if len(corpus_hashes) != 1:
        errors.append("result files use different corpus hashes")
        print(json.dumps({"valid": False, "errors": errors}, indent=2, sort_keys=True))
        raise SystemExit(2)

    engines = {}
    operations = sorted({
        operation
        for records in rows
        for row in records.values()
        for operation in row["operations"]
    })
    reference = results[args.reference]
    for engine, records in results.items():
        metadata = {(row["engine_version"], row.get("source_commit")) for row in records.values()}
        versions = sorted({version for version, _ in metadata})
        commits = sorted({commit for _, commit in metadata})
        engines[engine] = {
            "engine_version": versions[0],
            "source_commits": commits,
            "records": len(records),
        }

    scorecard_operations = {}
    for operation in operations:
        by_engine = {}
        for engine, records in results.items():
            counts = Counter(status_for(records.get(ident), operation)["status"]
                             for ident in sorted(set(records) | set(reference)))
            by_engine[engine] = dict(sorted(counts.items()))

        comparisons = {}
        for engine, records in results.items():
            if engine == args.reference:
                continue
            counts = Counter()
            for ident in sorted(set(records) | set(reference)):
                left = status_for(reference.get(ident), operation)
                right = status_for(records.get(ident), operation)
                if left["status"] == "ok" and right["status"] == "ok":
                    counts["match" if left.get("value") == right.get("value") else "mismatch"] += 1
                else:
                    counts["uncomparable"] += 1
            comparisons[engine] = dict(sorted(counts.items()))
        scorecard_operations[operation] = {"status_counts": by_engine, "against_reference": comparisons}

    report = {
        "valid": True,
        "schema_version": 1,
        "corpus_sha256": next(iter(corpus_hashes)),
        "corpus_records": len(reference),
        "reference_engine": args.reference,
        "engines": engines,
        "operations": scorecard_operations,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.fail_on_mismatch and any(
        counts.get("mismatch", 0)
        for operation in scorecard_operations.values()
        for counts in operation["against_reference"].values()
    ):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
