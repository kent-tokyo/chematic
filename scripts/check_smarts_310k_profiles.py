#!/usr/bin/env python3
"""Compare source default and opt-in SMARTS match sets and booleans on 310k cells.

The archived RDKit oracle, published Python rows and adjudicated 200-cell
match-set list are immutable inputs. Only 43 of those 200 cells change the
hit/no-hit Boolean. The source default must reproduce those 43 exact cell
identities, and the source default must reproduce all 200 match-set cell
identities, before any opt-in effect can be attributed to it.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ORACLE = ROOT / "validation/results/v1.0.30-rdkit-smarts-all-cells.jsonl.gz"
QUERIES = ROOT / "validation/rdkit_rebaseline_smarts_queries.json"
CLASSIFICATION = ROOT / "validation/results/rdkit-rebaseline-residual-classification-v1.0.27-issue634-vs-2026.03.6-2026-09-28.json"
PUBLISHED_ROWS = ROOT / "validation/results/v1.0.30-published-python-chemistry-rows.jsonl.gz"
PUBLISHED_ROWS_SHA256 = "4c429609828c193d84bde7200e0e8b36eaa91e7c01a8aaf4f6ebe73839412452"
PROFILES = ("parity", "shared_symmetrized")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_jsonl(path: Path):
    with path.open(encoding="utf-8") as stream:
        yield from (json.loads(line) for line in stream if line.strip())


def normalized_match_sets(value: object) -> frozenset[tuple[int, ...]]:
    if not isinstance(value, list) or any(
        not isinstance(match, list) or any(not isinstance(atom, int) for atom in match)
        for match in value
    ):
        raise ValueError("SMARTS match sets must be lists of atom-index lists")
    return frozenset(tuple(sorted(match)) for match in value)


def profile_matches(entry: dict, profile: str) -> list | None:
    if profile == "parity":
        if "parity_error" in entry or entry.get("parity_budget_exhausted"):
            return None
        matches = entry.get("parity")
    else:
        shared = entry.get("shared_symmetrized")
        if not isinstance(shared, dict):
            raise ValueError("shared SMARTS profile is missing")
        if "error" in shared or shared.get("budget_exhausted"):
            return None
        matches = shared.get("matches")
    if not isinstance(matches, list):
        raise ValueError(f"{profile} SMARTS profile has neither matches nor a refusal")
    return matches


def published_boolean_cells(path: Path, expected_match_set_cells: set[tuple[int, str]]) -> set[tuple[int, str]]:
    if sha256(path) != PUBLISHED_ROWS_SHA256:
        raise ValueError("published v1.0.30 row archive changed")
    match_set_cells = set()
    boolean_cells = set()
    row_count = 0
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for index, line in enumerate(stream):
            row_count += 1
            row = json.loads(line)
            if row["input_index"] != index:
                raise ValueError("published row order changed")
            for difference in row["smarts"]["differences"]:
                cell = (index, difference["query"])
                if cell in match_set_cells:
                    raise ValueError(f"duplicate published SMARTS difference: {cell}")
                match_set_cells.add(cell)
                if bool(difference["rdkit"]) != bool(difference["chematic"]):
                    boolean_cells.add(cell)
    if row_count != 10_000 or match_set_cells != expected_match_set_cells or len(boolean_cells) != 43:
        raise ValueError("published match-set/Boolean residual accounting changed")
    return boolean_cells


def audit(
    dump_path: Path,
    oracle_path: Path,
    queries_path: Path,
    classification_path: Path,
    published_rows_path: Path,
) -> dict:
    queries_document = json.loads(queries_path.read_text(encoding="utf-8"))
    queries = queries_document["queries"]
    if len(queries) != 31 or len(set(queries)) != 31:
        raise ValueError("expected 31 unique pinned SMARTS queries")
    classified = json.loads(classification_path.read_text(encoding="utf-8"))
    expected_cells = {
        (row["input_index"], cell["query"])
        for row in classified["smarts"]["rows"]
        for cell in row["cells"]
    }
    if len(expected_cells) != 200 or classified["smarts"]["differing_cells"] != 200:
        raise ValueError("adjudicated baseline does not contain 200 unique cells")
    published_boolean = published_boolean_cells(published_rows_path, expected_cells)

    with gzip.open(oracle_path, "rt", encoding="utf-8") as stream:
        oracle_records = (json.loads(line) for line in stream if line.strip())
        manifest = next(oracle_records)
        if (manifest.get("rdkit_version") != "2026.03.6"
                or manifest.get("input_count") != 10_000
                or manifest.get("query_count") != len(queries)
                or manifest.get("queries_sha256") != sha256(queries_path)):
            raise ValueError("oracle manifest does not pin this comparator and query set")

        dump_records = read_jsonl(dump_path)
        default_match_differences: set[tuple[int, str]] = set()
        default_boolean_differences: set[tuple[int, str]] = set()
        profile_match_differences: dict[str, set[tuple[int, str]]] = {
            profile: set() for profile in PROFILES
        }
        profile_boolean_differences: dict[str, set[tuple[int, str]]] = {
            profile: set() for profile in PROFILES
        }
        refused: dict[str, set[tuple[int, str]]] = {
            profile: set() for profile in PROFILES
        }
        source_parse_refusals = []
        hand_rows = 0
        compared_rows = 0
        for row in dump_records:
            if row.get("record_type") == "footer":
                footer = row
                break
            if not row["id"].startswith("corpus_"):
                hand_rows += 1
                continue
            index = compared_rows
            if row["id"] != f"corpus_{index}":
                raise ValueError(f"source row order differs at {index}")
            oracle = next(oracle_records)
            if (oracle["input_index"] != index or oracle["smiles"] != row["smiles"]
                    or not isinstance(oracle["matches"], list)
                    or len(oracle["matches"]) != len(queries)):
                raise ValueError(f"source/oracle rows differ at {index}")
            if "parse_error" in row:
                source_parse_refusals.append(index)
                compared_rows += 1
                continue
            if set(row["patterns"]) != set(queries):
                raise ValueError(f"query coverage differs at row {index}")
            for position, query in enumerate(queries):
                cell = (index, query)
                expected_matches = normalized_match_sets(oracle["matches"][position])
                expected_boolean = bool(expected_matches)
                entry = row["patterns"][query]
                default = entry.get("default")
                if not isinstance(default, list):
                    raise ValueError(f"default SMARTS result absent at {cell}")
                default_matches = normalized_match_sets(default)
                if default_matches != expected_matches:
                    default_match_differences.add(cell)
                if bool(default_matches) != expected_boolean:
                    default_boolean_differences.add(cell)
                for profile in PROFILES:
                    matches = profile_matches(entry, profile)
                    if matches is None:
                        refused[profile].add(cell)
                    else:
                        profile_sets = normalized_match_sets(matches)
                        if profile_sets != expected_matches:
                            profile_match_differences[profile].add(cell)
                        if bool(profile_sets) != expected_boolean:
                            profile_boolean_differences[profile].add(cell)
            compared_rows += 1
        else:
            raise ValueError("source dump has no completed footer")
        if (not footer.get("completed") or footer.get("input_rows") != 10_000
                or footer.get("hand_corpus_rows") != hand_rows
                or footer.get("emitted_molecule_rows") != 10_000 + hand_rows
                or compared_rows != 10_000 or next(dump_records, None) is not None
                or next(oracle_records, None) is not None):
            raise ValueError("source/oracle dump completion or row accounting failed")
        source_path = Path(footer["source_path"])
        if not source_path.is_absolute():
            source_path = ROOT / source_path
        if sha256(source_path) != manifest["corpus_sha256"]:
            raise ValueError("source dump and oracle do not pin the same corpus bytes")
        if source_parse_refusals:
            raise ValueError(f"source parser refused rows: {source_parse_refusals[:5]}")
        if default_match_differences != expected_cells:
            raise ValueError(
                "source default match-set residuals differ from published v1.0.30: "
                f"missing={list(expected_cells - default_match_differences)[:5]}, "
                f"extra={list(default_match_differences - expected_cells)[:5]}"
            )
        if default_boolean_differences != published_boolean:
            raise ValueError(
                "source default Boolean residuals differ from published v1.0.30: "
                f"missing={list(published_boolean - default_boolean_differences)[:5]}, "
                f"extra={list(default_boolean_differences - published_boolean)[:5]}"
            )

    profiles = {}
    for profile in PROFILES:
        remaining_match_sets = profile_match_differences[profile]
        remaining_booleans = profile_boolean_differences[profile]
        profile_refused = refused[profile]
        fixed_match_sets = default_match_differences - remaining_match_sets - profile_refused
        fixed_booleans = default_boolean_differences - remaining_booleans - profile_refused
        profiles[profile] = {
            "typed_refusal_cells": len(profile_refused),
            "refused_by_query": dict(sorted(Counter(q for _, q in profile_refused).items())),
            "refused_cells": [[index, query] for index, query in sorted(profile_refused)],
            "refused_rows": sorted({index for index, _ in profile_refused}),
            "match_sets": {
                "wrong_confident_cells": len(remaining_match_sets),
                "fixed_source_default_cells": len(fixed_match_sets),
                "refused_source_default_wrong_cells": len(profile_refused & default_match_differences),
                "refused_source_default_correct_cells": len(profile_refused - default_match_differences),
                "fixed_by_query": dict(sorted(Counter(q for _, q in fixed_match_sets).items())),
                "regressed_cells": len(remaining_match_sets - default_match_differences),
                "remaining_by_query": dict(sorted(Counter(q for _, q in remaining_match_sets).items())),
                "remaining_cells": [[index, query] for index, query in sorted(remaining_match_sets)],
                "remaining_rows": sorted({index for index, _ in remaining_match_sets}),
            },
            "booleans": {
                "wrong_confident_cells": len(remaining_booleans),
                "fixed_source_default_cells": len(fixed_booleans),
                "refused_source_default_wrong_cells": len(profile_refused & default_boolean_differences),
                "refused_source_default_correct_cells": len(profile_refused - default_boolean_differences),
                "fixed_by_query": dict(sorted(Counter(q for _, q in fixed_booleans).items())),
                "regressed_cells": len(remaining_booleans - default_boolean_differences),
                "remaining_by_query": dict(sorted(Counter(q for _, q in remaining_booleans).items())),
                "remaining_cells": [[index, query] for index, query in sorted(remaining_booleans)],
            },
        }
    return {
        "schema": "smarts-match-set-and-boolean-profile-310k/v1",
        "scope": "source-only opt-in diagnostic, not published-package evidence",
        "rdkit_version": "2026.03.6",
        "inputs": {
            "dump_sha256": sha256(dump_path),
            "oracle_sha256": sha256(oracle_path),
            "queries_sha256": sha256(queries_path),
            "classification_sha256": sha256(classification_path),
            "published_rows_sha256": sha256(published_rows_path),
        },
        "accounting": {
            "molecules": compared_rows,
            "queries": len(queries),
            "cells": compared_rows * len(queries),
            "published_match_set_wrong_cells": len(expected_cells),
            "published_match_set_wrong_rows": len({index for index, _ in expected_cells}),
            "published_boolean_wrong_cells": len(published_boolean),
            "source_default_match_set_wrong_cells": len(default_match_differences),
            "source_default_boolean_wrong_cells": len(default_boolean_differences),
            "source_default_boolean_wrong_by_query": dict(
                sorted(Counter(q for _, q in default_boolean_differences).items())
            ),
        },
        "profiles": profiles,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--oracle", type=Path, default=ORACLE)
    parser.add_argument("--queries", type=Path, default=QUERIES)
    parser.add_argument("--classification", type=Path, default=CLASSIFICATION)
    parser.add_argument("--published-rows", type=Path, default=PUBLISHED_ROWS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.dump, args.oracle, args.queries, args.classification, args.published_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    compact = {
        "accounting": report["accounting"],
        "profiles": {
            name: {
                "typed_refusal_cells": profile["typed_refusal_cells"],
                "match_sets": {
                    key: value for key, value in profile["match_sets"].items()
                    if key not in {"remaining_cells", "remaining_rows"}
                },
                "booleans": {
                    key: value for key, value in profile["booleans"].items()
                    if key != "remaining_cells"
                },
            }
            for name, profile in report["profiles"].items()
        },
    }
    print(json.dumps(compact, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
