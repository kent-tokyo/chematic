#!/usr/bin/env python3
"""Check the published-npm 2026.03.6 -> 2026.09.1 rebaseline packet.

This does not turn the missing 2026.09.1 Python/native lanes into successes.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation" / "results"
DAY = "2026-10-04"
CORPUS = ROOT / "validation" / "benchmark_corpora" / "rdkit-js-browser-10k-v1.smi"
EXECUTION = ROOT / "validation" / "rdkit_rebaseline_execution.json"
QUERIES = ROOT / "validation" / "rdkit_rebaseline_smarts_queries.json"
VERSIONS = ("2026.03.6", "2026.09.1")


def load(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def check_python_baseline(corpus_hash: str) -> None:
    """Verify the distributed 2026.03.6 Python baseline without importing a wheel."""
    prefix = "rdkit-rebaseline-python"
    provenance = load(f"{prefix}-provenance-v1.0.33-vs-2026.03.6-{DAY}.json")
    require(
        provenance["source"]["source_commit"]
        == "35f96e67a899992fec53e006d96f4532a4019da2"
        and provenance["source"]["dirty_worktree"] is False,
        "Python baseline source was not the corrected clean CIP comparator",
    )
    require(
        provenance["inputs"]["sealed_accuracy_cohort_reused"] is False,
        "Python baseline used sealed rows",
    )
    require(
        provenance["inputs"]["corpora"][0]["sha256"] == corpus_hash,
        "Python baseline corpus changed",
    )
    require(
        provenance["inputs"]["operation_configurations"][0]["sha256"]
        == sha256(EXECUTION.read_bytes()),
        "Python baseline command catalog changed",
    )
    require(len(provenance["lanes"]) == 1, "Python baseline has multiple oracle lanes")
    lane = provenance["lanes"][0]
    require(lane["missing_dimensions"] == [], "Python baseline provenance incomplete")
    require(
        lane["runtime"]["version"] == "2026.03.6", "Python baseline runtime changed"
    )
    require(
        lane["wrapper_backend"]["value"] == "boost_python",
        "Python baseline wrapper changed",
    )
    wheel_sha = lane["package"]["artifact"]["sha256"]
    require(
        wheel_sha == "e16c467cb254a223e59a0cf81358c6b39da15a99d2909170d693e95778fddb41",
        "Python baseline wheel changed",
    )

    contract = load(f"{prefix}-contract-vs-2026.03.6-{DAY}.json")
    require(contract["gate_passed"] is True, "Python baseline contract gate failed")
    require(
        contract["corpus"]["sha256"] == corpus_hash,
        "Python baseline contract corpus changed",
    )
    require(
        contract["package"]["artifact"]["sha256"] == wheel_sha,
        "Python baseline contract wheel changed",
    )
    require(
        contract["runtime"]["backend"]["value"] == "boost_python",
        "Python baseline contract backend changed",
    )
    require(
        contract["row_accounting"]["success_count"] == 10000,
        "Python baseline contract incomplete",
    )
    for operation in ("scalar_parse", "scalar_morgan", "batch_morgan"):
        timing = contract["timing"][operation]
        require(
            timing["repetitions"] == len(timing["samples_ns"]) == 7,
            f"Python baseline {operation} timing incomplete",
        )
        require(
            timing["result_count"] == 10000,
            f"Python baseline {operation} returned fewer rows",
        )

    chemistry = load(f"{prefix}-chemistry-v1.0.33-vs-2026.03.6-{DAY}.json")
    require(chemistry["gate_passed"] is True, "Python baseline chemistry gate failed")
    require(
        chemistry["rdkit_version"] == "2026.03.6"
        and chemistry["chematic_version"] == "1.0.33",
        "Python baseline chemistry versions changed",
    )
    require(
        chemistry["corpus"]["sha256"] == corpus_hash,
        "Python baseline chemistry corpus changed",
    )
    require(
        chemistry["queries"]["sha256"] == sha256(QUERIES.read_bytes()),
        "Python baseline queries changed",
    )
    require(
        chemistry["row_accounting"]
        == {"input_count": 10000, "completed_count": 10000, "parse_failure_count": 0},
        "Python baseline row accounting changed",
    )
    rows_path = ROOT / chemistry["rows"]["path"]
    compressed_rows = rows_path.read_bytes()
    require(
        sha256(compressed_rows) == chemistry["rows"]["compressed_sha256"],
        "Python baseline compressed rows checksum changed",
    )
    rows_body = gzip.decompress(compressed_rows)
    require(
        sha256(rows_body) == chemistry["rows"]["sha256"],
        "Python baseline rows checksum changed",
    )
    counts: Counter[str] = Counter()
    cip_abstentions: list[int] = []
    cip_mismatches: list[int] = []
    morgan_refusals: list[int] = []
    for index, line in enumerate(rows_body.splitlines()):
        row = json.loads(line)
        require(
            row["input_index"] == index and row["status"] == "completed",
            f"Python baseline row {index} incomplete",
        )
        require(
            row["smiles_parse_write"]["semantic_roundtrip"] is True,
            f"Python baseline row {index} lost structure",
        )
        require(
            row["cip"]["index_correspondence"]
            == {"atom_order": True, "bond_endpoints": True},
            f"Python baseline row {index} CIP index unproven",
        )
        require(
            row["smarts"]["query_count"] == 31,
            f"Python baseline row {index} SMARTS incomplete",
        )
        counts["smarts_differences"] += row["smarts"]["difference_count"]
        if row["cip"]["exact"]:
            counts["cip_exact"] += 1
        elif row["cip"]["chematic_unresolved"]:
            cip_abstentions.append(index)
        else:
            cip_mismatches.append(index)
            require(
                row["cip"]["rdkit_atoms"] == row["cip"]["chematic_atoms"]
                and row["cip"]["rdkit_bonds"] != row["cip"]["chematic_bonds"],
                f"Python baseline row {index} has an unexpected CIP mismatch",
            )
        if row["morgan"]["exact"]:
            counts["morgan_exact"] += 1
        else:
            morgan_refusals.append(index)
            require(
                "unsupported" in row["morgan"]["chematic_error"].lower(),
                f"Python baseline row {index} Morgan mismatch not typed",
            )
        if row["smiles_parse_write"]["exact_spelling"]:
            counts["smiles_exact"] += 1
    require(index + 1 == 10000, "Python baseline rows truncated")
    require(
        counts
        == {
            "smarts_differences": 200,
            "cip_exact": 9989,
            "morgan_exact": 9999,
            "smiles_exact": 64,
        },
        "Python baseline derived counts changed",
    )
    require(
        cip_abstentions == [2960, 4419, 4439, 4643, 4644],
        "Python baseline CIP abstentions changed",
    )
    require(
        cip_mismatches == [1206, 1213, 1214, 1287, 1370, 1371],
        "Python baseline imine E/Z mismatches changed",
    )
    require(morgan_refusals == [8341], "Python baseline Morgan refusal changed")
    require(
        chemistry["counts"]["smarts_differences"] == 200
        and chemistry["counts"]["cip_exact"] == 9989
        and chemistry["counts"]["cip_difference"] == 11
        and chemistry["counts"]["morgan_exact"] == 9999,
        "Python baseline summary and rows disagree",
    )


def check() -> None:
    corpus_hash = sha256(CORPUS.read_bytes())
    check_python_baseline(corpus_hash)
    availability = load(f"rdkit-2026-09-1-artifact-availability-{DAY}.json")
    require(availability["github_release"]["http_status"] == 200, "release missing")
    require(
        availability["pypi"]["http_status"] == 404,
        "Python availability snapshot changed",
    )
    require(
        availability["npm"]["http_status"] == 200, "npm availability snapshot changed"
    )

    artifacts: dict[str, dict] = {}
    cip_rows_digests: dict[str, str] = {}
    for version in VERSIONS:
        provenance = load(
            f"rdkit-rebaseline-npm-provenance-v1.0.33-vs-{version}-{DAY}.json"
        )
        require(
            provenance["inputs"]["sealed_accuracy_cohort_reused"] is False,
            "sealed corpus used",
        )
        require(
            provenance["inputs"]["corpora"][0]["sha256"] == corpus_hash,
            "provenance corpus changed",
        )
        require(
            provenance["inputs"]["operation_configurations"][0]["sha256"]
            == sha256(EXECUTION.read_bytes()),
            f"{version}: execution catalog changed",
        )
        require(len(provenance["lanes"]) == 1, "expected one provenance lane")
        lane = provenance["lanes"][0]
        require(lane["runtime"]["version"] == version, f"{version}: runtime mismatch")
        require(lane["missing_dimensions"] == [], f"{version}: missing provenance")
        require(
            lane["package"]["registry"]["status"] == "recorded",
            f"{version}: registry lock missing",
        )
        require(
            lane["package"]["registry"]["integrity"].startswith("sha512-"),
            f"{version}: no SRI",
        )
        artifact = lane["package"]["artifact"]
        require(
            artifact["sha256"] == lane["package"]["artifact_archive_sha256"],
            f"{version}: tarball hash mismatch",
        )
        artifacts[version] = lane

        parity = load(
            f"rdkit-rebaseline-npm-morgan-parity-v1.0.33-vs-{version}-{DAY}.json"
        )
        result = parity["result"]
        require(
            result["rdkit_version"] == version, f"{version}: parity runtime mismatch"
        )
        require(result["compared_rows"] == 10000, f"{version}: partial parity corpus")
        require(
            (
                result["exact_matches"],
                result["supported_rows"],
                result["unsupported_rows"],
            )
            == (9999, 9999, 1),
            f"{version}: Morgan counts changed",
        )
        require(
            result["unsupported_sample"][0]["index"] == 8341,
            f"{version}: refusal row changed",
        )
        require(result["first_mismatch"] is None, f"{version}: Morgan mismatch")

        browser = load(
            f"rdkit-rebaseline-npm-browser-runtime-v1.0.33-vs-{version}-{DAY}.json"
        )
        require(
            browser["corpus"]["sha256"] == corpus_hash,
            f"{version}: browser corpus changed",
        )
        require(
            browser["configuration"]["rows"] == 10000,
            f"{version}: partial browser corpus",
        )
        require(
            browser["configuration"]["repetitions"] == 20,
            f"{version}: browser repetitions changed",
        )
        require(
            browser["configuration"]["prepared_mode"] == "split",
            f"{version}: prepared operation changed",
        )
        require(
            len(browser["raw_runs"]["rdkit"]) == 20,
            f"{version}: RDKit raw runs missing",
        )
        require(
            len(browser["raw_runs"]["chematic"]) == 20,
            f"{version}: CheMatic raw runs missing",
        )
        wasm = next(
            asset
            for asset in lane["package"]["assets"]
            if asset["path"] == "dist/RDKit_minimal.wasm"
        )
        require(
            wasm["sha256"] == browser["artifacts"]["rdkit_wasm"]["sha256"],
            f"{version}: WASM artifact changed",
        )
        require(
            browser["artifacts"]["schematic_package"]["version"] == "1.0.33",
            f"{version}: CheMatic package changed",
        )
        cip = load(f"rdkit-rebaseline-npm-cip-parity-v1.0.33-vs-{version}-{DAY}.json")
        require(cip["rdkit"]["version"] == version, f"{version}: CIP runtime mismatch")
        require(
            cip["chematic"]["package"] == "1.0.33",
            f"{version}: CIP CheMatic package mismatch",
        )
        require(
            cip["rdkit"]["wasm_sha256"] == wasm["sha256"],
            f"{version}: CIP RDKit WASM mismatch",
        )
        require(
            cip["chematic"]["wasm_sha256"]
            == browser["artifacts"]["schematic_wasm"]["sha256"],
            f"{version}: CIP CheMatic WASM mismatch",
        )
        require(
            cip["corpus"]["sha256"] == corpus_hash, f"{version}: CIP corpus changed"
        )
        require(
            cip["counts"]
            == {
                "rows": 10000,
                "graph_correspondence": 9999,
                "exact_rows": 9988,
                "typed_abstention_rows": 5,
                "mismatch_rows": 6,
                "unproven_rows": 1,
            },
            f"{version}: CIP counts changed",
        )
        cip_compressed = (ROOT / cip["rows_output"]["path"]).read_bytes()
        require(
            sha256(cip_compressed) == cip["rows_output"]["compressed_sha256"],
            f"{version}: CIP compressed rows changed",
        )
        cip_body = gzip.decompress(cip_compressed)
        cip_digest = sha256(cip_body)
        require(
            cip_digest == cip["rows_output"]["uncompressed_sha256"],
            f"{version}: CIP rows changed",
        )
        cip_rows_digests[version] = cip_digest
        cip_status: Counter[str] = Counter()
        cip_nonexact: dict[str, list[int]] = {}
        for row_index, line in enumerate(cip_body.splitlines()):
            row = json.loads(line)
            require(
                row["input_index"] == row_index, f"{version}: CIP row order changed"
            )
            cip_status[row["status"]] += 1
            cip_nonexact.setdefault(row["status"], []).append(row_index)
        require(row_index + 1 == 10000, f"{version}: CIP rows incomplete")
        require(
            cip_status
            == {"exact": 9988, "typed_abstention": 5, "mismatch": 6, "unproven": 1},
            f"{version}: CIP row status count changed",
        )
        require(
            cip_nonexact["mismatch"] == [1206, 1213, 1214, 1287, 1370, 1371],
            f"{version}: CIP mismatch identities changed",
        )
        require(
            cip_nonexact["typed_abstention"] == [2960, 4419, 4439, 4643, 4644],
            f"{version}: CIP abstention identities changed",
        )
        require(
            cip_nonexact["unproven"] == [8341],
            f"{version}: CIP unproven identity changed",
        )

        smarts = load(
            f"rdkit-rebaseline-npm-smarts-parity-v1.0.33-vs-{version}-{DAY}.json"
        )
        require(
            smarts["rdkit"]["version"] == version, f"{version}: SMARTS runtime mismatch"
        )
        require(
            smarts["rdkit"]["wasm_sha256"] == wasm["sha256"],
            f"{version}: SMARTS RDKit WASM mismatch",
        )
        require(
            smarts["chematic"]["wasm_sha256"]
            == browser["artifacts"]["schematic_wasm"]["sha256"],
            f"{version}: SMARTS CheMatic WASM mismatch",
        )
        require(
            smarts["corpus"]["sha256"] == corpus_hash,
            f"{version}: SMARTS corpus changed",
        )
        require(
            smarts["queries"]["sha256"] == sha256(QUERIES.read_bytes()),
            f"{version}: SMARTS queries changed",
        )
        require(
            smarts["correspondence"]["rows_sha256"] == cip_digest,
            f"{version}: SMARTS correspondence proof changed",
        )
        expected_mismatches = 194 if version == VERSIONS[0] else 195
        require(
            smarts["counts"]
            == {
                "rows": 10000,
                "queries": 31,
                "total_cells": 310000,
                "exact_cells": 309969 - expected_mismatches,
                "mismatch_cells": expected_mismatches,
                "typed_refusal_cells": 0,
                "error_cells": 0,
                "unproven_cells": 31,
            },
            f"{version}: SMARTS counts changed",
        )
        smarts_compressed = (ROOT / smarts["rows_output"]["path"]).read_bytes()
        require(
            sha256(smarts_compressed) == smarts["rows_output"]["compressed_sha256"],
            f"{version}: SMARTS compressed rows changed",
        )
        smarts_body = gzip.decompress(smarts_compressed)
        require(
            sha256(smarts_body) == smarts["rows_output"]["uncompressed_sha256"],
            f"{version}: SMARTS rows changed",
        )
        smarts_mismatch_queries: Counter[str] = Counter()
        unproven_rows = []
        for row_index, line in enumerate(smarts_body.splitlines()):
            row = json.loads(line)
            require(
                row["input_index"] == row_index, f"{version}: SMARTS row order changed"
            )
            if row["status"] == "unproven_index_correspondence":
                unproven_rows.append(row_index)
            for difference in row["differences"]:
                require(
                    difference["status"] == "mismatch",
                    f"{version}: unexpected SMARTS cell outcome",
                )
                smarts_mismatch_queries[difference["query"]] += 1
        require(row_index + 1 == 10000, f"{version}: SMARTS rows incomplete")
        require(unproven_rows == [8341], f"{version}: SMARTS unproven row changed")
        expected_r3 = 62 if version == VERSIONS[0] else 63
        require(
            smarts_mismatch_queries
            == {
                "[R1]": 58,
                "[R2]": 63,
                "[R3]": expected_r3,
                "[k5]": 1,
                "[k6]": 10,
            },
            f"{version}: SMARTS residual families changed",
        )
    old_browser = load(
        f"rdkit-rebaseline-npm-browser-runtime-v1.0.33-vs-{VERSIONS[0]}-{DAY}.json"
    )
    new_browser = load(
        f"rdkit-rebaseline-npm-browser-runtime-v1.0.33-vs-{VERSIONS[1]}-{DAY}.json"
    )
    require(
        old_browser["artifacts"]["schematic_wasm"]["sha256"]
        == new_browser["artifacts"]["schematic_wasm"]["sha256"],
        "CheMatic WASM differs between arms",
    )
    require(
        cip_rows_digests[VERSIONS[0]] == cip_rows_digests[VERSIONS[1]],
        "CheMatic/RDKit CIP row outcomes differ by RDKit version",
    )

    comparison = load(
        f"rdkit-rebaseline-npm-oracle-delta-{VERSIONS[0]}-to-{VERSIONS[1]}-{DAY}.json"
    )
    counts = comparison["counts"]
    require(
        comparison["old"]["runtime_version"] == VERSIONS[0],
        "old oracle runtime mismatch",
    )
    require(
        comparison["new"]["runtime_version"] == VERSIONS[1],
        "new oracle runtime mismatch",
    )
    require(comparison["corpus"]["sha256"] == corpus_hash, "oracle corpus changed")
    require(
        (counts["rows"], counts["query_count"], counts["smarts_cells_compared"])
        == (10000, 31, 310000),
        "oracle comparison incomplete",
    )
    require(
        (
            counts["parse_status_differences"],
            counts["canonical_spelling_differences"],
            counts["morgan_bit_differences"],
            counts["smarts_cell_differences"],
        )
        == (0, 0, 0, 12),
        "oracle delta counts changed",
    )
    require(
        (
            counts["graph_json_differences"],
            counts["cip_tag_differences"],
            counts["cip_comparable_rows"],
        )
        == (0, 0, 10000),
        "CIP/index-correspondence evidence changed",
    )
    require(
        comparison["old"]["query_parse_failures"]
        == comparison["new"]["query_parse_failures"]
        == [],
        "query parser failures",
    )
    for version, name in zip(VERSIONS, ("old", "new")):
        package = artifacts[version]["package"]
        require(
            comparison[name]["package_json"]["sha256"]
            == package["package_json"]["sha256"],
            f"{version}: oracle package changed",
        )
        wasm = next(
            asset
            for asset in package["assets"]
            if asset["path"] == "dist/RDKit_minimal.wasm"
        )
        require(
            comparison[name]["wasm"]["sha256"] == wasm["sha256"],
            f"{version}: oracle WASM changed",
        )

    compressed = (ROOT / comparison["rows_output"]["path"]).read_bytes()
    require(
        sha256(compressed) == comparison["rows_output"]["compressed_sha256"],
        "compressed rows hash changed",
    )
    body = gzip.decompress(compressed)
    require(
        sha256(body) == comparison["rows_output"]["uncompressed_sha256"],
        "rows content hash changed",
    )
    status: Counter[tuple[str, str]] = Counter()
    changed: Counter[str] = Counter()
    for index, line in enumerate(body.splitlines()):
        row = json.loads(line)
        require(row["input_index"] == index, f"row {index}: index mismatch")
        require(
            row["old_graph_sha256"] == row["new_graph_sha256"],
            f"row {index}: graph correspondence mismatch",
        )
        require(
            row["old_cip_sha256"] == row["new_cip_sha256"],
            f"row {index}: CIP tag mismatch",
        )
        status[(row["old_status"], row["new_status"])] += 1
        for difference in row["differences"]:
            require(
                difference["operation"] == "smarts",
                f"row {index}: unexpected difference",
            )
            changed[difference["query"]] += 1
    require(index + 1 == 10000, "row file incomplete")
    require(status == {("ok", "ok"): 10000}, "row parse statuses changed")
    require(changed == {"[R2]": 6, "[R3]": 6}, "SMARTS difference families changed")
    print(
        "RDKit rebaseline evidence OK: published Python 2026.03.6 baseline; npm 2026.03.6/2026.09.1 10k, CIP, SMARTS, and browser lanes"
    )


if __name__ == "__main__":
    try:
        check()
    except (OSError, KeyError, IndexError, TypeError, ValueError) as exc:
        raise SystemExit(f"RDKit 2026.09.1 rebaseline evidence invalid: {exc}") from exc
