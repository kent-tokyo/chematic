#!/usr/bin/env python3
"""Verify the exposed v1.0.30 published-artifact packet without RDKit installed.

Integrity success is not P0/P1 acceptance: the 83-reaction extension has
confident differences, and binding-wide speed/memory gates remain incomplete.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from collections import Counter
from itertools import zip_longest
from pathlib import Path

if __package__:
    from .check_published_python_version_outputs import compare
    from .check_v1030_published_chemistry_residuals import check as check_chemistry_residuals
    from .check_published_rust_63op_outputs import build_report as build_rust_operation_report
    from .check_published_wasm_paired import check as check_published_wasm_paired
    from .check_published_browser_paired import check as check_published_browser_paired
    from .check_v1030_published_browser_morgan_rows import check as check_browser_morgan_rows
else:
    from check_published_python_version_outputs import compare
    from check_v1030_published_chemistry_residuals import check as check_chemistry_residuals
    from check_published_rust_63op_outputs import build_report as build_rust_operation_report
    from check_published_wasm_paired import check as check_published_wasm_paired
    from check_published_browser_paired import check as check_published_browser_paired
    from check_v1030_published_browser_morgan_rows import check as check_browser_morgan_rows

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation/results"
WHEEL_SHA256 = "f6b1b22dd898cf3f100001156d9504b3753cb6b392143cbc0e6f3be78ec723d5"
NPM_TARBALL_SHA256 = "fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201"
NPM_V1029_TARBALL_SHA256 = "d0af78a8d6b711a13b63985d4b36079e245441e99f1712532555dc69f11569fa"
CRATE_SHA256 = "b1a897689ab5b4983e56325dbd361e8ff0625dd752984af7fe75a10286f84a82"
NPM_REACTION_ROWS_SHA256 = "ce3d801ed55f610412553861382040356ac17a8004b0bec23223bc14ddbdee9b"
ALL_REACTION_REPORT_SHA256 = "a43cf9eec38b04664dc8ee269aa69f51a9ae3e310fa5a781177d07f389e85785"
REACTION_TEMPLATE_MAP_REPORT_SHA256 = "087f38ffa4ea766698ed14c7867238bdb5fc22cdfcab3bfceb4b0ae55f861ad8"
REACTION_NEW_MAP_REPORT_SHA256 = "98b469b3d3fbb6900031313ce5bfa457285968f2df6a5b03f49d4836cb60df82"
REACTION_MULTIPLICITY_REPORT_SHA256 = "aecd38403e092f3518aa121167077dbfb4ef7637d16d42b763460a2542b10a3d"
BROWSER_PAIRED_RAW_SHA256 = "de0e7f1f8ba7ad208aeb906acb435b10ca210c05bf4ddc860a4217618c99275d"
BROWSER_PAIRED_SUMMARY_SHA256 = "3f27f1c17824aabfee7b52ae9c122974736290eb00da7249c29f5053823d2233"
PREPARED_SPLIT_RAW_SHA256 = "feaa0a80067e18e34e82c3048ef0aed626bfbe1ae4f843b2a95f6235403e2af7"
PREPARED_SPLIT_SUMMARY_SHA256 = "16cc70762572b0b953418826fadb90172529bb0fb2dee72a138df3d51cd33bd8"
PUBLISHED_BROWSER_MORGAN_ROWS = {
    (250, "direct"): "3d1428b094bdcb756566a12b27690eb369c5e6ef6a8cbdac7264ca86a02cc105",
    (250, "prepared"): "2ffcff1c190881d45a26f62ba96d001afdc68c00c2b8549099d0bb20d09892a2",
    (10_000, "direct"): "f7c404c1e32d5f4b136bff17c8b4a974f6e3bfd2057309ba68d2c557e2bfb2f6",
    (10_000, "prepared"): "96f1eb2640243e6adcf71cbf330956f8fd7e4a4a40642415195475ca9df38b58",
}
PAIRED_63OP_SHA256 = "5448fd34190d5358227663a77e415e31dca71718cfd38c2d01e156ab539719da"


def read_json(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def raw_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fail_if(condition: bool, message: str) -> None:
    if condition:
        raise ValueError(message)


def formula_composition(value: object) -> dict[str, int] | None:
    if not isinstance(value, str):
        return None
    counts: dict[str, int] = {}
    at = 0
    for token in re.finditer(r"([A-Z][a-z]?)(\d*)", value):
        if token.start() != at:
            return None
        at = token.end()
        counts[token.group(1)] = counts.get(token.group(1), 0) + int(token.group(2) or 1)
    return counts if at == len(value) and at > 0 else None


def same_rounded_coordinates(actual: object, expected: object) -> bool:
    if not isinstance(actual, list) or not isinstance(expected, list) or len(actual) != len(expected):
        return False
    return all(isinstance(point, list) and isinstance(reference, list) and
               len(point) == len(reference) == 3 and
               all(isinstance(a, (int, float)) and isinstance(b, (int, float)) and
                   math.isfinite(a) and math.isfinite(b) and abs(a - b) <= 0.000051
                   for a, b in zip(point, reference, strict=True))
               for point, reference in zip(actual, expected, strict=True))


def npm_python_outcome(op: str, actual: dict, expected: dict) -> str:
    if actual == expected:
        return "exact"
    if actual.get("input_index") != expected.get("input_index"):
        return "different"
    av, ev = actual.get("value"), expected.get("value")
    if op in {"qed", "chi1v"} and isinstance(av, float) and isinstance(ev, float) and \
            math.isfinite(av) and math.isfinite(ev) and math.isclose(av, ev, rel_tol=1e-12, abs_tol=1e-12):
        return "numeric_roundoff"
    if op == "formula" and formula_composition(av) is not None and \
            formula_composition(av) == formula_composition(ev):
        return "formula_spelling_only"
    if op == "embed_3d" and same_rounded_coordinates(av, ev):
        return "coordinate_roundoff"
    return "different"


def check_npm_operation_slice(version: str) -> list[dict]:
    summary = read_json(f"v{version}-published-npm-python-63op-slice.json")
    raw_path = RESULTS / f"v{version}-published-npm-python-63op-slice.jsonl.gz"
    python_path = RESULTS / f"v{version}-published-python-63op-outputs.jsonl.gz"
    corpus_path = ROOT / "scripts/chembl_accuracy_corpus_4999.smi"
    raw_bytes = raw_path.read_bytes()
    python_bytes = python_path.read_bytes()
    expected_tarball = NPM_TARBALL_SHA256 if version == "1.0.30" else NPM_V1029_TARBALL_SHA256
    fail_if(summary["schema"] != "published-npm-python-63op-output-slice/v2" or
            summary["npm"] != {"version": version, "tarball_sha256": expected_tarball} or
            summary["python_archive_sha256"] != hashlib.sha256(python_bytes).hexdigest() or
            summary["corpus"] != {"sha256": hashlib.sha256(corpus_path.read_bytes()).hexdigest(), "input_count": 5000} or
            summary["output"]["compressed_sha256"] != hashlib.sha256(raw_bytes).hexdigest(),
            "npm operation artifact/corpus/archive digest")
    raw = gzip.decompress(raw_bytes)
    fail_if(summary["output"]["rows_sha256"] != hashlib.sha256(raw).hexdigest(),
            "npm operation uncompressed digest")
    rows = [json.loads(line) for line in raw.splitlines()]
    baseline = [json.loads(line) for line in gzip.decompress(python_bytes).splitlines()]
    by_name = {item["op"]: item["rows"] for item in baseline}
    names = [item["op"] for item in rows]
    coverage = summary["coverage"]
    absent = coverage["no_equivalent_public_api"]
    fail_if(len(baseline) != 63 or len(by_name) != 63 or len(rows) != 59 or
            coverage["total_operations"] != 63 or coverage["adapted"] != 59 or
            sorted(absent) != sorted(("rdkit_tpsa", "atom_pair(rdkit-compatible)",
                                     "pattern_fp(rdkit-compatible)", "embed+minimize_mmff94")) or
            any(not reason for reason in absent.values()) or
            sorted(names + list(absent)) != sorted(by_name) or
            names != [item["op"] for item in summary["comparison"]],
            "npm operation names or pending accounting")
    exact_operations = 0
    non_equivalent_differences = {"mol_block_write": 2000, "svg_depiction": 20, "2d_layout": 500}
    for item, recorded in zip(rows, summary["comparison"], strict=True):
        expected_rows = by_name[item["op"]]
        fail_if(len(item["rows"]) != len(expected_rows) or not 0 < len(expected_rows) <= 5000,
                f"npm {item['op']}: row denominator")
        outcomes = Counter()
        samples = []
        for index, (actual, expected) in enumerate(zip(item["rows"], expected_rows, strict=True)):
            fail_if(actual.get("input_index") != index, f"npm {item['op']}: row order {index}")
            result = npm_python_outcome(item["op"], actual, expected)
            outcomes[result] += 1
            if result != "exact" and len(samples) < 3:
                samples.append({"input_index": index, "actual": actual, "expected": expected})
        fail_if(recorded["outcomes"] != {key: outcomes[key] for key in
                ("exact", "numeric_roundoff", "formula_spelling_only", "coordinate_roundoff", "different")} or
                recorded["samples"] != samples or recorded["input_count"] != len(expected_rows),
                f"npm {item['op']}: archived comparison")
        fail_if(outcomes["different"] != non_equivalent_differences.get(item["op"], 0),
                f"npm {item['op']}: unexpected value mismatch")
        fail_if(outcomes["coordinate_roundoff"] != (100 if item["op"] == "embed_3d" else 0),
                f"npm {item['op']}: unexpected coordinate rounding")
        if item["op"] == "tanimoto_1xN(per target, compatible Morgan)":
            fail_if(len(item["rows"]) != 50 or
                    any(len(row.get("value", [])) != 5000 for row in item["rows"]),
                    "npm compatible-Morgan similarity 50 x 5000 accounting")
        exact_operations += outcomes["exact"] == len(expected_rows)
    fail_if(exact_operations != 52, "npm operation exact-count regression")
    return rows


def check_npm_version_diff() -> None:
    before = check_npm_operation_slice("1.0.29")
    after = check_npm_operation_slice("1.0.30")
    fail_if([op["op"] for op in before] != [op["op"] for op in after],
            "npm v1.0.29/v1.0.30 operation order")
    expected = {"hba": 1359, "lipinski_bundle(mw,logp,hbd,hba)": 1359}
    python_paths = {
        version: RESULTS / f"v{version}-published-python-63op-outputs.jsonl.gz"
        for version in ("1.0.29", "1.0.30")
    }
    python_rows = {
        version: {item["op"]: item["rows"] for item in
                  (json.loads(line) for line in gzip.open(path, "rt", encoding="utf-8"))}
        for version, path in python_paths.items()
    }
    observed = {}
    for old, new in zip(before, after, strict=True):
        name = old["op"]
        fail_if(len(old["rows"]) != len(new["rows"]), f"npm {name}: version denominator")
        npm_changed = [i for i, (a, b) in enumerate(zip(old["rows"], new["rows"], strict=True)) if a != b]
        python_changed = [i for i, (a, b) in enumerate(zip(
            python_rows["1.0.29"][name], python_rows["1.0.30"][name], strict=True)) if a != b]
        fail_if(npm_changed != python_changed, f"npm {name}: change indices differ from Python")
        if npm_changed:
            observed[name] = len(npm_changed)
    fail_if(observed != expected, f"npm cross-version operation deltas: {observed}")


def main() -> int:
    summaries = {
        channel: read_json(f"v1.0.30-published-{channel}-chemistry-summary.json")
        for channel in ("python", "npm", "rust")
    }
    paths = {
        channel: RESULTS / f"v1.0.30-published-{channel}-chemistry-rows.jsonl.gz"
        for channel in summaries
    }
    oracle_path = RESULTS / "v1.0.30-rdkit-smarts-all-cells.jsonl.gz"
    for channel, summary in summaries.items():
        fail_if(raw_digest(paths[channel]) != summary["rows"]["sha256"], f"{channel}: raw digest")
        counts = summary["counts"]
        for key, expected in {"input": 10000, "completed": 10000, "parse_failure": 0,
                              "cip_exact": 9995, "morgan_exact": 9999,
                              "smarts_cells": 310000, "smarts_differences": 200}.items():
            # The Python reference uses its historical summary key names.
            python_key = {"input": "completed", "parse_failure": None}.get(key, key)
            actual = counts.get(python_key) if channel == "python" and python_key else counts.get(key)
            if key == "parse_failure" and channel == "python":
                actual = summary["row_accounting"]["parse_failure_count"]
            fail_if(actual != expected, f"{channel}: {key} = {actual}")
    fail_if(summaries["python"]["rdkit_version"] != "2026.03.6", "RDKit oracle version")
    fail_if(summaries["python"]["chematic_version"] != "1.0.30", "Python version")
    fail_if(summaries["npm"]["package"]["version"] != "1.0.30" or
            summaries["npm"]["package"]["tarball_sha256"] != NPM_TARBALL_SHA256,
            "npm artifact identity")
    fail_if(summaries["rust"]["crate"]["version"] != "1.0.30", "Rust crate version")
    lock = (ROOT / "tools/published_rust_gate/Cargo.lock").read_text(encoding="utf-8")
    fail_if(f'name = "chematic"\nversion = "1.0.30"\nsource = "registry+https://github.com/rust-lang/crates.io-index"\nchecksum = "{CRATE_SHA256}"' not in lock,
            "crates.io lockfile checksum")

    differences = Counter()
    with gzip.open(paths["python"], "rt", encoding="utf-8") as py_rows, \
         gzip.open(paths["npm"], "rt", encoding="utf-8") as npm_rows, \
         gzip.open(paths["rust"], "rt", encoding="utf-8") as rust_rows, \
         gzip.open(oracle_path, "rt", encoding="utf-8") as oracle_rows:
        manifest = json.loads(next(oracle_rows))
        fail_if(not manifest.get("_manifest") or manifest["rdkit_version"] != "2026.03.6", "SMARTS manifest")
        for index, lines in enumerate(zip_longest(py_rows, npm_rows, rust_rows, oracle_rows)):
            fail_if(any(line is None for line in lines), f"row count mismatch at {index}")
            py, npm, rust, oracle = (json.loads(line) for line in lines)
            fail_if(any(row["input_index"] != index or row["smiles"] != py["smiles"]
                        for row in (py, npm, rust, oracle)), f"row correspondence {index}")
            fail_if(py["status"] != npm["status"] or py["status"] != rust["status"], f"status {index}")
            fail_if(py["smiles_parse_write"]["chematic_canonical"] != npm["canonical"]
                    or npm["canonical"] != rust["canonical"], f"canonical {index}")
            py_cip = py["cip"]
            for channel, row in (("npm", npm), ("rust", rust)):
                unresolved = row["cip"]["unresolved"]
                if channel == "npm":
                    reason_names = {"lonePairCenter": "lone_pair_center", "oracleUnstable": "oracle_unstable"}
                    unresolved = {atom: reason_names.get(reason, reason) for atom, reason in unresolved.items()}
                expected_cip = {"atoms": py_cip["chematic_atoms"],
                                "bonds": py_cip["chematic_bonds"],
                                "unresolved": py_cip["chematic_unresolved"]}
                fail_if({**row["cip"], "unresolved": unresolved} != expected_cip,
                        f"{channel}: CIP {index}")
                if py["morgan"]["chematic_sha256"] is None:
                    fail_if("morgan_sha256" in row or
                            row.get("morgan_error") not in py["morgan"]["chematic_error"],
                            f"{channel}: unaccounted Morgan refusal {index}")
                else:
                    fail_if(row.get("morgan_sha256") != py["morgan"]["chematic_sha256"],
                            f"{channel}: Morgan {index}")
                fail_if(row["smarts_matches"] != npm["smarts_matches"] or
                        row["smarts_matches"] != rust["smarts_matches"], f"{channel}: SMARTS {index}")
            differences["smarts_cells"] += len(oracle["matches"])
            differences["smarts_differences"] += sum(a != b for a, b in zip(npm["smarts_matches"], oracle["matches"], strict=True))
            differences["rows"] += 1
    fail_if(dict(differences) != {"rows": 10000, "smarts_cells": 310000, "smarts_differences": 200},
            f"cross-binding SMARTS accounting: {dict(differences)}")

    old_path = RESULTS / "v1.0.29-published-python-63op-outputs.jsonl.gz"
    new_path = RESULTS / "v1.0.30-published-python-63op-outputs.jsonl.gz"
    diff = compare(old_path, new_path)
    archived = read_json("v1.0.29-to-v1.0.30-python-63op-diff.json")
    fail_if(diff != archived or diff["error_rows_by_operation"], "63-operation archived diff or errors")
    fail_if(diff["operation_count"] != 63 or diff["output_count"] != 210410 or
            {key: value["count"] for key, value in diff["changed_rows_by_operation"].items()} !=
            {"hba": 1359, "lipinski_bundle(mw,logp,hbd,hba)": 1359}, "63-operation accounting")
    fail_if(build_rust_operation_report() != read_json("v1.0.29-to-v1.0.30-published-rust-63op-diff.json"),
            "published Rust 63-operation cross-binding or version differential")

    reaction = read_json("v1.0.30-published-reaction-83-strata.json")
    matrix = json.loads((ROOT / "benchmarks/2026-10-02-v1.0.30-python-op-matrix.json").read_text(encoding="utf-8"))
    fail_if(matrix["chematic"]["version"] != "1.0.30" or
            matrix["chematic"]["artifact_sha256"] != WHEEL_SHA256 or
            reaction["artifact"]["sha256"] != WHEEL_SHA256,
            "published Python wheel identity")
    fail_if(reaction["accounting"]["input"] != 83 or
            reaction["accounting"]["outcomes"] !=
            {"semantic_match": 75, "joint_invalid_input": 3, "wrong_confident": 5} or
            reaction["accounting"]["by_stratum"]["legacy"] != {"semantic_match": 57},
            "reaction stratification accounting")
    npm_reaction = read_json("v1.0.30-published-npm-reaction-83-rows.json")
    all_reactions = read_json("v1.0.30-published-reaction-83-all-bindings.json")
    fail_if(hashlib.sha256((RESULTS / "v1.0.30-published-npm-reaction-83-rows.json").read_bytes()).hexdigest() != NPM_REACTION_ROWS_SHA256 or
            hashlib.sha256((RESULTS / "v1.0.30-published-reaction-83-all-bindings.json").read_bytes()).hexdigest() != ALL_REACTION_REPORT_SHA256,
            "published reaction raw record digest")
    fail_if(npm_reaction["package"]["tarball_sha256"] != NPM_TARBALL_SHA256 or
            all_reactions["npm_artifact"] != npm_reaction["package"] or
            len(npm_reaction["rows"]) != 83 or
            [row["id"] for row in npm_reaction["rows"]] != [row["id"] for row in all_reactions["rows"]] or
            all_reactions["npm_outcomes"] !=
            {"semantic_match": 74, "joint_invalid_input": 3, "wrong_confident": 5,
             "invalid_serialization": 1},
            "published npm reaction accounting")
    multiplicity_path = RESULTS / "v1.0.30-published-reaction-83-multiplicity.json"
    multiplicity = json.loads(multiplicity_path.read_bytes())
    fail_if(hashlib.sha256(multiplicity_path.read_bytes()).hexdigest() != REACTION_MULTIPLICITY_REPORT_SHA256 or
            multiplicity["schema"] != "published-reaction-distinct-vs-raw-multiplicity/v1" or
            multiplicity["version"] != "1.0.30" or multiplicity["rdkit_version"] != "2026.03.6" or
            multiplicity["source_sha256"] != ALL_REACTION_REPORT_SHA256 or
            multiplicity["fixtures"] != all_reactions["fixtures"] or
            [row["id"] for row in multiplicity["rows"]] != [row["id"] for row in all_reactions["rows"]] or
            multiplicity["accounting"] != {"input": 83, "by_binding": {
                "python": {"distinct_and_raw_counts_match": 62, "distinct_match_raw_count_lower": 13,
                           "distinct_product_sets_differ": 5, "oracle_invalid_or_unavailable": 3},
                "rust": {"candidate_refusal_or_invalid": 1, "distinct_and_raw_counts_match": 62,
                         "distinct_match_raw_count_lower": 13, "distinct_product_sets_differ": 4,
                         "oracle_invalid_or_unavailable": 3},
                "npm": {"candidate_refusal_or_invalid": 1, "distinct_and_raw_counts_match": 62,
                        "distinct_match_raw_count_lower": 12, "distinct_product_sets_differ": 5,
                        "oracle_invalid_or_unavailable": 3}}},
            "published reaction distinct/raw multiplicity identity")
    bad_json = next(row for row in npm_reaction["rows"] if row["id"] == "v2_stereo_e_alkene_preserve")
    fail_if(bad_json["status"] != "invalid_serialization", "published npm E/Z failure classification")
    try:
        json.loads(bad_json["payload"])
    except json.JSONDecodeError:
        pass
    else:
        raise ValueError("published npm E/Z payload unexpectedly parses as JSON")
    provenance_rows_path = RESULTS / "v1.0.30-published-rust-reaction-83-provenance-rows.json"
    provenance_summary = read_json("v1.0.30-published-rust-reaction-83-provenance-summary.json")
    provenance_report = read_json("v1.0.30-published-rust-reaction-83-provenance-rdkit.json")
    provenance_rows = json.loads(provenance_rows_path.read_text(encoding="utf-8"))
    old_rust_rows = read_json("v1.0.30-published-rust-reaction-83-rows.json")
    fail_if(provenance_summary["schema"] != "published-rust-reactions/v2" or
            provenance_summary["crate"] != "chematic 1.0.30 from crates.io" or
            provenance_summary["rows_sha256"] != hashlib.sha256(provenance_rows_path.read_bytes()).hexdigest() or
            provenance_report["rust_rows_sha256"] != provenance_summary["rows_sha256"] or
            provenance_report["rdkit_version"] != "2026.03.6" or
            provenance_report["rust_summary_sha256"] != hashlib.sha256(
                (RESULTS / "v1.0.30-published-rust-reaction-83-provenance-summary.json").read_bytes()
            ).hexdigest(),
            "published Rust provenance input identity")
    fail_if(len(provenance_rows) != 83 or len(provenance_report["rows"]) != 83 or
            [row["id"] for row in provenance_rows] != [row["id"] for row in provenance_report["rows"]] or
            dict(Counter(row["outcome"] for row in provenance_report["rows"])) !=
            provenance_report["accounting"]["outcomes"] or
            any({key: value for key, value in tagged.items() if key != "sets_with_atom_sources"} != plain
                for tagged, plain in zip(provenance_rows, old_rust_rows, strict=True)) or
            provenance_report["accounting"] != {"input": 83, "outcomes": {
                "graph_and_provenance_match": 74, "graph_mismatch": 4,
                "joint_invalid_input": 3, "provenance_mismatch": 1,
                "typed_or_diagnosed_refusal": 1}},
            "published Rust provenance row accounting or untagged graph identity")
    map_rows_path = RESULTS / "v1.0.30-published-rust-reaction-83-template-map-rows.json"
    map_summary_path = RESULTS / "v1.0.30-published-rust-reaction-83-template-map-summary.json"
    map_report_path = RESULTS / "v1.0.30-published-rust-reaction-83-template-map-rdkit.json"
    map_rows = json.loads(map_rows_path.read_text(encoding="utf-8"))
    map_summary = json.loads(map_summary_path.read_text(encoding="utf-8"))
    map_report = json.loads(map_report_path.read_text(encoding="utf-8"))
    fail_if(hashlib.sha256(map_report_path.read_bytes()).hexdigest() != REACTION_TEMPLATE_MAP_REPORT_SHA256 or
            map_summary["schema"] != "published-rust-reaction-template-maps/v1" or
            map_summary["crate"] != "chematic 1.0.30 from crates.io" or
            map_summary["input_count"] != 83 or
            map_summary["rows_sha256"] != hashlib.sha256(map_rows_path.read_bytes()).hexdigest() or
            map_summary["base_cases_sha256"] != provenance_summary["base_cases_sha256"] or
            map_summary["strata_sha256"] != provenance_summary["strata_sha256"] or
            map_report["rdkit_version"] != "2026.03.6" or
            map_report["crate_version"] != "1.0.30" or
            map_report["rust_rows_sha256"] != map_summary["rows_sha256"] or
            map_report["rust_summary_sha256"] != hashlib.sha256(map_summary_path.read_bytes()).hexdigest() or
            map_report["baseline_rows_sha256"] != provenance_summary["rows_sha256"] or
            map_report["fixtures"] != {
                "base_sha256": map_summary["base_cases_sha256"],
                "strata_sha256": map_summary["strata_sha256"]},
            "published Rust template-map evidence identity")
    fail_if(len(map_rows) != 83 or len(map_report["rows"]) != 83 or
            [row["id"] for row in map_rows] != [row["id"] for row in map_report["rows"]] or
            [row["id"] for row in map_rows] != [row["id"] for row in provenance_rows] or
            any(row["status"] != baseline["status"]
                for row, baseline in zip(map_rows, provenance_rows, strict=True)) or
            dict(Counter(row["outcome"] for row in map_report["rows"])) !=
            map_report["accounting"]["outcomes"] or
            map_report["accounting"] != {"input": 83, "outcomes": {
                "graph_origin_map_match": 73, "graph_mismatch": 4,
                "map_label_mismatch": 1, "provenance_mismatch": 1,
                "typed_or_diagnosed_refusal": 1, "joint_invalid_input": 3}} or
            [row["id"] for row in map_report["rows"] if row["outcome"] == "map_label_mismatch"] !=
            ["alkene_hydrogenation_styrene"] or
            any(len(product["atom_sources"]) != len(product["template_map_numbers"])
                for row in map_rows if row["status"] == "products"
                for products in row["sets"] for product in products),
            "published Rust template-map outcome and row accounting")
    supplement_path = ROOT / "validation/reaction_product_new_maps_v3.json"
    supplement = json.loads(supplement_path.read_bytes())
    new_rows_path = RESULTS / "v1.0.30-published-rust-reaction-86-new-map-rows.json"
    new_summary_path = RESULTS / "v1.0.30-published-rust-reaction-86-new-map-summary.json"
    new_report_path = RESULTS / "v1.0.30-published-rust-reaction-new-map-parity.json"
    new_rows = json.loads(new_rows_path.read_bytes())
    new_summary = json.loads(new_summary_path.read_bytes())
    new_report = json.loads(new_report_path.read_bytes())
    fail_if(hashlib.sha256(new_report_path.read_bytes()).hexdigest() != REACTION_NEW_MAP_REPORT_SHA256 or
            supplement["schema_version"] != 3 or len(supplement["cases"]) != 3 or
            supplement["base_fixture_sha256"] != map_summary["base_cases_sha256"] or
            supplement["strata_fixture_sha256"] != map_summary["strata_sha256"] or
            new_summary["schema"] != "published-rust-reaction-template-maps/v1" or
            new_summary["crate"] != "chematic 1.0.30 from crates.io" or
            new_summary["input_count"] != 86 or
            new_summary["supplement_sha256"] != hashlib.sha256(supplement_path.read_bytes()).hexdigest() or
            new_summary["rows_sha256"] != hashlib.sha256(new_rows_path.read_bytes()).hexdigest() or
            new_rows[:83] != map_rows or
            [row["id"] for row in new_rows[83:]] != [case["id"] for case in supplement["cases"]] or
            new_report["schema"] != "published-rust-reaction-new-product-maps/v1" or
            new_report["crate_version"] != "1.0.30" or
            new_report["rdkit_version"] != "2026.03.6" or
            new_report["rust_rows_sha256"] != new_summary["rows_sha256"] or
            new_report["rust_summary_sha256"] != hashlib.sha256(new_summary_path.read_bytes()).hexdigest() or
            new_report["baseline_rows_sha256"] != hashlib.sha256(map_rows_path.read_bytes()).hexdigest() or
            new_report["fixtures"] != {
                "base_sha256": map_summary["base_cases_sha256"],
                "strata_sha256": map_summary["strata_sha256"],
                "supplement_sha256": new_summary["supplement_sha256"]} or
            new_report["accounting"] != {"input": 3, "outcomes": {"graph_origin_map_match": 3}} or
            [row["id"] for row in new_report["rows"]] != [case["id"] for case in supplement["cases"]] or
            any(not any(source is None and label is not None
                        for products in row["sets"] for product in products
                        for source, label in zip(product["atom_sources"], product["template_map_numbers"], strict=True))
                for row, case in zip(new_rows[83:], supplement["cases"], strict=True)
                if "mapped" in case["strata"]),
            "published Rust new-product-map evidence identity")
    paired_path = ROOT / "benchmarks/2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json"
    fail_if(hashlib.sha256(paired_path.read_bytes()).hexdigest() != PAIRED_63OP_SHA256,
            "published Python paired 63-operation record digest")
    paired = json.loads(paired_path.read_text(encoding="utf-8"))
    fail_if(paired["chematic"]["artifact_sha256"] != WHEEL_SHA256 or
            paired["rdkit"]["version"] != "2026.03.6" or
            paired["protocol"]["repeats"] != 20 or
            paired["summary"]["operations"] != 63 or
            paired["summary"]["equivalent_output_speed_wins_with_ci"] != 20 or
            any(len(op["paired_blocks"]) != 20 for op in paired["operations"]),
            "published Python paired 63-operation matrix accounting")
    check_npm_version_diff()
    check_chemistry_residuals()
    check_published_wasm_paired()
    browser_raw = ROOT / "benchmarks/2026-10-03-v1030-rdkitjs-published-chromium-paired20.json"
    browser_summary = ROOT / "benchmarks/2026-10-03-v1030-rdkitjs-published-chromium-paired20-summary.json"
    fail_if(hashlib.sha256(browser_raw.read_bytes()).hexdigest() != BROWSER_PAIRED_RAW_SHA256 or
            hashlib.sha256(browser_summary.read_bytes()).hexdigest() != BROWSER_PAIRED_SUMMARY_SHA256 or
            check_published_browser_paired(browser_raw) != json.loads(browser_summary.read_bytes()),
            "published Chromium paired timing/output/memory record identity")
    prepared_raw = ROOT / "benchmarks/2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.json"
    prepared_summary = ROOT / "benchmarks/2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20-summary.json"
    fail_if(hashlib.sha256(prepared_raw.read_bytes()).hexdigest() != PREPARED_SPLIT_RAW_SHA256 or
            hashlib.sha256(prepared_summary.read_bytes()).hexdigest() != PREPARED_SPLIT_SUMMARY_SHA256 or
            check_published_browser_paired(prepared_raw) != json.loads(prepared_summary.read_bytes()),
            "published Chromium prepared first-use/reused record identity")
    for (rows, operation), digest in PUBLISHED_BROWSER_MORGAN_ROWS.items():
        prefix = "row" if rows == 250 else "10k"
        path = ROOT / "benchmarks" / f"2026-10-03-v1030-published-chromium-morgan-{prefix}-{operation}.json"
        fail_if(hashlib.sha256(path.read_bytes()).hexdigest() != digest or
                check_browser_morgan_rows(path, operation, rows)["exact"] != rows - (rows == 10_000),
                f"published Chromium {rows}-row {operation} Morgan bit parity")
    print("v1.0.30 packet integrity OK: 3 published artifacts x 10k/310k; "
          "Python/Rust 63 operations x 210,410 rows each per version (only HBA/bundle changed); "
          "npm 59/63 output adapters (52 exact, ETKDG rounded, 3 representation lanes differ), "
          "v1.0.29 to v1.0.30 only HBA/bundle changed on the same 1,359 rows; "
          "200 SMARTS residuals and 5 CIP abstentions classified on published rows; "
          "20-block Python speed matrix, 8 output-gated Node/WASM lanes and "
          "3 output-gated Chromium Morgan lanes with prepared first-use/reused split and "
          "direct browser row-level Morgan checks (250/250 plus 9999/10000 and one typed refusal); "
          "npm reaction 74 match, 5 confident "
          "differences, 1 invalid JSON; Rust reaction origins 74/83 matched; "
          "73/83 graph-origin-map matches, 1 map-only residual; "
          "3/3 exposed new-product-map supplements match; distinct/raw reaction counts classified. "
          "P0.1 output audit complete with explicit failures/gaps; "
          "strict parity and P0.2/P1 acceptance remain OPEN.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
