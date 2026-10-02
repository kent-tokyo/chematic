#!/usr/bin/env python3
"""Verify the exposed v1.0.30 published-artifact packet without RDKit installed.

Integrity success is not P0/P1 acceptance: chemistry has documented
unadjudicated rows, the 83-reaction extension has confident differences, and
the binding-wide 63-operation/speed matrices remain incomplete.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from itertools import zip_longest
from pathlib import Path

from check_published_python_version_outputs import compare

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation/results"
WHEEL_SHA256 = "f6b1b22dd898cf3f100001156d9504b3753cb6b392143cbc0e6f3be78ec723d5"
NPM_TARBALL_SHA256 = "fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201"
CRATE_SHA256 = "b1a897689ab5b4983e56325dbd361e8ff0625dd752984af7fe75a10286f84a82"
NPM_REACTION_ROWS_SHA256 = "ce3d801ed55f610412553861382040356ac17a8004b0bec23223bc14ddbdee9b"
ALL_REACTION_REPORT_SHA256 = "a43cf9eec38b04664dc8ee269aa69f51a9ae3e310fa5a781177d07f389e85785"
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
    print("v1.0.30 packet integrity OK: 3 published artifacts x 10k/310k; "
          "Python 63 operations x 210,410 rows (only HBA/bundle changed); "
          "20-block Python speed matrix; npm reaction 74 match, 5 confident "
          "differences, 1 invalid JSON; Rust reaction origins 74/83 matched. "
          "P0/P1 acceptance remains OPEN.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
