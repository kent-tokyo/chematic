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
    print("v1.0.30 packet integrity OK: 3 published artifacts x 10k/310k; "
          "Python 63 operations x 210,410 rows (only HBA/bundle changed); "
          "reaction extension 75 match, 3 jointly invalid, 5 confident differences. "
          "P0/P1 acceptance remains OPEN.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
