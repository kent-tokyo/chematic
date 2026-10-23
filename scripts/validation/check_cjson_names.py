#!/usr/bin/env python3
"""Authored name-preservation development packet for CJSON ⇄ MOL2.

python3 -m scripts.validation.check_cjson_names --cli CLI --obabel OBABEL
--output NEW.json. Standard root names, legacy aliases, conflicting names
and unrepresentable MOL2 titles have separate outcomes; refusal is not coverage.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
from scripts.validation.check_cjson_mol2_formal_charges import mol2_snapshot


def packet():
    specs = [
        ("root-ascii", "helium standard", None, False),
        ("root-unicode", "ヘリウム β", None, False),
        ("legacy", None, "helium legacy", False),
        ("both-equal", "same title", "same title", False),
        ("absent", None, None, False),
        ("conflict", "root title", "legacy title", True),
        ("empty", "", None, True),
        ("line-break", "head\n@<TRIPOS>ATOM", None, True),
        ("boundary-space", " title ", None, True),
    ]
    cases = []
    for name, root, legacy, refuse in specs:
        source = {"chemicalJson": 1, "atoms": {"elements": {"number": [2]},
                  "formalCharges": [0], "coords": {"3d": [0, 0, 0]}},
                  "bonds": {"connections": {"index": []}, "order": []}}
        if root is not None:
            source["name"] = root
        if legacy is not None:
            source["properties"] = {"name": legacy}
        cases.append({"id": name, "input_format": "cjson", "output_format": "mol2",
                      "source": json.dumps(source, ensure_ascii=False, sort_keys=True,
                                           separators=(",", ":")),
                      "expected": {"name": root if root is not None else legacy,
                                   "refuse": refuse}})
    for name in ["helium standard", "ヘリウム β"]:
        cases.append({"id": "reverse-" + ("ascii" if name.isascii() else "unicode"),
                      "input_format": "mol2", "output_format": "cjson",
                      "source": "@<TRIPOS>MOLECULE\n" + name
                      + "\n1 0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 He1 0 0 0 He\n@<TRIPOS>BOND\n",
                      "expected": {"name": name, "refuse": False}})
    return cases


def inspect(text, target, expected):
    if target == "mol2":
        snapshot = mol2_snapshot(text)
        name = text.splitlines()[1]
    else:
        value = strict_json(text)
        # A legacy location alone is not semantic data loss in this profile.
        name = value.get("name", value.get("properties", {}).get("name"))
        if "name" in value and "name" in value.get("properties", {}):
            require(value["name"] == value["properties"]["name"], "conflicting output names")
        snapshot = {"elements": value["atoms"]["elements"]["number"],
                    "formal_charges": value["atoms"].get("formalCharges", [0]),
                    "coordinates": value["atoms"]["coords"]["3d"],
                    "bonds": value.get("bonds", {}).get("connections", {}).get("index", [])}
    return {"name": expected["name"] is None or name == expected["name"],
            "explicit_graph_and_coordinates": snapshot["elements"] == [2]
            and snapshot["formal_charges"] == [0] and snapshot["coordinates"] == [0, 0, 0]
            and snapshot["bonds"] == []}, name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--obabel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "output already exists")
    version = bounded_run([str(args.obabel), "-V"])
    require(version["exit_code"] == 0 and b"Open Babel 3.2.1" in version["stdout"],
            "expected pinned Open Babel 3.2.1")
    cases, rows = packet(), []
    for case in cases:
        data = case["source"].encode()
        case["sha256"] = hashlib.sha256(data).hexdigest()
        for engine, command in [
            ("chematic", [str(args.cli), "convert-with-report", "--input-format",
                          case["input_format"], "--output-format", case["output_format"],
                          "--policy", "allow_loss"]),
            ("openbabel", [str(args.obabel), "-i" + case["input_format"],
                           "-o" + case["output_format"]]),
        ]:
            run = bounded_run(command, data, timeout=10, max_output=1048576)
            row = {"case": case["id"], "engine": engine, "command": command,
                   **{key: value.decode("utf8", errors="replace") if isinstance(value, bytes)
                      else value for key, value in run.items()}}
            row["terminal"] = "process_failure"
            if run["resource_status"] != "completed":
                row["terminal"] = run["resource_status"]
            else:
                try:
                    report = strict_json(row["stdout"]) if engine == "chematic" else None
                    if report and report["status"] == "refused" and report["output"] is None:
                        row["terminal"] = "refused"
                    elif run["exit_code"] == 0:
                        text = report["output"] if report else row["stdout"]
                        row["fields"], row["observed_name"] = inspect(
                            text, case["output_format"], case["expected"])
                        row["terminal"] = ("output_on_ambiguous_or_unrepresentable_source"
                                           if case["expected"]["refuse"] else
                                           "packet_fields_preserved" if all(row["fields"].values())
                                           else "packet_field_loss")
                except (ValueError, KeyError, IndexError, StopIteration, TypeError) as error:
                    row["terminal"], row["error"] = "unresolved_output", str(error)
            rows.append(row)
    result = {"schema": "chematic.cjson-name-development.v1",
              "scope": "Original development cases for names and one explicit helium graph; not independent holdout, chemistry perception, full-format preservation or toolkit superiority. Refusal is separate from successful conversion.",
              "provenance": "Original synthetic inputs, 2026-10-21; recipe packet() in this script; no third-party or sealed corpus input.",
              "license": "MIT OR Apache-2.0", "redistribution": "allowed: original inputs",
              "minimization": "One helium atom isolates title handling from typing/perception; source coordinates are supplied in angstrom, exact integers with no tolerance.",
              "oracle": {"name": "independent bounded packet inspection", "version": 1,
                         "script": identity(__file__), "expected": "Names compared verbatim; no title assigned when both declarations conflict; no MOL2 candidate required for names outside the declared title profile."},
              "engines": {"chematic": identity(args.cli), "openbabel": identity(args.obabel)},
              "openbabel_version": version["stdout"].decode().strip(), "cases": cases, "rows": rows,
              "counts": {e: dict(Counter(r["terminal"] for r in rows if r["engine"] == e))
                         for e in ["chematic", "openbabel"]}}
    with args.output.open("x") as output:
        output.write(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
