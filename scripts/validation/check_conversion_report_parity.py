#!/usr/bin/env python3
"""Check full checked-conversion reports across explicit local CLI/Python/WASM.

Development-fixture parity only; this is not an independent chemistry oracle
or an Open Babel comparison. Select a source-built wheel via PYTHONPATH.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]

PYTHON_CALL = 'import chematic,sys; print(chematic.convert_with_report_json(sys.stdin.buffer.read().decode("utf-8"), *sys.argv[1:]))'
NODE_CALL = """
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm = await import(pathToFileURL(process.argv[1]).href);
console.log(wasm.convert_with_report_json(readFileSync(0, 'utf8'), ...process.argv.slice(2)));
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--node", default="node")
    parser.add_argument("--wasm-module", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=ROOT / "validation/mol2_preservation/amide.mol2")
    args = parser.parse_args()
    source = args.source.read_text()
    cases = []
    for variant, text in [("lf", source), ("crlf", source.replace("\n", "\r\n"))]:
        for target, policies in [("mol2", ["strict", "allow_loss", "preserve_source"]), ("smiles", ["strict", "allow_loss"]), ("cjson", ["strict", "allow_loss"])]:
            for policy in policies:
                cases.append((f"{variant}-{target}-{policy}", text, "mol2", target, policy, False))
    cases.append(("zero-bond-mismatch", source.replace("12 7 3 am", "12 7 3 nc"), "mol2", "smiles", "allow_loss", False))
    for token in ["nc", "du", "un"]:
        cases.append((f"cjson-bond-{token}", source.replace("12 7 3 am", f"12 7 3 {token}"), "mol2", "cjson", "allow_loss", False))
    for value in ["-1.4175225513503917", "-0.0", "5e-324", "1.7976931348623157e308"]:
        cases.append((f"cjson-coordinate-{value}", source.replace("0.123456789", value), "mol2", "cjson", "allow_loss", False))
    cases.append(("cjson-N.4", source.replace("0.123456789 0 0 C.3", "0.123456789 2.5 -3.75 N.4"), "mol2", "cjson", "allow_loss", False))
    for name, text, source_format, target, policy in [
        ("invalid", "invalid", "mol2", "mol2", "strict"),
        ("bad-charge", source.replace("ACM -0.70", "ACM NaN"), "mol2", "mol2", "strict"),
        ("unsupported-source", source, "sdf", "smiles", "strict"),
        ("unsupported-target", source, "mol2", "cml", "allow_loss"),
        ("unsupported-policy", source, "mol2", "smiles", "typo"),
        ("cross-source-retention", source, "mol2", "smiles", "preserve_source"),
        ("cjson-source-retention", source, "mol2", "cjson", "preserve_source"),
    ]:
        cases.append((name, text, source_format, target, policy, True))
    for value in ["-0.0", "-1.4175225513503917", "5e-324", "1.7976931348623157e308"]:
        cases.append((f"cjson-partial-charge-{value}", source.replace("ACM 0.70", f"ACM {value}"), "mol2", "cjson", "allow_loss", False))
    for name, text in [("missing", source.replace("ACM 0.70", "ACM")), ("undeclared", source.replace("USER_CHARGES", "NO_CHARGES"))]:
        cases.append((f"cjson-partial-charge-{name}", text, "mol2", "cjson", "allow_loss", False))
    cjson = (ROOT / "validation/cjson_preservation/charged.cjson").read_text()
    for variant, text in [("lf", cjson), ("crlf", cjson.replace("\n", "\r\n"))]:
        for target, policies in [("cjson", ["strict", "allow_loss", "preserve_source"]), ("smiles", ["strict", "allow_loss"])]:
            for policy in policies:
                cases.append((f"cjson-input-{variant}-{target}-{policy}", text, "cjson", target, policy, False))
    for name, text in [
        ("minimal", '{"atoms":{"elements":{"number":[6]}}}'),
        ("branch", '{"atoms":{"elements":{"number":[6,8,7,6]},"isotopes":[13,18,15,12],"formalCharges":[0,0,0,0]},"bonds":{"connections":{"index":[0,1,0,2,1,3]},"order":[1,1,1]}}'),
        ("aromatic", '{"atoms":{"elements":{"number":[6,6]}},"bonds":{"connections":{"index":[0,1]},"order":[4]}}'),
    ]:
        cases.append((f"cjson-input-{name}", text, "cjson", "smiles", "allow_loss", False))
    for name, text, target, policy in [
        ("charge-length", cjson.replace("[0, 1, -1, 0]", "[]"), "cjson", "strict"),
        ("isotope-type", cjson.replace("[13, 15, 0, 2]", "null"), "cjson", "strict"),
        ("unknown-order", cjson.replace('"order": [1, 1]', '"order": [2.2, 1]'), "smiles", "allow_loss"),
        ("duplicate-edge", cjson.replace("[0, 1, 1, 3]", "[0, 1, 1, 0]"), "smiles", "allow_loss"),
        ("self-edge", cjson.replace("[0, 1, 1, 3]", "[0, 0, 1, 3]"), "cjson", "strict"),
        ("dangling-edge", cjson.replace("[0, 1, 1, 3]", "[0, 4, 1, 3]"), "cjson", "strict"),
        ("duplicate-opaque-key", cjson.replace('"stereo": "opaque"', '"stereo": "opaque", "stereo": "other"'), "cjson", "strict"),
        ("unsupported-version", cjson.replace('"chemicalJson": 1', '"chemicalJson": 0'), "cjson", "strict"),
        ("unsupported-route", cjson, "cml", "allow_loss"),
        ("cross-preserve", cjson, "smiles", "preserve_source"),
    ]:
        cases.append((f"cjson-input-{name}", text, "cjson", target, policy, True))
    partial = (ROOT / "validation/cjson_preservation/partial-charges.cjson").read_text()
    for variant, text in [("lf", partial), ("crlf", partial.replace("\n", "\r\n"))]:
        for policy in ["strict", "allow_loss"]:
            cases.append((f"cjson-mol2-{variant}-{policy}", text, "cjson", "mol2", policy, False))
    for n in [-0.0, -1.4175225513503917, 5e-324, 1.7976931348623157e308]:
        value = json.loads(partial)
        value["partialCharges"]["USER_CHARGES"][0] = n
        value["atoms"]["coords"]["3d"][0] = n
        cases.append((f"cjson-mol2-float-{n}", json.dumps(value), "cjson", "mol2", "allow_loss", False))
    for mode in ["multiple", "missing-coordinates", "name-newline", "method-whitespace", "absent-dataset", "isotope-loss", "signed-formal-charge"]:
        value = json.loads(partial)
        if mode == "multiple": value["partialCharges"]["SECOND"] = [0] * 4
        elif mode == "missing-coordinates": del value["atoms"]["coords"]
        elif mode == "name-newline": value["properties"]["name"] = "a\n@<TRIPOS>ATOM"
        elif mode == "method-whitespace": value["partialCharges"] = {"two methods": [0] * 4}
        elif mode == "absent-dataset": del value["partialCharges"]
        elif mode == "isotope-loss": value["atoms"]["isotopes"][0] = 13
        else: value["atoms"]["formalCharges"][0] = -1
        cases.append((f"cjson-mol2-{mode}", json.dumps(value), "cjson", "mol2", "allow_loss", False))
    for i, bad in enumerate([None, [], {"method": []}, {"method": [0, 0, 0]}, {"method": [0, 0, 0, "NaN"]}]):
        value = json.loads(partial)
        value["partialCharges"] = bad
        cases.append((f"cjson-mol2-malformed-dataset-{i}", json.dumps(value), "cjson", "mol2", "allow_loss", True))
    for root, legacy in [("ヘリウム β", None), ("root", "legacy")]:
        value = {"atoms": {"elements": {"number": [2]}, "coords": {"3d": [0, 0, 0]}}, "name": root}
        if legacy is not None:
            value["properties"] = {"name": legacy}
        for target, policy in [("cjson", "strict"), ("mol2", "allow_loss"), ("smiles", "allow_loss")]:
            cases.append((f"cjson-name-{root}-{target}", json.dumps(value), "cjson", target, policy, False))
    results = []
    for name, text, source_format, target, policy, expect_error in cases:
        commands = {
            "cli": [str(args.cli.resolve()), "convert-with-report", "--input-format", source_format, "--output-format", target, "--policy", policy],
            "python": [args.python, "-c", PYTHON_CALL, source_format, target, policy],
            "wasm": [args.node, "--input-type=module", "-e", NODE_CALL, str(args.wasm_module.resolve()), source_format, target, policy],
        }
        reports = {}
        for binding, command in commands.items():
            run = subprocess.run(command, input=text, text=True, capture_output=True, timeout=30)
            if expect_error:
                if run.returncode == 0 or run.stdout.strip():
                    raise RuntimeError(f"{name}/{binding}: expected error without report")
                continue
            try:
                report = json.loads(run.stdout)
            except json.JSONDecodeError as error:
                raise RuntimeError(f"{name}/{binding}: invalid report: {run.stderr}") from error
            expected_exit = 2 if binding == "cli" and report["status"] == "refused" else 0
            if run.returncode != expected_exit:
                raise RuntimeError(f"{name}/{binding}: exit {run.returncode}: {run.stderr}")
            reports[binding] = report
        if not expect_error and not (reports["cli"] == reports["python"] == reports["wasm"]):
            raise RuntimeError(f"{name}: reports differ across bindings")
        results.append({"case": name, "result": "shared_error" if expect_error else "identical_report"})
    print(json.dumps({"scope": "source-built development-fixture binding parity", "passed": len(results), "results": results}, indent=2))


if __name__ == "__main__":
    main()
