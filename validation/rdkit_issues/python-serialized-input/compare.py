#!/usr/bin/env python3
"""Compare actual baseline/candidate extensions in isolated Python processes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def probe():
    import chematic
    import chematic.chematic as native

    def snapshot(value):
        if isinstance(value, chematic.Mol):
            return {"smiles": value.smiles, "formula": value.formula,
                    "atoms": value.atom_count, "bonds": value.bond_table}
        if isinstance(value, (tuple, list)):
            return [snapshot(item) for item in value]
        if isinstance(value, dict):
            return {key: snapshot(item) for key, item in value.items()}
        return value

    def call(reader, value):
        try:
            return {"result": snapshot(reader(value))}
        except Exception as error:
            return {"error": type(error).__name__, "message": str(error),
                    "attributes": {key: getattr(error, key, None)
                                   for key in ("category", "code", "format")}}

    rows = []
    for case in json.loads((HERE / "cases.json").read_text())["cases"]:
        reader = getattr(chematic, case["function"])
        rows.append({"function": case["function"],
                     "str": call(reader, case["text"]),
                     "bytes": call(reader, case["text"].encode()),
                     "malformed_str": [call(reader, text) for text in ("", "not a format", "C\0(")],
                     "malformed_bytes": [call(reader, text.encode()) for text in ("", "not a format", "C\0(")],
                     "surrogate": call(reader, "\ud800")})
    return {"module_sha256": digest(native.__file__), "rows": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline")
    parser.add_argument("--candidate")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--probe", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.probe:
        print(json.dumps(probe(), ensure_ascii=True))
        return
    if not all((args.baseline, args.candidate, args.output)):
        parser.error("--baseline, --candidate and --output are required")
    runs = {}
    for name in ("baseline", "candidate"):
        env = {**os.environ, "PYTHONPATH": getattr(args, name)}
        runs[name] = json.loads(subprocess.check_output(
            [sys.executable, str(Path(__file__).resolve()), "--probe"], env=env, text=True))
    baseline, candidate = runs["baseline"]["rows"], runs["candidate"]["rows"]
    assert len(baseline) == len(candidate) == 13
    for before, after in zip(baseline, candidate):
        assert before["function"] == after["function"]
        assert before["bytes"]["error"] == "TypeError", before
        assert "result" in before["str"], before
        assert before["str"] == after["str"] == after["bytes"], after
        assert before["malformed_str"] == after["malformed_str"] == after["malformed_bytes"], after
        assert before["surrogate"] == after["surrogate"], after
    report = {"schema": "chematic.serialized-input-comparison.v1",
              "issue": "https://github.com/rdkit/rdkit/issues/9601",
              "scope": "Thirteen Python text inputs; local extensions, not released artifacts or a chemistry/speed benchmark.",
              "python": sys.version, "fixture_sha256": digest(HERE / "cases.json"),
              "checks": {"baseline_bytes_rejections": 13, "preserved_str_outputs": 13,
                         "candidate_bytes_matches": 13, "preserved_malformed_inputs": 39,
                         "preserved_surrogate_errors": 13}, "runs": runs}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n")
    print(json.dumps(report["checks"]))


if __name__ == "__main__":
    main()
