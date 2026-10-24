#!/usr/bin/env python3
"""Measure numeric XYZ reads separately from refusals and OB XYZ rewrite output."""
from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import importlib
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require

DECIMAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")
SOURCES = ["crates/chematic-3d/src/xyz.rs", "crates/chematic-mol/src/xyz.rs",
           "crates/chematic-mol/src/xyz/extxyz.rs", "crates/chematic-py/src/formats/xyz.rs"]


def observed(value):
    return dict(status="accepted", finite=math.isfinite(value),
                bits=struct.pack(">d", value).hex())


def output_scalar(text):
    """Only observe a single ordinary element/XYZ row, including nonfinite text."""
    lines = text.splitlines()
    require(len(lines) == 3 and lines[0].strip() == "1", "not exactly one XYZ record")
    row = lines[2].split()
    require(len(row) == 4 and row[0] == "C", "unexpected atom/row layout")
    values = list(map(float, row[1:]))
    require(values[1:] == [0.0, 0.0], "other axes changed")
    return observed(values[0])


def summarize(records, lane):
    valid = [r[lane] for r in records if r["expected"] == "accept"]
    invalid = [r[lane] for r in records if r["expected"] == "refuse"]
    return dict(valid_count=len(valid),
                exact_finite=sum(r.get("exact", False) for r in valid),
                refusal_count=len(invalid),
                explicit_refusals=sum(r["status"] == "refused" for r in invalid),
                unsafe_acceptances=sum(r["status"] == "accepted" for r in invalid),
                no_record_or_unobservable=sum(r["status"] not in {"accepted", "refused"} for r in invalid))


def source_receipts(reference, parent):
    if reference is None:
        return dict(reference="working tree", files=[identity(ROOT / p) for p in SOURCES])
    commit = subprocess.check_output(["git", "rev-parse", reference], cwd=ROOT, text=True).strip()
    files = []
    for path in SOURCES:
        data = subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)
        require(data == subprocess.check_output(["git", "show", f"{parent}:{path}"], cwd=ROOT),
                f"baseline source differs from comparison parent: {path}")
        files.append(dict(path=path, sha256=hashlib.sha256(data).hexdigest(), size_bytes=len(data)))
    return dict(reference=commit, matches_parent_sources=True, files=files)


def validate_oracle(case, accept):
    token = case["token"]
    if accept:
        require(DECIMAL.fullmatch(token), "expected decimal grammar")
        value = float(token)
        require(math.isfinite(value) and struct.pack(">d", value).hex() == case["bits"], "fixture bit oracle")
        require(value != 0.0 or Decimal(token).is_zero(), "nonzero decimal lost to zero in valid oracle")
    elif case["reason"] == "underflow":
        require(DECIMAL.fullmatch(token) and not Decimal(token).is_zero() and float(token) == 0.0,
                "fixture underflow oracle")
    elif case["reason"] == "nonfinite":
        require(not math.isfinite(float(token)), "fixture nonfinite oracle")
    else:
        require(not DECIMAL.fullmatch(token), "fixture grammar refusal oracle")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--baseline-ref")
    p.add_argument("--parent-ref", default="HEAD")
    p.add_argument("--obabel", type=Path, default=Path("/opt/homebrew/bin/obabel"))
    args = p.parse_args()
    import chematic

    parent = subprocess.check_output(["git", "rev-parse", args.parent_ref], cwd=ROOT, text=True).strip()
    fixture = ROOT / "validation/xyz_numeric_read/cases.json"
    packet = json.loads(fixture.read_text(encoding="utf-8"))
    records = []
    for group in ["valid", "invalid"]:
        for case in packet[group]:
            accept = group == "valid"
            validate_oracle(case, accept)
            text = f"1\nnumeric read\nC {case['token']} 0 0\n"
            record = dict(case, expected="accept" if accept else "refuse", input=text)
            for lane, reader in [("python_xyz", lambda: chematic.from_xyz(text)[1][0][0]),
                                 ("python_extxyz", lambda: chematic.from_extxyz(text)["coords"][0][0])]:
                try:
                    record[lane] = observed(reader())
                except ValueError as error:
                    record[lane] = dict(status="refused", error=str(error))
            command = [str(args.obabel), "-ixyz", "-oxyz", "-ab"]
            proc = bounded_run(command, text.encode("utf-8"))
            proc = {k: v.decode("utf-8", errors="replace") if isinstance(v, bytes) else v for k, v in proc.items()}
            observation = dict(status="no_record_or_unobservable")
            if proc["resource_status"] == "completed" and proc["exit_code"] == 0:
                try:
                    observation = output_scalar(proc["stdout"])
                except (ValueError, IndexError) as error:
                    observation["error"] = str(error)
            record["openbabel_xyz_rewrite"] = dict(observation, command=command, process=proc)
            for lane in ["python_xyz", "python_extxyz", "openbabel_xyz_rewrite"]:
                r = record[lane]
                r["exact"] = bool(accept and r["status"] == "accepted" and r["finite"] and r["bits"] == case["bits"])
            records.append(record)
    summary = {lane: summarize(records, lane) for lane in ["python_xyz", "python_extxyz", "openbabel_xyz_rewrite"]}
    extension = importlib.import_module("chematic.chematic")
    result = dict(schema="chematic.xyz-numeric-read-observation.v1", baseline=bool(args.baseline_ref),
                  scope="Authored binary64 numeric read/refusal packet. Refusals are separate from preserved values. "
                        "Open Babel is observed through plain XYZ rewrite output; parser and writer are not isolated. "
                        "No speed, arbitrary precision, general format or toolkit superiority claim.",
                  parent_commit=parent, fixture=identity(fixture), python_version=sys.version,
                  python_runtime=identity(extension.__file__), binding_sources=source_receipts(args.baseline_ref, parent),
                  harness=identity(Path(__file__)), process_runner=identity(ROOT / "scripts/checks/check_interchange_preservation.py"),
                  openbabel=identity(args.obabel), openbabel_version=subprocess.check_output([str(args.obabel), "-V"], text=True).strip(),
                  summary=summary, records=records)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary))
    if not args.baseline_ref:
        for lane in ["python_xyz", "python_extxyz"]:
            require(summary[lane]["exact_finite"] == len(packet["valid"]), "candidate finite bit mismatch")
            require(summary[lane]["explicit_refusals"] == len(packet["invalid"]), "candidate accepted unsafe scalar")


if __name__ == "__main__":
    main()
