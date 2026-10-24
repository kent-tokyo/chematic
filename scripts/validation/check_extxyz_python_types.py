#!/usr/bin/env python3
"""Observe bounded typed ExtXYZ Python round trips; OB uses its plain XYZ route."""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import shlex
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require


def bits(value):
    value = float(value)
    require(math.isfinite(value), "nonfinite real")
    return struct.pack(">d", value).hex()


def scalar(token, kind):
    if kind == "R":
        return bits(token)
    if kind == "I":
        value = int(token)
        require(-(1 << 63) <= value < 1 << 63, "integer outside i64")
        return value
    if kind == "L":
        require(token in {"T", "True", "true", "1", "F", "False", "false", "0"},
                "invalid logical token")
        return token in {"T", "True", "true", "1"}
    require(kind == "S", "unknown property kind")
    return token


def snapshot(text):
    """Independent observer for this authored token/quoted-header packet only."""
    lines = text.splitlines()
    atoms = int(lines[0])
    require(atoms >= 0 and len(lines) == atoms + 2, "frame extent")
    header = lines[1]
    info = {}
    if "=" in header:
        for token in shlex.split(header):
            key, sep, value = token.partition("=")
            require(key not in info, "duplicate header key")
            info[key] = value if sep else "T"
    layout = info.pop("Properties", "species:S:1:pos:R:3").split(":")
    require(len(layout) % 3 == 0, "schema triplets")
    columns = [(layout[i], layout[i + 1], int(layout[i + 2]))
               for i in range(0, len(layout), 3)]
    require(columns[:2] == [("species", "S", 1), ("pos", "R", 3)], "schema prefix")
    require(all(c > 0 and k in {"R", "I", "S", "L"} for _, k, c in columns), "column type/width")
    lattice = info.pop("Lattice", None)
    if lattice is not None:
        lattice = [bits(v) for v in lattice.split()]
        require(len(lattice) == 9, "lattice width")
    rows = []
    for line in lines[2:]:
        tokens = line.split()
        require(len(tokens) == sum(c for _, _, c in columns), "atom row width")
        row, offset = [], 0
        for _, kind, count in columns:
            row.append([scalar(v, kind) for v in tokens[offset:offset + count]])
            offset += count
        rows.append(row)
    plain_comment = header if lattice is None and len(columns) == 2 and not info else None
    return dict(schema=columns, rows=rows, lattice=lattice,
                info=list(info.items()), plain_comment=plain_comment)


def observe(text, expected):
    try:
        result = snapshot(text)
        return dict(preserved=result == expected, snapshot=result)
    except (ValueError, IndexError) as error:
        return dict(preserved=False, error=f"{type(error).__name__}: {error}")


def binding_sources(baseline_ref, parent_ref):
    paths = ["crates/chematic-py/src/formats/xyz.rs", "crates/chematic-mol/src/xyz.rs",
             "crates/chematic-mol/src/xyz/extxyz.rs",
             "crates/chematic-mol/src/xyz/extxyz_validation.rs"]
    if baseline_ref is None:
        paths.append("crates/chematic-py/src/formats/xyz/properties.rs")
        return dict(reference="working tree", files=[identity(ROOT / p) for p in paths])
    commit = subprocess.check_output(["git", "rev-parse", baseline_ref], cwd=ROOT, text=True).strip()
    files = []
    for path in paths:
        data = subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)
        parent = subprocess.check_output(["git", "show", f"{parent_ref}:{path}"], cwd=ROOT)
        require(data == parent, f"baseline source differs from parent: {path}")
        files.append(dict(path=path, sha256=hashlib.sha256(data).hexdigest(), size_bytes=len(data)))
    return dict(reference=commit, matches_parent_sources=True, files=files)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline", action="store_true")
    parser.add_argument("--parent-ref", default="HEAD", help="comparison parent commit")
    parser.add_argument("--baseline-ref", help="commit used for the preserved local baseline build")
    parser.add_argument("--obabel", type=Path, default=Path("/opt/homebrew/bin/obabel"))
    args = parser.parse_args()
    parent_commit = subprocess.check_output(
        ["git", "rev-parse", args.parent_ref], cwd=ROOT, text=True).strip()
    require(args.baseline == bool(args.baseline_ref), "--baseline requires --baseline-ref")
    import chematic

    fixture = ROOT / "validation/extxyz_python_types/cases.json"
    cases = json.loads(fixture.read_text(encoding="utf-8"))["valid"]
    records = []
    for case in cases:
        expected = snapshot(case["text"])
        record = dict(case, expected=expected)
        try:
            frame = chematic.from_extxyz(case["text"])
            output = chematic.to_extxyz(**frame)
            record["python"] = dict(output=output, **observe(output, expected))
        except (ValueError, TypeError) as error:
            record["python"] = dict(preserved=False, error=f"{type(error).__name__}: {error}")
        if expected["rows"]:
            command = [str(args.obabel), "-ixyz", "-oxyz", "-ab"]
            process = bounded_run(command, case["text"].encode("utf-8"))
            process = {k: v.decode("utf-8", errors="replace") if isinstance(v, bytes) else v
                       for k, v in process.items()}
            output = process["stdout"] if (process["resource_status"] == "completed"
                                          and process["exit_code"] == 0) else ""
            record["openbabel_xyz"] = dict(command=command, process=process,
                                          **observe(output, expected))
        else:
            record["openbabel_xyz"] = dict(skipped="empty frame excluded from OB workflow")
        records.append(record)
    extension = importlib.import_module("chematic.chematic")
    summary = {
        "python": dict(preserved=sum(r["python"]["preserved"] for r in records), total=len(records)),
        "openbabel_xyz": dict(preserved=sum(r["openbabel_xyz"].get("preserved", False) for r in records),
                              total=sum("skipped" not in r["openbabel_xyz"] for r in records)),
    }
    result = dict(schema="chematic.extxyz-python-types-observation.v1", baseline=args.baseline,
                  scope="Synthetic typed Python read/write packet with an independent bounded text observer; "
                        "Open Babel is measured through plain XYZ, not claimed ASE ExtXYZ support; "
                        "no speed, chemical graph or general format ranking",
                  parent_commit=parent_commit,
                  fixture=identity(fixture), python_runtime=identity(extension.__file__),
                  python_version=sys.version, openbabel=identity(args.obabel),
                  openbabel_version=subprocess.check_output([str(args.obabel), "-V"], text=True).strip(),
                  binding_sources=binding_sources(args.baseline_ref, parent_commit),
                  harness_sources=[identity(Path(__file__)), identity(
                      ROOT / "scripts/checks/check_interchange_preservation.py")],
                  records=records, summary=summary)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))
    if not args.baseline:
        require(summary["python"]["preserved"] == len(records), "candidate lost declared state")


if __name__ == "__main__":
    main()
