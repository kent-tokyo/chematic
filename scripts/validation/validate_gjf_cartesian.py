#!/usr/bin/env python3
"""Measure Cartesian GJF row integrity, keeping Open Babel read/write lanes separate."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import importlib
import json
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from check_interchange_preservation import bounded_run, identity, require, strict_json


def serial_run(command, text=""):
    run = bounded_run([str(arg) for arg in command], text.encode(), timeout=10, max_output=1_048_576)
    require(run["resource_status"] == "completed", f"process budget exceeded: {command}")
    return {**run, "stdout": run["stdout"].decode("utf-8", errors="replace"),
            "stderr": run["stderr"].decode("utf-8", errors="replace")}


def geometry_matches(observed, expected):
    return (observed["coords"] == expected["coords"] and
            observed["atom_count"] == len(expected["elements"]) and
            observed["charge"] == expected["charge"] and
            observed["multiplicity"] == expected["multiplicity"])


def probe_result(run, expected):
    if not run["stdout"].strip():
        return "refused" if run["exit_code"] == 2 else "process_error"
    if run["exit_code"] != 0:
        return "process_error"
    try:
        observed = strict_json(run["stdout"])
        return "match" if expected is not None and geometry_matches(observed, expected) else "emitted_mismatch"
    except (ValueError, KeyError, TypeError):
        return "unreadable_output"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-probe", type=Path, required=True)
    parser.add_argument("--candidate-probe", type=Path, required=True)
    parser.add_argument("--obabel", type=Path, default=Path("/opt/homebrew/bin/obabel"))
    parser.add_argument("--output", type=Path, default=ROOT / "validation/results/gjf-cartesian-2026-10-22.json")
    args = parser.parse_args()
    import chematic
    fixture = ROOT / "validation/gjf-cartesian-2026-10-22/cases.json"
    cases = strict_json(fixture.read_bytes())["cases"]
    counts, records = Counter(), []
    for case in cases:
        expected = case["expected"]
        for newline in ("LF", "CRLF"):
            text = case["input"].replace("\n", "\r\n") if newline == "CRLF" else case["input"]
            record = {"id": case["id"], "category": case["category"], "newline": newline}
            for engine, probe in (("baseline", args.baseline_probe), ("candidate_native", args.candidate_probe)):
                run = serial_run([probe], text)
                result = probe_result(run, expected)
                counts[f'{engine}.{case["category"]}.{result}'] += 1
                record[engine] = {"result": result, **run}
                if engine == "candidate_native":
                    require(result == ("match" if expected is not None else "refused"), str(record))
            try:
                g = chematic.parse_gjf_geometry(text)
            except ValueError as error:
                require(expected is None, f"candidate Python refused positive: {record}: {error}")
                record["candidate_python"] = {"result": "refused", "error": str(error)}
                counts[f'candidate_python.{case["category"]}.refused'] += 1
            else:
                require(expected is not None, f"candidate Python emitted outside profile: {record}")
                # Ordered species observed through a second format; coordinates and job data
                # are compared directly to independently authored fixture values.
                elements = [line.split()[0] for line in g["mol"].to_xyz(g["coords"]).splitlines()[2:]]
                actual = {key: g[key] for key in ("coords", "freeze_codes", "charge", "multiplicity")}
                actual["elements"] = elements
                require(actual == expected, f"candidate Python semantic mismatch: {record}: {actual}")
                record["candidate_python"] = {"result": "match", **actual}
                counts[f'candidate_python.{case["category"]}.match'] += 1
            records.append(record)

    # Open Babel advertises GJF write-only. Its invalid-format error is capability evidence,
    # not evidence that its parser validates or refuses malformed GJF rows.
    catalog = serial_run([args.obabel, "-L", "gjf"])
    read_control = serial_run([args.obabel, "-igjf", "-oxyz"], cases[0]["input"])
    require("Write-only" in catalog["stdout"] + catalog["stderr"], "unexpected GJF catalog: review comparison scope")
    require(read_control["exit_code"] != 0 and "Invalid input format" in read_control["stderr"],
            "Open Babel GJF read support changed: review comparison")
    writer_controls = []
    max_error = 0.0
    for case in cases:
        if case["category"] != "cartesian":
            continue
        expected = case["expected"]
        xyz = str(len(expected["elements"])) + "\noriginal geometry\n" + "".join(
            f'{element} {x!r} {y!r} {z!r}\n' for element, (x, y, z) in zip(expected["elements"], expected["coords"]))
        run = serial_run([args.obabel, "-ixyz", "-ogjf", "-xk", "# HF/STO-3G"], xyz)
        require(run["exit_code"] == 0 and run["stdout"], "Open Babel writer failed")
        g = chematic.parse_gjf_geometry(run["stdout"])
        elements = [line.split()[0] for line in g["mol"].to_xyz(g["coords"]).splitlines()[2:]]
        require(elements == expected["elements"], "writer control changed species/order")
        error = max(abs(a-b) for got, want in zip(g["coords"], expected["coords"]) for a, b in zip(got, want))
        require(len(g["coords"]) == len(expected["coords"]) and error <= 5.1e-6, "writer control outside five-decimal rounding envelope")
        max_error = max(max_error, error)
        writer_controls.append({"id": case["id"], "source_xyz": xyz, "max_coordinate_error": error, **run})
    packet = {"cases": records, "openbabel_catalog": catalog, "openbabel_read_control": read_control,
              "openbabel_writer_controls": writer_controls}
    packet_path = args.output.with_suffix(".packet.json.gz")
    packet_bytes = gzip.compress(json.dumps(packet, indent=2, allow_nan=False).encode(), mtime=0)
    packet_path.write_bytes(packet_bytes)
    native = importlib.import_module("chematic.chematic")
    report = {
        "schema": "chematic.gjf-cartesian-validation.v1", "utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Local source candidate, Cartesian geometry projection. Not complete job retention, Gaussian execution, or toolkit-wide superiority.",
        "baseline_commit": "7c34941f52c5e795fd3ca93b075aa31bc9fa748f", "platform": platform.platform(),
        "baseline_gaussian_source_sha256": "ffbf8c89870c5598fc853a647a8f23841e5736e6f7a7ab4e7ea08d4317c6562e",
        "probe_build": "cargo build --offline --locked -p chematic-mol, CARGO_PROFILE_DEV_DEBUG=0 CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=2; rustc --edition 2024 probe.rs --extern chematic_mol=<target>/debug/libchematic_mol.rlib -L dependency=<target>/debug/deps -o <probe>",
        "rustc": serial_run(["rustc", "-V"]),
        "fixture": identity(fixture), "packet": identity(packet_path), "counts": dict(sorted(counts.items())),
        "cases": len(cases), "newline_variants": ["LF", "CRLF"],
        "artifacts": {"baseline_probe": identity(args.baseline_probe), "candidate_probe": identity(args.candidate_probe),
                      "probe_source": identity(ROOT / "validation/gjf-cartesian-2026-10-22/probe.rs"),
                      "python_extension": identity(native.__file__), "openbabel": identity(args.obabel.resolve())},
        "source": [identity(ROOT / path) for path in ["crates/chematic-mol/src/gaussian.rs", "crates/chematic-mol/src/lib.rs", "crates/chematic-py/src/formats/gaussian.rs", "crates/chematic-py/src/formats.rs", "scripts/validation/validate_gjf_cartesian.py"]],
        "openbabel_version": serial_run([args.obabel, "-V"]),
        "openbabel_gjf_read": "write-only catalog and invalid-input-format error confirmed; stdout contains usage text, not molecular output; not scored as a validation lane",
        "openbabel_writer_controls": {"cases": len(writer_controls), "max_coordinate_error": max_error,
                                      "scope": "XYZ -> GJF -> candidate geometry reader; species/order and five-decimal writer rounding only. No source freeze/charge/spin retention claim."},
        "measurement_notes": ["Pilot incorrectly expected empty stdout for an unsupported Open Babel input format. It failed before writing a report; inspection found exit 1, Invalid input format on stderr, and usage text on stdout. Final comparison checks the explicit capability error and archives both pipes."],
        "process_limits": {"timeout_seconds": 10, "stdout_plus_stderr_bytes": 1_048_576, "memory_limit": "not imposed"},
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"counts": report["counts"], "writer_controls": report["openbabel_writer_controls"]}, indent=2))


if __name__ == "__main__":
    main()
