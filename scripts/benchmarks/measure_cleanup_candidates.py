#!/usr/bin/env python3
"""Run frozen-loop equivalence checks and same-binary prepared-stage timings."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    args.output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    settings = {
        "CARGO_TARGET_DIR": str(args.target_dir.resolve()),
        "CARGO_BUILD_JOBS": "2",
        "CARGO_PROFILE_RELEASE_DEBUG": "0",
        "CARGO_PROFILE_RELEASE_LTO": "false",
        "CARGO_PROFILE_RELEASE_CODEGEN_UNITS": "16",
        "CARGO_INCREMENTAL": "0",
    }
    env.update(settings)
    base = ["cargo", "test", "-p", "chematic-smiles", "--offline", "--release", "--lib"]
    logs = {}
    # Timings run alone after validation, not concurrently with another build/test.
    for name, command in [
        ("differential", base + ["cleanup_candidate", "--", "--nocapture"]),
        (
            "timing",
            base + ["cleanup_candidates_timing", "--", "--ignored", "--nocapture"],
        ),
    ]:
        run = subprocess.run(
            command,
            cwd=root,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=300,
        )
        (args.output / f"{name}.log").write_text(run.stdout)
        if run.returncode:
            print(run.stdout)
            raise SystemExit(run.returncode)
        logs[name] = run.stdout
    rows = [
        dict(
            group=group,
            round=int(i),
            reference_seconds=float(a),
            candidate_seconds=float(b),
        )
        for group, i, a, b in re.findall(
            r"cleanup_timing\t(\w+)\t(\d+)\t([0-9.]+)\t([0-9.]+)", logs["timing"]
        )
    ]
    assert len(rows) == 18
    counts = dict(
        (group, int(count))
        for group, count in re.findall(r"cleanup_group\t(\w+)\t(\d+)", logs["timing"])
    )
    prepared, refused = map(
        int,
        re.search(
            r"cleanup differential: (\d+) prepared, (\d+) refused",
            logs["differential"],
        ).groups(),
    )
    assert counts["fixture"] == prepared
    assert counts["controls"] == 22
    summaries = {}
    for group in counts:
        selected = [row for row in rows if row["group"] == group]
        assert len(selected) == 9
        reference = statistics.median(row["reference_seconds"] for row in selected)
        candidate = statistics.median(row["candidate_seconds"] for row in selected)
        summaries[group] = dict(
            models=counts[group],
            median_reference_seconds=reference,
            median_candidate_seconds=candidate,
            median_reduction_percent=(1 - candidate / reference) * 100,
        )
    files = [
        "crates/chematic-smiles/src/rdkit/sanitize.rs",
        "crates/chematic-smiles/src/rdkit/sanitize/cleanup_candidate_tests.rs",
        "crates/chematic-smiles/tests/data/rdkit_canonical_smiles.tsv",
        "Cargo.lock",
    ]
    result = {
        "scope": "clean_up on raw internal models only; not parsing/sanitization end-to-end or a release benchmark",
        "baseline": "frozen pre-change loop in the same test binary",
        "platform": platform.platform(),
        "machine": platform.machine(),
        "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
        "cargo_settings": settings,
        "source_sha256": {
            file: hashlib.sha256((root / file).read_bytes()).hexdigest()
            for file in files
        },
        "checkout_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        "prepared_models": prepared,
        "refused_models": refused,
        "passes_per_round": 300,
        "rounds": rows,
        "summaries": summaries,
    }
    (args.output / "timing.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "rounds"}))


if __name__ == "__main__":
    main()
