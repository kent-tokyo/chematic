#!/usr/bin/env python3
"""Append a veridict drift-detection section to docs/validation.md.

Compares per-molecule |chematic - RDKit| error between a previous and the
current scripts/bench5k.py --json run (see the "deltas" field added there),
using veridict's mean-diff bootstrap. This flags cases where individual
molecules stay within their assert tolerance but the *average* error has
drifted since the last recorded run — not a per-PR gate (see the plan for
why: both sides are deterministic, so a PR that doesn't touch descriptor
code always diffs to zero and the corpus can vary between manual/scheduled
runs). Report-only: never changes this script's own exit code based on the
veridict verdict.

Usage:
    python3 scripts/veridict_accuracy_report.py <previous.json> <current.json> [--out docs/validation.md]

If <previous.json> doesn't exist (no prior run yet), a short "no baseline"
note is appended instead of running veridict.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

# (descriptor key in bench5k.py's "deltas" dict, --fail-below starting point)
# Thresholds are starting points, not tuned promises.
DESCRIPTORS = [
    ("tpsa", -0.02),
    ("logp", -0.005),
]


def run_veridict(pairs, fail_below):
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".jsonl", delete=False) as f:
        for smi, (prev_err, cur_err) in pairs.items():
            f.write(json.dumps({"id": smi, "baseline": -abs(prev_err), "candidate": -abs(cur_err)}) + "\n")
        jsonl_path = f.name

    report_path = jsonl_path.replace(".jsonl", "_report.json")
    subprocess.run(
        [
            "veridict", "compare", jsonl_path,
            "--metric", "mean-diff",
            "--confidence", "0.95",
            "--pass-above", "0",
            f"--fail-below={fail_below}",
            "--report-json", report_path,
        ],
        capture_output=True,
    )
    return json.loads(Path(report_path).read_text())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("previous_json")
    ap.add_argument("current_json")
    ap.add_argument("--out", default="docs/validation.md")
    args = ap.parse_args()

    lines = ["", "## Accuracy drift vs previous run (veridict)", ""]

    prev_path = Path(args.previous_json)
    if not prev_path.exists():
        lines.append("_No previous run recorded yet — drift detection starts next run._")
        Path(args.out).open("a").write("\n".join(lines) + "\n")
        return

    previous = json.loads(prev_path.read_text())
    current = json.loads(Path(args.current_json).read_text())

    for key, fail_below in DESCRIPTORS:
        prev_deltas = previous.get("deltas", {}).get(key, {})
        cur_deltas = current.get("deltas", {}).get(key, {})
        shared = set(prev_deltas) & set(cur_deltas)

        if not shared:
            lines.append(f"- **{key.upper()}**: no overlapping molecules with the previous run, skipped.")
            continue

        pairs = {smi: (prev_deltas[smi], cur_deltas[smi]) for smi in shared}
        report = run_veridict(pairs, fail_below)
        verdict = report.get("verdict", "unknown")
        effect = report.get("effect")
        ci_low = report.get("ci_low")
        ci_high = report.get("ci_high")
        lines.append(
            f"- **{key.upper()}** ({len(pairs)} molecules): **{verdict}** "
            f"(mean |error| change {effect:+.4f}, 95% CI [{ci_low:+.4f}, {ci_high:+.4f}])"
        )

    Path(args.out).open("a").write("\n".join(lines) + "\n")
    print(f"Appended accuracy drift section to {args.out}")


if __name__ == "__main__":
    main()
