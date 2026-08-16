#!/usr/bin/env python3
"""Issue #227 Phase 2: single-state `pipeline_v2_mmff94_strict` quality
summary. Reads one state's `chematic_pipeline_v2_vs_rdkit_dump.rs` output
(chematic_rows.jsonl) plus the paired RMSD/TFD JSONL files produced against
it (via `pipeline_v2_vs_rdkit_common_scorer.rs --pair` and
`pipeline_v2_vs_rdkit_tfd_227.py`), and reports every number named by its
producing entry point per the issue #227 Phase 2 directive: success rate,
typed-failure/timeout/crash breakdown, RMSD/TFD percentiles,
coverage@0.5/1.0/2.0A, declared/satisfied/violated/unevaluable stereo
counts, wall time.

`pipeline_v2_mmff94_strict` names the arm this whole report is about
throughout -- never conflated with the coverage-gate-only tooling
(`mmff94_strict_gate_remeasure_227.rs`) Phase 1 used for its OWN headline
numbers.

Usage:
    .venv/bin/python scripts/mmff94_bci_phase2_report_227.py \\
        --label STATE1_v0_16_0 \\
        --chematic-rows <path> \\
        [--rmsd <paired_rmsd.jsonl>] [--tfd <paired_tfd.jsonl>] \\
        [--subset <molecule_names.txt>] \\
        [--arm chematic_pipeline_v2_mmff94_strict]
"""

import argparse
import json
import re
import statistics
import sys


def load_jsonl(path):
    if path is None:
        return []
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


_FAILURE_CAUSE_RE = re.compile(r"^(\w+)(?:\((\w+)\b.*)?$", re.DOTALL)


def failure_cause_bucket(raw):
    """Collapses a full `{:?}`-formatted PipelineV2FailureCause (which
    embeds an entire Mmff94CoverageReport / error message inline, making
    each row's cause near-unique) down to its outer two variant names --
    e.g. `ForceField(MissingParameters(Mmff94CoverageReport { ... }))` ->
    `ForceField::MissingParameters`, `Timeout` -> `Timeout` -- readable
    aggregate counts without losing the raw string (still in the row
    JSONL itself for anyone who needs the full detail)."""
    if raw is None:
        return "unknown"
    m = _FAILURE_CAUSE_RE.match(raw)
    if not m:
        return raw
    outer, inner = m.group(1), m.group(2)
    return f"{outer}::{inner}" if inner else outer


def percentiles(values, ps=(0.5, 0.75, 0.9, 0.95)):
    if not values:
        return {f"p{int(p * 100)}": None for p in ps}
    s = sorted(values)
    n = len(s)
    out = {}
    for p in ps:
        idx = min(n - 1, int(round(p * (n - 1))))
        out[f"p{int(p * 100)}"] = s[idx]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--chematic-rows", required=True)
    ap.add_argument("--rmsd", default=None)
    ap.add_argument("--tfd", default=None)
    ap.add_argument("--subset", default=None, help="file, one molecule name per line, restricts the report to this subset")
    ap.add_argument("--arm", default="chematic_pipeline_v2_mmff94_strict")
    args = ap.parse_args()

    subset = None
    if args.subset:
        with open(args.subset) as f:
            subset = {line.strip() for line in f if line.strip()}

    rows = [r for r in load_jsonl(args.chematic_rows) if r.get("arm") == args.arm]
    if subset is not None:
        rows = [r for r in rows if r.get("name") in subset]

    total = len(rows)
    status_counts = {}
    failure_cause_counts = {}
    elapsed = []
    for r in rows:
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1
        if r["status"] in ("timeout", "typed_failure"):
            fc = failure_cause_bucket(r.get("failure_cause"))
            failure_cause_counts[fc] = failure_cause_counts.get(fc, 0) + 1
        if "elapsed_ms" in r:
            elapsed.append(r["elapsed_ms"])

    n_success = status_counts.get("success", 0)
    n_typed_failure = status_counts.get("typed_failure", 0)
    n_timeout = status_counts.get("timeout", 0)
    n_internal_error = status_counts.get("internal_error", 0)

    # Stereo, direct from the dump rows' own final_stereo_* fields (success
    # rows only -- failure rows may carry partial pre-failure stereo
    # evidence too, reported separately, never blended into "final").
    stereo_declared = stereo_satisfied = stereo_violated = stereo_unevaluable = 0
    for r in rows:
        if r["status"] == "success":
            stereo_declared += r.get("final_stereo_declared", 0) or 0
            stereo_satisfied += r.get("final_stereo_satisfied", 0) or 0
            stereo_violated += r.get("final_stereo_violations", 0) or 0
            stereo_unevaluable += r.get("final_stereo_unevaluable", 0) or 0

    rmsd_rows = [r for r in load_jsonl(args.rmsd) if r.get("status") == "paired_rmsd"]
    if subset is not None:
        rmsd_rows = [r for r in rmsd_rows if r.get("name") in subset]
    rmsd_values = [r["rmsd_symmetric_angstrom"] for r in rmsd_rows]

    tfd_rows = [r for r in load_jsonl(args.tfd) if r.get("status") == "paired_tfd"]
    if subset is not None:
        tfd_rows = [r for r in tfd_rows if r.get("name") in subset]
    tfd_values = [r["tfd"] for r in tfd_rows]

    coverage = {}
    if rmsd_values:
        for thresh in (0.5, 1.0, 2.0):
            coverage[f"coverage_at_{thresh}A"] = sum(1 for v in rmsd_values if v <= thresh) / total if total else None

    summary = {
        "label": args.label,
        "arm": args.arm,
        "subset_restricted": subset is not None,
        "subset_size_requested": len(subset) if subset is not None else None,
        "total_rows": total,
        "status_counts": status_counts,
        "success_rate": n_success / total if total else None,
        "failure_cause_counts": failure_cause_counts,
        "n_typed_failure": n_typed_failure,
        "n_timeout": n_timeout,
        "n_internal_error_crash": n_internal_error,
        "wall_time_ms": {
            "n": len(elapsed),
            "mean": statistics.mean(elapsed) if elapsed else None,
            **percentiles(elapsed),
        },
        "rmsd_symmetric_angstrom": {
            "n_paired": len(rmsd_values),
            "mean": statistics.mean(rmsd_values) if rmsd_values else None,
            **percentiles(rmsd_values),
        },
        "tfd": {
            "n_paired": len(tfd_values),
            "mean": statistics.mean(tfd_values) if tfd_values else None,
            **percentiles(tfd_values),
        },
        "coverage": coverage,
        "stereo": {
            "declared": stereo_declared,
            "satisfied": stereo_satisfied,
            "violated": stereo_violated,
            "unevaluable": stereo_unevaluable,
            "satisfaction_rate": stereo_satisfied / stereo_declared if stereo_declared else None,
        },
        "best_of_n": {
            "best_of_1": "measured (this report)",
            "best_of_10": "NOT MEASURABLE with existing tooling -- pipeline_v2_vs_rdkit_dump.rs runs exactly one embedding attempt sequence per (molecule, arm); MAX_ATTEMPTS=8 is retry-on-failure, not a best-of-N conformer selection. Building best-of-N tooling is out of scope for this measurement PR per the roadmap's 'use what's feasible' clause.",
        },
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
