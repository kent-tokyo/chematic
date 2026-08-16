#!/usr/bin/env python3
"""Issue #227 Phase 2, Step 5: full-corpus per-molecule transition table
between two states' `pipeline_v2_mmff94_strict` dumps (State 1 vs State 3
final, per the directive). Genuine per-molecule-ID join (Phase 1's own
`per_molecule_join_regressions`/`_improvements` naming and methodology,
`mmff94_torsion_gap_227_phase1_summary.json`'s `join_methodology` field) --
never aggregate-count arithmetic.

Pre-registered wall-clock-timeout-boundary jitter set (from
`validation/results/pipeline_v2_vs_rdkit_3point_paired_diff_summary.json`'s
`known_jitter_molecules`/`byte_identical_verification`, which shows even the
SAME commit re-run twice differs on 8/3181 rows, always on this exact
molecule set) -- flagged, not excluded, on any status flip so a jitter flip
isn't misread as caused by this PR's own changes.

Usage:
    .venv/bin/python scripts/mmff94_bci_phase2_transition_227.py \\
        --before <state1_chematic_rows.jsonl> --before-label STATE1 \\
        --after <state3_chematic_rows.jsonl> --after-label STATE3 \\
        [--before-rmsd ...] [--after-rmsd ...] \\
        [--before-tfd ...] [--after-tfd ...] \\
        --arm chematic_pipeline_v2_mmff94_strict
"""

import argparse
import json

KNOWN_JITTER_MOLECULES = {
    "chembl_tier_b_0166",
    "chembl_tier_b_0114",
    "chembl_tier_b_0117",
    "atorvastatin_fragment",
    "cholesterol",
}


def load_jsonl(path):
    if path is None:
        return []
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


def index_by_name(rows, arm):
    return {r["name"]: r for r in rows if r.get("arm") == arm}


def index_pairs_by_name(rows, key_field):
    return {r["name"]: r[key_field] for r in rows if key_field in r}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", required=True)
    ap.add_argument("--before-label", required=True)
    ap.add_argument("--after", required=True)
    ap.add_argument("--after-label", required=True)
    ap.add_argument("--before-rmsd", default=None)
    ap.add_argument("--after-rmsd", default=None)
    ap.add_argument("--before-tfd", default=None)
    ap.add_argument("--after-tfd", default=None)
    ap.add_argument("--arm", default="chematic_pipeline_v2_mmff94_strict")
    ap.add_argument("--subset", default=None, help="file, one molecule name per line, restricts the join to this subset")
    args = ap.parse_args()

    subset = None
    if args.subset:
        with open(args.subset) as f:
            subset = {line.strip() for line in f if line.strip()}

    before = index_by_name(load_jsonl(args.before), args.arm)
    after = index_by_name(load_jsonl(args.after), args.arm)
    if subset is not None:
        before = {k: v for k, v in before.items() if k in subset}
        after = {k: v for k, v in after.items() if k in subset}

    before_rmsd = index_pairs_by_name(
        [r for r in load_jsonl(args.before_rmsd) if r.get("status") == "paired_rmsd"],
        "rmsd_symmetric_angstrom",
    )
    after_rmsd = index_pairs_by_name(
        [r for r in load_jsonl(args.after_rmsd) if r.get("status") == "paired_rmsd"],
        "rmsd_symmetric_angstrom",
    )
    before_tfd = index_pairs_by_name(
        [r for r in load_jsonl(args.before_tfd) if r.get("status") == "paired_tfd"], "tfd"
    )
    after_tfd = index_pairs_by_name(
        [r for r in load_jsonl(args.after_tfd) if r.get("status") == "paired_tfd"], "tfd"
    )

    names = sorted(set(before) & set(after))
    transitions = {}
    per_molecule_join_regressions = []
    per_molecule_join_improvements = []
    rmsd_improved = []
    rmsd_worsened = []
    tfd_improved = []
    tfd_worsened = []
    stereo_newly_violated = []
    jitter_flips = []

    for name in names:
        b, a = before[name], after[name]
        bs, as_ = b["status"], a["status"]
        key = f"{bs}->{as_}"
        transitions[key] = transitions.get(key, 0) + 1

        if name in KNOWN_JITTER_MOLECULES and bs != as_:
            jitter_flips.append({"molecule": name, "before": bs, "after": as_})

        if bs == "success" and as_ != "success":
            per_molecule_join_regressions.append({"molecule": name, "before": bs, "after": as_})
        elif bs != "success" and as_ == "success":
            per_molecule_join_improvements.append({"molecule": name, "before": bs, "after": as_})

        if name in before_rmsd and name in after_rmsd:
            d = after_rmsd[name] - before_rmsd[name]
            if d < -1e-6:
                rmsd_improved.append({"molecule": name, "before": before_rmsd[name], "after": after_rmsd[name]})
            elif d > 1e-6:
                rmsd_worsened.append({"molecule": name, "before": before_rmsd[name], "after": after_rmsd[name]})

        if name in before_tfd and name in after_tfd:
            d = after_tfd[name] - before_tfd[name]
            if d < -1e-6:
                tfd_improved.append({"molecule": name, "before": before_tfd[name], "after": after_tfd[name]})
            elif d > 1e-6:
                tfd_worsened.append({"molecule": name, "before": before_tfd[name], "after": after_tfd[name]})

        if bs == "success" and as_ == "success":
            bv = b.get("final_stereo_violations", 0) or 0
            av = a.get("final_stereo_violations", 0) or 0
            if av > bv:
                stereo_newly_violated.append({"molecule": name, "before_violations": bv, "after_violations": av})

    result = {
        "before_label": args.before_label,
        "after_label": args.after_label,
        "arm": args.arm,
        "common_rows": len(names),
        "status_transitions": transitions,
        "per_molecule_join_regressions": per_molecule_join_regressions,
        "per_molecule_join_regressions_count": len(per_molecule_join_regressions),
        "per_molecule_join_improvements": per_molecule_join_improvements,
        "per_molecule_join_improvements_count": len(per_molecule_join_improvements),
        "join_methodology": "before/after rows keyed by molecule `name` into a dict per state, regression = name present in both with before-status success and after-status != success, improvement = before-status != success and after-status == success -- a genuine per-molecule-ID join, matching mmff94_torsion_gap_227_phase1_summary.json's own methodology field, not a before-count/after-count comparison.",
        "known_jitter_flips": jitter_flips,
        "known_jitter_flips_note": "molecules pre-registered (before this diff was computed) as wall-clock-timeout-boundary-sensitive per pipeline_v2_vs_rdkit_3point_paired_diff_summary.json's own byte_identical_verification (same commit re-run twice already flips these) -- any flip here is flagged, not silently folded into the regression/improvement counts' causal story.",
        "rmsd_improved_count": len(rmsd_improved),
        "rmsd_worsened_count": len(rmsd_worsened),
        "rmsd_worsened": rmsd_worsened,
        "tfd_improved_count": len(tfd_improved),
        "tfd_worsened_count": len(tfd_worsened),
        "tfd_worsened": tfd_worsened,
        "stereo_newly_violated_count": len(stereo_newly_violated),
        "stereo_newly_violated": stereo_newly_violated,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
