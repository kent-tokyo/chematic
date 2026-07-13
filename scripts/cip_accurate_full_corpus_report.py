#!/usr/bin/env python3
"""
Full-corpus accounting for `chematic-cip`'s experimental accurate-CIP engine
(`assign_cip_accurate_experimental` vs its stable `..._without_mancude` reference point),
formalizing the ad hoc verification Milestone 3B-1b used to confirm the MANCUDE live-wiring
switch had zero regressions. Every future full-corpus CIP gate (Milestone 4 included) needs
the same rigor again, so this is a checked-in, reusable script rather than a one-off.

Consumes two TSV snapshots produced by
`cargo run -p chematic-cip --release --example corpus_snapshot -- --baseline|--candidate
<SMILES.csv> <out.tsv>` -- one row per candidate stereocenter
(`smiles\tatom_idx\tvalue`, value is R/S/E/Z, `skip:*`, or `ERR\t<message>`).

Prints exactly the accounting fields a numeric-discrepancy audit needs (see
docs/cip_accurate_rfc.md's Milestone 3B closeout entry for how these resolved a genuine
4055-vs-4047 / 4188-vs-4186 gap between this session's quick verification and the project's
canonical `cip_ground_truth.py`-based numbers): row counts at every filtering stage, both
engines' correctness against a freshly-regenerated RDKit oracle, and **two independently
computed** regression counts (not one number reported twice) -- Method A only inspects rows
where the two snapshots differ textually (a pure diff, classified via the oracle); Method B
independently re-derives an oracle label for *every* row and counts baseline-correct/
candidate-incorrect from scratch, without assuming or reusing Method A's diff set at all.

Usage:
    .venv/bin/python scripts/cip_accurate_full_corpus_report.py \\
        <baseline.tsv> <candidate.tsv> [SMILES.csv]
"""

import hashlib
import sys
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

BASELINE_ENGINE_MODE = "assign_cip_accurate_experimental_without_mancude"
CANDIDATE_ENGINE_MODE = "assign_cip_accurate_experimental"


def load_snapshot(path):
    rows = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            smi, idx = parts[0], int(parts[1])
            value = parts[2] if len(parts) > 2 else parts[1]
            rows[(smi, idx)] = value
    return rows


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print(__doc__)
        sys.exit(1)
    baseline_path, candidate_path = args[0], args[1]
    csv_path = args[2] if len(args) > 2 else f"{sys.argv[0]}/../../SMILES.csv"

    baseline = load_snapshot(baseline_path)
    candidate = load_snapshot(candidate_path)

    all_keys = set(baseline) | set(candidate)
    baseline_assigned_rows = sum(1 for v in baseline.values() if v in ("R", "S", "E", "Z"))
    candidate_assigned_rows = sum(1 for v in candidate.values() if v in ("R", "S", "E", "Z"))

    # Fresh oracle labels, one RDKit parse+CIP-label pass per distinct SMILES (shared
    # across every atom_idx row for that molecule).
    smis = {smi for smi, _ in all_keys}
    oracle_by_smiles = {}
    for smi in smis:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            oracle_by_smiles[smi] = None
            continue
        try:
            rdCIPLabeler.AssignCIPLabels(rd)
        except Exception:
            oracle_by_smiles[smi] = None
            continue
        oracle_by_smiles[smi] = {a.GetIdx(): a.GetProp("_CIPCode") for a in rd.GetAtoms() if a.HasProp("_CIPCode")}

    oracle_assigned_rows = 0
    oracle_unassigned_rows = 0
    excluded_rows = 0
    baseline_correct = 0
    candidate_correct = 0
    newly_correct = 0
    regressions_from_full_recompute = 0
    regression_examples = []

    for key in all_keys:
        smi, idx = key
        labels = oracle_by_smiles.get(smi)
        modern = labels.get(idx) if labels else None
        if modern is None:
            oracle_unassigned_rows += 1
            excluded_rows += 1
            continue
        oracle_assigned_rows += 1

        b_val = baseline.get(key)
        c_val = candidate.get(key)
        b_ok = b_val == modern
        c_ok = c_val == modern
        if b_ok:
            baseline_correct += 1
        if c_ok:
            candidate_correct += 1
        if c_ok and not b_ok:
            newly_correct += 1
        if b_ok and not c_ok:
            regressions_from_full_recompute += 1
            regression_examples.append((smi, idx, b_val, c_val, modern))

    # Method A: only rows where the two snapshots differ textually -- a pure diff,
    # independent of Method B's full-set recompute above.
    changed = [k for k in all_keys if baseline.get(k) != candidate.get(k)]
    regressions_from_diff = 0
    fixes_from_diff = 0
    neutral_from_diff = 0
    for key in changed:
        smi, idx = key
        labels = oracle_by_smiles.get(smi)
        modern = labels.get(idx) if labels else None
        if modern is None:
            continue
        b_ok = baseline.get(key) == modern
        c_ok = candidate.get(key) == modern
        if b_ok and not c_ok:
            regressions_from_diff += 1
        elif c_ok and not b_ok:
            fixes_from_diff += 1
        else:
            neutral_from_diff += 1

    corpus_sha256 = None
    try:
        with open(csv_path, "rb") as f:
            corpus_sha256 = hashlib.sha256(f.read()).hexdigest()
    except OSError:
        pass

    print("=== cip_accurate_full_corpus_report ===")
    print(f"input_rows (snapshot union):        {len(all_keys)}")
    print(f"baseline_assigned_rows:              {baseline_assigned_rows}")
    print(f"candidate_assigned_rows:             {candidate_assigned_rows}")
    print(f"oracle_assigned_rows:                {oracle_assigned_rows}")
    print(f"oracle_unassigned_rows:              {oracle_unassigned_rows}")
    print(f"excluded_rows (== oracle_unassigned): {excluded_rows}")
    print(f"baseline_correct:                    {baseline_correct}/{oracle_assigned_rows} "
          f"({100 * baseline_correct / oracle_assigned_rows:.2f}%)")
    print(f"candidate_correct:                   {candidate_correct}/{oracle_assigned_rows} "
          f"({100 * candidate_correct / oracle_assigned_rows:.2f}%)")
    print(f"newly_correct:                       {newly_correct}")
    print(f"regressions_from_full_recompute:     {regressions_from_full_recompute}  (Method B -- independent full-set oracle recompute)")
    print(f"regressions_from_diff:               {regressions_from_diff}  (Method A -- diff-set only, {len(changed)} changed rows: "
          f"{fixes_from_diff} fixes, {regressions_from_diff} regressions, {neutral_from_diff} neutral)")
    if regression_examples:
        print("regression examples (baseline correct, candidate wrong):")
        for smi, idx, b, c, m in regression_examples[:10]:
            print(f"  {smi} atom {idx}: baseline={b} candidate={c} modern={m}")
    print(f"corpus_sha256:                       {corpus_sha256}")
    print(f"rdkit_version:                       {Chem.rdBase.rdkitVersion}")
    print(f"baseline_engine_mode:                {BASELINE_ENGINE_MODE}")
    print(f"candidate_engine_mode:               {CANDIDATE_ENGINE_MODE}")


if __name__ == "__main__":
    main()
