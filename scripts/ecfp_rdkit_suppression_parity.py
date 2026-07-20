#!/usr/bin/env python3
"""Phase B validation: does chematic's RDKit-equivalent redundant-environment
suppression (`ecfp_with_bitinfo_rdkit_environment_experimental`, via the new
`crates/chematic-fp/src/morgan_environment.rs`) emit the same set of
`(atom_idx, radius)` environments as RDKit's own default
(`includeRedundantEnvironments=False`) Morgan generator?

This asks a different question than PR #120's `ecfp_rdkit_environment_parity.py`,
which classifies where EmitAll-chematic first diverges from RDKit and whose
`redundant_environment_mismatch` bucket is purely RDKit-internal (rd.full vs
rd.default) -- it never references chematic's own emitted set, so re-running
it after adding suppression would not measure this fix at all (would keep
showing ~98% "divergence" regardless of whether chematic's suppression is
correct). This script instead directly compares chematic's *suppressed*
emitted-pair set against RDKit's *already-suppressed* (`default` variant)
emitted-pair set -- ignoring raw hash values throughout (same
partition/set-membership-only discipline as every prior script in this
project; FNV-1a never numerically matches RDKit's hash by construction).

Inputs:
  - chematic side: `crates/chematic-fp/examples/morgan_suppression_dump.rs`
    JSONL, one row per molecule: {row_id, smiles, parse_ok, atom_count,
    emitted: [[atom_idx, radius], ...]}.
  - RDKit side: `scripts/gen_ecfp_rdkit_environment_oracle.py`'s existing,
    UNMODIFIED oracle JSONL (same script PR #120 used) -- this script reads
    `row["default"]["sparse_bit_info"]` (raw-hash -> [[atom,radius],...]) and
    flattens it to get RDKit's real suppressed emitted-pair set.

Usage:
    python scripts/ecfp_rdkit_suppression_parity.py \
        --chematic <morgan_suppression_dump.jsonl> \
        --rdkit-oracle <gen_ecfp_rdkit_environment_oracle.py --rows-out output> \
        --summary-out <out.json>
    python scripts/ecfp_rdkit_suppression_parity.py --self-test
"""

from __future__ import annotations

import argparse
import json
import sys


def load_jsonl(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def chematic_pairs(row):
    return {(a, r) for a, r in row["emitted"]}


def rdkit_pairs(row):
    pairs = set()
    for envs in row["default"]["sparse_bit_info"].values():
        for a, r in envs:
            pairs.add((a, r))
    return pairs


def classify(chem, rd):
    """Returns (bucket, extra_chematic, missing_chematic) for one molecule.

    bucket is one of: parse_fail_chematic, parse_fail_rdkit, exact_match,
    over_emission (chematic emits pairs RDKit doesn't -- a suppression bug:
    chematic failed to suppress something RDKit does), under_emission
    (RDKit emits pairs chematic doesn't -- either a suppression bug in the
    other direction, or the expected hash-dependent representative-tie
    residual), both (both directions differ simultaneously).
    """
    if not chem.get("parse_ok", True):
        return "parse_fail_chematic", set(), set()
    if not rd.get("parse_ok", True):
        return "parse_fail_rdkit", set(), set()

    chem_set = chematic_pairs(chem)
    rd_set = rdkit_pairs(rd)
    extra = chem_set - rd_set
    missing = rd_set - chem_set

    if not extra and not missing:
        return "exact_match", extra, missing
    if extra and missing:
        return "both", extra, missing
    if extra:
        return "over_emission", extra, missing
    return "under_emission", extra, missing


def sparse_counts_match(chem_row, rd_row):
    """Secondary check (plan item 4): do the 8 single-atom/degree-0
    `sparse_count_mismatch` fixtures from PR #120 now also match RDKit's
    real `sparse_counts` shape? Compares counts-per-radius derived from the
    emitted-pair sets (both sides), not raw hash values."""
    if not chem_row.get("parse_ok", True) or not rd_row.get("parse_ok", True):
        return None
    chem_radii = sorted(r for _, r in chematic_pairs(chem_row))
    rd_radii = sorted(r for _, r in rdkit_pairs(rd_row))
    return chem_radii == rd_radii


def run(chematic_rows, rdkit_rows):
    if len(chematic_rows) != len(rdkit_rows):
        print(
            f"PIPELINE ERROR: row count mismatch chematic={len(chematic_rows)} "
            f"rdkit={len(rdkit_rows)}",
            file=sys.stderr,
        )
        sys.exit(1)

    chem_ids = [r["row_id"] for r in chematic_rows]
    rd_ids = [r["row_id"] for r in rdkit_rows]
    duplicate_input_ids = {
        "chematic": len(chem_ids) - len(set(chem_ids)),
        "rdkit_oracle": len(rd_ids) - len(set(rd_ids)),
    }

    buckets = {
        "exact_match": 0,
        "over_emission": 0,
        "under_emission": 0,
        "both": 0,
        "parse_fail_chematic": 0,
        "parse_fail_rdkit": 0,
    }
    mismatches = []

    for idx, (chem, rd) in enumerate(zip(chematic_rows, rdkit_rows)):
        if chem.get("row_id") != idx or rd.get("row_id") != idx:
            print(f"PIPELINE ERROR at position {idx}: row_id out of sync", file=sys.stderr)
            sys.exit(1)
        if chem.get("smiles") != rd.get("smiles"):
            print(
                f"PIPELINE ERROR at row {idx}: chematic smiles={chem.get('smiles')!r} "
                f"!= rdkit smiles={rd.get('smiles')!r}",
                file=sys.stderr,
            )
            sys.exit(1)

        bucket, extra, missing = classify(chem, rd)
        buckets[bucket] += 1
        if bucket != "exact_match":
            mismatches.append(
                {
                    "row_id": idx,
                    "smiles": chem.get("smiles"),
                    "bucket": bucket,
                    "extra_chematic_pairs": sorted(extra),
                    "missing_chematic_pairs": sorted(missing),
                }
            )

    total = len(chematic_rows)
    bucket_sum = sum(buckets.values())
    if bucket_sum != total:
        print(f"ACCOUNTING ERROR: bucket_sum={bucket_sum} != total_inputs={total}", file=sys.stderr)
        sys.exit(1)

    sparse_count_fixtures = [
        "C",
        "[CH4]",
        "[13CH4]",
        "[12CH4]",
        "[15NH4+]",
        "[15OH2].[16OH2]",
        "[Cl-]",
        "[CH3]",
    ]
    sparse_count_results = {}
    for chem, rd in zip(chematic_rows, rdkit_rows):
        smi = chem.get("smiles")
        if smi in sparse_count_fixtures:
            sparse_count_results[smi] = sparse_counts_match(chem, rd)

    summary = {
        "total_inputs": total,
        "buckets": buckets,
        "exact_match_pct": round(100.0 * buckets["exact_match"] / total, 4) if total else None,
        "duplicate_input_ids": duplicate_input_ids,
        "sparse_count_mismatch_fixtures_now_resolved": sparse_count_results,
        "mismatch_count": len(mismatches),
        "mismatches_sample": mismatches[:50],
    }
    return summary, mismatches


def _self_test_cases():
    """One hand-crafted (chematic, rdkit) row pair per bucket label, run
    through the real `classify()` -- not a reimplementation."""
    return [
        (
            "exact_match",
            {"parse_ok": True, "emitted": [[0, 0], [1, 0], [0, 1]]},
            {"parse_ok": True, "default": {"sparse_bit_info": {"111": [[0, 0]], "222": [[1, 0]], "333": [[0, 1]]}}},
        ),
        (
            "over_emission",
            {"parse_ok": True, "emitted": [[0, 0], [1, 0], [0, 1], [1, 1]]},
            {"parse_ok": True, "default": {"sparse_bit_info": {"111": [[0, 0]], "222": [[1, 0]], "333": [[0, 1]]}}},
        ),
        (
            "under_emission",
            {"parse_ok": True, "emitted": [[0, 0], [1, 0]]},
            {"parse_ok": True, "default": {"sparse_bit_info": {"111": [[0, 0]], "222": [[1, 0]], "333": [[0, 1]]}}},
        ),
        (
            "both",
            {"parse_ok": True, "emitted": [[0, 0], [1, 0], [1, 1]]},
            {"parse_ok": True, "default": {"sparse_bit_info": {"111": [[0, 0]], "222": [[1, 0]], "333": [[0, 1]]}}},
        ),
        (
            "parse_fail_chematic",
            {"parse_ok": False},
            {"parse_ok": True, "default": {"sparse_bit_info": {}}},
        ),
        (
            "parse_fail_rdkit",
            {"parse_ok": True, "emitted": [[0, 0]]},
            {"parse_ok": False},
        ),
    ]


def run_self_test():
    failures = []
    for expected_bucket, chem, rd in _self_test_cases():
        bucket, _, _ = classify(chem, rd)
        if bucket != expected_bucket:
            failures.append(f"expected {expected_bucket!r}, got {bucket!r} for {chem!r} vs {rd!r}")
    if failures:
        print("SELF-TEST FAILED:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        sys.exit(1)
    print(f"self-test OK: {len(_self_test_cases())} bucket labels all reachable via classify()")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chematic", default=None)
    p.add_argument("--rdkit-oracle", default=None)
    p.add_argument("--summary-out", default=None)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()

    if args.self_test:
        run_self_test()
        return

    if not args.chematic or not args.rdkit_oracle:
        print("--chematic and --rdkit-oracle are required (or pass --self-test)", file=sys.stderr)
        sys.exit(1)

    chematic_rows = load_jsonl(args.chematic)
    rdkit_rows = load_jsonl(args.rdkit_oracle)
    summary, mismatches = run(chematic_rows, rdkit_rows)

    print(json.dumps(summary["buckets"], indent=2))
    print(f"exact_match: {summary['buckets']['exact_match']}/{summary['total_inputs']} "
          f"({summary['exact_match_pct']}%)")
    print("sparse_count_mismatch fixtures now resolved:", summary["sparse_count_mismatch_fixtures_now_resolved"])

    if args.summary_out:
        with open(args.summary_out, "w") as f:
            json.dump({"summary": summary, "mismatches": mismatches}, f, indent=2, sort_keys=True)
        print(f"summary written to {args.summary_out}")


if __name__ == "__main__":
    main()
