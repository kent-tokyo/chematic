#!/usr/bin/env python3
"""Phase B validation: does chematic's RDKit-equivalent redundant-environment
suppression (`crates/chematic-fp/src/morgan_environment.rs`,
`SuppressRdkitRedundant` mode) emit the same set of `(atom_idx, radius)`
environments -- and the same sparse-count *shape* -- as RDKit's own default
(`includeRedundantEnvironments=False`) Morgan generator?

This asks a different question than PR #120's `ecfp_rdkit_environment_parity.py`,
which classifies where chematic's *pre-suppression* (`includeRedundantEnvironments
=True`-equivalent) expansion first diverges from RDKit, and whose
`redundant_environment_mismatch` bucket is purely RDKit-internal (rd.full vs
rd.default) -- it never references chematic's own emitted set, so re-running
it after adding suppression would not measure this fix at all (would keep
showing ~98% "divergence" regardless of whether chematic's suppression is
correct). This script instead directly compares chematic's *suppressed*
emitted-pair set, AND its raw-identifier sparse-count multiset shape, against
RDKit's *already-suppressed* (`default` variant) equivalents -- ignoring raw
hash VALUES throughout (same partition/set/multiset-shape-only discipline as
every prior script in this project; FNV-1a never numerically matches RDKit's
hash by construction).

This is a real GATE, not just a report: exits non-zero if any of the
following holds --
  - `duplicate_input_ids` nonzero on either side (pipeline desync)
  - any `parse_fail_chematic` / `parse_fail_rdkit` / `over_emission` /
    `under_emission` row (chematic parsing something RDKit can't or vice
    versa, or a one-directional pair mismatch, are never an expected shape)
  - any "both"-bucket mismatch, ANYWHERE in the input (not just the pinned
    fixtures), that is not *at least* a clean single-pair swap at the same
    radius -- the known representative-selection-tie residual shape. A
    multi-pair or cross-radius "both" mismatch would mean a real regression
    hiding behind the same bucket label, and must fail loudly.
  - any of the 4 pinned representative-swap fixtures NOT in one of the two
    acceptable states: (A) an outright exact pair-set match (the hash-
    matching milestone closing this gap in the future), or (B) the full
    well-shaped swap (pair-swap shape AND total-emitted-count AND
    sparse-count-shape AND unique-raw-identifier-count all preserved).
  - any of the 8 `sparse_count_mismatch` fixtures (from the Phase A
    diagnostic) with `sparse_count_shape_match` not True.
Only a REGRESSION (a worse-shaped mismatch, or a previously-passing fixture
degrading) is rejected -- a future improvement (e.g. exact match once a
hash-matching milestone lands) is explicitly accepted, not required.

Inputs:
  - chematic side: `crates/chematic-fp/examples/morgan_suppression_dump.rs`
    JSONL (requires the `diagnostics` feature), one row per molecule:
    {row_id, smiles, parse_ok, atom_count,
    emitted: [[atom_idx, radius, raw_environment_id], ...]}.
  - RDKit side: `scripts/gen_ecfp_rdkit_environment_oracle.py`'s existing,
    UNMODIFIED oracle JSONL (same script PR #120 used) -- this script reads
    `row["default"]["sparse_bit_info"]` (raw-hash -> [[atom,radius],...]) for
    the pair-set comparison, and `row["default"]["sparse_counts"]`
    (raw-hash -> count) for the count-shape comparison.

Usage:
    python scripts/ecfp_rdkit_suppression_parity.py \
        --chematic <morgan_suppression_dump.jsonl> \
        --rdkit-oracle <gen_ecfp_rdkit_environment_oracle.py --rows-out output> \
        --summary-out <out.json>
    python scripts/ecfp_rdkit_suppression_parity.py --self-test
"""

from __future__ import annotations

import argparse
import collections
import json
import sys

# The 9 representative-selection-tie residuals found in the full 5,041-input
# run are all a single-pair swap at the same radius (see mismatches in the
# summary). These 4 are pinned as permanent fixtures -- structural sanity
# checked (and GATED, see module docstring) on every run, not just observed
# once -- so a future hash-matching milestone can re-run this exact set and
# see the mismatch collapse to zero. "steroid-like" and "large/complex" per
# the residual's own shape, not a claim that either is acyclic -- both still
# contain rings, chosen as the largest/most topologically complex of the 9.
REPRESENTATIVE_SWAP_FIXTURES = {
    "CC(=O)NO": "small, non-ring, atom 1 vs atom 3 swap",
    "[10CH3][11CH3]": "isotope-labeled symmetric methyl pair swap",
    "C[C@]12CCC3C(CCC4=CC(=O)CC[C@@]43C3CO3)C1CCC2=O": "steroid-like fused-ring epoxide swap",
    "CCSc1ccnc(CSc2nc3ccccc3n2CC2CO2)c1C": "large polycyclic aromatic swap",
}

SPARSE_COUNT_MISMATCH_FIXTURES = [
    "C",
    "[CH4]",
    "[13CH4]",
    "[12CH4]",
    "[15NH4+]",
    "[15OH2].[16OH2]",
    "[Cl-]",
    "[CH3]",
]


def load_jsonl(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def chematic_pairs(row):
    return {(a, r) for a, r, _rid in row["emitted"]}


def chematic_counts(row):
    """raw_environment_id -> emission count, from chematic's suppressed dump."""
    counts = collections.Counter()
    for _a, _r, rid in row["emitted"]:
        counts[rid] += 1
    return counts


def rdkit_pairs(row):
    pairs = set()
    for envs in row["default"]["sparse_bit_info"].values():
        for a, r in envs:
            pairs.add((a, r))
    return pairs


def rdkit_counts(row):
    """raw_hash (str key) -> count, straight from RDKit's own sparse fingerprint."""
    return {k: v for k, v in row["default"]["sparse_counts"].items()}


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


def sparse_count_shape_match(chem_row, rd_row):
    """Real sparse-count-shape parity: does the MULTISET of per-raw-identifier
    emission counts match between chematic and RDKit? Compares
    `sorted(counts.values())` on both sides -- never the raw hash keys
    themselves (FNV-1a vs RDKit's hash never match by construction), only
    the shape of how many identifiers repeat how many times. This is the
    same "multiset shape" comparison PR #120 used for its own sparse-count
    check."""
    if not chem_row.get("parse_ok", True) or not rd_row.get("parse_ok", True):
        return None
    chem_shape = sorted(chematic_counts(chem_row).values())
    rd_shape = sorted(rdkit_counts(rd_row).values())
    return chem_shape == rd_shape


def swap_shape_check(chem_row, rd_row):
    """Characterizes how well a mismatch (if any) fits the known
    "representative-selection tie" shape: a clean single-pair swap at the
    same radius. Used both as the GLOBAL safety property every "both"-bucket
    mismatch in the full run must satisfy (see `is_acceptable_representative_swap`)
    and, more strictly, for the 4 pinned fixtures (see
    `is_acceptable_fixture_state`).

    NOTE on naming: `*_unique_raw_identifier_count` counts DISTINCT raw
    hash/environment IDENTIFIERS (`len({raw_environment_id})` on chematic's
    side, `len(sparse_bit_info)` on RDKit's), not literally-verified-distinct
    *cumulative bond environments* -- those are two different things
    whenever two structurally different environments happen to hash to the
    same raw identifier (a real, observed case: Kekule pyridine's
    `sparse_count_shape_match` mismatch traces to exactly this kind of
    accidental cross-radius hash collision). This check does not diagnose
    bond-index-sets directly, so it cannot claim bond-environment identity --
    only that the raw-identifier-count, total-emission-count, and
    count-multiset shape are each independently preserved.
    """
    bucket, extra, missing = classify(chem_row, rd_row)

    chem_unique_rids = len({rid for _a, _r, rid in chem_row["emitted"]}) if chem_row.get("parse_ok", True) else None
    rd_unique_rids = (
        len(rd_row["default"]["sparse_bit_info"]) if rd_row.get("parse_ok", True) else None
    )
    chem_total_emitted = len(chem_row["emitted"]) if chem_row.get("parse_ok", True) else None
    rd_total_emitted = (
        sum(len(v) for v in rd_row["default"]["sparse_bit_info"].values())
        if rd_row.get("parse_ok", True)
        else None
    )

    extra_radii = {r for _a, r in extra}
    missing_radii = {r for _a, r in missing}

    return {
        "bucket": bucket,
        "unique_raw_identifier_count_match": chem_unique_rids == rd_unique_rids,
        "chematic_unique_raw_identifier_count": chem_unique_rids,
        "rdkit_unique_raw_identifier_count": rd_unique_rids,
        "total_emitted_count_match": chem_total_emitted == rd_total_emitted,
        "chematic_total_emitted_count": chem_total_emitted,
        "rdkit_total_emitted_count": rd_total_emitted,
        "sparse_count_shape_match": sparse_count_shape_match(chem_row, rd_row),
        "is_exactly_one_pair_swap": len(extra) == 1 and len(missing) == 1,
        "swap_same_radius": bool(extra_radii) and bool(missing_radii) and extra_radii == missing_radii,
        "extra_pairs": sorted(extra),
        "missing_pairs": sorted(missing),
    }


def is_acceptable_representative_swap(result):
    """GLOBAL, weaker safety property every "both"-bucket mismatch in the
    full corpus must satisfy: exactly one pair replaced by another at the
    SAME radius -- never a larger, multi-pair, or cross-radius divergence
    hiding behind the same bucket label."""
    return result["bucket"] == "both" and result["is_exactly_one_pair_swap"] and result["swap_same_radius"]


def is_acceptable_fixture_state(result):
    """STRONGER property the 4 pinned representative-swap fixtures must
    satisfy: either an outright exact pair-set match (acceptable -- a future
    hash-matching milestone closing the gap), or the full well-shaped 1-1
    swap (pair-swap shape AND total-emitted-count AND sparse-count-shape AND
    unique-raw-identifier-count all preserved). Anything else is a
    regression."""
    if result["bucket"] == "exact_match":
        return True
    return (
        is_acceptable_representative_swap(result)
        and result["total_emitted_count_match"]
        and result["sparse_count_shape_match"]
        and result["unique_raw_identifier_count_match"]
    )


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
    sparse_shape_evaluated = 0
    sparse_shape_match = 0
    gate_failures = []

    representative_swap_results = {}
    sparse_count_fixture_results = {}

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

        if bucket in ("parse_fail_chematic", "parse_fail_rdkit", "over_emission", "under_emission"):
            gate_failures.append(
                f"row {idx} ({chem.get('smiles')!r}): bucket={bucket!r} is never an acceptable "
                f"shape (one-directional or parse-failure mismatch)"
            )
        elif bucket == "both":
            swap_result = swap_shape_check(chem, rd)
            if not is_acceptable_representative_swap(swap_result):
                gate_failures.append(
                    f"row {idx} ({chem.get('smiles')!r}): 'both'-bucket mismatch is NOT a clean "
                    f"single-pair same-radius swap (extra={swap_result['extra_pairs']}, "
                    f"missing={swap_result['missing_pairs']}) -- possible regression"
                )

        shape_match = sparse_count_shape_match(chem, rd)
        if shape_match is not None:
            sparse_shape_evaluated += 1
            if shape_match:
                sparse_shape_match += 1

        smi = chem.get("smiles")
        if smi in SPARSE_COUNT_MISMATCH_FIXTURES:
            sparse_count_fixture_results[smi] = shape_match
            if shape_match is not True:
                gate_failures.append(
                    f"sparse_count_mismatch fixture {smi!r}: sparse_count_shape_match={shape_match!r}, expected True"
                )
        if smi in REPRESENTATIVE_SWAP_FIXTURES:
            swap_result = swap_shape_check(chem, rd)
            representative_swap_results[smi] = swap_result
            if not is_acceptable_fixture_state(swap_result):
                gate_failures.append(
                    f"representative-swap fixture {smi!r}: not in an acceptable state "
                    f"(bucket={swap_result['bucket']!r}, "
                    f"is_exactly_one_pair_swap={swap_result['is_exactly_one_pair_swap']}, "
                    f"swap_same_radius={swap_result['swap_same_radius']}, "
                    f"total_emitted_count_match={swap_result['total_emitted_count_match']}, "
                    f"sparse_count_shape_match={swap_result['sparse_count_shape_match']}, "
                    f"unique_raw_identifier_count_match={swap_result['unique_raw_identifier_count_match']})"
                )

    total = len(chematic_rows)
    bucket_sum = sum(buckets.values())
    if bucket_sum != total:
        print(f"ACCOUNTING ERROR: bucket_sum={bucket_sum} != total_inputs={total}", file=sys.stderr)
        sys.exit(1)

    missing_swap_fixtures = set(REPRESENTATIVE_SWAP_FIXTURES) - set(representative_swap_results)
    missing_count_fixtures = set(SPARSE_COUNT_MISMATCH_FIXTURES) - set(sparse_count_fixture_results)
    if missing_swap_fixtures or missing_count_fixtures:
        print(
            f"FIXTURE COVERAGE ERROR: representative-swap fixtures not found in input: "
            f"{missing_swap_fixtures}; sparse-count fixtures not found: {missing_count_fixtures}",
            file=sys.stderr,
        )
        sys.exit(1)

    if duplicate_input_ids["chematic"] != 0 or duplicate_input_ids["rdkit_oracle"] != 0:
        gate_failures.append(f"duplicate_input_ids nonzero: {duplicate_input_ids}")

    sparse_count_shape = {
        "exact_match": sparse_shape_match,
        "mismatch": sparse_shape_evaluated - sparse_shape_match,
        "evaluated": sparse_shape_evaluated,
        "total_inputs": total,
        "exact_match_pct": (
            round(100.0 * sparse_shape_match / sparse_shape_evaluated, 4) if sparse_shape_evaluated else None
        ),
    }

    summary = {
        "total_inputs": total,
        "buckets": buckets,
        "exact_match_pct": round(100.0 * buckets["exact_match"] / total, 4) if total else None,
        "duplicate_input_ids": duplicate_input_ids,
        "sparse_count_shape": sparse_count_shape,
        "sparse_count_mismatch_fixtures_resolved": sparse_count_fixture_results,
        "representative_swap_fixtures": representative_swap_results,
        "mismatch_count": len(mismatches),
        "mismatches_sample": mismatches[:50],
        "gate_failures": gate_failures,
    }

    if gate_failures:
        print(f"GATE FAILED: {len(gate_failures)} violation(s):", file=sys.stderr)
        for f in gate_failures:
            print(f"  {f}", file=sys.stderr)
        return summary, mismatches, False

    return summary, mismatches, True


def _self_test_cases():
    """One hand-crafted (chematic, rdkit) row pair per bucket label, run
    through the real `classify()` -- not a reimplementation."""
    return [
        (
            "exact_match",
            {"parse_ok": True, "emitted": [[0, 0, 111], [1, 0, 222], [0, 1, 333]]},
            {"parse_ok": True, "default": {"sparse_bit_info": {"111": [[0, 0]], "222": [[1, 0]], "333": [[0, 1]]}}},
        ),
        (
            "over_emission",
            {"parse_ok": True, "emitted": [[0, 0, 111], [1, 0, 222], [0, 1, 333], [1, 1, 444]]},
            {"parse_ok": True, "default": {"sparse_bit_info": {"111": [[0, 0]], "222": [[1, 0]], "333": [[0, 1]]}}},
        ),
        (
            "under_emission",
            {"parse_ok": True, "emitted": [[0, 0, 111], [1, 0, 222]]},
            {"parse_ok": True, "default": {"sparse_bit_info": {"111": [[0, 0]], "222": [[1, 0]], "333": [[0, 1]]}}},
        ),
        (
            "both",
            {"parse_ok": True, "emitted": [[0, 0, 111], [1, 0, 222], [1, 1, 555]]},
            {"parse_ok": True, "default": {"sparse_bit_info": {"111": [[0, 0]], "222": [[1, 0]], "333": [[0, 1]]}}},
        ),
        (
            "parse_fail_chematic",
            {"parse_ok": False},
            {"parse_ok": True, "default": {"sparse_bit_info": {}}},
        ),
        (
            "parse_fail_rdkit",
            {"parse_ok": True, "emitted": [[0, 0, 111]]},
            {"parse_ok": False},
        ),
    ]


def run_self_test():
    failures = []
    for expected_bucket, chem, rd in _self_test_cases():
        bucket, _, _ = classify(chem, rd)
        if bucket != expected_bucket:
            failures.append(f"expected {expected_bucket!r}, got {bucket!r} for {chem!r} vs {rd!r}")

    # sparse_count_shape_match self-test: same shape, different raw ids/keys.
    chem_same_shape = {"parse_ok": True, "emitted": [[0, 0, 111], [0, 1, 222], [0, 2, 222]]}
    rd_same_shape = {
        "parse_ok": True,
        "default": {"sparse_counts": {"999": 1, "888": 2}},
    }
    if sparse_count_shape_match(chem_same_shape, rd_same_shape) is not True:
        failures.append("sparse_count_shape_match: expected True for matching shapes [1,2] vs [1,2]")

    chem_diff_shape = {"parse_ok": True, "emitted": [[0, 0, 111], [0, 1, 222]]}
    rd_diff_shape = {"parse_ok": True, "default": {"sparse_counts": {"999": 1, "888": 2}}}
    if sparse_count_shape_match(chem_diff_shape, rd_diff_shape) is not False:
        failures.append("sparse_count_shape_match: expected False for shapes [1,1] vs [1,2]")

    # swap_shape_check + is_acceptable_fixture_state self-test: a clean 1-1
    # same-radius swap with everything else preserved must be acceptable.
    chem_swap = {"emitted": [[0, 0, 111], [1, 1, 222]], "parse_ok": True}
    rd_swap = {
        "parse_ok": True,
        "default": {
            "sparse_bit_info": {"111": [[0, 0]], "222": [[2, 1]]},
            "sparse_counts": {"111": 1, "222": 1},
        },
    }
    good_result = swap_shape_check(chem_swap, rd_swap)
    if not is_acceptable_fixture_state(good_result):
        failures.append(f"is_acceptable_fixture_state: expected True for a clean 1-1 swap, got {good_result}")
    if not is_acceptable_representative_swap(good_result):
        failures.append(
            f"is_acceptable_representative_swap: expected True for a clean 1-1 swap, got {good_result}"
        )

    # Positive control: a "both" mismatch that is NOT a clean 1-1 same-radius
    # swap (two extra pairs, one at a DIFFERENT radius than the missing one)
    # must be REJECTED by both gate predicates -- proving the gate actually
    # catches a worse-shaped regression, not just labels every "both" bucket
    # as acceptable.
    chem_bad_swap = {"emitted": [[0, 0, 111], [1, 1, 222], [2, 2, 333]], "parse_ok": True}
    rd_bad_swap = {
        "parse_ok": True,
        "default": {
            "sparse_bit_info": {"111": [[0, 0]], "444": [[5, 1]]},
            "sparse_counts": {"111": 1, "444": 1},
        },
    }
    bad_result = swap_shape_check(chem_bad_swap, rd_bad_swap)
    if is_acceptable_representative_swap(bad_result):
        failures.append(
            f"is_acceptable_representative_swap: expected False for a 2-extra/1-missing/"
            f"cross-radius mismatch, got True ({bad_result})"
        )
    if is_acceptable_fixture_state(bad_result):
        failures.append(
            f"is_acceptable_fixture_state: expected False for a 2-extra/1-missing/"
            f"cross-radius mismatch, got True ({bad_result})"
        )

    if failures:
        print("SELF-TEST FAILED:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        sys.exit(1)
    print(
        f"self-test OK: {len(_self_test_cases())} bucket labels + shape/swap checks + "
        f"positive-control (bad-shape rejection) all reachable"
    )


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
    summary, mismatches, gate_passed = run(chematic_rows, rdkit_rows)

    print(json.dumps(summary["buckets"], indent=2))
    print(
        f"exact_match: {summary['buckets']['exact_match']}/{summary['total_inputs']} "
        f"({summary['exact_match_pct']}%)"
    )
    print("sparse_count_shape:", summary["sparse_count_shape"])
    print("sparse_count_mismatch fixtures resolved:", summary["sparse_count_mismatch_fixtures_resolved"])
    print("representative_swap fixtures:", json.dumps(summary["representative_swap_fixtures"], indent=2))

    if args.summary_out:
        with open(args.summary_out, "w") as f:
            json.dump({"summary": summary, "mismatches": mismatches}, f, indent=2, sort_keys=True)
        print(f"summary written to {args.summary_out}")

    if not gate_passed:
        sys.exit(1)


if __name__ == "__main__":
    main()
