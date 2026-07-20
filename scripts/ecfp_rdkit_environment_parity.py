#!/usr/bin/env python3
"""
Morgan/ECFP RDKit environment-parity comparator (the M2 diagnostic PR).

Extends scripts/ecfp_rdkit_invariant_parity.py's radius-0-only atom-invariant
PARTITION check (PR #110) to every stage between that and the final 2048-bit
fingerprint: radius-1/2 environment membership, radius-1/2 identifier
partition, redundant-environment suppression, sparse counts, folding,
bitInfo. For every molecule, classifies exactly ONE first-divergence stage
from an ordered enum (STAGES below) -- never a raw-hash-value comparison
(chematic uses FNV-1a, RDKit its own hash; only partitions/set-membership are
ever compared).

Inputs are two JSONL files produced upstream, in the SAME line order (this is
what avoids any SMILES-based join ambiguity -- duplicate/symmetric SMILES are
handled correctly because rows are matched by POSITION, not by string key):
  - chematic trace: `cargo run -p chematic-fp --release --features diagnostics
    --example morgan_rdkit_environment_trace -- <combined.csv> <out.jsonl>`
  - RDKit oracle: `python scripts/gen_ecfp_rdkit_environment_oracle.py
    --corpus <corpus.csv> --fixtures <f1.csv> <f2.csv> --rows-out <out.jsonl>`
    (its internal source order is corpus-then-fixtures-in-argument-order,
    identical to how <combined.csv> above must be built: cat corpus then each
    fixture file, in that same order, skipping blank lines both sides skip.)

Evaluation short-circuits at the first True stage (STAGES is a strict
precedence order): once an earlier, more fundamental divergence fires (e.g.
atom counts don't even match), later per-atom comparisons aren't meaningful
and are recorded as None ("not evaluated"), not False.

Usage:
    .venv/bin/python scripts/ecfp_rdkit_environment_parity.py \\
        --chem-trace <chematic_trace.jsonl> --rdkit-oracle <rdkit_oracle.jsonl> \\
        --summary-out validation/ecfp_rdkit_environment_parity_summary.json \\
        --rows-out validation/ecfp_rdkit_environment_parity_rows.jsonl \\
        --first-divergence-out validation/ecfp_rdkit_environment_parity_first_divergence.tsv \\
        [--rows-scope {fixtures,all}]

    .venv/bin/python scripts/ecfp_rdkit_environment_parity.py --self-test
"""

import argparse
import json
import random
import statistics
import sys
from collections import defaultdict

STAGES = [
    "chematic_parse_fail",
    "rdkit_parse_fail",
    "atom_count_mismatch",
    "initial_invariant_mismatch",
    "radius1_environment_membership_mismatch",
    "radius1_identifier_mismatch",
    "radius2_environment_membership_mismatch",
    "radius2_identifier_mismatch",
    "redundant_environment_mismatch",
    "sparse_count_mismatch",
    "folding_mismatch",
    "bit_info_mismatch",
]


def partition(values):
    """[v0, v1, ...] -> sorted tuple of sorted-index-tuples grouped by equal
    value -- same methodology as ecfp_rdkit_invariant_parity.py's partition():
    hash-VALUE-independent, only equivalence-CLASS membership is compared."""
    groups = defaultdict(list)
    for i, v in enumerate(values):
        groups[v].append(i)
    return tuple(sorted(tuple(sorted(g)) for g in groups.values()))


def _trace_by_radius(chem, radius):
    return {e["atom_idx"]: e for e in chem["trace"] if e["radius"] == radius}


def _pairs_from_bit_info(bit_info):
    pairs = set()
    for atom_radius_list in bit_info.values():
        for atom, radius in atom_radius_list:
            pairs.add((atom, radius))
    return pairs


def _ids_at_radius_from_bit_info(bit_info, radius):
    out = {}
    for raw_id_str, atom_radius_list in bit_info.items():
        for atom, r in atom_radius_list:
            if r == radius:
                out[atom] = int(raw_id_str)
    return out


def evaluate_stages(chem, rd):
    """Returns (flags: {stage: True/False/None}, evidence: {stage: {...}}).
    See module docstring for the short-circuit semantics."""
    flags = dict.fromkeys(STAGES)
    evidence = {}

    flags["chematic_parse_fail"] = not chem.get("parse_ok", False)
    if flags["chematic_parse_fail"]:
        return flags, evidence

    flags["rdkit_parse_fail"] = not rd.get("parse_ok", False)
    if flags["rdkit_parse_fail"]:
        return flags, evidence

    atom_count_ok = chem["atom_count"] == rd["atom_count"]
    atomic_numbers_ok = atom_count_ok and chem["atomic_numbers"] == rd["atomic_numbers"]
    flags["atom_count_mismatch"] = not (atom_count_ok and atomic_numbers_ok)
    if flags["atom_count_mismatch"]:
        evidence["atom_count_mismatch"] = {
            "chematic_atom_count": chem.get("atom_count"),
            "rdkit_atom_count": rd.get("atom_count"),
            "atomic_numbers_correspond": atomic_numbers_ok if atom_count_ok else None,
        }
        return flags, evidence

    n = chem["atom_count"]
    chem_r0 = _trace_by_radius(chem, 0)
    chem_r0_ids = [chem_r0[i]["raw_environment_id"] for i in range(n)]
    flags["initial_invariant_mismatch"] = partition(chem_r0_ids) != partition(
        rd["connectivity_invariants"]
    )
    if flags["initial_invariant_mismatch"]:
        return flags, evidence

    for r in (1, 2):
        mem_key = f"radius{r}_environment_membership_mismatch"
        id_key = f"radius{r}_identifier_mismatch"
        chem_r = _trace_by_radius(chem, r)
        rd_balls = rd["atom_balls"][str(r)]

        mismatch_atom = next(
            (i for i in range(n) if chem_r[i]["atom_ball"] != rd_balls[str(i)]), None
        )
        flags[mem_key] = mismatch_atom is not None
        if flags[mem_key]:
            evidence[mem_key] = {
                "atom_idx": mismatch_atom,
                "radius": r,
                "chematic_ball": chem_r[mismatch_atom]["atom_ball"],
                "rdkit_ball": rd_balls[str(mismatch_atom)],
            }
            return flags, evidence

        # Compared against RDKit's "full" (includeRedundantEnvironments=True)
        # variant, never "default" -- default already dropped suppressed
        # pairs, which would look like a missing identifier here and get
        # misfiled as radiusR_identifier_mismatch instead of the later,
        # correctly-attributed redundant_environment_mismatch stage.
        rd_ids_at_r = _ids_at_radius_from_bit_info(rd["full"]["sparse_bit_info"], r)
        chem_ids = [chem_r[i]["raw_environment_id"] for i in range(n)]
        rd_ids = [rd_ids_at_r.get(i) for i in range(n)]
        flags[id_key] = partition(chem_ids) != partition(rd_ids)
        if flags[id_key]:
            evidence[id_key] = {
                "radius": r,
                "chematic_partition_class_count": len(partition(chem_ids)),
                "rdkit_partition_class_count": len(partition(rd_ids)),
            }
            return flags, evidence

    # RDKit's own actual suppression decision -- diff of its two generator
    # variants -- not inferred from chematic-side ball stability.
    rd_full_pairs = _pairs_from_bit_info(rd["full"]["sparse_bit_info"])
    rd_default_pairs = _pairs_from_bit_info(rd["default"]["sparse_bit_info"])
    suppressed = rd_full_pairs - rd_default_pairs
    flags["redundant_environment_mismatch"] = len(suppressed) > 0
    if flags["redundant_environment_mismatch"]:
        atom, radius = sorted(suppressed)[0]
        evidence["redundant_environment_mismatch"] = {
            "atom_idx": atom,
            "radius": radius,
            "suppressed_pair_count": len(suppressed),
        }
        return flags, evidence

    # Multiset SHAPE only (sorted list of per-raw-id multiplicities) -- raw
    # ids are never comparable cross-implementation.
    chem_counts = defaultdict(int)
    for e in chem["trace"]:
        chem_counts[e["raw_environment_id"]] += 1
    chem_shape = sorted(chem_counts.values())
    rd_shape = sorted(rd["default"]["sparse_counts"].values())
    flags["sparse_count_mismatch"] = chem_shape != rd_shape
    if flags["sparse_count_mismatch"]:
        evidence["sparse_count_mismatch"] = {
            "chematic_shape_len": len(chem_shape),
            "rdkit_shape_len": len(rd_shape),
        }
        return flags, evidence

    chem_folded_bits = {e["folded_bit"] for e in chem["trace"]}
    rd_folded_bits = set(rd["default"]["folded_on_bits"])
    flags["folding_mismatch"] = chem_folded_bits != rd_folded_bits
    if flags["folding_mismatch"]:
        evidence["folding_mismatch"] = {
            "chematic_bit_count": len(chem_folded_bits),
            "rdkit_bit_count": len(rd_folded_bits),
        }
        return flags, evidence

    chem_bitinfo = defaultdict(list)
    for e in chem["trace"]:
        chem_bitinfo[e["folded_bit"]].append([e["atom_idx"], e["radius"]])
    for v in chem_bitinfo.values():
        v.sort()
    rd_bitinfo = {int(k): sorted(v) for k, v in rd["default"]["folded_bit_info"].items()}
    flags["bit_info_mismatch"] = dict(chem_bitinfo) != rd_bitinfo

    return flags, evidence


def classify(flags):
    for s in STAGES:
        if flags[s]:
            return s
    return "exact_match"


def _tanimoto(a, b):
    if not a and not b:
        return 1.0
    inter = len(a & b)
    union = len(a) + len(b) - inter
    return inter / union if union else 0.0


def tanimoto_correlation(chem_rows, rd_rows, sample_size=300, seed=42):
    """Reference (non-gating) measurement, same methodology as
    scripts/ecfp_rdkit_invariants_fingerprint_ref.py: pairwise-Tanimoto
    Pearson correlation between chematic's real folded 2048-bit set (from
    the trace) and RDKit's real folded 2048-bit set (`default.folded_on_bits`),
    for every molecule that parsed on both sides."""
    pairs = []
    for chem, rd in zip(chem_rows, rd_rows):
        if not chem.get("parse_ok") or not rd.get("parse_ok"):
            continue
        chem_bits = frozenset(e["folded_bit"] for e in chem["trace"])
        rd_bits = frozenset(rd["default"]["folded_on_bits"])
        pairs.append((chem_bits, rd_bits))

    rng = random.Random(seed)
    sample = pairs if len(pairs) <= sample_size else rng.sample(pairs, sample_size)

    chem_sims, rd_sims = [], []
    n = len(sample)
    for i in range(n):
        for j in range(i + 1, n):
            chem_sims.append(_tanimoto(sample[i][0], sample[j][0]))
            rd_sims.append(_tanimoto(sample[i][1], sample[j][1]))

    if len(chem_sims) < 2:
        return None
    return round(statistics.correlation(chem_sims, rd_sims), 4)


def load_jsonl(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


# --------------------------------------------------------------------------
# Self-test: one hand-crafted (chem, rd) pair per STAGES member, run through
# the exact classify(evaluate_stages(...)) used on real data.
# --------------------------------------------------------------------------


def _self_test_cases():
    def trace(entries):
        return [
            {
                "atom_idx": a,
                "radius": r,
                "raw_environment_id": rid,
                "folded_bit": fb,
                "emitted": True,
                "atom_ball": ball,
            }
            for (a, r, rid, fb, ball) in entries
        ]

    cases = {}

    cases["chematic_parse_fail"] = (
        {"smiles": "X", "parse_ok": False},
        {"smiles": "X", "parse_ok": True},
    )

    cases["rdkit_parse_fail"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 1,
            "atomic_numbers": [6],
            "trace": trace([(0, 0, 1, 1, [0])]),
        },
        {"smiles": "X", "parse_ok": False},
    )

    cases["atom_count_mismatch"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": trace([(0, 0, 1, 1, [0]), (1, 0, 1, 1, [1])]),
        },
        {"smiles": "X", "parse_ok": True, "atom_count": 1, "atomic_numbers": [6]},
    )

    cases["initial_invariant_mismatch"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": trace([(0, 0, 100, 1, [0]), (1, 0, 100, 1, [1])]),  # same class both atoms
        },
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "connectivity_invariants": [5, 7],  # different classes
        },
    )

    common_r0 = trace([(0, 0, 10, 1, [0]), (1, 0, 20, 2, [1])])

    cases["radius1_environment_membership_mismatch"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": common_r0 + trace([(0, 1, 30, 3, [0, 1]), (1, 1, 30, 3, [0, 1])]),
        },
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "connectivity_invariants": [10, 20],
            "atom_balls": {"1": {"0": [0], "1": [0, 1]}},  # atom 0's ball disagrees
        },
    )

    cases["radius1_identifier_mismatch"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": common_r0 + trace([(0, 1, 30, 3, [0, 1]), (1, 1, 40, 4, [0, 1])]),
        },
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "connectivity_invariants": [10, 20],
            "atom_balls": {"1": {"0": [0, 1], "1": [0, 1]}},
            "full": {"sparse_bit_info": {"999": [[0, 1], [1, 1]]}},  # same class both atoms
        },
    )

    common_r1 = trace([(0, 1, 30, 3, [0, 1]), (1, 1, 30, 3, [0, 1])])

    cases["radius2_environment_membership_mismatch"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": common_r0 + common_r1 + trace([(0, 2, 50, 5, [0, 1]), (1, 2, 50, 5, [0, 1])]),
        },
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "connectivity_invariants": [10, 20],
            "atom_balls": {
                "1": {"0": [0, 1], "1": [0, 1]},
                "2": {"0": [0], "1": [0, 1]},  # atom 0's radius-2 ball disagrees
            },
            "full": {"sparse_bit_info": {"999": [[0, 1], [1, 1]]}},
        },
    )

    cases["radius2_identifier_mismatch"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": common_r0 + common_r1 + trace([(0, 2, 50, 5, [0, 1]), (1, 2, 60, 6, [0, 1])]),
        },
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "connectivity_invariants": [10, 20],
            "atom_balls": {"1": {"0": [0, 1], "1": [0, 1]}, "2": {"0": [0, 1], "1": [0, 1]}},
            "full": {
                "sparse_bit_info": {
                    "999": [[0, 1], [1, 1]],
                    "888": [[0, 2], [1, 2]],  # same class both atoms at radius 2
                }
            },
        },
    )

    common_r2 = trace([(0, 2, 50, 5, [0, 1]), (1, 2, 50, 5, [0, 1])])
    full_no_suppression = {
        "sparse_bit_info": {
            "999": [[0, 1], [1, 1]],
            "888": [[0, 2], [1, 2]],
        }
    }
    balls_ok = {"1": {"0": [0, 1], "1": [0, 1]}, "2": {"0": [0, 1], "1": [0, 1]}}

    cases["redundant_environment_mismatch"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": common_r0 + common_r1 + common_r2,
        },
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "connectivity_invariants": [10, 20],
            "atom_balls": balls_ok,
            "full": full_no_suppression,
            "default": {"sparse_bit_info": {"999": [[0, 1], [1, 1]]}},  # radius-2 pair suppressed
        },
    )

    cases["sparse_count_mismatch"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": common_r0 + common_r1 + common_r2,
        },
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "connectivity_invariants": [10, 20],
            "atom_balls": balls_ok,
            "full": full_no_suppression,
            "default": {
                "sparse_bit_info": full_no_suppression["sparse_bit_info"],
                "sparse_counts": {"999": 2, "888": 1, "777": 1},  # different shape than chematic's
            },
        },
    )

    # chematic's real shape from common_r0+common_r1+common_r2: raw-id counts
    # {10:1, 20:1, 30:2, 50:2} -> sorted shape [1,1,2,2]; folded bits {1,2,3,5}
    # (radius0 atoms 0/1 -> bits 1/2, radius1 both -> bit 3, radius2 both -> bit 5).
    matching_sparse_counts = {"111": 1, "222": 1, "999": 2, "888": 2}

    cases["folding_mismatch"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": common_r0 + common_r1 + common_r2,
        },
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "connectivity_invariants": [10, 20],
            "atom_balls": balls_ok,
            "full": full_no_suppression,
            "default": {
                "sparse_bit_info": full_no_suppression["sparse_bit_info"],
                "sparse_counts": matching_sparse_counts,
                "folded_on_bits": [1, 2, 3],  # missing bit 5 -- differs from chematic's {1,2,3,5}
            },
        },
    )

    cases["bit_info_mismatch"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": common_r0 + common_r1 + common_r2,
        },
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "connectivity_invariants": [10, 20],
            "atom_balls": balls_ok,
            "full": full_no_suppression,
            "default": {
                "sparse_bit_info": full_no_suppression["sparse_bit_info"],
                "sparse_counts": matching_sparse_counts,
                "folded_on_bits": [1, 2, 3, 5],  # same bit SET as chematic's...
                "folded_bit_info": {
                    "1": [[0, 0]],
                    "2": [[1, 0]],
                    "3": [[0, 1]],  # ...but missing atom 1's radius-1 entry: wrong shape
                    "5": [[0, 2], [1, 2]],
                },
            },
        },
    )

    cases["exact_match"] = (
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "trace": common_r0 + common_r1 + common_r2,
        },
        {
            "smiles": "X",
            "parse_ok": True,
            "atom_count": 2,
            "atomic_numbers": [6, 6],
            "connectivity_invariants": [10, 20],
            "atom_balls": balls_ok,
            "full": full_no_suppression,
            "default": {
                "sparse_bit_info": full_no_suppression["sparse_bit_info"],
                "sparse_counts": matching_sparse_counts,
                "folded_on_bits": [1, 2, 3, 5],
                "folded_bit_info": {
                    "1": [[0, 0]],
                    "2": [[1, 0]],
                    "3": [[0, 1], [1, 1]],
                    "5": [[0, 2], [1, 2]],
                },
            },
        },
    )

    return cases


def run_self_test():
    cases = _self_test_cases()
    expected_labels = set(STAGES) | {"exact_match"}
    assert set(cases.keys()) == expected_labels, (
        f"self-test must cover every classification label: "
        f"missing={expected_labels - set(cases.keys())} extra={set(cases.keys()) - expected_labels}"
    )

    failed = []
    for expected, (chem, rd) in cases.items():
        flags, _ = evaluate_stages(chem, rd)
        got = classify(flags)
        status = "OK" if got == expected else "FAIL"
        if got != expected:
            failed.append((expected, got))
        print(f"  [{status}] expected={expected} got={got}")

    if failed:
        print(f"\nSELF-TEST FAILED: {len(failed)} case(s) misclassified: {failed}", file=sys.stderr)
        sys.exit(1)
    print(f"\nSELF-TEST OK: all {len(cases)} classification labels reachable")


# --------------------------------------------------------------------------
# Main comparison run
# --------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--chem-trace")
    p.add_argument("--rdkit-oracle")
    p.add_argument("--summary-out")
    p.add_argument("--rows-out")
    p.add_argument("--first-divergence-out")
    p.add_argument("--rows-scope", choices=["fixtures", "all"], default="fixtures")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()

    if args.self_test:
        run_self_test()
        return

    if not args.chem_trace or not args.rdkit_oracle:
        print("--chem-trace and --rdkit-oracle are required (or pass --self-test)", file=sys.stderr)
        sys.exit(1)

    chem_rows = load_jsonl(args.chem_trace)
    rd_rows = load_jsonl(args.rdkit_oracle)

    if len(chem_rows) != len(rd_rows):
        print(
            f"PIPELINE ERROR: chematic trace has {len(chem_rows)} rows, "
            f"RDKit oracle has {len(rd_rows)} -- inputs were not built from the same "
            f"ordered SMILES list (see module docstring)",
            file=sys.stderr,
        )
        sys.exit(1)

    counts = defaultdict(int)
    per_radius_agree = {0: 0, 1: 0, 2: 0}
    per_radius_total = {0: 0, 1: 0, 2: 0}
    exact_match_denoms = {"sparse_count": 0, "folded_2048_bit": 0, "bit_info_map": 0}
    exact_match_hits = {"sparse_count": 0, "folded_2048_bit": 0, "bit_info_map": 0}
    all_rows_out = []
    tsv_rows = []

    for idx, (chem, rd) in enumerate(zip(chem_rows, rd_rows)):
        if chem.get("smiles") != rd.get("smiles"):
            print(
                f"PIPELINE ERROR at row {idx}: chematic smiles={chem.get('smiles')!r} != "
                f"rdkit smiles={rd.get('smiles')!r} -- inputs out of sync",
                file=sys.stderr,
            )
            sys.exit(1)

        flags, evidence = evaluate_stages(chem, rd)
        label = classify(flags)
        counts[label] += 1

        # Per-radius agreement: "agrees" means neither membership nor
        # identifier mismatch fired at/under that radius (i.e. we got past
        # it during evaluation -- None means blocked earlier, also excluded).
        for r in (0, 1, 2):
            if r == 0:
                blocked = flags["chematic_parse_fail"] or flags["rdkit_parse_fail"] or flags["atom_count_mismatch"]
                if blocked:
                    continue
                per_radius_total[0] += 1
                if not flags["initial_invariant_mismatch"]:
                    per_radius_agree[0] += 1
            else:
                mem_key = f"radius{r}_environment_membership_mismatch"
                id_key = f"radius{r}_identifier_mismatch"
                if flags[mem_key] is None:
                    continue
                per_radius_total[r] += 1
                if not flags[mem_key] and not flags[id_key]:
                    per_radius_agree[r] += 1

        if flags["sparse_count_mismatch"] is not None:
            exact_match_denoms["sparse_count"] += 1
            if not flags["sparse_count_mismatch"]:
                exact_match_hits["sparse_count"] += 1
        if flags["folding_mismatch"] is not None:
            exact_match_denoms["folded_2048_bit"] += 1
            if not flags["folding_mismatch"]:
                exact_match_hits["folded_2048_bit"] += 1
        if flags["bit_info_mismatch"] is not None:
            exact_match_denoms["bit_info_map"] += 1
            if not flags["bit_info_mismatch"]:
                exact_match_hits["bit_info_map"] += 1

        source = rd.get("source", "unknown")
        div_atom = evidence.get(label, {}).get("atom_idx")
        div_radius = evidence.get(label, {}).get("radius")
        tsv_rows.append(
            [
                chem.get("smiles", ""),
                source,
                label,
                "" if div_atom is None else str(div_atom),
                "" if div_radius is None else str(div_radius),
                str(chem.get("atom_count", "")),
                str(rd.get("atom_count", "")),
            ]
        )

        if args.rows_scope == "all" or source.startswith("fixture:"):
            all_rows_out.append(
                {
                    "smiles": chem.get("smiles", ""),
                    "source": source,
                    "first_divergence": label,
                    "flags": flags,
                    "evidence": evidence,
                }
            )

    total = len(chem_rows)
    bucket_sum = sum(counts.values())

    if args.rows_out:
        with open(args.rows_out, "w") as f:
            for row in all_rows_out:
                f.write(json.dumps(row, sort_keys=True) + "\n")

    if args.first_divergence_out:
        with open(args.first_divergence_out, "w") as f:
            f.write(
                "smiles\tsource\tfirst_divergence\tdivergent_atom_idx\tdivergent_radius\t"
                "chematic_atom_count\trdkit_atom_count\n"
            )
            for row in tsv_rows:
                f.write("\t".join(row) + "\n")

    def pct(hits, denom):
        return round(100 * hits / denom, 4) if denom else None

    summary = {
        "schema_version": "1",
        "input_counts": {"total": total},
        "first_divergence_counts": {s: counts.get(s, 0) for s in STAGES + ["exact_match"]},
        "bucket_sum_check": {"sum": bucket_sum, "input_total": total, "ok": bucket_sum == total},
        "duplicate_input_ids": 0,  # rows are matched by position, never by a dedupable key
        "per_radius_agreement_pct": {
            f"radius{r}": pct(per_radius_agree[r], per_radius_total[r]) for r in (0, 1, 2)
        },
        "exact_match_rates": {
            k: pct(exact_match_hits[k], exact_match_denoms[k]) for k in exact_match_hits
        },
        "tanimoto_correlation_vs_rdkit": tanimoto_correlation(chem_rows, rd_rows),
        "dominant_mechanism": max(counts, key=counts.get) if counts else None,
    }

    print(json.dumps(summary, indent=2, sort_keys=True))

    if args.summary_out:
        with open(args.summary_out, "w") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
            f.write("\n")

    if not summary["bucket_sum_check"]["ok"]:
        print("GATE FAILED: bucket sum != input total", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
