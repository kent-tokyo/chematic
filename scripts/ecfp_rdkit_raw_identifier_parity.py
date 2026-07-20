#!/usr/bin/env python3
"""Morgan M4-A0: numeric raw-identifier parity between chematic's RDKit-exact
hash port (`crates/chematic-fp/src/rdkit_morgan_hash.rs`,
`examples/rdkit_morgan_hash_dump.rs`) and a real RDKit oracle
(`scripts/gen_ecfp_rdkit_environment_oracle.py`).

Unlike PR #120/#123's scripts, which compare *partitions* (which atoms are
chemically equivalent) and *emission lifecycle* (who wins/dies), this script
compares actual 32-bit hash VALUES -- the thing PR #120/#123 explicitly
never claimed. See `rdkit_morgan_hash.rs`'s module doc for the ported
formulas and the pinned RDKit commit.

Per molecule, classifies the FIRST layer at which chematic's port and RDKit
diverge, in this fixed order (each layer assumes every earlier layer already
matched):

    parse_mismatch            -- one side parses, the other doesn't
    both_parse_fail           -- both fail (not a divergence signal, still
                                  accounted for separately from the other
                                  buckets)
    radius0_numeric_mismatch  -- some atom's radius-0 identifier differs
    radius1_numeric_mismatch  -- radius 0 matches; some radius-1 (atom,id)
                                  pair differs or is present on only one side
    radius2_numeric_mismatch  -- radius 0-1 match; likewise for radius 2
    representative_selection_mismatch -- every raw id at every radius
                                  matches, but *which* atom emits under
                                  suppression (RDKit's "default" lifecycle)
                                  differs
    sparse_count_mismatch     -- representative selection matches, but the
                                  raw-id count multiset differs
    folded_bit_mismatch       -- sparse counts match, folded 2048-bit
                                  on-bit set differs
    bitinfo_mismatch          -- folded bits match, bit->(atom,radius)
                                  attribution differs
    exact_match               -- every layer matches

Usage:
    python scripts/ecfp_rdkit_raw_identifier_parity.py \\
        --chematic <rdkit_morgan_hash_dump.jsonl> \\
        --rdkit-oracle <gen_ecfp_rdkit_environment_oracle.py --rows-out output> \\
        --summary-out <out.json> [--mismatches-out <out.jsonl>]

    python scripts/ecfp_rdkit_raw_identifier_parity.py --self-test
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter

BUCKET_ORDER = [
    "parse_mismatch",
    "both_parse_fail",
    "radius0_numeric_mismatch",
    "radius1_numeric_mismatch",
    "radius2_numeric_mismatch",
    "representative_selection_mismatch",
    "sparse_count_mismatch",
    "folded_bit_mismatch",
    "bitinfo_mismatch",
    "exact_match",
]

NBITS = 2048

# Pinned ceiling on the 5,048-input M4-A0 corpus (`combined_input.csv` =
# 5,000-molecule corpus + PR #120's 41 fixtures + PR #123's 4 fixtures +
# `ecfp_rdkit_m4a0_hash_fixtures.csv`'s 3), measured 2026-07-20 with
# production Hueckel aromaticity applied before hashing (the main
# `rdkit_morgan_hash_dump` path). radius0 is pinned hard at 0 -- fully
# achieved, any regression is a real hash-port defect. radius1 accepts up to
# the documented residual (all 59 traced to Hueckel-vs-RDKit aromaticity
# *perception* disagreement on fused/macrocyclic rings, not a hash defect --
# resolves to 0/5048 under `apply_aromaticity_rdkit_parity_experimental`,
# see `rdkit_morgan_hash_dump_aromaticity_variant`). Every other bucket is
# pinned at 0: already fully achieved, no known residual to accept. This is
# a ceiling, not a target -- a run with FEWER mismatches than the ceiling
# still passes (matching PR #123's "reject regressions, accept
# improvements" gate philosophy); only exceeding it fails.
MAX_ACCEPTED_MISMATCHES = {
    "parse_mismatch": 0,
    "radius0_numeric_mismatch": 0,
    "radius1_numeric_mismatch": 59,
    "radius2_numeric_mismatch": 0,
    "representative_selection_mismatch": 0,
    "sparse_count_mismatch": 0,
    "folded_bit_mismatch": 0,
    "bitinfo_mismatch": 0,
}


def load_jsonl(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def chem_by_radius_full(row):
    """(atom, radius) -> raw_id under RDKit's `includeRedundantEnvironments=
    true` ("full") lifecycle. The two identifiers are tracked independently
    in `entries` (see `RdkitMorganRawTraceEntry`'s doc comment) because they
    can genuinely differ -- do not assume `raw_identifier_full ==
    raw_identifier_default` when both are present."""
    return {(a, r): rid_full for a, r, rid_full, _rid_default in row["entries"] if rid_full is not None}


def chem_default_pairs(row):
    """(atom, radius) -> raw_id, only entries chematic's port emits under
    the suppressed (RDKit `default`) lifecycle."""
    return {
        (a, r): rid_default for a, r, _rid_full, rid_default in row["entries"] if rid_default is not None
    }


def rd_full_by_radius(row):
    inv = {}
    for raw_id_str, entries in row["full"]["sparse_bit_info"].items():
        raw_id = int(raw_id_str) & 0xFFFFFFFF
        for atom, radius in entries:
            inv[(atom, radius)] = raw_id
    return inv


def rd_default_pairs(row):
    inv = {}
    for raw_id_str, entries in row["default"]["sparse_bit_info"].items():
        raw_id = int(raw_id_str) & 0xFFFFFFFF
        for atom, radius in entries:
            inv[(atom, radius)] = raw_id
    return inv


def rd_sparse_counts(row):
    return {int(k) & 0xFFFFFFFF: v for k, v in row["default"]["sparse_counts"].items()}


def chem_sparse_counts(row):
    return Counter(rid for (_a, _r), rid in chem_default_pairs(row).items())


def chem_folded_bits(row):
    return {rid % NBITS for (_a, _r), rid in chem_default_pairs(row).items()}


def chem_folded_bit_info(row):
    info = {}
    for (a, r), rid in chem_default_pairs(row).items():
        info.setdefault(rid % NBITS, set()).add((a, r))
    return info


def rd_folded_bit_info(row):
    info = {}
    for bit_str, entries in row["default"]["folded_bit_info"].items():
        info[int(bit_str)] = {(a, r) for a, r in entries}
    return info


def _radius_layer_mismatch(chem_full, rd_full, radius, atom_count):
    """True if any (atom, radius) pair at this exact radius differs in
    presence or value between the two sides, for atoms 0..atom_count."""
    for atom in range(atom_count):
        key = (atom, radius)
        chem_has = key in chem_full
        rd_has = key in rd_full
        if chem_has != rd_has:
            return True
        if chem_has and chem_full[key] != rd_full[key]:
            return True
    return False


def classify(chem, rd):
    chem_ok = chem.get("parse_ok", False)
    rd_ok = rd.get("parse_ok", False)
    if chem_ok != rd_ok:
        return "parse_mismatch"
    if not chem_ok and not rd_ok:
        return "both_parse_fail"

    atom_count = chem["atom_count"]
    chem_full = chem_by_radius_full(chem)
    rd_full = rd_full_by_radius(rd)

    if _radius_layer_mismatch(chem_full, rd_full, 0, atom_count):
        return "radius0_numeric_mismatch"
    if _radius_layer_mismatch(chem_full, rd_full, 1, atom_count):
        return "radius1_numeric_mismatch"
    if _radius_layer_mismatch(chem_full, rd_full, 2, atom_count):
        return "radius2_numeric_mismatch"

    chem_def = chem_default_pairs(chem)
    rd_def = rd_default_pairs(rd)
    if set(chem_def.keys()) != set(rd_def.keys()):
        return "representative_selection_mismatch"
    for key in chem_def:
        if chem_def[key] != rd_def[key]:
            return "representative_selection_mismatch"

    if sorted(chem_sparse_counts(chem).values()) != sorted(rd_sparse_counts(rd).values()):
        return "sparse_count_mismatch"
    if dict(chem_sparse_counts(chem)) != rd_sparse_counts(rd):
        return "sparse_count_mismatch"

    if chem_folded_bits(chem) != set(rd["default"]["folded_on_bits"]):
        return "folded_bit_mismatch"

    if chem_folded_bit_info(chem) != rd_folded_bit_info(rd):
        return "bitinfo_mismatch"

    return "exact_match"


def run(chematic_rows, rdkit_rows):
    gate_failures = []

    if len(chematic_rows) != len(rdkit_rows):
        print(
            f"PIPELINE ERROR: row count mismatch chematic={len(chematic_rows)} "
            f"rdkit={len(rdkit_rows)}",
            file=sys.stderr,
        )
        sys.exit(1)

    # A duplicate/missing/reordered row_id is caught here unconditionally,
    # as a hard pipeline failure (exit 1) rather than a soft gate_failures
    # entry -- matching ecfp_rdkit_suppression_parity.py's established
    # convention: row identity desync is a *tooling* bug (the two dumps got
    # out of sync), not a measurement finding, and self-test's
    # `duplicate_row_id_exits_nonzero` exercises this path directly.
    counts = {b: 0 for b in BUCKET_ORDER}
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

        bucket = classify(chem, rd)
        counts[bucket] += 1
        if bucket not in ("exact_match", "both_parse_fail"):
            mismatches.append({"row_id": idx, "smiles": chem.get("smiles"), "bucket": bucket})

    total = len(chematic_rows)
    accounted = sum(counts.values())
    if accounted != total:
        gate_failures.append(f"input accounting incomplete: {accounted} classified of {total} total rows")

    if counts["parse_mismatch"] > 0:
        gate_failures.append(f"{counts['parse_mismatch']} parse_mismatch row(s) (asymmetric parse failure)")

    for bucket, ceiling in MAX_ACCEPTED_MISMATCHES.items():
        if bucket == "parse_mismatch":
            continue  # already gated above with its own message
        if counts[bucket] > ceiling:
            gate_failures.append(
                f"{bucket}: {counts[bucket]} exceeds the pinned ceiling of {ceiling} "
                f"(regression against the documented M4-A0 baseline)"
            )

    summary = {
        "total_rows": total,
        "accounted_rows": accounted,
        "bucket_counts": counts,
        "radius0_exact_match_rate": round(
            (total - counts["radius0_numeric_mismatch"] - counts["parse_mismatch"] - counts["both_parse_fail"])
            / max(1, total - counts["parse_mismatch"] - counts["both_parse_fail"]),
            4,
        ),
        "full_exact_match_rate": round(
            counts["exact_match"] / max(1, total - counts["parse_mismatch"] - counts["both_parse_fail"]), 4
        ),
        "gate_failures": gate_failures,
        "gate_passed": len(gate_failures) == 0,
    }

    if gate_failures:
        print(f"GATE FAILED: {len(gate_failures)} violation(s):", file=sys.stderr)
        for msg in gate_failures:
            print(f"  - {msg}", file=sys.stderr)

    return summary, mismatches, summary["gate_passed"]


def _self_test_row(entries, rd_full_bit_info, rd_default_bit_info, rd_counts=None, rd_folded_bits=None, rd_folded_info=None, atom_count=2, parse_ok=True):
    rd_counts = rd_counts if rd_counts is not None else {}
    for raw_id_str, pairs in rd_default_bit_info.items():
        rd_counts.setdefault(raw_id_str, len(pairs))
    if rd_folded_bits is None:
        rd_folded_bits = sorted({int(k) % NBITS for k in rd_default_bit_info})
    if rd_folded_info is None:
        rd_folded_info = {}
        for raw_id_str, pairs in rd_default_bit_info.items():
            bit = str(int(raw_id_str) % NBITS)
            rd_folded_info.setdefault(bit, []).extend(pairs)
    chem = {"row_id": 0, "smiles": "X", "parse_ok": parse_ok, "atom_count": atom_count, "entries": entries}
    rd = {
        "row_id": 0,
        "smiles": "X",
        "parse_ok": parse_ok,
        "full": {"sparse_bit_info": rd_full_bit_info},
        "default": {
            "sparse_bit_info": rd_default_bit_info,
            "sparse_counts": rd_counts,
            "folded_on_bits": rd_folded_bits,
            "folded_bit_info": rd_folded_info,
        },
    }
    return chem, rd


def run_self_test():
    checks = []

    # exact_match: everything lines up.
    chem, rd = _self_test_row(
        entries=[[0, 0, 100, 100], [1, 0, 100, 100]],
        rd_full_bit_info={"100": [[0, 0], [1, 0]]},
        rd_default_bit_info={"100": [[0, 0], [1, 0]]},
    )
    checks.append(("exact_match", classify(chem, rd) == "exact_match"))

    # radius0_numeric_mismatch: chematic's atom 0 radius-0 id differs.
    chem, rd = _self_test_row(
        entries=[[0, 0, 999, 999], [1, 0, 100, 100]],
        rd_full_bit_info={"100": [[0, 0], [1, 0]]},
        rd_default_bit_info={"100": [[0, 0], [1, 0]]},
    )
    checks.append(("radius0_numeric_mismatch", classify(chem, rd) == "radius0_numeric_mismatch"))

    # radius1_numeric_mismatch: radius 0 matches, radius 1 value differs.
    chem, rd = _self_test_row(
        entries=[
            [0, 0, 100, 100],
            [1, 0, 100, 100],
            [0, 1, 555, 555],
        ],
        rd_full_bit_info={"100": [[0, 0], [1, 0]], "200": [[0, 1]]},
        rd_default_bit_info={"100": [[0, 0], [1, 0]], "200": [[0, 1]]},
        atom_count=2,
    )
    checks.append(("radius1_numeric_mismatch", classify(chem, rd) == "radius1_numeric_mismatch"))

    # representative_selection_mismatch: raw ids all match, but chematic
    # picks atom 1 to emit at radius 1 where RDKit picks atom 0.
    chem, rd = _self_test_row(
        entries=[
            [0, 0, 100, 100],
            [1, 0, 100, 100],
            [0, 1, 200, None],
            [1, 1, 200, 200],
        ],
        rd_full_bit_info={"100": [[0, 0], [1, 0]], "200": [[0, 1], [1, 1]]},
        rd_default_bit_info={"100": [[0, 0], [1, 0]], "200": [[0, 1]]},
    )
    checks.append(
        ("representative_selection_mismatch", classify(chem, rd) == "representative_selection_mismatch")
    )

    # sparse_count_mismatch: same emitted pairs/ids, but RDKit's oracle
    # reports a different count multiplicity for the shared raw id (can't
    # happen from a consistent oracle in practice, but the classifier must
    # still catch a data-level inconsistency rather than silently passing).
    chem, rd = _self_test_row(
        entries=[[0, 0, 100, 100], [1, 0, 100, 100]],
        rd_full_bit_info={"100": [[0, 0], [1, 0]]},
        rd_default_bit_info={"100": [[0, 0], [1, 0]]},
        rd_counts={"100": 5},
    )
    checks.append(("sparse_count_mismatch", classify(chem, rd) == "sparse_count_mismatch"))

    # folded_bit_mismatch: counts match but RDKit's folded_on_bits omits a
    # bit chematic's own raw-id%2048 computation would produce.
    chem, rd = _self_test_row(
        entries=[[0, 0, 100, 100], [1, 0, 100, 100]],
        rd_full_bit_info={"100": [[0, 0], [1, 0]]},
        rd_default_bit_info={"100": [[0, 0], [1, 0]]},
        rd_folded_bits=[999999],
    )
    checks.append(("folded_bit_mismatch", classify(chem, rd) == "folded_bit_mismatch"))

    # bitinfo_mismatch: folded bits match but attribution set differs.
    chem, rd = _self_test_row(
        entries=[[0, 0, 100, 100], [1, 0, 100, 100]],
        rd_full_bit_info={"100": [[0, 0], [1, 0]]},
        rd_default_bit_info={"100": [[0, 0], [1, 0]]},
        rd_folded_info={"100": [[0, 0]]},
    )
    checks.append(("bitinfo_mismatch", classify(chem, rd) == "bitinfo_mismatch"))

    # parse_mismatch positive control.
    chem = {"row_id": 0, "smiles": "X", "parse_ok": False}
    rd = {"row_id": 0, "smiles": "X", "parse_ok": True, "full": {"sparse_bit_info": {}}, "default": {"sparse_bit_info": {}, "sparse_counts": {}, "folded_on_bits": [], "folded_bit_info": {}}, "atom_count": 0}
    checks.append(("parse_mismatch", classify(chem, rd) == "parse_mismatch"))

    ok = True
    for name, passed in checks:
        status = "OK" if passed else "FAIL"
        print(f"  self-test {name}: {status}")
        ok = ok and passed

    # run()-level gate positive control 1: an asymmetric parse failure
    # (chematic fails, RDKit doesn't) must make gate_passed False.
    chem_rows = [{"row_id": 0, "smiles": "A", "parse_ok": False, "atom_count": 0, "entries": []}]
    rd_rows = [
        {
            "row_id": 0,
            "smiles": "A",
            "parse_ok": True,
            "full": {"sparse_bit_info": {}},
            "default": {"sparse_bit_info": {}, "sparse_counts": {}, "folded_on_bits": [], "folded_bit_info": {}},
        }
    ]
    _summary, _mismatches, gate_passed = run(chem_rows, rd_rows)
    live_control_ok = gate_passed is False
    print(f"  live gate control (asymmetric parse fail -> gate fails): {'OK' if live_control_ok else 'FAIL'}")
    ok = ok and live_control_ok

    # run()-level gate positive control 2: a duplicate/desynced row_id must
    # be a hard pipeline failure (SystemExit(1)), not a silently-passing row.
    chem_rows_dup = [
        {"row_id": 0, "smiles": "A", "parse_ok": True, "atom_count": 1, "entries": [[0, 0, 1, 1]]},
        {"row_id": 0, "smiles": "B", "parse_ok": True, "atom_count": 1, "entries": [[0, 0, 1, 1]]},
    ]
    rd_rows_dup = [
        {
            "row_id": 0,
            "smiles": "A",
            "parse_ok": True,
            "full": {"sparse_bit_info": {"1": [[0, 0]]}},
            "default": {
                "sparse_bit_info": {"1": [[0, 0]]},
                "sparse_counts": {"1": 1},
                "folded_on_bits": [1],
                "folded_bit_info": {"1": [[0, 0]]},
            },
        },
        {
            "row_id": 1,
            "smiles": "B",
            "parse_ok": True,
            "full": {"sparse_bit_info": {"1": [[0, 0]]}},
            "default": {
                "sparse_bit_info": {"1": [[0, 0]]},
                "sparse_counts": {"1": 1},
                "folded_on_bits": [1],
                "folded_bit_info": {"1": [[0, 0]]},
            },
        },
    ]
    try:
        run(chem_rows_dup, rd_rows_dup)
        duplicate_control_ok = False
    except SystemExit as e:
        duplicate_control_ok = e.code != 0
    print(f"  live gate control (duplicate row_id -> SystemExit(1)): {'OK' if duplicate_control_ok else 'FAIL'}")
    ok = ok and duplicate_control_ok

    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chematic")
    p.add_argument("--rdkit-oracle")
    p.add_argument("--summary-out", default=None)
    p.add_argument("--mismatches-out", default=None)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()

    if args.self_test:
        ok = run_self_test()
        sys.exit(0 if ok else 1)

    if not args.chematic or not args.rdkit_oracle:
        p.error("--chematic and --rdkit-oracle are required unless --self-test")

    chematic_rows = load_jsonl(args.chematic)
    rdkit_rows = load_jsonl(args.rdkit_oracle)
    summary, mismatches, gate_passed = run(chematic_rows, rdkit_rows)

    print(json.dumps(summary, indent=2))
    if args.summary_out:
        with open(args.summary_out, "w") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
        print(f"summary written to {args.summary_out}")
    if args.mismatches_out:
        with open(args.mismatches_out, "w") as f:
            for m in mismatches:
                f.write(json.dumps(m, sort_keys=True) + "\n")
        print(f"{len(mismatches)} mismatches written to {args.mismatches_out}")

    sys.exit(0 if gate_passed else 1)


if __name__ == "__main__":
    main()
