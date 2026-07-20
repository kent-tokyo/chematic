#!/usr/bin/env python3
"""Phase B acceptance gate: `rdkit_morgan_ecfp4_experimental` (the production,
fallible, RDKit-bit-exact ECFP4 API in `crates/chematic-fp/src/rdkit_morgan_ecfp4.rs`)
vs the real RDKit oracle (`scripts/gen_ecfp_rdkit_environment_oracle.py`).

Unlike `ecfp_rdkit_raw_identifier_parity.py`'s `classify()`, this script never
touches RDKit's `includeRedundantEnvironments=True` ("full") lifecycle --
the production API only ever computes the suppressed ("default") lifecycle,
which is its entire claim surface. Reuses the oracle-side helpers
(`rd_default_pairs`/`rd_sparse_counts`/`rd_folded_bit_info`) directly from
`ecfp_rdkit_raw_identifier_parity.py` so the oracle-reading logic can never
silently drift between the two scripts.

Denominator discipline (same as `ecfp_rdkit_raw_identifier_parity_aromaticity_variant.py`):
the exact-match rate is computed ONLY over rows where the chematic dump's own
`status == "success"`. Error rows are reported separately, never pooled in.

Usage:
    python scripts/ecfp_rdkit_morgan_ecfp4_parity.py \\
        --chematic <rdkit_morgan_ecfp4_dump.jsonl> \\
        --rdkit-oracle <gen_ecfp_rdkit_environment_oracle.py --rows-out output> \\
        --summary-out <out.json> [--mismatches-out <out.jsonl>]

    python scripts/ecfp_rdkit_morgan_ecfp4_parity.py --self-test
"""

from __future__ import annotations

import argparse
import json
import sys

from ecfp_rdkit_raw_identifier_parity import rd_default_pairs, rd_folded_bit_info, rd_sparse_counts

NBITS = 2048


def load_jsonl(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def chem_default_pairs(row):
    return {(a, r): rid for a, r, rid in row["default_pairs"]}


def chem_sparse_counts(row):
    return {rid: count for rid, count in row["sparse_counts"]}


def chem_folded_bit_info(row):
    return {bit: {(a, r) for a, r in envs} for bit, envs in row["folded_bit_info"]}


def is_exact_match(chem, rd):
    """True only if every layer this production API claims (default-lifecycle
    raw pairs, sparse counts, folded on-bits, folded bitInfo) matches. Only
    valid for `status == "success"` rows."""
    if chem_default_pairs(chem) != rd_default_pairs(rd):
        return False
    if chem_sparse_counts(chem) != rd_sparse_counts(rd):
        return False
    if set(chem["folded_on_bits"]) != set(rd["default"]["folded_on_bits"]):
        return False
    if chem_folded_bit_info(chem) != rd_folded_bit_info(rd):
        return False
    return True


def run(chematic_rows, rdkit_rows):
    if len(chematic_rows) != len(rdkit_rows):
        print(
            f"PIPELINE ERROR: row count mismatch chematic={len(chematic_rows)} rdkit={len(rdkit_rows)}",
            file=sys.stderr,
        )
        sys.exit(1)

    total_inputs = len(chematic_rows)
    success_rows = []
    error_rows = []

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

        if chem.get("status") == "success":
            success_rows.append((idx, chem, rd))
        else:
            error_rows.append({"row_id": idx, "smiles": chem.get("smiles"), "status": chem.get("status")})

    exact = 0
    non_exact_row_ids = []
    for idx, chem, rd in success_rows:
        if is_exact_match(chem, rd):
            exact += 1
        else:
            non_exact_row_ids.append(idx)

    success_count = len(success_rows)
    exact_pct = round(exact / success_count, 4) if success_count else None

    gate_failures = []
    if non_exact_row_ids:
        gate_failures.append(f"{len(non_exact_row_ids)} success row(s) not exact_match: {non_exact_row_ids[:10]}")

    summary = {
        "total_inputs": total_inputs,
        "success_count": success_count,
        "error_count": len(error_rows),
        "exact_match_among_success": exact,
        "exact_match_pct_among_success": exact_pct,
        "non_exact_success_row_ids": non_exact_row_ids,
        "error_rows": error_rows,
        "gate_failures": gate_failures,
        "gate_passed": len(gate_failures) == 0,
    }

    if gate_failures:
        print(f"GATE FAILED: {len(gate_failures)} violation(s):", file=sys.stderr)
        for msg in gate_failures:
            print(f"  - {msg}", file=sys.stderr)

    return summary, summary["gate_passed"]


def run_self_test():
    checks = []

    def row(default_pairs, sparse_counts=None, folded_bit_info=None, status="success"):
        sparse_counts = sparse_counts if sparse_counts is not None else {}
        folded_bit_info = folded_bit_info if folded_bit_info is not None else {}
        for (_a, _r), rid in {(a, r): rid for a, r, rid in default_pairs}.items():
            sparse_counts.setdefault(rid, 0)
            sparse_counts[rid] += 1
        chem = {
            "row_id": 0,
            "smiles": "X",
            "status": status,
            "default_pairs": default_pairs,
            "sparse_counts": list(sparse_counts.items()),
            "folded_on_bits": sorted({rid % NBITS for _a, _r, rid in default_pairs}),
            "folded_bit_info": [
                (bit, sorted(envs)) for bit, envs in (folded_bit_info or _derive_folded_info(default_pairs)).items()
            ],
        }
        rd = {
            "row_id": 0,
            "smiles": "X",
            "default": {
                "sparse_bit_info": {str(rid): [[a, r] for a, r, rid2 in default_pairs if rid2 == rid] for _a, _r, rid in default_pairs},
                "sparse_counts": {str(k): v for k, v in sparse_counts.items()},
                "folded_on_bits": sorted({rid % NBITS for _a, _r, rid in default_pairs}),
                "folded_bit_info": {
                    str(bit): [list(e) for e in envs]
                    for bit, envs in (folded_bit_info or _derive_folded_info(default_pairs)).items()
                },
            },
        }
        return chem, rd

    def _derive_folded_info(default_pairs):
        info = {}
        for a, r, rid in default_pairs:
            info.setdefault(rid % NBITS, set()).add((a, r))
        return info

    chem, rd = row([(0, 0, 100), (1, 0, 100)])
    checks.append(("exact_match", is_exact_match(chem, rd) is True))

    chem, rd = row([(0, 0, 100), (1, 0, 100)])
    chem["default_pairs"] = [(0, 0, 999), (1, 0, 100)]
    checks.append(("raw_id_mismatch", is_exact_match(chem, rd) is False))

    chem, rd = row([(0, 0, 100), (1, 0, 100)])
    chem["sparse_counts"] = [(100, 5)]
    checks.append(("sparse_count_mismatch", is_exact_match(chem, rd) is False))

    chem, rd = row([(0, 0, 100), (1, 0, 100)])
    chem["folded_on_bits"] = [999999]
    checks.append(("folded_bit_mismatch", is_exact_match(chem, rd) is False))

    chem, rd = row([(0, 0, 100), (1, 0, 100)])
    chem["folded_bit_info"] = [(100, [[0, 0]])]
    checks.append(("bitinfo_mismatch", is_exact_match(chem, rd) is False))

    ok = True
    for name, passed in checks:
        status = "OK" if passed else "FAIL"
        print(f"  self-test {name}: {status}")
        ok = ok and passed

    # live gate control: an error row must be excluded from the denominator,
    # not counted as a mismatch and not counted as a match.
    chem_success, rd_success = row([(0, 0, 100)])
    chem_error = {"row_id": 1, "smiles": "Y", "status": "rdkit_parity_kekulization_failed"}
    rd_error = {"row_id": 1, "smiles": "Y", "default": {"sparse_bit_info": {}, "sparse_counts": {}, "folded_on_bits": [], "folded_bit_info": {}}}
    summary, gate_passed = run([chem_success, chem_error], [rd_success, rd_error])
    live_control_ok = (
        summary["success_count"] == 1
        and summary["error_count"] == 1
        and summary["exact_match_among_success"] == 1
        and gate_passed is True
    )
    print(f"  live gate control (error row excluded from denominator): {'OK' if live_control_ok else 'FAIL'}")
    ok = ok and live_control_ok

    # live gate control: a genuine mismatch in a success row must fail the gate.
    chem_bad, rd_bad = row([(0, 0, 100), (1, 0, 100)])
    chem_bad["default_pairs"] = [(0, 0, 999), (1, 0, 100)]
    _summary2, gate_passed2 = run([chem_bad], [rd_bad])
    mismatch_gate_ok = gate_passed2 is False
    print(f"  live gate control (success-row mismatch -> gate fails): {'OK' if mismatch_gate_ok else 'FAIL'}")
    ok = ok and mismatch_gate_ok

    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chematic")
    p.add_argument("--rdkit-oracle")
    p.add_argument("--summary-out", default=None)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()

    if args.self_test:
        ok = run_self_test()
        sys.exit(0 if ok else 1)

    if not args.chematic or not args.rdkit_oracle:
        p.error("--chematic and --rdkit-oracle are required unless --self-test")

    chematic_rows = load_jsonl(args.chematic)
    rdkit_rows = load_jsonl(args.rdkit_oracle)
    summary, gate_passed = run(chematic_rows, rdkit_rows)

    print(json.dumps(summary, indent=2))
    if args.summary_out:
        with open(args.summary_out, "w") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
        print(f"summary written to {args.summary_out}")

    sys.exit(0 if gate_passed else 1)


if __name__ == "__main__":
    main()
