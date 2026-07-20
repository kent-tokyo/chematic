#!/usr/bin/env python3
"""Morgan M4-A0 residual-mechanism confirmation, measured honestly.

Compares `rdkit_morgan_hash_dump_aromaticity_variant`'s output (RDKit-parity
aromaticity, `Result`-fallible, NO Hueckel fallback -- see that example's
doc comment for why an earlier version's silent fallback was wrong) against
the real RDKit oracle, with denominators kept strictly separate:

  - Exact-match rate is computed ONLY over rows where
    `aromaticity_status == "rdkit_parity_success"`.
  - Rows where RDKit-parity aromaticity preprocessing itself failed
    (`rdkit_parity_kekulization_failed` / `rdkit_parity_internal_error`) are
    reported as errors, never silently folded into the exact-match
    denominator, and never silently "resolved" by any other algorithm.
  - A Hueckel-engine control on JUST the error rows is reported separately,
    explicitly marked `"gating": false` -- it answers "does the OLD
    (Hueckel) path at least agree with RDKit on these two inputs", not
    "does the RDKit-parity engine work here" (it explicitly doesn't, that's
    why they're error rows).

Also cross-references every row `ecfp_rdkit_raw_identifier_parity.py`
classified as `radius1_numeric_mismatch` under Hueckel aromaticity, and
reports precisely how many of those actually got a chance to resolve
(RDKit-parity succeeded AND became exact_match) vs how many could not be
evaluated (RDKit-parity itself failed on that input) vs any that remain
mismatching despite RDKit-parity succeeding.

Usage:
    python scripts/ecfp_rdkit_raw_identifier_parity_aromaticity_variant.py \\
        --chematic-variant <rdkit_morgan_hash_dump_aromaticity_variant.jsonl> \\
        --chematic-hueckel <rdkit_morgan_hash_dump.jsonl> \\
        --rdkit-oracle <gen_ecfp_rdkit_environment_oracle.py --rows-out output> \\
        --hueckel-mismatches <ecfp_rdkit_raw_identifier_parity.py --mismatches-out output> \\
        --summary-out <out.json>
"""

from __future__ import annotations

import argparse
import json
import sys

from ecfp_rdkit_raw_identifier_parity import classify


def load_jsonl(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def is_exact_match(chem_row, rd_row):
    """True only if every layer (radius 0-2, representative selection,
    sparse counts, folded bits, bitInfo) matches -- delegates to the SAME
    `classify()` used by `ecfp_rdkit_raw_identifier_parity.py`'s main run,
    so the two scripts can never silently drift on what "exact match"
    means. Only valid for `rdkit_parity_success` rows (a real `entries`
    array); callers must not invoke this on an error row (`entries: null`)."""
    return classify(chem_row, rd_row) == "exact_match"


def run(variant_rows, oracle_rows, hueckel_rows=None, hueckel_mismatches=None):
    if len(variant_rows) != len(oracle_rows):
        print(
            f"PIPELINE ERROR: row count mismatch variant={len(variant_rows)} "
            f"oracle={len(oracle_rows)}",
            file=sys.stderr,
        )
        sys.exit(1)

    total_inputs = len(variant_rows)
    success_rows = []
    error_rows = []

    for idx, (chem, rd) in enumerate(zip(variant_rows, oracle_rows)):
        if chem.get("row_id") != idx or rd.get("row_id") != idx:
            print(f"PIPELINE ERROR at position {idx}: row_id out of sync", file=sys.stderr)
            sys.exit(1)
        if chem.get("smiles") != rd.get("smiles"):
            print(
                f"PIPELINE ERROR at row {idx}: variant smiles={chem.get('smiles')!r} "
                f"!= rdkit smiles={rd.get('smiles')!r}",
                file=sys.stderr,
            )
            sys.exit(1)

        status = chem.get("aromaticity_status")
        if status == "rdkit_parity_success":
            success_rows.append((idx, chem, rd))
        else:
            error_rows.append(
                {
                    "row_id": idx,
                    "smiles": chem.get("smiles"),
                    "aromaticity_status": status,
                    "aromaticity_error": chem.get("aromaticity_error"),
                }
            )

    exact_among_success = 0
    non_exact_success_rows = []
    for idx, chem, rd in success_rows:
        if is_exact_match(chem, rd):
            exact_among_success += 1
        else:
            non_exact_success_rows.append(idx)

    rdkit_parity_success = len(success_rows)
    rdkit_parity_error = len(error_rows)
    exact_match_pct = round(exact_among_success / rdkit_parity_success, 4) if rdkit_parity_success else None

    summary = {
        "total_inputs": total_inputs,
        "rdkit_parity_success": rdkit_parity_success,
        "rdkit_parity_error": rdkit_parity_error,
        "exact_match_among_rdkit_parity_success": exact_among_success,
        "exact_match_pct_among_rdkit_parity_success": exact_match_pct,
        "non_exact_rdkit_parity_success_row_ids": non_exact_success_rows,
        "error_rows": error_rows,
    }

    # Optional, explicitly non-gating Hueckel control restricted to the
    # error rows: does the OLD path at least agree with RDKit here.
    if hueckel_rows is not None and error_rows:
        if len(hueckel_rows) != total_inputs:
            print(
                f"PIPELINE ERROR: hueckel row count {len(hueckel_rows)} != total_inputs {total_inputs}",
                file=sys.stderr,
            )
            sys.exit(1)
        error_row_ids = {e["row_id"] for e in error_rows}
        evaluated = 0
        exact = 0
        for idx in sorted(error_row_ids):
            hchem = hueckel_rows[idx]
            rd = oracle_rows[idx]
            if hchem.get("row_id") != idx:
                print(f"PIPELINE ERROR: hueckel row_id desync at {idx}", file=sys.stderr)
                sys.exit(1)
            if not hchem.get("parse_ok"):
                continue
            evaluated += 1
            if is_exact_match(hchem, rd):
                exact += 1
        summary["hueckel_control_on_parity_errors"] = {
            "evaluated": evaluated,
            "exact_match": exact,
            "gating": False,
        }

    # Cross-reference the 59 Hueckel radius1_numeric_mismatch rows: how many
    # actually got a chance to resolve under RDKit-parity vs couldn't be
    # evaluated vs still mismatch.
    if hueckel_mismatches is not None:
        radius1_residual_ids = [
            m["row_id"] for m in hueckel_mismatches if m["bucket"] == "radius1_numeric_mismatch"
        ]
        status_by_id = {idx: chem.get("aromaticity_status") for idx, chem, _ in success_rows}
        for e in error_rows:
            status_by_id[e["row_id"]] = e["aromaticity_status"]
        exact_by_id = {idx for idx, chem, rd in success_rows if is_exact_match(chem, rd)}

        resolved = 0
        not_evaluable = 0
        still_mismatching = 0
        for rid in radius1_residual_ids:
            status = status_by_id.get(rid)
            if status == "rdkit_parity_success":
                if rid in exact_by_id:
                    resolved += 1
                else:
                    still_mismatching += 1
            else:
                not_evaluable += 1

        summary["hueckel_radius1_residuals"] = len(radius1_residual_ids)
        summary["resolved_by_rdkit_parity"] = resolved
        summary["not_evaluable_due_to_aromaticity_error"] = not_evaluable
        summary["still_mismatching"] = still_mismatching

    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chematic-variant", required=True)
    p.add_argument("--rdkit-oracle", required=True)
    p.add_argument("--chematic-hueckel", default=None)
    p.add_argument("--hueckel-mismatches", default=None)
    p.add_argument("--summary-out", default=None)
    args = p.parse_args()

    variant_rows = load_jsonl(args.chematic_variant)
    oracle_rows = load_jsonl(args.rdkit_oracle)
    hueckel_rows = load_jsonl(args.chematic_hueckel) if args.chematic_hueckel else None
    hueckel_mismatches = load_jsonl(args.hueckel_mismatches) if args.hueckel_mismatches else None

    summary = run(variant_rows, oracle_rows, hueckel_rows, hueckel_mismatches)

    print(json.dumps(summary, indent=2))
    if args.summary_out:
        with open(args.summary_out, "w") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
        print(f"summary written to {args.summary_out}")


if __name__ == "__main__":
    main()
