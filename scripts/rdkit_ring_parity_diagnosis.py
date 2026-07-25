#!/usr/bin/env python3
"""SMARTS-R2: cross-reference chematic-smarts's opt-in RDKit-parity matching
mode (`find_matches_rdkit_parity`, `crates/chematic-smarts/src/
rdkit_ring_model.rs`) against a live RDKit oracle.

Reads validation/results/rdkit_ring_parity_dump.jsonl (one row per molecule,
emitted by `cargo run -p chematic-smarts --release --example
rdkit_parity_dump`), which itself contains chematic's DEFAULT match set and
the opt-in RDKit-PARITY match set for every (molecule, pattern) cell. This
script independently recomputes RDKit's own match set for every cell (never
trusts a comparison chematic already computed) and classifies each cell into
one of a fixed set of named buckets. Every mismatch must land in a named
bucket with a real mechanism -- an unrecognized bucket name is a hard error
(fail-closed), not a silent skip.

Also verifies chematic-atom-i <-> RDKit-atom-i alignment (element-by-element)
before trusting any index-based match-set comparison for a given molecule --
added specifically because an unverified alignment assumption would silently
invalidate every downstream comparison for that molecule.

Run:
    /tmp/chematic-smartsC-venv/bin/python scripts/rdkit_ring_parity_diagnosis.py
"""

import json
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
DUMP_PATH = ROOT / "validation" / "results" / "rdkit_ring_parity_dump.jsonl"
SUMMARY_PATH = ROOT / "validation" / "results" / "rdkit_ring_parity_diagnosis_summary.json"

RDKIT_PINNED_COMMIT = "8afba32ec539dcb2369bc84549d802aca3f7eb39"
EXPECTED_RDKIT_VERSION = "2026.03.3"

# `[RN]`-family patterns: the one primitive this mode's ring model actually
# changes (see rdkit_ring_model.rs's module doc comment for the proof that
# every other ring-shaped primitive is basis-invariant).
RING_COUNT_PATTERNS = {"[R]", "[R0]", "[R1]", "[R2]", "[R3]"}

# Fixed bucket whitelist -- an unrecognized bucket name aborts the run.
EXPECTED_BUCKETS = {
    "agree_all",  # default == parity == rdkit
    "parity_fixes_default",  # default != rdkit, parity == rdkit (the fix working)
    "parity_regresses",  # default == rdkit, parity != rdkit (would be a bug)
    "parity_worsens_agreement",  # default != rdkit, parity != rdkit, and parity's
    # set differs from default's (candidate-generation mismatch, see below)
    "both_disagree_same_as_default",  # default != rdkit, parity != rdkit, parity ==
    # default (opt-in mode made no difference -- residual carried over unchanged)
    "chematic_parse_error",
    "rdkit_smarts_parse_error",
    "chematic_parity_error",
}


def match_set(matches):
    return frozenset(frozenset(m) for m in matches)


def main():
    if rdBase.rdkitVersion != EXPECTED_RDKIT_VERSION:
        sys.exit(
            f"FATAL: expected rdkit=={EXPECTED_RDKIT_VERSION}, got "
            f"{rdBase.rdkitVersion}. Re-run in /tmp/chematic-smartsC-venv."
        )

    if not DUMP_PATH.exists():
        sys.exit(
            f"FATAL: {DUMP_PATH} not found. Run:\n"
            "  cargo run -p chematic-smarts --release --example rdkit_parity_dump "
            f"-- ~/Downloads/SMILES.csv > {DUMP_PATH}"
        )

    rows = [json.loads(line) for line in open(DUMP_PATH) if line.strip()]
    if not rows:
        sys.exit(f"FATAL: {DUMP_PATH} is empty")

    ids_seen = set()
    for r in rows:
        if r["id"] in ids_seen:
            sys.exit(f"FATAL: duplicate id {r['id']!r} in dump -- aborting (fail-closed)")
        ids_seen.add(r["id"])

    # Pre-compile every distinct pattern once with RDKit.
    all_patterns = sorted({p for r in rows for p in r["patterns"]})
    rd_queries = {}
    for p in all_patterns:
        rd_queries[p] = Chem.MolFromSmarts(p)

    bucket_counts = Counter()
    ring_count_bucket_counts = Counter()
    alignment_failures = []
    mismatch_examples = {b: [] for b in EXPECTED_BUCKETS}
    n_molecules = 0
    n_alignment_checked = 0

    for row in rows:
        smi = row["smiles"]
        rm = Chem.MolFromSmiles(smi)
        if rm is None:
            # chematic parsed it (it's in the dump) but RDKit didn't --
            # exclude from comparison, record separately, not silently dropped.
            alignment_failures.append({"id": row["id"], "smiles": smi, "reason": "rdkit_parse_failed"})
            continue

        # --- Alignment check: chematic atom i <-> RDKit atom i, element-by-element ---
        chematic_elems = row["atom_elements"]
        rdkit_elems = [a.GetSymbol() for a in rm.GetAtoms()]
        n_alignment_checked += 1
        if chematic_elems != rdkit_elems:
            alignment_failures.append(
                {
                    "id": row["id"],
                    "smiles": smi,
                    "reason": "element_mismatch",
                    "chematic": chematic_elems,
                    "rdkit": rdkit_elems,
                }
            )
            continue

        n_molecules += 1

        for pat, entry in row["patterns"].items():
            if entry.get("parse_error"):
                bucket_counts["chematic_parse_error"] += 1
                continue
            rq = rd_queries.get(pat)
            if rq is None:
                bucket_counts["rdkit_smarts_parse_error"] += 1
                continue
            rdkit_set = match_set(rm.GetSubstructMatches(rq, uniquify=True))

            default_set = frozenset(frozenset(m) for m in entry["default"])

            if "parity_error" in entry:
                bucket_counts["chematic_parity_error"] += 1
                if len(mismatch_examples["chematic_parity_error"]) < 10:
                    mismatch_examples["chematic_parity_error"].append(
                        {"id": row["id"], "smiles": smi, "pattern": pat, "error": entry["parity_error"]}
                    )
                continue

            parity_set = frozenset(frozenset(m) for m in entry["parity"])

            default_ok = default_set == rdkit_set
            parity_ok = parity_set == rdkit_set

            if default_ok and parity_ok:
                bucket = "agree_all"
            elif (not default_ok) and parity_ok:
                bucket = "parity_fixes_default"
            elif default_ok and (not parity_ok):
                bucket = "parity_regresses"
            else:
                # both disagree with RDKit
                if parity_set == default_set:
                    bucket = "both_disagree_same_as_default"
                else:
                    bucket = "parity_worsens_agreement"

            bucket_counts[bucket] += 1
            if pat in RING_COUNT_PATTERNS:
                ring_count_bucket_counts[bucket] += 1
            if bucket not in EXPECTED_BUCKETS:
                sys.exit(f"FATAL: unrecognized bucket {bucket!r} -- fail-closed abort")
            if bucket != "agree_all" and len(mismatch_examples[bucket]) < 15:
                mismatch_examples[bucket].append(
                    {
                        "id": row["id"],
                        "smiles": smi,
                        "pattern": pat,
                        "default": sorted(sorted(s) for s in default_set),
                        "parity": sorted(sorted(s) for s in parity_set),
                        "rdkit": sorted(sorted(s) for s in rdkit_set),
                    }
                )

    total_cells = sum(bucket_counts.values())

    summary = {
        "rdkit_version": rdBase.rdkitVersion,
        "rdkit_pinned_source_commit": RDKIT_PINNED_COMMIT,
        "n_rows_in_dump": len(rows),
        "n_molecules_compared": n_molecules,
        "n_alignment_checked": n_alignment_checked,
        "n_alignment_failures": len(alignment_failures),
        "alignment_failures": alignment_failures[:20],
        "total_cells": total_cells,
        "bucket_counts": dict(bucket_counts),
        "ring_count_pattern_bucket_counts": dict(ring_count_bucket_counts),
        "mismatch_examples": mismatch_examples,
    }
    SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(SUMMARY_PATH, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"rdkit version: {rdBase.rdkitVersion}  (expected {EXPECTED_RDKIT_VERSION})")
    print(f"molecules in dump: {len(rows)}, RDKit-parseable+aligned: {n_molecules}")
    print(f"alignment failures: {len(alignment_failures)}")
    print(f"total (molecule, pattern) cells compared: {total_cells}")
    print()
    for b in sorted(bucket_counts):
        pct = 100 * bucket_counts[b] / total_cells if total_cells else 0
        print(f"  {b:32s} {bucket_counts[b]:8d}  ({pct:.4f}%)")
    print()
    print("[RN] family only (the ring-count model's actual target):")
    ring_total = sum(ring_count_bucket_counts.values())
    for b in sorted(ring_count_bucket_counts):
        pct = 100 * ring_count_bucket_counts[b] / ring_total if ring_total else 0
        print(f"  {b:32s} {ring_count_bucket_counts[b]:8d}  ({pct:.4f}%)  / {ring_total}")
    print()
    agree_all = bucket_counts.get("agree_all", 0)
    print(f"overall agreement (default AND parity both == rdkit): {agree_all}/{total_cells} "
          f"({100*agree_all/total_cells:.4f}%)")
    default_agree = agree_all + bucket_counts.get("parity_regresses", 0)
    parity_agree = agree_all + bucket_counts.get("parity_fixes_default", 0)
    print(f"default-matcher agreement with rdkit:  {default_agree}/{total_cells} "
          f"({100*default_agree/total_cells:.4f}%)")
    print(f"parity-matcher agreement with rdkit:    {parity_agree}/{total_cells} "
          f"({100*parity_agree/total_cells:.4f}%)")
    print()
    print(f"wrote {SUMMARY_PATH.relative_to(ROOT)}")

    # Fail-closed: parity_regresses must be zero or explicitly investigated --
    # this run does not assert exit(1) on it (a real regression may need a
    # design decision, not just a script abort) but it prints loudly.
    if bucket_counts.get("parity_regresses", 0) > 0:
        print(
            f"\nWARNING: {bucket_counts['parity_regresses']} cells where the opt-in "
            "mode REGRESSES vs the default matcher's agreement with RDKit -- "
            "investigate before shipping.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
