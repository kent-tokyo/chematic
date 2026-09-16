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
import os
import sys
import tempfile
import argparse
import hashlib
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
DUMP_PATH = ROOT / "validation" / "results" / "rdkit_ring_parity_dump.jsonl"
SUMMARY_PATH = ROOT / "validation" / "results" / "rdkit_ring_parity_diagnosis_summary.json"

RDKIT_PINNED_COMMIT = "8afba32ec539dcb2369bc84549d802aca3f7eb39"
EXPECTED_RDKIT_VERSION = "2026.03.6"

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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump", type=Path, default=DUMP_PATH)
    parser.add_argument("--summary", type=Path, default=SUMMARY_PATH)
    parser.add_argument("--expected-rdkit-version", default=EXPECTED_RDKIT_VERSION)
    parser.add_argument(
        "--rdkit-pinned-source",
        default=RDKIT_PINNED_COMMIT,
        help="commit, package identifier, or other immutable comparator provenance",
    )
    args = parser.parse_args()
    dump_path = args.dump.resolve()
    summary_path = args.summary.resolve()

    if rdBase.rdkitVersion != args.expected_rdkit_version:
        sys.exit(
            f"FATAL: expected rdkit=={args.expected_rdkit_version}, got "
            f"{rdBase.rdkitVersion}. Re-run in /tmp/chematic-smartsC-venv."
        )

    if not dump_path.exists():
        sys.exit(
            f"FATAL: {dump_path} not found. Run:\n"
            "  cargo run -p chematic-smarts --release --example rdkit_parity_dump "
            f"-- ~/Downloads/SMILES.csv > {dump_path}"
        )

    records = [json.loads(line) for line in dump_path.open(encoding="utf-8") if line.strip()]
    if not records:
        sys.exit(f"FATAL: {dump_path} is empty")

    footers = [record for record in records if record.get("record_type") == "footer"]
    rows = [record for record in records if record.get("record_type", "molecule") == "molecule"]
    if len(footers) != 1 or records[-1].get("record_type") != "footer":
        sys.exit(
            "FATAL: dump must end with exactly one completed footer; "
            "a missing footer means the producer did not reach a clean exit"
        )
    footer = footers[0]
    if footer.get("completed") is not True:
        sys.exit("FATAL: dump footer is not marked completed")
    if footer.get("emitted_molecule_rows") != len(rows):
        sys.exit("FATAL: footer row count does not match molecule records")
    expected_input_rows = footer.get("input_rows")
    if not isinstance(expected_input_rows, int) or expected_input_rows < 0:
        sys.exit("FATAL: footer input_rows is missing or invalid")
    if footer.get("hand_corpus_rows") + expected_input_rows != len(rows):
        sys.exit("FATAL: footer source accounting does not match molecule records")
    source_path = Path(str(footer.get("source_path", "")))
    if not source_path.is_absolute():
        source_path = ROOT / source_path
    source_corpus_sha256 = sha256(source_path) if source_path.is_file() else None

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
    # Examples are deliberately capped for routine reports, but unresolved
    # parity cells must remain fully auditable: a truncated sample can make a
    # single apparent scaffold hide a second mechanism. Keep this inventory
    # only for cells where both matchers disagree with RDKit; parser errors and
    # candidate regressions retain their existing separate accounting.
    unresolved_residuals = []
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
            if bucket == "both_disagree_same_as_default":
                unresolved_residuals.append(
                    {
                        "id": row["id"],
                        "smiles": smi,
                        "pattern": pat,
                        "default": sorted(sorted(s) for s in default_set),
                        "rdkit": sorted(sorted(s) for s in rdkit_set),
                    }
                )

    total_cells = sum(bucket_counts.values())

    def display_path(path: Path) -> str:
        try:
            return path.relative_to(ROOT).as_posix()
        except ValueError:
            return path.name

    summary = {
        "rdkit_version": rdBase.rdkitVersion,
        "rdkit_pinned_source": args.rdkit_pinned_source,
        "dump_path": display_path(dump_path),
        "dump_sha256": sha256(dump_path),
        "source_corpus_path": display_path(source_path),
        "source_corpus_sha256": source_corpus_sha256,
        "source_corpus_rows": expected_input_rows,
        "comparison_config": {
            "match_set": "sorted atom-index sets",
            "rdkit_uniquify": True,
            "unsupported_policy": "count and report; never remove from denominator",
            "requires_completed_footer": True,
        },
        "n_rows_in_dump": len(rows),
        "dump_footer": footer,
        "n_molecules_compared": n_molecules,
        "n_alignment_checked": n_alignment_checked,
        "n_alignment_failures": len(alignment_failures),
        "alignment_failures": alignment_failures[:20],
        "total_cells": total_cells,
        "bucket_counts": dict(bucket_counts),
        "ring_count_pattern_bucket_counts": dict(ring_count_bucket_counts),
        "mismatch_examples": mismatch_examples,
        "unresolved_residuals": sorted(
            unresolved_residuals,
            key=lambda row: (row["id"], row["pattern"]),
        ),
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    # Do not truncate the last good diagnosis if serialization or the process
    # is interrupted. A failed comparator run must not leave an empty artifact
    # that can be mistaken for a completed measurement.
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=summary_path.parent, prefix=f".{summary_path.name}.", delete=False
    ) as f:
        json.dump(summary, f, indent=2)
        f.write("\n")
        temporary_path = Path(f.name)
    os.replace(temporary_path, summary_path)

    print(f"rdkit version: {rdBase.rdkitVersion}  (expected {args.expected_rdkit_version})")
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
    try:
        display_summary = summary_path.relative_to(ROOT)
    except ValueError:
        display_summary = summary_path
    print(f"wrote {display_summary}")

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
