#!/usr/bin/env python3
"""Validate the checked-in footer-verified SMARTS evidence packet."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fail(message: str) -> None:
    raise SystemExit(f"SMARTS evidence invalid: {message}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "validation/results/rdkit-smarts-direct-chembl-5000-footer-verified-v1.0.13.json",
    )
    args = parser.parse_args()
    report_path = args.report if args.report.is_absolute() else ROOT / args.report
    if not report_path.is_file():
        fail(f"report not found: {report_path}")
    report = json.loads(report_path.read_text(encoding="utf-8"))

    required = {
        "rdkit_version",
        "rdkit_pinned_source",
        "source_corpus_path",
        "source_corpus_sha256",
        "source_corpus_rows",
        "comparison_config",
        "n_rows_in_dump",
        "dump_footer",
        "n_molecules_compared",
        "n_alignment_checked",
        "n_alignment_failures",
        "total_cells",
        "bucket_counts",
    }
    missing = sorted(required - report.keys())
    if missing:
        fail(f"missing fields: {', '.join(missing)}")

    config = report["comparison_config"]
    if config.get("requires_completed_footer") is not True:
        fail("completion footer is not required")
    footer = report["dump_footer"]
    if footer.get("record_type") != "footer" or footer.get("completed") is not True:
        fail("dump footer is absent or incomplete")
    if footer.get("input_rows") != report["source_corpus_rows"]:
        fail("footer input row count does not match source corpus rows")
    if footer.get("emitted_molecule_rows") != report["n_molecules_compared"]:
        fail("footer molecule count does not match report")
    if footer.get("parse_failures") != 0:
        fail("source dump contains parse failures")
    if report["n_alignment_failures"] != 0 or report["n_alignment_checked"] != report["n_molecules_compared"]:
        fail("molecule alignment is incomplete")

    source = ROOT / report["source_corpus_path"]
    if not source.is_file():
        fail(f"source corpus not found: {source}")
    if sha256(source) != report["source_corpus_sha256"]:
        fail("source corpus SHA-256 mismatch")
    actual_rows = sum(1 for line in source.read_text(encoding="utf-8").splitlines() if line.strip())
    if actual_rows != report["source_corpus_rows"]:
        fail("source corpus row count mismatch")

    buckets = report["bucket_counts"]
    if sum(buckets.values()) != report["total_cells"]:
        fail("bucket counts do not sum to total cells")
    if report["n_rows_in_dump"] != report["n_molecules_compared"]:
        fail("dump row accounting mismatch")
    if report["total_cells"] <= 0:
        fail("empty comparison")

    print(
        "RDKit SMARTS evidence OK: "
        f"{report['n_molecules_compared']} molecules, "
        f"{report['total_cells']} cells, "
        f"{buckets.get('both_disagree_same_as_default', 0)} residuals, "
        f"{buckets.get('rdkit_smarts_parse_error', 0)} RDKit parse errors"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
