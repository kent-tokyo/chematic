#!/usr/bin/env python3
"""Build one independently attributable SMARTS/RDKit evidence arm.

Run this once in each clean baseline/candidate checkout, then compare the two
summaries with ``check_smarts_baseline_candidate.py``. A successful process is
not enough: the runner requires the producer footer and records the executable,
source snapshot, corpus, and query-set digests needed for an auditable pair.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(checkout: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=checkout, check=True, capture_output=True, text=True
    ).stdout.strip()


def atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False) as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary_path = Path(handle.name)
    os.replace(temporary_path, path)


def query_digest(dump_path: Path) -> str:
    patterns: set[str] = set()
    with dump_path.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            if record.get("record_type", "molecule") == "molecule":
                patterns.update(record.get("patterns", {}))
    if not patterns:
        raise ValueError("dump contains no query patterns")
    return sha256_bytes(json.dumps(sorted(patterns), separators=(",", ":")).encode("utf-8"))


def corpus_rows(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line.strip())


def promote_completed_dump(raw_path: Path, dump_path: Path, corpus_path: Path) -> None:
    """Atomically add a verifier footer when an older producer lacks one.

    The wrapper only performs this after the child exits successfully and every
    emitted molecule row is parseable, unique, and accounted for against the
    requested corpus. It records the footer origin so an old producer is never
    presented as if it implemented the newer footer protocol itself.
    """
    records = [json.loads(line) for line in raw_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    footers = [record for record in records if record.get("record_type") == "footer"]
    rows = [record for record in records if record.get("record_type", "molecule") == "molecule"]
    if not rows or len(footers) > 1 or (footers and records[-1].get("record_type") != "footer"):
        raise ValueError("producer rows are empty or footer placement is invalid")
    ids = [record.get("id") for record in rows]
    if any(not isinstance(identifier, str) or not identifier for identifier in ids) or len(set(ids)) != len(ids):
        raise ValueError("producer emitted missing or duplicate molecule IDs")
    if any(not isinstance(record.get("smiles"), str) or not isinstance(record.get("patterns"), dict) for record in rows):
        raise ValueError("producer emitted malformed molecule rows")
    expected_input_rows = corpus_rows(corpus_path)
    actual_input_rows = sum(identifier.startswith("corpus_") for identifier in ids)
    if actual_input_rows != expected_input_rows:
        raise ValueError(f"producer corpus row accounting {actual_input_rows} != expected {expected_input_rows}")
    parse_failures = sum("parse_error" in record for record in rows)
    if footers:
        footer = footers[0]
        if footer.get("completed") is not True or footer.get("emitted_molecule_rows") != len(rows):
            raise ValueError("producer footer is incomplete or has wrong row count")
        if footer.get("input_rows") != expected_input_rows:
            raise ValueError("producer footer input row count is wrong")
        if footer.get("parse_failures") != parse_failures:
            raise ValueError("producer footer parse failure count is wrong")
        footer = {**footer, "footer_origin": "producer"}
    else:
        footer = {
            "record_type": "footer",
            "completed": True,
            "footer_origin": "runner_verified",
            "producer_exit_code": 0,
            "source_path": str(corpus_path),
            "hand_corpus_rows": len(rows) - actual_input_rows,
            "input_rows": actual_input_rows,
            "emitted_molecule_rows": len(rows),
            "parse_failures": parse_failures,
        }
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=dump_path.parent, prefix=f".{dump_path.name}.", delete=False) as handle:
        for record in rows:
            handle.write(json.dumps(record, separators=(",", ":")) + "\n")
        handle.write(json.dumps(footer, separators=(",", ":")) + "\n")
        temporary_path = Path(handle.name)
    os.replace(temporary_path, dump_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--role", choices=("baseline", "candidate"), required=True)
    parser.add_argument("--checkout", type=Path, default=ROOT)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-rdkit-version", default="2025.09.3")
    parser.add_argument("--rdkit-pinned-source", default="installed-python-package-2025.09.3")
    args = parser.parse_args()

    checkout = args.checkout.resolve()
    corpus = args.corpus.resolve()
    output_dir = args.output_dir.resolve()
    if not corpus.is_file():
        parser.error(f"corpus not found: {corpus}")
    if git(checkout, "status", "--porcelain"):
        parser.error("checkout must be clean; use a separate committed baseline/candidate worktree")

    try:
        commit = git(checkout, "rev-parse", "HEAD")
        tracked_snapshot = git(checkout, "ls-files", "-s").encode("utf-8")
        source_diff = git(checkout, "diff", "--no-ext-diff", "HEAD").encode("utf-8")
        subprocess.run(
            ["cargo", "build", "-p", "chematic-smarts", "--release", "--example", "rdkit_parity_dump"],
            cwd=checkout,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        print(f"SMARTS evidence build failed: {exc}", file=sys.stderr)
        return exc.returncode or 1

    binary = checkout / "target" / "release" / "examples" / "rdkit_parity_dump"
    if not binary.is_file():
        print(f"SMARTS evidence build produced no example binary: {binary}", file=sys.stderr)
        return 1
    output_dir.mkdir(parents=True, exist_ok=True)
    dump_path = output_dir / f"rdkit-smarts-{args.role}.jsonl"
    raw_dump_path = output_dir / f".rdkit-smarts-{args.role}.raw.jsonl"
    summary_path = output_dir / f"rdkit-smarts-{args.role}.json"
    with raw_dump_path.open("w", encoding="utf-8") as dump_handle:
        run = subprocess.run([str(binary), str(corpus)], cwd=checkout, stdout=dump_handle)
    if run.returncode != 0:
        print(f"SMARTS evidence producer exited {run.returncode}; no summary was promoted", file=sys.stderr)
        return run.returncode
    try:
        promote_completed_dump(raw_dump_path, dump_path, corpus)
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "rdkit_ring_parity_diagnosis.py"),
                "--dump", str(dump_path),
                "--summary", str(summary_path),
                "--expected-rdkit-version", args.expected_rdkit_version,
                "--rdkit-pinned-source", args.rdkit_pinned_source,
            ],
            cwd=checkout,
            check=True,
        )
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary.update(
            {
                "schema_version": 1,
                "role": args.role,
                "provenance": {
                    "source_commit": commit,
                    "source_tree_sha256": sha256_bytes(tracked_snapshot),
                    "source_diff_sha256": sha256_bytes(source_diff),
                    "binary_sha256": sha256(binary),
                    "corpus_sha256": sha256(corpus),
                    "query_sha256": query_digest(dump_path),
                },
            }
        )
        atomic_json(summary_path, summary)
    except (OSError, ValueError, json.JSONDecodeError, subprocess.CalledProcessError) as exc:
        print(f"SMARTS evidence diagnosis failed: {exc}", file=sys.stderr)
        return 1
    print(f"SMARTS {args.role} evidence: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
