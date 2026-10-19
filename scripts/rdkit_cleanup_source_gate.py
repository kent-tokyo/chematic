#!/usr/bin/env python3
"""Compare source-built Cleanup JSONL with RDKit on the stress corpus.

Generate the JSONL with ``cargo run -p chematic-chem --example
cleanup_rdkit_dump --offline < validation/cosmolkit_comparison/corpora/stress-v1.smi``.
The shared compatible subset also requires a successful RDKit Cleanup oracle.
The separate raw Reionize count scores every source row.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

from rdkit_standardization_classification import (
    classify_cleanup_mismatch,
    classify_normalize_mismatch,
    classify_reionize_mismatch,
)

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "validation/cosmolkit_comparison/corpora/stress-v1.smi"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def value(result: dict) -> tuple[str | None, str | None]:
    if "Ok" in result:
        return result["Ok"], None
    return None, result.get("Err", "missing result")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--previous-source", type=Path)
    args = parser.parse_args()
    if args.baseline and args.previous_source:
        parser.error("choose one comparison report")
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1; found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")

    corpus = [
        line.split()[0] if line.split() else None
        for line in CORPUS.read_text(encoding="utf-8").splitlines()
    ]
    records = [json.loads(line) for line in args.dump.read_text().splitlines()]
    if len(records) != sum(value is not None for value in corpus):
        raise ValueError("source dump does not cover every nonblank corpus row")
    rows = {}
    for record in records:
        index = record["index"]
        if index in rows or index >= len(corpus) or record["input"] != corpus[index]:
            raise ValueError(f"source dump index/input mismatch at {index}")
        rows[index] = record

    compatible, exact, changed = set(), set(), set()
    normalize_compatible, normalize_exact = set(), set()
    reionize_compatible, reionize_exact, reionize_all_exact = set(), set(), set()
    excluded, errors, mismatches = [], [], []
    normalize_oracle_errors, reionize_oracle_errors = [], []
    reionize_all_oracle_errors = []
    normalize_mismatches, reionize_mismatches = [], []
    cleanup_classes: Counter[str] = Counter()
    cleanup_statuses: Counter[str] = Counter()
    normalize_classes: Counter[str] = Counter()
    normalize_statuses: Counter[str] = Counter()
    reionize_classes: Counter[str] = Counter()
    reionize_statuses: Counter[str] = Counter()
    for index, source in sorted(rows.items()):
        input_text = source["input"]
        rdkit_mol = Chem.MolFromSmiles(input_text)
        rdkit_input = (
            Chem.MolToSmiles(rdkit_mol, canonical=True, isomericSmiles=True)
            if rdkit_mol is not None else None
        )
        try:
            expected_reionize_all = Chem.MolToSmiles(
                rdMolStandardize.Reionize(rdkit_mol), canonical=True, isomericSmiles=True
            ) if rdkit_mol is not None else None
        except Exception as exc:
            expected_reionize_all = None
            reionize_all_oracle_errors.append({
                "index": index, "input": input_text,
                "oracle_error": f"{type(exc).__name__}: {exc}",
            })
        actual_reionize_all, _ = value(source.get("reionize", {}))
        if expected_reionize_all is not None and actual_reionize_all == expected_reionize_all:
            reionize_all_exact.add(index)
        source_input, parse_error = value(source.get("rdkit_input", {}))
        if parse_error is not None or source_input != rdkit_input:
            excluded.append({"index": index, "input": input_text,
                             "source_input": source_input, "rdkit_input": rdkit_input})
            continue
        try:
            expected = Chem.MolToSmiles(
                rdMolStandardize.Cleanup(rdkit_mol), canonical=True, isomericSmiles=True
            )
        except Exception as exc:
            errors.append({"index": index, "input": input_text,
                           "oracle_error": f"{type(exc).__name__}: {exc}"})
            continue
        compatible.add(index)
        try:
            expected_normalize = Chem.MolToSmiles(
                rdMolStandardize.Normalize(rdkit_mol), canonical=True, isomericSmiles=True
            )
            normalize_compatible.add(index)
        except Exception as exc:
            expected_normalize = None
            normalize_oracle_errors.append({
                "index": index,
                "input": input_text,
                "oracle_error": f"{type(exc).__name__}: {exc}",
            })
        if expected_reionize_all is not None:
            expected_reionize = expected_reionize_all
            reionize_compatible.add(index)
        else:
            expected_reionize = None
            reionize_oracle_errors.append({
                "index": index,
                "input": input_text,
                "oracle_error": "Reionize oracle failed or input could not be parsed",
            })
        actual_normalize, normalize_error = value(source.get("normalize", {}))
        actual, source_error = value(source.get("cleanup", {}))
        actual_reionize, reionize_error = value(source.get("reionize", {}))
        if index not in normalize_compatible:
            pass
        elif actual_normalize == expected_normalize:
            normalize_exact.add(index)
        else:
            classification, support_status = classify_normalize_mismatch(
                expected_normalize, actual_normalize, normalize_error
            )
            normalize_classes[classification] += 1
            normalize_statuses[support_status] += 1
            normalize_mismatches.append({
                "index": index,
                "input": input_text,
                "rdkit_normalize": expected_normalize,
                "chematic_normalize": actual_normalize,
                "source_error": normalize_error,
                "classification": classification,
                "support_status": support_status,
            })
        if actual == expected:
            exact.add(index)
        else:
            classification, support_status = classify_cleanup_mismatch(
                rdkit_mol, expected_normalize, expected, actual, source_error
            )
            cleanup_classes[classification] += 1
            cleanup_statuses[support_status] += 1
            mismatches.append({
                "index": index,
                "input": input_text,
                "rdkit_cleanup": expected,
                "chematic_cleanup": actual,
                "source_error": source_error,
                "classification": classification,
                "support_status": support_status,
            })
        if index not in reionize_compatible:
            pass
        elif actual_reionize == expected_reionize:
            reionize_exact.add(index)
        else:
            classification, support_status = classify_reionize_mismatch(
                rdkit_input, expected_reionize, actual_reionize, reionize_error
            )
            reionize_classes[classification] += 1
            reionize_statuses[support_status] += 1
            reionize_mismatches.append({
                "index": index,
                "input": input_text,
                "rdkit_reionize": expected_reionize,
                "chematic_reionize": actual_reionize,
                "source_error": reionize_error,
                "classification": classification,
                "support_status": support_status,
            })
        if expected != rdkit_input:
            changed.add(index)

    comparison = None
    if args.baseline:
        baseline = json.loads(args.baseline.read_text())
        if baseline["sampled_rows"] != len(corpus):
            raise ValueError("baseline and source sample counts differ")
        old_scored = set(range(len(corpus))) - {
            row["sample_index"] for field in ("excluded", "oracle_errors")
            for row in baseline[field]
        }
        if old_scored != compatible:
            raise ValueError("baseline and source input compatibility sets differ")
        old_mismatches = {row["sample_index"] for row in baseline["mismatches"]}
        comparison = {
            "baseline_exact": baseline["cleanup_exact"],
            "fixed_indices": sorted(old_mismatches & exact),
            "regressed_indices": sorted((old_scored - old_mismatches) - exact),
        }
    if args.previous_source:
        previous = json.loads(args.previous_source.read_text())
        if (previous["source_rows"] != len(corpus)
                or previous["compatible_rows"] != len(compatible)
                or previous["corpus_sha256"] != sha256(CORPUS)
                or previous["rdkit_version"] != rdBase.rdkitVersion
                or {row["index"] for row in previous["oracle_errors"]}
                    != {row["index"] for row in errors}
                or {row["index"] for row in previous["excluded"]}
                    != {row["index"] for row in excluded}):
            raise ValueError("previous source report uses a different corpus or oracle")
        previous_mismatches = {row["index"] for row in previous["mismatches"]}
        current_mismatches = {row["index"] for row in mismatches}
        comparison = {
            "previous_source": str(args.previous_source),
            "previous_exact": previous["exact_rows"],
            "fixed_indices": sorted(previous_mismatches - current_mismatches),
            "regressed_indices": sorted(current_mismatches - previous_mismatches),
        }

    report = {
        "profile": "source_cleanup_vs_rdkit_2026_09_1",
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                ROOT / "crates/chematic-chem/src/standardize.rs",
                ROOT / "crates/chematic-chem/src/standardize/reionize.rs",
                ROOT / "crates/chematic-chem/src/standardize/cleanup.rs",
                ROOT / "crates/chematic-chem/src/standardize/groups.rs",
                ROOT / "crates/chematic-core/src/molecule/builder.rs",
                ROOT / "crates/chematic-chem/examples/cleanup_rdkit_dump.rs",
                ROOT / "scripts/rdkit_cleanup_source_gate.py",
                ROOT / "crates/chematic-smiles/src/lib.rs",
                ROOT / "crates/chematic-smiles/src/rdkit/mod.rs",
            )
        },
        "corpus": str(CORPUS.relative_to(ROOT)),
        "corpus_sha256": sha256(CORPUS),
        "source_dump_sha256": sha256(args.dump),
        "rdkit_version": rdBase.rdkitVersion,
        "source_rows": len(corpus),
        "compatible_rows": len(compatible),
        "normalize_compatible_rows": len(normalize_compatible),
        "normalize_exact_rows": len(normalize_exact),
        "normalize_oracle_errors": normalize_oracle_errors,
        "normalize_classification_counts": dict(sorted(normalize_classes.items())),
        "normalize_support_status_counts": dict(sorted(normalize_statuses.items())),
        "normalize_mismatches": normalize_mismatches,
        "exact_rows": len(exact),
        "classification_counts": dict(sorted(cleanup_classes.items())),
        "support_status_counts": dict(sorted(cleanup_statuses.items())),
        "reionize_compatible_rows": len(reionize_compatible),
        "reionize_exact_rows": len(reionize_exact),
        "reionize_all_rows_exact": len(reionize_all_exact),
        "reionize_all_rows_oracle_errors": reionize_all_oracle_errors,
        "reionize_oracle_errors": reionize_oracle_errors,
        "reionize_classification_counts": dict(sorted(reionize_classes.items())),
        "reionize_support_status_counts": dict(sorted(reionize_statuses.items())),
        "reionize_mismatches": reionize_mismatches,
        "rdkit_changed_rows": len(changed),
        "excluded": excluded,
        "oracle_errors": errors,
        "mismatches": mismatches,
        "baseline_comparison": comparison,
        "limits": [
            "shared compatible metrics require the same parsed input and a successful Cleanup oracle",
            "raw all-row Reionize count includes rows where Cleanup fails or input interpretations differ",
            "RDKit agreement is not independent chemical truth",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"Normalize: {len(normalize_exact)}/{len(normalize_compatible)} exact")
    print(f"Cleanup: {len(exact)}/{len(compatible)} exact")
    print(f"Reionize (Cleanup-eligible rows): {len(reionize_exact)}/{len(reionize_compatible)} exact")
    print(f"Reionize raw all rows: {len(reionize_all_exact)}/{len(corpus)} exact")
    if comparison:
        print(f"Fixed: {len(comparison['fixed_indices'])}; "
              f"regressed: {len(comparison['regressed_indices'])}")


if __name__ == "__main__":
    main()
