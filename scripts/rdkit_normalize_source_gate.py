#!/usr/bin/env python3
"""Compare source-built Normalize output with pinned RDKit on fixed corpora.

Generate each JSONL with ``cargo run -p chematic-chem --example
cleanup_rdkit_dump --offline < CORPUS``. Both versions are scored directly on
all corpus rows; the Cleanup source gate has a separate input-compatibility
filter that is not part of this Normalize measurement.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Optional

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[1]
CORPORA = (
    ("stress", "validation/cosmolkit_comparison/corpora/stress-v1.smi", 4072),
    ("exposed_10k", "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi", 10000),
)
FIXED_STRESS_INDICES = [1147, 1459, 3584, 3801]
FIXED_TENK_INDICES = [7529, 8150, 9473, 9680, 9703]
FOLLOWUP_TENK_INDICES = [1230, 8256]
PREVIOUS_REPORT = "validation/results/rdkit-normalize-source-2026-10-20.json"
SOURCE_FILES = (
    "crates/chematic-chem/src/standardize/groups.rs",
    "crates/chematic-chem/src/standardize.rs",
    "crates/chematic-chem/examples/cleanup_rdkit_dump.rs",
    "scripts/rdkit_normalize_source_gate.py",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_dump(path: Path, corpus: list[str]) -> list[dict]:
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    if len(rows) != len(corpus):
        raise ValueError(f"{path}: expected {len(corpus)} rows, found {len(rows)}")
    for index, row in enumerate(rows):
        if row["index"] != index or row["input"] != corpus[index]:
            raise ValueError(f"{path}: row {index} differs from the corpus")
    return rows


def ok_value(row: dict, field: str) -> Optional[str]:
    return row.get(field, {}).get("Ok")


def score_corpus(
    name: str, relative_path: str, expected_rows: int,
    baseline_path: Path, source_path: Path, previous_path: Optional[Path],
) -> dict:
    corpus_path = ROOT / relative_path
    smiles_rows = [
        line.split()[0] for line in corpus_path.read_text(encoding="utf-8").splitlines()
    ]
    if len(smiles_rows) != expected_rows:
        raise ValueError(f"{relative_path}: unexpected corpus length")
    baseline = read_dump(baseline_path, smiles_rows)
    source = read_dump(source_path, smiles_rows)
    previous = read_dump(previous_path, smiles_rows) if previous_path else None
    old_mismatches, new_mismatches, input_mismatches = [], [], []
    changed, cleanup_changed, reionize_changed = [], [], []
    previous_mismatches, followup_changed = [], []
    followup_cleanup_changed, followup_reionize_changed = [], []
    fixed_cases = []
    fixed_indices = FIXED_STRESS_INDICES if name == "stress" else FIXED_TENK_INDICES
    if name == "exposed_10k" and previous is not None:
        fixed_indices = sorted(fixed_indices + FOLLOWUP_TENK_INDICES)
    for index, smiles in enumerate(smiles_rows):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"RDKit cannot parse {name} row {index}")
        rdkit_input = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        if baseline[index].get("rdkit_input") != source[index].get("rdkit_input"):
            raise ValueError(f"{name} row {index}: source input changed")
        if previous and previous[index].get("rdkit_input") != source[index].get("rdkit_input"):
            raise ValueError(f"{name} row {index}: previous input changed")
        if ok_value(source[index], "rdkit_input") != rdkit_input:
            input_mismatches.append(index)
        expected = Chem.MolToSmiles(
            rdMolStandardize.Normalize(mol), canonical=True, isomericSmiles=True
        )
        old = ok_value(baseline[index], "normalize")
        new = ok_value(source[index], "normalize")
        prior = ok_value(previous[index], "normalize") if previous else None
        if old != expected:
            old_mismatches.append(index)
        if new != expected:
            new_mismatches.append(index)
        if previous and prior != expected:
            previous_mismatches.append(index)
        if old != new:
            changed.append(index)
        if baseline[index].get("cleanup") != source[index].get("cleanup"):
            cleanup_changed.append(index)
        if baseline[index].get("reionize") != source[index].get("reionize"):
            reionize_changed.append(index)
        if previous:
            if prior != new:
                followup_changed.append(index)
            if previous[index].get("cleanup") != source[index].get("cleanup"):
                followup_cleanup_changed.append(index)
            if previous[index].get("reionize") != source[index].get("reionize"):
                followup_reionize_changed.append(index)
        if index in fixed_indices:
            fixed_cases.append({
                "index": index, "input": smiles, "baseline": old,
                **({"previous": prior} if previous else {}),
                "source": new, "rdkit": expected,
            })
    old_set, new_set = set(old_mismatches), set(new_mismatches)
    result = {
        "path": relative_path,
        "corpus_sha256": sha256(corpus_path),
        "baseline_dump_sha256": sha256(baseline_path),
        "source_dump_sha256": sha256(source_path),
        "rows": len(smiles_rows),
        "input_canonical_mismatch_indices": input_mismatches,
        "baseline_exact": len(smiles_rows) - len(old_mismatches),
        "source_exact": len(smiles_rows) - len(new_mismatches),
        "baseline_mismatch_indices": old_mismatches,
        "source_mismatch_indices": new_mismatches,
        "fixed_indices": sorted(old_set - new_set),
        "regressed_indices": sorted(new_set - old_set),
        "normalize_changed_indices": changed,
        "cleanup_changed_indices": cleanup_changed,
        "reionize_changed_indices": reionize_changed,
        "fixed_cases": fixed_cases,
    }
    if previous_path:
        prior_set = set(previous_mismatches)
        result["previous"] = {
            "source_dump_sha256": sha256(previous_path),
            "exact": len(smiles_rows) - len(previous_mismatches),
            "mismatch_indices": previous_mismatches,
            "fixed_indices": sorted(prior_set - new_set),
            "regressed_indices": sorted(new_set - prior_set),
            "normalize_changed_indices": followup_changed,
            "cleanup_changed_indices": followup_cleanup_changed,
            "reionize_changed_indices": followup_reionize_changed,
        }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-stress-jsonl", type=Path, required=True)
    parser.add_argument("--source-stress-jsonl", type=Path, required=True)
    parser.add_argument("--baseline-tenk-jsonl", type=Path, required=True)
    parser.add_argument("--source-tenk-jsonl", type=Path, required=True)
    parser.add_argument("--baseline-source-commit", required=True)
    parser.add_argument("--previous-stress-jsonl", type=Path)
    parser.add_argument("--previous-tenk-jsonl", type=Path)
    parser.add_argument("--previous-source-commit")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    previous_args = (
        args.previous_stress_jsonl,
        args.previous_tenk_jsonl,
        args.previous_source_commit,
    )
    if any(value is not None for value in previous_args) and not all(
        value is not None for value in previous_args
    ):
        parser.error("all three --previous-* arguments must be supplied together")
    has_previous = all(value is not None for value in previous_args)
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1; found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")
    subprocess.run(
        ["git", "rev-parse", "--verify", args.baseline_source_commit + "^{commit}"],
        cwd=ROOT, check=True, stdout=subprocess.DEVNULL,
    )
    if has_previous:
        subprocess.run(
            ["git", "rev-parse", "--verify", args.previous_source_commit + "^{commit}"],
            cwd=ROOT, check=True, stdout=subprocess.DEVNULL,
        )
        previous_report = json.loads(subprocess.check_output(
            ["git", "show", f"HEAD:{PREVIOUS_REPORT}"], cwd=ROOT, text=True,
        ))
        if previous_report["source_commit"] != args.previous_source_commit:
            parser.error("previous source commit does not match the recorded report")
        if previous_report["baseline_source_commit"] != args.baseline_source_commit:
            parser.error("baseline source commit does not match the recorded report")
        for name, path in (("stress", args.previous_stress_jsonl),
                           ("exposed_10k", args.previous_tenk_jsonl)):
            if previous_report["corpora"][name]["source_dump_sha256"] != sha256(path):
                parser.error(f"{name}: previous dump does not match the recorded report")
        for name, path in (("stress", args.baseline_stress_jsonl),
                           ("exposed_10k", args.baseline_tenk_jsonl)):
            if previous_report["corpora"][name]["baseline_dump_sha256"] != sha256(path):
                parser.error(f"{name}: baseline dump does not match the recorded report")
    subprocess.run(
        ["git", "diff", "--quiet", "HEAD", "--", *SOURCE_FILES], cwd=ROOT,
        check=True,
    )
    dump_paths = (
        (args.baseline_stress_jsonl, args.source_stress_jsonl,
         args.previous_stress_jsonl),
        (args.baseline_tenk_jsonl, args.source_tenk_jsonl,
         args.previous_tenk_jsonl),
    )
    corpora = {
        name: score_corpus(name, path, size, baseline, source, previous)
        for (name, path, size), (baseline, source, previous) in zip(CORPORA, dump_paths)
    }
    report = {
        "schema_version": 2 if has_previous else 1,
        "profile": "source_rdkit_2026_09_1_normalize",
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "baseline_source_commit": args.baseline_source_commit,
        **({"previous_source_commit": args.previous_source_commit} if has_previous else {}),
        "rdkit_version": rdBase.rdkitVersion,
        "source_file_sha256": {path: sha256(ROOT / path) for path in SOURCE_FILES},
        "corpora": corpora,
        "limits": [
            "Exact canonical isomeric SMILES agreement on two fixed corpora only.",
            "All input SMILES are parsed independently by RDKit and CheMatic.",
            "Input canonical forms agree on both corpora in this source comparison.",
            "The CheMatic results come from source builds, not published packages.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print({name: f"{row['source_exact']}/{row['rows']}" for name, row in corpora.items()})
    stress, tenk = corpora["stress"], corpora["exposed_10k"]
    expected_tenk_fixed = (
        sorted(FIXED_TENK_INDICES + FOLLOWUP_TENK_INDICES)
        if has_previous else FIXED_TENK_INDICES
    )
    return 0 if (
        stress["baseline_exact"] == 4068
        and stress["source_exact"] == 4072
        and stress["fixed_indices"] == FIXED_STRESS_INDICES
        and not stress["input_canonical_mismatch_indices"]
        and tenk["baseline_exact"] == 9993
        and tenk["source_exact"] == (10000 if has_previous else 9998)
        and tenk["fixed_indices"] == expected_tenk_fixed
        and not tenk["input_canonical_mismatch_indices"]
        and all(not row["regressed_indices"] for row in corpora.values())
        and all(row["normalize_changed_indices"] == row["cleanup_changed_indices"]
                and not row["reionize_changed_indices"] for row in corpora.values())
        and (not has_previous or all(
            row["previous"]["exact"] == (4072 if name == "stress" else 9998)
            and row["previous"]["fixed_indices"] == (
                [] if name == "stress" else FOLLOWUP_TENK_INDICES
            )
            and not row["previous"]["regressed_indices"]
            and row["previous"]["normalize_changed_indices"]
                == row["previous"]["cleanup_changed_indices"]
            and not row["previous"]["reionize_changed_indices"]
            for name, row in corpora.items()
        ))
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
