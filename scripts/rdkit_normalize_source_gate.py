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
    baseline_path: Path, source_path: Path,
) -> dict:
    corpus_path = ROOT / relative_path
    smiles_rows = [
        line.split()[0] for line in corpus_path.read_text(encoding="utf-8").splitlines()
    ]
    if len(smiles_rows) != expected_rows:
        raise ValueError(f"{relative_path}: unexpected corpus length")
    baseline = read_dump(baseline_path, smiles_rows)
    source = read_dump(source_path, smiles_rows)
    old_mismatches, new_mismatches, input_mismatches = [], [], []
    changed, cleanup_changed, reionize_changed = [], [], []
    fixed_cases = []
    fixed_indices = FIXED_STRESS_INDICES if name == "stress" else FIXED_TENK_INDICES
    for index, smiles in enumerate(smiles_rows):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"RDKit cannot parse {name} row {index}")
        rdkit_input = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        if baseline[index].get("rdkit_input") != source[index].get("rdkit_input"):
            raise ValueError(f"{name} row {index}: source input changed")
        if ok_value(source[index], "rdkit_input") != rdkit_input:
            input_mismatches.append(index)
        expected = Chem.MolToSmiles(
            rdMolStandardize.Normalize(mol), canonical=True, isomericSmiles=True
        )
        old = ok_value(baseline[index], "normalize")
        new = ok_value(source[index], "normalize")
        if old != expected:
            old_mismatches.append(index)
        if new != expected:
            new_mismatches.append(index)
        if old != new:
            changed.append(index)
        if baseline[index].get("cleanup") != source[index].get("cleanup"):
            cleanup_changed.append(index)
        if baseline[index].get("reionize") != source[index].get("reionize"):
            reionize_changed.append(index)
        if index in fixed_indices:
            fixed_cases.append({
                "index": index, "input": smiles, "baseline": old,
                "source": new, "rdkit": expected,
            })
    old_set, new_set = set(old_mismatches), set(new_mismatches)
    return {
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-stress-jsonl", type=Path, required=True)
    parser.add_argument("--source-stress-jsonl", type=Path, required=True)
    parser.add_argument("--baseline-tenk-jsonl", type=Path, required=True)
    parser.add_argument("--source-tenk-jsonl", type=Path, required=True)
    parser.add_argument("--baseline-source-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1; found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")
    subprocess.run(
        ["git", "rev-parse", "--verify", args.baseline_source_commit + "^{commit}"],
        cwd=ROOT, check=True, stdout=subprocess.DEVNULL,
    )
    subprocess.run(
        ["git", "diff", "--quiet", "HEAD", "--", *SOURCE_FILES], cwd=ROOT,
        check=True,
    )
    dump_paths = (
        (args.baseline_stress_jsonl, args.source_stress_jsonl),
        (args.baseline_tenk_jsonl, args.source_tenk_jsonl),
    )
    corpora = {
        name: score_corpus(name, path, size, baseline, source)
        for (name, path, size), (baseline, source) in zip(CORPORA, dump_paths)
    }
    report = {
        "schema_version": 1,
        "profile": "source_rdkit_2026_09_1_normalize",
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "baseline_source_commit": args.baseline_source_commit,
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
    return 0 if (
        stress["baseline_exact"] == 4068
        and stress["source_exact"] == 4072
        and stress["fixed_indices"] == FIXED_STRESS_INDICES
        and not stress["input_canonical_mismatch_indices"]
        and tenk["baseline_exact"] == 9993
        and tenk["source_exact"] == 9998
        and tenk["fixed_indices"] == FIXED_TENK_INDICES
        and not tenk["input_canonical_mismatch_indices"]
        and all(not row["regressed_indices"] for row in corpora.values())
        and all(row["normalize_changed_indices"] == row["cleanup_changed_indices"]
                and not row["reionize_changed_indices"] for row in corpora.values())
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
