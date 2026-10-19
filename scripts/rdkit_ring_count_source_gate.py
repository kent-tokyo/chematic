#!/usr/bin/env python3
"""Compare source-built ring counts with pinned RDKit on two fixed corpora."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
RDKIT_VERSION = "2026.09.1"
CORPORA = (
    ("stress", "validation/cosmolkit_comparison/corpora/stress-v1.smi", 4072),
    ("exposed_10k", "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi", 10000),
)
FIXED_STRESS_INDICES = (605, 608, 1014, 3774)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_dump(path: Path, corpus: list[str]) -> list[dict]:
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    if len(rows) != len(corpus):
        raise ValueError(f"{path}: expected {len(corpus)} rows, found {len(rows)}")
    for index, (row, smiles) in enumerate(zip(rows, corpus, strict=True)):
        if row["index"] != index or row["smiles"] != smiles:
            raise ValueError(f"{path}: source row {index} differs from corpus")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stress-source-jsonl", type=Path, required=True)
    parser.add_argument("--tenk-source-jsonl", type=Path, required=True)
    parser.add_argument("--baseline-stress-jsonl", type=Path, required=True)
    parser.add_argument("--baseline-source-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != RDKIT_VERSION:
        parser.error(f"expected RDKit {RDKIT_VERSION}; found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")

    source_paths = (args.stress_source_jsonl, args.tenk_source_jsonl)
    corpora = {}
    for (name, relative_path, expected_rows), source_path in zip(
        CORPORA, source_paths, strict=True
    ):
        corpus_path = ROOT / relative_path
        smiles_rows = [
            line.split()[0] for line in corpus_path.read_text(encoding="utf-8").splitlines()
        ]
        if len(smiles_rows) != expected_rows:
            raise ValueError(f"{relative_path}: unexpected corpus size")
        source_rows = read_dump(source_path, smiles_rows)
        mismatches = []
        oracle = []
        for index, (smiles, source_row) in enumerate(
            zip(smiles_rows, source_rows, strict=True)
        ):
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                raise ValueError(f"RDKit cannot parse {name} row {index}")
            expected = rdMolDescriptors.CalcNumRings(molecule)
            oracle.append(expected)
            if source_row["num_rings"] != expected:
                mismatches.append({
                    "index": index,
                    "smiles": smiles,
                    "source": source_row["num_rings"],
                    "rdkit": expected,
                })
        corpora[name] = {
            "path": relative_path,
            "corpus_sha256": sha256(corpus_path),
            "source_dump_sha256": sha256(source_path),
            "rows": len(source_rows),
            "exact": len(source_rows) - len(mismatches),
            "mismatches": mismatches,
        }
        if name == "stress":
            baseline_rows = read_dump(args.baseline_stress_jsonl, smiles_rows)
            baseline_mismatches = [
                index for index, (row, expected) in enumerate(
                    zip(baseline_rows, oracle, strict=True)
                ) if row["num_rings"] != expected
            ]
            corpora[name]["baseline"] = {
                "source_commit": args.baseline_source_commit,
                "source_dump_sha256": sha256(args.baseline_stress_jsonl),
                "exact": len(baseline_rows) - len(baseline_mismatches),
                "mismatch_indices": baseline_mismatches,
            }
            corpora[name]["fixed_cases"] = [
                {
                    "index": index,
                    "smiles": smiles_rows[index],
                    "baseline": baseline_rows[index]["num_rings"],
                    "source": source_rows[index]["num_rings"],
                    "rdkit": oracle[index],
                }
                for index in FIXED_STRESS_INDICES
            ]

    source_files = [
        "crates/chematic-chem/src/descriptors/rings_and_filters.rs",
        "crates/chematic-chem/src/descriptors/ring_component_count_tests.rs",
        "crates/chematic-chem/src/descriptors_tests.rs",
        "crates/chematic-chem/examples/ring_count_stress_dump.rs",
        "scripts/rdkit_ring_count_source_gate.py",
    ]
    report = {
        "schema_version": 1,
        "profile": "source_rdkit_2026_09_1_num_rings",
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "rdkit_version": rdBase.rdkitVersion,
        "source_file_sha256": {
            path: sha256(ROOT / path) for path in source_files
        },
        "corpora": corpora,
        "limits": [
            "Counts compare two fixed SMILES corpora with RDKit 2026.09.1 CalcNumRings.",
            "Exact agreement on these corpora does not prove agreement for every molecular graph.",
            "The source-built Rust example, not a published package, produced CheMatic rows.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print({name: f"{row['exact']}/{row['rows']}" for name, row in corpora.items()})
    return 0 if (
        not any(row["mismatches"] for row in corpora.values())
        and corpora["stress"]["baseline"]["mismatch_indices"]
        == list(FIXED_STRESS_INDICES)
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
