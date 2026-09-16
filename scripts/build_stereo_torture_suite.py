#!/usr/bin/env python3
"""Build the checked-in *development* stereo torture suite.

The suite intentionally combines already-exposed regression structures with a
fixed, RDKit-filtered sample from the committed descriptor census corpus. It
is not a sealed or independent accuracy evaluation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from rdkit import Chem, rdBase


ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "validation" / "stereo_torture_suite_development.jsonl"
SAMPLE_SOURCE = ROOT / "scripts" / "descriptor_census_corpus.smi"
SOURCES = (
    ROOT / "validation" / "cip_label_corpus.jsonl",
    ROOT / "validation" / "cip_residual_classification_corpus.jsonl",
    ROOT / "validation" / "cip_rule4b_discriminating_corpus.jsonl",
    ROOT / "validation" / "cip_negative_resonance_order_invariance.jsonl",
    ROOT / "validation" / "cip_oracle_instability.jsonl",
    ROOT / "validation" / "canonical_residual_fixtures.jsonl",
)
TARGET_CASES = 300
SEED = "chematic-stereo-torture-development-v1"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows_from_jsonl(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        item = json.loads(line)
        if item.get("_manifest"):
            continue
        smiles = item.get("smiles") or item.get("input")
        if isinstance(smiles, str) and smiles:
            rows.append({"smiles": smiles, "source": path.relative_to(ROOT).as_posix(), "metadata": item})
    return rows


def known_cases() -> list[dict[str, Any]]:
    by_smiles: dict[str, dict[str, Any]] = {}
    for path in SOURCES:
        for row in rows_from_jsonl(path):
            existing = by_smiles.get(row["smiles"])
            if existing is None:
                row["source_paths"] = [row.pop("source")]
                by_smiles[row["smiles"]] = row
            else:
                existing["source_paths"].append(row["source"])
    return [by_smiles[key] for key in sorted(by_smiles)]


def stereo_types(smiles: str) -> list[str]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return []
    return sorted({str(info.type) for info in Chem.FindPotentialStereo(molecule)})


def sampled_cases(excluded: set[str], count: int) -> list[dict[str, Any]]:
    candidates = []
    for line_number, line in enumerate(SAMPLE_SOURCE.read_text(encoding="utf-8").splitlines(), start=1):
        smiles = line.strip()
        if not smiles or smiles in excluded or not stereo_types(smiles):
            continue
        key = hashlib.sha256(f"{SEED}\0{smiles}".encode()).hexdigest()
        candidates.append((key, line_number, smiles, stereo_types(smiles)))
    candidates.sort()
    if len(candidates) < count:
        raise ValueError(f"need {count} stereo candidates, found {len(candidates)}")
    return [
        {
            "smiles": smiles,
            "source_paths": [SAMPLE_SOURCE.relative_to(ROOT).as_posix()],
            "sample_line": line_number,
            "rdkit_stereo_types": types,
        }
        for _, line_number, smiles, types in candidates[:count]
    ]


def build() -> list[dict[str, Any]]:
    known = known_cases()
    sampled = sampled_cases({row["smiles"] for row in known}, TARGET_CASES - len(known))
    cases = []
    for index, row in enumerate(known + sampled, start=1):
        metadata = row.pop("metadata", {})
        cases.append({
            "id": f"stereo-dev-{index:03d}",
            "smiles": row["smiles"],
            "tier": "existing_regression" if metadata else "fixed_sampled_development",
            "source_paths": row["source_paths"],
            "source_bucket": metadata.get("bucket") or metadata.get("category") or metadata.get("case"),
            "rdkit_stereo_types": row.get("rdkit_stereo_types") or stereo_types(row["smiles"]),
            "sample_line": row.get("sample_line"),
        })
    if len(cases) != TARGET_CASES or len({row["smiles"] for row in cases}) != TARGET_CASES:
        raise ValueError("suite must contain 300 unique structures")
    return cases


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=TARGET)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    cases = build()
    type_counts = Counter(kind for row in cases for kind in row["rdkit_stereo_types"])
    manifest = {
        "_manifest": True,
        "schema_version": 1,
        "suite": "stereo_torture_development_v1",
        "scope": "development regression only; not sealed, independent, or an accuracy claim",
        "rdkit_version": rdBase.rdkitVersion,
        "seed": SEED,
        "unique_structures": len(cases),
        "existing_regression": sum(row["tier"] == "existing_regression" for row in cases),
        "fixed_sampled_development": sum(row["tier"] == "fixed_sampled_development" for row in cases),
        "rdkit_stereo_type_counts": dict(sorted(type_counts.items())),
        "source_sha256": {path.relative_to(ROOT).as_posix(): digest(path) for path in (*SOURCES, SAMPLE_SOURCE)},
    }
    text = "\n".join(json.dumps(row, sort_keys=True) for row in [manifest, *cases]) + "\n"
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if args.check:
        if not output.is_file() or output.read_text(encoding="utf-8") != text:
            raise SystemExit("suite differs; rerun without --check")
    else:
        output.write_text(text, encoding="utf-8")
    print(json.dumps({"output": str(output), **manifest}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
