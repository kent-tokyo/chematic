#!/usr/bin/env python3
"""Fail closed on RDKit 2026.03.6 HBA parity for an installed source wheel.

The exposed ChEMBL corpus is development data, not the sealed A0 cohort.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors


ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "scripts/chembl_accuracy_corpus_4999.smi"
CORPUS_SHA256 = "1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e"
ORACLE_VERSION = "2026.03.6"
EXPECTED_ROWS = 5000


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wheel", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--build-profile", required=True, choices=("dev", "release"))
    args = parser.parse_args()

    if rdBase.rdkitVersion != ORACLE_VERSION:
        raise SystemExit(f"RDKit version changed: {rdBase.rdkitVersion}")
    corpus_bytes = CORPUS.read_bytes()
    corpus_sha256 = hashlib.sha256(corpus_bytes).hexdigest()
    if corpus_sha256 != CORPUS_SHA256:
        raise SystemExit(f"HBA corpus changed: {corpus_sha256}")
    smiles_rows = [line.split()[0] for line in corpus_bytes.decode("utf-8").splitlines() if line.strip()]
    if len(smiles_rows) != EXPECTED_ROWS:
        raise SystemExit(f"HBA corpus row count changed: {len(smiles_rows)}")

    counts = {"input": len(smiles_rows), "compared": 0, "exact_named": 0, "exact_native": 0,
              "parse_failed": 0}
    differences = []
    for index, smiles in enumerate(smiles_rows):
        try:
            native = chematic.from_smiles(smiles)
            reference = Chem.MolFromSmiles(smiles)
        except Exception as exc:  # retained as a failed row, never an agreement
            counts["parse_failed"] += 1
            differences.append({"index": index, "status": "parse_failed", "error": str(exc)})
            continue
        if reference is None:
            counts["parse_failed"] += 1
            differences.append({"index": index, "status": "rdkit_parse_failed"})
            continue
        expected = rdMolDescriptors.CalcNumHBA(reference)
        named = native.rdkit_hba
        ordinary = native.hba
        counts["compared"] += 1
        counts["exact_named"] += named == expected
        counts["exact_native"] += ordinary == expected
        if named != expected:
            differences.append({"index": index, "smiles": smiles,
                                "expected": expected, "named": named, "native": ordinary})

    report = {
        "schema_version": 1,
        "profile": "rdkit_2026_03_6_hba_source_wheel",
        "chematic_version": importlib.metadata.version("chematic"),
        "chematic_module": str(Path(chematic.__file__).resolve()),
        "wheel": args.wheel.name,
        "wheel_sha256": hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
        "build_profile": args.build_profile,
        "rdkit_version": rdBase.rdkitVersion,
        "corpus": {"path": str(CORPUS.relative_to(ROOT)), "sha256": corpus_sha256},
        "counts": counts,
        "differences": differences,
        "gate_passed": (counts["compared"] == EXPECTED_ROWS
                        and counts["exact_named"] == EXPECTED_ROWS
                        and counts["parse_failed"] == 0),
        "not_claimed": ["sealed-cohort accuracy", "other descriptor parity", "WASM parity"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("HBA source-wheel parity:", counts, "gate_passed=", report["gate_passed"])
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
