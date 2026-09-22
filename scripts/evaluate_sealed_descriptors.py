#!/usr/bin/env python3
"""Evaluate a frozen candidate's eight-field descriptor profile once.

Raw input and per-row output stay outside the repository.  The commit-safe
summary contains aggregate accounting and cryptographic provenance only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


FIELD_TOLERANCES = {
    "molecular_weight": 1e-6,
    "hba": 0.0,
    "hbd": 0.0,
    "tpsa": 1e-6,
    "logp": 1e-6,
    "molar_refractivity": 1e-6,
    "fsp3": 1e-6,
    "aromatic_ring_count": 0.0,
}
INTEGER_FIELDS = {"hba", "hbd", "aromatic_ring_count"}


def sha256(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def resolve_commit(worktree: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(worktree), "rev-parse", "HEAD^{commit}"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def read_smiles(path: Path) -> list[str]:
    return [line.split()[0] for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def descriptor_values(rd_mol: Any) -> dict[str, float | int]:
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    return {
        "molecular_weight": Descriptors.MolWt(rd_mol),
        "hba": Lipinski.NumHAcceptors(rd_mol),
        "hbd": Lipinski.NumHDonors(rd_mol),
        "tpsa": rdMolDescriptors.CalcTPSA(rd_mol, includeSandP=True),
        "logp": Crippen.MolLogP(rd_mol),
        "molar_refractivity": Crippen.MolMR(rd_mol),
        "fsp3": rdMolDescriptors.CalcFractionCSP3(rd_mol),
        "aromatic_ring_count": rdMolDescriptors.CalcNumAromaticRings(rd_mol),
    }


def compare_fields(
    actual: object, expected: dict[str, float | int]
) -> dict[str, dict[str, object]]:
    if not isinstance(actual, dict):
        actual = {}
    compared: dict[str, dict[str, object]] = {}
    for field, tolerance in FIELD_TOLERANCES.items():
        observed = actual.get(field)
        supported = (
            isinstance(observed, (int, float))
            and not isinstance(observed, bool)
            and math.isfinite(float(observed))
        )
        if not supported:
            compared[field] = {
                "status": "unsupported",
                "absolute_error": None,
                "strict_match": False,
            }
            continue
        expected_value = expected[field]
        if field in INTEGER_FIELDS and not float(observed).is_integer():
            compared[field] = {
                "status": "mismatch",
                "absolute_error": abs(float(observed) - float(expected_value)),
                "strict_match": False,
            }
            continue
        error = abs(float(observed) - float(expected_value))
        matched = error <= tolerance
        compared[field] = {
            "status": "match" if matched else "mismatch",
            "absolute_error": error,
            "strict_match": matched,
        }
    return compared


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--split", choices=("development", "sealed_holdout"), required=True)
    parser.add_argument("--candidate-worktree", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--raw-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    args = parser.parse_args()

    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if args.split == "sealed_holdout" and not (
        manifest.get("status") == "sealed" and manifest.get("sealed") is True
    ):
        parser.error("sealed_holdout evaluation requires a sealed cohort manifest")
    freeze = manifest.get("candidate_freeze") or {}
    candidate_commit = freeze.get("candidate_commit")
    candidate_tag = freeze.get("candidate_tag")
    worktree = args.candidate_worktree.resolve()
    if resolve_commit(worktree) != candidate_commit:
        parser.error("candidate worktree HEAD does not match the manifest freeze")
    binary = args.binary.resolve()
    if not binary.is_file() or not binary.is_relative_to(worktree):
        parser.error("evaluation binary must exist inside the frozen candidate worktree")

    split = manifest["splits"][args.split]
    split_path = Path(split["path"]).resolve()
    if sha256(split_path) != split["sha256"]:
        parser.error("cohort split digest does not match the manifest")
    smiles = read_smiles(split_path)
    if len(smiles) != split["rows"]:
        parser.error("cohort split row count does not match the manifest")

    completed = subprocess.run(
        [str(binary)],
        input="\n".join(smiles) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    records = [json.loads(line) for line in completed.stdout.splitlines() if line.strip()]
    if len(records) != len(smiles):
        parser.error("candidate binary output row count does not match the cohort")

    import rdkit
    from rdkit import Chem

    field_results = {
        field: {"strict_matches": 0, "mismatches": 0, "unsupported": 0, "max_abs_error": 0.0}
        for field in FIELD_TOLERANCES
    }
    parse_failures = 0
    oracle_failures = 0
    raw_path = args.raw_output.resolve()
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    with raw_path.open("w", encoding="utf-8") as raw:
        for index, (source, record) in enumerate(zip(smiles, records, strict=True)):
            rd_mol = Chem.MolFromSmiles(source)
            if rd_mol is None:
                oracle_failures += 1
                raw.write(json.dumps({"index": index, "status": "oracle_parse_failure"}) + "\n")
                continue
            if not isinstance(record, dict) or record.get("index") != index or record.get("status") != "ok":
                parse_failures += 1
                raw.write(json.dumps({"index": index, "status": "candidate_parse_failure"}) + "\n")
                continue
            compared = compare_fields(record.get("descriptors"), descriptor_values(rd_mol))
            for field, result in compared.items():
                totals = field_results[field]
                if result["status"] == "unsupported":
                    totals["unsupported"] += 1
                elif result["strict_match"]:
                    totals["strict_matches"] += 1
                else:
                    totals["mismatches"] += 1
                error = result["absolute_error"]
                if isinstance(error, float):
                    totals["max_abs_error"] = max(float(totals["max_abs_error"]), error)
            raw.write(json.dumps({"index": index, "status": "evaluated", "fields": compared}, separators=(",", ":")) + "\n")

    rows = len(smiles)
    all_fields_pass = all(
        result["strict_matches"] == rows
        and result["mismatches"] == 0
        and result["unsupported"] == 0
        for result in field_results.values()
    )
    accepted = parse_failures == 0 and oracle_failures == 0 and all_fields_pass
    attestation = manifest.get("unused_data_attestation") or {}
    summary = {
        "schema_version": 1,
        "profile": "rdkit_core_eight_frozen_candidate_evaluation_v1",
        "status": (
            "accepted" if accepted else "rejected"
        ) if args.split == "sealed_holdout" else ("development_passed" if accepted else "development_failed"),
        "scope": (
            "One-time sealed eight-field descriptor evaluation; raw input and per-row output remain local-only and exposed."
            if args.split == "sealed_holdout"
            else "Development split rehearsal; not sealed evidence."
        ),
        "evaluated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "candidate": {"commit": candidate_commit, "tag": candidate_tag},
        "oracle": {
            "engine": "RDKit",
            "version": rdkit.__version__,
            "profile": "rdkit_compat_v2",
        },
        "cohort": {
            "split": args.split,
            "rows": rows,
            "sha256": split["sha256"],
            "manifest_sha256": sha256(manifest_path),
            "source_sha256": manifest["source"]["source_sha256"],
        },
        "provenance": {
            "candidate_binary_sha256": sha256(binary),
            "unused_data_attestation_sha256": attestation.get("sha256"),
            "raw_result_sha256": sha256(raw_path),
        },
        "accounting": {
            "rows": rows,
            "parsed": rows - parse_failures - oracle_failures,
            "candidate_parse_failures": parse_failures,
            "oracle_parse_failures": oracle_failures,
            "all_rows_accounted": True,
        },
        "strict_tolerances": FIELD_TOLERANCES,
        "fields": field_results,
        "decision": {
            "accepted": accepted,
            "requirement": "every declared field must strictly match on every row with no parse failure or unsupported value",
        },
    }
    summary_path = args.summary_output.resolve()
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "rows": rows, "fields": field_results}, sort_keys=True))
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
