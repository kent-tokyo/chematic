#!/usr/bin/env python3
"""Gate stereo canonicalization against 32 deterministic RDKit SMILES spellings.

This checks semantic identity and canonical reparse stability, not RDKit
canonical-string equality or independent CIP-label correctness.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from rdkit import Chem, rdBase


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = ROOT / "target" / "debug" / "chematic"
CORPUS = ROOT / "validation" / "cip_label_corpus.jsonl"
VARIANTS_PER_CASE = 32
EXPECTED_CASES = 155
RANDOM_SEED = 0x051E_ED26


def rdkit_identity(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit rejected SMILES: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def load_cases(path: Path) -> list[dict[str, Any]]:
    cases = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("_manifest"):
            continue
        smiles = row.get("smiles")
        if not isinstance(smiles, str) or not smiles:
            raise ValueError(f"corpus line {line_number} lacks a SMILES string")
        cases.append({"line": line_number, "smiles": smiles, "rdkit_identity": rdkit_identity(smiles)})
    if len(cases) != EXPECTED_CASES:
        raise ValueError(f"expected {EXPECTED_CASES} corpus rows, got {len(cases)}")
    return cases


def run_batch(cli: Path, smiles: list[str]) -> list[dict[str, Any]]:
    completed = subprocess.run(
        [str(cli), "batch-canonicalize", "--max-records", str(len(smiles))],
        cwd=ROOT, input="\n".join(smiles) + "\n", text=True, capture_output=True,
    )
    if completed.returncode:
        raise RuntimeError(f"batch-canonicalize failed: {completed.stderr.strip()}")
    payload = json.loads(completed.stdout)
    records = payload.get("records")
    if payload.get("status") != "complete" or not isinstance(records, list) or len(records) != len(smiles):
        raise RuntimeError("batch-canonicalize did not retain one complete result per input")
    return records


def canonical_from_record(record: dict[str, Any]) -> str | None:
    value = record.get("canonical_smiles")
    return value if record.get("error") is None and isinstance(value, str) else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli, corpus = args.cli.resolve(), args.corpus.resolve()
    if not cli.is_file():
        parser.error(f"chematic CLI not found: {cli}")
    cases = load_cases(corpus)
    rdBase.SeedRandomNumberGenerator(RANDOM_SEED)
    inputs: list[dict[str, Any]] = []
    for case_index, case in enumerate(cases):
        inputs.append({"case_index": case_index, "variant": "source", "smiles": case["smiles"]})
        molecule = Chem.MolFromSmiles(case["smiles"])
        assert molecule is not None
        for variant in range(VARIANTS_PER_CASE):
            inputs.append({
                "case_index": case_index, "variant": variant,
                "smiles": Chem.MolToSmiles(molecule, canonical=False, doRandom=True, isomericSmiles=True),
            })

    first_pass = run_batch(cli, [entry["smiles"] for entry in inputs])
    baseline: dict[int, str] = {}
    canonical_values: list[str] = []
    failures: list[dict[str, Any]] = []
    for entry, record in zip(inputs, first_pass, strict=True):
        canonical = canonical_from_record(record)
        case = cases[entry["case_index"]]
        if canonical is None:
            failures.append({"case_line": case["line"], "variant": entry["variant"], "kind": "chematic_error", "detail": record.get("error")})
            continue
        canonical_values.append(canonical)
        try:
            semantic_equal = rdkit_identity(canonical) == case["rdkit_identity"]
        except ValueError as error:
            semantic_equal, detail = False, str(error)
        else:
            detail = None
        if not semantic_equal:
            failures.append({"case_line": case["line"], "variant": entry["variant"], "kind": "semantic_mismatch", "detail": detail})
        elif entry["variant"] == "source":
            baseline[entry["case_index"]] = canonical
        elif canonical != baseline.get(entry["case_index"]):
            failures.append({"case_line": case["line"], "variant": entry["variant"], "kind": "canonical_spelling_mismatch", "detail": canonical})

    reparse_failures: list[dict[str, Any]] = []
    if len(canonical_values) == len(inputs):
        for index, (canonical, record) in enumerate(zip(canonical_values, run_batch(cli, canonical_values), strict=True)):
            if canonical_from_record(record) != canonical:
                reparse_failures.append({"input_index": index, "kind": "canonical_not_idempotent", "detail": canonical_from_record(record)})
    else:
        reparse_failures.append({"kind": "not_run_after_first_pass_failure"})

    duplicate_spellings = sum(
        VARIANTS_PER_CASE - len({entry["smiles"] for entry in inputs if entry["case_index"] == index and entry["variant"] != "source"})
        for index in range(len(cases))
    )
    result = {
        "schema_version": 1, "profile": "stereo_spelling_invariance_v1", "rdkit_version": rdBase.rdkitVersion,
        "cli": cli.name,
        "corpus": {"path": corpus.relative_to(ROOT).as_posix(), "sha256": hashlib.sha256(corpus.read_bytes()).hexdigest(), "cases": len(cases)},
        "randomized_spellings": {"variants_per_case": VARIANTS_PER_CASE, "seed": RANDOM_SEED, "duplicate_spellings": duplicate_spellings, "inputs_including_source": len(inputs)},
        "first_pass": {"records": len(first_pass), "failures": len(failures), "failure_samples": failures[:20]},
        "canonical_reparse": {"records": len(canonical_values), "failures": len(reparse_failures), "failure_samples": reparse_failures[:20]},
        "gate_passed": not failures and not reparse_failures,
        "scope": "RDKit randomized isomeric SMILES spelling invariance plus chematic canonical reparse idempotency on fixed CIP oracle inputs",
        "not_claimed": ["RDKit canonical string equality", "independent CIP-label correctness", "MOL V2000 or V3000 round-trip preservation", "300-structure Stereo Torture Suite completion"],
    }
    output = args.output if args.output.is_absolute() else ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("gate_passed", "first_pass", "canonical_reparse")}, sort_keys=True))
    return 0 if result["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
