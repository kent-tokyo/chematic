#!/usr/bin/env python3
"""Gate stereo-sensitive canonicalization against randomized RDKit spellings.

The source corpus supplies oracle-labelled stereochemical targets.  This gate
does not reuse those labels as a canonicalization oracle: it only requires
that each RDKit-generated isomeric spelling has the same RDKit semantic
identity as its source, and that chematic produces one stable, idempotent
canonical spelling for the source plus 32 deterministic randomized spellings.
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
# RDKit's public Python RNG setter accepts an unsigned 32-bit seed.
RANDOM_SEED = 0x051E_ED26


def rdkit_identity(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit rejected SMILES: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def load_cases(path: Path) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("_manifest") is True:
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
        cwd=ROOT,
        input="\n".join(smiles) + "\n",
        text=True,
        capture_output=True,
    )
    if completed.returncode:
        raise RuntimeError(f"batch-canonicalize failed: {completed.stderr.strip()}")
    payload = json.loads(completed.stdout)
    if payload.get("status") != "complete":
        raise RuntimeError(f"batch-canonicalize status was {payload.get('status')!r}")
    records = payload.get("records")
    if not isinstance(records, list) or len(records) != len(smiles):
        raise RuntimeError("batch-canonicalize did not retain one result per input")
    return records


def canonical_from_record(record: dict[str, Any]) -> str | None:
    if record.get("error") is not None:
        return None
    value = record.get("canonical_smiles")
    return value if isinstance(value, str) else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        parser.error(f"chematic CLI not found: {cli}")

    cases = load_cases(args.corpus.resolve())
    rdBase.SeedRandomNumberGenerator(RANDOM_SEED)
    inputs: list[dict[str, Any]] = []
    for case_index, case in enumerate(cases):
        inputs.append({"case_index": case_index, "variant": "source", "smiles": case["smiles"]})
        molecule = Chem.MolFromSmiles(case["smiles"])
        assert molecule is not None  # checked by rdkit_identity above
        for variant in range(VARIANTS_PER_CASE):
            inputs.append(
                {
                    "case_index": case_index,
                    "variant": variant,
                    "smiles": Chem.MolToSmiles(
                        molecule, canonical=False, doRandom=True, isomericSmiles=True
                    ),
                }
            )

    first_pass = run_batch(cli, [entry["smiles"] for entry in inputs])
    baseline: dict[int, str] = {}
    first_pass_failures: list[dict[str, Any]] = []
    canonical_values: list[str] = []
    for entry, record in zip(inputs, first_pass, strict=True):
        canonical = canonical_from_record(record)
        case = cases[entry["case_index"]]
        if canonical is None:
            first_pass_failures.append({"case_line": case["line"], "variant": entry["variant"], "kind": "chematic_error", "detail": record.get("error")})
            continue
        canonical_values.append(canonical)
        try:
            semantic_equal = rdkit_identity(canonical) == case["rdkit_identity"]
        except ValueError as error:
            semantic_equal = False
            detail = str(error)
        else:
            detail = None
        if not semantic_equal:
            first_pass_failures.append({"case_line": case["line"], "variant": entry["variant"], "kind": "semantic_mismatch", "detail": detail})
            continue
        if entry["variant"] == "source":
            baseline[entry["case_index"]] = canonical
        elif canonical != baseline.get(entry["case_index"]):
            first_pass_failures.append({"case_line": case["line"], "variant": entry["variant"], "kind": "canonical_spelling_mismatch", "detail": canonical})

    if len(canonical_values) != len(inputs):
        second_pass_failures = [{"kind": "not_run_after_first_pass_failure"}]
    else:
        second_pass = run_batch(cli, canonical_values)
        second_pass_failures = []
        for input_index, (canonical, record) in enumerate(zip(canonical_values, second_pass, strict=True)):
            reparsed = canonical_from_record(record)
            if reparsed != canonical:
                second_pass_failures.append({"input_index": input_index, "kind": "canonical_not_idempotent", "detail": reparsed})

    duplicate_spellings = sum(
        VARIANTS_PER_CASE - len({entry["smiles"] for entry in inputs if entry["case_index"] == case_index and entry["variant"] != "source"})
        for case_index in range(len(cases))
    )
    result = {
        "schema_version": 1,
        "profile": "stereo_spelling_invariance_v1",
        "rdkit_version": rdBase.rdkitVersion,
        "cli": str(cli),
        "corpus": {
            "path": str(args.corpus.resolve().relative_to(ROOT)),
            "sha256": hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
            "cases": len(cases),
        },
        "randomized_spellings": {
            "variants_per_case": VARIANTS_PER_CASE,
            "seed": RANDOM_SEED,
            "duplicate_spellings": duplicate_spellings,
            "inputs_including_source": len(inputs),
        },
        "first_pass": {
            "records": len(first_pass),
            "failures": len(first_pass_failures),
            "failure_samples": first_pass_failures[:20],
        },
        "canonical_reparse": {
            "records": len(canonical_values),
            "failures": len(second_pass_failures),
            "failure_samples": second_pass_failures[:20],
        },
        "gate_passed": not first_pass_failures and not second_pass_failures,
        "scope": "RDKit randomized isomeric SMILES spelling invariance plus chematic canonical reparse idempotency on the fixed CIP oracle inputs",
        "not_claimed": [
            "RDKit canonical string equality",
            "independent CIP-label correctness",
            "MOL V2000 or V3000 round-trip preservation",
            "300-structure Stereo Torture Suite completion",
        ],
    }
    output = args.output if args.output.is_absolute() else ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("gate_passed", "first_pass", "canonical_reparse")}, sort_keys=True))
    return 0 if result["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
