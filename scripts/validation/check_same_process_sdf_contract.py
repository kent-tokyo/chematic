#!/usr/bin/env python3
"""Compare chematic and RDKit SDF parsing in one Python process.

This is a semantic contract gate, not a throughput ranking. Both libraries
receive the same bytes in the same process. Parser-specific failure behavior
is retained as counts, while canonical SMILES provides a small structural
comparison for every accepted record.
"""

from __future__ import annotations

import argparse
import io
import json
import tempfile
import time
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib.same_process_report import make_report, write_report


ROOT = Path(__file__).resolve().parents[2]
MALFORMED_CASES = {
    "bad_atom_line": "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\nthis is not an atom line\nM  END\n$$$$\n",
    "bad_coordinate": "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n  nope 0.0 0.0 C  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n",
}


def chematic_records(path: Path) -> tuple[int, int, list[str]]:
    import chematic

    records = 0
    failures = 0
    smiles: list[str] = []
    try:
        stream = chematic.iter_sdf_batched(str(path), batch_size=100)
        for batch in stream:
            records += len(batch)
            smiles.extend(record.smiles for record in batch)
        manifest = json.loads(stream.manifest_json())
        failures = int(manifest["rejected_records"])
    except Exception:
        return records, failures + 1, smiles
    return records, failures, smiles


def rdkit_records(payload: bytes) -> tuple[int, int, list[str]]:
    from rdkit import Chem

    records = 0
    failures = 0
    smiles: list[str] = []
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(payload), sanitize=True, removeHs=False)
    for molecule in supplier:
        if molecule is None:
            failures += 1
            continue
        records += 1
        smiles.append(Chem.MolToSmiles(molecule, canonical=True))
    return records, failures, smiles


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, default=Path("benchmarks/fixtures/streaming.sdf"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = args.sdf if args.sdf.is_absolute() else ROOT / args.sdf
    payload = path.read_bytes()
    text = payload.decode("utf-8")

    with tempfile.TemporaryDirectory(prefix="chematic-same-process-sdf-") as directory:
        input_path = Path(directory) / "input.sdf"
        input_path.write_bytes(payload)
        started = time.perf_counter()
        chematic = [chematic_records(input_path) for _ in range(args.repeats)]
        chematic_seconds = time.perf_counter() - started
        started = time.perf_counter()
        rdkit = [rdkit_records(payload) for _ in range(args.repeats)]
        rdkit_seconds = time.perf_counter() - started

    chematic_records_total = sum(row[0] for row in chematic)
    chematic_failures = sum(row[1] for row in chematic)
    rdkit_records_total = sum(row[0] for row in rdkit)
    rdkit_failures = sum(row[1] for row in rdkit)
    smiles_mismatches = sum(
        1
        for chematic_row, rdkit_row in zip(chematic, rdkit)
        if chematic_row[2] != rdkit_row[2]
    )
    expected_records = 2 * args.repeats
    errors = []
    if chematic_records_total != expected_records or chematic_failures != 0:
        errors.append("chematic record/failure count mismatch")
    if rdkit_records_total != expected_records or rdkit_failures != 0:
        errors.append("RDKit record/failure count mismatch")
    if smiles_mismatches:
        errors.append(f"canonical SMILES mismatch in {smiles_mismatches} repetitions")

    malformed = {}
    malformed_failure_count_mismatches = []
    with tempfile.TemporaryDirectory(prefix="chematic-same-process-sdf-malformed-") as directory:
        for case_id, case in MALFORMED_CASES.items():
            case_path = Path(directory) / f"{case_id}.sdf"
            case_path.write_text(case, encoding="utf-8")
            chematic_result = chematic_records(case_path)
            rdkit_result = rdkit_records(case.encode("utf-8"))
            malformed[case_id] = {
                "chematic": {"records": chematic_result[0], "failures": chematic_result[1]},
                "rdkit": {"records": rdkit_result[0], "failures": rdkit_result[1]},
            }
            for engine, result in (("chematic", chematic_result), ("rdkit", rdkit_result)):
                if result[0] != 0:
                    errors.append(f"{engine} malformed case {case_id} produced a record")
            if chematic_result[1] != rdkit_result[1]:
                malformed_failure_count_mismatches.append(case_id)

    report = make_report(
        ROOT, path, payload, args.repeats,
        "same_process_sdf_semantic_contract",
        rows={
            "chematic": {
                "records": chematic_records_total,
                "failures": chematic_failures,
                "seconds": round(chematic_seconds, 6),
                "boundary": "schematic Python extension in the current process",
            },
            "rdkit": {
                "records": rdkit_records_total,
                "failures": rdkit_failures,
                "seconds": round(rdkit_seconds, 6),
                "boundary": "RDKit ForwardSDMolSupplier in the current process",
            },
        },
        comparison={
            "expected_records": expected_records,
            "canonical_smiles_mismatch_repetitions": smiles_mismatches,
            "malformed_cases": malformed,
            "malformed_failure_count_mismatches": malformed_failure_count_mismatches,
            "malformed_rejection_contract": "zero accepted records; failure-count semantics are reported, not normalized",
            "speed_is_non_ranking_context": True,
        },
        errors=errors,
    )
    return write_report(ROOT, args.output, report)


if __name__ == "__main__":
    raise SystemExit(main())
