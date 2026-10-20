#!/usr/bin/env python3
"""Compare schematic and RDKit CDXML parsing in one Python process.

This is a bounded semantic contract.  The fixture is a page/fragment CDXML
document accepted by both parsers; timing is recorded only as non-ranking
context because the two APIs expose different document abstractions.
"""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib.same_process_report import make_report, write_report


ROOT = Path(__file__).resolve().parents[2]
MALFORMED = {
    "unknown_element": '<CDXML><page><fragment><n id="1" p="0 0" Element="999"/></fragment></page></CDXML>',
    "truncated_fragment": "<CDXML><page><fragment><n id=\"1\"",
}


def schematic(text: str) -> tuple[int, tuple[int, int, str] | None]:
    import chematic
    from rdkit import Chem

    try:
        molecule = chematic.from_cdxml(text)
        parsed = Chem.MolFromSmiles(molecule.smiles)
        canonical = Chem.MolToSmiles(parsed, canonical=True) if parsed is not None else molecule.smiles
        return 1, (len(molecule.atom_table), molecule.heavy_atoms, canonical)
    except Exception:
        return 0, None


def rdkit(text: str) -> tuple[int, tuple[int, int, str] | None]:
    from rdkit import Chem

    try:
        molecules = Chem.MolsFromCDXML(text.encode("utf-8"))
    except Exception:
        return 0, None
    if not molecules:
        return 0, None
    molecule = molecules[0]
    if molecule.GetNumConformers() == 0:
        return 0, None
    conformer = molecule.GetConformer()
    finite = all(
        math.isfinite(conformer.GetAtomPosition(i).x)
        and math.isfinite(conformer.GetAtomPosition(i).y)
        and math.isfinite(conformer.GetAtomPosition(i).z)
        for i in range(molecule.GetNumAtoms())
    )
    if not finite:
        return 0, None
    return 1, (
        molecule.GetNumAtoms(),
        molecule.GetNumHeavyAtoms(),
        Chem.MolToSmiles(molecule, canonical=True),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cdxml", type=Path, default=Path("benchmarks/fixtures/ethanol.cdxml"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = args.cdxml if args.cdxml.is_absolute() else ROOT / args.cdxml
    payload = path.read_bytes()
    text = payload.decode("utf-8")

    started = time.perf_counter()
    schematic_rows = [schematic(text) for _ in range(args.repeats)]
    schematic_seconds = time.perf_counter() - started
    started = time.perf_counter()
    rdkit_rows = [rdkit(text) for _ in range(args.repeats)]
    rdkit_seconds = time.perf_counter() - started

    schematic_records = sum(row[0] for row in schematic_rows)
    rdkit_records = sum(row[0] for row in rdkit_rows)
    signature_mismatches = sum(left[1] != right[1] for left, right in zip(schematic_rows, rdkit_rows))
    errors = []
    if schematic_records != args.repeats:
        errors.append("schematic record count mismatch")
    if rdkit_records != args.repeats:
        errors.append("RDKit record count mismatch")
    if signature_mismatches:
        errors.append(f"CDXML signature mismatch in {signature_mismatches} repetitions")

    malformed = {}
    malformed_rejection_mismatches = []
    for case_id, case in MALFORMED.items():
        left, right = schematic(case), rdkit(case)
        malformed[case_id] = {
            "chematic": {"records": left[0], "failures": 1 - left[0]},
            "rdkit": {"records": right[0], "failures": 1 - right[0]},
        }
        if left[0] != right[0]:
            malformed_rejection_mismatches.append(case_id)

    report = make_report(
        ROOT, path, payload, args.repeats,
        "same_process_cdxml_semantic_contract",
        rows={
            "chematic": {"records": schematic_records, "failures": args.repeats - schematic_records, "seconds": round(schematic_seconds, 6), "boundary": "chematic.from_cdxml in the current process"},
            "rdkit": {"records": rdkit_records, "failures": args.repeats - rdkit_records, "seconds": round(rdkit_seconds, 6), "boundary": "RDKit MolsFromCDXML in the current process"},
        },
        comparison={
            "expected_records": args.repeats,
            "signature": "total atoms, heavy atoms, canonical SMILES; RDKit coordinates finite",
            "signature_mismatch_repetitions": signature_mismatches,
            "malformed_cases": malformed,
            "malformed_rejection_mismatches": malformed_rejection_mismatches,
            "malformed_contract": "reported only; schematic CDXML parser is intentionally lenient",
            "speed_is_non_ranking_context": True,
        },
        errors=errors,
    )
    return write_report(ROOT, args.output, report)


if __name__ == "__main__":
    raise SystemExit(main())
