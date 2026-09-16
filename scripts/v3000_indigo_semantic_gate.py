#!/usr/bin/env python3
"""Check bounded Indigo -> chematic -> RDKit V3000 interoperability.

The gate deliberately covers only ordinary molecular V3000 records generated
by a pinned Indigo Python binding.  It checks chematic's output through both
RDKit and Indigo, while keeping SGROUP/COLLECTION semantics, ENDPTS/ATTACH
interpretation, haptic chemistry, query, and polymer semantics out of scope.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

from indigo import Indigo
from rdkit import Chem

from v3000_rdkit_semantic_gate import CASES, cli_provenance, signature


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = ROOT / "target" / "debug" / "chematic"
EXPECTED_COMPARABLE_CASES = 18
EXPECTED_SOURCE_CONTRACT_EXCLUSIONS = {"alanine", "lactic_acid"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        parser.error(f"chematic CLI not found: {cli}")

    indigo = Indigo()
    indigo.setOption("molfile-saving-mode", "3000")
    indigo_version = indigo.version()
    rows = []
    with tempfile.TemporaryDirectory(prefix="chematic-v3000-indigo-") as raw_tmp:
        tmp = Path(raw_tmp)
        for case_id, smiles in CASES.items():
            source_molecule = indigo.loadMolecule(smiles)
            source = source_molecule.molfile()
            row = {"id": case_id, "smiles": smiles}
            # Indigo 1.46.0 can emit a flat V3000 record with atom CFG for a
            # tetrahedral input which it then rejects itself. Such a record is
            # not a valid Indigo-to-other-engine interoperability input, so
            # retain it as an explicit source-contract exclusion rather than
            # attributing the rejection to chematic.
            try:
                source_reopened = indigo.loadMolecule(source)
            except Exception as exc:  # see comment above about exception type
                row.update(
                    {
                        "source_indigo_reopen": "rejected",
                        "source_indigo_error": str(exc),
                        "excluded_by_source_contract": True,
                    }
                )
                rows.append(row)
                continue
            source_path = tmp / f"{case_id}.indigo.v3000"
            output_path = tmp / f"{case_id}.chematic.v3000"
            source_path.write_text(source, encoding="utf-8")
            run = subprocess.run(
                [
                    str(cli), "convert", "--input-format", "mol_v3000", "--output-format",
                    "mol_v3000", "--input", str(source_path), "--output", str(output_path),
                ],
                text=True,
                capture_output=True,
            )
            row["cli_exit_code"] = run.returncode
            if run.returncode:
                row["error"] = run.stderr.strip()
                rows.append(row)
                continue
            written = output_path.read_text(encoding="utf-8")
            rdkit_reopened = Chem.MolFromMolBlock(written, sanitize=True, removeHs=False)
            if rdkit_reopened is None:
                row["error"] = "RDKit rejected chematic V3000 output"
                rows.append(row)
                continue
            try:
                indigo_reopened = indigo.loadMolecule(written)
            except Exception as exc:  # Indigo exceptions do not share a stable public base class.
                row["error"] = f"Indigo rejected chematic V3000 output: {exc}"
                rows.append(row)
                continue

            # Compare against the actual Indigo V3000 payload, not the
            # original SMILES. Indigo may deliberately omit information it
            # cannot represent in this export (for example E/Z in a flat,
            # bond-direction-free V3000 record); treating that omission as a
            # chematic regression would make this an invalid gate.
            expected = Chem.MolFromMolBlock(source, sanitize=True, removeHs=False)
            if expected is None:
                raise RuntimeError(f"RDKit could not parse Indigo V3000 fixture {case_id}: {smiles}")
            before = signature(expected)
            after = signature(rdkit_reopened)
            indigo_before = source_reopened.canonicalSmiles()
            indigo_after = indigo_reopened.canonicalSmiles()
            row.update(
                {
                    "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
                    "output_sha256": hashlib.sha256(written.encode()).hexdigest(),
                    "rdkit_before": before,
                    "rdkit_after": after,
                    "rdkit_semantic_equal": before == after,
                    "indigo_canonical_before": indigo_before,
                    "indigo_canonical_after": indigo_after,
                    "indigo_semantic_equal": indigo_before == indigo_after,
                }
            )
            rows.append(row)

    comparable_rows = [row for row in rows if not row.get("excluded_by_source_contract")]
    failures = [
        row["id"]
        for row in comparable_rows
        if not row.get("rdkit_semantic_equal", False)
        or not row.get("indigo_semantic_equal", False)
    ]
    exclusion_ids = {
        row["id"] for row in rows if row.get("excluded_by_source_contract")
    }
    gate_passed = (
        not failures
        and len(comparable_rows) == EXPECTED_COMPARABLE_CASES
        and exclusion_ids == EXPECTED_SOURCE_CONTRACT_EXCLUSIONS
    )
    result = {
        "schema_version": 1,
        "profile": "v3000_indigo_ordinary_semantic_roundtrip_v1",
        "indigo_version": indigo_version,
        "rdkit_version": Chem.rdBase.rdkitVersion,
        **cli_provenance(cli),
        "cases": len(rows),
        "comparable_cases": len(comparable_rows),
        "semantic_equal": len(comparable_rows) - len(failures),
        "failures": failures,
        "expected_comparable_cases": EXPECTED_COMPARABLE_CASES,
        "expected_source_contract_exclusions": sorted(EXPECTED_SOURCE_CONTRACT_EXCLUSIONS),
        "gate_passed": gate_passed,
        "source_contract_exclusions": [
            {"id": row["id"], "reason": row["source_indigo_error"]}
            for row in rows
            if row.get("excluded_by_source_contract")
        ],
        "scope": "Indigo-written and Indigo-reopenable ordinary V3000: RDKit topology/atom number/formal charge/isotope/bond direction/isomeric identity plus Indigo canonical-SMILES reopen",
        "not_claimed": [
            "SGROUP/COLLECTION typed semantic interoperability",
            "ENDPTS/ATTACH interpretation or haptic chemistry",
            "query or polymer expansion semantics",
            "Indigo version or platform portability beyond this pinned local binding",
        ],
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {key: result[key] for key in ("cases", "comparable_cases", "semantic_equal", "failures")}
        )
    )
    return 0 if gate_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
