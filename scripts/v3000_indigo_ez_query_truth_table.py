#!/usr/bin/env python3
"""Record a V3000 E/Z query truth table through chematic.

RDKit supplies E and Z V3000 inputs with 2D coordinates.  For each query,
this runner checks both E and Z targets in Indigo before and after the bounded
V3000 -> chematic -> V3000 conversion.  RDKit also reopens both payloads to
prove chematic retained the source isomeric identity.

This is deliberately an observation gate, not a claim that Indigo preserves
V3000 query E/Z semantics.  A known Indigo behaviour can yield all four
matches; the result records that discrepancy explicitly while requiring that
chematic does not change the observed table or the RDKit identity.
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


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = ROOT / "target" / "debug" / "chematic"
# Keep both the symmetric 2-butene pair and a mixed-substituent pair. The
# latter prevents a table pass from being an artefact of swapping identical
# methyl substituents around the double bond.
CASES = {
    "e": "C/C=C/C",
    "z": "C/C=C\\C",
    "e_mixed": "Cl/C=C/Br",
    "z_mixed": "Cl/C=C\\Br",
}


def canonical_isomeric(molblock: str) -> str | None:
    mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def rdkit_truth_table(query_molblock: str) -> dict[str, bool]:
    """Check RDKit's chiral substructure predicate for the V3000 query itself."""
    query = Chem.MolFromMolBlock(query_molblock, sanitize=True, removeHs=False)
    if query is None:
        raise RuntimeError("RDKit rejected V3000 query payload")
    return {
        target_id: Chem.MolFromSmiles(target).HasSubstructMatch(query, useChirality=True)
        for target_id, target in CASES.items()
    }


def indigo_truth_table(query_molblock: str) -> dict[str, bool]:
    indigo = Indigo()
    query = indigo.loadQueryMolecule(query_molblock)
    return {
        target_id: bool(indigo.substructureMatcher(indigo.loadMolecule(target)).match(query))
        for target_id, target in CASES.items()
    }


def cli_sha256(cli: Path) -> str:
    with cli.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        parser.error(f"chematic CLI not found: {cli}")

    rows = []
    with tempfile.TemporaryDirectory(prefix="chematic-v3000-indigo-ez-") as raw_tmp:
        tmp = Path(raw_tmp)
        for query_id, smiles in CASES.items():
            source_mol = Chem.MolFromSmiles(smiles)
            if source_mol is None:
                raise RuntimeError(f"RDKit rejected fixture: {smiles}")
            source = Chem.MolToMolBlock(source_mol, forceV3000=True)
            source_path = tmp / f"{query_id}.source.v3000"
            output_path = tmp / f"{query_id}.chematic.v3000"
            source_path.write_text(source, encoding="utf-8")
            run = subprocess.run(
                [
                    str(cli), "convert", "--input-format", "mol_v3000", "--output-format",
                    "mol_v3000", "--input", str(source_path), "--output", str(output_path),
                ],
                capture_output=True,
                text=True,
            )
            row: dict[str, object] = {"query": query_id, "smiles": smiles, "cli_exit_code": run.returncode}
            if run.returncode:
                row["error"] = run.stderr.strip()
                rows.append(row)
                continue
            written = output_path.read_text(encoding="utf-8")
            source_identity = canonical_isomeric(source)
            output_identity = canonical_isomeric(written)
            rdkit_source_table = rdkit_truth_table(source)
            rdkit_output_table = rdkit_truth_table(written)
            source_table = indigo_truth_table(source)
            output_table = indigo_truth_table(written)
            expected = {target_id: target_id == query_id for target_id in CASES}
            row.update(
                {
                    "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
                    "output_sha256": hashlib.sha256(written.encode()).hexdigest(),
                    "rdkit_source_identity": source_identity,
                    "rdkit_output_identity": output_identity,
                    "rdkit_identity_preserved": source_identity == output_identity,
                    "expected_stereo_truth_table": expected,
                    "rdkit_source_truth_table": rdkit_source_table,
                    "rdkit_output_truth_table": rdkit_output_table,
                    "rdkit_table_preserved_by_chematic": rdkit_source_table
                    == rdkit_output_table,
                    "rdkit_matches_expected_stereo": rdkit_source_table == expected,
                    "indigo_source_truth_table": source_table,
                    "indigo_output_truth_table": output_table,
                    "indigo_table_preserved_by_chematic": source_table == output_table,
                    "indigo_matches_expected_stereo": source_table == expected,
                }
            )
            rows.append(row)

    conversion_failures = [row["query"] for row in rows if row.get("cli_exit_code") != 0]
    identity_failures = [row["query"] for row in rows if not row.get("rdkit_identity_preserved", False)]
    rdkit_query_failures = [
        row["query"]
        for row in rows
        if not row.get("rdkit_table_preserved_by_chematic", False)
        or not row.get("rdkit_matches_expected_stereo", False)
    ]
    table_failures = [row["query"] for row in rows if not row.get("indigo_table_preserved_by_chematic", False)]
    indigo_semantic_loss = [row["query"] for row in rows if not row.get("indigo_matches_expected_stereo", False)]
    result = {
        "schema_version": 1,
        "profile": "v3000_indigo_ez_query_truth_table_v2",
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "indigo_version": Indigo().version(),
        "cli_name": cli.name,
        "cli_sha256": cli_sha256(cli),
        "cases": len(rows),
        "conversion_failures": conversion_failures,
        "rdkit_identity_failures": identity_failures,
        "rdkit_query_truth_table_failures": rdkit_query_failures,
        "indigo_table_preservation_failures": table_failures,
        "indigo_stereo_semantic_loss_queries": indigo_semantic_loss,
        "chematic_roundtrip_gate_passed": not (
            conversion_failures or identity_failures or rdkit_query_failures or table_failures
        ),
        "scope": "Four E/Z V3000 queries (symmetric 2-butene plus mixed-substituent chloro/bromo alkene) and four E/Z targets, checked before and after bounded chematic V3000 round trip with RDKit and Indigo query predicates.",
        "not_claimed": [
            "Indigo query E/Z semantic correctness",
            "typed V3000 query interpretation by chematic",
            "general V3000 query interoperability",
        ],
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "cases", "conversion_failures", "rdkit_identity_failures",
        "rdkit_query_truth_table_failures",
        "indigo_table_preservation_failures", "indigo_stereo_semantic_loss_queries",
        "chematic_roundtrip_gate_passed",
    )}))
    return 0 if result["chematic_roundtrip_gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
