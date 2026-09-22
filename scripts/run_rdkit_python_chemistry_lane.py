#!/usr/bin/env python3
"""Run complete-row SMILES/CIP/SMARTS/Morgan evidence against one RDKit wheel.

The corpus must be exposed development data.  Every input row is written once
to JSONL, including parser failures, typed CheMatic refusals, and mismatches.
This runner never uses accuracy-sealed cohorts and never changes product code.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_smiles(path: Path, limit: int | None) -> list[str]:
    rows: list[str] = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), 1
    ):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("{"):
            record = json.loads(stripped)
            if record.get("_manifest") is True:
                continue
            value = record.get("smiles")
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{path}:{line_number}: missing smiles")
            rows.append(value.strip())
        else:
            rows.append(stripped.split()[0])
        if limit is not None and len(rows) == limit:
            break
    if not rows:
        raise ValueError(f"no SMILES rows in {path}")
    return rows


def load_queries(path: Path) -> list[str]:
    document = json.loads(path.read_text(encoding="utf-8"))
    queries = document.get("queries")
    if (
        not isinstance(queries, list)
        or not queries
        or not all(isinstance(item, str) and item for item in queries)
    ):
        raise ValueError(f"{path}: queries must be a non-empty string array")
    if len(queries) != len(set(queries)):
        raise ValueError(f"{path}: duplicate SMARTS queries")
    return queries


def bond_endpoint_key(atom1: int, atom2: int) -> str:
    left, right = sorted((atom1, atom2))
    return f"{left}-{right}"


def rdkit_cip(Chem, rdCIPLabeler, molecule) -> tuple[dict[int, str], dict[str, str]]:
    rdCIPLabeler.AssignCIPLabels(molecule)
    atoms = {
        atom.GetIdx(): atom.GetProp("_CIPCode")
        for atom in molecule.GetAtoms()
        if atom.HasProp("_CIPCode")
    }
    bonds: dict[str, str] = {}
    for bond in molecule.GetBonds():
        endpoint = bond_endpoint_key(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if bond.GetStereo() == Chem.BondStereo.STEREOTRANS:
            bonds[endpoint] = "E"
        elif bond.GetStereo() == Chem.BondStereo.STEREOCIS:
            bonds[endpoint] = "Z"
    return atoms, bonds


def chematic_cip(molecule) -> tuple[dict[int, str], dict[str, str], dict[int, str]]:
    atoms: dict[int, str] = {}
    bonds: dict[str, str] = {}
    bond_table = molecule.bond_table
    for item in molecule.cip_stereo(mode="accurate"):
        index = int(item["atom_idx"])
        descriptor = str(item["descriptor"])
        if descriptor in {"R", "S", "r", "s"}:
            atoms[index] = descriptor
        elif descriptor in {"E", "Z"}:
            if index >= len(bond_table):
                bonds[f"invalid-bond-index-{index}"] = descriptor
            else:
                atom1, atom2, _bond_type, _is_aromatic = bond_table[index]
                bonds[bond_endpoint_key(int(atom1), int(atom2))] = descriptor
    unresolved = {
        int(item["atom_idx"]): str(item["reason"])
        for item in molecule.cip_stereo_unresolved()
    }
    return atoms, bonds, unresolved


def index_correspondence(rdkit_molecule, chematic_molecule) -> dict[str, bool]:
    rdkit_atoms = [atom.GetAtomicNum() for atom in rdkit_molecule.GetAtoms()]
    chematic_atoms = [int(row[1]) for row in chematic_molecule.atom_table]
    rdkit_bonds = {
        bond_endpoint_key(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for bond in rdkit_molecule.GetBonds()
    }
    chematic_bonds = {
        bond_endpoint_key(int(row[0]), int(row[1]))
        for row in chematic_molecule.bond_table
    }
    return {
        "atom_order": rdkit_atoms == chematic_atoms,
        "bond_endpoints": rdkit_bonds == chematic_bonds,
    }


def normalized_match_sets(matches) -> list[list[int]]:
    normalized = sorted(
        {tuple(sorted(int(index) for index in match)) for match in matches}
    )
    return [list(match) for match in normalized]


def difference(
    operation: str,
    detail: str,
    classification: str = "unresolved",
) -> dict[str, str]:
    return {
        "operation": operation,
        "classification": classification,
        "detail": detail,
    }


def count_difference_classes(
    differences: list[dict[str, str]], counts: Counter[str]
) -> None:
    for item in differences:
        counts[f"difference_class_{item['classification']}"] += 1


def example_record(row: dict[str, object]) -> dict[str, object]:
    return {
        "input_index": row["input_index"],
        "smiles": row["smiles"],
        "status": row["status"],
        "differences": row["differences"],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rows-output", type=Path, required=True)
    parser.add_argument("--expected-rdkit", required=True)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    for path in (args.corpus, args.queries):
        if not path.is_file():
            parser.error(f"input not found: {path}")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    return args


def main() -> int:
    args = parse_args()
    import chematic
    from rdkit import Chem, DataStructs, RDLogger, rdBase
    from rdkit.Chem import rdCIPLabeler, rdFingerprintGenerator

    RDLogger.DisableLog("rdApp.*")
    smiles_rows = load_smiles(args.corpus, args.limit)
    queries = load_queries(args.queries)
    query_molecules = {query: Chem.MolFromSmarts(query) for query in queries}
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    counts: Counter[str] = Counter()
    examples: list[dict[str, object]] = []

    args.rows_output.parent.mkdir(parents=True, exist_ok=True)
    with args.rows_output.open("w", encoding="utf-8") as rows_handle:
        for index, smiles in enumerate(smiles_rows):
            row: dict[str, object] = {
                "input_index": index,
                "smiles": smiles,
                "status": "completed",
                "differences": [],
            }
            differences: list[dict[str, str]] = row["differences"]  # type: ignore[assignment]
            rd_mol = Chem.MolFromSmiles(smiles)
            try:
                candidate = chematic.from_smiles(smiles)
                candidate_error = None
            except Exception as exc:
                candidate = None
                candidate_error = f"{type(exc).__name__}: {exc}"
            if rd_mol is None or candidate is None:
                row["status"] = "parse_failure"
                row["rdkit_parse"] = "success" if rd_mol is not None else "failed"
                row["chematic_parse"] = (
                    "success" if candidate is not None else {"error": candidate_error}
                )
                differences.append(
                    difference(
                        "smiles_parse_write", "parse outcome differs or both failed"
                    )
                )
                counts["parse_failure"] += 1
                count_difference_classes(differences, counts)
                rows_handle.write(json.dumps(row, sort_keys=True) + "\n")
                if len(examples) < 50:
                    examples.append(example_record(row))
                continue

            oracle_smiles = Chem.MolToSmiles(
                rd_mol, canonical=True, isomericSmiles=True
            )
            candidate_smiles = candidate.smiles
            candidate_roundtrip = Chem.MolFromSmiles(candidate_smiles)
            semantic_smiles = (
                Chem.MolToSmiles(
                    candidate_roundtrip, canonical=True, isomericSmiles=True
                )
                if candidate_roundtrip is not None
                else None
            )
            smiles_exact = candidate_smiles == oracle_smiles
            smiles_semantic = semantic_smiles == oracle_smiles
            nonisomeric_roundtrip = (
                candidate_roundtrip is not None
                and Chem.MolToSmiles(rd_mol, canonical=True, isomericSmiles=False)
                == Chem.MolToSmiles(
                    candidate_roundtrip, canonical=True, isomericSmiles=False
                )
            )
            row["smiles_parse_write"] = {
                "rdkit_canonical": oracle_smiles,
                "chematic_canonical": candidate_smiles,
                "exact_spelling": smiles_exact,
                "semantic_roundtrip": smiles_semantic,
                "nonisomeric_roundtrip": nonisomeric_roundtrip,
            }
            counts[
                "smiles_exact" if smiles_exact else "smiles_spelling_difference"
            ] += 1
            if not smiles_exact and smiles_semantic:
                differences.append(
                    difference(
                        "smiles_parse_write",
                        "canonical spelling differs while RDKit semantic identity is preserved",
                        "contract_difference",
                    )
                )
            if not smiles_semantic:
                if nonisomeric_roundtrip:
                    differences.append(
                        difference(
                            "smiles_parse_write",
                            "stereo identity is lost or altered while non-isomeric graph identity is preserved",
                            "chematic_regression",
                        )
                    )
                else:
                    differences.append(
                        difference(
                            "smiles_parse_write",
                            "canonical output changes RDKit graph or semantic identity",
                        )
                    )
                counts["smiles_semantic_difference"] += 1

            oracle_atoms, oracle_bonds = rdkit_cip(Chem, rdCIPLabeler, rd_mol)
            candidate_atoms, candidate_bonds, unresolved = chematic_cip(candidate)
            correspondence = index_correspondence(rd_mol, candidate)
            cip_exact = (
                all(correspondence.values())
                and oracle_atoms == candidate_atoms
                and oracle_bonds == candidate_bonds
            )
            row["cip"] = {
                "rdkit_atoms": oracle_atoms,
                "chematic_atoms": candidate_atoms,
                "rdkit_bonds": oracle_bonds,
                "chematic_bonds": candidate_bonds,
                "chematic_unresolved": unresolved,
                "index_correspondence": correspondence,
                "exact": cip_exact,
            }
            counts["cip_exact" if cip_exact else "cip_difference"] += 1
            if not cip_exact:
                detail = (
                    "atom-index correspondence or bond-endpoint topology differs"
                    if not all(correspondence.values())
                    else "atom or bond-endpoint label map differs"
                )
                differences.append(difference("cip", detail))

            oracle_fp = DataStructs.BitVectToBinaryText(
                generator.GetFingerprint(rd_mol)
            )
            try:
                candidate_fp = bytes(candidate.rdkit_ecfp4())
                morgan_error = None
            except Exception as exc:
                candidate_fp = None
                morgan_error = f"{type(exc).__name__}: {exc}"
            morgan_exact = candidate_fp == oracle_fp
            row["morgan"] = {
                "exact": morgan_exact,
                "rdkit_sha256": hashlib.sha256(oracle_fp).hexdigest(),
                "chematic_sha256": (
                    hashlib.sha256(candidate_fp).hexdigest()
                    if candidate_fp is not None
                    else None
                ),
                "chematic_error": morgan_error,
            }
            counts["morgan_exact" if morgan_exact else "morgan_difference"] += 1
            if not morgan_exact:
                if morgan_error is not None and "unsupported" in morgan_error.lower():
                    differences.append(
                        difference(
                            "morgan",
                            "typed unsupported outcome; no packed fingerprint was produced",
                            "contract_difference",
                        )
                    )
                else:
                    differences.append(
                        difference(
                            "morgan",
                            "packed radius-2/2048 fingerprint differs or failed without a typed unsupported outcome",
                        )
                    )

            smarts_differences: list[dict[str, object]] = []
            for query in queries:
                rd_query = query_molecules[query]
                if rd_query is None:
                    rd_matches = None
                else:
                    rd_matches = normalized_match_sets(
                        rd_mol.GetSubstructMatches(rd_query, uniquify=True)
                    )
                try:
                    candidate_matches = normalized_match_sets(
                        chematic.smarts_find(query, candidate)
                    )
                    candidate_query_error = None
                except Exception as exc:
                    candidate_matches = None
                    candidate_query_error = f"{type(exc).__name__}: {exc}"
                if rd_matches != candidate_matches:
                    smarts_differences.append(
                        {
                            "query": query,
                            "rdkit": rd_matches,
                            "chematic": candidate_matches,
                            "chematic_error": candidate_query_error,
                        }
                    )
            row["smarts"] = {
                "query_count": len(queries),
                "difference_count": len(smarts_differences),
                "differences": smarts_differences,
            }
            counts["smarts_cells"] += len(queries)
            counts["smarts_differences"] += len(smarts_differences)
            if smarts_differences:
                differences.append(
                    difference("smarts", "one or more target-atom-set cells differ")
                )
            counts["completed"] += 1
            counts["difference_rows" if differences else "exact_rows"] += 1
            count_difference_classes(differences, counts)
            rows_handle.write(json.dumps(row, sort_keys=True) + "\n")
            if differences and len(examples) < 50:
                examples.append(example_record(row))

    summary = {
        "schema_version": 1,
        "profile": "rdkit_python_chemistry_rebaseline_v1",
        "rdkit_version": rdBase.rdkitVersion,
        "chematic_version": getattr(chematic, "__version__", None),
        "corpus": {"path": str(args.corpus.resolve()), "sha256": sha256(args.corpus)},
        "queries": {
            "path": str(args.queries.resolve()),
            "sha256": sha256(args.queries),
            "count": len(queries),
        },
        "rows": {
            "path": str(args.rows_output.resolve()),
            "sha256": sha256(args.rows_output),
        },
        "row_accounting": {
            "input_count": len(smiles_rows),
            "completed_count": counts["completed"],
            "parse_failure_count": counts["parse_failure"],
        },
        "counts": dict(sorted(counts.items())),
        "difference_classification": {
            "allowed": [
                "chematic_regression",
                "oracle_change",
                "contract_difference",
                "unresolved",
            ],
            "default_for_unadjudicated_rows": "unresolved",
        },
        "difference_examples": examples,
        "sealed_accuracy_cohort_reused": False,
    }
    summary["gate"] = {
        "runtime_matches": rdBase.rdkitVersion == args.expected_rdkit,
        "row_accounting_complete": len(smiles_rows)
        == counts["completed"] + counts["parse_failure"],
        "rows_hash_recorded": bool(summary["rows"]["sha256"]),
    }
    summary["gate_passed"] = all(summary["gate"].values())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["gate_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
