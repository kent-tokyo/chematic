#!/usr/bin/env python3
"""Audit accurate CIP atom and E/Z bond labels against a pinned RDKit oracle.

The input is a provenance-carrying JSONL corpus with one ``smiles`` field per
record.  Atom R/S labels and bond E/Z labels are compared as separate maps;
accurate-mode unresolved entries are retained as a third outcome rather than
being counted as silently correct or incorrect assignments.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdCIPLabeler

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]
EXPECTED_RDKIT = "2025.09.3"


def load_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        record = json.loads(line)
        if not isinstance(record, dict):
            raise ValueError(f"{path}:{line_number}: row must be an object")
        if record.get("_manifest") is True:
            continue
        if not isinstance(record.get("smiles"), str) or not record["smiles"].strip():
            raise ValueError(f"{path}:{line_number}: missing non-empty smiles")
        rows.append(record)
    return rows


def rdkit_labels(smiles: str) -> tuple[dict[int, str], dict[int, str]] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rdCIPLabeler.AssignCIPLabels(mol)
    atoms = {
        atom.GetIdx(): atom.GetProp("_CIPCode")
        for atom in mol.GetAtoms()
        if atom.HasProp("_CIPCode")
    }
    bonds: dict[int, str] = {}
    for bond in mol.GetBonds():
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREOTRANS:
            bonds[bond.GetIdx()] = "E"
        elif stereo == Chem.BondStereo.STEREOCIS:
            bonds[bond.GetIdx()] = "Z"
    return atoms, bonds


def chematic_labels(mol) -> tuple[dict[int, str], dict[int, str], dict[int, str]]:
    assignments = mol.cip_stereo(mode="accurate")
    atoms: dict[int, str] = {}
    bonds: dict[int, str] = {}
    for item in assignments:
        index = int(item["atom_idx"])
        descriptor = str(item["descriptor"])
        if descriptor in {"R", "S", "r", "s"}:
            atoms[index] = descriptor
        elif descriptor in {"E", "Z"}:
            bonds[index] = descriptor
    unresolved = {
        int(item["atom_idx"]): str(item["reason"])
        for item in mol.cip_stereo_unresolved()
    }
    return atoms, bonds, unresolved


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rdkit", default=EXPECTED_RDKIT)
    args = parser.parse_args()

    if rdBase.rdkitVersion != args.expected_rdkit:
        print(
            f"FATAL: expected rdkit=={args.expected_rdkit}, got {rdBase.rdkitVersion}",
            file=sys.stderr,
        )
        return 2

    try:
        records = load_rows(args.corpus)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 2

    summary = {
        "schema_version": 1,
        "profile": "A2-cip-label-parity",
        "rdkit_version": rdBase.rdkitVersion,
        "corpus": str(args.corpus.resolve().relative_to(ROOT))
        if args.corpus.resolve().is_relative_to(ROOT)
        else str(args.corpus.resolve()),
        "corpus_sha256": hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
        "input_rows": len(records),
        "compared_rows": 0,
        "parse_failures": 0,
        "engine_failures": 0,
        "atom": {"oracle": 0, "assigned": 0, "exact": 0, "mismatches": 0},
        "bond": {"oracle": 0, "assigned": 0, "exact": 0, "mismatches": 0},
        "unresolved_reasons": Counter(),
        "atom_mismatch_classes": Counter(),
        "atom_mismatch_elements": Counter(),
        "atom_mismatch_oracle_labels": Counter(),
        "mismatch_examples": [],
    }

    import chematic

    for row_index, record in enumerate(records):
        smiles = record["smiles"].strip()
        oracle = rdkit_labels(smiles)
        if oracle is None:
            summary["parse_failures"] += 1
            continue
        rd_mol = Chem.MolFromSmiles(smiles)
        assert rd_mol is not None
        element_by_index = {atom.GetIdx(): atom.GetSymbol() for atom in rd_mol.GetAtoms()}
        try:
            candidate = chematic.from_smiles(smiles)
            candidate_atoms, candidate_bonds, unresolved = chematic_labels(candidate)
        except Exception as exc:  # pragma: no cover - environment-specific
            summary["engine_failures"] += 1
            if len(summary["mismatch_examples"]) < 50:
                summary["mismatch_examples"].append(
                    {"row": row_index, "smiles": smiles, "kind": "engine_error", "error": repr(exc)}
                )
            continue

        summary["compared_rows"] += 1
        oracle_atoms, oracle_bonds = oracle
        for reason in unresolved.values():
            summary["unresolved_reasons"][reason] += 1

        for kind, expected, actual in (
            ("atom", oracle_atoms, candidate_atoms),
            ("bond", oracle_bonds, candidate_bonds),
        ):
            stats = summary[kind]
            stats["oracle"] += len(expected)
            stats["assigned"] += len(actual)
            if expected == actual:
                stats["exact"] += 1
            else:
                stats["mismatches"] += 1
                if kind == "atom":
                    unresolved_for_expected = {
                        unresolved[index]
                        for index in expected
                        if index in unresolved
                    }
                    if unresolved_for_expected:
                        mismatch_class = "+".join(sorted(unresolved_for_expected))
                    else:
                        mismatch_class = "assigned_label_mismatch"
                    summary["atom_mismatch_classes"][mismatch_class] += 1
                    for index in expected:
                        summary["atom_mismatch_elements"][element_by_index[index]] += 1
                        summary["atom_mismatch_oracle_labels"][expected[index]] += 1
                if len(summary["mismatch_examples"]) < 50:
                    summary["mismatch_examples"].append(
                        {
                            "row": row_index,
                            "smiles": smiles,
                            "kind": kind,
                            "oracle": {str(k): v for k, v in sorted(expected.items())},
                            "candidate": {str(k): v for k, v in sorted(actual.items())},
                            "unresolved": {str(k): v for k, v in sorted(unresolved.items())},
                        }
                    )

    summary["unresolved_reasons"] = dict(sorted(summary["unresolved_reasons"].items()))
    summary["atom_mismatch_classes"] = dict(sorted(summary["atom_mismatch_classes"].items()))
    summary["atom_mismatch_elements"] = dict(sorted(summary["atom_mismatch_elements"].items()))
    summary["atom_mismatch_oracle_labels"] = dict(sorted(summary["atom_mismatch_oracle_labels"].items()))
    summary["gate_passed"] = (
        summary["compared_rows"] == summary["input_rows"]
        and summary["parse_failures"] == 0
        and summary["engine_failures"] == 0
        and summary["atom"]["mismatches"] == 0
        and summary["bond"]["mismatches"] == 0
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
