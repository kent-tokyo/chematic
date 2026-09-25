#!/usr/bin/env python3
"""Census of MMFF94 atom-type agreement between CheMatic and RDKit.

Every SMILES in the corpus is parsed by both engines, hydrogens are made
explicit (RDKit ``AddHs``, CheMatic ``add_hydrogens``; both append H atoms
after the heavy atoms in the same order), and each atom's MMFF94 numeric type
is compared with RDKit's ``MMFFGetMMFFAtomType``. Hydrogen disagreements are
keyed by (parent element, RDKit parent type, CheMatic parent type, RDKit H
type, CheMatic H type) so a hydrogen that only disagrees because its parent
does is visible as such. Heavy-atom disagreements are keyed by (element,
RDKit type, CheMatic type).

Rows are never dropped silently: RDKit parse failures, RDKit MMFF setup
failures, CheMatic errors and atom-count mismatches are counted by status.
This measures typing only, not parameters or energies (see
``mmff94_same_explicit_h_energy.py`` for those).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from collections import Counter
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem


def call(value):
    return value() if callable(value) else value


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--examples", type=int, default=1, help="example rows kept per key")
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")

    status: Counter[str] = Counter()
    h_diff: Counter[tuple] = Counter()
    heavy_diff: Counter[tuple] = Counter()
    examples: dict[tuple, list[int]] = {}
    atoms = Counter()
    rows_with_h_diff = 0
    rows_with_heavy_diff = 0
    lines = [line.split()[0] for line in args.corpus.read_text(encoding="utf-8").splitlines() if line.strip()]
    for index, smiles in enumerate(lines):
        rdkit_base = Chem.MolFromSmiles(smiles)
        if rdkit_base is None:
            status["rdkit_parse_failure"] += 1
            continue
        rdkit_mol = Chem.AddHs(rdkit_base)
        properties = AllChem.MMFFGetMoleculeProperties(rdkit_mol, mmffVariant="MMFF94")
        if properties is None:
            status["rdkit_mmff_unsupported"] += 1
            continue
        try:
            types = list(call(chematic.from_smiles(smiles).add_hydrogens().mmff94_numeric_atom_types))
        except Exception:  # noqa: BLE001 - typed refusal is a counted status
            status["chematic_error"] += 1
            continue
        if len(types) != rdkit_mol.GetNumAtoms():
            status["atom_count_mismatch"] += 1
            continue
        status["compared"] += 1
        row_h = row_heavy = False
        for atom in rdkit_mol.GetAtoms():
            i = atom.GetIdx()
            rdkit_type = properties.GetMMFFAtomType(i)
            if atom.GetAtomicNum() == 1:
                atoms["hydrogen"] += 1
                if rdkit_type != types[i]:
                    parent = atom.GetNeighbors()[0]
                    j = parent.GetIdx()
                    key = ("H", parent.GetSymbol(), properties.GetMMFFAtomType(j), types[j], rdkit_type, types[i])
                    h_diff[key] += 1
                    row_h = True
                    examples.setdefault(key, [])
                    if len(examples[key]) < args.examples and index not in examples[key]:
                        examples[key].append(index)
            else:
                atoms["heavy"] += 1
                if rdkit_type != types[i]:
                    key = ("heavy", atom.GetSymbol(), rdkit_type, types[i])
                    heavy_diff[key] += 1
                    row_heavy = True
                    examples.setdefault(key, [])
                    if len(examples[key]) < args.examples and index not in examples[key]:
                        examples[key].append(index)
        rows_with_h_diff += row_h
        rows_with_heavy_diff += row_heavy

    def listing(counter: Counter[tuple], fields: list[str]) -> list[dict]:
        return [
            dict(zip(fields, key[1:]), count=count, example_input_indices=examples[key])
            for key, count in counter.most_common()
        ]

    hydrogen_parent_agrees = sum(
        count for key, count in h_diff.items() if key[2] == key[3]
    )
    summary = {
        "schema_version": 1,
        "profile": "mmff94_atom_type_census_v1",
        "corpus": {"path": str(args.corpus), "sha256": sha256(args.corpus), "rows": len(lines)},
        "chematic_version": getattr(chematic, "__version__", None),
        "rdkit_version": rdBase.rdkitVersion,
        "host": {"platform": platform.platform(), "machine": platform.machine()},
        "row_status": dict(sorted(status.items())),
        "atoms_compared": dict(atoms),
        "hydrogen": {
            "differing_atoms": sum(h_diff.values()),
            "differing_atoms_with_agreeing_parent_type": hydrogen_parent_agrees,
            "rows_with_difference": rows_with_h_diff,
            "by_key": listing(h_diff, ["parent_element", "rdkit_parent_type", "chematic_parent_type", "rdkit_type", "chematic_type"]),
        },
        "heavy": {
            "differing_atoms": sum(heavy_diff.values()),
            "rows_with_difference": rows_with_heavy_diff,
            "by_key": listing(heavy_diff, ["element", "rdkit_type", "chematic_type"]),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("row_status", "atoms_compared")} | {
        "hydrogen_differing": summary["hydrogen"]["differing_atoms"],
        "hydrogen_differing_with_agreeing_parent": hydrogen_parent_agrees,
        "heavy_differing": summary["heavy"]["differing_atoms"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
