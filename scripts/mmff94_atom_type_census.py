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
    parser.add_argument("--wheel", type=Path, help="source wheel used for an artifact gate")
    parser.add_argument("--module-root", type=Path, help="installation root for the source wheel")
    parser.add_argument(
        "--source-wheel-gate", action="store_true",
        help="require the predeclared v1.0.31-source MMFF94 typing profile",
    )
    parser.add_argument("--examples", type=int, default=1, help="example rows kept per key")
    parser.add_argument(
        "--type37-context",
        action="store_true",
        help="diagnose the RDKit 37 / CheMatic 2 carbon bucket without changing the base census",
    )
    args = parser.parse_args()
    if args.source_wheel_gate:
        if not (args.wheel and args.module_root and args.type37_context):
            parser.error("source-wheel gate requires --wheel, --module-root and --type37-context")
        module_path = Path(chematic.__file__).resolve()
        if not module_path.is_relative_to(args.module_root.resolve()):
            parser.error(f"imported {module_path}, not the isolated wheel")
        if rdBase.rdkitVersion != "2026.03.6":
            parser.error(f"expected RDKit 2026.03.6, got {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")

    status: Counter[str] = Counter()
    h_diff: Counter[tuple] = Counter()
    heavy_diff: Counter[tuple] = Counter()
    examples: dict[tuple, list[int]] = {}
    atoms = Counter()
    rows_with_h_diff = 0
    rows_with_heavy_diff = 0
    type37_context: Counter[str] = Counter()
    type37_examples: list[dict] = []
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
            chematic_mol = chematic.from_smiles(smiles).add_hydrogens()
            types = list(call(chematic_mol.mmff94_numeric_atom_types))
        except Exception:  # noqa: BLE001 - typed refusal is a counted status
            status["chematic_error"] += 1
            continue
        if len(types) != rdkit_mol.GetNumAtoms():
            status["atom_count_mismatch"] += 1
            continue
        status["compared"] += 1
        row_h = row_heavy = False
        type37_atoms: list[int] = []
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
                    if atom.GetAtomicNum() == 6 and rdkit_type == 37 and types[i] == 2:
                        type37_atoms.append(i)
                    examples.setdefault(key, [])
                    if len(examples[key]) < args.examples and index not in examples[key]:
                        examples[key].append(index)
        rows_with_h_diff += row_h
        rows_with_heavy_diff += row_heavy
        if args.type37_context and type37_atoms:
            # These two aromaticity flags are not equivalent definitions:
            # RDKit's flag is the sanitized molecule's general aromaticity;
            # CheMatic's is its force-field-specific re-perception. Record
            # both as context, not as an MMFF aromaticity parity percentage.
            flags = list(call(chematic_mol.mmff94_numeric_aromatic_flags))
            if len(flags) != len(types):
                raise ValueError(f"MMFF aromatic flag count mismatch at row {index}")
            rdkit_rings = {tuple(sorted(ring)) for ring in rdkit_mol.GetRingInfo().AtomRings()}
            chematic_rings = {
                tuple(sorted(ring)) for ring in call(chematic_mol.symmetrized_sssr_atom_rings)
            }
            missing_rings = rdkit_rings - chematic_rings
            extra_rings = chematic_rings - rdkit_rings
            type37_context["rows"] += 1
            type37_context["atoms"] += len(type37_atoms)
            type37_context["rdkit_default_aromatic_atoms"] += sum(
                rdkit_mol.GetAtomWithIdx(i).GetIsAromatic() for i in type37_atoms
            )
            type37_context["chematic_mmff_aromatic_atoms"] += sum(flags[i] for i in type37_atoms)
            type37_context["rows_missing_rdkit_rings"] += bool(missing_rings)
            type37_context["rows_with_extra_chematic_rings"] += bool(extra_rings)
            type37_context["missing_rdkit_rings"] += len(missing_rings)
            type37_context["extra_chematic_rings"] += len(extra_rings)
            if len(type37_examples) < 10:
                type37_examples.append(
                    {
                        "input_index": index,
                        "mismatched_atom_indices": type37_atoms,
                        "rdkit_ring_count": len(rdkit_rings),
                        "chematic_ring_count": len(chematic_rings),
                        "missing_rdkit_rings": [list(ring) for ring in sorted(missing_rings)],
                    }
                )

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
    if args.type37_context:
        summary["type37_context"] = {
            "scope": "RDKit numeric type 37 / CheMatic numeric type 2 carbon atoms only",
            "note": "RDKit general aromatic and CheMatic MMFF aromatic flags are diagnostic context, not comparable MMFF flags",
            "counts": dict(sorted(type37_context.items())),
            "examples": type37_examples,
        }
    if args.source_wheel_gate:
        summary["scope"] = "source-built Python wheel; not a published registry artifact"
        summary["wheel_sha256"] = sha256(args.wheel)
        summary["module_path_within_wheel"] = str(module_path.relative_to(args.module_root.resolve()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    if args.source_wheel_gate:
        type37_bucket = next(
            (item["count"] for item in summary["heavy"]["by_key"]
             if item["element"] == "C" and item["rdkit_type"] == 37 and item["chematic_type"] == 2),
            0,
        )
        expected_corpus_sha = "f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f"
        if (summary["corpus"]["sha256"] != expected_corpus_sha
                or summary["corpus"]["rows"] != 10_000
                or summary["row_status"] != {"compared": 9_774, "rdkit_mmff_unsupported": 204, "chematic_error": 22}
                or summary["atoms_compared"] != {"heavy": 214_990, "hydrogen": 190_737}
                or summary["heavy"]["differing_atoms"] != 89
                or summary["hydrogen"]["differing_atoms"] != 1
                or type37_bucket != 16):
            raise ValueError("source-wheel MMFF94 type gate failed; inspect the written census report")
    print(json.dumps({k: summary[k] for k in ("row_status", "atoms_compared")} | {
        "hydrogen_differing": summary["hydrogen"]["differing_atoms"],
        "hydrogen_differing_with_agreeing_parent": hydrogen_parent_agrees,
        "heavy_differing": summary["heavy"]["differing_atoms"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
