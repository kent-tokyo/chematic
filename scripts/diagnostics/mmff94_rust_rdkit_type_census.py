#!/usr/bin/env python3
"""Compare source-built CheMatic MMFF94 types with RDKit at atom indices.

Build the ``mmff94_type_probe`` Rust example first. The probe writes atomic
number alongside each type so mismatched atom ordering cannot be mistaken for
a typing difference. RDKit MMFF-unsupported rows are counted separately.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-exact", action="store_true")
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")

    corpus_text = args.corpus.read_text(encoding="utf-8")
    smiles = [line.split()[0] for line in corpus_text.splitlines()]
    probe = subprocess.run(
        [str(args.binary.resolve())], input=corpus_text, text=True,
        capture_output=True, check=True,
    )
    records = probe.stdout.splitlines()
    if len(records) != len(smiles):
        raise ValueError(f"probe returned {len(records)} rows for {len(smiles)} inputs")

    status: Counter[str] = Counter()
    mismatch_kinds: Counter[tuple[str, int, int]] = Counter()
    mismatch_rows: list[dict] = []
    compared_atoms: Counter[str] = Counter()
    for index, (smi, record) in enumerate(zip(smiles, records)):
        returned_index, payload = record.split("\t", 1)
        if int(returned_index) != index:
            raise ValueError(f"probe row out of order at {index}")
        rdkit_base = Chem.MolFromSmiles(smi)
        if rdkit_base is None:
            status["rdkit_parse_failure"] += 1
            continue
        rdkit_mol = Chem.AddHs(rdkit_base)
        properties = AllChem.MMFFGetMoleculeProperties(rdkit_mol, mmffVariant="MMFF94")
        if properties is None:
            status["rdkit_mmff_unsupported"] += 1
            continue
        if payload == "E":
            status["chematic_error"] += 1
            continue
        atoms = [tuple(map(int, item.split(":"))) for item in payload.split(",")]
        if len(atoms) != rdkit_mol.GetNumAtoms():
            status["atom_count_mismatch"] += 1
            continue
        if any(a.GetAtomicNum() != atomic_number for a, (atomic_number, _) in zip(rdkit_mol.GetAtoms(), atoms)):
            status["atom_order_or_element_mismatch"] += 1
            continue
        status["compared"] += 1
        row_diffs = 0
        for atom, (_, chematic_type) in zip(rdkit_mol.GetAtoms(), atoms):
            bucket = "hydrogen" if atom.GetAtomicNum() == 1 else "heavy"
            compared_atoms[bucket] += 1
            rdkit_type = properties.GetMMFFAtomType(atom.GetIdx())
            if rdkit_type != chematic_type:
                mismatch_kinds[(atom.GetSymbol(), rdkit_type, chematic_type)] += 1
                row_diffs += 1
        if row_diffs:
            mismatch_rows.append({"input_index": index, "differing_atoms": row_diffs})

    result = {
        "schema_version": 1,
        "profile": "mmff94_source_rust_vs_rdkit_atom_types",
        "corpus": {"path": str(args.corpus), "sha256": sha256(args.corpus), "rows": len(smiles)},
        "binary": {"path": str(args.binary), "sha256": sha256(args.binary)},
        "rdkit_version": rdBase.rdkitVersion,
        "host": {"platform": platform.platform(), "machine": platform.machine()},
        "row_status": dict(sorted(status.items())),
        "atoms_compared": dict(sorted(compared_atoms.items())),
        "differing_rows": mismatch_rows,
        "differing_atoms": sum(mismatch_kinds.values()),
        "difference_categories": [
            {"element": element, "rdkit_type": rdkit_type,
             "chematic_type": chematic_type, "count": count}
            for (element, rdkit_type, chematic_type), count in sorted(mismatch_kinds.items())
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"row_status": result["row_status"],
                      "differing_rows": len(mismatch_rows),
                      "differing_atoms": result["differing_atoms"]}, sort_keys=True))
    return int(args.require_exact and (
        mismatch_rows or status["chematic_error"] or status["atom_count_mismatch"]
        or status["atom_order_or_element_mismatch"]
    ))


if __name__ == "__main__":
    raise SystemExit(main())
