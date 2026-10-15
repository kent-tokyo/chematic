#!/usr/bin/env python3
"""Compare source-preserving RDKit V3000 SDF reads across three engines.

Build the current CheMatic source first:
  cargo build -p chematic-mol --example sdf_v3000_accuracy_dump --offline

Then run this script with Python containing RDKit and COSMolKit:
  python check_sdf_v3000_accuracy.py --corpus corpora/stress-v1.smi \
      --chematic ../../target/debug/examples/sdf_v3000_accuracy_dump

Only RDKit MOL blocks that RDKit itself reads back with identical canonical
isomeric SMILES are eligible. Explicit H atoms remain explicit during the
SMILES normalization used for all three engine outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

import cosmolkit
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_with_explicit_h(smiles: str) -> str | None:
    params = Chem.SmilesParserParams()
    params.removeHs = False
    molecule = Chem.MolFromSmiles(smiles, params)
    return Chem.MolToSmiles(molecule, isomericSmiles=True) if molecule is not None else None


def prepare_corpus(corpus: Path) -> tuple[list[dict], str, dict]:
    rows: list[dict] = []
    blocks: list[str] = []
    skipped = {"invalid_source": 0, "write_error": 0, "source_changed": 0}
    for line_number, line in enumerate(corpus.read_text().splitlines()):
        if not line.strip():
            continue
        try:
            mol = Chem.MolFromSmiles(line.split()[0])
            if mol is None:
                skipped["invalid_source"] += 1
                continue
            rdDepictor.Compute2DCoords(mol)
            name = f"row-{line_number}"
            mol.SetProp("_Name", name)
            block = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=True)
            reread = Chem.MolFromMolBlock(block, removeHs=False)
            if reread is None:
                skipped["write_error"] += 1
                continue
            expected = Chem.MolToSmiles(mol, isomericSmiles=True)
            if Chem.MolToSmiles(reread, isomericSmiles=True) != expected:
                skipped["source_changed"] += 1
                continue
        except Exception:
            skipped["write_error"] += 1
            continue
        rows.append({"name": name, "expected": expected, "atoms": reread.GetNumAtoms()})
        blocks.append(block + "$$$$\n")
    return rows, "".join(blocks), skipped


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", required=True, type=Path)
    parser.add_argument("--chematic", required=True, type=Path)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")
    corpus = args.corpus.resolve()
    binary = args.chematic.resolve()
    rows, sdf, skipped = prepare_corpus(corpus)
    with tempfile.TemporaryDirectory(prefix="chematic-v3000-accuracy-") as tmpdir:
        sdf_path = Path(tmpdir) / "source-gated.sdf"
        sdf_path.write_text(sdf)
        completed = subprocess.run(
            [str(binary), str(sdf_path)], capture_output=True, text=True, check=True
        )
        chematic_rows = json.loads(completed.stdout)
        cosmolkit_rows = cosmolkit.SdfDataset.open(str(sdf_path))
        if len(chematic_rows) != len(rows) or len(cosmolkit_rows) != len(rows):
            raise AssertionError("engine record count differs from source-gated count")
        misses = {"chematic": [], "cosmolkit": []}
        for index, expected in enumerate(rows):
            chematic_row = chematic_rows[index]
            try:
                cosmolkit_row = cosmolkit_rows[index]
                cosmolkit_smiles = cosmolkit_row.molecule().to_smiles()
                cosmolkit_atoms = cosmolkit_row.molecule().num_atoms()
            except Exception as error:
                cosmolkit_smiles = None
                cosmolkit_atoms = None
                cosmolkit_error = str(error)
            else:
                cosmolkit_error = None
            for engine, smiles, atoms, error in [
                ("chematic", chematic_row.get("smiles"), chematic_row.get("atoms"), chematic_row.get("error")),
                ("cosmolkit", cosmolkit_smiles, cosmolkit_atoms, cosmolkit_error),
            ]:
                actual = canonical_with_explicit_h(smiles) if smiles is not None else None
                if actual != expected["expected"] or atoms != expected["atoms"]:
                    misses[engine].append({
                        "index": index, "name": expected["name"],
                        "expected": expected["expected"], "actual": actual,
                        "expected_atoms": expected["atoms"], "actual_atoms": atoms,
                        "error": error,
                    })
    output = {
        "corpus": str(corpus),
        "corpus_sha256": sha256(corpus.read_bytes()),
        "generated_sdf_sha256": sha256(sdf.encode()),
        "chematic_binary_sha256": sha256(binary.read_bytes()),
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "cosmolkit_version": cosmolkit.__version__,
        "eligible": len(rows), "skipped": skipped,
        "matched": {engine: len(rows) - len(failures) for engine, failures in misses.items()},
        "mismatches": {engine: failures[:20] for engine, failures in misses.items()},
    }
    print(json.dumps(output, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
