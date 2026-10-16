#!/usr/bin/env python3
"""Measure SMILES -> V3000 -> RDKit read-back with three independent writers.

The JSONL corpus contains ``id`` and ``smiles`` fields. Atom maps are removed
from identities for all engines. This is a development comparison, not an
independent chemical truth or an unseen holdout.
"""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import cosmolkit
from rdkit import Chem, RDLogger, rdBase


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CORPUS = ROOT / "validation/openbabel_comparison/v3000-stereo-writer-stress-316.jsonl"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(molecule: Chem.Mol | None) -> str | None:
    if molecule is None:
        return None
    molecule = Chem.Mol(molecule)
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def readback(block: str | None) -> str | None:
    if block is None or "V3000" not in block or "M  V30 BEGIN CTAB" not in block:
        return None
    try:
        return identity(Chem.MolFromMolBlock(block, removeHs=False))
    except (RuntimeError, ValueError):
        return None


def chematic_block(cli: Path, smiles: str) -> str | None:
    try:
        process = subprocess.run(
            [str(cli), "convert", "--input-format", "smiles", "--output-format", "mol_v3000"],
            input=smiles,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return process.stdout if process.returncode == 0 else None


def compare(corpus: Path, cli: Path) -> dict:
    params = cosmolkit.MolBlockWriteParams(format=cosmolkit.SdfFormat.V3000)
    rows = []
    for line in corpus.read_text().splitlines():
        if not line:
            continue
        case = json.loads(line)
        source = case["smiles"]
        molecule = Chem.MolFromSmiles(source)
        if molecule is None:
            raise ValueError(f"RDKit cannot parse {case['id']}")
        expected = identity(molecule)
        try:
            cosmolkit_block = cosmolkit.Molecule.from_smiles(source).to_mol_with_params(params)
        except Exception:
            cosmolkit_block = None
        blocks = {
            "rdkit": Chem.MolToMolBlock(molecule, forceV3000=True),
            "cosmolkit": cosmolkit_block,
            "chematic": chematic_block(cli, source),
        }
        rows.append(
            {
                "id": case["id"],
                "source_smiles": source,
                "expected": expected,
                "source_has_stereo_token": any(token in expected for token in "@/\\"),
                "identities": {engine: readback(block) for engine, block in blocks.items()},
                "block_sha256": {
                    engine: hashlib.sha256(block.encode()).hexdigest() if block else None
                    for engine, block in blocks.items()
                },
            }
        )
    summary = {"cases": len(rows)}
    for engine in ("rdkit", "cosmolkit", "chematic"):
        summary[engine] = {
            "exact": sum(row["identities"][engine] == row["expected"] for row in rows),
            "stereo_exact": sum(
                row["source_has_stereo_token"]
                and row["identities"][engine] == row["expected"]
                for row in rows
            ),
            "failed_ids": [
                row["id"] for row in rows if row["identities"][engine] != row["expected"]
            ],
        }
    return {
        "scope": "development SMILES to V3000 to RDKit read-back; atom maps ignored",
        "corpus": str(corpus),
        "corpus_sha256": digest(corpus),
        "cli": str(cli),
        "cli_sha256": digest(cli),
        "rdkit_version": rdBase.rdkitVersion,
        "cosmolkit_version": cosmolkit.__version__,
        "summary": summary,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")
    report = compare(args.corpus.resolve(), args.cli.resolve())
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
