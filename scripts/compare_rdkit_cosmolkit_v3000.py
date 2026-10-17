#!/usr/bin/env python3
"""Measure SMILES -> V3000 -> RDKit read-back with three independent writers.

JSONL corpora contain ``id`` and ``smiles`` fields; SMILES corpora contain one
source molecule per line. Atom maps are removed
from identities for all engines. Each writer receives the same RDKit-normalized
SMILES so that source spellings RDKit rewrites (such as ``[H]=C``) cannot make
the writers operate on different molecular graphs. This is a development
comparison, not an independent chemical truth or an unseen holdout.
COSMolKit also reads CheMatic's block as a separate interoperability measure.
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


def cosmolkit_readback(block: str | None) -> tuple[str | None, str | None]:
    if block is None:
        return None, "writer produced no block"
    try:
        smiles = cosmolkit.Molecule.from_mol(block).to_smiles()
        actual = identity(Chem.MolFromSmiles(smiles))
    except Exception as error:
        return None, f"{type(error).__name__}: {error}"
    return (actual, None) if actual is not None else (None, "unreadable SMILES")


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


def rdkit_block(molecule: Chem.Mol) -> str | None:
    try:
        return Chem.MolToMolBlock(molecule, forceV3000=True)
    except (RuntimeError, ValueError):
        return None


def cases(corpus: Path):
    for line_number, line in enumerate(corpus.read_text().splitlines(), 1):
        if not line.strip():
            continue
        if corpus.suffix == ".smi":
            yield {"id": f"line-{line_number}", "smiles": line.split()[0]}
        elif corpus.suffix == ".jsonl":
            case = json.loads(line)
            if case.get("_manifest") is True:
                continue
            yield case
        else:
            raise ValueError("corpus must be .smi or .jsonl")


def compare(corpus: Path, cli: Path) -> dict:
    params = cosmolkit.MolBlockWriteParams(format=cosmolkit.SdfFormat.V3000)
    rows = []
    for case in cases(corpus):
        source = case["smiles"]
        molecule = Chem.MolFromSmiles(source)
        if molecule is None:
            rows.append({"id": case["id"], "source_smiles": source,
                         "excluded": "RDKit cannot parse source"})
            continue
        expected = identity(molecule)
        assert expected is not None
        # The original text can denote a different graph under different
        # parsers. Give all three writers the same RDKit-judged molecule.
        comparison_molecule = Chem.MolFromSmiles(expected)
        if comparison_molecule is None or identity(comparison_molecule) != expected:
            rows.append({"id": case["id"], "source_smiles": source,
                         "comparison_smiles": expected,
                         "excluded": "RDKit-normalized SMILES does not preserve identity"})
            continue
        try:
            cosmolkit_block = cosmolkit.Molecule.from_smiles(expected).to_mol_with_params(params)
        except Exception:
            cosmolkit_block = None
        blocks = {
            "rdkit": rdkit_block(comparison_molecule),
            "cosmolkit": cosmolkit_block,
            "chematic": chematic_block(cli, expected),
        }
        cosmolkit_read, cosmolkit_error = cosmolkit_readback(blocks["chematic"])
        rows.append(
            {
                "id": case["id"],
                "source_smiles": source,
                "comparison_smiles": expected,
                "expected": expected,
                "source_has_stereo_token": any(token in expected for token in "@/\\"),
                "identities": {
                    **{engine: readback(block) for engine, block in blocks.items()},
                    "chematic_cosmolkit_reader": cosmolkit_read,
                },
                "readback_errors": {"chematic_cosmolkit_reader": cosmolkit_error},
                "block_sha256": {
                    engine: hashlib.sha256(block.encode()).hexdigest() if block else None
                    for engine, block in blocks.items()
                },
            }
        )
    comparable = [row for row in rows if "identities" in row]
    summary = {"cases": len(rows), "comparable": len(comparable),
               "excluded": len(rows) - len(comparable)}
    for engine in ("rdkit", "cosmolkit", "chematic", "chematic_cosmolkit_reader"):
        summary[engine] = {
            "exact": sum(row["identities"][engine] == row["expected"] for row in comparable),
            "stereo_exact": sum(
                row["source_has_stereo_token"]
                and row["identities"][engine] == row["expected"]
                for row in comparable
            ),
            "failed_ids": [
                row["id"] for row in comparable if row["identities"][engine] != row["expected"]
            ],
        }
    return {
        "scope": "development RDKit-normalized SMILES to V3000 to RDKit read-back; atom maps ignored",
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
