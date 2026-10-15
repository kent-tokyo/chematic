#!/usr/bin/env python3
"""Compare default V3000 stereo round trips using an independent RDKit reader.

Both writers receive the same SMILES. Open Babel is measured with its default
settings and with --gen2D; CheMatic uses its default V3000 conversion. Atom
map numbers are removed before comparing chemical identities, because map
retention is a separate format-metadata question.
"""

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase


ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "validation/stereo_torture_suite_development.jsonl"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command(args: list[str], source: str | None = None) -> str | None:
    try:
        result = subprocess.run(
            args, input=source, capture_output=True, text=True, timeout=30, check=False
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return result.stdout if result.returncode == 0 else None


def identity(molecule: Chem.Mol | None) -> str | None:
    if molecule is None:
        return None
    molecule = Chem.Mol(molecule)
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def mol_identity(block: str | None) -> str | None:
    if block is None:
        return None
    try:
        return identity(Chem.MolFromMolBlock(block, removeHs=False))
    except (ValueError, RuntimeError):
        return None


def source_rows(corpus: Path) -> list[dict]:
    rows = [json.loads(line) for line in corpus.read_text().splitlines() if line]
    return [row for row in rows if row.get("smiles")]


def compare(corpus: Path, cli: Path) -> dict:
    rows = []
    with tempfile.TemporaryDirectory() as directory:
        mol_path = Path(directory) / "case.mol"
        for case in source_rows(corpus):
            smiles = case["smiles"]
            source = identity(Chem.MolFromSmiles(smiles))
            if source is None:
                raise ValueError(f"RDKit cannot parse {case.get('id', smiles)}")
            chematic_block = command(
                [
                    str(cli),
                    "convert",
                    "--input-format",
                    "smiles",
                    "--output-format",
                    "mol_v3000",
                ],
                smiles,
            )
            babel_identities = []
            babel_hashes = []
            for options in ([], ["--gen2D"]):
                mol_path.unlink(missing_ok=True)
                result = command(
                    [
                        "obabel",
                        "-:" + smiles,
                        "-omol",
                        "-O",
                        str(mol_path),
                        "-x3",
                        *options,
                    ]
                )
                block = (
                    mol_path.read_text()
                    if result is not None and mol_path.is_file()
                    else None
                )
                babel_identities.append(mol_identity(block))
                babel_hashes.append(
                    hashlib.sha256(block.encode()).hexdigest() if block else None
                )
            rows.append(
                {
                    "id": case.get("id"),
                    "source_smiles": smiles,
                    "source_identity": source,
                    "source_has_stereo_token": any(token in source for token in "@/\\"),
                    "chematic_identity": mol_identity(chematic_block),
                    "openbabel_default_identity": babel_identities[0],
                    "openbabel_gen2d_identity": babel_identities[1],
                    "chematic_mol_sha256": (
                        hashlib.sha256(chematic_block.encode()).hexdigest()
                        if chematic_block
                        else None
                    ),
                    "openbabel_default_mol_sha256": babel_hashes[0],
                    "openbabel_gen2d_mol_sha256": babel_hashes[1],
                }
            )
    summary = {
        "cases": len(rows),
        "source_with_stereo_token": sum(row["source_has_stereo_token"] for row in rows),
    }
    for name in ("chematic", "openbabel_default", "openbabel_gen2d"):
        summary[f"{name}_identity_matches"] = sum(
            row[f"{name}_identity"] == row["source_identity"] for row in rows
        )
        summary[f"{name}_stereo_identity_matches"] = sum(
            row["source_has_stereo_token"]
            and row[f"{name}_identity"] == row["source_identity"]
            for row in rows
        )
    return {
        "schema_version": 1,
        "scope": "development corpus; default SMILES to V3000 writers, RDKit MOL reread, atom-map-independent isomeric identity",
        "tools": {
            "rdkit": rdBase.rdkitVersion,
            "openbabel": (command(["obabel", "-V"]) or "").strip(),
            "os": platform.system(),
            "architecture": platform.machine(),
            "python": sys.version.split()[0],
            "chematic_source_revision": (
                command(["git", "-C", str(ROOT), "rev-parse", "HEAD"]) or ""
            ).strip(),
            "chematic_cli": str(cli),
            "chematic_cli_sha256": digest(cli),
            "corpus": str(corpus),
            "corpus_sha256": digest(corpus),
        },
        "summary": summary,
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")
    result = compare(args.corpus, args.cli)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["summary"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
