#!/usr/bin/env python3
"""Compare 3D MOL2 identity from identical Open Babel-written input bytes.

RDKit embeds and source-gates each molecule before Open Babel writes MOL2.
The resulting MOL2 is read by both CheMatic and Open Babel. This is a
development comparison with a common RDKit identity judge, not chemical truth.
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CORPUS = ROOT / "validation/cosmolkit_comparison/corpora/stress-v1.smi"
BRIDGE_SOURCES = [
    "C1C[C@H]2CC[C@@H]1C2",
    "[C@@H]12CC[C@@H](CC1)O2",
    "C1=C[C@@H]2C[C@H]1C=C2",
]
ALLOWED_ATOMS = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command(args, timeout=20):
    try:
        return subprocess.run(args, capture_output=True, text=True, timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return None


def identity(smiles):
    mol = Chem.MolFromSmiles(smiles.split()[0]) if smiles.strip() else None
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol is not None else None


def stress_cases(corpus, limit=160):
    seen = set()
    selected = []
    for line_number, line in enumerate(corpus.read_text().splitlines(), 1):
        source = line.split()[0] if line.strip() else ""
        if "@" not in source:
            continue
        mol = Chem.MolFromSmiles(source)
        if mol is None or not 4 <= mol.GetNumHeavyAtoms() <= 20:
            continue
        if any(atom.GetAtomicNum() not in ALLOWED_ATOMS for atom in mol.GetAtoms()):
            continue
        expected = Chem.MolToSmiles(mol, isomericSmiles=True)
        if expected in seen:
            continue
        seen.add(expected)
        selected.append((f"stress-line-{line_number}", source, 20261017 + len(selected)))
        if limit is not None and len(selected) == limit:
            break
    return selected


def bridge_holdout_cases():
    return [
        (f"bridge-{index}-seed-{seed}", source, seed)
        for index, source in enumerate(BRIDGE_SOURCES)
        for seed in range(500000 + index * 100, 500010 + index * 100)
    ]


def compare_case(case_id, source, seed, cli, mol_path, mol2_path):
    expected = identity(source)
    molecule = Chem.AddHs(Chem.MolFromSmiles(source))
    try:
        embedded = AllChem.EmbedMolecule(molecule, randomSeed=seed, maxAttempts=20)
    except (RuntimeError, ValueError):
        embedded = -1
    row = {"id": case_id, "source": source, "seed": seed, "expected": expected}
    if embedded != 0:
        row["excluded"] = "RDKit embedding failed"
        return row
    mol_path.write_text(Chem.MolToMolBlock(molecule, forceV3000=False))
    reread = Chem.MolFromMolFile(str(mol_path), removeHs=True)
    if reread is None or Chem.MolToSmiles(reread, isomericSmiles=True) != expected:
        row["excluded"] = "RDKit V2000 read-back changed source identity"
        return row
    mol2_path.unlink(missing_ok=True)
    converted = command(["obabel", "-imol", str(mol_path), "-omol2", "-O", str(mol2_path)])
    if converted is None or converted.returncode != 0 or not mol2_path.is_file():
        row["excluded"] = "Open Babel MOL2 writer failed"
        return row
    chematic = command([
        str(cli), "convert", "--input-format", "mol2", "--output-format", "smiles",
        "--input", str(mol2_path),
    ])
    openbabel = command(["obabel", "-imol2", str(mol2_path), "-osmi"])
    row.update({
        "v2000_sha256": sha256(mol_path),
        "mol2_sha256": sha256(mol2_path),
        "identities": {
            "chematic": identity(chematic.stdout) if chematic and chematic.returncode == 0 else None,
            "openbabel": identity(openbabel.stdout) if openbabel and openbabel.returncode == 0 else None,
        },
    })
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--bridge-holdout", action="store_true")
    parser.add_argument(
        "--all-stereo", action="store_true",
        help="include every eligible distinct stereo identity instead of the first 160",
    )
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")
    obabel = shutil.which("obabel")
    if obabel is None:
        parser.error("Open Babel executable `obabel` is required")
    cli = args.cli.resolve(strict=True)
    corpus = args.corpus.resolve(strict=True)
    cases = bridge_holdout_cases() if args.bridge_holdout else stress_cases(
        corpus, None if args.all_stereo else 160)
    with tempfile.TemporaryDirectory() as directory:
        directory = Path(directory)
        rows = [
            compare_case(case_id, source, seed, cli, directory / "source.mol", directory / "source.mol2")
            for case_id, source, seed in cases
        ]
    comparable = [row for row in rows if "identities" in row]
    summary = {"selected": len(rows), "comparable": len(comparable),
               "excluded": len(rows) - len(comparable)}
    for engine in ("chematic", "openbabel"):
        summary[engine] = sum(row["identities"][engine] == row["expected"] for row in comparable)
    summary["openbabel_only"] = sum(
        row["identities"]["openbabel"] == row["expected"]
        and row["identities"]["chematic"] != row["expected"]
        for row in comparable
    )
    report = {
        "scope": "development 3D MOL2 identity; RDKit V2000 source gate; same MOL2 bytes",
        "rdkit_version": rdBase.rdkitVersion,
        "openbabel_version": command([obabel, "-V"]).stdout.strip(),
        "openbabel_sha256": sha256(Path(obabel)),
        "cli_sha256": sha256(cli),
        "corpus_sha256": sha256(corpus) if not args.bridge_holdout else None,
        "bridge_holdout": args.bridge_holdout,
        "all_stereo": args.all_stereo,
        "summary": summary,
        "rows": rows,
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
