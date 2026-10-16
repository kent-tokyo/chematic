#!/usr/bin/env python3
"""Check preservation of RDKit-enumerated aromatic sulfur stereo SMILES.

Build `cargo build -p chematic-smiles --example canonical_roundtrip_dump --offline`
and run this script with RDKit 2025.09.3 and COSMolKit 0.5.0 installed.
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
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

ROWS = (4385, 4471, 4498, 4503, 4528)
CORPUS_SHA256 = "f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f"


def chematic_roundtrip(binary: Path, smiles: list[str]) -> list[str]:
    with tempfile.TemporaryDirectory(prefix="chematic-aromatic-s-stereo-") as directory:
        corpus = Path(directory) / "stereo.smi"
        corpus.write_text("\n".join(smiles) + "\n")
        completed = subprocess.run(
            [str(binary), str(corpus)], capture_output=True, text=True, check=True
        )
    rows = [line.split("\t", 1) for line in completed.stdout.splitlines()]
    if [int(index) for index, _ in rows] != list(range(len(smiles))):
        raise AssertionError("CheMatic output indices or row count differ")
    if any(value.startswith("ERROR:") for _, value in rows):
        raise AssertionError("CheMatic rejected an enumerated stereo SMILES")
    return [value for _, value in rows]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--corpus", required=True, type=Path)
    parser.add_argument("--chematic", required=True, type=Path)
    args = parser.parse_args()
    if Chem.rdBase.rdkitVersion != "2025.09.3" or cosmolkit.__version__ != "0.5.0":
        raise SystemExit("requires RDKit 2025.09.3 and COSMolKit 0.5.0")
    RDLogger.DisableLog("rdApp.*")
    corpus = args.corpus.resolve()
    binary = args.chematic.resolve()
    corpus_hash = hashlib.sha256(corpus.read_bytes()).hexdigest()
    if corpus_hash != CORPUS_SHA256:
        raise SystemExit(f"corpus SHA-256 differs: {corpus_hash}")
    lines = corpus.read_text().splitlines()
    variants: list[str] = []
    sources: list[str] = []
    for index in ROWS:
        source = lines[index].split()[0]
        enumerated = sorted(
            Chem.MolToSmiles(mol)
            for mol in EnumerateStereoisomers(Chem.MolFromSmiles(source))
        )
        if len(enumerated) != 2 or not all("[s@" in s for s in enumerated):
            raise AssertionError(f"row {index}: expected two tagged aromatic sulfur isomers")
        sources.append(source)
        variants.extend(enumerated)

    chematic = chematic_roundtrip(binary, variants)
    chematic_again = chematic_roundtrip(binary, chematic)
    rdkit = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in variants]
    cosmolkit_result = [cosmolkit.Molecule.from_smiles(s).to_smiles() for s in variants]
    results = {}
    for name, output in (
        ("chematic", chematic), ("rdkit", rdkit), ("cosmolkit", cosmolkit_result)
    ):
        results[name] = sum(
            output[2 * i] != output[2 * i + 1] for i in range(len(ROWS))
        )
    if results["chematic"] != len(ROWS) or chematic != chematic_again:
        raise AssertionError("CheMatic failed to preserve distinct, stable stereo variants")
    print(json.dumps({
        "corpus_sha256": corpus_hash,
        "chematic_binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "cosmolkit_version": cosmolkit.__version__,
        "source_indices": ROWS,
        "source_smiles": sources,
        "pairs": len(ROWS),
        "distinct_pairs_after_read_write": results,
        "chematic_idempotent": True,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
