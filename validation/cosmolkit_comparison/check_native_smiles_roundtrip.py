#!/usr/bin/env python3
"""Score native SMILES read-back against pinned RDKit and COSMolKit.

Build the current source first:
  cargo build -p chematic-smiles --example canonical_roundtrip_dump --offline

Then run with RDKit 2025.09.3 and COSMolKit 0.5.0:
  python check_native_smiles_roundtrip.py --corpus corpora/stress-v1.smi \
      --chematic ../../target/debug/examples/canonical_roundtrip_dump

The oracle is RDKit's canonical isomeric SMILES of each source molecule.
Each engine's native output is parsed by RDKit before comparison. This is a
strict string comparison after read-back, including non-tetrahedral tags;
symmetry-equivalent ligand permutations can therefore count as differences.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import cosmolkit
from rdkit import Chem, RDLogger


def normalized(smiles: str | None) -> str | None:
    if smiles is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol is not None else None


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
    lines = corpus.read_text().splitlines()
    completed = subprocess.run([str(binary), str(corpus)], text=True, capture_output=True, check=True)
    output = {}
    for line in completed.stdout.splitlines():
        index, smiles = line.split("\t", 1)
        output[int(index)] = smiles

    eligible = 0
    skipped = 0
    matched = {"chematic": 0, "cosmolkit": 0}
    misses = {"chematic": [], "cosmolkit": []}
    for index, line in enumerate(lines):
        tokens = line.split()
        if not tokens:
            continue
        source = tokens[0]
        expected = normalized(source)
        if expected is None:
            skipped += 1
            continue
        eligible += 1
        try:
            cosmolkit_smiles = cosmolkit.Molecule.from_smiles(source).to_smiles()
        except Exception:
            cosmolkit_smiles = None
        for engine, smiles in (
            ("chematic", output.get(index)),
            ("cosmolkit", cosmolkit_smiles),
        ):
            actual = normalized(smiles) if smiles and not smiles.startswith("ERROR:") else None
            if actual == expected:
                matched[engine] += 1
            else:
                misses[engine].append({"index": index, "expected": expected, "actual": actual})
    print(json.dumps({
        "corpus_sha256": hashlib.sha256(corpus.read_bytes()).hexdigest(),
        "chematic_binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "cosmolkit_version": cosmolkit.__version__,
        "eligible": eligible,
        "skipped": skipped,
        "matched": matched,
        "mismatches": {engine: rows[:20] for engine, rows in misses.items()},
    }, sort_keys=True))


if __name__ == "__main__":
    main()
