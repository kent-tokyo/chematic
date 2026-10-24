"""Compare the CLI stereo-search and named fingerprint profiles with RDKit.

Usage: python scripts/validation/rdkit_cli_stereo_fingerprint_gate.py \
    --binary target/debug/chematic \
    --stereo-self-corpus scripts/fixtures/descriptor_census_corpus.smi

The optional corpus checks stereo-bearing SMILES against themselves as SMARTS.
This is a regression gate, not a toolkit-wide accuracy measurement. Run after
building the current CLI binary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator


STEREO_SMILES = (
    "F/C=C/F",
    "F/C=C\\F",
    "C/C=C/C",
    "C/C=C\\C",
    "N[C@@H](C)C(=O)O",
    "N[C@H](C)C(=O)O",
    "F[C@H](Cl)Br",
    "F[C@@H](Cl)Br",
    "C1[C@H](F)CC1",
    "C1[C@@H](F)CC1",
    "F/C(Cl)=C/F",
    "F/C(Cl)=C\\F",
    "C/C=C(/C)C",
    "C/C=C(\\C)C",
    "C/C=C/C=C/C",
    "C/C=C\\C=C/C",
    "C1CCCCC/C=C/1",
    "C1CCCCC/C=C\\1",
)

FINGERPRINT_SMILES = (
    "CCO",
    "CCN",
    "c1ccccc1",
    "CC(=O)O",
    "N[C@@H](C)C(=O)O",
    "F/C=C/F",
    "F/C=C\\F",
    "C1CCCCC1",
    "O=C(O)c1ccccc1",
    "C#N",
)


def cli_json(binary: Path, args: list[str], stdin: str) -> dict:
    result = subprocess.run(
        [str(binary), *args],
        input=stdin,
        text=True,
        capture_output=True,
        check=True,
    )
    return json.loads(result.stdout)


def run(binary: Path, stereo_self_corpus: Path | None = None) -> dict:
    stereo_pairs = [
        (target, query)
        for target in STEREO_SMILES
        for query in STEREO_SMILES
    ]
    stereo_input = "".join(f"{target}\t{query}\n" for target, query in stereo_pairs)
    stereo_rows = cli_json(binary, ["batch-substructure", "--use-chirality"], stereo_input)[
        "records"
    ]
    stereo_mismatches = []
    for (target, query), row in zip(stereo_pairs, stereo_rows, strict=True):
        expected = Chem.MolFromSmiles(target).HasSubstructMatch(
            Chem.MolFromSmarts(query), useChirality=True
        )
        actual = (
            bool(row["substructure"]["match_count"])
            if row["substructure"] is not None
            else None
        )
        if actual != expected or row["error"] is not None:
            stereo_mismatches.append(
                {"target": target, "query": query, "rdkit": expected, "chematic": actual,
                 "error": row["error"]}
            )

    fp_mismatches = []
    morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    for algorithm in ("rdkit_morgan", "rdkit_pattern"):
        rows = cli_json(
            binary,
            ["batch-fingerprints", "--algorithm", algorithm],
            "".join(f"{smiles}\n" for smiles in FINGERPRINT_SMILES),
        )["records"]
        for smiles, row in zip(FINGERPRINT_SMILES, rows, strict=True):
            molecule = Chem.MolFromSmiles(smiles)
            expected = list(
                (morgan.GetFingerprint(molecule) if algorithm == "rdkit_morgan"
                 else Chem.PatternFingerprint(molecule, fpSize=2048)).GetOnBits()
            )
            actual = (
                row["fingerprint"]["set_bits"] if row["fingerprint"] is not None else None
            )
            if actual != expected or row["error"] is not None:
                fp_mismatches.append(
                    {"smiles": smiles, "algorithm": algorithm, "rdkit_bits": expected,
                     "chematic_bits": actual, "error": row["error"]}
                )

    self_cases = []
    self_mismatches = []
    corpus_sha256 = None
    if stereo_self_corpus is not None:
        corpus_bytes = stereo_self_corpus.read_bytes()
        corpus_sha256 = hashlib.sha256(corpus_bytes).hexdigest()
        self_cases = [
            line for line in corpus_bytes.decode("utf-8").splitlines()
            if line and not line.startswith("#")
            and ("@" in line or "/" in line or "\\" in line)
        ]
        rows = cli_json(
            binary,
            ["batch-substructure", "--use-chirality"],
            "".join(f"{smiles}\t{smiles}\n" for smiles in self_cases),
        )["records"]
        for index, (smiles, row) in enumerate(zip(self_cases, rows, strict=True)):
            target = Chem.MolFromSmiles(smiles)
            query = Chem.MolFromSmarts(smiles)
            if target is None or query is None:
                raise ValueError(f"RDKit could not parse corpus row {index}: {smiles}")
            expected = target.HasSubstructMatch(query, useChirality=True)
            actual = (
                bool(row["substructure"]["match_count"])
                if row["substructure"] is not None else None
            )
            if actual != expected or row["error"] is not None:
                self_mismatches.append({
                    "index": index, "smiles": smiles, "rdkit": expected,
                    "chematic": actual, "error": row["error"],
                })

    return {
        "rdkit_version": rdBase.rdkitVersion,
        "stereo_search_pairs": len(stereo_pairs),
        "stereo_search_mismatches": stereo_mismatches,
        "fingerprint_cases": len(FINGERPRINT_SMILES) * 2,
        "fingerprint_mismatches": fp_mismatches,
        "stereo_self_corpus_sha256": corpus_sha256,
        "stereo_self_cases": len(self_cases),
        "stereo_self_mismatches": self_mismatches,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--stereo-self-corpus", type=Path)
    args = parser.parse_args()
    result = run(args.binary.resolve(), args.stereo_self_corpus)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return int(bool(result["stereo_search_mismatches"] or result["fingerprint_mismatches"]
                    or result["stereo_self_mismatches"]))


if __name__ == "__main__":
    raise SystemExit(main())
