#!/usr/bin/env python3
"""Emit common-schema smoke results for chematic or RDKit."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORPUS = Path(__file__).with_name("smoke_corpus.jsonl")


def corpus_hash() -> str:
    return hashlib.sha256(CORPUS.read_bytes()).hexdigest()


def source_commit() -> str | None:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def result_base(engine: str, version: str, row: dict) -> dict:
    return {"schema_version": 1, "engine": engine, "engine_version": version,
            "source_commit": source_commit(), "corpus_sha256": corpus_hash(),
            "id": row["id"], "smiles": row["smiles"], "status": "ok", "operations": {}}


def ok(value):
    return {"status": "ok", "value": value}


def run_rdkit(rows):
    from rdkit import Chem, rdBase
    from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMolDescriptors
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    for row in rows:
        out = result_base("rdkit", rdBase.rdkitVersion, row)
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            out["status"] = "parse_error"
        else:
            out["operations"] = {
                "canonical_smiles": ok(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)),
                "formula": ok(rdMolDescriptors.CalcMolFormula(mol)),
                "molecular_weight": ok(Descriptors.MolWt(mol)),
                "rdkit_morgan_bits": ok(sorted(generator.GetFingerprint(mol).GetOnBits())),
            }
        yield out


def run_chematic(rows):
    import chematic
    version = getattr(chematic, "__version__", "unknown")
    for row in rows:
        out = result_base("chematic", version, row)
        try:
            mol = chematic.from_smiles(row["smiles"])
        except Exception as exc:  # the contract preserves the failure category
            out["status"] = "parse_error"
            out["operations"] = {"parse": {"status": "parse_error", "error": str(exc)}}
        else:
            out["operations"] = {
                "canonical_smiles": ok(mol.smiles),
                "formula": ok(mol.formula),
                "molecular_weight": ok(mol.mw),
            }
            if hasattr(mol, "rdkit_ecfp_config"):
                fp = mol.rdkit_ecfp_config(2, 2048)
                out["operations"]["rdkit_morgan_bits"] = ok(
                    [i for i in range(len(fp) * 8) if fp[i // 8] & (1 << (i % 8))]
                )
            else:
                out["operations"]["rdkit_morgan_bits"] = {
                    "status": "unsupported",
                    "error": "requires chematic v0.26.0 or newer",
                }
        yield out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", choices=["rdkit", "chematic"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = [json.loads(line) for line in CORPUS.read_text().splitlines() if line.strip()]
    records = run_rdkit(rows) if args.engine == "rdkit" else run_chematic(rows)
    args.output.write_text("".join(json.dumps(record, sort_keys=True) + "\n" for record in records))


if __name__ == "__main__":
    main()
