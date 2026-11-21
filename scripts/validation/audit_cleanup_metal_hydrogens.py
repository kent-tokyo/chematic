#!/usr/bin/env python3
"""Verify formula preservation on the 47 Cleanup rows and metal-H controls."""

import argparse
import hashlib
import json
from pathlib import Path
import zipfile

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[2]
RECORD = ROOT / "validation/results/rdkit-cleanup-v1.2.2-source-stress-2026-11-21.json"
CORPUS = ROOT / "validation/cosmolkit_comparison/corpora/stress-v1.smi"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wheel", required=True, type=Path)
    parser.add_argument("--candidate-commit", required=True)
    args = parser.parse_args()
    if chematic.__version__ != "1.2.2" or rdBase.rdkitVersion != "2026.09.1":
        raise ValueError("expected chematic 1.2.2 and RDKit 2026.09.1")
    RDLogger.DisableLog("rdApp.*")
    package_root = Path(chematic.__file__).parent.parent
    with zipfile.ZipFile(args.wheel) as archive:
        files = [name for name in archive.namelist() if name.startswith("chematic/") and not name.endswith("/")]
        if not files:
            raise ValueError("wheel contains no chematic package files")
        for name in files:
            path = package_root / name
            if not path.is_file() or path.read_bytes() != archive.read(name):
                raise ValueError(f"installed package differs from wheel: {path}")
    cases = {row["index"] for row in json.loads(RECORD.read_text())["mismatches"]}
    cases.update((1294, 3075, 3402))
    corpus = CORPUS.read_text().splitlines()
    rows = []
    for index in sorted(cases):
        text = corpus[index].split()[0]
        source = chematic.from_smiles(text)
        cleaned = source.cleanup()
        output = cleaned.rdkit_smiles
        reparsed = chematic.from_smiles(output)
        repeated = cleaned.cleanup()
        rdkit_input = Chem.MolFromSmiles(text)
        rdkit_output = rdMolStandardize.Cleanup(rdkit_input)
        charges = [sum(atom[2] for atom in mol.atom_table)
                   for mol in (source, cleaned, reparsed, repeated)]
        rows.append({"index": index, "input_formula": source.formula,
                     "output_formula": cleaned.formula, "reparsed_formula": reparsed.formula,
                     "repeated_formula": repeated.formula, "formal_charges": charges,
                     "rdkit_input_formula": rdMolDescriptors.CalcMolFormula(rdkit_input),
                     "rdkit_output_formula": rdMolDescriptors.CalcMolFormula(rdkit_output),
                     "output": output})
    report = {"candidate_commit": args.candidate_commit, "wheel": args.wheel.name,
              "rdkit_version": rdBase.rdkitVersion,
              "wheel_sha256": hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
              "installed_package_files_verified": len(files),
              "corpus_sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
              "source_record": str(RECORD.relative_to(ROOT)), "rows": rows}
    path = ROOT / "validation/results/cleanup-metal-hydrogen-wheel-2026-11-21.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    failures = [row["index"] for row in rows if len({row["input_formula"], row["output_formula"],
                                                      row["reparsed_formula"], row["repeated_formula"]}) != 1
                or len(set(row["formal_charges"])) != 1]
    if failures:
        raise ValueError(f"formula differs on rows {failures}")
    print(f"Formula and formal charge preserved across input, Cleanup, round trip and repeat for {len(rows)} rows")


if __name__ == "__main__":
    main()
