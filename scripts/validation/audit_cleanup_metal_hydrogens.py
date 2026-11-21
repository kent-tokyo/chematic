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
    parser.add_argument("--output", type=Path, required=True)
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
        try:
            cleaned = source.cleanup()
        except chematic.ChematicInputError as error:
            if error.category != "ambiguous" or error.code != "cleanup_metal_hydrogen":
                raise
            rows.append({"index": index, "input_formula": source.formula,
                         "refused": {"category": error.category, "code": error.code}})
            continue
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
    spelling_pairs = []
    for label, inputs, expect_error in (
        ("ambiguous_tin", ("C/N=C(\\S)N/N=C/C(C)=[N]/[SnH]=[CH]/N",
                           "C/N=C(\\S)N/N=C/C(C)=[N]/[Sn]([H])=[CH]/N"), True),
        ("resolved_platinum", ("C[PtH]C", "C[Pt]([H])C"), False),
    ):
        variants = []
        for text in inputs:
            source = chematic.from_smiles(text)
            entry = {"input": text, "input_smiles": source.rdkit_smiles,
                     "input_formula": source.formula}
            try:
                cleaned = source.cleanup()
                entry.update(output=cleaned.rdkit_smiles, output_formula=cleaned.formula)
            except chematic.ChematicInputError as error:
                entry["error"] = {"category": error.category, "code": error.code}
            variants.append(entry)
        if (len({row["input_smiles"] for row in variants}) != 1
                or len({row["input_formula"] for row in variants}) != 1):
            raise ValueError(f"inputs are not equivalent for {label}")
        if expect_error:
            if any(row.get("error") != {"category": "ambiguous", "code": "cleanup_metal_hydrogen"}
                   for row in variants):
                raise ValueError(f"inconsistent refusal for {label}")
        elif (any("error" in row or row["output_formula"] != row["input_formula"] for row in variants)
              or len({row["output"] for row in variants}) != 1):
            raise ValueError(f"inconsistent accepted output for {label}")
        spelling_pairs.append({"case": label, "variants": variants})
    report = {"candidate_commit": args.candidate_commit, "wheel": args.wheel.name,
              "rdkit_version": rdBase.rdkitVersion,
              "wheel_sha256": hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
              "installed_package_files_verified": len(files),
              "corpus_sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
              "source_record": str(RECORD.relative_to(ROOT)), "rows": rows,
              "explicit_hydrogen_spelling_pairs": spelling_pairs}
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    failures = [row["index"] for row in rows if "refused" not in row and (len({row["input_formula"], row["output_formula"],
                                                      row["reparsed_formula"], row["repeated_formula"]}) != 1
                or len(set(row["formal_charges"])) != 1)]
    if failures:
        raise ValueError(f"formula differs on rows {failures}")
    print(f"Formula and charge preserved for {sum('refused' not in row for row in rows)} accepted rows; refused {[row['index'] for row in rows if 'refused' in row]}")


if __name__ == "__main__":
    main()
