#!/usr/bin/env python3
"""Compare Cleanup on recorded inputs and equivalent explicit-H graphs."""

import argparse
import hashlib
import json
from pathlib import Path
import zipfile

import chematic
from rdkit import Chem, RDLogger, rdBase

ROOT = Path(__file__).resolve().parents[2]
RECORD = ROOT / "validation/results/rdkit-cleanup-v1.2.2-source-stress-2026-11-21.json"
CORPUS = ROOT / "validation/cosmolkit_comparison/corpora/stress-v1.smi"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wheel", required=True, type=Path)
    parser.add_argument("--candidate-commit", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if chematic.__version__ != "1.2.2" or rdBase.rdkitVersion != "2026.09.1":
        raise ValueError("expected chematic 1.2.2 and RDKit 2026.09.1")
    package_root = Path(chematic.__file__).parent.parent
    with zipfile.ZipFile(args.wheel) as archive:
        files = [n for n in archive.namelist() if n.startswith("chematic/") and not n.endswith("/")]
        if not files or any(not (package_root / n).is_file() or
                            (package_root / n).read_bytes() != archive.read(n) for n in files):
            raise ValueError("installed package differs from wheel")
    RDLogger.DisableLog("rdApp.*")
    cases = {row["index"] for row in json.loads(RECORD.read_text())["mismatches"]}
    cases.update((1294, 3075, 3402))
    corpus = CORPUS.read_text().splitlines()
    rows = []
    for index in sorted(cases):
        original = corpus[index].split()[0]
        reference = Chem.MolFromSmiles(original)
        if reference is None:
            raise ValueError(f"RDKit could not parse row {index}")
        explicit = Chem.MolToSmiles(Chem.AddHs(reference), canonical=False, isomericSmiles=True)
        reparsed = Chem.MolFromSmiles(explicit)
        if reparsed is None or Chem.MolToSmiles(reparsed, isomericSmiles=True) != Chem.MolToSmiles(reference, isomericSmiles=True):
            raise ValueError(f"explicit H graph is not equivalent on row {index}")
        variants = []
        for smiles in (original, explicit):
            source = chematic.from_smiles(smiles)
            row = {"input": smiles, "input_formula": source.formula}
            try:
                cleaned = source.cleanup()
                row.update(output=cleaned.rdkit_smiles, output_formula=cleaned.formula)
            except chematic.ChematicInputError as error:
                row["error"] = {"category": error.category, "code": error.code}
            variants.append(row)
        accepted = all("error" not in x for x in variants)
        refused = all(x.get("error") == {"category": "ambiguous", "code": "cleanup_metal_hydrogen"} for x in variants)
        matches = (variants[0]["input_formula"] == variants[1]["input_formula"] and
                   ((accepted and variants[0]["output"] == variants[1]["output"] and
                     all(x["output_formula"] == x["input_formula"] for x in variants)) or refused))
        rows.append({"index": index, "matches": matches, "variants": variants})
    report = {"candidate_commit": args.candidate_commit, "wheel": args.wheel.name,
              "wheel_sha256": digest(args.wheel), "rdkit_version": rdBase.rdkitVersion,
              "installed_package_files_verified": len(files), "corpus_sha256": digest(CORPUS),
              "source_record_sha256": digest(RECORD), "method": "RDKit AddHs and noncanonical isomeric SMILES; verify RDKit reparse, compare chematic formulas, output and typed refusals", "rows": rows}
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    failures = [row["index"] for row in rows if not row["matches"]]
    if failures:
        raise ValueError(f"explicit-H variants differ on rows {failures}")
    print(f"{len(rows)} explicit-H pairs have matching formulas and Cleanup outcomes")


if __name__ == "__main__":
    main()
