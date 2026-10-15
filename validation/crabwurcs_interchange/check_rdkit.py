"""Compare generated CTfiles with source SMILES using an independent reader."""

import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase


TABLES = ("glycoshape_notations.tsv", "glycoshape_derived_notations.tsv")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", type=Path, help="crabWURCS data directory")
    parser.add_argument("output_dir", type=Path, help="directory containing generated SDFs")
    parser.add_argument("--summary", type=Path, help="write a JSON summary here")
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")

    source = []
    for table in TABLES:
        for line_number, line in enumerate((args.data_dir / table).read_text().splitlines(), 1):
            fields = line.split("\t")
            source.append((f"{table}:{line_number}", fields[4] if len(fields) > 4 else ""))

    summary = {"rdkit": rdBase.rdkitVersion, "source_rows": len(source), "formats": {}}
    failed = False
    for format_name in ("v2000", "v3000"):
        path = args.output_dir / f"{format_name}.sdf"
        supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=True)
        counts = {
            "records": len(supplier),
            "source_invalid": 0,
            "mol_invalid": 0,
            "connectivity_diff": 0,
            "stereo_diff": 0,
            "same": 0,
        }
        examples = []
        if len(supplier) != len(source):
            failed = True
        for (label, smiles), mol in zip(source, supplier):
            original = Chem.MolFromSmiles(smiles)
            if original is None:
                kind = "source_invalid"
            elif mol is None:
                kind = "mol_invalid"
            else:
                before = Chem.MolToSmiles(original, isomericSmiles=True)
                after = Chem.MolToSmiles(mol, isomericSmiles=True)
                if before == after:
                    counts["same"] += 1
                    continue
                original_no_stereo = Chem.Mol(original)
                mol_no_stereo = Chem.Mol(mol)
                Chem.RemoveStereochemistry(original_no_stereo)
                Chem.RemoveStereochemistry(mol_no_stereo)
                kind = (
                    "stereo_diff"
                    if Chem.MolToSmiles(original_no_stereo) == Chem.MolToSmiles(mol_no_stereo)
                    else "connectivity_diff"
                )
            counts[kind] += 1
            if len(examples) < 20:
                examples.append({"row": label, "kind": kind})
        if any(counts[key] for key in ("source_invalid", "mol_invalid", "connectivity_diff", "stereo_diff")):
            failed = True
        summary["formats"][format_name] = {
            "counts": counts,
            "examples": examples,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }

    result = json.dumps(summary, indent=2) + "\n"
    if args.summary:
        args.summary.write_text(result)
    print(result, end="")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
