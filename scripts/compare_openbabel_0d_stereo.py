#!/usr/bin/env python3
"""Compare source-identity retention for Open Babel 0D V2000 MOL files.

The source SMILES is independent of the MOL re-read. Only cases where Open
Babel's own re-read agrees with the RDKit source identity are scored. Both
CheMatic binaries receive the identical Open Babel-generated MOL bytes.
"""

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem, rdBase


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(*args: str) -> str | None:
    process = subprocess.run(args, capture_output=True, text=True, timeout=30)
    return process.stdout.strip() if process.returncode == 0 else None


def identity(smiles: str | None) -> str | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles.split()[0])
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol else None


def compare(corpus: Path, before: Path, after: Path) -> dict:
    source_smiles = list(
        dict.fromkeys(
            row["smiles"]
            for line in corpus.read_text().splitlines()
            if (row := json.loads(line)).get("smiles")
        )
    )
    rows = []
    with tempfile.TemporaryDirectory() as directory:
        mol_path = Path(directory) / "case.mol"
        for smiles in source_smiles:
            source = identity(smiles)
            if source is None:
                continue
            mol_path.unlink(missing_ok=True)
            if run("obabel", "-:" + smiles, "-omol", "-O", str(mol_path)) is None:
                continue
            if not mol_path.is_file():
                continue
            lines = mol_path.read_text().splitlines()
            atom_count = int(lines[3][:3])
            parity_count = sum(
                line[39:42].strip() in ("1", "2")
                for line in lines[4 : 4 + atom_count]
            )
            if not parity_count:
                continue
            babel = identity(run("obabel", str(mol_path), "-osmi"))

            def converted(binary: Path) -> str | None:
                return identity(
                    run(
                        str(binary),
                        "convert",
                        "--input",
                        str(mol_path),
                        "--input-format",
                        "mol",
                        "--output-format",
                        "smiles",
                    )
                )

            rows.append(
                {
                    "source_smiles": smiles,
                    "input_sha256": sha256(mol_path),
                    "source_identity": source,
                    "openbabel_identity": babel,
                    "before_identity": converted(before),
                    "after_identity": converted(after),
                    "atom_parity_centers": parity_count,
                }
            )
    comparable = [row for row in rows if row["openbabel_identity"] == row["source_identity"]]
    summary = {
        "source_unique": len(source_smiles),
        "with_atom_parity": len(rows),
        "comparable": len(comparable),
        "before_source_identity_matches": sum(
            row["before_identity"] == row["source_identity"] for row in comparable
        ),
        "after_source_identity_matches": sum(
            row["after_identity"] == row["source_identity"] for row in comparable
        ),
        "excluded_by_source_gate": len(rows) - len(comparable),
    }
    return {
        "tools": {
            "rdkit": rdBase.rdkitVersion,
            "openbabel": run("obabel", "-V"),
            "before": str(before),
            "after": str(after),
            "corpus_sha256": sha256(corpus),
            "before_sha256": sha256(before),
            "after_sha256": sha256(after),
        },
        "summary": summary,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=Path("validation/cip_label_corpus.jsonl"))
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = compare(args.corpus, args.before, args.after)
    header = json.dumps(
        {"tools": report["tools"], "summary": report["summary"]},
        indent=2,
        ensure_ascii=False,
    )
    rows = ",\n".join(
        "    " + json.dumps(row, ensure_ascii=False, separators=(", ", ": "))
        for row in report["rows"]
    )
    args.report.write_text(header[:-2] + ',\n  "rows": [\n' + rows + "\n  ]\n}\n")
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
