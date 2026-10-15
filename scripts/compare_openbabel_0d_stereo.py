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
    try:
        process = subprocess.run(args, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return process.stdout.strip() if process.returncode == 0 else None


def identity(smiles: str | None) -> str | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles.split()[0])
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol else None


def compare(corpus: Path, before: Path, after: Path, stereo_only: bool = False) -> dict:
    if corpus.suffix == ".smi":
        candidates = [
            line.split()[0]
            for line in corpus.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]
    else:
        candidates = [
            row["smiles"]
            for line in corpus.read_text().splitlines()
            if (row := json.loads(line)).get("smiles")
        ]
    source_smiles = list(
        dict.fromkeys(
            smiles
            for smiles in candidates
            if not stereo_only or any(token in smiles for token in ("@", "/", "\\"))
        )
    )
    rows = []
    invalid_source = 0
    invalid_mol_output = 0
    non_v2000_output = 0
    without_stereo_markers = 0
    skipped = []
    with tempfile.TemporaryDirectory() as directory:
        mol_path = Path(directory) / "case.mol"
        for smiles in source_smiles:
            source = identity(smiles)
            if source is None:
                invalid_source += 1
                skipped.append({"source_smiles": smiles, "reason": "invalid_source"})
                continue
            lines = []
            for _ in range(3):
                mol_path.unlink(missing_ok=True)
                if run("obabel", "-:" + smiles, "-omol", "-O", str(mol_path)) is None:
                    continue
                if mol_path.is_file():
                    lines = mol_path.read_text().splitlines()
                if len(lines) >= 4:
                    break
            if len(lines) < 4:
                invalid_mol_output += 1
                skipped.append(
                    {"source_smiles": smiles, "reason": "invalid_mol_output"}
                )
                continue
            if not lines[3].endswith("V2000"):
                non_v2000_output += 1
                skipped.append({"source_smiles": smiles, "reason": "non_v2000_output"})
                continue
            try:
                atom_count = int(lines[3][:3])
                bond_count = int(lines[3][3:6])
            except ValueError:
                invalid_mol_output += 1
                skipped.append(
                    {"source_smiles": smiles, "reason": "invalid_mol_counts"}
                )
                continue
            if len(lines) < 4 + atom_count + bond_count:
                invalid_mol_output += 1
                skipped.append(
                    {"source_smiles": smiles, "reason": "truncated_mol_output"}
                )
                continue
            parity_count = sum(
                line[39:42].strip() in ("1", "2") for line in lines[4 : 4 + atom_count]
            )
            directional_count = sum(
                line[9:12].strip() in ("1", "6")
                for line in lines[4 + atom_count : 4 + atom_count + bond_count]
            )
            if not parity_count and not directional_count:
                without_stereo_markers += 1
                skipped.append({"source_smiles": smiles, "reason": "no_stereo_markers"})
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
                    "directional_bonds": directional_count,
                }
            )
    assert len(source_smiles) == (
        invalid_source
        + invalid_mol_output
        + non_v2000_output
        + without_stereo_markers
        + len(rows)
    )
    comparable = [
        row for row in rows if row["openbabel_identity"] == row["source_identity"]
    ]
    summary = {
        "source_unique": len(source_smiles),
        "invalid_source": invalid_source,
        "invalid_mol_output": invalid_mol_output,
        "non_v2000_output": non_v2000_output,
        "without_stereo_markers": without_stereo_markers,
        "with_atom_parity": sum(row["atom_parity_centers"] > 0 for row in rows),
        "with_directional_bonds": sum(row["directional_bonds"] > 0 for row in rows),
        "with_stereo_markers": len(rows),
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
            "stereo_only": stereo_only,
            "mol_generation_attempts": 3,
            "before_sha256": sha256(before),
            "after_sha256": sha256(after),
        },
        "summary": summary,
        "rows": rows,
        "skipped": skipped,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--corpus", type=Path, default=Path("validation/cip_label_corpus.jsonl")
    )
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--stereo-only", action="store_true")
    args = parser.parse_args()
    report = compare(args.corpus, args.before, args.after, args.stereo_only)
    header = json.dumps(
        {"tools": report["tools"], "summary": report["summary"]},
        indent=2,
        ensure_ascii=False,
    )
    rows = ",\n".join(
        "    " + json.dumps(row, ensure_ascii=False, separators=(", ", ": "))
        for row in report["rows"]
    )
    skipped = ",\n".join(
        "    " + json.dumps(row, ensure_ascii=False, separators=(", ", ": "))
        for row in report["skipped"]
    )
    args.report.write_text(
        header[:-2]
        + ',\n  "rows": [\n'
        + rows
        + '\n  ],\n  "skipped": [\n'
        + skipped
        + "\n  ]\n}\n"
    )
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
