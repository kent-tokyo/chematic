#!/usr/bin/env python3
"""Paired 3D MOL2 tetrahedral identity check against Open Babel 3.2.1.

RDKit embeds source SMILES, Open Babel writes MOL2, and RDKit rereads that
exact MOL2 before either converter is scored. Source-gate failures remain in
the report. This is a bounded generated panel, not independent chemical truth.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

SEEDS = [
    "C[C@H](F)Cl",
    "C[C@H](F)Br",
    "C[C@H](Cl)Br",
    "CC[C@H](F)Cl",
    "CC[C@H](Cl)Br",
    "N[C@H](C)C(=O)O",
    "O=C(O)[C@H](O)C",
    "C[C@H](O)N",
    "C[C@H](F)O",
    "C[C@H]1CCCO1",
    "F[C@](Cl)(Br)I",
    "C[C@H](F)[C@H](Cl)Br",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(mol: Chem.Mol | None) -> str | None:
    if mol is None:
        return None
    return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)


def output_identity(smiles: str) -> str | None:
    return identity(Chem.MolFromSmiles(smiles)) if smiles else None


def run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, capture_output=True, text=True, check=False, timeout=30)


def converted(tool: Path, mol2: Path, *, openbabel: bool) -> dict:
    if openbabel:
        result = run(str(tool), "-imol2", str(mol2), "-ocan", "-xn")
        smiles = result.stdout.splitlines()[0].split("\t", 1)[0].strip() if result.stdout else ""
    else:
        result = run(str(tool), "convert", "--input-format", "mol2",
                     "--output-format", "smiles", "--input", str(mol2))
        smiles = result.stdout.strip()
    return {"exit_code": result.returncode, "smiles": smiles,
            "identity": output_identity(smiles) if result.returncode == 0 else None,
            "stderr": result.stderr.strip()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chematic", type=Path, required=True)
    parser.add_argument("--before", type=Path)
    parser.add_argument("--obabel", type=Path, default=Path("/opt/homebrew/bin/obabel"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    tools = {"chematic": args.chematic.resolve(), "openbabel": args.obabel.resolve()}
    if args.before:
        tools["before"] = args.before.resolve()
    for executable in tools.values():
        if not executable.is_file():
            parser.error(f"missing executable: {executable}")

    rows = []
    with tempfile.TemporaryDirectory(prefix="chematic-mol2-3d-") as directory:
        temporary = Path(directory)
        for index, seed in enumerate(SEEDS):
            for invert in (False, True):
                source = seed.replace("@", "@@", 1) if invert else seed
                row = {"source": source, "seed_index": index, "inverted": invert}
                original = Chem.MolFromSmiles(source)
                row["expected"] = identity(original)
                embedded = Chem.AddHs(original)
                if AllChem.EmbedMolecule(embedded, randomSeed=20261016 + 2 * index + invert) != 0:
                    row["gate"] = "embedding_failed"
                    rows.append(row)
                    continue
                mol = temporary / "source.mol"
                mol2 = temporary / "source.mol2"
                mol.write_text(Chem.MolToMolBlock(embedded))
                mol2.unlink(missing_ok=True)
                written = run(str(tools["openbabel"]), "-imol", str(mol), "-omol2", "-O", str(mol2))
                if written.returncode != 0 or not mol2.is_file():
                    row["gate"] = "openbabel_writer_failed"
                    row["writer_stderr"] = written.stderr.strip()
                    rows.append(row)
                    continue
                row["mol2_sha256"] = sha256(mol2)
                reread = Chem.MolFromMol2File(str(mol2), sanitize=True, removeHs=False)
                row["source_reread"] = identity(reread)
                if row["source_reread"] != row["expected"]:
                    row["gate"] = "source_reread_mismatch"
                    rows.append(row)
                    continue
                row["gate"] = "comparable"
                for name, executable in tools.items():
                    row[name] = converted(executable, mol2, openbabel=name == "openbabel")
                    row[name]["matches_source"] = row[name]["identity"] == row["expected"]
                rows.append(row)

    comparable = [row for row in rows if row["gate"] == "comparable"]
    outcomes = Counter(
        (row["chematic"]["matches_source"], row["openbabel"]["matches_source"])
        for row in comparable
    )
    summary = {"attempted": len(rows), "comparable": len(comparable),
               "source_exclusions": len(rows) - len(comparable),
               "chematic_correct": sum(row["chematic"]["matches_source"] for row in comparable),
               "openbabel_correct": sum(row["openbabel"]["matches_source"] for row in comparable),
               "chematic_only": outcomes[(True, False)],
               "openbabel_only": outcomes[(False, True)],
               "both_wrong": outcomes[(False, False)]}
    if args.before:
        summary["before_correct"] = sum(row["before"]["matches_source"] for row in comparable)
    root = Path(__file__).resolve().parents[1]
    source_files = (
        "crates/chematic-mol/src/mol2_tripos.rs",
        "crates/chematic-mol/src/mol2_source.rs",
        "crates/chematic-mol/src/mol2000.rs",
    )
    report = {"schema": "openbabel_mol2_3d_tetrahedral_paired_v1",
              "scope": "generated RDKit 3D chiral molecules with explicit H, Open Babel MOL2 writer, RDKit source/re-read gate",
              "rdkit_version": rdBase.rdkitVersion,
              "script_sha256": sha256(Path(__file__)),
              "source_sha256": {name: sha256(root / name) for name in source_files},
              "tools": {name: {"path": str(path), "sha256": sha256(path)} for name, path in tools.items()},
              "summary": summary, "rows": rows,
              "not_claimed": ["all MOL2 dialects", "independent experimental stereochemical truth",
                              "toolkit-wide superiority"]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
