#!/usr/bin/env python3
"""Compare CheMatic and Open Babel stereo preservation on identical MOL inputs.

RDKit writes each 2D MOL block and reopens it before either converter runs.
RDKit then reads both output SMILES and checks semantic identity. Rows whose
source MOL already changed identity are retained as exclusions, never wins.
This measures the two CLI MOL-to-SMILES workflows on one declared corpus.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

from rdkit import Chem, rdBase

from stereo_format_roundtrip_gate import (
    CORPUS,
    DEFAULT_CLI,
    FORMATS,
    ROOT,
    classify_stereo_difference,
    identity,
    load_cases,
    mol_block,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def openbabel_convert(executable: Path, path: Path) -> tuple[int, str, str]:
    result = subprocess.run(
        [str(executable), "-imol", str(path), "-osmi"],
        text=True,
        capture_output=True,
        timeout=30,
        check=False,
    )
    first_line = result.stdout.splitlines()[0] if result.stdout else ""
    smiles = first_line.split("\t", 1)[0].strip()
    return result.returncode, smiles, result.stderr.strip()


def chematic_convert(executable: Path, input_format: str, path: Path) -> tuple[int, str, str]:
    result = subprocess.run(
        [
            str(executable), "convert", "--input-format", input_format,
            "--output-format", "smiles", "--input", str(path),
        ],
        text=True,
        capture_output=True,
        timeout=30,
        check=False,
    )
    return result.returncode, result.stdout.strip(), result.stderr.strip()


def tool_version(executable: Path, flag: str) -> str:
    result = subprocess.run(
        [str(executable), flag], text=True, capture_output=True, check=True
    )
    return (result.stdout or result.stderr).strip()


def evaluate(expected: Chem.Mol, expected_identity: str, exit_code: int, smiles: str) -> dict:
    parsed = Chem.MolFromSmiles(smiles) if exit_code == 0 and smiles else None
    actual_identity = identity(parsed)
    equal = actual_identity == expected_identity
    classification = None
    if not equal:
        classification, _, _ = classify_stereo_difference(expected, parsed)
    return {
        "exit_code": exit_code,
        "smiles": smiles,
        "identity": actual_identity,
        "semantic_equal": equal,
        "difference_class": classification,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--obabel", type=Path, default=Path(shutil.which("obabel") or "obabel"))
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    obabel = args.obabel.resolve()
    corpus = args.corpus.resolve()
    if not cli.is_file() or not obabel.is_file():
        parser.error("both --cli and --obabel must name executable files")

    cases = load_cases(corpus)
    rows = []
    with tempfile.TemporaryDirectory(prefix="chematic-openbabel-stereo-") as directory:
        temporary = Path(directory)
        for case_index, case in enumerate(cases):
            expected = Chem.MolFromSmiles(case["smiles"])
            assert expected is not None
            for profile, cli_format in FORMATS:
                source = mol_block(case["smiles"], v3000=profile == "mol_v3000")
                source_hash = hashlib.sha256(source.encode()).hexdigest()
                reopened = Chem.MolFromMolBlock(source, sanitize=True, removeHs=False)
                row = {
                    "case_line": case["line"],
                    "bucket": case["bucket"],
                    "profile": profile,
                    "source_smiles": case["smiles"],
                    "expected_identity": case["identity"],
                    "source_sha256": source_hash,
                }
                if identity(reopened) != case["identity"]:
                    row["source_contract"] = "excluded"
                    row["source_error"] = "RDKit writer/reopen identity mismatch"
                    rows.append(row)
                    continue

                path = temporary / f"{case_index}-{profile}.mol"
                path.write_text(source, encoding="utf-8")
                chematic_exit, chematic_smiles, chematic_error = chematic_convert(
                    cli, cli_format, path
                )
                obabel_exit, obabel_smiles, obabel_error = openbabel_convert(
                    obabel, path
                )
                row["source_contract"] = "comparable"
                row["chematic"] = evaluate(
                    expected, case["identity"], chematic_exit, chematic_smiles
                )
                row["openbabel"] = evaluate(
                    expected, case["identity"], obabel_exit, obabel_smiles
                )
                row["chematic"]["stderr"] = chematic_error
                row["openbabel"]["stderr"] = obabel_error
                rows.append(row)

    profiles = {}
    for profile, _ in FORMATS:
        profile_rows = [row for row in rows if row["profile"] == profile]
        comparable = [row for row in profile_rows if row["source_contract"] == "comparable"]
        outcomes = Counter(
            (
                "both_correct"
                if row["chematic"]["semantic_equal"] and row["openbabel"]["semantic_equal"]
                else "chematic_only_correct"
                if row["chematic"]["semantic_equal"]
                else "openbabel_only_correct"
                if row["openbabel"]["semantic_equal"]
                else "both_incorrect"
            )
            for row in comparable
        )
        profiles[profile] = {
            "source_exclusions": len(profile_rows) - len(comparable),
            "comparable": len(comparable),
            "distinct_source_identities": len({row["expected_identity"] for row in comparable}),
            "both_correct": outcomes["both_correct"],
            "chematic_only_correct": outcomes["chematic_only_correct"],
            "distinct_chematic_only_correct_identities": len({
                row["expected_identity"]
                for row in comparable
                if row["chematic"]["semantic_equal"]
                and not row["openbabel"]["semantic_equal"]
            }),
            "openbabel_only_correct": outcomes["openbabel_only_correct"],
            "both_incorrect": outcomes["both_incorrect"],
        }

    result = {
        "schema": "openbabel_stereo_format_paired_v1",
        "scope": "One-host RDKit-written 2D MOL V2000/V3000 to CLI SMILES, adjudicated by RDKit identity",
        "host": {"system": platform.system(), "machine": platform.machine()},
        "script_sha256": sha256(Path(__file__).resolve()),
        "source_revision": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
        ).stdout.strip(),
        "rdkit_version": rdBase.rdkitVersion,
        "chematic_cli": {"path": str(cli), "sha256": sha256(cli), "version": tool_version(cli, "--version")},
        "openbabel_cli": {"path": str(obabel), "sha256": sha256(obabel), "version": tool_version(obabel, "-V")},
        "corpus": {
            "path": str(corpus.relative_to(ROOT)),
            "sha256": sha256(corpus),
            "cases": len(cases),
        },
        "profiles": profiles,
        "rows": rows,
        "not_claimed": [
            "toolkit-wide superiority",
            "independent stereochemical truth beyond the RDKit source/reopen contract",
            "published-package behavior",
            "cross-platform reproducibility",
            "speed or memory advantage",
        ],
    }
    output = args.output if args.output.is_absolute() else ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(profiles, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
