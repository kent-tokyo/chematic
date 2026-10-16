#!/usr/bin/env python3
"""Compare CheMatic and Open Babel stereo preservation on identical MOL inputs.

RDKit writes each 2D MOL block and reopens it before either converter runs.
RDKit then reads both output SMILES and checks semantic identity. Rows that
cannot be written in the requested MOL version or lose identity before either
converter runs are retained as exclusions, never wins.
This measures the two CLI MOL-to-SMILES workflows on one declared corpus.
The optional bounded stress selection covers additional tetrahedral/E/Z cases.
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


def nonstereogenic_explicit_h_equivalent(expected: Chem.Mol, actual: Chem.Mol) -> bool:
    """Normalize a plain explicit H only when neither molecule specifies stereo."""
    if not any(atom.GetAtomicNum() == 1 for atom in expected.GetAtoms()):
        return False
    for mol in (expected, actual):
        if any(
            atom.GetAtomicNum() == 1 and (atom.GetIsotope() or atom.GetAtomMapNum())
            for atom in mol.GetAtoms()
        ):
            return False
        if any(
            str(info.specified) == "Specified"
            for info in Chem.FindPotentialStereo(Chem.Mol(mol), cleanIt=True, flagPossible=True)
        ):
            return False
    return identity(Chem.RemoveHs(expected)) == identity(Chem.RemoveHs(actual))


def evaluate(expected: Chem.Mol, expected_identity: str, exit_code: int, smiles: str) -> dict:
    parsed = Chem.MolFromSmiles(smiles) if exit_code == 0 and smiles else None
    actual_identity = identity(parsed)
    strict_equal = actual_identity == expected_identity
    hydrogen_equivalent = (
        not strict_equal
        and parsed is not None
        and nonstereogenic_explicit_h_equivalent(expected, parsed)
    )
    equal = strict_equal or hydrogen_equivalent
    classification = None
    if not equal:
        classification, _, _ = classify_stereo_difference(expected, parsed)
    return {
        "exit_code": exit_code,
        "smiles": smiles,
        "identity": actual_identity,
        "strict_equal": strict_equal,
        "hydrogen_representation_equivalent": hydrogen_equivalent,
        "semantic_equal": equal,
        "difference_class": classification,
    }


def load_stress_cases(path: Path, limit: int) -> list[dict]:
    """Select the first unique RDKit identities with tetrahedral or E/Z syntax."""
    seen = set()
    cases = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        smiles = line.split()[0]
        if not any(tag in smiles for tag in ("@", "/", "\\")):
            continue
        if any(tag in smiles for tag in ("@SP", "@TB", "@OH")):
            continue
        mol = Chem.MolFromSmiles(smiles)
        expected_identity = identity(mol)
        if expected_identity is None or expected_identity in seen:
            continue
        seen.add(expected_identity)
        cases.append({
            "line": line_number,
            "bucket": "stress-first-unique",
            "smiles": smiles,
            "identity": expected_identity,
        })
        if len(cases) == limit:
            break
    return cases


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--obabel", type=Path, default=Path(shutil.which("obabel") or "obabel"))
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument(
        "--stress-first-unique",
        type=int,
        help="read the first N unique tetrahedral/E/Z identities from an SMI corpus",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    obabel = args.obabel.resolve()
    corpus = args.corpus.resolve()
    if not cli.is_file() or not obabel.is_file():
        parser.error("both --cli and --obabel must name executable files")

    if args.stress_first_unique is not None and args.stress_first_unique < 1:
        parser.error("--stress-first-unique must be positive")
    cases = (
        load_stress_cases(corpus, args.stress_first_unique)
        if args.stress_first_unique is not None
        else load_cases(corpus)
    )
    rows = []
    with tempfile.TemporaryDirectory(prefix="chematic-openbabel-stereo-") as directory:
        temporary = Path(directory)
        for case_index, case in enumerate(cases):
            expected = Chem.MolFromSmiles(case["smiles"])
            assert expected is not None
            for profile, cli_format in FORMATS:
                row = {
                    "case_line": case["line"],
                    "bucket": case["bucket"],
                    "profile": profile,
                    "source_smiles": case["smiles"],
                    "expected_identity": case["identity"],
                }
                try:
                    source = mol_block(case["smiles"], v3000=profile == "mol_v3000")
                    reopened = Chem.MolFromMolBlock(source, sanitize=True, removeHs=False)
                except Exception as exc:
                    row["source_contract"] = "excluded"
                    row["source_error"] = f"RDKit writer/reopen failed: {exc}"
                    rows.append(row)
                    continue
                row["source_sha256"] = hashlib.sha256(source.encode()).hexdigest()
                version = "V3000" if profile == "mol_v3000" else "V2000"
                if version not in source.splitlines()[3]:
                    row["source_contract"] = "excluded"
                    row["source_error"] = f"RDKit writer emitted a different MOL version than {version}"
                    rows.append(row)
                    continue
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
                "both_match_source"
                if row["chematic"]["semantic_equal"] and row["openbabel"]["semantic_equal"]
                else "chematic_only_matches_source"
                if row["chematic"]["semantic_equal"]
                else "openbabel_only_matches_source"
                if row["openbabel"]["semantic_equal"]
                else "neither_matches_source"
            )
            for row in comparable
        )
        profiles[profile] = {
            "source_exclusions": len(profile_rows) - len(comparable),
            "comparable": len(comparable),
            "distinct_source_identities": len({row["expected_identity"] for row in comparable}),
            "both_match_source": outcomes["both_match_source"],
            "chematic_only_matches_source": outcomes["chematic_only_matches_source"],
            "distinct_chematic_only_matching_identities": len({
                row["expected_identity"]
                for row in comparable
                if row["chematic"]["semantic_equal"]
                and not row["openbabel"]["semantic_equal"]
            }),
            "openbabel_only_matches_source": outcomes["openbabel_only_matches_source"],
            "neither_matches_source": outcomes["neither_matches_source"],
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
