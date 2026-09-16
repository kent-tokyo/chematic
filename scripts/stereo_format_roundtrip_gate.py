#!/usr/bin/env python3
"""Check stereo-sensitive RDKit MOL V2000/V3000 inputs through chematic.

Each source SMILES is written by a pinned RDKit binding after 2D depiction.
The exact RDKit MOL source must first reopen with the same RDKit semantic
identity; records that fail that source contract are retained as exclusions.
For comparable records, chematic converts the MOL text to SMILES and RDKit
checks the semantic identity of that result.  This tests reader-to-SMILES
stereo preservation, not byte-identical MOL writing or cross-engine canonical
string equality.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = ROOT / "target" / "debug" / "chematic"
CORPUS = ROOT / "validation" / "cip_label_corpus.jsonl"
EXPECTED_CASES = 155
FORMATS = (("mol_v2000", "mol"), ("mol_v3000", "mol_v3000"))
MAX_TOPOLOGY_MATCHES = 256


def identity(molecule: Chem.Mol | None) -> str | None:
    if molecule is None:
        return None
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def stereo_summary(molecule: Chem.Mol | None) -> dict[str, Any]:
    """Return an RDKit-native, atom-map-independent stereo inventory.

    The reader-to-SMILES comparison deliberately uses RDKit semantic identity
    as its oracle.  A mismatch alone does not say whether chematic discarded
    a declared stereofact, introduced one, or retained the count while
    changing its assignment.  `FindPotentialStereo` gives a stable inventory
    without pretending that the two independently parsed molecules share atom
    indices.  It is a classification aid, not an atom-by-atom adjudication.
    """
    if molecule is None:
        return {
            "specified_total": 0,
            "specified_by_type": {},
            "specified_descriptor_counts": {},
        }
    infos = Chem.FindPotentialStereo(molecule, cleanIt=True, flagPossible=True)
    specified = [info for info in infos if str(info.specified) == "Specified"]
    by_type = Counter(str(info.type) for info in specified)
    descriptors = Counter(f"{info.type}:{info.descriptor}" for info in specified)
    return {
        "specified_total": len(specified),
        "specified_by_type": dict(sorted(by_type.items())),
        "specified_descriptor_counts": dict(sorted(descriptors.items())),
    }


def stereo_records(
    molecule: Chem.Mol | None, output_to_expected: dict[int, int] | None = None
) -> dict[str, str] | None:
    """Describe specified stereo facts in source-atom coordinates.

    `output_to_expected` is a topology-only atom correspondence from a
    separately parsed output molecule back to the original source.  It is
    deliberately independent of atom maps: adding maps to the RDKit CTAB
    changed a few source depictions and invalidated the comparison itself.
    """
    if molecule is None:
        return {}
    records: dict[str, str] = {}
    for info in Chem.FindPotentialStereo(molecule, cleanIt=True, flagPossible=True):
        if str(info.specified) != "Specified":
            continue
        if str(info.type) == "Atom_Tetrahedral":
            center = info.centeredOn
            if output_to_expected is not None:
                center = output_to_expected.get(center, -1)
            if center < 0:
                return None
            key = f"atom:{center}"
        elif str(info.type) == "Bond_Double":
            bond = molecule.GetBondWithIdx(info.centeredOn)
            endpoints = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            if output_to_expected is not None:
                endpoints = tuple(output_to_expected.get(atom, -1) for atom in endpoints)
            if min(endpoints) < 0:
                return None
            first, second = sorted(endpoints)
            key = f"bond:{first}-{second}"
        else:
            return None
        records[key] = str(info.descriptor)
    return dict(sorted(records.items()))


def topology_only(molecule: Chem.Mol) -> Chem.Mol:
    """Copy a graph while dropping all stereo state for isomorphism matching."""
    query = Chem.Mol(molecule)
    Chem.RemoveStereochemistry(query)
    return query


def topology_aligned_stereo(
    expected: Chem.Mol | None, output: Chem.Mol | None
) -> dict[str, Any] | None:
    """Align stereo facts through bounded topology isomorphisms.

    This uses RDKit only to align the same non-stereo molecular graph.  It
    does not decide the desired configuration.  Symmetric graphs can have
    several valid mappings, so the result retains the number searched and
    chooses the mapping with the smallest remaining mismatch.  A capped
    search is disclosed as non-adjudicating rather than treated as exact.
    """
    if expected is None or output is None:
        return None
    expected_records = stereo_records(expected)
    if expected_records is None:
        return None
    matches = output.GetSubstructMatches(
        topology_only(expected),
        uniquify=True,
        useChirality=False,
        maxMatches=MAX_TOPOLOGY_MATCHES,
    )
    if not matches:
        return None

    candidates = []
    for match in matches:
        output_to_expected = {output_index: expected_index for expected_index, output_index in enumerate(match)}
        output_records = stereo_records(output, output_to_expected)
        if output_records is None:
            continue
        missing = sorted(set(expected_records) - set(output_records))
        extra = sorted(set(output_records) - set(expected_records))
        changed = sorted(
            key
            for key in set(expected_records) & set(output_records)
            if expected_records[key] != output_records[key]
        )
        categories = []
        if missing:
            categories.append("stereo_information_loss")
        if changed:
            categories.append("stereo_descriptor_assignment_mismatch")
        if extra:
            categories.append("spurious_stereo_assignment")
        candidates.append((
            (len(missing) + len(changed) + len(extra), len(categories), tuple(categories)),
            output_records,
            missing,
            changed,
            extra,
            categories,
        ))
    if not candidates:
        return None
    _, output_records, missing, changed, extra, categories = min(candidates, key=lambda item: item[0])
    return {
        "expected": expected_records,
        "output": output_records,
        "missing": missing,
        "descriptor_changed": changed,
        "extra": extra,
        "categories": categories,
        "topology_isomorphisms_searched": len(matches),
        "topology_isomorphism_search_capped": len(matches) == MAX_TOPOLOGY_MATCHES,
    }


def classify_stereo_difference(
    expected: Chem.Mol | None, output: Chem.Mol | None
) -> tuple[str, dict[str, Any], dict[str, Any]]:
    """Classify a semantic mismatch without conflating loss and reassignment.

    Exact atom/bond correspondence requires a separate mapped oracle lane.
    This gate instead records conservative inventory-level categories.  They
    make E/Z loss and tetrahedral reassignment independently actionable while
    avoiding a false claim that symmetric atoms have been aligned.
    """
    expected_summary = stereo_summary(expected)
    output_summary = stereo_summary(output)
    if output is None:
        return "rdkit_rejected_chematic_smiles", expected_summary, output_summary

    aligned = topology_aligned_stereo(expected, output)
    if aligned and aligned["categories"]:
        return aligned["categories"][0], expected_summary, output_summary

    expected_total = expected_summary["specified_total"]
    output_total = output_summary["specified_total"]
    if output_total < expected_total:
        return "stereo_information_loss", expected_summary, output_summary
    if output_total > expected_total:
        return "spurious_stereo_assignment", expected_summary, output_summary
    if expected_summary["specified_by_type"] != output_summary["specified_by_type"]:
        return "stereo_type_count_mismatch", expected_summary, output_summary
    if expected_summary["specified_descriptor_counts"] != output_summary["specified_descriptor_counts"]:
        return "stereo_descriptor_assignment_mismatch", expected_summary, output_summary
    return "nonstereo_or_unaligned_stereo_mismatch", expected_summary, output_summary


def load_cases(path: Path) -> list[dict[str, Any]]:
    cases = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("_manifest") is True:
            continue
        smiles = row.get("smiles")
        if not isinstance(smiles, str) or not smiles:
            raise ValueError(f"corpus line {line_number} lacks a SMILES string")
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"RDKit rejected corpus line {line_number}")
        cases.append({
            "line": line_number,
            "bucket": row.get("bucket"),
            "smiles": smiles,
            "identity": identity(molecule),
        })
    if len(cases) != EXPECTED_CASES:
        raise ValueError(f"expected {EXPECTED_CASES} corpus rows, got {len(cases)}")
    return cases


def mol_block(smiles: str, v3000: bool) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    rdDepictor.Compute2DCoords(molecule)
    return Chem.MolToMolBlock(molecule, forceV3000=v3000)


def convert(cli: Path, input_format: str, source_path: Path) -> tuple[int, str, str]:
    completed = subprocess.run(
        [
            str(cli), "convert", "--input-format", input_format,
            "--output-format", "smiles", "--input", str(source_path),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    return completed.returncode, completed.stdout.strip(), completed.stderr.strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        parser.error(f"chematic CLI not found: {cli}")
    corpus_path = args.corpus.resolve()
    cases = load_cases(corpus_path)
    rows: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="chematic-stereo-format-") as directory:
        tmp = Path(directory)
        for case_index, case in enumerate(cases):
            for profile, cli_format in FORMATS:
                row: dict[str, Any] = {
                    "case_line": case["line"],
                    "bucket": case["bucket"],
                    "profile": profile,
                    "source_smiles": case["smiles"],
                }
                try:
                    source = mol_block(case["smiles"], v3000=profile == "mol_v3000")
                except Exception as error:
                    row.update({"source_contract": "excluded", "source_error": str(error)})
                    rows.append(row)
                    continue
                source_reopened = Chem.MolFromMolBlock(source, sanitize=True, removeHs=False)
                if identity(source_reopened) != case["identity"]:
                    row.update({
                        "source_contract": "excluded",
                        "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
                        "source_error": "RDKit writer/reopen semantic identity mismatch",
                    })
                    rows.append(row)
                    continue
                path = tmp / f"{case_index}-{profile}.mol"
                path.write_text(source, encoding="utf-8")
                exit_code, written_smiles, error = convert(cli, cli_format, path)
                row.update({
                    "source_contract": "comparable",
                    "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
                    "cli_exit_code": exit_code,
                })
                output_molecule = Chem.MolFromSmiles(written_smiles) if exit_code == 0 else None
                semantic_equal = identity(output_molecule) == case["identity"]
                expected_molecule = Chem.MolFromSmiles(case["smiles"])
                classification, expected_stereo, output_stereo = classify_stereo_difference(
                    expected_molecule, output_molecule
                )
                topology_alignment = topology_aligned_stereo(expected_molecule, output_molecule)
                row.update({
                    "semantic_equal": semantic_equal,
                    "failure_kind": (
                        None
                        if semantic_equal
                        else "rdkit_rejected_chematic_smiles"
                        if output_molecule is None
                        else "semantic_identity_mismatch"
                    ),
                    "stereo_difference_class": None if semantic_equal else classification,
                    "expected_stereo": expected_stereo if not semantic_equal else None,
                    "output_stereo": output_stereo if not semantic_equal else None,
                    "topology_aligned_stereo": topology_alignment if not semantic_equal else None,
                    "expected_identity": case["identity"] if not semantic_equal else None,
                    "output_smiles": written_smiles if not semantic_equal else None,
                    "output_identity": identity(output_molecule) if not semantic_equal else None,
                    "error": error or None,
                })
                rows.append(row)

    comparable = [row for row in rows if row["source_contract"] == "comparable"]
    failures = [row for row in comparable if not row.get("semantic_equal", False)]
    by_profile = {
        profile: {
            "source_exclusions": sum(1 for row in rows if row["profile"] == profile and row["source_contract"] == "excluded"),
            "comparable": sum(1 for row in comparable if row["profile"] == profile),
            "semantic_equal": sum(1 for row in comparable if row["profile"] == profile and row.get("semantic_equal")),
            "failures": sum(1 for row in failures if row["profile"] == profile),
        }
        for profile, _ in FORMATS
    }
    by_bucket = {
        (bucket if bucket is not None else "unbucketed"): {
            "comparable": sum(1 for row in comparable if row["bucket"] == bucket),
            "semantic_equal": sum(
                1 for row in comparable if row["bucket"] == bucket and row.get("semantic_equal")
            ),
            "failures": sum(1 for row in failures if row["bucket"] == bucket),
        }
        for bucket in sorted({row["bucket"] for row in rows}, key=lambda value: str(value))
    }
    failure_kinds = {
        kind: sum(1 for row in failures if row["failure_kind"] == kind)
        for kind in sorted({row["failure_kind"] for row in failures})
    }
    stereo_difference_classes = {
        kind: sum(1 for row in failures if row.get("stereo_difference_class") == kind)
        for kind in sorted({row.get("stereo_difference_class") for row in failures})
    }
    topology_aligned_categories = {
        category: sum(
            1
            for row in failures
            if category in (row.get("topology_aligned_stereo") or {}).get("categories", [])
        )
        for category in (
            "stereo_information_loss",
            "stereo_descriptor_assignment_mismatch",
            "spurious_stereo_assignment",
        )
    }
    unique_failure_sources = {
        (row["profile"], row["expected_identity"])
        for row in failures
    }
    residual_structures: dict[str, dict[str, Any]] = {}
    for row in failures:
        identity_key = row["expected_identity"]
        entry = residual_structures.setdefault(
            identity_key,
            {
                "source_smiles": row["source_smiles"],
                "source_identity": identity_key,
                "corpus_lines": [],
                "profiles": {},
            },
        )
        entry["corpus_lines"].append(row["case_line"])
        profile_entry = entry["profiles"].setdefault(
            row["profile"],
            {"count": 0, "classes": [], "topology_categories": []},
        )
        profile_entry["count"] += 1
        profile_entry["classes"].append(row["stereo_difference_class"])
        profile_entry["topology_categories"].extend(
            (row.get("topology_aligned_stereo") or {}).get("categories", [])
        )
    for entry in residual_structures.values():
        entry["corpus_lines"] = sorted(set(entry["corpus_lines"]))
        for profile_entry in entry["profiles"].values():
            profile_entry["classes"] = sorted(set(profile_entry["classes"]))
            profile_entry["topology_categories"] = sorted(
                set(profile_entry["topology_categories"])
            )
    result = {
        "schema_version": 2,
        "profile": "stereo_mol_format_roundtrip_v1",
        "rdkit_version": rdBase.rdkitVersion,
        "cli": str(cli),
        "corpus": {
            "path": str(corpus_path.relative_to(ROOT)),
            "sha256": hashlib.sha256(corpus_path.read_bytes()).hexdigest(),
            "cases": len(cases),
        },
        "profiles": by_profile,
        "buckets": by_bucket,
        "failure_count": len(failures),
        "failure_summary": {
            "unique_source_structures": len({identity for _, identity in unique_failure_sources}),
            "unique_source_profiles": len(unique_failure_sources),
            "kinds": failure_kinds,
            "stereo_difference_classes": stereo_difference_classes,
            "topology_aligned_stereo_categories": topology_aligned_categories,
        },
        "failures": failures,
        "residual_structures": list(residual_structures.values()),
        "failure_samples": failures[:20],
        "source_contract_exclusions": [row for row in rows if row["source_contract"] == "excluded"],
        "gate_passed": not failures,
        "scope": "RDKit 2D MOL V2000/V3000 source to chematic reader to SMILES, checked by RDKit semantic identity",
        "not_claimed": [
            "byte-identical MOL output",
            "chematic MOL writer round-trip",
            "cross-engine canonical string equality",
            "300-structure Stereo Torture Suite completion",
            "independent CIP-label correctness",
            "unbounded graph-isomorphism enumeration",
        ],
    }
    output = args.output if args.output.is_absolute() else ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"profiles": by_profile, "failure_count": len(failures), "gate_passed": result["gate_passed"]}, sort_keys=True))
    return 0 if result["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
