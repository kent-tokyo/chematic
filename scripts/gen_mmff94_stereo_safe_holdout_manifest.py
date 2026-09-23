#!/usr/bin/env python3
"""Freeze an outcome-blind MMFF94 stereo-safe post-candidate holdout.

The candidate commit is supplied before selection. Molecules are selected only
from source and RDKit parse metadata; CheMatic is never executed here. The
lowest SHA-256 ranks provide a deterministic random-like sample independent of
source-file order. Half the rows have declared atom or bond stereochemistry so
the stereo-safe contract is exercised directly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).parent.parent
DEFAULT_SOURCE = ROOT / "scripts" / "descriptor_census_corpus.smi"
DEFAULT_OUTPUT = (
    ROOT
    / "validation"
    / "manifests"
    / "mmff94_stereo_safe_post_freeze_holdout_v1.json"
)
BASE_EXCLUDED_MANIFESTS = (
    ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
    ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
)
ALLOWED_ELEMENTS = {"C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
MIN_HEAVY_ATOMS = 8
MAX_HEAVY_ATOMS = 60
STEREO_COUNT = 50
GENERAL_COUNT = 50
DEFAULT_HASH_SALT = "chematic-mmff94-stereo-safe-holdout-v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate-commit",
        required=True,
        help="full frozen source-candidate commit evaluated after manifest creation",
    )
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--salt", default=DEFAULT_HASH_SALT)
    parser.add_argument(
        "--exclude-manifest",
        type=Path,
        action="append",
        default=[],
        help="additional manifest whose canonical molecules must be excluded",
    )
    return parser.parse_args()


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    args = parse_args()
    from rdkit import Chem, __version__ as rdkit_version

    source = args.source if args.source.is_absolute() else ROOT / args.source
    output = args.output if args.output.is_absolute() else ROOT / args.output
    extra_exclusions = [
        path if path.is_absolute() else ROOT / path for path in args.exclude_manifest
    ]
    excluded_manifests = (*BASE_EXCLUDED_MANIFESTS, *extra_exclusions)

    excluded: set[str] = set()
    for path in excluded_manifests:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        for row in manifest["molecules"]:
            molecule = Chem.MolFromSmiles(row["smiles"])
            if molecule is not None:
                excluded.add(Chem.MolToSmiles(molecule, isomericSmiles=True))

    candidates: dict[str, dict[str, object]] = {}
    parse_failures = 0
    filtered = 0
    for line_number, line in enumerate(
        source.read_text(encoding="utf-8").splitlines(), start=1
    ):
        smiles = line.strip()
        if not smiles or "." in smiles:
            filtered += bool(smiles)
            continue
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            parse_failures += 1
            continue
        elements = {atom.GetSymbol() for atom in molecule.GetAtoms()}
        heavy_atoms = molecule.GetNumAtoms()
        if (
            not elements.issubset(ALLOWED_ELEMENTS)
            or not MIN_HEAVY_ATOMS <= heavy_atoms <= MAX_HEAVY_ATOMS
        ):
            filtered += 1
            continue
        canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
        if canonical in excluded or canonical in candidates:
            continue
        has_declared_stereo = any(
            atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            for atom in molecule.GetAtoms()
        ) or any(
            bond.GetStereo() != Chem.BondStereo.STEREONONE
            for bond in molecule.GetBonds()
        )
        candidates[canonical] = {
            "smiles": canonical,
            "heavy_atom_count": heavy_atoms,
            "has_declared_stereo": has_declared_stereo,
            "source_line": line_number,
            "rank": sha256(f"{args.salt}\0{canonical}".encode()),
        }

    stereo = sorted(
        (row for row in candidates.values() if row["has_declared_stereo"]),
        key=lambda row: row["rank"],
    )[:STEREO_COUNT]
    general = sorted(
        (row for row in candidates.values() if not row["has_declared_stereo"]),
        key=lambda row: row["rank"],
    )[:GENERAL_COUNT]
    if len(stereo) != STEREO_COUNT or len(general) != GENERAL_COUNT:
        raise SystemExit(
            f"insufficient candidates: stereo={len(stereo)} general={len(general)}"
        )

    molecules: list[dict[str, object]] = []
    for cohort, rows in (("stereo", stereo), ("general", general)):
        for index, row in enumerate(rows):
            molecules.append(
                {
                    "name": f"mmff94_holdout_{cohort}_{index:03d}",
                    "smiles": row["smiles"],
                    "heavy_atom_count": row["heavy_atom_count"],
                    "primary_category": f"post_freeze_{cohort}",
                    "has_declared_stereo": row["has_declared_stereo"],
                    "source_line": row["source_line"],
                    "selection_sha256": row["rank"],
                }
            )

    manifest = {
        "schema_version": 1,
        "tier": "H",
        "description": (
            "Post-candidate MMFF94 stereo-safe holdout: 50 declared-stereo and "
            "50 general drug-like molecules selected without running CheMatic."
        ),
        "candidate_commit_frozen_before_selection": args.candidate_commit,
        "generator": "scripts/gen_mmff94_stereo_safe_holdout_manifest.py",
        "generator_rdkit_version": rdkit_version,
        "source_file": str(source.relative_to(ROOT)),
        "source_sha256": sha256(source.read_bytes()),
        "excluded_manifests": [str(path.relative_to(ROOT)) for path in excluded_manifests],
        "selection_rule": [
            "single-fragment RDKit-parseable molecules only",
            f"elements restricted to {sorted(ALLOWED_ELEMENTS)}",
            f"heavy-atom count in [{MIN_HEAVY_ATOMS}, {MAX_HEAVY_ATOMS}]",
            "exclude canonical isomeric SMILES already in Tier A or Tier B",
            "deduplicate by RDKit canonical isomeric SMILES",
            f"rank by SHA-256 of salt {args.salt!r} plus canonical SMILES",
            f"take lowest {STEREO_COUNT} declared-stereo and {GENERAL_COUNT} general ranks",
            "never execute or inspect CheMatic outcomes during selection",
        ],
        "parse_failure_count": parse_failures,
        "filtered_count": filtered,
        "eligible_unique_count": len(candidates),
        "molecule_count": len(molecules),
        "corpus_sha256": sha256(
            json.dumps(molecules, sort_keys=True).encode("utf-8")
        ),
        "molecules": molecules,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {output.relative_to(ROOT)} with {len(molecules)} molecules")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
