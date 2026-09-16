#!/usr/bin/env python3
"""Check RDKit OR enhanced stereo through chematic V3000 output.

RDKit constructs the source record, including a `MDLV30/STEREL1` COLLECTION
entry for its `STEREO_OR` group. The gate verifies that chematic preserves the group for RDKit
and emits a V3000 record accepted by Indigo. It does not compare rendering,
CIP assignment, or stereo-group editing semantics.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

from indigo import Indigo
from rdkit import Chem

from v3000_rdkit_semantic_gate import cli_provenance


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = ROOT / "target" / "debug" / "chematic"


def relative_signature(mol: Chem.Mol) -> list[tuple[str, tuple[int, ...]]]:
    return sorted(
        (str(group.GetGroupType()), tuple(atom.GetIdx() for atom in group.GetAtoms()))
        for group in mol.GetStereoGroups()
    )


def source_record() -> tuple[str, list[tuple[str, tuple[int, ...]]]]:
    molecule = Chem.RWMol(Chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
    # RDKit's Python STEREO_OR writer emits the observed `MDLV30/STEREL1`
    # token; use its actual API/token pairing rather than inferred names.
    group = Chem.CreateStereoGroup(Chem.StereoGroupType.STEREO_OR, molecule, [1])
    molecule.SetStereoGroups([group])
    source = Chem.MolToMolBlock(molecule, forceV3000=True)
    expected = relative_signature(molecule)
    if "MDLV30/STEREL1" not in source:
        raise RuntimeError("RDKit did not generate the expected relative stereo V3000 fixture")
    return source, expected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        parser.error(f"chematic CLI not found: {cli}")

    source, expected_groups = source_record()
    indigo = Indigo()
    indigo.setOption("molfile-saving-mode", "3000")
    with tempfile.TemporaryDirectory(prefix="chematic-v3000-relative-") as raw_tmp:
        tmp = Path(raw_tmp)
        input_path = tmp / "rdkit-relative.v3000"
        output_path = tmp / "chematic-relative.v3000"
        input_path.write_text(source, encoding="utf-8")
        run = subprocess.run(
            [
                str(cli), "convert", "--input-format", "mol_v3000", "--output-format",
                "mol_v3000", "--input", str(input_path), "--output", str(output_path),
            ],
            text=True,
            capture_output=True,
        )
        written = output_path.read_text(encoding="utf-8") if output_path.exists() else ""

    rdkit_after = Chem.MolFromMolBlock(written, sanitize=True, removeHs=False) if not run.returncode else None
    groups_after = relative_signature(rdkit_after) if rdkit_after is not None else []
    try:
        indigo.loadMolecule(written) if not run.returncode else None
        indigo_accepted = not run.returncode
    except Exception as exc:  # stable public exception class is not exposed
        indigo_accepted = False
        indigo_error = str(exc)
    else:
        indigo_error = None

    result = {
        "schema_version": 1,
        "profile": "v3000_rdkit_or_stereo_v1",
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "indigo_version": indigo.version(),
        **cli_provenance(cli),
        "cli_exit_code": run.returncode,
        "or_collection_retained": "MDLV30/STEREL1" in written,
        "rdkit_expected_groups": expected_groups,
        "rdkit_after_groups": groups_after,
        "rdkit_group_equal": expected_groups == groups_after,
        "indigo_accepted": indigo_accepted,
        "indigo_error": indigo_error,
        "gate_passed": (
            not run.returncode
            and "MDLV30/STEREL1" in written
            and expected_groups == groups_after
            and indigo_accepted
        ),
        "scope": "one RDKit-generated OR enhanced-stereo group survives chematic V3000 round trip",
        "not_claimed": [
            "absolute CIP label equivalence",
            "OR/AND/relative group editing semantics",
            "multiple group or bond-member coverage",
            "Indigo group-type semantic equivalence",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("gate_passed", "rdkit_group_equal", "indigo_accepted")}))
    return 0 if result["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
