#!/usr/bin/env python3
"""Check an Indigo-created and edited V3000 SGROUP record through chematic.

The gate uses Indigo's public SGroup API to add an empty SUP group, set a
one-atom membership, and replace it with a two-atom membership before writing
V3000.  chematic therefore receives an actual Indigo API-created and edited
record, not a hand-written seed merely reserialized by Indigo.  This remains a
syntax/preservation gate, not typed SGROUP semantic equivalence: it proves
that this bounded SUP line remains present with a synchronized CTAB SGROUP
count and can be reopened by pinned RDKit and Indigo bindings.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

from indigo import Indigo
from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = ROOT / "target" / "debug" / "chematic"
INITIAL_ATOM_INDICES = [0]
FINAL_ATOM_INDICES = [0, 1]
SGROUP_NAME = "core"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        parser.error(f"chematic CLI not found: {cli}")

    indigo = Indigo()
    indigo.setOption("molfile-saving-mode", "3000")
    try:
        indigo_source = indigo.loadMolecule("CC")
        sgroup = indigo_source.addSGroup("SUP", 0)
        if sgroup.setSGroupAtoms(INITIAL_ATOM_INDICES) != 1:
            raise RuntimeError("Indigo did not accept the initial SUP atom membership")
        if sgroup.setSGroupAtoms(FINAL_ATOM_INDICES) != 1:
            raise RuntimeError("Indigo did not accept the edited SUP atom membership")
        if sgroup.setSGroupName(SGROUP_NAME) != 1:
            raise RuntimeError("Indigo did not accept the SUP label")
        indigo_written_source = indigo_source.molfile()
    except Exception as exc:  # stable public exception class is not exposed
        parser.error(f"Indigo could not create and write the SUP group: {exc}")
    source_sgroup_lines = [
        line for line in indigo_written_source.splitlines() if " SUP " in line
    ]
    if not source_sgroup_lines:
        parser.error("Indigo API-created SUP source contains no SUP SGROUP line")
    expected_membership = "ATOMS=(2 1 2)"
    expected_label = f"LABEL={SGROUP_NAME}"
    source_edit_reflected = all(
        expected_membership in line and expected_label in line
        for line in source_sgroup_lines
    )

    with tempfile.TemporaryDirectory(prefix="chematic-v3000-sgroup-") as raw_tmp:
        tmp = Path(raw_tmp)
        source = tmp / "source.v3000"
        output = tmp / "chematic.v3000"
        source.write_text(indigo_written_source, encoding="utf-8")
        run = subprocess.run(
            [
                str(cli), "convert", "--input-format", "mol_v3000", "--output-format",
                "mol_v3000", "--input", str(source), "--output", str(output),
            ],
            text=True,
            capture_output=True,
        )
        written = output.read_text(encoding="utf-8") if output.exists() else ""

    rdkit_mol = Chem.MolFromMolBlock(written, sanitize=True, removeHs=False) if not run.returncode else None
    try:
        indigo.loadMolecule(written) if not run.returncode else None
        indigo_accepted = not run.returncode
    except Exception as exc:  # stable public exception class is not exposed
        indigo_accepted = False
        indigo_error = str(exc)
    else:
        indigo_error = None

    retained = all(line in written for line in source_sgroup_lines)
    count_synchronized = "M  V30 COUNTS 2 1 1 0 0" in written
    result = {
        "schema_version": 1,
        "profile": "v3000_sgroup_external_reader_v1",
        "indigo_version": indigo.version(),
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "cli": str(cli),
        "cli_exit_code": run.returncode,
        "indigo_written_source_sha256": hashlib.sha256(
            indigo_written_source.encode("utf-8")
        ).hexdigest(),
        "indigo_written_sgroup_lines": source_sgroup_lines,
        "source_is_indigo_api_created": True,
        "source_sgroup_api": {
            "type": "SUP",
            "initial_atom_indices": INITIAL_ATOM_INDICES,
            "final_atom_indices": FINAL_ATOM_INDICES,
            "name": SGROUP_NAME,
        },
        "source_edit_reflected": source_edit_reflected,
        "sgroup_line_retained": retained,
        "sgroup_count_synchronized": count_synchronized,
        "rdkit_accepted": rdkit_mol is not None,
        "indigo_accepted": indigo_accepted,
        "indigo_error": indigo_error,
        "gate_passed": (
            not run.returncode
            and source_edit_reflected
            and retained
            and count_synchronized
            and rdkit_mol is not None
            and indigo_accepted
        ),
        "scope": (
            "one Indigo API-created and edited V3000 SUP SGROUP is preserved "
            "and accepted by external readers"
        ),
        "not_claimed": [
            "typed SGROUP semantic equivalence",
            "SGROUP editing or atom-deletion reference updates",
            "COLLECTION enhanced-stereo interoperability",
            "ENDPTS/ATTACH or haptic chemistry semantics",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("gate_passed", "rdkit_accepted", "indigo_accepted")}))
    return 0 if result["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
