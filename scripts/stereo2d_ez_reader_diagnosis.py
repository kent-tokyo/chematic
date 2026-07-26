#!/usr/bin/env python3
"""Differential validation for P1-S2 (E/Z direction perception wired into the
MOL V2000/V3000 readers via `chematic_perception::apply_ez_directions_from_2d_ex`).

Companion to `scripts/stereo2d_reader_diagnosis.py` (the tetrahedral
reader-integration diagnosis) -- same structure, this time for double-bond
E/Z direction.

Reads validation/results/stereo2d_ez_reader_fixture_dump.jsonl (produced by
`cargo run -p chematic-mol --example stereo2d_ez_reader_fixture_dump`) and,
for each fixture, records:

1. RDKit's own E/Z verdict on the ORIGINAL MOL block (`Chem.MolFromMolBlock`).
2. chematic's `write()` output.
3. chematic's `canonical_smiles()` output.
4. RDKit's E/Z verdict after re-parsing (2) and (3).
5. Standard InChI /b layer for all three parses.

The comparison uses the InChI /b layer (atom-order-independent) as the
primary "same semantic E/Z" oracle, not exact SMILES text or a specific
physical bond -- per the task spec, only the semantics need to match.

Acceptance: verdict(original) == verdict(after chematic write) ==
verdict(after chematic canonical_smiles) for every fixture.

As a negative control, one fixture's expected verdict is deliberately
inverted before comparison to confirm this harness actually detects a
mismatch (exits 1) -- run with --negative-control to exercise this path
without corrupting the real summary file.

Fail-closed: exits 1 on an unexpected fixture-ID set, a duplicate ID, a
parse failure, or an InChI/b mismatch across the three parses.

Run:
    /path/to/venv/bin/python scripts/stereo2d_ez_reader_diagnosis.py
"""

import argparse
import io
import json
import sys
from contextlib import redirect_stderr
from pathlib import Path

from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parent.parent
DUMP_PATH = ROOT / "validation" / "results" / "stereo2d_ez_reader_fixture_dump.jsonl"
SUMMARY_PATH = ROOT / "validation" / "results" / "stereo2d_ez_reader_diagnosis_summary.json"

RDKIT_PINNED_COMMIT = "8afba32ec539dcb2369bc84549d802aca3f7eb39"
EXPECTED_RDKIT_VERSION = "2026.03.3"

EXPECTED_FIXTURE_IDS = {
    "z_2butene_v2000",
    "e_2butene_v2000",
    "z_2butene_v3000",
    "e_2butene_v3000",
    "trisubstituted_alkene_v2000",
    "tetrasubstituted_alkene_v2000",
    "conjugated_diene_v2000",
    "exocyclic_double_bond_v2000",
    "isotopic_substituent_v3000",
    "wedge_and_ez_coexist_v2000",
    "wedge_adjacent_to_double_bond_v2000",
    "atom_renumbered_v2000",
    "rotated_v2000",
    "mirrored_v2000",
}


def inchi_b_layer(mol):
    """Return the `/b...` layer of `mol`'s standard InChI, or None if the
    InChI has no b-layer (RDKit's own verdict: no resolvable E/Z stereo)."""
    if mol is None:
        return None
    stderr_buf = io.StringIO()
    with redirect_stderr(stderr_buf):
        inchi = Chem.MolToInchi(mol)
    if not inchi:
        return None
    for part in inchi.split("/"):
        if part.startswith("b"):
            return part
    return None


def parse_and_verdict_molblock(mol_block):
    stderr_buf = io.StringIO()
    with redirect_stderr(stderr_buf):
        mol = Chem.MolFromMolBlock(mol_block)
    return mol, inchi_b_layer(mol), stderr_buf.getvalue()


def parse_and_verdict_smiles(smiles):
    if not smiles:
        return None, None, ""
    stderr_buf = io.StringIO()
    with redirect_stderr(stderr_buf):
        mol = Chem.MolFromSmiles(smiles)
    return mol, inchi_b_layer(mol), stderr_buf.getvalue()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--negative-control",
        action="store_true",
        help="Invert the first fixture's expected verdict to confirm the harness detects a mismatch, then exit without writing the summary file.",
    )
    args = parser.parse_args()

    installed_version = rdBase.rdkitVersion
    version_mismatch = installed_version != EXPECTED_RDKIT_VERSION
    if version_mismatch:
        print(
            f"WARNING: installed rdkit=={installed_version} != pinned reference "
            f"{EXPECTED_RDKIT_VERSION}",
            file=sys.stderr,
        )

    if not DUMP_PATH.exists():
        print(f"FAIL: {DUMP_PATH} does not exist -- run the fixture dump example first", file=sys.stderr)
        sys.exit(1)

    rows = [json.loads(line) for line in DUMP_PATH.read_text().splitlines() if line.strip()]
    ids = [r["id"] for r in rows]
    if len(ids) != len(set(ids)):
        print("FAIL: duplicate fixture IDs in dump", file=sys.stderr)
        sys.exit(1)
    id_set = set(ids)
    if id_set != EXPECTED_FIXTURE_IDS:
        missing = EXPECTED_FIXTURE_IDS - id_set
        extra = id_set - EXPECTED_FIXTURE_IDS
        print(f"FAIL: fixture ID set mismatch. missing={missing} extra={extra}", file=sys.stderr)
        sys.exit(1)

    results = []
    unexplained = []
    for row in rows:
        if "parse_error" in row:
            print(f"FAIL: fixture {row['id']} failed to parse: {row['parse_error']}", file=sys.stderr)
            sys.exit(1)

        _, verdict_original, stderr_original = parse_and_verdict_molblock(row["mol_block"])
        _, verdict_write, stderr_write = parse_and_verdict_smiles(row["write"])
        _, verdict_canon, stderr_canon = parse_and_verdict_smiles(row["canonical_smiles"])

        if args.negative_control and row is rows[0]:
            # Deliberately corrupt the recorded "original" verdict so the
            # comparison below must fail -- proves this harness is not
            # vacuously green.
            verdict_original = (verdict_original or "") + "_CORRUPTED_FOR_NEGATIVE_CONTROL"

        agree = verdict_original == verdict_write == verdict_canon
        entry = {
            "id": row["id"],
            "description": row["description"],
            "mol_block_verdict": verdict_original,
            "write_output": row["write"],
            "write_verdict": verdict_write,
            "canonical_smiles": row["canonical_smiles"],
            "canonical_verdict": verdict_canon,
            "agree": agree,
            "ez_diagnostics": row.get("ez_diagnostics", []),
            "rdkit_stderr": (stderr_original + stderr_write + stderr_canon).strip(),
        }
        results.append(entry)
        if not agree:
            unexplained.append(entry)

    if args.negative_control:
        if unexplained:
            print(
                f"OK (negative control): harness correctly detected {len(unexplained)} "
                "deliberately-corrupted mismatch(es); exiting 1 as expected."
            )
            sys.exit(1)
        print(
            "FAIL (negative control): harness did NOT detect the deliberately-corrupted "
            "mismatch -- the comparison logic itself is broken.",
            file=sys.stderr,
        )
        sys.exit(1)

    summary = {
        "rdkit_version": installed_version,
        "rdkit_version_mismatch": version_mismatch,
        "rdkit_pinned_commit": RDKIT_PINNED_COMMIT,
        "total_fixtures": len(rows),
        "unexplained_count": len(unexplained),
        "unexplained": unexplained,
        "results": results,
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2))

    agree_count = sum(1 for r in results if r["agree"])
    print(f"E/Z semantic agreement (InChI /b layer): {agree_count}/{len(results)}")

    if unexplained:
        print(f"FAIL: {len(unexplained)} fixture(s) disagree across original/write/canonical:", file=sys.stderr)
        for u in unexplained:
            print(f"  {u['id']}: original={u['mol_block_verdict']} write={u['write_verdict']} canonical={u['canonical_verdict']}", file=sys.stderr)
        sys.exit(1)

    print(f"OK: {len(results)}/{len(results)} fixtures agree, 0 unexplained.")
    sys.exit(0)


if __name__ == "__main__":
    main()
