#!/usr/bin/env python3
"""P1-A0 diagnostic: cross-reference chematic's current 2D-stereo pipeline
against RDKit on the frozen fixtures emitted by
`cargo run -p chematic-mol --example stereo2d_fixture_dump`.

Diagnostic only. Does not change chematic's production behavior. Reads
validation/results/stereo2d_fixture_dump.jsonl (one row per fixture, produced
by the Rust example) and, for each fixture's raw MOL block, re-parses it with
RDKit (pinned commit 8afba32ec539dcb2369bc84549d802aca3f7eb39 / rdkit==2026.03.3
in this project's venv) to get RDKit's own chiral-tag/CIP/E-Z assignment and
canonical SMILES, then classifies the (chematic, RDKit) pair into a failure
bucket. Writes validation/results/stereo2d_diagnosis_summary.json.

Run:
    .venv/bin/python scripts/stereo2d_diagnosis.py
"""

import json
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

ROOT = Path(__file__).resolve().parent.parent
DUMP_PATH = ROOT / "validation" / "results" / "stereo2d_fixture_dump.jsonl"
SUMMARY_PATH = ROOT / "validation" / "results" / "stereo2d_diagnosis_summary.json"

RDKIT_PINNED_COMMIT = "8afba32ec539dcb2369bc84549d802aca3f7eb39"


def rdkit_read(mol_block):
    """Parse `mol_block` with RDKit (default sanitize=True) and report its
    final stereo assignment.

    Note: RDKit logs parser warnings (e.g. "ambiguous stereochemistry -
    opposing bonds have opposite wedging") straight to process stderr via its
    C++-side log handler, not as Python exceptions or a return value -- they
    are NOT captured here. This script's fixtures were deliberately built to
    avoid triggering that specific warning (see the comment on fixture #2 in
    stereo2d_fixture_dump.rs); if you add new fixtures, watch stderr when
    running this script.
    """
    mol = Chem.MolFromMolBlock(mol_block)
    if mol is None:
        return {"parsed": False}

    rdCIPLabeler.AssignCIPLabels(mol)

    atoms = []
    for a in mol.GetAtoms():
        atoms.append(
            {
                "idx": a.GetIdx(),
                "symbol": a.GetSymbol(),
                "chiral_tag": str(a.GetChiralTag()),
                "cip_code": a.GetPropsAsDict().get("_CIPCode"),
            }
        )
    bonds = []
    for b in mol.GetBonds():
        bonds.append(
            {
                "a1": b.GetBeginAtomIdx(),
                "a2": b.GetEndAtomIdx(),
                "bond_type": str(b.GetBondType()),
                "bond_dir": str(b.GetBondDir()),
                "stereo": str(b.GetStereo()),
            }
        )
    return {
        "parsed": True,
        "atoms": atoms,
        "bonds": bonds,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "any_chiral_tag": any(a["chiral_tag"] != "CHI_UNSPECIFIED" for a in atoms),
        "any_cip_code": any(a["cip_code"] for a in atoms),
        "any_bond_stereo": any(b["stereo"] not in ("STEREONONE",) for b in bonds),
    }


def classify(fixture, rdkit_result):
    """Assign one failure bucket per fixture, based on what chematic
    computed/wrote vs what RDKit's default parse produces."""

    mech = fixture["mechanism"]
    chematic_rs = fixture.get("assign_stereo_from_2d_result", [])
    chematic_ez = fixture.get("assign_ez_from_2d_result", [])
    chirality_reached_writer = fixture.get("chirality_reached_writer", False)
    naive_smiles = fixture.get("naive_smiles_write", "") or ""
    has_naive_direction_token = "/" in naive_smiles or "\\" in naive_smiles

    rdkit_found_stereo = rdkit_result.get("any_chiral_tag") or rdkit_result.get(
        "any_bond_stereo"
    )

    if mech in ("negative_control", "terminal_alkene", "cip_priority_tie"):
        if not chematic_rs and not chematic_ez and not rdkit_found_stereo:
            return "correctly_no_stereo_both_agree"
        return "unexpected_disagreement_on_no_stereo_case"

    if mech == "degenerate_coordinates":
        if not chematic_rs and rdkit_result.get("parsed"):
            return "degenerate_coords_correctly_yields_no_stereo"
        return "unexpected_result_on_degenerate_coords"

    if mech == "coord_atom_count_mismatch":
        return "silent_result_from_corrupted_fallback_positions_not_error"

    if mech == "tetrahedral_3heavy_implicit_h":
        if not chematic_rs and rdkit_result.get("any_chiral_tag"):
            return "rs_not_computed_3heavy_implicit_h_gap"
        return "unexpected_3heavy_result"

    if mech in ("tetrahedral_4heavy", "solid_wedge", "dashed_wedge"):
        if chematic_rs and not chirality_reached_writer:
            if has_naive_direction_token:
                return "rs_computed_but_writer_emits_meaningless_bond_direction_token"
            return "rs_computed_but_not_written_to_chirality"
        if not chematic_rs and rdkit_result.get("any_chiral_tag"):
            return "rs_not_computed_despite_rdkit_success"
        return "unexpected_tetrahedral_result"

    if mech == "wedge_atom_order_reversed":
        # This fixture's own file is non-standard (atom1=substituent for the
        # wedge bond), so "agrees with RDKit's reading of the SAME file" is
        # the right bar -- not "agrees with the standard-order fixture",
        # which is a separate cross-check done in the report by comparing
        # against tetrahedral_4heavy_explicit_h's row directly.
        rdkit_cip = {a["cip_code"] for a in rdkit_result.get("atoms", []) if a.get("cip_code")}
        chematic_cip = {a["cip_code"] for a in chematic_rs}
        if chematic_rs and rdkit_cip and chematic_cip == rdkit_cip:
            return "wedge_atom_order_reversed_agrees_with_rdkit_on_same_file"
        if chematic_rs and rdkit_cip and chematic_cip != rdkit_cip:
            return "wedge_atom_order_reversed_disagrees_with_rdkit_on_same_file"
        if chematic_rs and not rdkit_cip:
            return "wedge_atom_order_reversed_chematic_only"
        return "unexpected_wedge_atom_order_result"

    if mech == "multiple_stereocenters":
        if not chematic_rs and rdkit_result.get("any_chiral_tag"):
            return "rs_not_computed_3heavy_implicit_h_gap"
        return "unexpected_multi_stereocenter_result"

    if mech == "ez_geometry":
        if chematic_ez and not has_naive_direction_token:
            return "ez_computed_but_no_bond_direction_for_writer"
        return "unexpected_ez_result"

    if mech == "contradictory_wedges":
        return "no_consistency_check_both_wedges_silently_tokenized"

    return "unclassified"


def main():
    if not DUMP_PATH.exists():
        print(
            f"missing {DUMP_PATH} -- run: cargo run -p chematic-mol --example stereo2d_fixture_dump "
            f"> {DUMP_PATH}",
            file=sys.stderr,
        )
        sys.exit(1)

    fixtures = [json.loads(line) for line in DUMP_PATH.read_text().splitlines() if line.strip()]

    rows = []
    bucket_counts = {}
    for fx in fixtures:
        mol_block = fx.get("mol_block", "")
        rdkit_result = rdkit_read(mol_block) if mol_block else {"parsed": False}
        bucket = classify(fx, rdkit_result)
        bucket_counts[bucket] = bucket_counts.get(bucket, 0) + 1
        rows.append(
            {
                "id": fx["id"],
                "mechanism": fx["mechanism"],
                "description": fx["description"],
                "chematic": {
                    k: fx.get(k)
                    for k in (
                        "assign_stereo_from_2d_result",
                        "assign_ez_from_2d_result",
                        "post_apply_stereo_from_2d_fields",
                        "naive_smiles_write",
                        "naive_canonical_smiles",
                        "chirality_reached_writer",
                        "note",
                    )
                },
                "rdkit": rdkit_result,
                "failure_bucket": bucket,
            }
        )

    unexpected = [r for r in rows if r["failure_bucket"].startswith("unexpected")]

    summary = {
        "rdkit_pinned_commit": RDKIT_PINNED_COMMIT,
        "rdkit_python_version": __import__("rdkit").rdBase.rdkitVersion,
        "fixture_count": len(rows),
        "bucket_counts": bucket_counts,
        "unexplained_count": len(unexpected),
        "rows": rows,
    }

    SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2))

    print(f"fixtures: {len(rows)}")
    print("bucket counts:")
    for bucket, count in sorted(bucket_counts.items()):
        print(f"  {bucket}: {count}")
    print(f"unexplained (should be 0): {len(unexpected)}")
    if unexpected:
        for r in unexpected:
            print(f"  UNEXPECTED: {r['id']} -> {r['failure_bucket']}")
    print(f"wrote {SUMMARY_PATH}")


if __name__ == "__main__":
    main()
