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

Fail-closed by design: any fixture whose classification isn't one of the
explicit EXPECTED_BUCKETS below (a whitelist, not a "startswith unexpected"
convention) makes this script exit(1). A mismatched fixture-ID set, a
duplicate fixture ID, or a failed self-test also exits(1). This means a
future contributor who adds a 15th fixture, renames a mechanism, or changes
chematic's runtime behavior in a way this script doesn't already have a
bucket for gets a hard failure here rather than a silently-still-green
"0 unexplained" result.

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
EXPECTED_RDKIT_VERSION = "2026.03.3"

# The full fixture-ID set this script knows how to classify. Adding, renaming,
# or removing a fixture in stereo2d_fixture_dump.rs without updating this set
# is caught below (missing/extra IDs), not silently ignored.
EXPECTED_FIXTURE_IDS = {
    "tetrahedral_3heavy_implicit_h",
    "tetrahedral_4neighbors_explicit_h",
    "tetrahedral_4heavy_no_h",
    "solid_wedge_only",
    "dashed_wedge_only",
    "wedge_atom_order_reversed",
    "multiple_stereocenters",
    "no_wedge_negative_control",
    "cip_priority_tie",
    "degenerate_2d_coordinates",
    "ez_geometry_2butene",
    "terminal_alkene_propene",
    "contradictory_wedge_annotations",
    "coord_atom_count_mismatch",
}

# Whitelist of buckets classify() may legitimately return. Anything else
# (including "unclassified", any "unexpected_*" bucket, or "rdkit_parse_failed")
# is, by construction, not in this set and therefore fails the run.
EXPECTED_BUCKETS = {
    "correctly_no_stereo_both_agree",
    "degenerate_coords_correctly_yields_no_stereo",
    "silent_result_from_corrupted_fallback_positions_not_error",
    "rs_not_computed_3heavy_implicit_h_gap",
    "rs_computed_but_writer_emits_meaningless_bond_direction_token",
    "rs_not_computed_despite_rdkit_success",
    "wedge_atom_order_reversed_agrees_with_rdkit_on_same_file",
    "wedge_atom_order_reversed_disagrees_with_rdkit_on_same_file",
    "wedge_atom_order_reversed_chematic_only",
    "ez_computed_but_no_bond_direction_for_writer",
    "no_consistency_check_both_wedges_silently_tokenized",
}


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
    computed/wrote vs what RDKit's default parse produces.

    Every branch that returns an EXPECTED_BUCKETS member first checks the
    concrete evidence that bucket claims to describe (not just the mechanism
    name) -- e.g. "no_consistency_check_both_wedges_silently_tokenized" only
    fires if the naive SMILES actually contains two direction tokens, not
    just because the fixture's mechanism string says "contradictory_wedges".
    A mismatch falls through to an "unexpected_*" bucket, which is NOT in
    EXPECTED_BUCKETS and therefore fails the run in main().
    """

    mech = fixture["mechanism"]
    chematic_rs = fixture.get("assign_stereo_from_2d_result", [])
    chematic_ez = fixture.get("assign_ez_from_2d_result", [])
    chirality_reached_writer = fixture.get("chirality_reached_writer", False)
    naive_smiles = fixture.get("naive_smiles_write", "") or ""
    direction_token_count = naive_smiles.count("/") + naive_smiles.count("\\")

    if not rdkit_result.get("parsed"):
        return "rdkit_parse_failed"

    rdkit_found_stereo = rdkit_result.get("any_chiral_tag") or rdkit_result.get(
        "any_bond_stereo"
    )

    if mech in ("negative_control", "terminal_alkene", "cip_priority_tie"):
        if not chematic_rs and not chematic_ez and not rdkit_found_stereo:
            return "correctly_no_stereo_both_agree"
        return "unexpected_disagreement_on_no_stereo_case"

    if mech == "degenerate_coordinates":
        if not chematic_rs and not rdkit_result.get("any_chiral_tag"):
            return "degenerate_coords_correctly_yields_no_stereo"
        return "unexpected_result_on_degenerate_coords"

    if mech == "coord_atom_count_mismatch":
        # This bucket's whole claim is "a nonempty, non-error assignment came
        # out of corrupted-fallback geometry" -- if that's not what actually
        # happened (e.g. a future fix makes this correctly return empty),
        # this row must NOT silently pass as the same known-bad bucket.
        if chematic_rs:
            return "silent_result_from_corrupted_fallback_positions_not_error"
        return "unexpected_coord_mismatch_result"

    if mech == "tetrahedral_3heavy_implicit_h":
        if not chematic_rs and rdkit_result.get("any_chiral_tag"):
            return "rs_not_computed_3heavy_implicit_h_gap"
        return "unexpected_3heavy_result"

    if mech in ("tetrahedral_4neighbors", "tetrahedral_4heavy_no_h", "solid_wedge", "dashed_wedge"):
        if chematic_rs and not chirality_reached_writer:
            if direction_token_count > 0:
                return "rs_computed_but_writer_emits_meaningless_bond_direction_token"
            return "unexpected_tetrahedral_result"
        if not chematic_rs and rdkit_result.get("any_chiral_tag"):
            return "rs_not_computed_despite_rdkit_success"
        return "unexpected_tetrahedral_result"

    if mech == "wedge_atom_order_reversed":
        # This fixture's own file is non-standard (atom1=substituent for the
        # wedge bond), so "agrees with RDKit's reading of the SAME file" is
        # the right bar -- not "agrees with the standard-order fixture",
        # which is a separate cross-check done in the report by comparing
        # against tetrahedral_4neighbors_explicit_h's row directly.
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
        if chematic_ez and direction_token_count == 0:
            return "ez_computed_but_no_bond_direction_for_writer"
        return "unexpected_ez_result"

    if mech == "contradictory_wedges":
        # The bucket's claim is specifically "both wedges got tokenized
        # independently" -- verify two direction tokens actually appear,
        # not just that the mechanism ran.
        if direction_token_count >= 2:
            return "no_consistency_check_both_wedges_silently_tokenized"
        return "unexpected_contradictory_wedge_result"

    return "unclassified"


def _self_test():
    """Positive controls for the fail-closed machinery itself: prove that an
    unrecognized mechanism and a duplicate-ID input are actually caught,
    rather than trusting that they would be. Runs before any RDKit call, no
    network/subprocess needed -- a plain assertion failure here means the
    fail-closed logic itself is broken and must hard-crash immediately."""

    # Control A: an unknown mechanism must classify as "unclassified", which
    # must NOT be in the bucket whitelist (i.e. it would fail the run).
    bogus_fixture = {
        "mechanism": "some_mechanism_that_does_not_exist",
        "assign_stereo_from_2d_result": [],
        "assign_ez_from_2d_result": [],
        "chirality_reached_writer": False,
        "naive_smiles_write": "C",
    }
    bogus_bucket = classify(bogus_fixture, {"parsed": True, "any_chiral_tag": False, "any_bond_stereo": False})
    assert bogus_bucket == "unclassified", f"self-test A: expected 'unclassified', got {bogus_bucket!r}"
    assert bogus_bucket not in EXPECTED_BUCKETS, "self-test A: 'unclassified' leaked into EXPECTED_BUCKETS"

    # Control B: an "unexpected_*" bucket produced when a claimed bucket's
    # evidence doesn't hold (e.g. contradictory-wedge claimed but only one
    # direction token actually present) must also not be in the whitelist.
    weak_evidence_fixture = {
        "mechanism": "contradictory_wedges",
        "assign_stereo_from_2d_result": [],
        "assign_ez_from_2d_result": [],
        "chirality_reached_writer": False,
        "naive_smiles_write": "C(/F)Cl",  # only ONE direction token, not two
    }
    weak_bucket = classify(weak_evidence_fixture, {"parsed": True, "any_chiral_tag": False, "any_bond_stereo": False})
    assert weak_bucket == "unexpected_contradictory_wedge_result", (
        f"self-test B: expected the evidence-check to reject weak input, got {weak_bucket!r}"
    )
    assert weak_bucket not in EXPECTED_BUCKETS, "self-test B: an 'unexpected_*' bucket leaked into EXPECTED_BUCKETS"

    # Control C: duplicate fixture IDs must be detected by the same check
    # main() runs on the real dump.
    dup_ids = ["a", "b", "a"]
    assert len(set(dup_ids)) != len(dup_ids), "self-test C: duplicate-ID fixture itself has no duplicates"


def main():
    _self_test()

    if not DUMP_PATH.exists():
        print(
            f"missing {DUMP_PATH} -- run: cargo run -p chematic-mol --example stereo2d_fixture_dump "
            f"> {DUMP_PATH}",
            file=sys.stderr,
        )
        sys.exit(1)

    fixtures = [json.loads(line) for line in DUMP_PATH.read_text().splitlines() if line.strip()]

    fail = False

    ids = [fx["id"] for fx in fixtures]
    if len(set(ids)) != len(ids):
        seen = set()
        dupes = [i for i in ids if i in seen or seen.add(i)]
        print(f"FATAL: duplicate fixture IDs in {DUMP_PATH}: {dupes}", file=sys.stderr)
        sys.exit(1)

    actual_ids = set(ids)
    missing = EXPECTED_FIXTURE_IDS - actual_ids
    extra = actual_ids - EXPECTED_FIXTURE_IDS
    if missing or extra:
        if missing:
            print(f"FATAL: expected fixture IDs missing from the dump: {sorted(missing)}", file=sys.stderr)
        if extra:
            print(
                f"FATAL: dump contains fixture IDs this script doesn't know how to classify "
                f"(update EXPECTED_FIXTURE_IDS and classify() in this script): {sorted(extra)}",
                file=sys.stderr,
            )
        sys.exit(1)

    rdkit_version = __import__("rdkit").rdBase.rdkitVersion
    version_mismatch = rdkit_version != EXPECTED_RDKIT_VERSION
    if version_mismatch:
        print(
            f"WARNING: installed rdkit=={rdkit_version} != pinned reference {EXPECTED_RDKIT_VERSION} "
            f"(pinned commit {RDKIT_PINNED_COMMIT}) -- results below may not reflect the pinned behavior.",
            file=sys.stderr,
        )

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

    unexplained = [r for r in rows if r["failure_bucket"] not in EXPECTED_BUCKETS]
    if unexplained:
        fail = True

    summary = {
        "rdkit_pinned_commit": RDKIT_PINNED_COMMIT,
        "rdkit_python_version": rdkit_version,
        "rdkit_version_mismatch": version_mismatch,
        "fixture_count": len(rows),
        "bucket_counts": bucket_counts,
        "unexplained_count": len(unexplained),
        "rows": rows,
    }

    SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2))

    print(f"fixtures: {len(rows)}")
    print("bucket counts:")
    for bucket, count in sorted(bucket_counts.items()):
        print(f"  {bucket}: {count}")
    print(f"unexplained (must be 0): {len(unexplained)}")
    if unexplained:
        for r in unexplained:
            print(f"  UNEXPLAINED: {r['id']} -> {r['failure_bucket']}")
    print(f"wrote {SUMMARY_PATH}")

    if fail:
        print("FATAL: unexplained fixtures present -- see above.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
