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

# Frozen baseline: the exact bucket each fixture ID must land in. This is
# stricter than "some bucket in the EXPECTED_BUCKETS whitelist below" --
# several buckets in that whitelist describe mutually exclusive outcomes for
# the SAME mechanism (e.g. wedge_atom_order_reversed's "agrees with RDKit" vs
# "disagrees with RDKit" vs "chematic only") precisely because classify()
# has to be able to express all three depending on runtime evidence. Without
# this per-ID map, a fixture whose bucket silently drifts from one whitelist
# member to another (e.g. wedge_atom_order_reversed flipping from "agrees"
# to "disagrees", or tetrahedral_4heavy_no_h regressing from "computed" to
# "not computed despite RDKit success") would still report 0 unexplained.
# Adding, renaming, or removing a fixture in stereo2d_fixture_dump.rs without
# updating this map is caught below (missing/extra IDs), not silently ignored.
# 2026-07-25 update: three rows below changed bucket after the P1-A0
# reader-integration PR shipped (`chematic_mol::read_mol_with_diagnostics`
# now calls `apply_local_parity_from_wedges_with_diagnostics` from inside
# `parse_mol_with_coords` -- exactly the fix this diagnosis's own §7 called
# for). This example's `main()` still calls `parse_mol_with_coords` on each
# fixture, so `chirality_reached_writer` (computed from the SAME parsed
# `Molecule`, before this script's own subsequent `apply_stereo_from_2d`
# call) now flips from `False` to `True` for any fixture whose wedge(s)
# resolve to a valid local parity. Re-verified against the regenerated
# `stereo2d_fixture_dump.jsonl`, not assumed:
#   - tetrahedral_4neighbors_explicit_h / tetrahedral_4heavy_no_h: chematic
#     already computed R/S via the old CIP engine; the only thing that
#     changed is `chirality` now ALSO reaches the writer correctly (was the
#     literal bug this diagnosis reported -- "computed but writer emits a
#     meaningless bond-direction token" -- now fixed).
#   - contradictory_wedge_annotations: its two "up" wedges (tripod_atoms
#     geometry, F and Cl) turn out to AGREE in isolation under the
#     calibrated per-wedge consistency check (`wedges_agree_3` --
#     "same direction" is not itself the discriminator, see
#     docs/stereo2d_local_parity_calibration.md's "Multi-wedge consistency"
#     section) -- so this specific fixture was never a genuine contradiction,
#     and chematic now correctly resolves a definite local parity for it
#     instead of the old "no consistency check at all" bug. A geometry that
#     IS a genuine per-wedge disagreement (e.g. the 4-heavy quad_positions
#     layout) is separately covered by the new reader-integration corpus's
#     own `contradictory_wedge_v2000` fixture, which chematic still
#     correctly rejects.
EXPECTED_BUCKET_BY_ID = {
    "tetrahedral_3heavy_implicit_h": "rs_not_computed_3heavy_implicit_h_gap",
    "tetrahedral_4neighbors_explicit_h": "rs_computed_and_chirality_now_reaches_writer",
    "tetrahedral_4heavy_no_h": "rs_computed_and_chirality_now_reaches_writer",
    "solid_wedge_only": "rs_not_computed_despite_rdkit_success",
    "dashed_wedge_only": "rs_not_computed_despite_rdkit_success",
    "wedge_atom_order_reversed": "wedge_atom_order_reversed_agrees_with_rdkit_on_same_file",
    "multiple_stereocenters": "rs_not_computed_3heavy_implicit_h_gap",
    "no_wedge_negative_control": "correctly_no_stereo_both_agree",
    "cip_priority_tie": "correctly_no_stereo_both_agree",
    "degenerate_2d_coordinates": "degenerate_coords_correctly_yields_no_stereo",
    "ez_geometry_2butene": "ez_computed_but_no_bond_direction_for_writer",
    "terminal_alkene_propene": "correctly_no_stereo_both_agree",
    "contradictory_wedge_annotations": "wedges_agree_in_isolation_chirality_now_resolved",
    "coord_atom_count_mismatch": "silent_result_from_corrupted_fallback_positions_not_error",
}

EXPECTED_FIXTURE_IDS = set(EXPECTED_BUCKET_BY_ID)

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
    "rs_computed_and_chirality_now_reaches_writer",
    "wedges_agree_in_isolation_chirality_now_resolved",
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
        if chematic_rs and chirality_reached_writer:
            # The reader-integration fix: chirality (from the new
            # apply_local_parity_from_wedges_with_diagnostics wiring) now
            # correctly reaches the writer alongside the old engine's R/S.
            return "rs_computed_and_chirality_now_reaches_writer"
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
        # not just that the mechanism ran. This is pre-PR#140/pre-reader-
        # integration behavior; kept as a reachable bucket in case a future
        # fixture reproduces it, but this specific corpus's fixture no
        # longer lands here (see the two buckets below).
        if direction_token_count >= 2:
            return "no_consistency_check_both_wedges_silently_tokenized"
        if chirality_reached_writer:
            # This fixture's two "up" wedges agree in isolation under the
            # calibrated per-wedge consistency check (wedges_agree_3) -- not
            # a genuine contradiction, so chematic now resolves a definite
            # local parity instead of silently doing nothing.
            return "wedges_agree_in_isolation_chirality_now_resolved"
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

    # Control D: the frozen per-ID baseline must actually pin ONE bucket, not
    # just any whitelisted one -- prove that wedge_atom_order_reversed's
    # three mutually-exclusive possible outcomes are all in EXPECTED_BUCKETS
    # (so the whitelist alone could not tell them apart) while
    # EXPECTED_BUCKET_BY_ID pins exactly one of them. This is what actually
    # catches a result silently drifting between "agrees with RDKit" and
    # "disagrees with RDKit" across runs.
    reversed_wedge_outcomes = {
        "wedge_atom_order_reversed_agrees_with_rdkit_on_same_file",
        "wedge_atom_order_reversed_disagrees_with_rdkit_on_same_file",
        "wedge_atom_order_reversed_chematic_only",
    }
    assert reversed_wedge_outcomes <= EXPECTED_BUCKETS, (
        "self-test D: expected all three wedge_atom_order_reversed outcomes to be individually valid buckets"
    )
    pinned = EXPECTED_BUCKET_BY_ID["wedge_atom_order_reversed"]
    assert pinned in reversed_wedge_outcomes, "self-test D: frozen baseline for wedge_atom_order_reversed is not one of its own possible outcomes"
    other_outcomes = reversed_wedge_outcomes - {pinned}
    assert other_outcomes, "self-test D: no alternative outcome to distinguish from the pinned one"
    for other in other_outcomes:
        assert other != pinned, "self-test D: pinned baseline must differ from the alternative outcomes it's meant to catch"


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
    baseline_drift = []
    for fx in fixtures:
        mol_block = fx.get("mol_block", "")
        rdkit_result = rdkit_read(mol_block) if mol_block else {"parsed": False}
        bucket = classify(fx, rdkit_result)
        bucket_counts[bucket] = bucket_counts.get(bucket, 0) + 1
        expected_bucket = EXPECTED_BUCKET_BY_ID[fx["id"]]
        if bucket != expected_bucket:
            baseline_drift.append((fx["id"], expected_bucket, bucket))
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
                "expected_bucket": expected_bucket,
            }
        )

    # Two independent fail-closed checks: (1) is the bucket even a
    # recognized/whitelisted outcome at all, and (2) -- stricter -- does it
    # match the EXACT bucket this fixture was frozen at, not just some
    # whitelisted one. (2) catches drift between mutually-exclusive
    # whitelisted outcomes for the same fixture (e.g. wedge_atom_order_reversed
    # silently flipping from "agrees with RDKit" to "disagrees") that (1) alone
    # would miss.
    unexplained = [r for r in rows if r["failure_bucket"] not in EXPECTED_BUCKETS]
    if unexplained or baseline_drift:
        fail = True

    summary = {
        "rdkit_pinned_commit": RDKIT_PINNED_COMMIT,
        "rdkit_python_version": rdkit_version,
        "rdkit_version_mismatch": version_mismatch,
        "fixture_count": len(rows),
        "bucket_counts": bucket_counts,
        "unexplained_count": len(unexplained),
        "baseline_drift_count": len(baseline_drift),
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
    print(f"baseline drift vs frozen EXPECTED_BUCKET_BY_ID (must be 0): {len(baseline_drift)}")
    if baseline_drift:
        for fx_id, expected, actual in baseline_drift:
            print(f"  DRIFT: {fx_id} -> expected {expected!r}, got {actual!r}")
    print(f"wrote {SUMMARY_PATH}")

    if fail:
        print("FATAL: unexplained and/or baseline-drifted fixtures present -- see above.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
