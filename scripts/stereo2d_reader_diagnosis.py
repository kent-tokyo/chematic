#!/usr/bin/env python3
"""Differential validation for the 2D wedge/hash *reader-integration* work
(`chematic_mol::read_mol_with_diagnostics` / `read_mol_v3000_with_diagnostics`
wiring `chematic_perception::apply_local_parity_from_wedges_with_diagnostics`
into the V2000/V3000 MOL readers).

Companion to, not a replacement for, `scripts/stereo2d_diagnosis.py` (the
original P1-A0 diagnosis of the parity math itself, still frozen and
accurate -- its sign convention was already calibrated against RDKit in
docs/rfcs/stereo2d_local_parity_calibration.md). This script instead validates
that wiring that math into the production reader path produces the right
end-to-end results on the NEW fixture categories this integration adds:
V2000/V3000 agreement, renumbering/reflection/rotation invariance,
round-trip losslessness, ring/charged/isotopic stereocenters, and the new
structured StereoDiagnostic API.

Reads validation/results/stereo2d_reader_fixture_dump.jsonl (produced by
`cargo run -p chematic-mol --example stereo2d_reader_integration_fixture_dump`)
and, for each fixture's raw MOL block, re-parses it with RDKit to get two
SEPARATE, never-conflated comparisons:

1. Local-parity agreement: chematic's `Atom.chirality` (CW/CCW/None) vs
   RDKit's raw chiral tag from `Chem.MolFromMolBlock(block, sanitize=False)`
   -- unconditional, no CIP ranking involved on either side (matches RDKit's
   own `assignChiralTypesFromBondDirs`, which runs regardless of `sanitize`).
2. Final CIP R/S agreement: chematic's accurate-CIP-mode label
   (`chematic_chem::assign_cip_with_mode(mol, CipMode::Accurate)`) vs
   RDKit's `rdCIPLabeler.AssignCIPLabels` on a sanitized parse.

Fail-closed: exits 1 on an unexpected fixture-ID set, a duplicate ID, or a
result outside this script's own expected-outcome table for that fixture.

Run:
    .venv/bin/python scripts/stereo2d_reader_diagnosis.py
"""

import io
import json
import sys
from contextlib import redirect_stderr
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdCIPLabeler

ROOT = Path(__file__).resolve().parent.parent
DUMP_PATH = ROOT / "validation" / "results" / "stereo2d_reader_fixture_dump.jsonl"
SUMMARY_PATH = ROOT / "validation" / "results" / "stereo2d_reader_diagnosis_summary.json"

RDKIT_PINNED_COMMIT = "8afba32ec539dcb2369bc84549d802aca3f7eb39"
EXPECTED_RDKIT_VERSION = "2026.03.3"

# Frozen baseline: for each fixture ID, whether RDKit's raw chiral tag
# ("resolved") or the accurate-CIP label ("resolved") is expected to agree
# with chematic, or whether the fixture is intentionally one RDKit itself
# leaves unresolved (e.g. a contradictory-wedge center, or a bond that isn't
# a stereocenter candidate at all). Values: "agree" (chematic and RDKit both
# assign the same non-null tag/label), "both_none" (neither assigns
# anything), "rdkit_none_chematic_none" is the same as "both_none".
#   "chematic_rejects_rdkit_unreliable_fallback": chematic conservatively
#   refuses (a StereoDiagnostic, per wedges_agree_4's tri-state isolated-
#   volume check) while RDKit's own dual-volume fallback silently produces
#   SOME tag/label with no warning -- measured, not assumed: this is the
#   exact same divergence already characterized in
#   docs/rfcs/stereo2d_local_parity_calibration.md's "Multi-wedge consistency"
#   section (RDKit's fallback tag doesn't even agree with itself across
#   similar disagreeing-wedge fixtures there, so it isn't a target to match).
EXPECTED_LOCAL_PARITY_OUTCOME = {
    "valid_wedge_v2000": "agree",
    "valid_wedge_v3000": "agree",
    "contradictory_wedge_v2000": "chematic_rejects_rdkit_unreliable_fallback",
    "atom_renumbered_v2000": "agree",
    "reflected_v2000": "agree",
    "rotated_translated_v2000": "agree",
    "bond_order_reversed_v2000": "agree",
    "charged_n_center_v2000": "agree",
    "isotopic_stereocenter_v2000": "agree",
    "ring_stereocenter_v2000": "agree",
    "multi_stereocenter_v2000": "agree",
    "achiral_negative_control_v2000": "both_none",
    "non_tetrahedral_wedge_v2000": "both_none",
}

# ring_stereocenter_v2000's expected accurate-CIP outcome is "both_none", NOT
# "agree": bromocyclopropane's substituted carbon has a well-defined,
# wedge-derived LOCAL parity (both engines' raw/geometric read agrees, see
# EXPECTED_LOCAL_PARITY_OUTCOME above) but is NOT a genuine CIP stereocenter
# at all -- the ring's own local mirror symmetry (swapping the two
# unsubstituted ring carbons gives the identical molecule) makes it
# achiral. Verified directly: RDKit's `FindMolChiralCenters(mol,
# includeUnassigned=True)` returns `[]` for this molecule, and a freshly
# canonicalized parse (no synthetic wedge) shows `CHI_UNSPECIFIED` on every
# atom. chematic's accurate-CIP engine correctly agrees (also assigns
# nothing) -- this is exactly the CIP-independence local parity is designed
# to demonstrate: a drawing can have a definite geometric parity without
# corresponding to a real stereocenter.
EXPECTED_ACCURATE_CIP_OUTCOME = {
    "valid_wedge_v2000": "agree",
    "valid_wedge_v3000": "agree",
    "contradictory_wedge_v2000": "chematic_rejects_rdkit_unreliable_fallback",
    "atom_renumbered_v2000": "agree",
    "reflected_v2000": "agree",
    "rotated_translated_v2000": "agree",
    "bond_order_reversed_v2000": "agree",
    "charged_n_center_v2000": "agree",
    "isotopic_stereocenter_v2000": "agree",
    "ring_stereocenter_v2000": "both_none",
    "multi_stereocenter_v2000": "agree",
    "achiral_negative_control_v2000": "both_none",
    "non_tetrahedral_wedge_v2000": "both_none",
}

EXPECTED_FIXTURE_IDS = set(EXPECTED_LOCAL_PARITY_OUTCOME)

CHIRAL_TAG_TO_CHEMATIC = {
    Chem.ChiralType.CHI_TETRAHEDRAL_CW: "Clockwise",
    Chem.ChiralType.CHI_TETRAHEDRAL_CCW: "CounterClockwise",
}


def rdkit_raw_chiral_tags(mol_block):
    """Return {atom_idx: "Clockwise"|"CounterClockwise"} from RDKit's raw,
    unconditional chiral-tag assignment (sanitize=False -- per RFC section
    5a, RDKit's chiral-tag step runs regardless of the sanitize flag; only
    CIP labeling is gated on it)."""
    stderr_buf = io.StringIO()
    with redirect_stderr(stderr_buf):
        mol = Chem.MolFromMolBlock(mol_block, sanitize=False)
    if mol is None:
        return {}, stderr_buf.getvalue()
    tags = {}
    for atom in mol.GetAtoms():
        tag = atom.GetChiralTag()
        if tag in CHIRAL_TAG_TO_CHEMATIC:
            tags[atom.GetIdx()] = CHIRAL_TAG_TO_CHEMATIC[tag]
    return tags, stderr_buf.getvalue()


def rdkit_accurate_cip_labels(mol_block):
    """Return {atom_idx: "R"|"S"} from a sanitized parse + rdCIPLabeler."""
    stderr_buf = io.StringIO()
    with redirect_stderr(stderr_buf):
        mol = Chem.MolFromMolBlock(mol_block, sanitize=True)
        if mol is None:
            return {}, stderr_buf.getvalue()
        rdCIPLabeler.AssignCIPLabels(mol)
    labels = {}
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            labels[atom.GetIdx()] = atom.GetProp("_CIPCode")
    return labels, stderr_buf.getvalue()


def classify_local_parity(row, rdkit_tags):
    chematic_tags = {c["atom"]: c["chirality"] for c in row["chirality"]}
    if not chematic_tags and not rdkit_tags:
        return "both_none"
    if chematic_tags == rdkit_tags:
        return "agree"
    if not chematic_tags and rdkit_tags:
        return "chematic_rejects_rdkit_unreliable_fallback"
    return f"disagree(chematic={chematic_tags}, rdkit={rdkit_tags})"


def classify_accurate_cip(row, rdkit_labels):
    cip = row.get("accurate_cip_labels") or []
    chematic_labels = {c["atom"]: c["cip_code"] for c in cip}
    if not chematic_labels and not rdkit_labels:
        return "both_none"
    if chematic_labels == rdkit_labels:
        return "agree"
    if not chematic_labels and rdkit_labels:
        return "chematic_rejects_rdkit_unreliable_fallback"
    return f"disagree(chematic={chematic_labels}, rdkit={rdkit_labels})"


def main():
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

        rdkit_tags, tag_stderr = rdkit_raw_chiral_tags(row["mol_block"])
        rdkit_labels, cip_stderr = rdkit_accurate_cip_labels(row["mol_block"])

        local_parity_outcome = classify_local_parity(row, rdkit_tags)
        accurate_cip_outcome = classify_accurate_cip(row, rdkit_labels)

        expected_local = EXPECTED_LOCAL_PARITY_OUTCOME[row["id"]]
        expected_cip = EXPECTED_ACCURATE_CIP_OUTCOME[row["id"]]

        local_ok = local_parity_outcome == expected_local
        cip_ok = accurate_cip_outcome == expected_cip
        if not local_ok or not cip_ok:
            unexplained.append(
                {
                    "id": row["id"],
                    "expected_local_parity": expected_local,
                    "actual_local_parity": local_parity_outcome,
                    "expected_accurate_cip": expected_cip,
                    "actual_accurate_cip": accurate_cip_outcome,
                    "rdkit_stderr": (tag_stderr + cip_stderr).strip(),
                }
            )

        results.append(
            {
                "id": row["id"],
                "description": row["description"],
                "chematic_chirality": row["chirality"],
                "rdkit_raw_chiral_tags": rdkit_tags,
                "local_parity_outcome": local_parity_outcome,
                "chematic_accurate_cip": row.get("accurate_cip_labels"),
                "rdkit_cip_labels": rdkit_labels,
                "accurate_cip_outcome": accurate_cip_outcome,
                "stereo_diagnostics": row.get("stereo_diagnostics", []),
                "rdkit_stderr": (tag_stderr + cip_stderr).strip(),
            }
        )

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

    local_agree = sum(1 for r in results if r["local_parity_outcome"] in ("agree", "both_none"))
    cip_agree = sum(1 for r in results if r["accurate_cip_outcome"] in ("agree", "both_none"))
    print(
        f"Local-parity agreement: {local_agree}/{len(results)}  "
        f"Accurate-CIP agreement: {cip_agree}/{len(results)}"
    )

    if unexplained:
        print(f"FAIL: {len(unexplained)} fixture(s) did not match their expected outcome:", file=sys.stderr)
        for u in unexplained:
            print(f"  {u}", file=sys.stderr)
        sys.exit(1)

    print(f"OK: {len(results)}/{len(results)} fixtures matched their expected outcome, 0 unexplained.")
    sys.exit(0)


if __name__ == "__main__":
    main()
