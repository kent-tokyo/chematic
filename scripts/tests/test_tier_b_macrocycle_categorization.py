#!/usr/bin/env python3
"""
Regression test for the Tier B manifest's macrocycle categorization fix
(A1's Finding 2: `gen_pipeline_v2_vs_rdkit_tier_b_manifest.py` used to
hardcode every accepted molecule as "drug_like" with no ring-size check at
all, silently mis-tagging 6 large bis-pyridinium macrocycles that a later
benchmark round found by auditing ring sizes directly).

Two checks, not one:

1. POSITIVE: the 6 specific molecules this round identified
   (`chembl_tier_b_0009/0023/0028/0029/0030/0034`, 28-32 atom rings) are
   tagged "macrocycle", not "drug_like".
2. NEGATIVE (the one that matters -- catches "moved everything" or "moved
   the wrong ones" mistakes a hand-edit could make): independently
   re-derive macrocycle status for ALL 200 molecules directly from their
   SMILES via RDKit ring info (the same rule
   `gen_pipeline_v2_vs_rdkit_tier_b_manifest.py` now applies), and assert
   it matches the manifest's stored `primary_category` for every single
   molecule -- not just the 6 known ones. A fix that accidentally
   recategorized an unrelated drug_like molecule, or missed a 7th
   macrocycle, fails here even though it would pass a check that only
   looks at the 6 named molecules.

Usage: python scripts/tests/test_tier_b_macrocycle_categorization.py
"""

import json
import os
import sys

MACROCYCLE_MIN = 9  # must match chematic_3d::etkdg_knowledge::classify::MACROCYCLE_MIN
KNOWN_MACROCYCLES = {
    "chembl_tier_b_0009",
    "chembl_tier_b_0023",
    "chembl_tier_b_0028",
    "chembl_tier_b_0029",
    "chembl_tier_b_0030",
    "chembl_tier_b_0034",
}

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
MANIFEST_PATH = os.path.join(
    ROOT, "validation", "manifests", "pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"
)


def positive_control(molecules):
    for m in molecules:
        if m["name"] in KNOWN_MACROCYCLES:
            assert m["primary_category"] == "macrocycle", (
                f"{m['name']} (a known 28-32-atom-ring macrocycle) is tagged "
                f"{m['primary_category']!r}, not 'macrocycle'"
            )
    tagged = {m["name"] for m in molecules if m["primary_category"] == "macrocycle"}
    assert tagged == KNOWN_MACROCYCLES, (
        f"exactly the 6 known macrocycles should be tagged 'macrocycle', got {sorted(tagged)}"
    )


def negative_control(molecules):
    """Every molecule's stored category must match an independent,
    freshly-computed ring-size check -- not just the 6 known cases."""
    from rdkit import Chem

    mismatches = []
    for m in molecules:
        mol = Chem.MolFromSmiles(m["smiles"])
        assert mol is not None, f"{m['name']}: stored SMILES failed to re-parse"
        ring_sizes = [len(r) for r in mol.GetRingInfo().AtomRings()]
        is_macrocycle = any(s >= MACROCYCLE_MIN for s in ring_sizes)
        expected = "macrocycle" if is_macrocycle else "drug_like"
        if m["primary_category"] != expected:
            mismatches.append((m["name"], m["primary_category"], expected))

    assert not mismatches, (
        f"{len(mismatches)} molecule(s) have a stored primary_category that disagrees "
        f"with an independent ring-size recomputation: {mismatches}"
    )


def main():
    with open(MANIFEST_PATH) as f:
        manifest = json.load(f)
    molecules = manifest["molecules"]
    assert len(molecules) == 200, f"expected 200 molecules, got {len(molecules)}"

    positive_control(molecules)
    negative_control(molecules)

    print(f"OK: {len(molecules)} molecules, 6 correctly tagged macrocycle, "
          "194 correctly tagged drug_like, independently re-verified via RDKit ring info.")


if __name__ == "__main__":
    sys.exit(main())
