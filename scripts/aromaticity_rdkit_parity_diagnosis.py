#!/usr/bin/env python3
"""Diag/aromaticity-rdkit-parity: cross-reference chematic's aromaticity
perception and kekulization against RDKit on a frozen, deliberately-built
SMILES fixture corpus emitted by
`cargo run -p chematic-perception --example aromaticity_rdkit_parity_dump`.

Diagnostic only. Does not change chematic's production behavior. Reads
validation/results/aromaticity_rdkit_parity_fixture_dump.jsonl (one row per
fixture) and, for each fixture's SMILES, independently parses it with RDKit
(pinned commit 8afba32ec539dcb2369bc84549d802aca3f7eb39 / rdkit==2026.03.3 in
this project's venv) to get RDKit's own per-atom/per-bond aromatic flags,
ring info, and Kekule bond assignment, then classifies the (chematic, RDKit)
pair into a failure bucket -- per fixture, but every bucket's branch inspects
concrete per-atom/per-bond evidence, not just the fixture's category label.
Writes validation/results/aromaticity_rdkit_parity_diagnosis_summary.json.

Fail-closed by design, matching scripts/stereo2d_diagnosis.py's style: a
bucket whitelist (not a "startswith unexpected" convention), an
EXPECTED_BUCKET_BY_ID frozen per-fixture baseline (stricter than the
whitelist -- several buckets describe mutually-exclusive outcomes of the
same mechanism), a fixture-ID-set check (missing/extra IDs abort), duplicate-
ID detection, and self-tests of the fail-closed logic itself (hand-injected
failure modes) run before touching any real data. sys.exit(1) on any of
these.

Run:
    .venv/bin/python scripts/aromaticity_rdkit_parity_diagnosis.py
"""

import json
import sys
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parent.parent
DUMP_PATH = ROOT / "validation" / "results" / "aromaticity_rdkit_parity_fixture_dump.jsonl"
SUMMARY_PATH = ROOT / "validation" / "results" / "aromaticity_rdkit_parity_diagnosis_summary.json"

RDKIT_PINNED_COMMIT = "8afba32ec539dcb2369bc84549d802aca3f7eb39"
EXPECTED_RDKIT_VERSION = "2026.03.3"

_BOND_TYPE_STR = {
    Chem.BondType.SINGLE: "Single",
    Chem.BondType.DOUBLE: "Double",
    Chem.BondType.TRIPLE: "Triple",
    Chem.BondType.AROMATIC: "Aromatic",
}

# Frozen baseline: the exact bucket each fixture ID must land in. Stricter
# than "some bucket in the EXPECTED_BUCKETS whitelist below" for the same
# reason as scripts/stereo2d_diagnosis.py's own EXPECTED_BUCKET_BY_ID: several
# buckets describe mutually exclusive outcomes of the SAME mechanism (e.g.
# "matches RDKit because both correctly disallow it" vs "matches RDKit only
# because an unverified flag survived"), so the whitelist alone cannot tell a
# silent drift between them apart from a stable result.
EXPECTED_BUCKET_BY_ID = {
    "benzene": "matches_rdkit_kekule_valid_alternate",
    "cyclohexane": "both_correctly_nonaromatic",
    "cyclohexadiene_1_3": "both_correctly_nonaromatic",
    "cyclooctatetraene": "both_correctly_nonaromatic",
    "pyridine": "matches_rdkit_exact_kekule",
    "pyrimidine": "matches_rdkit_exact_kekule",
    "furan": "matches_rdkit_exact_kekule",
    "thiophene": "matches_rdkit_exact_kekule",
    "pyrrole": "matches_rdkit_exact_kekule",
    "n_methylpyrrole": "matches_rdkit_exact_kekule",
    "imidazole": "matches_rdkit_exact_kekule",
    "pyrazole": "matches_rdkit_exact_kekule",
    "oxazole": "matches_rdkit_exact_kekule",
    "thiazole": "matches_rdkit_exact_kekule",
    "isoxazole": "matches_rdkit_exact_kekule",
    "triazole_1_2_3": "matches_rdkit_exact_kekule",
    "tetrazole": "matches_rdkit_exact_kekule",
    "naphthalene": "matches_rdkit_exact_kekule",
    "anthracene": "matches_rdkit_exact_kekule",
    "quinoline": "matches_rdkit_exact_kekule",
    "isoquinoline": "matches_rdkit_exact_kekule",
    "indole": "matches_rdkit_exact_kekule",
    "indolizine": "sssr_bridge_artifact_not_reproduced_docs_stale",
    "purine": "matches_rdkit_exact_kekule",
    "azulene": "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent",
    # kekulize/K1 fix (fix/kekulize-charge-aware-k1): these 6 fixtures used to
    # hard-fail kekulize() outright (bucket
    # "kekulize_fails_atom_bond_flags_survive_coincidentally") -- see
    # docs/aromaticity_rdkit_parity_rfc.md §1 root causes A-D. `kekulize()`
    # now succeeds and its bond-by-bond Kekule assignment is verified
    # byte-identical to RDKit's own choice (kekule_bond_mismatch_pairs == []
    # for all 6, checked directly against the diagnosis summary JSON, not
    # assumed). What moves them into
    # "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent"
    # instead of "matches_rdkit_exact_kekule" is unrelated to the K1 fix
    # itself: `build_molecule_from_model`'s atom-flag rebuild loop only ever
    # *promotes* `atom.aromatic` to true, never demotes a stale pre-existing
    # true when the Huckel model disagrees (RFC §1b) -- deliberately out of
    # scope here, tracked separately as "K2". Same bucket/shape as
    # selenophene/azulene above, for the same deferred reason.
    "tropylium_cation": "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent",
    "cyclopentadienyl_anion": "matches_rdkit_exact_kekule",
    "imidazolium": "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent",
    "pyridinium": "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent",
    "pyrylium": "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent",
    "selenophene": "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent",
    "tellurophene": "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent",
    "phosphole": "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent",
    "borole": "both_correctly_nonaromatic",
    "borazine": "both_correctly_nonaromatic",
    "pyridone_2": "matches_rdkit_exact_kekule_exocyclic_bond_excluded_correctly",
    "tropone": "matches_rdkit_exact_kekule_exocyclic_bond_excluded_correctly",
    "benzoquinone_1_4": "both_correctly_nonaromatic",
    "cyclopentadienone": "both_correctly_nonaromatic",
    "thiophene_1_oxide": "both_correctly_nonaromatic",
}

EXPECTED_FIXTURE_IDS = set(EXPECTED_BUCKET_BY_ID)

# Whitelist of buckets classify() may legitimately return. Anything else
# (including "unclassified" or any "unexpected_*" bucket) is, by
# construction, not in this set and therefore fails the run.
EXPECTED_BUCKETS = {
    "matches_rdkit_exact_kekule",
    "matches_rdkit_exact_kekule_exocyclic_bond_excluded_correctly",
    "matches_rdkit_kekule_valid_alternate",
    "both_correctly_nonaromatic",
    "sssr_bridge_artifact_not_reproduced_docs_stale",
    "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent",
    "kekulize_fails_atom_bond_flags_survive_coincidentally",
}


def rdkit_read(smiles):
    """Parse `smiles` with RDKit (default sanitize=True) and report its
    per-atom/per-bond aromaticity, ring membership, and Kekule bond types.
    Atom indices are EXPECTED to match chematic's (both preserve SMILES input
    token order; RDKit never reorders atoms on parse, only on request for a
    *canonical* SMILES output, which this function does not ask for) -- but
    this is an assumption the whole per-index comparison in classify() rests
    on, not a guarantee, so `atom_element` is returned here specifically so
    classify() can verify it index-by-index against chematic's own per-atom
    element field before trusting any other index-aligned comparison.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"parsed": False}

    atom_aromatic = [a.GetIsAromatic() for a in mol.GetAtoms()]
    atom_element = [a.GetSymbol() for a in mol.GetAtoms()]
    bond_aromatic_by_pair = {}
    for b in mol.GetBonds():
        pair = frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        bond_aromatic_by_pair[pair] = b.GetIsAromatic()

    ring_atom_sets = [frozenset(r) for r in mol.GetRingInfo().AtomRings()]
    aromatic_ring_count = sum(
        1 for ring in ring_atom_sets if all(atom_aromatic[i] for i in ring)
    )

    # Kekulize a copy to get RDKit's own chosen Kekule bond assignment,
    # keyed by atom-index pair (order-independent) for direct comparison
    # against chematic's own kekulized_order per bond.
    kek = Chem.Mol(mol)
    had_aromatic_bonds = any(bond_aromatic_by_pair.values())
    kekule_bond_type_by_pair = {}
    if had_aromatic_bonds:
        Chem.Kekulize(kek, clearAromaticFlags=True)
        for b in kek.GetBonds():
            pair = frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
            kekule_bond_type_by_pair[pair] = _BOND_TYPE_STR.get(b.GetBondType(), str(b.GetBondType()))

    return {
        "parsed": True,
        "atom_aromatic": atom_aromatic,
        "atom_element": atom_element,
        "bond_aromatic_by_pair": bond_aromatic_by_pair,
        "ring_atom_sets": ring_atom_sets,
        "aromatic_ring_count": aromatic_ring_count,
        "had_aromatic_bonds": had_aromatic_bonds,
        "kekule_bond_type_by_pair": kekule_bond_type_by_pair,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "_sanitized_mol": mol,
    }


def kekule_structure_is_valid_alternate(rdkit_mol, chematic_kekule_by_pair):
    """Independently verify (via RDKit itself, not by trusting either
    engine's own bond-alternation bookkeeping) that chematic's kekulized
    bond orders describe a chemically valid resonance structure of the SAME
    molecule RDKit parsed -- even if it differs bond-by-bond from the
    specific Kekule structure RDKit's own Kekulize() happened to choose.
    Kekule structures are not unique (naphthalene has 2); a raw bond-by-bond
    diff alone cannot distinguish "different but equally valid" from
    "genuinely wrong."

    Rebuilds every ring bond's type from chematic's assignment on a working
    copy, clears all aromatic flags, and asks RDKit to re-sanitize (which
    re-derives aromaticity from scratch and raises on invalid valence/
    unkekulizable input) and re-canonicalize. Both must succeed and the
    canonical SMILES must match the original for this to return True.
    """
    rw = Chem.RWMol(rdkit_mol)
    for b in rw.GetBonds():
        pair = frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        order = chematic_kekule_by_pair.get(pair)
        if order == "Single":
            b.SetBondType(Chem.BondType.SINGLE)
        elif order == "Double":
            b.SetBondType(Chem.BondType.DOUBLE)
        elif order is not None:
            return False, f"chematic kekulized_order {order!r} is not Single/Double"
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False)
    try:
        Chem.SanitizeMol(rw)
    except Exception as e:  # noqa: BLE001 -- RDKit raises several distinct C++-backed types here
        return False, f"RDKit rejects chematic's kekule structure: {e}"
    try:
        new_canon = Chem.MolToSmiles(rw)
        orig_canon = Chem.MolToSmiles(rdkit_mol)
    except Exception as e:  # noqa: BLE001
        return False, f"failed to canonicalize: {e}"
    if new_canon != orig_canon:
        return False, f"valid structure but a DIFFERENT molecule: {new_canon!r} != {orig_canon!r}"
    return True, "confirmed valid alternate (or identical) Kekule structure by RDKit re-sanitize + canonical match"


def classify(fixture, rdkit_result):
    """Assign one failure/agreement bucket per fixture. Every branch checks
    the concrete per-atom/per-bond/per-ring evidence a bucket claims to
    describe, not just the fixture's category label -- a mismatch falls
    through to an "unexpected_*" bucket, which is NOT in EXPECTED_BUCKETS and
    therefore fails the run in main()."""

    if not rdkit_result.get("parsed"):
        return "rdkit_parse_failed", {}

    n = len(fixture["atoms"])
    chem_atom_arom = [a["default_aromatic"] for a in fixture["atoms"]]
    rdkit_atom_arom = rdkit_result["atom_aromatic"]
    if len(rdkit_atom_arom) != n:
        return "unexpected_atom_count_mismatch", {"chematic_n": n, "rdkit_n": len(rdkit_atom_arom)}

    # Every per-index comparison below (atoms AND bonds-by-pair) assumes
    # chematic's atom i is RDKit's atom i for the SAME input SMILES. That
    # assumption is normally sound (neither parser reorders atoms), but it
    # is exactly the kind of silent precondition that must be verified, not
    # trusted -- a counts-equal-but-order-differs case would otherwise
    # produce a plausible-looking but wrong bucket instead of an honest
    # failure. Check it explicitly, element-by-element, before using the
    # index alignment for anything else.
    chem_elements = [a["element"] for a in fixture["atoms"]]
    rdkit_elements = rdkit_result["atom_element"]
    element_mismatch = [
        i for i in range(n) if chem_elements[i].upper() != rdkit_elements[i].upper()
    ]
    if element_mismatch:
        return "unexpected_atom_order_or_element_mismatch", {
            "element_mismatch_idxs": element_mismatch,
            "chematic_elements": chem_elements,
            "rdkit_elements": rdkit_elements,
        }

    atom_mismatch = [i for i in range(n) if chem_atom_arom[i] != rdkit_atom_arom[i]]

    bond_mismatch = []
    for b in fixture["bonds"]:
        pair = frozenset((b["a1"], b["a2"]))
        rdkit_b_arom = rdkit_result["bond_aromatic_by_pair"].get(pair)
        if rdkit_b_arom is None:
            return "unexpected_bond_not_found_in_rdkit", {"pair": sorted(pair)}
        if b["default_aromatic"] != rdkit_b_arom:
            bond_mismatch.append(list(pair))

    ring_count_match = fixture["raw_sssr_ring_count"] == len(rdkit_result["ring_atom_sets"])
    aromatic_ring_count_match = fixture["count_aromatic_rings"] == rdkit_result["aromatic_ring_count"]

    evidence = {
        "atom_mismatch_idxs": atom_mismatch,
        "bond_mismatch_pairs": bond_mismatch,
        "ring_count_match": ring_count_match,
        "aromatic_ring_count_match": aromatic_ring_count_match,
        "rdkit_aromatic_ring_count": rdkit_result["aromatic_ring_count"],
    }

    # --- kekulization direction (only meaningful if RDKit itself found
    # aromatic bonds to kekulize AND chematic's own kekulization succeeded) ---
    kekule_bond_mismatch = []
    kekule_checked = False
    kekule_valid_alternate = None
    if rdkit_result["had_aromatic_bonds"] and fixture["kekulize_ok"]:
        kekule_checked = True
        chematic_kekule_by_pair = {}
        for b in fixture["bonds"]:
            pair = frozenset((b["a1"], b["a2"]))
            chematic_kekule_by_pair[pair] = b["kekulized_order"]
            rdkit_kek = rdkit_result["kekule_bond_type_by_pair"].get(pair)
            if rdkit_kek is not None and rdkit_kek != b["kekulized_order"]:
                kekule_bond_mismatch.append(list(pair))
        if kekule_bond_mismatch:
            kekule_valid_alternate, detail = kekule_structure_is_valid_alternate(
                rdkit_result["_sanitized_mol"], chematic_kekule_by_pair
            )
            evidence["kekule_valid_alternate_detail"] = detail
    evidence["kekule_checked"] = kekule_checked
    evidence["kekule_bond_mismatch_pairs"] = kekule_bond_mismatch
    evidence["kekule_valid_alternate"] = kekule_valid_alternate

    model_says_aromatic = fixture["huckel_model_aromatic_atom_count"] > 0
    evidence["huckel_model_aromatic_atom_count"] = fixture["huckel_model_aromatic_atom_count"]
    evidence["model_says_aromatic"] = model_says_aromatic

    # --- indolizine gets one extra, additional check layered on top of
    # whatever bucket it lands in below: does the raw (pre-augmentation)
    # SSSR already return the correct [5, 6] decomposition, or does the
    # 9-ring fundamental-cycle artifact CLAUDE.md describes still reproduce? ---
    if fixture["id"] == "indolizine":
        artifact_reproduced = sorted(len(r) for r in fixture["raw_sssr_rings"]) != [5, 6]
        evidence["artifact_reproduced"] = artifact_reproduced
        evidence["raw_sssr_ring_sizes"] = sorted(len(r) for r in fixture["raw_sssr_rings"])
        if artifact_reproduced:
            return "unexpected_indolizine_sssr_artifact_reproduced", evidence

    both_nonaromatic = not any(rdkit_atom_arom) and not any(chem_atom_arom)
    # "Fully agree" is only the genuine, trustworthy bucket when the raw
    # Huckel model ITSELF independently confirms aromaticity -- not merely
    # when the final flags happen to match RDKit's. Without this,
    # tellurophene/phosphole (flags coincidentally match RDKit only because
    # kekulize() failed and left the original Aromatic bond order untouched)
    # would wrongly land here instead of in case 3 below.
    both_fully_aromatic_agree = (
        not atom_mismatch
        and not bond_mismatch
        and aromatic_ring_count_match
        and any(chem_atom_arom)
        and model_says_aromatic
    )

    # --- case 1: both engines fully agree, atom AND bond flags consistent,
    # AND the raw model genuinely confirms it (not just the final flags) ---
    if both_nonaromatic and not atom_mismatch and not bond_mismatch and aromatic_ring_count_match:
        bucket = "both_correctly_nonaromatic"
    elif both_fully_aromatic_agree:
        if not kekule_checked or not kekule_bond_mismatch:
            if fixture["category"] == "exocyclic_multiple_bond":
                bucket = "matches_rdkit_exact_kekule_exocyclic_bond_excluded_correctly"
            else:
                bucket = "matches_rdkit_exact_kekule"
        elif kekule_valid_alternate:
            bucket = "matches_rdkit_kekule_valid_alternate"
        else:
            return "unexpected_invalid_kekule_structure", evidence
    else:
        bucket = None

    if bucket is not None:
        # indolizine additionally requires the SSSR-bridge-artifact check
        # above to have found no artifact for this specific, otherwise-
        # unremarkable "matches RDKit" outcome to earn its own distinct,
        # more informative bucket name (see CLAUDE.md's 9-ring claim).
        if fixture["id"] == "indolizine" and bucket == "matches_rdkit_exact_kekule":
            return "sssr_bridge_artifact_not_reproduced_docs_stale", evidence
        return bucket, evidence

    # --- case 2: atoms agree (both True) but bonds disagree, and the raw
    # Huckel model itself does NOT confirm aromaticity while chematic's own
    # kekulization succeeded -- the model's "no" got silently overridden by
    # a stale pre-existing atom.aromatic=true flag that `apply_aromaticity`
    # never clears (see RFC finding on `build_molecule_from_model`'s
    # promote-only atom loop), while the bond loop correctly reflects the
    # real (non-aromatic) Kekule structure kekulize() computed -- the two
    # loops disagree with EACH OTHER, not just with RDKit. ---
    if (
        not atom_mismatch
        and bond_mismatch
        and not model_says_aromatic
        and fixture["kekulize_ok"]
    ):
        return "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent", evidence

    # --- case 3: atoms AND bonds both coincidentally match RDKit (both
    # True), the raw model does NOT confirm aromaticity, but this time
    # chematic's OWN kekulization failed outright -- so the "Aromatic" bond
    # order never got converted to Single/Double in the first place, and the
    # promote-only bug's fallback (`bond.order` when the model disagrees)
    # just echoes that already-Aromatic order back, matching RDKit's bond
    # flag by coincidence rather than by any successful computation. ---
    if (
        not atom_mismatch
        and not bond_mismatch
        and not model_says_aromatic
        and not fixture["kekulize_ok"]
    ):
        return "kekulize_fails_atom_bond_flags_survive_coincidentally", evidence

    return "unexpected_mismatch", evidence


def _self_test():
    """Positive controls for the fail-closed machinery itself, run before any
    RDKit call. A plain assertion failure here means the fail-closed logic
    itself is broken and must hard-crash immediately -- mirrors
    scripts/stereo2d_diagnosis.py's own four self-tests."""

    # Control A: an atom-count mismatch (chematic vs RDKit parsed a different
    # number of atoms for the "same" SMILES) must be caught explicitly, not
    # silently zipped/truncated.
    bogus_fixture = {
        "id": "self_test_bogus",
        "category": "heteroaromatic",
        "atoms": [{"default_aromatic": True}],
        "bonds": [],
        "kekulize_ok": True,
        "raw_sssr_ring_count": 0,
        "count_aromatic_rings": 0,
        "huckel_model_aromatic_atom_count": 0,
    }
    bucket, _ = classify(bogus_fixture, {"parsed": True, "atom_aromatic": [True, True], "bond_aromatic_by_pair": {}, "ring_atom_sets": [], "aromatic_ring_count": 0, "had_aromatic_bonds": False, "kekule_bond_type_by_pair": {}})
    assert bucket == "unexpected_atom_count_mismatch", f"self-test A: expected atom-count-mismatch bucket, got {bucket!r}"
    assert bucket not in EXPECTED_BUCKETS, "self-test A: an 'unexpected_*' bucket leaked into EXPECTED_BUCKETS"

    # Control B: a genuine atom-aromaticity disagreement (not explained by
    # any known mechanism) must fall through to "unexpected_mismatch", not
    # get silently absorbed into "both_correctly_nonaromatic" or similar.
    weak_fixture = {
        "id": "self_test_weak",
        "category": "heteroaromatic",
        "atoms": [
            {"default_aromatic": True, "element": "C"},
            {"default_aromatic": False, "element": "N"},
        ],
        "bonds": [{"a1": 0, "a2": 1, "default_aromatic": False, "kekulized_order": "Single"}],
        "kekulize_ok": True,
        "raw_sssr_ring_count": 0,
        "count_aromatic_rings": 0,
        "huckel_model_aromatic_atom_count": 2,
    }
    weak_rdkit = {
        "parsed": True,
        "atom_aromatic": [False, False],
        "atom_element": ["C", "N"],
        "bond_aromatic_by_pair": {frozenset((0, 1)): False},
        "ring_atom_sets": [],
        "aromatic_ring_count": 0,
        "had_aromatic_bonds": False,
        "kekule_bond_type_by_pair": {},
    }
    bucket, evidence = classify(weak_fixture, weak_rdkit)
    assert bucket == "unexpected_mismatch", f"self-test B: expected unexpected_mismatch, got {bucket!r}"
    assert evidence["atom_mismatch_idxs"] == [0], "self-test B: evidence did not actually record the mismatching atom index"
    assert bucket not in EXPECTED_BUCKETS, "self-test B: an 'unexpected_*' bucket leaked into EXPECTED_BUCKETS"

    # Control B2: an atom-order/element mismatch (chematic's atom i is NOT
    # RDKit's atom i, even though both parsed the same atom count) must be
    # caught explicitly and BEFORE any aromaticity comparison runs on the
    # wrongly-aligned indices -- this is the guard the advisor flagged as
    # missing, since every other check in this script silently assumes
    # index alignment holds.
    misaligned_rdkit = dict(weak_rdkit, atom_element=["N", "C"])  # swapped vs weak_fixture's C, N
    bucket_b2, evidence_b2 = classify(weak_fixture, misaligned_rdkit)
    assert bucket_b2 == "unexpected_atom_order_or_element_mismatch", (
        f"self-test B2: expected the order/element-mismatch bucket, got {bucket_b2!r}"
    )
    assert evidence_b2["element_mismatch_idxs"] == [0, 1], "self-test B2: did not record which indices disagree"
    assert bucket_b2 not in EXPECTED_BUCKETS, "self-test B2: an 'unexpected_*' bucket leaked into EXPECTED_BUCKETS"

    # Control C: duplicate fixture IDs must be detected by the same check
    # main() runs on the real dump.
    dup_ids = ["a", "b", "a"]
    assert len(set(dup_ids)) != len(dup_ids), "self-test C: duplicate-ID fixture itself has no duplicates"

    # Control D: "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent"'s
    # whole evidentiary claim is (a) atom flags coincidentally match RDKit
    # while bond flags do NOT, AND (b) the raw Huckel model independently
    # says non-aromatic. Prove classify() actually checks the model count
    # (field `huckel_model_aromatic_atom_count`), not just the atom/bond
    # mismatch shape -- an otherwise-identical case where the model DID
    # confirm aromaticity (>0) must NOT land in this bucket, since bonds
    # still disagreeing with RDKit despite the model agreeing would be a
    # genuinely different, unexplained problem, not this known one.
    inconsistent_shape = {
        "id": "self_test_inconsistent",
        "category": "fused_polycyclic",
        "atoms": [{"default_aromatic": True, "element": "C"}] * 5,
        "bonds": [
            {"a1": i, "a2": (i + 1) % 5, "default_aromatic": False, "kekulized_order": "Single"}
            for i in range(5)
        ],
        "kekulize_ok": True,
        "raw_sssr_ring_count": 1,
        "count_aromatic_rings": 1,
        "huckel_model_aromatic_atom_count": 0,
    }
    inconsistent_rdkit = {
        "parsed": True,
        "atom_aromatic": [True] * 5,
        "atom_element": ["C"] * 5,
        "bond_aromatic_by_pair": {frozenset((i, (i + 1) % 5)): True for i in range(5)},
        "ring_atom_sets": [frozenset(range(5))],
        "aromatic_ring_count": 1,
        "had_aromatic_bonds": False,
        "kekule_bond_type_by_pair": {},
    }
    bucket, evidence = classify(inconsistent_shape, inconsistent_rdkit)
    assert bucket == "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent", (
        f"self-test D1: expected the inconsistent-flags bucket, got {bucket!r}"
    )
    assert evidence["model_says_aromatic"] is False, "self-test D1: evidence did not record the model verdict"

    model_agrees_shape = dict(inconsistent_shape, huckel_model_aromatic_atom_count=5)
    bucket2, evidence2 = classify(model_agrees_shape, inconsistent_rdkit)
    assert bucket2 != "kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent", (
        "self-test D2: classify() ignored the model count -- landed in the known bucket even though the model agreed"
    )
    assert evidence2["model_says_aromatic"] is True, "self-test D2: evidence did not record the model-agrees case"


def main():
    _self_test()

    if not DUMP_PATH.exists():
        print(
            f"missing {DUMP_PATH} -- run: cargo run -p chematic-perception "
            f"--example aromaticity_rdkit_parity_dump > {DUMP_PATH}",
            file=sys.stderr,
        )
        sys.exit(1)

    fixtures = [json.loads(line) for line in DUMP_PATH.read_text().splitlines() if line.strip()]

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
                f"(update EXPECTED_BUCKET_BY_ID and classify() in this script): {sorted(extra)}",
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
        rdkit_result = rdkit_read(fx["smiles"])
        bucket, evidence = classify(fx, rdkit_result)
        bucket_counts[bucket] = bucket_counts.get(bucket, 0) + 1
        expected_bucket = EXPECTED_BUCKET_BY_ID[fx["id"]]
        if bucket != expected_bucket:
            baseline_drift.append((fx["id"], expected_bucket, bucket))
        rows.append(
            {
                "id": fx["id"],
                "category": fx["category"],
                "description": fx["description"],
                "smiles": fx["smiles"],
                "chematic": {
                    k: fx.get(k)
                    for k in (
                        "kekulize_ok",
                        "kekulize_error",
                        "experimental_ok",
                        "experimental_error",
                        "raw_sssr_ring_count",
                        "augmented_ring_count",
                        "count_aromatic_rings",
                        "huckel_model_aromatic_atom_count",
                        "rdkitlike_model_aromatic_atom_count",
                    )
                },
                "rdkit": {
                    "parsed": rdkit_result.get("parsed"),
                    "aromatic_ring_count": rdkit_result.get("aromatic_ring_count"),
                    "ring_count": len(rdkit_result.get("ring_atom_sets", [])),
                    "canonical_smiles": rdkit_result.get("canonical_smiles"),
                },
                "failure_bucket": bucket,
                "expected_bucket": expected_bucket,
                "evidence": evidence,
            }
        )

    unexplained = [r for r in rows if r["failure_bucket"] not in EXPECTED_BUCKETS]
    fail = bool(unexplained or baseline_drift)

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
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, default=str))

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
