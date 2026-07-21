#!/usr/bin/env python3
"""Canonical SMILES residual diagnosis: chematic vs RDKit (diag/canonical-smiles-residual).

Builds on the existing `scripts/canonical_diff.py` (round-trip + idempotency,
committed results in `validation/results/canonical_diff.jsonl`) and
`scripts/canonical_structural_correctness.py` (per-molecule structural
correctness under RDKit-generated variants). This script adds the four checks
requested for this diagnosis round, measured and reported SEPARATELY (never
pooled into one "agreement %"):

  1. RDKit exact canonical string parity  -- chematic's canonical SMILES
     byte-for-byte equal to RDKit's, same input. NOT expected to be 100% even
     for structurally-correct output (different canonicalization algorithms) --
     a mismatch here is only an entry point into check 4, not proof of a bug.
  2. Permutation invariance -- K reproducibly (seeded) randomly-relabeled
     RDKit spellings of the SAME parsed molecule (`RenumberAtoms` + a fixed
     `random.Random(seed)`, not RDKit's own unseeded `doRandom=True` --
     picked deliberately so this script's output is byte-reproducible run
     to run), fed through chematic; chematic's own output must be identical
     across all K+1 spellings. A failure here IS always a real chematic bug
     (chematic-internal self-consistency, RDKit is only used to generate
     alternate valid spellings of one molecule). NOTE: this is a lower
     bound at the tested K -- passing K relabelings is not a proof of
     invariance under all possible relabelings, only evidence against it
     not being found in K samples.
  3. Idempotence -- canonical(canonical(s)) == canonical(s).
  4. Semantic structure parity -- for every check-1 mismatch, reparse
     chematic's canonical output through RDKit and compare ACTUAL STRUCTURE
     (formula / heavy-atom multiset / aromatic atom+bond counts / ring-size
     multiset / CIP stereocenter labels / bond E-Z multiset / isotope+charge+
     atom-map multiset) against RDKit's own canonicalization of the original
     input -- never string comparison alone.

Every check-1/check-4 residual lands in exactly one named bucket (see
BUCKETS below) or "unclassified" -- no silent drops.

Usage:
    python scripts/canonical_residual_diagnosis.py [SMILES.csv] [--limit N] [--k N]
    python scripts/canonical_residual_diagnosis.py --self-test

Writes:
    validation/results/canonical_residual_diagnosis.jsonl   (per-mismatch rows)
    validation/results/canonical_residual_diagnosis_summary.json (headline numbers)
"""
import argparse
import itertools
import json
import os
import sys
from dataclasses import dataclass, asdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_JSONL = os.path.join(ROOT, "validation", "results", "canonical_residual_diagnosis.jsonl")
OUT_SUMMARY = os.path.join(ROOT, "validation", "results", "canonical_residual_diagnosis_summary.json")

BUCKETS = [
    "aromaticity_kekulization",
    "tetrahedral_parity",
    "ez_direction",
    "ring_closure_ordering",
    "bridged_fused_spiro",
    "disconnected_fragment_ordering",
    "isotope_charge_atommap",
    "symmetry_tie_breaking",
    "writer_token_bug",
    "unclassified",
]


# ---------------------------------------------------------------------------
# Pure-data feature extraction + classification (RDKit-object-free below the
# extraction boundary, so the classification decision tree is directly
# unit-testable without constructing tricky RDKit Mol objects -- see
# --self-test).
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MolFeatures:
    formula: str
    heavy_atoms: tuple      # sorted atomic numbers
    aromatic_counts: tuple  # (n_aromatic_atoms, n_aromatic_bonds)
    ring_sizes: tuple       # sorted SSSR ring sizes
    stereocenters: tuple    # sorted CIP labels (incl. "?" for unassigned)
    ez_bonds: tuple         # sorted bond-stereo labels (E/Z only)
    isotope_charge_map: tuple  # sorted (isotope, charge, atom-map) triples


def extract_features(mol) -> MolFeatures:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    heavy = tuple(sorted(a.GetAtomicNum() for a in mol.GetAtoms()))
    arom_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    arom_bonds = sum(1 for b in mol.GetBonds() if b.GetIsAromatic())
    rings = tuple(sorted(len(r) for r in mol.GetRingInfo().AtomRings()))
    centers = tuple(sorted(
        lbl for _, lbl in Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, useLegacy=False)
    ))
    ez = tuple(sorted(
        str(b.GetStereo()) for b in mol.GetBonds()
        if b.GetStereo() not in (Chem.BondStereo.STEREONONE,)
    ))
    icm = tuple(sorted(
        (a.GetIsotope(), a.GetFormalCharge(), a.GetAtomMapNum())
        for a in mol.GetAtoms()
    ))
    return MolFeatures(
        formula=rdMolDescriptors.CalcMolFormula(mol),
        heavy_atoms=heavy,
        aromatic_counts=(arom_atoms, arom_bonds),
        ring_sizes=rings,
        stereocenters=centers,
        ez_bonds=ez,
        isotope_charge_map=icm,
    )


def classify_real_diff(a: MolFeatures, b: MolFeatures):
    """a = features of RDKit's own parse of the ORIGINAL input.
    b = features of RDKit's parse of chematic's canonical output.
    Precondition: their RDKit-recanonicalized strings differ (a genuine
    semantic residual, not just a spelling difference)."""
    if a.heavy_atoms != b.heavy_atoms:
        return "writer_token_bug", f"heavy-atom multiset differs: {a.heavy_atoms} vs {b.heavy_atoms}"
    if a.formula != b.formula:
        return "aromaticity_kekulization", f"formula differs (H count / valence): {a.formula} vs {b.formula}"
    if a.aromatic_counts != b.aromatic_counts:
        return "aromaticity_kekulization", f"aromatic (atom,bond) counts differ: {a.aromatic_counts} vs {b.aromatic_counts}"
    if a.ring_sizes != b.ring_sizes:
        return "bridged_fused_spiro", f"SSSR ring-size multiset differs: {a.ring_sizes} vs {b.ring_sizes}"
    if a.stereocenters != b.stereocenters:
        return "tetrahedral_parity", f"CIP stereocenter label multiset differs: {a.stereocenters} vs {b.stereocenters}"
    if a.ez_bonds != b.ez_bonds:
        return "ez_direction", f"bond E/Z multiset differs: {a.ez_bonds} vs {b.ez_bonds}"
    if a.isotope_charge_map != b.isotope_charge_map:
        return "isotope_charge_atommap", "isotope/charge/atom-map multiset differs"
    return "unclassified", "no structural probe found a difference, yet RDKit re-canonicalization strings differ"


def classify_cosmetic(is_multi_fragment: bool, has_symmetry_tie: bool, ring_relationship: str):
    """Classify a check-1 mismatch that check-4 found to be semantically
    IDENTICAL (RDKit re-canonicalizes both sides to the same string) -- i.e.
    a real algorithm-choice difference, not a chematic bug.
    ring_relationship in {"none", "simple", "fused", "spiro", "bridged"}.

    NOTE on priority order (most-specific-structural-fact first, symmetry
    last): almost every check-1 mismatch in this "semantically identical"
    branch is, at bottom, just "two different valid canonicalization
    algorithms disagreeing" -- not a bug. This sub-classification is a
    best-effort DESCRIPTIVE heuristic over the molecule's topology (does it
    have a ring / a bridged-fused-spiro system / multiple fragments / a
    tied canonical rank somewhere), reported to characterize what KINDS of
    molecules dominate the cosmetic-mismatch population -- it is NOT a
    proof that the specific string divergence was CAUSED by that trait.
    Ring/fragment topology is checked before the symmetry-tie catch-all
    because it is closer to a literal, well-precedented mechanism (see
    `bridged_bicyclic_canonical_gap_documentation` in
    crates/chematic-smiles/tests/canonical_robustness.rs); tied canonical
    ranks are near-ubiquitous in drug-like molecules (any monosubstituted
    phenyl ring has one) and are checked last so they don't swallow the
    more specific ring-shape buckets."""
    if is_multi_fragment:
        return "disconnected_fragment_ordering"
    if ring_relationship in ("spiro", "bridged"):
        return "bridged_fused_spiro"
    if ring_relationship in ("simple", "fused"):
        return "ring_closure_ordering"
    if has_symmetry_tie:
        return "symmetry_tie_breaking"
    return "unclassified"


def ring_relationship(mol) -> str:
    atom_rings = mol.GetRingInfo().AtomRings()
    if not atom_rings:
        return "none"
    if len(atom_rings) == 1:
        return "simple"
    bonds = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds()}
    worst = "fused"
    for ra, rb in itertools.combinations(atom_rings, 2):
        shared = set(ra) & set(rb)
        if len(shared) == 1:
            return "spiro"
        if len(shared) >= 2:
            pair_is_edge = any(frozenset(p) in bonds for p in itertools.combinations(shared, 2))
            if not pair_is_edge or len(shared) > 2:
                return "bridged"
    return worst


def has_symmetry_tie(mol) -> bool:
    from rdkit import Chem
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    return len(set(ranks)) < mol.GetNumAtoms()


# ---------------------------------------------------------------------------
# Main diagnosis loop
# ---------------------------------------------------------------------------

def random_relabeling(rm, rng):
    """A reproducible (seeded) random atom-relabeling of `rm` as a valid,
    non-canonical SMILES string of the SAME molecule. Deliberately not
    RDKit's `doRandom=True` (which draws from RDKit's own unseeded global
    RNG state -- every run picks different variants, so a committed JSONL/
    JSON result and a re-run of this script would disagree even with no
    code change). `RenumberAtoms` + `canonical=False` gives a real
    randomly-relabeled spelling (RDKit just serializes the atoms in the
    given order, no re-canonicalization), fully determined by `rng`."""
    from rdkit import Chem
    order = list(range(rm.GetNumAtoms()))
    rng.shuffle(order)
    return Chem.MolToSmiles(Chem.RenumberAtoms(rm, order), canonical=False)


def run(corpus_path, limit, k_variants, out_jsonl=OUT_JSONL, out_summary=OUT_SUMMARY, seed=0):
    import random
    import chematic
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    rng = random.Random(seed)

    smis = [l.strip() for l in open(os.path.expanduser(corpus_path)) if l.strip()]
    if limit:
        smis = smis[:limit]

    n_total = 0
    n_chematic_parse_fail = 0

    exact_match = 0
    exact_mismatch = 0

    idem_ok = 0
    idem_fail = 0

    perm_invariant = 0
    perm_fail = 0
    perm_fail_examples = []

    bucket_counts = {b: 0 for b in BUCKETS}
    bucket_examples = {b: [] for b in BUCKETS}
    rows = []

    for s in smis:
        rm = Chem.MolFromSmiles(s)
        if rm is None:
            continue
        try:
            cm = chematic.from_smiles(s).smiles
        except Exception:
            n_chematic_parse_fail += 1
            continue
        n_total += 1

        # --- Check 3: idempotence ---
        try:
            cm2 = chematic.from_smiles(cm).smiles
        except Exception:
            cm2 = None
        if cm2 == cm:
            idem_ok += 1
        else:
            idem_fail += 1

        # --- Check 2: permutation invariance (chematic self-consistency) ---
        outputs = {cm}
        perm_error = False
        for _ in range(k_variants):
            variant = random_relabeling(rm, rng)
            try:
                outputs.add(chematic.from_smiles(variant).smiles)
            except Exception:
                perm_error = True
        if len(outputs) == 1 and not perm_error:
            perm_invariant += 1
        else:
            perm_fail += 1
            if len(perm_fail_examples) < 10000:  # effectively unbounded; capped only against pathological input
                # Distinguish "writer emits >1 spelling for the SAME
                # molecule" (still a real canonicalization-stability bug,
                # but chemistry preserved) from "the outputs are actually
                # DIFFERENT molecules" (a much more serious correctness bug)
                # -- classified via RDKit re-canonicalization of each output,
                # never by string comparison of the outputs themselves.
                canon_forms = set()
                for o in outputs:
                    om = Chem.MolFromSmiles(o)
                    canon_forms.add(Chem.MolToSmiles(om) if om is not None else None)
                perm_fail_examples.append({
                    "smiles": s,
                    "distinct_outputs": sorted(outputs),
                    "chematic_variant_parse_error": perm_error,
                    "outputs_semantically_identical": len(canon_forms) == 1,
                    "has_ez_marker": any("/" in o or "\\" in o for o in outputs),
                })

        # --- Check 1: exact RDKit string parity ---
        rd_native = Chem.MolToSmiles(rm)
        if cm == rd_native:
            exact_match += 1
            continue
        exact_mismatch += 1

        # --- Check 4: semantic structure parity for this check-1 mismatch ---
        mol_cm = Chem.MolFromSmiles(cm)
        if mol_cm is None:
            bucket = "writer_token_bug"
            detail = "chematic's canonical output is not RDKit-parseable"
            bucket_counts[bucket] += 1
            if len(bucket_examples[bucket]) < 10:
                bucket_examples[bucket].append({"smiles": s, "chematic": cm, "detail": detail})
            rows.append({"smiles": s, "chematic": cm, "rdkit_native": rd_native,
                         "bucket": bucket, "detail": detail, "semantically_identical": False})
            continue

        rd_of_cm = Chem.MolToSmiles(mol_cm)
        if rd_of_cm == rd_native:
            bucket = classify_cosmetic(
                is_multi_fragment=len(Chem.GetMolFrags(rm)) > 1,
                has_symmetry_tie=has_symmetry_tie(rm),
                ring_relationship=ring_relationship(rm),
            )
            detail = "structurally identical molecule; different valid canonicalization spelling"
            sem_identical = True
        else:
            feat_a = extract_features(rm)
            feat_b = extract_features(mol_cm)
            bucket, detail = classify_real_diff(feat_a, feat_b)
            sem_identical = False

        bucket_counts[bucket] += 1
        if len(bucket_examples[bucket]) < 10:
            bucket_examples[bucket].append({"smiles": s, "chematic": cm, "detail": detail})
        rows.append({"smiles": s, "chematic": cm, "rdkit_native": rd_native,
                     "rdkit_of_chematic": rd_of_cm, "bucket": bucket, "detail": detail,
                     "semantically_identical": sem_identical})

    os.makedirs(os.path.dirname(out_jsonl), exist_ok=True)
    with open(out_jsonl, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")

    summary = {
        "corpus": corpus_path,
        "n_input_lines": len(smis),
        "n_total_parsed_by_both": n_total,
        "chematic_parse_failures": n_chematic_parse_fail,
        "k_variants_per_molecule": k_variants,
        "permutation_seed": seed,
        "check1_exact_string_parity": {
            "match": exact_match, "mismatch": exact_mismatch, "n": n_total,
            "pct_match": round(100 * exact_match / n_total, 2) if n_total else None,
        },
        "check2_permutation_invariance": {
            "invariant": perm_invariant, "fail": perm_fail, "n": n_total,
            "pct_invariant": round(100 * perm_invariant / n_total, 2) if n_total else None,
            "fail_semantically_identical_outputs": sum(
                1 for e in perm_fail_examples if e["outputs_semantically_identical"]),
            "fail_semantically_DIFFERENT_outputs": sum(
                1 for e in perm_fail_examples if not e["outputs_semantically_identical"]),
            "fail_examples_captured": len(perm_fail_examples),
        },
        "check3_idempotence": {
            "ok": idem_ok, "fail": idem_fail, "n": n_total,
            "pct_ok": round(100 * idem_ok / n_total, 2) if n_total else None,
        },
        "check4_bucket_breakdown_of_check1_mismatches": bucket_counts,
        "check4_bucket_totals_check": {
            "sum_of_buckets": sum(bucket_counts.values()),
            "check1_mismatch_count": exact_mismatch,
            "accounted_for": sum(bucket_counts.values()) == exact_mismatch,
        },
    }
    with open(out_summary, "w") as f:
        json.dump({"summary": summary, "bucket_examples": bucket_examples,
                    "permutation_invariance_failures_sample": perm_fail_examples}, f, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"\nwrote {len(rows)} check-1-mismatch rows -> {os.path.relpath(out_jsonl, ROOT)}")
    print(f"wrote summary -> {os.path.relpath(out_summary, ROOT)}")
    print(f"\npermutation-invariance FAILURES (real chematic bugs): {perm_fail}/{n_total}")
    for ex in perm_fail_examples[:5]:
        print(f"  {ex['smiles']!r} -> {len(ex['distinct_outputs'])} distinct outputs")
    return summary


# ---------------------------------------------------------------------------
# Self-test: verifies the classification decision tree actually discriminates
# every bucket (fail-closed check, not just "ran without crashing") --
# per this repo's measurement-harness-controls convention (positive AND
# negative controls, not just positive).
# ---------------------------------------------------------------------------

def self_test():
    base = MolFeatures(
        formula="C6H6", heavy_atoms=(6,) * 6, aromatic_counts=(6, 6),
        ring_sizes=(6,), stereocenters=(), ez_bonds=(), isotope_charge_map=((0, 0, 0),) * 6,
    )

    def with_(**kw):
        d = asdict(base)
        d.update(kw)
        return MolFeatures(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in d.items()})

    checks = []

    b, _ = classify_real_diff(base, with_(heavy_atoms=(6,) * 5 + (7,)))
    checks.append(("writer_token_bug", b))
    assert b == "writer_token_bug"

    b, _ = classify_real_diff(base, with_(formula="C6H8"))
    checks.append(("aromaticity_kekulization (formula)", b))
    assert b == "aromaticity_kekulization"

    b, _ = classify_real_diff(base, with_(aromatic_counts=(0, 0)))
    checks.append(("aromaticity_kekulization (arom counts)", b))
    assert b == "aromaticity_kekulization"

    b, _ = classify_real_diff(base, with_(ring_sizes=(3, 3)))
    checks.append(("bridged_fused_spiro (ring sizes)", b))
    assert b == "bridged_fused_spiro"

    b, _ = classify_real_diff(base, with_(stereocenters=("R",)))
    checks.append(("tetrahedral_parity", b))
    assert b == "tetrahedral_parity"

    b, _ = classify_real_diff(base, with_(ez_bonds=("STEREOE",)))
    checks.append(("ez_direction", b))
    assert b == "ez_direction"

    b, _ = classify_real_diff(base, with_(isotope_charge_map=((0, 0, 1),) * 6))
    checks.append(("isotope_charge_atommap", b))
    assert b == "isotope_charge_atommap"

    b, _ = classify_real_diff(base, base)
    checks.append(("unclassified (identical features)", b))
    assert b == "unclassified"

    assert classify_cosmetic(True, False, "none") == "disconnected_fragment_ordering"
    assert classify_cosmetic(False, True, "none") == "symmetry_tie_breaking"
    assert classify_cosmetic(False, False, "spiro") == "bridged_fused_spiro"
    assert classify_cosmetic(False, False, "bridged") == "bridged_fused_spiro"
    assert classify_cosmetic(False, False, "fused") == "ring_closure_ordering"
    assert classify_cosmetic(False, False, "simple") == "ring_closure_ordering"
    assert classify_cosmetic(False, False, "none") == "unclassified"
    # priority: fragment ordering and ring-shape facts win over the symmetry
    # catch-all (see classify_cosmetic's docstring for why)
    assert classify_cosmetic(True, True, "bridged") == "disconnected_fragment_ordering"
    assert classify_cosmetic(False, True, "bridged") == "bridged_fused_spiro"

    # ring_relationship on real RDKit mols -- spiro/fused/bridged smoke test
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    assert ring_relationship(Chem.MolFromSmiles("C1CCCCC1")) == "simple"
    assert ring_relationship(Chem.MolFromSmiles("c1ccc2ccccc2c1")) == "fused"
    assert ring_relationship(Chem.MolFromSmiles("C1CCC2(CC1)CCCC2")) == "spiro"
    assert ring_relationship(Chem.MolFromSmiles("C1CC2CCC1CC2")) == "bridged"  # bicyclo[2.2.2]octane

    # Positive control: a KNOWN real permutation-invariance bug must be
    # detected end-to-end (not just at the pure-function level). This is the
    # bridged-bicyclic ring-closure-ordering gap already documented in
    # crates/chematic-smiles/tests/canonical_robustness.rs
    # (`bridged_bicyclic_canonical_gap_documentation`) -- verifies the harness
    # is fail-closed, not just fail-silent-passing.
    import chematic
    a = chematic.from_smiles("C1CC2CCC1CC2").smiles
    b2 = chematic.from_smiles("C1CCC2CC1CC2").smiles
    assert a != b2, (
        "expected the known bicyclo[2.2.2]octane permutation-invariance bug "
        "to still reproduce (positive control) -- if this now passes, the "
        "bug may have been fixed upstream; update the RFC and this control"
    )

    print("self-test OK:")
    for label, bucket in checks:
        print(f"  {label:45s} -> {bucket}")
    print("  ring_relationship simple/fused/spiro/bridged smoke test -> OK")
    print("  positive control: bicyclo[2.2.2]octane permutation-invariance bug -> reproduces (expected)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("corpus", nargs="?", default="~/Downloads/SMILES.csv")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--k", type=int, default=8, help="permutation variants per molecule (check 2)")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed for reproducible permutation relabeling")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out-jsonl", default=OUT_JSONL)
    ap.add_argument("--out-summary", default=OUT_SUMMARY)
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    run(args.corpus, args.limit, args.k, args.out_jsonl, args.out_summary, args.seed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
