#!/usr/bin/env python3
"""
ECFP4 vs RDKit agreement -- the "Round 1" migration-decision metric, never
measured before this script existed.

Headline result (5,000-mol ChEMBL corpus): chematic's ECFP4 is NOT the
standard Rogers-Hahn/RDKit ECFP definition -- its radius-0 invariant
includes atom.aromatic, which RDKit's default invariant set omits
(source-confirmed, crates/chematic-fp/src/ecfp.rs:initial_atom_id). This is
a legitimate design choice (cf. FCFP), not a hash-folding artifact, and it
is why only ~77% of molecules produce an identical invariant-equivalence
structure to RDKit's (tier 2) and pairwise similarity correlates at r=0.94,
not ~1.0 (tier 3). Practical consequence: bit vectors are not
RDKit-compatible (expected -- different hash), and similarity is *close but
not identical* to RDKit's -- RDKit-trained thresholds/models should not be
assumed to transfer without re-validation.

Separately, tier 5 found a real, previously-unknown self-consistency defect:
because atom.aromatic is not auto-perceived for Kekule-written SMILES
(chematic requires an explicit apply_aromaticity() call), ecfp4() gives a
DIFFERENT fingerprint for two equally valid spellings of the identical
molecule in 92% of cases if the caller skips that call (a footgun, but
arguably working as documented -- see CLAUDE.md/README on the
aromaticity-perception contract). More seriously, ~13% of molecules STILL
mismatch even *after* calling apply_aromaticity() as documented -- of
those, roughly a third (41/130) also disagree on the full aromatic-atom/bond
assignment multiset (same class of bug as the known aromatic_context
perception issue, literal-same-code-path unverified) and two-thirds (89/130)
have an IDENTICAL assignment multiset yet still a different fingerprint --
an unattributed, real defect.

Tier 6 is the follow-up to that residual: is it ECFP4-specific, and is it
really about aromaticity? Answer to both: no, and not entirely. The 130
residual molecules are the SAME ones canonical_smiles() and InChI diverge on
(100% overlap each way) -- one shared, systemic root cause across at least
three core output functions, not an ECFP4 quirk. SSSR ring-size
decomposition is identical for all 130 (rules out ring-finding as the
cause). And a discriminating control that varies ONLY atom traversal order
(never touching Kekule/aromatic origin) shows ECFP4 barely moves (0.0%),
ruling out plain order-dependence as ECFP4's explanation and pointing back
at the Kekule-vs-aromatic origin itself as the operative variable --
though canonical_smiles DOES carry its own separate, smaller
order-sensitivity (3.8%) on the same control, worth tracking independently.
The exact mechanism (what specifically differs between the two origins,
given flags/rings/order are all ruled out) is NOT yet identified -- would
require reading apply_aromaticity()'s and canonical.rs's bond-order handling
directly, out of scope for this measurement-only script.

Raw bit-vector equality is NOT a meaningful cross-implementation metric here
and is intentionally NOT the headline number: chematic hashes atom
environments with FNV-1a (see crates/chematic-fp/src/ecfp.rs, whose own
doc comment says bit positions are not meant to match RDKit's). Two
independent hash functions landing on the same bit index for the same
chemistry is a ~1/2048 coincidence per environment, so bit equality would
report a manufactured near-random number regardless of whether the
underlying chemistry is correct -- and because both fingerprints are sparse,
raw per-position agreement is actually dominated by 0/0 non-matches and can
look misleadingly HIGH, not low. It is printed once (tier 0) purely so the
number exists on the record and nobody re-derives the same false alarm in a
future round.

Tiers, all hash-independent, all built on introspection that already existed
on both sides before this script (chematic's
Mol.ecfp_bitinfo/morgan_fp_counts/bond_table, RDKit's bitInfo /
GetSparseCountFingerprint):

  1. Coverage parity   -- does chematic generate an environment at every
                           (atom, radius) RDKit does, full corpus. Tests
                           *emission slots*, not what each slot encodes.
  2. Invariant partition agreement -- the real chemistry check, and the one
                           that actually found a difference. Within each
                           implementation, group that molecule's raw
                           (unfolded) environment hashes by equality -- two
                           environments landing on the same raw hash means
                           THAT implementation considers them chemically
                           identical. Compare the resulting partition SHAPE
                           (sorted multiset of group sizes) between chematic
                           and RDKit -- hash-*value*-independent, so this
                           isolates genuine invariant-encoding disagreement.
                           Result: chematic's radius-0 invariant includes
                           `atom.aromatic` (crates/chematic-fp/src/ecfp.rs:
                           initial_atom_id); RDKit's default invariant does
                           not. Root-caused, not assumed -- see tier 5 for
                           the practical consequence this has.
  3. Similarity-structure preservation -- for a sample of molecule pairs,
                           does chematic's Tanimoto(A,B) correlate with
                           RDKit's Tanimoto(A,B)? This is the practical
                           "is it a valid drop-in for similarity search /
                           clustering / QSAR" answer. The residual gap here
                           is explained by tier 2's invariant difference,
                           not purely by hash-fold collisions.
  4. Connectivity sanity check (auxiliary, NOT a fingerprint test) -- does
                           independently-run BFS adjacency agree between the
                           two libraries' parsed graphs? This only checks
                           that both parsers agree on which atoms are bonded
                           to which; it never touches fingerprint code and
                           cannot by itself detect an invariant-encoding
                           difference (that's tier 2's job). Kept because a
                           parser-level disagreement would invalidate every
                           other tier's atom-index correspondence.
  5. Aromaticity representation-dependence -- the practical consequence of
                           tier 2's finding: because `atom.aromatic` feeds
                           the invariant and is NOT auto-perceived for
                           Kekule-written SMILES (chematic requires an
                           explicit apply_aromaticity() call), ecfp4() can
                           give a DIFFERENT fingerprint for two equally
                           valid spellings of the identical molecule. This
                           is chematic-internal self-consistency, not an
                           RDKit comparison, and confirms whether
                           apply_aromaticity() is a complete fix at scale.
  6. Layer 2 shared-mechanism check -- is tier 5's post-apply_aromaticity()
                           residual ECFP4-specific, or does it hit
                           canonical_smiles()/InChI on the SAME molecules?
                           Cross-consumer overlap, SSSR ring-size multiset
                           (rules out ring decomposition), and an
                           order-only control that isolates plain
                           atom-order-dependence from Kekule-vs-aromatic
                           origin -- a confound present in tier 5's own
                           naive/residual comparison (RDKit's non-canonical
                           Kekule respelling changes atom order too, not
                           just origin).

Usage:
    .venv/bin/python scripts/ecfp4_agreement.py [SMILES.csv] [--limit N]
        [--partition-sample N] [--pairs-sample N] [--json out.json]
"""
import argparse
import json
import os
import random
import statistics
from collections import Counter
import sys


def bfs_ball(adj, start, radius):
    seen = {start}
    frontier = {start}
    for _ in range(radius):
        nxt = set()
        for u in frontier:
            for v in adj.get(u, ()):
                if v not in seen:
                    nxt.add(v)
        seen |= nxt
        frontier = nxt
        if not frontier:
            break
    return seen


def chematic_adjacency(mol):
    adj = {}
    for a1, a2, _btype, _arom in mol.bond_table:
        adj.setdefault(a1, set()).add(a2)
        adj.setdefault(a2, set()).add(a1)
    return adj


def rdkit_adjacency(rd):
    adj = {}
    for bond in rd.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adj.setdefault(a1, set()).add(a2)
        adj.setdefault(a2, set()).add(a1)
    return adj


def tier1_coverage_parity(smis, chematic, Chem, rdFingerprintGenerator, limit):
    # RDKit's default GetMorganFingerprintAsBitVect / GetMorganGenerator prune
    # "redundant" environments (includeRedundantEnvironments=False) -- an
    # atom's environment stops growing once RDKit's internal duplicate-
    # detection decides a larger radius wouldn't add new discriminating
    # structure. This is a real RDKit optimization, not a hash-folding
    # artifact (confirmed via the unfolded GetMorganFingerprint too) -- but
    # it means comparing chematic's (complete) per-atom-radius coverage
    # against RDKit's *default*-pruned coverage measures "did we replicate
    # RDKit's pruning heuristic," not "did we implement the same chemistry."
    # includeRedundantEnvironments=True disables the pruning so this tier
    # compares the actual chemistry, not RDKit's default trimming.
    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeRedundantEnvironments=True
    )
    n_mol = 0
    n_match = 0
    n_mismatch = 0
    examples = []
    for smi in smis[:limit] if limit else smis:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        try:
            cm = chematic.from_smiles(smi)
        except Exception:
            continue
        n_mol += 1

        _fp_c, info_c = cm.ecfp_bitinfo(2)
        chem_pairs = {p for lst in info_c.values() for p in lst}

        ao = rdFingerprintGenerator.AdditionalOutput()
        ao.AllocateBitInfoMap()
        gen.GetFingerprint(rd, additionalOutput=ao)
        rd_pairs = {p for lst in ao.GetBitInfoMap().values() for p in lst}

        if chem_pairs == rd_pairs:
            n_match += 1
        else:
            n_mismatch += 1
            if len(examples) < 10:
                examples.append(
                    {
                        "smiles": smi,
                        "chematic_only": sorted(chem_pairs - rd_pairs)[:5],
                        "rdkit_only": sorted(rd_pairs - chem_pairs)[:5],
                    }
                )
    return {
        "n_molecules": n_mol,
        "exact_coverage_match": n_match,
        "coverage_mismatch": n_mismatch,
        "agreement_pct": round(100.0 * n_match / n_mol, 2) if n_mol else None,
        "examples": examples,
    }


def tier2_invariant_partition_agreement(smis, chematic, Chem, rdFingerprintGenerator, limit):
    # Same fairness fix as tier 1: RDKit's default prunes "redundant"
    # environments, which would make its partition trivially different from
    # chematic's for reasons unrelated to invariant correctness.
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, includeRedundantEnvironments=True)

    sample = smis[:limit] if limit else smis

    n_mol = 0
    n_exact_profile_match = 0
    examples = []
    for smi in sample:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        try:
            cm = chematic.from_smiles(smi)
        except Exception:
            continue
        n_mol += 1

        chem_profile = sorted(cm.morgan_fp_counts(2).values())

        rd_sparse = gen.GetSparseCountFingerprint(rd)
        rd_profile = sorted(rd_sparse.GetNonzeroElements().values())

        if chem_profile == rd_profile:
            n_exact_profile_match += 1
        elif len(examples) < 10:
            examples.append(
                {
                    "smiles": smi,
                    "chematic_n_groups": len(chem_profile),
                    "chematic_total_envs": sum(chem_profile),
                    "rdkit_n_groups": len(rd_profile),
                    "rdkit_total_envs": sum(rd_profile),
                }
            )
    return {
        "n_molecules": n_mol,
        "exact_profile_match": n_exact_profile_match,
        "agreement_pct": round(100.0 * n_exact_profile_match / n_mol, 2) if n_mol else None,
        "examples": examples,
    }


def _atom_multiset(mol):
    return Counter((row[1], row[3], row[5]) for row in mol.atom_table)


def _bond_multiset(mol):
    at = mol.atom_table
    ms = Counter()
    for a1, a2, btype, _arom in mol.bond_table:
        e1, e2 = at[a1][1], at[a2][1]
        ms[(min(e1, e2), max(e1, e2), btype)] += 1
    return ms


def _ring_size_multiset(mol):
    return Counter(len(r) for r in mol.sssr_atom_rings)


def tier5_aromaticity_representation_dependence(smis, chematic, Chem, sample_n, seed):
    # Root-caused via crates/chematic-fp/src/ecfp.rs:initial_atom_id -- the
    # radius-0 invariant explicitly includes `atom.aromatic as u8`. Since
    # `atom.aromatic` is not auto-perceived for Kekule-written SMILES
    # (chematic requires an explicit apply_aromaticity() call, unlike
    # RDKit's default sanitize-on-parse), ecfp4() is representation-
    # dependent for any caller who fingerprints Kekule-spelled input without
    # perceiving aromaticity first -- the same molecule gets a different
    # fingerprint depending on which valid SMILES spelling was used. This
    # tier quantifies how often that actually triggers on a real corpus, and
    # confirms apply_aromaticity() is a complete fix (not just for benzene).
    rng = random.Random(seed)
    sample = smis if len(smis) <= sample_n else rng.sample(smis, sample_n)

    n_checked = 0
    n_naive_mismatch = 0
    n_perceived_still_mismatch = 0
    # Of the still-mismatching cases: do the two spellings agree on the
    # FULL aromaticity ASSIGNMENT, not just aggregate counts? Two earlier,
    # coarser checks were tried and rejected as insufficient to prove a
    # "new defect" label: (1) aromatic_ring_count -- same ring count can
    # hide a different SET of perceived-aromatic atoms in fused systems;
    # (2) plain aromatic-atom/bond COUNTS -- count-equality doesn't imply
    # assignment-equality (two spellings could assign aromaticity to a
    # DIFFERENT set of atoms/bonds while preserving the totals, which the
    # invariant -- computed per-atom/per-bond, not from aggregate counts --
    # would still see as different). The decisive, order-independent check:
    # multiset of (atomic_number, aromatic, degree) per atom, and multiset
    # of (min_element, max_element, bond_type) per bond -- finest-grained
    # check that avoids needing atom-index correspondence between the two
    # spellings (kek_smi is a different string, so atom order can differ).
    # If either multiset differs, the residual is still explained by
    # aromaticity perception (folds into the known aromatic_context bug or
    # an extension of it). If BOTH multisets are identical AND the
    # fingerprint still differs, this is very likely NOT an
    # aromaticity-perception issue -- a genuine, separate defect. NOTE: the
    # multiset is not a per-atom-correspondence proof -- it's blind to a
    # symmetric swap between two atoms of identical (atomic_number, degree)
    # that trade aromaticity flags between spellings (same multiset, but a
    # genuine per-atom perception difference). This reproduces the identical
    # split as the two coarser checks tried before it, so treat the 41/89
    # split as a strong convergent estimate from three independent checks,
    # not an airtight proof -- a rare symmetric swap would land in the
    # perception bucket, slightly deflating the "genuine defect" count.
    n_perception_disagrees = 0
    n_perception_agrees_but_fp_differs = 0
    examples_perception = []
    examples_unattributed = []

    for smi in sample:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        rd_kek = Chem.MolFromSmiles(smi)
        try:
            Chem.Kekulize(rd_kek, clearAromaticFlags=True)
        except Exception:
            continue
        kek_smi = Chem.MolToSmiles(rd_kek, kekuleSmiles=True, canonical=False)
        try:
            cm_arom = chematic.from_smiles(smi)
            cm_kek_naive = chematic.from_smiles(kek_smi)
        except Exception:
            continue
        n_checked += 1

        naive_mismatch = cm_arom.ecfp4() != cm_kek_naive.ecfp4()
        if naive_mismatch:
            n_naive_mismatch += 1

        cm_kek_perceived = cm_kek_naive.apply_aromaticity()
        if cm_arom.ecfp4() != cm_kek_perceived.ecfp4():
            n_perceived_still_mismatch += 1

            same_assignment = (
                _atom_multiset(cm_arom) == _atom_multiset(cm_kek_perceived)
                and _bond_multiset(cm_arom) == _bond_multiset(cm_kek_perceived)
            )
            if not same_assignment:
                n_perception_disagrees += 1
                if len(examples_perception) < 5:
                    examples_perception.append({"smiles": smi})
            else:
                n_perception_agrees_but_fp_differs += 1
                if len(examples_unattributed) < 5:
                    examples_unattributed.append({"smiles": smi, "kekule_smiles": kek_smi})

    return {
        "n_molecules": n_checked,
        "naive_mismatch": n_naive_mismatch,
        "naive_mismatch_pct": round(100.0 * n_naive_mismatch / n_checked, 2) if n_checked else None,
        "apply_aromaticity_mitigated_mismatch": n_perceived_still_mismatch,
        "residual_assignment_multiset_disagrees": n_perception_disagrees,
        "residual_assignment_multiset_agrees_but_fp_differs": n_perception_agrees_but_fp_differs,
        "examples_perception_disagrees": examples_perception,
        "examples_unattributed": examples_unattributed,
    }


def tier6_layer2_shared_mechanism(smis, chematic, Chem, sample_n, seed):
    # Follow-up to tier 5's ~13% post-apply_aromaticity residual: is that
    # residual an ECFP4-specific defect (as an earlier round's writeup
    # assumed), or does it hit other apply_aromaticity()-consuming functions
    # on the SAME molecules -- i.e. one shared root cause, not three
    # coincidentally-similar ones? Uses the SAME (seed, sample_n) draw as
    # tier 5 so the residual-mismatch sets are directly comparable.
    #
    # Two sub-checks against the ECFP4-residual set:
    #   (a) cross-consumer overlap with canonical_smiles/InChI's own
    #       post-apply_aromaticity residual.
    #   (b) SSSR ring-size multiset -- does the residual correspond to a
    #       DIFFERENT ring decomposition, or an identical one (ruling out
    #       ring-finding as the mechanism)?
    #
    # Plus one discriminating check that is NOT restricted to the residual
    # set, because it targets a confound in tier 5's own naive/residual
    # comparison: RDKit's non-canonical Kekule respelling
    # (MolToSmiles(canonical=False)) changes BOTH the aromaticity origin
    # AND the atom traversal order at once. Holding origin fixed (both
    # spellings aromatic-written, neither ever Kekulized) and varying ONLY
    # atom order via RDKit's doRandom respelling isolates plain
    # atom-order-dependence from a Kekule-vs-aromatic-origin-specific
    # effect. If order-only mismatch is near the residual rate, the
    # "aromaticity" framing is likely a red herring for at least part of
    # the residual; if it's near zero, order-dependence is ruled out as an
    # explanation and the residual really does track Kekule-vs-aromatic
    # origin specifically.
    rng = random.Random(seed)
    sample = smis if len(smis) <= sample_n else rng.sample(smis, sample_n)

    n_checked = 0
    ecfp4_residual = set()
    csmi_residual = set()
    inchi_residual = set()
    ring_multiset_differs = 0
    ring_multiset_same = 0

    for smi in sample:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        rd_kek = Chem.MolFromSmiles(smi)
        try:
            Chem.Kekulize(rd_kek, clearAromaticFlags=True)
        except Exception:
            continue
        kek_smi = Chem.MolToSmiles(rd_kek, kekuleSmiles=True, canonical=False)
        try:
            cm_arom = chematic.from_smiles(smi)
            cm_kek_perceived = chematic.from_smiles(kek_smi).apply_aromaticity()
        except Exception:
            continue
        n_checked += 1

        if cm_arom.ecfp4() != cm_kek_perceived.ecfp4():
            ecfp4_residual.add(smi)
            if _ring_size_multiset(cm_arom) != _ring_size_multiset(cm_kek_perceived):
                ring_multiset_differs += 1
            else:
                ring_multiset_same += 1
        if cm_arom.canonical_smiles_mode("normal") != cm_kek_perceived.canonical_smiles_mode("normal"):
            csmi_residual.add(smi)
        try:
            if cm_arom.inchi != cm_kek_perceived.inchi:
                inchi_residual.add(smi)
        except Exception:
            pass

    both_csmi = ecfp4_residual & csmi_residual
    both_inchi = ecfp4_residual & inchi_residual

    # Order-only discriminator, seeded (MolToRandomSmilesVect's randomSeed, not
    # MolToSmiles(doRandom=True) which has no Python-exposed seed) so the
    # reported rates are exactly reproducible, not "usually around this".
    # Also tracks *which* molecules mismatch (not just counts) so we can check
    # whether a consumer's order-sensitive set is disjoint from or a subset of
    # its own apply_aromaticity residual set -- "3.8% order-only" only means
    # "a separate, additive defect" if those molecules are mostly NOT already
    # in the residual set. Assuming disjointness without checking is exactly
    # the kind of unverified-independence claim this project's whole
    # measurement discipline exists to catch.
    n_order_checked = 0
    order_mismatch_ecfp4 = set()
    order_mismatch_csmi = set()
    order_mismatch_inchi = set()
    for smi in sample:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        try:
            respelled = Chem.MolToRandomSmilesVect(rd, 1, randomSeed=seed)[0]
            cm_a = chematic.from_smiles(smi)
            cm_b = chematic.from_smiles(respelled)
        except Exception:
            continue
        n_order_checked += 1
        if cm_a.ecfp4() != cm_b.ecfp4():
            order_mismatch_ecfp4.add(smi)
        if cm_a.canonical_smiles_mode("normal") != cm_b.canonical_smiles_mode("normal"):
            order_mismatch_csmi.add(smi)
        try:
            if cm_a.inchi != cm_b.inchi:
                order_mismatch_inchi.add(smi)
        except Exception:
            pass

    csmi_order_inside_residual = order_mismatch_csmi & csmi_residual
    csmi_order_outside_residual = order_mismatch_csmi - csmi_residual
    inchi_order_inside_residual = order_mismatch_inchi & inchi_residual
    inchi_order_outside_residual = order_mismatch_inchi - inchi_residual
    # InChI-specific excess: order-only mismatches InChI has that canonical_smiles
    # doesn't, on the same molecules -- isolates a defect unique to InChI's own code
    # path from ordering-layer churn already known/accounted for via canonical_smiles.
    inchi_specific_excess = order_mismatch_inchi - order_mismatch_csmi

    return {
        "n_molecules": n_checked,
        "ecfp4_residual_count": len(ecfp4_residual),
        "canonical_smiles_residual_count": len(csmi_residual),
        "inchi_residual_count": len(inchi_residual),
        "ecfp4_and_canonical_smiles_overlap": len(both_csmi),
        "ecfp4_and_canonical_smiles_overlap_pct_of_ecfp4": (
            round(100.0 * len(both_csmi) / len(ecfp4_residual), 1) if ecfp4_residual else None
        ),
        "ecfp4_and_inchi_overlap": len(both_inchi),
        "ecfp4_and_inchi_overlap_pct_of_ecfp4": (
            round(100.0 * len(both_inchi) / len(ecfp4_residual), 1) if ecfp4_residual else None
        ),
        "ring_size_multiset_differs_in_ecfp4_residual": ring_multiset_differs,
        "ring_size_multiset_same_in_ecfp4_residual": ring_multiset_same,
        "order_only_n_molecules": n_order_checked,
        "order_only_ecfp4_mismatch": len(order_mismatch_ecfp4),
        "order_only_ecfp4_mismatch_pct": (
            round(100.0 * len(order_mismatch_ecfp4) / n_order_checked, 2) if n_order_checked else None
        ),
        "order_only_canonical_smiles_mismatch": len(order_mismatch_csmi),
        "order_only_canonical_smiles_mismatch_pct": (
            round(100.0 * len(order_mismatch_csmi) / n_order_checked, 2) if n_order_checked else None
        ),
        "order_only_inchi_mismatch": len(order_mismatch_inchi),
        "order_only_inchi_mismatch_pct": (
            round(100.0 * len(order_mismatch_inchi) / n_order_checked, 2) if n_order_checked else None
        ),
        "canonical_smiles_order_mismatch_inside_its_own_residual": len(csmi_order_inside_residual),
        "canonical_smiles_order_mismatch_outside_its_own_residual": len(csmi_order_outside_residual),
        "inchi_order_mismatch_inside_its_own_residual": len(inchi_order_inside_residual),
        "inchi_order_mismatch_outside_its_own_residual": len(inchi_order_outside_residual),
        "inchi_specific_excess_over_canonical_smiles": len(inchi_specific_excess),
        "inchi_specific_excess_pct_of_inchi_order_mismatch": (
            round(100.0 * len(inchi_specific_excess) / len(order_mismatch_inchi), 1)
            if order_mismatch_inchi
            else None
        ),
    }


def tier4_connectivity_sanity_check(smis, chematic, Chem, sample_n, seed):
    rng = random.Random(seed)
    sample = smis if len(smis) <= sample_n else rng.sample(smis, sample_n)

    n_mol = 0
    n_atom_radius_checks = 0
    n_match = 0
    examples = []
    for smi in sample:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        try:
            cm = chematic.from_smiles(smi)
        except Exception:
            continue
        if len(cm.atom_table) != rd.GetNumAtoms():
            continue  # heavy-atom parse mismatch -- not what this tier measures
        n_mol += 1

        adj_c = chematic_adjacency(cm)
        adj_r = rdkit_adjacency(rd)

        for atom_idx in range(len(cm.atom_table)):
            for radius in (1, 2):
                ball_c = bfs_ball(adj_c, atom_idx, radius)
                ball_r = bfs_ball(adj_r, atom_idx, radius)
                n_atom_radius_checks += 1
                if ball_c == ball_r:
                    n_match += 1
                elif len(examples) < 10:
                    examples.append(
                        {
                            "smiles": smi,
                            "atom_idx": atom_idx,
                            "radius": radius,
                            "chematic_ball": sorted(ball_c),
                            "rdkit_ball": sorted(ball_r),
                        }
                    )
    return {
        "n_molecules": n_mol,
        "n_atom_radius_checks": n_atom_radius_checks,
        "n_match": n_match,
        "agreement_pct": round(100.0 * n_match / n_atom_radius_checks, 4)
        if n_atom_radius_checks
        else None,
        "examples": examples,
    }


def tier3_similarity_correlation(smis, chematic, Chem, AllChem, DataStructs, sample_n, seed):
    rng = random.Random(seed)
    sample = smis if len(smis) <= sample_n else rng.sample(smis, sample_n)

    chem_fps = []
    rd_fps = []
    for smi in sample:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        try:
            cm = chematic.from_smiles(smi)
        except Exception:
            continue
        chem_fps.append(cm.ecfp4())
        rd_fps.append(AllChem.GetMorganFingerprintAsBitVect(rd, 2, 2048))

    chem_sims = []
    rd_sims = []
    n = len(chem_fps)
    for i in range(n):
        for j in range(i + 1, n):
            chem_sims.append(chematic.tanimoto(chem_fps[i], chem_fps[j]))
            rd_sims.append(DataStructs.TanimotoSimilarity(rd_fps[i], rd_fps[j]))

    corr = statistics.correlation(chem_sims, rd_sims) if len(chem_sims) > 1 else None
    mean_abs_diff = (
        sum(abs(a - b) for a, b in zip(chem_sims, rd_sims)) / len(chem_sims)
        if chem_sims
        else None
    )
    return {
        "n_molecules": n,
        "n_pairs": len(chem_sims),
        "pearson_correlation": round(corr, 4) if corr is not None else None,
        "mean_abs_tanimoto_diff": round(mean_abs_diff, 4) if mean_abs_diff is not None else None,
    }


def tier0_raw_bit_equality(smis, chematic, Chem, AllChem, sample_n, seed):
    rng = random.Random(seed)
    sample = smis if len(smis) <= sample_n else rng.sample(smis, sample_n)

    agree_fracs = []
    for smi in sample:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        try:
            cm = chematic.from_smiles(smi)
        except Exception:
            continue
        chem_fp = cm.ecfp4()
        chem_bits = "".join(f"{b:08b}"[::-1] for b in chem_fp)  # LSB-first per byte
        rd_bv = AllChem.GetMorganFingerprintAsBitVect(rd, 2, 2048)
        rd_bits = rd_bv.ToBitString()
        agree = sum(1 for a, b in zip(chem_bits, rd_bits) if a == b) / len(rd_bits)
        agree_fracs.append(agree)

    return {
        "n_molecules": len(agree_fracs),
        "mean_per_position_agreement_pct": round(100.0 * sum(agree_fracs) / len(agree_fracs), 2)
        if agree_fracs
        else None,
        "note": (
            "NOT a correctness signal -- expected to look high due to sparse-vector "
            "0/0 matches dominating, or near-random, depending on density; hash "
            "functions differ by design (FNV-1a vs RDKit), so bit POSITIONS were "
            "never meant to align. See tiers 1-3 for the real metrics."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("smiles_csv", nargs="?", default="~/Downloads/SMILES.csv")
    parser.add_argument("--limit", type=int, default=None, help="Cap corpus size (tiers 1/2, cheap)")
    parser.add_argument("--bit-sample", type=int, default=300, help="Molecules for tier 0")
    parser.add_argument("--pairs-sample", type=int, default=300, help="Molecules for tier 3 (pairwise, O(n^2))")
    parser.add_argument("--connectivity-sample", type=int, default=300, help="Molecules for tier 4")
    parser.add_argument("--aromaticity-sample", type=int, default=300, help="Molecules for tier 5")
    parser.add_argument("--layer2-sample", type=int, default=300,
                         help="Molecules for tier 6 (uses the same seed as tier 5, so pass the "
                              "same value as --aromaticity-sample to compare identical residual sets)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--json", default=None)
    args = parser.parse_args()

    try:
        from rdkit import Chem, RDLogger, DataStructs
        from rdkit.Chem import AllChem, rdFingerprintGenerator

        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        sys.exit("rdkit not installed. pip install rdkit")
    import chematic

    path = os.path.expanduser(args.smiles_csv)
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    smis = [l.split(",")[0].strip() for l in lines if l.split(",")[0].strip().lower() != "smiles"]

    print(f"Corpus: {path} ({len(smis)} SMILES)")
    print()

    print("Tier 0 -- raw bit-vector equality (NOT a correctness signal, see note)...")
    t0 = tier0_raw_bit_equality(smis, chematic, Chem, AllChem, args.bit_sample, args.seed)
    print(f"  mean per-position agreement: {t0['mean_per_position_agreement_pct']}% "
          f"(n={t0['n_molecules']})")
    print(f"  {t0['note']}")
    print()

    print("Tier 1 -- coverage parity (per-(atom,radius) environment existence, "
          "RDKit includeRedundantEnvironments=True for a fair chemistry-only comparison)...")
    t1 = tier1_coverage_parity(smis, chematic, Chem, rdFingerprintGenerator, args.limit)
    print(f"  {t1['exact_coverage_match']}/{t1['n_molecules']} molecules "
          f"({t1['agreement_pct']}%) have identical (atom,radius) coverage sets")
    if t1["examples"]:
        print(f"  example mismatch: {t1['examples'][0]}")
    print()

    print("Tier 2 -- invariant partition agreement (the real chemistry check: do the two "
          "implementations group environments into identical-size equivalence classes)...")
    t2 = tier2_invariant_partition_agreement(smis, chematic, Chem, rdFingerprintGenerator, args.limit)
    print(f"  {t2['exact_profile_match']}/{t2['n_molecules']} molecules "
          f"({t2['agreement_pct']}%) have an identical invariant-equivalence-class profile")
    if t2["examples"]:
        print(f"  example mismatch: {t2['examples'][0]}")
    print()

    print("Tier 3 -- similarity-structure preservation (pairwise Tanimoto correlation)...")
    t3 = tier3_similarity_correlation(smis, chematic, Chem, AllChem, DataStructs, args.pairs_sample, args.seed)
    print(f"  Pearson r = {t3['pearson_correlation']}, mean |Δ Tanimoto| = "
          f"{t3['mean_abs_tanimoto_diff']} over {t3['n_pairs']} pairs "
          f"({t3['n_molecules']} molecules)")
    print()

    print("Tier 4 -- connectivity sanity check (auxiliary; parser agreement, NOT a fingerprint test)...")
    t4 = tier4_connectivity_sanity_check(smis, chematic, Chem, args.connectivity_sample, args.seed)
    print(f"  {t4['n_match']}/{t4['n_atom_radius_checks']} atom-radius checks match "
          f"({t4['agreement_pct']}%) across {t4['n_molecules']} molecules")
    if t4["examples"]:
        print(f"  example mismatch: {t4['examples'][0]}")
    print()

    print("Tier 5 -- aromaticity representation-dependence (chematic self-consistency, "
          "practical consequence of tier 2)...")
    t5 = tier5_aromaticity_representation_dependence(smis, chematic, Chem, args.aromaticity_sample, args.seed)
    print(f"  naive (no apply_aromaticity): {t5['naive_mismatch']}/{t5['n_molecules']} "
          f"({t5['naive_mismatch_pct']}%) molecules get a DIFFERENT ecfp4() for the Kekule "
          f"spelling vs. the aromatic spelling of the same molecule")
    print(f"  after apply_aromaticity(): {t5['apply_aromaticity_mitigated_mismatch']}/"
          f"{t5['n_molecules']} still mismatching, of which:")
    print(f"    {t5['residual_assignment_multiset_disagrees']} also disagree on the full "
          f"aromatic-atom/bond assignment (still aromaticity perception -- known "
          f"aromatic_context bug or an extension)")
    print(f"    {t5['residual_assignment_multiset_agrees_but_fp_differs']} have an IDENTICAL "
          f"aromatic-atom/bond assignment but still a different fingerprint (very likely not "
          f"perception -- genuine separate defect; see tier 5's code comment on the multiset's "
          f"one known blind spot)")
    if t5["examples_unattributed"]:
        print(f"  unattributed example: {t5['examples_unattributed'][0]}")
    print()

    print("Tier 6 -- Layer 2 shared-mechanism check (is tier 5's residual ECFP4-specific, "
          "or does it hit canonical_smiles/InChI on the SAME molecules)...")
    t6 = tier6_layer2_shared_mechanism(smis, chematic, Chem, args.layer2_sample, args.seed)
    print(f"  residual counts (n={t6['n_molecules']}): ecfp4={t6['ecfp4_residual_count']}, "
          f"canonical_smiles={t6['canonical_smiles_residual_count']}, "
          f"inchi={t6['inchi_residual_count']}")
    print(f"  ecfp4 ∩ canonical_smiles residual = {t6['ecfp4_and_canonical_smiles_overlap']} "
          f"({t6['ecfp4_and_canonical_smiles_overlap_pct_of_ecfp4']}% of ecfp4's residual set)")
    print(f"  ecfp4 ∩ inchi residual = {t6['ecfp4_and_inchi_overlap']} "
          f"({t6['ecfp4_and_inchi_overlap_pct_of_ecfp4']}% of ecfp4's residual set)")
    print(f"  SSSR ring-size multiset, within ecfp4's residual set: "
          f"{t6['ring_size_multiset_differs_in_ecfp4_residual']} differ, "
          f"{t6['ring_size_multiset_same_in_ecfp4_residual']} identical "
          f"(identical yet ecfp4() still differs -- rules out ring decomposition)")
    print(f"  order-only discriminator (two aromatic-preserving respellings, no Kekule "
          f"involved, seeded, n={t6['order_only_n_molecules']}):")
    print(f"    ecfp4 mismatch: {t6['order_only_ecfp4_mismatch']} "
          f"({t6['order_only_ecfp4_mismatch_pct']}%)")
    print(f"    canonical_smiles mismatch: {t6['order_only_canonical_smiles_mismatch']} "
          f"({t6['order_only_canonical_smiles_mismatch_pct']}%)")
    print(f"    inchi mismatch: {t6['order_only_inchi_mismatch']} "
          f"({t6['order_only_inchi_mismatch_pct']}%)")
    print(f"  of canonical_smiles's order-only mismatches, "
          f"{t6['canonical_smiles_order_mismatch_inside_its_own_residual']} fall INSIDE its "
          f"own apply_aromaticity residual set and "
          f"{t6['canonical_smiles_order_mismatch_outside_its_own_residual']} fall OUTSIDE it "
          f"-- tests whether order-sensitivity is a separate defect (mostly outside) or "
          f"entangled with the shared residual (mostly inside)")
    print(f"  of inchi's order-only mismatches, "
          f"{t6['inchi_order_mismatch_inside_its_own_residual']} fall INSIDE its own "
          f"apply_aromaticity residual set and "
          f"{t6['inchi_order_mismatch_outside_its_own_residual']} fall OUTSIDE it")
    print(f"  inchi-specific excess (order-only mismatches inchi has that "
          f"canonical_smiles doesn't, same molecules): "
          f"{t6['inchi_specific_excess_over_canonical_smiles']} "
          f"({t6['inchi_specific_excess_pct_of_inchi_order_mismatch']}% of inchi's own "
          f"order-only mismatch set) -- isolates a defect unique to inchi's code path "
          f"from ordering-layer churn already shared with canonical_smiles")
    print()

    result = {"tier0_raw_bit_equality": t0, "tier1_coverage_parity": t1,
              "tier2_invariant_partition_agreement": t2, "tier3_similarity_correlation": t3,
              "tier4_connectivity_sanity_check": t4,
              "tier5_aromaticity_representation_dependence": t5,
              "tier6_layer2_shared_mechanism": t6}
    if args.json:
        with open(args.json, "w") as f:
            json.dump(result, f, indent=2)
        print(f"Wrote {args.json}")


if __name__ == "__main__":
    main()
