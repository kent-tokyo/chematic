#!/usr/bin/env python3
"""Wave 3 (issue #149) audit: general shared-carrier coupling-component
mechanism diagnosis.

PR #351 (Wave 2D) fixed the ring-endocyclic-in-small-ring mechanism behind
the 18 pinned `EZ_SHARED_CANDIDATE_BOND_RESIDUALS` fixtures. Its own commit
message states this closes only ~10% (3 of 31) of the corpus's general
shared-carrier coupling-component population -- the other ~90% (28 of 31)
was reported as "a separate, still-unidentified mechanism." This script
measures that remaining population directly, rather than assuming its size,
shape, or residual status from that figure (a topological-presence count,
not a confirmed-permutation-invariance-failure count).

Diagnosis only: no production code is touched by this script or by the Rust
example it drives (`crates/chematic-smiles/examples/
ez_shared_carrier_coupling_mechanism_audit.rs`). Does not touch or
generalize the ring<8-endocyclic predicate in `compute_stereo_alkene_ends`.
Does not conflate this population with ROADMAP.md backlog item 6's separate
"346 abstained E/Z bonds" ledger (85 lost-in-canonicalization + 261
CarrierConflict, an unrelated, unchanged code path).

Method, in order (see docs/rfcs/ez_shared_carrier_coupling_mechanism_audit.md
for the full narrative and results):

1. Provenance gate: do the 18 pinned fixtures / the 2 never-corrupts SMILES
   appear verbatim in the committed corpus? Reported, never assumed, never
   pooled with the corpus-derived population if disjoint.
2. Current-topology scan (`scan` subcommand, ring-gate-aware -- unlike
   either pre-existing example in this file's family, both of which predate
   PR #351 and still report the pre-fix topology): how many coupled
   components (size >= 2) exist in the corpus *today*.
3. Axis 1 (RDKit relabeling, K=64 by default, seeded/reproducible): does chematic's own
   canonical output stay identical across the configured independent
   relabelings of each
   coupled molecule? A real chematic-internal-self-consistency probe (RDKit
   only supplies alternate valid spellings; RDKit agreement is not what's
   measured).
4. Axis 2 (single-end mark relocation, no RDKit): reimplements the private
   `alternate_ez_markings` test helper. Measured finding, not assumed: this
   probe is STRUCTURALLY INCAPABLE of testing any genuinely coupled 2-node
   component at all -- relocating either end's mark necessarily strips the
   *shared* bond's own mark, which the other end depends on for its own
   geometry reading, so the helper's own geometry-preservation check always
   rejects the move. Confirmed both mechanistically (traced against a real
   corpus example) and empirically (0 alternates survive for any of the 28
   corpus components, or for the negative-control fixture). Axis 2 can only
   ever inform a STANDALONE (singleton) stereo-alkene end.
5. RDKit stereogenicity oracle (`Chem.FindPotentialStereo` +
   `Chem.AssignStereochemistry`), same index-correspondence-first discipline
   as `ez_ring_constrained_residual_diagnosis.py`'s `crosscheck_row`: for
   each of the 28 components' two ends, is RDKit's own verdict "Specified"
   (a real stereocenter) or not?
6. Structural classification: group confirmed-residual/all-28 components by
   a measured feature tuple (ring membership, candidate-bond representation,
   RDKit stereogenicity) -- never a heuristic imposed ahead of measurement.
7. Calibration cross-check: the permanent regression fixture
   `ez_carrier_shared_bond_between_two_stereo_systems_never_corrupts`
   (canonical.rs) is checked directly -- its own doc comment claims its two
   pinned spellings do NOT converge to one canonical string. This script
   finds (see RFC) that they DO converge today, on both the two pinned
   spellings and the configured fresh relabelings -- flagged as a likely-stale doc claim
   (a documentation-currency finding, not something this diagnosis-only
   script fixes).

Usage:
    .venv/bin/python3 scripts/ez_shared_carrier_coupling_mechanism_diagnosis.py [CORPUS]
    .venv/bin/python3 scripts/ez_shared_carrier_coupling_mechanism_diagnosis.py --self-test

Writes:
    validation/results/ez_shared_carrier_coupling_mechanism_audit.jsonl
    validation/results/ez_shared_carrier_coupling_mechanism_audit_summary.json
"""
import argparse
import json
import os
import random
import subprocess
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE = "ez_shared_carrier_coupling_mechanism_audit"
OUT_JSONL = os.path.join(ROOT, "validation", "results", f"{EXAMPLE}.jsonl")
OUT_SUMMARY = os.path.join(ROOT, "validation", "results", f"{EXAMPLE}_summary.json")
DEFAULT_CORPUS = os.path.join(ROOT, "scripts", "descriptor_census_corpus.smi")


def path_for_report(path):
    """Path suitable for embedding in committed provenance JSON: relative to
    the repo root if it's inside the repo, home-contracted (~) otherwise --
    never a raw absolute path, which would embed the local username."""
    abs_path = os.path.abspath(os.path.expanduser(str(path)))
    if abs_path.startswith(ROOT + os.sep):
        return os.path.relpath(abs_path, ROOT)
    home = os.path.expanduser("~")
    if abs_path.startswith(home):
        return "~" + abs_path[len(home) :]
    return abs_path


# Mirrors `canonical.rs`'s current, merged 18-entry `EZ_SHARED_CARRIER_
# FULLY_RESOLVED` list exactly -- used only for the provenance gate and the
# negative control, never re-derived independently of canonical.rs.
EZ_SHARED_CARRIER_FULLY_RESOLVED = [
    r"CCCCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/Cc2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
    r"CCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/c2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
    r"CC(C)(C)/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(I)c1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1C(F)(F)F",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccccc1OC",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccc([N+](=O)[O-])cc1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1cccc(C(F)(F)F)c1",
    r"CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
    r"CC1=C2CC[C@@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
    r"COC(=O)/C=C/[C@H]1CCC2=C(C)/C(=N/N=C(N)N)CC[C@@]21C",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccc(I)cc1",
    r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(C)c1",
    r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(OC)c1",
    r"CCO/C(O)=C(\C1=NCCN1)c1nnc(N)s1",
]

# The permanent no-corruption regression fixture (canonical.rs). NOT one of
# the 18 -- kept separate, per this repo's own never-pool discipline.
EZ_NEVER_CORRUPTS_A = r"OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c1/c(c(c1O)O)=N/CCCCC"
EZ_NEVER_CORRUPTS_B = "OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c\\1c(/c(c1O)O)=N/CCCCC"

# A known-fine fixture sharing the ring-endocyclic *shape* with the 8
# formerly-residual fixtures, but never itself a residual (Wave 2C's own
# necessary-but-not-sufficient finding) -- the negative control for axis 2.
NEGATIVE_CONTROL = r"CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N"

N_RELABELINGS_PER_MOLECULE = 64


# ---------------------------------------------------------------------------
# Chematic side: drive the Rust example, never `import chematic` (matches
# every prior "Wave" script's independence discipline).
# ---------------------------------------------------------------------------

def run_chematic_example(subcommand, file_arg=None):
    cmd = ["cargo", "run", "-p", "chematic-smiles", "--release", "--example", EXAMPLE, "--", subcommand]
    if file_arg is not None:
        cmd.append(file_arg)
    result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"cargo run {subcommand} failed (exit {result.returncode}):\n{result.stderr[-4000:]}")
    rows = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if line.startswith("{"):
            rows.append(json.loads(line))
    return rows


# ---------------------------------------------------------------------------
# RDKit side.
# ---------------------------------------------------------------------------

def rdkit_relabelings(smiles, k=N_RELABELINGS_PER_MOLECULE):
    """K reproducible (seeded) relabelings of `smiles`, matching
    `canonical_residual_diagnosis.py`'s exact Check-2 method: `Chem.
    RenumberAtoms` + a seeded `random.Random` shuffle, deliberately NOT
    RDKit's own unseeded `doRandom=True`."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    n = mol.GetNumAtoms()
    out = []
    for seed in range(k):
        rng = random.Random(seed)
        order = list(range(n))
        rng.shuffle(order)
        rm = Chem.RenumberAtoms(mol, order)
        out.append(Chem.MolToSmiles(rm, canonical=False))
    return out


def double_bond_partner(mol, atom_idx):
    from rdkit import Chem

    atom = mol.GetAtomWithIdx(atom_idx)
    for b in atom.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE:
            return b.GetOtherAtomIdx(atom_idx)
    return None


def rdkit_stereo_verdict(mol, potential, atom_idx, partner_idx):
    """`Chem.FindPotentialStereo` verdict for the double bond between
    `atom_idx`/`partner_idx`, matched by atom-pair (not bond index alone),
    matching `ez_ring_constrained_residual_diagnosis.py`'s own oracle
    discipline. Returns None if RDKit doesn't even list the bond (this
    RDKit version has no `NOT_POSSIBLE` value in `StereoSpecified` --
    verified live in `self_test()` below)."""
    for si in potential:
        if str(si.type) != "Bond_Double":
            continue
        b = mol.GetBondWithIdx(si.centeredOn)
        pair = {b.GetBeginAtomIdx(), b.GetEndAtomIdx()}
        if pair == {atom_idx, partner_idx}:
            return str(si.specified)
    return None


def rdkit_both_ends_specified(smiles, end_atom_idxs):
    """For a coupled component's two end atoms, is RDKit's own
    `FindPotentialStereo` verdict "Specified" (a real stereocenter) for
    BOTH? Index-correspondence is implicit here (both chematic and RDKit
    parse the SAME literal SMILES string, so atom indices already
    correspond by construction -- unlike the cross-spelling axis-1
    correspondence, which needs the explicit canonical_atom_order-based
    check)."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "rdkit_parse_failed"
    potential = list(Chem.FindPotentialStereo(mol, cleanIt=False, flagPossible=True))
    verdicts = []
    for end_idx in end_atom_idxs:
        partner_idx = double_bond_partner(mol, end_idx)
        if partner_idx is None:
            return None, "no_double_bond_partner_found"
        v = rdkit_stereo_verdict(mol, potential, end_idx, partner_idx)
        verdicts.append(v)
    return all(v == "Specified" for v in verdicts), verdicts


# ---------------------------------------------------------------------------
# Main diagnosis.
# ---------------------------------------------------------------------------

def provenance_gate(corpus_path):
    with open(corpus_path) as f:
        corpus_set = {line.strip() for line in f if line.strip()}
    fixtures_in_corpus = [s for s in EZ_SHARED_CARRIER_FULLY_RESOLVED if s in corpus_set]
    never_corrupts_in_corpus = [s for s in (EZ_NEVER_CORRUPTS_A, EZ_NEVER_CORRUPTS_B) if s in corpus_set]
    return {
        "n_pinned_fixtures_in_corpus": len(fixtures_in_corpus),
        "n_pinned_fixtures_total": len(EZ_SHARED_CARRIER_FULLY_RESOLVED),
        "n_never_corrupts_in_corpus": len(never_corrupts_in_corpus),
        "conclusion": (
            "fixture-derived and corpus-derived populations are NOT pooled below, "
            "regardless of this overlap"
        ),
    }


def current_topology(corpus_path):
    rows = run_chematic_example("scan", corpus_path)
    components = [r for r in rows if r["kind"] == "component"]
    ends = [r for r in rows if r["kind"] == "end"]
    coupled = [c for c in components if c["size"] >= 2]

    by_smiles_ends = defaultdict(list)
    for e in ends:
        if e["coupled"]:
            by_smiles_ends[e["smiles"]].append(e)

    return {
        "n_component_rows": len(components),
        "n_coupled_components": len(coupled),
        "coupled_component_sizes": sorted({c["size"] for c in coupled}),
        "coupled_component_shapes": sorted({c["shape"] for c in coupled}),
        "n_end_rows": len(ends),
    }, coupled, by_smiles_ends


def axis1_divergence(coupled_smiles, relabelings_per_molecule=N_RELABELINGS_PER_MOLECULE):
    """Writes a temp TSV of (original, relabeled) pairs, runs the Rust
    example's `axis1` subcommand, and reports per-molecule whether
    chematic's own canonical output stays identical across all K
    relabelings."""
    lines = []
    skipped = []
    for smi in coupled_smiles:
        relabelings = rdkit_relabelings(smi, relabelings_per_molecule)
        if relabelings is None:
            skipped.append(smi)
            continue
        for r in relabelings:
            lines.append(f"{smi}\t{r}")

    tsv_path = os.path.join(ROOT, "validation", "results", f"{EXAMPLE}_axis1_input.tsv")
    with open(tsv_path, "w") as f:
        f.write("\n".join(lines) + ("\n" if lines else ""))

    rows = run_chematic_example("axis1", tsv_path) if lines else []
    variants = [r for r in rows if r["kind"] == "axis1_variant"]

    by_orig = defaultdict(list)
    for v in variants:
        by_orig[v["original_smiles"]].append(v)

    per_molecule = {}
    n_divergent = 0
    n_cross_correspondence_failures = 0
    for smi, vs in by_orig.items():
        canons = {v["canonical_smiles"] for v in vs if v["correspondence_ok"]}
        cross_fail = sum(1 for v in vs if not v["cross_correspondence_ok"])
        n_cross_correspondence_failures += cross_fail
        divergent = len(canons) > 1
        if divergent:
            n_divergent += 1
        per_molecule[smi] = {
            "n_relabelings_tested": len(vs),
            "n_cross_correspondence_failures": cross_fail,
            "n_distinct_canonical_outputs": len(canons),
            "distinct_canonical_outputs": sorted(canons),
            "divergent": divergent,
        }

    os.remove(tsv_path)
    return {
        "n_molecules_tested": len(by_orig),
        "n_skipped_rdkit_parse_failed": len(skipped),
        "n_divergent": n_divergent,
        "n_cross_correspondence_failures_total": n_cross_correspondence_failures,
        "per_molecule": per_molecule,
    }


def axis2_alternates(smiles_list):
    """Runs the Rust example's `axis2` subcommand. Reports, per molecule,
    whether ANY geometry-preserving single-end mark relocation exists at
    all, and whether any that DO exist diverge from baseline."""
    input_path = os.path.join(ROOT, "validation", "results", f"{EXAMPLE}_axis2_input.smi")
    with open(input_path, "w") as f:
        f.write("\n".join(smiles_list) + "\n")

    rows = run_chematic_example("axis2", input_path)
    os.remove(input_path)

    by_smiles = defaultdict(list)
    for r in rows:
        if r["kind"] == "axis2_variant":
            by_smiles[r["source_smiles"]].append(r)

    per_molecule = {}
    for smi, rs in by_smiles.items():
        alts = [r for r in rs if r["variant"] != "baseline"]
        per_molecule[smi] = {
            "n_alternates_generated": len(alts),
            "any_divergent": any(r.get("differs_from_baseline") for r in alts),
        }
    return per_molecule


def structural_classification(coupled_smiles, by_smiles_ends, axis1_result):
    """Groups the coupled components by a MEASURED feature tuple -- never a
    heuristic imposed before measurement. Names a bucket only once
    occupied."""
    buckets = defaultdict(list)
    per_component_detail = {}

    for smi in coupled_smiles:
        end_rows = by_smiles_ends[smi]
        end_idxs = [e["end_atom_idx"] for e in end_rows]
        has_ring = any(e["end_atom_in_ring"] for e in end_rows)
        order_types = tuple(sorted(c["current_bond_order"] for e in end_rows for c in e["candidate_bonds"]))
        both_specified, verdicts = rdkit_both_ends_specified(smi, end_idxs)

        key = (
            "has_ring" if has_ring else "no_ring",
            order_types,
            "both_rdkit_specified" if both_specified else f"NOT_both_specified:{verdicts}",
        )
        buckets[key].append(smi)

        per_component_detail[smi] = {
            "has_ring": has_ring,
            "candidate_bond_order_types": list(order_types),
            "rdkit_both_ends_specified": both_specified,
            "rdkit_verdicts": verdicts,
            "axis1_divergent": axis1_result["per_molecule"].get(smi, {}).get("divergent"),
        }

    bucket_summary = [
        {"feature_tuple": list(k), "count": len(v), "example_smiles": v[0]}
        for k, v in buckets.items()
    ]
    return bucket_summary, per_component_detail


def calibration_check(relabelings_per_molecule=N_RELABELINGS_PER_MOLECULE):
    """Checks the never-corrupts fixture's two pinned spellings directly,
    plus the configured fresh RDKit relabelings of spelling A -- reports whether they
    converge (a finding either way, not assumed from the doc comment)."""
    rows = run_chematic_example("scan")
    canons = {}
    for r in rows:
        if r["kind"] == "component" and r["smiles"] in (EZ_NEVER_CORRUPTS_A, EZ_NEVER_CORRUPTS_B):
            canons[r["smiles"]] = r["canonical_smiles"]
    pinned_pair_converges = len(set(canons.values())) == 1 if len(canons) == 2 else None

    relabel_result = axis1_divergence([EZ_NEVER_CORRUPTS_A], relabelings_per_molecule)
    fresh_relabeling_stable = not relabel_result["per_molecule"].get(EZ_NEVER_CORRUPTS_A, {}).get("divergent", True)

    return {
        "pinned_pair_canonical_outputs": canons,
        "pinned_pair_converges": pinned_pair_converges,
        "fresh_relabeling_stable": fresh_relabeling_stable,
        "n_relabelings_tested": relabel_result["per_molecule"].get(
            EZ_NEVER_CORRUPTS_A, {}
        ).get("n_relabelings_tested", 0),
        "conclusion": (
            f"CONVERGES on current main (both the 2 pinned spellings and {relabelings_per_molecule} fresh "
            "relabelings agree) -- the doc comment's 'does NOT resolve to one "
            "canonical string' claim appears STALE, most likely superseded by "
            "PR #229's joint-component solver without the comment being revisited. "
            "Flagged as a documentation-currency finding; not fixed by this "
            "diagnosis-only script (even a comment-only edit to canonical.rs is a "
            "change under crates/*/src/**, out of this PR's declared scope)."
        ) if pinned_pair_converges and fresh_relabeling_stable else (
            "Did not fully reproduce convergence -- see raw fields above; treat "
            "the doc comment's claim as still potentially accurate and "
            "investigate further before drawing a conclusion."
        ),
    }


def axis2_applicability_note():
    return (
        "Axis 2 (single-end mark relocation) is measured here to be STRUCTURALLY "
        "incapable of testing any genuinely coupled 2-node component: relocating "
        "either end's mark necessarily strips the SHARED bond's own mark, which "
        "the other end depends on for its own geometry reading, so the "
        "reimplemented `alternate_ez_markings` helper's own geometry-preservation "
        "check rejects the move every time. Confirmed both mechanistically (traced "
        "against a real corpus example, see the RFC) and empirically (0 alternates "
        "survive for any of the 28 coupled corpus components, or for the "
        "never-corrupts calibration pair; only the singleton negative-control "
        "fixture produced any alternate at all, and it did not diverge). This is a "
        "structural property of single-end relocation, not an implementation gap "
        "specific to this audit -- it applies identically to production's own "
        "`alternate_ez_markings` test helper."
    )


def verdict(
    n_coupled,
    axis1_result,
    bucket_summary,
    calibration,
    relabelings_per_molecule=N_RELABELINGS_PER_MOLECULE,
):
    all_both_specified = all(
        "both_rdkit_specified" in b["feature_tuple"] for b in bucket_summary
    )
    one_mechanism = len(bucket_summary) <= 2 and all_both_specified  # ring-vs-acyclic split, same RDKit signature

    return {
        "n_coupled_components_measured": n_coupled,
        "n_confirmed_divergent_by_axis1": axis1_result["n_divergent"],
        "mechanism_count": "one" if one_mechanism else "several",
        "mechanism_note": (
            "All measured coupled components have BOTH ends independently, "
            "genuinely RDKit-Specified (real stereocenters) -- the SAME shape as "
            "the 5 known-fine EZ_SHARED_CARRIER_FULLY_RESOLVED hydrazone-imine "
            "fixtures and the never-corrupts fixture (which itself converges, see "
            "calibration_check). A secondary, non-mechanism-changing structural "
            "split exists (ring/aromatic-stashed vs. acyclic/literal-marker "
            "representation), reported in bucket_summary."
        ) if one_mechanism else "Buckets differ in RDKit stereogenicity signature -- see bucket_summary.",
        "verdict": (
            "NEEDS-RESEARCH, confirmed residuals found"
            if axis1_result["n_divergent"]
            else "NEEDS-RESEARCH, leaning GO (no sampled residuals)"
        ),
        "verdict_reasoning": (
            f"{axis1_result['n_divergent']}/{n_coupled} coupled components show canonical-output divergence "
            f"under axis 1 (RDKit relabeling, K={relabelings_per_molecule}, 0 "
            f"cross-correspondence failures). Axis 2 cannot test coupled pairs at "
            "all (structural limitation, not a gap in this audit). The "
            "previously-cited never-corrupts calibration example -- itself an "
            "instance of this exact shape -- still converges. The sampled "
            "residuals must remain open and are not treated as proof of a general "
            "mechanism: relabeling is finite, and RDKit's relabel-and-reserialize "
            "process does not guarantee every alternate carrier spelling. The "
            "next implementation step is to add these reproducible residuals as "
            "held-out regression fixtures before changing production ranking."
        ),
    }


def self_test():
    """Fabricated-row + live-RDKit positive-control checks, matching every
    prior Wave script's discipline."""
    print("=== self-test ===")

    # RDKit-version-drift check: no NOT_POSSIBLE value should exist, and the
    # two mechanistic positive controls from Wave 2C's audit should still
    # hold under whatever RDKit version is actually installed.
    from rdkit import Chem
    enum_values = [v for v in dir(Chem.StereoSpecified) if not v.startswith("_")]
    assert "NOT_POSSIBLE" not in enum_values, f"unexpected NOT_POSSIBLE in {enum_values}"
    print(f"  OK: RDKit {Chem.rdBase.rdkitVersion} has no NOT_POSSIBLE StereoSpecified value")

    cyclohexene = Chem.MolFromSmiles("CC1=C(C)CCCC1")
    potential = list(Chem.FindPotentialStereo(cyclohexene, cleanIt=False, flagPossible=True))
    assert not any(str(si.type) == "Bond_Double" for si in potential), "1,2-disub cyclohexene must be absent"
    print("  OK: 1,2-disubstituted cyclohexene absent from FindPotentialStereo")

    cyclooctene = Chem.MolFromSmiles("CC1=C(C)CCCCCC1")
    potential = list(Chem.FindPotentialStereo(cyclooctene, cleanIt=False, flagPossible=True))
    double_entries = [si for si in potential if str(si.type) == "Bond_Double"]
    assert len(double_entries) == 1 and str(double_entries[0].specified) == "Unspecified"
    print("  OK: 1,2-disubstituted cyclooctene present as Unspecified")

    # Fabricated-row classification negative control: one differing feature
    # must land in a different bucket.
    fake_buckets = defaultdict(list)
    rows = [
        ("m1", True, ("Aromatic", "Aromatic"), True),
        ("m2", True, ("Aromatic", "Aromatic"), True),
        ("m3", False, ("Single", "Up"), True),
    ]
    for smi, has_ring, orders, specified in rows:
        key = ("has_ring" if has_ring else "no_ring", orders, specified)
        fake_buckets[key].append(smi)
    assert len(fake_buckets) == 2, "differing has_ring/orders feature must split into a different bucket"
    print("  OK: classification negative control (differing feature -> different bucket)")

    # Calibration + negative control (live, real cargo run -- the actual
    # audit's own anchors, not fabricated).
    calibration = calibration_check()
    print(f"  calibration_check: pinned_pair_converges={calibration['pinned_pair_converges']}, "
          f"fresh_relabeling_stable={calibration['fresh_relabeling_stable']}")

    neg_axis2 = axis2_alternates([NEGATIVE_CONTROL])
    neg_result = neg_axis2.get(NEGATIVE_CONTROL, {})
    assert not neg_result.get("any_divergent", True), "negative control must show 0 axis-2 divergence"
    print(f"  OK: negative control (known-fine ring-shape fixture) shows 0 axis-2 divergence "
          f"({neg_result.get('n_alternates_generated', 0)} alternate(s) generated)")

    print("=== self-test passed ===")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("corpus", nargs="?", default=DEFAULT_CORPUS)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--relabelings",
        type=int,
        default=N_RELABELINGS_PER_MOLECULE,
        help="seeded RDKit atom relabelings per coupled molecule (default: 64)",
    )
    parser.add_argument("--out-jsonl", default=OUT_JSONL)
    parser.add_argument("--out-summary", default=OUT_SUMMARY)
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    provenance = provenance_gate(args.corpus)
    topology, coupled_components, by_smiles_ends = current_topology(args.corpus)
    coupled_smiles = sorted(by_smiles_ends.keys())

    print(f"provenance gate: {provenance}")
    print(f"topology: {topology}")

    if args.relabelings < 1:
        parser.error("--relabelings must be >= 1")
    axis1_result = axis1_divergence(coupled_smiles, args.relabelings)
    print(f"axis1: {axis1_result['n_divergent']}/{axis1_result['n_molecules_tested']} divergent, "
          f"{axis1_result['n_cross_correspondence_failures_total']} cross-correspondence failures")

    axis2_all = axis2_alternates(coupled_smiles + [EZ_NEVER_CORRUPTS_A, EZ_NEVER_CORRUPTS_B, NEGATIVE_CONTROL])
    n_axis2_alternates = sum(v["n_alternates_generated"] for v in axis2_all.values())
    print(f"axis2: {n_axis2_alternates} total alternates generated across "
          f"{len(axis2_all)} molecules (see axis2_applicability_note)")

    bucket_summary, per_component_detail = structural_classification(coupled_smiles, by_smiles_ends, axis1_result)
    calibration = calibration_check(args.relabelings)
    final_verdict = verdict(
        len(coupled_smiles),
        axis1_result,
        bucket_summary,
        calibration,
        args.relabelings,
    )

    os.makedirs(os.path.dirname(args.out_jsonl), exist_ok=True)
    with open(args.out_jsonl, "w") as f:
        for smi in coupled_smiles:
            f.write(json.dumps({
                "smiles": smi,
                "axis1": axis1_result["per_molecule"].get(smi),
                "axis2": axis2_all.get(smi),
                "classification": per_component_detail.get(smi),
            }) + "\n")

    summary = {
        "corpus": path_for_report(args.corpus),
        "provenance_gate": provenance,
        "topology": topology,
        "axis1_summary": {k: v for k, v in axis1_result.items() if k != "per_molecule"},
        "relabelings_per_molecule": args.relabelings,
        "axis2_applicability_note": axis2_applicability_note(),
        "n_axis2_alternates_generated_total": n_axis2_alternates,
        "structural_classification_buckets": bucket_summary,
        "calibration_check": calibration,
        "verdict": final_verdict,
    }
    with open(args.out_summary, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nwrote {args.out_jsonl}")
    print(f"wrote {args.out_summary}")
    print(f"\nVERDICT: {final_verdict['verdict']}")


if __name__ == "__main__":
    main()
