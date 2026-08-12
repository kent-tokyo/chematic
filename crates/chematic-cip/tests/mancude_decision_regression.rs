//! MANCUDE-Decision-A0 regression fixtures.
//!
//! Freezes the 3 stereocenters found (full-corpus diagnostic,
//! `examples/mancude_decision_diagnosis.rs`, `SMILES.csv` sha256 `1c47371d...`) whose
//! resolved R/S label differs between the un-Kekulized baseline
//! (`assign_cip_accurate_experimental_without_mancude` on the original, aromatic-form
//! molecule) and the live production path -- see `compare.rs`'s
//! `CompareContext::fractional_decisions` doc comment and
//! `docs/rfcs/cip_accurate_rfc.md`'s MANCUDE-Decision-A0 tripwire closeout entry. All three
//! were confirmed to match both RDKit's modern (`rdCIPLabeler`) and legacy
//! (`_CIPCode`) oracles at diagnosis time (`.venv/bin/python3`, RDKit 2026.03.3) --
//! frozen here as fixture expectations, not re-verified live against RDKit on every
//! test run.
//!
//! **Classification, isolating structure from fraction (the naive with/without-
//! `MancudeContext` contrast conflates both -- see the doc comment above for how the
//! first attempt at this got it wrong):** holding the Kekule-respelled structure fixed
//! and only toggling the attached `MancudeContext` leaves the root-child partition
//! (and hence the resolved label) **unchanged** for all 3 -- the fix from the baseline's
//! `R` to the correct `S` comes entirely from Kekule-respelling (real double-bond
//! duplicate nodes an aromatic-bond digraph never has), not from the fractional atomic
//! numbers. These are classification **D** ("fraction locally load-bearing --
//! `fractional_decisions > 0` -- but final-label-inert"), the same bucket as the other
//! 33 affected stereocenters this diagnostic found; **zero** stereocenters in this
//! corpus are classification E ("fraction changes the final label").
//!
//! Deliberately does **not** assert that
//! `assign_cip_accurate_experimental_without_mancude` gives `R` on the *original*
//! (un-Kekulized) molecule for these atoms: today it does (that's the whole finding),
//! but locking that in as an expected value would make a future improvement to the
//! plain, un-Kekulized path that also reaches `S` *fail* this suite for getting more
//! correct. Only the production path's correctness is frozen here, plus the structural
//! facts (a load-bearing fractional decision occurs somewhere in Pass 1's recursive
//! comparison, and the resolved label is structure-driven, not fraction-driven) that
//! explain why `fractional_decisions > 0` here is expected and benign.

use chematic_cip::digraph_diff::renumber_molecule;
use chematic_cip::{
    CipBudget, CipDigraph, CompareContext, NodeId, assign_cip_accurate_experimental,
    assign_cip_accurate_experimental_without_mancude, prepare_kekule_form, rank_children,
};
use chematic_core::{AtomIdx, Chirality, CipCode, Molecule};

struct Fixture {
    smiles: &'static str,
    atom_idx: u32,
    expected: CipCode,
}

const FIXTURES: &[Fixture] = &[
    Fixture {
        smiles: "Oc1cc(C(F)(F)F)c2cc3c(cc2n1)NC[C@@H]1CCCC[C@H]31",
        atom_idx: 22,
        expected: CipCode::S,
    },
    Fixture {
        smiles: "Oc1cc(C(F)(F)F)c2cc3c(cc2n1)NC[C@@H]1CCC[C@H]31",
        atom_idx: 21,
        expected: CipCode::S,
    },
    Fixture {
        smiles: "Oc1cc(C(F)(F)F)c2cc3c(cc2n1)NC[C@H]1CCCC[C@H]31",
        atom_idx: 22,
        expected: CipCode::S,
    },
];

#[test]
fn mancude_decision_a0_fixtures() {
    for f in FIXTURES {
        check_fixture(f);
    }
}

#[track_caller]
fn check_fixture(f: &Fixture) {
    let budget = CipBudget::default_budget();
    let mol = chematic_smiles::parse(f.smiles).expect("valid SMILES");
    let idx = AtomIdx(f.atom_idx);

    // 1. The live, production (mancude-aware) engine gives the expected, RDKit-agreeing
    // label.
    assert_eq!(
        live_code(&mol, idx, budget),
        Some(f.expected),
        "{}: live engine assignment at atom {}",
        f.smiles,
        f.atom_idx
    );

    // 2. A MANCUDE fraction is load-bearing *somewhere* in this center's Pass-1
    // comparison (fractional_decisions > 0) -- confirming this fixture actually
    // exercises the fractional-key path, not dead code.
    let (kekule_mol, ctx) = prepare_kekule_form(&mol).expect("kekulizable fixture");
    let fractional_decisions = fractional_decisions_at_root(&kekule_mol, idx, &ctx, budget);
    assert!(
        fractional_decisions > 0,
        "{}: expected a load-bearing fractional decision at atom {}",
        f.smiles,
        f.atom_idx
    );

    // 3. Classification D, not E: holding the Kekule-respelled *structure* fixed and
    // only toggling the attached MancudeContext leaves the root-child partition (and
    // therefore the resolved label) unchanged -- the fraction is locally load-bearing
    // (item 2) but never decides this center's final ranking. Uses the public
    // `_without_mancude` entry point directly on `kekule_mol` (not the original
    // aromatic `mol`) for the no-fraction side, since `apply_kekule` preserves
    // `AtomIdx` (assign.rs's own note) -- this is the properly isolated
    // "integer-collapsed control", not a with/without-mancude-on-different-structures
    // contrast.
    let structure_only_label =
        assign_cip_accurate_experimental_without_mancude(&kekule_mol, budget)
            .expect("budget not exceeded")
            .assignments
            .iter()
            .find(|(a, _)| *a == idx)
            .map(|(_, c)| *c);
    assert_eq!(
        structure_only_label,
        Some(f.expected),
        "{}: expected the Kekule-respelled structure ALONE (no MancudeContext) to \
         already resolve atom {} to {:?} -- if this fails, the fraction has become \
         load-bearing for the final label and this fixture is actually \
         classification E, not D; re-classify before trusting this test",
        f.smiles,
        f.atom_idx,
        f.expected
    );

    // 4. Atom-renumbering invariance.
    let perm: Vec<usize> = (0..mol.atom_count()).rev().collect();
    let (renumbered, old_to_new) = renumber_molecule(&mol, &perm);
    let new_idx = AtomIdx(old_to_new[f.atom_idx as usize]);
    assert_eq!(
        live_code(&renumbered, new_idx, budget),
        Some(f.expected),
        "{}: renumbered assignment (reversed atom order)",
        f.smiles
    );

    // 5. Canonical-SMILES round-trip invariance.
    let canonical = chematic_smiles::canonical_smiles(&mol);
    let reparsed = chematic_smiles::parse(&canonical).expect("canonical SMILES reparses");
    let target = locate_by_branch_signature(&mol, idx, &reparsed, budget).unwrap_or_else(|| {
        panic!(
            "{}: no unambiguous atom match after canonical round-trip (-> {canonical})",
            f.smiles
        )
    });
    assert_eq!(
        live_code(&reparsed, target, budget),
        Some(f.expected),
        "{}: canonical round-trip assignment (-> {canonical})",
        f.smiles
    );
}

fn live_code(mol: &Molecule, idx: AtomIdx, budget: CipBudget) -> Option<CipCode> {
    let result = assign_cip_accurate_experimental(mol, budget)
        .expect("budget not exceeded for these fixtures");
    result
        .assignments
        .iter()
        .find(|(a, _)| *a == idx)
        .map(|(_, c)| *c)
}

/// `fractional_decisions` for one stereocenter's root-children ranking, on the given
/// (already Kekule-respelled) molecule and its `MancudeContext`.
fn fractional_decisions_at_root(
    kekule_mol: &Molecule,
    idx: AtomIdx,
    ctx: &chematic_cip::MancudeContext,
    budget: CipBudget,
) -> u64 {
    let mut graph =
        CipDigraph::new_with_mancude(kekule_mol, idx, budget, ctx).expect("digraph builds");
    let root = graph.root();
    let children = graph.expand_children(root).expect("root expands");
    let mut cmp_ctx = CompareContext::new();
    let _groups: Vec<Vec<NodeId>> =
        rank_children(&mut graph, &children, &mut cmp_ctx).expect("ranking succeeds");
    cmp_ctx.fractional_decisions
}

/// Locate "the same" stereocenter in `mol2` via [`CipDigraph::branch_signature`] (a
/// chirality-independent structural subtree hash -- see its own doc comment). These
/// fixtures each have 2 structurally-distinct-but-locally-similar stereocenters, so
/// `digraph_diff::find_atom_by_signature`'s coarser (element/degree/chirality-tagged)
/// signature is ambiguous here (verified empirically); `branch_signature` disambiguates
/// uniquely on all 3 fixtures.
fn locate_by_branch_signature(
    mol1: &Molecule,
    atom1: AtomIdx,
    mol2: &Molecule,
    budget: CipBudget,
) -> Option<AtomIdx> {
    let sig = |m: &Molecule, a: AtomIdx| -> Option<u64> {
        let mut g = CipDigraph::new(m, a, budget).ok()?;
        let root = g.root();
        g.branch_signature(root).ok()
    };
    let want = sig(mol1, atom1)?;
    let mut matches = (0..mol2.atom_count())
        .map(|i| AtomIdx(i as u32))
        .filter(|&idx| mol2.atom(idx).chirality != Chirality::None)
        .filter(|&idx| sig(mol2, idx) == Some(want));
    let first = matches.next()?;
    if matches.next().is_some() {
        None
    } else {
        Some(first)
    }
}
