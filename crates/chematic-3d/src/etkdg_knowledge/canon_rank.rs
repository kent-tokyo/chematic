//! Minimal Morgan-style (extended connectivity) atom-rank invariant, used
//! only by `matcher.rs`'s [`super::matcher`]-internal
//! [`canonical_atoms`](super::matcher) tie-break for choosing a
//! deterministic outer-atom quadruple when multiple same-tier candidates
//! match one central bond (see that function's doc comment for the bug this
//! fixes).
//!
//! `chematic-smiles::canonical::morgan_ranks` already implements this exact
//! algorithm (same neighbor-hash-refinement-to-fixpoint idea), but
//! `chematic-3d` does not depend on `chematic-smiles` in production code
//! (only as a dev-dependency, per that crate's own Cargo.toml comment and
//! this workspace's crate-layering convention in the root `CLAUDE.md`) --
//! reimplemented locally rather than adding a new cross-crate dependency for
//! one function, using only `chematic_core` types this crate already
//! depends on. If a third consumer needs this, hoisting it to
//! `chematic-perception` (which both `chematic-smiles` and `chematic-3d`
//! already depend on) would be the principled de-duplication, not attempted
//! here (out of this PR's file-ownership scope, spec §14).
//!
//! ponytail: intentionally a smaller feature set than the original --
//! element/degree/charge/aromaticity/bond-order only, no isotope/wildcard/
//! explicit-H-count refinement (irrelevant to distinguishing a torsion
//! rule's outer-atom choice in practice). Upgrade by porting the extra
//! fields from `chematic-smiles::canonical::initial_invariant` if a case
//! ever needs them.

use chematic_core::{AtomIdx, BondOrder, Molecule};

fn initial_invariant(mol: &Molecule, idx: AtomIdx) -> u64 {
    let atom = mol.atom(idx);
    let an = atom.element.atomic_number() as u64;
    let degree = mol.degree(idx) as u64;
    let charge = (atom.charge as i64 + 128) as u64;
    let arom = atom.aromatic as u64;
    (an << 24) | (degree << 16) | (charge << 8) | arom
}

fn bond_order_value(order: BondOrder) -> u64 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Quadruple => 5,
        _ => 0,
    }
}

fn fnv_hash_sequence(base: u64, values: &[u64]) -> u64 {
    const FNV_PRIME: u64 = 0x0000_0100_0000_01B3;
    const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    let mut h = FNV_OFFSET ^ base.wrapping_mul(FNV_PRIME);
    for &v in values {
        h ^= v;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

fn normalize_ranks(ranks: &[u64]) -> Vec<u64> {
    let mut sorted: Vec<u64> = ranks.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    ranks
        .iter()
        .map(|r| sorted.partition_point(|&s| s < *r) as u64)
        .collect()
}

/// Topological (structure-only) per-atom rank: two atoms get the same rank
/// iff no graph feature (element/degree/charge/aromaticity/neighbor
/// environment) distinguishes them -- i.e. iff they are in the same
/// automorphism-orbit-or-larger refinement cell. Does NOT resolve remaining
/// ties beyond that (no individualize-refine); genuine automorphism orbits
/// keep equal ranks by design, same as `chematic-smiles`'s `morgan_ranks`.
pub(super) fn morgan_ranks(mol: &Molecule) -> Vec<u64> {
    let n = mol.atom_count();
    let mut ranks: Vec<u64> = (0..n)
        .map(|i| initial_invariant(mol, AtomIdx(i as u32)))
        .collect();

    let max_iter = n + 2;
    for _ in 0..max_iter {
        let old_distinct = {
            let mut s = ranks.clone();
            s.sort_unstable();
            s.dedup();
            s.len()
        };

        let new_ranks: Vec<u64> = (0..n)
            .map(|i| {
                let idx = AtomIdx(i as u32);
                let mut contributions: Vec<u64> = mol
                    .neighbors(idx)
                    .map(|(nb, bidx)| {
                        let bond_val = bond_order_value(mol.bond(bidx).order);
                        fnv_hash_sequence(ranks[nb.0 as usize], &[bond_val])
                    })
                    .collect();
                contributions.sort_unstable();
                fnv_hash_sequence(ranks[i], &contributions)
            })
            .collect();

        let new_distinct = {
            let mut s = new_ranks.clone();
            s.sort_unstable();
            s.dedup();
            s.len()
        };
        ranks = new_ranks;
        if new_distinct <= old_distinct {
            break;
        }
    }

    normalize_ranks(&ranks)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn benzene_all_carbons_tie() {
        let mol = parse("c1ccccc1").unwrap();
        let ranks = morgan_ranks(&mol);
        assert!(ranks.iter().all(|&r| r == ranks[0]), "{ranks:?}");
    }

    #[test]
    fn toluene_distinguishes_ipso_ortho_meta_para_and_methyl() {
        // Methyl-C, ipso-C, ortho x2, meta x2, para -- 5 distinct classes.
        let mol = parse("Cc1ccccc1").unwrap();
        let ranks = morgan_ranks(&mol);
        let distinct: std::collections::HashSet<u64> = ranks.iter().copied().collect();
        assert_eq!(distinct.len(), 5, "{ranks:?}");
    }

    #[test]
    fn acetic_acid_distinguishes_the_two_oxygens() {
        // CC(=O)O: carbonyl O (degree 1, double bond) vs hydroxyl O (degree
        // 1, single bond) must NOT tie.
        let mol = parse("CC(=O)O").unwrap();
        let ranks = morgan_ranks(&mol);
        // atoms: C0, C1, O2(=O), O3(-OH)
        assert_ne!(ranks[2], ranks[3], "{ranks:?}");
    }

    #[test]
    fn is_invariant_to_atom_relabeling() {
        // Same molecule, atoms written in reverse order -- the MULTISET of
        // ranks (not the per-index assignment, which is numbering-dependent
        // by construction) must match.
        let forward = parse("CC(=O)Nc1ccc(O)cc1").unwrap(); // paracetamol-ish
        let mut fwd_ranks = morgan_ranks(&forward);
        fwd_ranks.sort_unstable();

        let backward = parse("Oc1ccc(NC(C)=O)cc1").unwrap();
        let mut bwd_ranks = morgan_ranks(&backward);
        bwd_ranks.sort_unstable();

        assert_eq!(fwd_ranks, bwd_ranks);
    }

    /// Regression evidence for the residual atom-order-energy cases the
    /// gap-check example's `atom_order_energy_invariance` still reports as a
    /// disclosed (not unexplained) FAIL after the `uniquify: false` +
    /// `canonical_atoms` fix (see `matcher.rs`'s doc comment for the full
    /// writeup): biphenyl, adamantane, norbornane, and cubane. In each, the
    /// specific outer atoms the torsion match picks tie in Morgan rank --
    /// but **rank-tying is necessary, not sufficient, for the substitution
    /// to be a genuine graph automorphism**, and this test only checks the
    /// necessary half. An earlier draft of this doc claimed all 4 (then 3)
    /// cases were confirmed true automorphisms; independent review found
    /// that was wrong for cubane specifically, via a stricter, constrained
    /// check this test does NOT perform (see below) -- corrected here rather
    /// than left standing.
    ///
    /// - `biphenyl`, `adamantane`, `norbornane`: independently confirmed (via
    ///   `mol.GetSubstructMatches(mol, uniquify=False)`, constrained to hold
    ///   the central bond AND the other outer atom fixed -- the condition
    ///   that actually matters, not just "does *some* automorphism map atom
    ///   X to atom Y" unconstrained) to be genuine automorphisms: these
    ///   atoms truly are chemically interchangeable, and no purely-
    ///   topological rule can order-break them.
    /// - `cubane`: the SAME constrained check finds this is **NOT** a
    ///   genuine automorphism, even though atoms 2 and 5 really do tie in
    ///   `morgan_ranks` -- and that tie is *correct*, not a refinement
    ///   artefact: `morgan_ranks`'s stable partition matches cubane's real
    ///   automorphism-orbit partition exactly (independently verified by
    ///   enumerating `Aut(G)` directly). Atoms 2 and 5 genuinely are in one
    ///   global orbit. The problem is that global-orbit membership is the
    ///   wrong equivalence for this use: `canonical_atoms` (in `matcher.rs`)
    ///   needs equivalence under the *stabilizer of the central bond*, and
    ///   cubane's only non-trivial automorphism maps atom 2 to atom 5 while
    ///   also moving the central bond's own endpoints elsewhere -- so no
    ///   automorphism fixes the central bond and sends 2 to 5 at once. No
    ///   per-atom rank can distinguish this, however precisely it computes
    ///   global orbits; the ambiguity is a property of the quadruple as a
    ///   whole. This fixture's residual is a live, unresolved instance of the
    ///   same tie-break bug fixed elsewhere in this pass (menthol/
    ///   testosterone/cholesterol/penicillin_core), disclosed rather than
    ///   mislabeled. See `matcher.rs`'s `canonical_atoms` doc comment for the
    ///   full mechanism.
    ///
    /// The rank-tie assertions below remain correct and useful regardless:
    /// they are a true, falsifiable fact about `morgan_ranks`'s own output
    /// (falsifiable if a future change accidentally starts or stops
    /// distinguishing these atoms), just not, on their own, proof of a real
    /// automorphism -- that stronger claim needs the constrained check
    /// above, done once via an external oracle, not re-derived by this test.
    ///
    /// penicillin_core's gem-dimethyl pair (asserted below too, also a
    /// genuine automorphism) ties in rank but is NOT one of the 4 named
    /// residuals: `atom_in_ring_size_range` (round-4 fix, `matcher.rs`)
    /// constrains a ring-torsion rule's outer atoms to ring members, and the
    /// gem-dimethyl carbons aren't in any ring, so neither ever gets picked
    /// as the outer atom in the first place -- the tie exists in the graph
    /// but the matcher never has to break it. Kept here as a standing
    /// regression check on `morgan_ranks` itself, independent of whether the
    /// matcher currently reaches this tie.
    #[test]
    fn known_rank_ties_behind_the_residual_atom_order_cases() {
        // biphenyl: the two ortho carbons flanking each ipso carbon (a
        // genuine local mirror symmetry) tie -- atoms 2 and 10 both flank
        // atom 3 (excluding the central atom 4). Confirmed a genuine
        // constrained automorphism (see doc comment above).
        let biphenyl = parse("c1ccc(-c2ccccc2)cc1").unwrap();
        let r = morgan_ranks(&biphenyl);
        assert_eq!(r[2], r[10], "biphenyl ortho pair must tie: {r:?}");

        // adamantane: atom 8's two non-atom-10 neighbors (atoms 7 and 9,
        // both bridgehead-adjacent methylenes) tie. Confirmed a genuine
        // constrained automorphism (see doc comment above).
        let adamantane = parse("C1CC2CC3CC1CC(C2)C3").unwrap();
        let r = morgan_ranks(&adamantane);
        assert_eq!(
            r[7], r[9],
            "adamantane bridge-neighbor pair must tie: {r:?}"
        );

        // norbornane: the two real candidate substitutions independent
        // review's methodology traces (atoms 0/4 and atoms 1/3) both tie.
        // Confirmed a genuine constrained automorphism (see doc comment
        // above) -- this fixture's small residual is an embedding-geometry
        // artifact, not a tie-break defect.
        let norbornane = parse("C1CC2CCC1C2").unwrap();
        let r = morgan_ranks(&norbornane);
        assert_eq!(r[0], r[4], "norbornane bridgehead pair must tie: {r:?}");
        assert_eq!(r[1], r[3], "norbornane bridgehead pair must tie: {r:?}");

        // cubane: atom 1's two non-atom-0 neighbors (atoms 2 and 5) tie in
        // rank, and the two atoms genuinely ARE chemically interchangeable
        // (a real global automorphism swaps them) -- but that is the wrong
        // equivalence for `canonical_atoms` to use here: the constrained
        // check in the doc comment above confirms no automorphism realizes
        // that swap while also fixing the central bond, which is the
        // relation that actually matters for picking one quadruple. Kept as
        // a standing regression check that `morgan_ranks` still produces
        // this specific (correct, but not bond-stabilizer-aware) tie, not
        // evidence the resulting quadruple choice is harmless.
        let cubane = parse("C1C2C3C1C4C2C3C4").unwrap();
        let r = morgan_ranks(&cubane);
        assert_eq!(r[2], r[5], "cubane symmetric-neighbor pair must tie: {r:?}");

        // penicillin_core: the gem-dimethyl group's two methyl carbons
        // (atoms 0 and 2, both bonded only to the same quaternary ring
        // carbon 1) tie -- the one non-trivial tie in an otherwise fully
        // discrete-ranked molecule. Confirmed a genuine automorphism, but
        // never reached by the matcher (see doc comment above).
        let penicillin = parse("CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O").unwrap();
        let r = morgan_ranks(&penicillin);
        assert_eq!(r[0], r[2], "penicillin gem-dimethyl pair must tie: {r:?}");
    }
}
