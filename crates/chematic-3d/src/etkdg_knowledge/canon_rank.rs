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

    /// Regression evidence for the residual 4/72 corpus fixtures the
    /// gap-check example's `atom_order_energy_invariance` still reports as
    /// FAIL after the `uniquify: false` + `canonical_atoms` fix (see
    /// `matcher.rs`'s doc comment): in every one of these 4, the specific
    /// outer atoms the torsion match picks tie in Morgan rank -- i.e. they
    /// are GENUINELY, chemically interchangeable substituents (a true graph
    /// automorphism, not two different atoms an insufficiently-refined
    /// invariant merely failed to distinguish). No purely-topological rule
    /// can pick one over the other in a way that survives adversarial
    /// relabeling, since by definition nothing about the graph distinguishes
    /// them -- confirmed here directly rather than assumed, so this is
    /// falsifiable if a future change to `morgan_ranks` accidentally starts
    /// distinguishing (or stops distinguishing) the wrong atoms.
    #[test]
    fn known_true_symmetry_ties_behind_the_residual_atom_order_cases() {
        // biphenyl: the two ortho carbons flanking each ipso carbon (a
        // genuine local mirror symmetry) tie -- atoms 2 and 10 both flank
        // atom 3 (excluding the central atom 4).
        let biphenyl = parse("c1ccc(-c2ccccc2)cc1").unwrap();
        let r = morgan_ranks(&biphenyl);
        assert_eq!(r[2], r[10], "biphenyl ortho pair must tie: {r:?}");

        // adamantane: atom 8's two non-atom-10 neighbors (atoms 7 and 9,
        // both bridgehead-adjacent methylenes) tie.
        let adamantane = parse("C1CC2CC3CC1CC(C2)C3").unwrap();
        let r = morgan_ranks(&adamantane);
        assert_eq!(
            r[7], r[9],
            "adamantane bridge-neighbor pair must tie: {r:?}"
        );

        // cubane: atom 1's two non-atom-0 neighbors (atoms 2 and 5, related
        // by the cube's own symmetry) tie.
        let cubane = parse("C1C2C3C1C4C2C3C4").unwrap();
        let r = morgan_ranks(&cubane);
        assert_eq!(r[2], r[5], "cubane symmetric-neighbor pair must tie: {r:?}");

        // penicillin_core: the gem-dimethyl group's two methyl carbons
        // (atoms 0 and 2, both bonded only to the same quaternary ring
        // carbon 1) tie -- the one non-trivial tie in an otherwise fully
        // discrete-ranked molecule.
        let penicillin = parse("CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O").unwrap();
        let r = morgan_ranks(&penicillin);
        assert_eq!(r[0], r[2], "penicillin gem-dimethyl pair must tie: {r:?}");
    }
}
