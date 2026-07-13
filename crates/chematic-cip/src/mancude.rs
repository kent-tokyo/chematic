//! A bounded, brute-force MANCUDE (maximum non-cumulated double bonds) oracle.
//!
//! `chematic_core::kekulization::kekulize` returns *one* canonical Kekulé structure per
//! molecule (a single maximum matching). IUPAC's actual CIP treatment of a mancude ring
//! system needs every valid perfect matching, since a ring atom's duplicate contribution
//! is the *mean* atomic number of whichever neighbor it's double-bonded to, averaged
//! across all of them -- not whichever one a single matching happens to pick.
//! [`enumerate_kekule_matchings`] fills that gap by exhaustive backtracking, reusing the
//! exact same must-match/lone-pair-donor classification `kekulize()` itself uses (via
//! [`chematic_core::kekulization::atom_must_be_matched`], widened to `pub` for this
//! purpose) so the two can never silently disagree about what counts as a valid
//! placement.
//!
//! **Deliberately not a production algorithm.** Full enumeration of every perfect
//! matching is combinatorially explosive on large fused ring systems; [`MancudeBudget`]
//! bounds it and *errors* rather than silently truncating when exceeded (same discipline
//! as [`crate::budget::CipBudget`]). This module exists as a small-ring test oracle to
//! design and verify Milestone 3B-1's production representation against -- M3B-2 will
//! need a counting/matching-based approach for scale, not this brute force.
//!
//! **Design note for Milestone 3B-1 (recorded here, not resolved this round)**: this
//! oracle keys off `BondOrder::Aromatic` on bonds, exactly like `kekulize()` does. A
//! `Molecule` produced by `apply_kekule` has `Single`/`Double` bonds and *no* aromatic
//! bonds left (only atoms keep their `aromatic` flag) -- so this oracle can only compute
//! a signature from aromatic-notation input, never from an already-Kekulé-respelled one.
//! That's fine for this milestone's own tests, but it means M3B-1/2's production digraph
//! construction must detect a MANCUDE component from **atom aromatic flags**, not bond
//! order, or Kekulé-respelled input would silently produce no `MancudeDuplicate` at all.

use std::collections::HashMap;

use chematic_core::kekulization::{KekuleResult, atom_must_be_matched, build_kekule_result};
use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

use crate::rational::RationalAtomicNumber;

/// Bounds on [`enumerate_kekule_matchings`]'s search, to keep it a small-ring oracle
/// rather than an accidental production path. Both bounds *error* the whole call rather
/// than returning a truncated/partial result -- an incomplete enumeration would silently
/// corrupt any mean computed from it.
#[derive(Debug, Clone, Copy)]
pub struct MancudeBudget {
    /// Maximum number of must-match atoms in one connected component.
    pub max_atoms: usize,
    /// Maximum number of complete matchings collected.
    pub max_matchings: usize,
    /// Maximum number of backtracking search steps (guards pathological search trees
    /// that explore heavily before finding few, or zero, complete matchings).
    pub max_search_steps: usize,
}

impl Default for MancudeBudget {
    fn default() -> Self {
        Self {
            max_atoms: 24,
            max_matchings: 64,
            max_search_steps: 100_000,
        }
    }
}

/// Errors from [`enumerate_kekule_matchings`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MancudeError {
    /// The must-match component is larger than [`MancudeBudget::max_atoms`].
    TooManyAtoms { count: usize, max: usize },
    /// More valid matchings exist than [`MancudeBudget::max_matchings`] allows.
    TooManyMatchings { max: usize },
    /// The backtracking search exceeded [`MancudeBudget::max_search_steps`].
    SearchBudgetExceeded { max: usize },
}

impl core::fmt::Display for MancudeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            MancudeError::TooManyAtoms { count, max } => {
                write!(
                    f,
                    "mancude component has {count} must-match atoms, over budget {max}"
                )
            }
            MancudeError::TooManyMatchings { max } => {
                write!(
                    f,
                    "mancude component has more than {max} valid Kekulé matchings"
                )
            }
            MancudeError::SearchBudgetExceeded { max } => {
                write!(f, "mancude matching search exceeded {max} steps")
            }
        }
    }
}

impl std::error::Error for MancudeError {}

/// Enumerate *every* valid perfect matching of `mol`'s aromatic must-match subgraph,
/// bounded by `budget`. Returns one [`KekuleResult`] per valid matching, in the exact
/// same format `chematic_core::kekulization::kekulize` produces (every aromatic bond
/// mapped to `Single` or `Double`) -- so `kekulize(mol)`'s own single result is always
/// expected to appear as a member of this function's output (checked in this module's
/// own tests).
///
/// A molecule with no aromatic bonds has exactly one (empty) valid placement, mirroring
/// `kekulize()`'s own no-op convention.
pub fn enumerate_kekule_matchings(
    mol: &Molecule,
    budget: MancudeBudget,
) -> Result<Vec<KekuleResult>, MancudeError> {
    let mut aromatic_bonds: Vec<BondIdx> = Vec::new();
    let mut aromatic_atoms: Vec<AtomIdx> = Vec::new();
    for (bidx, bond) in mol.bonds() {
        if bond.order == BondOrder::Aromatic {
            aromatic_bonds.push(bidx);
            if !aromatic_atoms.contains(&bond.atom1) {
                aromatic_atoms.push(bond.atom1);
            }
            if !aromatic_atoms.contains(&bond.atom2) {
                aromatic_atoms.push(bond.atom2);
            }
        }
    }
    if aromatic_bonds.is_empty() {
        return Ok(vec![KekuleResult::new()]);
    }

    let mut must_match: Vec<AtomIdx> = aromatic_atoms
        .into_iter()
        .filter(|&idx| atom_must_be_matched(mol, idx))
        .collect();
    must_match.sort();

    if must_match.len() > budget.max_atoms {
        return Err(MancudeError::TooManyAtoms {
            count: must_match.len(),
            max: budget.max_atoms,
        });
    }

    let mut adj: HashMap<AtomIdx, Vec<AtomIdx>> = HashMap::new();
    for &bidx in &aromatic_bonds {
        let bond = mol.bond(bidx);
        if must_match.contains(&bond.atom1) && must_match.contains(&bond.atom2) {
            adj.entry(bond.atom1).or_default().push(bond.atom2);
            adj.entry(bond.atom2).or_default().push(bond.atom1);
        }
    }

    let mut results: Vec<HashMap<AtomIdx, AtomIdx>> = Vec::new();
    let mut current: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let mut steps = 0usize;
    backtrack(
        &must_match,
        &adj,
        &mut current,
        &mut results,
        &budget,
        &mut steps,
    )?;

    Ok(results
        .into_iter()
        .map(|matching| build_kekule_result(&aromatic_bonds, mol, &matching))
        .collect())
}

/// Enumerate every perfect matching of `remaining` under `adj` by always extending the
/// smallest still-unmatched atom -- this visits each complete matching exactly once (the
/// atom's partner in any given valid matching is unique, so branching over its candidate
/// partners partitions the search space without overlap or omission).
fn backtrack(
    remaining: &[AtomIdx],
    adj: &HashMap<AtomIdx, Vec<AtomIdx>>,
    current: &mut HashMap<AtomIdx, AtomIdx>,
    results: &mut Vec<HashMap<AtomIdx, AtomIdx>>,
    budget: &MancudeBudget,
    steps: &mut usize,
) -> Result<(), MancudeError> {
    *steps += 1;
    if *steps > budget.max_search_steps {
        return Err(MancudeError::SearchBudgetExceeded {
            max: budget.max_search_steps,
        });
    }

    let still_free: Vec<AtomIdx> = remaining
        .iter()
        .copied()
        .filter(|a| !current.contains_key(a))
        .collect();
    let Some(&v) = still_free.first() else {
        if results.len() >= budget.max_matchings {
            return Err(MancudeError::TooManyMatchings {
                max: budget.max_matchings,
            });
        }
        results.push(current.clone());
        return Ok(());
    };

    let candidates: Vec<AtomIdx> = adj
        .get(&v)
        .into_iter()
        .flatten()
        .copied()
        .filter(|u| !current.contains_key(u))
        .collect();
    for u in candidates {
        current.insert(v, u);
        current.insert(u, v);
        backtrack(remaining, adj, current, results, budget, steps)?;
        current.remove(&v);
        current.remove(&u);
    }
    Ok(())
}

/// For one atom, the mean atomic number of whichever neighbor it is double-bonded to,
/// across every matching in `matchings` -- the quantity Milestone 3B-1's mancude
/// duplicate node will store. `None` if `atom` is never double-bonded in any of the
/// matchings (a lone-pair donor, e.g. furan's O or pyrrole's `[nH]`, or an atom outside
/// any mancude system) -- such atoms contribute no extra duplicate at all, matching
/// `kekulize()`'s existing single/double-only bond model.
pub fn effective_atomic_number(
    mol: &Molecule,
    atom: AtomIdx,
    matchings: &[KekuleResult],
) -> Option<RationalAtomicNumber> {
    let mut partners = Vec::new();
    for matching in matchings {
        for (nb, bidx) in mol.neighbors(atom) {
            if matching.get(&bidx) == Some(&BondOrder::Double) {
                partners.push(mol.atom(nb).element.atomic_number() as u32);
                break;
            }
        }
    }
    if partners.is_empty() {
        None
    } else {
        Some(RationalAtomicNumber::mean(&partners))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{Atom, Element, MoleculeBuilder};

    fn benzene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// N1=C2-C3=C4-C5=C6(-N1) around the ring; C2 and C6 are the carbons adjacent to N.
    fn pyridine() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(Atom::aromatic(Element::N));
        let cs: Vec<_> = (0..5)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        let atoms = [n, cs[0], cs[1], cs[2], cs[3], cs[4]];
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    fn pyridine_n_and_adjacent_carbon() -> (Molecule, AtomIdx, AtomIdx) {
        let mol = pyridine();
        // atom 0 = N, atom 1 = the carbon adjacent to it (per construction above).
        (mol, AtomIdx(0), AtomIdx(1))
    }

    /// O-C1=C2-C3=C4(-O) around the ring; O is a lone-pair donor, excluded from must_match.
    fn furan() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::aromatic(Element::O));
        let cs: Vec<_> = (0..4)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        let atoms = [o, cs[0], cs[1], cs[2], cs[3]];
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    #[test]
    fn no_aromatic_bonds_is_one_empty_matching() {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mol = b.build();
        let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
        assert_eq!(all, vec![KekuleResult::new()]);
    }

    /// Consistency check: kekulize()'s own single result must be a member of the full
    /// enumeration, for every fixture -- if it weren't, the two would be using different
    /// notions of "valid matching" and the oracle couldn't be trusted to design against.
    #[test]
    fn kekulize_result_is_a_member_of_the_full_enumeration() {
        for mol in [benzene(), pyridine(), furan()] {
            let single = chematic_core::kekulization::kekulize(&mol).unwrap();
            let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
            assert!(
                all.contains(&single),
                "kekulize()'s own result must appear in the full enumeration"
            );
        }
    }

    #[test]
    fn benzene_has_exactly_two_kekule_forms() {
        let all = enumerate_kekule_matchings(&benzene(), MancudeBudget::default()).unwrap();
        assert_eq!(all.len(), 2);
    }

    #[test]
    fn furan_has_exactly_one_kekule_form() {
        // O is a lone-pair donor (excluded from must-match); the remaining 4 ring
        // carbons form a PATH (not a cycle, since neither C-O bond is a matching
        // candidate), and a 4-atom path has exactly one perfect matching: {C1-C2, C3-C4}.
        let all = enumerate_kekule_matchings(&furan(), MancudeBudget::default()).unwrap();
        assert_eq!(all.len(), 1);
    }

    /// Hand-derived IUPAC example, matching the design conversation's own "6½" value:
    /// pyridine's ring carbon adjacent to N is double-bonded to N in one Kekulé form and
    /// to its other (carbon) ring neighbor in the other -- effective atomic number
    /// (7 + 6) / 2 = 13/2 = 6½.
    #[test]
    fn pyridine_adjacent_carbon_is_six_and_a_half() {
        let (mol, _n, adjacent_c) = pyridine_n_and_adjacent_carbon();
        let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
        assert_eq!(all.len(), 2);
        let signature = effective_atomic_number(&mol, adjacent_c, &all).unwrap();
        assert_eq!(signature.numerator(), 13);
        assert_eq!(signature.denominator(), 2);
    }

    /// Furan's O is never double-bonded (lone-pair donor) -- no duplicate contribution.
    #[test]
    fn furan_oxygen_has_no_effective_atomic_number() {
        let mol = furan();
        let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
        assert_eq!(effective_atomic_number(&mol, AtomIdx(0), &all), None);
    }

    /// The common-signature property M3B-0 exists to demonstrate: two individually-valid,
    /// genuinely different Kekulé forms of the same molecule disagree on a given atom's
    /// immediate double-bond partner, yet both are members of the one enumeration whose
    /// mean is the single MANCUDE signature -- checked on a hydrocarbon fixture (trivial,
    /// can't fail: both forms agree since every ring atom is carbon) and a hetero fixture
    /// (genuinely fractional: the two forms disagree, and the signature averages them).
    #[test]
    fn kekule_form_a_and_b_share_one_common_signature_hydrocarbon() {
        let mol = benzene();
        let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
        let (form_a, form_b) = (&all[0], &all[1]);
        assert_ne!(
            form_a, form_b,
            "must be two genuinely different resonance structures"
        );

        let atom = AtomIdx(0);
        let partner_in = |form: &KekuleResult| {
            mol.neighbors(atom)
                .find(|&(_, bidx)| form.get(&bidx) == Some(&BondOrder::Double))
                .map(|(nb, _)| mol.atom(nb).element.atomic_number() as u32)
                .unwrap()
        };
        let (pa, pb) = (partner_in(form_a), partner_in(form_b));
        // Hydrocarbon: both partners are carbon, so the two forms happen to agree here --
        // that's expected, not a bug; the fraction only becomes visible with a heteroatom
        // (see the hetero variant of this test below).
        assert_eq!(pa, 6);
        assert_eq!(pb, 6);

        let signature = effective_atomic_number(&mol, atom, &all).unwrap();
        assert_eq!(signature, RationalAtomicNumber::mean(&[pa, pb]));
        assert_eq!(signature, RationalAtomicNumber::integer(6));
    }

    #[test]
    fn kekule_form_a_and_b_share_one_common_signature_hetero() {
        let (mol, _n, adjacent_c) = pyridine_n_and_adjacent_carbon();
        let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
        let (form_a, form_b) = (&all[0], &all[1]);
        assert_ne!(
            form_a, form_b,
            "must be two genuinely different resonance structures"
        );

        let partner_in = |form: &KekuleResult| {
            mol.neighbors(adjacent_c)
                .find(|&(_, bidx)| form.get(&bidx) == Some(&BondOrder::Double))
                .map(|(nb, _)| mol.atom(nb).element.atomic_number() as u32)
                .unwrap()
        };
        let (pa, pb) = (partner_in(form_a), partner_in(form_b));
        // The two forms genuinely disagree (one has this carbon double-bonded to N, the
        // other to its carbon neighbor) -- this divergence is exactly what a single-form
        // representation (today's digraph) can't average away, and what the common
        // signature below reconciles into one value.
        assert_ne!(
            pa, pb,
            "hetero fixture must show the two forms actually disagreeing"
        );

        let signature = effective_atomic_number(&mol, adjacent_c, &all).unwrap();
        assert_eq!(signature, RationalAtomicNumber::mean(&[pa, pb]));
        assert_eq!(signature.numerator(), 13);
        assert_eq!(signature.denominator(), 2);
    }
}
