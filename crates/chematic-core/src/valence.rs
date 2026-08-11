//! Valence model: implicit hydrogen count for organic-subset atoms.
//!
//! Reference: OpenSMILES specification, section 3.4 (Implicit hydrogen)
//! <http://opensmiles.org/opensmiles-spec.html>

use crate::bond::BondOrder;
use crate::molecule::{AtomIdx, Molecule};
use std::fmt;

/// Compute the implicit hydrogen count for atom `idx`.
///
/// - Bracket atoms (`hydrogen_count.is_some()`): return the stored value directly.
/// - Wildcard atoms: return 0.
/// - Organic-subset atoms: derive from the normal-valence table.
/// - All other atoms: return 0 (no implicit H rule defined).
///
/// # Algorithm
/// 1. Sum the integer bond orders of all bonds on the atom.
/// 2. Find the smallest normal valence >= bond_sum (adjusted for formal charge).
/// 3. implicit_H = adjusted_valence - bond_sum.
///
/// Charge adjustment:
/// - Positive charge: increases target valence (e.g. [NH4]+ has 4 bonds → valence 4).
/// - Negative charge: decreases target valence.
pub fn implicit_hcount(mol: &Molecule, idx: AtomIdx) -> u8 {
    let atom = mol.atom(idx);

    // Wildcards have no defined implicit H.
    if atom.wildcard {
        return 0;
    }

    // Bracket atoms store the explicit H count.
    if let Some(h) = atom.hydrogen_count {
        return h;
    }

    valence_inferred_hcount(mol, idx)
}

/// Compute the hydrogen count that organic-subset (unbracketed) valence
/// inference would give for atom `idx`, **ignoring** any stored explicit
/// `hydrogen_count` -- unlike [`implicit_hcount`], which returns the stored
/// value directly for bracket atoms.
///
/// Used to decide whether a bracket atom's explicit H count is genuinely
/// disambiguating information (differs from what organic-subset spelling
/// would infer) or merely repeats it, in which case the atom can be
/// canonically re-spelled without brackets for that reason -- see
/// `chematic-smiles`'s `emit_atom`/`initial_invariant`, which both need
/// "would this atom's H count survive being unbracketed" independent of
/// whatever notation the atom happened to be parsed from.
pub fn valence_inferred_hcount(mol: &Molecule, idx: AtomIdx) -> u8 {
    let atom = mol.atom(idx);

    if atom.wildcard {
        return 0;
    }

    // Only the organic subset gets implicit H.
    if !atom.element.is_organic_subset() {
        return 0;
    }

    let normal_valences = atom.element.normal_valences();
    if normal_valences.is_empty() {
        return 0;
    }

    let charge = atom.charge as i32;

    // Separate aromatic bonds from non-aromatic bonds.
    let mut aromatic_count: usize = 0;
    let mut non_aromatic_sum: i32 = 0;
    for (_, bidx) in mol.neighbors(idx) {
        let bond = mol.bond(bidx);
        let order = bond.order;
        if order == BondOrder::Aromatic {
            aromatic_count += 1;
        } else if order == BondOrder::Dative && bond.atom1 == idx {
            // Donor side of a dative (coordinate) bond: per `BondOrder::Dative`'s
            // own documented `atom1 (donor) -> atom2 (acceptor)` convention, the
            // donor shares a lone pair without spending any of its own normal
            // covalent valence -- e.g. `N->[Pt]` must still imply NH3 (3 implicit
            // H, valence fully intact), not NH2. Matches RDKit's identical
            // treatment of the same SMILES dative-arrow syntax. Contributes 0,
            // unlike a plain covalent bond of the same `order_int()==1`.
            // The acceptor side (and any other atom for which `idx` is `atom2`)
            // is intentionally left counted as before -- out of scope here,
            // since every acceptor in the motivating corpus (Pt, Fe, Co, ...)
            // is already outside the organic subset and short-circuits above.
        } else {
            non_aromatic_sum += order.order_int() as i32;
        }
    }

    if aromatic_count > 0 {
        // Aromatic molecule (pre-Kekulization): each aromatic bond contributes 1.5
        // to the effective bond order (OpenSMILES convention).
        //
        // floor(1.5 × n) gives the contribution from n aromatic bonds:
        //   n=2 → 3  benzene CH:   4−3=1H ✓   pyridine N: 3−3=0H ✓
        //   n=3 → 4  junction C:   4−4=0H ✓
        //
        // Combined with non-aromatic substituents (e.g. N−CH₃) this correctly yields
        // 0 H for all substituted aromatic atoms without needing Kekulization.
        // Always use the lowest normal valence; aromatic atoms cannot be hypervalent.
        let aromatic_contribution = (aromatic_count as f64 * 1.5).floor() as i32;
        let effective_sum = aromatic_contribution.saturating_add(non_aromatic_sum);
        let v = normal_valences[0] as i32 + charge;
        if v <= 0 || effective_sum >= v {
            return 0;
        }
        return (v - effective_sum) as u8;
    }

    // Non-aromatic path (or post-Kekulization molecule where all bonds are explicit).
    let bond_sum = non_aromatic_sum;

    // For atoms that carry the aromatic flag but reside in a kekulized molecule
    // (bonds are Single/Double, not Aromatic), use only the lowest normal valence.
    // Rationale: after Kekulization, a substituted aromatic N (e.g. N−CH₃ in caffeine
    // with one ring double bond) has bond_sum=4, which would select valence 5 and
    // give 1 implicit H.  Capping at the primary valence (3) returns 0 H instead.
    let valences_to_check: &[u8] = if atom.aromatic {
        &normal_valences[..1]
    } else {
        normal_valences
    };

    // Iterate through valences (ascending) and pick the smallest ≥ bond_sum.
    for &v in valences_to_check {
        let target = v as i32 + charge;
        if target < 0 {
            continue;
        }
        if target >= bond_sum {
            return (target - bond_sum) as u8;
        }
    }

    // bond_sum exceeds all consulted valences → 0 implicit H.
    0
}

#[deprecated(
    since = "0.1.95",
    note = "use `implicit_hcount` directly — the two functions are identical"
)]
/// Alias for [`implicit_hcount`]; kept for API compatibility.
pub fn total_hcount(mol: &Molecule, idx: AtomIdx) -> u8 {
    implicit_hcount(mol, idx)
}

/// Sum of integer bond orders for heavy-atom bonds on `idx`.
/// Aromatic bonds count as 1 (pre-Kekulization representation).
pub fn bond_order_sum(mol: &Molecule, idx: AtomIdx) -> u8 {
    mol.neighbors(idx)
        .map(|(_, bidx)| mol.bond(bidx).order.order_int())
        .fold(0u8, |acc, x| acc.saturating_add(x))
}

/// Returns true if the bond is counted as a "double bond equivalent" in valence sums.
pub fn is_pi_bond(order: BondOrder) -> bool {
    matches!(
        order,
        BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
    )
}

// ---------------------------------------------------------------------------
// Valence validation
// ---------------------------------------------------------------------------

/// A valence violation on a specific atom.
///
/// Returned by [`validate_valence`] for each atom whose observed bond-order sum
/// exceeds all allowed normal valences (after formal-charge adjustment).
#[derive(Debug, Clone)]
pub struct ValenceError {
    /// Index of the over-valenced atom.
    pub atom: AtomIdx,
    /// Observed bond-order sum (+ explicit bracket H count).
    pub actual: u8,
    /// Allowed normal valences for the element (from [`crate::Element::normal_valences`]).
    pub allowed: &'static [u8],
}

impl fmt::Display for ValenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let valences_str = self
            .allowed
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(", ");
        write!(
            f,
            "atom {} has valence {} (allowed: [{}])",
            self.atom.0, self.actual, valences_str
        )
    }
}

impl std::error::Error for ValenceError {}

/// Check every atom in `mol` for valence violations.
///
/// Returns one [`ValenceError`] per over-valenced atom; an empty `Vec` means
/// all atoms have valid valence.
///
/// Atoms without defined normal valences (transition metals, etc.) are skipped.
/// Formal charge shifts the effective maximum: each unit of positive charge
/// adds one to the allowed ceiling (e.g. `[NH4+]` with 4 bonds is valid).
///
/// Aromatic bonds are counted as 1 each (`order_int()`).  Molecules still
/// written with `BondOrder::Aromatic` are handled correctly; fully kekulized
/// molecules are also supported.
pub fn validate_valence(mol: &Molecule) -> Vec<ValenceError> {
    let mut errors = Vec::new();
    for (idx, atom) in mol.atoms() {
        if atom.wildcard {
            continue;
        }
        let valences = atom.element.normal_valences();
        if valences.is_empty() {
            continue;
        }

        let bos = bond_order_sum(mol, idx);
        let explicit_h = atom.hydrogen_count.unwrap_or(0);
        let used = bos.saturating_add(explicit_h);
        let charge = atom.charge as i16;

        let has_valid = valences.iter().any(|&v| {
            let effective = (v as i16 + charge).max(0) as u8;
            effective >= used
        });

        if !has_valid {
            errors.push(ValenceError {
                atom: idx,
                actual: used,
                allowed: valences,
            });
        }
    }
    errors
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::BondOrder;
    use crate::element::Element;
    use crate::molecule::MoleculeBuilder;

    fn single_atom(elem: Element) -> Molecule {
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::organic(elem));
        b.build()
    }

    fn two_atoms(e1: Element, e2: Element, order: BondOrder) -> Molecule {
        let mut b = MoleculeBuilder::new();
        let a = b.add_atom(Atom::organic(e1));
        let c = b.add_atom(Atom::organic(e2));
        b.add_bond(a, c, order).unwrap();
        b.build()
    }

    #[test]
    fn test_methane() {
        // C alone: 0 bonds, valence 4 → 4 implicit H
        let mol = single_atom(Element::C);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 4);
    }

    #[test]
    fn test_ethane_c() {
        // CC: each C has 1 single bond → valence 4 → 3 implicit H
        let mol = two_atoms(Element::C, Element::C, BondOrder::Single);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 3);
        assert_eq!(implicit_hcount(&mol, AtomIdx(1)), 3);
    }

    #[test]
    fn test_ethylene_c() {
        // C=C: double bond → bond_sum=2 → 4-2=2 implicit H
        let mol = two_atoms(Element::C, Element::C, BondOrder::Double);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 2);
    }

    #[test]
    fn test_acetylene_c() {
        // C#C: triple bond → bond_sum=3 → 4-3=1 implicit H
        let mol = two_atoms(Element::C, Element::C, BondOrder::Triple);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 1);
    }

    #[test]
    fn test_nitrogen_amine() {
        // N alone: 0 bonds, first normal valence=3 → 3 implicit H (NH3)
        let mol = single_atom(Element::N);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 3);
    }

    #[test]
    fn test_nitrogen_triple() {
        // N#C: N has triple bond → bond_sum=3 → 3-3=0 (nitrile N)
        let mol = two_atoms(Element::N, Element::C, BondOrder::Triple);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 0);
    }

    #[test]
    fn test_oxygen_ether() {
        // O alone: 0 bonds, valence 2 → 2 implicit H (water)
        let mol = single_atom(Element::O);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 2);
    }

    #[test]
    fn test_fluorine() {
        // F alone: valence 1 → 1 implicit H (HF)
        let mol = single_atom(Element::F);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 1);
    }

    #[test]
    fn test_bracket_atom_explicit_h() {
        // [NH4+] — bracket atom: explicit H=4 returned directly
        let mut b = MoleculeBuilder::new();
        let atom = Atom::bracket(Element::N, None, Default::default(), 4, 1, None);
        b.add_atom(atom);
        let mol = b.build();
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 4);
    }

    #[test]
    fn test_hypervalent_sulfur() {
        // S with four single bonds: bond_sum=4, S valences=[2,4,6] → target=4 → 0 H
        let mut b = MoleculeBuilder::new();
        let s = b.add_atom(Atom::organic(Element::S));
        for _ in 0..4 {
            let c = b.add_atom(Atom::organic(Element::C));
            b.add_bond(s, c, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 0);
    }

    // ---------------------------------------------------------------------------
    // validate_valence tests
    // ---------------------------------------------------------------------------

    #[test]
    fn test_validate_valence_valid_molecules() {
        // All normal molecules should produce no errors.
        // methane (C, 0 bonds): valid
        let mol = single_atom(Element::C);
        assert!(
            validate_valence(&mol).is_empty(),
            "isolated C must be valid"
        );

        // water (O, 0 bonds): valid
        let mol = single_atom(Element::O);
        assert!(
            validate_valence(&mol).is_empty(),
            "isolated O must be valid"
        );

        // ethane (C–C): C has bond_sum=1, max valence 4 → valid
        let mol = two_atoms(Element::C, Element::C, BondOrder::Single);
        assert!(validate_valence(&mol).is_empty(), "ethane must be valid");

        // formaldehyde (C=O): C bond_sum=2, O bond_sum=2 → both valid
        let mol = two_atoms(Element::C, Element::O, BondOrder::Double);
        assert!(
            validate_valence(&mol).is_empty(),
            "formaldehyde must be valid"
        );
    }

    #[test]
    fn test_validate_valence_pentavalent_carbon() {
        // C with 5 single bonds: bond_sum=5 > max(C valences)=4 → error
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::organic(Element::C));
        for _ in 0..5 {
            let h = b.add_atom(Atom::new(Element::C));
            b.add_bond(c, h, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        let errors = validate_valence(&mol);
        assert_eq!(
            errors.len(),
            1,
            "C with 5 bonds must produce exactly 1 error"
        );
        assert_eq!(errors[0].atom, AtomIdx(0));
        assert_eq!(errors[0].actual, 5);
    }

    #[test]
    fn test_validate_valence_trivalent_oxygen() {
        // O with 3 single bonds: bond_sum=3 > max(O valences)=2 → error
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::organic(Element::O));
        for _ in 0..3 {
            let c = b.add_atom(Atom::organic(Element::C));
            b.add_bond(o, c, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        let errors = validate_valence(&mol);
        assert!(
            !errors.is_empty(),
            "O with 3 bonds must be flagged as over-valenced"
        );
        assert_eq!(errors[0].atom, AtomIdx(0));
    }

    #[test]
    fn test_validate_valence_ammonium_valid() {
        // [NH4+]: N with charge +1 and 4 bonds: effective max = 3+1=4 → valid
        let mut b = MoleculeBuilder::new();
        let mut n_atom = Atom::organic(Element::N);
        n_atom.charge = 1;
        let n = b.add_atom(n_atom);
        for _ in 0..4 {
            let c = b.add_atom(Atom::organic(Element::C));
            b.add_bond(n, c, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        assert!(
            validate_valence(&mol).is_empty(),
            "N+ with 4 bonds must be valid (ammonium-like)"
        );
    }

    // ---------------------------------------------------------------------------
    // Te (tellurium) tests — element.rs now has a real normal_valences() entry
    // for atomic number 52 ([2, 4, 6], source-verified against RDKit; see
    // element.rs::test_te_valence_source_verified). These pin the resulting
    // implicit_hcount()/validate_valence() behavior against RDKit's own output.
    // ---------------------------------------------------------------------------

    #[test]
    fn test_te_implicit_hcount_no_explicit_h_stays_zero() {
        // Te is outside the OpenSMILES organic subset, so implicit_hcount() must
        // still return 0 for a bare (non-bracket) Te atom, unchanged by the new
        // valence entry (guarded by the is_organic_subset() short-circuit).
        let mol = single_atom(Element::TE);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 0);
    }

    #[test]
    fn test_validate_valence_te_divalent_neutral() {
        // C[Te]C (dimethyl telluride analog). RDKit: explicitValence=2,
        // totalValence=2, sanitizes OK.
        let mut b = MoleculeBuilder::new();
        let te = b.add_atom(Atom::organic(Element::TE));
        for _ in 0..2 {
            let c = b.add_atom(Atom::organic(Element::C));
            b.add_bond(te, c, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        assert!(
            validate_valence(&mol).is_empty(),
            "divalent Te must be valid (RDKit sanitizes C[Te]C OK)"
        );
    }

    #[test]
    fn test_validate_valence_te_tetravalent_neutral() {
        // Cl[Te](Cl)(Cl)Cl (TeCl4 analog). RDKit: explicitValence=4, sanitizes OK.
        let mut b = MoleculeBuilder::new();
        let te = b.add_atom(Atom::organic(Element::TE));
        for _ in 0..4 {
            let c = b.add_atom(Atom::organic(Element::C));
            b.add_bond(te, c, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        assert!(
            validate_valence(&mol).is_empty(),
            "tetravalent Te must be valid (RDKit sanitizes TeCl4-analog OK)"
        );
    }

    #[test]
    fn test_validate_valence_te_hexavalent_neutral() {
        // F[Te](F)(F)(F)(F)F (TeF6 analog). RDKit: explicitValence=6, sanitizes OK.
        let mut b = MoleculeBuilder::new();
        let te = b.add_atom(Atom::organic(Element::TE));
        for _ in 0..6 {
            let c = b.add_atom(Atom::organic(Element::C));
            b.add_bond(te, c, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        assert!(
            validate_valence(&mol).is_empty(),
            "hexavalent Te must be valid (RDKit sanitizes TeF6-analog OK)"
        );
    }

    #[test]
    fn test_validate_valence_te_overvalent_invalid() {
        // 8 single bonds on Te (Cl x8 analog). RDKit rejects during sanitization:
        // "Explicit valence for atom # 1 Te, 8, is greater than permitted".
        let mut b = MoleculeBuilder::new();
        let te = b.add_atom(Atom::organic(Element::TE));
        for _ in 0..8 {
            let c = b.add_atom(Atom::organic(Element::C));
            b.add_bond(te, c, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        let errors = validate_valence(&mol);
        assert_eq!(
            errors.len(),
            1,
            "8-bonded Te must be flagged invalid, matching RDKit's AtomValenceException"
        );
        assert_eq!(errors[0].actual, 8);
        assert_eq!(errors[0].allowed, &[2, 4, 6]);
    }

    #[test]
    fn test_validate_valence_te_cation_telluronium() {
        // C[Te+](C)C (trimethyltelluronium). RDKit: charge=+1, degree=3,
        // explicitValence=3 (effective valence 2+1=3), sanitizes OK.
        let mut b = MoleculeBuilder::new();
        let mut te_atom = Atom::organic(Element::TE);
        te_atom.charge = 1;
        let te = b.add_atom(te_atom);
        for _ in 0..3 {
            let c = b.add_atom(Atom::organic(Element::C));
            b.add_bond(te, c, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        assert!(
            validate_valence(&mol).is_empty(),
            "Te+ telluronium (3 bonds) must be valid, matching RDKit"
        );
    }

    #[test]
    fn test_validate_valence_te_anion_telluride() {
        // [Te-2] (isolated telluride dianion). RDKit: charge=-2, degree=0,
        // explicitValence=0, sanitizes OK.
        let mut b = MoleculeBuilder::new();
        let te_atom = Atom::bracket(Element::TE, None, Default::default(), 0, -2, None);
        b.add_atom(te_atom);
        let mol = b.build();
        assert!(
            validate_valence(&mol).is_empty(),
            "[Te-2] must be valid, matching RDKit"
        );
    }

    #[test]
    fn test_validate_valence_te_anion_hydrotelluride() {
        // [TeH-]. RDKit: charge=-1, explicit H=1, totalValence=1, sanitizes OK.
        let mut b = MoleculeBuilder::new();
        let te_atom = Atom::bracket(Element::TE, None, Default::default(), 1, -1, None);
        b.add_atom(te_atom);
        let mol = b.build();
        assert!(
            validate_valence(&mol).is_empty(),
            "[TeH-] must be valid, matching RDKit"
        );
    }

    #[test]
    fn test_te_bracket_h2_implicit_and_valid() {
        // [TeH2]. RDKit: charge=0, explicit H=2, totalValence=2, sanitizes OK.
        // Bracket atoms return their stored H count directly from implicit_hcount().
        let mut b = MoleculeBuilder::new();
        let te_atom = Atom::bracket(Element::TE, None, Default::default(), 2, 0, None);
        let te = b.add_atom(te_atom);
        let mol = b.build();
        assert_eq!(implicit_hcount(&mol, te), 2);
        assert!(
            validate_valence(&mol).is_empty(),
            "[TeH2] must be valid, matching RDKit"
        );
    }

    #[test]
    fn test_validate_valence_transition_metal_skipped() {
        // Fe has no normal_valences → always valid regardless of bonds
        let mut b = MoleculeBuilder::new();
        let fe = b.add_atom(Atom::new(Element::FE));
        for _ in 0..6 {
            let c = b.add_atom(Atom::organic(Element::C));
            b.add_bond(fe, c, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        assert!(
            validate_valence(&mol).is_empty(),
            "Fe with 6 bonds must be skipped"
        );
    }

    // -------------------------------------------------------------------
    // Donor-side dative-bond implicit H -- regression tests for the
    // platinum coordination-chemistry benchmark
    // (validation/platinum/FEASIBILITY.md). A bare, un-bracketed donor
    // atom's own normal covalent valence must not be spent by its own
    // dative bond: `N->[Pt]` must still mean NH3, matching RDKit's
    // identical treatment of the same `->` SMILES syntax. Before this fix,
    // `order_int()`'s `1` for `Dative` was summed exactly like a real
    // covalent bond, giving NH2 instead.
    // -------------------------------------------------------------------

    #[test]
    fn test_dative_donor_n_keeps_full_valence() {
        // N->Pt : N is the donor (atom1). Bare N, no other substituents.
        let mol = two_atoms(Element::N, Element::PT, BondOrder::Dative);
        assert_eq!(
            implicit_hcount(&mol, AtomIdx(0)),
            3,
            "a dative donor's own lone pair, not one of its 3 normal covalent \
             slots, is shared with the acceptor -- N->[Pt] must mean NH3"
        );
    }

    #[test]
    fn test_dative_donor_o_keeps_full_valence() {
        // O->Fe : O is the donor. Matches the corpus's water/DMSO-oxygen
        // style dative ligands (see cisplatin_diaqua_activation_product in
        // pt_corpus.jsonl).
        let mol = two_atoms(Element::O, Element::FE, BondOrder::Dative);
        assert_eq!(
            implicit_hcount(&mol, AtomIdx(0)),
            2,
            "O->[Fe] must mean H2O"
        );
    }

    #[test]
    fn test_dative_acceptor_side_unaffected() {
        // Pt->N (reversed: Pt is the donor per this specific bond's stored
        // direction, N is the acceptor) -- the fix is donor-side-only by
        // design (see valence_inferred_hcount's doc comment), so N here
        // gets no special treatment; this pins that scope boundary rather
        // than leaving it implicit.
        let mol = two_atoms(Element::PT, Element::N, BondOrder::Dative);
        // N is atom2 (the acceptor) here: bond_sum=1 (unchanged, counted
        // like a normal covalent bond), giving valence 3 - 1 = 2H, not 3.
        assert_eq!(
            implicit_hcount(&mol, AtomIdx(1)),
            2,
            "acceptor-side implicit H is intentionally untouched by this fix"
        );
    }

    #[test]
    fn test_generalization_not_platinum_specific() {
        // Same fix, non-Pt acceptors (Fe, Co) -- confirms this is a general
        // dative-bond fix, not special-cased to Pt (task's generalization
        // gate, see FEASIBILITY.md section 16/6).
        for acceptor in [Element::FE, Element::CO, Element::PD, Element::RU] {
            let mol = two_atoms(Element::N, acceptor, BondOrder::Dative);
            assert_eq!(
                implicit_hcount(&mol, AtomIdx(0)),
                3,
                "N->{acceptor:?} must mean NH3 regardless of which metal accepts"
            );
        }
    }

    #[test]
    fn test_known_divergence_bracketed_dative_donor_can_still_false_positive_in_validate_valence() {
        // KNOWN, DOCUMENTED, NOT FIXED HERE (see FEASIBILITY.md's residual
        // limitations): `validate_valence` uses `bond_order_sum`, a
        // separate, PUBLIC function (also used by chematic-cip/tautomer/
        // chematic-ff) that was deliberately left untouched -- widening it
        // to exempt donor-side dative bonds too would change its meaning
        // for those other consumers, a much larger blast radius than this
        // benchmark measured or scoped. The practical consequence: a
        // BRACKETED dative donor whose only listed normal valence is
        // already exactly met by its own explicit H count (e.g. O, whose
        // only valence is 2) still gets a spurious `ValenceError` from
        // `validate_valence`, even though `implicit_hcount`/
        // `valence_inferred_hcount` correctly agree the atom is valid.
        // Nitrogen is NOT affected in practice ([NH3]->[Pt] does not
        // trigger this) only because N's valence list [3, 5] happens to
        // have a second, higher tier that absorbs the extra count -- an
        // element-specific coincidence, not a general exemption. This
        // corpus never hits this path (it uses bare, un-bracketed donor
        // atoms throughout), which is why `pt_corpus.jsonl`'s
        // `valence_errors` fields are all empty.
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::bracket(
            Element::O,
            None,
            Default::default(),
            2,
            0,
            None,
        ));
        let pt = b.add_atom(Atom::new(Element::PT));
        b.add_bond(o, pt, BondOrder::Dative).unwrap();
        let mol = b.build();

        assert_eq!(
            implicit_hcount(&mol, o),
            2,
            "bracket H is stored directly regardless of bond wiring"
        );
        assert_eq!(
            validate_valence(&mol).len(),
            1,
            "known false positive: [OH2]->[Pt] is valid water-donor chemistry \
             but validate_valence still flags it via bond_order_sum's \
             unchanged (donor-side-inclusive) count -- if this assertion \
             ever starts failing because it's empty, bond_order_sum was \
             fixed too and this whole test (and its FEASIBILITY.md note) \
             should be deleted, not \"corrected\" back to failing"
        );
    }
}
