//! Atropisomer detection and axial stereochemistry assignment.
//!
//! Detects rotationally constrained bonds (biaryl, allene) and assigns
//! M/P stereochemistry based on CIP rules adapted for axial centers.

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};
use chematic_perception::find_sssr;

/// Type of atropisomeric bond/center.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtropisomerType {
    /// Biaryl: aromatic C - aromatic C with ortho substituents
    Biaryl,
    /// Allene: C=C=C linear system
    Allene,
    /// Constrained: bond in fused ring or strained system
    Constrained,
}

/// Given an SSSR ring (ordered atom cycle) and one of its members, return the
/// two atoms adjacent to it *within that ring* (its ortho positions relative
/// to `atom`'s exocyclic bond).
fn ring_ortho_atoms(ring: &[AtomIdx], atom: AtomIdx) -> Option<(AtomIdx, AtomIdx)> {
    let idx = ring.iter().position(|&a| a == atom)?;
    let n = ring.len();
    Some((ring[(idx + n - 1) % n], ring[(idx + 1) % n]))
}

/// Detect atropisomeric (rotationally constrained) bonds in a molecule.
///
/// Returns list of bond indices and their atropisomer type.
/// Detects: biaryl systems with ortho substitution, allenes, and other constrained bonds.
///
/// ## Scope / simplifications
///
/// Biaryl detection uses an intentionally simplified heuristic, not full
/// IUPAC-grade atropisomer chemistry (e.g. Cahn's "3 of 4 ortho substituents"
/// bulk-threshold rule, or a real steric-bulk/A-value model): a bond is
/// flagged as `AtropisomerType::Biaryl` when it connects two aromatic carbons
/// in two *different* SSSR rings (i.e. a genuine inter-ring connector, not a
/// fused-ring shared edge) and **at least one** ring atom adjacent to each
/// ipso carbon (an "ortho position") has degree > 2 (i.e. carries some
/// non-H substituent) on **both** rings. This does not weigh substituent
/// size/bulk, doesn't require *both* ortho positions per ring, and doesn't
/// model 5-membered heteroaryl axes. It is structural (SSSR-based) rather
/// than `BondOrder`-based, so it is invariant to whether the SMILES wrote an
/// explicit `-` between the rings or left the bond implicit.
///
/// Two further known imprecisions of the "degree > 2" ortho test: a
/// ring-fusion atom (e.g. the shared carbon of a 1-arylnaphthalene's fused
/// ring) has degree 3 from its own ring bonds alone and reads as
/// "substituted" even with no exocyclic group there — arguably reasonable
/// (a fused ring is real steric bulk too) but not distinguished from an
/// actual substituent. And for biaryls where one ring is itself part of a
/// fused polycyclic system, which SSSR cycle the ortho lookup walks can be
/// affected by the same envelope-ring/fundamental-cycle SSSR decomposition
/// artifacts documented on [`chematic_perception::count_aromatic_rings`];
/// this function does not apply that correction.
pub fn detect_atropisomers(mol: &Molecule) -> Vec<(BondIdx, AtropisomerType)> {
    let mut result = Vec::new();
    let rings = find_sssr(mol);
    let rings = rings.rings();

    for (bidx, bond) in mol.bonds() {
        let a1 = mol.atom(bond.atom1);
        let a2 = mol.atom(bond.atom2);

        // Check for biaryl: both aromatic C, each in its own aromatic SSSR
        // ring, connected by a genuine inter-ring bond (no shared ring).
        // Deliberately independent of `bond.order`: whether the SMILES wrote
        // an explicit single bond or left it implicit must not change the
        // answer for the same real molecule (issue #262).
        if a1.aromatic
            && a2.aromatic
            && a1.element.atomic_number() == 6
            && a2.element.atomic_number() == 6
        {
            let ring1 = rings
                .iter()
                .find(|r| r.contains(&bond.atom1) && r.iter().all(|&a| mol.atom(a).aromatic));
            let ring2 = rings
                .iter()
                .find(|r| r.contains(&bond.atom2) && r.iter().all(|&a| mol.atom(a).aromatic));

            let shares_ring = rings
                .iter()
                .any(|r| r.contains(&bond.atom1) && r.contains(&bond.atom2));

            if let (Some(ring1), Some(ring2)) = (ring1, ring2)
                && !shares_ring
            {
                // Ortho check: at least one ring-adjacent neighbor of each
                // ipso carbon must carry a non-ring substituent (degree > 2).
                let ring1_ortho_substituted =
                    ring_ortho_atoms(ring1, bond.atom1).is_some_and(|(p, n)| {
                        mol.neighbors(p).count() > 2 || mol.neighbors(n).count() > 2
                    });
                let ring2_ortho_substituted =
                    ring_ortho_atoms(ring2, bond.atom2).is_some_and(|(p, n)| {
                        mol.neighbors(p).count() > 2 || mol.neighbors(n).count() > 2
                    });

                if ring1_ortho_substituted && ring2_ortho_substituted {
                    result.push((bidx, AtropisomerType::Biaryl));
                }
            }
        }

        // Check for allene: C=C=C
        if a1.element.atomic_number() == 6
            && a2.element.atomic_number() == 6
            && bond.order == BondOrder::Double
        {
            // Check if both neighbors have another double bond (allene pattern)
            let a1_has_double = mol
                .neighbors(bond.atom1)
                .filter(|(n, _)| n != &bond.atom2)
                .any(|(_, nb)| mol.bond(nb).order == BondOrder::Double);

            let a2_has_double = mol
                .neighbors(bond.atom2)
                .filter(|(n, _)| n != &bond.atom1)
                .any(|(_, nb)| mol.bond(nb).order == BondOrder::Double);

            if a1_has_double && a2_has_double {
                result.push((bidx, AtropisomerType::Allene));
            }
        }
    }

    result
}

/// Assign M/P stereochemistry to atropisomeric bonds based on CIP priorities.
///
/// Returns molecule with wedge/dash bonds (Up/Down) annotated for atropisomers.
/// For each atropisomeric bond, assigns chirality direction based on:
/// - Priority of substituents at each stereogenic atom
/// - CIP rule (atomic number, mass, connectivity, recursion)
pub fn assign_atropisomer_chirality(mol: &Molecule) -> Molecule {
    use chematic_core::MoleculeBuilder;

    let atropisomers = detect_atropisomers(mol);

    // Always rebuild to ensure consistent output type
    let mut builder = MoleculeBuilder::new();
    let mut remap = std::collections::HashMap::new();

    // Copy all atoms
    for (idx, atom) in mol.atoms() {
        let new_idx = builder.add_atom(atom.clone());
        remap.insert(idx, new_idx);
    }

    // Copy bonds, applying stereochemistry to atropisomeric bonds
    for (bidx, bond) in mol.bonds() {
        let mut new_bond_order = bond.order;

        // Check if this bond is atropisomeric and apply stereochemistry.
        // Gate on `AtropisomerType::Biaryl` (detect_atropisomers's own
        // structural classification), not `bond.order == BondOrder::Single`:
        // the biaryl inter-ring bond parses as `BondOrder::Aromatic` when
        // the SMILES leaves it implicit (two aromatic-carbon atoms with no
        // explicit bond symbol between them), so the old `Single`-only
        // check silently skipped chirality assignment for that notation
        // even though `detect_atropisomers` (fixed in #271) already
        // correctly reports the bond as atropisomeric independent of
        // `bond.order` (issue #276). Restricting to `Biaryl` here (rather
        // than dropping the type check altogether) preserves the existing,
        // separate no-op behavior for `AtropisomerType::Allene`, whose
        // reported bond is the allene's own central *double* bond -- that
        // one must never be overwritten with `Up`/`Down`, which would
        // silently destroy the double bond.
        if let Some((_, atrop_type)) = atropisomers.iter().find(|(b, _)| b == &bidx)
            && *atrop_type == AtropisomerType::Biaryl
        {
            let a1_neighbors: Vec<_> = mol.neighbors(bond.atom1).collect();
            let a2_neighbors: Vec<_> = mol.neighbors(bond.atom2).collect();

            let a1_max_an = a1_neighbors
                .iter()
                .filter(|(n, _)| n != &bond.atom2)
                .map(|(n, _)| mol.atom(*n).element.atomic_number())
                .max()
                .unwrap_or(0);

            let a2_max_an = a2_neighbors
                .iter()
                .filter(|(n, _)| n != &bond.atom1)
                .map(|(n, _)| mol.atom(*n).element.atomic_number())
                .max()
                .unwrap_or(0);

            if a1_max_an > a2_max_an {
                new_bond_order = BondOrder::Up;
            } else if a2_max_an > a1_max_an {
                new_bond_order = BondOrder::Down;
            }
        }

        if let (Some(&a), Some(&b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(a, b, new_bond_order);
        }
    }

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn detect_atropisomers_biaryl() {
        // Unsubstituted biphenyl has no ortho substituents on either ring,
        // so it must not be flagged as an atropisomer (implicit inter-ring
        // bond notation).
        let m = mol("c1ccccc1c2ccccc2");
        let atrops = detect_atropisomers(&m);
        assert!(
            atrops.is_empty(),
            "biphenyl without ortho subs should have no atropisomers"
        );
    }

    #[test]
    fn detect_atropisomers_biphenyl_notation_invariant() {
        // issue #262 bug 1: the same real molecule (biphenyl) must return the
        // same result regardless of whether the SMILES spells the inter-ring
        // bond explicitly (`-`) or leaves it implicit.
        let explicit = mol("c1ccccc1-c2ccccc2");
        let implicit = mol("c1ccccc1c2ccccc2");
        assert_eq!(
            detect_atropisomers(&explicit).len(),
            detect_atropisomers(&implicit).len(),
            "biphenyl must give the same result under both notations"
        );
        assert!(detect_atropisomers(&explicit).is_empty());
        assert!(detect_atropisomers(&implicit).is_empty());
    }

    #[test]
    fn detect_atropisomers_ortho_dimethylbiphenyl() {
        // 2,2'-dimethylbiphenyl: bulky methyl at the ortho position on both
        // rings -> genuine atropisomer.
        let m = mol("Cc1ccccc1-c1ccccc1C");
        let atrops = detect_atropisomers(&m);
        assert_eq!(
            atrops.len(),
            1,
            "2,2'-dimethylbiphenyl should have exactly 1 atropisomeric bond"
        );
        assert_eq!(atrops[0].1, AtropisomerType::Biaryl);
    }

    #[test]
    fn detect_atropisomers_ortho_dimethylbiphenyl_notation_invariant() {
        // Same as the biphenyl notation-invariance test above, but on a
        // *positive* case: an empty-vs-empty match alone can't distinguish
        // "notation-invariant" from "the biaryl branch never fires" (which
        // was exactly how the old `detect_atropisomers_biaryl` test passed
        // for the wrong reason, see issue #262). This asserts the same
        // nonzero result under both notations.
        let explicit = mol("Cc1ccccc1-c2ccccc2C");
        let implicit = mol("Cc1ccccc1c2ccccc2C");
        assert_eq!(
            detect_atropisomers(&explicit).len(),
            detect_atropisomers(&implicit).len(),
            "2,2'-dimethylbiphenyl must give the same result under both notations"
        );
        assert_eq!(detect_atropisomers(&explicit).len(), 1);
        assert_eq!(detect_atropisomers(&implicit).len(), 1);
    }

    #[test]
    fn detect_atropisomers_para_dimethylbiphenyl_not_flagged() {
        // issue #262 bug 2: 4,4'-dimethylbiphenyl has substituents at the
        // *para* position, not ortho -> no rotational hindrance, must not be
        // flagged.
        let m = mol("Cc1ccc(cc1)-c1ccc(C)cc1");
        let atrops = detect_atropisomers(&m);
        assert!(
            atrops.is_empty(),
            "4,4'-dimethylbiphenyl (para-substituted) should have no atropisomers"
        );
    }

    #[test]
    fn detect_atropisomers_naphthalene_fusion_bond_not_biaryl() {
        // Regression/robustness check on a real fused polycyclic aromatic:
        // both bridgehead atoms of naphthalene have degree 3 (like any
        // biaryl ipso carbon), so the old `a1_degree >= 3 && a2_degree >= 3`
        // heuristic was one `bond.order` accident away from misclassifying
        // a ring-fusion bond as a biaryl axis. It happened not to fire here
        // because `apply_aromaticity()` (verified below) assigns
        // `BondOrder::Aromatic` to ring bonds, which the old code's
        // `bond.order == BondOrder::Single` gate excluded — but that made
        // the old code's correctness on fused systems an accident of
        // bond-order bookkeeping, not a real check. The new SSSR-based
        // `shares_ring` exclusion is not order-dependent and rejects this
        // bond on structural grounds: naphthalene's two SSSR rings both
        // contain the bridgehead pair, so the fusion bond is correctly
        // recognized as a shared/fused-ring edge, not a genuine inter-ring
        // connector, regardless of what bond order it carries.
        let m = mol("C1=CC2=CC=CC=C2C=C1");
        let m = chematic_perception::apply_aromaticity(&m);

        // Pin down the premise this test claims to exercise, rather than
        // asserting only the end result (which could pass for unrelated
        // reasons, e.g. apply_aromaticity not flagging anything aromatic).
        assert!(
            m.atoms().all(|(_, a)| a.aromatic),
            "apply_aromaticity must flag every naphthalene atom aromatic, else this test is vacuous"
        );
        let bridgeheads: Vec<_> = m
            .atoms()
            .filter(|(idx, _)| m.neighbors(*idx).count() == 3)
            .map(|(idx, _)| idx)
            .collect();
        assert_eq!(
            bridgeheads.len(),
            2,
            "naphthalene should have exactly 2 degree-3 bridgehead atoms"
        );
        let (_, fusion_bond) = m
            .bond_between(bridgeheads[0], bridgeheads[1])
            .expect("bridgeheads must be directly bonded (the fusion bond)");
        assert_eq!(
            fusion_bond.order,
            BondOrder::Aromatic,
            "apply_aromaticity is expected to assign BondOrder::Aromatic to ring bonds; if this \
             changes, the comment above (and the PR description) needs updating to match"
        );

        let atrops = detect_atropisomers(&m);
        assert!(
            atrops.is_empty(),
            "naphthalene's ring-fusion bond must not be flagged as a biaryl axis, got {atrops:?}"
        );
    }

    #[test]
    fn detect_atropisomers_none() {
        // No atropisomers
        let m = mol("CC");
        let atrops = detect_atropisomers(&m);
        assert_eq!(atrops.len(), 0, "ethane should have no atropisomers");
    }

    #[test]
    fn assign_atropisomer_chirality_preserves_atoms() {
        let m = mol("c1ccccc1c2ccccc2");
        let result = assign_atropisomer_chirality(&m);
        assert_eq!(
            result.atom_count(),
            m.atom_count(),
            "atom count should match"
        );
    }

    #[test]
    fn assign_atropisomer_chirality_preserves_bonds() {
        let m = mol("c1ccccc1c2ccccc2");
        let result = assign_atropisomer_chirality(&m);
        assert_eq!(
            result.bond_count(),
            m.bond_count(),
            "bond count should match"
        );
    }

    /// Returns the `BondOrder` `assign_atropisomer_chirality` wrote onto the
    /// molecule's single `Biaryl` atropisomeric bond. Panics if there isn't
    /// exactly one, since these tests need that premise to hold to mean
    /// anything.
    fn biaryl_result_order(smi: &str) -> BondOrder {
        let m = mol(smi);
        let atrops = detect_atropisomers(&m);
        let biaryl: Vec<_> = atrops
            .iter()
            .filter(|(_, t)| *t == AtropisomerType::Biaryl)
            .collect();
        assert_eq!(
            biaryl.len(),
            1,
            "'{smi}' should have exactly 1 biaryl atropisomeric bond, got {atrops:?}"
        );
        let (bidx, _) = biaryl[0];
        // Bond indices are preserved 1:1 by assign_atropisomer_chirality's
        // rebuild (it iterates `mol.bonds()` in order and appends to the
        // builder in the same order), so `*bidx` still identifies the same
        // bond in the rebuilt molecule.
        assign_atropisomer_chirality(&m).bond(*bidx).order
    }

    #[test]
    fn assign_atropisomer_chirality_notation_invariant() {
        // issue #276: detect_atropisomers (fixed in #271) is notation-
        // invariant, but assign_atropisomer_chirality had its own separate
        // `bond.order == BondOrder::Single` gate that silently skipped
        // chirality assignment for the same bond when the SMILES left the
        // inter-ring bond implicit (which parses as `BondOrder::Aromatic`,
        // not `Single` -- confirmed empirically, see PR description).
        //
        // A plain 2,2'-disubstituted biphenyl (e.g. `Cc1ccccc1-c1ccccc1C`)
        // doesn't actually exercise this: this function's CIP-priority
        // heuristic only compares the *immediate* ring neighbors of each
        // ipso carbon, which are both plain aromatic carbons on either side
        // of an all-carbocyclic biphenyl axis, so the comparison ties and
        // no Up/Down is assigned either way (masking the bug rather than
        // demonstrating it). Using a biaryl with a ring nitrogen ortho to
        // one ipso carbon (2-methylpyridin-3-yl vs 2-methylphenyl) breaks
        // the tie and produces a real, distinguishing assignment.
        let explicit = "Cc1cccnc1-c1ccccc1C";
        let implicit = "Cc1cccnc1c1ccccc1C";
        let explicit_order = biaryl_result_order(explicit);
        let implicit_order = biaryl_result_order(implicit);

        assert!(
            matches!(explicit_order, BondOrder::Up | BondOrder::Down),
            "explicit-notation biaryl should get a real Up/Down chirality assignment, got {explicit_order:?}"
        );
        assert_eq!(
            explicit_order, implicit_order,
            "assign_atropisomer_chirality must give the same chirality assignment for the same \
             real molecule regardless of whether the inter-ring bond is written explicitly or \
             left implicit"
        );
    }

    #[test]
    fn assign_atropisomer_chirality_leaves_allene_double_bond_untouched() {
        // Regression guard for the *other* way this fix could have gone
        // wrong: naively deleting the old `bond.order == BondOrder::Single`
        // check entirely (rather than replacing it with the
        // `AtropisomerType::Biaryl` type gate) would have let this loop
        // also try to stamp Up/Down onto whatever bond
        // `detect_atropisomers` reports as `AtropisomerType::Allene` --
        // which is the allene's own central *double* bond, not a single
        // bond. Overwriting that would silently corrupt the bond order
        // (Double -> Up/Down). The `Biaryl`-only gate this fix actually
        // uses excludes `Allene` entries, preserving the pre-fix behavior
        // that allene bonds are never rewritten here (this function's
        // actual M/P assignment is only implemented for the biaryl case).
        //
        // Note: `detect_atropisomers`'s `Allene` branch, as currently
        // written, only fires for a cumulated *triene* (4+ consecutive
        // double bonds) and never for a plain 3-carbon allene (`C=C=C`) --
        // a separate, pre-existing gap outside #276's scope (see PR
        // description). `CC=C=C=CC` (a butatriene) is used here because it
        // is the shortest input that actually produces an `Allene` entry
        // to exercise this guard against.
        let m = mol("CC=C=C=CC");
        let atrops = detect_atropisomers(&m);
        let (bidx, _) = atrops
            .iter()
            .find(|(_, t)| *t == AtropisomerType::Allene)
            .expect("'CC=C=C=CC' should have an Allene atropisomer entry");
        assert_eq!(
            assign_atropisomer_chirality(&m).bond(*bidx).order,
            BondOrder::Double,
            "allene central double bond must not be overwritten with Up/Down"
        );
    }
}
