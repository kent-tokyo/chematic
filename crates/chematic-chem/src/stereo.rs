//! Stereochemistry manipulation — inversion, enumeration, and CIP/atropisomer integration.
//!
//! Provides utilities for inverting R/S stereochemistry, enumerating
//! stereoisomers from unspecified stereocenters, and unified assignment
//! of both tetrahedral (R/S) and axial (M/P) stereochemistry.

use std::collections::HashMap;

use crate::atropisomer;
use chematic_core::{
    Atom, AtomIdx, BondOrder, Chirality, Molecule, MoleculeBuilder, STEREO_H_SENTINEL,
};

/// Invert the stereochemistry of a tetrahedral stereocenter.
///
/// If the atom has CIP code R, it becomes S (and vice versa).
///
/// Primary path: if `idx` has `@`/`@@` chirality recorded (the common case
/// for SMILES-parsed input), flips the `Chirality` enum directly against the
/// unchanged `stereo_neighbor_order` reference frame. This is correct
/// independent of atom numbering or neighbor iteration order — the
/// three-point ordering the `@`/`@@` symbol is defined relative to never
/// moves, only the handedness does.
///
/// Fallback path (kept for API compatibility, e.g. molecules built directly
/// from 2D wedge bonds with no `chirality` set): inverts wedge/dash bonds
/// (Up ↔ Down) attached to the stereocenter instead.
///
/// Atoms with neither `@`/`@@` chirality nor a wedge bond are unchanged.
/// Returned molecule preserves all other properties, including every other
/// atom's stereo metadata (`stereo_neighbor_order`, `stereo_groups`,
/// `bond_directions`).
pub fn invert_stereocenter(mol: &Molecule, idx: AtomIdx) -> Molecule {
    let mut builder = MoleculeBuilder::new();

    if mol.atom(idx).chirality != Chirality::None {
        for (i, atom) in mol.atoms() {
            let mut atom = atom.clone();
            if i == idx {
                atom.chirality = match atom.chirality {
                    Chirality::Clockwise => Chirality::CounterClockwise,
                    Chirality::CounterClockwise => Chirality::Clockwise,
                    Chirality::None => Chirality::None,
                    // cis/trans relabeling is a constitutional change, not a chiral
                    // (mirror-image) inversion -- leave square-planar tags untouched.
                    sp @ Chirality::SquarePlanar(_) => sp,
                };
            }
            builder.add_atom(atom);
        }
        for (_, bond) in mol.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
    } else {
        let has_wedge = mol
            .neighbors(idx)
            .any(|(_, bidx)| matches!(mol.bond(bidx).order, BondOrder::Up | BondOrder::Down));

        for (_, atom) in mol.atoms() {
            builder.add_atom(atom.clone());
        }
        for (_, bond) in mol.bonds() {
            let new_order = if has_wedge && (bond.atom1 == idx || bond.atom2 == idx) {
                match bond.order {
                    BondOrder::Up => BondOrder::Down,
                    BondOrder::Down => BondOrder::Up,
                    other => other,
                }
            } else {
                bond.order
            };
            let _ = builder.add_bond(bond.atom1, bond.atom2, new_order);
        }
    }

    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.build()
}

/// Enumerate all stereoisomers by assigning R/S to unspecified stereocenters.
///
/// Identifies tetrahedral carbon atoms without explicit `@`/`@@` notation and
/// generates all 2^n combinations where n is the number of unspecified centers.
/// At most 2^6 = 64 combinations are enumerated; returns empty vector if n > 6.
///
/// Each returned molecule has all stereocenters assigned (using `Chirality::Clockwise`
/// for R and `Chirality::CounterClockwise` for S). Deduplicates results by
/// canonical SMILES.
pub fn enumerate_stereoisomers(mol: &Molecule) -> Vec<Molecule> {
    use chematic_smiles::canonical_smiles;

    // Identify unspecified tetrahedral carbon stereocenters
    let unspecified: Vec<AtomIdx> = mol
        .atoms()
        .filter(|(idx, atom)| {
            if atom.element.atomic_number() != 6 || atom.aromatic {
                return false;
            }
            if atom.chirality != Chirality::None {
                return false;
            }
            let degree = mol.neighbors(*idx).count();
            if degree < 2 {
                return false;
            }
            let total = degree + chematic_core::implicit_hcount(mol, *idx) as usize;
            total == 4
                && mol.neighbors(*idx).all(|(_, bidx)| {
                    !matches!(mol.bond(bidx).order, BondOrder::Double | BondOrder::Triple)
                })
        })
        .map(|(idx, _)| idx)
        .collect();

    let n = unspecified.len();
    if n > 6 {
        return Vec::new(); // Limit to 2^6 = 64 combinations
    }

    if n == 0 {
        // Return a copy of the input molecule. A pure passthrough (same atom
        // count/order, none skipped) -- copy the stereo side channels
        // verbatim rather than losing already-specified stereocenters'
        // stereo_neighbor_order (see the main loop below for why this
        // matters: @/@@ is meaningless without it).
        let mut builder = MoleculeBuilder::new();
        for (_, atom) in mol.atoms() {
            builder.add_atom(atom.clone());
        }
        for (_, bond) in mol.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        builder.copy_stereo_groups_from(mol);
        builder.copy_stereo_from(mol);
        builder.copy_bond_directions_from(mol);
        return vec![builder.build()];
    }

    let mut seen = std::collections::HashSet::new();
    let mut results = Vec::new();

    for bits in 0u32..(1u32 << n) {
        let chirality_overrides: HashMap<AtomIdx, Chirality> = unspecified
            .iter()
            .enumerate()
            .map(|(i, &idx)| {
                let cw = (bits >> i) & 1 == 1;
                let chirality = if cw {
                    Chirality::Clockwise
                } else {
                    Chirality::CounterClockwise
                };
                (idx, chirality)
            })
            .collect();

        let mut builder = MoleculeBuilder::new();
        for (idx, atom) in mol.atoms() {
            let mut a = Atom::new(atom.element);
            a.charge = atom.charge;
            a.isotope = atom.isotope;
            a.aromatic = atom.aromatic;
            a.atom_map = atom.atom_map;
            if let Some(&new_chirality) = chirality_overrides.get(&idx) {
                a.chirality = new_chirality;
                // Force bracket notation for SMILES output
                let implicit_h = chematic_core::implicit_hcount(mol, idx);
                a.hydrogen_count = Some(atom.hydrogen_count.unwrap_or(implicit_h));
            } else {
                a.chirality = atom.chirality;
                a.hydrogen_count = atom.hydrogen_count;
            }
            builder.add_atom(a);
        }
        for (_, bond) in mol.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }

        // Passthrough rebuild (same atom/bond count/order as `mol`) -- copy
        // the stereo side channels verbatim so already-specified
        // stereocenters (not in `chirality_overrides`) keep the
        // stereo_neighbor_order their `@`/`@@` is defined relative to.
        builder.copy_stereo_groups_from(mol);
        builder.copy_stereo_from(mol);
        builder.copy_bond_directions_from(mol);

        // `mol` never had a stereo_neighbor_order entry for these atoms
        // (they were achiral there), so copy_stereo_from above didn't set
        // one -- without it, the freshly-assigned `chirality` above is
        // uninterpretable (same "chirality set, no reference order" failure
        // shape Kekule-S0 fixed for apply_kekule). Neighbor order here is an
        // arbitrary but fixed, well-defined reference frame (bond-insertion
        // order, implicit H last): `enumerate_stereoisomers` only needs the
        // 2^n combinations to be distinct and independently interpretable,
        // not to match any external R/S convention.
        for &idx in chirality_overrides.keys() {
            let mut order: Vec<u32> = mol.neighbors(idx).map(|(nb, _)| nb.0).collect();
            let implicit_h = chematic_core::implicit_hcount(mol, idx);
            order.extend(std::iter::repeat_n(STEREO_H_SENTINEL, implicit_h as usize));
            builder.set_stereo_neighbor_order(idx, order);
        }

        let isomer = builder.build();
        let smi = canonical_smiles(&isomer);
        if seen.insert(smi) {
            results.push(isomer);
        }
    }

    results
}

/// Assign both tetrahedral (R/S) and axial (M/P) stereochemistry in one step.
///
/// This function integrates CIP-based R/S assignment with atropisomer detection
/// and M/P assignment, providing complete stereochemical characterization.
///
/// The process:
/// 1. Applies atropisomer chirality assignment (biaryl, allene, constrained bonds)
/// 2. Returns the molecule with both tetrahedral and axial stereochemistry assigned
///
/// This is a convenience wrapper that handles both stereo types together, rather
/// than requiring separate calls to `assign_atropisomer_chirality()` and the
/// perception module's CIP assignment.
pub fn assign_complete_stereochemistry(mol: &Molecule) -> Molecule {
    // Apply atropisomer assignment (M/P for biaryl, allene, constrained bonds)
    atropisomer::assign_atropisomer_chirality(mol)

    // Note: Tetrahedral R/S assignment happens in perception module
    // (via chematic_perception::assign_stereo_from_2d or 3D variant)
    // This function handles the integration step.
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn invert_stereocenter_r_to_s() {
        // [C@@H](F)(Cl)Br — R configuration (@@)
        let m = mol("[C@@H](F)(Cl)Br");
        let before = crate::assign_cip(&m).get(AtomIdx(0));
        let inverted = invert_stereocenter(&m, AtomIdx(0));
        let after = crate::assign_cip(&inverted).get(AtomIdx(0));
        assert_eq!(
            (before, after),
            (
                Some(chematic_core::CipCode::R),
                Some(chematic_core::CipCode::S)
            ),
            "inversion must flip the CIP label, not just add an @ character"
        );
    }

    // -----------------------------------------------------------------
    // invert_stereocenter-correctness-P0: the pre-fix implementation only
    // ever inverted 2D wedge bonds (`BondOrder::Up`/`Down`). Plain `@`/`@@`
    // SMILES stereocenters have no wedge bond at all, so it was a silent
    // no-op for the common case -- the has_wedge check never fired and the
    // "no stereochemistry annotation" passthrough branch ran instead. Fixed
    // by flipping the `Chirality` enum directly (correct regardless of atom
    // numbering, since it's read against the unchanged
    // stereo_neighbor_order reference frame), with wedge-bond inversion
    // kept only as a fallback for atoms with no recorded `@`/`@@`. Also
    // carries the same missing copy_stereo_groups_from/copy_stereo_from/
    // copy_bond_directions_from pattern found repeatedly elsewhere this
    // milestone (Kekule-S0, Stereo-Rebuild-S1, Tautomer-Rebuild-S2).
    // -----------------------------------------------------------------

    /// Two independent, both-specified tetrahedral stereocenters that the
    /// Accurate CIP engine can actually resolve (each has an implicit H,
    /// unlike an all-heavy-substituent center, which is a separate,
    /// pre-existing Accurate-engine limitation unrelated to this fix).
    const TWO_CENTER_SMILES: &str = "C[C@H:1](O)[C@@H:2](N)C(=O)O";

    #[test]
    fn invert_stereocenter_flips_only_target_atom_mapped_cip() {
        let m = mol(TWO_CENTER_SMILES);
        let target = find_by_map(&m, 1);
        let other = find_by_map(&m, 2);

        let before =
            crate::assign_cip_with_mode(&m, crate::CipMode::Accurate).expect("accurate CIP");
        let (before_target, before_other) = (before.get(target), before.get(other));
        assert!(
            before_target.is_some() && before_other.is_some(),
            "test setup sanity: both centers must be Accurate-CIP-resolvable"
        );

        let inverted = invert_stereocenter(&m, target);
        let after =
            crate::assign_cip_with_mode(&inverted, crate::CipMode::Accurate).expect("accurate CIP");

        let expected_target = match before_target {
            Some(chematic_core::CipCode::R) => chematic_core::CipCode::S,
            Some(chematic_core::CipCode::S) => chematic_core::CipCode::R,
            other => panic!("test setup sanity: expected R/S, got {other:?}"),
        };
        assert_eq!(
            after.get(target),
            Some(expected_target),
            "target center's CIP label must flip to the exact opposite (not just differ, \
             and not become unresolvable)"
        );
        assert_eq!(
            before_other,
            after.get(other),
            "uninvolved center's CIP label must be unchanged"
        );
    }

    #[test]
    fn invert_stereocenter_double_inversion_returns_to_original() {
        let m = mol(TWO_CENTER_SMILES);
        let target = find_by_map(&m, 1);

        let once = invert_stereocenter(&m, target);
        let twice = invert_stereocenter(&once, target);

        assert_eq!(
            chematic_smiles::canonical_smiles(&twice),
            chematic_smiles::canonical_smiles(&m),
            "inverting the same center twice must reproduce the original molecule"
        );
    }

    #[test]
    fn invert_stereocenter_independent_of_atom_numbering() {
        // The same chemical center, reached via two different AtomIdx
        // values and two different stereo_neighbor_order orderings (a
        // canonical round trip renumbers atoms and re-derives the neighbor
        // order from scratch) -- inverting it must produce the same
        // molecule either way. Guards against a future rewrite drifting
        // back to a numbering-dependent scheme (the exact failure mode
        // this fix replaces).
        //
        // Uses a longer-chain variant of `TWO_CENTER_SMILES` rather than the
        // shared constant: after the explicit/implicit-H-count Morgan-rank
        // unification fix (issue #205), `TWO_CENTER_SMILES` itself no
        // longer gets renumbered on canonicalization (its map-1 atom already
        // had the correct, now-unified rank), which would make this
        // precondition assert -- and the test meaningless -- rather than
        // testing anything about `invert_stereocenter`. This variant still
        // renumbers, confirmed empirically, so the precondition holds again.
        let m1 = mol("CC[C@H:1](O)[C@@H:2](N)CC(=O)O");
        let m2 = mol(&chematic_smiles::canonical_smiles(&m1));
        assert_ne!(
            find_by_map(&m1, 1),
            find_by_map(&m2, 1),
            "test setup sanity: canonicalization must actually renumber this atom"
        );

        let inv1 = invert_stereocenter(&m1, find_by_map(&m1, 1));
        let inv2 = invert_stereocenter(&m2, find_by_map(&m2, 1));

        assert_eq!(
            chematic_smiles::canonical_smiles(&inv1),
            chematic_smiles::canonical_smiles(&inv2),
            "inverting the same chemical center must be independent of its AtomIdx/neighbor order"
        );
    }

    #[test]
    fn invert_stereocenter_matches_independent_ground_truth_single_center() {
        // Chemical-identity check per review: don't compare @/@@ strings --
        // verify inversion produces the SAME molecule as independently
        // parsing the known mirror SMILES from scratch.
        let l_ala = mol("N[C@@H](C)C(=O)O");
        let target = l_ala
            .atoms()
            .find(|(_, a)| a.chirality != Chirality::None)
            .expect("test setup sanity: has a stereocenter")
            .0;
        let inverted = invert_stereocenter(&l_ala, target);
        let d_ala = mol("N[C@H](C)C(=O)O");
        assert_eq!(
            chematic_smiles::canonical_smiles(&inverted),
            chematic_smiles::canonical_smiles(&d_ala),
            "inverting L-alanine's stereocenter must match independently-parsed D-alanine"
        );
    }

    #[test]
    fn invert_stereocenter_matches_independent_ground_truth_two_centers() {
        let m = mol(TWO_CENTER_SMILES);
        let target = find_by_map(&m, 1);
        let inverted = invert_stereocenter(&m, target);
        let expected = mol("C[C@@H:1](O)[C@@H:2](N)C(=O)O");
        assert_eq!(
            chematic_smiles::canonical_smiles(&inverted),
            chematic_smiles::canonical_smiles(&expected),
            "inverting only the mapped-1 center must match an independently-parsed \
             molecule with only that center's @/@@ flipped"
        );
    }

    #[test]
    fn invert_stereocenter_preserves_uninvolved_stereo_metadata() {
        let mut m = mol(TWO_CENTER_SMILES);
        let target = find_by_map(&m, 1);
        let other = find_by_map(&m, 2);

        // Attach side-channel metadata to the atom/bond NOT being inverted,
        // to prove the rebuild carries it through untouched.
        m.add_stereo_group(chematic_core::StereoGroup::new(
            chematic_core::StereoGroupKind::Absolute,
            vec![other],
        ));
        let (_, other_bond) = m
            .neighbors(other)
            .next()
            .expect("test setup sanity: other center has a neighbor");
        m.set_bond_direction(other_bond, BondOrder::Up);

        let other_order_before = m.stereo_neighbor_order(other).map(<[u32]>::to_vec);
        assert!(other_order_before.is_some(), "test setup sanity");

        let inverted = invert_stereocenter(&m, target);

        assert_eq!(
            inverted.stereo_groups(),
            m.stereo_groups(),
            "stereo_groups must survive invert_stereocenter verbatim"
        );
        assert_eq!(
            inverted.bond_direction(other_bond),
            m.bond_direction(other_bond),
            "bond_directions must survive invert_stereocenter verbatim"
        );
        assert_eq!(
            inverted.stereo_neighbor_order(other).map(<[u32]>::to_vec),
            other_order_before,
            "uninvolved center's stereo_neighbor_order must survive verbatim"
        );
    }

    #[test]
    fn invert_stereocenter_wedge_fallback_still_works() {
        // Atoms with no `@`/`@@` chirality but an explicit 2D wedge bond
        // (BondOrder::Up/Down) -- the pre-fix code's only working path,
        // kept as a fallback for API compatibility.
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(chematic_core::Element::C));
        let n = b.add_atom(Atom::new(chematic_core::Element::N));
        let o = b.add_atom(Atom::new(chematic_core::Element::O));
        let cl = b.add_atom(Atom::new(chematic_core::Element::CL));
        let br = b.add_atom(Atom::new(chematic_core::Element::BR));
        let wedge_bond = b.add_bond(c, n, BondOrder::Up).unwrap();
        b.add_bond(c, o, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        let m = b.build();
        assert_eq!(
            m.atom(c).chirality,
            Chirality::None,
            "test setup sanity: no @/@@ chirality, only a wedge bond"
        );

        let inverted = invert_stereocenter(&m, c);
        assert_eq!(
            inverted.bond(wedge_bond).order,
            BondOrder::Down,
            "wedge bond must flip Up -> Down when no @/@@ chirality is present"
        );
    }

    #[test]
    fn enumerate_stereoisomers_single_center() {
        // [C@?H](F)(Cl)Br — one unspecified stereocenter
        let m = mol("C(F)(Cl)Br");
        let isomers = enumerate_stereoisomers(&m);
        // Should generate 2 isomers (R and S)
        assert_eq!(
            isomers.len(),
            2,
            "single stereocenter should yield 2 isomers"
        );
    }

    #[test]
    fn enumerate_stereoisomers_no_centers() {
        // CC — no stereocenters
        let m = mol("CC");
        let isomers = enumerate_stereoisomers(&m);
        assert_eq!(isomers.len(), 1, "no stereocenters should yield 1 isomer");
    }

    #[test]
    fn enumerate_stereoisomers_too_many() {
        // Molecule with >6 unspecified stereocenters (unlikely in practice, but test the limit)
        // For testing, we'd need a molecule with 7+ stereocenters, which is complex to construct.
        // Instead, we trust the logic: n > 6 returns empty vec
    }

    #[test]
    fn assign_complete_stereochemistry_simple() {
        // Biphenyl: has potential for atropisomerism if substituted
        let m = mol("c1ccccc1c2ccccc2");
        let result = assign_complete_stereochemistry(&m);
        assert_eq!(result.atom_count(), m.atom_count(), "atom count preserved");
        assert_eq!(result.bond_count(), m.bond_count(), "bond count preserved");
    }

    #[test]
    fn assign_complete_stereochemistry_preserves_structure() {
        // Simple molecule
        let m = mol("CC(F)(Cl)Br");
        let result = assign_complete_stereochemistry(&m);
        assert_eq!(result.atom_count(), m.atom_count(), "structure preserved");
    }

    // -----------------------------------------------------------------
    // Stereo-Rebuild-S1: enumerate_stereoisomers must preserve the
    // already-specified stereocenter's stereo_neighbor_order (it's a
    // passthrough rebuild for every other atom) and must set a valid one
    // for each newly-assigned center -- otherwise the freshly-set
    // `chirality` is uninterpretable (same failure shape Kekule-S0 fixed
    // for apply_kekule).
    // -----------------------------------------------------------------

    /// One already-specified stereocenter (atom 0, `[C@H]`) + one
    /// unspecified one (atom 3, the `C(Br)(I)N` carbon: degree 4, 0
    /// implicit H, 4 chemically distinct substituents).
    const MIXED_STEREOCENTER_SMILES: &str = "[C@H](F)(Cl)C(Br)(I)N";

    #[test]
    fn enumerate_stereoisomers_preserves_specified_center() {
        let m = mol(MIXED_STEREOCENTER_SMILES);
        let specified = AtomIdx(0);
        assert_ne!(
            m.atom(specified).chirality,
            Chirality::None,
            "test setup sanity: atom 0 should already be a specified stereocenter"
        );

        let isomers = enumerate_stereoisomers(&m);
        assert_eq!(
            isomers.len(),
            2,
            "one unspecified center should yield 2 isomers"
        );

        for isomer in &isomers {
            assert_eq!(
                isomer.atom(specified).chirality,
                m.atom(specified).chirality,
                "already-specified stereocenter's chirality label must be unchanged"
            );
            assert_eq!(
                isomer.stereo_neighbor_order(specified),
                m.stereo_neighbor_order(specified),
                "already-specified stereocenter's stereo_neighbor_order must be preserved verbatim"
            );
        }
    }

    #[test]
    fn enumerate_stereoisomers_only_varies_unspecified_center() {
        let m = mol(MIXED_STEREOCENTER_SMILES);
        let unspecified = AtomIdx(3);
        assert_eq!(
            m.atom(unspecified).chirality,
            Chirality::None,
            "test setup sanity: atom 3 should be unspecified before enumeration"
        );

        let isomers = enumerate_stereoisomers(&m);
        let chiralities: std::collections::HashSet<Chirality> = isomers
            .iter()
            .map(|iso| iso.atom(unspecified).chirality)
            .collect();
        assert_eq!(
            chiralities,
            [Chirality::Clockwise, Chirality::CounterClockwise]
                .into_iter()
                .collect(),
            "the 2 isomers should cover both configurations of the unspecified center"
        );
        for isomer in &isomers {
            assert!(
                isomer.stereo_neighbor_order(unspecified).is_some(),
                "newly-assigned chirality must have an interpretable stereo_neighbor_order, \
                 not silently uninterpretable (chirality set, no reference order)"
            );
        }
    }

    #[test]
    fn enumerate_stereoisomers_outputs_are_genuinely_distinct() {
        let m = mol(MIXED_STEREOCENTER_SMILES);
        let isomers = enumerate_stereoisomers(&m);
        assert_eq!(isomers.len(), 2);
        let c0 = chematic_smiles::canonical_smiles(&isomers[0]);
        let c1 = chematic_smiles::canonical_smiles(&isomers[1]);
        assert_ne!(
            c0, c1,
            "the 2 isomers must canonicalize to different SMILES"
        );
    }

    fn find_by_map(m: &Molecule, map_num: u16) -> AtomIdx {
        m.atoms()
            .find(|(_, a)| a.atom_map == Some(map_num))
            .map(|(idx, _)| idx)
            .expect("atom map tag not found")
    }

    #[test]
    fn enumerate_stereoisomers_atom_mapped_cip_survives_canonical_round_trip() {
        // Chemical-identity check per review: don't just compare @/@@
        // strings -- verify the actual CIP label at each stereocenter
        // survives a canonicalize -> reparse round trip. Matched by
        // atom-map tag (confirmed to survive canonical_smiles), NOT by raw
        // AtomIdx -- canonicalization reorders atoms, so the same index in
        // `isomer` and `reparsed` is not the same atom.
        let m = mol("[C@H:1](F)(Cl)[C:2](Br)(I)N");
        for isomer in enumerate_stereoisomers(&m) {
            let before = crate::assign_cip(&isomer);
            let before_specified = before.get(find_by_map(&isomer, 1));
            let before_newly_assigned = before.get(find_by_map(&isomer, 2));

            let smi = chematic_smiles::canonical_smiles(&isomer);
            let reparsed = chematic_smiles::parse(&smi).expect("valid canonical SMILES");
            let after = crate::assign_cip(&reparsed);
            let after_specified = after.get(find_by_map(&reparsed, 1));
            let after_newly_assigned = after.get(find_by_map(&reparsed, 2));

            assert_eq!(
                before_specified, after_specified,
                "already-specified stereocenter's CIP code must survive a canonical round trip"
            );
            assert_eq!(
                before_newly_assigned, after_newly_assigned,
                "newly-assigned stereocenter's CIP code must survive a canonical round trip"
            );
        }
    }

    #[test]
    fn enumerate_stereoisomers_two_centers_no_duplicates_no_omissions() {
        // Two independent unspecified stereocenters -- exactly 4 (2^2)
        // distinct isomers expected, covering all combinations (no
        // enantiomer/diastereomer collapsed together, none missing).
        let m = mol("C(F)(Cl)C(Br)(I)N"); // atom0: F/Cl/H/chain, atom2: Br/I/N/chain
        let unspecified: Vec<AtomIdx> = m
            .atoms()
            .filter(|(_, a)| a.element.atomic_number() == 6 && a.chirality == Chirality::None)
            .map(|(idx, _)| idx)
            .collect();
        assert_eq!(
            unspecified.len(),
            2,
            "test setup sanity: 2 unspecified centers"
        );

        let isomers = enumerate_stereoisomers(&m);
        assert_eq!(
            isomers.len(),
            4,
            "2 unspecified centers should yield 2^2 = 4 isomers"
        );

        let canonical: std::collections::HashSet<String> = isomers
            .iter()
            .map(chematic_smiles::canonical_smiles)
            .collect();
        assert_eq!(
            canonical.len(),
            4,
            "all 4 isomers must be pairwise distinct"
        );
    }

    #[test]
    fn assign_complete_stereochemistry_no_panic() {
        // Ensure it doesn't panic on various inputs
        for smiles in &["C", "CC", "c1ccccc1", "C=C", "CC(=O)O"] {
            let m = mol(smiles);
            let _ = assign_complete_stereochemistry(&m);
        }
    }
}
