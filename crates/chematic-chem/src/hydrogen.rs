//! Explicit hydrogen management.
//!
//! Converts between the compact hydrogen-implicit representation used
//! throughout chematic and a fully-explicit hydrogen graph where each
//! hydrogen is an atom node.

use std::collections::HashMap;

use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Chirality, Element, Molecule, MoleculeBuilder,
    STEREO_H_SENTINEL, implicit_hcount,
};

/// The declared stereo neighbor order for `idx` on `mol`, including the
/// `STEREO_H_SENTINEL` marker for an implicit H if one is declared.
///
/// Own copy of `chematic-3d`'s (private) `stereo_constraints::
/// declared_neighbor_order` -- identical two-branch logic (prefer the
/// parser-recorded `stereo_neighbor_order`, else reconstruct from raw
/// adjacency + the bracket-H insertion heuristic) -- reimplemented here
/// rather than depending on `chematic-3d`, which itself depends on this
/// crate (`chematic-chem`), so the reverse dependency isn't available.
fn declared_neighbor_order(mol: &Molecule, idx: AtomIdx) -> Option<Vec<u32>> {
    if let Some(order) = mol.stereo_neighbor_order(idx) {
        return Some(order.to_vec());
    }
    let atom = mol.atom(idx);
    if atom.chirality == Chirality::None {
        return None;
    }
    let mut neighbors: Vec<u32> = mol.neighbors(idx).map(|(nb, _)| nb.0).collect();
    let has_bracket_h = atom.hydrogen_count.is_some_and(|h| h > 0);
    if has_bracket_h {
        let has_preceding = neighbors.first().map(|&nb| nb < idx.0).unwrap_or(false);
        let h_pos = if has_preceding { 1 } else { 0 };
        neighbors.insert(h_pos, STEREO_H_SENTINEL);
    }
    Some(neighbors)
}

/// Return a new molecule in which every implicit hydrogen is converted to an
/// explicit H atom node bonded to its parent heavy atom.
///
/// The heavy atoms in the returned molecule have `hydrogen_count = Some(0)`,
/// preventing further implicit-H generation.  All original bonds and atom
/// properties are preserved.
///
/// Declared tetrahedral chirality's neighbor order is also preserved: for
/// every stereocenter whose declared order records an implicit H (the
/// `STEREO_H_SENTINEL` marker), that entry is remapped to the newly-added
/// real H atom's index. Without this, the returned molecule's
/// `stereo_neighbor_order` would simply be missing for these atoms (it's a
/// `Molecule`-level side table, not part of `Atom`, so a fresh
/// `MoleculeBuilder` rebuild loses it by default) -- and any downstream
/// consumer falling back to a bracket-H heuristic would silently reconstruct
/// the WRONG order, since this function also sets `hydrogen_count =
/// Some(0)` on every converted atom, defeating the standard "has an implicit
/// H" test that heuristic relies on. Confirmed as a real, previously-latent
/// bug via direct testosterone/cholesterol embedding: without this fix,
/// `chematic-3d`'s stereo repair mechanism can silently accept a
/// wrong-chirality geometry as satisfied after this function runs.
pub fn add_hydrogens(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Copy heavy atoms with hydrogen_count set to Some(0).
    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();
        atom.hydrogen_count = Some(0);
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }

    // Copy all original bonds.
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        if let (Some(&na), Some(&nb)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(na, nb, bond.order);
        }
    }

    // Heavy-atom indices are unchanged by this function (copied 1:1, above),
    // so every existing stereo_neighbor_order entry referencing only real
    // heavy atoms is already correct verbatim. Entries with a sentinel are
    // patched below once each atom's new explicit H index is known.
    builder.copy_stereo_from(mol);

    // Add explicit H atoms for each implicit hydrogen, and fix up declared
    // stereo order for any stereocenter that gains one.
    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let h_count = implicit_hcount(mol, old_idx);
        if h_count == 0 {
            continue;
        }
        let heavy_new = remap[&old_idx];
        let mut new_h_atoms: Vec<AtomIdx> = Vec::with_capacity(h_count as usize);
        for _ in 0..h_count {
            let h_atom = Atom::new(Element::H);
            let h_new = builder.add_atom(h_atom);
            let _ = builder.add_bond(heavy_new, h_new, BondOrder::Single);
            new_h_atoms.push(h_new);
        }

        let atom = mol.atom(old_idx);
        if atom.chirality == Chirality::None {
            continue;
        }
        let Some(order) = declared_neighbor_order(mol, old_idx) else {
            continue;
        };
        // A declared tetrahedral center carries at most one sentinel, and
        // only when it has exactly one implicit H (`TetrahedralConstraint`'s
        // own invariant in chematic-3d) -- anything else here means this
        // atom isn't actually a simple tetrahedral stereocenter as declared;
        // leave its (already bulk-copied) order alone rather than guess.
        if new_h_atoms.len() != 1 {
            continue;
        }
        let new_h_idx = new_h_atoms[0].0;
        let new_order: Vec<u32> = order
            .into_iter()
            .map(|v| {
                if v == STEREO_H_SENTINEL {
                    new_h_idx
                } else {
                    remap[&AtomIdx(v)].0
                }
            })
            .collect();
        builder.set_stereo_neighbor_order(heavy_new, new_order);
    }

    builder.build()
}

/// Return a new molecule in which explicit H atom nodes are removed and their
/// bonds are converted back to implicit hydrogens.
///
/// Only explicit hydrogen atoms (nodes with `element == H`) are removed.
/// Chirality annotations and other atom properties are preserved.
///
/// Heavy atoms that had explicit H neighbors will have `hydrogen_count` set
/// to `None` so that implicit H is recomputed from valence.
pub fn remove_hydrogens(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        if mol.atom(old_idx).element == Element::H {
            continue;
        }
        let mut atom = mol.atom(old_idx).clone();
        // Restore implicit H computation for atoms that had explicit H set.
        if atom.hydrogen_count == Some(0) {
            atom.hydrogen_count = None;
        }
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }

    // Copy heavy–heavy bonds only.
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        let a1_is_h = mol.atom(bond.atom1).element == Element::H;
        let a2_is_h = mol.atom(bond.atom2).element == Element::H;
        if a1_is_h || a2_is_h {
            continue;
        }
        if let (Some(&na), Some(&nb)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(na, nb, bond.order);
        }
    }

    builder.build()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn add_h_methane_atom_count() {
        // C → 1 C + 4 H = 5 atoms
        let m = add_hydrogens(&mol("C"));
        assert_eq!(m.atom_count(), 5, "methane + H should have 5 atoms");
    }

    #[test]
    fn add_h_methane_bond_count() {
        // 4 C-H bonds
        let m = add_hydrogens(&mol("C"));
        assert_eq!(m.bond_count(), 4, "methane + H should have 4 bonds");
    }

    #[test]
    fn add_h_ethane() {
        // CC → 2 C + 6 H = 8 atoms, 1 C-C + 6 C-H = 7 bonds
        let m = add_hydrogens(&mol("CC"));
        assert_eq!(m.atom_count(), 8, "ethane + H atoms");
        assert_eq!(m.bond_count(), 7, "ethane + H bonds");
    }

    #[test]
    fn add_h_benzene() {
        // c1ccccc1 → 6 C + 6 H = 12 atoms, 6 ring + 6 C-H = 12 bonds
        let m = add_hydrogens(&mol("c1ccccc1"));
        assert_eq!(m.atom_count(), 12, "benzene + H atoms");
        assert_eq!(m.bond_count(), 12, "benzene + H bonds");
    }

    #[test]
    fn add_remove_roundtrip_ethanol() {
        let orig = mol("CCO");
        let with_h = add_hydrogens(&orig);
        let restored = remove_hydrogens(&with_h);
        // Heavy-atom count and bond count should match original.
        assert_eq!(
            restored.atom_count(),
            orig.atom_count(),
            "roundtrip atom count"
        );
        assert_eq!(
            restored.bond_count(),
            orig.bond_count(),
            "roundtrip bond count"
        );
    }

    #[test]
    fn add_remove_roundtrip_aspirin() {
        let orig = mol("CC(=O)Oc1ccccc1C(=O)O");
        let with_h = add_hydrogens(&orig);
        let restored = remove_hydrogens(&with_h);
        assert_eq!(restored.atom_count(), orig.atom_count());
        assert_eq!(restored.bond_count(), orig.bond_count());
    }

    #[test]
    fn remove_h_no_h_atoms_unchanged() {
        // A molecule with no explicit H nodes: remove_hydrogens should be a no-op.
        let orig = mol("CC");
        let result = remove_hydrogens(&orig);
        assert_eq!(result.atom_count(), 2);
        assert_eq!(result.bond_count(), 1);
    }

    #[test]
    fn add_h_water() {
        // O → 1 O + 2 H = 3 atoms, 2 bonds
        let m = add_hydrogens(&mol("O"));
        assert_eq!(m.atom_count(), 3);
        assert_eq!(m.bond_count(), 2);
    }

    #[test]
    fn add_h_preserves_element_distribution() {
        // Aspirin: 9 C + 4 O = 13 heavy; 8 H added → 21 total
        let orig = mol("CC(=O)Oc1ccccc1C(=O)O");
        let with_h = add_hydrogens(&orig);
        let h_count = with_h
            .atoms()
            .filter(|(_, a)| a.element == Element::H)
            .count();
        assert_eq!(h_count, 8, "aspirin should gain 8 H atoms (C9H8O4)");
    }

    // ─── Declared-chirality preservation (issue #291) ──────────────────────

    #[test]
    fn add_h_implicit_h_stereocenter_remaps_sentinel_to_new_h_atom() {
        // N[C@@H](C)C(=O)O (L-alanine): atom 1 is the stereocenter, declared
        // order [N(0), STEREO_H_SENTINEL, C(2), C(3)] at parse time.
        let orig = mol("N[C@@H](C)C(=O)O");
        let stereocenter = AtomIdx(1);
        let orig_order = orig
            .stereo_neighbor_order(stereocenter)
            .expect("parser must record stereo order for a declared @@ center")
            .to_vec();
        assert!(
            orig_order.contains(&STEREO_H_SENTINEL),
            "original order must record the implicit H: {orig_order:?}"
        );

        let with_h = add_hydrogens(&orig);
        assert_eq!(
            with_h.atom(stereocenter).hydrogen_count,
            Some(0),
            "heavy atom index must be unchanged by add_hydrogens"
        );

        let new_order = with_h
            .stereo_neighbor_order(stereocenter)
            .expect("stereo order must survive add_hydrogens, not just get dropped")
            .to_vec();
        assert!(
            !new_order.contains(&STEREO_H_SENTINEL),
            "sentinel must be replaced by a real atom index: {new_order:?}"
        );
        assert_eq!(
            new_order.len(),
            orig_order.len(),
            "remapping must not change the neighbor count"
        );

        // The sentinel's replacement must be the new H atom actually bonded
        // to the stereocenter -- not just any new atom.
        let sentinel_pos = orig_order
            .iter()
            .position(|&v| v == STEREO_H_SENTINEL)
            .unwrap();
        let new_h_idx = AtomIdx(new_order[sentinel_pos]);
        assert_eq!(
            with_h.atom(new_h_idx).element,
            Element::H,
            "sentinel must be replaced by an H atom, got {:?}",
            with_h.atom(new_h_idx).element
        );
        assert!(
            with_h
                .neighbors(stereocenter)
                .any(|(nb, _)| nb == new_h_idx),
            "the substituted H atom must actually be bonded to the stereocenter"
        );

        // Every non-sentinel slot must be untouched (same real neighbor
        // index -- heavy atoms don't move).
        for (i, &v) in orig_order.iter().enumerate() {
            if v != STEREO_H_SENTINEL {
                assert_eq!(
                    new_order[i], v,
                    "non-H neighbor at slot {i} must be unchanged"
                );
            }
        }
    }

    #[test]
    fn add_h_quaternary_stereocenter_order_unchanged() {
        // [C@](F)(Cl)(Br)I: no implicit H, no sentinel -- add_hydrogens adds
        // nothing to this atom, so its declared order must be preserved
        // verbatim (already correct via the bulk copy, not the H-remap path).
        let orig = mol("[C@](F)(Cl)(Br)I");
        let stereocenter = AtomIdx(0);
        let orig_order = orig
            .stereo_neighbor_order(stereocenter)
            .expect("quaternary center must have a declared order")
            .to_vec();
        assert!(!orig_order.contains(&STEREO_H_SENTINEL));

        let with_h = add_hydrogens(&orig);
        let new_order = with_h
            .stereo_neighbor_order(stereocenter)
            .expect("order must survive add_hydrogens even with zero implicit H")
            .to_vec();
        assert_eq!(
            new_order, orig_order,
            "no-implicit-H center must be untouched"
        );
    }

    #[test]
    fn add_h_multi_stereocenter_molecule_all_orders_correct() {
        // L-threonine: two implicit-H stereocenters in the same molecule --
        // confirms the fix handles more than one sentinel-bearing atom
        // independently and correctly in a single call.
        let orig = mol("C[C@H](O)[C@@H](N)C(=O)O");
        let centers: Vec<AtomIdx> = (0..orig.atom_count() as u32)
            .map(AtomIdx)
            .filter(|&idx| orig.atom(idx).chirality != Chirality::None)
            .collect();
        assert_eq!(centers.len(), 2, "threonine has 2 declared stereocenters");

        let with_h = add_hydrogens(&orig);
        for &center in &centers {
            let orig_order = orig.stereo_neighbor_order(center).unwrap().to_vec();
            let new_order = with_h
                .stereo_neighbor_order(center)
                .unwrap_or_else(|| panic!("order for atom {center:?} must survive"))
                .to_vec();
            assert!(
                !new_order.contains(&STEREO_H_SENTINEL),
                "atom {center:?}: sentinel must be resolved, got {new_order:?}"
            );
            let sentinel_pos = orig_order.iter().position(|&v| v == STEREO_H_SENTINEL);
            if let Some(pos) = sentinel_pos {
                let new_h_idx = AtomIdx(new_order[pos]);
                assert_eq!(with_h.atom(new_h_idx).element, Element::H);
                assert!(with_h.neighbors(center).any(|(nb, _)| nb == new_h_idx));
            } else {
                assert_eq!(new_order, orig_order);
            }
        }
    }
}
