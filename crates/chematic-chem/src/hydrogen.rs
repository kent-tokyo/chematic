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

/// Whether atom `a` is a *removable* explicit hydrogen: element `H` with no
/// isotope specified. An isotope-labeled hydrogen (`[2H]` deuterium, `[3H]`
/// tritium, ...) is real, distinguishing chemical information -- akin to
/// charge or element identity, not a notation choice -- so it is never
/// removed by this function regardless of `element == H`. See
/// [`remove_hydrogens`]'s own doc comment for the full rationale and the
/// regression this guards (issue: isotope labels silently destroyed by
/// unconditional H-node removal).
fn is_removable_explicit_h(a: &Atom) -> bool {
    a.element == Element::H && a.isotope.is_none()
}

/// Return a new molecule in which removable explicit H atom nodes (see
/// [`is_removable_explicit_h`]) are removed and their bonds are converted
/// back to implicit hydrogens.
///
/// Only *non-isotopic* explicit hydrogen atoms are removed -- an
/// isotope-labeled H (`[2H]`, `[3H]`, ...) is kept as an explicit atom node,
/// exactly like any other heavy atom, since collapsing it into an ordinary
/// atom's opaque `hydrogen_count` would silently discard the isotope label
/// (there is no way to record "N implicit hydrogens, one of which is
/// deuterium" in that compact representation). A heavy atom that retains an
/// isotopic-H neighbor keeps that bond untouched; only bonds to a removed
/// (non-isotopic) H are dropped. Chirality annotations and other atom
/// properties are preserved.
///
/// Heavy atoms that had *removable* explicit H neighbors will have
/// `hydrogen_count` reset to `None` so implicit H is recomputed from
/// valence -- correct even when the atom also keeps an explicit isotopic-H
/// neighbor, since valence inference counts every bonded neighbor (isotopic
/// or not) via its bond order, not by element identity (see
/// `chematic_core::valence::valence_inferred_hcount`).
///
/// Declared tetrahedral/square-planar chirality's neighbor order, and E/Z
/// bond direction, are also restored: this is the exact inverse of what
/// [`add_hydrogens`] already does for the opposite direction (see its own
/// doc comment for why this matters -- `stereo_neighbor_order`/
/// `bond_directions` are `Molecule`-level side tables, not part of `Atom`/
/// `BondEntry`, so a fresh `MoleculeBuilder` rebuild loses them by default,
/// even on a call that removes nothing at all). For every stereocenter
/// that survives, this atom's *original* declared order (parser-recorded,
/// or the same bracket-H heuristic reconstruction [`add_hydrogens`] uses
/// when nothing was recorded) is remapped: an entry pointing at a removed
/// (non-isotopic) H becomes the `STEREO_H_SENTINEL` marker again (mirroring
/// the sentinel this function's own H-removal just re-created conceptually
/// -- an entry pointing at a *kept* isotopic H, or at any surviving heavy
/// atom, is remapped to its new index instead), and every surviving bond
/// that carried a declared E/Z direction keeps it at its new bond index.
///
/// Without this, `chematic-smiles`'s canonical writer -- which requires
/// `stereo_neighbor_order` to safely reinterpret a stored `@`/`@@` tag
/// against a *different* (e.g. canonically-reordered) neighbor sequence --
/// has no order to consult for any stereocenter that survived a
/// `remove_hydrogens` call, and silently passes the raw stored tag through
/// unchanged against whatever new order it picks (`corrected_chirality`'s
/// `stereo_neighbor_order(atom).is_none()` fallback). That tag's meaning is
/// only valid relative to the order it was declared against; reusing it
/// against an unrelated order can silently encode the *wrong*
/// configuration. This was a real, confirmed bug: re-canonicalizing an
/// already-canonical SMILES (produced via `standardize` with
/// `remove_explicit_h: true`, this crate's own default) could flip a
/// declared stereocenter to its mirror image on some symmetric-ranking-
/// ambiguous molecules, independently confirmed via RDKit InChIKey
/// divergence on a real-world eMolecules corpus (9.47M compounds; 289 of
/// 290 confirmed InChIKey mismatches are resolved by this fix -- the one
/// residual case is a coupled/shared-bond E/Z system with a confirmed
/// *different* root cause, independent of `remove_hydrogens` entirely;
/// see issue #390).
pub fn remove_hydrogens(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        if is_removable_explicit_h(mol.atom(old_idx)) {
            continue;
        }
        let mut atom = mol.atom(old_idx).clone();
        // Restore implicit H computation, but only for atoms that actually
        // had a removable explicit H *atom* neighbor this call is dropping --
        // `add_hydrogens` sets `hydrogen_count = Some(0)` on every atom it
        // converts specifically to pair with the new explicit H atom
        // neighbors it adds, so "had a removable-H neighbor" is the correct
        // signal to undo that, not "stored H count happens to be 0". An atom
        // whose H count is genuinely, deliberately 0 (e.g. a dative-bonded
        // `[O+]` disconnected from a metal by an earlier standardize stage,
        // never had, and isn't gaining, any H atom neighbor) must not have
        // that 0 reinterpreted via valence-based inference: doing so
        // silently invents an implicit hydrogen from a valence change that
        // happened elsewhere in the pipeline, and does so too late for
        // `neutralize_charges` (which already ran) to neutralize the
        // resulting charge -- a real, confirmed idempotency bug (issue #403).
        if atom.hydrogen_count == Some(0)
            && mol
                .neighbors(old_idx)
                .any(|(nb, _)| is_removable_explicit_h(mol.atom(nb)))
        {
            atom.hydrogen_count = None;
        }
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }

    // Copy bonds, dropping only those touching a removed (non-isotopic) H,
    // and tracking each surviving bond's new index so `bond_directions`
    // (E/Z `/`/`\` markers, a separate Molecule-level side table keyed by
    // bond index) can be remapped below -- the same "fresh MoleculeBuilder
    // loses any side table nothing explicitly re-populates" issue
    // `stereo_neighbor_order` has, just for bonds instead of atoms.
    // `Molecule::with_atom_removed`/`with_bond_removed` already carry an
    // equivalent remap for exactly this reason (see their own doc comments
    // on `bond_directions` "silent loss" -- the same failure category,
    // just not yet closed for this function).
    let mut bond_remap: HashMap<BondIdx, BondIdx> = HashMap::new();
    for i in 0..mol.bond_count() {
        let old_bidx = BondIdx(i as u32);
        let bond = mol.bond(old_bidx);
        let a1_removed = is_removable_explicit_h(mol.atom(bond.atom1));
        let a2_removed = is_removable_explicit_h(mol.atom(bond.atom2));
        if a1_removed || a2_removed {
            continue;
        }
        if let (Some(&na), Some(&nb)) = (remap.get(&bond.atom1), remap.get(&bond.atom2))
            && let Ok(new_bidx) = builder.add_bond(na, nb, bond.order)
        {
            bond_remap.insert(old_bidx, new_bidx);
        }
    }

    // Restore bond_directions (E/Z `/`/`\` markers) for every surviving
    // bond that carried one -- without this, `chematic-smiles`'s canonical
    // writer has no declared direction to reinterpret against a
    // canonically-reordered traversal for a double bond whose ends are
    // organic-subset atoms (no bracket, so nothing lives on `Atom` itself
    // to carry this), and silently drops or guesses E/Z geometry the same
    // way `corrected_chirality` does for tetrahedral centers without
    // `stereo_neighbor_order` (see that fix above).
    for (&old_bidx, &new_bidx) in &bond_remap {
        if let Some(direction) = mol.bond_direction(old_bidx) {
            builder.set_bond_direction(new_bidx, direction);
        }
    }

    // Restore stereo_neighbor_order for every surviving stereocenter --
    // the exact inverse of add_hydrogens' own sentinel-remap (above).
    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let Some(&new_idx) = remap.get(&old_idx) else {
            continue; // this atom itself was removed (an H can't be chiral anyway)
        };
        if mol.atom(old_idx).chirality == Chirality::None {
            continue;
        }
        let Some(order) = declared_neighbor_order(mol, old_idx) else {
            continue;
        };
        let new_order: Vec<u32> = order
            .into_iter()
            .map(|v| {
                if v == STEREO_H_SENTINEL {
                    // Already an implicit-H slot in the original -- stays one.
                    STEREO_H_SENTINEL
                } else {
                    match remap.get(&AtomIdx(v)) {
                        // Surviving neighbor (heavy atom, or a kept
                        // isotopic H): follow it to its new index.
                        Some(&nb_new) => nb_new.0,
                        // This neighbor was itself a removed explicit H --
                        // it's now an implicit hydrogen on this atom again.
                        None => STEREO_H_SENTINEL,
                    }
                }
            })
            .collect();
        builder.set_stereo_neighbor_order(new_idx, new_order);
    }

    builder.build()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::{canonical_smiles, parse};

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    // ─── Isotopic-hydrogen preservation ────────────────────────────────────

    fn h_isotopes(m: &Molecule) -> Vec<Option<u16>> {
        let mut v: Vec<Option<u16>> = m
            .atoms()
            .filter(|(_, a)| a.element == Element::H)
            .map(|(_, a)| a.isotope)
            .collect();
        v.sort();
        v
    }

    #[test]
    fn remove_hydrogens_keeps_fully_deuterated_methane() {
        // CD4: 4 explicit D atom nodes, no plain H at all.
        let orig = mol("[2H]C([2H])([2H])[2H]");
        assert_eq!(orig.atom_count(), 5);
        let result = remove_hydrogens(&orig);
        assert_eq!(
            result.atom_count(),
            5,
            "all 4 deuteriums must be kept as explicit atom nodes, not folded away"
        );
        assert_eq!(
            h_isotopes(&result),
            vec![Some(2), Some(2), Some(2), Some(2)]
        );
        let canon = canonical_smiles(&result);
        assert_eq!(
            canon.matches("[2H]").count(),
            4,
            "canonical SMILES must show all 4 deuteriums: {canon}"
        );
        assert!(
            !canon.contains("[H]"),
            "must not have degraded to plain explicit H: {canon}"
        );
    }

    #[test]
    fn remove_hydrogens_keeps_tritium() {
        // Tritiated water analog: O bonded to one T and (implicitly) one
        // ordinary H it never had explicitly -- OT.
        let orig = mol("[3H]O");
        let result = remove_hydrogens(&orig);
        assert_eq!(
            result.atom_count(),
            2,
            "the tritium atom must be kept, not removed"
        );
        assert_eq!(h_isotopes(&result), vec![Some(3)]);
        let canon = canonical_smiles(&result);
        assert!(canon.contains("[3H]"), "canonical SMILES lost T: {canon}");
    }

    #[test]
    fn remove_hydrogens_mixed_deuterium_and_plain_hydrogen() {
        // One explicit D + three explicit plain H on the same carbon
        // (CHD3-style methane, fully explicit form). The plain H's must be
        // removed (folded into hydrogen_count); the D must survive as an
        // explicit neighbor; the carbon's total valence (4 bonds either way)
        // must stay correct.
        let orig = mol("[2H]C([H])([H])[H]");
        assert_eq!(orig.atom_count(), 5);
        let result = remove_hydrogens(&orig);
        assert_eq!(
            result.atom_count(),
            2,
            "C + the single kept deuterium; the 3 plain H must be removed"
        );
        assert_eq!(h_isotopes(&result), vec![Some(2)]);

        let carbon = result
            .atoms()
            .find(|(_, a)| a.element == Element::C)
            .map(|(idx, _)| idx)
            .expect("carbon must survive");
        // Valence must be satisfied by 1 explicit D bond + 3 recomputed
        // implicit H -- not 4 implicit H (which would silently duplicate
        // the D's own bond) and not 0 (which would under-saturate carbon).
        assert_eq!(
            chematic_core::implicit_hcount(&result, carbon),
            3,
            "carbon must recompute exactly 3 implicit H, correctly accounting for \
             the kept deuterium's own bond -- not 4 (double-counting D) or fewer"
        );

        let canon = canonical_smiles(&result);
        assert_eq!(canon.matches("[2H]").count(), 1, "{canon}");
        assert!(!canon.contains("[H]"), "plain H must not survive: {canon}");
        // total_formula() is deliberately isotope-blind (it keys purely by
        // element symbol, per its own doc comment) -- 4 H-family atoms
        // total either way (1 kept D + 3 recomputed implicit protium),
        // matching the original explicit-everything form's own H count.
        assert_eq!(result.total_formula(), orig.total_formula());
        // formula_with_isotopes() *is* isotope-aware and must show the kept D.
        assert!(
            result.formula_with_isotopes().contains("2H"),
            "isotope-aware formula must show the surviving deuterium: {}",
            result.formula_with_isotopes()
        );
    }

    #[test]
    fn remove_hydrogens_heavy_isotopes_untouched_13c_14c_15n_18o() {
        // None of these carry an isotope-tagged H at all -- confirms the
        // fix's H-specific guard doesn't regress the already-correct
        // heavy-atom-isotope behavior (unaffected before and after, since
        // `element == H` was always false for these atoms).
        for smi in ["[13CH4]", "[14CH4]", "[15NH3]", "[18OH2]"] {
            let orig = mol(smi);
            let before_isotope = orig
                .atoms()
                .find(|(_, a)| a.element != Element::H)
                .and_then(|(_, a)| a.isotope);
            let result = remove_hydrogens(&orig);
            let after_isotope = result
                .atoms()
                .find(|(_, a)| a.element != Element::H)
                .and_then(|(_, a)| a.isotope);
            assert_eq!(
                before_isotope, after_isotope,
                "{smi}: heavy-atom isotope must be untouched by remove_hydrogens"
            );
            assert!(
                result.atoms().all(|(_, a)| a.element != Element::H),
                "{smi}: no explicit H atom node exists in bracket-implicit-H \
                 form to begin with, so nothing should remain either way"
            );
        }
    }

    #[test]
    fn remove_hydrogens_preserves_isotope_on_a_tetrahedral_stereocenter() {
        // Alpha-deuterated alanine analog: a declared tetrahedral
        // stereocenter with an explicit deuterium substituent (a real,
        // common labeling pattern in metabolic-stability medicinal
        // chemistry). Scope of this fix: the D atom and the stereocenter's
        // atom-level `chirality` tag must both survive structurally sound
        // (right atom/bond count, right element set) -- this test does NOT
        // assert the exact `@`/`@@` parity is textually preserved, since
        // `remove_hydrogens` does not restore the parser-recorded
        // `stereo_neighbor_order` side table (a separate, pre-existing gap,
        // out of scope for this isotope-only fix -- see the module doc
        // comment on `remove_hydrogens`).
        let orig = mol("N[C@@H]([2H])C(=O)O");
        let stereocenter = orig
            .atoms()
            .find(|(_, a)| a.chirality != Chirality::None)
            .map(|(idx, _)| idx)
            .expect("a declared stereocenter must exist");
        assert_ne!(orig.atom(stereocenter).chirality, Chirality::None);

        let result = remove_hydrogens(&orig);
        assert_eq!(h_isotopes(&result), vec![Some(2)], "the D must survive");
        assert_eq!(
            result.atom(stereocenter).chirality,
            orig.atom(stereocenter).chirality,
            "the stereocenter's own chirality tag (an Atom-level field, not \
             the separate neighbor-order side table) must be copied verbatim"
        );
        assert!(
            result
                .neighbors(stereocenter)
                .any(|(nb, _)| result.atom(nb).element == Element::H
                    && result.atom(nb).isotope == Some(2)),
            "the kept deuterium must still be bonded to the stereocenter"
        );
        // Structural soundness: same heavy-atom formula, same total H-family
        // count (1 explicit D + whatever protium the other atoms need).
        assert_eq!(result.total_formula(), orig.total_formula());
    }

    #[test]
    fn standardize_round_trip_preserves_isotope() {
        // The exact policy RENKIN's stock-identity pipeline uses
        // (`remove_explicit_h: true`), run through the real, full
        // `standardize()` entry point -- not just the bare `remove_hydrogens`
        // helper -- to confirm the fix actually reaches consumers through
        // that path, not only when called directly.
        use crate::standardize::{StandardizeOptions, ZwitterionHandling, standardize};
        let opts = StandardizeOptions {
            canonical_tautomer: false,
            neutralize_charges: false,
            remove_explicit_h: true,
            largest_fragment_only: false,
            zwitterion_handling: ZwitterionHandling::Keep,
        };
        let orig = mol("[2H]C([2H])([2H])[2H]");
        let result = standardize(&orig, &opts);
        assert_eq!(
            h_isotopes(&result),
            vec![Some(2), Some(2), Some(2), Some(2)]
        );
        let canon = canonical_smiles(&result);
        assert_eq!(canon.matches("[2H]").count(), 4, "{canon}");
    }

    #[test]
    fn canonical_round_trip_preserves_isotope() {
        // canon(parse(canon(parse(s)))) must still carry the isotope --
        // the exact re-canonicalization scenario (RENKIN's `renkin doctor
        // stock` reimport_idempotency check) that originally surfaced this
        // bug: a real stock file's already-canonical line, canonicalized
        // again, silently losing D/T on the very first pass, let alone a
        // second one.
        use crate::standardize::{StandardizeOptions, ZwitterionHandling, standardize};
        let opts = StandardizeOptions {
            canonical_tautomer: false,
            neutralize_charges: false,
            remove_explicit_h: true,
            largest_fragment_only: false,
            zwitterion_handling: ZwitterionHandling::Keep,
        };
        let stock_identity =
            |s: &str| -> String { canonical_smiles(&standardize(&parse(s).unwrap(), &opts)) };

        for smi in [
            "[2H]C([2H])([2H])[2H]",
            "[3H]O",
            "[2H]C([2H])([2H])C(=O)N[C@H]1CCc2cc(OC)c(OC)c(OC)c2-c2ccc(OC)c(=O)cc21",
        ] {
            let once = stock_identity(smi);
            let twice = stock_identity(&once);
            assert_eq!(
                once, twice,
                "re-canonicalizing an already-canonical isotope-labeled SMILES \
                 must be a no-op: {smi}"
            );
            assert!(
                once.contains("[2H]") || once.contains("[3H]"),
                "isotope must survive even one full standardize+canonicalize pass: \
                 {smi} -> {once}"
            );
        }
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

    // ─── remove_hydrogens: stereo_neighbor_order/bond_directions restoration
    // (canonical round-trip non-idempotency fix) ───────────────────────────

    use crate::cip::assign_cip;
    use crate::standardize::{StandardizeOptions, ZwitterionHandling, standardize};

    /// RENKIN's own stock-identity policy (`remove_explicit_h: true`) --
    /// the exact real pipeline that surfaced this bug via a `renkin doctor
    /// stock reimport_idempotency` FAIL on a 9.48M-compound corpus.
    fn stock_identity_opts() -> StandardizeOptions {
        StandardizeOptions {
            canonical_tautomer: false,
            neutralize_charges: false,
            remove_explicit_h: true,
            largest_fragment_only: false,
            zwitterion_handling: ZwitterionHandling::Keep,
        }
    }

    fn stock_identity(s: &str) -> String {
        canonical_smiles(&standardize(&mol(s), &stock_identity_opts()))
    }

    #[test]
    fn remove_hydrogens_restores_stereo_neighbor_order_with_nothing_removed() {
        // A stereocenter written with bracket-implicit H ([C@H]) -- there is
        // no separate explicit H atom node here at all, so remove_hydrogens
        // has literally nothing to remove. Before this fix, the function
        // still built a brand-new MoleculeBuilder without ever copying
        // stereo_neighbor_order, so this table was silently wiped even in
        // this complete no-op case.
        let orig = mol("N[C@@H](C)C(=O)O");
        let stereocenter = orig
            .atoms()
            .find(|(_, a)| a.chirality != Chirality::None)
            .map(|(idx, _)| idx)
            .unwrap();
        let result = remove_hydrogens(&orig);
        assert!(
            result.stereo_neighbor_order(stereocenter).is_some(),
            "stereo_neighbor_order must survive a remove_hydrogens call that \
             removed nothing at all"
        );
    }

    #[test]
    fn remove_hydrogens_converts_removed_h_neighbor_back_to_sentinel() {
        // A stereocenter whose declared order records a REAL (non-sentinel)
        // explicit-H atom index, not a bracket-implicit sentinel -- the
        // parser itself always collapses a *written* `[C@@H]([H])`-style
        // explicit H straight into the sentinel form (confirmed directly:
        // `mol("N[C@@H]([H])C(=O)O")`'s own `stereo_neighbor_order` already
        // contains `STEREO_H_SENTINEL`, not a real index), so the only way
        // to reach a REAL index in this table is via `add_hydrogens`'s own
        // sentinel -> new-real-H-atom substitution. Chaining
        // `remove_hydrogens(add_hydrogens(orig))` is therefore the actual
        // round trip this fix must get right -- the exact inverse of
        // add_h_implicit_h_stereocenter_remaps_sentinel_to_new_h_atom above.
        let orig = mol("N[C@@H](C)C(=O)O");
        let stereocenter = orig
            .atoms()
            .find(|(_, a)| a.chirality != Chirality::None)
            .map(|(idx, _)| idx)
            .unwrap();
        let orig_order = orig.stereo_neighbor_order(stereocenter).unwrap().to_vec();
        assert!(
            orig_order.contains(&STEREO_H_SENTINEL),
            "sanity: the original bracket-implicit form must start as a sentinel: {orig_order:?}"
        );

        let with_h = add_hydrogens(&orig);
        let with_h_order = with_h.stereo_neighbor_order(stereocenter).unwrap().to_vec();
        assert!(
            !with_h_order.contains(&STEREO_H_SENTINEL),
            "sanity: add_hydrogens must have substituted a real H atom index: {with_h_order:?}"
        );

        let result = remove_hydrogens(&with_h);
        let new_order = result
            .stereo_neighbor_order(stereocenter)
            .expect("order must survive the round trip")
            .to_vec();
        assert!(
            new_order.contains(&STEREO_H_SENTINEL),
            "the removed H's slot must become a sentinel again: {new_order:?}"
        );
        assert_eq!(
            new_order, orig_order,
            "round trip must reproduce the original order exactly"
        );
    }

    #[test]
    fn canonical_round_trip_tetrahedral_witness() {
        // Minimized tetrahedral witness (a fluorenylmethyl-carbamate-type
        // stereocenter) that flipped to its mirror image on a second
        // canonicalization pass before this fix.
        let witness = "c1c2c(C(c3c2cccc3)[C@H](OC(Cl)=O)C)ccc1";
        let once = stock_identity(witness);
        let twice = stock_identity(&once);
        assert_eq!(once, twice, "canon(x) must equal canon(canon(x))");
    }

    #[test]
    fn canonical_round_trip_ez_witness() {
        // Minimized E/Z witness (fumaric acid, a genuine textbook E-alkene).
        let witness = "OC(=O)/C=C/C(=O)O";
        let once = stock_identity(witness);
        let twice = stock_identity(&once);
        assert_eq!(once, twice, "canon(x) must equal canon(canon(x))");
    }

    #[test]
    fn cip_assignment_preserved_across_canonical_round_trip() {
        let witness = "c1c2c(C(c3c2cccc3)[C@H](OC(Cl)=O)C)ccc1";
        let once = mol(&stock_identity(witness));
        let twice = mol(&stock_identity(&stock_identity(witness)));

        let mut once_codes: Vec<_> = assign_cip(&once)
            .assignments
            .iter()
            .map(|(_, c)| *c)
            .collect();
        let mut twice_codes: Vec<_> = assign_cip(&twice)
            .assignments
            .iter()
            .map(|(_, c)| *c)
            .collect();
        once_codes.sort_by_key(|c| format!("{c:?}"));
        twice_codes.sort_by_key(|c| format!("{c:?}"));
        assert_eq!(
            once_codes, twice_codes,
            "the same set of CIP descriptors must survive a second \
             canonicalization pass, not just the raw @/@@ token"
        );
        assert!(
            !once_codes.is_empty(),
            "witness must have an assignable CIP center"
        );
    }

    #[test]
    fn mirror_image_gets_a_different_canonical_identity() {
        // Sanity check in the opposite direction from idempotency: fixing
        // "don't accidentally flip stereo" must not degrade into "don't
        // distinguish stereo at all". @ and @@ on the same skeleton must
        // still canonicalize to two different strings.
        let r = stock_identity("N[C@H](C)C(=O)O");
        let s = stock_identity("N[C@@H](C)C(=O)O");
        assert_ne!(
            r, s,
            "enantiomers must still get distinct canonical identities"
        );
    }

    #[test]
    fn atom_order_permutation_invariance_tetrahedral() {
        // Two textually different, independently-written SMILES for the
        // exact same real molecule and the exact same real configuration
        // (verified by hand-checking the substituent priority order, not
        // just asserting a chosen answer) -- both must canonicalize to the
        // same string, not just be self-idempotent individually.
        let a = stock_identity("N[C@@H](C)C(=O)O"); // L-alanine, written from N
        let b = stock_identity("C[C@H](N)C(=O)O"); // same molecule, written from the methyl
        assert_eq!(
            a, b,
            "two valid spellings of the same molecule/configuration must converge"
        );
    }

    #[test]
    fn regression_boc_protected_amine_tetrahedral_stable() {
        // Boc (tert-butyloxycarbonyl) protecting group next to a
        // stereocenter -- a very common real-world substructure, deliberately
        // included per this fix's own regression-test requirement.
        let witness = "CC(C)(C)OC(=O)N[C@@H](C)C(=O)O";
        let once = stock_identity(witness);
        let twice = stock_identity(&once);
        assert_eq!(once, twice);
    }

    #[test]
    fn regression_fused_ring_system_tetrahedral_stable() {
        // A fused bicyclic (indane-type) system with an adjacent
        // stereocenter -- the ring-fusion/symmetric-ranking-ambiguous shape
        // this whole investigation's corpus scan found the bug concentrated in.
        let witness = "c1ccc2c(c1)CC[C@H]2N";
        let once = stock_identity(witness);
        let twice = stock_identity(&once);
        assert_eq!(once, twice);
    }

    #[test]
    fn simple_isolated_ez_bonds_are_correct_and_stable() {
        // A spread of simple, non-coupled E/Z cases -- all confirmed via
        // this fix's own investigation to already round-trip correctly
        // (both idempotent AND matching an independent InChIKey check),
        // pinned here as an explicit regression net.
        for (e_form, z_form) in [
            ("C/C=C/C", "C/C=C\\C"),
            ("OC(=O)/C=C/C(=O)O", "OC(=O)/C=C\\C(=O)O"),
        ] {
            let e_once = stock_identity(e_form);
            let z_once = stock_identity(z_form);
            assert_eq!(e_once, stock_identity(&e_once), "{e_form} must be stable");
            assert_eq!(z_once, stock_identity(&z_once), "{z_form} must be stable");
            assert_ne!(
                e_once, z_once,
                "E and Z isomers must not collapse to the same canonical form"
            );
        }
    }
}
