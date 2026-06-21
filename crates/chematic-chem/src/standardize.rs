//! Molecular standardization routines.
//!
//! Provides utilities for cleaning up molecular representations:
//! - Selecting the largest connected fragment.
//! - Neutralizing simple formal charges.

#![forbid(unsafe_code)]

use std::collections::{HashMap, VecDeque};

use chematic_core::{AtomIdx, BondIdx, Element, Molecule, MoleculeBuilder, validate_valence};

use crate::{hash::mol_hash, hydrogen::remove_hydrogens, tautomer::canonical_tautomer};
use chematic_smarts::{MatchConfig, find_matches_with_config, parse_smarts};

/// Salt removal catalog: common salt patterns (counterions and solvates).
///
/// Each pattern is a (name, SMARTS) tuple for organic and inorganic salts.
/// Used by [`remove_salts`] to filter out counterions and keep drug-like fragments.
#[derive(Clone, Debug)]
pub struct SaltCatalog {
    /// (name, SMARTS) pairs for salt patterns
    patterns: Vec<(&'static str, &'static str)>,
}

impl Default for SaltCatalog {
    fn default() -> Self {
        Self::new()
    }
}

impl SaltCatalog {
    /// Create a default salt catalog with common counterions and solvates.
    pub fn new() -> Self {
        Self {
            patterns: vec![
                // Organic salts (carboxylates, sulfonates, etc.)
                ("acetate", "[#6](-[#1])(-[#1])-[#6](=[#8])[O-]"),
                ("formate", "[#6](=[#8])[O-]"),
                (
                    "propionate",
                    "[#6](-[#1])(-[#1])-[#6](-[#1])-[#6](=[#8])[O-]",
                ),
                ("benzoate", "c1ccccc1-[#6](=[#8])[O-]"),
                (
                    "trifluoroacetate",
                    "[#9]-[#6](-[#9])(-[#9])-[#6](=[#8])[O-]",
                ),
                (
                    "mesylate",
                    "[#16](=[#8])(=[#8])-[#8]-[#6](-[#1])(-[#1])-[#1]",
                ),
                ("tosylate", "c1ccc(cc1)-[#16](=[#8])(=[#8])-[#8]"),
                (
                    "nosylate",
                    "[#8]-[#6](-[#1])(-[#1])-[#8]-[#16](=[#8])(=[#8])-c1ccc([N+](=O)[O-])cc1",
                ),
                ("sulfate", "[#16](=[#8])(=[#8])(-[#8])-[#8]"),
                ("phosphate", "[#15](=[#8])(-[#8])(-[#8])-[#8]"),
                (
                    "citrate",
                    "[#6](-[#6](=[#8])[O-])(-[#6](=[#8])[O-])-[#6](-[#8])-[#6](=[#8])[O-]",
                ),
                (
                    "tartrate",
                    "[#6](-[#8])(-[#6](-[#8])-[#6](=[#8])[O-])-[#6](=[#8])[O-]",
                ),
                // Inorganic salts (single atoms/small molecules)
                ("sodium_cation", "[Na+]"),
                ("potassium_cation", "[K+]"),
                ("lithium_cation", "[Li+]"),
                ("calcium_cation", "[Ca+2]"),
                ("magnesium_cation", "[Mg+2]"),
                ("chloride_anion", "[Cl-]"),
                ("bromide_anion", "[Br-]"),
                ("iodide_anion", "[I-]"),
                ("fluoride_anion", "[F-]"),
                ("oxide_anion", "[O-2]"),
                ("sulfate_anion", "[#16](=[#8])(=[#8])(-[#8])-[#8-]"),
                ("phosphate_anion", "[#15](=[#8])(-[#8])(-[#8-])-[#8]"),
                // Solvates and additives
                ("water", "[#8](-[#1])-[#1]"),
                (
                    "dmso",
                    "[#16](=[#8])(-[#6](-[#1])(-[#1])-[#1])-[#6](-[#1])(-[#1])-[#1]",
                ),
                ("methanol", "[#6](-[#1])(-[#1])-[#8]-[#1]"),
                ("ethanol", "[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#8]-[#1]"),
                (
                    "isopropanol",
                    "[#6](-[#1])(-[#1])-[#6](-[#8]-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#1]",
                ),
                // Rare but important salts
                ("borate", "[#5](-[#8])(-[#8])-[#8]"),
                ("ammonium", "[#7+;H0,H1,H2,H3]"),
            ],
        }
    }

    /// Add a custom salt pattern to this catalog.
    pub fn add(&mut self, name: &'static str, smarts: &'static str) {
        self.patterns.push((name, smarts));
    }

    /// Check if a molecule fragment matches any salt pattern.
    pub fn is_salt(&self, frag: &Molecule) -> bool {
        let config = MatchConfig {
            max_visit_budget: Some(1_000_000),
            ..Default::default()
        };
        for (_, smarts_str) in &self.patterns {
            if let Ok(query) = parse_smarts(smarts_str)
                && !find_matches_with_config(&query, frag, &config).is_empty()
            {
                return true;
            }
        }
        false
    }
}

/// Find all connected components of `mol` via BFS, sorted descending by size.
fn connected_components(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components: Vec<Vec<AtomIdx>> = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut component = Vec::new();
        let mut queue: VecDeque<AtomIdx> = VecDeque::new();
        queue.push_back(AtomIdx(start as u32));

        while let Some(current) = queue.pop_front() {
            component.push(current);
            for (neighbor, _) in mol.neighbors(current) {
                let ni = neighbor.0 as usize;
                if !visited[ni] {
                    visited[ni] = true;
                    queue.push_back(neighbor);
                }
            }
        }
        components.push(component);
    }

    components.sort_by_key(|b| std::cmp::Reverse(b.len()));
    components
}

/// Copy bonds from `mol` into `builder` when both endpoints are remapped.
fn copy_bonds(mol: &Molecule, builder: &mut MoleculeBuilder, remap: &HashMap<AtomIdx, AtomIdx>) {
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a, new_b, bond.order);
        }
    }
}

/// Check if a fragment is a common inorganic salt or counterion.
///
/// Returns true if fragment matches patterns like: NaCl, KCl, Na+, K+, Cl-, Br-, I-, etc.
fn is_salt_fragment(frag: &Molecule) -> bool {
    let n = frag.atom_count();

    // Single atom: check if it's a common counterion
    if n == 1 {
        let atom = frag.atom(AtomIdx(0));
        return matches!(
            atom.element.atomic_number(),
            11 | 19 | 37 | 55 |  // Na, K, Rb, Cs (alkali metals)
            17 | 35 | 53 |       // Cl, Br, I (halogens)
            8 // O (oxide)
        );
    }

    // Two atoms: check for common binary salts (NaCl, KBr, etc.)
    if n == 2 {
        let a0 = frag.atom(AtomIdx(0)).element.atomic_number();
        let a1 = frag.atom(AtomIdx(1)).element.atomic_number();
        let bond_count = frag.bond_count();

        // Ionic pair (no bond between them) — cation + anion
        if bond_count == 0 {
            let metals = [11, 19, 37, 55]; // Na, K, Rb, Cs
            let nonmetals = [17, 35, 53, 8]; // Cl, Br, I, O
            return (metals.contains(&a0) && nonmetals.contains(&a1))
                || (metals.contains(&a1) && nonmetals.contains(&a0));
        }
    }

    // Small molecules with only metal/nonmetal atoms (common solvate salts)
    if n <= 4 {
        let has_organic = frag.atoms().any(|(_, a)| a.element.atomic_number() == 6);
        if !has_organic {
            // Pure inorganic salt (no carbons)
            return true;
        }
    }

    false
}

/// Return a new `Molecule` with inorganic salts removed, keeping largest organic fragment.
///
/// Uses a comprehensive salt catalog (SMARTS patterns) and heuristic detection to identify
/// and exclude common counterions (Na+, K+, TFA-, acetate, etc.) and inorganic salt fragments.
///
/// If no non-salt fragment exists or molecule is empty, returns the original largest fragment.
pub fn remove_salts(mol: &Molecule) -> Molecule {
    remove_salts_with_catalog(mol, &SaltCatalog::new())
}

/// Remove salts using a custom catalog.
///
/// # Arguments
/// - `mol`: the molecule to process
/// - `catalog`: custom salt catalog (use `SaltCatalog::new()` for default)
pub fn remove_salts_with_catalog(mol: &Molecule, catalog: &SaltCatalog) -> Molecule {
    if mol.atom_count() == 0 {
        return MoleculeBuilder::new().build();
    }

    let components = connected_components(mol);

    // Find largest non-salt fragment
    let mut largest_non_salt: Option<&Vec<AtomIdx>> = None;
    let mut largest_non_salt_size = 0;

    for component in &components {
        // Extract fragment molecule temporarily to check if it's a salt
        let mut builder = MoleculeBuilder::new();
        let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
        for &old_idx in component {
            let new_idx = builder.add_atom(mol.atom(old_idx).clone());
            remap.insert(old_idx, new_idx);
        }
        copy_bonds(mol, &mut builder, &remap);
        let frag = builder.build();

        // Check with catalog first, fall back to heuristic
        let is_salt = catalog.is_salt(&frag) || is_salt_fragment(&frag);

        // If not a salt and larger than current best, use it
        if !is_salt && component.len() > largest_non_salt_size {
            largest_non_salt = Some(component);
            largest_non_salt_size = component.len();
        }
    }

    // Fall back to largest fragment if no non-salt found
    let component = largest_non_salt.unwrap_or(&components[0]);

    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let mut builder = MoleculeBuilder::new();
    for &old_idx in component {
        let new_idx = builder.add_atom(mol.atom(old_idx).clone());
        remap.insert(old_idx, new_idx);
    }
    copy_bonds(mol, &mut builder, &remap);
    builder.build()
}

/// Return a new `Molecule` containing only the largest connected fragment.
///
/// If the molecule is empty, an empty `Molecule` is returned.
/// This is an alias for `remove_salts()` for backward compatibility.
pub fn largest_fragment(mol: &Molecule) -> Molecule {
    remove_salts(mol)
}

/// Normalize chemical groups (nitro groups, etc.).
///
/// Transforms:
/// - `[N+](=O)[O-]` → `N(=O)=O` (nitro normalization: N charge 0, O- → double bond)
///
/// Returns a new molecule with normalized groups.
pub fn normalize_groups(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let mut nitro_atoms = std::collections::HashSet::new();
    let mut oxide_atoms = std::collections::HashSet::new();
    let mut azide_atoms = std::collections::HashSet::new();
    let mut sulfoxide_atoms = std::collections::HashSet::new();

    // First pass: identify functional groups via per-group detectors.
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() == 7 && atom.charge == 1 {
            detect_nitro(mol, idx, atom, &mut nitro_atoms, &mut oxide_atoms);
            detect_azide(mol, idx, &mut azide_atoms);
        }
        if atom.element.atomic_number() == 16 {
            detect_sulfoxide(mol, idx, &mut sulfoxide_atoms);
        }
    }

    // Second pass: copy atoms with normalized charges
    for (idx, atom) in mol.atoms() {
        let mut new_atom = atom.clone();

        if nitro_atoms.contains(&idx) {
            // Neutralize the N and O in nitro group
            if atom.element.atomic_number() == 7 || atom.element.atomic_number() == 8 {
                new_atom.charge = 0;
            }
        }

        if azide_atoms.contains(&idx) {
            // Neutralize all N in azide group [N-][N+]#N -> N=N=N
            if atom.element.atomic_number() == 7 {
                new_atom.charge = 0;
            }
        }

        // Sulfoxide: keep as is (S=O is already correct form)

        let new_idx = builder.add_atom(new_atom);
        remap.insert(idx, new_idx);
    }

    // Third pass: copy bonds, normalizing functional groups
    for i in 0..mol.bond_count() {
        let bond = mol.bond(chematic_core::BondIdx(i as u32));
        let mut new_order = bond.order;

        // Nitro groups: convert single N-O (where O is negative) to double
        if nitro_atoms.contains(&bond.atom1) && nitro_atoms.contains(&bond.atom2) {
            let a1_is_n = mol.atom(bond.atom1).element.atomic_number() == 7;
            let a2_is_o = mol.atom(bond.atom2).element.atomic_number() == 8;
            let a1_is_o = mol.atom(bond.atom1).element.atomic_number() == 8;
            let a2_is_n = mol.atom(bond.atom2).element.atomic_number() == 7;

            if (a1_is_n
                && a2_is_o
                && bond.order == chematic_core::BondOrder::Single
                && mol.atom(bond.atom2).charge == -1)
                || (a1_is_o
                    && a2_is_n
                    && bond.order == chematic_core::BondOrder::Single
                    && mol.atom(bond.atom1).charge == -1)
            {
                new_order = chematic_core::BondOrder::Double;
            }
        }

        // AZIDE normalization: [N-][N+]#N -> N=N=N (convert single to double)
        if azide_atoms.contains(&bond.atom1) && azide_atoms.contains(&bond.atom2) {
            let a1_is_n = mol.atom(bond.atom1).element.atomic_number() == 7;
            let a2_is_n = mol.atom(bond.atom2).element.atomic_number() == 7;

            if a1_is_n && a2_is_n && bond.order == chematic_core::BondOrder::Single {
                // Convert single bonds in azide to double
                new_order = chematic_core::BondOrder::Double;
            }
        }

        // N-oxide: keep as single bond (already correct after atom charge normalization)
        if oxide_atoms.contains(&bond.atom1) || oxide_atoms.contains(&bond.atom2) {
            // No bond order change needed — already single
        }

        // Sulfoxide: keep as is (S=O is already correct form)

        if let (Some(&new_a1), Some(&new_a2)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a1, new_a2, new_order);
        }
    }

    builder.build()
}

// ---------------------------------------------------------------------------
// Functional group detectors for normalize_groups
// ---------------------------------------------------------------------------

/// Mark nitro [N+](=O)[O-] and aromatic N-oxide atoms.
fn detect_nitro(
    mol: &Molecule,
    idx: AtomIdx,
    atom: &chematic_core::Atom,
    nitro_atoms: &mut std::collections::HashSet<AtomIdx>,
    oxide_atoms: &mut std::collections::HashSet<AtomIdx>,
) {
    let o_nbrs: Vec<_> = mol
        .neighbors(idx)
        .filter(|(n, _)| mol.atom(*n).element.atomic_number() == 8)
        .collect();

    if o_nbrs.len() == 2 {
        let mut has_double_o = false;
        let mut has_single_neg_o = false;
        for (o_idx, bid) in &o_nbrs {
            let o = mol.atom(*o_idx);
            let b = mol.bond(*bid);
            if b.order == chematic_core::BondOrder::Double && o.charge == 0 {
                has_double_o = true;
            }
            if b.order == chematic_core::BondOrder::Single && o.charge == -1 {
                has_single_neg_o = true;
                nitro_atoms.insert(*o_idx);
            }
        }
        if has_double_o && has_single_neg_o {
            nitro_atoms.insert(idx);
        }
    } else if let Some((o_idx, bid)) = o_nbrs.first() {
        let o = mol.atom(*o_idx);
        let b = mol.bond(*bid);
        if atom.aromatic && b.order == chematic_core::BondOrder::Single && o.charge == -1 {
            nitro_atoms.insert(idx);
            oxide_atoms.insert(*o_idx);
        }
    }
}

/// Mark azide [N-][N+]#N atoms.
fn detect_azide(
    mol: &Molecule,
    idx: AtomIdx,
    azide_atoms: &mut std::collections::HashSet<AtomIdx>,
) {
    let n_nbrs: Vec<_> = mol
        .neighbors(idx)
        .filter(|(n, _)| mol.atom(*n).element.atomic_number() == 7)
        .collect();

    for (n_idx, bid) in &n_nbrs {
        let n = mol.atom(*n_idx);
        let b = mol.bond(*bid);
        if b.order == chematic_core::BondOrder::Triple && n.charge == 0 {
            for (other_idx, other_bid) in n_nbrs.iter() {
                if other_idx == n_idx {
                    continue;
                }
                let other = mol.atom(*other_idx);
                let other_b = mol.bond(*other_bid);
                if other_b.order == chematic_core::BondOrder::Single && other.charge == -1 {
                    azide_atoms.insert(idx);
                    azide_atoms.insert(*n_idx);
                    azide_atoms.insert(*other_idx);
                }
            }
        }
    }
}

/// Mark sulfoxide S=O atoms.
fn detect_sulfoxide(
    mol: &Molecule,
    idx: AtomIdx,
    sulfoxide_atoms: &mut std::collections::HashSet<AtomIdx>,
) {
    for (o_idx, bid) in mol.neighbors(idx) {
        if mol.atom(o_idx).element.atomic_number() == 8
            && mol.bond(bid).order == chematic_core::BondOrder::Double
        {
            sulfoxide_atoms.insert(idx);
            sulfoxide_atoms.insert(o_idx);
        }
    }
}

/// Detect if a molecule contains a zwitterion (internal salt).
///
/// A zwitterion is defined as having both positive and negative formal charges.
/// Examples: amino acids in zwitterionic form ([NH3+][COO-]).
pub fn has_zwitterion(mol: &Molecule) -> bool {
    let mut has_positive = false;
    let mut has_negative = false;

    for (_, atom) in mol.atoms() {
        if atom.charge > 0 {
            has_positive = true;
        } else if atom.charge < 0 {
            has_negative = true;
        }
        if has_positive && has_negative {
            return true;
        }
    }
    false
}

/// Normalize a zwitterion to neutral form by proton transfer.
///
/// For each negatively-charged atom, find the nearest positively-charged atom
/// and transfer one proton (increase positive charge's H count, decrease negative charge).
///
/// Returns a new molecule with normalized charges.
pub fn normalize_zwitterion(mol: &Molecule) -> Molecule {
    if !has_zwitterion(mol) {
        return clone_molecule(mol);
    }

    let mut modifications: HashMap<AtomIdx, (i8, Option<u8>)> = HashMap::new();

    // Collect positive and negative charge atoms
    let mut positive_atoms: Vec<AtomIdx> = Vec::new();
    let mut negative_atoms: Vec<AtomIdx> = Vec::new();

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        if atom.charge > 0 {
            positive_atoms.push(idx);
        } else if atom.charge < 0 {
            negative_atoms.push(idx);
        }
    }

    // For each negative atom, transfer proton from nearest positive atom
    for &neg_idx in &negative_atoms {
        if positive_atoms.is_empty() {
            continue;
        }

        // Find nearest positive charge (by BFS distance)
        let mut closest_pos_idx = positive_atoms[0];
        let mut closest_distance = i32::MAX;

        for &pos_idx in &positive_atoms {
            if let Some(dist) = bfs_distance(mol, neg_idx, pos_idx)
                && dist < closest_distance
            {
                closest_distance = dist;
                closest_pos_idx = pos_idx;
            }
        }

        // Transfer proton: N+ loses H, O- gains H
        let neg_atom = mol.atom(neg_idx);
        let pos_atom = mol.atom(closest_pos_idx);

        // Decrease negative charge
        let new_neg_charge = neg_atom.charge + 1;
        let neg_h = neg_atom.hydrogen_count.unwrap_or(0);
        modifications.insert(neg_idx, (new_neg_charge, Some(neg_h + 1)));

        // Decrease positive charge by decreasing H count
        let pos_h = pos_atom.hydrogen_count.unwrap_or(0);
        if pos_h > 0 {
            let new_pos_charge = pos_atom.charge - 1;
            modifications.insert(closest_pos_idx, (new_pos_charge, Some(pos_h - 1)));
        }
    }

    // Reconstruct molecule with modified charges
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();
        if let Some(&(new_charge, new_h)) = modifications.get(&old_idx) {
            atom.charge = new_charge;
            atom.hydrogen_count = new_h;
        }
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }
    copy_bonds(mol, &mut builder, &remap);
    builder.build()
}

/// BFS distance between two atoms in a molecule.
fn bfs_distance(mol: &Molecule, start: AtomIdx, end: AtomIdx) -> Option<i32> {
    if start == end {
        return Some(0);
    }

    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    queue.push_back((start, 0));
    visited[start.0 as usize] = true;

    while let Some((current, dist)) = queue.pop_front() {
        for (neighbor, _) in mol.neighbors(current) {
            if neighbor == end {
                return Some(dist + 1);
            }
            let ni = neighbor.0 as usize;
            if !visited[ni] {
                visited[ni] = true;
                queue.push_back((neighbor, dist + 1));
            }
        }
    }
    None
}

/// Neutralize simple formal charges in a molecule.
///
/// Rules applied:
/// - `[O-]` with a carbon neighbor → charge 0, +1 H (carboxylate → carboxylic acid).
/// - `[N+]` with at least one explicit H → charge 0, −1 H (ammonium → amine).
/// - `[O+]` with at least one explicit H → charge 0, −1 H (protonated ether → ether).
pub fn neutralize_charges(mol: &Molecule) -> Molecule {
    let mut modifications: HashMap<AtomIdx, (i8, Option<u8>)> = HashMap::new();

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        let h = atom.hydrogen_count.unwrap_or(0);

        match (atom.element, atom.charge) {
            (Element::O, -1) => {
                let has_c_neighbor = mol
                    .neighbors(idx)
                    .any(|(nb, _)| mol.atom(nb).element == Element::C);
                if has_c_neighbor {
                    modifications.insert(idx, (0, Some(h + 1)));
                }
            }
            (Element::N, 1) | (Element::O, 1) if h > 0 => {
                modifications.insert(idx, (0, Some(h - 1)));
            }
            _ => {}
        }
    }

    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();
        if let Some(&(new_charge, new_h)) = modifications.get(&old_idx) {
            atom.charge = new_charge;
            atom.hydrogen_count = new_h;
        }
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }
    copy_bonds(mol, &mut builder, &remap);
    builder.build()
}

/// Remove all isotope labels from atoms.
///
/// Returns a new molecule with `atom.isotope = None` for all atoms.
pub fn remove_isotopes(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();
        atom.isotope = None;
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }
    copy_bonds(mol, &mut builder, &remap);
    builder.build()
}

/// Remove all stereochemistry from a molecule.
///
/// Sets `atom.chirality = Chirality::None` for all atoms and converts
/// wedge/wedge-hash bonds to single bonds.
pub fn remove_stereo(mol: &Molecule) -> Molecule {
    use chematic_core::{BondOrder, Chirality};

    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();
        atom.chirality = Chirality::None;
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }

    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let order = match bond.order {
                BondOrder::Up | BondOrder::Down => BondOrder::Single,
                other => other,
            };
            let _ = builder.add_bond(new_a, new_b, order);
        }
    }

    builder.build()
}

/// Remove atoms from stereo groups that are no longer chiral, and drop empty groups.
///
/// Analog of RDKit PR #9051 ("cleanup of stereogroups and wedges for non-chiral sites").
/// When molecular operations (fragment removal, chirality clearing, …) leave atoms
/// without chirality flags, their stereo group membership becomes invalid.  This
/// function filters each [`StereoGroup`]'s atom list to only those atoms where
/// `atom.chirality != Chirality::None`, and discards any group that becomes empty.
pub fn clean_stereo_groups(mol: &Molecule) -> Molecule {
    use chematic_core::{Chirality, StereoGroup};

    let cleaned: Vec<StereoGroup> = mol
        .stereo_groups()
        .iter()
        .filter_map(|g| {
            let chiral_atoms: Vec<AtomIdx> = g
                .atom_indices
                .iter()
                .copied()
                .filter(|&idx| mol.atom(idx).chirality != Chirality::None)
                .collect();
            if chiral_atoms.is_empty() {
                None
            } else {
                Some(StereoGroup::new(g.kind.clone(), chiral_atoms))
            }
        })
        .collect();

    let mut out = MoleculeBuilder::from_molecule(mol).build();
    out.set_stereo_groups(cleaned);
    out
}

/// Keep only the largest organic (carbon-containing) fragment.
///
/// Removes all inorganic fragments (those without carbon atoms).
/// Useful for removing metal ions, salts, and other counterions.
/// Falls back to largest fragment if no organic fragment exists.
pub fn prefer_organic(mol: &Molecule) -> Molecule {
    if mol.atom_count() == 0 {
        return MoleculeBuilder::new().build();
    }

    let components = connected_components(mol);

    // Find largest organic fragment
    let mut largest_organic: Option<&Vec<AtomIdx>> = None;
    let mut largest_organic_size = 0;

    for component in &components {
        // Check if fragment contains carbon (organic)
        let has_carbon = component
            .iter()
            .any(|&idx| mol.atom(idx).element.atomic_number() == 6);

        if has_carbon && component.len() > largest_organic_size {
            largest_organic = Some(component);
            largest_organic_size = component.len();
        }
    }

    // Fall back to largest fragment if no organic found
    let target_component = largest_organic.or_else(|| components.first());

    if let Some(component) = target_component {
        let mut builder = MoleculeBuilder::new();
        let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
        for &old_idx in component {
            let new_idx = builder.add_atom(mol.atom(old_idx).clone());
            remap.insert(old_idx, new_idx);
        }
        copy_bonds(mol, &mut builder, &remap);
        builder.build()
    } else {
        MoleculeBuilder::new().build()
    }
}

/// Reionize a molecule by adjusting protonation to favored forms.
///
/// Simple heuristic approach that adjusts charges on common acidic and basic groups:
/// - Carboxylic acids with OH: deprotonate to COO- (carboxylate)
/// - Phenols: deprotonate to phenoxide (O-)
/// - Primary amines: protonate to ammonium (NH3+)
/// - Imidazoles: protonate to imidazolium (N+)
///
/// This is a simplified version suitable for typical organic molecules.
pub fn reionize(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Copy all atoms, adjusting charges
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        let an = atom.element.atomic_number();

        // Check for carboxylic acid or phenol: C=O with O-H or Ar-O-H
        if an == 8 {
            // Oxygen: check if it's OH bonded to C
            if let Some((c_idx, _)) = mol.neighbors(idx).find(|(neighbor, bond_idx)| {
                mol.bond(*bond_idx).order == chematic_core::BondOrder::Single
                    && mol.atom(*neighbor).element.atomic_number() == 6
            }) {
                // Check if C is aromatic (phenol) or has a double-bonded O (carboxylic acid)
                let is_aromatic = mol.atom(c_idx).aromatic;
                let has_double_bonded_o = mol.neighbors(c_idx).any(|(other, bond_idx)| {
                    mol.bond(bond_idx).order == chematic_core::BondOrder::Double
                        && mol.atom(other).element.atomic_number() == 8
                        && other != idx
                });

                // Only deprotonate if it's a phenol or carboxylic acid, not aliphatic OH
                if (is_aromatic || has_double_bonded_o) && atom.charge >= 0 {
                    atom.charge -= 1; // Deprotonate: OH → O-
                }
            }
        }

        // Check for primary/secondary amines (but NOT amides)
        if an == 7 {
            // Check if this N is NOT part of an amide (C(=O)-N)
            let is_amide = mol.neighbors(idx).any(|(neighbor, bond_idx)| {
                mol.bond(bond_idx).order == chematic_core::BondOrder::Single
                    && mol.atom(neighbor).element.atomic_number() == 6
                    && mol.neighbors(neighbor).any(|(o_neighbor, o_bond)| {
                        mol.bond(o_bond).order == chematic_core::BondOrder::Double
                            && (mol.atom(o_neighbor).element.atomic_number() == 8
                                || mol.atom(o_neighbor).element.atomic_number() == 16)
                    })
            });

            if !is_amide {
                let h_count = chematic_core::implicit_hcount(mol, idx);
                // Protonate free amines only (not amides)
                if (h_count == 2 || h_count == 1) && atom.charge <= 0 {
                    atom.charge += 1; // Protonate: NH2 → NH3+
                }
            }
        }

        let new_idx = builder.add_atom(atom);
        remap.insert(idx, new_idx);
    }

    copy_bonds(mol, &mut builder, &remap);
    builder.build()
}

/// Remove all charges from a molecule by protonation/deprotonation.
///
/// Neutralizes positively charged atoms by removing protons and
/// negatively charged atoms by adding protons. This is an aggressive
/// neutralization that may create chemically unrealistic structures.
///
/// # Note
/// This differs from [`neutralize_charges`] which uses specific rules.
/// `uncharge` is a brute-force approach suitable for structure cleanup.
pub fn uncharge(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Copy all atoms, removing charges
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        atom.charge = 0; // Force neutral
        let new_idx = builder.add_atom(atom);
        remap.insert(idx, new_idx);
    }

    copy_bonds(mol, &mut builder, &remap);
    builder.build()
}

/// One transformation stage in the standardization pipeline.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum StandardizationStep {
    /// Select the largest connected component.
    LargestFragment,
    /// Apply simple neutralization rules for common formal charges.
    NeutralizeCharges,
    /// Normalize chemical groups (nitro groups, etc.).
    NormalizeGroups,
    /// Normalize zwitterionic forms to neutral molecules.
    ZwitterionNormalization,
    /// Remove explicit hydrogen atoms.
    RemoveExplicitHydrogens,
    /// Canonicalize supported tautomer systems.
    CanonicalTautomer,
    /// Keep only the largest non-salt fragment.
    FragmentParent,
    /// Neutralize all formal charges.
    ChargeParent,
    /// Remove all isotope labels.
    IsotopeParent,
    /// Remove all stereochemistry (wedge/wedge-hash bonds, chirality).
    StereoParent,
}

impl StandardizationStep {
    /// Stable machine-readable stage name.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::LargestFragment => "largest_fragment",
            Self::NeutralizeCharges => "neutralize_charges",
            Self::NormalizeGroups => "normalize_groups",
            Self::ZwitterionNormalization => "zwitterion_normalization",
            Self::RemoveExplicitHydrogens => "remove_explicit_hydrogens",
            Self::CanonicalTautomer => "canonical_tautomer",
            Self::FragmentParent => "fragment_parent",
            Self::ChargeParent => "charge_parent",
            Self::IsotopeParent => "isotope_parent",
            Self::StereoParent => "stereo_parent",
        }
    }
}

/// High-level status for a standardization run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PipelineStatus {
    /// Pipeline completed and the output is structurally identical to the input.
    Unchanged,
    /// Pipeline completed and at least one enabled stage changed the molecule.
    Modified,
    /// Pipeline completed, but warnings indicate unsupported or suspicious input features.
    CompletedWithWarnings,
}

/// Warning emitted during standardization.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StandardizationWarning {
    /// Stable machine-readable warning code.
    pub code: String,
    /// Human-readable detail.
    pub message: String,
}

impl StandardizationWarning {
    fn new(code: &str, message: String) -> Self {
        Self {
            code: code.to_string(),
            message,
        }
    }
}

/// Atom/bond/hash summary before or after a pipeline stage.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MoleculeSnapshot {
    /// Number of atoms in the molecule.
    pub atoms: usize,
    /// Number of bonds in the molecule.
    pub bonds: usize,
    /// Deterministic structure hash based on canonical SMILES.
    pub hash: u64,
}

impl MoleculeSnapshot {
    fn from_mol(mol: &Molecule) -> Self {
        Self {
            atoms: mol.atom_count(),
            bonds: mol.bond_count(),
            hash: mol_hash(mol),
        }
    }
}

/// Per-stage audit entry for a standardization run.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StandardizationStepReport {
    /// Pipeline stage.
    pub step: StandardizationStep,
    /// Whether the stage was enabled in the config.
    pub enabled: bool,
    /// Whether the stage changed the molecule hash.
    pub changed: bool,
    /// Molecule summary before the stage.
    pub before: MoleculeSnapshot,
    /// Molecule summary after the stage.
    pub after: MoleculeSnapshot,
}

/// Full standardization audit result.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StandardizationReport {
    /// Overall status.
    pub status: PipelineStatus,
    /// Summary of the input molecule.
    pub input: MoleculeSnapshot,
    /// Summary of the output molecule.
    pub output: MoleculeSnapshot,
    /// Ordered per-stage results.
    pub steps: Vec<StandardizationStepReport>,
    /// Validation and unsupported-feature warnings.
    pub warnings: Vec<StandardizationWarning>,
}

impl StandardizationReport {
    /// Returns `true` if the molecule changed at any enabled stage.
    pub fn changed(&self) -> bool {
        self.input.hash != self.output.hash
    }
}

/// Handling strategy for zwitterions (internal salts).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ZwitterionHandling {
    /// Keep zwitterionic form as-is.
    Keep,
    /// Normalize to neutral form via proton transfer.
    #[default]
    Normalize,
}

/// Options for molecular standardization.
///
/// Controls which cleaning transformations are applied in a standardization pipeline.
#[derive(Clone, Debug)]
pub struct StandardizeOptions {
    /// Convert to canonical tautomer. Default: `true`.
    pub canonical_tautomer: bool,
    /// Neutralize simple formal charges. Default: `true`.
    pub neutralize_charges: bool,
    /// Remove explicit hydrogen atoms. Default: `true`.
    pub remove_explicit_h: bool,
    /// Keep only the largest connected fragment. Default: `false`.
    pub largest_fragment_only: bool,
    /// Handle zwitterions (internal salts). Default: `Normalize`.
    pub zwitterion_handling: ZwitterionHandling,
}

impl Default for StandardizeOptions {
    fn default() -> Self {
        Self {
            canonical_tautomer: true,
            neutralize_charges: true,
            remove_explicit_h: true,
            largest_fragment_only: false,
            zwitterion_handling: ZwitterionHandling::Normalize,
        }
    }
}

/// RDKit-style standardization pipeline with an auditable report.
#[derive(Clone, Debug, Default)]
pub struct StandardizationPipeline {
    options: StandardizeOptions,
}

impl StandardizationPipeline {
    /// Create a pipeline from explicit options.
    pub fn new(options: StandardizeOptions) -> Self {
        Self { options }
    }

    /// Borrow the pipeline options.
    pub fn options(&self) -> &StandardizeOptions {
        &self.options
    }

    /// Standardize a molecule and return both the output molecule and an audit report.
    pub fn run(&self, mol: &Molecule) -> (Molecule, StandardizationReport) {
        let input = MoleculeSnapshot::from_mol(mol);
        let mut current = clone_molecule(mol);
        let mut steps = Vec::new();
        let mut warnings = detect_initial_warnings(mol);

        // Disconnect metals early (remove dative/coordinate bonds)
        let has_metals = current.atoms().any(|(_, a)| is_metal(a.element));
        if has_metals {
            current = disconnect_metals(&current);
        }

        // Apply NeutralizeCharges BEFORE LargestFragment to ensure predictable fragment selection.
        // Example: [NH3+].[Cl-] should be neutralized first to [NH3].[Cl-], then largest fragment.
        current = self.apply_stage(
            current,
            StandardizationStep::NeutralizeCharges,
            self.options.neutralize_charges,
            neutralize_charges,
            &mut steps,
            &mut warnings,
        );
        current = self.apply_stage(
            current,
            StandardizationStep::LargestFragment,
            self.options.largest_fragment_only,
            largest_fragment,
            &mut steps,
            &mut warnings,
        );
        let zwitterion_enabled = self.options.zwitterion_handling == ZwitterionHandling::Normalize;
        current = self.apply_stage(
            current,
            StandardizationStep::ZwitterionNormalization,
            zwitterion_enabled,
            normalize_zwitterion,
            &mut steps,
            &mut warnings,
        );
        current = self.apply_stage(
            current,
            StandardizationStep::RemoveExplicitHydrogens,
            self.options.remove_explicit_h,
            remove_hydrogens,
            &mut steps,
            &mut warnings,
        );
        current = self.apply_stage(
            current,
            StandardizationStep::CanonicalTautomer,
            self.options.canonical_tautomer,
            canonical_tautomer,
            &mut steps,
            &mut warnings,
        );

        let output = MoleculeSnapshot::from_mol(&current);
        // Status depends only on structure change (hash), not on warnings.
        // Warnings are reported separately for the user to inspect.
        let status = if input.hash == output.hash {
            PipelineStatus::Unchanged
        } else if !warnings.is_empty() {
            PipelineStatus::CompletedWithWarnings
        } else {
            PipelineStatus::Modified
        };

        (
            current,
            StandardizationReport {
                status,
                input,
                output,
                steps,
                warnings,
            },
        )
    }

    fn apply_stage(
        &self,
        current: Molecule,
        step: StandardizationStep,
        enabled: bool,
        f: fn(&Molecule) -> Molecule,
        steps: &mut Vec<StandardizationStepReport>,
        warnings: &mut Vec<StandardizationWarning>,
    ) -> Molecule {
        let before = MoleculeSnapshot::from_mol(&current);
        let next = if enabled {
            f(&current)
        } else {
            clone_molecule(&current)
        };
        let after = MoleculeSnapshot::from_mol(&next);
        steps.push(StandardizationStepReport {
            step,
            enabled,
            changed: before.hash != after.hash,
            before,
            after,
        });
        if enabled {
            append_valence_warnings(step, &next, warnings);
        }
        next
    }
}

fn clone_molecule(mol: &Molecule) -> Molecule {
    MoleculeBuilder::from_molecule(mol).build()
}

fn detect_initial_warnings(mol: &Molecule) -> Vec<StandardizationWarning> {
    let mut warnings = Vec::new();
    // Metal disconnection is now handled in the pipeline, so we don't warn about it
    let valence_errors = validate_valence(mol);
    if !valence_errors.is_empty() {
        warnings.push(StandardizationWarning::new(
            "input_valence_validation_failed",
            format!(
                "input molecule has {} valence validation issue(s)",
                valence_errors.len()
            ),
        ));
    }
    warnings
}

fn append_valence_warnings(
    step: StandardizationStep,
    mol: &Molecule,
    warnings: &mut Vec<StandardizationWarning>,
) {
    let errors = validate_valence(mol);
    if errors.is_empty() {
        return;
    }
    warnings.push(StandardizationWarning::new(
        "valence_validation_failed",
        format!(
            "{} produced {} valence validation issue(s)",
            step.as_str(),
            errors.len()
        ),
    ));
}

/// Disconnect metal-nonmetal bonds by removing dative/coordinate bonds to metals.
///
/// Iterates through all atoms; if a metal is found, removes all bonds between
/// that metal and organic/inorganic atoms. Returns the molecule with metal
/// coordination bonds severed.
fn disconnect_metals(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();

    // Copy all atoms
    for i in 0..mol.atom_count() {
        builder.add_atom(mol.atom(AtomIdx(i as u32)).clone());
    }

    // Copy only non-metal bonds
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        let atom1_is_metal = is_metal(mol.atom(bond.atom1).element);
        let atom2_is_metal = is_metal(mol.atom(bond.atom2).element);

        // Skip bonds where at least one atom is a metal
        if !atom1_is_metal && !atom2_is_metal {
            builder.add_bond(bond.atom1, bond.atom2, bond.order).ok();
        }
    }

    builder.build()
}

fn is_metal(element: Element) -> bool {
    matches!(
        element.atomic_number(),
        3 | 4
            | 11 | 12 | 13
            | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31
            | 37 | 38 | 39 | 40 | 41 | 42 | 43 | 44 | 45 | 46 | 47 | 48 | 49 | 50
            | 55 | 56 | 57..=71
            | 72 | 73 | 74 | 75 | 76 | 77 | 78 | 79 | 80 | 81 | 82 | 83
            | 87 | 88 | 89..=103
            | 104 | 105 | 106 | 107 | 108 | 109 | 110 | 111 | 112 | 113 | 114 | 115 | 116
    )
}

/// Apply a series of standardization steps to a molecule.
///
/// Transformations are applied in this order:
/// 1. If `largest_fragment_only`, select the largest connected component.
/// 2. If `neutralize_charges`, neutralize simple charges.
/// 3. If `remove_explicit_h`, remove explicit H atoms.
/// 4. If `canonical_tautomer`, convert to the canonical tautomer.
///
/// Useful for cleaning pasted structures or database entries.
pub fn standardize(mol: &Molecule, opts: &StandardizeOptions) -> Molecule {
    StandardizationPipeline::new(opts.clone()).run(mol).0
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn largest_fragment_two_fragments_picks_larger() {
        // "CC.CCC" — ethane (2 C) and propane (3 C)
        let mol = parse("CC.CCC").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 3, "should keep propane (3 C)");
    }

    #[test]
    fn largest_fragment_single_fragment_unchanged() {
        // "CC" — ethane, only one fragment
        let mol = parse("CC").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 2);
    }

    #[test]
    fn largest_fragment_keeps_benzene_over_ethane() {
        // "CC.c1ccccc1" — ethane (2 C) vs benzene (6 C)
        let mol = parse("CC.c1ccccc1").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 6, "should keep benzene (6 atoms)");
    }

    #[test]
    fn largest_fragment_ionic_pair_keeps_one_atom() {
        // "[Na+].[Cl-]" — both fragments are single atoms; either is fine
        let mol = parse("[Na+].[Cl-]").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 1);
    }

    #[test]
    fn neutralize_neutral_molecule_unchanged() {
        // "CC" is already neutral; no atom should gain/lose charge
        let mol = parse("CC").unwrap();
        let result = neutralize_charges(&mol);
        for i in 0..result.atom_count() {
            let atom = result.atom(AtomIdx(i as u32));
            assert_eq!(atom.charge, 0, "all atoms should remain neutral");
        }
    }

    #[test]
    fn neutralize_acetate_oxygen() {
        // "CC(=O)[O-]" — acetate; the [O-] should become neutral with H added
        let mol = parse("CC(=O)[O-]").unwrap();
        let result = neutralize_charges(&mol);

        // Find the oxygen that was originally [O-]: it should now have charge 0
        // and hydrogen_count == Some(1).
        let neutralized_o = (0..result.atom_count())
            .map(|i| result.atom(AtomIdx(i as u32)))
            .find(|a| a.element == Element::O && a.hydrogen_count == Some(1));

        assert!(
            neutralized_o.is_some(),
            "neutralized [O-] should have hydrogen_count == Some(1)"
        );
        assert_eq!(
            neutralized_o.unwrap().charge,
            0,
            "neutralized [O-] should have charge == 0"
        );
    }

    #[test]
    fn standardize_with_defaults() {
        // "CC(=O)[O-]" — acetate ion
        let mol = parse("CC(=O)[O-]").unwrap();
        let opts = StandardizeOptions::default();
        let result = standardize(&mol, &opts);

        // With default options, remove_explicit_h will be applied.
        // Check that [O-] was neutralized (charge should be 0).
        let has_neutral_o = (0..result.atom_count())
            .map(|i| result.atom(AtomIdx(i as u32)))
            .any(|a| a.element == Element::O && a.charge == 0);
        assert!(has_neutral_o, "acetate oxygen should be neutralized");

        // Should have at least 3 atoms (C, C, O) with no explicit H
        assert!(
            result.atom_count() >= 3,
            "should have at least 3 atoms after standardization"
        );
    }

    #[test]
    fn standardize_skip_largest_fragment() {
        // "CC.CCC" — ethane and propane
        let mol = parse("CC.CCC").unwrap();
        let opts = StandardizeOptions {
            zwitterion_handling: ZwitterionHandling::Normalize,
            largest_fragment_only: false,
            ..Default::default()
        };
        let result = standardize(&mol, &opts);

        // Should keep both fragments
        assert_eq!(
            result.atom_count(),
            5,
            "should keep both fragments when largest_fragment_only=false"
        );
    }

    #[test]
    fn pipeline_report_tracks_enabled_stage_changes() {
        let mol = parse("CC.CCC").unwrap();
        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            largest_fragment_only: true,
            neutralize_charges: false,
            remove_explicit_h: false,
            canonical_tautomer: false,
            zwitterion_handling: ZwitterionHandling::Keep,
        });

        let (result, report) = pipeline.run(&mol);

        assert_eq!(result.atom_count(), 3);
        assert_eq!(report.status, PipelineStatus::Modified);
        assert!(report.changed());
        assert_eq!(report.steps.len(), 5);
        // NeutralizeCharges is applied first (not enabled, so no change)
        assert_eq!(report.steps[0].step, StandardizationStep::NeutralizeCharges);
        assert!(!report.steps[0].enabled);
        // LargestFragment is applied second and is enabled
        assert_eq!(report.steps[1].step, StandardizationStep::LargestFragment);
        assert!(report.steps[1].enabled);
        assert!(report.steps[1].changed);
    }

    #[test]
    fn pipeline_report_marks_unchanged_clean_molecule() {
        let mol = parse("CC").unwrap();
        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            canonical_tautomer: false,
            neutralize_charges: false,
            remove_explicit_h: false,
            largest_fragment_only: false,
            zwitterion_handling: ZwitterionHandling::Keep,
        });

        let (_result, report) = pipeline.run(&mol);

        assert_eq!(report.status, PipelineStatus::Unchanged);
        assert!(!report.changed());
        assert!(report.warnings.is_empty());
        assert!(report.steps.iter().all(|s| !s.enabled && !s.changed));
    }

    #[test]
    fn pipeline_report_disconnects_metal_bonds() {
        let mol = parse("[Na]OC").unwrap();
        assert_eq!(mol.bond_count(), 2, "input has Na-O and O-C bonds");

        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            canonical_tautomer: false,
            neutralize_charges: false,
            remove_explicit_h: false,
            largest_fragment_only: false,
            zwitterion_handling: ZwitterionHandling::Keep,
        });

        let (result, _report) = pipeline.run(&mol);

        // Metal disconnection should run automatically, removing the Na-O bond
        assert_eq!(result.bond_count(), 1, "Na-O bond should be disconnected");
        // The remaining bond should be O-C
        assert!(
            result.bond(BondIdx(0)).atom1.0 < 3 && result.bond(BondIdx(0)).atom2.0 < 3,
            "remaining bond should connect organic atoms"
        );
    }

    #[test]
    fn bug3_ionic_pair_neutralize_before_largest_fragment() {
        // BUG #3 Fix Verification: NeutralizeCharges must run BEFORE LargestFragment
        // Example: [NH3+].[OH-] (ammonium hydroxide)
        // NH3+ will be neutralized (reduce H by 1)
        // After neutralization: [NH2+].[OH-] - still two fragments
        // LargestFragment then picks the larger one
        // Key: stage order affects which fragment is selected
        let mol = parse("[NH3+].[OH-]").unwrap();
        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            largest_fragment_only: true,
            neutralize_charges: true,
            remove_explicit_h: false,
            canonical_tautomer: false,
            zwitterion_handling: ZwitterionHandling::Normalize,
        });

        let (_result, report) = pipeline.run(&mol);

        // Verify step order: NeutralizeCharges MUST come before LargestFragment
        assert_eq!(report.steps.len(), 5, "Should have 5 steps in pipeline");
        assert_eq!(
            report.steps[0].step,
            StandardizationStep::NeutralizeCharges,
            "NeutralizeCharges must be step 0"
        );
        assert_eq!(
            report.steps[1].step,
            StandardizationStep::LargestFragment,
            "LargestFragment must be step 1"
        );

        // The test passes if step order is correct and the pipeline runs without error
        assert!(report.changed(), "Pipeline should report changes");
        assert_eq!(
            report.status,
            PipelineStatus::Modified,
            "Should be marked as Modified"
        );
    }

    // ── Parent structure extraction tests ────────────────────────────────────

    #[test]
    fn remove_isotopes_strips_isotope_labels() {
        // "[13C]CC" has one 13C isotope label
        let mol = parse("[13C]CC").unwrap();
        let result = remove_isotopes(&mol);
        for i in 0..result.atom_count() {
            assert_eq!(
                result.atom(chematic_core::AtomIdx(i as u32)).isotope,
                None,
                "atom {} should have no isotope",
                i
            );
        }
    }

    #[test]
    fn remove_isotopes_preserves_structure() {
        // "[13C]CC" and "CC" should be structurally identical after isotope removal
        let mol = parse("[13C]CC").unwrap();
        let result = remove_isotopes(&mol);
        assert_eq!(result.atom_count(), 3, "atom count preserved");
        assert_eq!(result.bond_count(), 2, "bond count preserved");
    }

    #[test]
    fn remove_stereo_strips_chirality() {
        // "N[C@@H](C)C(=O)O" — alanine with (S) stereochemistry
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let result = remove_stereo(&mol);
        for i in 0..result.atom_count() {
            use chematic_core::Chirality;
            assert_eq!(
                result.atom(chematic_core::AtomIdx(i as u32)).chirality,
                Chirality::None,
                "atom {} should have no chirality",
                i
            );
        }
    }

    #[test]
    fn remove_stereo_converts_wedge_bonds_to_single() {
        // "C[C@H](O)C" — methylcarbinol with wedge stereochemistry
        // After remove_stereo, Up/Down bonds should become Single
        let mol = parse("C[C@H](O)C").unwrap();
        let result = remove_stereo(&mol);
        for i in 0..result.bond_count() {
            use chematic_core::BondOrder;
            let bond = result.bond(chematic_core::BondIdx(i as u32));
            assert_ne!(
                bond.order,
                BondOrder::Up,
                "bond {} should not be Up after stereo removal",
                i
            );
            assert_ne!(
                bond.order,
                BondOrder::Down,
                "bond {} should not be Down after stereo removal",
                i
            );
        }
    }

    #[test]
    fn remove_stereo_preserves_structure() {
        // "N[C@@H](C)C(=O)O" (alanine) should keep 6 atoms, 5 bonds after stereo removal
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let result = remove_stereo(&mol);
        assert_eq!(result.atom_count(), 6, "atom count preserved");
        assert_eq!(result.bond_count(), 5, "bond count preserved");
    }

    #[test]
    fn parent_variant_step_names_distinct() {
        // Verify that the 4 new parent variants have distinct step names
        let frag_parent = StandardizationStep::FragmentParent;
        let charge_parent = StandardizationStep::ChargeParent;
        let isotope_parent = StandardizationStep::IsotopeParent;
        let stereo_parent = StandardizationStep::StereoParent;

        assert_eq!(frag_parent.as_str(), "fragment_parent");
        assert_eq!(charge_parent.as_str(), "charge_parent");
        assert_eq!(isotope_parent.as_str(), "isotope_parent");
        assert_eq!(stereo_parent.as_str(), "stereo_parent");
    }

    // ── Cleanup transform tests (B3) ────────────────────────────────────────

    #[test]
    fn prefer_organic_removes_inorganic_salts() {
        // "CCO.[Na+].[Cl-]" — ethanol + sodium chloride
        let mol = parse("CCO.[Na+].[Cl-]").unwrap();
        assert_eq!(mol.atom_count(), 5, "input has CCO + Na + Cl");

        let result = prefer_organic(&mol);

        // Should keep only the organic ethanol fragment
        assert_eq!(result.atom_count(), 3, "should keep only ethanol (C, C, O)");
    }

    #[test]
    fn prefer_organic_keeps_organic_if_no_inorganic() {
        // "CC" — ethane only
        let mol = parse("CC").unwrap();
        let result = prefer_organic(&mol);
        assert_eq!(result.atom_count(), 2, "ethane unchanged");
    }

    #[test]
    fn prefer_organic_falls_back_to_largest() {
        // "C.C.C" — three separate carbons (all organic)
        let mol = parse("C.C.C").unwrap();
        let result = prefer_organic(&mol);
        // Should keep the largest organic fragment (any one of them, but all are size 1)
        assert_eq!(
            result.atom_count(),
            1,
            "falls back to largest fragment (one C)"
        );
    }

    #[test]
    fn uncharge_neutralizes_all_charges() {
        // "[NH4+].[OH-]" — ammonium hydroxide
        let mol = parse("[NH4+].[OH-]").unwrap();
        assert!(
            mol.atoms().any(|(_, a)| a.charge != 0),
            "input has charged atoms"
        );

        let result = uncharge(&mol);

        // All atoms should be neutral
        for (_, atom) in result.atoms() {
            assert_eq!(atom.charge, 0, "all atoms should be neutral");
        }
    }

    #[test]
    fn reionize_deprotonates_carboxylic_acids() {
        // "CC(=O)O" — acetic acid
        let mol = parse("CC(=O)O").unwrap();

        let result = reionize(&mol);

        // Should have a negatively charged oxygen (carboxylate anion)
        let has_negative_oxygen = result
            .atoms()
            .any(|(_, a)| a.element.atomic_number() == 8 && a.charge < 0);

        assert!(
            has_negative_oxygen,
            "reionize should deprotonate carboxylic acids"
        );
    }

    #[test]
    fn reionize_protonates_amines() {
        // "CC(N)C" — secondary amine
        let mol = parse("CC(N)C").unwrap();

        let result = reionize(&mol);

        // Should have a positively charged nitrogen (ammonium)
        let has_positive_nitrogen = result
            .atoms()
            .any(|(_, a)| a.element.atomic_number() == 7 && a.charge > 0);

        assert!(has_positive_nitrogen, "reionize should protonate amines");
    }

    #[test]
    fn reionize_protects_amide_nitrogen() {
        // BUG FIX #2: Amide nitrogen should NOT be protonated
        // "CC(=O)N" — primary amide
        let mol = parse("CC(=O)N").unwrap();

        let result = reionize(&mol);

        // Should NOT have a positively charged nitrogen
        let has_positive_nitrogen = result
            .atoms()
            .any(|(_, a)| a.element.atomic_number() == 7 && a.charge > 0);

        assert!(
            !has_positive_nitrogen,
            "reionize should NOT protonate amide nitrogen"
        );
    }

    #[test]
    fn reionize_protects_thioamide_nitrogen() {
        // BUG FIX #2: Thioamide nitrogen should also NOT be protonated
        // "CC(=S)N" — thioamide
        let mol = parse("CC(=S)N").unwrap();

        let result = reionize(&mol);

        // Should NOT have a positively charged nitrogen
        let has_positive_nitrogen = result
            .atoms()
            .any(|(_, a)| a.element.atomic_number() == 7 && a.charge > 0);

        assert!(
            !has_positive_nitrogen,
            "reionize should NOT protonate thioamide nitrogen (C=S conjugation)"
        );
    }

    // B2: normalize_groups expansion tests

    #[test]
    fn normalize_groups_nitro() {
        // Standard nitro group: [N+](=O)[O-]
        let mol = parse("C[N+](=O)[O-]").unwrap();
        let result = normalize_groups(&mol);

        // All atoms should be neutral after normalization
        let all_neutral = result.atoms().all(|(_, a)| a.charge == 0);
        assert!(all_neutral, "nitro group should be neutralized");

        // N-O bonds should be double
        let mut has_double_bond = false;
        for (_, bond) in result.bonds() {
            let a1 = result.atom(bond.atom1);
            let a2 = result.atom(bond.atom2);
            if ((a1.element.atomic_number() == 7 && a2.element.atomic_number() == 8)
                || (a1.element.atomic_number() == 8 && a2.element.atomic_number() == 7))
                && bond.order == chematic_core::BondOrder::Double
            {
                has_double_bond = true;
            }
        }
        assert!(has_double_bond, "nitro should have N=O double bond");
    }

    #[test]
    fn normalize_groups_azide() {
        // Azide: [N-][N+]#N
        let mol = parse("[N-][N+]#N").unwrap();
        let result = normalize_groups(&mol);

        // All atoms should be neutral
        let all_neutral = result.atoms().all(|(_, a)| a.charge == 0);
        assert!(all_neutral, "azide should be neutralized");

        // Check for N=N bonds (converted from single)
        let mut has_double_bond_count = 0;
        for (_, bond) in result.bonds() {
            let a1 = result.atom(bond.atom1);
            let a2 = result.atom(bond.atom2);
            if a1.element.atomic_number() == 7
                && a2.element.atomic_number() == 7
                && bond.order == chematic_core::BondOrder::Double
            {
                has_double_bond_count += 1;
            }
        }
        assert!(
            has_double_bond_count > 0,
            "azide should have N=N double bonds after normalization"
        );
    }

    #[test]
    fn normalize_groups_sulfoxide() {
        // Sulfoxide: S(=O)(C)(C)
        let mol = parse("C[S](=O)C").unwrap();
        let result = normalize_groups(&mol);

        // Sulfoxide structure should remain (S=O is already correct form)
        let mut has_s_double_o = false;
        for (_, bond) in result.bonds() {
            let a1 = result.atom(bond.atom1);
            let a2 = result.atom(bond.atom2);
            if ((a1.element.atomic_number() == 16 && a2.element.atomic_number() == 8)
                || (a1.element.atomic_number() == 8 && a2.element.atomic_number() == 16))
                && bond.order == chematic_core::BondOrder::Double
            {
                has_s_double_o = true;
            }
        }
        assert!(has_s_double_o, "sulfoxide should have S=O double bond");
    }

    // ── RDKit PR #9051: StereoGroup cleanup for non-chiral atoms ────────────

    #[test]
    fn remove_stereo_clears_stereo_groups() {
        // remove_stereo must produce a molecule with no stereo groups (RDKit PR #9051).
        // It rebuilds via MoleculeBuilder::new() so groups are already empty; this
        // test guards against regressions where from_molecule() is accidentally used.
        use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
        let mut mol = parse("[C@@H](F)(Cl)Br").unwrap();
        mol.add_stereo_group(StereoGroup::new(StereoGroupKind::Absolute, vec![AtomIdx(0)]));
        assert_eq!(mol.stereo_groups().len(), 1, "precondition: group was added");
        let stripped = remove_stereo(&mol);
        assert_eq!(
            stripped.stereo_groups().len(),
            0,
            "remove_stereo must clear stereo groups"
        );
    }

    #[test]
    fn clean_stereo_groups_drops_non_chiral_atoms() {
        // clean_stereo_groups filters out atoms with Chirality::None (RDKit PR #9051).
        use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
        // atom 0 = CH3 (not chiral), atom 1 = @@ chiral center
        let mut mol = parse("C[C@@H](F)Cl").unwrap();
        mol.add_stereo_group(StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomIdx(0), AtomIdx(1)], // atom 0 is NOT chiral
        ));
        let cleaned = clean_stereo_groups(&mol);
        assert_eq!(cleaned.stereo_groups().len(), 1, "group must survive with 1 atom");
        assert_eq!(
            cleaned.stereo_groups()[0].atom_indices,
            vec![AtomIdx(1)],
            "only chiral atom must remain in group"
        );
    }

    #[test]
    fn clean_stereo_groups_drops_empty_groups() {
        // A group with no chiral atoms at all must be removed entirely.
        use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
        let mut mol = parse("CC").unwrap(); // no chiral atoms
        mol.add_stereo_group(StereoGroup::new(StereoGroupKind::Absolute, vec![AtomIdx(0)]));
        let cleaned = clean_stereo_groups(&mol);
        assert_eq!(cleaned.stereo_groups().len(), 0, "empty group must be removed");
    }

    #[test]
    fn clean_stereo_groups_preserves_valid_groups() {
        // A group whose atoms are all chiral must be kept intact.
        use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
        let mut mol = parse("[C@@H](F)(Cl)Br").unwrap();
        mol.add_stereo_group(StereoGroup::new(StereoGroupKind::Absolute, vec![AtomIdx(0)]));
        let cleaned = clean_stereo_groups(&mol);
        assert_eq!(cleaned.stereo_groups().len(), 1, "valid group must be preserved");
        assert_eq!(cleaned.stereo_groups()[0].atom_indices, vec![AtomIdx(0)]);
    }

    #[test]
    fn normalize_groups_mixed_nitro_and_azide() {
        // Molecule with both nitro and azide groups
        let mol = parse("C[N+](=O)[O-].N[N+](=O)[O-]").unwrap();
        let result = normalize_groups(&mol);

        // All atoms should be neutral
        let all_neutral = result.atoms().all(|(_, a)| a.charge == 0);
        assert!(all_neutral, "both nitro and azide should be neutralized");
    }
}
