//! Hückel aromaticity perception with antiaromaticity detection.
//!
//! Works on kekulized molecules (no `Aromatic` bond orders) **or** on molecules
//! that retain `Aromatic` bond orders from the SMILES parser (pre-kekulization).
//! Call `kekulize` + `apply_kekule` from `chematic-core` before calling
//! `assign_aromaticity` if you need the explicit double-bond form.
//!
//! Algorithm:
//! 1. Find all SSSR rings via `find_sssr`.
//! 2. **Pass 1**: evaluate each ring independently using Hückel electron counting.
//!    Aromatic (`BondOrder::Aromatic`) bonds are treated equivalently to double bonds
//!    so that pre-kekulization input is handled correctly.
//!    A special "bridgehead N" rule covers fused-ring N atoms whose entire valence
//!    is satisfied by single σ-bonds (like indolizine's junction nitrogen).
//! 3. **Pass 2**: iterative propagation. Rings that were `NonAromatic` or
//!    indeterminate in Pass 1 are re-evaluated using the already-aromatic atom set
//!    as context: confirmed-aromatic atoms contribute 1π unconditionally, allowing
//!    fused rings to be recognised bottom-up (e.g. the 6-ring of indolizine).
//! 4. Classify rings by electron count:
//!    - 4n+2 electrons (n >= 0): aromatic (favorable)
//!    - 4n electrons (n > 0): antiaromatic (unfavorable, strongly disfavored)
//!    - Other: non-aromatic
//! 5. Record all aromatic atoms, bonds, and antiaromatic rings in an `AromaticityModel`.

// ---------------------------------------------------------------------------
// Algorithm selector
// ---------------------------------------------------------------------------

/// Algorithm used to classify ring aromaticity.
///
/// Passed to [`assign_aromaticity_ex`] and [`apply_aromaticity_ex`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum AromaticityAlgorithm {
    /// Strict Hückel 4n+2 rule (default). Supports C, N, O, S.
    #[default]
    Huckel,
    /// RDKit-compatible extension. Adds Se (34) and Te (52) as chalcogen lone-pair
    /// donors (2π), matching the RDKit DEFAULT aromaticity model for common
    /// organic and chalcogen heteroaromatics.
    ///
    /// P-containing aromatic rings are NOT supported in this mode (separate sprint).
    /// Keto-lactam aromaticity is NOT included (TautomerMode, separate sprint).
    RdkitLike,
}

use rustc_hash::{FxHashMap, FxHashSet};

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, implicit_hcount};

use crate::sssr::find_sssr;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Ring aromaticity classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingAromaticity {
    /// 4n+2 electrons: aromatic (favorable)
    Aromatic,
    /// 4n electrons (n > 0): antiaromatic (unfavorable)
    Antiaromatic,
    /// Any other electron count: non-aromatic
    NonAromatic,
}

/// Aromaticity assignment for a molecule.
///
/// Records which atoms and bonds belong to aromatic rings according to
/// the Hückel 4n+2 rule applied to SSSR rings (with fused-ring propagation).
/// Also tracks antiaromatic rings (4n electrons) for chemical accuracy.
#[derive(Debug, Clone)]
pub struct AromaticityModel {
    aromatic_atoms: FxHashSet<AtomIdx>,
    aromatic_bonds: FxHashSet<BondIdx>,
    antiaromatic_rings: Vec<Vec<AtomIdx>>,
    ring_classifications: Vec<(Vec<AtomIdx>, RingAromaticity, u32)>,
}

impl AromaticityModel {
    /// Whether atom `idx` is part of an aromatic ring.
    pub fn is_atom_aromatic(&self, idx: AtomIdx) -> bool {
        self.aromatic_atoms.contains(&idx)
    }

    /// Whether bond `idx` is part of an aromatic ring.
    pub fn is_bond_aromatic(&self, idx: BondIdx) -> bool {
        self.aromatic_bonds.contains(&idx)
    }

    /// Total number of atoms flagged as aromatic.
    pub fn aromatic_atom_count(&self) -> usize {
        self.aromatic_atoms.len()
    }

    /// Get all rings and their classification with electron counts.
    ///
    /// Each entry is `(ring_atoms, classification, π_electron_count)`.
    /// Rings that could not be evaluated (sp3 atoms, unsupported elements) are omitted.
    pub fn ring_classifications(&self) -> &[(Vec<AtomIdx>, RingAromaticity, u32)] {
        &self.ring_classifications
    }

    /// Get all antiaromatic rings (4n electrons, n > 0).
    pub fn antiaromatic_rings(&self) -> &[Vec<AtomIdx>] {
        &self.antiaromatic_rings
    }

    /// Check if any atom belongs to an antiaromatic ring.
    pub fn has_antiaromaticity(&self) -> bool {
        !self.antiaromatic_rings.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Main entry points
// ---------------------------------------------------------------------------

/// Classify a ring by its pi electron count using Hückel and antiaromaticity rules.
#[allow(clippy::manual_is_multiple_of)]
fn classify_ring_aromaticity(pi_electrons: u32) -> (RingAromaticity, u32) {
    if pi_electrons >= 2 && (pi_electrons - 2) % 4 == 0 {
        (RingAromaticity::Aromatic, pi_electrons)
    } else if pi_electrons > 0 && pi_electrons % 4 == 0 {
        (RingAromaticity::Antiaromatic, pi_electrons)
    } else {
        (RingAromaticity::NonAromatic, pi_electrons)
    }
}

/// Mark all atoms and bonds in `ring` as aromatic in the provided sets.
fn mark_ring_aromatic(
    mol: &Molecule,
    ring: &[AtomIdx],
    aromatic_atoms: &mut FxHashSet<AtomIdx>,
    aromatic_bonds: &mut FxHashSet<BondIdx>,
) {
    for &atom in ring {
        aromatic_atoms.insert(atom);
    }
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        if let Some((bidx, _)) = mol.bond_between(a, b) {
            aromatic_bonds.insert(bidx);
        }
    }
}

/// Assign aromaticity to a molecule using the Hückel 4n+2 rule with fused-ring
/// propagation (Pass 2) and antiaromaticity detection (4n electrons).
///
/// The molecule may be kekulized (`Single`/`Double` bonds) **or** may retain
/// `BondOrder::Aromatic` bonds from the SMILES parser.  In the latter case,
/// aromatic bonds are treated as equivalent to double bonds for electron
/// counting, allowing correct detection without an explicit kekulization step.
///
/// For kekulized input from aromatic SMILES, call `chematic_core::kekulize`
/// then `chematic_core::apply_kekule` first.
///
/// Uses [`AromaticityAlgorithm::Huckel`] (default). See [`assign_aromaticity_ex`]
/// for the RdkitLike variant.
pub fn assign_aromaticity(mol: &Molecule) -> AromaticityModel {
    assign_aromaticity_ex(mol, AromaticityAlgorithm::Huckel)
}

/// Assign aromaticity using the specified algorithm.
///
/// The default ([`assign_aromaticity`]) uses [`AromaticityAlgorithm::Huckel`].
/// Pass [`AromaticityAlgorithm::RdkitLike`] to additionally recognise Se/Te
/// as lone-pair donors in aromatic rings.
pub fn assign_aromaticity_ex(mol: &Molecule, algo: AromaticityAlgorithm) -> AromaticityModel {
    let ring_set = find_sssr(mol);
    let sssr_rings = ring_set.rings();

    // Augment SSSR rings with smaller XOR sub-rings (GF(2) differences between pairs).
    // This corrects the case where the SSSR algorithm stores a large fundamental cycle
    // instead of its smaller GF(2)-reduced equivalent (e.g. the 5-ring of indolizine).
    let rings: Vec<Vec<AtomIdx>> = augmented_ring_set(mol, sssr_rings);

    let mut aromatic_atoms: FxHashSet<AtomIdx> = FxHashSet::default();
    let mut aromatic_bonds: FxHashSet<BondIdx> = FxHashSet::default();
    let mut antiaromatic_rings: Vec<Vec<AtomIdx>> = Vec::new();

    // Per-ring classification: None means "not yet evaluated / indeterminate".
    let mut classifications: Vec<Option<(RingAromaticity, u32)>> = vec![None; rings.len()];

    // Indices of rings that are candidates for Pass 2 re-evaluation
    // (returned None or NonAromatic in Pass 1).
    let mut pass2_candidates: Vec<usize> = Vec::new();

    // ----- Pass 1: independent Hückel per ring -----
    let empty_context = FxHashSet::default();
    for (ring_idx, ring) in rings.iter().enumerate() {
        match ring_pi_electrons(mol, ring, &empty_context, algo) {
            Some(pi) => {
                let (cls, count) = classify_ring_aromaticity(pi);
                classifications[ring_idx] = Some((cls, count));
                match cls {
                    RingAromaticity::Aromatic => {
                        mark_ring_aromatic(mol, ring, &mut aromatic_atoms, &mut aromatic_bonds);
                    }
                    RingAromaticity::Antiaromatic => {
                        antiaromatic_rings.push(ring.to_vec());
                        // Antiaromatic is definitive — do not retry in Pass 2.
                    }
                    RingAromaticity::NonAromatic => {
                        pass2_candidates.push(ring_idx);
                    }
                }
            }
            None => {
                // Indeterminate (sp3 atoms, unsupported elements, etc.).
                pass2_candidates.push(ring_idx);
            }
        }
    }

    // ----- Pass 2: propagate through fused ring systems -----
    // Re-evaluate rings adjacent to already-aromatic rings.  Repeat until
    // convergence (no newly aromatic ring found in the last iteration).
    loop {
        let mut any_new = false;
        let mut still_pending: Vec<usize> = Vec::new();

        for ring_idx in pass2_candidates {
            let ring = &rings[ring_idx];
            // Only rings that share an atom with an already-aromatic ring qualify.
            if !ring.iter().any(|a| aromatic_atoms.contains(a)) {
                still_pending.push(ring_idx);
                continue;
            }
            match ring_pi_electrons(mol, ring, &aromatic_atoms, algo) {
                Some(pi) => {
                    let (cls, count) = classify_ring_aromaticity(pi);
                    classifications[ring_idx] = Some((cls, count));
                    if matches!(cls, RingAromaticity::Aromatic) {
                        mark_ring_aromatic(mol, ring, &mut aromatic_atoms, &mut aromatic_bonds);
                        any_new = true;
                    }
                    // NonAromatic even in Pass 2 context: do not retry further.
                }
                None => {
                    still_pending.push(ring_idx);
                }
            }
        }

        pass2_candidates = still_pending;
        if !any_new {
            break;
        }
    }

    // Build the public ring_classifications list (SSSR rings only, omitting augmented/indeterminate).
    let ring_classifications: Vec<(Vec<AtomIdx>, RingAromaticity, u32)> = rings
        .iter()
        .take(sssr_rings.len()) // only expose SSSR rings in the public API
        .enumerate()
        .filter_map(|(i, ring)| classifications[i].map(|(cls, count)| (ring.to_vec(), cls, count)))
        .collect();

    AromaticityModel {
        aromatic_atoms,
        aromatic_bonds,
        antiaromatic_rings,
        ring_classifications,
    }
}

/// Apply aromaticity perception to a molecule.
///
/// Returns a new [`Molecule`] where atoms in Hückel-aromatic rings have
/// `atom.aromatic = true` and their bonds carry [`BondOrder::Aromatic`].
/// Non-aromatic atoms and bonds are unchanged.
///
/// The input may be kekulized (no `Aromatic` bond orders) or may retain
/// aromatic bond orders from the SMILES parser.
///
/// Uses [`AromaticityAlgorithm::Huckel`] (default). See [`apply_aromaticity_ex`]
/// for the RdkitLike variant.
pub fn apply_aromaticity(mol: &Molecule) -> Molecule {
    apply_aromaticity_ex(mol, AromaticityAlgorithm::Huckel)
}

/// Apply aromaticity using the specified algorithm.
///
/// Returns a new [`Molecule`] with aromatic flags set according to `algo`.
pub fn apply_aromaticity_ex(mol: &Molecule, algo: AromaticityAlgorithm) -> Molecule {
    use chematic_core::{BondOrder, MoleculeBuilder};

    let model = assign_aromaticity_ex(mol, algo);
    let mut builder = MoleculeBuilder::new();

    for (idx, atom) in mol.atoms() {
        let mut a = atom.clone();
        if model.is_atom_aromatic(idx) {
            a.aromatic = true;
        }
        builder.add_atom(a);
    }
    for (bidx, bond) in mol.bonds() {
        let order = if model.is_bond_aromatic(bidx) {
            BondOrder::Aromatic
        } else {
            bond.order
        };
        if let Ok(new_bidx) = builder.add_bond(bond.atom1, bond.atom2, order)
            && order == BondOrder::Aromatic
            && matches!(bond.order, BondOrder::Up | BondOrder::Down)
        {
            // Kekule input promoted to Aromatic here loses its E/Z direction
            // the same way the SMILES parser's aromatic-aromatic coercion
            // does — stash it so an exocyclic double bond anchored on this
            // ring bond still round-trips through the canonical writer.
            builder.set_bond_direction(new_bidx, bond.order);
        }
    }
    // Atoms/bonds above are re-added in `mol`'s own enumeration order with
    // none skipped, so indices line up 1:1 — safe to copy side-channel
    // metadata wholesale. (This rebuild previously dropped stereo_groups and
    // stereo_neighbor_order silently; closing that here too.)
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.build()
}

// ---------------------------------------------------------------------------
// Ring augmentation (XOR sub-rings)
// ---------------------------------------------------------------------------

/// Return the sorted set of bond indices that form `ring`.
fn ring_bond_set(mol: &Molecule, ring: &[AtomIdx]) -> Vec<BondIdx> {
    let n = ring.len();
    let mut bonds: Vec<BondIdx> = (0..n)
        .filter_map(|i| {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            mol.bond_between(a, b).map(|(bidx, _)| bidx)
        })
        .collect();
    bonds.sort();
    bonds
}

/// Sorted symmetric difference of two sorted slices.
fn bond_sym_diff(a: &[BondIdx], b: &[BondIdx]) -> Vec<BondIdx> {
    let mut result: Vec<BondIdx> = Vec::new();
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                result.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                result.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    result.extend_from_slice(&a[i..]);
    result.extend_from_slice(&b[j..]);
    result
}

/// Reconstruct an ordered atom sequence from a set of bond indices forming a simple cycle.
/// Returns `None` if the bonds do not form a valid simple cycle.
fn ring_atoms_from_bond_set(mol: &Molecule, bonds: &[BondIdx]) -> Option<Vec<AtomIdx>> {
    if bonds.is_empty() {
        return None;
    }
    let mut adj: FxHashMap<AtomIdx, [Option<AtomIdx>; 2]> = FxHashMap::default();
    for &bidx in bonds {
        let bond = mol.bond(bidx);
        for (a, b) in [(bond.atom1, bond.atom2), (bond.atom2, bond.atom1)] {
            let e = adj.entry(a).or_insert([None; 2]);
            if e[0].is_none() {
                e[0] = Some(b);
            } else if e[1].is_none() {
                e[1] = Some(b);
            } else {
                return None; // degree > 2 — not a simple ring
            }
        }
    }
    // All atoms must have exactly 2 neighbours.
    if adj.values().any(|e| e[1].is_none()) {
        return None;
    }
    let start = *adj.keys().next()?;
    let mut path = vec![start];
    let mut prev = start;
    let mut current = adj[&start][0]?;
    while current != start {
        path.push(current);
        let [n0, n1] = adj[&current];
        let next = if n0 == Some(prev) { n1? } else { n0? };
        prev = current;
        current = next;
    }
    if path.len() != bonds.len() {
        return None;
    }
    Some(path)
}

/// Augment the SSSR ring list with smaller XOR sub-rings found by pairwise GF(2)
/// differences between SSSR rings that share atoms.
///
/// The standard SSSR algorithm sometimes stores a large fundamental cycle rather
/// than its smaller GF(2)-reduced equivalent (e.g. the 5-ring of indolizine is
/// the XOR of the 6-ring and the 9-ring the algorithm reports).
/// This augmentation adds such missing smaller rings so that aromaticity
/// perception works on the correct smallest rings without modifying the SSSR.
///
/// The returned `Vec` starts with all SSSR rings in their original order; any
/// additional sub-rings derived by GF(2) pairwise XOR follow.  The function
/// only adds a ring if it is strictly smaller than *both* parents, ensuring
/// that envelope rings (e.g. the 10-membered perimeter of naphthalene) are
/// never introduced.
pub fn augmented_ring_set(mol: &Molecule, sssr_rings: &[Vec<AtomIdx>]) -> Vec<Vec<AtomIdx>> {
    let mut rings: Vec<Vec<AtomIdx>> = sssr_rings.to_vec();

    // Track which atom-sets we already have (as sorted atom lists).
    let mut known: FxHashSet<Vec<AtomIdx>> = sssr_rings
        .iter()
        .map(|r| {
            let mut s = r.clone();
            s.sort();
            s
        })
        .collect();

    // Iterative pairwise XOR until convergence.
    //
    // A single pass only finds rings that are the XOR of two SSSR rings.
    // Iterating also finds rings that require XOR of 3+ SSSR rings
    // (e.g. the inner hexagon of coronene, or sub-rings in multi-step
    // fused PAHs where the SSSR chose large perimeter cycles).
    // Termination is guaranteed because each new ring is strictly smaller
    // than both of its parents, so ring size can only decrease.
    loop {
        let mut changed = false;
        let n = rings.len();
        let bond_sets: Vec<Vec<BondIdx>> = rings.iter().map(|r| ring_bond_set(mol, r)).collect();

        for i in 0..n {
            for j in (i + 1)..n {
                // Only consider pairs that share atoms (fused rings).
                let shares_atom = rings[i].iter().any(|a| rings[j].contains(a));
                if !shares_atom {
                    continue;
                }
                let xor_bonds = bond_sym_diff(&bond_sets[i], &bond_sets[j]);
                if xor_bonds.is_empty() {
                    continue;
                }
                // Only interesting if the XOR ring is not larger than the larger
                // parent.  Using max() recovers cases where SSSR chose a large
                // cycle (e.g. 10-ring macro vs 6-ring benzene twin).
                // Using `>` (not `>=`) also allows same-size XOR rings, which
                // handles bridged bicyclics (e.g. tropane or dioxolane spirocycles)
                // where both parent rings are 6-membered and the missing bridge
                // ring is also 6-membered.  Termination is still guaranteed:
                // the `known` set prevents duplicates, and a finite molecule has
                // finitely many valid cycles.
                if xor_bonds.len() > rings[i].len().max(rings[j].len()) {
                    continue;
                }
                if let Some(new_ring) = ring_atoms_from_bond_set(mol, &xor_bonds) {
                    let mut key = new_ring.clone();
                    key.sort();
                    if known.insert(key) {
                        rings.push(new_ring);
                        changed = true;
                    }
                }
            }
        }

        // 3-ring XOR: catches small rings that require XOR of 3 SSSR rings
        // when no intermediate 2-ring XOR produces a valid smaller ring.
        for i in 0..n {
            for j in (i + 1)..n {
                let shares_ij = rings[i].iter().any(|a| rings[j].contains(a));
                if !shares_ij {
                    continue;
                }
                let xor_ij = bond_sym_diff(&bond_sets[i], &bond_sets[j]);
                if xor_ij.is_empty() {
                    continue;
                }
                for k in (j + 1)..n {
                    let shares_k = rings[k]
                        .iter()
                        .any(|a| rings[i].contains(a) || rings[j].contains(a));
                    if !shares_k {
                        continue;
                    }
                    let xor_ijk = bond_sym_diff(&xor_ij, &bond_sets[k]);
                    let max_size = rings[i].len().max(rings[j].len()).max(rings[k].len());
                    if xor_ijk.is_empty() || xor_ijk.len() > max_size {
                        continue;
                    }
                    if let Some(new_ring) = ring_atoms_from_bond_set(mol, &xor_ijk) {
                        let mut key = new_ring.clone();
                        key.sort();
                        if known.insert(key) {
                            rings.push(new_ring);
                            changed = true;
                        }
                    }
                }
            }
        }

        if !changed {
            break;
        }
    }

    rings
}

/// Shared inner: SSSR → augmented_ring_set → strip_envelope_rings, no aromaticity filter.
fn all_ring_list_inner(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let sssr = crate::sssr::find_sssr(mol);
    let aug = augmented_ring_set(mol, sssr.rings());
    if aug.len() <= 1 {
        return aug;
    }
    let bond_sets: Vec<Vec<BondIdx>> = aug.iter().map(|r| ring_bond_set(mol, r)).collect();
    let mut is_envelope = vec![false; aug.len()];
    strip_envelope_rings(&aug, &bond_sets, &mut is_envelope);
    aug.into_iter()
        .zip(is_envelope)
        .filter(|(_, e)| !e)
        .map(|(r, _)| r)
        .collect()
}

/// Return all rings after augmented-ring-set expansion and envelope stripping.
///
/// Same pipeline as [`aromatic_ring_list`] but with no aromaticity filter — useful
/// for aliphatic/saturated ring counting and bridgehead detection where SSSR
/// envelope rings cause over-counting.
pub fn all_ring_list(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    all_ring_list_inner(mol)
}

/// True when all ring bonds between ring atoms are `BondOrder::Aromatic`.
///
/// Rings written with aromatic-SMILES notation but containing an explicit single
/// bond (`c-n`, `nc-2`, etc.) are NOT truly aromatic.  RDKit canonicalises such
/// SMILES with lowercase atoms and a `-` bond, which the parser stores as
/// `BondOrder::Single` between two aromatic-flagged atoms.  Returning `false`
/// here lets callers exclude them from the aromatic ring count.
pub fn ring_bonds_all_aromatic(mol: &Molecule, ring: &[AtomIdx]) -> bool {
    let n = ring.len();
    (0..n).all(|i| {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        mol.bond_between(a, b)
            .map(|(bidx, _)| mol.bond(bidx).order == BondOrder::Aromatic)
            .unwrap_or(true)
    })
}

/// Return the de-duplicated list of aromatic rings after augmented-ring-set expansion
/// and envelope stripping.  Useful for filtering (e.g. counting only aromatic heterocycles).
pub fn aromatic_ring_list(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let mol_with_arom;
    let mol = if mol.atoms().any(|(_, a)| a.aromatic) {
        mol
    } else {
        mol_with_arom = apply_aromaticity(mol);
        &mol_with_arom
    };
    all_ring_list_inner(mol)
        .into_iter()
        .filter(|ring| {
            ring.iter().all(|&idx| mol.atom(idx).aromatic) && ring_bonds_all_aromatic(mol, ring)
        })
        .collect()
}

/// Mark which rings in `aromatic` are GF(2) sums (bond-XOR) of 2–4 smaller rings.
fn strip_envelope_rings(
    aromatic: &[Vec<AtomIdx>],
    bond_sets: &[Vec<BondIdx>],
    is_envelope: &mut [bool],
) {
    let n = aromatic.len();
    for i in 0..n {
        let si = aromatic[i].len();
        'jk: for j in 0..n {
            if j == i || aromatic[j].len() >= si {
                continue;
            }
            for k in (j + 1)..n {
                if k == i || aromatic[k].len() >= si {
                    continue;
                }
                if bond_sym_diff(&bond_sets[j], &bond_sets[k]) == bond_sets[i] {
                    is_envelope[i] = true;
                    break 'jk;
                }
            }
        }
        if !is_envelope[i] {
            'jkl: for j in 0..n {
                if j == i || aromatic[j].len() >= si {
                    continue;
                }
                for k in (j + 1)..n {
                    if k == i || aromatic[k].len() >= si {
                        continue;
                    }
                    let xor_jk = bond_sym_diff(&bond_sets[j], &bond_sets[k]);
                    for l in (k + 1)..n {
                        if l == i || aromatic[l].len() >= si {
                            continue;
                        }
                        if bond_sym_diff(&xor_jk, &bond_sets[l]) == bond_sets[i] {
                            is_envelope[i] = true;
                            break 'jkl;
                        }
                    }
                }
            }
        }
        if !is_envelope[i] {
            'jklm: for j in 0..n {
                if j == i || aromatic[j].len() >= si {
                    continue;
                }
                for k in (j + 1)..n {
                    if k == i || aromatic[k].len() >= si {
                        continue;
                    }
                    let xor_jk = bond_sym_diff(&bond_sets[j], &bond_sets[k]);
                    for l in (k + 1)..n {
                        if l == i || aromatic[l].len() >= si {
                            continue;
                        }
                        let xor_jkl = bond_sym_diff(&xor_jk, &bond_sets[l]);
                        for m in (l + 1)..n {
                            if m == i || aromatic[m].len() >= si {
                                continue;
                            }
                            if bond_sym_diff(&xor_jkl, &bond_sets[m]) == bond_sets[i] {
                                is_envelope[i] = true;
                                break 'jklm;
                            }
                        }
                    }
                }
            }
        }
    }
}

pub fn count_aromatic_rings(mol: &Molecule) -> usize {
    // For Kekulé-form input (uppercase atoms, no aromatic flags yet), run Hückel
    // perception first so ring detection works correctly (RDKit #9271).
    let mol_with_arom;
    let mol = if mol.atoms().any(|(_, a)| a.aromatic) {
        mol // aromatic SMILES — flags already set during parsing
    } else {
        mol_with_arom = apply_aromaticity(mol);
        &mol_with_arom
    };

    let sssr = crate::sssr::find_sssr(mol);
    let aug = augmented_ring_set(mol, sssr.rings());

    // Keep only rings where every atom carries the aromatic flag.
    let aromatic: Vec<Vec<AtomIdx>> = aug
        .into_iter()
        .filter(|ring| ring.iter().all(|&idx| mol.atom(idx).aromatic))
        .collect();

    if aromatic.len() <= 1 {
        return aromatic.len();
    }

    // Build sorted bond-index sets for each aromatic ring.
    let bond_sets: Vec<Vec<BondIdx>> = aromatic.iter().map(|r| ring_bond_set(mol, r)).collect();

    // Mark rings that are the GF(2) sum (bond-XOR) of 2, 3, or 4 strictly
    // smaller aromatic rings.  Such rings are "envelope" cycles introduced
    // when the SSSR chose a large fundamental cycle instead of its smaller
    // GF(2) components.
    // 2-ring XOR: handles linear/angular fused systems (naphthalene, indolizine…).
    // 3-ring XOR: handles compact PAHs like pyrene.
    // 4-ring XOR: handles coronene-class PAHs where the outer perimeter is the
    //   GF(2) sum of four inner hexagons.
    let n = aromatic.len();
    let mut is_envelope = vec![false; n];
    strip_envelope_rings(&aromatic, &bond_sets, &mut is_envelope);
    is_envelope.iter().filter(|&&e| !e).count()
}

// ---------------------------------------------------------------------------
// Per-ring pi electron count
// ---------------------------------------------------------------------------

/// Count pi electrons for a ring atom, returning `None` if the atom is
/// incompatible with aromaticity (e.g. sp3 carbon).
///
/// `aromatic_context`: atoms already confirmed aromatic (from Pass 1 or a
/// previous Pass 2 iteration).  Such atoms contribute 1π unconditionally,
/// without requiring an explicit double bond.
///
/// Rules:
/// - **C**: if already in `aromatic_context` → 1π (confirmed sp2).
///   1. No double bond anywhere: carbanion (`charge == -1`) → 2π (lone pair,
///      e.g. cyclopentadienyl anion); otherwise sp3 → None.
///   2. Has a double bond, but only exocyclic and to a more electronegative
///      atom (O/N/S) → 0π (its p-orbital electrons are in the exocyclic π
///      bond, e.g. the carbonyl carbon in tropone/pyridone/pyranone).
///   3. Otherwise (has an endocyclic Double/Aromatic bond) → 1π.
/// - **N**:
///   1. Has H → 2π (pyrrole-type lone pair).
///   2. Has an explicit `Double` bond → 1π (pyridine-type).
///   3. total_degree == 3 AND ring_degree < total_degree AND no explicit
///      double bond → 2π (lone pair in p orbital): covers both a bridgehead
///      N shared by two fused rings (indolizine) and a substituted
///      pyrrole-type N (N-methylpyrrole, N-glycosylated purine); the overall
///      4n+2 sum, not the substituent, decides ring aromaticity.
///   4. Has in-ring `Aromatic` bond → 1π (pyridine-like aromatic N).
///   5. Already in `aromatic_context` → 1π.
///   6. Otherwise → None.
/// - **O/S**: ring_degree must be 2; contributes 2π (lone pair).
/// - **Se (34) / Te (52)**: analogous to S; only in [`AromaticityAlgorithm::RdkitLike`] mode.
/// - **Other elements**: None (unsupported).
fn ring_pi_electrons(
    mol: &Molecule,
    ring: &[AtomIdx],
    aromatic_context: &FxHashSet<AtomIdx>,
    algo: AromaticityAlgorithm,
) -> Option<u32> {
    let ring_atom_set: FxHashSet<AtomIdx> = ring.iter().copied().collect();
    let mut total_pi: u32 = 0;

    for &atom_idx in ring {
        // Atoms already confirmed aromatic in an adjacent ring contribute 1π.
        if aromatic_context.contains(&atom_idx) {
            total_pi += 1;
            continue;
        }

        let atom = mol.atom(atom_idx);
        let an = atom.element.atomic_number();

        let ring_degree = mol
            .neighbors(atom_idx)
            .filter(|(nb, _)| ring_atom_set.contains(nb))
            .count();

        let total_degree = mol.degree(atom_idx);

        // Explicit Double bond anywhere (not counting Aromatic).
        let has_explicit_double = mol
            .neighbors(atom_idx)
            .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Double);

        // Double OR Aromatic bond anywhere (for C sp2 check).
        let has_double_any = has_explicit_double
            || mol
                .neighbors(atom_idx)
                .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Aromatic);

        // Aromatic bond within the ring (for pyridine-like N in aromatic SMILES).
        let has_aromatic_in_ring = mol
            .neighbors(atom_idx)
            .filter(|(nb, _)| ring_atom_set.contains(nb))
            .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Aromatic);

        let pi = match an {
            // Carbon: must be sp2 (has a double or aromatic bond somewhere).
            6 => {
                if !has_double_any {
                    // No double bond: a ring carbanion still donates its lone
                    // pair (e.g. cyclopentadienyl anion), otherwise sp3.
                    if atom.charge == -1 {
                        2
                    } else {
                        return None; // sp3 carbon — ring cannot be aromatic
                    }
                } else if has_explicit_double
                    && !has_aromatic_in_ring
                    && !mol.neighbors(atom_idx).any(|(nb, bidx)| {
                        ring_atom_set.contains(&nb) && mol.bond(bidx).order == BondOrder::Double
                    })
                    && mol.neighbors(atom_idx).any(|(nb, bidx)| {
                        !ring_atom_set.contains(&nb)
                            && mol.bond(bidx).order == BondOrder::Double
                            && matches!(mol.atom(nb).element.atomic_number(), 7 | 8 | 16)
                    })
                {
                    // Only double bond is exocyclic, to a more electronegative
                    // atom (O/N/S): p-orbital electrons sit in that exocyclic π
                    // bond, contributing 0π to the ring (e.g. carbonyl carbon
                    // in tropone/pyridone/pyranone).
                    0
                } else {
                    1
                }
            }

            // Nitrogen
            7 => {
                if implicit_hcount(mol, atom_idx) > 0 {
                    // Pyrrole-type N with H: lone pair → 2π.
                    2
                } else if has_explicit_double {
                    // Pyridine-type N with explicit double bond → 1π.
                    1
                } else if total_degree == 3 && ring_degree < total_degree {
                    // N with no H, no explicit double bond, and all three σ-bonds
                    // exactly filling its valence (3): a bridgehead N shared by two
                    // fused rings (e.g. indolizine) and a substituted pyrrole-type N
                    // (e.g. N-methylpyrrole, N-glycosylated purine/pyrimidine) have
                    // the identical local shape — the lone pair occupies the p
                    // orbital → 2π either way. Whether the ring this atom sits in is
                    // actually aromatic is decided by the overall 4n+2 sum below, not
                    // by inspecting the substituent: an imide N (phthalimide) still
                    // correctly comes out non-aromatic because its ring's carbonyl
                    // carbons contribute 0π each (exocyclic C=O rule above), giving
                    // 4π total, not 4n+2.
                    2
                } else if has_aromatic_in_ring {
                    // N in an aromatic ring (pre-kekulization input) without an
                    // explicit double bond and not a bridgehead → pyridine-like → 1π.
                    1
                } else {
                    // Cannot determine pi contribution.
                    return None;
                }
            }

            // Oxygen / sulfur: lone-pair donor, must be 2-connected in the ring.
            8 | 16 => {
                if ring_degree != 2 {
                    return None;
                }
                // Sulfoxide/sulfone: exocyclic S=O ties up the lone pair; cannot donate 2π
                if an == 16
                    && mol.neighbors(atom_idx).any(|(nb, bidx)| {
                        !ring_atom_set.contains(&nb) && mol.bond(bidx).order == BondOrder::Double
                    })
                {
                    return None;
                }
                2
            }

            // Se (34) / Te (52): chalcogen lone-pair donors (2π), analogous to S.
            // Only recognised in RdkitLike mode.
            34 | 52 => {
                if algo != AromaticityAlgorithm::RdkitLike {
                    return None;
                }
                if ring_degree != 2 {
                    return None;
                }
                // Exocyclic Se=O / Te=O ties up the lone pair.
                if mol.neighbors(atom_idx).any(|(nb, bidx)| {
                    !ring_atom_set.contains(&nb) && mol.bond(bidx).order == BondOrder::Double
                }) {
                    return None;
                }
                2
            }

            // Unsupported element.
            _ => return None,
        };

        total_pi += pi;
    }

    Some(total_pi)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

    // =========================================================================
    // Molecule builder helpers (kekulized, manually constructed)
    // =========================================================================

    fn benzene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(atoms[i], atoms[(i + 1) % 6], order).unwrap();
        }
        b.build()
    }

    fn cyclohexane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Single)
                .unwrap();
        }
        b.build()
    }

    fn pyridine_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(Atom::new(Element::N));
        let atoms_c: Vec<_> = (0..5).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        let ring = [
            n, atoms_c[0], atoms_c[1], atoms_c[2], atoms_c[3], atoms_c[4],
        ];
        for i in 0..6 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(ring[i], ring[(i + 1) % 6], order).unwrap();
        }
        b.build()
    }

    fn furan_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::new(Element::O));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [o, c1, c2, c3, c4];
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        b.build()
    }

    fn pyrrole_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let mut n_atom = Atom::new(Element::N);
        n_atom.hydrogen_count = Some(1);
        let n = b.add_atom(n_atom);
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [n, c1, c2, c3, c4];
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        b.build()
    }

    fn naphthalene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        let ring1 = [0usize, 1, 2, 3, 4, 9];
        let orders1 = [
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ];
        for i in 0..6 {
            b.add_bond(atoms[ring1[i]], atoms[ring1[(i + 1) % 6]], orders1[i])
                .unwrap();
        }
        let ring2_extra = [(4, 5), (5, 6), (6, 7), (7, 8), (8, 9)];
        let orders2 = [
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ];
        for (i, &(a, bb)) in ring2_extra.iter().enumerate() {
            b.add_bond(atoms[a], atoms[bb], orders2[i]).unwrap();
        }
        b.build()
    }

    fn cyclobutadiene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..4 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(atoms[i], atoms[(i + 1) % 4], order).unwrap();
        }
        b.build()
    }

    fn cyclooctatetraene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..8).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..8 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(atoms[i], atoms[(i + 1) % 8], order).unwrap();
        }
        b.build()
    }

    /// Helper: parse an aromatic SMILES and return the molecule with aromatic bonds
    /// (no kekulization).  Use for compounds where kekulization is unsupported.
    #[cfg(test)]
    fn mol_aromatic(smiles: &str) -> chematic_core::Molecule {
        chematic_smiles::parse(smiles).expect("valid SMILES")
    }

    /// Helper: parse SMILES and kekulize.  Panics if kekulization fails.
    #[cfg(test)]
    fn mol_kekulized(smiles: &str) -> chematic_core::Molecule {
        let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
        let k = chematic_core::kekulize(&mol).expect("kekulizable");
        chematic_core::apply_kekule(&mol, &k)
    }

    // =========================================================================
    // Regression: kekulized single-ring aromatics (Pass 1 only, no context)
    // =========================================================================

    #[test]
    fn test_benzene_is_aromatic() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "all 6 benzene atoms aromatic"
        );
        for i in 0..6u32 {
            assert!(model.is_atom_aromatic(AtomIdx(i)));
        }
    }

    #[test]
    fn test_cyclohexane_not_aromatic() {
        let mol = cyclohexane();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 0, "cyclohexane not aromatic");
    }

    #[test]
    fn test_pyridine_is_aromatic() {
        let mol = pyridine_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6);
    }

    #[test]
    fn test_furan_is_aromatic() {
        let mol = furan_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5);
    }

    #[test]
    fn test_pyrrole_is_aromatic() {
        let mol = pyrrole_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5);
    }

    #[test]
    fn test_naphthalene_both_rings_aromatic() {
        let mol = naphthalene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "all 10 naphthalene atoms aromatic"
        );
    }

    #[test]
    fn test_bond_aromaticity_benzene() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        let count = mol
            .bonds()
            .filter(|(b, _)| model.is_bond_aromatic(*b))
            .count();
        assert_eq!(count, 6);
    }

    #[test]
    fn test_apply_aromaticity_benzene() {
        let mol = benzene_kekule();
        let aromatic = apply_aromaticity(&mol);
        for (_, atom) in aromatic.atoms() {
            assert!(atom.aromatic, "every benzene carbon should be aromatic");
        }
        let aromatic_bond_count = aromatic
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Aromatic)
            .count();
        assert_eq!(aromatic_bond_count, 6);
    }

    #[test]
    fn test_apply_aromaticity_cyclohexane_unchanged() {
        let mol = cyclohexane();
        let result = apply_aromaticity(&mol);
        for (_, atom) in result.atoms() {
            assert!(!atom.aromatic);
        }
        for (_, bond) in result.bonds() {
            assert_ne!(bond.order, BondOrder::Aromatic);
        }
    }

    // =========================================================================
    // Antiaromaticity
    // =========================================================================

    #[test]
    fn test_cyclobutadiene_antiaromatic() {
        let mol = cyclobutadiene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            0,
            "cyclobutadiene not aromatic"
        );
        assert!(model.has_antiaromaticity(), "cyclobutadiene antiaromatic");
        assert_eq!(model.antiaromatic_rings().len(), 1);
        let classifications = model.ring_classifications();
        assert_eq!(classifications.len(), 1);
        assert_eq!(classifications[0].1, RingAromaticity::Antiaromatic);
        assert_eq!(classifications[0].2, 4);
    }

    #[test]
    fn test_cyclooctatetraene_antiaromatic() {
        let mol = cyclooctatetraene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 0, "COT not aromatic");
        assert!(model.has_antiaromaticity(), "COT antiaromatic");
        assert_eq!(model.antiaromatic_rings().len(), 1);
        let cls = &model.ring_classifications()[0];
        assert_eq!(cls.1, RingAromaticity::Antiaromatic);
        assert_eq!(cls.2, 8);
    }

    // =========================================================================
    // Ring classifications
    // =========================================================================

    #[test]
    fn test_ring_classifications_benzene() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        let classifications = model.ring_classifications();
        assert_eq!(classifications.len(), 1);
        assert_eq!(classifications[0].1, RingAromaticity::Aromatic);
        assert_eq!(classifications[0].2, 6);
    }

    #[test]
    fn test_ring_classifications_naphthalene() {
        let mol = naphthalene_kekule();
        let model = assign_aromaticity(&mol);
        let classifications = model.ring_classifications();
        assert_eq!(classifications.len(), 2, "naphthalene has two rings");
        for (_, classification, count) in classifications {
            assert_eq!(*classification, RingAromaticity::Aromatic);
            assert_eq!(*count, 6);
        }
    }

    #[test]
    fn test_non_aromatic_cyclohexane() {
        let mol = cyclohexane();
        let model = assign_aromaticity(&mol);
        for (_, classification, _) in model.ring_classifications() {
            assert_ne!(*classification, RingAromaticity::Aromatic);
            assert_ne!(*classification, RingAromaticity::Antiaromatic);
        }
    }

    // =========================================================================
    // Electron distribution
    // =========================================================================

    #[test]
    fn test_thiophene_aromatic() {
        let mut b = MoleculeBuilder::new();
        let s = b.add_atom(Atom::new(Element::S));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [s, c1, c2, c3, c4];
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        let mol = b.build();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5);
        assert_eq!(model.ring_classifications()[0].2, 6);
    }

    #[test]
    fn test_electron_distribution_tracking() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.ring_classifications()[0].2, 6, "benzene: 6 × 1π = 6");

        let mol = pyrrole_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.ring_classifications()[0].2,
            6,
            "pyrrole: N(2π) + 4C(1π) = 6"
        );

        let mol = furan_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.ring_classifications()[0].2,
            6,
            "furan: O(2π) + 4C(1π) = 6"
        );
    }

    // =========================================================================
    // Aromatic-SMILES input (BondOrder::Aromatic, no kekulization)
    // Verifies that assign_aromaticity works on pre-kekulization molecules.
    // =========================================================================

    #[test]
    fn test_benzene_aromatic_smiles() {
        // c1ccccc1 — parsed with BondOrder::Aromatic bonds
        let mol = mol_aromatic("c1ccccc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "benzene from aromatic SMILES"
        );
    }

    #[test]
    fn test_naphthalene_aromatic_smiles() {
        let mol = mol_aromatic("c1ccc2ccccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "naphthalene from aromatic SMILES"
        );
    }

    #[test]
    fn test_pyridine_aromatic_smiles() {
        let mol = mol_aromatic("c1ccncc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "pyridine from aromatic SMILES"
        );
    }

    #[test]
    fn test_furan_aromatic_smiles() {
        let mol = mol_aromatic("c1ccoc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5, "furan from aromatic SMILES");
    }

    #[test]
    fn test_pyrrole_aromatic_smiles() {
        // [nH] bracket atom: hydrogen_count = Some(1)
        let mol = mol_aromatic("c1cc[nH]c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            5,
            "pyrrole from aromatic SMILES"
        );
    }

    #[test]
    fn test_thiophene_aromatic_smiles() {
        let mol = mol_aromatic("c1ccsc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            5,
            "thiophene from aromatic SMILES"
        );
    }

    // =========================================================================
    // Fused-ring kekulized systems (Pass 2 propagation)
    // =========================================================================

    #[test]
    fn test_indole_aromatic() {
        // c1ccc2[nH]ccc2c1 — indole (9 atoms, 5-ring + 6-ring fused)
        let mol = mol_kekulized("c1ccc2[nH]ccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 indole atoms aromatic"
        );
    }

    #[test]
    fn test_benzimidazole_aromatic() {
        // Two N atoms in fused 5+6 ring system
        let mol = mol_kekulized("c1ccc2[nH]cnc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 9, "all 9 benzimidazole atoms");
    }

    #[test]
    fn test_quinoline_aromatic() {
        let mol = mol_kekulized("c1ccc2ncccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 10, "all 10 quinoline atoms");
    }

    #[test]
    fn test_acridine_aromatic() {
        // 3 fused 6-membered rings, central N: 13 atoms
        let mol = mol_kekulized("c1ccc2nc3ccccc3cc2c1");
        let model = assign_aromaticity(&mol);
        // acridine is C13H9N → 14 heavy atoms (13 C + 1 N), all aromatic
        assert_eq!(model.aromatic_atom_count(), 14, "all 14 acridine atoms");
    }

    // =========================================================================
    // Fused-ring aromatic-SMILES input (BondOrder::Aromatic, kekulize fails)
    // =========================================================================

    #[test]
    fn test_indolizine_aromatic() {
        // c1ccn2cccc2c1 — indolizine: bridgehead N, kekulization unsupported.
        // The SSSR finds a 6-ring and a 9-ring; the 5-ring is recovered via
        // augmentation (XOR of 6- and 9-ring).
        // Pass 1: 5-ring (augmented) detected via bridgehead-N rule → 6π.
        // Pass 2: 6-ring detected using N already aromatic from 5-ring → 6π.
        // The 9-ring (SSSR artifact) is NonAromatic (9π ≠ 4n+2), but all
        // 9 atoms are correctly flagged aromatic via the 5- and 6-ring.
        let mol = mol_aromatic("c1ccn2cccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 indolizine atoms aromatic"
        );
        // At least the 6-ring should be classified as Aromatic in the SSSR set.
        let has_aromatic_ring = model
            .ring_classifications()
            .iter()
            .any(|(_, cls, _)| *cls == RingAromaticity::Aromatic);
        assert!(has_aromatic_ring, "at least one SSSR ring aromatic");
    }

    #[test]
    fn test_purine_aromatic() {
        // c1cnc2[nH]cnc2n1 — purine: 9 atoms, kekulizable
        let mol = mol_kekulized("c1cnc2[nH]cnc2n1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 purine atoms aromatic"
        );
    }

    #[test]
    fn test_purine_aromatic_from_aromatic_smiles() {
        let mol = mol_aromatic("c1cnc2[nH]cnc2n1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "purine from aromatic SMILES"
        );
    }

    #[test]
    fn test_2_pyridinone_aromatic() {
        // O=c1ccncc1 — 2-pyridinone (aromatic SMILES, N without H, exo C=O).
        // Kekulization fails; tested on the aromatic-bond form directly.
        // The exo C=O gives the C atom has_double_any=true → 1π.
        // N has Aromatic bonds in ring → 1π (pyridine-like).
        // Total: 6 × 1π = 6π → aromatic.
        let mol = mol_aromatic("O=c1ccncc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "all 6 ring atoms of 2-pyridinone aromatic"
        );
    }

    #[test]
    fn test_quinolone_aromatic() {
        // O=c1ccc2ncccc2c1 — quinolone: fused 6+6 with exo C=O, kekulize fails
        let mol = mol_aromatic("O=c1ccc2ncccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "all 10 quinolone ring atoms aromatic"
        );
        assert_eq!(
            model.ring_classifications().len(),
            2,
            "two rings classified"
        );
    }

    #[test]
    fn test_indole_aromatic_smiles() {
        let mol = mol_aromatic("c1ccc2[nH]ccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "indole from aromatic SMILES"
        );
    }

    // =========================================================================
    // Bridgehead N rule: specifically test that the rule fires correctly
    // =========================================================================

    #[test]
    fn test_bridgehead_n_contributes_lone_pair() {
        // Indolizine: the bridgehead N (degree 3, no H, no explicit double bond)
        // must be detected as a 2π contributor for the 5-membered ring.
        // We verify by checking the 5-ring classification (if accessible).
        let mol = mol_aromatic("c1ccn2cccc2c1");
        let model = assign_aromaticity(&mol);
        // All 9 atoms aromatic: both rings must be aromatic.
        assert_eq!(model.aromatic_atom_count(), 9);
        // The bridgehead N itself must be in the aromatic set.
        // In the SMILES c1ccn2cccc2c1, n is atom index 3.
        assert!(
            model.is_atom_aromatic(AtomIdx(3)),
            "bridgehead N must be aromatic"
        );
    }

    #[test]
    fn test_non_bridgehead_n_no_false_positive() {
        // Pyrimidine: two N atoms in a 6-membered ring, no bridgehead.
        // Both N have ring_degree == total_degree == 2.
        // Should be detected as aromatic via has_aromatic_in_ring (Aromatic bonds).
        let mol = mol_aromatic("c1ccncn1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6, "pyrimidine is aromatic");
    }

    #[test]
    fn test_imidazole_aromatic() {
        // c1cn[nH]c1 / c1c[nH]cn1 — imidazole: one pyridine-type N, one pyrrole-type N
        let mol = mol_aromatic("c1cn[nH]c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5, "imidazole is aromatic");
    }

    // =========================================================================
    // Pass 2 specifically: rings that need fused-ring context
    // =========================================================================

    #[test]
    fn test_pass2_needed_for_indolizine_6ring() {
        // The augmented 5-ring (XOR of SSSR 6-ring and 9-ring) is detected aromatic in Pass 1.
        // The SSSR 6-ring is then detected aromatic in Pass 2 (N already aromatic → 1π).
        // The SSSR 9-ring (9π) remains NonAromatic per Hückel.
        // Key assertion: all 9 atoms are aromatic (correct overall perception).
        let mol = mol_aromatic("c1ccn2cccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 indolizine atoms aromatic"
        );
        // The bridgehead N must be aromatic.
        assert!(
            model.is_atom_aromatic(AtomIdx(3)),
            "bridgehead N is aromatic"
        );
        // The 6-ring (SSSR ring, improved by Pass 2) should be classified Aromatic.
        let aromatic_count = model
            .ring_classifications()
            .iter()
            .filter(|(_, cls, _)| *cls == RingAromaticity::Aromatic)
            .count();
        assert!(aromatic_count >= 1, "at least one SSSR ring is aromatic");
    }

    #[test]
    fn test_no_pass2_needed_for_naphthalene() {
        // Naphthalene: both rings pass independently in Pass 1.
        // Verifies Pass 2 doesn't break things that already work.
        let mol = naphthalene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 10);
        let classes = model.ring_classifications();
        assert_eq!(classes.len(), 2);
        for (_, cls, _) in classes {
            assert_eq!(*cls, RingAromaticity::Aromatic);
        }
    }

    #[test]
    fn test_anthracene_aromatic() {
        // c1ccc2cc3ccccc3cc2c1 — anthracene: 3 linearly fused 6-rings, 14 atoms
        let mol = mol_kekulized("c1ccc2cc3ccccc3cc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 14, "all 14 anthracene atoms");
    }

    // =========================================================================
    // Regression: aromatic-bond path must not perturb kekulized correctness
    // =========================================================================

    #[test]
    fn test_kekulized_path_unaffected_by_aromatic_bond_changes() {
        // Kekulized benzene: bonds are Double/Single, not Aromatic.
        // The new Aromatic-bond branches must stay dormant.
        let mol = benzene_kekule();
        // Verify no aromatic bonds in input.
        for (_, bond) in mol.bonds() {
            assert_ne!(bond.order, BondOrder::Aromatic, "input must be kekulized");
        }
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6);
        // All 6 bonds in benzene ring should be aromatic.
        let aromatic_bonds = mol
            .bonds()
            .filter(|(b, _)| model.is_bond_aromatic(*b))
            .count();
        assert_eq!(aromatic_bonds, 6);
    }

    #[test]
    fn test_keto_pyridinone_aromatic() {
        // O=C1NC=CC=C1 — 2-pyridinone keto form with N-H.
        // π count: C(=O)(0π, exocyclic-only double bond to O) + N-H(2π) +
        // 4×C in 2 ring C=C (1π each) = 6π → aromatic. Matches RDKit, which
        // marks all 6 ring atoms aromatic (exocyclic O stays non-aromatic).
        let mol = mol_kekulized("O=C1NC=CC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "keto pyridinone ring is Hückel aromatic (6π = 4n+2)"
        );
    }

    #[test]
    fn test_tropone_aromatic() {
        // O=C1C=CC=CC=C1 — tropone (cycloheptatrienone), Kekulized input.
        // Carbonyl C contributes 0π (exocyclic-only double bond to O); the
        // other 6 ring carbons contribute 1π each from 3 endocyclic C=C.
        // Total 6π → aromatic, matching RDKit (all 7 ring atoms aromatic).
        let mol = mol_kekulized("O=C1C=CC=CC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            7,
            "all 7 tropone ring atoms aromatic"
        );
    }

    #[test]
    fn test_4_pyridone_aromatic() {
        // O=C1C=CNC=C1 — 4-pyridone, Kekulized input. Same 6π accounting as
        // 2-pyridone, just with N para to the carbonyl. Matches RDKit.
        let mol = mol_kekulized("O=C1C=CNC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "all 6 4-pyridone ring atoms aromatic"
        );
    }

    #[test]
    fn test_pyranone_aromatic() {
        // O=C1C=COC=C1 — 4H-pyran-4-one, Kekulized input. Ring O contributes
        // 2π (lone pair), carbonyl C contributes 0π, remaining 4 ring carbons
        // contribute 1π each from 2 endocyclic C=C. Total 6π. Matches RDKit.
        let mol = mol_kekulized("O=C1C=COC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "all 6 pyranone ring atoms aromatic"
        );
    }

    #[test]
    fn test_cyclopentadienyl_anion_aromatic() {
        // [CH-]1C=CC=C1 — cyclopentadienyl anion. The carbanion carbon has no
        // double bond but contributes 2π (lone pair); the other 4 carbons
        // contribute 1π each from 2 endocyclic C=C. Total 6π. Matches RDKit
        // (all 5 atoms aromatic).
        let mol = mol_kekulized("[CH-]1C=CC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            5,
            "all 5 cyclopentadienyl anion atoms aromatic"
        );
    }

    // ── N-substituted pyrrole-type N: bridgehead-branch guard removal ────────
    //
    // The bridgehead-N branch used to require the exocyclic substituent to be
    // sp2, to defensively block imide N (phthalimide). That guard also
    // blocked the much more common case of a plain alkyl/aryl/sugar
    // substituent on an otherwise-aromatic pyrrole-type N. It was removed;
    // these tests cover both the newly-fixed cases and the phthalimide
    // regression it was guarding against (which stays correct via the
    // overall 4n+2 sum, not the substituent).

    #[test]
    fn test_n_methylpyrrole_aromatic() {
        let mol = mol_kekulized("CN1C=CC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            5,
            "all 5 N-methylpyrrole ring atoms aromatic"
        );
    }

    #[test]
    fn test_n_methylimidazole_aromatic() {
        let mol = mol_kekulized("CN1C=CN=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            5,
            "all 5 N-methylimidazole ring atoms aromatic"
        );
    }

    #[test]
    fn test_n_methylindole_aromatic() {
        let mol = mol_kekulized("CN1C=CC2=CC=CC=C21");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 N-methylindole ring atoms aromatic"
        );
    }

    #[test]
    fn test_9_methylpurine_aromatic() {
        let mol = mol_kekulized("CN1C=NC2=NC=NC=C21");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 9-methylpurine ring atoms aromatic"
        );
    }

    #[test]
    fn test_phthalimide_5ring_not_aromatic() {
        // O=C1NC(=O)c2ccccc21 — only the fused benzo ring is aromatic (6
        // atoms); the imide 5-ring (2 carbonyl C + N) is not: carbonyl
        // carbons contribute 0π each (exocyclic C=O rule), N contributes 2π,
        // the two ring-fusion carbons contribute 1π each — 4π total, not
        // 4n+2. Regression guard for the bridgehead-N guard removal above.
        let mol = mol_kekulized("O=C1NC(=O)c2ccccc21");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "only the 6 benzo atoms of phthalimide are aromatic"
        );
    }

    #[test]
    fn test_n_methylphthalimide_5ring_not_aromatic() {
        // O=C1N(C)C(=O)c2ccccc21 — same as phthalimide but N-methylated;
        // same accounting applies (N still contributes 2π regardless of
        // substituent), 5-ring still non-aromatic.
        let mol = mol_kekulized("O=C1N(C)C(=O)c2ccccc21");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "only the 6 benzo atoms of N-methylphthalimide are aromatic"
        );
    }

    #[test]
    fn test_azulene_kekulized_aromatic() {
        // C1=CC2=CC=CC=CC2=C1 — non-alternant fused bicyclic, all 10 atoms
        // aromatic per RDKit. Regression coverage: this was previously
        // (incorrectly) believed to need a ring-system rewrite, based on a
        // test that never called apply_aromaticity() on Kekulized input.
        let mol = mol_kekulized("C1=CC2=CC=CC=CC2=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "all 10 azulene atoms aromatic"
        );
    }

    // ── RDKit #9271: charged / zwitterionic aromatic systems ─────────────────

    #[test]
    fn test_fluorescein_dianion_aromatic() {
        // Fluorescein dianion: RDKit #9271 incorrectly marked xanthene bonds as
        // single instead of aromatic. Verify chematic parses and identifies
        // aromatic atoms correctly (two benzene rings + xanthene O-bridge ring).
        // Kekulé-form SMILES: all atoms uppercase.
        let smi = "C1=CC=C(C(=C1)C2=C3C=CC(=O)C=C3OC4=C2C=CC(=C4)[O-])C(=O)[O-]";
        let mol = chematic_smiles::parse(smi).expect("fluorescein dianion should parse");
        // The molecule should parse without panic. Verify aromatic ring count:
        // fluorescein has 3 aromatic rings (2 benzene + xanthene core).
        let arc = count_aromatic_rings(&mol);
        assert!(
            arc >= 2,
            "fluorescein dianion: expected ≥2 aromatic rings, got {arc} \
             (RDKit #9271: charged aromatics may be misclassified)"
        );
    }

    #[test]
    fn test_rhodamine_zwitterion_parses() {
        // Rhodamine-type zwitterion with N+ and bridging O (RDKit #9271).
        // Must parse cleanly and produce a valid aromatic ring count.
        let smi = "CCN(CC)c1ccc2c(-c3ccccc3C(=O)O)c3ccc(=[N+](CC)CC)cc-3oc2c1";
        let mol = chematic_smiles::parse(smi).expect("rhodamine zwitterion should parse");
        let arc = count_aromatic_rings(&mol);
        assert!(arc >= 3, "rhodamine: expected ≥3 aromatic rings, got {arc}");
    }

    #[test]
    fn test_cyclopentadienyl_not_aromatic_kekulized() {
        // C1=CC=CC1 — cyclopentadiene (4 C with doubles + 1 sp3 CH2): not aromatic.
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C)); // sp3
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Double).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        b.add_bond(c3, c4, BondOrder::Double).unwrap();
        b.add_bond(c4, c0, BondOrder::Single).unwrap();
        let mol = b.build();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            0,
            "cyclopentadiene not aromatic"
        );
    }

    // =========================================================================
    // RdkitLike mode: Se/Te chalcogen heteroaromatics
    // =========================================================================

    #[test]
    fn test_selenophene_huckel_not_aromatic() {
        // c1cc[se]c1 — in strict Hückel mode, Se is unsupported → 0 aromatic atoms
        // (assign_aromaticity_ex re-derives from scratch, ignoring parser's aromatic flags)
        let mol = mol_aromatic("c1cc[se]c1");
        let m = assign_aromaticity(&mol); // default Hückel
        assert_eq!(
            m.aromatic_atom_count(),
            0,
            "selenophene: Se not aromatic in Hückel mode"
        );
    }

    #[test]
    fn test_selenophene_rdkit_aromatic() {
        // c1cc[se]c1 — in RdkitLike mode, Se donates 2π → 6π total → aromatic
        let mol = mol_aromatic("c1cc[se]c1");
        let m = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m.aromatic_atom_count(),
            5,
            "selenophene: all 5 atoms aromatic in RdkitLike"
        );
    }

    #[test]
    fn test_tellurophene_rdkit_aromatic() {
        // c1cc[te]c1 — Te analogous to Se (2π donor)
        let mol = mol_aromatic("c1cc[te]c1");
        let m = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m.aromatic_atom_count(),
            5,
            "tellurophene: all 5 atoms aromatic in RdkitLike"
        );
    }

    #[test]
    fn test_benzoselenophene_rdkit() {
        // Fused benzene + selenophene
        let mol = mol_aromatic("c1ccc2[se]ccc2c1");
        let m = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m.aromatic_atom_count(),
            9,
            "benzoselenophene: 9 atoms aromatic"
        );
    }

    #[test]
    fn test_rdkit_mode_does_not_break_benzene() {
        // Benzene must give same result in both modes
        let mol = mol_aromatic("c1ccccc1");
        let m_h = assign_aromaticity(&mol);
        let m_r = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(m_h.aromatic_atom_count(), m_r.aromatic_atom_count());
    }

    #[test]
    fn test_rdkit_mode_does_not_break_thiophene() {
        let mol = mol_aromatic("c1ccsc1");
        let m_h = assign_aromaticity(&mol);
        let m_r = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m_h.aromatic_atom_count(),
            m_r.aromatic_atom_count(),
            "thiophene same in both modes"
        );
    }

    // ── Known regressions from fix #2 (bridgehead-N guard removal) ──────────
    //
    // PROVISIONAL: measured pre-SSSR-fix. find_sssr is itself non-deterministic
    // and non-minimal for ~50% of sampled molecules (single spanning-tree
    // fundamental-cycle basis with zero candidate redundancy -- see project
    // plan). These counts are only valid for the current find_sssr; expect
    // them to change once ring perception is fixed, independent of anything
    // below being "resolved."
    //
    // These 32 molecules share one root cause: a "fake bridgehead" N (same
    // local shape as a genuine bridgehead or N-substituted azole) feeds a
    // central ring that only closes via the `aromatic_context` bypass reusing
    // an unrelated ring's atoms. Fixing this requires removing the bypass in
    // favor of proper ring-system candidate enumeration (see project plan/
    // issue tracker). Pinned here as *known-wrong* so the eventual fix is
    // measurable by how many of these flip from this assertion to correct,
    // not just by an aggregate corpus percentage.
    #[test]
    fn test_known_regressions_from_bridgehead_n_fix() {
        // (kekulized SMILES, current chematic aromatic_atom_count(), RDKit's correct count)
        let cases: &[(&str, usize, usize)] = &[
            ("C[Si](C)(C)C1=CC=C(C2=CC3=CC=CC=C3C3=NCCCN23)C=C1", 16, 12),
            (
                "C1=C(C2=CC=C(CCC3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
                22,
                18,
            ),
            ("ClC1=CC=C(OCC2=CC3=CC=CC=C3C3=NCCCN23)C=C1", 16, 12),
            ("N[C@@H](CC1=CC=CC=C1)C1=CC2=CC=CC=C2C2=NCCCN12", 16, 12),
            (
                "CC(C)(C)C1=CC=C(C2=C(CC3=CC=CC=C3)C3=CC=CC=C3C3=NCCCN32)C=C1",
                22,
                18,
            ),
            (
                "C[Si](C)(C)C1=CC=C(C2=C(CC3=CC=CC=C3)C3=CC=CC=C3C3=NCCCN32)C=C1",
                22,
                18,
            ),
            (
                "C1=C(C2=CC=C(C3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
                22,
                18,
            ),
            (
                "C1=C(C2=CC=C(OCC3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
                22,
                18,
            ),
            ("COC1=C(OC)C(OC)=CC(C2=CC3=CC=CC=C3C3=NCCCN23)=C1", 16, 12),
            ("CC1=CC2=CC=CC=C2C2=NCCCN12", 10, 6),
            (
                "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(NC(=O)NC4CCCCC4)C=C3)C3=NCCCN23)C=C1",
                16,
                12,
            ),
            (
                "C1=CC=C(CCC2=CC=C(C3=C(CC4=CC=CC=C4)C4=CC=CC=C4C4=NCCCN43)C=C2)C=C1",
                28,
                24,
            ),
            (
                "CCCCC1=C(C2=CC=C(CCC3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
                22,
                18,
            ),
            (
                "CCCCC1=C(C2=CC=C(C(C)(C)C)C=C2)N2CCCN=C2C2=CC=CC=C12",
                16,
                12,
            ),
            ("CCCCCCC1=CC2=CC=CC=C2C2=NCCCN12", 10, 6),
            (
                "CCOC1=CC=C(CC2=C(CCCC3=CC=CC4=CC=CC=C34)N3CCCN=C3C3=CC=CC=C23)C=C1",
                26,
                22,
            ),
            (
                "CCOC1=CC=C(CC2=C(C3=CC=C(CCC4=CC=CC=C4)C=C3)N3CCCN=C3C3=CC=CC=C23)C=C1",
                28,
                24,
            ),
            (
                "CN(C)CCC1=C(C2=CC=C(C(C)(C)C)C=C2)N2CCCN=C2C2=CC=CC=C12",
                16,
                12,
            ),
            (
                "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(N/C(S)=N/C4CCCCC4)C=C3)C3=NCCCN23)C=C1",
                16,
                12,
            ),
            ("C1=C(/C=C/C2=CC=CC=C2)N2CCCN=C2C2=CC=CC=C12", 16, 12),
            ("CC(C)(C)C1=CC=C(C2=CC3=CC=CC=C3C3=NCCCN23)C=C1", 16, 12),
            (
                "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(NC(=O)CC4=CC=CC=N4)C=C3)C3=NCCCN23)C=C1",
                22,
                18,
            ),
            (
                "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(NC(=O)NC4=C(Cl)C=C(Cl)C=C4)C=C3)C3=NCCCN23)C=C1",
                22,
                18,
            ),
            ("C1=C(CC2=CC=CC=C2)C2=CC=CC=C2C2=NCCCN12", 16, 12),
            ("ClC1=CC=C(C2=CC3=CC=CC=C3C3=NCCCN23)C=C1", 16, 12),
            ("C1=C(C2=CC=CC=C2)N2CCCN=C2C2=CC=CC=C12", 16, 12),
            (
                "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(N(CC4=CC=CC=C4)CC4=CC=CC=C4)C=C3)C3=NCCCN23)C=C1",
                28,
                24,
            ),
            (
                "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(N)C=C3)C3=NCCCN23)C=C1",
                16,
                12,
            ),
            ("CC1=C2C(=NC=C1)N(C1CC1)C1=NC=CC=C1C(=O)N2C", 15, 12),
            ("CC(=O)N1C2=NC=CC=C2C(=O)N(C)C2=CC=CN=C21", 15, 12),
            ("CN1C(=O)C2=CC=CN=C2N(C(C)(C)C)C2=NC=CC=C21", 15, 12),
            ("CCCN1C2=NC=CC=C2C(=O)N(C)C2=CC=CN=C21", 15, 12),
        ];
        for (smi, expected_wrong, rdkit_correct) in cases {
            let mol = mol_kekulized(smi);
            let model = assign_aromaticity(&mol);
            assert_eq!(
                model.aromatic_atom_count(),
                *expected_wrong,
                "{smi}: expected current (wrong) count {expected_wrong} (RDKit correct: {rdkit_correct})"
            );
        }
    }

    // ── Known order-dependence: same molecule, different Kekulized traversal ─
    //
    // PROVISIONAL: measured pre-SSSR-fix (see note above).
    //
    // These 3 molecules pass with RDKit's canonical Kekulized SMILES but fail
    // with at least one other valid Kekulized ordering of the identical
    // structure -- confirmed via atom-map-number alignment (no substructure
    // matching). Root cause is NOT Pass 1/Pass 2 (verified order-invariant by
    // construction: Pass 1 evaluates each ring independently, Pass 2 loops to
    // a fixed point) -- it's `find_sssr` itself, which builds a single BFS
    // spanning tree rooted at the lowest atom index and generates exactly one
    // fundamental cycle per non-tree edge (zero candidate redundancy), so it
    // can return a non-minimal ring (e.g. a 10-membered ring in place of two
    // 6-membered ones) depending on traversal order. Each entry below is one
    // concrete failing traversal, pinned as a fixed regression.
    #[test]
    fn test_known_order_dependent_regressions() {
        let cases: &[(&str, usize, usize)] = &[
            (
                "N1=C2C(N(CC(O)=O)C(=O)N=C2N(C2C=C(C(F)(F)F)C=C(C=2)C(F)(F)F)C2C1=CC=CC=2)=O",
                16,
                20,
            ),
            (
                "[C@H]12N(C([C@H](NC(=O)[C@H]([C@H](OC(=O)[C@@H](N(C)C(CN(C)C1=O)=O)C(C)C)C)NC(=O)C1C=C(OC)C(C)=C3OC4=C(C)C(=O)C(=C(C4=NC=13)C(=O)N[C@H]1C(=O)N[C@@H](C(C)C)C(N3[C@H](C(=O)N(CC(N([C@H](C(C)C)C(O[C@H]1C)=O)C)=O)C)CCC3)=O)N)C(C)C)=O)CCC2",
                6,
                14,
            ),
            ("C12N(C3C=CC=CC=3)C3=NC(=O)N(C)C(C3=NC1=CC=CC=2)=O", 16, 20),
        ];
        for (smi, expected_wrong, rdkit_correct) in cases {
            let mol = mol_kekulized(smi);
            let model = assign_aromaticity(&mol);
            assert_eq!(
                model.aromatic_atom_count(),
                *expected_wrong,
                "{smi}: expected current (wrong) count {expected_wrong} (RDKit correct: {rdkit_correct})"
            );
        }
    }
}
