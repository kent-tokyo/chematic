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
    /// RDKit-compatible extension. Adds P (15), Se (34), and Te (52) as
    /// heteroatom lone-pair donors (2π), matching the RDKit DEFAULT
    /// aromaticity model for common organic heteroaromatics.
    ///
    /// Keto-lactam aromaticity is NOT included (TautomerMode, separate sprint).
    RdkitLike,
}

use rustc_hash::{FxHashMap, FxHashSet};

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, implicit_hcount};

use crate::ring_family::RingFamily;
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
/// The default model also has a deliberately narrow whole-envelope fallback
/// for all-carbon odd/odd fused systems such as azulene.
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

    /// Build a model directly from an aromatic atom/bond set, with no ring
    /// classification or antiaromaticity data.
    ///
    /// Used by engines (e.g. `rdkit_parity`'s experimental production API)
    /// that determine an aromatic atom/bond set directly rather than via
    /// this module's own per-ring Hückel passes -- `ring_classifications()`
    /// and `antiaromatic_rings()` are empty on the result.
    pub(crate) fn from_atom_bond_sets(
        aromatic_atoms: FxHashSet<AtomIdx>,
        aromatic_bonds: FxHashSet<BondIdx>,
    ) -> Self {
        AromaticityModel {
            aromatic_atoms,
            aromatic_bonds,
            antiaromatic_rings: Vec::new(),
            ring_classifications: Vec::new(),
        }
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
/// Pass [`AromaticityAlgorithm::RdkitLike`] to additionally recognise P/Se/Te
/// as lone-pair donors in aromatic rings.
///
/// Byte-identical to this function's behavior before the K2b
/// authoritative-demotion work started (`ring_pi_electrons`'s carbon rule
/// does not get the ring-fusion-aware fix here -- see
/// [`assign_aromaticity_authoritative_experimental`] for the opt-in variant
/// that does).
pub fn assign_aromaticity_ex(mol: &Molecule, algo: AromaticityAlgorithm) -> AromaticityModel {
    assign_aromaticity_ex_impl(mol, algo, false)
}

/// Opt-in variant of [`assign_aromaticity_ex`] with the K2b fused-diazine
/// ring-fusion fix enabled in `ring_pi_electrons`'s carbon rule (see its doc
/// comment) -- a ring-fusion bond into an adjacent ring's heteroatom is no
/// longer wrongly treated as a genuine exocyclic substituent. Always uses
/// [`AromaticityAlgorithm::Huckel`], matching [`assign_aromaticity`]'s own
/// default (this mechanism is orthogonal to the `RdkitLike` Se/Te
/// extension; the ordinary `RdkitLike` path now uses the verified fused-ring
/// parity engine when its pre-kekulized-input precondition can be met).
///
/// **Known limitation, honestly documented, not a blocker to using this**:
/// resolves 29/33 of the corpus cluster this fix targets
/// (`fused_diazine_quinazoline_quinoxaline_purine`, see
/// `validation/results/aromaticity_flag_demotion_k2b_fused_diazine_fix_summary.json`)
/// but does NOT fix two other, architecturally distinct, still-open gaps in
/// the underlying per-ring Pass 1/Pass 2 Hückel model: non-alternant
/// whole-perimeter systems like azulene (49 corpus molecules; see
/// `validation/results/aromaticity_flag_demotion_k2b_azulene_cluster_finding.json`
/// for why this is not boundable by a rule-level fix) and 2 large fused
/// polycyclic cage molecules with a similar odd-π-count blind spot (plus 4
/// molecules that combine both the now-fixed and the still-open mechanism in
/// the same molecule). Real, verified improvement over the promote-only
/// default nonetheless -- see `test_authoritative_experimental_*` below.
pub fn assign_aromaticity_authoritative_experimental(mol: &Molecule) -> AromaticityModel {
    assign_aromaticity_ex_impl(mol, AromaticityAlgorithm::Huckel, true)
}

fn assign_aromaticity_ex_impl(
    mol: &Molecule,
    algo: AromaticityAlgorithm,
    ring_fusion_aware: bool,
) -> AromaticityModel {
    // The RDKit-compatible mode uses the independently verified parity engine
    // as its production path. Unlike this module's historical per-ring Hückel
    // pass, that engine evaluates connected fused-ring subsets and therefore
    // handles non-alternant whole-perimeter systems such as azulene. Keep the
    // old infallible implementation as a defensive fallback for molecules the
    // parity engine cannot kekulize; callers needing to distinguish that case
    // can use the fallible `assign_aromaticity_rdkit_parity_experimental` API.
    if algo == AromaticityAlgorithm::RdkitLike
        && let Ok(model) = crate::rdkit_parity::assign_aromaticity_rdkit_parity_experimental(mol)
    {
        return model;
    }

    let ring_set = find_sssr(mol);
    let sssr_rings = ring_set.rings();

    // Augment SSSR rings with smaller XOR sub-rings (GF(2) differences between pairs).
    // This corrects the case where the SSSR algorithm stores a large fundamental cycle
    // instead of its smaller GF(2)-reduced equivalent (e.g. the 5-ring of indolizine).
    let rings: Vec<Vec<AtomIdx>> = augmented_ring_set(mol, sssr_rings);

    // K2b fused-diazine fix, opt-in only (`ring_fusion_aware`): the
    // whole-molecule set of bonds that lie on ANY ring (not just the one
    // ring currently being evaluated). Computed once here (cheap:
    // proportional to total ring length, reusing the existing
    // `ring_bond_set` helper) and threaded into `ring_pi_electrons` so its
    // carbon "genuine exocyclic double bond" rule can tell a real substituent
    // (tropone's C=O, whose far atom is on no ring at all) apart from a
    // ring-fusion bond whose far atom just happens to lie in a DIFFERENT ring
    // than the one under evaluation (see `ring_pi_electrons`'s doc comment).
    // Deliberately not recomputed per-atom inside the hot loop -- an O(V+E)
    // ring-bond check per query there previously caused a real 10-14x perf
    // regression (SSSR misused as a boolean ring-bond check); this set is the
    // same for every ring in this call, so it is built exactly once.
    //
    // `assign_aromaticity_ex` (the default, byte-identical-to-pre-K2b entry
    // point) passes `ring_fusion_aware = false` here, which keeps this set
    // EMPTY -- `ring_pi_electrons`'s `!all_ring_bonds.contains(&bidx)` check
    // is then unconditionally true, exactly reproducing the pre-fix
    // `!ring_atom_set.contains(&nb)` check it replaced (that check was
    // itself already guaranteed true by this point: the preceding sibling
    // condition already established no Double-bonded neighbor is in
    // `ring_atom_set`, so a real neighbor reaching this check was never in
    // it either way). Verified byte-identical against `main` pre-K2b via the
    // full 5000-molecule corpus (both calling conventions), not just
    // reasoned about -- see the authoritative-experimental test module.
    let all_ring_bonds: FxHashSet<BondIdx> = if ring_fusion_aware {
        rings.iter().flat_map(|r| ring_bond_set(mol, r)).collect()
    } else {
        FxHashSet::default()
    };

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
        match ring_pi_electrons(mol, ring, &empty_context, algo, &all_ring_bonds) {
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
            match ring_pi_electrons(mol, ring, &aromatic_atoms, algo, &all_ring_bonds) {
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
        // Once every atom in the candidate ring set is already aromatic, no
        // pending ring can add information to the aromatic context. This is
        // RDKit's `aromRingsAllSet` fixed-point short circuit; in particular,
        // it prevents a later indeterminate ring from reopening a converged
        // fused-ring component.
        let arom_rings_all_set = rings
            .iter()
            .flatten()
            .all(|atom| aromatic_atoms.contains(atom));
        if !any_new || arom_rings_all_set {
            break;
        }
    }

    // A strict per-ring pass cannot seed azulene's 5+7 fused system because
    // both constituent rings have an odd local count. Apply only the narrow,
    // fail-closed all-carbon odd/odd envelope rule here; this does not route
    // the default model through the broader RdkitLike implementation.
    if algo == AromaticityAlgorithm::Huckel {
        apply_huckel_nonalternant_fused_fallback(
            mol,
            &rings,
            &mut aromatic_atoms,
            &mut aromatic_bonds,
        );
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

fn apply_huckel_nonalternant_fused_fallback(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    aromatic_atoms: &mut FxHashSet<AtomIdx>,
    aromatic_bonds: &mut FxHashSet<BondIdx>,
) {
    let families = crate::ring_family::find_ring_families_over(mol, rings);
    let all_ring_bonds: FxHashSet<BondIdx> = rings
        .iter()
        .flat_map(|ring| ring_bond_set(mol, ring))
        .collect();

    for candidate in build_conjugated_components(
        mol,
        rings,
        &families,
        AromaticityAlgorithm::Huckel,
        &all_ring_bonds,
    ) {
        if candidate.source_rings.len() < 2
            || candidate
                .source_rings
                .iter()
                .any(|&ring_idx| rings[ring_idx].len().is_multiple_of(2))
            || !candidate.atoms.len().wrapping_sub(2).is_multiple_of(4)
            || candidate
                .atoms
                .iter()
                .any(|&atom_idx| mol.atom(atom_idx).element.atomic_number() != 6)
        {
            continue;
        }

        // Every eligible carbon contributes one electron in this narrow
        // envelope. The size check above is therefore the 4n+2 test.
        for &atom_idx in &candidate.atoms {
            aromatic_atoms.insert(atom_idx);
        }
        for &atom_idx in &candidate.atoms {
            for (neighbor, bond_idx) in mol.neighbors(atom_idx) {
                if candidate.atoms.contains(&neighbor)
                    && matches!(
                        mol.bond(bond_idx).order,
                        BondOrder::Double | BondOrder::Aromatic
                    )
                {
                    aromatic_bonds.insert(bond_idx);
                }
            }
        }
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
///
/// Byte-identical to this function's behavior before the K2b
/// authoritative-demotion work started -- promote-only, matching `main`
/// pre-K2b (see [`build_molecule_from_model`]'s doc comment). See
/// [`apply_aromaticity_authoritative_experimental`] for the opt-in variant.
pub fn apply_aromaticity_ex(mol: &Molecule, algo: AromaticityAlgorithm) -> Molecule {
    let model = assign_aromaticity_ex(mol, algo);
    build_molecule_from_model(mol, &model)
}

/// Apply aromaticity using the opt-in, authoritative-demotion engine (see
/// [`assign_aromaticity_authoritative_experimental`] for the mechanism and
/// its documented, still-open limitations).
///
/// Returns a new [`Molecule`] where an atom's aromatic flag reflects the
/// model's verdict in BOTH directions -- promoted when the model confirms
/// it, DEMOTED when a stale parser-set `aromatic: true` the model does not
/// independently confirm survived from the input. [`apply_aromaticity_ex`]
/// (the default) only ever promotes.
///
/// Explicitly opt-in and separate from [`apply_aromaticity`]/
/// [`apply_aromaticity_ex`] -- those remain unchanged, matching this
/// codebase's existing pattern for `_experimental` production surfaces (see
/// `apply_aromaticity_rdkit_parity_experimental`). Infallible: unlike the
/// `rdkit_parity` engine, this one does not perform its own internal
/// kekulization, so it has no failure mode `apply_aromaticity_ex` doesn't
/// already have.
pub fn apply_aromaticity_authoritative_experimental(mol: &Molecule) -> Molecule {
    let model = assign_aromaticity_authoritative_experimental(mol);
    build_molecule_from_model_authoritative(mol, &model)
}

/// Build a new [`Molecule`] from `mol` with atom/bond aromaticity flags set
/// according to an already-computed `model`, using the model's verdict to
/// only ever PROMOTE an atom to aromatic, never demote a stale parser-set
/// `aromatic: true` the model doesn't independently confirm.
///
/// This is the original, pre-K2b behavior -- unchanged since before the
/// authoritative-demotion work started, and what [`apply_aromaticity_ex`]
/// (the default entry point) still uses. Bond orders ARE fully authoritative
/// (a bond's order always reflects the model's verdict; there is no
/// "promote-only" ambiguity for bonds, since `bond.order` is unconditionally
/// either the model's `Aromatic` or its own already-Kekulized value) -- only
/// the ATOM flag is promote-only. See [`build_molecule_from_model_authoritative`]
/// for the opt-in, fully bidirectional variant
/// ([`apply_aromaticity_authoritative_experimental`]) that also demotes atom
/// flags, backing `apply_aromaticity_rdkit_parity_experimental` too (a no-op
/// distinction for that caller, since its input is always freshly
/// re-Kekulized with every atom's `aromatic` flag already reset to `false`
/// beforehand -- there is nothing to demote FROM).
pub(crate) fn build_molecule_from_model(mol: &Molecule, model: &AromaticityModel) -> Molecule {
    let bond_orders = compute_bond_orders(mol, model);
    // Promote-only: an atom ends up aromatic if the model confirms it OR it
    // was ALREADY aromatic on `mol` to begin with (`atom.aromatic`) --
    // never demoted. This is NOT the same as "assign from the model's set
    // alone" (that would be a silent demotion of every atom the model
    // doesn't confirm, which is exactly the authoritative variant's job,
    // not this one's) -- the `|| atom.aromatic` term is what makes this
    // function promote-only rather than fully authoritative.
    let atom_aromatic: FxHashSet<AtomIdx> = mol
        .atoms()
        .filter_map(|(idx, atom)| (model.is_atom_aromatic(idx) || atom.aromatic).then_some(idx))
        .collect();
    finish_molecule_with_flags(mol, &atom_aromatic, &bond_orders)
}

/// Authoritative variant of [`build_molecule_from_model`]: the model is
/// authoritative in BOTH directions -- promote AND demote -- instead of only
/// ever promoting (see docs/rfcs/aromaticity_rdkit_parity_rfc.md section 1b/6). A
/// stale parser-set `aromatic: true` the model does not independently
/// confirm does not survive.
///
/// The one deliberate exception: an atom incident to a bond that ends up
/// `Aromatic` in `bond_orders` is always kept aromatic too, even if the
/// model itself didn't confirm it. This is not a reintroduction of the
/// promote-only bug -- it only ever fires when `bond.order` was itself
/// still `Aromatic` going in and the model gave no verdict to demote it
/// with. There is no independently-computed Kekule value to fall back to in
/// that case, so leaving both the atom and its bond flagged aromatic
/// together is the "clean, well-defined fallback state" for that molecule.
/// This can never mask a genuine demotion: for every already-Kekulized
/// input, `bond.order` is a real Single/Double value and this fallback
/// never triggers.
///
/// Backs [`apply_aromaticity_authoritative_experimental`] (opt-in, general
/// mechanism including the fused-diazine ring-fusion fix -- see
/// `assign_aromaticity_authoritative_experimental`) and
/// `apply_aromaticity_rdkit_parity_experimental` (already relies on this
/// behavior; a no-op distinction for it, since its input molecule is always
/// a fresh re-Kekulized clone with every atom's `aromatic` flag reset to
/// `false` first -- there is no stale flag to demote).
pub(crate) fn build_molecule_from_model_authoritative(
    mol: &Molecule,
    model: &AromaticityModel,
) -> Molecule {
    let bond_orders = compute_bond_orders(mol, model);
    let mut atom_aromatic: FxHashSet<AtomIdx> = mol
        .atoms()
        .filter_map(|(idx, _)| model.is_atom_aromatic(idx).then_some(idx))
        .collect();
    for (bidx, bond) in mol.bonds() {
        if bond_orders[&bidx] == BondOrder::Aromatic {
            atom_aromatic.insert(bond.atom1);
            atom_aromatic.insert(bond.atom2);
        }
    }
    finish_molecule_with_flags(mol, &atom_aromatic, &bond_orders)
}

/// The model's per-bond verdict: `Aromatic` when the model confirms it,
/// `bond.order` otherwise (either already a genuine Kekule value, or, for a
/// caller that never Kekulized an unsupported/gap ring first, still
/// `Aromatic`). Shared by both [`build_molecule_from_model`] and
/// [`build_molecule_from_model_authoritative`] -- this part of the
/// computation never differed between the two; only the ATOM flag's
/// promote-only-vs-authoritative decision does.
fn compute_bond_orders(mol: &Molecule, model: &AromaticityModel) -> FxHashMap<BondIdx, BondOrder> {
    mol.bonds()
        .map(|(bidx, bond)| {
            let order = if model.is_bond_aromatic(bidx) {
                BondOrder::Aromatic
            } else {
                bond.order
            };
            (bidx, order)
        })
        .collect()
}

/// Shared "finish" step for [`build_molecule_from_model`] and
/// [`build_molecule_from_model_authoritative`]: given final per-atom
/// aromatic flags and per-bond orders already decided (the only place the
/// two variants differ), builds the normalized [`Molecule`] -- implicit-H
/// preservation, bond-direction stashing, and stereo-metadata copying are
/// identical either way.
fn finish_molecule_with_flags(
    mol: &Molecule,
    atom_aromatic: &FxHashSet<AtomIdx>,
    bond_orders: &FxHashMap<BondIdx, BondOrder>,
) -> Molecule {
    use chematic_core::{MoleculeBuilder, implicit_hcount};

    // Implicit-H counts computed BEFORE bond orders are normalized below, for
    // organic-subset atoms without an explicit bracket H count. Needed because
    // normalizing every aromatic-model bond to `BondOrder::Aromatic` (below)
    // discards the Kekule Single/Double pattern that distinguishes a
    // lone-pair-donating "pyrrole-type" heteroatom (2 ring single bonds pre-
    // normalization, needs 1 implicit H) from a "pyridine-type" one (1 ring
    // single + 1 ring double, needs 0) -- post-normalization both look
    // identical (aromatic, 2 aromatic-order ring bonds, no substituent), so
    // `implicit_hcount`'s aromatic-path heuristic (correct for SMILES that
    // was aromatic-written from the start, per OpenSMILES convention: bare
    // aromatic `n` is pyridine-type, pyrrole-type is always `[nH]`) silently
    // returns the wrong value for atoms that reach this function via
    // Kekule-then-perceive instead. This under-counts molecular weight and
    // formula, not just fingerprints/canonical SMILES.
    let pre_h: Vec<Option<u8>> = mol
        .atoms()
        .map(|(idx, atom)| {
            if atom.hydrogen_count.is_some() {
                None // already explicit; nothing to preserve
            } else {
                Some(implicit_hcount(mol, idx))
            }
        })
        .collect();

    let mut builder = MoleculeBuilder::new();
    for (idx, atom) in mol.atoms() {
        let mut a = atom.clone();
        a.aromatic = atom_aromatic.contains(&idx);
        builder.add_atom(a);
    }
    for (bidx, bond) in mol.bonds() {
        let order = bond_orders[&bidx];
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
    let normalized = builder.build();

    // Compare the pre-normalization implicit H against what the same
    // (already-tested, unmodified) `implicit_hcount` computes on the
    // normalized bonds; only atoms where normalization actually changed the
    // answer get an explicit H frozen in. Benzene CH and pyridine-type N
    // (heuristic already agrees) are left untouched -- no spurious bracket
    // notation for atoms that didn't need it.
    let needs_patch: Vec<(chematic_core::AtomIdx, u8)> = normalized
        .atoms()
        .filter_map(|(idx, _)| {
            let pre = pre_h[idx.0 as usize]?;
            let post = implicit_hcount(&normalized, idx);
            (pre != post).then_some((idx, pre))
        })
        .collect();
    if needs_patch.is_empty() {
        return normalized;
    }

    let mut patched = MoleculeBuilder::new();
    for (idx, atom) in normalized.atoms() {
        let mut a = atom.clone();
        if let Some(&(_, h)) = needs_patch.iter().find(|(pidx, _)| *pidx == idx) {
            a.hydrogen_count = Some(h);
        }
        patched.add_atom(a);
    }
    for (_bond_idx, bond) in normalized.bonds() {
        let _ = patched.add_bond(bond.atom1, bond.atom2, bond.order);
    }
    patched.copy_stereo_groups_from(&normalized);
    patched.copy_stereo_from(&normalized);
    patched.copy_bond_directions_from(&normalized);
    patched.build()
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
///   2. Has a double bond whose far atom is on NO ring at all (a genuine
///      exocyclic substituent, not a ring-fusion bond into a different ring)
///      and is a more electronegative atom (O/N/S) → 0π (its p-orbital
///      electrons are in the exocyclic π bond, e.g. the carbonyl carbon in
///      tropone/pyridone/pyranone). A double bond whose far atom lies in a
///      DIFFERENT ring (e.g. a fusion carbon whose own Kekule double bond
///      happens to point into the other ring of a fused bicyclic, as in
///      quinazoline/quinoxaline) is a ring bond, not a substituent, and
///      falls through to rule 3 instead — see `all_ring_bonds` below.
///   3. Otherwise (has an endocyclic Double/Aromatic bond, or a double bond
///      into another ring) → 1π.
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
/// - **P (15) / Se (34) / Te (52)**: analogous lone-pair donors; only in
///   [`AromaticityAlgorithm::RdkitLike`] mode.
/// - **Other elements**: None (unsupported).
fn ring_pi_electrons(
    mol: &Molecule,
    ring: &[AtomIdx],
    aromatic_context: &FxHashSet<AtomIdx>,
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
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
                if atom.charge > 0 {
                    // Cationic ring carbon (tropylium's `[cH+]`): empty
                    // p-orbital electron acceptor, 0π, regardless of
                    // representation -- mirrors RDKit's carbon-specific
                    // charge-sign flip (see `kekulization.rs`'s
                    // `atom_must_be_matched` doc comment for the same rule
                    // in the Kekule-matching layer) and this function's own
                    // symmetric anion rule below (charge == -1 => 2π).
                    0
                } else if !has_double_any {
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
                        !all_ring_bonds.contains(&bidx)
                            && mol.bond(bidx).order == BondOrder::Double
                            && matches!(mol.atom(nb).element.atomic_number(), 7 | 8 | 16)
                    })
                {
                    // Only double bond is a genuine exocyclic substituent (its
                    // bond is on NO ring at all, not merely "not in the ring
                    // currently being evaluated") to a more electronegative
                    // atom (O/N/S): p-orbital electrons sit in that exocyclic π
                    // bond, contributing 0π to the ring (e.g. carbonyl carbon
                    // in tropone/pyridone/pyranone). A double bond into a
                    // DIFFERENT ring (a ring-fusion bond, e.g. a quinazoline
                    // fusion carbon whose own Kekule double bond points at the
                    // other ring's N) is excluded by the `all_ring_bonds`
                    // check and falls through to the sp2 default below instead
                    // of being wrongly zeroed (K2b fused-diazine fix).
                    0
                } else {
                    1
                }
            }

            // Nitrogen
            7 => {
                if implicit_hcount(mol, atom_idx) > 0 && atom.charge <= 0 {
                    // Pyrrole-type N with H, neutral or anionic: lone pair → 2π.
                    2
                } else if has_explicit_double {
                    // Pyridine-type N with an explicit double bond → 1π. Also
                    // catches a protonated ring N (pyridinium's `[nH+]`): the
                    // added proton consumes the lone pair the H-count check
                    // above would otherwise have claimed, and
                    // `chematic_core::kekulize` (charge-aware per K1) routes
                    // such an atom to a real Kekule double bond, exactly like
                    // neutral pyridine's bare N -- so this branch is reached
                    // instead of the one above once `atom.charge <= 0` fails.
                    1
                } else if total_degree == 3 && ring_degree < total_degree && atom.charge <= 0 {
                    // N with no H, no explicit double bond, all three σ-bonds
                    // exactly filling its valence (3), and neutral/anionic: a
                    // bridgehead N shared by two fused rings (e.g. indolizine)
                    // and a substituted pyrrole-type N (e.g. N-methylpyrrole,
                    // N-glycosylated purine/pyrimidine) have the identical
                    // local shape — the lone pair occupies the p orbital → 2π
                    // either way. Whether the ring this atom sits in is
                    // actually aromatic is decided by the overall 4n+2 sum below, not
                    // by inspecting the substituent: an imide N (phthalimide) still
                    // correctly comes out non-aromatic because its ring's carbonyl
                    // carbons contribute 0π each (exocyclic C=O rule above), giving
                    // 4π total, not 4n+2. The `charge <= 0` guard keeps a charged
                    // N with an H (pyridinium's `[nH+]`, degree 3 = 2 ring + 1 H)
                    // from being wrongly routed here in the aromatic-bond
                    // (pre-Kekulization) representation, where it has no
                    // explicit double bond to be caught by the branch above —
                    // it falls through to the pyridine-type branch below instead.
                    2
                } else if has_aromatic_in_ring {
                    // N in an aromatic ring (pre-kekulization input) without an
                    // explicit double bond and not a bridgehead → pyridine-like
                    // → 1π. Also the protonated-N fallback for the aromatic-bond
                    // representation (see the guards above).
                    1
                } else {
                    // Cannot determine pi contribution.
                    return None;
                }
            }

            // Oxygen / sulfur: lone-pair donor, must be 2-connected in the ring
            // -- *unless* a positive charge (pyrylium's `[o+]`) has consumed
            // the lone pair, in which case it needs pyridine-type treatment
            // (1π via its own ring double/aromatic bond) instead, mirroring
            // `kekulization.rs`'s charge-aware donor-exemption rule (K1).
            8 | 16 => {
                if atom.charge > 0 {
                    if has_explicit_double || has_aromatic_in_ring {
                        1
                    } else {
                        return None;
                    }
                } else {
                    if ring_degree != 2 {
                        return None;
                    }
                    // Sulfoxide/sulfone: exocyclic S=O ties up the lone pair; cannot donate 2π
                    if an == 16
                        && mol.neighbors(atom_idx).any(|(nb, bidx)| {
                            !ring_atom_set.contains(&nb)
                                && mol.bond(bidx).order == BondOrder::Double
                        })
                    {
                        return None;
                    }
                    2
                }
            }

            // P (15) / Se (34) / Te (52): heteroatom lone-pair donors (2π),
            // analogous to S. Only recognised in RdkitLike mode. P-H and
            // substituted P in a five-membered ring are the phosphole
            // counterparts of pyrrole; the ring-degree and exocyclic-double
            // guards keep hypervalent/exocyclic forms fail-closed.
            15 | 34 | 52 => {
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
// Diagnostic trace (Aromaticity-A1-0) — observational only, no production
// behavior change. `ring_pi_electrons` above is untouched and remains the
// single source of truth for `assign_aromaticity_ex`'s actual decisions;
// this is a parallel, read-only explanation layer for `component/atom/reason`
// tracing, used by `aromaticity_a1_0_report` and the corpus diagnostics in
// `validation/aromaticity_a1_0_corpus.jsonl`. See `docs/rfcs/aromaticity_a1_rfc.md`.
// ---------------------------------------------------------------------------

/// Reason a ring atom contributes (or fails to contribute) pi electrons,
/// mirroring `ring_pi_electrons`'s branches one-to-one. Purely diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContributionReason {
    /// Already aromatic from a previous Pass 1/Pass 2 ring: contributes 1π unconditionally.
    AlreadyAromaticContext,
    /// Carbon with an endocyclic double/aromatic bond: 1π.
    CarbonEndocyclicDouble,
    /// Carbon whose only double bond is exocyclic to O/N/S: 0π (e.g. a carbonyl carbon).
    CarbonExocyclicHeteroatomDouble,
    /// Carbanion with no double bond: 2π (lone pair).
    CarbonCarbanionLonePair,
    /// Cationic ring carbon (e.g. tropylium's `[cH+]`): empty p-orbital
    /// electron acceptor, 0π, regardless of representation (Kekule or
    /// aromatic-bond) -- mirrors `CarbonCarbanionLonePair`'s anion rule at
    /// the opposite electron-count extreme.
    CarbonCationVacant,
    /// sp3 carbon (no double bond, not a carbanion): ineligible.
    CarbonSp3Ineligible,
    /// Pyrrole-type N with an H, neutral or anionic: 2π.
    NitrogenPyrroleTypeH,
    /// Pyridine-type N with an explicit double bond (bare, or protonated
    /// N-H+ once it has a Kekule double bond): 1π.
    NitrogenPyridineTypeExplicitDouble,
    /// Bridgehead N (or N-substituted azole N), neutral or anionic:
    /// all-sigma valence, lone pair in p orbital: 2π.
    NitrogenBridgeheadOrSubstitutedLonePair,
    /// N with an in-ring aromatic bond, not a bridgehead (pyridine-type
    /// notation, or a charged N-H+ in aromatic-bond representation): 1π.
    NitrogenAromaticInRing,
    /// N matching none of the above rules: ineligible.
    NitrogenIneligible,
    /// O/S/Se/Te lone-pair donor, neutral or anionic, ring-degree 2: 2π.
    ChalcogenLonePair,
    /// P lone-pair donor in the opt-in RDKit-compatible model: 2π.
    PnictogenOrChalcogenLonePair,
    /// Charged O/S (e.g. pyrylium's `[o+]`): the positive charge consumes
    /// the lone pair, so this atom needs pyridine-type treatment (1π via
    /// its own ring double/aromatic bond) instead of donating 2π.
    ChalcogenCationPyridineType,
    /// O/S/Se/Te with the wrong ring degree, an exocyclic X=O, (Se/Te)
    /// non-RdkitLike mode, or a charged O/S with no ring double/aromatic
    /// bond to fall back on: ineligible.
    ChalcogenIneligible,
    /// Element not supported by the model: ineligible.
    UnsupportedElement,
}

impl ContributionReason {
    /// Whether this reason is an eligible contribution (matches
    /// `ring_pi_electrons` returning `Some`) rather than one that disqualifies
    /// the whole ring (matches it returning `None`).
    pub fn is_eligible(self) -> bool {
        !matches!(
            self,
            ContributionReason::CarbonSp3Ineligible
                | ContributionReason::NitrogenIneligible
                | ContributionReason::ChalcogenIneligible
                | ContributionReason::UnsupportedElement
        )
    }

    /// Coarse `PiEligibility` bucket for this fine-grained reason
    /// (Aromaticity-A1-1a). `AlreadyAromaticContext` has no single fixed
    /// bucket -- it always carries exactly 1π, so it maps to `OneElectron`.
    pub fn eligibility(self) -> PiEligibility {
        use ContributionReason::*;
        match self {
            AlreadyAromaticContext
            | CarbonEndocyclicDouble
            | NitrogenPyridineTypeExplicitDouble
            | NitrogenAromaticInRing
            | ChalcogenCationPyridineType => PiEligibility::OneElectron,
            CarbonCarbanionLonePair
            | NitrogenPyrroleTypeH
            | NitrogenBridgeheadOrSubstitutedLonePair
            | ChalcogenLonePair
            | PnictogenOrChalcogenLonePair => PiEligibility::LonePairDonor,
            CarbonExocyclicHeteroatomDouble | CarbonCationVacant => PiEligibility::ZeroElectron,
            CarbonSp3Ineligible | NitrogenIneligible | ChalcogenIneligible | UnsupportedElement => {
                PiEligibility::Ineligible
            }
        }
    }
}

/// Coarse per-atom pi-eligibility bucket (Aromaticity-A1-1a). A summary view
/// over [`ContributionReason`]'s finer-grained rules -- `electrons()` gives
/// the electron count implied by each bucket.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PiEligibility {
    /// Contributes exactly 1π (e.g. an endocyclic double/aromatic bond).
    OneElectron,
    /// Contributes 2π (a lone pair: pyrrole-type N, chalcogen, bridgehead N, carbanion).
    LonePairDonor,
    /// Contributes 0π but is still sp2 (p-orbital spent on an exocyclic multiple bond).
    ZeroElectron,
    /// Not eligible to be part of any conjugated system (e.g. sp3).
    Ineligible,
}

impl PiEligibility {
    /// Electron count implied by this bucket, or `None` for `Ineligible`.
    pub fn electrons(self) -> Option<u8> {
        match self {
            PiEligibility::OneElectron => Some(1),
            PiEligibility::LonePairDonor => Some(2),
            PiEligibility::ZeroElectron => Some(0),
            PiEligibility::Ineligible => None,
        }
    }
}

/// A candidate conjugated system: some atoms/bonds evaluated together as one
/// pi-electron-counting problem (Aromaticity-A1-1a). Two distinct uses:
/// - a single SSSR/augmented ring, reinterpreted as a trivial one-ring
///   candidate (what `trace_ring_pi_electrons` builds today);
/// - a genuine multi-ring fused envelope, built by
///   [`build_conjugated_components`] as a connected component of the
///   "conjugation graph" (double/aromatic-bonded atoms, plus lone-pair-donor
///   atoms bridging across single bonds) -- the azulene-class candidate
///   `augmented_ring_set`'s own docstring already named as future work
///   ("candidate rings = SSSR ∪ fused envelopes").
#[derive(Debug, Clone)]
pub struct ConjugatedComponent {
    pub atoms: Vec<AtomIdx>,
    pub bonds: Vec<BondIdx>,
    /// Ring indices (into whatever ring list the caller built this from) this
    /// candidate derives from -- one entry for a plain single-ring candidate,
    /// 2+ for a fused envelope spanning multiple rings.
    pub source_rings: Vec<usize>,
}

impl ConjugatedComponent {
    /// Build a trivial single-ring candidate from one ring's atom list (no
    /// bond list needed by [`evaluate_atom_pi_contribution`], which only
    /// consults `atoms` membership).
    fn from_ring(ring: &[AtomIdx], ring_idx: usize) -> Self {
        ConjugatedComponent {
            atoms: ring.to_vec(),
            bonds: Vec::new(),
            source_rings: vec![ring_idx],
        }
    }
}

/// The full per-atom decision from [`evaluate_atom_pi_contribution`]: the
/// coarse eligibility bucket plus the specific rule that produced it.
#[derive(Debug, Clone, Copy)]
pub struct ContributionDecision {
    pub eligibility: PiEligibility,
    pub reason: ContributionReason,
}

impl ContributionDecision {
    pub fn electrons(&self) -> Option<u8> {
        self.eligibility.electrons()
    }
}

/// Per-atom trace entry from [`trace_ring_pi_electrons`].
#[derive(Debug, Clone, Copy)]
pub struct AtomElectronTrace {
    pub atom_idx: AtomIdx,
    /// `None` iff `reason.is_eligible()` is false.
    pub contribution: Option<u8>,
    pub reason: ContributionReason,
}

/// Full per-atom pi-electron trace for one ring — the diagnostic twin of
/// [`ring_pi_electrons`]. Unlike `ring_pi_electrons` (which returns `None` at
/// the first ineligible atom), this always scans every atom so a caller can
/// see exactly which atom(s) disqualify a ring, not just that one did.
#[derive(Debug, Clone)]
pub struct RingElectronTrace {
    pub atoms: Vec<AtomElectronTrace>,
    /// `Some(sum)` iff every atom was eligible — must equal
    /// `ring_pi_electrons(mol, ring, aromatic_context, algo, all_ring_bonds)`
    /// for the same inputs (checked by
    /// `trace_matches_ring_pi_electrons_on_corpus` below).
    pub total: Option<u32>,
}

/// Diagnostic twin of [`ring_pi_electrons`]: identical per-atom rules
/// (delegating to [`evaluate_atom_pi_contribution`], the single source of
/// truth for both this trace and any future experimental production path —
/// see `docs/rfcs/aromaticity_a1_rfc.md`'s A1-1a section), but returns a full
/// trace instead of a single early-exiting `Option<u32>`. Does not call,
/// wrap, or change `ring_pi_electrons` itself — zero effect on
/// `assign_aromaticity_ex`'s behavior. `trace_matches_ring_pi_electrons_on_corpus`
/// is the anti-drift guard that keeps this and `ring_pi_electrons` in sync.
pub fn trace_ring_pi_electrons(
    mol: &Molecule,
    ring: &[AtomIdx],
    aromatic_context: &FxHashSet<AtomIdx>,
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> RingElectronTrace {
    let component = ConjugatedComponent::from_ring(ring, 0);
    let mut atoms = Vec::with_capacity(ring.len());
    let mut total: Option<u32> = Some(0);

    for &atom_idx in ring {
        let (contribution, reason) = if aromatic_context.contains(&atom_idx) {
            (Some(1u8), ContributionReason::AlreadyAromaticContext)
        } else {
            let decision =
                evaluate_atom_pi_contribution(mol, atom_idx, &component, algo, all_ring_bonds);
            (decision.electrons(), decision.reason)
        };

        total = match (total, contribution) {
            (Some(t), Some(c)) => Some(t + c as u32),
            _ => None,
        };

        atoms.push(AtomElectronTrace {
            atom_idx,
            contribution,
            reason,
        });
    }

    RingElectronTrace { atoms, total }
}

/// Single source of truth for per-atom pi-electron contribution
/// (Aromaticity-A1-1a): identical rules to `ring_pi_electrons`'s match arms,
/// condition-for-condition, parameterized by an arbitrary candidate
/// [`ConjugatedComponent`] instead of one fixed SSSR ring — the same
/// function evaluates a plain single-ring candidate (via
/// `ConjugatedComponent::from_ring`) or a genuine multi-ring fused envelope
/// (via `build_conjugated_components`) identically. Currently called by
/// `trace_ring_pi_electrons` only — NOT wired into `ring_pi_electrons` or
/// `assign_aromaticity_ex` (that wiring, behind a new opt-in
/// `AromaticityAlgorithm` variant, is Aromaticity-A1-1b, not this round).
pub fn evaluate_atom_pi_contribution(
    mol: &Molecule,
    atom_idx: AtomIdx,
    component: &ConjugatedComponent,
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> ContributionDecision {
    let component_atoms: FxHashSet<AtomIdx> = component.atoms.iter().copied().collect();
    let (_electrons, reason) =
        evaluate_atom_pi_contribution_inner(mol, atom_idx, &component_atoms, algo, all_ring_bonds);
    // `reason.eligibility().electrons()` is asserted equal to `_electrons`
    // for every branch by `contribution_decision_electrons_match_inner_on_corpus`.
    ContributionDecision {
        eligibility: reason.eligibility(),
        reason,
    }
}

/// Per-atom contribution logic, mirroring `ring_pi_electrons`'s match arms
/// condition-for-condition, but returning a reason alongside the
/// contribution instead of returning early on `None`.
fn evaluate_atom_pi_contribution_inner(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_atom_set: &FxHashSet<AtomIdx>,
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> (Option<u8>, ContributionReason) {
    let atom = mol.atom(atom_idx);
    let an = atom.element.atomic_number();

    let ring_degree = mol
        .neighbors(atom_idx)
        .filter(|(nb, _)| ring_atom_set.contains(nb))
        .count();
    let total_degree = mol.degree(atom_idx);

    let has_explicit_double = mol
        .neighbors(atom_idx)
        .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Double);
    let has_double_any = has_explicit_double
        || mol
            .neighbors(atom_idx)
            .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Aromatic);
    let has_aromatic_in_ring = mol
        .neighbors(atom_idx)
        .filter(|(nb, _)| ring_atom_set.contains(nb))
        .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Aromatic);

    match an {
        6 => {
            if atom.charge > 0 {
                (Some(0), ContributionReason::CarbonCationVacant)
            } else if !has_double_any {
                if atom.charge == -1 {
                    (Some(2), ContributionReason::CarbonCarbanionLonePair)
                } else {
                    (None, ContributionReason::CarbonSp3Ineligible)
                }
            } else if has_explicit_double
                && !has_aromatic_in_ring
                && !mol.neighbors(atom_idx).any(|(nb, bidx)| {
                    ring_atom_set.contains(&nb) && mol.bond(bidx).order == BondOrder::Double
                })
                && mol.neighbors(atom_idx).any(|(nb, bidx)| {
                    !all_ring_bonds.contains(&bidx)
                        && mol.bond(bidx).order == BondOrder::Double
                        && matches!(mol.atom(nb).element.atomic_number(), 7 | 8 | 16)
                })
            {
                // See `ring_pi_electrons`'s identical rule (K2b fused-diazine
                // fix): a double bond into a DIFFERENT ring is a ring-fusion
                // bond, not a genuine exocyclic substituent, and must not be
                // zeroed here either -- this function must stay in lockstep
                // with `ring_pi_electrons` (checked by
                // `trace_matches_ring_pi_electrons_on_corpus`).
                (Some(0), ContributionReason::CarbonExocyclicHeteroatomDouble)
            } else {
                (Some(1), ContributionReason::CarbonEndocyclicDouble)
            }
        }
        7 => {
            if implicit_hcount(mol, atom_idx) > 0 && atom.charge <= 0 {
                (Some(2), ContributionReason::NitrogenPyrroleTypeH)
            } else if has_explicit_double {
                (
                    Some(1),
                    ContributionReason::NitrogenPyridineTypeExplicitDouble,
                )
            } else if total_degree == 3 && ring_degree < total_degree && atom.charge <= 0 {
                (
                    Some(2),
                    ContributionReason::NitrogenBridgeheadOrSubstitutedLonePair,
                )
            } else if has_aromatic_in_ring {
                (Some(1), ContributionReason::NitrogenAromaticInRing)
            } else {
                (None, ContributionReason::NitrogenIneligible)
            }
        }
        8 | 16 => {
            if atom.charge > 0 {
                if has_explicit_double || has_aromatic_in_ring {
                    (Some(1), ContributionReason::ChalcogenCationPyridineType)
                } else {
                    (None, ContributionReason::ChalcogenIneligible)
                }
            } else {
                let exocyclic_double = an == 16
                    && mol.neighbors(atom_idx).any(|(nb, bidx)| {
                        !ring_atom_set.contains(&nb) && mol.bond(bidx).order == BondOrder::Double
                    });
                if ring_degree != 2 || exocyclic_double {
                    (None, ContributionReason::ChalcogenIneligible)
                } else {
                    (Some(2), ContributionReason::ChalcogenLonePair)
                }
            }
        }
        15 | 34 | 52 => {
            let exocyclic_double = mol.neighbors(atom_idx).any(|(nb, bidx)| {
                !ring_atom_set.contains(&nb) && mol.bond(bidx).order == BondOrder::Double
            });
            if algo != AromaticityAlgorithm::RdkitLike || ring_degree != 2 || exocyclic_double {
                (None, ContributionReason::ChalcogenIneligible)
            } else {
                (Some(2), ContributionReason::PnictogenOrChalcogenLonePair)
            }
        }
        _ => (None, ContributionReason::UnsupportedElement),
    }
}

/// Evaluate an atom's pi contribution using its "home ring" within a
/// (possibly multi-ring) candidate, instead of the candidate's flattened
/// atom set directly: tries each of `candidate.source_rings` that actually
/// contains the atom, evaluating against *that one ring's own* atom set, and
/// returns the first eligible result found. Falls back to evaluating
/// directly against the flattened `candidate` if `source_rings` is empty or
/// none of them contain the atom (shouldn't happen for well-formed
/// candidates, but keeps this total rather than panicking).
///
/// Needed because degree-sensitive rules (the N bridgehead/substituted-azole
/// rule, `total_degree == 3 && ring_degree < total_degree`) test "does this
/// atom have a bond that points outside THIS ring" -- a genuine multi-ring
/// bridgehead's every bond is "in-family" once the evaluation context is the
/// flattened whole envelope (every neighbor is, by construction, some other
/// family member), which silently defeats that test and makes a real
/// bridgehead N (e.g. indolizine's) look `Ineligible`. Evaluating against
/// one constituent ring at a time preserves the rule's original, correct,
/// per-ring meaning even when the candidate spans multiple rings. This does
/// **not** attempt to resolve whether a bridgehead's lone-pair credit is
/// *legitimately shared* between two rings that are both otherwise valid vs.
/// wrongly borrowed by one ring from another that's actually broken (e.g.
/// by an sp3 atom) -- that is a distinct, harder, open question, deliberately
/// left to Aromaticity-A1-1b (see `docs/rfcs/aromaticity_a1_rfc.md`).
fn evaluate_atom_via_home_ring(
    mol: &Molecule,
    atom_idx: AtomIdx,
    candidate: &ConjugatedComponent,
    rings: &[Vec<AtomIdx>],
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> ContributionDecision {
    let mut last = None;
    for &ri in &candidate.source_rings {
        if !rings[ri].contains(&atom_idx) {
            continue;
        }
        let home = ConjugatedComponent::from_ring(&rings[ri], ri);
        let decision = evaluate_atom_pi_contribution(mol, atom_idx, &home, algo, all_ring_bonds);
        if decision.electrons().is_some() {
            return decision;
        }
        last = Some(decision);
    }
    last.unwrap_or_else(|| {
        evaluate_atom_pi_contribution(mol, atom_idx, candidate, algo, all_ring_bonds)
    })
}

/// Build genuine multi-ring conjugated-system candidates (Aromaticity-A1-1a):
/// connected components of the "conjugation graph" over each ring family's
/// atoms -- nodes are atoms whose eligibility (evaluated per-atom against its
/// own home ring, via `evaluate_atom_via_home_ring` -- not the flattened
/// family) is not `Ineligible`; edges are any bond (single, double, or
/// aromatic) between two independently-eligible family atoms: ordinary
/// carbon-carbon single-bond conjugation (butadiene's C=C-C=C middle bond,
/// styrene's vinyl-to-phenyl bond) connects just as directly as a
/// lone-pair-donor heteroatom bridging a sigma bond.
///
/// A pure candidate *generator* -- callers (currently only
/// `exhaustive_aromaticity_oracle`) still run full 4n+2 electron counting on
/// each result. Only components spanning 2+ of a family's rings are
/// returned: a single unfused ring is already covered by
/// `ConjugatedComponent::from_ring`, so this only adds the fused-envelope
/// candidates `augmented_ring_set`'s docstring named as future work
/// ("candidate rings = SSSR ∪ fused envelopes").
pub fn build_conjugated_components(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    ring_families: &[RingFamily],
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> Vec<ConjugatedComponent> {
    let mut out = Vec::new();

    for family in ring_families {
        if family.ring_indices.len() < 2 {
            continue; // single-ring families add nothing beyond from_ring.
        }
        let family_component = ConjugatedComponent {
            atoms: family.atoms.clone(),
            bonds: Vec::new(),
            source_rings: family.ring_indices.clone(),
        };

        // Eligibility per atom, evaluated against its *home* constituent
        // ring (not the flattened family) -- see `evaluate_atom_via_home_ring`'s
        // doc comment for why the flattened version breaks degree-sensitive
        // rules (bridgehead N) for any atom whose every bond happens to be
        // "in-family" once the family itself is the context.
        let eligible: FxHashMap<AtomIdx, bool> = family
            .atoms
            .iter()
            .map(|&a| {
                let decision = evaluate_atom_via_home_ring(
                    mol,
                    a,
                    &family_component,
                    rings,
                    algo,
                    all_ring_bonds,
                );
                (a, decision.electrons().is_some())
            })
            .collect();
        // Union-find over eligible family atoms, connected by conjugation edges.
        let atoms: Vec<AtomIdx> = family.atoms.clone();
        let index_of: FxHashMap<AtomIdx, usize> =
            atoms.iter().enumerate().map(|(i, &a)| (a, i)).collect();
        let mut parent: Vec<usize> = (0..atoms.len()).collect();
        fn find(parent: &mut [usize], x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }
        fn union(parent: &mut [usize], x: usize, y: usize) {
            let (px, py) = (find(parent, x), find(parent, y));
            if px != py {
                parent[px] = py;
            }
        }

        // Any bond (single, double, or aromatic) between two independently
        // eligible atoms conjugation-connects them: two sp2 atoms bridge
        // across a single bond exactly like butadiene's C=C-C=C middle bond
        // or styrene's vinyl-to-phenyl bond -- ordinary carbon-carbon
        // conjugation, not just lone-pair-donor bridging. (First version of
        // this rule only bridged single bonds via a `LonePairDonor`
        // endpoint, which is too narrow: it left azulene's all-carbon
        // alternating single/double perimeter as 5 disconnected 2-atom
        // pairs, never forming the one 10-atom fused-envelope candidate it
        // needs -- caught by `exhaustive_aromaticity_oracle` returning an
        // empty set for azulene instead of the whole ring.) The
        // `is_lone_pair_donor` check is now unused for connectivity, kept
        // only where a NON-eligible atom's neighbor still needs distinguishing
        // (none currently) -- eligibility alone (both endpoints not
        // `Ineligible`) is the connectivity condition; bond order still fully
        // determines each atom's *electron count* via
        // `evaluate_atom_pi_contribution`, just not graph connectivity.
        let family_atom_set: FxHashSet<AtomIdx> = family.atoms.iter().copied().collect();
        let mut conjugation_bonds: Vec<BondIdx> = Vec::new();
        for &a in &atoms {
            if !eligible[&a] {
                continue;
            }
            for (nb, bidx) in mol.neighbors(a) {
                if !family_atom_set.contains(&nb) || !eligible.get(&nb).copied().unwrap_or(false) {
                    continue;
                }
                // Both endpoints eligible -> connected (see comment above).
                union(&mut parent, index_of[&a], index_of[&nb]);
                conjugation_bonds.push(bidx);
            }
        }

        let mut groups: FxHashMap<usize, Vec<AtomIdx>> = FxHashMap::default();
        for &a in &atoms {
            if !eligible[&a] {
                continue;
            }
            let root = find(&mut parent, index_of[&a]);
            groups.entry(root).or_default().push(a);
        }

        for group_atoms in groups.into_values() {
            let group_set: FxHashSet<AtomIdx> = group_atoms.iter().copied().collect();
            let source_rings: Vec<usize> = family
                .ring_indices
                .iter()
                .copied()
                .filter(|&ri| rings[ri].iter().all(|a| group_set.contains(a)))
                .collect();
            if source_rings.len() < 2 {
                continue; // doesn't actually span multiple full rings.
            }
            let group_bonds: Vec<BondIdx> = conjugation_bonds
                .iter()
                .copied()
                .filter(|&bidx| {
                    let b = mol.bond(bidx);
                    group_set.contains(&b.atom1) && group_set.contains(&b.atom2)
                })
                .collect();
            out.push(ConjugatedComponent {
                atoms: group_atoms,
                bonds: group_bonds,
                source_rings,
            });
        }
    }

    out
}

/// Test/diagnostic-only exhaustive-candidate reference oracle
/// (Aromaticity-A1-1a) — **not** used by production or by
/// `trace_ring_pi_electrons`. Evaluates every SSSR/augmented ring AND every
/// multi-ring fused-envelope candidate from `build_conjugated_components`,
/// marking an atom/bond aromatic if ANY candidate containing it
/// independently satisfies 4n+2 via `evaluate_atom_pi_contribution`'s
/// per-atom rules — every candidate is evaluated from a clean slate, with NO
/// `aromatic_context` bootstrapping at all (unlike `assign_aromaticity_ex`'s
/// production Pass 1/Pass 2). Exists to cross-check hypotheses about which
/// per-atom rule needs to change, per the MANCUDE-style bounded-enumeration
/// precedent — see `docs/rfcs/aromaticity_a1_rfc.md`'s A1-1a section.
/// Deliberately simple/slow: O(rings + fused envelopes) candidates, no
/// attempt at Pass-2-style iteration, memoization, or performance tuning.
pub fn exhaustive_aromaticity_oracle(
    mol: &Molecule,
    algo: AromaticityAlgorithm,
) -> (FxHashSet<AtomIdx>, FxHashSet<BondIdx>) {
    let sssr = find_sssr(mol);
    let rings = augmented_ring_set(mol, sssr.rings());
    let families = crate::ring_family::find_ring_families_over(mol, &rings);
    let all_ring_bonds: FxHashSet<BondIdx> =
        rings.iter().flat_map(|r| ring_bond_set(mol, r)).collect();

    let mut candidates: Vec<ConjugatedComponent> = rings
        .iter()
        .enumerate()
        .map(|(i, r)| ConjugatedComponent::from_ring(r, i))
        .collect();
    candidates.extend(build_conjugated_components(
        mol,
        &rings,
        &families,
        algo,
        &all_ring_bonds,
    ));

    let mut aromatic_atoms: FxHashSet<AtomIdx> = FxHashSet::default();
    let mut aromatic_bonds: FxHashSet<BondIdx> = FxHashSet::default();

    for candidate in &candidates {
        let mut total: Option<u32> = Some(0);
        for &atom_idx in &candidate.atoms {
            // Multi-ring candidates evaluate each atom against its home ring
            // (see `evaluate_atom_via_home_ring`'s doc comment); single-ring
            // candidates fall through to the same code path with exactly one
            // source ring, unchanged from evaluating against `candidate` directly.
            let decision = evaluate_atom_via_home_ring(
                mol,
                atom_idx,
                candidate,
                &rings,
                algo,
                &all_ring_bonds,
            );
            total = match (total, decision.electrons()) {
                (Some(t), Some(e)) => Some(t + e as u32),
                _ => None,
            };
        }
        let Some(pi) = total else { continue };
        let (cls, _) = classify_ring_aromaticity(pi);
        if !matches!(cls, RingAromaticity::Aromatic) {
            continue;
        }
        for &a in &candidate.atoms {
            aromatic_atoms.insert(a);
        }
        for &a in &candidate.atoms {
            for (nb, bidx) in mol.neighbors(a) {
                if candidate.atoms.contains(&nb)
                    && matches!(
                        mol.bond(bidx).order,
                        BondOrder::Double | BondOrder::Aromatic
                    )
                {
                    aromatic_bonds.insert(bidx);
                }
            }
        }
    }

    (aromatic_atoms, aromatic_bonds)
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

    /// Same ring as `pyrrole_kekule()`, but the N has NO explicit
    /// `hydrogen_count` — matching how the SMILES parser actually builds a
    /// bare, non-bracket `N` (e.g. from `Chem.Kekulize` + non-canonical
    /// `MolToSmiles(kekuleSmiles=True)` round-tripping an `[nH]`-written
    /// pyrrole/imidazole/purine nitrogen). `pyrrole_kekule()` above sidesteps
    /// the bug this reproduces by setting `hydrogen_count` manually.
    fn pyrrole_kekule_implicit_h() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(Atom::new(Element::N));
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
    fn test_apply_aromaticity_preserves_pyrrole_nh_implicit_hydrogen() {
        // Regression test: apply_aromaticity_ex() normalizes all aromatic-
        // model ring bonds to BondOrder::Aromatic, which discards the
        // Kekule Single/Double pattern that distinguishes a pyrrole-type N
        // (needs 1 implicit H) from a pyridine-type N (needs 0) once both
        // have exactly 2 aromatic-order ring bonds and no explicit bracket
        // H count. Without preserving the pre-normalization value,
        // implicit_hcount() on the perceived molecule silently returns 0
        // instead of 1 for the unsubstituted pyrrole N -- wrong molecular
        // formula/weight, and a representation-dependent divergence from
        // the same molecule parsed directly from aromatic-written SMILES
        // (where `[nH]`'s bracket H count is correct by construction).
        let mol = pyrrole_kekule_implicit_h();
        let n_idx = AtomIdx(0);
        assert_eq!(
            implicit_hcount(&mol, n_idx),
            1,
            "pre-normalization: bare N with 2 single ring bonds must show 1 implicit H"
        );

        let perceived = apply_aromaticity(&mol);
        assert!(perceived.atom(n_idx).aromatic, "ring N must be aromatic");
        assert_eq!(
            implicit_hcount(&perceived, n_idx),
            1,
            "post-apply_aromaticity: pyrrole N must still show 1 implicit H, not 0"
        );
    }

    #[test]
    fn test_apply_aromaticity_does_not_add_h_to_pyridine_type_n() {
        // Sibling check to the pyrrole regression above: a pyridine-type
        // ring N (1 ring single + 1 ring double pre-normalization, no H)
        // must NOT gain a spurious implicit H from the preservation logic --
        // its pre- and post-normalization implicit_hcount already agree
        // (both 0), so it must be left untouched.
        let mol = pyridine_kekule();
        let n_idx = AtomIdx(0);
        assert_eq!(implicit_hcount(&mol, n_idx), 0);

        let perceived = apply_aromaticity(&mol);
        assert!(perceived.atom(n_idx).aromatic);
        assert_eq!(implicit_hcount(&perceived, n_idx), 0);
        assert_eq!(
            perceived.atom(n_idx).hydrogen_count,
            None,
            "pyridine N must not gain an explicit hydrogen_count -- would force spurious bracket notation"
        );
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
    #[ignore = "PROVISIONAL: regressed by the Horton SSSR fix, see comment below"]
    fn test_purine_aromatic() {
        // c1cnc2[nH]cnc2n1 — purine: 9 atoms, kekulizable
        //
        // Regressed by the Horton SSSR rewrite (confirmed passing on the old
        // single-spanning-tree find_sssr, failing only after Horton; see
        // debug dump captured during diagnosis). Root cause, empirically
        // confirmed: the 6-membered ring (pyrimidine-type) passes Pass 1
        // alone (6π) and marks its atoms aromatic. The 5-membered ring
        // (imidazole-type) evaluates to 4π in isolation — its two fusion
        // carbons each have their only double bond exocyclic to a ring N,
        // which the exocyclic-to-heteroatom rule scores as 0π — and 4π trips
        // `classify_ring_aromaticity`'s "4n → Antiaromatic" branch. Pass 1
        // treats Antiaromatic as definitive and never retries it in Pass 2,
        // even though the fusion carbons would each contribute 1π (not 0π)
        // once `aromatic_context` recognizes them as already-aromatic — that
        // recount gives 6π (aromatic). The old, non-minimal SSSR never hit
        // this path because it fed a different (structurally wrong) ring set
        // into Pass 1 in the first place.
        //
        // Fix belongs in the aromatic_context-removal PR (see
        // greedy-hopping-crescent.md step 5), not here: retrying
        // Antiaromatic rings in Pass 2 is a real fix, but must not be
        // bundled into the SSSR PR per the "measure free recoveries with
        // zero aromaticity.rs changes" staging requirement.
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

    // ── K2a: charge-aware ring_pi_electrons -- tropylium/imidazolium/
    // pyridinium/pyrylium now genuinely confirmed aromatic by the raw
    // Huckel model itself (not just a stale parser flag surviving), under
    // BOTH documented calling conventions (`apply_aromaticity`'s own doc
    // comment: "may be kekulized... or may retain Aromatic bond orders from
    // the SMILES parser"). RDKit-verified: all four are aromatic cations,
    // all-atom/all-bond, per rdkit==2026.03.3 (see
    // docs/rfcs/aromaticity_rdkit_parity_rfc.md and the K2a PR description for
    // the full 40-fixture oracle re-run against a live RDKit).
    //
    // K1 (fix/kekulize-charge-aware-k1, already merged) made
    // chematic_core::kekulize() succeed for all four; this fix is the
    // separate, independent charge-blindness bug in the Huckel
    // pi-electron-counting layer (`ring_pi_electrons`) that K1 explicitly
    // did not touch. Deliberately does NOT touch `build_molecule_from_model`
    // (that promote-only-vs-demote question is tracked separately as K2b) --
    // these four fixtures need no demotion at all: their atom flags were
    // already `true` from the aromatic-notation parse, and once the model
    // itself confirms the ring, the EXISTING promote-only bond loop already
    // correctly promotes their bonds to `Aromatic` for the first time. That
    // is what actually fixes the pre-existing atom/bond flag inconsistency
    // for these four -- no demotion capability required.
    fn assert_fully_aromatic(mol: &Molecule, n: usize, label: &str) {
        let applied = apply_aromaticity(mol);
        for (idx, atom) in applied.atoms() {
            assert!(atom.aromatic, "{label}: atom {idx:?} should be aromatic");
        }
        assert_eq!(applied.atom_count(), n, "{label}: unexpected atom count");
        for (_, bond) in applied.bonds() {
            assert_eq!(
                bond.order,
                BondOrder::Aromatic,
                "{label}: every ring bond should end up Aromatic order"
            );
        }
    }

    #[test]
    fn test_tropylium_cation_aromatic_raw_and_kekulized() {
        let raw = chematic_smiles::parse("c1ccc[cH+]cc1").expect("valid SMILES");
        assert_fully_aromatic(&raw, 7, "tropylium (raw)");
        let kek = mol_kekulized("c1ccc[cH+]cc1");
        assert_fully_aromatic(&kek, 7, "tropylium (kekulized)");
        assert_eq!(
            assign_aromaticity(&raw).aromatic_atom_count(),
            7,
            "tropylium: raw model itself must confirm all 7 atoms, not rely on a stale flag"
        );
        assert_eq!(
            assign_aromaticity(&kek).aromatic_atom_count(),
            7,
            "tropylium: kekulized model itself must confirm all 7 atoms"
        );
    }

    #[test]
    fn test_imidazolium_aromatic_raw_and_kekulized() {
        let raw = chematic_smiles::parse("c1c[nH+]c[nH]1").expect("valid SMILES");
        assert_fully_aromatic(&raw, 5, "imidazolium (raw)");
        let kek = mol_kekulized("c1c[nH+]c[nH]1");
        assert_fully_aromatic(&kek, 5, "imidazolium (kekulized)");
        assert_eq!(assign_aromaticity(&raw).aromatic_atom_count(), 5);
        assert_eq!(assign_aromaticity(&kek).aromatic_atom_count(), 5);
    }

    #[test]
    fn test_pyridinium_aromatic_raw_and_kekulized() {
        let raw = chematic_smiles::parse("c1cc[nH+]cc1").expect("valid SMILES");
        assert_fully_aromatic(&raw, 6, "pyridinium (raw)");
        let kek = mol_kekulized("c1cc[nH+]cc1");
        assert_fully_aromatic(&kek, 6, "pyridinium (kekulized)");
        assert_eq!(assign_aromaticity(&raw).aromatic_atom_count(), 6);
        assert_eq!(assign_aromaticity(&kek).aromatic_atom_count(), 6);
    }

    #[test]
    fn test_pyrylium_aromatic_raw_and_kekulized() {
        let raw = chematic_smiles::parse("c1cc[o+]cc1").expect("valid SMILES");
        assert_fully_aromatic(&raw, 6, "pyrylium (raw)");
        let kek = mol_kekulized("c1cc[o+]cc1");
        assert_fully_aromatic(&kek, 6, "pyrylium (kekulized)");
        assert_eq!(assign_aromaticity(&raw).aromatic_atom_count(), 6);
        assert_eq!(assign_aromaticity(&kek).aromatic_atom_count(), 6);
    }

    // ── K2a scope guard: tellurophene/phosphole are explicitly NOT fixed by
    // the charge-aware change above (they need real Se/Te/P electron-donor
    // support in the default Huckel engine, out of scope -- see the K2a/K2b
    // PR descriptions). Pin the current (still-gap) count so a future
    // change to this area doesn't silently start claiming these are fixed
    // without an explicit, source-grounded review.
    #[test]
    fn test_tellurophene_and_phosphole_still_unsupported_under_default_huckel() {
        let te = mol_kekulized("c1cc[te]c1");
        assert_eq!(
            assign_aromaticity(&te).aromatic_atom_count(),
            0,
            "tellurophene: still unsupported under default Huckel (K2a does not add Te support)"
        );
        let p = mol_kekulized("c1cc[pH]c1");
        assert_eq!(
            assign_aromaticity(&p).aromatic_atom_count(),
            0,
            "phosphole: still unsupported under default Huckel (K2a does not add P support)"
        );
    }

    // ── K2b fused-diazine fix (fix/aromaticity-flag-demotion-k2b follow-up) ─
    //
    // Opt-in only, via `assign_aromaticity_authoritative_experimental` --
    // per coordinator decision, `apply_aromaticity`/`apply_aromaticity_ex`
    // (and the plain `assign_aromaticity`/`assign_aromaticity_ex` they call)
    // stay byte-identical to their pre-K2b behavior. The
    // `test_known_gap_fused_diazine_exocyclic_misfire_antiaromatic` pin that
    // used to live here (asserting the DEFAULT engine's wrong 6/10 count) is
    // superseded by `test_default_engine_unaffected_by_fused_diazine_fix`
    // below (same assertion, renamed for clarity: this is now a permanent
    // "default stays reverted" guard, not a "known gap" pin -- the gap is
    // only closed for the opt-in engine, not fixed in the default at all).
    // The azulene pin further below is untouched either way (separate,
    // still-open, out-of-scope mechanism, never affected by this fix in
    // ANY engine).

    #[test]
    fn test_authoritative_experimental_fixes_fused_diazine_ring_fusion() {
        // c1cnc2ccccc2n1 -- a bare, unsubstituted naphthyridine isomer (15
        // chars, no substituents). RDKit: fully aromatic, all 10 atoms/bonds
        // (verified live against rdkit==2026.03.3). Under the DEFAULT engine
        // (`assign_aromaticity`), chematic confirms only 6/10 (the
        // pyridine-type ring) -- see
        // `test_default_engine_unaffected_by_fused_diazine_fix` below: the
        // benzo ring's Pass 1 evaluation wrongly zeroes out BOTH its fusion
        // carbons via `CarbonExocyclicHeteroatomDouble` (each fusion
        // carbon's own Kekule double bond points into the OTHER (pyridine)
        // ring, toward a nitrogen there -- from the benzo ring's own,
        // single-ring-only perspective using only `ring_atom_set`, that bond
        // looks exactly like a genuine exocyclic C=O/C=N substituent
        // (tropone's shape), which is what that rule is actually meant to
        // catch). Landing on EXACTLY pi=4 classifies the ring `Antiaromatic`,
        // which Pass 2 never retries ("definitive, do not retry").
        //
        // Fixed under the OPT-IN `assign_aromaticity_authoritative_experimental`
        // engine by making the rule bond-level (`all_ring_bonds`, built once
        // from every SSSR/augmented ring): a double bond whose far atom sits
        // on a DIFFERENT ring is a ring-fusion bond, not a substituent, so it
        // no longer zeroes the atom -- both fusion carbons now fall through
        // to the ordinary sp2 default (1π each), the benzo ring lands on
        // pi=6 (Aromatic) directly in Pass 1, and Pass 2 promotes the
        // pyridine-type ring via `AlreadyAromaticContext` as before.
        // Confirmed Kekule-choice-dependent, not shape-dependent: plain
        // quinoxaline and quinazoline (`c1ccc2nccnc2c1`, `c1ccc2ncncc2c1`)
        // never reproduced this in the first place (only ONE fusion carbon
        // was affected for those, landing on the retryable odd pi=5
        // NonAromatic case). This molecule was constructed as a minimal
        // repro for the dominant pattern seen in 33/84 corpus regressions
        // K2b's demotion fix surfaced (fused quinazoline/quinoxaline/
        // purine-shaped bicyclics with an N-substituent elsewhere in the
        // molecule); it is not itself one of the 84 (it is unsubstituted).
        let mol = mol_kekulized("c1cnc2ccccc2n1");
        let model = assign_aromaticity_authoritative_experimental(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "all 10 atoms should be aromatic under the opt-in engine, matching RDKit"
        );
        assert!(
            mol.atoms().all(|(idx, _)| model.is_atom_aromatic(idx)),
            "every atom should be aromatic"
        );
    }

    #[test]
    fn test_default_engine_unaffected_by_fused_diazine_fix() {
        // Same molecule as above, through the DEFAULT engine
        // (`assign_aromaticity`) -- must stay exactly as it was before the
        // K2b fused-diazine follow-up fix existed (6/10, still wrong vs
        // RDKit), confirming `apply_aromaticity`/`apply_aromaticity_ex`
        // remain byte-identical to pre-K2b behavior per the coordinator
        // decision to ship this as opt-in only.
        let mol = mol_kekulized("c1cnc2ccccc2n1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "default engine must stay unaffected: only the pyridine-type ring \
             (6/10 atoms) confirmed, matching pre-K2b behavior"
        );
    }

    /// A handful of the 33-molecule `fused_diazine_quinazoline_quinoxaline_purine`
    /// corpus cluster (K2b's own diagnosis; see the PR description), pinned
    /// as permanent regression tests against the OPT-IN engine now that this
    /// fix resolves them there. Not exhaustive -- the fixed corpus-vs-RDKit
    /// comparison (`scripts/aromaticity_atom_parity.py` equivalent run
    /// against `scripts/descriptor_census_corpus.smi`) is the authoritative
    /// check; these are a stable, minimal sample.
    #[test]
    fn test_authoritative_experimental_fused_diazine_cluster_sample_matches_rdkit() {
        // (smiles, expected RDKit-aromatic atom count, all-aromatic?)
        let cases: &[(&str, usize)] = &[
            ("COc1cccc2nc(N3CCNCC3)cnc12", 10),
            ("Fc1cccc2nc(N3CCNCC3)cnc12", 10),
            ("Clc1cccc2nc(N3CCNCC3)cnc12", 10),
            ("CN1CCN(c2cnc3cc(Cl)ccc3n2)CC1", 10),
            ("Clc1cc2ncc(N3CCNCC3)nc2cc1Cl", 10),
            ("O=C(O)C1CN(c2cnc3ccccc3n2)CCN1", 10),
        ];
        for (smi, expected) in cases {
            let mol = mol_kekulized(smi);
            let model = assign_aromaticity_authoritative_experimental(&mol);
            assert_eq!(
                model.aromatic_atom_count(),
                *expected,
                "{smi}: expected {expected} aromatic atoms (the fused \
                 quinoxaline/naphthyridine core) under the opt-in engine, matching RDKit"
            );
        }
    }

    #[test]
    fn test_known_gap_azulene_nonalternant_odd_odd_split() {
        // c1ccc2cccc-2cc1 -- azulene itself (already the canonical example
        // in this codebase and in docs/rfcs/aromaticity_a1_rfc.md). RDKit: fully
        // aromatic, all 10 atoms, 9/10 bonds (the explicit fusion bond the
        // SMILES itself writes non-aromatic, `-2`, stays a formal single
        // bond even in RDKit's own answer). chematic: 0/10 -- both the
        // 5-ring and 7-ring independently get an ODD pi count (5 and 7) in
        // Pass 1, so neither is Aromatic nor Antiaromatic (both
        // `NonAromatic`), and Pass 2 never seeds because seeding requires
        // an ALREADY-aromatic adjacent ring, which neither ring is able to
        // become on its own. The default model now applies a deliberately
        // narrow all-carbon odd/odd fused-envelope fallback for this case.
        let mol = mol_kekulized("c1ccc2cccc-2cc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "default Hückel's bounded fused-envelope fallback recognizes azulene"
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
        //
        // Regressed by the Horton SSSR rewrite (confirmed passing on the old
        // single-spanning-tree find_sssr, failing only after Horton). Root
        // cause, empirically confirmed via debug dump: Horton's correct,
        // minimal SSSR is exactly the 5-ring + 7-ring (matches RDKit). Each
        // evaluated standalone has an ODD pi-electron count (5-ring: 5pi,
        // 7-ring: 7pi — every ring atom contributes 1pi via a double bond,
        // whether the double bond is endo- or exocyclic-to-a-carbon), so
        // neither passes Pass 1 and neither can seed Pass 2's
        // aromatic_context bootstrap. Azulene's aromaticity is a genuinely
        // non-alternant, whole-perimeter (10-atom, 10pi) delocalized system
        // — it needs the full-ring-system envelope as a Hückel candidate,
        // which `augmented_ring_set` deliberately excludes (its docstring
        // names naphthalene's spurious 10-ring as the exact case to avoid).
        // The old, non-minimal SSSR happened to hand a large fundamental
        // cycle straight to Pass 1 that included the whole perimeter,
        // papering over this gap by coincidence.
        //
        // The bounded all-carbon odd/odd fused-envelope fallback now handles
        // this case without changing the broader default model.
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
    // RdkitLike mode: P/Se/Te heteroaromatics
    // =========================================================================

    #[test]
    fn test_phosphole_rdkit_aromatic() {
        // c1cc[pH]c1 — P donates its lone pair in the RDKit-compatible mode.
        let mol = mol_aromatic("c1cc[pH]c1");
        let m = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m.aromatic_atom_count(),
            5,
            "phosphole: all 5 atoms aromatic in RdkitLike"
        );
    }

    #[test]
    fn test_azulene_rdkit_like_uses_whole_perimeter() {
        // The strict per-ring Hückel pass sees azulene as an odd/odd fused
        // split. RDKit evaluates the connected 10π perimeter instead.
        let mol = mol_kekulized("C1=CC2=CC=CC=CC2=C1");
        let m = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m.aromatic_atom_count(),
            10,
            "azulene: whole perimeter must be aromatic in RdkitLike"
        );
    }

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
    // Re-measured after the Horton SSSR rewrite landed (find_sssr is now
    // minimal and deterministic, 0% self-instability on the 5000-molecule
    // corpus): all 32 counts below are UNCHANGED under the DEFAULT engine.
    // Zero free recoveries there.
    //
    // These 32 molecules share one root cause: a "fake bridgehead" N (same
    // local shape as a genuine bridgehead or N-substituted azole) feeds a
    // central ring that only closes via the `aromatic_context` bypass reusing
    // an unrelated ring's atoms. Fixing this requires removing the bypass in
    // favor of proper ring-system candidate enumeration (see project plan/
    // issue tracker).
    //
    // RESOLVED, but only under the OPT-IN `assign_aromaticity_authoritative_experimental`
    // engine (K2b fused-diazine follow-up fix; see
    // `test_authoritative_experimental_fixes_bridgehead_n_false_positives`
    // below): all 32 of these benzo-fused bridgehead-N tricyclics
    // (`...C3=NCCCN23`-shaped) ALSO have a fusion carbon whose own Kekule
    // double bond points into the adjacent ring at a heteroatom -- the exact
    // same misclassification the fused-diazine fix targets, just in a
    // three-ring rather than two-ring shape. Spot-checked live against
    // rdkit==2026.03.3 for 4 of the 32 (the shortest, a 15/12 case, and two
    // of the 28/24 cases), all matching. The originally-suspected root cause
    // above (the `aromatic_context`/`AlreadyAromaticContext` bypass) was
    // evidently either wrong or not the operative mechanism for this
    // specific molecule class -- not re-investigated further, since the fix
    // that resolved it was general (scoped to the fused-diazine cluster) and
    // not bridgehead-N-specific. This is opt-in only: the DEFAULT engine
    // (`assign_aromaticity`) is unaffected and still shows the original
    // `expected_wrong` counts below (see the coordinator decision requiring
    // `apply_aromaticity`/`apply_aromaticity_ex` to stay byte-identical to
    // pre-K2b behavior).
    // (kekulized SMILES, current chematic aromatic_atom_count() under the
    // default engine, RDKit's correct count).
    // Named at module level (not a local in the test below) so
    // Aromaticity-A1-0's corpus tests, further down this module, can reuse
    // the identical pinned data instead of re-deriving a copy that could
    // silently drift out of sync with it.
    const KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES: &[(&str, usize, usize)] = &[
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

    #[test]
    fn test_known_regressions_from_bridgehead_n_fix() {
        for (smi, expected_wrong, rdkit_correct) in KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES {
            let mol = mol_kekulized(smi);
            let model = assign_aromaticity(&mol);
            assert_eq!(
                model.aromatic_atom_count(),
                *expected_wrong,
                "{smi}: expected current (wrong) count {expected_wrong} under the default \
                 engine (RDKit correct: {rdkit_correct})"
            );
        }
    }

    #[test]
    fn test_authoritative_experimental_fixes_bridgehead_n_false_positives() {
        // Beneficial, unattempted side effect of the K2b fused-diazine
        // follow-up fix, now reachable only via the opt-in engine -- see
        // this const's preceding doc comment.
        for (smi, _expected_wrong, rdkit_correct) in KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES {
            let mol = mol_kekulized(smi);
            let model = assign_aromaticity_authoritative_experimental(&mol);
            assert_eq!(
                model.aromatic_atom_count(),
                *rdkit_correct,
                "{smi}: expected {rdkit_correct} aromatic atoms under the opt-in \
                 authoritative-experimental engine, matching RDKit"
            );
        }
    }

    // ── Known order-dependence: same molecule, different Kekulized traversal ─
    //
    // Originally found because these 3 molecules passed with RDKit's
    // canonical Kekulized SMILES but failed with at least one other valid
    // Kekulized ordering of the identical structure -- confirmed via
    // atom-map-number alignment (no substructure matching). Root cause was
    // NOT Pass 1/Pass 2 (verified order-invariant by construction) -- it was
    // `find_sssr` itself, non-deterministic and non-minimal.
    //
    // Re-measured after the Horton SSSR rewrite (find_sssr is now
    // deterministic and minimal, 0% self-instability on the 5000-molecule
    // corpus): the 3 pinned failing-traversal counts below are UNCHANGED.
    // The original order-dependence *mechanism* (find_sssr picking a
    // different non-minimal ring depending on traversal) is resolved -- but
    // these 3 specific SMILES still disagree with RDKit's count, so at least
    // one more bug (likely `aromatic_context`, same as the 32-molecule
    // corpus above) also affects this molecule class. Not re-diagnosed here;
    // a fresh worst-of-N run against the full corpus would confirm whether
    // order-dependence itself (canonical vs. this pinned variant disagreeing
    // with each other) is now fully gone, separate from RDKit agreement.
    //
    // The K2b fused-diazine follow-up fix (`assign_aromaticity_authoritative_experimental`)
    // does shift 2 of these 3 counts (16->12) when run through the OPT-IN
    // engine -- confirmed unrelated to and not fixing this bucket (still
    // wrong, by a different amount, a separate multi-causal bug). This test
    // asserts the DEFAULT (`assign_aromaticity`) engine only, which is
    // unaffected by that opt-in fix, so the pinned values below stay as
    // originally measured.
    // Named at module level for the same reason as
    // `KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES` above -- Aromaticity-A1-0's corpus
    // tests reuse this exact pinned data instead of a second copy.
    const KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES: &[(&str, usize, usize)] = &[
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

    #[test]
    fn test_known_order_dependent_regressions() {
        for (smi, expected_wrong, rdkit_correct) in KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES {
            let mol = mol_kekulized(smi);
            let model = assign_aromaticity(&mol);
            assert_eq!(
                model.aromatic_atom_count(),
                *expected_wrong,
                "{smi}: expected current (wrong) count {expected_wrong} (RDKit correct: {rdkit_correct})"
            );
        }
    }

    // ── Aromaticity-A1-0: anti-drift guard for `trace_ring_pi_electrons` ────
    //
    // `trace_ring_pi_electrons` is a deliberately separate implementation
    // from `ring_pi_electrons` (see the doc comment above it) so it can
    // report *why* each atom scored what it did. That separateness is a
    // drift risk: nothing stops the two from silently diverging as either
    // one is edited. This test is the guard -- for every ring in every
    // molecule of the known false-positive/false-negative/negative-control
    // corpus (the same molecules `docs/rfcs/aromaticity_a1_rfc.md`'s diagnostic
    // corpus uses), both functions must agree exactly, in both an empty
    // context (Pass-1-equivalent) and the model's final converged context
    // (an upper-bound Pass-2-equivalent). This does not assert anything
    // about correctness vs RDKit -- only that the trace and the real engine
    // never disagree with each other.
    #[test]
    fn trace_matches_ring_pi_electrons_on_corpus() {
        let smiles: Vec<&str> = KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES
            .iter()
            .map(|(smi, _, _)| *smi)
            .chain(
                KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES
                    .iter()
                    .map(|(smi, _, _)| *smi),
            )
            .chain([
                "C1=CC2=CC=CC=CC2=C1",    // azulene (Kekulized) -- known false negative
                "c1cnc2[nH]cnc2n1",       // purine -- known false negative
                "C1=Cc2ccccc2C2=NCCCN12", // PR #86 minimal false-positive reproducer
                "C1=Cc2ccccc2C2=CCCC12",  // negative control: no bridgehead N
                "C1=Cc2ccccc2C2=CCNC12",  // negative control: N not at bridgehead
                "C1Cc2ccccc2C2=NCCCN12",  // negative control: bridgehead N, no exocyclic C=C
                "c1ccc2[nH]ccc2c1",       // indole -- must stay correct
                "c1ccc2ncccc2c1",         // quinoline -- must stay correct
                "c1ccc2ccccc2c1",         // naphthalene -- must stay correct
            ])
            .collect();

        for algo in [
            AromaticityAlgorithm::Huckel,
            AromaticityAlgorithm::RdkitLike,
        ] {
            for smi in &smiles {
                let mol = mol_kekulized(smi);
                let model = assign_aromaticity_ex(&mol, algo);
                let final_context: FxHashSet<AtomIdx> = mol
                    .atoms()
                    .map(|(idx, _)| idx)
                    .filter(|&idx| model.is_atom_aromatic(idx))
                    .collect();

                let sssr = find_sssr(&mol);
                let rings = augmented_ring_set(&mol, sssr.rings());
                let empty_context: FxHashSet<AtomIdx> = FxHashSet::default();
                let all_ring_bonds: FxHashSet<BondIdx> =
                    rings.iter().flat_map(|r| ring_bond_set(&mol, r)).collect();

                for ring in &rings {
                    for ctx in [&empty_context, &final_context] {
                        let expected = ring_pi_electrons(&mol, ring, ctx, algo, &all_ring_bonds);
                        let traced =
                            trace_ring_pi_electrons(&mol, ring, ctx, algo, &all_ring_bonds);
                        assert_eq!(
                            traced.total,
                            expected,
                            "{smi} (algo={algo:?}, ring={ring:?}, ctx_len={}): \
                             trace_ring_pi_electrons diverged from ring_pi_electrons",
                            ctx.len()
                        );
                        // Cross-check the per-atom eligibility bookkeeping too.
                        for a in &traced.atoms {
                            assert_eq!(
                                a.contribution.is_some(),
                                a.reason.is_eligible(),
                                "{smi}: atom {:?} contribution/reason eligibility mismatch",
                                a.atom_idx
                            );
                        }
                    }
                }
            }
        }
    }

    // ── Aromaticity-A1-0: false-positive/false-negative polarity sanity ────
    //
    // These are cheap, structural sanity checks that the corpus buckets are
    // labeled the direction they claim -- not a re-measurement of the full
    // corpus (that's `aromaticity_a1_0_report` + the Python RDKit join, see
    // `docs/rfcs/aromaticity_a1_rfc.md`). Catches an accidental swap or a stale
    // pinned count silently going the other way.
    #[test]
    fn false_positive_corpus_over_counts_vs_rdkit() {
        for (smi, expected_wrong, rdkit_correct) in KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES {
            assert!(
                expected_wrong > rdkit_correct,
                "{smi}: false-positive bucket entry should over-count \
                 (chematic={expected_wrong} should be > rdkit={rdkit_correct})"
            );
        }
    }

    #[test]
    fn false_negative_corpus_under_counts_vs_rdkit() {
        for (smi, expected_wrong, rdkit_correct) in KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES {
            assert!(
                expected_wrong < rdkit_correct,
                "{smi}: false-negative bucket entry should under-count \
                 (chematic={expected_wrong} should be < rdkit={rdkit_correct})"
            );
        }
    }

    // ── Aromaticity-A1-1a: exhaustive_aromaticity_oracle pinned cases ──────
    //
    // The oracle is a discovery tool, not a correct-answer generator: its
    // candidates are built from the SAME per-atom local rules
    // (`evaluate_atom_pi_contribution`) that are wrong for the false-positive
    // family, so it can't independently arbitrate that family. This test
    // pins what the oracle DOES get right (RDKit-atom-index-verified, not
    // guessed) after two real fixes made during this milestone:
    //
    // 1. Connectivity: `build_conjugated_components`'s conjugation graph
    //    originally only bridged single bonds via a `LonePairDonor` endpoint,
    //    leaving azulene's all-carbon alternating perimeter as 5 disconnected
    //    2-atom pairs (oracle returned an empty set). Fixed: any bond between
    //    two independently-eligible atoms connects (ordinary carbon-carbon
    //    single-bond conjugation, ordinary organic chemistry).
    // 2. Home-ring evaluation: evaluating a multi-ring candidate's electron
    //    sum against its own *flattened* atom set broke the N
    //    bridgehead/substituted-azole rule for any TRUE bridgehead (every
    //    bond looks "in-family" once the family itself is the context) --
    //    indolizine's own bridgehead N came out `Ineligible`, an oracle bug,
    //    not a chematic bug. Fixed via `evaluate_atom_via_home_ring`.
    //
    // Both fixes were originally confirmed correct AND confirmed NOT to
    // silently "fix" the false-positive family by accident.
    //
    // UPDATE (K2b fused-diazine follow-up fix): both the false-positive
    // reproducer AND purine are now RDKit-exact too, as a side effect of the
    // same general `CarbonExocyclicHeteroatomDouble` ring-fusion fix
    // (`evaluate_atom_pi_contribution_inner` mirrors `ring_pi_electrons`'s
    // rule exactly -- see its doc comment). The false-positive reproducer's
    // own fusion carbon (whose double bond points into the bridgehead-N
    // ring's own nitrogen) no longer gets wrongly zeroed, so the oracle
    // stops over-aromatizing into the bridgehead ring and correctly confirms
    // only the plain benzo ring. Purine's 5-ring fusion carbons no longer
    // get wrongly zeroed by the same rule either, so the oracle now confirms
    // all 9 atoms, matching RDKit -- resolving the open finding below.
    // Verified live against rdkit==2026.03.3 for both (not assumed from the
    // fix's general mechanism alone).
    #[test]
    fn exhaustive_oracle_pinned_cases() {
        let algo = AromaticityAlgorithm::RdkitLike;

        // (name, smiles, expected oracle-aromatic atom indices, sorted)
        let matches_rdkit: &[(&str, &str, &[u32])] = &[
            (
                "azulene",
                "C1=CC2=CC=CC=CC2=C1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "naphthalene",
                "c1ccc2ccccc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "anthracene",
                "c1ccc2cc3ccccc3cc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
            ),
            ("indole", "c1ccc2[nH]ccc2c1", &[0, 1, 2, 3, 4, 5, 6, 7, 8]),
            (
                "quinoline",
                "c1ccc2ncccc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "indolizine (bridgehead N, both rings valid)",
                "c1ccn2ccccc12",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            ),
            ("tropone", "O=c1cccccc1", &[1, 2, 3, 4, 5, 6, 7]),
            ("2-pyridone", "O=c1cccc[nH]1", &[1, 2, 3, 4, 5, 6]),
        ];
        for (name, smi, expected) in matches_rdkit {
            let mol = mol_kekulized(smi);
            let (atoms, _bonds) = exhaustive_aromaticity_oracle(&mol, algo);
            let mut got: Vec<u32> = atoms.iter().map(|a| a.0).collect();
            got.sort();
            assert_eq!(&got, expected, "{name} ({smi}): oracle should match RDKit");
        }

        // Now RDKit-exact -- see this test's doc comment (K2b fused-diazine
        // follow-up fix). RDKit: only the plain benzo ring (6 atoms) is
        // aromatic; the bridgehead-N ring is not (verified live).
        let (fp_atoms, _) =
            exhaustive_aromaticity_oracle(&mol_kekulized("C1=Cc2ccccc2C2=NCCCN12"), algo);
        let mut fp_got: Vec<u32> = fp_atoms.iter().map(|a| a.0).collect();
        fp_got.sort();
        assert_eq!(
            fp_got,
            vec![2, 3, 4, 5, 6, 7],
            "false-positive reproducer: oracle now matches RDKit exactly"
        );

        // Now RDKit-exact -- see this test's doc comment (K2b fused-diazine
        // follow-up fix). RDKit: all 9 atoms aromatic (verified live).
        let (purine_atoms, _) =
            exhaustive_aromaticity_oracle(&mol_kekulized("c1cnc2[nH]cnc2n1"), algo);
        let mut purine_got: Vec<u32> = purine_atoms.iter().map(|a| a.0).collect();
        purine_got.sort();
        assert_eq!(
            purine_got,
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
            "purine: oracle now matches RDKit exactly (all 9 atoms)"
        );
    }
}
