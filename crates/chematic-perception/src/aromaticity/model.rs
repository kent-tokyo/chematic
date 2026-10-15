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
    let sssr_count = ring_set.rings().len();
    let rings = augmented_ring_set(mol, ring_set.rings());
    let all_ring_bonds = aromaticity_ring_bonds(mol, &rings, ring_fusion_aware);
    let mut state = AromaticityPassState::new(rings.len());
    state.first_pass(mol, &rings, algo, &all_ring_bonds);
    state.fused_pass(mol, &rings, algo, &all_ring_bonds);

    // A strict per-ring pass cannot seed azulene's 5+7 fused system because
    // both constituent rings have an odd local count. Apply only the narrow,
    // fail-closed all-carbon odd/odd envelope rule here; this does not route
    // the default model through the broader RdkitLike implementation.
    if algo == AromaticityAlgorithm::Huckel {
        apply_huckel_nonalternant_fused_fallback(
            mol,
            &rings,
            &mut state.aromatic_atoms,
            &mut state.aromatic_bonds,
        );
    }
    AromaticityModel {
        aromatic_atoms: state.aromatic_atoms,
        aromatic_bonds: state.aromatic_bonds,
        antiaromatic_rings: state.antiaromatic_rings,
        ring_classifications: public_ring_classifications(
            &rings,
            sssr_count,
            &state.classifications,
        ),
    }
}

fn aromaticity_ring_bonds(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    ring_fusion_aware: bool,
) -> FxHashSet<BondIdx> {
    if ring_fusion_aware {
        rings
            .iter()
            .flat_map(|ring| ring_bond_set(mol, ring))
            .collect()
    } else {
        FxHashSet::default()
    }
}

struct AromaticityPassState {
    aromatic_atoms: FxHashSet<AtomIdx>,
    aromatic_bonds: FxHashSet<BondIdx>,
    antiaromatic_rings: Vec<Vec<AtomIdx>>,
    classifications: Vec<Option<(RingAromaticity, u32)>>,
    pending: Vec<usize>,
}

impl AromaticityPassState {
    fn new(ring_count: usize) -> Self {
        Self {
            aromatic_atoms: FxHashSet::default(),
            aromatic_bonds: FxHashSet::default(),
            antiaromatic_rings: Vec::new(),
            classifications: vec![None; ring_count],
            pending: Vec::new(),
        }
    }

    fn first_pass(
        &mut self,
        mol: &Molecule,
        rings: &[Vec<AtomIdx>],
        algo: AromaticityAlgorithm,
        all_ring_bonds: &FxHashSet<BondIdx>,
    ) {
        let empty_context = FxHashSet::default();
        for (index, ring) in rings.iter().enumerate() {
            let Some(electrons) =
                ring_pi_electrons(mol, ring, &empty_context, algo, all_ring_bonds)
            else {
                self.pending.push(index);
                continue;
            };
            let classification = classify_ring_aromaticity(electrons);
            self.classifications[index] = Some(classification);
            match classification.0 {
                RingAromaticity::Aromatic => mark_ring_aromatic(
                    mol,
                    ring,
                    &mut self.aromatic_atoms,
                    &mut self.aromatic_bonds,
                ),
                RingAromaticity::Antiaromatic => self.antiaromatic_rings.push(ring.clone()),
                RingAromaticity::NonAromatic => self.pending.push(index),
            }
        }
    }

    fn fused_pass(
        &mut self,
        mol: &Molecule,
        rings: &[Vec<AtomIdx>],
        algo: AromaticityAlgorithm,
        all_ring_bonds: &FxHashSet<BondIdx>,
    ) {
        while !self.pending.is_empty() {
            let (pending, added) = self.fused_iteration(mol, rings, algo, all_ring_bonds);
            self.pending = pending;
            let all_set = rings
                .iter()
                .flatten()
                .all(|atom| self.aromatic_atoms.contains(atom));
            if !added || all_set {
                break;
            }
        }
    }

    fn fused_iteration(
        &mut self,
        mol: &Molecule,
        rings: &[Vec<AtomIdx>],
        algo: AromaticityAlgorithm,
        all_ring_bonds: &FxHashSet<BondIdx>,
    ) -> (Vec<usize>, bool) {
        let mut pending = Vec::new();
        let mut added = false;
        for index in std::mem::take(&mut self.pending) {
            let ring = &rings[index];
            if !ring.iter().any(|atom| self.aromatic_atoms.contains(atom)) {
                pending.push(index);
                continue;
            }
            let Some(electrons) =
                ring_pi_electrons(mol, ring, &self.aromatic_atoms, algo, all_ring_bonds)
            else {
                pending.push(index);
                continue;
            };
            let classification = classify_ring_aromaticity(electrons);
            self.classifications[index] = Some(classification);
            if classification.0 == RingAromaticity::Aromatic {
                mark_ring_aromatic(
                    mol,
                    ring,
                    &mut self.aromatic_atoms,
                    &mut self.aromatic_bonds,
                );
                added = true;
            }
        }
        (pending, added)
    }
}

fn public_ring_classifications(
    rings: &[Vec<AtomIdx>],
    sssr_count: usize,
    classifications: &[Option<(RingAromaticity, u32)>],
) -> Vec<(Vec<AtomIdx>, RingAromaticity, u32)> {
    rings
        .iter()
        .take(sssr_count)
        .enumerate()
        .filter_map(|(index, ring)| {
            classifications[index]
                .map(|(classification, count)| (ring.clone(), classification, count))
        })
        .collect()
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
    // Lowercase aromatic input is normally preserved for compatibility, but
    // a neutral group-16 atom with more than two coordinated neighbours
    // cannot be a lone-pair donor.  Rebuild that narrow invalid case from a
    // Kekule form so stale parser flags do not survive the corrected #767
    // donor verdict.  If the input cannot be Kekulized, retain the historical
    // infallible fallback rather than inventing a partial bond assignment.
    let has_invalid_explicit_chalcogen = mol.atoms().any(|(idx, atom)| {
        atom.aromatic
            && atom.charge == 0
            && matches!(atom.element.atomic_number(), 8 | 16 | 34 | 52)
            && mol.degree(idx) + implicit_hcount(mol, idx) as usize > 2
    });
    if has_invalid_explicit_chalcogen && let Ok(kekule) = chematic_core::kekulize(mol) {
        let kekulized = chematic_core::apply_kekule(mol, &kekule);
        let model = assign_aromaticity_ex(&kekulized, algo);
        return build_molecule_from_model_authoritative(&kekulized, &model);
    }

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
    let mut pre_h: Vec<Option<u8>> = mol
        .atoms()
        .map(|(idx, atom)| {
            if atom.hydrogen_count.is_some() {
                None // already explicit; nothing to preserve
            } else {
                Some(implicit_hcount(mol, idx))
            }
        })
        .collect();

    let mut builder = MoleculeBuilder::with_capacity(mol.atom_count(), mol.bond_count());
    for (idx, atom) in mol.atoms() {
        let mut a = atom.clone();
        a.aromatic = atom_aromatic.contains(&idx);
        builder.add_atom(a);
    }
    let mut stash = Vec::new();
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
            stash.push((new_bidx, bond.order, bond.atom1));
        }
    }
    // Atoms/bonds above are re-added in `mol`'s own enumeration order with
    // none skipped, so indices line up 1:1 — safe to copy side-channel
    // metadata wholesale. (This rebuild previously dropped stereo_groups and
    // stereo_neighbor_order silently; closing that here too.)
    builder.copy_r_groups_from(mol);
    builder.copy_atom_tags_from(mol);
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    // After the wholesale copy, which would otherwise erase them.
    for (bidx, dir, anchor) in stash {
        builder.set_bond_direction(bidx, dir);
        builder.set_bond_direction_anchor(bidx, anchor);
    }
    let normalized = builder.build();

    // Compare the pre-normalization implicit H against what the same
    // (already-tested, unmodified) `implicit_hcount` computes on the
    // normalized bonds; only atoms where normalization actually changed the
    // answer get an explicit H frozen in. Benzene CH and pyridine-type N
    // (heuristic already agrees) are left untouched -- no spurious bracket
    // notation for atoms that didn't need it.
    // Reuse the pre-normalization storage as an indexed patch map. Clearing
    // unchanged entries avoids a second allocation and per-atom linear scans.
    let mut needs_patch = false;
    for (index, h) in pre_h.iter_mut().enumerate() {
        *h = h.filter(|&pre| pre != implicit_hcount(&normalized, AtomIdx(index as u32)));
        needs_patch |= h.is_some();
    }
    if !needs_patch {
        return normalized;
    }

    let mut patched =
        MoleculeBuilder::with_capacity(normalized.atom_count(), normalized.bond_count());
    for (idx, atom) in normalized.atoms() {
        let mut a = atom.clone();
        if let Some(h) = pre_h[idx.0 as usize] {
            a.hydrogen_count = Some(h);
        }
        patched.add_atom(a);
    }
    for (_bond_idx, bond) in normalized.bonds() {
        let _ = patched.add_bond(bond.atom1, bond.atom2, bond.order);
    }
    patched.copy_stereo_groups_from(&normalized);
    patched.copy_stereo_from(&normalized);
    patched.copy_atom_tags_from(&normalized);
    patched.copy_bond_directions_from(&normalized);
    patched.build()
}

#[cfg(test)]
#[path = "normalization_storage_tests.rs"]
mod normalization_storage_tests;
