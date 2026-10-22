//! Auditable standardization pipeline and stage diagnostics.

use std::collections::HashSet;

use chematic_core::{
    AtomIdx, BondIdx, Molecule, MoleculeBuilder, valence_inferred_hcount, validate_valence,
};

use crate::{hash::mol_hash, hydrogen::remove_hydrogens, tautomer::canonical_tautomer};

use super::{clone_molecule, is_metal, largest_fragment, neutralize_charges, normalize_zwitterion};

/// One transformation stage in the standardization pipeline.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum StandardizationStep {
    /// Remove metal coordination bonds before charge and fragment processing.
    DisconnectMetals,
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
            Self::DisconnectMetals => "disconnect_metals",
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
    pub(super) fn new(code: &str, message: String) -> Self {
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
    pub(crate) fn from_mol(mol: &Molecule) -> Self {
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

        // Disconnect metals early (remove dative/coordinate bonds). Keep this
        // implicit chemistry step in the report: otherwise the input snapshot
        // and the first reported stage describe different molecules.
        let has_metals = current.atoms().any(|(_, a)| is_metal(a.element));
        current = self.apply_stage(
            current,
            StandardizationStep::DisconnectMetals,
            has_metals,
            disconnect_metals,
            &mut steps,
            &mut warnings,
        );

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

    /// The molecule [`Self::run`] returns, without the audit report.
    ///
    /// `run` snapshots the molecule (a canonical-SMILES hash) before and after
    /// every stage and validates valences for warnings; none of that feeds
    /// back into the structure. This applies the identical stage sequence —
    /// including the `clone_molecule` pass-through for disabled stages — and
    /// skips the bookkeeping.
    pub fn run_output_only(&self, mol: &Molecule) -> Molecule {
        let stage = |current: Molecule, enabled: bool, f: fn(&Molecule) -> Molecule| {
            if enabled {
                f(&current)
            } else {
                clone_molecule(&current)
            }
        };
        let mut current = clone_molecule(mol);
        let has_metals = current.atoms().any(|(_, a)| is_metal(a.element));
        current = stage(current, has_metals, disconnect_metals);
        current = stage(current, self.options.neutralize_charges, neutralize_charges);
        current = stage(
            current,
            self.options.largest_fragment_only,
            largest_fragment,
        );
        let zwitterion_enabled = self.options.zwitterion_handling == ZwitterionHandling::Normalize;
        current = stage(current, zwitterion_enabled, normalize_zwitterion);
        current = stage(current, self.options.remove_explicit_h, remove_hydrogens);
        stage(current, self.options.canonical_tautomer, canonical_tautomer)
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
///
/// For every non-metal atom losing a bond this way, an explicitly-stored
/// `hydrogen_count` is recomputed by valence inference from the *surviving*
/// bonds -- a dative M-O/M-N bond is commonly written with a formal charge on
/// the non-metal atom that exactly balances the bond (e.g. `[O+]`
/// single-bonded to the metal, satisfying O+'s valence-3), and leaving the
/// stale, now-too-low H count in place would let the very next pipeline
/// stage, `neutralize_charges`, see `h == 0` (its guard is `h > 0`), skip
/// neutralizing it, and leave a dangling formal charge with no bond left to
/// justify it -- a real, confirmed idempotency bug (issue #403): the charge
/// only got neutralized on a *second* standardize pass, once a fresh parse of
/// the first pass's (incorrectly charged) output SMILES stored the H count
/// explicitly instead of inferring it. Recomputing here (not merely resetting
/// to `None` -- `neutralize_charges` reads the raw stored field, not
/// `implicit_hcount`, so an unresolved `None` would look identical to a
/// stored `Some(0)` to its guard) means `neutralize_charges` sees the true
/// post-disconnection valence on the very first pass and can correctly
/// neutralize (or correctly decide not to) right away, converging to the
/// same result every time.
fn disconnect_metals(mol: &Molecule) -> Molecule {
    // Atoms adjacent to a metal, on the *non-metal* side of the bond -- their
    // hydrogen_count needs recomputing once that bond is dropped.
    let mut affected: HashSet<AtomIdx> = HashSet::new();
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        let atom1_is_metal = is_metal(mol.atom(bond.atom1).element);
        let atom2_is_metal = is_metal(mol.atom(bond.atom2).element);
        if atom1_is_metal && !atom2_is_metal {
            affected.insert(bond.atom2);
        } else if atom2_is_metal && !atom1_is_metal {
            affected.insert(bond.atom1);
        }
    }

    // Pass 1: copy atoms unchanged, drop metal bonds -- builds the
    // post-disconnection *topology* so valence inference below has the
    // surviving bonds (not the original, metal-including ones) to work from.
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        builder.add_atom(mol.atom(AtomIdx(i as u32)).clone());
    }
    for i in 0..mol.bond_count() {
        let old_bidx = BondIdx(i as u32);
        let bond = mol.bond(old_bidx);
        let atom1_is_metal = is_metal(mol.atom(bond.atom1).element);
        let atom2_is_metal = is_metal(mol.atom(bond.atom2).element);
        if !atom1_is_metal
            && !atom2_is_metal
            && let Ok(new_bidx) = builder.add_bond(bond.atom1, bond.atom2, bond.order)
            && let Some(direction) = mol.bond_direction(old_bidx)
        {
            builder.set_bond_direction(new_bidx, direction);
            if let Some(anchor) = mol.bond_direction_anchor(old_bidx) {
                builder.set_bond_direction_anchor(new_bidx, anchor);
            }
        }
    }
    let disconnected = builder.build();

    // Pass 2: recompute the affected atoms' hydrogen_count against the
    // now-disconnected topology (charge unchanged; `neutralize_charges` runs
    // next and decides whether the charge itself should move).
    let mut builder = MoleculeBuilder::new();
    for i in 0..disconnected.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = disconnected.atom(idx).clone();
        if affected.contains(&idx) && atom.hydrogen_count.is_some() {
            atom.hydrogen_count = Some(valence_inferred_hcount(&disconnected, idx));
        }
        builder.add_atom(atom);
    }
    for i in 0..disconnected.bond_count() {
        let bond = disconnected.bond(BondIdx(i as u32));
        if let Ok(new_bidx) = builder.add_bond(bond.atom1, bond.atom2, bond.order)
            && let Some(direction) = disconnected.bond_direction(BondIdx(i as u32))
        {
            builder.set_bond_direction(new_bidx, direction);
            if let Some(anchor) = disconnected.bond_direction_anchor(BondIdx(i as u32)) {
                builder.set_bond_direction_anchor(new_bidx, anchor);
            }
        }
    }
    builder.copy_stereo_from(&disconnected);
    builder.copy_stereo_groups_from(&disconnected);
    builder.build()
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
///
/// With `canonical_tautomer` enabled (the default), a stereocenter's CIP
/// (R/S) label can legitimately change even though its real spatial
/// configuration never moves -- see [`crate::tautomer::canonical_tautomer`]'s
/// own doc comment (issue #402) before treating a pre- vs. post-standardize
/// CIP-label difference as a bug on its own.
pub fn standardize(mol: &Molecule, opts: &StandardizeOptions) -> Molecule {
    StandardizationPipeline::new(opts.clone()).run_output_only(mol)
}
