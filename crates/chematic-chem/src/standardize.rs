//! Molecular standardization routines.
//!
//! Provides utilities for cleaning up molecular representations:
//! - Selecting the largest connected fragment.
//! - Neutralizing simple formal charges.

#![forbid(unsafe_code)]

use std::collections::{HashMap, VecDeque};

use chematic_core::{AtomIdx, BondIdx, Element, Molecule, MoleculeBuilder, validate_valence};
use serde::{Deserialize, Serialize};

use crate::{hash::mol_hash, hydrogen::remove_hydrogens, tautomer::canonical_tautomer};

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

    components.sort_by(|a, b| b.len().cmp(&a.len()));
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

/// Return a new `Molecule` containing only the largest connected fragment.
///
/// If the molecule is empty, an empty `Molecule` is returned.
pub fn largest_fragment(mol: &Molecule) -> Molecule {
    if mol.atom_count() == 0 {
        return MoleculeBuilder::new().build();
    }

    let components = connected_components(mol);
    let largest = &components[0];

    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let mut builder = MoleculeBuilder::new();
    for &old_idx in largest {
        let new_idx = builder.add_atom(mol.atom(old_idx).clone());
        remap.insert(old_idx, new_idx);
    }
    copy_bonds(mol, &mut builder, &remap);
    builder.build()
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
            (Element::N, 1) | (Element::O, 1) => {
                if h > 0 {
                    modifications.insert(idx, (0, Some(h - 1)));
                }
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

/// One transformation stage in the standardization pipeline.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum StandardizationStep {
    /// Select the largest connected component.
    LargestFragment,
    /// Apply simple neutralization rules for common formal charges.
    NeutralizeCharges,
    /// Remove explicit hydrogen atoms.
    RemoveExplicitHydrogens,
    /// Canonicalize supported tautomer systems.
    CanonicalTautomer,
}

impl StandardizationStep {
    /// Stable machine-readable stage name.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::LargestFragment => "largest_fragment",
            Self::NeutralizeCharges => "neutralize_charges",
            Self::RemoveExplicitHydrogens => "remove_explicit_hydrogens",
            Self::CanonicalTautomer => "canonical_tautomer",
        }
    }
}

/// High-level status for a standardization run.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PipelineStatus {
    /// Pipeline completed and the output is structurally identical to the input.
    Unchanged,
    /// Pipeline completed and at least one enabled stage changed the molecule.
    Modified,
    /// Pipeline completed, but warnings indicate unsupported or suspicious input features.
    CompletedWithWarnings,
}

/// Warning emitted during standardization.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
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
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
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
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
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
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
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
}

impl Default for StandardizeOptions {
    fn default() -> Self {
        Self {
            canonical_tautomer: true,
            neutralize_charges: true,
            remove_explicit_h: true,
            largest_fragment_only: false,
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

        current = self.apply_stage(
            current,
            StandardizationStep::LargestFragment,
            self.options.largest_fragment_only,
            largest_fragment,
            &mut steps,
            &mut warnings,
        );
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
        let status = if !warnings.is_empty() {
            PipelineStatus::CompletedWithWarnings
        } else if input.hash == output.hash {
            PipelineStatus::Unchanged
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
    let mut metal_bonds = 0usize;
    for (idx, atom) in mol.atoms() {
        if !is_metal(atom.element) {
            continue;
        }
        metal_bonds += mol
            .neighbors(idx)
            .filter(|(nb, _)| !is_metal(mol.atom(*nb).element))
            .count();
    }
    if metal_bonds > 0 {
        warnings.push(StandardizationWarning::new(
            "metal_disconnection_not_applied",
            format!(
                "found {metal_bonds} metal-to-nonmetal bond(s); metal disconnection is not implemented yet"
            ),
        ));
    }
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
        });

        let (result, report) = pipeline.run(&mol);

        assert_eq!(result.atom_count(), 3);
        assert_eq!(report.status, PipelineStatus::Modified);
        assert!(report.changed());
        assert_eq!(report.steps.len(), 4);
        assert_eq!(report.steps[0].step, StandardizationStep::LargestFragment);
        assert!(report.steps[0].enabled);
        assert!(report.steps[0].changed);
        assert_eq!(report.steps[1].step, StandardizationStep::NeutralizeCharges);
        assert!(!report.steps[1].enabled);
    }

    #[test]
    fn pipeline_report_marks_unchanged_clean_molecule() {
        let mol = parse("CC").unwrap();
        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            canonical_tautomer: false,
            neutralize_charges: false,
            remove_explicit_h: false,
            largest_fragment_only: false,
        });

        let (_result, report) = pipeline.run(&mol);

        assert_eq!(report.status, PipelineStatus::Unchanged);
        assert!(!report.changed());
        assert!(report.warnings.is_empty());
        assert!(report.steps.iter().all(|s| !s.enabled && !s.changed));
    }

    #[test]
    fn pipeline_report_warns_about_metal_disconnection_gap() {
        let mol = parse("[Na]OC").unwrap();
        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            canonical_tautomer: false,
            neutralize_charges: false,
            remove_explicit_h: false,
            largest_fragment_only: false,
        });

        let (_result, report) = pipeline.run(&mol);

        assert_eq!(report.status, PipelineStatus::CompletedWithWarnings);
        assert!(
            report
                .warnings
                .iter()
                .any(|w| w.code == "metal_disconnection_not_applied")
        );
    }
}
