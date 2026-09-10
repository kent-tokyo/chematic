//! Evidence-scoped reaction stoichiometry diagnostics.
//!
//! This module checks explicit atom/isotope inventories and formal charges.
//! It deliberately does not claim chemical completeness, mechanism validity,
//! or product prediction.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::document::{ComponentRole, ReactionDocument, ReactionStep};

const EPSILON: f64 = 1e-9;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StoichiometryEvidenceScope {
    ExplicitAtomInventory,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ChemicalCompleteness {
    NotEvaluated,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StoichiometryStatus {
    Balanced,
    Unbalanced,
    UnderSpecified,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DiagnosticSeverity {
    Info,
    Warning,
    Error,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StoichiometryIssueCode {
    AtomInventoryImbalance,
    IsotopeInventoryImbalance,
    ChargeImbalance,
    MissingReactants,
    MissingProducts,
    AgentOnly,
    FractionalCoefficient,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StoichiometryDiagnostic {
    pub path: String,
    pub code: StoichiometryIssueCode,
    pub severity: DiagnosticSeverity,
    pub message: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AtomInventory {
    pub atomic_number: u8,
    pub isotope: Option<u16>,
    pub count: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ComponentEvidence {
    pub path: String,
    pub id: String,
    pub role: ComponentRole,
    pub coefficient: f64,
    pub atoms: Vec<AtomInventory>,
    pub formal_charge: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StepStoichiometryReport {
    pub step_id: String,
    pub status: StoichiometryStatus,
    pub components: Vec<ComponentEvidence>,
    pub reactant_atoms: Vec<AtomInventory>,
    pub product_atoms: Vec<AtomInventory>,
    pub reactant_charge: f64,
    pub product_charge: f64,
    pub diagnostics: Vec<StoichiometryDiagnostic>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StoichiometryReport {
    pub evidence_scope: StoichiometryEvidenceScope,
    pub chemical_completeness: ChemicalCompleteness,
    pub status: StoichiometryStatus,
    pub steps: Vec<StepStoichiometryReport>,
}

/// Input component for evidence analysis. Coefficients may be positive
/// integers or positive fractional values; they are never rounded.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StoichiometryComponent {
    pub id: String,
    pub role: ComponentRole,
    pub smiles: String,
    pub coefficient: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StoichiometryStep {
    pub id: String,
    pub components: Vec<StoichiometryComponent>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum StoichiometryError {
    InvalidCoefficient { path: String, coefficient: f64 },
    EmptyComponent { path: String },
    Parse { path: String, message: String },
}

impl core::fmt::Display for StoichiometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidCoefficient { path, coefficient } => {
                write!(f, "invalid coefficient at {path}: {coefficient}")
            }
            Self::EmptyComponent { path } => write!(f, "empty component at {path}"),
            Self::Parse { path, message } => {
                write!(f, "failed to parse component at {path}: {message}")
            }
        }
    }
}

impl std::error::Error for StoichiometryError {}

/// Analyze all document steps while preserving their stable IDs and paths.
pub fn analyze_reaction_document(
    document: &ReactionDocument,
) -> Result<StoichiometryReport, StoichiometryError> {
    document.validate().map_err(|e| StoichiometryError::Parse {
        path: "/document".into(),
        message: e.to_string(),
    })?;
    let mut steps = Vec::with_capacity(document.steps.len());
    for step in &document.steps {
        steps.push(analyze_reaction_step(step)?);
    }
    let status = if steps
        .iter()
        .any(|s| s.status == StoichiometryStatus::Unbalanced)
    {
        StoichiometryStatus::Unbalanced
    } else if steps
        .iter()
        .any(|s| s.status == StoichiometryStatus::UnderSpecified)
    {
        StoichiometryStatus::UnderSpecified
    } else {
        StoichiometryStatus::Balanced
    };
    Ok(StoichiometryReport {
        evidence_scope: StoichiometryEvidenceScope::ExplicitAtomInventory,
        chemical_completeness: ChemicalCompleteness::NotEvaluated,
        status,
        steps,
    })
}

/// Analyze one rich document step using integer component coefficients.
pub fn analyze_reaction_step(
    step: &ReactionStep,
) -> Result<StepStoichiometryReport, StoichiometryError> {
    analyze_components(&StoichiometryStep {
        id: step.id.clone(),
        components: step
            .components
            .iter()
            .map(|component| StoichiometryComponent {
                id: component.id.clone(),
                role: component.role,
                smiles: component.smiles.clone(),
                coefficient: component.coefficient as f64,
            })
            .collect(),
    })
}

/// Analyze one step with integer or fractional coefficients.
pub fn analyze_components(
    step: &StoichiometryStep,
) -> Result<StepStoichiometryReport, StoichiometryError> {
    let mut components = Vec::with_capacity(step.components.len());
    let mut reactants = BTreeMap::new();
    let mut products = BTreeMap::new();
    let mut reactant_charge = 0.0;
    let mut product_charge = 0.0;
    let mut has_reactants = false;
    let mut has_products = false;
    let mut has_agents = false;
    let mut diagnostics = Vec::new();
    for component in &step.components {
        let path = format!("/steps/{}/components/{}", step.id, component.id);
        if !component.coefficient.is_finite() || component.coefficient <= 0.0 {
            return Err(StoichiometryError::InvalidCoefficient {
                path,
                coefficient: component.coefficient,
            });
        }
        if component.smiles.is_empty() {
            return Err(StoichiometryError::EmptyComponent { path });
        }
        let molecule =
            chematic_smiles::parse(&component.smiles).map_err(|e| StoichiometryError::Parse {
                path: path.clone(),
                message: e.to_string(),
            })?;
        let (atoms, charge) = inventory(&molecule, component.coefficient);
        if component.coefficient.fract().abs() > EPSILON {
            diagnostics.push(StoichiometryDiagnostic {
                path: path.clone(),
                code: StoichiometryIssueCode::FractionalCoefficient,
                severity: DiagnosticSeverity::Info,
                message: "fractional coefficient was evaluated exactly as supplied".into(),
            });
        }
        match component.role {
            ComponentRole::Reactant => {
                has_reactants = true;
                merge(&mut reactants, &atoms);
                reactant_charge += charge;
            }
            ComponentRole::Product => {
                has_products = true;
                merge(&mut products, &atoms);
                product_charge += charge;
            }
            ComponentRole::Agent => {
                has_agents = true;
            }
        }
        components.push(ComponentEvidence {
            path,
            id: component.id.clone(),
            role: component.role,
            coefficient: component.coefficient,
            atoms: atoms
                .into_iter()
                .map(|((atomic_number, isotope), count)| AtomInventory {
                    atomic_number,
                    isotope,
                    count,
                })
                .collect(),
            formal_charge: charge,
        });
    }
    if !has_reactants {
        diagnostics.push(diag(
            "/steps",
            StoichiometryIssueCode::MissingReactants,
            "reaction has no reactant components",
            DiagnosticSeverity::Warning,
        ));
    }
    if !has_products {
        diagnostics.push(diag(
            "/steps",
            StoichiometryIssueCode::MissingProducts,
            "reaction has no product components",
            DiagnosticSeverity::Warning,
        ));
    }
    if has_agents && !has_reactants && !has_products {
        diagnostics.push(diag(
            "/steps",
            StoichiometryIssueCode::AgentOnly,
            "agents are present without reactant or product components",
            DiagnosticSeverity::Warning,
        ));
    }
    let reactant_atoms = to_inventory(reactants);
    let product_atoms = to_inventory(products);
    let reactant_elements = element_totals(&reactant_atoms);
    let product_elements = element_totals(&product_atoms);
    if reactant_elements != product_elements {
        diagnostics.push(diag(
            "/steps",
            StoichiometryIssueCode::AtomInventoryImbalance,
            "explicit atom inventories differ",
            DiagnosticSeverity::Error,
        ));
    }
    if reactant_elements == product_elements && reactant_atoms != product_atoms {
        diagnostics.push(diag(
            "/steps",
            StoichiometryIssueCode::IsotopeInventoryImbalance,
            "isotope-resolved atom inventories differ",
            DiagnosticSeverity::Error,
        ));
    }
    if (reactant_charge - product_charge).abs() > EPSILON {
        diagnostics.push(diag(
            "/steps",
            StoichiometryIssueCode::ChargeImbalance,
            "formal charges differ",
            DiagnosticSeverity::Error,
        ));
    }
    let status = if !has_reactants || !has_products {
        StoichiometryStatus::UnderSpecified
    } else if diagnostics
        .iter()
        .any(|d| d.severity == DiagnosticSeverity::Error)
    {
        StoichiometryStatus::Unbalanced
    } else {
        StoichiometryStatus::Balanced
    };
    Ok(StepStoichiometryReport {
        step_id: step.id.clone(),
        status,
        components,
        reactant_atoms,
        product_atoms,
        reactant_charge,
        product_charge,
        diagnostics,
    })
}

fn inventory(
    molecule: &chematic_core::Molecule,
    coefficient: f64,
) -> (BTreeMap<(u8, Option<u16>), f64>, f64) {
    let mut atoms = BTreeMap::new();
    let mut charge = 0.0;
    for (_, atom) in molecule.atoms() {
        *atoms
            .entry((atom.element.atomic_number(), atom.isotope))
            .or_insert(0.0) += coefficient;
        charge += f64::from(atom.charge) * coefficient;
    }
    (atoms, charge)
}

fn merge(target: &mut BTreeMap<(u8, Option<u16>), f64>, source: &BTreeMap<(u8, Option<u16>), f64>) {
    for (key, value) in source {
        *target.entry(*key).or_insert(0.0) += value;
    }
}
fn to_inventory(map: BTreeMap<(u8, Option<u16>), f64>) -> Vec<AtomInventory> {
    map.into_iter()
        .map(|((atomic_number, isotope), count)| AtomInventory {
            atomic_number,
            isotope,
            count,
        })
        .collect()
}
fn element_totals(atoms: &[AtomInventory]) -> BTreeMap<u8, f64> {
    let mut out = BTreeMap::new();
    for atom in atoms {
        *out.entry(atom.atomic_number).or_insert(0.0) += atom.count;
    }
    out
}
fn diag(
    path: &str,
    code: StoichiometryIssueCode,
    message: &str,
    severity: DiagnosticSeverity,
) -> StoichiometryDiagnostic {
    StoichiometryDiagnostic {
        path: path.into(),
        code,
        severity,
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn component(
        id: &str,
        role: ComponentRole,
        smiles: &str,
        coefficient: f64,
    ) -> StoichiometryComponent {
        StoichiometryComponent {
            id: id.into(),
            role,
            smiles: smiles.into(),
            coefficient,
        }
    }

    #[test]
    fn balanced_fractional_and_isotope_inventory_are_reported() {
        let report = analyze_components(&StoichiometryStep {
            id: "s1".into(),
            components: vec![
                component("r", ComponentRole::Reactant, "[13CH4]", 0.5),
                component("p", ComponentRole::Product, "[13CH4]", 0.5),
            ],
        })
        .unwrap();
        assert_eq!(report.status, StoichiometryStatus::Balanced);
        assert_eq!(report.components[0].atoms[0].isotope, Some(13));
        assert!(
            report
                .diagnostics
                .iter()
                .any(|d| d.code == StoichiometryIssueCode::FractionalCoefficient)
        );
        assert_eq!(
            serde_json::from_str::<StepStoichiometryReport>(
                &serde_json::to_string(&report).unwrap()
            )
            .unwrap(),
            report
        );
    }

    #[test]
    fn distinguishes_atom_isotope_and_charge_imbalances() {
        let report = analyze_components(&StoichiometryStep {
            id: "s1".into(),
            components: vec![
                component("r", ComponentRole::Reactant, "[13CH4+]", 1.0),
                component("p", ComponentRole::Product, "[CH4]", 1.0),
            ],
        })
        .unwrap();
        assert_eq!(report.status, StoichiometryStatus::Unbalanced);
        assert!(
            report
                .diagnostics
                .iter()
                .any(|d| d.code == StoichiometryIssueCode::IsotopeInventoryImbalance)
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|d| d.code == StoichiometryIssueCode::ChargeImbalance)
        );
    }

    #[test]
    fn agent_only_is_under_specified() {
        let report = analyze_components(&StoichiometryStep {
            id: "s1".into(),
            components: vec![component("a", ComponentRole::Agent, "O", 1.0)],
        })
        .unwrap();
        assert_eq!(report.status, StoichiometryStatus::UnderSpecified);
        assert!(
            report
                .diagnostics
                .iter()
                .any(|d| d.code == StoichiometryIssueCode::AgentOnly)
        );
    }
}
