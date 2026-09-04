//! Typed, loss-aware reaction-document model.
//!
//! This is the foundation for rich reaction exchange. It deliberately keeps
//! molecule payloads as serialized SMILES so the document model can be
//! serialized without making `Molecule` part of the public wire contract.
//! Legacy reaction SMILES conversion is checked and rejects information that
//! the legacy three-section format cannot represent.

use serde::{Deserialize, Serialize};

use crate::reaction::{Reaction, parse_reaction, write_reaction};

/// Semantic role of a reaction component.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ComponentRole {
    Reactant,
    Agent,
    Product,
}

/// Whether content was supplied by the author or derived by an adapter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ContentOrigin {
    Authored,
    Derived,
}

/// A stable, typed reaction component.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReactionComponent {
    /// Stable within the containing document.
    pub id: String,
    pub role: ComponentRole,
    /// Canonical or source SMILES for the component.
    pub smiles: String,
    /// Positive stoichiometric coefficient.
    #[serde(default = "one")]
    pub coefficient: u32,
    #[serde(default = "default_origin")]
    pub origin: ContentOrigin,
}

/// A named reaction condition, kept ordered for deterministic serialization.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReactionCondition {
    pub key: String,
    pub value: String,
}

/// Provenance attached to authored or derived document content.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProvenanceRecord {
    pub source: String,
    pub kind: String,
    #[serde(default)]
    pub note: Option<String>,
}

/// One ordered reaction step.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReactionStep {
    pub id: String,
    pub components: Vec<ReactionComponent>,
    #[serde(default)]
    pub conditions: Vec<ReactionCondition>,
    #[serde(default)]
    pub provenance: Vec<ProvenanceRecord>,
    #[serde(default = "default_origin")]
    pub origin: ContentOrigin,
}

/// A rich reaction document. Steps are ordered and are never implicitly
/// flattened by a legacy adapter.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReactionDocument {
    pub id: String,
    pub steps: Vec<ReactionStep>,
    #[serde(default)]
    pub provenance: Vec<ProvenanceRecord>,
}

/// Information that prevents a lossless conversion to a legacy format.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReactionLoss {
    pub field: String,
    pub detail: String,
}

/// Error returned when a rich document cannot be represented losslessly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReactionDocumentError {
    Parse(RxnErrorMessage),
    InvalidDocument(String),
    Losses(Vec<ReactionLoss>),
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RxnErrorMessage(pub String);

impl core::fmt::Display for ReactionDocumentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Parse(e) => write!(f, "reaction document parse error: {}", e.0),
            Self::InvalidDocument(e) => write!(f, "invalid reaction document: {e}"),
            Self::Losses(losses) => {
                write!(f, "legacy conversion would lose {} field(s)", losses.len())
            }
        }
    }
}

impl std::error::Error for ReactionDocumentError {}

impl ReactionDocumentError {
    /// Build a typed parse error for format adapters that cannot expose their
    /// parser's concrete error type through this crate's dependency graph.
    pub fn parse_message(message: impl Into<String>) -> Self {
        Self::Parse(RxnErrorMessage(message.into()))
    }
}

impl ReactionDocument {
    /// Parse a legacy reaction SMILES into a one-step derived document.
    pub fn from_reaction_smiles(s: &str) -> Result<Self, ReactionDocumentError> {
        let reaction = parse_reaction(s)
            .map_err(|e| ReactionDocumentError::Parse(RxnErrorMessage(e.to_string())))?;
        Ok(Self::from_reaction(&reaction))
    }

    /// Convert the existing reaction model into a derived one-step document.
    pub fn from_reaction(reaction: &Reaction) -> Self {
        let mut components = Vec::new();
        for (role, molecules) in [
            (ComponentRole::Reactant, &reaction.reactants),
            (ComponentRole::Agent, &reaction.agents),
            (ComponentRole::Product, &reaction.products),
        ] {
            for (index, molecule) in molecules.iter().enumerate() {
                components.push(ReactionComponent {
                    id: format!("{}-{}", role_name(role), index + 1),
                    role,
                    smiles: chematic_smiles::write(molecule),
                    coefficient: 1,
                    origin: ContentOrigin::Derived,
                });
            }
        }
        Self {
            id: "step-1".to_string(),
            steps: vec![ReactionStep {
                id: "step-1".to_string(),
                components,
                conditions: Vec::new(),
                provenance: Vec::new(),
                origin: ContentOrigin::Derived,
            }],
            provenance: Vec::new(),
        }
    }

    /// Validate IDs, roles, coefficients, and non-empty SMILES payloads.
    pub fn validate(&self) -> Result<(), ReactionDocumentError> {
        if self.id.is_empty() || self.steps.is_empty() {
            return Err(ReactionDocumentError::InvalidDocument(
                "document id and at least one step are required".to_string(),
            ));
        }
        let mut ids = std::collections::HashSet::new();
        for step in &self.steps {
            if step.id.is_empty() || !ids.insert(step.id.clone()) {
                return Err(ReactionDocumentError::InvalidDocument(
                    "step IDs must be non-empty and unique".to_string(),
                ));
            }
            for component in &step.components {
                if component.id.is_empty() || !ids.insert(component.id.clone()) {
                    return Err(ReactionDocumentError::InvalidDocument(
                        "component IDs must be non-empty and globally unique".to_string(),
                    ));
                }
                if component.coefficient == 0 || component.smiles.is_empty() {
                    return Err(ReactionDocumentError::InvalidDocument(
                        "component SMILES and positive coefficient are required".to_string(),
                    ));
                }
                parse_reaction(&format!("{}>>{}", component.smiles, component.smiles))
                    .map_err(|e| ReactionDocumentError::Parse(RxnErrorMessage(e.to_string())))?;
            }
        }
        Ok(())
    }

    /// Convert to legacy reaction SMILES, rejecting every unsupported loss.
    pub fn to_reaction_smiles(&self) -> Result<String, ReactionDocumentError> {
        self.validate()?;
        let mut losses = Vec::new();
        if self.steps.len() != 1 {
            losses.push(ReactionLoss {
                field: "steps".to_string(),
                detail: "legacy reaction SMILES has one reaction boundary".to_string(),
            });
        }
        if !self.provenance.is_empty() {
            losses.push(ReactionLoss {
                field: "provenance".to_string(),
                detail: "legacy reaction SMILES has no provenance channel".to_string(),
            });
        }
        let step = &self.steps[0];
        for condition in &step.conditions {
            losses.push(ReactionLoss {
                field: format!("conditions.{}", condition.key),
                detail: "legacy reaction SMILES has no conditions channel".to_string(),
            });
        }
        if !step.provenance.is_empty() {
            losses.push(ReactionLoss {
                field: "step.provenance".to_string(),
                detail: "legacy reaction SMILES has no provenance channel".to_string(),
            });
        }
        for component in &step.components {
            if component.coefficient != 1 {
                losses.push(ReactionLoss {
                    field: format!("{}.coefficient", component.id),
                    detail: "legacy reaction SMILES has no stoichiometric coefficient".to_string(),
                });
            }
        }
        if !losses.is_empty() {
            return Err(ReactionDocumentError::Losses(losses));
        }
        let mut reaction = Reaction {
            reactants: Vec::new(),
            agents: Vec::new(),
            products: Vec::new(),
        };
        for component in &step.components {
            let molecule = chematic_smiles::parse(&component.smiles)
                .map_err(|e| ReactionDocumentError::Parse(RxnErrorMessage(e.to_string())))?;
            match component.role {
                ComponentRole::Reactant => reaction.reactants.push(molecule),
                ComponentRole::Agent => reaction.agents.push(molecule),
                ComponentRole::Product => reaction.products.push(molecule),
            }
        }
        Ok(write_reaction(&reaction))
    }
}

fn one() -> u32 {
    1
}

fn default_origin() -> ContentOrigin {
    ContentOrigin::Derived
}

fn role_name(role: ComponentRole) -> &'static str {
    match role {
        ComponentRole::Reactant => "reactant",
        ComponentRole::Agent => "agent",
        ComponentRole::Product => "product",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_reaction_round_trips_through_derived_document() {
        let document = ReactionDocument::from_reaction_smiles("CCO>O>[CH3:1][OH:2]").unwrap();
        assert_eq!(
            document.to_reaction_smiles().unwrap(),
            "CCO>O>[CH3:1][OH:2]"
        );
    }

    #[test]
    fn rich_fields_reject_legacy_flattening() {
        let mut document = ReactionDocument::from_reaction_smiles("CC>>CC").unwrap();
        document.steps[0].components[0].coefficient = 2;
        let err = document.to_reaction_smiles().unwrap_err();
        assert!(matches!(err, ReactionDocumentError::Losses(_)));
    }

    #[test]
    fn serde_round_trip_preserves_typed_document() {
        let document = ReactionDocument::from_reaction_smiles("CC>>CC").unwrap();
        let json = serde_json::to_string(&document).unwrap();
        let decoded: ReactionDocument = serde_json::from_str(&json).unwrap();
        assert_eq!(decoded, document);
    }
}
