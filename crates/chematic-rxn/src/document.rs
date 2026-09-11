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
    /// Explicit atom-map identities in the component's serialized atom order.
    ///
    /// This is additive metadata for document consumers that need stable
    /// references without reparsing SMILES. An omitted/empty list preserves
    /// compatibility with older authored JSON; adapters that populate it must
    /// agree with the map annotations in `smiles`.
    #[serde(default)]
    pub atom_maps: Vec<ReactionAtomMap>,
    /// Positive stoichiometric coefficient.
    #[serde(default = "one")]
    pub coefficient: u32,
    #[serde(default = "default_origin")]
    pub origin: ContentOrigin,
}

/// An atom-map number and its stable zero-based position within one component.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReactionAtomMap {
    pub map_number: u16,
    pub atom_index: u32,
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

/// Bounded, ID-addressed edits exposed by document bindings.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind")]
pub enum ReactionDocumentEdit {
    SetDocumentId {
        id: String,
    },
    SetStepCondition {
        step_id: String,
        key: String,
        value: String,
    },
    SetComponentCoefficient {
        component_id: String,
        coefficient: u32,
    },
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
    /// Deserialize and validate a typed reaction document from JSON.
    ///
    /// `serde_json::from_str::<ReactionDocument>` remains available for
    /// callers that need raw serde behaviour, but format and binding
    /// boundaries should use this constructor so malformed documents cannot
    /// bypass the model's invariants.
    pub fn from_json_str(s: &str) -> Result<Self, ReactionDocumentError> {
        let document: Self = serde_json::from_str(s).map_err(|error| {
            ReactionDocumentError::InvalidDocument(format!("invalid JSON: {error}"))
        })?;
        document.validate()?;
        Ok(document)
    }

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
                let smiles = chematic_smiles::write(molecule);
                components.push(ReactionComponent {
                    id: format!("{}-{}", role_name(role), index + 1),
                    role,
                    atom_maps: atom_maps_for_smiles(&smiles),
                    smiles,
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

    /// Apply one bounded edit, addressed by stable document IDs, and revalidate.
    pub fn apply_json_edit(&self, edit_json: &str) -> Result<Self, ReactionDocumentError> {
        let edit: ReactionDocumentEdit = serde_json::from_str(edit_json).map_err(|error| {
            ReactionDocumentError::InvalidDocument(format!("invalid reaction edit JSON: {error}"))
        })?;
        let mut next = self.clone();
        match edit {
            ReactionDocumentEdit::SetDocumentId { id } => next.id = id,
            ReactionDocumentEdit::SetStepCondition {
                step_id,
                key,
                value,
            } => {
                let step = next
                    .steps
                    .iter_mut()
                    .find(|step| step.id == step_id)
                    .ok_or_else(|| {
                        ReactionDocumentError::InvalidDocument(format!(
                            "unknown step id '{step_id}'"
                        ))
                    })?;
                if let Some(condition) = step
                    .conditions
                    .iter_mut()
                    .find(|condition| condition.key == key)
                {
                    condition.value = value;
                } else {
                    step.conditions.push(ReactionCondition { key, value });
                }
            }
            ReactionDocumentEdit::SetComponentCoefficient {
                component_id,
                coefficient,
            } => {
                let component = next
                    .steps
                    .iter_mut()
                    .flat_map(|step| &mut step.components)
                    .find(|component| component.id == component_id)
                    .ok_or_else(|| {
                        ReactionDocumentError::InvalidDocument(format!(
                            "unknown component id '{component_id}'"
                        ))
                    })?;
                component.coefficient = coefficient;
            }
        }
        next.validate()?;
        Ok(next)
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
            let mut condition_keys = std::collections::HashSet::new();
            for condition in &step.conditions {
                if condition.key.is_empty() {
                    return Err(ReactionDocumentError::InvalidDocument(
                        "reaction condition keys must be non-empty".to_string(),
                    ));
                }
                if !condition_keys.insert(&condition.key) {
                    return Err(ReactionDocumentError::InvalidDocument(
                        "reaction condition keys must be unique within a step".to_string(),
                    ));
                }
            }
            for provenance in &step.provenance {
                if provenance.source.is_empty() || provenance.kind.is_empty() {
                    return Err(ReactionDocumentError::InvalidDocument(
                        "step provenance source and kind must be non-empty".to_string(),
                    ));
                }
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
                // A component is one molecule, not a reaction section.  Going
                // through `parse_reaction` here could accidentally accept a
                // payload containing `>` when the surrounding text happened
                // to split into three parseable sections.  Keep this boundary
                // loss-aware by validating the payload with the molecule
                // parser directly.
                let molecule = chematic_smiles::parse(&component.smiles).map_err(|e| {
                    ReactionDocumentError::Parse(RxnErrorMessage(format!(
                        "invalid component '{}': {e}",
                        component.id
                    )))
                })?;
                if !component.atom_maps.is_empty() {
                    let mut seen = std::collections::HashSet::new();
                    for atom_map in &component.atom_maps {
                        if atom_map.atom_index as usize >= molecule.atom_count() {
                            return Err(ReactionDocumentError::InvalidDocument(format!(
                                "atom-map index {} is outside component '{}'",
                                atom_map.atom_index, component.id
                            )));
                        }
                        if !seen.insert(atom_map.map_number) {
                            return Err(ReactionDocumentError::InvalidDocument(format!(
                                "atom-map number {} is duplicated in component '{}'",
                                atom_map.map_number, component.id
                            )));
                        }
                        if molecule
                            .atom(chematic_core::AtomIdx(atom_map.atom_index))
                            .atom_map
                            != Some(atom_map.map_number)
                        {
                            return Err(ReactionDocumentError::InvalidDocument(format!(
                                "atom-map identity {}:{} does not match component '{}' SMILES",
                                atom_map.map_number, atom_map.atom_index, component.id
                            )));
                        }
                    }
                    let expected = atom_maps_for_molecule(&molecule);
                    if component.atom_maps != expected {
                        return Err(ReactionDocumentError::InvalidDocument(format!(
                            "atom-map identities do not cover component '{}' SMILES",
                            component.id
                        )));
                    }
                }
            }
        }
        for provenance in &self.provenance {
            if provenance.source.is_empty() || provenance.kind.is_empty() {
                return Err(ReactionDocumentError::InvalidDocument(
                    "document provenance source and kind must be non-empty".to_string(),
                ));
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

fn atom_maps_for_molecule(molecule: &chematic_core::Molecule) -> Vec<ReactionAtomMap> {
    molecule
        .atoms()
        .map(|(index, atom)| (index.0, atom))
        .filter_map(|(atom_index, atom)| {
            atom.atom_map.map(|map_number| ReactionAtomMap {
                map_number,
                atom_index,
            })
        })
        .collect()
}

fn atom_maps_for_smiles(smiles: &str) -> Vec<ReactionAtomMap> {
    chematic_smiles::parse(smiles)
        .map(|molecule| atom_maps_for_molecule(&molecule))
        .unwrap_or_default()
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
        assert_eq!(ReactionDocument::from_json_str(&json).unwrap(), document);
    }

    #[test]
    fn checked_json_constructor_rejects_invalid_component() {
        let document = ReactionDocument::from_reaction_smiles("CC>>CC").unwrap();
        let mut value = serde_json::to_value(document).unwrap();
        value["steps"][0]["components"][0]["smiles"] = serde_json::json!("C>C");

        let error =
            ReactionDocument::from_json_str(&serde_json::to_string(&value).unwrap()).unwrap_err();
        assert!(matches!(error, ReactionDocumentError::Parse(_)));
        assert!(error.to_string().contains("invalid component"));
    }

    #[test]
    fn derived_document_exposes_atom_map_identities() {
        let document = ReactionDocument::from_reaction_smiles("[CH3:7]>>[CH3:7]").unwrap();
        assert_eq!(
            document.steps[0].components[0].atom_maps,
            vec![ReactionAtomMap {
                map_number: 7,
                atom_index: 0,
            }]
        );
        assert_eq!(
            document.steps[0].components[1].atom_maps,
            vec![ReactionAtomMap {
                map_number: 7,
                atom_index: 0,
            }]
        );
        document.validate().unwrap();
    }

    #[test]
    fn validation_rejects_atom_map_identity_mismatch() {
        let mut document = ReactionDocument::from_reaction_smiles("[CH3:7]>>[CH3:7]").unwrap();
        document.steps[0].components[0].atom_maps[0].map_number = 8;
        let error = document.validate().unwrap_err();
        assert!(
            matches!(error, ReactionDocumentError::InvalidDocument(message) if message.contains("atom-map identity"))
        );
    }

    #[test]
    fn component_validation_rejects_reaction_payloads() {
        let mut document = ReactionDocument::from_reaction_smiles("CC>>CC").unwrap();
        document.steps[0].components[0].smiles = "C>C".to_string();

        let error = document.validate().unwrap_err();
        assert!(matches!(error, ReactionDocumentError::Parse(_)));
        assert!(error.to_string().contains("invalid component"));
    }

    #[test]
    fn bounded_json_edit_preserves_stable_ids_and_metadata() {
        let mut document = ReactionDocument::from_reaction_smiles("CC>>CC").unwrap();
        document.provenance.push(ProvenanceRecord {
            source: "fixture".into(),
            kind: "authored".into(),
            note: Some("keep".into()),
        });
        let edited = document
            .apply_json_edit(r#"{"kind":"set_step_condition","step_id":"step-1","key":"temperature","value":"25 C"}"#)
            .unwrap();
        assert_eq!(edited.id, document.id);
        assert_eq!(
            edited.steps[0].components[0].id,
            document.steps[0].components[0].id
        );
        assert_eq!(edited.provenance, document.provenance);
        assert_eq!(edited.steps[0].conditions[0].value, "25 C");
    }

    #[test]
    fn bounded_json_edit_rejects_unknown_ids_and_zero_coefficients() {
        let document = ReactionDocument::from_reaction_smiles("CC>>CC").unwrap();
        let unknown = document.apply_json_edit(
            r#"{"kind":"set_component_coefficient","component_id":"missing","coefficient":2}"#,
        );
        assert!(
            unknown
                .unwrap_err()
                .to_string()
                .contains("unknown component id")
        );
        let zero = document.apply_json_edit(
            r#"{"kind":"set_component_coefficient","component_id":"reactant-1","coefficient":0}"#,
        );
        assert!(
            zero.unwrap_err()
                .to_string()
                .contains("positive coefficient")
        );
    }

    #[test]
    fn validation_rejects_unidentifiable_metadata() {
        let mut document = ReactionDocument::from_reaction_smiles("CC>>CC").unwrap();
        document.steps[0].conditions.push(ReactionCondition {
            key: String::new(),
            value: "25 C".into(),
        });
        assert!(matches!(
            document.validate(),
            Err(ReactionDocumentError::InvalidDocument(_))
        ));

        let mut document = ReactionDocument::from_reaction_smiles("CC>>CC").unwrap();
        document.provenance.push(ProvenanceRecord {
            source: "source-id".into(),
            kind: String::new(),
            note: None,
        });
        assert!(matches!(
            document.validate(),
            Err(ReactionDocumentError::InvalidDocument(_))
        ));
    }

    #[test]
    fn validation_rejects_duplicate_condition_keys() {
        let mut document = ReactionDocument::from_reaction_smiles("CC>>CC").unwrap();
        document.steps[0].conditions = vec![
            ReactionCondition {
                key: "temperature".into(),
                value: "20 C".into(),
            },
            ReactionCondition {
                key: "temperature".into(),
                value: "25 C".into(),
            },
        ];
        let error = document.validate().unwrap_err();
        assert!(
            matches!(error, ReactionDocumentError::InvalidDocument(message) if message.contains("unique"))
        );
    }
}
