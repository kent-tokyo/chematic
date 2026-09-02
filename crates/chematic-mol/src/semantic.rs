//! Explicit semantic constructs that cannot be represented by an ordinary molecule.
//!
//! This module is intentionally conservative: ambiguous references and unsafe
//! expansions are rejected instead of being guessed into a molecular graph.

use std::collections::BTreeMap;

use serde_json::{Map, Value};

/// Stable identifier used by source-level semantic objects.
pub type SemanticId = String;

/// A source atom reference used by Markush and polymer constructs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomRef {
    pub atom_id: SemanticId,
}

/// A variable substituent definition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RGroupDefinition {
    pub id: SemanticId,
    pub attachment_atoms: Vec<AtomRef>,
    pub alternatives: Vec<String>,
    pub selected_alternative: Option<usize>,
}

/// A polymer repeat unit with explicit linkage and end-group references.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolymerRepeatUnit {
    pub id: SemanticId,
    pub attachment_atoms: Vec<AtomRef>,
    pub end_groups: Vec<String>,
    pub repeat_count: Option<u32>,
}

/// Loss/unsupported reason returned by validation or expansion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SemanticError {
    DuplicateId(String),
    MissingAtom(String),
    MissingAlternative(String),
    InvalidAlternative { id: String, reason: String },
    AmbiguousAttachment(String),
    Unsupported { construct: String, reason: String },
    InvalidJson(String),
}

impl std::fmt::Display for SemanticError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::DuplicateId(id) => write!(f, "duplicate semantic id: {id}"),
            Self::MissingAtom(id) => write!(f, "semantic reference targets missing atom: {id}"),
            Self::MissingAlternative(id) => write!(f, "missing R-group alternative: {id}"),
            Self::InvalidAlternative { id, reason } => {
                write!(f, "invalid alternative {id}: {reason}")
            }
            Self::AmbiguousAttachment(id) => write!(f, "ambiguous attachment topology: {id}"),
            Self::Unsupported { construct, reason } => {
                write!(f, "unsupported {construct}: {reason}")
            }
            Self::InvalidJson(reason) => write!(f, "invalid semantic JSON: {reason}"),
        }
    }
}

impl std::error::Error for SemanticError {}

/// Explicit source-level semantic model. It does not silently flatten into a molecule.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SemanticModel {
    pub atom_ids: Vec<SemanticId>,
    pub bond_ids: Vec<SemanticId>,
    pub r_groups: Vec<RGroupDefinition>,
    pub polymer_units: Vec<PolymerRepeatUnit>,
    pub extensions: BTreeMap<String, Value>,
}

impl SemanticModel {
    /// Validate IDs, references, attachment arity, and SMARTS-like alternatives.
    pub fn validate(&self) -> Result<(), SemanticError> {
        let mut ids = std::collections::BTreeSet::new();
        for id in self.atom_ids.iter().chain(self.bond_ids.iter()) {
            if !ids.insert(id.clone()) {
                return Err(SemanticError::DuplicateId(id.clone()));
            }
        }
        for group in &self.r_groups {
            if !ids.insert(group.id.clone()) {
                return Err(SemanticError::DuplicateId(group.id.clone()));
            }
            if group.attachment_atoms.is_empty() {
                return Err(SemanticError::AmbiguousAttachment(group.id.clone()));
            }
            for a in &group.attachment_atoms {
                if !self.atom_ids.contains(&a.atom_id) {
                    return Err(SemanticError::MissingAtom(a.atom_id.clone()));
                }
            }
            for pattern in &group.alternatives {
                if pattern.trim().is_empty() {
                    return Err(SemanticError::InvalidAlternative {
                        id: group.id.clone(),
                        reason: "empty query".into(),
                    });
                }
            }
            if let Some(i) = group.selected_alternative
                && i >= group.alternatives.len()
            {
                return Err(SemanticError::MissingAlternative(group.id.clone()));
            }
        }
        for unit in &self.polymer_units {
            if !ids.insert(unit.id.clone()) {
                return Err(SemanticError::DuplicateId(unit.id.clone()));
            }
            if unit.attachment_atoms.len() != 2 {
                return Err(SemanticError::AmbiguousAttachment(unit.id.clone()));
            }
            for a in &unit.attachment_atoms {
                if !self.atom_ids.contains(&a.atom_id) {
                    return Err(SemanticError::MissingAtom(a.atom_id.clone()));
                }
            }
        }
        Ok(())
    }

    /// Serialize the model to a stable JSON object for WASM/editor exchange.
    pub fn to_json(&self) -> Value {
        let mut root = Map::new();
        root.insert(
            "schema".into(),
            Value::String("chematic.semantic.v1".into()),
        );
        root.insert(
            "atom_ids".into(),
            Value::Array(self.atom_ids.iter().cloned().map(Value::String).collect()),
        );
        root.insert(
            "bond_ids".into(),
            Value::Array(self.bond_ids.iter().cloned().map(Value::String).collect()),
        );
        root.insert("r_groups".into(), Value::Array(self.r_groups.iter().map(|g| serde_json::json!({
            "id": g.id, "attachment_atoms": g.attachment_atoms.iter().map(|a| &a.atom_id).collect::<Vec<_>>(),
            "alternatives": g.alternatives, "selected_alternative": g.selected_alternative
        })).collect()));
        root.insert("polymer_units".into(), Value::Array(self.polymer_units.iter().map(|u| serde_json::json!({
            "id": u.id, "attachment_atoms": u.attachment_atoms.iter().map(|a| &a.atom_id).collect::<Vec<_>>(),
            "end_groups": u.end_groups, "repeat_count": u.repeat_count
        })).collect()));
        root.insert(
            "extensions".into(),
            Value::Object(self.extensions.clone().into_iter().collect()),
        );
        Value::Object(root)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validates_r_group_references_and_serializes_schema() {
        let model = SemanticModel {
            atom_ids: vec!["a1".into(), "a2".into()],
            r_groups: vec![RGroupDefinition {
                id: "r1".into(),
                attachment_atoms: vec![AtomRef {
                    atom_id: "a1".into(),
                }],
                alternatives: vec!["[*]C".into()],
                selected_alternative: Some(0),
            }],
            ..Default::default()
        };
        model.validate().unwrap();
        assert_eq!(model.to_json()["schema"], "chematic.semantic.v1");
    }

    #[test]
    fn rejects_dangling_and_ambiguous_constructs() {
        let model = SemanticModel {
            atom_ids: vec!["a1".into()],
            polymer_units: vec![PolymerRepeatUnit {
                id: "p1".into(),
                attachment_atoms: vec![AtomRef {
                    atom_id: "missing".into(),
                }],
                end_groups: vec![],
                repeat_count: None,
            }],
            ..Default::default()
        };
        assert!(matches!(
            model.validate(),
            Err(SemanticError::AmbiguousAttachment(_))
        ));
    }
}
