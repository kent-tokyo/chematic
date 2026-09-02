//! Explicit semantic constructs that cannot be represented by an ordinary molecule.
//!
//! This module is intentionally conservative: ambiguous references and unsafe
//! expansions are rejected instead of being guessed into a molecular graph.

use std::collections::BTreeMap;

use serde_json::{Map, Value};

use chematic_core::{AtomIdx, Molecule};

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
    /// Optional SMILES repeat fragment. It may use `[*]` at both ends, or use
    /// `repeat_endpoint_atoms` to identify two explicit endpoint atoms.
    pub repeat_smiles: Option<String>,
    /// Zero-based atom indices in `repeat_smiles` used when no `[*]` markers
    /// are present. Exactly two distinct endpoints are required.
    pub repeat_endpoint_atoms: Option<[u32; 2]>,
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
    InvalidExpansion { id: String, reason: String },
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
            Self::InvalidExpansion { id, reason } => write!(f, "cannot expand {id}: {reason}"),
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

/// Result of a checked semantic expansion. The mapping is required for undo/edit flows.
#[derive(Clone)]
pub struct ExpandedSemantic {
    pub molecule: Molecule,
    pub source_to_expanded: BTreeMap<SemanticId, Vec<AtomIdx>>,
}

/// Immutable, auditable edits to a semantic model.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SemanticCommand {
    SelectRGroupAlternative {
        group_id: SemanticId,
        alternative: usize,
    },
}

impl SemanticModel {
    /// Apply one command while preserving stable IDs and rejecting ambiguity.
    pub fn apply(&self, command: &SemanticCommand) -> Result<Self, SemanticError> {
        let mut next = self.clone();
        match command {
            SemanticCommand::SelectRGroupAlternative {
                group_id,
                alternative,
            } => {
                let group = next
                    .r_groups
                    .iter_mut()
                    .find(|g| &g.id == group_id)
                    .ok_or_else(|| SemanticError::InvalidExpansion {
                        id: group_id.clone(),
                        reason: "unknown R-group".into(),
                    })?;
                if *alternative >= group.alternatives.len() {
                    return Err(SemanticError::MissingAlternative(group_id.clone()));
                }
                group.selected_alternative = Some(*alternative);
            }
        }
        next.validate()?;
        Ok(next)
    }

    /// Expand only explicitly selected R-groups. No alternative is guessed.
    /// Supported alternatives use a leading `[*]` attachment placeholder.
    pub fn expand(&self, base: &Molecule) -> Result<ExpandedSemantic, SemanticError> {
        self.validate()?;
        if self.atom_ids.len() != base.atom_count() {
            return Err(SemanticError::InvalidExpansion {
                id: "model".into(),
                reason: "atom_ids must match base molecule atom count".into(),
            });
        }
        let mut molecule = base.clone();
        let mut mapping: BTreeMap<SemanticId, Vec<AtomIdx>> = self
            .atom_ids
            .iter()
            .cloned()
            .zip((0..base.atom_count()).map(|i| vec![AtomIdx(i as u32)]))
            .collect();
        for group in &self.r_groups {
            let Some(selected) = group.selected_alternative else {
                return Err(SemanticError::Unsupported {
                    construct: group.id.clone(),
                    reason: "no explicit alternative selected".into(),
                });
            };
            let pattern = &group.alternatives[selected];
            let fragment =
                chematic_smiles::parse(pattern).map_err(|e| SemanticError::InvalidExpansion {
                    id: group.id.clone(),
                    reason: e.to_string(),
                })?;
            if fragment.atom_count() < 2 || group.attachment_atoms.len() != 1 {
                return Err(SemanticError::Unsupported {
                    construct: group.id.clone(),
                    reason: "requires one attachment and a fragment with a leading placeholder"
                        .into(),
                });
            }
            let base_atom = AtomIdx(
                self.atom_ids
                    .iter()
                    .position(|id| id == &group.attachment_atoms[0].atom_id)
                    .ok_or_else(|| {
                        SemanticError::MissingAtom(group.attachment_atoms[0].atom_id.clone())
                    })? as u32,
            );
            if base_atom.0 as usize >= molecule.atom_count() {
                return Err(SemanticError::MissingAtom(
                    group.attachment_atoms[0].atom_id.clone(),
                ));
            }
            let mut remap = BTreeMap::new();
            for (idx, atom) in fragment.atoms() {
                if idx.0 == 0 {
                    continue;
                }
                remap.insert(idx, molecule.add_atom(atom.clone()));
            }
            for (_, bond) in fragment.bonds() {
                if bond.atom1.0 == 0 || bond.atom2.0 == 0 {
                    continue;
                }
                molecule
                    .add_bond(remap[&bond.atom1], remap[&bond.atom2], bond.order)
                    .map_err(|e| SemanticError::InvalidExpansion {
                        id: group.id.clone(),
                        reason: e.to_string(),
                    })?;
            }
            let attach =
                remap
                    .values()
                    .next()
                    .copied()
                    .ok_or_else(|| SemanticError::InvalidExpansion {
                        id: group.id.clone(),
                        reason: "empty expansion".into(),
                    })?;
            molecule
                .add_bond(base_atom, attach, chematic_core::BondOrder::Single)
                .map_err(|e| SemanticError::InvalidExpansion {
                    id: group.id.clone(),
                    reason: e.to_string(),
                })?;
            mapping.insert(group.id.clone(), remap.values().copied().collect());
        }
        for unit in &self.polymer_units {
            let pattern =
                unit.repeat_smiles
                    .as_deref()
                    .ok_or_else(|| SemanticError::Unsupported {
                        construct: unit.id.clone(),
                        reason: "repeat_smiles is not provided".into(),
                    })?;
            let placeholder_endpoints = pattern.trim().starts_with("[*]")
                && pattern.trim().ends_with("[*]")
                && pattern.matches("[*]").count() == 2;
            if !placeholder_endpoints && unit.repeat_endpoint_atoms.is_none() {
                return Err(SemanticError::InvalidExpansion {
                    id: unit.id.clone(),
                    reason: "repeat requires two [*] markers or explicit endpoint atoms".into(),
                });
            }
            let repeats = unit
                .repeat_count
                .ok_or_else(|| SemanticError::Unsupported {
                    construct: unit.id.clone(),
                    reason: "repeat count must be explicit".into(),
                })?;
            let left = AtomIdx(
                self.atom_ids
                    .iter()
                    .position(|id| id == &unit.attachment_atoms[0].atom_id)
                    .ok_or_else(|| {
                        SemanticError::MissingAtom(unit.attachment_atoms[0].atom_id.clone())
                    })? as u32,
            );
            let right = AtomIdx(
                self.atom_ids
                    .iter()
                    .position(|id| id == &unit.attachment_atoms[1].atom_id)
                    .ok_or_else(|| {
                        SemanticError::MissingAtom(unit.attachment_atoms[1].atom_id.clone())
                    })? as u32,
            );
            let mut unit_atoms = Vec::new();
            let mut previous_right: Option<AtomIdx> = None;
            for repeat_index in 0..repeats {
                let fragment = chematic_smiles::parse(pattern).map_err(|e| {
                    SemanticError::InvalidExpansion {
                        id: unit.id.clone(),
                        reason: e.to_string(),
                    }
                })?;
                let (endpoint_left, endpoint_right, excluded) = if placeholder_endpoints {
                    let last = AtomIdx(fragment.atom_count().saturating_sub(1) as u32);
                    let ln = fragment.neighbors(AtomIdx(0)).collect::<Vec<_>>();
                    let rn = fragment.neighbors(last).collect::<Vec<_>>();
                    if fragment.atom_count() < 4 || ln.len() != 1 || rn.len() != 1 {
                        return Err(SemanticError::Unsupported {
                            construct: unit.id.clone(),
                            reason: "repeat requires exactly one neighbor at each linkage".into(),
                        });
                    }
                    (ln[0].0, rn[0].0, Some((AtomIdx(0), last)))
                } else {
                    let [left, right] = unit.repeat_endpoint_atoms.ok_or_else(|| {
                        SemanticError::InvalidExpansion {
                            id: unit.id.clone(),
                            reason: "repeat endpoints are missing".into(),
                        }
                    })?;
                    if left == right
                        || left as usize >= fragment.atom_count()
                        || right as usize >= fragment.atom_count()
                    {
                        return Err(SemanticError::InvalidExpansion {
                            id: unit.id.clone(),
                            reason: "repeat endpoint atom index is invalid".into(),
                        });
                    }
                    (AtomIdx(left), AtomIdx(right), None)
                };
                let mut remap = BTreeMap::new();
                for (idx, atom) in fragment.atoms() {
                    if excluded.is_none_or(|(a, b)| idx != a && idx != b) {
                        let added = molecule.add_atom(atom.clone());
                        unit_atoms.push(added);
                        remap.insert(idx, added);
                    }
                }
                for (_, bond) in fragment.bonds() {
                    if remap.contains_key(&bond.atom1) && remap.contains_key(&bond.atom2) {
                        molecule
                            .add_bond(remap[&bond.atom1], remap[&bond.atom2], bond.order)
                            .map_err(|e| SemanticError::InvalidExpansion {
                                id: unit.id.clone(),
                                reason: e.to_string(),
                            })?;
                    }
                }
                let chain_left = previous_right.unwrap_or(left);
                molecule
                    .add_bond(
                        chain_left,
                        remap[&endpoint_left],
                        chematic_core::BondOrder::Single,
                    )
                    .map_err(|e| SemanticError::InvalidExpansion {
                        id: unit.id.clone(),
                        reason: e.to_string(),
                    })?;
                let current_right = remap[&endpoint_right];
                if repeat_index + 1 == repeats {
                    molecule
                        .add_bond(right, current_right, chematic_core::BondOrder::Single)
                        .map_err(|e| SemanticError::InvalidExpansion {
                            id: unit.id.clone(),
                            reason: e.to_string(),
                        })?;
                } else {
                    previous_right = Some(current_right);
                }
            }
            mapping.insert(unit.id.clone(), unit_atoms);
        }
        Ok(ExpandedSemantic {
            molecule,
            source_to_expanded: mapping,
        })
    }
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
            if unit.repeat_count.is_none() {
                return Err(SemanticError::Unsupported {
                    construct: unit.id.clone(),
                    reason: "repeat count must be explicit".into(),
                });
            }
            if let Some(pattern) = unit.repeat_smiles.as_deref() {
                let pattern = pattern.trim();
                let has_placeholders = pattern.starts_with("[*]")
                    && pattern.ends_with("[*]")
                    && pattern.matches("[*]").count() == 2;
                if !has_placeholders && unit.repeat_endpoint_atoms.is_none() {
                    return Err(SemanticError::InvalidExpansion {
                        id: unit.id.clone(),
                        reason: "repeat requires two [*] markers or explicit endpoint atoms".into(),
                    });
                }
                if has_placeholders && pattern.matches("[*]").count() != 2 {
                    return Err(SemanticError::InvalidExpansion {
                        id: unit.id.clone(),
                        reason: "repeat_smiles must contain exactly two [*] markers".into(),
                    });
                }
            }
            if let Some([left, right]) = unit.repeat_endpoint_atoms
                && left == right
            {
                return Err(SemanticError::InvalidExpansion {
                    id: unit.id.clone(),
                    reason: "repeat endpoint atoms must be distinct".into(),
                });
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
            "end_groups": u.end_groups, "repeat_count": u.repeat_count,
            "repeat_smiles": u.repeat_smiles, "repeat_endpoint_atoms": u.repeat_endpoint_atoms
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
                repeat_smiles: None,
                repeat_endpoint_atoms: None,
            }],
            ..Default::default()
        };
        assert!(matches!(
            model.validate(),
            Err(SemanticError::AmbiguousAttachment(_))
        ));
    }

    #[test]
    fn command_selects_explicit_r_group_and_expands_with_mapping() {
        let base = chematic_smiles::parse("CC").unwrap();
        let model = SemanticModel {
            atom_ids: vec!["a1".into(), "a2".into()],
            r_groups: vec![RGroupDefinition {
                id: "r1".into(),
                attachment_atoms: vec![AtomRef {
                    atom_id: "a2".into(),
                }],
                alternatives: vec!["[*]O".into()],
                selected_alternative: None,
            }],
            ..Default::default()
        };
        let selected = model
            .apply(&SemanticCommand::SelectRGroupAlternative {
                group_id: "r1".into(),
                alternative: 0,
            })
            .unwrap();
        let expanded = selected.expand(&base).unwrap();
        assert_eq!(expanded.molecule.atom_count(), 3);
        assert_eq!(expanded.source_to_expanded["r1"].len(), 1);
    }

    #[test]
    fn expands_explicit_two_ended_polymer_repeat() {
        let base = chematic_smiles::parse("CC").unwrap();
        let model = SemanticModel {
            atom_ids: vec!["a1".into(), "a2".into()],
            polymer_units: vec![PolymerRepeatUnit {
                id: "p1".into(),
                attachment_atoms: vec![
                    AtomRef {
                        atom_id: "a1".into(),
                    },
                    AtomRef {
                        atom_id: "a2".into(),
                    },
                ],
                end_groups: vec![],
                repeat_count: Some(2),
                repeat_smiles: Some("[*]CC[*]".into()),
                repeat_endpoint_atoms: None,
            }],
            ..Default::default()
        };
        let expanded = model.expand(&base).unwrap();
        assert_eq!(expanded.molecule.atom_count(), 6);
        assert_eq!(expanded.source_to_expanded["p1"].len(), 4);
    }

    #[test]
    fn expands_explicit_polymer_endpoints_without_placeholders() {
        let base = chematic_smiles::parse("CC").unwrap();
        let model = SemanticModel {
            atom_ids: vec!["a1".into(), "a2".into()],
            polymer_units: vec![PolymerRepeatUnit {
                id: "p1".into(),
                attachment_atoms: vec![
                    AtomRef {
                        atom_id: "a1".into(),
                    },
                    AtomRef {
                        atom_id: "a2".into(),
                    },
                ],
                end_groups: vec![],
                repeat_count: Some(2),
                repeat_smiles: Some("CCO".into()),
                repeat_endpoint_atoms: Some([0, 2]),
            }],
            ..Default::default()
        };
        let expanded = model.expand(&base).unwrap();
        assert_eq!(expanded.molecule.atom_count(), 8);
        assert_eq!(expanded.source_to_expanded["p1"].len(), 6);
    }
}
