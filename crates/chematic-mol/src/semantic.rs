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
    /// Decode the stable `chematic.semantic.v1` JSON representation.
    pub fn from_json(value: &Value) -> Result<Self, SemanticError> {
        if value.get("schema").and_then(Value::as_str) != Some("chematic.semantic.v1") {
            return Err(SemanticError::InvalidJson(
                "schema must be chematic.semantic.v1".into(),
            ));
        }
        let strings = |key: &str| -> Result<Vec<String>, SemanticError> {
            value
                .get(key)
                .and_then(Value::as_array)
                .ok_or_else(|| SemanticError::InvalidJson(format!("{key} must be an array")))?
                .iter()
                .map(|item| {
                    item.as_str().map(str::to_owned).ok_or_else(|| {
                        SemanticError::InvalidJson(format!("{key} entries must be strings"))
                    })
                })
                .collect()
        };
        let atom_ids = strings("atom_ids")?;
        let bond_ids = strings("bond_ids")?;
        let r_groups = value
            .get("r_groups")
            .and_then(Value::as_array)
            .ok_or_else(|| SemanticError::InvalidJson("r_groups must be an array".into()))?
            .iter()
            .map(|group| {
                let id = json_string(group, "id")?;
                let attachment_atoms = json_string_array(group, "attachment_atoms")?
                    .into_iter()
                    .map(|atom_id| AtomRef { atom_id })
                    .collect();
                let alternatives = json_string_array(group, "alternatives")?;
                let selected_alternative = match group.get("selected_alternative") {
                    None | Some(Value::Null) => None,
                    Some(value) => {
                        Some(value.as_u64().map(|index| index as usize).ok_or_else(|| {
                            SemanticError::InvalidJson(
                                "selected_alternative must be an integer or null".into(),
                            )
                        })?)
                    }
                };
                Ok(RGroupDefinition {
                    id,
                    attachment_atoms,
                    alternatives,
                    selected_alternative,
                })
            })
            .collect::<Result<Vec<_>, SemanticError>>()?;
        let polymer_units = value
            .get("polymer_units")
            .and_then(Value::as_array)
            .ok_or_else(|| SemanticError::InvalidJson("polymer_units must be an array".into()))?
            .iter()
            .map(|unit| {
                let id = json_string(unit, "id")?;
                let attachment_atoms = json_string_array(unit, "attachment_atoms")?
                    .into_iter()
                    .map(|atom_id| AtomRef { atom_id })
                    .collect();
                let end_groups = json_string_array(unit, "end_groups")?;
                let repeat_count = match unit.get("repeat_count") {
                    None | Some(Value::Null) => None,
                    Some(value) => {
                        Some(value.as_u64().map(|count| count as u32).ok_or_else(|| {
                            SemanticError::InvalidJson(
                                "repeat_count must be an integer or null".into(),
                            )
                        })?)
                    }
                };
                let repeat_smiles = unit
                    .get("repeat_smiles")
                    .and_then(Value::as_str)
                    .map(str::to_owned);
                let repeat_endpoint_atoms = match unit.get("repeat_endpoint_atoms") {
                    None | Some(Value::Null) => None,
                    Some(value) => Some({
                        let values = value.as_array().ok_or_else(|| {
                            SemanticError::InvalidJson(
                                "repeat_endpoint_atoms must be an array or null".into(),
                            )
                        })?;
                        if values.len() != 2 {
                            return Err(SemanticError::InvalidJson(
                                "repeat_endpoint_atoms must contain two indices".into(),
                            ));
                        }
                        Ok([
                            values[0].as_u64().ok_or_else(|| {
                                SemanticError::InvalidJson(
                                    "endpoint indices must be integers".into(),
                                )
                            })? as u32,
                            values[1].as_u64().ok_or_else(|| {
                                SemanticError::InvalidJson(
                                    "endpoint indices must be integers".into(),
                                )
                            })? as u32,
                        ])
                    }?),
                };
                Ok(PolymerRepeatUnit {
                    id,
                    attachment_atoms,
                    end_groups,
                    repeat_count,
                    repeat_smiles,
                    repeat_endpoint_atoms,
                })
            })
            .collect::<Result<Vec<_>, SemanticError>>()?;
        let extensions = value
            .get("extensions")
            .and_then(Value::as_object)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect();
        let model = Self {
            atom_ids,
            bond_ids,
            r_groups,
            polymer_units,
            extensions,
        };
        model.validate()?;
        Ok(model)
    }

    /// Decode and apply a JSON command using the same contract as the Rust API.
    pub fn apply_json_command(&self, value: &Value) -> Result<Self, SemanticError> {
        let group_id = value
            .get("group_id")
            .and_then(Value::as_str)
            .ok_or_else(|| SemanticError::InvalidJson("group_id must be a string".into()))?;
        let alternative = value
            .get("alternative")
            .and_then(Value::as_u64)
            .ok_or_else(|| SemanticError::InvalidJson("alternative must be an integer".into()))?
            as usize;
        self.apply(&SemanticCommand::SelectRGroupAlternative {
            group_id: group_id.into(),
            alternative,
        })
    }

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
            let wildcards = fragment
                .atoms()
                .filter_map(|(idx, atom)| atom.wildcard.then_some(idx))
                .collect::<Vec<_>>();
            if fragment.atom_count() < 2 || wildcards.len() != group.attachment_atoms.len() {
                return Err(SemanticError::InvalidExpansion {
                    id: group.id.clone(),
                    reason: format!(
                        "alternative has {} wildcard attachment markers for {} source attachments",
                        wildcards.len(),
                        group.attachment_atoms.len()
                    ),
                });
            }
            let base_atoms = group
                .attachment_atoms
                .iter()
                .map(|reference| {
                    self.atom_ids
                        .iter()
                        .position(|id| id == &reference.atom_id)
                        .map(|index| AtomIdx(index as u32))
                        .ok_or_else(|| SemanticError::MissingAtom(reference.atom_id.clone()))
                })
                .collect::<Result<Vec<_>, _>>()?;
            for wildcard in &wildcards {
                let neighbors = fragment.neighbors(*wildcard).collect::<Vec<_>>();
                if neighbors.len() != 1 {
                    return Err(SemanticError::AmbiguousAttachment(group.id.clone()));
                }
            }
            let mut remap = BTreeMap::new();
            // Wildcard atoms are linkage markers, never copied into the result.
            for (idx, atom) in fragment.atoms() {
                if !atom.wildcard {
                    remap.insert(idx, molecule.add_atom(atom.clone()));
                }
            }
            for (_, bond) in fragment.bonds() {
                if !remap.contains_key(&bond.atom1) || !remap.contains_key(&bond.atom2) {
                    continue;
                }
                molecule
                    .add_bond(remap[&bond.atom1], remap[&bond.atom2], bond.order)
                    .map_err(|e| SemanticError::InvalidExpansion {
                        id: group.id.clone(),
                        reason: e.to_string(),
                    })?;
            }
            for (wildcard, base_atom) in wildcards.iter().zip(base_atoms) {
                let attach = remap[&fragment.neighbors(*wildcard).next().unwrap().0];
                molecule
                    .add_bond(base_atom, attach, chematic_core::BondOrder::Single)
                    .map_err(|e| SemanticError::InvalidExpansion {
                        id: group.id.clone(),
                        reason: e.to_string(),
                    })?;
            }
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
            for (side, end_group) in unit.end_groups.iter().enumerate() {
                let fragment = chematic_smiles::parse(end_group).map_err(|e| {
                    SemanticError::InvalidExpansion {
                        id: unit.id.clone(),
                        reason: format!("invalid end-group {side}: {e}"),
                    }
                })?;
                let wildcards = fragment
                    .atoms()
                    .filter_map(|(idx, atom)| atom.wildcard.then_some(idx))
                    .collect::<Vec<_>>();
                if fragment.atom_count() < 2 || wildcards.len() != 1 {
                    return Err(SemanticError::InvalidExpansion {
                        id: unit.id.clone(),
                        reason: format!("end-group {side} requires exactly one [*] marker"),
                    });
                }
                let wildcard = wildcards[0];
                let neighbors = fragment.neighbors(wildcard).collect::<Vec<_>>();
                if neighbors.len() != 1 {
                    return Err(SemanticError::AmbiguousAttachment(unit.id.clone()));
                }
                let mut remap = BTreeMap::new();
                for (idx, atom) in fragment.atoms() {
                    if !atom.wildcard {
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
                let endpoint = if side == 0 { left } else { right };
                molecule
                    .add_bond(
                        endpoint,
                        remap[&neighbors[0].0],
                        chematic_core::BondOrder::Single,
                    )
                    .map_err(|e| SemanticError::InvalidExpansion {
                        id: unit.id.clone(),
                        reason: e.to_string(),
                    })?;
                mapping.insert(
                    format!(
                        "{}.end_group_{}",
                        unit.id,
                        if side == 0 { "left" } else { "right" }
                    ),
                    remap.values().copied().collect(),
                );
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
            let mut attachments = std::collections::BTreeSet::new();
            for a in &group.attachment_atoms {
                if !self.atom_ids.contains(&a.atom_id) {
                    return Err(SemanticError::MissingAtom(a.atom_id.clone()));
                }
                if !attachments.insert(&a.atom_id) {
                    return Err(SemanticError::AmbiguousAttachment(group.id.clone()));
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
            if unit.repeat_count == Some(0) {
                return Err(SemanticError::InvalidExpansion {
                    id: unit.id.clone(),
                    reason: "repeat count must be greater than zero".into(),
                });
            }
            if !unit.end_groups.is_empty() && unit.end_groups.len() != 2 {
                return Err(SemanticError::InvalidExpansion {
                    id: unit.id.clone(),
                    reason: "end_groups must contain exactly [left, right] when provided".into(),
                });
            }
            for end_group in &unit.end_groups {
                if end_group.trim().is_empty() {
                    return Err(SemanticError::InvalidExpansion {
                        id: unit.id.clone(),
                        reason: "end-group SMILES must not be empty".into(),
                    });
                }
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

impl ExpandedSemantic {
    /// Serialize the expanded graph and its source mapping for binding use.
    pub fn to_json(&self) -> Value {
        serde_json::json!({
            "schema": "chematic.semantic-expanded.v1",
            "smiles": chematic_smiles::write(&self.molecule),
            "source_to_expanded": self.source_to_expanded.iter().map(|(id, atoms)|
                (id.clone(), atoms.iter().map(|atom| atom.0).collect::<Vec<_>>())
            ).collect::<BTreeMap<_, _>>(),
        })
    }
}

fn json_string(value: &Value, key: &str) -> Result<String, SemanticError> {
    value
        .get(key)
        .and_then(Value::as_str)
        .map(str::to_owned)
        .ok_or_else(|| SemanticError::InvalidJson(format!("{key} must be a string")))
}

fn json_string_array(value: &Value, key: &str) -> Result<Vec<String>, SemanticError> {
    value
        .get(key)
        .and_then(Value::as_array)
        .ok_or_else(|| SemanticError::InvalidJson(format!("{key} must be an array")))?
        .iter()
        .map(|item| {
            item.as_str()
                .map(str::to_owned)
                .ok_or_else(|| SemanticError::InvalidJson(format!("{key} entries must be strings")))
        })
        .collect()
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

    #[test]
    fn json_contract_round_trips_and_expands_selected_markush() {
        let base = chematic_smiles::parse("CC").unwrap();
        let model = SemanticModel {
            atom_ids: vec!["a1".into(), "a2".into()],
            r_groups: vec![RGroupDefinition {
                id: "r1".into(),
                attachment_atoms: vec![AtomRef {
                    atom_id: "a2".into(),
                }],
                alternatives: vec!["[*]O".into()],
                selected_alternative: Some(0),
            }],
            ..Default::default()
        };
        let decoded = SemanticModel::from_json(&model.to_json()).unwrap();
        let expanded = decoded.expand(&base).unwrap();
        assert_eq!(expanded.molecule.atom_count(), 3);
        assert_eq!(expanded.to_json()["source_to_expanded"]["r1"][0], 2);
    }

    #[test]
    fn expands_markush_with_multiple_attachment_points_in_marker_order() {
        let base = chematic_smiles::parse("CCCC").unwrap();
        let model = SemanticModel {
            atom_ids: vec!["a1".into(), "a2".into(), "a3".into(), "a4".into()],
            r_groups: vec![RGroupDefinition {
                id: "r1".into(),
                attachment_atoms: vec![
                    AtomRef {
                        atom_id: "a1".into(),
                    },
                    AtomRef {
                        atom_id: "a4".into(),
                    },
                ],
                alternatives: vec!["[*]O[*]".into()],
                selected_alternative: Some(0),
            }],
            ..Default::default()
        };
        let expanded = model.expand(&base).unwrap();
        assert_eq!(expanded.molecule.atom_count(), 5);
        assert_eq!(expanded.source_to_expanded["r1"].len(), 1);
    }

    #[test]
    fn expands_polymer_end_groups_with_stable_mapping() {
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
                end_groups: vec!["[*]O".into(), "[*]N".into()],
                repeat_count: Some(1),
                repeat_smiles: Some("[*]CC[*]".into()),
                repeat_endpoint_atoms: None,
            }],
            ..Default::default()
        };
        let expanded = model.expand(&base).unwrap();
        assert_eq!(expanded.molecule.atom_count(), 6);
        assert_eq!(expanded.source_to_expanded["p1"].len(), 4);
        assert_eq!(expanded.source_to_expanded["p1.end_group_left"].len(), 1);
        assert_eq!(expanded.source_to_expanded["p1.end_group_right"].len(), 1);
    }

    #[test]
    fn rejects_zero_repeat_count() {
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
                repeat_count: Some(0),
                repeat_smiles: Some("[*]CC[*]".into()),
                repeat_endpoint_atoms: None,
            }],
            ..Default::default()
        };
        assert!(matches!(
            model.validate(),
            Err(SemanticError::InvalidExpansion { .. })
        ));
    }
}
