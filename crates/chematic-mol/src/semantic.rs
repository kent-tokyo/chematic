//! Explicit semantic constructs that cannot be represented by an ordinary molecule.
//!
//! This module is intentionally conservative: ambiguous references and unsafe
//! expansions are rejected instead of being guessed into a molecular graph.

use std::collections::{BTreeMap, BTreeSet};

use serde_json::{Map, Value};

use chematic_core::{AtomIdx, Molecule};

mod expansion;

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
    /// Nested choices owned by this group. They are kept explicit rather than
    /// flattened into an ordinary molecule.
    pub nested_groups: Vec<RGroupDefinition>,
}

/// A polymer repeat unit with explicit linkage and end-group references.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolymerRepeatUnit {
    pub id: SemanticId,
    pub attachment_atoms: Vec<AtomRef>,
    pub end_groups: Vec<String>,
    /// Optional typed end-group definitions. When present, exactly two
    /// entries are required in left/right order and their IDs are retained in
    /// expansion provenance.
    pub end_group_definitions: Vec<PolymerEndGroup>,
    pub repeat_count: Option<u32>,
    /// Optional SMILES repeat fragment. It may use `[*]` at both ends, or use
    /// `repeat_endpoint_atoms` to identify two explicit endpoint atoms.
    pub repeat_smiles: Option<String>,
    /// Zero-based atom indices in `repeat_smiles` used when no `[*]` markers
    /// are present. Exactly two distinct endpoints are required.
    pub repeat_endpoint_atoms: Option<[u32; 2]>,
}

/// A stable, typed polymer end-group definition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolymerEndGroup {
    pub id: SemanticId,
    pub smiles: String,
}

/// A typed linkage between two source atoms in a polymer/S-group definition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolymerLinkage {
    pub id: SemanticId,
    pub left_atom: AtomRef,
    pub right_atom: AtomRef,
    pub bond_order: u8,
}

/// A loss-aware S-group definition. Unsupported S-group chemistry remains
/// explicit in `kind` instead of being flattened into the base molecule.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SemanticSGroup {
    pub id: SemanticId,
    pub kind: String,
    pub member_atoms: Vec<AtomRef>,
    pub repeat_unit_id: Option<SemanticId>,
    pub linkage: Option<PolymerLinkage>,
}

/// Loss/unsupported reason returned by validation or expansion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SemanticError {
    DuplicateId(String),
    MissingAtom(String),
    MissingAlternative(String),
    InvalidAlternative {
        id: String,
        reason: String,
    },
    AmbiguousAttachment(String),
    Unsupported {
        construct: String,
        reason: String,
    },
    InvalidJson(String),
    InvalidExpansion {
        id: String,
        reason: String,
    },
    ExpansionLimit {
        id: String,
        resource: &'static str,
        requested: usize,
        limit: usize,
    },
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
            Self::ExpansionLimit {
                id,
                resource,
                requested,
                limit,
            } => write!(
                f,
                "expansion {id} needs {requested} {resource}, limit is {limit}"
            ),
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
    pub s_groups: Vec<SemanticSGroup>,
    pub extensions: BTreeMap<String, Value>,
}

/// Result of a checked semantic expansion. The mapping is required for undo/edit flows.
#[derive(Clone)]
pub struct ExpandedSemantic {
    pub molecule: Molecule,
    pub source_to_expanded: BTreeMap<SemanticId, Vec<AtomIdx>>,
    base_molecule: Molecule,
}

/// Resource limits for semantic expansion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SemanticExpansionLimits {
    /// Maximum atom count in the expanded molecule.
    pub max_atoms: usize,
    /// Maximum repeat count for an individual polymer unit.
    pub max_repeat_count: u32,
}

impl Default for SemanticExpansionLimits {
    fn default() -> Self {
        Self {
            max_atoms: 100_000,
            max_repeat_count: 10_000,
        }
    }
}

/// Immutable, auditable edits to a semantic model.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SemanticCommand {
    SelectRGroupAlternative {
        group_id: SemanticId,
        alternative: usize,
    },
    /// Replace the explicitly allowed substituent set for an R-group.
    ReplaceRGroupAlternatives {
        group_id: SemanticId,
        alternatives: Vec<String>,
    },
    /// Clear a previously selected alternative so the model returns to an
    /// explicit, non-expandable Markush state.
    ClearRGroupAlternative { group_id: SemanticId },
    /// Change a typed S-group kind while retaining its stable identity.
    SetSGroupKind { group_id: SemanticId, kind: String },
    /// Clear a polymer repeat count so the model returns to an explicit,
    /// non-expandable editing state.
    ClearPolymerRepeatCount { unit_id: SemanticId },
    SetPolymerRepeatCount {
        unit_id: SemanticId,
        repeat_count: u32,
    },
}

fn parse_end_group_definitions(value: &Value) -> Result<Vec<PolymerEndGroup>, SemanticError> {
    value
        .get("end_group_definitions")
        .map(|definitions| {
            definitions
                .as_array()
                .ok_or_else(|| {
                    SemanticError::InvalidJson("end_group_definitions must be an array".into())
                })?
                .iter()
                .map(|definition| {
                    Ok(PolymerEndGroup {
                        id: json_string(definition, "id")?,
                        smiles: json_string(definition, "smiles")?,
                    })
                })
                .collect()
        })
        .transpose()
        .map(Option::unwrap_or_default)
}

fn parse_optional_repeat_count(value: &Value) -> Result<Option<u32>, SemanticError> {
    match value.get("repeat_count") {
        None | Some(Value::Null) => Ok(None),
        Some(count) => count
            .as_u64()
            .and_then(|count| u32::try_from(count).ok())
            .map(Some)
            .ok_or_else(|| {
                SemanticError::InvalidJson("repeat_count must be an integer or null".into())
            }),
    }
}

fn parse_optional_endpoints(value: &Value) -> Result<Option<[u32; 2]>, SemanticError> {
    let Some(endpoints) = value.get("repeat_endpoint_atoms") else {
        return Ok(None);
    };
    if endpoints.is_null() {
        return Ok(None);
    }
    let endpoints = endpoints.as_array().ok_or_else(|| {
        SemanticError::InvalidJson("repeat_endpoint_atoms must be an array or null".into())
    })?;
    if endpoints.len() != 2 {
        return Err(SemanticError::InvalidJson(
            "repeat_endpoint_atoms must contain two indices".into(),
        ));
    }
    let endpoint = |value: &Value| {
        value
            .as_u64()
            .and_then(|index| u32::try_from(index).ok())
            .ok_or_else(|| {
                SemanticError::InvalidJson("endpoint indices must be u32 integers".into())
            })
    };
    Ok(Some([endpoint(&endpoints[0])?, endpoint(&endpoints[1])?]))
}

fn parse_polymer_unit(value: &Value) -> Result<PolymerRepeatUnit, SemanticError> {
    Ok(PolymerRepeatUnit {
        id: json_string(value, "id")?,
        attachment_atoms: json_string_array(value, "attachment_atoms")?
            .into_iter()
            .map(|atom_id| AtomRef { atom_id })
            .collect(),
        end_groups: json_string_array(value, "end_groups")?,
        end_group_definitions: parse_end_group_definitions(value)?,
        repeat_count: parse_optional_repeat_count(value)?,
        repeat_smiles: value
            .get("repeat_smiles")
            .and_then(Value::as_str)
            .map(str::to_owned),
        repeat_endpoint_atoms: parse_optional_endpoints(value)?,
    })
}

fn parse_array<T>(
    value: &Value,
    key: &str,
    parse: impl FnMut(&Value) -> Result<T, SemanticError>,
) -> Result<Vec<T>, SemanticError> {
    value
        .get(key)
        .and_then(Value::as_array)
        .ok_or_else(|| SemanticError::InvalidJson(format!("{key} must be an array")))?
        .iter()
        .map(parse)
        .collect()
}

fn parse_optional_s_groups(value: &Value) -> Result<Vec<SemanticSGroup>, SemanticError> {
    value
        .get("s_groups")
        .map(|groups| {
            groups
                .as_array()
                .ok_or_else(|| SemanticError::InvalidJson("s_groups must be an array".into()))?
                .iter()
                .map(parse_s_group)
                .collect()
        })
        .transpose()
        .map(Option::unwrap_or_default)
}

impl SemanticModel {
    /// Decode the stable `chematic.semantic.v1` JSON representation.
    pub fn from_json(value: &Value) -> Result<Self, SemanticError> {
        if value.get("schema").and_then(Value::as_str) != Some("chematic.semantic.v1") {
            return Err(SemanticError::InvalidJson(
                "schema must be chematic.semantic.v1".into(),
            ));
        }
        if value.get("biomolecules").is_some() {
            return Err(SemanticError::Unsupported {
                construct: "biomolecule".into(),
                reason: "typed biomolecule semantics are not implemented".into(),
            });
        }
        let model = Self {
            atom_ids: json_string_array(value, "atom_ids")?,
            bond_ids: json_string_array(value, "bond_ids")?,
            r_groups: parse_array(value, "r_groups", parse_r_group)?,
            polymer_units: parse_array(value, "polymer_units", parse_polymer_unit)?,
            s_groups: parse_optional_s_groups(value)?,
            extensions: value
                .get("extensions")
                .and_then(Value::as_object)
                .cloned()
                .unwrap_or_default()
                .into_iter()
                .collect(),
        };
        model.validate()?;
        Ok(model)
    }

    /// Decode and apply a JSON command using the same contract as the Rust API.
    pub fn apply_json_command(&self, value: &Value) -> Result<Self, SemanticError> {
        if let Some(group_id) = value.get("replace_group_id").and_then(Value::as_str) {
            let alternatives = value
                .get("alternatives")
                .and_then(Value::as_array)
                .ok_or_else(|| SemanticError::InvalidJson("alternatives must be an array".into()))?
                .iter()
                .map(|alternative| {
                    alternative.as_str().map(str::to_owned).ok_or_else(|| {
                        SemanticError::InvalidJson("alternatives entries must be strings".into())
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            return self.apply(&SemanticCommand::ReplaceRGroupAlternatives {
                group_id: group_id.into(),
                alternatives,
            });
        }
        if let Some(group_id) = value.get("s_group_id").and_then(Value::as_str) {
            let kind = value
                .get("kind")
                .and_then(Value::as_str)
                .ok_or_else(|| SemanticError::InvalidJson("kind must be a string".into()))?;
            return self.apply(&SemanticCommand::SetSGroupKind {
                group_id: group_id.into(),
                kind: kind.into(),
            });
        }
        if let Some(group_id) = value.get("clear_group_id").and_then(Value::as_str) {
            return self.apply(&SemanticCommand::ClearRGroupAlternative {
                group_id: group_id.into(),
            });
        }
        if let Some(unit_id) = value.get("clear_unit_id").and_then(Value::as_str) {
            return self.apply(&SemanticCommand::ClearPolymerRepeatCount {
                unit_id: unit_id.into(),
            });
        }
        if let Some(unit_id) = value.get("unit_id").and_then(Value::as_str) {
            let repeat_count = value
                .get("repeat_count")
                .and_then(Value::as_u64)
                .and_then(|count| u32::try_from(count).ok())
                .ok_or_else(|| {
                    SemanticError::InvalidJson("repeat_count must be a positive u32 integer".into())
                })?;
            return self.apply(&SemanticCommand::SetPolymerRepeatCount {
                unit_id: unit_id.into(),
                repeat_count,
            });
        }
        let group_id = value
            .get("group_id")
            .and_then(Value::as_str)
            .ok_or_else(|| SemanticError::InvalidJson("group_id must be a string".into()))?;
        let alternative = value
            .get("alternative")
            .and_then(Value::as_u64)
            .and_then(|index| usize::try_from(index).ok())
            .ok_or_else(|| SemanticError::InvalidJson("alternative must be an integer".into()))?;
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
                let group = find_r_group_mut(&mut next.r_groups, group_id).ok_or_else(|| {
                    SemanticError::InvalidExpansion {
                        id: group_id.clone(),
                        reason: "unknown R-group".into(),
                    }
                })?;
                if *alternative >= group.alternatives.len() {
                    return Err(SemanticError::MissingAlternative(group_id.clone()));
                }
                group.selected_alternative = Some(*alternative);
            }
            SemanticCommand::ReplaceRGroupAlternatives {
                group_id,
                alternatives,
            } => {
                let group = find_r_group_mut(&mut next.r_groups, group_id).ok_or_else(|| {
                    SemanticError::InvalidExpansion {
                        id: group_id.clone(),
                        reason: "unknown R-group".into(),
                    }
                })?;
                group.alternatives = alternatives.clone();
                group.selected_alternative = None;
            }
            SemanticCommand::ClearRGroupAlternative { group_id } => {
                let group = find_r_group_mut(&mut next.r_groups, group_id).ok_or_else(|| {
                    SemanticError::InvalidExpansion {
                        id: group_id.clone(),
                        reason: "unknown R-group".into(),
                    }
                })?;
                group.selected_alternative = None;
            }
            SemanticCommand::SetSGroupKind { group_id, kind } => {
                let group = next
                    .s_groups
                    .iter_mut()
                    .find(|group| &group.id == group_id)
                    .ok_or_else(|| SemanticError::InvalidExpansion {
                        id: group_id.clone(),
                        reason: "unknown S-group".into(),
                    })?;
                group.kind = kind.clone();
            }
            SemanticCommand::ClearPolymerRepeatCount { unit_id } => {
                let unit = next
                    .polymer_units
                    .iter_mut()
                    .find(|unit| &unit.id == unit_id)
                    .ok_or_else(|| SemanticError::InvalidExpansion {
                        id: unit_id.clone(),
                        reason: "unknown polymer unit".into(),
                    })?;
                unit.repeat_count = None;
            }
            SemanticCommand::SetPolymerRepeatCount {
                unit_id,
                repeat_count,
            } => {
                let unit = next
                    .polymer_units
                    .iter_mut()
                    .find(|unit| &unit.id == unit_id)
                    .ok_or_else(|| SemanticError::InvalidExpansion {
                        id: unit_id.clone(),
                        reason: "unknown polymer unit".into(),
                    })?;
                if *repeat_count == 0 {
                    return Err(SemanticError::InvalidExpansion {
                        id: unit_id.clone(),
                        reason: "repeat count must be greater than zero".into(),
                    });
                }
                unit.repeat_count = Some(*repeat_count);
            }
        }
        next.validate()?;
        Ok(next)
    }

    /// Expand only explicitly selected R-groups. No alternative is guessed.
    /// Supported alternatives use a leading `[*]` attachment placeholder.
    pub fn expand(&self, base: &Molecule) -> Result<ExpandedSemantic, SemanticError> {
        self.expand_with_limits(base, &SemanticExpansionLimits::default())
    }

    /// Expand with an explicit atom and repeat budget.
    pub fn expand_with_limits(
        &self,
        base: &Molecule,
        limits: &SemanticExpansionLimits,
    ) -> Result<ExpandedSemantic, SemanticError> {
        expansion::expand(self, base, limits)
    }
}

fn validate_end_groups(
    unit: &PolymerRepeatUnit,
    ids: &mut BTreeSet<SemanticId>,
) -> Result<(), SemanticError> {
    if !unit.end_groups.is_empty() && unit.end_groups.len() != 2 {
        return Err(SemanticError::InvalidExpansion {
            id: unit.id.clone(),
            reason: "end_groups must contain exactly [left, right] when provided".into(),
        });
    }
    if !unit.end_group_definitions.is_empty() {
        if !unit.end_groups.is_empty() || unit.end_group_definitions.len() != 2 {
            return Err(SemanticError::InvalidExpansion {
                id: unit.id.clone(),
                reason:
                    "use exactly two typed end_group_definitions or legacy end_groups, not both"
                        .into(),
            });
        }
        let mut local_ids = BTreeSet::new();
        for definition in &unit.end_group_definitions {
            if definition.id.trim().is_empty() || !local_ids.insert(&definition.id) {
                return Err(SemanticError::DuplicateId(definition.id.clone()));
            }
            if !ids.insert(definition.id.clone()) {
                return Err(SemanticError::DuplicateId(definition.id.clone()));
            }
            if definition.smiles.trim().is_empty() {
                return Err(SemanticError::InvalidExpansion {
                    id: definition.id.clone(),
                    reason: "end-group SMILES must not be empty".into(),
                });
            }
        }
    }
    if unit.end_groups.iter().any(|group| group.trim().is_empty()) {
        return Err(SemanticError::InvalidExpansion {
            id: unit.id.clone(),
            reason: "end-group SMILES must not be empty".into(),
        });
    }
    Ok(())
}

fn validate_repeat_pattern(unit: &PolymerRepeatUnit) -> Result<(), SemanticError> {
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
    Ok(())
}

fn validate_polymer_unit(
    unit: &PolymerRepeatUnit,
    atom_ids: &[SemanticId],
    ids: &mut BTreeSet<SemanticId>,
) -> Result<(), SemanticError> {
    if !ids.insert(unit.id.clone()) {
        return Err(SemanticError::DuplicateId(unit.id.clone()));
    }
    if unit.attachment_atoms.len() != 2 {
        return Err(SemanticError::AmbiguousAttachment(unit.id.clone()));
    }
    if unit.repeat_count == Some(0) {
        return Err(SemanticError::InvalidExpansion {
            id: unit.id.clone(),
            reason: "repeat count must be greater than zero".into(),
        });
    }
    validate_end_groups(unit, ids)?;
    validate_repeat_pattern(unit)?;
    for atom in &unit.attachment_atoms {
        if !atom_ids.contains(&atom.atom_id) {
            return Err(SemanticError::MissingAtom(atom.atom_id.clone()));
        }
    }
    Ok(())
}

fn validate_s_group(
    group: &SemanticSGroup,
    model: &SemanticModel,
    ids: &mut BTreeSet<SemanticId>,
) -> Result<(), SemanticError> {
    if !ids.insert(group.id.clone()) {
        return Err(SemanticError::DuplicateId(group.id.clone()));
    }
    if group.kind.trim().is_empty() {
        return Err(SemanticError::Unsupported {
            construct: "S-group".into(),
            reason: format!("S-group {} has an empty kind", group.id),
        });
    }
    if group.member_atoms.is_empty() {
        return Err(SemanticError::InvalidExpansion {
            id: group.id.clone(),
            reason: "S-group must contain at least one member atom".into(),
        });
    }
    for atom in &group.member_atoms {
        if !model.atom_ids.contains(&atom.atom_id) {
            return Err(SemanticError::MissingAtom(atom.atom_id.clone()));
        }
    }
    if let Some(unit_id) = &group.repeat_unit_id
        && !model.polymer_units.iter().any(|unit| &unit.id == unit_id)
    {
        return Err(SemanticError::InvalidExpansion {
            id: group.id.clone(),
            reason: format!("unknown repeat unit reference: {unit_id}"),
        });
    }
    let Some(linkage) = &group.linkage else {
        return Ok(());
    };
    if linkage.id.trim().is_empty() || !ids.insert(linkage.id.clone()) {
        return Err(SemanticError::DuplicateId(linkage.id.clone()));
    }
    for atom in [&linkage.left_atom, &linkage.right_atom] {
        if !model.atom_ids.contains(&atom.atom_id) {
            return Err(SemanticError::MissingAtom(atom.atom_id.clone()));
        }
    }
    if linkage.left_atom.atom_id == linkage.right_atom.atom_id {
        return Err(SemanticError::AmbiguousAttachment(linkage.id.clone()));
    }
    if !matches!(linkage.bond_order, 1..=4) {
        return Err(SemanticError::Unsupported {
            construct: "polymer linkage".into(),
            reason: format!("unsupported bond order: {}", linkage.bond_order),
        });
    }
    Ok(())
}

impl SemanticModel {
    /// Validate IDs, references, attachment arity, and SMARTS-like alternatives.
    pub fn validate(&self) -> Result<(), SemanticError> {
        let mut ids = BTreeSet::new();
        for id in self.atom_ids.iter().chain(self.bond_ids.iter()) {
            if !ids.insert(id.clone()) {
                return Err(SemanticError::DuplicateId(id.clone()));
            }
        }
        for group in &self.r_groups {
            validate_r_group(group, &self.atom_ids, &mut ids)?;
        }
        for unit in &self.polymer_units {
            validate_polymer_unit(unit, &self.atom_ids, &mut ids)?;
        }
        for group in &self.s_groups {
            validate_s_group(group, self, &mut ids)?;
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
        root.insert(
            "r_groups".into(),
            Value::Array(self.r_groups.iter().map(r_group_to_json).collect()),
        );
        root.insert("polymer_units".into(), Value::Array(self.polymer_units.iter().map(|u| serde_json::json!({
            "id": u.id, "attachment_atoms": u.attachment_atoms.iter().map(|a| &a.atom_id).collect::<Vec<_>>(),
            "end_groups": u.end_groups, "repeat_count": u.repeat_count,
            "end_group_definitions": u.end_group_definitions.iter().map(|definition| serde_json::json!({
                "id": definition.id, "smiles": definition.smiles
            })).collect::<Vec<_>>(),
            "repeat_smiles": u.repeat_smiles, "repeat_endpoint_atoms": u.repeat_endpoint_atoms
        })).collect()));
        root.insert(
            "s_groups".into(),
            Value::Array(self.s_groups.iter().map(s_group_to_json).collect()),
        );
        root.insert(
            "extensions".into(),
            Value::Object(self.extensions.clone().into_iter().collect()),
        );
        Value::Object(root)
    }
}

fn validate_r_group(
    group: &RGroupDefinition,
    atom_ids: &[SemanticId],
    ids: &mut std::collections::BTreeSet<SemanticId>,
) -> Result<(), SemanticError> {
    if !ids.insert(group.id.clone()) {
        return Err(SemanticError::DuplicateId(group.id.clone()));
    }
    if group.attachment_atoms.is_empty() {
        return Err(SemanticError::AmbiguousAttachment(group.id.clone()));
    }
    let mut attachments = std::collections::BTreeSet::new();
    for a in &group.attachment_atoms {
        if !atom_ids.contains(&a.atom_id) {
            return Err(SemanticError::MissingAtom(a.atom_id.clone()));
        }
        if !attachments.insert(&a.atom_id) {
            return Err(SemanticError::AmbiguousAttachment(group.id.clone()));
        }
    }
    if group.alternatives.is_empty() {
        return Err(SemanticError::InvalidAlternative {
            id: group.id.clone(),
            reason: "allowed substituent set must not be empty".into(),
        });
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
    for nested in &group.nested_groups {
        validate_r_group(nested, atom_ids, ids)?;
    }
    Ok(())
}

fn parse_s_group(group: &Value) -> Result<SemanticSGroup, SemanticError> {
    let member_atoms = json_string_array(group, "member_atoms")?
        .into_iter()
        .map(|atom_id| AtomRef { atom_id })
        .collect();
    let linkage = group
        .get("linkage")
        .map(|value| {
            let bond_order = value
                .get("bond_order")
                .and_then(Value::as_u64)
                .and_then(|order| u8::try_from(order).ok())
                .ok_or_else(|| {
                    SemanticError::InvalidJson("linkage bond_order must be a u8".into())
                })?;
            Ok(PolymerLinkage {
                id: json_string(value, "id")?,
                left_atom: AtomRef {
                    atom_id: json_string(value, "left_atom")?,
                },
                right_atom: AtomRef {
                    atom_id: json_string(value, "right_atom")?,
                },
                bond_order,
            })
        })
        .transpose()?;
    Ok(SemanticSGroup {
        id: json_string(group, "id")?,
        kind: json_string(group, "kind")?,
        member_atoms,
        repeat_unit_id: group
            .get("repeat_unit_id")
            .and_then(Value::as_str)
            .map(str::to_owned),
        linkage,
    })
}

fn s_group_to_json(group: &SemanticSGroup) -> Value {
    serde_json::json!({
        "id": group.id,
        "kind": group.kind,
        "member_atoms": group.member_atoms.iter().map(|atom| &atom.atom_id).collect::<Vec<_>>(),
        "repeat_unit_id": group.repeat_unit_id,
        "linkage": group.linkage.as_ref().map(|linkage| serde_json::json!({
            "id": linkage.id,
            "left_atom": linkage.left_atom.atom_id,
            "right_atom": linkage.right_atom.atom_id,
            "bond_order": linkage.bond_order,
        })),
    })
}

fn find_r_group_mut<'a>(
    groups: &'a mut [RGroupDefinition],
    id: &str,
) -> Option<&'a mut RGroupDefinition> {
    for group in groups {
        if group.id == id {
            return Some(group);
        }
        if let Some(found) = find_r_group_mut(&mut group.nested_groups, id) {
            return Some(found);
        }
    }
    None
}

fn collect_r_groups<'a>(groups: &'a [RGroupDefinition], output: &mut Vec<&'a RGroupDefinition>) {
    for group in groups {
        output.push(group);
        collect_r_groups(&group.nested_groups, output);
    }
}

fn r_group_to_json(group: &RGroupDefinition) -> Value {
    serde_json::json!({
        "id": group.id,
        "attachment_atoms": group.attachment_atoms.iter().map(|a| &a.atom_id).collect::<Vec<_>>(),
        "alternatives": group.alternatives,
        "selected_alternative": group.selected_alternative,
        "nested_groups": group.nested_groups.iter().map(r_group_to_json).collect::<Vec<_>>(),
    })
}

impl ExpandedSemantic {
    /// Contract the expansion back to the exact base graph supplied to
    /// [`SemanticModel::expand`]. This is deliberately provenance-backed rather
    /// than inferred by deleting atoms from an arbitrary edited graph.
    pub fn contract(&self) -> Molecule {
        self.base_molecule.clone()
    }

    /// Serialize the expanded graph and its source mapping for binding use.
    pub fn to_json(&self) -> Value {
        serde_json::json!({
            "schema": "chematic.semantic-expanded.v1",
            "smiles": chematic_smiles::write(&self.molecule),
            "contracted_smiles": chematic_smiles::write(&self.base_molecule),
            "source_to_expanded": self.source_to_expanded.iter().map(|(id, atoms)|
                (id.clone(), atoms.iter().map(|atom| atom.0).collect::<Vec<_>>())
            ).collect::<BTreeMap<_, _>>(),
        })
    }
}

fn ensure_atom_budget(
    molecule: &Molecule,
    additional: usize,
    limits: &SemanticExpansionLimits,
    id: &str,
) -> Result<(), SemanticError> {
    let requested = molecule
        .atom_count()
        .checked_add(additional)
        .ok_or_else(|| SemanticError::ExpansionLimit {
            id: id.into(),
            resource: "atoms",
            requested: usize::MAX,
            limit: limits.max_atoms,
        })?;
    if requested > limits.max_atoms {
        return Err(SemanticError::ExpansionLimit {
            id: id.into(),
            resource: "atoms",
            requested,
            limit: limits.max_atoms,
        });
    }
    Ok(())
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

fn parse_r_group(group: &Value) -> Result<RGroupDefinition, SemanticError> {
    let selected_alternative = match group.get("selected_alternative") {
        None | Some(Value::Null) => None,
        Some(value) => Some(
            value
                .as_u64()
                .and_then(|index| usize::try_from(index).ok())
                .ok_or_else(|| {
                    SemanticError::InvalidJson(
                        "selected_alternative must be an integer or null".into(),
                    )
                })?,
        ),
    };
    let nested_groups = group
        .get("nested_groups")
        .map(|value| {
            value
                .as_array()
                .ok_or_else(|| SemanticError::InvalidJson("nested_groups must be an array".into()))?
                .iter()
                .map(parse_r_group)
                .collect()
        })
        .transpose()?
        .unwrap_or_default();
    Ok(RGroupDefinition {
        id: json_string(group, "id")?,
        attachment_atoms: json_string_array(group, "attachment_atoms")?
            .into_iter()
            .map(|atom_id| AtomRef { atom_id })
            .collect(),
        alternatives: json_string_array(group, "alternatives")?,
        selected_alternative,
        nested_groups,
    })
}

#[cfg(test)]
#[path = "semantic/tests.rs"]
mod tests;

#[cfg(test)]
mod expansion_boundary_contract_tests {
    use super::*;
    use chematic_smiles::{canonical_smiles, parse};

    fn polymer(repeats: u32) -> SemanticModel {
        SemanticModel {
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
                end_group_definitions: vec![],
                repeat_count: Some(repeats),
                repeat_smiles: Some("[*]CC[*]".into()),
                repeat_endpoint_atoms: None,
            }],
            ..Default::default()
        }
    }

    fn rejected_without_mutation(model: &SemanticModel, context: &str) -> SemanticError {
        let base = parse("CC").unwrap();
        let before = model.to_json();
        let error = model
            .expand(&base)
            .err()
            .expect("invalid expansion must fail");
        assert!(error.to_string().contains(context), "{context}: {error}");
        assert_eq!(model.to_json(), before);
        assert_eq!(canonical_smiles(&base), "CC");
        assert_eq!(base.atom_count(), 2);
        assert_eq!(base.bond_count(), 1);
        error
    }

    #[test]
    fn capped_repeat_expansion_preserves_identity_mapping_at_exact_atom_limits() {
        let base = parse("CC").unwrap();
        for (repeats, expected) in [(1, "NC1CCC1O"), (2, "NC1CCCCC1O"), (3, "NC1CCCCCCC1O")] {
            for typed in [false, true] {
                let mut model = polymer(repeats);
                let ids = if typed {
                    model.polymer_units[0].end_group_definitions = vec![
                        PolymerEndGroup {
                            id: "left-cap".into(),
                            smiles: "[*]O".into(),
                        },
                        PolymerEndGroup {
                            id: "right-cap".into(),
                            smiles: "[*]N".into(),
                        },
                    ];
                    ["left-cap", "right-cap"]
                } else {
                    model.polymer_units[0].end_groups = vec!["[*]O".into(), "[*]N".into()];
                    ["p1.end_group_left", "p1.end_group_right"]
                };
                let before = model.to_json();
                let total = 2 * repeats as usize + 4;
                let limits = SemanticExpansionLimits {
                    max_atoms: total,
                    max_repeat_count: repeats,
                };
                let result = model.expand_with_limits(&base, &limits).unwrap();
                assert_eq!(result.molecule.atom_count(), total);
                assert_eq!(
                    canonical_smiles(&result.molecule),
                    canonical_smiles(&parse(expected).unwrap())
                );
                assert_eq!(result.source_to_expanded["a1"], vec![AtomIdx(0)]);
                assert_eq!(result.source_to_expanded["a2"], vec![AtomIdx(1)]);
                assert_eq!(result.source_to_expanded["p1"].len(), total - 2);
                assert_eq!(
                    result.source_to_expanded[ids[0]],
                    vec![AtomIdx((total - 2) as u32)]
                );
                assert_eq!(
                    result.source_to_expanded[ids[1]],
                    vec![AtomIdx((total - 1) as u32)]
                );
                let error = model
                    .expand_with_limits(
                        &base,
                        &SemanticExpansionLimits {
                            max_atoms: total - 1,
                            ..limits
                        },
                    )
                    .err()
                    .unwrap();
                assert_eq!(
                    error,
                    SemanticError::ExpansionLimit {
                        id: "p1".into(),
                        resource: "atoms",
                        requested: total,
                        limit: total - 1
                    }
                );
                assert_eq!(model.to_json(), before);
                assert_eq!(canonical_smiles(&base), "CC");
            }
        }
    }

    #[test]
    fn malformed_repeat_and_cap_topologies_fail_before_returning_a_graph() {
        for (pattern, endpoints, context) in [
            ("CC", None, "two [*] markers or explicit endpoint atoms"),
            ("CC", Some([0, 0]), "must be distinct"),
            ("CC", Some([0, 2]), "endpoint atom index is invalid"),
            ("CC", Some([3, 1]), "endpoint atom index is invalid"),
            ("[*]C1C[*]", None, "unmatched ring closure"),
            ("[*]C[*]", None, "exactly one neighbor at each linkage"),
            ("[*](C)CC[*]", None, "exactly one neighbor at each linkage"),
        ] {
            let mut model = polymer(2);
            model.polymer_units[0].repeat_smiles = Some(pattern.into());
            model.polymer_units[0].repeat_endpoint_atoms = endpoints;
            rejected_without_mutation(&model, context);
        }
        for (cap, context) in [
            ("O", "exactly one [*] marker"),
            ("[*]", "exactly one [*] marker"),
            ("[*]O[*]", "exactly one [*] marker"),
            ("[*]C1", "invalid end-group 0"),
            ("[*](C)O", "ambiguous attachment topology: p1"),
        ] {
            let mut model = polymer(2);
            model.polymer_units[0].end_groups = vec![cap.into(), "[*]N".into()];
            rejected_without_mutation(&model, context);
        }
        let mut model = polymer(2);
        model.polymer_units[0].repeat_smiles = None;
        assert!(matches!(
            rejected_without_mutation(&model, "repeat_smiles is not provided"),
            SemanticError::Unsupported { .. }
        ));
        let mut model = polymer(2);
        model.polymer_units[0].repeat_count = None;
        assert!(matches!(
            rejected_without_mutation(&model, "repeat count must be explicit"),
            SemanticError::Unsupported { .. }
        ));
    }

    #[test]
    fn selected_r_group_alternatives_require_an_unambiguous_linkage() {
        for (pattern, context) in [
            ("C", "0 wildcard attachment markers for 1"),
            ("[*]", "alternative has"),
            ("[*]C[*]", "2 wildcard attachment markers for 1"),
            ("[*](C)O", "ambiguous attachment topology: r1"),
            ("[*]C1", "unmatched ring closure"),
        ] {
            let model = SemanticModel {
                atom_ids: vec!["a1".into(), "a2".into()],
                r_groups: vec![RGroupDefinition {
                    id: "r1".into(),
                    attachment_atoms: vec![AtomRef {
                        atom_id: "a2".into(),
                    }],
                    alternatives: vec![pattern.into()],
                    selected_alternative: Some(0),
                    nested_groups: vec![],
                }],
                ..Default::default()
            };
            rejected_without_mutation(&model, context);
        }
    }
}
