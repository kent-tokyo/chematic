//! Bounded, lossless nucleic-acid document metadata.
//!
//! This model deliberately does not infer residues from a molecular graph and
//! does not flatten strands into an ordinary [`chematic_core::Molecule`]. It
//! preserves an explicit editor/document contract and rejects ambiguous atom
//! ownership or unsupported linkage topology.

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

pub const NUCLEIC_ACID_SCHEMA: &str = "chematic.nucleic-acid.v1";
pub const NUCLEIC_ACID_VALIDATION_SCHEMA: &str = "chematic.nucleic-acid-validation.v1";

/// Resource limits shared by Rust, Python, and WASM entry points.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NucleicAcidLimits {
    pub max_input_bytes: usize,
    pub max_strands: usize,
    pub max_residues: usize,
    pub max_atom_ids: usize,
    pub max_atom_refs_per_residue: usize,
    pub max_linkages: usize,
    pub max_annotations: usize,
    pub max_annotation_bytes: usize,
    pub max_id_bytes: usize,
}

impl Default for NucleicAcidLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 1_048_576,
            max_strands: 64,
            max_residues: 10_000,
            max_atom_ids: 100_000,
            max_atom_refs_per_residue: 512,
            max_linkages: 10_000,
            max_annotations: 10_000,
            max_annotation_bytes: 262_144,
            max_id_bytes: 256,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NucleicAcidStrandKind {
    Dna,
    Rna,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum NucleicAcidBase {
    A,
    C,
    G,
    T,
    U,
    #[serde(rename = "other")]
    Other,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NucleicAcidSugar {
    Deoxyribose,
    Ribose,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NucleicAcidResidue {
    pub id: String,
    pub base: NucleicAcidBase,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub modification: Option<String>,
    pub sugar: NucleicAcidSugar,
    pub atom_refs: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub five_prime_linkage: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub three_prime_linkage: Option<String>,
    #[serde(default)]
    pub annotations: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NucleicAcidStrand {
    pub id: String,
    pub kind: NucleicAcidStrandKind,
    pub residues: Vec<NucleicAcidResidue>,
    #[serde(default)]
    pub annotations: BTreeMap<String, Value>,
}

/// One supported linear 3'-to-5' phosphodiester linkage.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NucleicAcidLinkage {
    pub id: String,
    pub kind: String,
    pub three_prime_residue_id: String,
    pub five_prime_residue_id: String,
    pub three_prime_atom_ref: String,
    pub five_prime_atom_ref: String,
    #[serde(default)]
    pub annotations: BTreeMap<String, Value>,
}

/// A source-level nucleic-acid document. Atom IDs are opaque references; no
/// graph-to-sequence inference is performed.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NucleicAcidDocument {
    pub schema: String,
    pub atom_ids: Vec<String>,
    pub strands: Vec<NucleicAcidStrand>,
    pub linkages: Vec<NucleicAcidLinkage>,
    #[serde(default)]
    pub annotations: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum NucleicAcidAnnotationTarget {
    Document,
    Strand { id: String },
    Residue { id: String },
    Linkage { id: String },
}

/// Bounded metadata-only edits. Linkage topology and atom ownership cannot be
/// mutated through this first-slice command API.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum NucleicAcidCommand {
    SetAnnotation {
        target: NucleicAcidAnnotationTarget,
        key: String,
        value: Value,
    },
    RemoveAnnotation {
        target: NucleicAcidAnnotationTarget,
        key: String,
    },
    SetResidueIdentity {
        residue_id: String,
        base: NucleicAcidBase,
        #[serde(default)]
        modification: Option<String>,
        sugar: NucleicAcidSugar,
    },
}

/// Stable typed validation categories.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NucleicAcidError {
    InvalidJson {
        message: String,
    },
    InvalidSchema {
        found: Option<String>,
    },
    InvalidId {
        path: String,
        reason: String,
    },
    DuplicateId {
        path: String,
        id: String,
    },
    MissingReference {
        path: String,
        id: String,
    },
    AmbiguousAtomMapping {
        path: String,
        atom_id: String,
        reason: String,
    },
    UnsupportedLinkageTopology {
        path: String,
        reason: String,
    },
    UnknownModification {
        path: String,
        modification: String,
    },
    InvalidResidueIdentity {
        path: String,
        reason: String,
    },
    ResourceLimit {
        path: String,
        resource: String,
        requested: usize,
        limit: usize,
    },
    InvalidCommand {
        message: String,
    },
    MissingTarget {
        path: String,
        id: String,
    },
}

impl NucleicAcidError {
    pub fn code(&self) -> &'static str {
        match self {
            Self::InvalidJson { .. } => "invalid_json",
            Self::InvalidSchema { .. } => "invalid_schema",
            Self::InvalidId { .. } => "invalid_id",
            Self::DuplicateId { .. } => "duplicate_id",
            Self::MissingReference { .. } => "missing_reference",
            Self::AmbiguousAtomMapping { .. } => "ambiguous_atom_mapping",
            Self::UnsupportedLinkageTopology { .. } => "unsupported_linkage_topology",
            Self::UnknownModification { .. } => "unknown_modification",
            Self::InvalidResidueIdentity { .. } => "invalid_residue_identity",
            Self::ResourceLimit { .. } => "resource_limit",
            Self::InvalidCommand { .. } => "invalid_command",
            Self::MissingTarget { .. } => "missing_target",
        }
    }

    pub fn path(&self) -> &str {
        match self {
            Self::InvalidJson { .. } | Self::InvalidSchema { .. } | Self::InvalidCommand { .. } => {
                "/"
            }
            Self::InvalidId { path, .. }
            | Self::DuplicateId { path, .. }
            | Self::MissingReference { path, .. }
            | Self::AmbiguousAtomMapping { path, .. }
            | Self::UnsupportedLinkageTopology { path, .. }
            | Self::UnknownModification { path, .. }
            | Self::InvalidResidueIdentity { path, .. }
            | Self::ResourceLimit { path, .. }
            | Self::MissingTarget { path, .. } => path,
        }
    }

    pub fn to_json(&self) -> Value {
        json!({
            "code": self.code(),
            "path": self.path(),
            "message": self.to_string(),
        })
    }
}

impl std::fmt::Display for NucleicAcidError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidJson { message } => write!(f, "invalid nucleic-acid JSON: {message}"),
            Self::InvalidSchema { found } => write!(
                f,
                "schema must be {NUCLEIC_ACID_SCHEMA}, found {}",
                found.as_deref().unwrap_or("missing")
            ),
            Self::InvalidId { path, reason } => write!(f, "invalid id at {path}: {reason}"),
            Self::DuplicateId { path, id } => write!(f, "duplicate id at {path}: {id}"),
            Self::MissingReference { path, id } => {
                write!(f, "missing reference at {path}: {id}")
            }
            Self::AmbiguousAtomMapping {
                path,
                atom_id,
                reason,
            } => write!(
                f,
                "ambiguous atom mapping at {path} for {atom_id}: {reason}"
            ),
            Self::UnsupportedLinkageTopology { path, reason } => {
                write!(f, "unsupported linkage topology at {path}: {reason}")
            }
            Self::UnknownModification { path, modification } => {
                write!(f, "unknown modification at {path}: {modification}")
            }
            Self::InvalidResidueIdentity { path, reason } => {
                write!(f, "invalid residue identity at {path}: {reason}")
            }
            Self::ResourceLimit {
                path,
                resource,
                requested,
                limit,
            } => write!(
                f,
                "resource limit at {path}: {resource} requested {requested}, limit is {limit}"
            ),
            Self::InvalidCommand { message } => {
                write!(f, "invalid nucleic-acid command: {message}")
            }
            Self::MissingTarget { path, id } => write!(f, "missing edit target at {path}: {id}"),
        }
    }
}

impl std::error::Error for NucleicAcidError {}

impl NucleicAcidDocument {
    pub fn from_json(value: &Value) -> Result<Self, NucleicAcidError> {
        Self::from_json_with_limits(value, &NucleicAcidLimits::default())
    }

    pub fn from_json_with_limits(
        value: &Value,
        limits: &NucleicAcidLimits,
    ) -> Result<Self, NucleicAcidError> {
        let document: Self = serde_json::from_value(value.clone()).map_err(|error| {
            NucleicAcidError::InvalidJson {
                message: error.to_string(),
            }
        })?;
        document.validate_with_limits(limits)?;
        Ok(document)
    }

    pub fn from_json_str(input: &str) -> Result<Self, NucleicAcidError> {
        Self::from_json_str_with_limits(input, &NucleicAcidLimits::default())
    }

    pub fn from_json_str_with_limits(
        input: &str,
        limits: &NucleicAcidLimits,
    ) -> Result<Self, NucleicAcidError> {
        check_limit("/", "input_bytes", input.len(), limits.max_input_bytes)?;
        let value: Value =
            serde_json::from_str(input).map_err(|error| NucleicAcidError::InvalidJson {
                message: error.to_string(),
            })?;
        Self::from_json_with_limits(&value, limits)
    }

    pub fn to_json(&self) -> Value {
        serde_json::to_value(self).expect("nucleic-acid document is JSON serializable")
    }

    pub fn validate(&self) -> Result<(), NucleicAcidError> {
        self.validate_with_limits(&NucleicAcidLimits::default())
    }

    pub fn validate_with_limits(&self, limits: &NucleicAcidLimits) -> Result<(), NucleicAcidError> {
        if self.schema != NUCLEIC_ACID_SCHEMA {
            return Err(NucleicAcidError::InvalidSchema {
                found: Some(self.schema.clone()),
            });
        }
        check_limit(
            "/strands",
            "strands",
            self.strands.len(),
            limits.max_strands,
        )?;
        check_limit(
            "/atom_ids",
            "atom_ids",
            self.atom_ids.len(),
            limits.max_atom_ids,
        )?;
        check_limit(
            "/linkages",
            "linkages",
            self.linkages.len(),
            limits.max_linkages,
        )?;
        if self.strands.is_empty() {
            return Err(NucleicAcidError::InvalidResidueIdentity {
                path: "/strands".into(),
                reason: "at least one strand is required".into(),
            });
        }

        let total_residues = self
            .strands
            .iter()
            .map(|strand| strand.residues.len())
            .sum();
        check_limit("/strands", "residues", total_residues, limits.max_residues)?;

        let mut global_ids = BTreeSet::new();
        let mut atom_ids = BTreeSet::new();
        for (index, id) in self.atom_ids.iter().enumerate() {
            validate_id(id, &format!("/atom_ids/{index}"), limits)?;
            if !atom_ids.insert(id.clone()) {
                return Err(NucleicAcidError::DuplicateId {
                    path: format!("/atom_ids/{index}"),
                    id: id.clone(),
                });
            }
        }

        let mut residue_locations = BTreeMap::<String, (usize, usize)>::new();
        let mut atom_owner = BTreeMap::<String, String>::new();
        for (strand_index, strand) in self.strands.iter().enumerate() {
            let strand_path = format!("/strands/{strand_index}");
            validate_id(&strand.id, &format!("{strand_path}/id"), limits)?;
            if !global_ids.insert(strand.id.clone()) {
                return Err(NucleicAcidError::DuplicateId {
                    path: format!("{strand_path}/id"),
                    id: strand.id.clone(),
                });
            }
            if strand.residues.is_empty() {
                return Err(NucleicAcidError::InvalidResidueIdentity {
                    path: format!("{strand_path}/residues"),
                    reason: "a strand must contain at least one residue".into(),
                });
            }
            for (residue_index, residue) in strand.residues.iter().enumerate() {
                let residue_path = format!("{strand_path}/residues/{residue_index}");
                validate_id(&residue.id, &format!("{residue_path}/id"), limits)?;
                if !global_ids.insert(residue.id.clone()) {
                    return Err(NucleicAcidError::DuplicateId {
                        path: format!("{residue_path}/id"),
                        id: residue.id.clone(),
                    });
                }
                residue_locations.insert(residue.id.clone(), (strand_index, residue_index));
                check_limit(
                    &format!("{residue_path}/atom_refs"),
                    "atom_refs_per_residue",
                    residue.atom_refs.len(),
                    limits.max_atom_refs_per_residue,
                )?;
                if residue.atom_refs.is_empty() {
                    return Err(NucleicAcidError::AmbiguousAtomMapping {
                        path: format!("{residue_path}/atom_refs"),
                        atom_id: String::new(),
                        reason: "a residue must reference at least one atom".into(),
                    });
                }
                validate_residue_identity(strand.kind, residue, &residue_path)?;
                for (atom_index, atom_id) in residue.atom_refs.iter().enumerate() {
                    let path = format!("{residue_path}/atom_refs/{atom_index}");
                    if !atom_ids.contains(atom_id) {
                        return Err(NucleicAcidError::MissingReference {
                            path,
                            id: atom_id.clone(),
                        });
                    }
                    if let Some(owner) = atom_owner.insert(atom_id.clone(), residue.id.clone()) {
                        return Err(NucleicAcidError::AmbiguousAtomMapping {
                            path,
                            atom_id: atom_id.clone(),
                            reason: format!("referenced by both {owner} and {}", residue.id),
                        });
                    }
                }
            }
        }

        for atom_id in &self.atom_ids {
            if !atom_owner.contains_key(atom_id) {
                return Err(NucleicAcidError::AmbiguousAtomMapping {
                    path: "/atom_ids".into(),
                    atom_id: atom_id.clone(),
                    reason: "atom is not assigned to any residue".into(),
                });
            }
        }

        let mut linkage_by_id = BTreeMap::new();
        for (index, linkage) in self.linkages.iter().enumerate() {
            let path = format!("/linkages/{index}");
            validate_id(&linkage.id, &format!("{path}/id"), limits)?;
            if !global_ids.insert(linkage.id.clone()) {
                return Err(NucleicAcidError::DuplicateId {
                    path: format!("{path}/id"),
                    id: linkage.id.clone(),
                });
            }
            if linkage.kind != "phosphodiester" {
                return Err(NucleicAcidError::UnsupportedLinkageTopology {
                    path: format!("{path}/kind"),
                    reason: format!("unsupported linkage kind {}", linkage.kind),
                });
            }
            for (field, residue_id) in [
                ("three_prime_residue_id", &linkage.three_prime_residue_id),
                ("five_prime_residue_id", &linkage.five_prime_residue_id),
            ] {
                if !residue_locations.contains_key(residue_id) {
                    return Err(NucleicAcidError::MissingReference {
                        path: format!("{path}/{field}"),
                        id: residue_id.clone(),
                    });
                }
            }
            for (field, atom_id, residue_id) in [
                (
                    "three_prime_atom_ref",
                    &linkage.three_prime_atom_ref,
                    &linkage.three_prime_residue_id,
                ),
                (
                    "five_prime_atom_ref",
                    &linkage.five_prime_atom_ref,
                    &linkage.five_prime_residue_id,
                ),
            ] {
                if atom_owner.get(atom_id) != Some(residue_id) {
                    return Err(NucleicAcidError::AmbiguousAtomMapping {
                        path: format!("{path}/{field}"),
                        atom_id: atom_id.clone(),
                        reason: format!("atom is not owned by residue {residue_id}"),
                    });
                }
            }
            linkage_by_id.insert(linkage.id.clone(), (index, linkage));
        }

        let expected_linkages: usize = self
            .strands
            .iter()
            .map(|strand| strand.residues.len().saturating_sub(1))
            .sum();
        if self.linkages.len() != expected_linkages {
            return Err(NucleicAcidError::UnsupportedLinkageTopology {
                path: "/linkages".into(),
                reason: format!(
                    "linear strands require {expected_linkages} linkages, found {}",
                    self.linkages.len()
                ),
            });
        }

        let mut used_linkages = BTreeSet::new();
        for (strand_index, strand) in self.strands.iter().enumerate() {
            for (residue_index, residue) in strand.residues.iter().enumerate() {
                let path = format!("/strands/{strand_index}/residues/{residue_index}");
                if residue_index == 0 && residue.five_prime_linkage.is_some() {
                    return Err(NucleicAcidError::UnsupportedLinkageTopology {
                        path: format!("{path}/five_prime_linkage"),
                        reason: "the first residue must have a free 5' terminus".into(),
                    });
                }
                if residue_index + 1 == strand.residues.len()
                    && residue.three_prime_linkage.is_some()
                {
                    return Err(NucleicAcidError::UnsupportedLinkageTopology {
                        path: format!("{path}/three_prime_linkage"),
                        reason: "the last residue must have a free 3' terminus".into(),
                    });
                }
                if residue_index + 1 < strand.residues.len() {
                    let next = &strand.residues[residue_index + 1];
                    let Some(linkage_id) = residue.three_prime_linkage.as_ref() else {
                        return Err(NucleicAcidError::UnsupportedLinkageTopology {
                            path: format!("{path}/three_prime_linkage"),
                            reason: "adjacent residues require an explicit linkage".into(),
                        });
                    };
                    if next.five_prime_linkage.as_ref() != Some(linkage_id) {
                        return Err(NucleicAcidError::UnsupportedLinkageTopology {
                            path: format!(
                                "/strands/{strand_index}/residues/{}/five_prime_linkage",
                                residue_index + 1
                            ),
                            reason: "adjacent residues must reference the same linkage".into(),
                        });
                    }
                    let Some((_, linkage)) = linkage_by_id.get(linkage_id) else {
                        return Err(NucleicAcidError::MissingReference {
                            path: format!("{path}/three_prime_linkage"),
                            id: linkage_id.clone(),
                        });
                    };
                    if linkage.three_prime_residue_id != residue.id
                        || linkage.five_prime_residue_id != next.id
                    {
                        return Err(NucleicAcidError::UnsupportedLinkageTopology {
                            path: format!("/linkages/{linkage_id}"),
                            reason: "linkage endpoints do not match adjacent strand order".into(),
                        });
                    }
                    used_linkages.insert(linkage_id.clone());
                }
            }
        }
        if used_linkages.len() != self.linkages.len() {
            return Err(NucleicAcidError::UnsupportedLinkageTopology {
                path: "/linkages".into(),
                reason: "cross-strand, cyclic, branched, or unreferenced linkage".into(),
            });
        }

        validate_annotations(self, limits)
    }

    pub fn apply(&self, command: &NucleicAcidCommand) -> Result<Self, NucleicAcidError> {
        self.apply_with_limits(command, &NucleicAcidLimits::default())
    }

    pub fn apply_with_limits(
        &self,
        command: &NucleicAcidCommand,
        limits: &NucleicAcidLimits,
    ) -> Result<Self, NucleicAcidError> {
        let mut next = self.clone();
        match command {
            NucleicAcidCommand::SetAnnotation { target, key, value } => {
                validate_annotation_key(key, limits)?;
                next.annotations_for_target_mut(target)?
                    .insert(key.clone(), value.clone());
            }
            NucleicAcidCommand::RemoveAnnotation { target, key } => {
                validate_annotation_key(key, limits)?;
                next.annotations_for_target_mut(target)?.remove(key);
            }
            NucleicAcidCommand::SetResidueIdentity {
                residue_id,
                base,
                modification,
                sugar,
            } => {
                let Some(residue) = next
                    .strands
                    .iter_mut()
                    .flat_map(|strand| strand.residues.iter_mut())
                    .find(|residue| residue.id == *residue_id)
                else {
                    return Err(NucleicAcidError::MissingTarget {
                        path: "/command/residue_id".into(),
                        id: residue_id.clone(),
                    });
                };
                residue.base = *base;
                residue.modification.clone_from(modification);
                residue.sugar = *sugar;
            }
        }
        next.validate_with_limits(limits)?;
        Ok(next)
    }

    pub fn apply_json_command_with_limits(
        &self,
        command: &Value,
        limits: &NucleicAcidLimits,
    ) -> Result<Self, NucleicAcidError> {
        let command: NucleicAcidCommand = serde_json::from_value(command.clone()).map_err(|e| {
            NucleicAcidError::InvalidCommand {
                message: e.to_string(),
            }
        })?;
        self.apply_with_limits(&command, limits)
    }

    fn annotations_for_target_mut(
        &mut self,
        target: &NucleicAcidAnnotationTarget,
    ) -> Result<&mut BTreeMap<String, Value>, NucleicAcidError> {
        match target {
            NucleicAcidAnnotationTarget::Document => Ok(&mut self.annotations),
            NucleicAcidAnnotationTarget::Strand { id } => self
                .strands
                .iter_mut()
                .find(|strand| strand.id == *id)
                .map(|strand| &mut strand.annotations)
                .ok_or_else(|| NucleicAcidError::MissingTarget {
                    path: "/command/target/id".into(),
                    id: id.clone(),
                }),
            NucleicAcidAnnotationTarget::Residue { id } => self
                .strands
                .iter_mut()
                .flat_map(|strand| strand.residues.iter_mut())
                .find(|residue| residue.id == *id)
                .map(|residue| &mut residue.annotations)
                .ok_or_else(|| NucleicAcidError::MissingTarget {
                    path: "/command/target/id".into(),
                    id: id.clone(),
                }),
            NucleicAcidAnnotationTarget::Linkage { id } => self
                .linkages
                .iter_mut()
                .find(|linkage| linkage.id == *id)
                .map(|linkage| &mut linkage.annotations)
                .ok_or_else(|| NucleicAcidError::MissingTarget {
                    path: "/command/target/id".into(),
                    id: id.clone(),
                }),
        }
    }
}

pub fn validate_nucleic_acid_json(input: &str, limits: &NucleicAcidLimits) -> Value {
    validation_envelope(NucleicAcidDocument::from_json_str_with_limits(
        input, limits,
    ))
}

pub fn apply_nucleic_acid_json_command(
    input: &str,
    command_input: &str,
    limits: &NucleicAcidLimits,
) -> Value {
    let result = (|| {
        let document = NucleicAcidDocument::from_json_str_with_limits(input, limits)?;
        check_limit(
            "/command",
            "command_bytes",
            command_input.len(),
            limits.max_input_bytes,
        )?;
        let command: Value = serde_json::from_str(command_input).map_err(|error| {
            NucleicAcidError::InvalidCommand {
                message: error.to_string(),
            }
        })?;
        document.apply_json_command_with_limits(&command, limits)
    })();
    validation_envelope(result)
}

fn validation_envelope(result: Result<NucleicAcidDocument, NucleicAcidError>) -> Value {
    match result {
        Ok(document) => json!({
            "schema": NUCLEIC_ACID_VALIDATION_SCHEMA,
            "ok": true,
            "document": document,
        }),
        Err(error) => json!({
            "schema": NUCLEIC_ACID_VALIDATION_SCHEMA,
            "ok": false,
            "error": error.to_json(),
        }),
    }
}

fn validate_residue_identity(
    kind: NucleicAcidStrandKind,
    residue: &NucleicAcidResidue,
    path: &str,
) -> Result<(), NucleicAcidError> {
    let expected_sugar = match kind {
        NucleicAcidStrandKind::Dna => NucleicAcidSugar::Deoxyribose,
        NucleicAcidStrandKind::Rna => NucleicAcidSugar::Ribose,
    };
    if residue.sugar != expected_sugar {
        return Err(NucleicAcidError::InvalidResidueIdentity {
            path: format!("{path}/sugar"),
            reason: format!("{:?} strand requires {:?}", kind, expected_sugar),
        });
    }
    let base_allowed = match kind {
        NucleicAcidStrandKind::Dna => !matches!(residue.base, NucleicAcidBase::U),
        NucleicAcidStrandKind::Rna => !matches!(residue.base, NucleicAcidBase::T),
    };
    if !base_allowed {
        return Err(NucleicAcidError::InvalidResidueIdentity {
            path: format!("{path}/base"),
            reason: format!("base {:?} is incompatible with {:?}", residue.base, kind),
        });
    }
    match residue.modification.as_deref() {
        None if residue.base == NucleicAcidBase::Other => {
            Err(NucleicAcidError::InvalidResidueIdentity {
                path: format!("{path}/base"),
                reason: "base 'other' requires a known modification".into(),
            })
        }
        None => Ok(()),
        Some("5mC") if residue.base == NucleicAcidBase::C => Ok(()),
        Some("m6A") if residue.base == NucleicAcidBase::A => Ok(()),
        Some("pseudouridine")
            if kind == NucleicAcidStrandKind::Rna
                && matches!(residue.base, NucleicAcidBase::U | NucleicAcidBase::Other) =>
        {
            Ok(())
        }
        Some(known @ ("5mC" | "m6A" | "pseudouridine")) => {
            Err(NucleicAcidError::InvalidResidueIdentity {
                path: format!("{path}/modification"),
                reason: format!("modification {known} is incompatible with the selected base"),
            })
        }
        Some(modification) => Err(NucleicAcidError::UnknownModification {
            path: format!("{path}/modification"),
            modification: modification.to_owned(),
        }),
    }
}

fn validate_annotations(
    document: &NucleicAcidDocument,
    limits: &NucleicAcidLimits,
) -> Result<(), NucleicAcidError> {
    let mut count = document.annotations.len();
    for strand in &document.strands {
        count += strand.annotations.len();
        for residue in &strand.residues {
            count += residue.annotations.len();
        }
    }
    for linkage in &document.linkages {
        count += linkage.annotations.len();
    }
    check_limit("/annotations", "annotations", count, limits.max_annotations)?;
    let bytes = serde_json::to_vec(&document.annotations).map_or(usize::MAX, |v| v.len())
        + document
            .strands
            .iter()
            .map(|strand| {
                serde_json::to_vec(&strand.annotations).map_or(usize::MAX, |v| v.len())
                    + strand
                        .residues
                        .iter()
                        .map(|residue| {
                            serde_json::to_vec(&residue.annotations).map_or(usize::MAX, |v| v.len())
                        })
                        .sum::<usize>()
            })
            .sum::<usize>()
        + document
            .linkages
            .iter()
            .map(|linkage| serde_json::to_vec(&linkage.annotations).map_or(usize::MAX, |v| v.len()))
            .sum::<usize>();
    check_limit(
        "/annotations",
        "annotation_bytes",
        bytes,
        limits.max_annotation_bytes,
    )
}

fn validate_annotation_key(key: &str, limits: &NucleicAcidLimits) -> Result<(), NucleicAcidError> {
    if key.is_empty() {
        return Err(NucleicAcidError::InvalidCommand {
            message: "annotation key must not be empty".into(),
        });
    }
    check_limit(
        "/command/key",
        "annotation_key_bytes",
        key.len(),
        limits.max_id_bytes,
    )
}

fn validate_id(id: &str, path: &str, limits: &NucleicAcidLimits) -> Result<(), NucleicAcidError> {
    if id.is_empty() {
        return Err(NucleicAcidError::InvalidId {
            path: path.into(),
            reason: "id must not be empty".into(),
        });
    }
    check_limit(path, "id_bytes", id.len(), limits.max_id_bytes)
}

fn check_limit(
    path: &str,
    resource: &str,
    requested: usize,
    limit: usize,
) -> Result<(), NucleicAcidError> {
    if requested > limit {
        Err(NucleicAcidError::ResourceLimit {
            path: path.into(),
            resource: resource.into(),
            requested,
            limit,
        })
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna_document() -> Value {
        json!({
            "schema": NUCLEIC_ACID_SCHEMA,
            "atom_ids": ["r1-5", "r1-3", "r2-5", "r2-3"],
            "strands": [{
                "id": "strand-1",
                "kind": "dna",
                "residues": [
                    {"id":"r1","base":"A","sugar":"deoxyribose","atom_refs":["r1-5","r1-3"],"three_prime_linkage":"l1","annotations":{}},
                    {"id":"r2","base":"C","modification":"5mC","sugar":"deoxyribose","atom_refs":["r2-5","r2-3"],"five_prime_linkage":"l1","annotations":{"label":"5-methylcytosine"}}
                ],
                "annotations": {}
            }],
            "linkages": [{
                "id":"l1","kind":"phosphodiester","three_prime_residue_id":"r1","five_prime_residue_id":"r2","three_prime_atom_ref":"r1-3","five_prime_atom_ref":"r2-5","annotations":{}
            }],
            "annotations": {"title":"bounded DNA"}
        })
    }

    #[test]
    fn valid_document_roundtrips_without_loss() {
        let input = dna_document();
        let document = NucleicAcidDocument::from_json(&input).unwrap();
        assert_eq!(document.to_json(), input);
        let reparsed = NucleicAcidDocument::from_json(&document.to_json()).unwrap();
        assert_eq!(reparsed, document);
    }

    #[test]
    fn ambiguous_atom_mapping_is_typed() {
        let mut input = dna_document();
        input["strands"][0]["residues"][1]["atom_refs"] = json!(["r1-3", "r2-3"]);
        let error = NucleicAcidDocument::from_json(&input).unwrap_err();
        assert_eq!(error.code(), "ambiguous_atom_mapping");
    }

    #[test]
    fn unknown_modification_is_typed() {
        let mut input = dna_document();
        input["strands"][0]["residues"][1]["modification"] = json!("unknown-X");
        let error = NucleicAcidDocument::from_json(&input).unwrap_err();
        assert_eq!(error.code(), "unknown_modification");
    }

    #[test]
    fn annotation_edit_preserves_topology() {
        let document = NucleicAcidDocument::from_json(&dna_document()).unwrap();
        let next = document
            .apply(&NucleicAcidCommand::SetAnnotation {
                target: NucleicAcidAnnotationTarget::Residue { id: "r1".into() },
                key: "selected".into(),
                value: json!(true),
            })
            .unwrap();
        assert_eq!(next.strands[0].residues[0].annotations["selected"], true);
        assert_eq!(next.linkages, document.linkages);
        assert_eq!(next.atom_ids, document.atom_ids);
    }

    #[test]
    fn input_limit_is_typed_before_parsing() {
        let error = NucleicAcidDocument::from_json_str_with_limits(
            "{}",
            &NucleicAcidLimits {
                max_input_bytes: 1,
                ..NucleicAcidLimits::default()
            },
        )
        .unwrap_err();
        assert_eq!(error.code(), "resource_limit");
    }
}
