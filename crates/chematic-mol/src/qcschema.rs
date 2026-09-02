//! QCSchema (MolSSI Quantum Chemistry Schema) JSON support.
//!
//! QCSchema is MolSSI's JSON schema family for exchanging quantum-chemistry
//! molecules, inputs, and results between packages (Psi4, QCEngine,
//! QCFractal, xtb-python, ...). This module implements three of its
//! document types by hand against `serde_json::Value` (see "No serde derive"
//! below):
//!
//! - [`QcMolecule`] -- the `qcschema_molecule` geometry object.
//! - [`AtomicInput`] -- the `qcschema_input` object (molecule + driver +
//!   model + keywords).
//! - [`AtomicResult`] -- the `qcschema_output` object (an `AtomicInput` plus
//!   `return_result`/`properties`/`success`/`error`/provenance).
//!
//! `OptimizationInput`/`OptimizationResult` (trajectory-bearing multi-step
//! procedures) are explicitly out of scope for this module.
//!
//! ## Spec sources (public, official)
//!
//! Field names, types, and defaults below are grounded in:
//! - Molecule field reference: <https://molssi.github.io/QCElemental/model_molecule.html>
//! - Molecule JSON Schema (v1 series): <https://github.com/MolSSI/QCSchema/blob/master/qcschema/data/v1/qc_schema_molecule.schema>
//! - `schema_name`/`schema_version` definitions: <https://github.com/MolSSI/QCSchema/blob/master/qcschema/dev/molecule.py>
//! - `AtomicResult`/`AtomicInput`/`Model`/`Provenance`/`DriverEnum` field
//!   reference: <https://molssi.github.io/QCElemental/dev/api/qcelemental.models.AtomicResult.html>
//! - `QCEngine` (context for how `driver`/`model`/`return_result` are used
//!   in practice): <https://github.com/MolSSI/QCEngine>
//!
//! No code or comments were copied from `qcelemental`/`QCEngine`; only the
//! documented field shapes were used as a reference.
//!
//! ## Schema-version ambiguity, resolved
//!
//! The v1 JSON Schema (`data/v1/qc_schema_molecule.schema`) requires only
//! `symbols` and `geometry`; the newer `dev/molecule.py` model also treats
//! `schema_name`/`schema_version` as required fields. Real-world producers
//! (QCEngine, QCArchive) always emit `schema_name`/`schema_version`, but
//! some hand-written / minimal fixtures omit them. This module targets the
//! v1 object family (`schema_version` integer `1`) and splits the
//! difference: a missing `schema_name`/`schema_version` on read defaults to
//! the canonical v1 value (`"qcschema_molecule"`/`1` etc.), but a
//! *present-and-wrong* `schema_name` is rejected. The schema's own regex
//! (`^(qc_?schema_input)$`, from the same repo) also accepts the legacy
//! underscored spelling `qc_schema_*` alongside `qcschema_*`; both are
//! accepted here and the input's exact spelling is preserved verbatim
//! through round-trip.
//!
//! ## No serde derive
//!
//! `chematic-mol` depends on `serde_json` but not `serde` (no `derive`
//! feature available in this crate -- see `Cargo.toml`, which this module
//! must not edit). All (de)serialization here is hand-written against
//! `serde_json::Value` rather than `#[derive(Serialize, Deserialize)]`.
//!
//! ## Unrecognized-field handling (consistent across this module)
//!
//! QCSchema objects are intentionally extensible (`keywords`, `extras`,
//! `properties`, `protocols`, `native_files` are open bags by design, not a
//! gap in this implementation). Two distinct buckets are kept apart on every
//! object that has them:
//! - Fields the spec itself defines as an open map (`keywords`, `extras`,
//!   `protocols`, `native_files`, `properties`) are stored as
//!   `BTreeMap<String, serde_json::Value>` (alias [`JsonObject`]) --
//!   iteration order is always the sorted key order, so output is
//!   deterministic without any extra sorting step.
//! - Any *other* top-level key not part of the documented schema for that
//!   object is preserved verbatim in a separate `unknown_fields` bag and
//!   re-emitted at the top level (not nested under `extras`) on write, so a
//!   round trip never silently drops data.
//!
//! ## Serialization rule
//!
//! A struct field typed as a plain (non-`Option`) Rust value mirrors a
//! QCSchema field that has a spec-defined default (e.g. `molecular_charge`,
//! `fix_com`, `schema_version`) -- it is *always* emitted. A field typed
//! `Option<T>` mirrors a genuinely optional QCSchema field and is omitted
//! from output when `None`. Open-bag map fields (`extras`, `keywords`, ...)
//! are omitted when empty, since an absent key and an empty object are
//! semantically identical for an open map.
//!
//! ## NaN/Infinity
//!
//! Every `parse_*` entry point below walks the *entire* parsed
//! `serde_json::Value` tree once, up front, rejecting any numeric leaf that
//! is not finite ([`QcSchemaError::NonFinite`]). Structured field extraction
//! never has to re-check finiteness. In practice, on this workspace's
//! pinned `serde_json` (1.0.151, no `arbitrary_precision` feature anywhere
//! in the dependency graph), a literal like `1e400` is already rejected by
//! `serde_json::from_str` itself ("number out of range") before it ever
//! becomes a `Value` -- `serde_json::Number` cannot represent a non-finite
//! value via any safe public constructor in that configuration. The
//! `check_finite` pass below is kept anyway as defense-in-depth (a
//! different `serde_json` feature set could change that guarantee) and
//! because it is what the task asked for explicitly; either way, a
//! document containing such a literal is guaranteed to fail closed with a
//! typed `Err`, never a panic and never a silently-carried infinity.

use std::collections::{BTreeMap, BTreeSet};

use chematic_core::{
    Atom, AtomIdx, BondOrder, Coords3D, Element, Molecule, MoleculeBuilder, Point3,
};
use serde_json::{Map, Value, json};

/// An open, order-deterministic JSON object bag (spec-extensible fields,
/// and captured unrecognized fields). See the module docs' "Unrecognized-
/// field handling" section.
pub type JsonObject = BTreeMap<String, Value>;

/// `(atom1, atom2, bond_order)` triples, 0-indexed -- [`QcMolecule::connectivity`].
pub type Connectivity = Vec<(usize, usize, f64)>;

// ─── Error ──────────────────────────────────────────────────────────────────

/// Resource limits applied before a QCSchema JSON document is decoded into
/// domain objects.
///
/// The existing parsing functions use these defaults. Call the
/// `*_with_limits` variants when an application has a tighter input budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QcSchemaParseLimits {
    /// Maximum UTF-8 input size, in bytes.
    pub max_input_bytes: usize,
    /// Maximum nesting depth of JSON arrays and objects.
    pub max_json_depth: usize,
    /// Maximum number of entries in any JSON array.
    pub max_array_items: usize,
    /// Maximum size of any JSON string, in UTF-8 bytes.
    pub max_string_bytes: usize,
}

impl Default for QcSchemaParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_json_depth: 128,
            max_array_items: 10_000_000,
            max_string_bytes: 16 * 1024 * 1024,
        }
    }
}

/// Errors that can occur while parsing a QCSchema document. Never panics --
/// every malformed/incomplete input maps to one of these variants.
#[derive(Debug, Clone, PartialEq)]
pub enum QcSchemaError {
    /// The input was not valid JSON at all.
    InvalidJson(String),
    /// A required field was absent.
    MissingField(String),
    /// A field was present but had the wrong JSON type.
    WrongType {
        field: String,
        expected: &'static str,
    },
    /// A numeric field was NaN or +/-Infinity.
    NonFinite { field: String },
    /// `schema_name` was present but did not match any accepted spelling.
    InvalidSchemaName {
        object: &'static str,
        expected: &'static [&'static str],
        found: String,
    },
    /// A parallel-array field's length did not match the atom/fragment count.
    LengthMismatch { detail: String },
    /// A connectivity/fragment atom index was out of range.
    IndexOutOfRange { detail: String },
    /// `connectivity` listed the same atom pair more than once.
    DuplicateBond { detail: String },
    /// An enum-valued field (e.g. `driver`) had an unrecognized string.
    InvalidEnumValue {
        field: &'static str,
        expected: &'static [&'static str],
        found: String,
    },
    /// A cross-field invariant was violated (e.g. `success: true` with no
    /// `return_result`).
    Inconsistent { detail: String },
    /// The document exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        path: String,
        actual: usize,
        limit: usize,
    },
}

impl std::fmt::Display for QcSchemaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidJson(s) => write!(f, "QCSchema: invalid JSON: {s}"),
            Self::MissingField(field) => write!(f, "QCSchema: missing required field '{field}'"),
            Self::WrongType { field, expected } => {
                write!(f, "QCSchema: field '{field}' expected {expected}")
            }
            Self::NonFinite { field } => {
                write!(
                    f,
                    "QCSchema: field '{field}' is NaN or Infinity, which is not allowed"
                )
            }
            Self::InvalidSchemaName {
                object,
                expected,
                found,
            } => write!(
                f,
                "QCSchema: {object}.schema_name = '{found}' does not match any of {expected:?}"
            ),
            Self::LengthMismatch { detail } => write!(f, "QCSchema: length mismatch: {detail}"),
            Self::IndexOutOfRange { detail } => write!(f, "QCSchema: index out of range: {detail}"),
            Self::DuplicateBond { detail } => write!(f, "QCSchema: duplicate bond: {detail}"),
            Self::InvalidEnumValue {
                field,
                expected,
                found,
            } => write!(
                f,
                "QCSchema: field '{field}' = '{found}' is not one of {expected:?}"
            ),
            Self::Inconsistent { detail } => write!(f, "QCSchema: inconsistent document: {detail}"),
            Self::ResourceLimit {
                resource,
                path,
                actual,
                limit,
            } => write!(
                f,
                "QCSchema: {resource} at {path} has size {actual}, exceeding limit {limit}"
            ),
        }
    }
}

impl std::error::Error for QcSchemaError {}

// ─── Generic JSON-Value helpers (hand-rolled parsing, no serde derive) ──────

fn check_finite(v: &Value, path: &str) -> Result<(), QcSchemaError> {
    match v {
        Value::Number(n) => {
            if let Some(f) = n.as_f64()
                && !f.is_finite()
            {
                return Err(QcSchemaError::NonFinite {
                    field: path.to_string(),
                });
            }
            Ok(())
        }
        Value::Array(arr) => {
            for (i, item) in arr.iter().enumerate() {
                check_finite(item, &format!("{path}[{i}]"))?;
            }
            Ok(())
        }
        Value::Object(o) => {
            for (k, val) in o {
                check_finite(val, &format!("{path}.{k}"))?;
            }
            Ok(())
        }
        _ => Ok(()),
    }
}

fn check_resource_limits(
    v: &Value,
    path: &str,
    depth: usize,
    limits: &QcSchemaParseLimits,
) -> Result<(), QcSchemaError> {
    if depth > limits.max_json_depth {
        return Err(QcSchemaError::ResourceLimit {
            resource: "JSON depth",
            path: path.to_string(),
            actual: depth,
            limit: limits.max_json_depth,
        });
    }
    match v {
        Value::String(s) if s.len() > limits.max_string_bytes => {
            Err(QcSchemaError::ResourceLimit {
                resource: "string bytes",
                path: path.to_string(),
                actual: s.len(),
                limit: limits.max_string_bytes,
            })
        }
        Value::Array(arr) => {
            if arr.len() > limits.max_array_items {
                return Err(QcSchemaError::ResourceLimit {
                    resource: "array items",
                    path: path.to_string(),
                    actual: arr.len(),
                    limit: limits.max_array_items,
                });
            }
            for (i, item) in arr.iter().enumerate() {
                check_resource_limits(item, &format!("{path}[{i}]"), depth + 1, limits)?;
            }
            Ok(())
        }
        Value::Object(o) => {
            for (k, val) in o {
                if k.len() > limits.max_string_bytes {
                    return Err(QcSchemaError::ResourceLimit {
                        resource: "object key bytes",
                        path: format!("{path}.{k}"),
                        actual: k.len(),
                        limit: limits.max_string_bytes,
                    });
                }
                check_resource_limits(val, &format!("{path}.{k}"), depth + 1, limits)?;
            }
            Ok(())
        }
        _ => Ok(()),
    }
}

fn parse_json_with_limits(
    input: &str,
    limits: &QcSchemaParseLimits,
) -> Result<Value, QcSchemaError> {
    if input.len() > limits.max_input_bytes {
        return Err(QcSchemaError::ResourceLimit {
            resource: "input bytes",
            path: "$".to_string(),
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let root: Value =
        serde_json::from_str(input).map_err(|e| QcSchemaError::InvalidJson(e.to_string()))?;
    check_resource_limits(&root, "$", 0, limits)?;
    check_finite(&root, "$")?;
    Ok(root)
}

fn obj<'a>(v: &'a Value, ctx: &str) -> Result<&'a Map<String, Value>, QcSchemaError> {
    v.as_object().ok_or_else(|| QcSchemaError::WrongType {
        field: ctx.to_string(),
        expected: "object",
    })
}

fn req<'a>(o: &'a Map<String, Value>, key: &str) -> Result<&'a Value, QcSchemaError> {
    o.get(key)
        .ok_or_else(|| QcSchemaError::MissingField(key.to_string()))
}

fn get_str(o: &Map<String, Value>, key: &str) -> Result<Option<String>, QcSchemaError> {
    match o.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.clone())),
        Some(_) => Err(QcSchemaError::WrongType {
            field: key.to_string(),
            expected: "string",
        }),
    }
}

fn req_str(o: &Map<String, Value>, key: &str) -> Result<String, QcSchemaError> {
    get_str(o, key)?.ok_or_else(|| QcSchemaError::MissingField(key.to_string()))
}

fn get_bool_field(o: &Map<String, Value>, key: &str) -> Result<Option<bool>, QcSchemaError> {
    match o.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(_) => Err(QcSchemaError::WrongType {
            field: key.to_string(),
            expected: "boolean",
        }),
    }
}

fn req_bool(o: &Map<String, Value>, key: &str) -> Result<bool, QcSchemaError> {
    get_bool_field(o, key)?.ok_or_else(|| QcSchemaError::MissingField(key.to_string()))
}

fn get_f64_checked(v: &Value, field: &str) -> Result<f64, QcSchemaError> {
    let f = v.as_f64().ok_or_else(|| QcSchemaError::WrongType {
        field: field.to_string(),
        expected: "number",
    })?;
    if !f.is_finite() {
        return Err(QcSchemaError::NonFinite {
            field: field.to_string(),
        });
    }
    Ok(f)
}

fn get_f64_field(o: &Map<String, Value>, key: &str) -> Result<Option<f64>, QcSchemaError> {
    match o.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(v) => Ok(Some(get_f64_checked(v, key)?)),
    }
}

fn item_i64(v: &Value, field: &str) -> Result<i64, QcSchemaError> {
    if let Some(n) = v.as_i64() {
        return Ok(n);
    }
    if let Some(f) = v.as_f64()
        && f.is_finite()
        && f.fract() == 0.0
    {
        return Ok(f as i64);
    }
    Err(QcSchemaError::WrongType {
        field: field.to_string(),
        expected: "integer",
    })
}

fn item_bool(v: &Value, field: &str) -> Result<bool, QcSchemaError> {
    v.as_bool().ok_or_else(|| QcSchemaError::WrongType {
        field: field.to_string(),
        expected: "boolean",
    })
}

fn item_string(v: &Value, field: &str) -> Result<String, QcSchemaError> {
    v.as_str()
        .map(str::to_string)
        .ok_or_else(|| QcSchemaError::WrongType {
            field: field.to_string(),
            expected: "string",
        })
}

fn get_i64_field(o: &Map<String, Value>, key: &str) -> Result<Option<i64>, QcSchemaError> {
    match o.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(v) => Ok(Some(item_i64(v, key)?)),
    }
}

fn item_usize(v: &Value, field: &str) -> Result<usize, QcSchemaError> {
    let n = item_i64(v, field)?;
    if n < 0 {
        return Err(QcSchemaError::WrongType {
            field: field.to_string(),
            expected: "non-negative integer",
        });
    }
    Ok(n as usize)
}

/// Parse an optional JSON array field, applying `item` to each element with
/// a path-qualified field name for error messages.
fn get_array_opt<T>(
    o: &Map<String, Value>,
    key: &str,
    item: impl Fn(&Value, &str) -> Result<T, QcSchemaError>,
) -> Result<Option<Vec<T>>, QcSchemaError> {
    match o.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Array(arr)) => arr
            .iter()
            .enumerate()
            .map(|(i, v)| item(v, &format!("{key}[{i}]")))
            .collect::<Result<Vec<T>, _>>()
            .map(Some),
        Some(_) => Err(QcSchemaError::WrongType {
            field: key.to_string(),
            expected: "array",
        }),
    }
}

fn req_array<T>(
    o: &Map<String, Value>,
    key: &str,
    item: impl Fn(&Value, &str) -> Result<T, QcSchemaError>,
) -> Result<Vec<T>, QcSchemaError> {
    get_array_opt(o, key, item)?.ok_or_else(|| QcSchemaError::MissingField(key.to_string()))
}

fn get_object_opt(o: &Map<String, Value>, key: &str) -> Result<Option<JsonObject>, QcSchemaError> {
    match o.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Object(m)) => Ok(Some(
            m.iter().map(|(k, v)| (k.clone(), v.clone())).collect(),
        )),
        Some(_) => Err(QcSchemaError::WrongType {
            field: key.to_string(),
            expected: "object",
        }),
    }
}

/// Every key of `o` not in `known`, preserved verbatim (see module docs).
fn collect_unknown(o: &Map<String, Value>, known: &[&str]) -> JsonObject {
    o.iter()
        .filter(|(k, _)| !known.contains(&k.as_str()))
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect()
}

fn json_object_value(m: &JsonObject) -> Value {
    Value::Object(m.iter().map(|(k, v)| (k.clone(), v.clone())).collect())
}

fn insert_unknown(o: &mut Map<String, Value>, unknown: &JsonObject) {
    for (k, v) in unknown {
        o.insert(k.clone(), v.clone());
    }
}

// ─── Provenance ──────────────────────────────────────────────────────────────

/// QCSchema `Provenance` object: which program/library/person created the
/// enclosing document.
#[derive(Debug, Clone, PartialEq)]
pub struct Provenance {
    pub creator: String,
    /// Defaults to `""` per spec.
    pub version: String,
    /// Defaults to `""` per spec.
    pub routine: String,
    pub unknown_fields: JsonObject,
}

const PROVENANCE_KEYS: &[&str] = &["creator", "version", "routine"];

fn parse_provenance(v: &Value) -> Result<Provenance, QcSchemaError> {
    let o = obj(v, "provenance")?;
    Ok(Provenance {
        creator: req_str(o, "creator")?,
        version: get_str(o, "version")?.unwrap_or_default(),
        routine: get_str(o, "routine")?.unwrap_or_default(),
        unknown_fields: collect_unknown(o, PROVENANCE_KEYS),
    })
}

fn provenance_to_value(p: &Provenance) -> Value {
    let mut o = Map::new();
    o.insert("creator".into(), json!(p.creator));
    o.insert("version".into(), json!(p.version));
    o.insert("routine".into(), json!(p.routine));
    insert_unknown(&mut o, &p.unknown_fields);
    Value::Object(o)
}

// ─── Molecule ────────────────────────────────────────────────────────────────

/// The QCSchema `Molecule` object (`schema_name: "qcschema_molecule"`).
///
/// Named `QcMolecule` (not `Molecule`) to avoid colliding with
/// [`chematic_core::Molecule`] when both are imported together. See
/// [`qc_molecule_to_chematic`] / [`chematic_to_qc_molecule`] for conversion.
#[derive(Debug, Clone, PartialEq)]
pub struct QcMolecule {
    /// `"qcschema_molecule"` or the legacy `"qc_schema_molecule"` spelling,
    /// preserved verbatim.
    pub schema_name: String,
    pub schema_version: i64,
    /// Element symbols, title case, length `nat`.
    pub symbols: Vec<String>,
    /// Flat Cartesian coordinates, length `3 * nat`, **Bohr (a0)** -- see
    /// [`BOHR_TO_ANGSTROM`].
    pub geometry: Vec<f64>,
    pub molecular_charge: f64,
    pub molecular_multiplicity: i64,
    pub fix_com: bool,
    pub fix_orientation: bool,
    pub masses: Option<Vec<f64>>,
    pub real: Option<Vec<bool>>,
    pub atomic_numbers: Option<Vec<i64>>,
    /// Isotope mass numbers; `-1` means unspecified (per spec).
    pub mass_numbers: Option<Vec<i64>>,
    pub atom_labels: Option<Vec<String>>,
    pub name: Option<String>,
    pub comment: Option<String>,
    /// `(atom1, atom2, bond_order)` triples, 0-indexed.
    pub connectivity: Option<Connectivity>,
    pub fragments: Option<Vec<Vec<usize>>>,
    pub fragment_charges: Option<Vec<f64>>,
    pub fragment_multiplicities: Option<Vec<i64>>,
    pub fix_symmetry: Option<String>,
    pub provenance: Option<Provenance>,
    pub id: Option<String>,
    pub extras: JsonObject,
    pub unknown_fields: JsonObject,
}

const QCSCHEMA_MOLECULE_NAMES: &[&str] = &["qcschema_molecule", "qc_schema_molecule"];

const MOLECULE_KEYS: &[&str] = &[
    "schema_name",
    "schema_version",
    "symbols",
    "geometry",
    "molecular_charge",
    "molecular_multiplicity",
    "fix_com",
    "fix_orientation",
    "masses",
    "real",
    "atomic_numbers",
    "mass_numbers",
    "atom_labels",
    "name",
    "comment",
    "connectivity",
    "fragments",
    "fragment_charges",
    "fragment_multiplicities",
    "fix_symmetry",
    "provenance",
    "id",
    "extras",
];

fn get_connectivity_opt(
    o: &Map<String, Value>,
    key: &str,
    natoms: usize,
) -> Result<Option<Connectivity>, QcSchemaError> {
    let arr = match o.get(key) {
        None | Some(Value::Null) => return Ok(None),
        Some(Value::Array(a)) => a,
        Some(_) => {
            return Err(QcSchemaError::WrongType {
                field: key.to_string(),
                expected: "array",
            });
        }
    };
    let mut out = Vec::with_capacity(arr.len());
    let mut seen: BTreeSet<(usize, usize)> = BTreeSet::new();
    for (i, item) in arr.iter().enumerate() {
        let field = format!("{key}[{i}]");
        let triple = item.as_array().ok_or_else(|| QcSchemaError::WrongType {
            field: field.clone(),
            expected: "[atom1, atom2, order]",
        })?;
        if triple.len() != 3 {
            return Err(QcSchemaError::WrongType {
                field,
                expected: "array of length 3",
            });
        }
        let a = item_usize(&triple[0], &format!("{field}[0]"))?;
        let b = item_usize(&triple[1], &format!("{field}[1]"))?;
        let order = get_f64_checked(&triple[2], &format!("{field}[2]"))?;
        if a >= natoms || b >= natoms {
            return Err(QcSchemaError::IndexOutOfRange {
                detail: format!(
                    "{field} references atom index {} but there are only {natoms} atoms",
                    a.max(b)
                ),
            });
        }
        if a == b {
            return Err(QcSchemaError::Inconsistent {
                detail: format!("{field} is a self-bond (atom {a} to itself)"),
            });
        }
        let pair = (a.min(b), a.max(b));
        if !seen.insert(pair) {
            return Err(QcSchemaError::DuplicateBond {
                detail: format!(
                    "{field} duplicates an earlier bond between atoms {} and {}",
                    pair.0, pair.1
                ),
            });
        }
        out.push((a, b, order));
    }
    Ok(Some(out))
}

fn get_fragments_opt(
    o: &Map<String, Value>,
    key: &str,
    natoms: usize,
) -> Result<Option<Vec<Vec<usize>>>, QcSchemaError> {
    let arr = match o.get(key) {
        None | Some(Value::Null) => return Ok(None),
        Some(Value::Array(a)) => a,
        Some(_) => {
            return Err(QcSchemaError::WrongType {
                field: key.to_string(),
                expected: "array",
            });
        }
    };
    let mut out = Vec::with_capacity(arr.len());
    for (i, item) in arr.iter().enumerate() {
        let field = format!("{key}[{i}]");
        let inner = item.as_array().ok_or_else(|| QcSchemaError::WrongType {
            field: field.clone(),
            expected: "array of atom indices",
        })?;
        let idxs: Vec<usize> = inner
            .iter()
            .enumerate()
            .map(|(j, v)| item_usize(v, &format!("{field}[{j}]")))
            .collect::<Result<_, _>>()?;
        for &ix in &idxs {
            if ix >= natoms {
                return Err(QcSchemaError::IndexOutOfRange {
                    detail: format!(
                        "{field} references atom index {ix} but there are only {natoms} atoms"
                    ),
                });
            }
        }
        out.push(idxs);
    }
    Ok(Some(out))
}

/// Parse a QCSchema `Molecule` JSON document.
pub fn parse_qcschema_molecule(input: &str) -> Result<QcMolecule, QcSchemaError> {
    parse_qcschema_molecule_with_limits(input, &QcSchemaParseLimits::default())
}

/// Parse a QCSchema `Molecule` JSON document with explicit resource limits.
pub fn parse_qcschema_molecule_with_limits(
    input: &str,
    limits: &QcSchemaParseLimits,
) -> Result<QcMolecule, QcSchemaError> {
    let root = parse_json_with_limits(input, limits)?;
    parse_qc_molecule_value(&root)
}

fn parse_qc_molecule_value(v: &Value) -> Result<QcMolecule, QcSchemaError> {
    let o = obj(v, "molecule")?;

    let schema_name = get_str(o, "schema_name")?.unwrap_or_else(|| "qcschema_molecule".to_string());
    if !QCSCHEMA_MOLECULE_NAMES.contains(&schema_name.as_str()) {
        return Err(QcSchemaError::InvalidSchemaName {
            object: "Molecule",
            expected: QCSCHEMA_MOLECULE_NAMES,
            found: schema_name,
        });
    }
    let schema_version = get_i64_field(o, "schema_version")?.unwrap_or(1);

    let symbols = req_array(o, "symbols", item_string)?;
    let geometry = req_array(o, "geometry", get_f64_checked)?;
    if geometry.len() != symbols.len() * 3 {
        return Err(QcSchemaError::LengthMismatch {
            detail: format!(
                "geometry has {} values, expected {} (3 * {} atoms)",
                geometry.len(),
                symbols.len() * 3,
                symbols.len()
            ),
        });
    }

    let molecular_charge = get_f64_field(o, "molecular_charge")?.unwrap_or(0.0);
    let molecular_multiplicity = get_i64_field(o, "molecular_multiplicity")?.unwrap_or(1);
    let fix_com = get_bool_field(o, "fix_com")?.unwrap_or(false);
    let fix_orientation = get_bool_field(o, "fix_orientation")?.unwrap_or(false);

    let masses = get_array_opt(o, "masses", get_f64_checked)?;
    let real = get_array_opt(o, "real", item_bool)?;
    let atomic_numbers = get_array_opt(o, "atomic_numbers", item_i64)?;
    let mass_numbers = get_array_opt(o, "mass_numbers", item_i64)?;
    let atom_labels = get_array_opt(o, "atom_labels", item_string)?;

    for (field, len) in [
        ("masses", masses.as_ref().map(Vec::len)),
        ("real", real.as_ref().map(Vec::len)),
        ("atomic_numbers", atomic_numbers.as_ref().map(Vec::len)),
        ("mass_numbers", mass_numbers.as_ref().map(Vec::len)),
        ("atom_labels", atom_labels.as_ref().map(Vec::len)),
    ] {
        if let Some(l) = len
            && l != symbols.len()
        {
            return Err(QcSchemaError::LengthMismatch {
                detail: format!(
                    "{field} has {l} entries, expected {} (one per atom)",
                    symbols.len()
                ),
            });
        }
    }

    let name = get_str(o, "name")?;
    let comment = get_str(o, "comment")?;
    let connectivity = get_connectivity_opt(o, "connectivity", symbols.len())?;
    let fragments = get_fragments_opt(o, "fragments", symbols.len())?;
    let fragment_charges = get_array_opt(o, "fragment_charges", get_f64_checked)?;
    let fragment_multiplicities = get_array_opt(o, "fragment_multiplicities", item_i64)?;
    let fix_symmetry = get_str(o, "fix_symmetry")?;
    let provenance = match o.get("provenance") {
        None | Some(Value::Null) => None,
        Some(pv) => Some(parse_provenance(pv)?),
    };
    let id = get_str(o, "id")?;
    let extras = get_object_opt(o, "extras")?.unwrap_or_default();
    let unknown_fields = collect_unknown(o, MOLECULE_KEYS);

    Ok(QcMolecule {
        schema_name,
        schema_version,
        symbols,
        geometry,
        molecular_charge,
        molecular_multiplicity,
        fix_com,
        fix_orientation,
        masses,
        real,
        atomic_numbers,
        mass_numbers,
        atom_labels,
        name,
        comment,
        connectivity,
        fragments,
        fragment_charges,
        fragment_multiplicities,
        fix_symmetry,
        provenance,
        id,
        extras,
        unknown_fields,
    })
}

fn qc_molecule_to_value(m: &QcMolecule) -> Value {
    let mut o = Map::new();
    o.insert("schema_name".into(), json!(m.schema_name));
    o.insert("schema_version".into(), json!(m.schema_version));
    o.insert("symbols".into(), json!(m.symbols));
    o.insert("geometry".into(), json!(m.geometry));
    o.insert("molecular_charge".into(), json!(m.molecular_charge));
    o.insert(
        "molecular_multiplicity".into(),
        json!(m.molecular_multiplicity),
    );
    o.insert("fix_com".into(), json!(m.fix_com));
    o.insert("fix_orientation".into(), json!(m.fix_orientation));
    if let Some(v) = &m.masses {
        o.insert("masses".into(), json!(v));
    }
    if let Some(v) = &m.real {
        o.insert("real".into(), json!(v));
    }
    if let Some(v) = &m.atomic_numbers {
        o.insert("atomic_numbers".into(), json!(v));
    }
    if let Some(v) = &m.mass_numbers {
        o.insert("mass_numbers".into(), json!(v));
    }
    if let Some(v) = &m.atom_labels {
        o.insert("atom_labels".into(), json!(v));
    }
    if let Some(v) = &m.name {
        o.insert("name".into(), json!(v));
    }
    if let Some(v) = &m.comment {
        o.insert("comment".into(), json!(v));
    }
    if let Some(c) = &m.connectivity {
        let arr: Vec<Value> = c.iter().map(|(a, b, ord)| json!([a, b, ord])).collect();
        o.insert("connectivity".into(), Value::Array(arr));
    }
    if let Some(v) = &m.fragments {
        o.insert("fragments".into(), json!(v));
    }
    if let Some(v) = &m.fragment_charges {
        o.insert("fragment_charges".into(), json!(v));
    }
    if let Some(v) = &m.fragment_multiplicities {
        o.insert("fragment_multiplicities".into(), json!(v));
    }
    if let Some(v) = &m.fix_symmetry {
        o.insert("fix_symmetry".into(), json!(v));
    }
    if let Some(p) = &m.provenance {
        o.insert("provenance".into(), provenance_to_value(p));
    }
    if let Some(v) = &m.id {
        o.insert("id".into(), json!(v));
    }
    if !m.extras.is_empty() {
        o.insert("extras".into(), json_object_value(&m.extras));
    }
    insert_unknown(&mut o, &m.unknown_fields);
    Value::Object(o)
}

/// Serialize a [`QcMolecule`] to a pretty-printed QCSchema JSON document.
pub fn write_qcschema_molecule(m: &QcMolecule) -> String {
    serde_json::to_string_pretty(&qc_molecule_to_value(m)).unwrap_or_default()
}

// ─── Model / Driver ──────────────────────────────────────────────────────────

/// QCSchema `Model.basis`: either a named basis set string (`"6-31G"`) or a
/// full inline `BasisSet` object. The `BasisSet` schema itself is not
/// independently typed here (out of scope) -- its JSON is preserved
/// verbatim.
#[derive(Debug, Clone, PartialEq)]
pub enum Basis {
    Name(String),
    Object(Value),
}

/// QCSchema `Model`: the computational method (and optional basis) to run.
#[derive(Debug, Clone, PartialEq)]
pub struct QcModel {
    pub method: String,
    pub basis: Option<Basis>,
    pub unknown_fields: JsonObject,
}

fn parse_model(v: &Value) -> Result<QcModel, QcSchemaError> {
    let o = obj(v, "model")?;
    let method = req_str(o, "method")?;
    let basis = match o.get("basis") {
        None | Some(Value::Null) => None,
        Some(Value::String(s)) => Some(Basis::Name(s.clone())),
        Some(other @ Value::Object(_)) => Some(Basis::Object(other.clone())),
        Some(_) => {
            return Err(QcSchemaError::WrongType {
                field: "model.basis".to_string(),
                expected: "string or object",
            });
        }
    };
    let unknown_fields = collect_unknown(o, &["method", "basis"]);
    Ok(QcModel {
        method,
        basis,
        unknown_fields,
    })
}

fn model_to_value(m: &QcModel) -> Value {
    let mut o = Map::new();
    o.insert("method".into(), json!(m.method));
    if let Some(b) = &m.basis {
        let bv = match b {
            Basis::Name(s) => json!(s),
            Basis::Object(v) => v.clone(),
        };
        o.insert("basis".into(), bv);
    }
    insert_unknown(&mut o, &m.unknown_fields);
    Value::Object(o)
}

/// QCSchema `DriverEnum`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Driver {
    Energy,
    Gradient,
    Hessian,
    Properties,
}

const DRIVER_VALUES: &[&str] = &["energy", "gradient", "hessian", "properties"];

impl Driver {
    pub fn as_str(self) -> &'static str {
        match self {
            Driver::Energy => "energy",
            Driver::Gradient => "gradient",
            Driver::Hessian => "hessian",
            Driver::Properties => "properties",
        }
    }

    fn parse(s: &str) -> Result<Self, QcSchemaError> {
        match s {
            "energy" => Ok(Driver::Energy),
            "gradient" => Ok(Driver::Gradient),
            "hessian" => Ok(Driver::Hessian),
            "properties" => Ok(Driver::Properties),
            _ => Err(QcSchemaError::InvalidEnumValue {
                field: "driver",
                expected: DRIVER_VALUES,
                found: s.to_string(),
            }),
        }
    }
}

// ─── AtomicInput ─────────────────────────────────────────────────────────────

/// The QCSchema `AtomicInput` object (`schema_name: "qcschema_input"`):
/// a molecule plus what to compute on it.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomicInput {
    pub schema_name: String,
    pub schema_version: i64,
    pub id: Option<String>,
    pub molecule: QcMolecule,
    pub driver: Driver,
    pub model: QcModel,
    /// Program-specific options; open bag by design.
    pub keywords: JsonObject,
    pub protocols: JsonObject,
    pub extras: JsonObject,
    pub provenance: Option<Provenance>,
    pub unknown_fields: JsonObject,
}

const QCSCHEMA_INPUT_NAMES: &[&str] = &["qcschema_input", "qc_schema_input"];

const ATOMIC_INPUT_KEYS: &[&str] = &[
    "schema_name",
    "schema_version",
    "id",
    "molecule",
    "driver",
    "model",
    "keywords",
    "protocols",
    "extras",
    "provenance",
];

/// Parse a QCSchema `AtomicInput` JSON document.
pub fn parse_atomic_input(input: &str) -> Result<AtomicInput, QcSchemaError> {
    parse_atomic_input_with_limits(input, &QcSchemaParseLimits::default())
}

/// Parse a QCSchema `AtomicInput` JSON document with explicit resource limits.
pub fn parse_atomic_input_with_limits(
    input: &str,
    limits: &QcSchemaParseLimits,
) -> Result<AtomicInput, QcSchemaError> {
    let root = parse_json_with_limits(input, limits)?;
    let o = obj(&root, "AtomicInput")?;

    let schema_name = get_str(o, "schema_name")?.unwrap_or_else(|| "qcschema_input".to_string());
    if !QCSCHEMA_INPUT_NAMES.contains(&schema_name.as_str()) {
        return Err(QcSchemaError::InvalidSchemaName {
            object: "AtomicInput",
            expected: QCSCHEMA_INPUT_NAMES,
            found: schema_name,
        });
    }
    let schema_version = get_i64_field(o, "schema_version")?.unwrap_or(1);
    let id = get_str(o, "id")?;
    let molecule = parse_qc_molecule_value(req(o, "molecule")?)?;
    let driver = Driver::parse(&req_str(o, "driver")?)?;
    let model = parse_model(req(o, "model")?)?;
    let keywords = get_object_opt(o, "keywords")?.unwrap_or_default();
    let protocols = get_object_opt(o, "protocols")?.unwrap_or_default();
    let extras = get_object_opt(o, "extras")?.unwrap_or_default();
    let provenance = match o.get("provenance") {
        None | Some(Value::Null) => None,
        Some(pv) => Some(parse_provenance(pv)?),
    };
    let unknown_fields = collect_unknown(o, ATOMIC_INPUT_KEYS);

    Ok(AtomicInput {
        schema_name,
        schema_version,
        id,
        molecule,
        driver,
        model,
        keywords,
        protocols,
        extras,
        provenance,
        unknown_fields,
    })
}

fn atomic_input_fields_to_map(a: &AtomicInput, schema_name: &str) -> Map<String, Value> {
    let mut o = Map::new();
    o.insert("schema_name".into(), json!(schema_name));
    o.insert("schema_version".into(), json!(a.schema_version));
    if let Some(v) = &a.id {
        o.insert("id".into(), json!(v));
    }
    o.insert("molecule".into(), qc_molecule_to_value(&a.molecule));
    o.insert("driver".into(), json!(a.driver.as_str()));
    o.insert("model".into(), model_to_value(&a.model));
    if !a.keywords.is_empty() {
        o.insert("keywords".into(), json_object_value(&a.keywords));
    }
    if !a.protocols.is_empty() {
        o.insert("protocols".into(), json_object_value(&a.protocols));
    }
    if !a.extras.is_empty() {
        o.insert("extras".into(), json_object_value(&a.extras));
    }
    if let Some(p) = &a.provenance {
        o.insert("provenance".into(), provenance_to_value(p));
    }
    insert_unknown(&mut o, &a.unknown_fields);
    o
}

/// Serialize an [`AtomicInput`] to a pretty-printed QCSchema JSON document.
pub fn write_atomic_input(a: &AtomicInput) -> String {
    let o = atomic_input_fields_to_map(a, &a.schema_name);
    serde_json::to_string_pretty(&Value::Object(o)).unwrap_or_default()
}

// ─── AtomicResult ────────────────────────────────────────────────────────────

/// QCSchema `ComputeError`: describes a failed computation.
#[derive(Debug, Clone, PartialEq)]
pub struct ComputeError {
    pub error_type: String,
    pub error_message: String,
    pub extras: JsonObject,
    pub unknown_fields: JsonObject,
}

const COMPUTE_ERROR_KEYS: &[&str] = &["error_type", "error_message", "extras"];

fn parse_compute_error(v: &Value) -> Result<ComputeError, QcSchemaError> {
    let o = obj(v, "error")?;
    Ok(ComputeError {
        error_type: req_str(o, "error_type")?,
        error_message: req_str(o, "error_message")?,
        extras: get_object_opt(o, "extras")?.unwrap_or_default(),
        unknown_fields: collect_unknown(o, COMPUTE_ERROR_KEYS),
    })
}

fn compute_error_to_value(e: &ComputeError) -> Value {
    let mut o = Map::new();
    o.insert("error_type".into(), json!(e.error_type));
    o.insert("error_message".into(), json!(e.error_message));
    if !e.extras.is_empty() {
        o.insert("extras".into(), json_object_value(&e.extras));
    }
    insert_unknown(&mut o, &e.unknown_fields);
    Value::Object(o)
}

/// The primary return of an `AtomicResult`, shaped by `driver`: a scalar for
/// `energy`, an array for `gradient`/`hessian`, or a property dict for
/// `properties`.
#[derive(Debug, Clone, PartialEq)]
pub enum ReturnResult {
    Scalar(f64),
    /// Always flat on this side (a nested 2-D input array, e.g. an `(nat,
    /// 3)` gradient, is flattened on read; always written back flat).
    Array(Vec<f64>),
    Properties(JsonObject),
}

fn flatten_numeric_array(v: &Value, field: &str, out: &mut Vec<f64>) -> Result<(), QcSchemaError> {
    match v {
        Value::Array(arr) => {
            for (i, item) in arr.iter().enumerate() {
                flatten_numeric_array(item, &format!("{field}[{i}]"), out)?;
            }
            Ok(())
        }
        Value::Number(_) => {
            out.push(get_f64_checked(v, field)?);
            Ok(())
        }
        _ => Err(QcSchemaError::WrongType {
            field: field.to_string(),
            expected: "number or nested array of numbers",
        }),
    }
}

fn parse_return_result(v: &Value) -> Result<ReturnResult, QcSchemaError> {
    match v {
        Value::Number(_) => Ok(ReturnResult::Scalar(get_f64_checked(v, "return_result")?)),
        Value::Array(_) => {
            let mut flat = Vec::new();
            flatten_numeric_array(v, "return_result", &mut flat)?;
            Ok(ReturnResult::Array(flat))
        }
        Value::Object(o) => Ok(ReturnResult::Properties(
            o.iter().map(|(k, v)| (k.clone(), v.clone())).collect(),
        )),
        _ => Err(QcSchemaError::WrongType {
            field: "return_result".to_string(),
            expected: "number, array, or object",
        }),
    }
}

fn return_result_to_value(r: &ReturnResult) -> Value {
    match r {
        ReturnResult::Scalar(f) => json!(f),
        ReturnResult::Array(v) => json!(v),
        ReturnResult::Properties(m) => json_object_value(m),
    }
}

/// The QCSchema `AtomicResult` object (`schema_name: "qcschema_output"`):
/// an [`AtomicInput`]'s fields plus the computed result.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomicResult {
    pub schema_name: String,
    pub schema_version: i64,
    pub id: Option<String>,
    pub molecule: QcMolecule,
    pub driver: Driver,
    pub model: QcModel,
    pub keywords: JsonObject,
    pub protocols: JsonObject,
    pub extras: JsonObject,
    /// Required (describes the program that produced this result), unlike
    /// `AtomicInput::provenance`.
    pub provenance: Provenance,
    /// Extensible property bag (energy, dipole, ...); spec-required but may
    /// be empty for a failed computation.
    pub properties: JsonObject,
    /// `Some` iff `success == true` (see [`QcSchemaError::Inconsistent`]).
    pub return_result: Option<ReturnResult>,
    pub success: bool,
    /// `Some` iff `success == false`.
    pub error: Option<ComputeError>,
    pub stdout: Option<String>,
    pub stderr: Option<String>,
    pub native_files: JsonObject,
    /// Raw `WavefunctionProperties` JSON, not independently typed (out of
    /// scope); preserved verbatim.
    pub wavefunction: Option<Value>,
    pub unknown_fields: JsonObject,
}

const QCSCHEMA_OUTPUT_NAMES: &[&str] = &["qcschema_output", "qc_schema_output"];

const ATOMIC_RESULT_EXTRA_KEYS: &[&str] = &[
    "provenance",
    "properties",
    "return_result",
    "success",
    "error",
    "stdout",
    "stderr",
    "native_files",
    "wavefunction",
];

/// Parse a QCSchema `AtomicResult` JSON document.
pub fn parse_atomic_result(input: &str) -> Result<AtomicResult, QcSchemaError> {
    parse_atomic_result_with_limits(input, &QcSchemaParseLimits::default())
}

/// Parse a QCSchema `AtomicResult` JSON document with explicit resource limits.
pub fn parse_atomic_result_with_limits(
    input: &str,
    limits: &QcSchemaParseLimits,
) -> Result<AtomicResult, QcSchemaError> {
    let root = parse_json_with_limits(input, limits)?;
    let o = obj(&root, "AtomicResult")?;

    let schema_name = get_str(o, "schema_name")?.unwrap_or_else(|| "qcschema_output".to_string());
    if !QCSCHEMA_OUTPUT_NAMES.contains(&schema_name.as_str()) {
        return Err(QcSchemaError::InvalidSchemaName {
            object: "AtomicResult",
            expected: QCSCHEMA_OUTPUT_NAMES,
            found: schema_name,
        });
    }
    let schema_version = get_i64_field(o, "schema_version")?.unwrap_or(1);
    let id = get_str(o, "id")?;
    let molecule = parse_qc_molecule_value(req(o, "molecule")?)?;
    let driver = Driver::parse(&req_str(o, "driver")?)?;
    let model = parse_model(req(o, "model")?)?;
    let keywords = get_object_opt(o, "keywords")?.unwrap_or_default();
    let protocols = get_object_opt(o, "protocols")?.unwrap_or_default();
    let extras = get_object_opt(o, "extras")?.unwrap_or_default();
    let provenance = parse_provenance(req(o, "provenance")?)?;
    let properties = get_object_opt(o, "properties")?.unwrap_or_default();
    let return_result = match o.get("return_result") {
        None | Some(Value::Null) => None,
        Some(rv) => Some(parse_return_result(rv)?),
    };
    let success = req_bool(o, "success")?;
    let error = match o.get("error") {
        None | Some(Value::Null) => None,
        Some(ev) => Some(parse_compute_error(ev)?),
    };
    let stdout = get_str(o, "stdout")?;
    let stderr = get_str(o, "stderr")?;
    let native_files = get_object_opt(o, "native_files")?.unwrap_or_default();
    let wavefunction = o.get("wavefunction").filter(|v| !v.is_null()).cloned();

    if success && return_result.is_none() {
        return Err(QcSchemaError::Inconsistent {
            detail: "success = true but return_result is missing".to_string(),
        });
    }
    if !success && error.is_none() {
        return Err(QcSchemaError::Inconsistent {
            detail: "success = false but error is missing".to_string(),
        });
    }

    let mut known: Vec<&str> = ATOMIC_INPUT_KEYS.to_vec();
    known.extend_from_slice(ATOMIC_RESULT_EXTRA_KEYS);
    let unknown_fields = collect_unknown(o, &known);

    Ok(AtomicResult {
        schema_name,
        schema_version,
        id,
        molecule,
        driver,
        model,
        keywords,
        protocols,
        extras,
        provenance,
        properties,
        return_result,
        success,
        error,
        stdout,
        stderr,
        native_files,
        wavefunction,
        unknown_fields,
    })
}

/// Serialize an [`AtomicResult`] to a pretty-printed QCSchema JSON document.
pub fn write_atomic_result(r: &AtomicResult) -> String {
    // Reuse the AtomicInput-shaped fields via a throwaway AtomicInput-like
    // map builder to avoid duplicating molecule/driver/model/keywords
    // serialization logic.
    let shared = AtomicInput {
        schema_name: r.schema_name.clone(),
        schema_version: r.schema_version,
        id: r.id.clone(),
        molecule: r.molecule.clone(),
        driver: r.driver,
        model: r.model.clone(),
        keywords: r.keywords.clone(),
        protocols: r.protocols.clone(),
        extras: r.extras.clone(),
        provenance: None, // AtomicResult's provenance is written separately below (required, not Option)
        unknown_fields: JsonObject::new(),
    };
    let mut o = atomic_input_fields_to_map(&shared, &r.schema_name);
    o.insert("provenance".into(), provenance_to_value(&r.provenance));
    o.insert("properties".into(), json_object_value(&r.properties));
    if let Some(rr) = &r.return_result {
        o.insert("return_result".into(), return_result_to_value(rr));
    }
    o.insert("success".into(), json!(r.success));
    if let Some(e) = &r.error {
        o.insert("error".into(), compute_error_to_value(e));
    }
    if let Some(v) = &r.stdout {
        o.insert("stdout".into(), json!(v));
    }
    if let Some(v) = &r.stderr {
        o.insert("stderr".into(), json!(v));
    }
    if !r.native_files.is_empty() {
        o.insert("native_files".into(), json_object_value(&r.native_files));
    }
    if let Some(v) = &r.wavefunction {
        o.insert("wavefunction".into(), v.clone());
    }
    insert_unknown(&mut o, &r.unknown_fields);
    serde_json::to_string_pretty(&Value::Object(o)).unwrap_or_default()
}

// ─── Conversion to/from chematic_core::Molecule ─────────────────────────────

/// CODATA 2018 Bohr radius in Angstrom (1 a0 = 0.529177210903 Å), matching
/// the constant set `qcelemental` uses as of this writing. QCSchema
/// `geometry` is always in Bohr; `chematic_core::Coords3D` is always in
/// Angstrom (see `chematic_core::coords3d` module docs) -- this constant is
/// the single conversion point between the two.
pub const BOHR_TO_ANGSTROM: f64 = 0.529177210903;

/// Errors converting between [`QcMolecule`] and `chematic_core`'s own
/// molecule/coordinate types.
#[derive(Debug, Clone, PartialEq)]
pub enum QcConvertError {
    UnknownElement(String),
    InvalidBondIndex {
        index: usize,
        atom_count: usize,
    },
    /// `coords.atom_count() != mol.atom_count()` on the chematic -> QCSchema
    /// direction.
    AtomCountMismatch {
        molecule: usize,
        coords: usize,
    },
}

impl std::fmt::Display for QcConvertError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnknownElement(s) => {
                write!(f, "QCSchema conversion: unknown element symbol '{s}'")
            }
            Self::InvalidBondIndex { index, atom_count } => write!(
                f,
                "QCSchema conversion: connectivity references atom index {index} but there are only {atom_count} atoms"
            ),
            Self::AtomCountMismatch { molecule, coords } => write!(
                f,
                "QCSchema conversion: molecule has {molecule} atoms but coords has {coords}"
            ),
        }
    }
}

impl std::error::Error for QcConvertError {}

/// Bundle returned by [`qc_molecule_to_chematic`]: `chematic_core::Molecule`
/// has no molecule-level charge/multiplicity field and no built-in
/// coordinate storage of its own (coordinates live in the separate
/// [`Coords3D`] type, same convention as this crate's XYZ reader), so both
/// are threaded through here rather than silently dropped.
// `chematic_core::Molecule` implements neither `Debug` nor `PartialEq`, so
// this bundle can't derive them either.
#[derive(Clone)]
pub struct ChematicMoleculeView {
    pub molecule: Molecule,
    /// Angstrom (converted from the source's Bohr geometry).
    pub coords: Coords3D,
    pub molecular_charge: f64,
    pub molecular_multiplicity: i64,
}

fn qc_bond_order_to_chematic(v: f64) -> BondOrder {
    if (v - 1.5).abs() < 1e-6 {
        return BondOrder::Aromatic;
    }
    match v.round() as i64 {
        0 => BondOrder::Zero,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Quadruple,
        _ => BondOrder::Single,
    }
}

fn chematic_bond_order_to_qc(o: BondOrder) -> f64 {
    match o {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Aromatic => 1.5,
        BondOrder::Zero => 0.0,
        BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => 1.0,
    }
}

/// Convert a QCSchema [`QcMolecule`] into chematic's own graph + coordinate
/// representation.
///
/// **Gained:** nothing beyond what QCSchema already carries.
/// **Lost:** `fix_com`/`fix_orientation`/`fragments`/`fragment_charges`/
/// `fragment_multiplicities`/`real` (ghost-atom flag)/`masses` (explicit
/// isotope masses)/`atom_labels`/`name`/`comment`/`fix_symmetry`/
/// `provenance` -- `chematic_core::Molecule` and [`Coords3D`] have no
/// equivalent slots for any of these; call sites that need them should keep
/// the original [`QcMolecule`] around. `connectivity` bond orders (plain
/// floats, e.g. `1.0`, `2.0`, `1.5`) are mapped onto the nearest
/// [`BondOrder`]; QCSchema has no stereo/aromaticity-flag concept distinct
/// from that numeric order. `mass_numbers` (when present and not `-1`) is
/// mapped onto `Atom::isotope`.
pub fn qc_molecule_to_chematic(qc: &QcMolecule) -> Result<ChematicMoleculeView, QcConvertError> {
    let mut builder = MoleculeBuilder::new();
    let mut coords = Coords3D::new_zeroed(qc.symbols.len());

    for (i, sym) in qc.symbols.iter().enumerate() {
        let element =
            Element::from_symbol(sym).ok_or_else(|| QcConvertError::UnknownElement(sym.clone()))?;
        let mut atom = Atom::new(element);
        if let Some(mn) = qc.mass_numbers.as_ref().and_then(|v| v.get(i))
            && *mn >= 0
        {
            atom.isotope = Some(*mn as u16);
        }
        let idx = builder.add_atom(atom);
        let p = Point3::new(
            qc.geometry[i * 3] * BOHR_TO_ANGSTROM,
            qc.geometry[i * 3 + 1] * BOHR_TO_ANGSTROM,
            qc.geometry[i * 3 + 2] * BOHR_TO_ANGSTROM,
        );
        coords.set(idx, p);
    }

    if let Some(conn) = &qc.connectivity {
        for (a, b, order) in conn {
            if *a >= qc.symbols.len() || *b >= qc.symbols.len() {
                return Err(QcConvertError::InvalidBondIndex {
                    index: (*a).max(*b),
                    atom_count: qc.symbols.len(),
                });
            }
            let bo = qc_bond_order_to_chematic(*order);
            // Length/duplicate/self-bond invariants were already enforced by
            // `parse_qc_molecule_value`; `add_bond` cannot fail here.
            let _ = builder.add_bond(AtomIdx(*a as u32), AtomIdx(*b as u32), bo);
        }
    }

    Ok(ChematicMoleculeView {
        molecule: builder.build(),
        coords,
        molecular_charge: qc.molecular_charge,
        molecular_multiplicity: qc.molecular_multiplicity,
    })
}

/// Convert a chematic `Molecule` + its `Coords3D` (Angstrom) + molecule-level
/// charge/multiplicity into a QCSchema [`QcMolecule`] (Bohr).
///
/// **Gained:** nothing. **Lost:** per-atom formal charge, isotope-vs-mass-
/// number distinction (only the isotope's mass number survives, as
/// `mass_numbers`; `masses` -- fractional atomic weight -- is left `None`,
/// letting a QCSchema consumer fall back to canonical weights), aromaticity,
/// stereochemistry (`Chirality`/CIP code), and wildcard/atom-map annotations
/// -- QCSchema's `Molecule` is a geometry+connectivity schema, not a full
/// cheminformatics model, and has no field for any of these.
pub fn chematic_to_qc_molecule(
    mol: &Molecule,
    coords: &Coords3D,
    molecular_charge: f64,
    molecular_multiplicity: i64,
) -> Result<QcMolecule, QcConvertError> {
    if coords.atom_count() != mol.atom_count() {
        return Err(QcConvertError::AtomCountMismatch {
            molecule: mol.atom_count(),
            coords: coords.atom_count(),
        });
    }

    let mut symbols = Vec::with_capacity(mol.atom_count());
    let mut geometry = Vec::with_capacity(mol.atom_count() * 3);
    let mut mass_numbers = Vec::with_capacity(mol.atom_count());
    let mut any_isotope = false;

    for (idx, atom) in mol.atoms() {
        symbols.push(atom.element.symbol().to_string());
        let p = coords.get(idx);
        geometry.push(p.x / BOHR_TO_ANGSTROM);
        geometry.push(p.y / BOHR_TO_ANGSTROM);
        geometry.push(p.z / BOHR_TO_ANGSTROM);
        match atom.isotope {
            Some(iso) => {
                mass_numbers.push(iso as i64);
                any_isotope = true;
            }
            None => mass_numbers.push(-1),
        }
    }

    let connectivity = if mol.bond_count() > 0 {
        Some(
            mol.bonds()
                .map(|(_, b)| {
                    (
                        b.atom1.0 as usize,
                        b.atom2.0 as usize,
                        chematic_bond_order_to_qc(b.order),
                    )
                })
                .collect(),
        )
    } else {
        None
    };

    Ok(QcMolecule {
        schema_name: "qcschema_molecule".to_string(),
        schema_version: 1,
        symbols,
        geometry,
        molecular_charge,
        molecular_multiplicity,
        fix_com: false,
        fix_orientation: false,
        masses: None,
        real: None,
        atomic_numbers: None,
        mass_numbers: if any_isotope {
            Some(mass_numbers)
        } else {
            None
        },
        atom_labels: None,
        name: None,
        comment: None,
        connectivity,
        fragments: None,
        fragment_charges: None,
        fragment_multiplicities: None,
        fix_symmetry: None,
        provenance: None,
        id: None,
        extras: JsonObject::new(),
        unknown_fields: JsonObject::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn minimal_molecule() -> &'static str {
        r#"{"symbols":["H"],"geometry":[0.0,0.0,0.0]}"#
    }

    #[test]
    fn bounded_molecule_parser_rejects_oversized_input() {
        let limits = QcSchemaParseLimits {
            max_input_bytes: 8,
            ..Default::default()
        };
        let err = parse_qcschema_molecule_with_limits(minimal_molecule(), &limits).unwrap_err();
        assert!(matches!(
            err,
            QcSchemaError::ResourceLimit {
                resource: "input bytes",
                ..
            }
        ));
    }

    #[test]
    fn bounded_molecule_parser_rejects_oversized_arrays() {
        let limits = QcSchemaParseLimits {
            max_array_items: 2,
            ..Default::default()
        };
        let err = parse_qcschema_molecule_with_limits(minimal_molecule(), &limits).unwrap_err();
        assert!(matches!(
            err,
            QcSchemaError::ResourceLimit {
                resource: "array items",
                ..
            }
        ));
    }

    #[test]
    fn bounded_atomic_parsers_apply_json_limits_before_schema_validation() {
        let limits = QcSchemaParseLimits {
            max_string_bytes: 8,
            ..Default::default()
        };
        let input = r#"{"x":"longvalue","symbols":["H"],"geometry":[0.0,0.0,0.0]}"#;
        assert!(matches!(
            parse_atomic_input_with_limits(input, &limits),
            Err(QcSchemaError::ResourceLimit {
                resource: "string bytes",
                ..
            })
        ));
        assert!(matches!(
            parse_atomic_result_with_limits(input, &limits),
            Err(QcSchemaError::ResourceLimit {
                resource: "string bytes",
                ..
            })
        ));
    }
}
