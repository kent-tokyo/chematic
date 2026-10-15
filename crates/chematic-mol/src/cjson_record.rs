//! Immutable ChemicalJSON records: strict graph validation and opaque retention.
use crate::cjson::{CjsonError, CjsonParseLimits, check_json_limits, parse_cjson_value};
use chematic_core::Molecule;
use serde::de::{MapAccess, SeqAccess, Visitor};
use serde::{Deserialize, Deserializer};
use serde_json::{Map, Value};
use std::collections::{BTreeMap, HashSet};

/// A validated ChemicalJSON source with a bounded graph projection.
///
/// The local profile accepts version 1 (or an absent version), atomic numbers,
/// formal charges, isotopes, 3D coordinates, bond indices/orders, properties.name
/// and root partialCharges datasets. Other fields are retained opaquely,
/// including stereo extensions. No edits or interpretation of opaque references are
/// promised. [`write_cjson_record`] returns the exact original UTF-8 source.
#[derive(Clone)]
pub struct CjsonRecord {
    source: String,
    molecule: Molecule,
    coordinates: Vec<(f64, f64, f64)>,
    opaque_fields: Vec<String>,
    declared_version: Option<u64>,
    coordinates_present: bool,
    formal_charges_explicit: bool,
    isotopes_explicit: bool,
    bond_orders_explicit: bool,
    name: Option<String>,
    partial_charge_datasets: BTreeMap<String, Vec<f64>>,
}

impl CjsonRecord {
    pub fn source(&self) -> &str {
        &self.source
    }
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
    pub fn coordinates(&self) -> &[(f64, f64, f64)] {
        &self.coordinates
    }
    /// JSON Pointer paths for whole fields not interpreted by the projection.
    /// Paths are sorted; an opaque subtree is represented by its root path.
    pub fn opaque_fields(&self) -> &[String] {
        &self.opaque_fields
    }
    pub fn declared_version(&self) -> Option<u64> {
        self.declared_version
    }
    pub fn coordinates_present(&self) -> bool {
        self.coordinates_present
    }
    pub fn formal_charges_explicit(&self) -> bool {
        self.formal_charges_explicit
    }
    pub fn isotopes_explicit(&self) -> bool {
        self.isotopes_explicit
    }
    pub fn bond_orders_explicit(&self) -> bool {
        self.bond_orders_explicit
    }
    /// Source properties.name, without trimming or chemical interpretation.
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }
    /// Complete finite per-atom datasets, keyed by the source method label.
    /// Labels are declarations only; no charge calculation is performed.
    pub fn partial_charge_datasets(&self) -> &BTreeMap<String, Vec<f64>> {
        &self.partial_charge_datasets
    }

    pub fn to_json(&self) -> Value {
        serde_json::json!({
            "schema": "chematic.cjson-record.v1",
            "source": self.source,
            "declared_version": self.declared_version,
            "coordinates": self.coordinates,
            "coordinates_present": self.coordinates_present,
            "formal_charges": self.molecule.atoms().map(|(_, a)| a.charge).collect::<Vec<_>>(),
            "isotopes": self.molecule.atoms().map(|(_, a)| a.isotope).collect::<Vec<_>>(),
            "opaque_fields": self.opaque_fields,
            "name": self.name,
            "partial_charge_datasets": self.partial_charge_datasets,
            "graph_projection": {
                "formal_charges_explicit": self.formal_charges_explicit,
                "isotopes_explicit": self.isotopes_explicit,
                "bond_orders_explicit": self.bond_orders_explicit,
                "aromatic_bond_4_policy": "local_alias_for_1.5",
                "atom_aromaticity": "not_perceived",
                "coordinate_stereo": "not_interpreted",
                "opaque_fields": "not_interpreted"
            }
        })
    }
}

pub fn parse_cjson_record(input: &str) -> Result<CjsonRecord, CjsonError> {
    parse_cjson_record_with_limits(input, &CjsonParseLimits::default())
}

pub fn parse_cjson_record_with_limits(
    input: &str,
    limits: &CjsonParseLimits,
) -> Result<CjsonRecord, CjsonError> {
    if input.len() > limits.max_input_bytes {
        return Err(CjsonError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    // serde_json's default tree reader replaces duplicate keys. Require
    // uniqueness recursively, including decoded/escaped spellings of keys.
    let value: UniqueValue =
        serde_json::from_str(input).map_err(|error| CjsonError::InvalidJson(error.to_string()))?;
    let value = value.0;
    check_json_limits(&value, 0, limits)?;
    validate_projection(&value)?;
    let (molecule, coordinates) = parse_cjson_value(&value, limits)?;
    let mut opaque_fields = Vec::new();
    collect_opaque(&value, "", &mut opaque_fields);
    opaque_fields.sort();
    Ok(CjsonRecord {
        source: input.to_owned(),
        molecule,
        coordinates,
        opaque_fields,
        declared_version: value.get("chemicalJson").and_then(Value::as_u64),
        coordinates_present: value.pointer("/atoms/coords/3d").is_some(),
        formal_charges_explicit: value.pointer("/atoms/formalCharges").is_some(),
        isotopes_explicit: value.pointer("/atoms/isotopes").is_some(),
        bond_orders_explicit: value.pointer("/bonds/order").is_some(),
        name: value
            .pointer("/properties/name")
            .and_then(Value::as_str)
            .map(str::to_owned),
        partial_charge_datasets: value
            .get("partialCharges")
            .and_then(Value::as_object)
            .map(|sets| {
                sets.iter()
                    .map(|(method, values)| {
                        (
                            method.clone(),
                            values
                                .as_array()
                                .unwrap()
                                .iter()
                                .map(|v| v.as_f64().unwrap())
                                .collect(),
                        )
                    })
                    .collect()
            })
            .unwrap_or_default(),
    })
}

pub fn write_cjson_record(record: &CjsonRecord) -> String {
    record.source.clone()
}

fn invalid(field: &'static str) -> CjsonError {
    CjsonError::InvalidValue { field }
}
fn object<'a>(value: &'a Value, field: &'static str) -> Result<&'a Map<String, Value>, CjsonError> {
    value.as_object().ok_or_else(|| invalid(field))
}
fn array<'a>(value: &'a Value, field: &'static str) -> Result<&'a Vec<Value>, CjsonError> {
    value.as_array().ok_or_else(|| invalid(field))
}
fn validate_projection(value: &Value) -> Result<(), CjsonError> {
    let root = object(value, "root")?;
    if let Some(version) = root.get("chemicalJson")
        && version.as_u64() != Some(1)
    {
        return Err(invalid("chemicalJson"));
    }
    let atoms = object(
        root.get("atoms").ok_or(CjsonError::MissingField("atoms"))?,
        "atoms",
    )?;
    let elements = object(
        atoms
            .get("elements")
            .ok_or(CjsonError::MissingField("atoms.elements"))?,
        "atoms.elements",
    )?;
    let numbers = array(
        elements
            .get("number")
            .ok_or(CjsonError::MissingField("atoms.elements.number"))?,
        "atoms.elements.number",
    )?;
    let n = numbers.len();
    if let Some(properties) = root.get("properties") {
        let properties = object(properties, "properties")?;
        if let Some(name) = properties.get("name")
            && !name.is_string()
        {
            return Err(invalid("properties.name"));
        }
    }
    if let Some(datasets) = root.get("partialCharges") {
        for values in object(datasets, "partialCharges")?.values() {
            let values = array(values, "partialCharges dataset")?;
            if values.len() != n
                || values
                    .iter()
                    .any(|v| !v.as_f64().is_some_and(f64::is_finite))
            {
                return Err(invalid("partialCharges dataset"));
            }
        }
    }
    for (key, field) in [
        ("formalCharges", "atoms.formalCharges"),
        ("isotopes", "atoms.isotopes"),
    ] {
        if let Some(value) = atoms.get(key)
            && array(value, field)?.len() != n
        {
            return Err(invalid(field));
        }
    }
    if let Some(coords) = atoms.get("coords") {
        let coords = object(coords, "atoms.coords")?;
        if let Some(value) = coords.get("3d")
            && array(value, "atoms.coords.3d")?.len() != n.saturating_mul(3)
        {
            return Err(invalid("atoms.coords.3d"));
        }
    }
    if let Some(bonds) = root.get("bonds") {
        let bonds = object(bonds, "bonds")?;
        let connections = object(
            bonds
                .get("connections")
                .ok_or(CjsonError::MissingField("bonds.connections"))?,
            "bonds.connections",
        )?;
        let indices = array(
            connections
                .get("index")
                .ok_or(CjsonError::MissingField("bonds.connections.index"))?,
            "bonds.connections.index",
        )?;
        if !indices.len().is_multiple_of(2) {
            return Err(CjsonError::OddBondIndexList);
        }
        if let Some(orders) = bonds.get("order") {
            let orders = array(orders, "bonds.order")?;
            if orders.len() != indices.len() / 2 {
                return Err(invalid("bonds.order"));
            }
            for value in orders {
                if !matches!(value.as_f64(), Some(1.0 | 1.5 | 2.0 | 3.0 | 4.0)) {
                    return Err(invalid("bonds.order"));
                }
            }
        }
        let mut edges = HashSet::new();
        for pair in indices.as_chunks::<2>().0 {
            let a = pair[0]
                .as_u64()
                .ok_or_else(|| invalid("bonds.connections.index"))?;
            let b = pair[1]
                .as_u64()
                .ok_or_else(|| invalid("bonds.connections.index"))?;
            if a == b || !edges.insert((a.min(b), a.max(b))) {
                return Err(invalid("bonds.connections.index"));
            }
        }
    }
    Ok(())
}

fn collect_opaque(value: &Value, path: &str, out: &mut Vec<String>) {
    let known: &[&str] = match path {
        "" => &[
            "chemicalJson",
            "atoms",
            "bonds",
            "properties",
            "partialCharges",
        ],
        "/properties" => &["name"],
        "/atoms" => &["elements", "formalCharges", "isotopes", "coords"],
        "/atoms/elements" => &["number"],
        "/atoms/coords" => &["3d"],
        "/bonds" => &["connections", "order"],
        "/bonds/connections" => &["index"],
        _ => return,
    };
    if let Some(fields) = value.as_object() {
        for (key, child) in fields {
            let pointer = format!("{path}/{}", key.replace('~', "~0").replace('/', "~1"));
            if known.contains(&key.as_str()) {
                collect_opaque(child, &pointer, out);
            } else {
                out.push(pointer);
            }
        }
    }
}

pub(crate) struct UniqueValue(pub(crate) Value);
impl<'de> Deserialize<'de> for UniqueValue {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(UniqueVisitor)
    }
}
struct UniqueVisitor;
impl<'de> Visitor<'de> for UniqueVisitor {
    type Value = UniqueValue;
    fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str("JSON with unique object keys")
    }
    fn visit_bool<E: serde::de::Error>(self, value: bool) -> Result<Self::Value, E> {
        Ok(UniqueValue(Value::Bool(value)))
    }
    fn visit_i64<E: serde::de::Error>(self, value: i64) -> Result<Self::Value, E> {
        Ok(UniqueValue(Value::Number(value.into())))
    }
    fn visit_u64<E: serde::de::Error>(self, value: u64) -> Result<Self::Value, E> {
        Ok(UniqueValue(Value::Number(value.into())))
    }
    fn visit_f64<E: serde::de::Error>(self, value: f64) -> Result<Self::Value, E> {
        serde_json::Number::from_f64(value)
            .map(|v| UniqueValue(Value::Number(v)))
            .ok_or_else(|| E::custom("nonfinite JSON number"))
    }
    fn visit_str<E: serde::de::Error>(self, value: &str) -> Result<Self::Value, E> {
        Ok(UniqueValue(Value::String(value.to_owned())))
    }
    fn visit_string<E: serde::de::Error>(self, value: String) -> Result<Self::Value, E> {
        Ok(UniqueValue(Value::String(value)))
    }
    fn visit_unit<E: serde::de::Error>(self) -> Result<Self::Value, E> {
        Ok(UniqueValue(Value::Null))
    }
    fn visit_none<E: serde::de::Error>(self) -> Result<Self::Value, E> {
        Ok(UniqueValue(Value::Null))
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
        let mut values = Vec::new();
        while let Some(value) = seq.next_element::<UniqueValue>()? {
            values.push(value.0);
        }
        Ok(UniqueValue(Value::Array(values)))
    }
    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        let mut fields = Map::new();
        while let Some(key) = map.next_key::<String>()? {
            if fields.contains_key(&key) {
                return Err(serde::de::Error::custom("duplicate object key"));
            }
            let value = map.next_value::<UniqueValue>()?;
            fields.insert(key, value.0);
        }
        Ok(UniqueValue(Value::Object(fields)))
    }
}
