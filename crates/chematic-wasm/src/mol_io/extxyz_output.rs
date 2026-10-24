//! Strict output-option decoding; never ignore malformed or duplicate fields.
use super::*;
use serde::de::{MapAccess, Visitor};
use serde::{Deserialize, Deserializer};
use std::collections::BTreeMap;
use std::marker::PhantomData;

struct UniqueMap<T>(BTreeMap<String, T>);
impl<'de, T: Deserialize<'de>> Deserialize<'de> for UniqueMap<T> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct MapVisitor<T>(PhantomData<T>);
        impl<'de, T: Deserialize<'de>> Visitor<'de> for MapVisitor<T> {
            type Value = UniqueMap<T>;
            fn expecting(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.write_str("an object with unique keys")
            }
            fn visit_map<M: MapAccess<'de>>(self, mut map: M) -> Result<Self::Value, M::Error> {
                let mut result = BTreeMap::new();
                while let Some((key, value)) = map.next_entry::<String, T>()? {
                    if result.insert(key.clone(), value).is_some() {
                        return Err(serde::de::Error::custom(format!("duplicate key {key:?}")));
                    }
                }
                Ok(UniqueMap(result))
            }
        }
        deserializer.deserialize_map(MapVisitor(PhantomData))
    }
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ExtxyzOptions {
    lattice: Option<[f64; 9]>,
    properties: Option<UniqueMap<Vec<Vec<f64>>>>,
    info: Option<UniqueMap<String>>,
}

/// Build the [`chematic_mol::XyzFrame`] for [`to_extxyz_json`] from its raw
/// string arguments, with a plain `String` error -- kept separate from the
/// `#[wasm_bindgen]` wrapper (which maps `Err` to `JsValue::from_str` at the
/// boundary) so tests can exercise the error paths directly, matching
/// `parse_pdb_molecule_and_coords`'s split: constructing a `JsValue` outside
/// a real wasm/JS runtime aborts the native test process ("JsValue
/// native-abort" -- see other tests in `tests.rs` for the same pattern).
pub(crate) fn extxyz_frame_from_json_args(
    mol: &chematic_core::Molecule,
    coords_json: &str,
    options_json: &str,
) -> Result<chematic_mol::XyzFrame, String> {
    if coords_json.len() > WASM_MAX_JSON_STRING_BYTES
        || options_json.len() > WASM_MAX_JSON_STRING_BYTES
    {
        return Err("input JSON too large".to_string());
    }
    let coords: Vec<[f64; 3]> =
        serde_json::from_str(coords_json).map_err(|e| format!("invalid coords_json: {e}"))?;
    let n = mol.atom_count();
    if coords.len() != n {
        return Err(format!(
            "coords_json has {} row(s), mol has {n} atom(s)",
            coords.len()
        ));
    }
    let atoms: Vec<chematic_mol::XyzAtom> = (0..n)
        .map(|i| {
            let element = mol.atom(chematic_core::AtomIdx(i as u32)).element;
            let [x, y, z] = coords[i];
            chematic_mol::XyzAtom { element, x, y, z }
        })
        .collect();

    if !options_json.trim_start().starts_with('{') {
        return Err("options_json must be an object".into());
    }
    let options: ExtxyzOptions =
        serde_json::from_str(options_json).map_err(|e| format!("invalid options_json: {e}"))?;
    let lattice = options.lattice;

    let mut properties = Vec::new();
    if let Some(UniqueMap(props)) = options.properties {
        for (name, rows) in props {
            if rows.len() != n {
                return Err(format!(
                    "properties['{name}'] has {} row(s), mol has {n} atom(s)",
                    rows.len()
                ));
            }
            let count = rows.first().map_or(0, |r| r.len());
            if rows.iter().any(|r| r.len() != count) {
                return Err(format!(
                    "properties['{name}'] rows have inconsistent lengths"
                ));
            }
            let values = rows
                .into_iter()
                .map(|r| r.into_iter().map(chematic_mol::XyzValue::Real).collect())
                .collect();
            properties.push(chematic_mol::XyzProperty {
                name,
                kind: chematic_mol::XyzPropertyKind::Real,
                count,
                values,
            });
        }
    }

    let mut info = Vec::new();
    if let Some(UniqueMap(info_obj)) = options.info {
        info.extend(info_obj);
    }

    Ok(chematic_mol::XyzFrame {
        atoms,
        comment: String::new(),
        lattice,
        properties,
        info,
    })
}

/// Write a molecule + coordinates as an Extended XYZ (extxyz) frame.
///
/// `coords_json`: `[[x,y,z],...]` (Å), same order and length as `mol`'s
/// atoms.
///
/// `options_json`: an optional JSON object,
/// `{"lattice":[9 numbers]|null,"properties":{"name":[[...]],...},"info":{"key":"value",...}}`
/// -- pass `"{}"` for a plain (non-extended) frame. Only real-valued
/// per-atom `properties` columns are supported from this binding (matches
/// the Python `to_extxyz` binding's scope); build a
/// `chematic_mol::XyzProperty` directly from Rust for integer/string/logical
/// columns.
///
/// Returns a JS error if `coords_json`/`options_json` are malformed, if
/// `coords_json`'s length doesn't match `mol`'s atom count, or if a
/// `properties` column's row count doesn't match it. Options must be an object
/// with only lattice/properties/info; duplicate keys and malformed values fail.
/// Nonfinite values, ambiguous/reserved names and invalid frame data fail before
/// output. Bare empty info values require explicit typed metadata encoding.
/// Molecular bonds, formal charges, isotope and stereo are not projected.
#[wasm_bindgen]
pub fn to_extxyz_json(
    mol: &MolHandle,
    coords_json: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let frame = extxyz_frame_from_json_args(&mol.inner, coords_json, options_json)
        .map_err(|e| JsValue::from_str(&e))?;
    chematic_mol::write_extxyz(&frame).map_err(|e| JsValue::from_str(&e.to_string()))
}
