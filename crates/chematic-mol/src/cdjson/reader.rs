use super::{CdjsonDocument, CdjsonError, CdjsonMolecule, CdjsonParseLimits, limit};
use crate::cjson::{CjsonParseLimits, check_json_limits};
use crate::cjson_record::UniqueValue;
use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
use serde_json::{Map, Value};
use std::collections::HashSet;

type Object = Map<String, Value>;

fn object<'a>(v: &'a Value, allowed: &[&str]) -> Result<&'a Object, CdjsonError> {
    let o = v
        .as_object()
        .ok_or_else(|| CdjsonError::Invalid("expected object".into()))?;
    for key in o.keys() {
        if !allowed.contains(&key.as_str()) {
            return Err(CdjsonError::Unsupported(format!("field {key}")));
        }
    }
    Ok(o)
}

fn array<'a>(o: &'a Object, key: &str, required: bool) -> Result<&'a [Value], CdjsonError> {
    match o.get(key) {
        None if !required => Ok(&[]),
        Some(Value::Array(values)) => Ok(values),
        _ => Err(CdjsonError::Invalid(format!("{key} must be an array"))),
    }
}

fn integer(o: &Object, key: &str, default: i64) -> Result<i64, CdjsonError> {
    match o.get(key) {
        None => Ok(default),
        Some(v) => v
            .as_i64()
            .ok_or_else(|| CdjsonError::Invalid(format!("{key} must be an integer"))),
    }
}

fn coordinate(o: &Object, key: &str, optional: bool) -> Result<f64, CdjsonError> {
    match o.get(key) {
        None if optional => Ok(0.),
        Some(v) => v
            .as_f64()
            .filter(|f| f.is_finite())
            .ok_or_else(|| CdjsonError::Invalid(format!("{key} must be a finite coordinate"))),
        _ => Err(CdjsonError::Invalid(format!("missing coordinate {key}"))),
    }
}

fn id(o: &Object, ids: &mut HashSet<String>) -> Result<Option<String>, CdjsonError> {
    match o.get("i") {
        None => Ok(None),
        Some(Value::String(s)) if !s.is_empty() => {
            if !ids.insert(s.clone()) {
                return Err(CdjsonError::Invalid(format!("duplicate ID {s}")));
            }
            Ok(Some(s.clone()))
        }
        _ => Err(CdjsonError::Invalid("i must be a nonempty string".into())),
    }
}

fn parse_atom(o: &Object) -> Result<Atom, CdjsonError> {
    let label = match o.get("l") {
        None => "C",
        Some(Value::String(s)) => s.as_str(),
        _ => return Err(CdjsonError::Invalid("l must be an element symbol".into())),
    };
    let element = Element::from_symbol(label)
        .ok_or_else(|| CdjsonError::Unsupported(format!("atom label {label}")))?;
    let mut atom = Atom::new(element);
    atom.charge = i8::try_from(integer(o, "c", 0)?)
        .map_err(|_| CdjsonError::Invalid("c outside i8 range".into()))?;
    let mass = integer(o, "m", -1)?;
    atom.isotope =
        if mass == -1 {
            None
        } else {
            Some(u16::try_from(mass).ok().filter(|&m| m > 0).ok_or_else(|| {
                CdjsonError::Invalid("m must be -1 or a positive u16 mass".into())
            })?)
        };
    let hydrogen = integer(o, "h", -1)?;
    atom.hydrogen_count = if hydrogen == -1 {
        None
    } else {
        Some(
            u8::try_from(hydrogen)
                .map_err(|_| CdjsonError::Invalid("h must be -1 or in [0,255]".into()))?,
        )
    };
    for key in ["r", "p"] {
        if integer(o, key, 0)? != 0 {
            return Err(CdjsonError::Unsupported(format!("nonzero {key}")));
        }
    }
    Ok(atom)
}

fn parse_bonds(
    bonds: &[Value],
    count: usize,
    builder: &mut MoleculeBuilder,
    ids: &mut HashSet<String>,
) -> Result<Vec<Option<String>>, CdjsonError> {
    let mut bond_ids = Vec::new();
    for bond in bonds {
        let o = object(bond, &["b", "e", "o", "s", "i"])?;
        let begin = integer(o, "b", -1)?;
        let end = integer(o, "e", -1)?;
        if begin < 0
            || end < 0
            || begin as u64 >= count as u64
            || end as u64 >= count as u64
            || begin == end
        {
            return Err(CdjsonError::Invalid(
                "bond endpoints must be distinct in-range atom indices".into(),
            ));
        }
        let order = match o.get("o").map(Value::as_f64).unwrap_or(Some(1.)) {
            Some(1.) => BondOrder::Single,
            Some(2.) => BondOrder::Double,
            Some(3.) => BondOrder::Triple,
            _ => return Err(CdjsonError::Unsupported("bond order outside 1/2/3".into())),
        };
        if o.get("s").is_some_and(|v| v.as_str() != Some("none")) {
            return Err(CdjsonError::Unsupported("bond stereochemistry".into()));
        }
        builder
            .add_bond(AtomIdx(begin as u32), AtomIdx(end as u32), order)
            .map_err(|e| CdjsonError::Invalid(e.to_string()))?;
        bond_ids.push(id(o, ids)?);
    }
    Ok(bond_ids)
}

fn parse_molecule(v: &Value, ids: &mut HashSet<String>) -> Result<CdjsonMolecule, CdjsonError> {
    let o = object(v, &["a", "b", "i"])?;
    let molecule_id = id(o, ids)?;
    let mut builder = MoleculeBuilder::new();
    let mut coordinates = Vec::new();
    let mut atom_ids = Vec::new();
    let atoms = array(o, "a", true)?;
    for value in atoms {
        let atom = object(value, &["x", "y", "z", "l", "c", "m", "h", "r", "p", "i"])?;
        builder.add_atom(parse_atom(atom)?);
        coordinates.push([
            coordinate(atom, "x", false)?,
            coordinate(atom, "y", false)?,
            coordinate(atom, "z", true)?,
        ]);
        atom_ids.push(id(atom, ids)?);
    }
    let bond_ids = parse_bonds(array(o, "b", false)?, atoms.len(), &mut builder, ids)?;
    Ok(CdjsonMolecule {
        molecule: builder.build(),
        coordinates,
        id: molecule_id,
        atom_ids,
        bond_ids,
    })
}

/// Read either one molecular root or a content root; never select only the first.
pub fn parse_cdjson(input: &str) -> Result<CdjsonDocument, CdjsonError> {
    parse_cdjson_with_limits(input, &CdjsonParseLimits::default())
}

pub fn parse_cdjson_with_limits(
    input: &str,
    limits: &CdjsonParseLimits,
) -> Result<CdjsonDocument, CdjsonError> {
    let hard = CdjsonParseLimits::default();
    limit(
        "input bytes",
        input.len(),
        limits.max_input_bytes,
        hard.max_input_bytes,
    )?;
    let value: UniqueValue =
        serde_json::from_str(input).map_err(|e| CdjsonError::Invalid(e.to_string()))?;
    let value = value.0;
    let json_limits = CjsonParseLimits {
        max_json_depth: 32,
        max_array_items: 131_072,
        max_string_bytes: 65_536,
        ..Default::default()
    };
    check_json_limits(&value, 0, &json_limits).map_err(|e| CdjsonError::Invalid(e.to_string()))?;
    let values: &[Value] = if value.get("a").is_some() {
        std::slice::from_ref(&value)
    } else {
        let o = object(&value, &["m", "s"])?;
        if !array(o, "s", false)?.is_empty() {
            return Err(CdjsonError::Unsupported("shapes".into()));
        }
        array(o, "m", false)?
    };
    limit(
        "molecules",
        values.len(),
        limits.max_molecules,
        hard.max_molecules,
    )?;
    let mut atoms = 0usize;
    let mut bonds = 0usize;
    for v in values {
        let o = object(v, &["a", "b", "i"])?;
        atoms += array(o, "a", true)?.len();
        bonds += array(o, "b", false)?.len();
        limit("atoms", atoms, limits.max_atoms, hard.max_atoms)?;
        limit("bonds", bonds, limits.max_bonds, hard.max_bonds)?;
    }
    let mut ids = HashSet::new();
    let molecules = values
        .iter()
        .map(|v| parse_molecule(v, &mut ids))
        .collect::<Result<_, _>>()?;
    Ok(CdjsonDocument {
        source: input.to_owned(),
        molecules,
    })
}
