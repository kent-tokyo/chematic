//! Complete preflight of hand-built frames, before the streaming sink is touched.
use super::{XyzFrame, XyzPropertyKind, XyzValue};
use std::collections::HashSet;

fn ambiguous_character(c: char) -> bool {
    c.is_whitespace() || matches!(c, '=' | '"' | '\'' | '\\' | '[' | ']' | '{' | '}')
}

fn valid_value(kind: XyzPropertyKind, value: &XyzValue) -> bool {
    match (kind, value) {
        (XyzPropertyKind::Real, XyzValue::Real(v)) => v.is_finite(),
        (XyzPropertyKind::Integer, XyzValue::Integer(_))
        | (XyzPropertyKind::Logical, XyzValue::Logical(_)) => true,
        (XyzPropertyKind::String, XyzValue::Str(v)) => {
            !v.is_empty() && !v.chars().any(|c| c.is_control() || ambiguous_character(c))
        }
        _ => false,
    }
}

fn validate_info(info: &[(String, String)]) -> Result<(), String> {
    let mut keys = HashSet::new();
    for (key, value) in info {
        if value.is_empty() {
            return Err(format!(
                "info value for {key:?} is empty; use an explicit typed metadata encoding"
            ));
        }
        if !keys.insert(key) {
            return Err(format!("duplicate info key {key:?}"));
        }
        if key
            .chars()
            .any(|c| c.is_control() || ambiguous_character(c))
        {
            return Err(format!("info key {key:?} has ambiguous header syntax"));
        }
        if value.chars().any(|c| c.is_control() && c != '\t') {
            return Err(format!(
                "info value for {key:?} has an unsupported control character"
            ));
        }
    }
    Ok(())
}

pub(super) fn validate(frame: &XyzFrame) -> Result<(), String> {
    for (i, atom) in frame.atoms.iter().enumerate() {
        if ![atom.x, atom.y, atom.z].into_iter().all(f64::is_finite) {
            return Err(format!("atom {i} has nonfinite coordinates"));
        }
    }
    if frame
        .lattice
        .is_some_and(|l| !l.into_iter().all(f64::is_finite))
    {
        return Err("lattice has nonfinite values".into());
    }
    validate_info(&frame.info)?;
    let mut names = HashSet::from(["species", "pos", "symbols", "positions", "numbers", "Z"]);
    let mut columns = 4usize;
    for prop in &frame.properties {
        // ASE maps charge to charges; also refuse aliases that would overwrite
        // the species/position columns rather than become independent data.
        let canonical = if prop.name == "charge" {
            "charges"
        } else {
            prop.name.as_str()
        };
        if !names.insert(canonical) {
            return Err(format!(
                "duplicate or reserved property name {:?}",
                prop.name
            ));
        }
        if prop
            .name
            .chars()
            .any(|c| c.is_control() || ambiguous_character(c))
        {
            return Err(format!(
                "property name {:?} has ambiguous header syntax",
                prop.name
            ));
        }
        columns = columns
            .checked_add(prop.count)
            .ok_or_else(|| "property column count overflow".to_owned())?;
        if prop.values.len() != frame.atoms.len() {
            return Err(format!(
                "property {:?} has {} rows; expected {}",
                prop.name,
                prop.values.len(),
                frame.atoms.len()
            ));
        }
        for (i, row) in prop.values.iter().enumerate() {
            if row.len() != prop.count || !row.iter().all(|v| valid_value(prop.kind, v)) {
                return Err(format!(
                    "property {:?} atom {i} has invalid component count, type or value",
                    prop.name
                ));
            }
        }
    }
    if frame.lattice.is_none()
        && frame.properties.is_empty()
        && frame.info.is_empty()
        && frame.comment.chars().any(|c| c.is_control() && c != '\t')
    {
        return Err("comment contains an unsupported control character".into());
    }
    Ok(())
}
