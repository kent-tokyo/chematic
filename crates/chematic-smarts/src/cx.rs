//! Minimal CXSMARTS wrapper.
//!
//! CXSMARTS uses the same trailing `|...|` metadata envelope as CXSMILES.
//! The query graph is parsed by the normal SMARTS parser and supported atom
//! metadata is preserved here instead of being discarded.

use crate::{QueryMolecule, SmartsError, parse_smarts};

/// One CXSMARTS atom property entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxQueryAtomProp {
    pub atom: usize,
    pub key: String,
    pub value: String,
}

/// Parsed SMARTS query plus CX metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct CxSmarts {
    pub query: QueryMolecule,
    pub atom_labels: Vec<Option<String>>,
    pub atom_props: Vec<CxQueryAtomProp>,
    pub atom_radicals: Vec<Option<u8>>,
}

impl CxSmarts {
    fn new(query: QueryMolecule) -> Self {
        let n = query.atom_count();
        Self {
            query,
            atom_labels: vec![None; n],
            atom_props: Vec::new(),
            atom_radicals: vec![None; n],
        }
    }
}

/// Parse CXSMARTS, preserving supported CX metadata.
pub fn parse_cxsmarts(input: &str) -> Result<CxSmarts, SmartsError> {
    let (base, cx) = split_cx(input);
    let mut out = CxSmarts::new(parse_smarts(base.trim())?);
    if let Some(cx) = cx {
        parse_cx_block(cx, &mut out);
    }
    Ok(out)
}

fn split_cx(input: &str) -> (&str, Option<&str>) {
    let trimmed = input.trim();
    if let Some(start) = trimmed.find('|')
        && let Some(end_rel) = trimmed[start + 1..].find('|')
    {
        let end = start + 1 + end_rel;
        return (&trimmed[..start], Some(&trimmed[start + 1..end]));
    }
    (trimmed, None)
}

fn parse_cx_block(cx: &str, out: &mut CxSmarts) {
    for field in split_cx_fields(cx) {
        if field.len() >= 2 && field.starts_with('$') && field.ends_with('$') {
            parse_labels(&field[1..field.len() - 1], out);
        } else if let Some(rest) = field.strip_prefix("atomProp:") {
            parse_atom_props(rest, out);
        } else if let Some(rest) = field.strip_prefix('^') {
            parse_radicals(rest, out);
        }
    }
}

fn split_cx_fields(cx: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut current = String::new();
    let mut in_labels = false;
    for ch in cx.chars() {
        match ch {
            '$' => {
                in_labels = !in_labels;
                current.push(ch);
            }
            ',' if !in_labels => {
                if !current.is_empty() {
                    fields.push(current.trim().to_string());
                    current.clear();
                }
            }
            _ => current.push(ch),
        }
    }
    if !current.is_empty() {
        fields.push(current.trim().to_string());
    }
    fields
}

fn parse_labels(labels: &str, out: &mut CxSmarts) {
    for (i, label) in labels.split(';').enumerate().take(out.atom_labels.len()) {
        if !label.is_empty() {
            out.atom_labels[i] = Some(unescape_cx_value(label));
        }
    }
}

fn parse_atom_props(props: &str, out: &mut CxSmarts) {
    for prop in props.split(':') {
        let mut parts = prop.splitn(3, '.');
        let Some(atom_raw) = parts.next() else {
            continue;
        };
        let Some(key) = parts.next() else { continue };
        let Some(value) = parts.next() else { continue };
        let Ok(atom) = atom_raw.parse::<usize>() else {
            continue;
        };
        if atom >= out.query.atom_count() {
            continue;
        }
        out.atom_props.push(CxQueryAtomProp {
            atom,
            key: unescape_cx_value(key),
            value: unescape_cx_value(value),
        });
    }
}

fn parse_radicals(rest: &str, out: &mut CxSmarts) {
    let Some((class_raw, atoms_raw)) = rest.split_once(':') else {
        return;
    };
    let Ok(class) = class_raw.parse::<u8>() else {
        return;
    };
    for atom_raw in atoms_raw.split(',') {
        let Ok(atom) = atom_raw.trim().parse::<usize>() else {
            continue;
        };
        if atom < out.atom_radicals.len() {
            out.atom_radicals[atom] = Some(class);
        }
    }
}

fn unescape_cx_value(value: &str) -> String {
    let mut out = String::new();
    let mut escape = false;
    for ch in value.chars() {
        if escape {
            out.push(ch);
            escape = false;
        } else if ch == '\\' {
            escape = true;
        } else {
            out.push(ch);
        }
    }
    if escape {
        out.push('\\');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondPrimitive, BondQuery};

    #[test]
    fn parse_cxsmarts_preserves_metadata() {
        let cx = parse_cxsmarts("[#6]~[#8] |$C1;O2$,atomProp:1.role.acceptor,^2:0|").unwrap();
        assert_eq!(cx.query.atom_count(), 2);
        assert_eq!(cx.atom_labels[0].as_deref(), Some("C1"));
        assert_eq!(cx.atom_labels[1].as_deref(), Some("O2"));
        assert_eq!(cx.atom_props[0].key, "role");
        assert_eq!(cx.atom_radicals[0], Some(2));
        assert!(matches!(
            cx.query.bonds[0].query,
            BondQuery::Primitive(BondPrimitive::Any) | BondQuery::Any
        ));
    }

    #[test]
    fn bug1_cxsmarts_trailing_backslash_preservation() {
        // BUG #1 Fix Verification: trailing backslashes in labels should be preserved
        // Input: "label\\" (escaped form) -> should parse to "label\" (single backslash)
        let cx = parse_cxsmarts(r#"[#6] |$label\\$|"#).unwrap();
        assert_eq!(
            cx.atom_labels[0].as_deref(),
            Some("label\\"),
            "Trailing backslash should be preserved after unescape"
        );
    }

    #[test]
    fn bug1_cxsmarts_double_trailing_backslash() {
        // More complex case: label with backslash in middle and at end
        // "C\\label\\" in escaped form -> "C\label\" after unescape
        let cx = parse_cxsmarts(r#"[#6] |$C\\label\\$|"#).unwrap();
        assert_eq!(
            cx.atom_labels[0].as_deref(),
            Some("C\\label\\"),
            "Both backslashes should be preserved"
        );
    }

    #[test]
    fn bug1_cxsmarts_escaped_comma_with_trailing_backslash() {
        // Test that escaped comma doesn't interfere with trailing backslash handling
        // Label with escaped comma and trailing backslash
        let cx = parse_cxsmarts(r#"[#6] |$label\,end\\$|"#).unwrap();
        assert_eq!(
            cx.atom_labels[0].as_deref(),
            Some("label,end\\"),
            "Escaped comma and trailing backslash should both be preserved"
        );
    }

    #[test]
    fn empty_cx_label_field_is_ignored_without_panicking() {
        let cx = parse_cxsmarts("C |$|").unwrap();
        assert_eq!(cx.query.atom_count(), 1);
        assert!(cx.atom_labels.iter().all(Option::is_none));
    }
}
