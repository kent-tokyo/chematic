//! Minimal CXSMILES/CXSMARTS metadata support.
//!
//! This module preserves the CX fields that are commonly lost during ordinary
//! SMILES round-trips: atom labels, atom properties, atom radicals, and
//! zero-order bonds (`Z:`). The molecular graph remains a normal [`Molecule`];
//! CX-only atom metadata is carried by [`CxSmiles`].

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

use crate::{SmilesError, parse, write};

/// One CX atom property entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxAtomProp {
    /// Atom index in molecule atom order.
    pub atom: AtomIdx,
    /// Property key.
    pub key: String,
    /// Property value.
    pub value: String,
}

/// A parsed CXSMILES molecule plus CX metadata.
pub struct CxSmiles {
    /// Parsed molecular graph.
    pub mol: Molecule,
    /// Optional atom labels indexed by atom order.
    pub atom_labels: Vec<Option<String>>,
    /// Atom property triples.
    pub atom_props: Vec<CxAtomProp>,
    /// Optional radical class indexed by atom order.
    pub atom_radicals: Vec<Option<u8>>,
    /// Bond indices (0-based) marked as "wavy" (STEREOANY / unknown E/Z stereo).
    ///
    /// Corresponds to the `|w:N,...|` field in CXSMILES (RDKit PR #9316).
    /// A wavy bond means the double bond has unspecified or unknown E/Z geometry
    /// (drawn as a wavy line in 2D depictions rather than a cis/trans wedge).
    pub wavy_bonds: Vec<BondIdx>,
}

impl CxSmiles {
    fn new(mol: Molecule) -> Self {
        let n = mol.atom_count();
        Self {
            mol,
            atom_labels: vec![None; n],
            atom_props: Vec::new(),
            atom_radicals: vec![None; n],
            wavy_bonds: Vec::new(),
        }
    }
}

/// Parse a CXSMILES string.
///
/// Unknown CX fields are ignored so that newer RDKit CX blocks remain readable.
pub fn parse_cxsmiles(input: &str) -> Result<CxSmiles, SmilesError> {
    let (base, cx) = split_cx(input);
    let mut cxmol = CxSmiles::new(parse(base.trim())?);
    if let Some(cx) = cx {
        parse_cx_block(cx, &mut cxmol);
    }
    Ok(cxmol)
}

/// Write a CXSMILES string, including supported CX fields when present.
pub fn write_cxsmiles(cx: &CxSmiles) -> String {
    let mut base = write(&cx.mol);
    let mut fields = Vec::new();

    if cx.atom_labels.iter().any(|label| label.is_some()) {
        let labels = (0..cx.mol.atom_count())
            .map(|i| {
                cx.atom_labels
                    .get(i)
                    .and_then(|v| v.as_deref())
                    .unwrap_or("")
            })
            .map(escape_cx_value)
            .collect::<Vec<_>>()
            .join(";");
        fields.push(format!("${labels}$"));
    }

    if !cx.atom_props.is_empty() {
        let props = cx
            .atom_props
            .iter()
            .map(|p| {
                format!(
                    "{}.{}.{}",
                    p.atom.0,
                    escape_cx_value(&p.key),
                    escape_cx_value(&p.value)
                )
            })
            .collect::<Vec<_>>()
            .join(":");
        fields.push(format!("atomProp:{props}"));
    }

    for class in 1..=7u8 {
        let atoms = cx
            .atom_radicals
            .iter()
            .enumerate()
            .filter_map(|(i, radical)| (*radical == Some(class)).then_some(i.to_string()))
            .collect::<Vec<_>>();
        if !atoms.is_empty() {
            fields.push(format!("^{class}:{}", atoms.join(",")));
        }
    }

    let zero_bonds = cx
        .mol
        .bonds()
        .filter_map(|(bidx, bond)| (bond.order == BondOrder::Zero).then_some(bidx.0.to_string()))
        .collect::<Vec<_>>();
    if !zero_bonds.is_empty() {
        fields.push(format!("Z:{}", zero_bonds.join(",")));
    }

    if !cx.wavy_bonds.is_empty() {
        let indices = cx
            .wavy_bonds
            .iter()
            .map(|bidx| bidx.0.to_string())
            .collect::<Vec<_>>()
            .join(",");
        fields.push(format!("w:{indices}"));
    }

    if !fields.is_empty() {
        base.push_str(" |");
        base.push_str(&fields.join(","));
        base.push('|');
    }
    base
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

fn parse_cx_block(cx: &str, out: &mut CxSmiles) {
    // Track the "active" multi-value prefix so that bare numeric tokens produced
    // when split_cx_fields splits `w:0,2` → `["w:0", "2"]` are attributed to the
    // correct handler (wavy bonds or zero bonds).
    let mut last_mv: Option<char> = None;

    for field in split_cx_fields(cx) {
        if field.starts_with('$') && field.ends_with('$') {
            parse_labels(&field[1..field.len() - 1], out);
            last_mv = None;
        } else if let Some(rest) = field.strip_prefix("atomProp:") {
            parse_atom_props(rest, out);
            last_mv = None;
        } else if let Some(rest) = field.strip_prefix('Z').and_then(|s| s.strip_prefix(':')) {
            parse_zero_bonds(rest, out);
            last_mv = Some('Z');
        } else if let Some(rest) = field.strip_prefix('w').and_then(|s| s.strip_prefix(':')) {
            parse_wavy_bonds(rest, out);
            last_mv = Some('w');
        } else if let Some(rest) = field.strip_prefix('^') {
            parse_radicals(rest, out);
            last_mv = None;
        } else if !field.is_empty() && field.trim().chars().all(|c| c.is_ascii_digit()) {
            // Bare digit token — continuation of a split multi-value field.
            // split_cx_fields splits `w:0,2` at the comma, yielding `["w:0","2"]`.
            match last_mv {
                Some('w') => parse_wavy_bonds(field.trim(), out),
                Some('Z') => parse_zero_bonds(field.trim(), out),
                _ => last_mv = None,
            }
        } else {
            last_mv = None;
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

fn parse_labels(labels: &str, out: &mut CxSmiles) {
    for (i, label) in labels.split(';').enumerate().take(out.atom_labels.len()) {
        if !label.is_empty() {
            out.atom_labels[i] = Some(unescape_cx_value(label));
        }
    }
}

fn parse_atom_props(props: &str, out: &mut CxSmiles) {
    for prop in props.split(':') {
        let mut parts = prop.splitn(3, '.');
        let Some(atom_raw) = parts.next() else {
            continue;
        };
        let Some(key) = parts.next() else { continue };
        let Some(value) = parts.next() else { continue };
        let Ok(atom) = atom_raw.parse::<u32>() else {
            continue;
        };
        if atom as usize >= out.mol.atom_count() {
            continue;
        }
        out.atom_props.push(CxAtomProp {
            atom: AtomIdx(atom),
            key: unescape_cx_value(key),
            value: unescape_cx_value(value),
        });
    }
}

fn parse_zero_bonds(rest: &str, out: &mut CxSmiles) {
    let mut mol = chematic_core::MoleculeBuilder::from_molecule(&out.mol).build();
    for item in rest.split(',') {
        let Ok(idx) = item.trim().parse::<u32>() else {
            continue;
        };
        if (idx as usize) < mol.bond_count() {
            mol = mol.with_bond_order(BondIdx(idx), BondOrder::Zero);
        }
    }
    out.mol = mol;
}

/// Parse `|w:N,M,...|` — bond indices of wavy (STEREOANY) double bonds.
fn parse_wavy_bonds(rest: &str, out: &mut CxSmiles) {
    for item in rest.split(',') {
        let Ok(idx) = item.trim().parse::<u32>() else {
            continue;
        };
        if (idx as usize) < out.mol.bond_count() {
            out.wavy_bonds.push(BondIdx(idx));
        }
    }
}

fn parse_radicals(rest: &str, out: &mut CxSmiles) {
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

fn escape_cx_value(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace(';', "\\;")
        .replace(',', "\\,")
        .replace('|', "\\|")
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

    #[test]
    fn parse_cxsmiles_atom_labels_and_props() {
        let cx = parse_cxsmiles("CO |$C1;O2$,atomProp:0.p1.5:1.note.acceptor|").unwrap();
        assert_eq!(cx.mol.atom_count(), 2);
        assert_eq!(cx.atom_labels[0].as_deref(), Some("C1"));
        assert_eq!(cx.atom_labels[1].as_deref(), Some("O2"));
        assert_eq!(cx.atom_props.len(), 2);
        assert_eq!(cx.atom_props[0].key, "p1");
        assert_eq!(cx.atom_props[0].value, "5");
    }

    #[test]
    fn parse_cxsmiles_zero_bond() {
        let cx = parse_cxsmiles("C~O |Z:0|").unwrap();
        assert_eq!(cx.mol.bond(BondIdx(0)).order, BondOrder::Zero);
        let out = write_cxsmiles(&cx);
        assert!(out.contains("Z:0"), "{out}");
    }

    #[test]
    fn parse_cxsmiles_radicals() {
        let cx = parse_cxsmiles("[CH3] |^2:0|").unwrap();
        assert_eq!(cx.atom_radicals[0], Some(2));
        assert!(write_cxsmiles(&cx).contains("^2:0"));
    }

    #[test]
    fn bug1_cxsmiles_trailing_backslash_preservation() {
        // BUG #1 Fix Verification: trailing backslashes in labels should be preserved
        // Input: "label\\" (escaped form) -> should parse to "label\" (single backslash)
        let cx = parse_cxsmiles(r#"CO |$label\\;O$|"#).unwrap();
        assert_eq!(
            cx.atom_labels[0].as_deref(),
            Some("label\\"),
            "Trailing backslash should be preserved after unescape"
        );

        // Round-trip test: parsed label should serialize back correctly
        let serialized = write_cxsmiles(&cx);
        let cx2 = parse_cxsmiles(&serialized).unwrap();
        assert_eq!(
            cx2.atom_labels[0], cx.atom_labels[0],
            "Trailing backslash should round-trip correctly"
        );
    }

    #[test]
    fn bug1_cxsmiles_double_trailing_backslash() {
        // More complex case: label with backslash in middle and at end
        // "C\\label\\" in escaped form -> "C\label\" after unescape
        let cx = parse_cxsmiles(r#"CO |$C\\label\\;O$|"#).unwrap();
        assert_eq!(
            cx.atom_labels[0].as_deref(),
            Some("C\\label\\"),
            "Both backslashes should be preserved"
        );

        // Verify round-trip
        let serialized = write_cxsmiles(&cx);
        let cx2 = parse_cxsmiles(&serialized).unwrap();
        assert_eq!(
            cx2.atom_labels[0], cx.atom_labels[0],
            "Double backslash pattern should round-trip"
        );
    }

    #[test]
    fn bug1_cxsmiles_escaped_comma_with_trailing_backslash() {
        // Test that escaped comma doesn't interfere with trailing backslash handling
        // Label with escaped comma and trailing backslash
        let cx = parse_cxsmiles(r#"CO |$label\,end\\;O$|"#).unwrap();
        assert_eq!(
            cx.atom_labels[0].as_deref(),
            Some("label,end\\"),
            "Escaped comma and trailing backslash should both be preserved"
        );
    }

    // ── RDKit PR #9316: CXSMILES STEREOANY / wavy bonds ─────────────────────

    #[test]
    fn parse_cxsmiles_wavy_bond_single() {
        // |w:1| marks bond index 1 as a wavy (STEREOANY) double bond.
        let cx = parse_cxsmiles("CC=CC |w:1|").unwrap();
        assert_eq!(cx.wavy_bonds.len(), 1, "should find one wavy bond");
        assert_eq!(cx.wavy_bonds[0], BondIdx(1), "wavy bond index must be 1");
    }

    #[test]
    fn parse_cxsmiles_wavy_bond_multiple() {
        // |w:0,2| marks two bonds as wavy.
        let cx = parse_cxsmiles("C/C=C/C=CC |w:0,2|").unwrap();
        assert_eq!(cx.wavy_bonds.len(), 2);
        assert!(cx.wavy_bonds.contains(&BondIdx(0)));
        assert!(cx.wavy_bonds.contains(&BondIdx(2)));
    }

    #[test]
    fn cxsmiles_wavy_bond_round_trip() {
        // Parsing |w:1| and writing back must preserve the marker.
        let input = "CC=CC |w:1|";
        let cx = parse_cxsmiles(input).unwrap();
        let out = write_cxsmiles(&cx);
        assert!(out.contains("w:1"), "wavy bond marker must round-trip: {out}");
    }

    #[test]
    fn cxsmiles_no_wavy_bond_no_w_field() {
        // If no wavy bonds are present, |w:| must not appear in output.
        let cx = parse_cxsmiles("CC=CC").unwrap();
        let out = write_cxsmiles(&cx);
        assert!(!out.contains("w:"), "no wavy bonds → no |w:| in output: {out}");
    }

    #[test]
    fn parse_cxsmiles_wavy_bond_out_of_range_ignored() {
        // Bond index beyond bond count must be silently ignored.
        let cx = parse_cxsmiles("CO |w:99|").unwrap();
        assert!(cx.wavy_bonds.is_empty(), "out-of-range bond index must be ignored");
    }

    #[test]
    fn cxsmiles_wavy_bond_combined_with_other_fields() {
        // Wavy bonds should coexist with atom labels and other CX fields.
        let input = "CC=CC |$a;b;c;d$,w:1|";
        let cx = parse_cxsmiles(input).unwrap();
        assert_eq!(cx.atom_labels[0].as_deref(), Some("a"));
        assert_eq!(cx.wavy_bonds, vec![BondIdx(1)]);
        let out = write_cxsmiles(&cx);
        assert!(out.contains("w:1"), "wavy bond preserved with other fields: {out}");
    }
}
