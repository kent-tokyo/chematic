//! ChemDraw XML (CDXML) parser — **read-only**.
//!
//! CDXML is a proprietary XML format produced by ChemDraw (PerkinElmer /
//! Revvity).  This parser handles the minimal subset needed to extract
//! molecular structure from a CDXML document.
//!
//! # Supported elements / attributes
//!
//! `<n>` (atom): `id`, `Element` (atomic number), `p` ("x y" 2D coords),
//!               `NumHydrogens`, `Charge`, `Isotope`
//!
//! `<b>` (bond): `B` (begin atom id), `E` (end atom id), `Order` (1/2/3,
//!               defaults to 1)
//!
//! # Limitations
//!
//! - Only the first `<fragment>` in the document is returned as a single
//!   molecule.  Multi-molecule documents (multiple fragments) are not yet
//!   supported.
//! - Write support is **not** implemented (CDXML is a proprietary format;
//!   writing files that ChemDraw will accept requires undocumented attributes).
//! - Stereochemistry attributes are ignored.
//! - Presentation-only nodes (text boxes, arrows, etc.) are silently skipped.

use std::collections::HashMap;

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

use crate::cml::parse_xml_attrs;

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Error returned when parsing a CDXML document fails.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CdxmlError {
    /// An atom `Element` attribute contained an unknown atomic number.
    UnknownAtomicNumber(u32),
    /// A bond referenced an atom id that was not defined.
    UnknownAtomRef(String),
    /// A `<b>` bond element is missing a `B` or `E` attribute.
    MissingBondEndpoint,
    /// The `p` coordinate attribute could not be parsed.
    InvalidCoords(String),
}

impl std::fmt::Display for CdxmlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CdxmlError::UnknownAtomicNumber(n) => write!(f, "unknown atomic number: {n}"),
            CdxmlError::UnknownAtomRef(s)      => write!(f, "unknown atom ref: {s}"),
            CdxmlError::MissingBondEndpoint    => write!(f, "bond missing B or E attribute"),
            CdxmlError::InvalidCoords(s)       => write!(f, "invalid p coords: {s}"),
        }
    }
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a CDXML document into a `(Molecule, 2D-coords)` pair.
///
/// The second element of the tuple contains one `(x, y)` entry per atom in
/// atom-insertion order.  Coordinates are in CDXML point units (1/72 inch);
/// no conversion to Ångströms is applied.
pub fn parse_cdxml(input: &str) -> Result<(Molecule, Vec<(f64, f64)>), CdxmlError> {
    struct AtomData {
        element: Element,
        charge: i8,
        isotope: Option<u16>,
        hydrogen_count: Option<u8>,
        x: f64,
        y: f64,
    }

    struct BondDataRaw { b: String, e: String, order: BondOrder }

    let mut atom_by_id: HashMap<String, (usize, AtomData)> = HashMap::new();
    let mut bonds_raw: Vec<BondDataRaw> = Vec::new();

    for raw_line in input.lines() {
        let line = raw_line.trim();
        if line.is_empty() { continue; }

        // Atom node <n ...>
        if is_n_tag(line) {
            let attrs = parse_xml_attrs(line);

            let id = match attrs.get("id") {
                Some(s) => s.clone(),
                None => continue, // skip nodes without id (non-atom nodes)
            };

            // Atomic number: default to Carbon (6) if missing.
            let element_num: u32 = attrs.get("Element")
                .and_then(|s| s.trim().parse::<u32>().ok())
                .unwrap_or(6);
            let element = Element::from_atomic_number(element_num as u8)
                .ok_or(CdxmlError::UnknownAtomicNumber(element_num))?;

            let charge = attrs.get("Charge")
                .and_then(|s| s.trim().parse::<i8>().ok())
                .unwrap_or(0);
            let isotope = attrs.get("Isotope")
                .and_then(|s| s.trim().parse::<u16>().ok())
                .filter(|&v| v > 0);
            let hydrogen_count = attrs.get("NumHydrogens")
                .and_then(|s| s.trim().parse::<u8>().ok());

            let (x, y) = if let Some(p) = attrs.get("p") {
                let parts: Vec<&str> = p.split_whitespace().collect();
                if parts.len() < 2 {
                    return Err(CdxmlError::InvalidCoords(p.clone()));
                }
                let x = parts[0].parse::<f64>().unwrap_or(0.0);
                let y = parts[1].parse::<f64>().unwrap_or(0.0);
                (x, y)
            } else {
                (0.0, 0.0)
            };

            let pos = atom_by_id.len();
            atom_by_id.insert(id, (pos, AtomData { element, charge, isotope, hydrogen_count, x, y }));
            continue;
        }

        // Bond element <b ...>
        if is_b_tag(line) {
            let attrs = parse_xml_attrs(line);
            let b = attrs.get("B").cloned().ok_or(CdxmlError::MissingBondEndpoint)?;
            let e = attrs.get("E").cloned().ok_or(CdxmlError::MissingBondEndpoint)?;
            let order: BondOrder = match attrs.get("Order").map(String::as_str) {
                Some("2") => BondOrder::Double,
                Some("3") => BondOrder::Triple,
                _         => BondOrder::Single,
            };
            bonds_raw.push(BondDataRaw { b, e, order });
        }
    }

    // Build molecule in insertion order (pos).
    let mut ordered: Vec<(&str, &AtomData)> = atom_by_id.iter()
        .map(|(id, (_, data))| (id.as_str(), data))
        .collect();
    // Sort by insertion order via position lookup.
    let mut pos_ordered: Vec<(&str, usize, &AtomData)> = atom_by_id.iter()
        .map(|(id, (pos, data))| (id.as_str(), *pos, data))
        .collect();
    pos_ordered.sort_by_key(|&(_, pos, _)| pos);

    let mut builder = MoleculeBuilder::new();
    let mut idx_by_id: HashMap<String, AtomIdx> = HashMap::new();
    let mut coords: Vec<(f64, f64)> = Vec::new();

    // Suppress unused variable warning for `ordered`
    let _ = ordered.len();

    for (id, _, data) in &pos_ordered {
        let mut atom = Atom::new(data.element);
        atom.charge = data.charge;
        atom.isotope = data.isotope;
        atom.hydrogen_count = data.hydrogen_count;
        let idx = builder.add_atom(atom);
        idx_by_id.insert(id.to_string(), idx);
        coords.push((data.x, data.y));
    }

    for bd in &bonds_raw {
        let a1 = *idx_by_id.get(&bd.b)
            .ok_or_else(|| CdxmlError::UnknownAtomRef(bd.b.clone()))?;
        let a2 = *idx_by_id.get(&bd.e)
            .ok_or_else(|| CdxmlError::UnknownAtomRef(bd.e.clone()))?;
        builder.add_bond(a1, a2, bd.order)
            .map_err(|_| CdxmlError::UnknownAtomRef(format!("{} {}", bd.b, bd.e)))?;
    }

    Ok((builder.build(), coords))
}

/// True if `line` starts a CDXML `<n` atom node tag.
fn is_n_tag(line: &str) -> bool {
    (line.starts_with("<n ") || line.starts_with("<n\t") || line == "<n>")
        && !line.starts_with("<node") // avoid accidental match on <node>
}

/// True if `line` starts a CDXML `<b` bond tag.
fn is_b_tag(line: &str) -> bool {
    line.starts_with("<b ") || line.starts_with("<b\t") || line == "<b>"
        || line.starts_with("<b/>") || line.starts_with("<b>")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Minimal hand-crafted CDXML for ethanol (C-C-O) with 2 bonds.
    const ETHANOL_CDXML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">
<CDXML>
<fragment>
<n id="1" p="10.0 20.0" Element="6" NumHydrogens="3"/>
<n id="2" p="25.0 20.0" Element="6" NumHydrogens="2"/>
<n id="3" p="40.0 20.0" Element="8" NumHydrogens="1"/>
<b B="1" E="2" Order="1"/>
<b B="2" E="3" Order="1"/>
</fragment>
</CDXML>"#;

    #[test]
    fn parse_cdxml_ethanol_atom_count() {
        let (mol, coords) = parse_cdxml(ETHANOL_CDXML).unwrap();
        assert_eq!(mol.atom_count(), 3, "ethanol: 3 heavy atoms");
        assert_eq!(mol.bond_count(), 2, "ethanol: 2 bonds");
        assert_eq!(coords.len(), 3, "one coord per atom");
    }

    #[test]
    fn parse_cdxml_ethanol_elements() {
        let (mol, _) = parse_cdxml(ETHANOL_CDXML).unwrap();
        let elems: Vec<&str> = mol.atoms().map(|(_, a)| a.element.symbol()).collect();
        assert!(elems.contains(&"C"), "should contain C");
        assert!(elems.contains(&"O"), "should contain O");
    }

    #[test]
    fn parse_cdxml_ethanol_coords() {
        let (_, coords) = parse_cdxml(ETHANOL_CDXML).unwrap();
        // Atom 1 (first C): p="10.0 20.0"
        assert!((coords[0].0 - 10.0).abs() < 0.01, "first atom x=10.0: {:?}", coords[0]);
        assert!((coords[0].1 - 20.0).abs() < 0.01, "first atom y=20.0: {:?}", coords[0]);
    }

    #[test]
    fn parse_cdxml_carbon_element_6() {
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
</fragment></CDXML>"#;
        let (mol, _) = parse_cdxml(cdxml).unwrap();
        assert_eq!(mol.atom_count(), 1);
        let atom = mol.atom(chematic_core::AtomIdx(0));
        assert_eq!(atom.element.symbol(), "C", "Element=6 → Carbon");
    }

    #[test]
    fn parse_cdxml_double_bond() {
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="8" p="10 0"/>
<b B="1" E="2" Order="2"/>
</fragment></CDXML>"#;
        let (mol, _) = parse_cdxml(cdxml).unwrap();
        let bond = mol.bond(chematic_core::BondIdx(0));
        assert_eq!(bond.order, BondOrder::Double, "Order=2 → Double");
    }

    #[test]
    fn parse_cdxml_charge() {
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="7" Charge="1" p="0 0"/>
</fragment></CDXML>"#;
        let (mol, _) = parse_cdxml(cdxml).unwrap();
        let atom = mol.atom(chematic_core::AtomIdx(0));
        assert_eq!(atom.charge, 1, "Charge=1 → N+");
    }

    #[test]
    fn parse_cdxml_unknown_atomic_number_returns_err() {
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="999" p="0 0"/>
</fragment></CDXML>"#;
        let result = parse_cdxml(cdxml);
        assert!(
            matches!(result, Err(CdxmlError::UnknownAtomicNumber(_))),
            "unknown atomic number should return Err"
        );
    }
}
