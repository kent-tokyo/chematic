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

/// Parse a CDXML document and return the first molecular fragment.
///
/// Convenience wrapper around [`parse_cdxml_all`].
pub fn parse_cdxml(input: &str) -> Result<(Molecule, Vec<(f64, f64)>), CdxmlError> {
    let mut all = parse_cdxml_all(input)?;
    if all.is_empty() {
        return Ok((MoleculeBuilder::new().build(), vec![]));
    }
    Ok(all.remove(0))
}

/// Parse all molecular fragments from a CDXML document.
///
/// Each `<fragment>` element in the document is parsed as a separate
/// molecule.  Returns a `Vec` of `(Molecule, 2D-coords)` pairs in document
/// order.  Coordinates are in CDXML point units (1/72 inch).
///
/// # Stereochemistry
///
/// Wedge bonds are derived from the `Display` attribute of `<b>` elements:
/// `"WedgeBegin"` / `"WedgedHashBegin"` → [`BondOrder::Up`];
/// `"Hash"` / `"Dash"` / `"WedgeEnd"` / `"WedgedHashEnd"` → [`BondOrder::Down`].
pub fn parse_cdxml_all(input: &str) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError> {
    // Per-fragment atom and bond accumulators.
    let mut atom_ids:   Vec<String>    = Vec::new();
    let mut atom_elems: Vec<Element>   = Vec::new();
    let mut atom_charges: Vec<i8>      = Vec::new();
    let mut atom_isotopes: Vec<Option<u16>> = Vec::new();
    let mut atom_h:     Vec<Option<u8>> = Vec::new();
    let mut atom_xs:    Vec<f64>       = Vec::new();
    let mut atom_ys:    Vec<f64>       = Vec::new();

    let mut bond_bs:    Vec<String>    = Vec::new();
    let mut bond_es:    Vec<String>    = Vec::new();
    let mut bond_ords:  Vec<BondOrder> = Vec::new();

    let mut results: Vec<(Molecule, Vec<(f64, f64)>)> = Vec::new();
    let mut in_fragment = false;

    let flush = |atom_ids: &mut Vec<String>, atom_elems: &mut Vec<Element>,
                 atom_charges: &mut Vec<i8>, atom_isotopes: &mut Vec<Option<u16>>,
                 atom_h: &mut Vec<Option<u8>>, atom_xs: &mut Vec<f64>, atom_ys: &mut Vec<f64>,
                 bond_bs: &mut Vec<String>, bond_es: &mut Vec<String>, bond_ords: &mut Vec<BondOrder>,
                 results: &mut Vec<(Molecule, Vec<(f64, f64)>)>|
        -> Result<(), CdxmlError>
    {
        if atom_ids.is_empty() && bond_bs.is_empty() { return Ok(()); }

        let mut id_to_pos: HashMap<&str, usize> = HashMap::new();
        for (i, id) in atom_ids.iter().enumerate() { id_to_pos.insert(id.as_str(), i); }

        let mut builder = MoleculeBuilder::new();
        let mut idx_map: HashMap<usize, AtomIdx> = HashMap::new();
        let mut coords: Vec<(f64, f64)> = Vec::new();

        for (i, _) in atom_ids.iter().enumerate() {
            let mut a = Atom::new(atom_elems[i]);
            a.charge = atom_charges[i];
            a.isotope = atom_isotopes[i];
            a.hydrogen_count = atom_h[i];
            let new_idx = builder.add_atom(a);
            idx_map.insert(i, new_idx);
            coords.push((atom_xs[i], atom_ys[i]));
        }

        for k in 0..bond_bs.len() {
            let pos_b = *id_to_pos.get(bond_bs[k].as_str())
                .ok_or_else(|| CdxmlError::UnknownAtomRef(bond_bs[k].clone()))?;
            let pos_e = *id_to_pos.get(bond_es[k].as_str())
                .ok_or_else(|| CdxmlError::UnknownAtomRef(bond_es[k].clone()))?;
            let a1 = idx_map[&pos_b];
            let a2 = idx_map[&pos_e];
            builder.add_bond(a1, a2, bond_ords[k])
                .map_err(|_| CdxmlError::UnknownAtomRef(format!("{} {}", bond_bs[k], bond_es[k])))?;
        }

        results.push((builder.build(), coords));
        atom_ids.clear(); atom_elems.clear(); atom_charges.clear(); atom_isotopes.clear();
        atom_h.clear(); atom_xs.clear(); atom_ys.clear();
        bond_bs.clear(); bond_es.clear(); bond_ords.clear();
        Ok(())
    };

    for raw_line in input.lines() {
        let line = raw_line.trim();
        if line.is_empty() { continue; }

        if line.starts_with("<fragment") {
            in_fragment = true;
            atom_ids.clear(); atom_elems.clear(); atom_charges.clear(); atom_isotopes.clear();
            atom_h.clear(); atom_xs.clear(); atom_ys.clear();
            bond_bs.clear(); bond_es.clear(); bond_ords.clear();
            continue;
        }

        if line.starts_with("</fragment>") {
            flush(&mut atom_ids, &mut atom_elems, &mut atom_charges, &mut atom_isotopes,
                  &mut atom_h, &mut atom_xs, &mut atom_ys,
                  &mut bond_bs, &mut bond_es, &mut bond_ords, &mut results)?;
            in_fragment = false;
            continue;
        }

        if is_n_tag(line) {
            let attrs = parse_xml_attrs(line);
            let id = match attrs.get("id") { Some(s) => s.clone(), None => continue };
            let element_num: u32 = attrs.get("Element")
                .and_then(|s| s.trim().parse().ok()).unwrap_or(6);
            let element = Element::from_atomic_number(element_num as u8)
                .ok_or(CdxmlError::UnknownAtomicNumber(element_num))?;
            let charge = attrs.get("Charge").and_then(|s| s.trim().parse().ok()).unwrap_or(0);
            let isotope = attrs.get("Isotope")
                .and_then(|s| s.trim().parse::<u16>().ok()).filter(|&v| v > 0);
            let hcount = attrs.get("NumHydrogens").and_then(|s| s.trim().parse().ok());
            let (x, y) = if let Some(p) = attrs.get("p") {
                let parts: Vec<&str> = p.split_whitespace().collect();
                if parts.len() < 2 { return Err(CdxmlError::InvalidCoords(p.clone())); }
                (parts[0].parse().unwrap_or(0.0), parts[1].parse().unwrap_or(0.0))
            } else { (0.0, 0.0) };
            atom_ids.push(id); atom_elems.push(element); atom_charges.push(charge);
            atom_isotopes.push(isotope); atom_h.push(hcount); atom_xs.push(x); atom_ys.push(y);
            continue;
        }

        if is_b_tag(line) {
            let attrs = parse_xml_attrs(line);
            let b = attrs.get("B").cloned().ok_or(CdxmlError::MissingBondEndpoint)?;
            let e = attrs.get("E").cloned().ok_or(CdxmlError::MissingBondEndpoint)?;
            let base: BondOrder = match attrs.get("Order").map(String::as_str) {
                Some("2") => BondOrder::Double,
                Some("3") => BondOrder::Triple,
                _         => BondOrder::Single,
            };
            let order = if base == BondOrder::Single {
                match attrs.get("Display").map(String::as_str) {
                    Some("WedgeBegin") | Some("WedgedHashBegin") => BondOrder::Up,
                    Some("Hash") | Some("Dash") | Some("WedgeEnd") | Some("WedgedHashEnd")
                        => BondOrder::Down,
                    _ => BondOrder::Single,
                }
            } else { base };
            bond_bs.push(b); bond_es.push(e); bond_ords.push(order);
        }
    }

    // Handle documents without explicit </fragment> closing tags.
    if !atom_ids.is_empty() || !bond_bs.is_empty() {
        flush(&mut atom_ids, &mut atom_elems, &mut atom_charges, &mut atom_isotopes,
              &mut atom_h, &mut atom_xs, &mut atom_ys,
              &mut bond_bs, &mut bond_es, &mut bond_ords, &mut results)?;
    }

    Ok(results)
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
