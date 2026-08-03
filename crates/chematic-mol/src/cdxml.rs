//! ChemDraw XML (CDXML) parser and writer.
//!
//! CDXML is a proprietary XML format produced by ChemDraw (PerkinElmer /
//! Revvity).  This module handles the minimal subset needed to read and
//! write molecular structure through a CDXML document.
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
//! - [`write_cdxml`] targets **self-round-trip** correctness (this parser
//!   can read what this writer produces), not full ChemDraw-application
//!   compatibility — CDXML is a proprietary format and writing files
//!   ChemDraw itself will accept requires undocumented attributes.
//! - Wedge/hash bond stereo (tetrahedral) is derived from `Display` attribute.
//! - E/Z double-bond stereo is derived from 2D coordinates via `assign_ez_from_2d`.
//! - Presentation-only nodes (text boxes, arrows, etc.) are silently skipped.

use std::collections::HashMap;

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_perception::assign_ez_from_2d;

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
    /// The document contains more atoms than the parser's safety limit.
    TooManyAtoms(usize),
}

impl std::fmt::Display for CdxmlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CdxmlError::UnknownAtomicNumber(n) => write!(f, "unknown atomic number: {n}"),
            CdxmlError::UnknownAtomRef(s) => write!(f, "unknown atom ref: {s}"),
            CdxmlError::MissingBondEndpoint => write!(f, "bond missing B or E attribute"),
            CdxmlError::InvalidCoords(s) => write!(f, "invalid p coords: {s}"),
            CdxmlError::TooManyAtoms(n) => write!(f, "CDXML document exceeds atom limit ({n})"),
        }
    }
}

impl std::error::Error for CdxmlError {}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a CDXML document and return the first molecular fragment.
///
/// Convenience wrapper around [`parse_cdxml_all`].
/// Parse a ChemDraw CDXML string into a molecule and 2D coordinates.
///
/// **Coordinate system:** The returned `coords` use **ChemDraw Y-down convention**
/// (Y increases downward, matching screen/SVG pixel space). No Y-axis conversion is required
/// for SVG rendering; coordinates can be used directly in [`crate::svg::render_svg`] or similar.
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
/// `"WedgeBegin"` / `"WedgedHashBegin"` / `"Bold"` → [`BondOrder::Up`];
/// `"Hash"` / `"Dash"` / `"WedgeEnd"` / `"WedgedHashEnd"` → [`BondOrder::Down`].
/// `"Bold"` is ChemDraw's simplified, non-directional "coming toward viewer"
/// convention (a plain thick line, used by some chemists in place of a real
/// wedge) -- it carries no begin/end distinction of its own, so it's
/// interpreted the same bond-direction-from-B-to-E way `"Hash"`/`"Dash"`
/// already are, not via a new direction-inference mechanism.
/// Accumulator for a single CDXML `<fragment>` being parsed.
#[derive(Default)]
struct FragAccum {
    atom_ids: Vec<String>,
    atom_elems: Vec<Element>,
    atom_charges: Vec<i8>,
    atom_isotopes: Vec<Option<u16>>,
    atom_h: Vec<Option<u8>>,
    atom_xs: Vec<f64>,
    atom_ys: Vec<f64>,
    bond_bs: Vec<String>,
    bond_es: Vec<String>,
    bond_ords: Vec<BondOrder>,
}

impl FragAccum {
    fn is_empty(&self) -> bool {
        self.atom_ids.is_empty() && self.bond_bs.is_empty()
    }

    fn flush(&mut self, results: &mut Vec<(Molecule, Vec<(f64, f64)>)>) -> Result<(), CdxmlError> {
        if self.is_empty() {
            return Ok(());
        }

        let mut id_to_pos: HashMap<&str, usize> = HashMap::new();
        for (i, id) in self.atom_ids.iter().enumerate() {
            id_to_pos.insert(id.as_str(), i);
        }

        let mut builder = MoleculeBuilder::new();
        let mut idx_map: HashMap<usize, AtomIdx> = HashMap::new();
        let mut coords: Vec<(f64, f64)> = Vec::new();

        for i in 0..self.atom_ids.len() {
            let mut a = Atom::new(self.atom_elems[i]);
            a.charge = self.atom_charges[i];
            a.isotope = self.atom_isotopes[i];
            a.hydrogen_count = self.atom_h[i];
            let new_idx = builder.add_atom(a);
            idx_map.insert(i, new_idx);
            coords.push((self.atom_xs[i], self.atom_ys[i]));
        }

        for k in 0..self.bond_bs.len() {
            let pos_b = *id_to_pos
                .get(self.bond_bs[k].as_str())
                .ok_or_else(|| CdxmlError::UnknownAtomRef(self.bond_bs[k].clone()))?;
            let pos_e = *id_to_pos
                .get(self.bond_es[k].as_str())
                .ok_or_else(|| CdxmlError::UnknownAtomRef(self.bond_es[k].clone()))?;
            let a1 = idx_map[&pos_b];
            let a2 = idx_map[&pos_e];
            builder.add_bond(a1, a2, self.bond_ords[k]).map_err(|_| {
                CdxmlError::UnknownAtomRef(format!("{} {}", self.bond_bs[k], self.bond_es[k]))
            })?;
        }

        let mut mol = builder.build();
        // Derive E/Z stereo from 2D atom positions (RDKit issue #9356: CDXML loses E/Z).
        assign_ez_from_2d(&mut mol, &coords);
        results.push((mol, coords));
        *self = FragAccum::default();
        Ok(())
    }
}

/// Maximum atoms allowed in a single CDXML fragment (DoS guard).
pub const CDXML_MAX_ATOMS: usize = 10_000;

#[allow(clippy::type_complexity)]
pub fn parse_cdxml_all(input: &str) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError> {
    let mut acc = FragAccum::default();
    let mut results: Vec<(Molecule, Vec<(f64, f64)>)> = Vec::new();

    for raw_line in input.lines() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with("<fragment") {
            acc = FragAccum::default();
            continue;
        }

        if line.starts_with("</fragment>") {
            acc.flush(&mut results)?;
            continue;
        }

        if is_n_tag(line) {
            let attrs = parse_xml_attrs(line);
            let id = match attrs.get("id") {
                Some(s) => s.clone(),
                None => continue,
            };
            let element_num: u32 = attrs
                .get("Element")
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(6);
            if element_num > 255 {
                return Err(CdxmlError::UnknownAtomicNumber(element_num));
            }
            let element = Element::from_atomic_number(element_num as u8)
                .ok_or(CdxmlError::UnknownAtomicNumber(element_num))?;
            let charge = attrs
                .get("Charge")
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0);
            let isotope = attrs
                .get("Isotope")
                .and_then(|s| s.trim().parse::<u16>().ok())
                .filter(|&v| v > 0);
            let hcount = attrs
                .get("NumHydrogens")
                .and_then(|s| s.trim().parse().ok());
            let (x, y) = if let Some(p) = attrs.get("p") {
                let parts: Vec<&str> = p.split_whitespace().collect();
                if parts.len() < 2 {
                    return Err(CdxmlError::InvalidCoords(p.clone()));
                }
                (
                    parts[0].parse().unwrap_or(0.0),
                    parts[1].parse().unwrap_or(0.0),
                )
            } else {
                (0.0, 0.0)
            };
            if acc.atom_ids.len() >= CDXML_MAX_ATOMS {
                return Err(CdxmlError::TooManyAtoms(CDXML_MAX_ATOMS));
            }
            acc.atom_ids.push(id);
            acc.atom_elems.push(element);
            acc.atom_charges.push(charge);
            acc.atom_isotopes.push(isotope);
            acc.atom_h.push(hcount);
            acc.atom_xs.push(x);
            acc.atom_ys.push(y);
            continue;
        }

        if is_b_tag(line) {
            let attrs = parse_xml_attrs(line);
            let b = attrs
                .get("B")
                .cloned()
                .ok_or(CdxmlError::MissingBondEndpoint)?;
            let e = attrs
                .get("E")
                .cloned()
                .ok_or(CdxmlError::MissingBondEndpoint)?;
            let base: BondOrder = match attrs.get("Order").map(String::as_str) {
                Some("2") => BondOrder::Double,
                Some("3") => BondOrder::Triple,
                // "1.5" is used by some CDXML writers (e.g. OpenBabel) for
                // aromatic bonds.  Store as Aromatic so that aromaticity
                // perception is not required to recover the correct bond type.
                Some("1.5") => BondOrder::Aromatic,
                _ => BondOrder::Single,
            };
            let order = if base == BondOrder::Single {
                match attrs.get("Display").map(String::as_str) {
                    Some("WedgeBegin") | Some("WedgedHashBegin") | Some("Bold") => BondOrder::Up,
                    Some("Hash") | Some("Dash") | Some("WedgeEnd") | Some("WedgedHashEnd") => {
                        BondOrder::Down
                    }
                    _ => BondOrder::Single,
                }
            } else {
                base
            };
            acc.bond_bs.push(b);
            acc.bond_es.push(e);
            acc.bond_ords.push(order);
        }
    }

    // Handle documents without explicit </fragment> closing tags.
    acc.flush(&mut results)?;

    Ok(results)
}

/// True if `line` starts a CDXML `<n` atom node tag.
fn is_n_tag(line: &str) -> bool {
    (line.starts_with("<n ") || line.starts_with("<n\t") || line == "<n>")
        && !line.starts_with("<node") // avoid accidental match on <node>
}

/// True if `line` starts a CDXML `<b` bond tag.
fn is_b_tag(line: &str) -> bool {
    line.starts_with("<b ")
        || line.starts_with("<b\t")
        || line == "<b>"
        || line.starts_with("<b/>")
        || line.starts_with("<b>")
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Serialise `mol` to a minimal CDXML document that [`parse_cdxml`] can read
/// back. See the module docs for the self-round-trip scope of this writer.
///
/// `coords[i]` is the `(x, y)` position for atom `i`, in the same
/// **ChemDraw Y-down convention** `parse_cdxml` produces. Atoms beyond
/// `coords.len()` receive `(0.0, 0.0)`.
pub fn write_cdxml(mol: &Molecule, coords: &[(f64, f64)]) -> String {
    let mut out =
        String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<CDXML>\n<page>\n<fragment>\n");

    for (i, (idx, atom)) in mol.atoms().enumerate() {
        let id = idx.0 + 1;
        let (x, y) = coords.get(i).copied().unwrap_or((0.0, 0.0));
        let mut parts = vec![
            format!("id=\"{id}\""),
            format!("p=\"{x} {y}\""),
            format!("Element=\"{}\"", atom.element.atomic_number()),
        ];
        if let Some(h) = atom.hydrogen_count {
            parts.push(format!("NumHydrogens=\"{h}\""));
        }
        if atom.charge != 0 {
            parts.push(format!("Charge=\"{}\"", atom.charge));
        }
        if let Some(iso) = atom.isotope {
            parts.push(format!("Isotope=\"{iso}\""));
        }
        out.push_str(&format!("<n {}/>\n", parts.join(" ")));
    }

    for (_, bond) in mol.bonds() {
        let b = bond.atom1.0 + 1;
        let e = bond.atom2.0 + 1;
        // Wedge/hash stereo is Single order + Display attribute (mirrors the
        // parser's decoding at is_b_tag handling above); other orders map
        // straight to the Order attribute and ignore Display.
        let (order_attr, display_attr) = match bond.order {
            BondOrder::Double => ("2", None),
            BondOrder::Triple => ("3", None),
            BondOrder::Aromatic => ("1.5", None),
            BondOrder::Up => ("1", Some("WedgeBegin")),
            BondOrder::Down => ("1", Some("Hash")),
            _ => ("1", None),
        };
        match display_attr {
            Some(d) => out.push_str(&format!(
                "<b B=\"{b}\" E=\"{e}\" Order=\"{order_attr}\" Display=\"{d}\"/>\n"
            )),
            None => out.push_str(&format!(
                "<b B=\"{b}\" E=\"{e}\" Order=\"{order_attr}\"/>\n"
            )),
        }
    }

    out.push_str("</fragment>\n</page>\n</CDXML>\n");
    out
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
        assert!(
            (coords[0].0 - 10.0).abs() < 0.01,
            "first atom x=10.0: {:?}",
            coords[0]
        );
        assert!(
            (coords[0].1 - 20.0).abs() < 0.01,
            "first atom y=20.0: {:?}",
            coords[0]
        );
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

    #[test]
    fn parse_cdxml_element_above_255_is_rejected() {
        // Element=300 would silently truncate to 44 (Ru) via u32 as u8.
        // Must be caught as UnknownAtomicNumber before the cast.
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="300" p="0 0"/>
</fragment></CDXML>"#;
        let result = parse_cdxml(cdxml);
        assert!(
            matches!(result, Err(CdxmlError::UnknownAtomicNumber(300))),
            "Element=300 must return UnknownAtomicNumber(300)"
        );
    }

    // B4: multi-fragment CDXML document tests

    const TWO_FRAGMENT_CDXML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<CDXML>
<fragment>
<n id="1" p="10.0 20.0" Element="6"/>
<n id="2" p="25.0 20.0" Element="8"/>
<b B="1" E="2" Order="1"/>
</fragment>
<fragment>
<n id="3" p="60.0 20.0" Element="7"/>
<n id="4" p="75.0 20.0" Element="6"/>
<n id="5" p="90.0 20.0" Element="6"/>
<b B="3" E="4" Order="1"/>
<b B="4" E="5" Order="1"/>
</fragment>
</CDXML>"#;

    #[test]
    fn parse_cdxml_all_two_fragments_count() {
        let mols = parse_cdxml_all(TWO_FRAGMENT_CDXML).unwrap();
        assert_eq!(mols.len(), 2, "two <fragment> elements → two molecules");
    }

    #[test]
    fn parse_cdxml_all_first_fragment_co() {
        let mols = parse_cdxml_all(TWO_FRAGMENT_CDXML).unwrap();
        let (mol, _) = &mols[0];
        assert_eq!(mol.atom_count(), 2, "first fragment: C + O");
        assert_eq!(mol.bond_count(), 1);
    }

    #[test]
    fn parse_cdxml_all_second_fragment_ncc() {
        let mols = parse_cdxml_all(TWO_FRAGMENT_CDXML).unwrap();
        let (mol, _) = &mols[1];
        assert_eq!(mol.atom_count(), 3, "second fragment: N + C + C");
        assert_eq!(mol.bond_count(), 2);
    }

    #[test]
    fn parse_cdxml_all_coords_independent() {
        let mols = parse_cdxml_all(TWO_FRAGMENT_CDXML).unwrap();
        let (_, coords0) = &mols[0];
        let (_, coords1) = &mols[1];
        // First fragment atom 0 at x=10
        assert!((coords0[0].0 - 10.0).abs() < 0.01);
        // Second fragment atom 0 at x=60
        assert!((coords1[0].0 - 60.0).abs() < 0.01);
    }

    #[test]
    fn parse_cdxml_empty_doc_returns_empty_vec() {
        let cdxml = r#"<?xml version="1.0"?><CDXML></CDXML>"#;
        let mols = parse_cdxml_all(cdxml).unwrap();
        assert!(mols.is_empty(), "empty CDXML → empty Vec");
    }

    // -----------------------------------------------------------------------
    // Order="1.5" aromatic bond (OpenBabel / some CDXML writers)
    // -----------------------------------------------------------------------

    #[test]
    fn parse_cdxml_aromatic_bond_order_1_5() {
        // Some CDXML producers (e.g. tools derived from OpenBabel) write
        // aromatic bonds as Order="1.5".  These must be stored as
        // BondOrder::Aromatic rather than falling through to Single.
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1.5"/>
</fragment></CDXML>"#;
        let (mol, _) = parse_cdxml(cdxml).unwrap();
        let bond = mol.bond(chematic_core::BondIdx(0));
        assert_eq!(bond.order, BondOrder::Aromatic, "Order=1.5 → Aromatic");
    }

    #[test]
    fn parse_cdxml_benzene_all_aromatic_bonds() {
        // Benzene written with all Order="1.5" bonds.
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<n id="3" Element="6" p="20 0"/>
<n id="4" Element="6" p="30 0"/>
<n id="5" Element="6" p="20 10"/>
<n id="6" Element="6" p="10 10"/>
<b B="1" E="2" Order="1.5"/>
<b B="2" E="3" Order="1.5"/>
<b B="3" E="4" Order="1.5"/>
<b B="4" E="5" Order="1.5"/>
<b B="5" E="6" Order="1.5"/>
<b B="6" E="1" Order="1.5"/>
</fragment></CDXML>"#;
        let (mol, _) = parse_cdxml(cdxml).unwrap();
        assert_eq!(mol.atom_count(), 6, "benzene has 6 atoms");
        assert_eq!(mol.bond_count(), 6, "benzene has 6 bonds");
        let all_aromatic = mol.bonds().all(|(_, b)| b.order == BondOrder::Aromatic);
        assert!(all_aromatic, "all bonds must be Aromatic");
    }

    // -----------------------------------------------------------------------
    // write_cdxml round-trip
    // -----------------------------------------------------------------------

    #[test]
    fn write_cdxml_roundtrip_ethanol() {
        let (mol, coords) = parse_cdxml(ETHANOL_CDXML).unwrap();
        let written = write_cdxml(&mol, &coords);
        let (mol2, coords2) = parse_cdxml(&written).unwrap();

        assert_eq!(mol2.atom_count(), mol.atom_count());
        assert_eq!(mol2.bond_count(), mol.bond_count());
        let elems: Vec<&str> = mol.atoms().map(|(_, a)| a.element.symbol()).collect();
        let elems2: Vec<&str> = mol2.atoms().map(|(_, a)| a.element.symbol()).collect();
        assert_eq!(elems, elems2);
        for (c1, c2) in coords.iter().zip(coords2.iter()) {
            assert!((c1.0 - c2.0).abs() < 0.01);
            assert!((c1.1 - c2.1).abs() < 0.01);
        }
    }

    #[test]
    fn write_cdxml_roundtrip_double_bond() {
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="8" p="10 0"/>
<b B="1" E="2" Order="2"/>
</fragment></CDXML>"#;
        let (mol, coords) = parse_cdxml(cdxml).unwrap();
        let written = write_cdxml(&mol, &coords);
        let (mol2, _) = parse_cdxml(&written).unwrap();
        let bond2 = mol2.bond(chematic_core::BondIdx(0));
        assert_eq!(bond2.order, BondOrder::Double);
    }

    #[test]
    fn write_cdxml_roundtrip_wedge_bond() {
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="WedgeBegin"/>
</fragment></CDXML>"#;
        let (mol, coords) = parse_cdxml(cdxml).unwrap();
        let written = write_cdxml(&mol, &coords);
        assert!(written.contains("Display=\"WedgeBegin\""));
        let (mol2, _) = parse_cdxml(&written).unwrap();
        let bond2 = mol2.bond(chematic_core::BondIdx(0));
        assert_eq!(bond2.order, BondOrder::Up);
    }

    /// Issue found while surveying RDKit's open issues (analogous to RDKit
    /// #9359, "CDXML reading doesn't use Bold or undirectional hash bonds
    /// for stereochemistry"): ChemDraw's simplified non-directional "Bold"
    /// bond display (a plain thick line, used by some chemists in place of
    /// a real wedge) must be interpreted as stereo, same as this reader
    /// already does for bare "Hash".
    #[test]
    fn parse_cdxml_bold_bond_is_wedge_up() {
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="Bold"/>
</fragment></CDXML>"#;
        let (mol, _) = parse_cdxml(cdxml).unwrap();
        let bond = mol.bond(chematic_core::BondIdx(0));
        assert_eq!(
            bond.order,
            BondOrder::Up,
            "Display=\"Bold\" must be read as a wedge-up stereo bond, not silently \
             dropped to a plain single bond"
        );
    }

    /// "Bold" and "WedgeBegin" both encode the same "coming toward viewer"
    /// stereo relative to a bond's B->E atom order -- for identical
    /// connectivity/geometry, both display attributes must parse to the
    /// same BondOrder.
    #[test]
    fn parse_cdxml_bold_bond_matches_wedge_begin() {
        let bold = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="Bold"/>
</fragment></CDXML>"#;
        let wedge = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="WedgeBegin"/>
</fragment></CDXML>"#;
        let (mol_bold, _) = parse_cdxml(bold).unwrap();
        let (mol_wedge, _) = parse_cdxml(wedge).unwrap();
        assert_eq!(
            mol_bold.bond(chematic_core::BondIdx(0)).order,
            mol_wedge.bond(chematic_core::BondIdx(0)).order,
        );
    }

    /// Negative control: bare "Hash"/"Dash" (already-handled non-directional
    /// stereo) and a plain undecorated bond must be unaffected by adding
    /// the "Bold" arm -- confirms no accidental widening of the match.
    #[test]
    fn parse_cdxml_bold_addition_does_not_change_other_display_values() {
        for (display, expected) in [
            ("Hash", BondOrder::Down),
            ("Dash", BondOrder::Down),
            ("WedgeEnd", BondOrder::Down),
            ("WedgedHashEnd", BondOrder::Down),
            ("WedgedHashBegin", BondOrder::Up),
        ] {
            let cdxml = format!(
                r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="{display}"/>
</fragment></CDXML>"#
            );
            let (mol, _) = parse_cdxml(&cdxml).unwrap();
            assert_eq!(
                mol.bond(chematic_core::BondIdx(0)).order,
                expected,
                "Display=\"{display}\" regressed"
            );
        }
        let plain = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1"/>
</fragment></CDXML>"#;
        let (mol, _) = parse_cdxml(plain).unwrap();
        assert_eq!(
            mol.bond(chematic_core::BondIdx(0)).order,
            BondOrder::Single,
            "a bond with no Display attribute at all must stay plain Single"
        );
    }

    #[test]
    fn write_cdxml_roundtrip_charge_isotope_hcount() {
        let cdxml = r#"<CDXML><fragment>
<n id="1" Element="7" Charge="1" Isotope="15" NumHydrogens="2" p="0 0"/>
</fragment></CDXML>"#;
        let (mol, coords) = parse_cdxml(cdxml).unwrap();
        let written = write_cdxml(&mol, &coords);
        let (mol2, _) = parse_cdxml(&written).unwrap();
        let atom2 = mol2.atom(chematic_core::AtomIdx(0));
        assert_eq!(atom2.charge, 1);
        assert_eq!(atom2.isotope, Some(15));
        assert_eq!(atom2.hydrogen_count, Some(2));
    }
}
