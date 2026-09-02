//! CML (Chemical Markup Language) parser and writer.
//!
//! CML is an open XML-based format for chemical structures.  This module
//! implements a minimal subset: `<molecule>`, `<atomArray>`, `<bondArray>`,
//! `<atom>`, and `<bond>` elements with the attributes needed to represent
//! 2D chemical structures.
//!
//! No external XML library is used; attributes are extracted with a simple
//! key="value" scanner.
//!
//! # Supported attributes
//!
//! `<atom>`: `id`, `elementType`, `x2`, `y2`, `formalCharge`, `isotope`,
//!           `hydrogenCount`, `isotopeNumber`
//!
//! `<bond>`: `atomRefs2` (two atom ids separated by a space), `order`
//!           ("1" / "2" / "3" / "A" / "S" / "D" / "T")
//!
//! # Limitations
//!
//! - 3D coordinates (x3/y3/z3) are parsed but not returned (only 2D x2/y2).
//! - Aromatic bonds ("A") are stored as alternating single/double bonds (no
//!   special aromatic bond order in the core model).
//! - Stereochemistry attributes are ignored.
//! - The CML namespace declaration is accepted but not validated.

use std::collections::HashMap;

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

/// Resource limits for the lightweight CML parser.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CmlParseLimits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_atoms: usize,
    pub max_bonds: usize,
    pub max_xml_elements: usize,
}

impl Default for CmlParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_line_bytes: 1024 * 1024,
            max_lines: 1_000_000,
            max_atoms: 1_000_000,
            max_bonds: 2_000_000,
            max_xml_elements: 2_000_000,
        }
    }
}

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Error returned when parsing a CML document fails.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CmlError {
    /// An atom referenced `elementType` that is not a known element symbol.
    UnknownElement(String),
    /// A `<bond>` element referenced an atom id that was not defined.
    UnknownAtomRef(String),
    /// A `<bond>` element's `atomRefs2` attribute was not two space-separated ids.
    InvalidAtomRefs2(String),
    /// A `<bond>` element's `order` attribute could not be interpreted.
    InvalidBondOrder(String),
    /// A coordinate attribute parsed as NaN or infinite.
    NonFiniteCoordinate(String),
    /// The input exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl std::fmt::Display for CmlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CmlError::UnknownElement(s) => write!(f, "unknown element symbol: {s}"),
            CmlError::UnknownAtomRef(s) => write!(f, "unknown atom ref: {s}"),
            CmlError::InvalidAtomRefs2(s) => write!(f, "invalid atomRefs2: {s}"),
            CmlError::InvalidBondOrder(s) => write!(f, "invalid bond order: {s}"),
            CmlError::NonFiniteCoordinate(s) => write!(f, "non-finite coordinate: {s}"),
            CmlError::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "CML: {resource} has size {actual}, exceeding limit {limit}"
            ),
        }
    }
}

impl std::error::Error for CmlError {}

// ---------------------------------------------------------------------------
// XML attribute scanner (no external dependency)
// ---------------------------------------------------------------------------

/// Extract all `key="value"` pairs from a single XML tag line.
#[doc(hidden)]
///
/// Handles the most common cases:
/// - Double-quoted values
/// - Basic XML entities: `&amp;` `&lt;` `&gt;` `&quot;` `&apos;`
/// - Ignores self-closing slash and angle brackets
pub(crate) fn parse_xml_attrs(line: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let bytes = line.as_bytes();
    let len = bytes.len();
    let mut i = 0usize;

    while i < len {
        // Skip whitespace and punctuation between attributes.
        while i < len && (bytes[i].is_ascii_whitespace() || bytes[i] == b'/' || bytes[i] == b'>') {
            i += 1;
        }
        if i >= len {
            break;
        }

        // Read attribute key: stop at '=' or whitespace.
        let key_start = i;
        while i < len && bytes[i] != b'=' && !bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        let key = line.get(key_start..i).unwrap_or("").trim().to_string();
        if key.is_empty() {
            i += 1;
            continue;
        }

        // Expect '='
        while i < len && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= len || bytes[i] != b'=' {
            continue;
        }
        i += 1;

        // Expect opening '"' or '\''
        while i < len && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= len {
            break;
        }
        let quote = bytes[i];
        if quote != b'"' && quote != b'\'' {
            i += 1;
            continue;
        }
        i += 1;

        // Read until closing quote.
        let value_start = i;
        while i < len && bytes[i] != quote {
            i += 1;
        }
        let raw_value = line.get(value_start..i).unwrap_or("");
        let value = unescape_xml(raw_value);
        if i < len {
            i += 1;
        } // consume closing quote

        if !key.is_empty() {
            map.insert(key, value);
        }
    }
    map
}

/// Decode basic XML character entities.
fn unescape_xml(s: &str) -> String {
    if !s.contains('&') {
        return s.to_string();
    }
    s.replace("&amp;", "&")
        .replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
}

// ---------------------------------------------------------------------------
// CML bond order decoding
// ---------------------------------------------------------------------------

fn parse_bond_order(s: &str) -> Result<BondOrder, CmlError> {
    match s.trim() {
        "0" | "zero" => Ok(BondOrder::Zero),
        "1" | "S" | "single" => Ok(BondOrder::Single),
        "2" | "D" | "double" => Ok(BondOrder::Double),
        "3" | "T" | "triple" => Ok(BondOrder::Triple),
        // Aromatic bonds: store as single (aromaticity is perceived from topology)
        "A" | "aromatic" => Ok(BondOrder::Aromatic),
        "any" => Ok(BondOrder::QueryAny),
        "S/D" => Ok(BondOrder::QuerySingleOrDouble),
        "S/A" => Ok(BondOrder::QuerySingleOrAromatic),
        "D/A" => Ok(BondOrder::QueryDoubleOrAromatic),
        other => Err(CmlError::InvalidBondOrder(other.to_string())),
    }
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Intermediate atom data before building the `Molecule`.
struct CmlAtomData {
    element: Element,
    charge: i8,
    isotope: Option<u16>,
    hydrogen_count: Option<u8>,
    x2: f64,
    y2: f64,
    aromatic: bool,
}

/// Parse a CML document into a `(Molecule, 2D-coords)` pair.
///
/// The second element of the tuple contains one `(x, y)` entry per atom in
/// atom-insertion order.  If a molecule has no `x2`/`y2` attributes, the
/// coordinate pairs default to `(0.0, 0.0)`.
/// Parse a Chemical Markup Language (CML) string into a molecule and 2D coordinates.
///
/// **Coordinate system:** The returned `coords` use **chemical Y-up convention**
/// (Y increases upward, matching mathematical axes). This differs from SVG/screen Y-down.
/// When rendering with an SVG renderer such as `render_svg` (see `chematic-depict`),
/// callers must negate Y coordinates to match SVG pixel space.
///
/// # Example
/// ```text
/// // CML molecule with chemical Y-up coords
/// let (mol, coords) = parse_cml(cml_str)?;
/// // coords[i].1 increases upward (chemical convention)
///
/// // To render in SVG (Y-down):
/// let svg_coords: Vec<(f64, f64)> = coords.iter()
///     .map(|(x, y)| (x, -y))  // Negate Y
///     .collect();
/// ```
pub fn parse_cml(input: &str) -> Result<(Molecule, Vec<(f64, f64)>), CmlError> {
    parse_cml_with_limits(input, &CmlParseLimits::default())
}

/// Parse CML with explicit resource limits.
pub fn parse_cml_with_limits(
    input: &str,
    limits: &CmlParseLimits,
) -> Result<(Molecule, Vec<(f64, f64)>), CmlError> {
    if input.len() > limits.max_input_bytes {
        return Err(CmlError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let lines = input.lines().collect::<Vec<_>>();
    if lines.len() > limits.max_lines {
        return Err(CmlError::ResourceLimit {
            resource: "lines",
            actual: lines.len(),
            limit: limits.max_lines,
        });
    }
    if let Some(line_bytes) = lines.iter().map(|line| line.len()).max()
        && line_bytes > limits.max_line_bytes
    {
        return Err(CmlError::ResourceLimit {
            resource: "line bytes",
            actual: line_bytes,
            limit: limits.max_line_bytes,
        });
    }
    let element_count = lines
        .iter()
        .filter(|line| line.contains('<') && line.contains('>'))
        .count();
    if element_count > limits.max_xml_elements {
        return Err(CmlError::ResourceLimit {
            resource: "XML elements",
            actual: element_count,
            limit: limits.max_xml_elements,
        });
    }
    // --- Pass 1: collect atoms and bonds from the CML -----------------------

    let mut atom_data: Vec<CmlAtomData> = Vec::new();
    let mut atom_id_to_pos: HashMap<String, usize> = HashMap::new();

    struct BondData {
        a1: usize,
        a2: usize,
        order: BondOrder,
        aromatic: bool,
    }
    let mut bond_data: Vec<BondData> = Vec::new();

    for raw_line in lines {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }

        if is_element_tag(line, "atom") {
            if atom_data.len() >= limits.max_atoms {
                return Err(CmlError::ResourceLimit {
                    resource: "atom elements",
                    actual: atom_data.len() + 1,
                    limit: limits.max_atoms,
                });
            }
            let attrs = parse_xml_attrs(line);
            let id = attrs.get("id").cloned().unwrap_or_default();

            let element_sym = attrs.get("elementType").map(String::as_str).unwrap_or("C");
            let element = Element::from_symbol(element_sym)
                .ok_or_else(|| CmlError::UnknownElement(element_sym.to_string()))?;

            let charge = attrs
                .get("formalCharge")
                .and_then(|s| s.trim().parse::<i8>().ok())
                .unwrap_or(0);
            let isotope = attrs
                .get("isotope")
                .or_else(|| attrs.get("isotopeNumber"))
                .and_then(|s| s.trim().parse::<u16>().ok())
                .filter(|&v| v > 0);
            let hydrogen_count = attrs
                .get("hydrogenCount")
                .and_then(|s| s.trim().parse::<u8>().ok());
            let x2 = attrs
                .get("x2")
                .and_then(|s| s.trim().parse::<f64>().ok())
                .unwrap_or(0.0);
            let y2 = attrs
                .get("y2")
                .and_then(|s| s.trim().parse::<f64>().ok())
                .unwrap_or(0.0);
            if !x2.is_finite() {
                return Err(CmlError::NonFiniteCoordinate(x2.to_string()));
            }
            if !y2.is_finite() {
                return Err(CmlError::NonFiniteCoordinate(y2.to_string()));
            }

            let pos = atom_data.len();
            if !id.is_empty() {
                atom_id_to_pos.insert(id.clone(), pos);
            }
            atom_data.push(CmlAtomData {
                element,
                charge,
                isotope,
                hydrogen_count,
                x2,
                y2,
                aromatic: false,
            });
            continue;
        }

        if is_element_tag(line, "bond") {
            if bond_data.len() >= limits.max_bonds {
                return Err(CmlError::ResourceLimit {
                    resource: "bond elements",
                    actual: bond_data.len() + 1,
                    limit: limits.max_bonds,
                });
            }
            let attrs = parse_xml_attrs(line);
            let refs = attrs
                .get("atomRefs2")
                .ok_or_else(|| CmlError::InvalidAtomRefs2("missing atomRefs2".to_string()))?;
            let parts: Vec<&str> = refs.split_whitespace().collect();
            if parts.len() < 2 {
                return Err(CmlError::InvalidAtomRefs2(refs.clone()));
            }
            let a1 = *atom_id_to_pos
                .get(parts[0])
                .ok_or_else(|| CmlError::UnknownAtomRef(parts[0].to_string()))?;
            let a2 = *atom_id_to_pos
                .get(parts[1])
                .ok_or_else(|| CmlError::UnknownAtomRef(parts[1].to_string()))?;
            let order_s = attrs.get("order").map(String::as_str).unwrap_or("1");
            let is_aromatic = matches!(order_s.trim(), "A" | "aromatic");
            let order = parse_bond_order(order_s)?;
            bond_data.push(BondData {
                a1,
                a2,
                order,
                aromatic: is_aromatic,
            });
        }
    }

    // --- Pass 2: mark aromatic atoms, then build Molecule -------------------

    for bd in &bond_data {
        if bd.aromatic {
            if let Some(a) = atom_data.get_mut(bd.a1) {
                a.aromatic = true;
            }
            if let Some(a) = atom_data.get_mut(bd.a2) {
                a.aromatic = true;
            }
        }
    }

    let mut builder = MoleculeBuilder::new();
    let mut idx_by_pos: Vec<AtomIdx> = Vec::new();
    let mut coords: Vec<(f64, f64)> = Vec::new();

    for ad in &atom_data {
        let mut atom = Atom::new(ad.element);
        atom.charge = ad.charge;
        atom.isotope = ad.isotope;
        atom.hydrogen_count = ad.hydrogen_count;
        atom.aromatic = ad.aromatic;
        let idx = builder.add_atom(atom);
        idx_by_pos.push(idx);
        coords.push((ad.x2, ad.y2));
    }

    for bd in &bond_data {
        let a1 = idx_by_pos[bd.a1];
        let a2 = idx_by_pos[bd.a2];
        builder
            .add_bond(a1, a2, bd.order)
            .map_err(|_| CmlError::InvalidAtomRefs2(format!("{} {}", bd.a1, bd.a2)))?;
    }

    Ok((builder.build(), coords))
}

/// True if `line` starts the opening tag for the given element name
/// (case-insensitive on the element name, as some writers vary case).
fn is_element_tag(line: &str, name: &str) -> bool {
    let check = |prefix: &str| -> bool {
        line.starts_with(prefix) && {
            let rest = &line[prefix.len()..];
            rest.is_empty()
                || rest.starts_with(' ')
                || rest.starts_with('/')
                || rest.starts_with('>')
                || rest.starts_with('\t')
        }
    };
    check(&format!("<{name}")) || check(&format!("<{}", name.to_uppercase()))
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Serialise a molecule to a CML string.
///
/// `coords` — optional slice of `(x, y)` pairs in atom-insertion order.
/// If `None`, no coordinate attributes are written.
pub fn write_cml(mol: &Molecule, coords: Option<&[(f64, f64)]>) -> String {
    let mut out = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <molecule xmlns=\"http://www.xml-cml.org/schema\">\n  <atomArray>\n",
    );

    for (i, (idx, atom)) in mol.atoms().enumerate() {
        let sym = atom.element.symbol();
        let mut parts = vec![
            format!("id=\"a{}\"", idx.0 + 1),
            format!("elementType=\"{sym}\""),
        ];
        if let Some(cs) = coords
            && let Some(&(x, y)) = cs.get(i)
        {
            parts.push(format!("x2=\"{x:.4}\""));
            parts.push(format!("y2=\"{y:.4}\""));
        }
        if atom.charge != 0 {
            parts.push(format!("formalCharge=\"{}\"", atom.charge));
        }
        if let Some(iso) = atom.isotope {
            parts.push(format!("isotopeNumber=\"{iso}\""));
        }
        if let Some(h) = atom.hydrogen_count {
            parts.push(format!("hydrogenCount=\"{h}\""));
        }
        out.push_str(&format!("    <atom {}/>\n", parts.join(" ")));
    }

    out.push_str("  </atomArray>\n  <bondArray>\n");

    for (_, bond) in mol.bonds() {
        let order_str = match bond.order {
            BondOrder::Zero => "0",
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => "1",
            BondOrder::Double => "2",
            BondOrder::Triple => "3",
            BondOrder::Aromatic => "A",
            BondOrder::Quadruple => "4",
            BondOrder::QueryAny => "any",
            BondOrder::QuerySingleOrDouble => "S/D",
            BondOrder::QuerySingleOrAromatic => "S/A",
            BondOrder::QueryDoubleOrAromatic => "D/A",
        };
        let a1_id = bond.atom1.0 + 1;
        let a2_id = bond.atom2.0 + 1;
        out.push_str(&format!(
            "    <bond atomRefs2=\"a{a1_id} a{a2_id}\" order=\"{order_str}\"/>\n"
        ));
    }

    out.push_str("  </bondArray>\n</molecule>\n");
    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const ETHANOL_CML: &str = r#"<?xml version="1.0"?>
<molecule xmlns="http://www.xml-cml.org/schema">
  <atomArray>
    <atom id="a1" elementType="C" x2="0.0" y2="0.0"/>
    <atom id="a2" elementType="C" x2="1.5" y2="0.0"/>
    <atom id="a3" elementType="O" x2="3.0" y2="0.0"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a2 a3" order="1"/>
  </bondArray>
</molecule>"#;

    #[test]
    fn parse_cml_ethanol_atom_count() {
        let (mol, coords) = parse_cml(ETHANOL_CML).unwrap();
        assert_eq!(mol.atom_count(), 3, "ethanol has 3 heavy atoms");
        assert_eq!(mol.bond_count(), 2, "ethanol has 2 bonds");
        assert_eq!(coords.len(), 3, "one coord per atom");
    }

    #[test]
    fn parse_cml_ethanol_elements() {
        let (mol, _) = parse_cml(ETHANOL_CML).unwrap();
        let elems: Vec<&str> = mol.atoms().map(|(_, a)| a.element.symbol()).collect();
        assert!(elems.contains(&"C"), "should have C");
        assert!(elems.contains(&"O"), "should have O");
    }

    #[test]
    fn parse_cml_ethanol_coords() {
        let (_, coords) = parse_cml(ETHANOL_CML).unwrap();
        assert!((coords[0].0 - 0.0).abs() < 0.001, "a1 x=0.0");
        assert!((coords[1].0 - 1.5).abs() < 0.001, "a2 x=1.5");
        assert!((coords[2].0 - 3.0).abs() < 0.001, "a3 x=3.0");
    }

    #[test]
    fn parse_cml_formal_charge() {
        let cml = r#"<molecule>
  <atomArray>
    <atom id="a1" elementType="N" formalCharge="1"/>
  </atomArray>
  <bondArray/>
</molecule>"#;
        let (mol, _) = parse_cml(cml).unwrap();
        let atom = mol.atom(chematic_core::AtomIdx(0));
        assert_eq!(atom.charge, 1, "N+ should have charge=1");
    }

    #[test]
    fn parse_cml_unknown_element_returns_err() {
        let cml = r#"<molecule>
  <atomArray>
    <atom id="a1" elementType="Xx"/>
  </atomArray>
</molecule>"#;
        let result = parse_cml(cml);
        assert!(
            matches!(result, Err(CmlError::UnknownElement(_))),
            "unknown element should return Err"
        );
    }

    #[test]
    fn write_cml_roundtrip() {
        let (mol, coords) = parse_cml(ETHANOL_CML).unwrap();
        let written = write_cml(&mol, Some(&coords));
        let (mol2, coords2) = parse_cml(&written).unwrap();

        assert_eq!(mol.atom_count(), mol2.atom_count(), "atom count preserved");
        assert_eq!(mol.bond_count(), mol2.bond_count(), "bond count preserved");
        assert!(
            (coords[1].0 - coords2[1].0).abs() < 0.001,
            "x coord preserved"
        );
    }

    #[test]
    fn write_cml_no_coords() {
        use chematic_core::MoleculeBuilder;
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        b.add_atom(Atom::new(Element::from_symbol("O").unwrap()));
        let mol = b.build();
        let cml = write_cml(&mol, None);
        assert!(cml.contains("elementType=\"C\""), "C atom in output");
        assert!(cml.contains("elementType=\"O\""), "O atom in output");
        assert!(!cml.contains("x2="), "no x2 when no coords");
    }

    #[test]
    fn parse_cml_double_bond() {
        let cml = r#"<molecule>
  <atomArray>
    <atom id="a1" elementType="C"/>
    <atom id="a2" elementType="O"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="2"/>
  </bondArray>
</molecule>"#;
        let (mol, _) = parse_cml(cml).unwrap();
        let bond = mol.bond(chematic_core::BondIdx(0));
        assert_eq!(bond.order, BondOrder::Double, "order=2 should give Double");
    }

    #[test]
    fn bounded_parser_rejects_input_and_line_limits() {
        assert!(matches!(
            parse_cml_with_limits(
                ETHANOL_CML,
                &CmlParseLimits {
                    max_input_bytes: 8,
                    ..Default::default()
                }
            ),
            Err(CmlError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));
        let long_line = format!("{}\n", "x".repeat(32));
        assert!(matches!(
            parse_cml_with_limits(
                &long_line,
                &CmlParseLimits {
                    max_line_bytes: 16,
                    ..Default::default()
                }
            ),
            Err(CmlError::ResourceLimit {
                resource: "line bytes",
                ..
            })
        ));
    }

    #[test]
    fn bounded_parser_rejects_atom_and_bond_limits() {
        assert!(matches!(
            parse_cml_with_limits(
                ETHANOL_CML,
                &CmlParseLimits {
                    max_atoms: 2,
                    ..Default::default()
                }
            ),
            Err(CmlError::ResourceLimit {
                resource: "atom elements",
                ..
            })
        ));
        assert!(matches!(
            parse_cml_with_limits(
                ETHANOL_CML,
                &CmlParseLimits {
                    max_bonds: 1,
                    ..Default::default()
                }
            ),
            Err(CmlError::ResourceLimit {
                resource: "bond elements",
                ..
            })
        ));
    }
}
