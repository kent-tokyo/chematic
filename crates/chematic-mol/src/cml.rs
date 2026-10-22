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
//! `hydrogenCount` is the total attached H count, including explicit H atom
//! neighbors (Murray-Rust & Rzepa, DOI 10.1021/ci990052b, Table 4). The reader
//! subtracts those neighbors before storing the core's non-node H count.
//! When the attribute is absent, explicit H neighbors suppress further
//! valence-based H inference (Open Babel #2655). The writer always emits the
//! total count, including zero, so partial explicit-H graphs round-trip.
//!
//! # Limitations
//!
//! - Legacy readers ignore 3D coordinates. `parse_cml_with_3d` returns complete
//!   finite x3/y3/z3 triplets through an additive, bounded single-molecule profile.
//! - Aromatic bonds ("A") are stored as alternating single/double bonds (no
//!   special aromatic bond order in the core model).
//! - Stereochemistry attributes are ignored.
//! - The CML namespace declaration is accepted but not validated.

use std::collections::HashMap;

use chematic_core::{Atom, BondOrder, Element, Molecule, MoleculeBuilder, implicit_hcount};

mod geometry;
pub use geometry::{CmlWith3d, parse_cml_with_3d, parse_cml_with_3d_with_limits};

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
    /// A strict parse did not contain a `<molecule>` element.
    MissingMolecule,
    /// A strict parse contained a molecule but no atom elements.
    EmptyMolecule,
    /// A strict parse found mismatched or unclosed XML elements.
    MalformedXml(String),
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
            CmlError::MissingMolecule => write!(f, "CML: missing molecule element"),
            CmlError::EmptyMolecule => write!(f, "CML: molecule contains no atom elements"),
            CmlError::MalformedXml(s) => write!(f, "malformed CML XML: {s}"),
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

mod attributes;
pub(crate) use attributes::parse_xml_attrs;

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
    hydrogen_count: Option<usize>,
    x2: f64,
    y2: f64,
    aromatic: bool,
}

struct CmlBondData {
    a1: usize,
    a2: usize,
    order: BondOrder,
    aromatic: bool,
}

fn cml_lines<'a>(input: &'a str, limits: &CmlParseLimits) -> Result<Vec<&'a str>, CmlError> {
    if input.len() > limits.max_input_bytes {
        return Err(CmlError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let lines = input.lines().collect::<Vec<_>>();
    let constraints = [
        ("lines", lines.len(), limits.max_lines),
        (
            "line bytes",
            lines.iter().map(|line| line.len()).max().unwrap_or(0),
            limits.max_line_bytes,
        ),
        (
            "XML elements",
            lines
                .iter()
                .filter(|line| line.contains('<') && line.contains('>'))
                .count(),
            limits.max_xml_elements,
        ),
    ];
    if let Some((resource, actual, limit)) = constraints
        .into_iter()
        .find(|(_, actual, limit)| actual > limit)
    {
        return Err(CmlError::ResourceLimit {
            resource,
            actual,
            limit,
        });
    }
    Ok(lines)
}

fn cml_value<T: std::str::FromStr>(
    attrs: &HashMap<String, String>,
    key: &str,
) -> Result<Option<T>, CmlError> {
    attrs
        .get(key)
        .map(|value| {
            value
                .trim()
                .parse()
                .map_err(|_| CmlError::MalformedXml(format!("invalid CML {key}: {value}")))
        })
        .transpose()
}

fn cml_isotope(attrs: &HashMap<String, String>) -> Result<Option<u16>, CmlError> {
    let legacy = cml_value::<u16>(attrs, "isotope")?;
    let number = cml_value::<u16>(attrs, "isotopeNumber")?;
    if matches!((legacy, number), (Some(a), Some(b)) if a != b) {
        return Err(CmlError::MalformedXml(
            "conflicting CML isotope and isotopeNumber".into(),
        ));
    }
    Ok(number.or(legacy).filter(|&value| value > 0))
}

fn parse_cml_atom(attrs: &HashMap<String, String>) -> Result<CmlAtomData, CmlError> {
    let symbol = attrs.get("elementType").map(String::as_str).unwrap_or("C");
    let element =
        Element::from_symbol(symbol).ok_or_else(|| CmlError::UnknownElement(symbol.to_owned()))?;
    let x2 = cml_value::<f64>(attrs, "x2")?.unwrap_or(0.0);
    let y2 = cml_value::<f64>(attrs, "y2")?.unwrap_or(0.0);
    if !x2.is_finite() {
        return Err(CmlError::NonFiniteCoordinate(x2.to_string()));
    }
    if !y2.is_finite() {
        return Err(CmlError::NonFiniteCoordinate(y2.to_string()));
    }
    Ok(CmlAtomData {
        element,
        charge: cml_value(attrs, "formalCharge")?.unwrap_or(0),
        isotope: cml_isotope(attrs)?,
        hydrogen_count: cml_value(attrs, "hydrogenCount")?,
        x2,
        y2,
        aromatic: false,
    })
}

fn parse_cml_bond(
    attrs: &HashMap<String, String>,
    positions: &HashMap<String, usize>,
) -> Result<CmlBondData, CmlError> {
    let refs = attrs
        .get("atomRefs2")
        .ok_or_else(|| CmlError::InvalidAtomRefs2("missing atomRefs2".to_owned()))?;
    let atoms = refs.split_whitespace().collect::<Vec<_>>();
    if atoms.len() < 2 {
        return Err(CmlError::InvalidAtomRefs2(refs.clone()));
    }
    let position = |id: &str| {
        positions
            .get(id)
            .copied()
            .ok_or_else(|| CmlError::UnknownAtomRef(id.to_owned()))
    };
    let order_text = attrs.get("order").map(String::as_str).unwrap_or("1");
    Ok(CmlBondData {
        a1: position(atoms[0])?,
        a2: position(atoms[1])?,
        order: parse_bond_order(order_text)?,
        aromatic: matches!(order_text.trim(), "A" | "aromatic"),
    })
}

fn collect_cml_data(
    lines: Vec<&str>,
    limits: &CmlParseLimits,
) -> Result<(Vec<CmlAtomData>, Vec<CmlBondData>), CmlError> {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut positions = HashMap::new();
    for line in lines
        .into_iter()
        .map(str::trim)
        .filter(|line| !line.is_empty())
    {
        if is_element_tag(line, "atom") {
            if atoms.len() >= limits.max_atoms {
                return Err(CmlError::ResourceLimit {
                    resource: "atom elements",
                    actual: atoms.len() + 1,
                    limit: limits.max_atoms,
                });
            }
            let attrs = parse_xml_attrs(line);
            if let Some(id) = attrs.get("id").filter(|id| !id.is_empty())
                && positions.insert(id.clone(), atoms.len()).is_some()
            {
                return Err(CmlError::MalformedXml(format!(
                    "duplicate CML atom id: {id}"
                )));
            }
            atoms.push(parse_cml_atom(&attrs)?);
        } else if is_element_tag(line, "bond") {
            if bonds.len() >= limits.max_bonds {
                return Err(CmlError::ResourceLimit {
                    resource: "bond elements",
                    actual: bonds.len() + 1,
                    limit: limits.max_bonds,
                });
            }
            bonds.push(parse_cml_bond(&parse_xml_attrs(line), &positions)?);
        }
    }
    Ok((atoms, bonds))
}

fn build_cml(
    mut atoms: Vec<CmlAtomData>,
    bonds: Vec<CmlBondData>,
) -> Result<(Molecule, Vec<(f64, f64)>), CmlError> {
    // Include H isotopes: they are still Element::H graph nodes.
    let mut explicit_h = vec![0usize; atoms.len()];
    for bond in &bonds {
        if atoms[bond.a1].element == Element::H {
            explicit_h[bond.a2] += 1;
        }
        if atoms[bond.a2].element == Element::H {
            explicit_h[bond.a1] += 1;
        }
    }
    for bond in bonds.iter().filter(|bond| bond.aromatic) {
        atoms[bond.a1].aromatic = true;
        atoms[bond.a2].aromatic = true;
    }
    let mut builder = MoleculeBuilder::new();
    let mut indices = Vec::with_capacity(atoms.len());
    let mut coords = Vec::with_capacity(atoms.len());
    for (index, source) in atoms.into_iter().enumerate() {
        let mut atom = Atom::new(source.element);
        atom.charge = source.charge;
        atom.isotope = source.isotope;
        atom.hydrogen_count = match source.hydrogen_count {
            Some(total) => {
                let remainder = total.checked_sub(explicit_h[index]).ok_or_else(|| {
                    CmlError::MalformedXml(format!(
                        "invalid CML hydrogenCount for atom {index}: total {total} is smaller than {} explicit H neighbors",
                        explicit_h[index]
                    ))
                })?;
                Some(u8::try_from(remainder).map_err(|_| {
                    CmlError::MalformedXml(format!(
                        "invalid CML hydrogenCount for atom {index}: {remainder} non-node hydrogens exceeds 255"
                    ))
                })?)
            }
            None if explicit_h[index] > 0 => Some(0),
            None => None,
        };
        atom.aromatic = source.aromatic;
        indices.push(builder.add_atom(atom));
        coords.push((source.x2, source.y2));
    }
    for bond in bonds {
        builder
            .add_bond(indices[bond.a1], indices[bond.a2], bond.order)
            .map_err(|_| CmlError::InvalidAtomRefs2(format!("{} {}", bond.a1, bond.a2)))?;
    }
    Ok((builder.build(), coords))
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

/// Parse a structurally valid, non-empty CML molecule.
///
/// This is the strict counterpart to [`parse_cml`]. The default parser keeps
/// its historical lenient behavior for compatibility (including an empty
/// molecule when no `<molecule>`/`<atom>` elements are found). Strict callers
/// opt into XML nesting checks and a non-empty molecule boundary.
pub fn parse_cml_strict(input: &str) -> Result<(Molecule, Vec<(f64, f64)>), CmlError> {
    parse_cml_strict_with_limits(input, &CmlParseLimits::default())
}

/// Parse CML with explicit limits and the strict structural boundary.
pub fn parse_cml_strict_with_limits(
    input: &str,
    limits: &CmlParseLimits,
) -> Result<(Molecule, Vec<(f64, f64)>), CmlError> {
    validate_strict_cml_structure(input)?;
    let normalized = split_cml_tags(input);
    let parsed = parse_cml_with_limits(&normalized, limits)?;
    if parsed.0.atom_count() == 0 {
        return Err(CmlError::EmptyMolecule);
    }
    Ok(parsed)
}

/// Make the lightweight line-oriented reader usable with compact XML too.
///
/// The legacy parser intentionally only inspects one tag per input line. The
/// strict entry point can normalize tag boundaries without changing the
/// legacy parser's compatibility behavior. A quote-aware scan avoids splitting
/// an attribute value that happens to contain `>`.
fn split_cml_tags(input: &str) -> String {
    let mut output = String::with_capacity(input.len() + input.len() / 8);
    let mut in_tag = false;
    let mut quote = None;
    let mut cursor = 0usize;
    while cursor < input.len() {
        // Comments, CDATA, and processing instructions may contain `>` (and
        // comments/CDATA may contain `<`) as ordinary content. Keep each
        // construct intact so the line-oriented compatibility reader never
        // interprets markup inside it as a molecule element. The strict
        // validator has already checked that the terminator exists.
        if !in_tag && input[cursor..].starts_with("<!--") {
            let end = input[cursor + 4..]
                .find("-->")
                .map(|offset| cursor + 4 + offset + 3)
                .unwrap_or(input.len());
            if !output.is_empty() && !output.ends_with('\n') {
                output.push('\n');
            }
            output.push_str(&input[cursor..end]);
            output.push('\n');
            cursor = end;
            continue;
        }
        if !in_tag && input[cursor..].starts_with("<![CDATA[") {
            let end = input[cursor + 9..]
                .find("]]>")
                .map(|offset| cursor + 9 + offset + 3)
                .unwrap_or(input.len());
            if !output.is_empty() && !output.ends_with('\n') {
                output.push('\n');
            }
            output.push_str(&input[cursor..end]);
            output.push('\n');
            cursor = end;
            continue;
        }
        if !in_tag && input[cursor..].starts_with("<?") {
            let end = input[cursor + 2..]
                .find("?>")
                .map(|offset| cursor + 2 + offset + 2)
                .unwrap_or(input.len());
            if !output.is_empty() && !output.ends_with('\n') {
                output.push('\n');
            }
            output.push_str(&input[cursor..end]);
            output.push('\n');
            cursor = end;
            continue;
        }

        let character = input[cursor..]
            .chars()
            .next()
            .expect("cursor is a valid non-empty UTF-8 boundary");
        if character == '<' && !in_tag {
            if !output.is_empty() && !output.ends_with('\n') {
                output.push('\n');
            }
            in_tag = true;
        }
        output.push(character);
        if in_tag {
            if let Some(active_quote) = quote {
                if character == active_quote {
                    quote = None;
                }
            } else if character == '\'' || character == '"' {
                quote = Some(character);
            } else if character == '>' {
                in_tag = false;
                output.push('\n');
            }
        }
        cursor += character.len_utf8();
    }
    output
}

/// Parse CML with explicit resource limits.
pub fn parse_cml_with_limits(
    input: &str,
    limits: &CmlParseLimits,
) -> Result<(Molecule, Vec<(f64, f64)>), CmlError> {
    let lines = cml_lines(input, limits)?;
    let (atoms, bonds) = collect_cml_data(lines, limits)?;
    build_cml(atoms, bonds)
}

fn validate_strict_cml_structure(input: &str) -> Result<(), CmlError> {
    let mut stack: Vec<String> = Vec::new();
    let mut saw_molecule = false;
    let mut saw_atom = false;
    let mut root_closed = false;
    let mut cursor = 0usize;
    while let Some(relative_start) = input[cursor..].find('<') {
        let start = cursor + relative_start;
        if let Some(end) = strict_cml_special_end(input, start)? {
            cursor = end;
            continue;
        }
        let relative_end = find_strict_tag_end(&input[start..]).ok_or_else(|| {
            CmlError::MalformedXml("unterminated tag or quoted attribute".to_string())
        })?;
        let end = start + relative_end;
        let raw = input[start + 1..end].trim();
        cursor = end + 1;
        if raw.is_empty() || raw.starts_with('!') || raw.starts_with('?') {
            continue;
        }
        if let Some(name) = raw.strip_prefix('/') {
            let name = tag_name(name);
            let Some(open) = stack.pop() else {
                return Err(CmlError::MalformedXml(format!(
                    "unexpected closing tag </{name}>"
                )));
            };
            if open != name {
                return Err(CmlError::MalformedXml(format!(
                    "closing tag </{name}> does not match <{open}>"
                )));
            }
            if stack.is_empty() {
                root_closed = true;
            }
            continue;
        }
        let self_closing = raw.ends_with('/');
        let name = tag_name(raw);
        if stack.is_empty() && root_closed {
            return Err(CmlError::MalformedXml("multiple root elements".to_string()));
        }
        if name == "molecule" {
            saw_molecule = true;
        }
        if name == "atom" {
            saw_atom = true;
        }
        if !self_closing {
            stack.push(name);
        } else if stack.is_empty() {
            root_closed = true;
        }
    }
    if !stack.is_empty() {
        return Err(CmlError::MalformedXml(format!(
            "unclosed tag <{}>",
            stack.last().expect("non-empty stack")
        )));
    }
    if !saw_molecule {
        return Err(CmlError::MissingMolecule);
    }
    if !saw_atom {
        return Err(CmlError::EmptyMolecule);
    }
    Ok(())
}

fn strict_cml_special_end(input: &str, start: usize) -> Result<Option<usize>, CmlError> {
    let suffix = &input[start..];
    let (prefix, terminator, description) = if suffix.starts_with("<!--") {
        (4, "-->", "comment")
    } else if suffix.starts_with("<![CDATA[") {
        (9, "]]>", "CDATA section")
    } else if suffix.starts_with("<?") {
        (2, "?>", "processing instruction")
    } else {
        return Ok(None);
    };
    input[start + prefix..]
        .find(terminator)
        .map(|offset| Some(start + prefix + offset + terminator.len()))
        .ok_or_else(|| CmlError::MalformedXml(format!("unterminated {description}")))
}

/// Find the end of one XML-ish tag without treating a `>` inside a quoted
/// attribute as the tag boundary.  The strict parser is intentionally a small
/// structural gate rather than a general XML implementation, but its boundary
/// scanner must agree with [`split_cml_tags`] on this basic lexical rule.
fn find_strict_tag_end(input: &str) -> Option<usize> {
    let mut quote = None;
    for (offset, character) in input.char_indices().skip(1) {
        if let Some(active_quote) = quote {
            if character == active_quote {
                quote = None;
            }
        } else if character == '\'' || character == '"' {
            quote = Some(character);
        } else if character == '>' {
            return Some(offset);
        }
    }
    None
}

fn tag_name(raw: &str) -> String {
    raw.trim_start_matches('/')
        .trim_end_matches('/')
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_ascii_lowercase()
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
/// `coords` — optional `(x, y)` pairs in atom order; finite values retain binary64 precision.
/// If `None`, no coordinate attributes are written.
/// Every atom receives a total `hydrogenCount`, combining non-node H with
/// explicit H neighbors, so the reader never adds or drops H through inference.
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
            parts.push(format!("x2=\"{x}\""));
            parts.push(format!("y2=\"{y}\""));
        }
        if atom.charge != 0 {
            parts.push(format!("formalCharge=\"{}\"", atom.charge));
        }
        if let Some(iso) = atom.isotope {
            parts.push(format!("isotopeNumber=\"{iso}\""));
        }
        let h = usize::from(implicit_hcount(mol, idx))
            + mol
                .neighbors(idx)
                .filter(|(neighbor, _)| mol.atom(*neighbor).element == Element::H)
                .count();
        parts.push(format!("hydrogenCount=\"{h}\""));
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
    fn strict_cml_accepts_writer_output() {
        let (mol, coords) = parse_cml_strict(ETHANOL_CML).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(coords.len(), 3);
    }

    #[test]
    fn strict_cml_accepts_compact_xml() {
        let compact = "<molecule><atomArray><atom id=\"a1\" elementType=\"C\"/></atomArray><bondArray/></molecule>";
        let (mol, coords) = parse_cml_strict(compact).unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(coords, vec![(0.0, 0.0)]);
    }

    #[test]
    fn strict_cml_scanner_keeps_greater_than_inside_quoted_attribute() {
        let compact = "<molecule><atomArray><atom id=\"a>1\" elementType=\"C\"/></atomArray><bondArray/></molecule>";
        let (mol, coords) = parse_cml_strict(compact).unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(coords.len(), 1);
    }

    #[test]
    fn strict_cml_scanner_rejects_unterminated_quoted_attribute() {
        let malformed =
            "<molecule><atomArray><atom id=\"a1 elementType=\"C\"/></atomArray></molecule>";
        assert!(matches!(
            parse_cml_strict(malformed),
            Err(CmlError::MalformedXml(_))
        ));
    }

    #[test]
    fn strict_cml_does_not_parse_markup_inside_comment_or_cdata() {
        let input = concat!(
            "<?xml version=\"1.0\"?>",
            "<molecule>",
            "<!-- <atom id=\"fake\" elementType=\"N\"/> > -->",
            "<atomArray><atom id=\"a1\" elementType=\"C\"/></atomArray>",
            "<metadata><![CDATA[<atom id=\"also-fake\"/> >]]></metadata>",
            "<bondArray/></molecule>"
        );
        let (mol, _) = parse_cml_strict(input).unwrap();
        assert_eq!(mol.atom_count(), 1);
    }

    #[test]
    fn strict_cml_rejects_unterminated_comment_and_cdata() {
        assert!(matches!(
            parse_cml_strict("<molecule><!-- <atom id=\"a1\"/>"),
            Err(CmlError::MalformedXml(_))
        ));
        assert!(matches!(
            parse_cml_strict("<molecule><![CDATA[unfinished"),
            Err(CmlError::MalformedXml(_))
        ));
    }

    #[test]
    fn strict_cml_rejects_empty_and_unbalanced_documents() {
        assert!(matches!(
            parse_cml_strict("<cml/>"),
            Err(CmlError::MissingMolecule)
        ));
        assert!(matches!(
            parse_cml_strict("<molecule><atomArray>"),
            Err(CmlError::MalformedXml(_))
        ));
        assert!(matches!(
            parse_cml_strict("<molecule><atomArray/></molecule>"),
            Err(CmlError::EmptyMolecule)
        ));
    }

    #[test]
    fn strict_cml_rejects_multiple_root_elements() {
        let input = concat!(
            "<molecule><atomArray><atom id=\"a1\" elementType=\"C\"/></atomArray></molecule>",
            "<molecule><atomArray><atom id=\"a2\" elementType=\"N\"/></atomArray></molecule>"
        );
        assert!(matches!(
            parse_cml_strict(input),
            Err(CmlError::MalformedXml(message)) if message == "multiple root elements"
        ));
    }

    #[test]
    fn lenient_cml_boundary_remains_compatible() {
        let (mol, coords) = parse_cml("<cml/>").unwrap();
        assert_eq!(mol.atom_count(), 0);
        assert!(coords.is_empty());
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
    fn shared_cml_roundtrip_contract_matches() {
        let document: serde_json::Value = serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/cross_binding_contract.json"
        )))
        .expect("contract JSON");
        let contract = &document["cml_contract"];
        let (mol, coords) = parse_cml(contract["input"].as_str().unwrap()).unwrap();
        let serialized = write_cml(&mol, Some(&coords));
        let (roundtripped, roundtripped_coords) = parse_cml(&serialized).unwrap();
        assert_eq!(
            roundtripped.atom_count(),
            contract["expected"]["atom_count"]
        );
        assert_eq!(roundtripped.bond_count(), 2);
        assert_eq!(roundtripped_coords.len(), coords.len());
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
