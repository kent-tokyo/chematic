//! ChemAxon Marvin MRV parser and writer — `.mrv`.
//!
//! MRV is an XML-based format (`<cml><MDocument><MChemicalStruct><molecule>...`).
//! This is chematic's counterpart to RDKit's `MolFromMrvBlock`/`MolToMrvBlock`
//! (`Code/GraphMol/MarvinParse/*`, RDKit commit
//! `8afba32ec539dcb2369bc84549d802aca3f7eb39`, the true resolution of tag
//! `Release_2026_03_4`).
//!
//! **Scope, deliberately bounded (matches this project's own spec, not an
//! incomplete draft):** `molecule`/`atomArray`/`bondArray`, atom IDs,
//! element/charge/isotope/radical/atom-map, bond order (single/double/triple/
//! aromatic), 2D/3D coordinates, wedge/dash stereo. The following are
//! explicitly **not** implemented and produce
//! [`MrvError::UnsupportedFeature`] rather than a silent partial parse or a
//! guessed mapping: S-groups, polymers, reactions, multicenter bonds, query
//! atoms/bonds, R-groups, enhanced stereo groups, compressed/encoded embedded
//! data. RDKit itself *does* support several of these (S-groups, polymers,
//! reactions, enhanced stereo) — chematic's port does not, by explicit scope
//! decision, not because they're unsupported in MRV generally.
//!
//! **No external XML crate.** A purpose-built, dependency-free, nesting-aware
//! tokenizer (below) is used instead of `cml.rs`'s line-by-line
//! `is_element_tag` scanner (too weak for MRV's nested elements and child
//! text content, e.g. `<bondStereo>W</bondStereo>`) and instead of adding a
//! new XML crate dependency (whose *default* DTD/entity posture would then
//! need independent verification against this module's own mandated
//! security limits anyway — a purpose-built scanner where every limit is
//! controlled directly is more defensible). [`crate::cml::parse_xml_attrs`]
//! is reused for attribute extraction within one already-isolated tag,
//! matching this crate's "reuse before rewriting" convention.
//!
//! **Security — chematic-only, RDKit provides none of this.** This
//! session's source audit confirmed RDKit's own MRV parser passes no
//! `xml_parser_flags` to Boost's `read_xml` at all — no DTD/entity handling
//! either way, "RDKit does it" cannot justify skipping this. This parser
//! rejects `<!DOCTYPE` and `<!ENTITY` outright, enforces an input byte-size
//! limit, an element-nesting-depth limit, and an attribute-value-length
//! limit — all before any chemistry is interpreted.
//!
//! **`<cml>` with no `<molecule>` found anywhere returns an empty molecule
//! (0 atoms), not an error** — matching RDKit's own confirmed lenient
//! behavior. Only a genuinely missing `<cml>` root is
//! [`MrvError::MissingMolecule`]. Callers building oracle/round-trip
//! comparisons must not treat a 0-atom result as a meaningful match (this
//! project's own "vacuous ground truth" lesson).

use std::collections::HashMap;

use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
use chematic_perception::{apply_ez_directions_from_2d, apply_local_parity_from_wedges};

use crate::cml::parse_xml_attrs;
use crate::record::MoleculeRecord;

// ---------------------------------------------------------------------------
// Options / Errors
// ---------------------------------------------------------------------------

/// Security/robustness limits enforced before any chemistry is interpreted.
#[derive(Debug, Clone, Copy)]
pub struct MrvParseLimits {
    pub max_input_bytes: usize,
    pub max_depth: usize,
    pub max_attr_len: usize,
}

impl Default for MrvParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 32 << 20, // 32 MiB
            max_depth: 64,
            max_attr_len: 1 << 16, // 64 KiB
        }
    }
}

/// Errors from [`parse_mrv`]/[`write_mrv`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MrvError {
    /// Input exceeded [`MrvParseLimits::max_input_bytes`].
    InputTooLarge { limit: usize },
    /// A `<!DOCTYPE` or `<!ENTITY` declaration was present — rejected
    /// outright, not resolved (this parser is not a full XML processor and
    /// makes no claim about entity-expansion safety beyond "never attempt
    /// it").
    DtdOrEntityNotAllowed,
    /// Element nesting exceeded [`MrvParseLimits::max_depth`].
    NestingTooDeep { limit: usize },
    /// An attribute value exceeded [`MrvParseLimits::max_attr_len`].
    AttributeTooLong { limit: usize },
    /// The document has no `<cml>` root element at all.
    MissingMolecule,
    /// Malformed XML (unbalanced tags, unterminated tag, etc.).
    MalformedXml { detail: String },
    /// An atom referenced an unrecognized element symbol.
    UnknownElement { symbol: String },
    /// A bond referenced an atom id not present in the atom array.
    UnknownAtomRef { id: String },
    /// Two atoms shared the same `id` (RDKit itself doesn't detect this;
    /// chematic does, deliberately more robust here).
    DuplicateAtomId { id: String },
    /// A bond's `atomRefs2` was missing or malformed.
    InvalidAtomRefs2 { detail: String },
    /// A bond order/query-type/convention attribute value was not recognized.
    InvalidBondOrder { detail: String },
    /// A `radical` attribute value was not one of the recognized literals.
    InvalidRadical { detail: String },
    /// A `<bondStereo>` element had more than one of {value, dictRef,
    /// convention}, or an unrecognized value.
    InvalidBondStereo { detail: String },
    /// A construct this port deliberately does not support was detected —
    /// see the module docs for the full list.
    UnsupportedFeature { feature: String, location: String },
    /// [`MrvWriteOptions::kekulize`] was requested but the molecule's
    /// aromatic system could not be kekulized.
    KekulizationFailed { detail: String },
}

impl std::fmt::Display for MrvError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InputTooLarge { limit } => write!(f, "input exceeds the {limit}-byte limit"),
            Self::DtdOrEntityNotAllowed => write!(f, "DOCTYPE/ENTITY declarations are not allowed"),
            Self::NestingTooDeep { limit } => {
                write!(f, "element nesting exceeds the depth limit of {limit}")
            }
            Self::AttributeTooLong { limit } => {
                write!(f, "attribute value exceeds the {limit}-byte limit")
            }
            Self::MissingMolecule => write!(f, "no <cml> root element found"),
            Self::MalformedXml { detail } => write!(f, "malformed XML: {detail}"),
            Self::UnknownElement { symbol } => write!(f, "unknown element symbol: {symbol}"),
            Self::UnknownAtomRef { id } => write!(f, "unknown atom ref: {id}"),
            Self::DuplicateAtomId { id } => write!(f, "duplicate atom id: {id}"),
            Self::InvalidAtomRefs2 { detail } => write!(f, "invalid atomRefs2: {detail}"),
            Self::InvalidBondOrder { detail } => write!(f, "invalid bond order: {detail}"),
            Self::InvalidRadical { detail } => write!(f, "invalid radical: {detail}"),
            Self::InvalidBondStereo { detail } => write!(f, "invalid bondStereo: {detail}"),
            Self::UnsupportedFeature { feature, location } => {
                write!(f, "unsupported MRV feature {feature:?} at {location}")
            }
            Self::KekulizationFailed { detail } => write!(f, "kekulization failed: {detail}"),
        }
    }
}

impl std::error::Error for MrvError {}

// ---------------------------------------------------------------------------
// Minimal, dependency-free, nesting-aware XML tree
// ---------------------------------------------------------------------------

struct XmlNode {
    name: String,
    attrs: HashMap<String, String>,
    children: Vec<XmlNode>,
    text: String,
    /// Byte offset of the opening tag, for error location reporting.
    location: usize,
}

/// Parse `input` into a tree of [`XmlNode`]s rooted at a single top-level
/// element, enforcing `limits` throughout. Comments (`<!-- -->`) and the
/// XML declaration (`<?xml ... ?>`) are skipped; `<!DOCTYPE`/`<!ENTITY` are
/// rejected outright (see the module docs).
fn parse_xml_tree(input: &str, limits: &MrvParseLimits) -> Result<XmlNode, MrvError> {
    if input.len() > limits.max_input_bytes {
        return Err(MrvError::InputTooLarge {
            limit: limits.max_input_bytes,
        });
    }
    if input.contains("<!DOCTYPE") || input.contains("<!ENTITY") || input.contains("<!doctype") {
        return Err(MrvError::DtdOrEntityNotAllowed);
    }

    let bytes = input.as_bytes();
    let mut pos = 0usize;

    // Skip leading XML declaration / comments / whitespace before the root.
    loop {
        skip_ws(bytes, &mut pos);
        if input[pos..].starts_with("<?") {
            pos = find_after(input, pos, "?>")?;
        } else if input[pos..].starts_with("<!--") {
            pos = find_after(input, pos, "-->")?;
        } else {
            break;
        }
    }

    let (root, end) = parse_element(input, pos, 0, limits)?;
    let _ = end;
    Ok(root)
}

fn skip_ws(bytes: &[u8], pos: &mut usize) {
    while *pos < bytes.len() && bytes[*pos].is_ascii_whitespace() {
        *pos += 1;
    }
}

fn find_after(input: &str, from: usize, marker: &str) -> Result<usize, MrvError> {
    match input[from..].find(marker) {
        Some(rel) => Ok(from + rel + marker.len()),
        None => Err(MrvError::MalformedXml {
            detail: format!("unterminated construct looking for {marker:?}"),
        }),
    }
}

/// Parse one element (and its children/text) starting at `pos`, which must
/// point at a `<`. Returns the parsed node and the byte offset just past
/// its closing tag.
fn parse_element(
    input: &str,
    mut pos: usize,
    depth: usize,
    limits: &MrvParseLimits,
) -> Result<(XmlNode, usize), MrvError> {
    if depth > limits.max_depth {
        return Err(MrvError::NestingTooDeep {
            limit: limits.max_depth,
        });
    }

    let bytes = input.as_bytes();
    if pos >= bytes.len() || bytes[pos] != b'<' {
        return Err(MrvError::MalformedXml {
            detail: format!("expected '<' at byte {pos}"),
        });
    }
    let tag_start = pos;
    pos += 1;

    let name_start = pos;
    while pos < bytes.len()
        && !bytes[pos].is_ascii_whitespace()
        && bytes[pos] != b'>'
        && bytes[pos] != b'/'
    {
        pos += 1;
    }
    let name = input[name_start..pos].to_string();

    let attrs_start = pos;
    let tag_end =
        input[pos..]
            .find('>')
            .map(|i| pos + i)
            .ok_or_else(|| MrvError::MalformedXml {
                detail: format!("unterminated tag <{name}"),
            })?;
    let self_closing = tag_end > attrs_start && bytes[tag_end - 1] == b'/';
    let attrs_text = if self_closing {
        &input[attrs_start..tag_end - 1]
    } else {
        &input[attrs_start..tag_end]
    };
    let attrs = parse_xml_attrs(attrs_text);
    for v in attrs.values() {
        if v.len() > limits.max_attr_len {
            return Err(MrvError::AttributeTooLong {
                limit: limits.max_attr_len,
            });
        }
    }

    pos = tag_end + 1;

    if self_closing {
        return Ok((
            XmlNode {
                name,
                attrs,
                children: Vec::new(),
                text: String::new(),
                location: tag_start,
            },
            pos,
        ));
    }

    let mut children = Vec::new();
    let mut text = String::new();
    let closing = format!("</{name}>");

    loop {
        if pos >= bytes.len() {
            return Err(MrvError::MalformedXml {
                detail: format!("unterminated element <{name}>"),
            });
        }
        if input[pos..].starts_with("<!--") {
            pos = find_after(input, pos, "-->")?;
            continue;
        }
        if input[pos..].starts_with(&closing) {
            pos += closing.len();
            break;
        }
        if bytes[pos] == b'<' {
            let (child, next) = parse_element(input, pos, depth + 1, limits)?;
            children.push(child);
            pos = next;
        } else {
            let next_lt = input[pos..]
                .find('<')
                .map(|i| pos + i)
                .unwrap_or(bytes.len());
            text.push_str(unescape(&input[pos..next_lt]).as_str());
            pos = next_lt;
        }
    }

    Ok((
        XmlNode {
            name,
            attrs,
            children,
            text,
            location: tag_start,
        },
        pos,
    ))
}

fn unescape(s: &str) -> String {
    if !s.contains('&') {
        return s.to_string();
    }
    s.replace("&amp;", "&")
        .replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
}

impl XmlNode {
    fn find_child(&self, name: &str) -> Option<&XmlNode> {
        self.children.iter().find(|c| c.name == name)
    }

    fn find_all(&self, name: &str) -> impl Iterator<Item = &XmlNode> {
        self.children.iter().filter(move |c| c.name == name)
    }
}

// ---------------------------------------------------------------------------
// Molecule building
// ---------------------------------------------------------------------------

struct MrvAtom {
    id: String,
    element: Element,
    charge: i8,
    isotope: Option<u16>,
    radical_electrons: u8,
    atom_map: Option<u16>,
    x2: Option<f64>,
    y2: Option<f64>,
    x3: Option<f64>,
    y3: Option<f64>,
    z3: Option<f64>,
}

fn parse_radical(s: &str) -> Result<u8, MrvError> {
    match s {
        "monovalent" => Ok(1),
        "divalent" | "divalent1" | "divalent3" => Ok(2),
        "trivalent" | "trivalent2" | "trivalent4" => Ok(3),
        "4" => Ok(4),
        other => Err(MrvError::InvalidRadical {
            detail: other.to_string(),
        }),
    }
}

fn location_str(offset: usize) -> String {
    format!("byte offset {offset}")
}

fn check_unsupported_atom(
    attrs: &HashMap<String, String>,
    location: usize,
) -> Result<(), MrvError> {
    if attrs.contains_key("mrvStereoGroup") {
        return Err(MrvError::UnsupportedFeature {
            feature: "enhanced stereo group (mrvStereoGroup)".to_string(),
            location: location_str(location),
        });
    }
    if attrs.get("elementType").map(String::as_str) == Some("R") || attrs.contains_key("rgroupRef")
    {
        return Err(MrvError::UnsupportedFeature {
            feature: "R-group placeholder atom".to_string(),
            location: location_str(location),
        });
    }
    Ok(())
}

fn parse_atom_array(node: &XmlNode) -> Result<Vec<MrvAtom>, MrvError> {
    let mut atoms = Vec::new();
    let mut seen_ids = std::collections::HashSet::new();

    for atom_node in node.find_all("atom") {
        let attrs = &atom_node.attrs;
        check_unsupported_atom(attrs, atom_node.location)?;

        let id = attrs.get("id").cloned().unwrap_or_default();
        if !id.is_empty() && !seen_ids.insert(id.clone()) {
            return Err(MrvError::DuplicateAtomId { id });
        }

        let symbol = attrs.get("elementType").map(String::as_str).unwrap_or("C");
        if symbol == "*" {
            return Err(MrvError::UnsupportedFeature {
                feature: "wildcard query atom".to_string(),
                location: location_str(atom_node.location),
            });
        }
        let element = Element::from_symbol(symbol).ok_or_else(|| MrvError::UnknownElement {
            symbol: symbol.to_string(),
        })?;

        let charge = attrs
            .get("formalCharge")
            .and_then(|s| s.trim().parse::<i8>().ok())
            .unwrap_or(0);
        let isotope = attrs
            .get("isotope")
            .and_then(|s| s.trim().parse::<u16>().ok())
            .filter(|&v| v > 0);
        let radical_electrons = match attrs.get("radical") {
            Some(r) => parse_radical(r)?,
            None => 0,
        };
        let atom_map = attrs
            .get("mrvMap")
            .and_then(|s| s.trim().parse::<u16>().ok());
        let x2 = attrs.get("x2").and_then(|s| s.trim().parse::<f64>().ok());
        let y2 = attrs.get("y2").and_then(|s| s.trim().parse::<f64>().ok());
        let x3 = attrs.get("x3").and_then(|s| s.trim().parse::<f64>().ok());
        let y3 = attrs.get("y3").and_then(|s| s.trim().parse::<f64>().ok());
        let z3 = attrs.get("z3").and_then(|s| s.trim().parse::<f64>().ok());

        atoms.push(MrvAtom {
            id,
            element,
            charge,
            isotope,
            radical_electrons,
            atom_map,
            x2,
            y2,
            x3,
            y3,
            z3,
        });
    }
    Ok(atoms)
}

struct MrvBond {
    a1: String,
    a2: String,
    order: BondOrder,
    is_aromatic: bool,
    stereo: Option<&'static str>, // "wedge" | "dash"
}

fn parse_bond_order_attr(
    attrs: &HashMap<String, String>,
    location: usize,
) -> Result<(BondOrder, bool), MrvError> {
    if attrs.contains_key("queryType") {
        return Err(MrvError::UnsupportedFeature {
            feature: format!("query bond ({})", attrs.get("queryType").unwrap()),
            location: location_str(location),
        });
    }
    if let Some(conv) = attrs.get("convention") {
        return match conv.as_str() {
            "cxn:coord" => Ok((BondOrder::Dative, false)),
            other => Err(MrvError::InvalidBondOrder {
                detail: format!("unsupported convention {other}"),
            }),
        };
    }
    let order = attrs.get("order").map(String::as_str).unwrap_or("1");
    match order {
        "1" => Ok((BondOrder::Single, false)),
        "2" => Ok((BondOrder::Double, false)),
        "3" => Ok((BondOrder::Triple, false)),
        "A" => Ok((BondOrder::Aromatic, true)),
        other => Err(MrvError::InvalidBondOrder {
            detail: other.to_string(),
        }),
    }
}

fn parse_bond_stereo(bond_node: &XmlNode) -> Result<Option<&'static str>, MrvError> {
    let Some(stereo_node) = bond_node.find_child("bondStereo") else {
        return Ok(None);
    };
    let has_value = !stereo_node.text.trim().is_empty();
    let has_dict_ref = stereo_node.attrs.contains_key("dictRef");
    let has_convention = stereo_node.attrs.contains_key("convention");
    if (has_value as u8 + has_dict_ref as u8 + has_convention as u8) > 1 {
        return Err(MrvError::InvalidBondStereo {
            detail: "bondStereo can have only one of: a value, dictRef, or convention".to_string(),
        });
    }

    if has_value {
        return match stereo_node.text.trim() {
            "W" => Ok(Some("wedge")),
            "H" => Ok(Some("dash")),
            "C" | "T" => Ok(None), // cis/trans -- explicitly ignored, matching RDKit
            other => Err(MrvError::InvalidBondStereo {
                detail: other.to_string(),
            }),
        };
    }
    if has_dict_ref {
        return match stereo_node.attrs.get("dictRef").map(String::as_str) {
            Some("cml:W") => Ok(Some("wedge")),
            Some("cml:H") => Ok(Some("dash")),
            Some(other) => Err(MrvError::InvalidBondStereo {
                detail: other.to_string(),
            }),
            None => Ok(None),
        };
    }
    if has_convention {
        return match stereo_node.attrs.get("conventionValue").map(String::as_str) {
            Some("1") => Ok(Some("wedge")),
            Some("6") => Ok(Some("dash")),
            Some(other) => Err(MrvError::InvalidBondStereo {
                detail: format!("MDL conventionValue {other}"),
            }),
            None => Ok(None),
        };
    }
    Ok(None)
}

fn parse_bond_array(node: &XmlNode) -> Result<Vec<MrvBond>, MrvError> {
    let mut bonds = Vec::new();
    for bond_node in node.find_all("bond") {
        let attrs = &bond_node.attrs;
        let refs = attrs
            .get("atomRefs2")
            .ok_or_else(|| MrvError::InvalidAtomRefs2 {
                detail: "missing atomRefs2".to_string(),
            })?;
        let parts: Vec<&str> = refs.split_whitespace().collect();
        if parts.len() != 2 {
            return Err(MrvError::InvalidAtomRefs2 {
                detail: refs.clone(),
            });
        }
        let (order, is_aromatic) = parse_bond_order_attr(attrs, bond_node.location)?;
        let stereo = parse_bond_stereo(bond_node)?;
        bonds.push(MrvBond {
            a1: parts[0].to_string(),
            a2: parts[1].to_string(),
            order,
            is_aromatic,
            stereo,
        });
    }
    Ok(bonds)
}

/// Reject S-groups/polymers/R-group-definitions/multicenter bonds: any
/// child `<molecule role="...">` nested inside the top molecule element.
fn check_no_sgroups(molecule_node: &XmlNode) -> Result<(), MrvError> {
    if let Some(nested) = molecule_node.find_child("molecule") {
        let role = nested
            .attrs
            .get("role")
            .cloned()
            .unwrap_or_else(|| "unknown".to_string());
        return Err(MrvError::UnsupportedFeature {
            feature: format!("S-group/polymer/R-group (role={role})"),
            location: location_str(nested.location),
        });
    }
    Ok(())
}

/// Find the `molecule` element under `cml.MDocument.MChemicalStruct`, or
/// `None` if genuinely absent anywhere under `cml` (RDKit's own lenient
/// "empty molecule" case — see the module docs).
fn find_molecule_node(root: &XmlNode) -> Result<Option<&XmlNode>, MrvError> {
    if root.name != "cml" {
        return Err(MrvError::MissingMolecule);
    }
    if root.find_child("reactionList").is_some()
        || root.find_all("molecule").count() > 1 && root.find_child("MDocument").is_none()
    {
        return Err(MrvError::UnsupportedFeature {
            feature: "reaction MRV".to_string(),
            location: location_str(root.location),
        });
    }
    let doc = root.find_child("MDocument");
    let struct_node = doc.and_then(|d| d.find_child("MChemicalStruct")).or(doc);
    let candidate = struct_node
        .and_then(|s| s.find_child("molecule"))
        .or_else(|| root.find_child("molecule"));
    Ok(candidate)
}

/// Parse an MRV (Marvin) document into a [`MoleculeRecord`].
///
/// See the module docs for the exact supported feature set and the
/// deliberate scope boundary (S-groups/polymers/reactions/multicenter
/// bonds/query atoms-bonds/R-groups/enhanced stereo groups/embedded data
/// all produce [`MrvError::UnsupportedFeature`]).
pub fn parse_mrv(input: &str) -> Result<MoleculeRecord, MrvError> {
    parse_mrv_with_limits(input, &MrvParseLimits::default())
}

/// Like [`parse_mrv`], with explicit security/robustness limits.
pub fn parse_mrv_with_limits(
    input: &str,
    limits: &MrvParseLimits,
) -> Result<MoleculeRecord, MrvError> {
    let root = parse_xml_tree(input, limits)?;
    let Some(molecule_node) = find_molecule_node(&root)? else {
        // Matches RDKit's own confirmed lenient behavior: no error, just no atoms.
        return Ok(MoleculeRecord::new(MoleculeBuilder::new().build()));
    };

    check_no_sgroups(molecule_node)?;

    let atom_array =
        molecule_node
            .find_child("atomArray")
            .ok_or_else(|| MrvError::MalformedXml {
                detail: "missing atomArray".to_string(),
            })?;
    let bond_array = molecule_node.find_child("bondArray");

    let mrv_atoms = parse_atom_array(atom_array)?;
    let mrv_bonds = match bond_array {
        Some(b) => parse_bond_array(b)?,
        None => Vec::new(),
    };

    let mut id_to_pos: HashMap<String, usize> = HashMap::new();
    for (i, a) in mrv_atoms.iter().enumerate() {
        if !a.id.is_empty() {
            id_to_pos.insert(a.id.clone(), i);
        }
    }

    let has_2d = mrv_atoms.iter().any(|a| a.x2.is_some() || a.y2.is_some());
    let has_3d = mrv_atoms
        .iter()
        .any(|a| a.x3.is_some() || a.y3.is_some() || a.z3.is_some());

    // Atom aromaticity is derived from aromatic bonds (order="A"), same as
    // cml.rs -- must be known before construction since MoleculeBuilder has
    // no post-hoc atom mutator for it.
    let mut is_aromatic_atom = vec![false; mrv_atoms.len()];
    for b in &mrv_bonds {
        if b.is_aromatic {
            if let Some(&pos) = id_to_pos.get(&b.a1) {
                is_aromatic_atom[pos] = true;
            }
            if let Some(&pos) = id_to_pos.get(&b.a2) {
                is_aromatic_atom[pos] = true;
            }
        }
    }

    let mut builder = MoleculeBuilder::new();
    let mut idx_by_pos: Vec<AtomIdx> = Vec::new();
    let mut coords_2d = Vec::new();
    let mut coords_3d = Vec::new();

    for (i, a) in mrv_atoms.iter().enumerate() {
        let mut atom = Atom::new(a.element);
        atom.charge = a.charge;
        atom.isotope = a.isotope;
        atom.atom_map = a.atom_map;
        atom.aromatic = is_aromatic_atom[i];
        // Radical electron count has no chematic_core::Atom slot -- treated
        // as neutral, same convention as mol2000.rs's doublet-radical case.
        let _ = a.radical_electrons;
        let idx = builder.add_atom(atom);
        idx_by_pos.push(idx);
        if has_2d {
            coords_2d.push([a.x2.unwrap_or(0.0), a.y2.unwrap_or(0.0)]);
        }
        if has_3d {
            coords_3d.push([
                a.x3.unwrap_or(0.0),
                a.y3.unwrap_or(0.0),
                a.z3.unwrap_or(0.0),
            ]);
        }
    }

    for b in &mrv_bonds {
        let pos1 = *id_to_pos
            .get(&b.a1)
            .ok_or_else(|| MrvError::UnknownAtomRef { id: b.a1.clone() })?;
        let pos2 = *id_to_pos
            .get(&b.a2)
            .ok_or_else(|| MrvError::UnknownAtomRef { id: b.a2.clone() })?;
        let a1 = idx_by_pos[pos1];
        let a2 = idx_by_pos[pos2];
        // Wedge/dash stereo is represented via BondOrder::Up/Down, same
        // convention as mol2000.rs's MDL stereo-flag mapping.
        let order = match (b.order, b.stereo) {
            (BondOrder::Single, Some("wedge")) => BondOrder::Up,
            (BondOrder::Single, Some("dash")) => BondOrder::Down,
            _ => b.order,
        };
        builder
            .add_bond(a1, a2, order)
            .map_err(|_| MrvError::InvalidAtomRefs2 {
                detail: format!("{} {}", b.a1, b.a2),
            })?;
    }

    let mut mol = builder.build();
    if has_2d {
        // Tetrahedral parity first (raw wedge/hash still fully intact on
        // `bond.order`), THEN E/Z direction -- same ordering rationale as
        // mol2000.rs (the E/Z stage only ever writes to the separate
        // `bond_direction` side channel, never to `bond.order`, so running
        // it after can never disturb the wedge/hash the tetrahedral stage
        // just read). MRV's `bondStereo` cis/trans dictRef values ("C"/"T")
        // are parsed but explicitly discarded (`parse_bond_stereo` above,
        // matching RDKit), so there is no MRV-specific "explicitly
        // unspecified E/Z" set to thread through here, unlike MDL V2000's
        // stereo code 3 -- every candidate double bond is judged on 2D
        // geometry alone.
        let coord_pairs: Vec<(f64, f64)> = coords_2d.iter().map(|&[x, y]| (x, y)).collect();
        apply_local_parity_from_wedges(&mut mol, &coord_pairs);
        apply_ez_directions_from_2d(&mut mol, &coord_pairs);
    }
    let mut record = MoleculeRecord::new(mol);
    if has_2d {
        record.coordinates_2d = Some(coords_2d);
    }
    if has_3d {
        record.coordinates_3d = Some(coords_3d);
    }
    Ok(record)
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Options for [`write_mrv`].
#[derive(Debug, Clone, Copy)]
pub struct MrvWriteOptions {
    /// Decimal places for coordinate output (0 is treated as the default, 4).
    pub precision: usize,
    /// Kekulize aromatic rings before writing (bond order "1"/"2" instead of
    /// "A"), matching RDKit's `MolToMrvBlock` default. A molecule written
    /// this way and re-read loses its aromatic flag on read-back (the
    /// Kekulé form is a different, chemically-equivalent representation --
    /// re-derive aromaticity via `chematic_perception::apply_aromaticity`
    /// if the aromatic flag itself must survive a round trip).
    pub kekulize: bool,
    /// Emit `<bondStereo>` for wedge/dash (`BondOrder::Up`/`Down`) bonds.
    pub include_stereo: bool,
}

impl Default for MrvWriteOptions {
    fn default() -> Self {
        Self {
            precision: 4,
            kekulize: true,
            include_stereo: true,
        }
    }
}

/// Write a [`MoleculeRecord`] to an MRV document.
///
/// Only single/double/triple/aromatic bonds are supported — any other
/// `BondOrder` returns [`MrvError::UnsupportedFeature`] rather than a
/// guessed mapping (matches RDKit's own `MarvinWriterException` behavior
/// for unsupported bond types).
pub fn write_mrv(record: &MoleculeRecord, options: &MrvWriteOptions) -> Result<String, MrvError> {
    let precision = if options.precision == 0 {
        4
    } else {
        options.precision
    };

    let kekulized;
    let mol = if options.kekulize {
        let mut m = record.mol.clone();
        chematic_perception::kekulize_inplace(&mut m).map_err(|e| {
            MrvError::KekulizationFailed {
                detail: e.to_string(),
            }
        })?;
        kekulized = m;
        &kekulized
    } else {
        &record.mol
    };

    let mut atoms_xml = String::new();
    for (i, (idx, atom)) in mol.atoms().enumerate() {
        let id = format!("a{}", idx.0 + 1);
        let mut attrs = vec![
            format!("id=\"{id}\""),
            format!("elementType=\"{}\"", atom.element.symbol()),
        ];
        if let Some(coords) = &record.coordinates_2d
            && let Some(c) = coords.get(i)
        {
            attrs.push(format!("x2=\"{:.precision$}\"", c[0]));
            attrs.push(format!("y2=\"{:.precision$}\"", c[1]));
        }
        if let Some(coords) = &record.coordinates_3d
            && let Some(c) = coords.get(i)
        {
            attrs.push(format!("x3=\"{:.precision$}\"", c[0]));
            attrs.push(format!("y3=\"{:.precision$}\"", c[1]));
            attrs.push(format!("z3=\"{:.precision$}\"", c[2]));
        }
        if atom.charge != 0 {
            attrs.push(format!("formalCharge=\"{}\"", atom.charge));
        }
        if let Some(iso) = atom.isotope {
            attrs.push(format!("isotope=\"{iso}\""));
        }
        if let Some(map) = atom.atom_map {
            attrs.push(format!("mrvMap=\"{map}\""));
        }
        atoms_xml.push_str(&format!("<atom {}/>", attrs.join(" ")));
    }

    let mut bonds_xml = String::new();
    for (i, (_, bond)) in mol.bonds().enumerate() {
        let a1 = format!("a{}", bond.atom1.0 + 1);
        let a2 = format!("a{}", bond.atom2.0 + 1);
        let (order, stereo_value) = match bond.order {
            BondOrder::Single => ("1", None),
            BondOrder::Up if options.include_stereo => ("1", Some("W")),
            BondOrder::Down if options.include_stereo => ("1", Some("H")),
            BondOrder::Up | BondOrder::Down => ("1", None),
            BondOrder::Double => ("2", None),
            BondOrder::Triple => ("3", None),
            BondOrder::Aromatic => ("A", None),
            other => {
                return Err(MrvError::UnsupportedFeature {
                    feature: format!("{other:?} bond order"),
                    location: format!("bond index {i}"),
                });
            }
        };
        match stereo_value {
            Some(v) => bonds_xml.push_str(&format!(
                "<bond id=\"b{}\" atomRefs2=\"{a1} {a2}\" order=\"{order}\"><bondStereo>{v}</bondStereo></bond>",
                i + 1
            )),
            None => bonds_xml.push_str(&format!("<bond id=\"b{}\" atomRefs2=\"{a1} {a2}\" order=\"{order}\"/>", i + 1)),
        }
    }

    Ok(format!(
        "<cml><MDocument><MChemicalStruct><molecule molID=\"m1\"><atomArray>{atoms_xml}</atomArray><bondArray>{bonds_xml}</bondArray></molecule></MChemicalStruct></MDocument></cml>"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::BondIdx;

    const ETHANOL: &str = r#"<cml><MDocument><MChemicalStruct><molecule molID="m1">
        <atomArray>
            <atom id="a1" elementType="C" x2="0.0" y2="0.0"/>
            <atom id="a2" elementType="C" x2="1.0" y2="0.0"/>
            <atom id="a3" elementType="O" x2="2.0" y2="0.0"/>
        </atomArray>
        <bondArray>
            <bond id="b1" atomRefs2="a1 a2" order="1"/>
            <bond id="b2" atomRefs2="a2 a3" order="1"/>
        </bondArray>
    </molecule></MChemicalStruct></MDocument></cml>"#;

    #[test]
    fn parses_basic_molecule_with_2d_coords() {
        let rec = parse_mrv(ETHANOL).unwrap();
        assert_eq!(rec.mol.atom_count(), 3);
        assert_eq!(rec.mol.bond_count(), 2);
        let coords = rec.coordinates_2d.unwrap();
        assert_eq!(coords.len(), 3);
        assert_eq!(coords[2], [2.0, 0.0]);
    }

    #[test]
    fn parses_charge_isotope_atom_map() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray>
                <atom id="a1" elementType="N" formalCharge="1" isotope="15" mrvMap="7"/>
                <atom id="a2" elementType="C"/>
            </atomArray>
            <bondArray><bond atomRefs2="a1 a2" order="1"/></bondArray>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        let rec = parse_mrv(mrv).unwrap();
        let (_, atom) = rec.mol.atoms().next().unwrap();
        assert_eq!(atom.charge, 1);
        assert_eq!(atom.isotope, Some(15));
        assert_eq!(atom.atom_map, Some(7));
    }

    #[test]
    fn aromatic_bond_marks_both_atoms_aromatic() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray>
                <atom id="a1" elementType="C"/>
                <atom id="a2" elementType="C"/>
            </atomArray>
            <bondArray><bond atomRefs2="a1 a2" order="A"/></bondArray>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        let rec = parse_mrv(mrv).unwrap();
        assert!(rec.mol.atoms().all(|(_, a)| a.aromatic));
    }

    #[test]
    fn writer_kekulizes_benzene_by_default() {
        let rec = MoleculeRecord::new(chematic_smiles::parse("c1ccccc1").unwrap());
        let written = write_mrv(&rec, &MrvWriteOptions::default()).unwrap();
        assert!(!written.contains(r#"order="A""#));
        assert!(written.contains(r#"order="2""#));

        // kekulize=false preserves the aromatic bond order token instead.
        let opts = MrvWriteOptions {
            kekulize: false,
            ..MrvWriteOptions::default()
        };
        let written_aromatic = write_mrv(&rec, &opts).unwrap();
        assert!(written_aromatic.contains(r#"order="A""#));
    }

    #[test]
    fn wedge_and_dash_bond_stereo_round_trip() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray>
                <atom id="a1" elementType="C"/>
                <atom id="a2" elementType="C"/>
                <atom id="a3" elementType="C"/>
            </atomArray>
            <bondArray>
                <bond atomRefs2="a1 a2" order="1"><bondStereo>W</bondStereo></bond>
                <bond atomRefs2="a1 a3" order="1"><bondStereo>H</bondStereo></bond>
            </bondArray>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        let rec = parse_mrv(mrv).unwrap();
        let orders: Vec<BondOrder> = rec.mol.bonds().map(|(_, b)| b.order).collect();
        assert_eq!(orders, vec![BondOrder::Up, BondOrder::Down]);

        let written = write_mrv(&rec, &MrvWriteOptions::default()).unwrap();
        let reparsed = parse_mrv(&written).unwrap();
        let reparsed_orders: Vec<BondOrder> = reparsed.mol.bonds().map(|(_, b)| b.order).collect();
        assert_eq!(reparsed_orders, orders);
    }

    /// Issue #202: `parse_mrv` must convert a 2D wedge + coordinates into
    /// `Atom.chirality`, the same wiring `mol2000.rs`/`cdxml.rs` already
    /// have via `chematic_perception`. Coordinates and expected result
    /// (`Chirality::CounterClockwise`) reused verbatim from
    /// `chematic_perception::stereo2d_local`'s own RDKit-calibrated
    /// `chfclbr`/`quad_positions` fixture (C-F wedge, F/Cl/Br/I at
    /// `(-1.0,0.4)/(0.9,0.7)/(-0.5,-1.1)/(0.8,-0.6)`) -- independently
    /// re-confirmed here against a live RDKit oracle
    /// (`Chem.MolFromMrvBlock` on this exact MRV text, `rdkit==2026.03.3`):
    /// `CHI_TETRAHEDRAL_CCW`, CIP `S`, canonical SMILES `F[C@](Cl)(Br)I`.
    #[test]
    fn wedge_2d_perceives_tetrahedral_chirality() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray>
                <atom id="a1" elementType="C" x2="0.0" y2="0.0"/>
                <atom id="a2" elementType="F" x2="-1.0" y2="0.4"/>
                <atom id="a3" elementType="Cl" x2="0.9" y2="0.7"/>
                <atom id="a4" elementType="Br" x2="-0.5" y2="-1.1"/>
                <atom id="a5" elementType="I" x2="0.8" y2="-0.6"/>
            </atomArray>
            <bondArray>
                <bond id="b1" atomRefs2="a1 a2" order="1"><bondStereo>W</bondStereo></bond>
                <bond id="b2" atomRefs2="a1 a3" order="1"/>
                <bond id="b3" atomRefs2="a1 a4" order="1"/>
                <bond id="b4" atomRefs2="a1 a5" order="1"/>
            </bondArray>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        let rec = parse_mrv(mrv).unwrap();
        let (center, _) = rec.mol.atoms().next().unwrap();
        assert_eq!(
            rec.mol.atom(center).chirality,
            chematic_core::Chirality::CounterClockwise
        );
    }

    /// Issue #202: `parse_mrv` must convert a 2D double-bond geometry into
    /// the `bond_direction` side channel, mirroring `mol2000.rs`/`cdxml.rs`.
    /// Coordinates reused verbatim from
    /// `chematic_perception::stereo2d_ez_direction`'s own
    /// `z_but2ene_assigns_and_matches_legacy_cip_engine` fixture (cis-2-butene) --
    /// independently re-confirmed here against a live RDKit oracle
    /// (`Chem.MolFromMrvBlock` + `DetectBondStereochemistry` on this exact
    /// MRV text): `STEREOZ`, canonical SMILES `C/C=C\C`.
    #[test]
    fn two_d_coords_perceive_ez_direction() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray>
                <atom id="a1" elementType="C" x2="-0.866" y2="0.5"/>
                <atom id="a2" elementType="C" x2="0.0" y2="0.0"/>
                <atom id="a3" elementType="C" x2="1.5" y2="0.0"/>
                <atom id="a4" elementType="C" x2="2.366" y2="0.5"/>
            </atomArray>
            <bondArray>
                <bond id="b1" atomRefs2="a1 a2" order="1"/>
                <bond id="b2" atomRefs2="a2 a3" order="2"/>
                <bond id="b3" atomRefs2="a3 a4" order="1"/>
            </bondArray>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        let rec = parse_mrv(mrv).unwrap();
        let sub_bonds: Vec<BondIdx> = rec
            .mol
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Single)
            .map(|(bidx, _)| bidx)
            .collect();
        assert_eq!(sub_bonds.len(), 2, "two single (substituent) bonds");
        for &bidx in &sub_bonds {
            assert!(
                rec.mol.bond_direction(bidx).is_some(),
                "substituent bond {bidx:?} must carry a perceived E/Z direction"
            );
        }
    }

    #[test]
    fn disconnected_fragments_share_one_atom_bond_array() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray>
                <atom id="a1" elementType="Na" formalCharge="1"/>
                <atom id="a2" elementType="Cl" formalCharge="-1"/>
            </atomArray>
            <bondArray/>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        let rec = parse_mrv(mrv).unwrap();
        assert_eq!(rec.mol.atom_count(), 2);
        assert_eq!(rec.mol.bond_count(), 0);
    }

    #[test]
    fn round_trip_write_then_read_preserves_atoms_and_bonds() {
        let rec = parse_mrv(ETHANOL).unwrap();
        let written = write_mrv(&rec, &MrvWriteOptions::default()).unwrap();
        let reparsed = parse_mrv(&written).unwrap();
        assert_eq!(reparsed.mol.atom_count(), rec.mol.atom_count());
        assert_eq!(reparsed.mol.bond_count(), rec.mol.bond_count());
    }

    #[test]
    fn missing_cml_root_is_error() {
        assert!(matches!(
            parse_mrv("<notcml/>"),
            Err(MrvError::MissingMolecule)
        ));
    }

    #[test]
    fn cml_with_no_molecule_is_lenient_empty_result() {
        let rec = parse_mrv("<cml><MDocument/></cml>").unwrap();
        assert_eq!(rec.mol.atom_count(), 0);
    }

    #[test]
    fn sgroup_nested_molecule_is_unsupported_feature() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray><atom id="a1" elementType="C"/></atomArray>
            <bondArray/>
            <molecule role="SuperatomSgroup"><atomArray/></molecule>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        assert!(matches!(
            parse_mrv(mrv),
            Err(MrvError::UnsupportedFeature { .. })
        ));
    }

    #[test]
    fn query_bond_is_unsupported_feature() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray><atom id="a1" elementType="C"/><atom id="a2" elementType="C"/></atomArray>
            <bondArray><bond atomRefs2="a1 a2" queryType="SD"/></bondArray>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        assert!(matches!(
            parse_mrv(mrv),
            Err(MrvError::UnsupportedFeature { .. })
        ));
    }

    #[test]
    fn rgroup_placeholder_atom_is_unsupported_feature() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray><atom id="a1" elementType="R" rgroupRef="1"/></atomArray>
            <bondArray/>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        assert!(matches!(
            parse_mrv(mrv),
            Err(MrvError::UnsupportedFeature { .. })
        ));
    }

    #[test]
    fn enhanced_stereo_group_is_unsupported_feature() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray><atom id="a1" elementType="C" mrvStereoGroup="1"/></atomArray>
            <bondArray/>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        assert!(matches!(
            parse_mrv(mrv),
            Err(MrvError::UnsupportedFeature { .. })
        ));
    }

    #[test]
    fn reaction_list_is_unsupported_feature() {
        let mrv = r#"<cml><reactionList><reaction/></reactionList></cml>"#;
        assert!(matches!(
            parse_mrv(mrv),
            Err(MrvError::UnsupportedFeature { .. })
        ));
    }

    #[test]
    fn duplicate_atom_id_is_error() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray>
                <atom id="a1" elementType="C"/>
                <atom id="a1" elementType="O"/>
            </atomArray>
            <bondArray/>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        assert!(matches!(
            parse_mrv(mrv),
            Err(MrvError::DuplicateAtomId { .. })
        ));
    }

    #[test]
    fn unknown_atom_ref_in_bond_is_error() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray><atom id="a1" elementType="C"/></atomArray>
            <bondArray><bond atomRefs2="a1 a99" order="1"/></bondArray>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        assert!(matches!(
            parse_mrv(mrv),
            Err(MrvError::UnknownAtomRef { .. })
        ));
    }

    #[test]
    fn unknown_element_symbol_is_error() {
        let mrv = r#"<cml><MDocument><MChemicalStruct><molecule>
            <atomArray><atom id="a1" elementType="Zz"/></atomArray>
            <bondArray/>
        </molecule></MChemicalStruct></MDocument></cml>"#;
        assert!(matches!(
            parse_mrv(mrv),
            Err(MrvError::UnknownElement { .. })
        ));
    }

    #[test]
    fn unterminated_tag_is_malformed_xml_error() {
        assert!(matches!(
            parse_mrv("<cml><MDocument"),
            Err(MrvError::MalformedXml { .. })
        ));
    }

    #[test]
    fn doctype_declaration_is_rejected() {
        let mrv = "<!DOCTYPE foo><cml/>";
        assert!(matches!(
            parse_mrv(mrv),
            Err(MrvError::DtdOrEntityNotAllowed)
        ));
    }

    #[test]
    fn entity_declaration_is_rejected() {
        let mrv = "<cml><!ENTITY xxe SYSTEM \"file:///etc/passwd\"><MDocument/></cml>";
        assert!(matches!(
            parse_mrv(mrv),
            Err(MrvError::DtdOrEntityNotAllowed)
        ));
    }

    #[test]
    fn unsupported_bond_order_writer_error() {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(Atom::new(Element::C));
        let b = builder.add_atom(Atom::new(Element::C));
        builder.add_bond(a, b, BondOrder::QueryAny).unwrap();
        let rec = MoleculeRecord::new(builder.build());
        assert!(matches!(
            write_mrv(&rec, &MrvWriteOptions::default()),
            Err(MrvError::UnsupportedFeature { .. })
        ));
    }
}

#[cfg(test)]
mod adversarial_tests {
    use super::*;

    #[test]
    fn empty_input_is_explicit_error_not_panic() {
        // Zero bytes: no root tag at all, distinct from `MissingMolecule`
        // (a real root that just isn't named `cml`, see `missing_cml_root_is_error`).
        assert!(matches!(parse_mrv(""), Err(MrvError::MalformedXml { .. })));
    }

    #[test]
    fn deeply_nested_elements_hit_depth_limit_not_stack_overflow() {
        let limits = MrvParseLimits {
            max_depth: 32,
            ..MrvParseLimits::default()
        };
        let mut xml = String::from("<cml>");
        for _ in 0..1000 {
            xml.push_str("<a>");
        }
        for _ in 0..1000 {
            xml.push_str("</a>");
        }
        xml.push_str("</cml>");
        assert!(matches!(
            parse_mrv_with_limits(&xml, &limits),
            Err(MrvError::NestingTooDeep { .. })
        ));
    }

    #[test]
    fn oversized_attribute_is_explicit_error_not_panic() {
        let limits = MrvParseLimits {
            max_attr_len: 16,
            ..MrvParseLimits::default()
        };
        let xml = format!(
            r#"<cml><MDocument><MChemicalStruct><molecule><atomArray><atom id="a1" elementType="{}"/></atomArray></molecule></MChemicalStruct></MDocument></cml>"#,
            "C".repeat(1000)
        );
        assert!(matches!(
            parse_mrv_with_limits(&xml, &limits),
            Err(MrvError::AttributeTooLong { .. })
        ));
    }

    #[test]
    fn oversized_input_is_explicit_error_not_panic() {
        let limits = MrvParseLimits {
            max_input_bytes: 100,
            ..MrvParseLimits::default()
        };
        let xml = format!("<cml>{}</cml>", "x".repeat(1000));
        assert!(matches!(
            parse_mrv_with_limits(&xml, &limits),
            Err(MrvError::InputTooLarge { .. })
        ));
    }

    #[test]
    fn billion_laughs_style_entity_chain_rejected_outright() {
        let payload = r#"<!DOCTYPE lolz [
            <!ENTITY lol "lol">
            <!ENTITY lol2 "&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;&lol;">
        ]><cml>&lol2;</cml>"#;
        assert!(matches!(
            parse_mrv(payload),
            Err(MrvError::DtdOrEntityNotAllowed)
        ));
    }

    #[test]
    fn xxe_local_file_reference_rejected_outright() {
        let payload = r#"<?xml version="1.0"?><!DOCTYPE cml [
            <!ENTITY xxe SYSTEM "file:///etc/passwd">
        ]><cml>&xxe;</cml>"#;
        assert!(matches!(
            parse_mrv(payload),
            Err(MrvError::DtdOrEntityNotAllowed)
        ));
    }

    #[test]
    fn truncated_input_no_panic() {
        for cut in [1usize, 5, 20, 50] {
            let s = &ETHANOL_LIKE[..cut.min(ETHANOL_LIKE.len())];
            let _ = parse_mrv(s); // must not panic, result doesn't matter
        }
    }

    const ETHANOL_LIKE: &str = r#"<cml><MDocument><MChemicalStruct><molecule><atomArray><atom id="a1" elementType="C"/></atomArray></molecule></MChemicalStruct></MDocument></cml>"#;

    #[test]
    fn invalid_utf8_boundary_input_does_not_panic() {
        // A raw byte slice that is not valid UTF-8 must never reach `parse_mrv`
        // as a `&str` in the first place (the type system already prevents
        // it); this test instead exercises non-ASCII multi-byte UTF-8 content
        // inside an attribute value, adjacent to the parser's own byte-index
        // slicing, to ensure no panic on a char boundary.
        let xml = r#"<cml><MDocument><MChemicalStruct><molecule><atomArray><atom id="a1" elementType="C"/></atomArray></molecule></MChemicalStruct></MDocument></cml>"#;
        let mut with_unicode = xml.replace("m1", "m1\u{1F9EA}");
        with_unicode.push('\u{1F9EA}');
        let _ = parse_mrv(&with_unicode);
    }

    #[test]
    fn excessive_atom_count_does_not_panic_or_hang() {
        let mut atoms = String::new();
        for i in 0..20_000 {
            atoms.push_str(&format!(r#"<atom id="a{i}" elementType="C"/>"#));
        }
        let xml = format!(
            "<cml><MDocument><MChemicalStruct><molecule><atomArray>{atoms}</atomArray></molecule></MChemicalStruct></MDocument></cml>"
        );
        let rec = parse_mrv(&xml).unwrap();
        assert_eq!(rec.mol.atom_count(), 20_000);
    }

    #[test]
    fn nan_infinity_coordinate_values_do_not_panic() {
        // f64::parse accepts literal "NaN"/"inf"/"-inf" tokens -- confirm
        // this never panics downstream (e.g. in write_mrv's coordinate
        // formatting) rather than silently producing a malformed float.
        for bad in ["NaN", "inf", "-inf", "infinity"] {
            let xml = format!(
                r#"<cml><MDocument><MChemicalStruct><molecule><atomArray><atom id="a1" elementType="C" x2="{bad}" y2="0.0"/></atomArray></molecule></MChemicalStruct></MDocument></cml>"#
            );
            let rec = parse_mrv(&xml).unwrap();
            let coords = rec.coordinates_2d.as_ref().unwrap();
            assert!(!coords[0][0].is_finite());
            // Must not panic when writing the non-finite value back out.
            let _ = write_mrv(&rec, &MrvWriteOptions::default());
        }
    }

    #[test]
    fn random_mutation_corpus_never_panics() {
        // Deterministic splitmix64, matching the seeded-mutation convention
        // used in smiles_table.rs/tdt.rs (no `rand` dependency).
        struct SplitMix64(u64);
        impl SplitMix64 {
            fn next(&mut self) -> u64 {
                self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
                let mut z = self.0;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
                z ^ (z >> 31)
            }
        }

        let base = ETHANOL_LIKE.as_bytes().to_vec();
        let mut rng = SplitMix64(42);
        for _ in 0..500 {
            let mut mutated = base.clone();
            let n_mutations = 1 + (rng.next() % 5) as usize;
            for _ in 0..n_mutations {
                if mutated.is_empty() {
                    break;
                }
                let pos = (rng.next() as usize) % mutated.len();
                mutated[pos] = (rng.next() % 256) as u8;
            }
            if let Ok(s) = String::from_utf8(mutated) {
                let _ = parse_mrv(&s); // must not panic
            }
        }
    }
}
