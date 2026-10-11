//! ChemDraw XML (CDXML) parser and writer.
//!
//! CDXML is a proprietary XML format produced by ChemDraw (PerkinElmer /
//! Revvity).  This module handles the minimal subset needed to read and
//! write molecular structure through a CDXML document.
//!
//! # Supported elements / attributes
//!
//! `<n>` (atom): `id`, `Element` (atomic number), `p` ("x y" 2D coords),
//!               `NumHydrogens`, `Charge`, `Isotope`, and bounded
//!               `GenericNickname` pseudoatoms (`*`, `R`, `R<n>`)
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
//! - Tetrahedral chirality (`Atom.chirality`) is perceived from wedge/hash
//!   `Display` bonds plus 2D coordinates via
//!   [`chematic_perception::apply_local_parity_from_wedges`] -- the same
//!   CIP-independent mechanism the MOL/MRV readers use, so it abstains
//!   (leaves `Chirality::None`) rather than guessing on contradictory
//!   wedges, missing coordinates, or degenerate/coplanar geometry.
//!   Directional wedges (`WedgeBegin`/`WedgeEnd`/`WedgedHashBegin`/
//!   `WedgedHashEnd`) are always perceived. Non-directional ones
//!   (`Bold`/`Hash`/`Dash` -- a plain thick or dashed line, no narrow/wide
//!   end) are perceived only when [`CdxmlParseOptions::infer_nondirectional_stereo`]
//!   is explicitly opted in (`parse_cdxml_with_options`/
//!   `parse_cdxml_all_with_options`), since `Bold` in particular is
//!   sometimes drawn for visual emphasis rather than stereo intent. When
//!   opted in, a non-directional bond whose BOTH endpoints independently
//!   qualify as stereocenter candidates (3-4 neighbours) still abstains --
//!   one such bond cannot mean "toward the viewer" from both ends at once
//!   -- while a directional wedge in the same situation is not ambiguous
//!   (the Begin atom is the reference by CDXML's own convention) and is
//!   unaffected either way.
//! - E/Z double-bond stereo is derived from 2D coordinates via `assign_ez_from_2d`.
//! - Presentation-only nodes (text boxes, arrows, etc.) are silently skipped.

use std::collections::HashMap;

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder, RGroupLabel};
use chematic_perception::{apply_local_parity_from_wedges, assign_ez_from_2d};

use crate::cml::parse_xml_attrs;

/// Resource limits for CDXML parsing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CdxmlParseLimits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_attribute_bytes: usize,
    pub max_atoms: usize,
    pub max_bonds: usize,
    pub max_fragments: usize,
}

impl Default for CdxmlParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_line_bytes: 16 * 1024 * 1024,
            max_lines: 5_000_000,
            max_attribute_bytes: 16 * 1024 * 1024,
            max_atoms: CDXML_MAX_ATOMS,
            max_bonds: 20_000,
            max_fragments: 10_000,
        }
    }
}

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Error returned when parsing a CDXML document fails.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CdxmlError {
    /// The input did not contain one well-formed CDXML document root.
    InvalidDocument(String),
    /// An atom `Element` attribute contained an unknown atomic number.
    UnknownAtomicNumber(u32),
    /// A bond referenced an atom id that was not defined.
    UnknownAtomRef(String),
    /// An edit referenced a page ID that occurs more than once.
    AmbiguousPageId(String),
    /// A `<b>` bond element is missing a `B` or `E` attribute.
    MissingBondEndpoint,
    /// The `p` coordinate attribute could not be parsed.
    InvalidCoords(String),
    /// A generic nickname is outside the lossless `*`/`R`/`R<n>` contract.
    UnsupportedPseudoatomLabel(String),
    /// The document contains more atoms than the parser's safety limit.
    TooManyAtoms(usize),
    /// The document exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl std::fmt::Display for CdxmlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CdxmlError::InvalidDocument(s) => write!(f, "invalid CDXML document: {s}"),
            CdxmlError::UnknownAtomicNumber(n) => write!(f, "unknown atomic number: {n}"),
            CdxmlError::UnknownAtomRef(s) => write!(f, "unknown atom ref: {s}"),
            CdxmlError::AmbiguousPageId(s) => write!(f, "ambiguous page id: {s}"),
            CdxmlError::MissingBondEndpoint => write!(f, "bond missing B or E attribute"),
            CdxmlError::InvalidCoords(s) => write!(f, "invalid p coords: {s}"),
            CdxmlError::UnsupportedPseudoatomLabel(label) => {
                write!(f, "unsupported CDXML pseudoatom label: {label}")
            }
            CdxmlError::TooManyAtoms(n) => write!(f, "CDXML document exceeds atom limit ({n})"),
            CdxmlError::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "CDXML: {resource} has size {actual}, exceeding limit {limit}"
            ),
        }
    }
}

impl std::error::Error for CdxmlError {}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Options controlling optional/heuristic CDXML parsing behavior.
///
/// Currently just the one knob (kept as a struct, not a bare `bool`
/// parameter, so a future option doesn't force every call site to change
/// signature again).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CdxmlParseOptions {
    /// When `true`, also perceive tetrahedral chirality from
    /// **non-directional** wedge displays (`Bold`/`Hash`/`Dash` -- a plain
    /// thick or dashed line, no narrow/wide end). Default **`false`**:
    /// `Bold` in particular is sometimes used by chemists as a plain thick
    /// line for visual emphasis, not stereo intent, so treating it as
    /// stereo unconditionally would manufacture false positives (RDKit
    /// issue #9359's own discussion leans the same way -- inferring stereo
    /// from these should be opt-in, not silently default).
    ///
    /// This flag has no effect on **directional** wedges
    /// (`WedgeBegin`/`WedgeEnd`/`WedgedHashBegin`/`WedgedHashEnd`), which
    /// are unambiguous stereo intent and are always perceived regardless.
    pub infer_nondirectional_stereo: bool,
}

/// Parse a CDXML document and return the first molecular fragment.
///
/// Convenience wrapper around [`parse_cdxml_all`], using
/// [`CdxmlParseOptions::default`]. Use [`parse_cdxml_with_options`] to
/// perceive chirality from non-directional (`Bold`/`Hash`/`Dash`) wedges
/// too.
///
/// **Coordinate system:** The returned `coords` use **ChemDraw Y-down convention**
/// (Y increases downward, matching screen/SVG pixel space). No Y-axis conversion is required
/// for SVG rendering; coordinates can be used directly with an SVG renderer such as
/// `render_svg` (see `chematic-depict`).
pub fn parse_cdxml(input: &str) -> Result<(Molecule, Vec<(f64, f64)>), CdxmlError> {
    parse_cdxml_with_options_and_limits(
        input,
        &CdxmlParseOptions::default(),
        &CdxmlParseLimits::default(),
    )
}

/// Parse a CDXML document with explicit resource limits.
pub fn parse_cdxml_with_limits(
    input: &str,
    limits: &CdxmlParseLimits,
) -> Result<(Molecule, Vec<(f64, f64)>), CdxmlError> {
    parse_cdxml_with_options_and_limits(input, &CdxmlParseOptions::default(), limits)
}

/// Same as [`parse_cdxml`], with explicit [`CdxmlParseOptions`].
pub fn parse_cdxml_with_options(
    input: &str,
    options: &CdxmlParseOptions,
) -> Result<(Molecule, Vec<(f64, f64)>), CdxmlError> {
    parse_cdxml_with_options_and_limits(input, options, &CdxmlParseLimits::default())
}

/// Same as [`parse_cdxml_with_options`], with explicit resource limits.
pub fn parse_cdxml_with_options_and_limits(
    input: &str,
    options: &CdxmlParseOptions,
    limits: &CdxmlParseLimits,
) -> Result<(Molecule, Vec<(f64, f64)>), CdxmlError> {
    let mut all = parse_cdxml_all_with_options_and_limits(input, options, limits)?;
    if all.is_empty() {
        return Ok((MoleculeBuilder::new().build(), vec![]));
    }
    Ok(all.remove(0))
}

/// Parse all molecular fragments from a CDXML document, using
/// [`CdxmlParseOptions::default`]. Use [`parse_cdxml_all_with_options`] to
/// perceive chirality from non-directional (`Bold`/`Hash`/`Dash`) wedges
/// too.
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
/// wedge), mapped to the same `BondOrder` bucket as `"WedgeBegin"` for
/// round-trip/structural purposes. `"Bold"`/`"Hash"`/`"Dash"` have no
/// Begin/End reference convention, unlike `"WedgeBegin"`/`"WedgeEnd"`.
/// When non-directional stereo inference is enabled
/// ([`CdxmlParseOptions::infer_nondirectional_stereo`]), the parser
/// identifies a unique stereocenter-candidate endpoint and normalizes a
/// temporary perception view so the result is invariant to CDXML B/E
/// ordering; if the reference endpoint is ambiguous (or there isn't one),
/// chirality perception abstains rather than guessing.
///
/// These wedge bonds, combined with 2D coordinates, are what actually
/// perceives `Atom.chirality` (see the module-level doc for the full
/// abstain-on-ambiguity behavior, and [`CdxmlParseOptions`] for the
/// non-directional opt-in).
#[allow(clippy::type_complexity)]
pub fn parse_cdxml_all(input: &str) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError> {
    parse_cdxml_all_with_options_and_limits(
        input,
        &CdxmlParseOptions::default(),
        &CdxmlParseLimits::default(),
    )
}

/// Parse all CDXML fragments with explicit resource limits.
#[allow(clippy::type_complexity)]
pub fn parse_cdxml_all_with_limits(
    input: &str,
    limits: &CdxmlParseLimits,
) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError> {
    parse_cdxml_all_with_options_and_limits(input, &CdxmlParseOptions::default(), limits)
}

/// Accumulator for a single CDXML `<fragment>` being parsed.
#[derive(Default)]
struct FragAccum {
    atom_ids: Vec<String>,
    atom_elems: Vec<Element>,
    atom_charges: Vec<i8>,
    atom_isotopes: Vec<Option<u16>>,
    atom_h: Vec<Option<u8>>,
    atom_wildcards: Vec<bool>,
    atom_rgroups: Vec<Option<RGroupLabel>>,
    atom_xs: Vec<f64>,
    atom_ys: Vec<f64>,
    bond_bs: Vec<String>,
    bond_es: Vec<String>,
    bond_ords: Vec<BondOrder>,
    /// Raw `Display` attribute value per bond, `None` if absent. Needed
    /// alongside `bond_ords` because `BondOrder::Up`/`Down` alone can't
    /// distinguish a directional wedge (`WedgeBegin`/`WedgeEnd`, which has
    /// an inherent narrow-end convention) from a non-directional one
    /// (`Bold`/`Hash`/`Dash`, which doesn't) -- see the ambiguity check in
    /// `flush()`.
    bond_display: Vec<Option<String>>,
}

struct ResolvedBonds {
    positions: Vec<(usize, usize)>,
    degree: Vec<usize>,
}

struct BuiltFragment {
    molecule: Molecule,
    indices: HashMap<usize, AtomIdx>,
    coords: Vec<(f64, f64)>,
}

impl FragAccum {
    fn is_empty(&self) -> bool {
        self.atom_ids.is_empty() && self.bond_bs.is_empty()
    }

    fn resolve_bond_positions(&self) -> Result<ResolvedBonds, CdxmlError> {
        let positions = self
            .atom_ids
            .iter()
            .enumerate()
            .map(|(index, id)| (id.as_str(), index))
            .collect::<HashMap<_, _>>();
        let mut bonds = Vec::with_capacity(self.bond_bs.len());
        let mut degree = vec![0; self.atom_ids.len()];
        for (begin, end) in self.bond_bs.iter().zip(&self.bond_es) {
            let begin_pos = *positions
                .get(begin.as_str())
                .ok_or_else(|| CdxmlError::UnknownAtomRef(begin.clone()))?;
            let end_pos = *positions
                .get(end.as_str())
                .ok_or_else(|| CdxmlError::UnknownAtomRef(end.clone()))?;
            degree[begin_pos] += 1;
            degree[end_pos] += 1;
            bonds.push((begin_pos, end_pos));
        }
        Ok(ResolvedBonds {
            positions: bonds,
            degree,
        })
    }

    fn build_molecule(
        &self,
        bond_positions: &[(usize, usize)],
    ) -> Result<BuiltFragment, CdxmlError> {
        let mut builder = MoleculeBuilder::new();
        let mut indices = HashMap::new();
        let mut coords = Vec::with_capacity(self.atom_ids.len());
        for index in 0..self.atom_ids.len() {
            let mut atom = Atom::new(self.atom_elems[index]);
            atom.wildcard = self.atom_wildcards[index];
            atom.charge = self.atom_charges[index];
            atom.isotope = self.atom_isotopes[index];
            atom.hydrogen_count = if atom.wildcard {
                Some(0)
            } else {
                self.atom_h[index]
            };
            let atom_index = builder.add_atom(atom);
            if let Some(label) = self.atom_rgroups[index] {
                builder.set_r_group(atom_index, label);
            }
            indices.insert(index, atom_index);
            coords.push((self.atom_xs[index], self.atom_ys[index]));
        }
        for (bond_index, &(begin, end)) in bond_positions.iter().enumerate() {
            builder
                .add_bond(indices[&begin], indices[&end], self.bond_ords[bond_index])
                .map_err(|_| {
                    CdxmlError::UnknownAtomRef(format!(
                        "{} {}",
                        self.bond_bs[bond_index], self.bond_es[bond_index]
                    ))
                })?;
        }
        Ok(BuiltFragment {
            molecule: builder.build(),
            indices,
            coords,
        })
    }

    fn has_nondirectional_wedge(&self) -> bool {
        self.bond_display.iter().any(|display| {
            matches!(
                display.as_deref(),
                Some("Bold") | Some("Hash") | Some("Dash")
            )
        })
    }

    fn normalize_nondirectional_wedges(
        &self,
        view: &mut Molecule,
        indices: &HashMap<usize, AtomIdx>,
        bond_positions: &[(usize, usize)],
        degree: &[usize],
        options: &CdxmlParseOptions,
    ) {
        for (bond_index, &(begin, end)) in bond_positions.iter().enumerate() {
            if !matches!(
                self.bond_display[bond_index].as_deref(),
                Some("Bold") | Some("Hash") | Some("Dash")
            ) {
                continue;
            }
            let (begin_atom, end_atom) = (indices[&begin], indices[&end]);
            let Some((index, _)) = view.bond_between(begin_atom, end_atom) else {
                continue;
            };
            if !options.infer_nondirectional_stereo {
                view.set_bond_order(index, BondOrder::Single);
                continue;
            }
            let eligible = |value: usize| (3..=4).contains(&value);
            let (begin_eligible, end_eligible) = (eligible(degree[begin]), eligible(degree[end]));
            if begin_eligible && end_eligible {
                view.set_bond_order(index, BondOrder::Single);
            } else if end_eligible && !begin_eligible {
                let flipped = match view.bond(index).order {
                    BondOrder::Up => BondOrder::Down,
                    BondOrder::Down => BondOrder::Up,
                    other => other,
                };
                view.set_bond_order(index, flipped);
            }
        }
    }

    fn apply_stereo(
        &self,
        molecule: &mut Molecule,
        indices: &HashMap<usize, AtomIdx>,
        bond_positions: &[(usize, usize)],
        degree: &[usize],
        coords: &[(f64, f64)],
        options: &CdxmlParseOptions,
    ) {
        if !self.has_nondirectional_wedge() {
            apply_local_parity_from_wedges(molecule, coords);
            return;
        }
        let mut view = molecule.clone();
        self.normalize_nondirectional_wedges(&mut view, indices, bond_positions, degree, options);
        apply_local_parity_from_wedges(&mut view, coords);
        for index in 0..self.atom_ids.len() {
            let atom = indices[&index];
            let chirality = view.atom(atom).chirality;
            if chirality != chematic_core::Chirality::None {
                molecule.set_chirality(atom, chirality);
                if let Some(order) = view.stereo_neighbor_order(atom) {
                    molecule.set_stereo_neighbor_order(atom, order.to_vec());
                }
            }
        }
    }

    fn flush(
        &mut self,
        results: &mut Vec<(Molecule, Vec<(f64, f64)>)>,
        options: &CdxmlParseOptions,
    ) -> Result<(), CdxmlError> {
        if self.is_empty() {
            return Ok(());
        }
        let bonds = self.resolve_bond_positions()?;
        let mut fragment = self.build_molecule(&bonds.positions)?;
        self.apply_stereo(
            &mut fragment.molecule,
            &fragment.indices,
            &bonds.positions,
            &bonds.degree,
            &fragment.coords,
            options,
        );
        assign_ez_from_2d(&mut fragment.molecule, &fragment.coords);
        results.push((fragment.molecule, fragment.coords));
        *self = Self::default();
        Ok(())
    }
}

/// Maximum atoms allowed in a single CDXML fragment (DoS guard).
pub const CDXML_MAX_ATOMS: usize = 10_000;

/// Same as [`parse_cdxml_all`], with explicit [`CdxmlParseOptions`].
#[allow(clippy::type_complexity)]
pub fn parse_cdxml_all_with_options(
    input: &str,
    options: &CdxmlParseOptions,
) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError> {
    parse_cdxml_all_with_options_and_limits(input, options, &CdxmlParseLimits::default())
}

fn check_cdxml_input_limits(input: &str, limits: &CdxmlParseLimits) -> Result<(), CdxmlError> {
    if input.len() > limits.max_input_bytes {
        return Err(CdxmlError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    for (line_index, line) in input.lines().enumerate() {
        if line_index >= limits.max_lines {
            return Err(CdxmlError::ResourceLimit {
                resource: "lines",
                actual: line_index + 1,
                limit: limits.max_lines,
            });
        }
        if line.len() > limits.max_line_bytes {
            return Err(CdxmlError::ResourceLimit {
                resource: "line bytes",
                actual: line.len(),
                limit: limits.max_line_bytes,
            });
        }
    }
    Ok(())
}

fn cdxml_atom_kind(
    attrs: &HashMap<String, String>,
) -> Result<(Element, bool, Option<RGroupLabel>), CdxmlError> {
    let pseudo = if attrs.get("NodeType").map(String::as_str) == Some("GenericNickname") {
        Some(attrs.get("GenericNickname").cloned().ok_or_else(|| {
            CdxmlError::UnsupportedPseudoatomLabel("missing GenericNickname attribute".into())
        })?)
    } else {
        None
    };
    match pseudo.as_deref() {
        Some("*") => Ok((Element::C, true, None)),
        Some("R") => Ok((Element::C, true, Some(RGroupLabel::unnumbered()))),
        Some(label) if label.starts_with('R') => label[1..]
            .parse::<u16>()
            .ok()
            .and_then(RGroupLabel::numbered)
            .map(|label| (Element::C, true, Some(label)))
            .ok_or_else(|| CdxmlError::UnsupportedPseudoatomLabel(label.into())),
        Some(label) => Err(CdxmlError::UnsupportedPseudoatomLabel(label.into())),
        None => {
            let number = attrs
                .get("Element")
                .and_then(|value| value.trim().parse::<u32>().ok())
                .unwrap_or(6);
            if number > 255 {
                return Err(CdxmlError::UnknownAtomicNumber(number));
            }
            Element::from_atomic_number(number as u8)
                .map(|element| (element, false, None))
                .ok_or(CdxmlError::UnknownAtomicNumber(number))
        }
    }
}

fn cdxml_coords(attrs: &HashMap<String, String>) -> Result<(f64, f64), CdxmlError> {
    let Some(value) = attrs.get("p") else {
        return Ok((0.0, 0.0));
    };
    let fields = value.split_whitespace().collect::<Vec<_>>();
    if fields.len() < 2 {
        return Err(CdxmlError::InvalidCoords(value.clone()));
    }
    Ok((
        fields[0].parse().unwrap_or(0.0),
        fields[1].parse().unwrap_or(0.0),
    ))
}

fn push_cdxml_atom(
    accumulator: &mut FragAccum,
    line: &str,
    limits: &CdxmlParseLimits,
) -> Result<(), CdxmlError> {
    let attrs = parse_xml_attrs(line);
    check_attribute_limits(&attrs, limits)?;
    let Some(id) = attrs.get("id").cloned() else {
        return Ok(());
    };
    if accumulator.atom_ids.len() >= limits.max_atoms {
        return Err(CdxmlError::ResourceLimit {
            resource: "atoms",
            actual: accumulator.atom_ids.len() + 1,
            limit: limits.max_atoms,
        });
    }
    let (element, wildcard, r_group) = cdxml_atom_kind(&attrs)?;
    let (x, y) = cdxml_coords(&attrs)?;
    accumulator.atom_ids.push(id);
    accumulator.atom_elems.push(element);
    accumulator.atom_charges.push(
        attrs
            .get("Charge")
            .and_then(|value| value.trim().parse().ok())
            .unwrap_or(0),
    );
    accumulator.atom_isotopes.push(
        attrs
            .get("Isotope")
            .and_then(|value| value.trim().parse::<u16>().ok())
            .filter(|&value| value > 0),
    );
    accumulator.atom_h.push(
        attrs
            .get("NumHydrogens")
            .and_then(|value| value.trim().parse().ok()),
    );
    accumulator.atom_wildcards.push(wildcard);
    accumulator.atom_rgroups.push(r_group);
    accumulator.atom_xs.push(x);
    accumulator.atom_ys.push(y);
    Ok(())
}

fn cdxml_bond_order(attrs: &HashMap<String, String>) -> BondOrder {
    let base = match attrs.get("Order").map(String::as_str) {
        Some("2") => BondOrder::Double,
        Some("3") => BondOrder::Triple,
        Some("1.5") => BondOrder::Aromatic,
        _ => BondOrder::Single,
    };
    if base != BondOrder::Single {
        return base;
    }
    match attrs.get("Display").map(String::as_str) {
        Some("WedgeBegin") | Some("WedgedHashBegin") | Some("Bold") => BondOrder::Up,
        Some("Hash") | Some("Dash") | Some("WedgeEnd") | Some("WedgedHashEnd") => BondOrder::Down,
        _ => BondOrder::Single,
    }
}

fn push_cdxml_bond(
    accumulator: &mut FragAccum,
    line: &str,
    limits: &CdxmlParseLimits,
) -> Result<(), CdxmlError> {
    let attrs = parse_xml_attrs(line);
    check_attribute_limits(&attrs, limits)?;
    if accumulator.bond_bs.len() >= limits.max_bonds {
        return Err(CdxmlError::ResourceLimit {
            resource: "bonds",
            actual: accumulator.bond_bs.len() + 1,
            limit: limits.max_bonds,
        });
    }
    accumulator.bond_bs.push(
        attrs
            .get("B")
            .cloned()
            .ok_or(CdxmlError::MissingBondEndpoint)?,
    );
    accumulator.bond_es.push(
        attrs
            .get("E")
            .cloned()
            .ok_or(CdxmlError::MissingBondEndpoint)?,
    );
    accumulator.bond_ords.push(cdxml_bond_order(&attrs));
    accumulator.bond_display.push(attrs.get("Display").cloned());
    Ok(())
}

/// Same as [`parse_cdxml_all_with_options`], with explicit resource limits.
#[allow(clippy::type_complexity)]
pub fn parse_cdxml_all_with_options_and_limits(
    input: &str,
    options: &CdxmlParseOptions,
    limits: &CdxmlParseLimits,
) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError> {
    check_cdxml_input_limits(input, limits)?;
    let mut accumulator = FragAccum::default();
    let mut results = Vec::new();
    let mut fragment_count = 0;
    for raw in crate::cdxml_document::logical_cdxml_lines(input) {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with("<fragment") {
            if fragment_count >= limits.max_fragments {
                return Err(CdxmlError::ResourceLimit {
                    resource: "fragments",
                    actual: fragment_count + 1,
                    limit: limits.max_fragments,
                });
            }
            fragment_count += 1;
            accumulator = FragAccum::default();
        } else if line.starts_with("</fragment>") {
            accumulator.flush(&mut results, options)?;
        } else if is_n_tag(line) {
            push_cdxml_atom(&mut accumulator, line, limits)?;
        } else if is_b_tag(line) {
            push_cdxml_bond(&mut accumulator, line, limits)?;
        }
    }
    accumulator.flush(&mut results, options)?;
    Ok(results)
}

fn check_attribute_limits(
    attrs: &HashMap<String, String>,
    limits: &CdxmlParseLimits,
) -> Result<(), CdxmlError> {
    for (key, value) in attrs {
        let bytes = key.len().saturating_add(value.len());
        if bytes > limits.max_attribute_bytes {
            return Err(CdxmlError::ResourceLimit {
                resource: "attribute bytes",
                actual: bytes,
                limit: limits.max_attribute_bytes,
            });
        }
    }
    Ok(())
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
        let mut parts = vec![format!("id=\"{id}\""), format!("p=\"{x} {y}\"")];
        if atom.wildcard {
            let label = match mol.r_group_label(idx).and_then(|label| label.number()) {
                Some(number) => format!("R{number}"),
                None if mol.r_group_label(idx).is_some() => "R".to_string(),
                None => "*".to_string(),
            };
            parts.push("NodeType=\"GenericNickname\"".to_string());
            parts.push(format!("GenericNickname=\"{label}\""));
        } else {
            parts.push(format!("Element=\"{}\"", atom.element.atomic_number()));
        }
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
#[path = "cdxml/tests.rs"]
mod tests;
