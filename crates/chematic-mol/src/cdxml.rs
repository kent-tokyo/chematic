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
//! `<n NodeType="Nickname">`/`<n NodeType="Fragment">`: ChemDraw
//!               abbreviations (Boc, Ph, template structures) are expanded
//!               from their nested `<fragment>`, attached through its
//!               `ExternalConnectionPoint` nodes (and the outer node's
//!               `BondOrdering` when there are several)
//!
//! `<b>` (bond): `B` (begin atom id), `E` (end atom id), `Order` (1/2/3,
//!               defaults to 1; also 4, 1.5/aromatic, directed dative, and
//!               the query orders `any`, `1 1.5`, `2 1.5`, `1 2`)
//! A dative bond retains B -> E regardless of wedge/hash display decoration.
//! Other declared orders return an [`CdxmlError::InvalidDocument`] classified
//! by [`CdxmlError::unsupported_bond_order`]; only an absent order defaults to
//! single. This is a bounded graph profile. [`parse_cdxml_all_with_report`]
//! keeps the rest of a document past generic query labels and
//! unrepresentable bond orders, reporting each as a diagnostic.
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
//!   [`chematic_perception::apply_ctab_local_parity_from_wedges_with_diagnostics`] -- the same
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
//!   (the narrow endpoint is the reference by CDXML's convention) and is
//!   unaffected either way.
//! - E/Z double-bond stereo is derived from 2D coordinates: bond directions
//!   via `apply_ez_directions_from_2d` (what the SMILES and MOL writers read,
//!   as for MOL and MRV input), and the E/Z label cache via `assign_ez_from_2d`.
//! - Presentation-only nodes (text boxes, arrows, etc.) are silently skipped.

use std::collections::HashMap;

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder, RGroupLabel};
use chematic_perception::{
    apply_ctab_local_parity_from_wedges_with_diagnostics, apply_ez_directions_from_2d_ex,
    assign_ez_from_2d,
};

use crate::cml::{parse_xml_attrs, unescape_xml};

#[path = "cdxml/abbreviation.rs"]
mod abbreviation;
pub use abbreviation::{CdxmlAbbreviation, CdxmlParseDiagnostic, CdxmlParseReport};

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
            CdxmlError::InvalidDocument(s)
                if s.starts_with(UNSUPPORTED_BOND_ORDER_PREFIX)
                    || CdxmlError::is_unexpandable_abbreviation(s) =>
            {
                f.write_str(s)
            }
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

const UNSUPPORTED_BOND_ORDER_PREFIX: &str = "unsupported CDXML bond order: ";

impl CdxmlError {
    /// Returns the authored bond-order token when this error is the bounded
    /// CDXML profile's unsupported-order refusal.
    ///
    /// This classifier keeps the public error enum source-compatible with
    /// earlier 1.x releases while allowing callers to distinguish this refusal
    /// without parsing [`std::fmt::Display`] output.
    pub fn unsupported_bond_order(&self) -> Option<&str> {
        match self {
            Self::InvalidDocument(message) => message.strip_prefix(UNSUPPORTED_BOND_ORDER_PREFIX),
            _ => None,
        }
    }
}

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
/// Solid wedges (`WedgeBegin`/`WedgeEnd`) and `Bold` use [`BondOrder::Up`];
/// hashed wedges (`WedgedHashBegin`/`WedgedHashEnd`) and `Hash`/`Dash` use
/// [`BondOrder::Down`]. End glyphs reverse the bond's reference endpoint.
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

/// One parsed `<n>` atom node.
#[derive(Clone, Debug)]
struct AtomRec {
    id: String,
    element: Element,
    charge: i8,
    isotope: Option<u16>,
    hydrogens: Option<u8>,
    wildcard: bool,
    r_group: Option<RGroupLabel>,
    x: f64,
    y: f64,
    /// `NodeType="Nickname"`/`"Fragment"`: an abbreviation placeholder that
    /// a nested `<fragment>` normally replaces (see `abbreviation.rs`).
    abbreviation: bool,
}

/// One parsed `<b>` bond.
#[derive(Clone, Debug)]
struct BondRec {
    id: Option<String>,
    begin: String,
    end: String,
    order: BondOrder,
    /// Raw `Display` attribute value, `None` if absent. Needed alongside
    /// `order` because `BondOrder::Up`/`Down` alone can't distinguish a
    /// directional wedge (`WedgeBegin`/`WedgeEnd`, which has an inherent
    /// narrow-end convention) from a non-directional one
    /// (`Bold`/`Hash`/`Dash`, which doesn't) -- see the ambiguity check in
    /// `normalize_nondirectional_wedges()`.
    display: Option<String>,
}

/// An `<n>` element whose children are still being read.
#[derive(Debug)]
struct OpenNode {
    id: Option<String>,
    abbreviation: bool,
    label: String,
    in_string: bool,
}

/// Accumulator for a single CDXML `<fragment>` being parsed.
#[derive(Default, Debug)]
struct FragAccum {
    atoms: Vec<AtomRec>,
    bonds: Vec<BondRec>,
    /// `NodeType="ExternalConnectionPoint"` nodes: `(id, ExternalConnectionNum)`.
    connection_points: Vec<(String, Option<u32>)>,
    /// `BondOrdering` per node id (bond ids, one per external connection).
    bond_orderings: HashMap<String, Vec<String>>,
    /// Nested fragments read inside abbreviation nodes, spliced at flush.
    expansions: Vec<abbreviation::Expansion>,
    /// Expanded abbreviations: `(node id, label, atom ids)`.
    abbreviations: Vec<(String, Option<String>, Vec<String>)>,
    open_node: Option<OpenNode>,
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

/// Parse-wide state: tolerance mode and the diagnostics it collects.
#[derive(Default)]
struct ParseCtx {
    tolerant: bool,
    diagnostics: Vec<CdxmlParseDiagnostic>,
    abbreviations: Vec<CdxmlAbbreviation>,
}

impl ParseCtx {
    fn diagnose(
        &mut self,
        molecule_index: usize,
        code: &'static str,
        object_id: Option<&str>,
        message: String,
    ) {
        self.diagnostics.push(CdxmlParseDiagnostic {
            molecule_index,
            code,
            object_id: object_id.map(str::to_string),
            message,
        });
    }
}

impl FragAccum {
    fn is_empty(&self) -> bool {
        self.atoms.is_empty() && self.bonds.is_empty()
    }

    fn resolve_bond_positions(&self) -> Result<ResolvedBonds, CdxmlError> {
        let positions = self
            .atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| (atom.id.as_str(), index))
            .collect::<HashMap<_, _>>();
        let mut bonds = Vec::with_capacity(self.bonds.len());
        let mut degree = vec![0; self.atoms.len()];
        for bond in &self.bonds {
            let begin_pos = *positions
                .get(bond.begin.as_str())
                .ok_or_else(|| CdxmlError::UnknownAtomRef(bond.begin.clone()))?;
            let end_pos = *positions
                .get(bond.end.as_str())
                .ok_or_else(|| CdxmlError::UnknownAtomRef(bond.end.clone()))?;
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
        let mut coords = Vec::with_capacity(self.atoms.len());
        for (index, rec) in self.atoms.iter().enumerate() {
            let mut atom = Atom::new(rec.element);
            atom.wildcard = rec.wildcard;
            atom.charge = rec.charge;
            atom.isotope = rec.isotope;
            atom.hydrogen_count = if atom.wildcard {
                Some(0)
            } else {
                rec.hydrogens
            };
            let atom_index = builder.add_atom(atom);
            if let Some(label) = rec.r_group {
                builder.set_r_group(atom_index, label);
            }
            indices.insert(index, atom_index);
            coords.push((rec.x, rec.y));
        }
        for (bond, &(begin, end)) in self.bonds.iter().zip(bond_positions) {
            // Single-bond End glyphs put the stereo reference end on E.
            // They never reverse the independent B -> E dative direction.
            let (begin, end) = if matches!(bond.order, BondOrder::Up | BondOrder::Down)
                && matches!(
                    bond.display.as_deref(),
                    Some("WedgeEnd") | Some("WedgedHashEnd")
                ) {
                (end, begin)
            } else {
                (begin, end)
            };
            builder
                .add_bond(indices[&begin], indices[&end], bond.order)
                .map_err(|_| CdxmlError::UnknownAtomRef(format!("{} {}", bond.begin, bond.end)))?;
        }
        Ok(BuiltFragment {
            molecule: builder.build(),
            indices,
            coords,
        })
    }

    fn has_nondirectional_wedge(&self) -> bool {
        self.bonds
            .iter()
            .any(|bond| is_nondirectional_display(bond.display.as_deref()))
    }

    fn normalize_nondirectional_wedges(
        &self,
        molecule: &Molecule,
        bond_positions: &[(usize, usize)],
        degree: &[usize],
        options: &CdxmlParseOptions,
    ) -> Result<Molecule, CdxmlError> {
        let mut builder =
            MoleculeBuilder::with_capacity(molecule.atom_count(), molecule.bond_count());
        for (_, atom) in molecule.atoms() {
            builder.add_atom(atom.clone());
        }
        for (bond_index, (_, bond)) in molecule.bonds().enumerate() {
            let (mut begin, mut end, mut order) = (bond.atom1, bond.atom2, bond.order);
            if is_nondirectional_display(self.bonds[bond_index].display.as_deref())
                && matches!(order, BondOrder::Up | BondOrder::Down)
            {
                let (b, e) = bond_positions[bond_index];
                let eligible = |value: usize| (3..=4).contains(&value);
                let (b_ok, e_ok) = (eligible(degree[b]), eligible(degree[e]));
                if !options.infer_nondirectional_stereo || b_ok == e_ok {
                    order = BondOrder::Single;
                } else if e_ok {
                    (begin, end) = (end, begin);
                }
            }
            builder.add_bond(begin, end, order).map_err(|error| {
                CdxmlError::InvalidDocument(format!("invalid perception bond: {error}"))
            })?;
        }
        Ok(builder.build())
    }

    fn apply_stereo(
        &self,
        molecule: &mut Molecule,
        indices: &HashMap<usize, AtomIdx>,
        bonds: &ResolvedBonds,
        coords: &[(f64, f64)],
        options: &CdxmlParseOptions,
    ) -> Result<(), CdxmlError> {
        // CDXML uses screen Y-down coordinates; perception expects Cartesian Y-up.
        // Keep the returned source coordinates unchanged.
        let cartesian: Vec<_> = coords.iter().map(|&(x, y)| (x, -y)).collect();
        if !self.has_nondirectional_wedge() {
            apply_ctab_local_parity_from_wedges_with_diagnostics(molecule, &cartesian);
            return Ok(());
        }
        let mut view = self.normalize_nondirectional_wedges(
            molecule,
            &bonds.positions,
            &bonds.degree,
            options,
        )?;
        apply_ctab_local_parity_from_wedges_with_diagnostics(&mut view, &cartesian);
        for index in 0..self.atoms.len() {
            let atom = indices[&index];
            let chirality = view.atom(atom).chirality;
            if chirality != chematic_core::Chirality::None {
                molecule.set_chirality(atom, chirality);
                if let Some(order) = view.stereo_neighbor_order(atom) {
                    molecule.set_stereo_neighbor_order(atom, order.to_vec());
                }
            }
        }
        Ok(())
    }

    /// A wavy carrier explicitly leaves each adjacent double bond unspecified.
    fn apply_wavy(&self, fragment: &mut BuiltFragment, positions: &[(usize, usize)]) {
        let mut unspecified = std::collections::HashSet::new();
        for (index, bond) in self.bonds.iter().enumerate() {
            if bond.display.as_deref() != Some("Wavy") {
                continue;
            }
            let (begin, end) = positions[index];
            for atom in [fragment.indices[&begin], fragment.indices[&end]] {
                for (_, bond_index) in fragment.molecule.neighbors(atom) {
                    if fragment.molecule.bond(bond_index).order == BondOrder::Double {
                        unspecified.insert(bond_index);
                    }
                }
            }
        }
        for &bond_index in &unspecified {
            let first = fragment.molecule.bond(bond_index).atom1;
            if matches!(
                fragment.molecule.atom(first).cip_code,
                Some(chematic_core::atom::CipCode::E | chematic_core::atom::CipCode::Z)
            ) {
                fragment.molecule.set_cip_code(first, None);
            }
            fragment.molecule.set_bond_stereo_unspecified(bond_index);
        }
        apply_ez_directions_from_2d_ex(&mut fragment.molecule, &fragment.coords, &unspecified);
    }

    fn flush(
        &mut self,
        results: &mut Vec<(Molecule, Vec<(f64, f64)>)>,
        options: &CdxmlParseOptions,
        limits: &CdxmlParseLimits,
        ctx: &mut ParseCtx,
    ) -> Result<(), CdxmlError> {
        let molecule_index = results.len();
        self.expand_abbreviations(limits)?;
        self.drop_top_level_connection_points(molecule_index, ctx);
        for atom in self.atoms.iter().filter(|atom| atom.abbreviation) {
            ctx.diagnose(
                molecule_index,
                "unexpanded_abbreviation",
                Some(&atom.id),
                "abbreviation node has no nested fragment; kept as a carbon placeholder".into(),
            );
        }
        if self.is_empty() {
            *self = Self::default();
            return Ok(());
        }
        let bonds = self.resolve_bond_positions()?;
        let mut fragment = self.build_molecule(&bonds.positions)?;
        self.apply_stereo(
            &mut fragment.molecule,
            &fragment.indices,
            &bonds,
            &fragment.coords,
            options,
        )?;
        assign_ez_from_2d(&mut fragment.molecule, &fragment.coords);
        self.apply_wavy(&mut fragment, &bonds.positions);
        self.record_abbreviations(molecule_index, ctx);
        results.push((fragment.molecule, fragment.coords));
        *self = Self::default();
        Ok(())
    }

    /// Text inside an open node's `<s>` element belongs to its label.
    fn push_text(&mut self, text: &str) {
        if let Some(node) = self.open_node.as_mut().filter(|node| node.in_string) {
            node.label.push_str(&unescape_xml(text));
        }
    }

    fn set_string_open(&mut self, open: bool) {
        if let Some(node) = self.open_node.as_mut() {
            node.in_string = open;
        }
    }
}

fn is_nondirectional_display(display: Option<&str>) -> bool {
    matches!(display, Some("Bold") | Some("Hash") | Some("Dash"))
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

/// ChemDraw generic query labels with a well-known meaning that a concrete
/// [`Molecule`] cannot express (no query-atom model): A (any atom except H),
/// Q (any heteroatom), X (halogen), M (metal), and their H-allowing variants.
const GENERIC_QUERY_LABELS: [&str; 8] = ["A", "AH", "Q", "QH", "X", "XH", "M", "MH"];

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

/// [`cdxml_atom_kind`], except that in tolerant mode an unsupported generic
/// label becomes a wildcard atom with a diagnostic instead of an error.
fn cdxml_atom_kind_in(
    attrs: &HashMap<String, String>,
    id: &str,
    molecule_index: usize,
    ctx: &mut ParseCtx,
) -> Result<(Element, bool, Option<RGroupLabel>), CdxmlError> {
    match cdxml_atom_kind(attrs) {
        Err(CdxmlError::UnsupportedPseudoatomLabel(label)) if ctx.tolerant => {
            let (code, meaning) = if GENERIC_QUERY_LABELS.contains(&label.as_str()) {
                ("generic_query_atom_as_wildcard", "generic query atom")
            } else {
                ("unrecognized_generic_label", "unrecognized generic label")
            };
            ctx.diagnose(
                molecule_index,
                code,
                Some(id),
                format!("{meaning} {label:?} read as a wildcard atom (*)"),
            );
            Ok((Element::C, true, None))
        }
        other => other,
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
    molecule_index: usize,
    ctx: &mut ParseCtx,
) -> Result<(), CdxmlError> {
    let attrs = parse_xml_attrs(line);
    check_attribute_limits(&attrs, limits)?;
    let node_type = attrs.get("NodeType").map(String::as_str);
    let abbreviation = matches!(node_type, Some("Nickname") | Some("Fragment"));
    if !line.ends_with("/>") {
        accumulator.open_node = Some(OpenNode {
            id: attrs.get("id").cloned(),
            abbreviation,
            label: String::new(),
            in_string: false,
        });
    }
    let Some(id) = attrs.get("id").cloned() else {
        return Ok(());
    };
    if node_type == Some("ExternalConnectionPoint") {
        let number = attrs
            .get("ExternalConnectionNum")
            .and_then(|value| value.trim().parse().ok());
        accumulator.connection_points.push((id, number));
        return Ok(());
    }
    if accumulator.atoms.len() >= limits.max_atoms {
        return Err(CdxmlError::ResourceLimit {
            resource: "atoms",
            actual: accumulator.atoms.len() + 1,
            limit: limits.max_atoms,
        });
    }
    if let Some(ordering) = attrs.get("BondOrdering") {
        let ordering = ordering.split_whitespace().map(str::to_string).collect();
        accumulator.bond_orderings.insert(id.clone(), ordering);
    }
    let (element, wildcard, r_group) = cdxml_atom_kind_in(&attrs, &id, molecule_index, ctx)?;
    let (x, y) = cdxml_coords(&attrs)?;
    let number = |key: &str| attrs.get(key).map(|value| value.trim());
    accumulator.atoms.push(AtomRec {
        id,
        element,
        charge: number("Charge").and_then(|v| v.parse().ok()).unwrap_or(0),
        isotope: number("Isotope")
            .and_then(|v| v.parse::<u16>().ok())
            .filter(|&value| value > 0),
        hydrogens: number("NumHydrogens").and_then(|v| v.parse().ok()),
        wildcard,
        r_group,
        x,
        y,
        abbreviation,
    });
    Ok(())
}

/// Decode a CDXML `Order` value. Space-separated multi-order values are
/// ChemDraw query bonds; the ones with a core query-bond type are accepted
/// in any token order.
fn cdxml_base_order(order: Option<&str>) -> Option<BondOrder> {
    let Some(order) = order else {
        return Some(BondOrder::Single);
    };
    let mut tokens = order.split_whitespace().collect::<Vec<_>>();
    tokens.sort_unstable();
    tokens.dedup();
    Some(match tokens.as_slice() {
        ["1"] => BondOrder::Single,
        ["2"] => BondOrder::Double,
        ["3"] => BondOrder::Triple,
        ["4"] => BondOrder::Quadruple,
        ["1.5"] => BondOrder::Aromatic,
        ["dative"] => BondOrder::Dative,
        ["any"] => BondOrder::QueryAny,
        ["1", "1.5"] => BondOrder::QuerySingleOrAromatic,
        ["1.5", "2"] => BondOrder::QueryDoubleOrAromatic,
        ["1", "2"] => BondOrder::QuerySingleOrDouble,
        _ => return None,
    })
}

fn cdxml_bond_order(attrs: &HashMap<String, String>) -> Result<BondOrder, CdxmlError> {
    let declared = attrs.get("Order").map(|value| value.trim());
    let base = cdxml_base_order(declared).ok_or_else(|| {
        CdxmlError::InvalidDocument(format!(
            "{UNSUPPORTED_BOND_ORDER_PREFIX}{}",
            declared.unwrap_or_default()
        ))
    })?;
    if base != BondOrder::Single {
        return Ok(base);
    }
    Ok(match attrs.get("Display").map(String::as_str) {
        Some("WedgeBegin") | Some("WedgeEnd") | Some("Bold") => BondOrder::Up,
        Some("Hash") | Some("Dash") | Some("WedgedHashBegin") | Some("WedgedHashEnd") => {
            BondOrder::Down
        }
        _ => BondOrder::Single,
    })
}

fn push_cdxml_bond(
    accumulator: &mut FragAccum,
    line: &str,
    limits: &CdxmlParseLimits,
    molecule_index: usize,
    ctx: &mut ParseCtx,
) -> Result<(), CdxmlError> {
    let attrs = parse_xml_attrs(line);
    check_attribute_limits(&attrs, limits)?;
    if accumulator.bonds.len() >= limits.max_bonds {
        return Err(CdxmlError::ResourceLimit {
            resource: "bonds",
            actual: accumulator.bonds.len() + 1,
            limit: limits.max_bonds,
        });
    }
    let begin = attrs
        .get("B")
        .cloned()
        .ok_or(CdxmlError::MissingBondEndpoint)?;
    let end = attrs
        .get("E")
        .cloned()
        .ok_or(CdxmlError::MissingBondEndpoint)?;
    let order = match cdxml_bond_order(&attrs) {
        Ok(order) => order,
        Err(error) if ctx.tolerant && error.unsupported_bond_order().is_some() => {
            ctx.diagnose(
                molecule_index,
                "unsupported_bond_order_dropped",
                attrs.get("id").map(String::as_str),
                format!("{error}; bond {begin}-{end} was dropped"),
            );
            return Ok(());
        }
        Err(error) => return Err(error),
    };
    accumulator.bonds.push(BondRec {
        id: attrs.get("id").cloned(),
        begin,
        end,
        order,
        display: attrs.get("Display").cloned(),
    });
    Ok(())
}

/// Same as [`parse_cdxml_all_with_options`], with explicit resource limits.
#[allow(clippy::type_complexity)]
pub fn parse_cdxml_all_with_options_and_limits(
    input: &str,
    options: &CdxmlParseOptions,
    limits: &CdxmlParseLimits,
) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError> {
    let mut ctx = ParseCtx::default();
    parse_cdxml_document(input, options, limits, &mut ctx)
}

/// Like [`parse_cdxml_all_with_options_and_limits`], but tolerant: a feature
/// that a concrete [`Molecule`] cannot represent no longer rejects the whole
/// document. Instead the rest of the document is kept and the feature is
/// reported in [`CdxmlParseReport::diagnostics`]:
///
/// - generic labels such as `A`, `Q`, `X`, `M` (and their `H` variants), and
///   unrecognized generic labels, become wildcard atoms (`*`), each with a
///   diagnostic naming the label (`generic_query_atom_as_wildcard` /
///   `unrecognized_generic_label`);
/// - bonds with an order that has no core bond type (for example
///   `hydrogen`, `ionic`) are dropped (`unsupported_bond_order_dropped`).
///
/// The report also lists the expanded ChemDraw abbreviations
/// ([`CdxmlParseReport::abbreviations`]), so an editor can show them
/// collapsed. Structural errors (malformed XML, unknown atom references,
/// resource limits, an abbreviation that cannot be expanded) still fail.
pub fn parse_cdxml_all_with_report(
    input: &str,
    options: &CdxmlParseOptions,
    limits: &CdxmlParseLimits,
) -> Result<CdxmlParseReport, CdxmlError> {
    let mut ctx = ParseCtx {
        tolerant: true,
        ..ParseCtx::default()
    };
    let molecules = parse_cdxml_document(input, options, limits, &mut ctx)?;
    Ok(CdxmlParseReport {
        molecules,
        abbreviations: ctx.abbreviations,
        diagnostics: ctx.diagnostics,
    })
}

#[allow(clippy::type_complexity)]
fn parse_cdxml_document(
    input: &str,
    options: &CdxmlParseOptions,
    limits: &CdxmlParseLimits,
    ctx: &mut ParseCtx,
) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError> {
    check_cdxml_input_limits(input, limits)?;
    // The bottom entry is the top-level fragment; entries above it are
    // fragments nested inside abbreviation nodes.
    let mut stack = vec![FragAccum::default()];
    let mut results = Vec::new();
    let mut fragment_count = 0;
    for raw in crate::cdxml_document::logical_cdxml_lines(input) {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        let index = results.len();
        let depth = stack.len();
        let top = stack.last_mut().expect("fragment stack is never empty");
        if line.starts_with("<fragment") {
            if fragment_count >= limits.max_fragments {
                return Err(CdxmlError::ResourceLimit {
                    resource: "fragments",
                    actual: fragment_count + 1,
                    limit: limits.max_fragments,
                });
            }
            fragment_count += 1;
            if top.open_node.is_none() {
                *top = FragAccum::default();
            } else if !line.ends_with("/>") {
                stack.push(FragAccum::default());
            }
        } else if line.starts_with("</fragment>") {
            if depth > 1 {
                let inner = stack.pop().expect("nested fragment");
                let parent = stack.last_mut().expect("parent fragment");
                parent.attach_nested(inner, limits, index, ctx)?;
            } else {
                top.flush(&mut results, options, limits, ctx)?;
            }
        } else if is_n_tag(line) {
            push_cdxml_atom(top, line, limits, index, ctx)?;
        } else if crate::cdxml_document::is_close_tag(line, "n") {
            top.close_node();
        } else if crate::cdxml_document::is_open_tag(line, "s") {
            top.set_string_open(!line.ends_with("/>"));
        } else if crate::cdxml_document::is_close_tag(line, "s") {
            top.set_string_open(false);
        } else if !line.starts_with('<') {
            top.push_text(line);
        } else if is_b_tag(line) {
            push_cdxml_bond(top, line, limits, index, ctx)?;
        }
    }
    stack.truncate(1);
    stack[0].flush(&mut results, options, limits, ctx)?;
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

/// True if `line` starts a CDXML `<n` atom node tag. ChemDraw puts a line
/// break right after the element name (`<n` then `id="..."` on the next line).
fn is_n_tag(line: &str) -> bool {
    crate::cdxml_document::is_open_tag(line, "n")
}

/// True if `line` starts a CDXML `<b` bond tag (see [`is_n_tag`]).
fn is_b_tag(line: &str) -> bool {
    crate::cdxml_document::is_open_tag(line, "b")
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

    let wavy = wavy_carrier_bonds(mol);
    for (bond_index, bond) in mol.bonds() {
        let b = bond.atom1.0 + 1;
        let e = bond.atom2.0 + 1;
        // Wedge/hash stereo is Single order + Display attribute (mirrors the
        // parser's decoding at is_b_tag handling above); other orders map
        // straight to the Order attribute and ignore Display. Both wedge
        // kinds are written directionally with the narrow end on `B`
        // (`atom1`, the stereo reference endpoint), so the parser never has
        // to treat them as ambiguous non-directional glyphs.
        let (order_attr, display_attr) = match bond.order {
            _ if wavy.contains(&bond_index) => (
                if bond.order == BondOrder::Double {
                    "2"
                } else {
                    "1"
                },
                Some("Wavy"),
            ),
            BondOrder::Double => ("2", None),
            BondOrder::Triple => ("3", None),
            BondOrder::Quadruple => ("4", None),
            BondOrder::Dative => ("dative", None),
            BondOrder::Aromatic => ("1.5", None),
            BondOrder::QueryAny => ("any", None),
            BondOrder::QuerySingleOrAromatic => ("1 1.5", None),
            BondOrder::QueryDoubleOrAromatic => ("2 1.5", None),
            BondOrder::QuerySingleOrDouble => ("1 2", None),
            BondOrder::Up => ("1", Some("WedgeBegin")),
            BondOrder::Down => ("1", Some("WedgedHashBegin")),
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

/// Bonds that should carry `Display="Wavy"` so every double bond whose
/// stereo the molecule records as explicitly unspecified reads back that way.
///
/// The parser marks every double bond touching a wavy bond's endpoints as
/// unspecified, so a single-bond carrier is only chosen when it cannot
/// spill onto a double bond with defined stereo; otherwise the double bond
/// itself carries the glyph.
fn wavy_carrier_bonds(mol: &Molecule) -> std::collections::HashSet<chematic_core::BondIdx> {
    let mut carriers = std::collections::HashSet::new();
    let safe_atom = |atom: AtomIdx| {
        mol.neighbors(atom).all(|(_, bond)| {
            mol.bond(bond).order != BondOrder::Double || mol.bond_stereo_is_unspecified(bond)
        })
    };
    for (double, bond) in mol.bonds() {
        if bond.order != BondOrder::Double || !mol.bond_stereo_is_unspecified(double) {
            continue;
        }
        let carrier = [bond.atom1, bond.atom2]
            .into_iter()
            .flat_map(|atom| mol.neighbors(atom))
            .map(|(_, candidate)| candidate)
            .find(|&candidate| {
                let candidate_bond = mol.bond(candidate);
                candidate_bond.order == BondOrder::Single
                    && safe_atom(candidate_bond.atom1)
                    && safe_atom(candidate_bond.atom2)
            });
        carriers.insert(carrier.unwrap_or(double));
    }
    carriers
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "cdxml/tests.rs"]
mod tests;
