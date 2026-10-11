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

impl FragAccum {
    fn is_empty(&self) -> bool {
        self.atom_ids.is_empty() && self.bond_bs.is_empty()
    }

    fn flush(
        &mut self,
        results: &mut Vec<(Molecule, Vec<(f64, f64)>)>,
        options: &CdxmlParseOptions,
    ) -> Result<(), CdxmlError> {
        if self.is_empty() {
            return Ok(());
        }

        let mut id_to_pos: HashMap<&str, usize> = HashMap::new();
        for (i, id) in self.atom_ids.iter().enumerate() {
            id_to_pos.insert(id.as_str(), i);
        }

        // Resolve every bond's endpoints to atom positions up front, both to
        // validate atom refs early and so degree (needed by the ambiguity
        // check below) can be computed before any bond is added.
        let mut bond_pos: Vec<(usize, usize)> = Vec::with_capacity(self.bond_bs.len());
        let mut degree: Vec<usize> = vec![0; self.atom_ids.len()];
        for k in 0..self.bond_bs.len() {
            let pos_b = *id_to_pos
                .get(self.bond_bs[k].as_str())
                .ok_or_else(|| CdxmlError::UnknownAtomRef(self.bond_bs[k].clone()))?;
            let pos_e = *id_to_pos
                .get(self.bond_es[k].as_str())
                .ok_or_else(|| CdxmlError::UnknownAtomRef(self.bond_es[k].clone()))?;
            degree[pos_b] += 1;
            degree[pos_e] += 1;
            bond_pos.push((pos_b, pos_e));
        }

        let mut builder = MoleculeBuilder::new();
        let mut idx_map: HashMap<usize, AtomIdx> = HashMap::new();
        let mut coords: Vec<(f64, f64)> = Vec::new();

        for i in 0..self.atom_ids.len() {
            let mut a = Atom::new(self.atom_elems[i]);
            a.wildcard = self.atom_wildcards[i];
            a.charge = self.atom_charges[i];
            a.isotope = self.atom_isotopes[i];
            a.hydrogen_count = if a.wildcard { Some(0) } else { self.atom_h[i] };
            let new_idx = builder.add_atom(a);
            if let Some(label) = self.atom_rgroups[i] {
                builder.set_r_group(new_idx, label);
            }
            idx_map.insert(i, new_idx);
            coords.push((self.atom_xs[i], self.atom_ys[i]));
        }

        // The REAL molecule keeps every wedge bond's Display-derived
        // BondOrder::Up/Down and its original B/E atom order exactly as
        // parsed, for ANY Display value -- structural/round-trip fidelity
        // is unconditional, independent of whether that bond ends up
        // feeding chirality perception below.
        for (k, &(pos_b, pos_e)) in bond_pos.iter().enumerate() {
            let a1 = idx_map[&pos_b];
            let a2 = idx_map[&pos_e];
            builder.add_bond(a1, a2, self.bond_ords[k]).map_err(|_| {
                CdxmlError::UnknownAtomRef(format!("{} {}", self.bond_bs[k], self.bond_es[k]))
            })?;
        }

        let mut mol = builder.build();

        // Non-directional wedge displays (`Bold`/`Hash`/`Dash` -- a plain
        // thick or dashed line, no narrow/wide end) carry NO Begin/End
        // convention at all in the CDXML spec, unlike `WedgeBegin`/
        // `WedgeEnd`/`WedgedHashBegin`/`WedgedHashEnd`. Which atom a `<b>`
        // element happens to list first as `B` is an arbitrary artifact of
        // however the file was written -- it is not a "narrow end" signal,
        // so `chematic_perception::wedge_z`'s `bond.atom1` == the tip/
        // reference-atom convention (correct for a directional wedge, whose
        // Begin atom genuinely IS the tip) must not be applied naively to
        // these using the raw B/E order. Concretely: swapping which atom a
        // CDXML author lists as B vs E for the identical drawing must NOT
        // change the perceived configuration for `Bold`/`Hash`/`Dash`,
        // while it correctly DOES for a real `WedgeBegin`/`WedgeEnd` pair
        // (the Begin atom is genuinely a different physical reference).
        //
        // Resolved by perceiving parity on a throwaway VIEW of the
        // molecule, not the real one -- so the real molecule's own bond
        // orders/atom order are never touched by any of this -- with each
        // non-directional bond's contribution to that view decided as:
        //   - `options.infer_nondirectional_stereo` is `false` (the
        //     default): downgraded to `Single` in the view unconditionally.
        //     `Bold` in particular is sometimes drawn for visual emphasis,
        //     not stereo intent, so this reader does not manufacture
        //     stereo from it unless explicitly asked to.
        //   - `true`, exactly one endpoint qualifies as a stereocenter
        //     candidate (3-4 neighbours, matching
        //     [`chematic_perception`]'s own gate): reordered in the view so
        //     that endpoint is atom1 -- then `wedge_z`'s existing
        //     atom1-is-the-reference convention becomes correct by
        //     construction, independent of the original B/E order.
        //   - `true`, both endpoints qualify: one non-directional bond
        //     cannot mean "toward the viewer" from both ends of the same
        //     bond at once -- downgraded to `Single` in the view so neither
        //     gets a chirality it can't justify (a directional wedge in the
        //     same situation is NOT ambiguous, since the Begin atom is
        //     unambiguously the reference, and is never altered here).
        //   - `true`, neither qualifies: no-op (the shared perception gate
        //     already skips both regardless).
        // Directional wedges are never altered in the view either way.
        let any_nondirectional_wedge = self
            .bond_display
            .iter()
            .any(|d| matches!(d.as_deref(), Some("Bold") | Some("Hash") | Some("Dash")));
        if !any_nondirectional_wedge {
            apply_local_parity_from_wedges(&mut mol, &coords);
        } else {
            let mut view = mol.clone();
            for (k, &(pos_b, pos_e)) in bond_pos.iter().enumerate() {
                let is_nondirectional_wedge = matches!(
                    self.bond_display[k].as_deref(),
                    Some("Bold") | Some("Hash") | Some("Dash")
                );
                if !is_nondirectional_wedge {
                    continue;
                }
                let (a1, a2) = (idx_map[&pos_b], idx_map[&pos_e]);
                let Some((bond_idx, _)) = view.bond_between(a1, a2) else {
                    continue;
                };
                if !options.infer_nondirectional_stereo {
                    view.set_bond_order(bond_idx, BondOrder::Single);
                    continue;
                }
                let eligible = |d: usize| (3..=4).contains(&d);
                let (b_elig, e_elig) = (eligible(degree[pos_b]), eligible(degree[pos_e]));
                if b_elig && e_elig {
                    view.set_bond_order(bond_idx, BondOrder::Single);
                } else if e_elig && !b_elig {
                    // The candidate (stereocenter) is E, not B. `wedge_z`
                    // reads `bond.atom1` as the reference atom -- flip
                    // Up<->Down (equivalent to "as seen from E instead of
                    // B") rather than swapping atom1/atom2 in place: a
                    // swap would require remove_bond+add_bond, which
                    // perturbs `neighbors()` iteration order for BOTH
                    // endpoints (adjacency is rebuilt from the bond list's
                    // new order), silently changing which neighbor
                    // `chematic_perception` picks as tetrahedral apex --
                    // confirmed empirically: an earlier swap-based version
                    // of this fix produced a real, non-obvious B/E-order
                    // dependence via exactly this side channel. Flipping
                    // just the order field has no such effect (checked in
                    // `set_bond_order`'s own doc).
                    let flipped = match view.bond(bond_idx).order {
                        BondOrder::Up => BondOrder::Down,
                        BondOrder::Down => BondOrder::Up,
                        other => other,
                    };
                    view.set_bond_order(bond_idx, flipped);
                }
            }
            apply_local_parity_from_wedges(&mut view, &coords);
            for i in 0..self.atom_ids.len() {
                let idx = idx_map[&i];
                let chirality = view.atom(idx).chirality;
                if chirality != chematic_core::Chirality::None {
                    mol.set_chirality(idx, chirality);
                    if let Some(order) = view.stereo_neighbor_order(idx) {
                        mol.set_stereo_neighbor_order(idx, order.to_vec());
                    }
                }
            }
        }

        // Tetrahedral parity (above) runs before E/Z direction (below) --
        // same ordering rationale as mol2000.rs/mol3000.rs/mrv.rs: E/Z
        // writes to a separate side channel (cip_code here) and never
        // touches bond.order, so it can't disturb the wedge/hash data
        // parity perception just consumed.
        // Derive E/Z stereo from 2D atom positions (RDKit issue #9356: CDXML loses E/Z).
        assign_ez_from_2d(&mut mol, &coords);
        results.push((mol, coords));
        *self = FragAccum::default();
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

/// Same as [`parse_cdxml_all_with_options`], with explicit resource limits.
#[allow(clippy::type_complexity)]
pub fn parse_cdxml_all_with_options_and_limits(
    input: &str,
    options: &CdxmlParseOptions,
    limits: &CdxmlParseLimits,
) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError> {
    if input.len() > limits.max_input_bytes {
        return Err(CdxmlError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let mut acc = FragAccum::default();
    let mut results: Vec<(Molecule, Vec<(f64, f64)>)> = Vec::new();
    let mut fragment_count = 0usize;

    for (line_index, raw_line) in input.lines().enumerate() {
        if line_index >= limits.max_lines {
            return Err(CdxmlError::ResourceLimit {
                resource: "lines",
                actual: line_index + 1,
                limit: limits.max_lines,
            });
        }
        if raw_line.len() > limits.max_line_bytes {
            return Err(CdxmlError::ResourceLimit {
                resource: "line bytes",
                actual: raw_line.len(),
                limit: limits.max_line_bytes,
            });
        }
    }

    // Physical-line budgets apply before tokenization. Producers may place
    // several tags on one line; reuse the document editor's quote-aware
    // splitter so minified fragments cannot silently become empty graphs.
    for raw_line in crate::cdxml_document::logical_cdxml_lines(input) {
        let line = raw_line.trim();
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
            acc = FragAccum::default();
            continue;
        }

        if line.starts_with("</fragment>") {
            acc.flush(&mut results, options)?;
            continue;
        }

        if is_n_tag(line) {
            let attrs = parse_xml_attrs(line);
            check_attribute_limits(&attrs, limits)?;
            let id = match attrs.get("id") {
                Some(s) => s.clone(),
                None => continue,
            };
            let pseudo_label =
                if attrs.get("NodeType").map(String::as_str) == Some("GenericNickname") {
                    Some(attrs.get("GenericNickname").cloned().ok_or_else(|| {
                        CdxmlError::UnsupportedPseudoatomLabel(
                            "missing GenericNickname attribute".to_string(),
                        )
                    })?)
                } else {
                    None
                };
            let (element, wildcard, r_group) = match pseudo_label.as_deref() {
                Some("*") => (Element::C, true, None),
                Some("R") => (Element::C, true, Some(RGroupLabel::unnumbered())),
                Some(label) if label.starts_with('R') => {
                    let number = label[1..]
                        .parse::<u16>()
                        .ok()
                        .and_then(RGroupLabel::numbered)
                        .ok_or_else(|| CdxmlError::UnsupportedPseudoatomLabel(label.to_string()))?;
                    (Element::C, true, Some(number))
                }
                Some(label) => {
                    return Err(CdxmlError::UnsupportedPseudoatomLabel(label.to_string()));
                }
                None => {
                    let element_num: u32 = attrs
                        .get("Element")
                        .and_then(|s| s.trim().parse().ok())
                        .unwrap_or(6);
                    if element_num > 255 {
                        return Err(CdxmlError::UnknownAtomicNumber(element_num));
                    }
                    let element = Element::from_atomic_number(element_num as u8)
                        .ok_or(CdxmlError::UnknownAtomicNumber(element_num))?;
                    (element, false, None)
                }
            };
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
            if acc.atom_ids.len() >= limits.max_atoms {
                return Err(CdxmlError::ResourceLimit {
                    resource: "atoms",
                    actual: acc.atom_ids.len() + 1,
                    limit: limits.max_atoms,
                });
            }
            acc.atom_ids.push(id);
            acc.atom_elems.push(element);
            acc.atom_charges.push(charge);
            acc.atom_isotopes.push(isotope);
            acc.atom_h.push(hcount);
            acc.atom_wildcards.push(wildcard);
            acc.atom_rgroups.push(r_group);
            acc.atom_xs.push(x);
            acc.atom_ys.push(y);
            continue;
        }

        if is_b_tag(line) {
            let attrs = parse_xml_attrs(line);
            check_attribute_limits(&attrs, limits)?;
            if acc.bond_bs.len() >= limits.max_bonds {
                return Err(CdxmlError::ResourceLimit {
                    resource: "bonds",
                    actual: acc.bond_bs.len() + 1,
                    limit: limits.max_bonds,
                });
            }
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
            acc.bond_display.push(attrs.get("Display").cloned());
        }
    }

    // Handle documents without explicit </fragment> closing tags.
    acc.flush(&mut results, options)?;

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
