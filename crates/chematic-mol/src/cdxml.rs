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

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};
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
            CdxmlError::UnknownAtomicNumber(n) => write!(f, "unknown atomic number: {n}"),
            CdxmlError::UnknownAtomRef(s) => write!(f, "unknown atom ref: {s}"),
            CdxmlError::MissingBondEndpoint => write!(f, "bond missing B or E attribute"),
            CdxmlError::InvalidCoords(s) => write!(f, "invalid p coords: {s}"),
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
            a.charge = self.atom_charges[i];
            a.isotope = self.atom_isotopes[i];
            a.hydrogen_count = self.atom_h[i];
            let new_idx = builder.add_atom(a);
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
    use chematic_core::Chirality;

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

    // ── RDKit issue #9359: actual tetrahedral chirality perception ────────
    //
    // The tests above only ever checked `bond.order` -- they confirm the
    // *bond* survives as a wedge, not that a *chirality* was ever perceived
    // from it. Before this fix, `flush()` never called any wedge-to-parity
    // mechanism at all, for ANY Display value including plain WedgeBegin --
    // so `Atom.chirality` was always `Chirality::None` for every CDXML
    // molecule, no matter how it was drawn. These tests assert the actual
    // RDKit #9359 expectation (`center_atom(mol).GetChiralTag() !=
    // CHI_UNSPECIFIED`, i.e. `Atom.chirality != Chirality::None` here), and
    // then go further to check *meaning*, not just non-None.

    /// RDKit #9359's exact repro fixture (`BOLD_BASE_CDXML`, atom/bond ids
    /// preserved), parameterised on the `Display` of the C1-F2 bond (`101`)
    /// so the no-wedge/Bold/Hash/WedgeBegin cases share one source of
    /// truth. Atom 1 (implicit Carbon, the default when `Element` is
    /// omitted) is the only degree-4 atom -- F(2), the ring-continuation
    /// C(3)->C(4), Cl(7), and CH3(6) -- matching RDKit's own
    /// `center_atom` (`atomic_number==6 and degree==4`) selector.
    fn rdkit_9359_cdxml(display: Option<&str>) -> String {
        let bond101 = match display {
            Some(d) => format!(r#"<b id="101" B="1" E="2" Display="{d}"/>"#),
            None => r#"<b id="101" B="1" E="2"/>"#.to_string(),
        };
        format!(
            r#"<CDXML>
<page>
<fragment>
<n id="1" p="270.73 355.40"/>
<n id="2" p="263.44 367.76" Element="9"/>
<n id="3" p="283.22 348.33"/>
<n id="4" p="295.58 355.62"/>
<n id="6" p="258.37 348.11"/>
<n id="7" p="277.79 367.89" Element="17" NumHydrogens="0"><t><s>Cl</s></t></n>
{bond101}
<b id="102" B="1" E="3"/>
<b id="103" B="1" E="7"/>
<b id="104" B="1" E="6"/>
<b id="105" B="3" E="4"/>
</fragment>
</page>
</CDXML>"#
        )
    }

    const RDKIT_9359_CENTER: chematic_core::AtomIdx = chematic_core::AtomIdx(0);

    /// Parses with [`CdxmlParseOptions::infer_nondirectional_stereo`] on --
    /// the opt-in path, used by tests specifically exercising Bold/Hash
    /// chirality perception (which is off by default, see the
    /// default-behavior tests below).
    fn parse_infer(cdxml: &str) -> (Molecule, Vec<(f64, f64)>) {
        parse_cdxml_with_options(
            cdxml,
            &CdxmlParseOptions {
                infer_nondirectional_stereo: true,
            },
        )
        .unwrap()
    }

    #[test]
    fn rdkit_9359_no_display_has_no_chirality() {
        // Flat drawing, no wedge at all -- must NOT invent stereo, with or
        // without the opt-in (nothing to infer from either way).
        let (mol, _) = parse_cdxml(&rdkit_9359_cdxml(None)).unwrap();
        assert_eq!(mol.atom(RDKIT_9359_CENTER).chirality, Chirality::None);
        let (mol, _) = parse_infer(&rdkit_9359_cdxml(None));
        assert_eq!(mol.atom(RDKIT_9359_CENTER).chirality, Chirality::None);
    }

    // ── Default behavior: Bold/Hash do NOT perceive chirality unopted ─────
    //
    // `Bold` in particular is sometimes drawn by chemists as a plain thick
    // line for visual emphasis, not stereo intent (RDKit #9359's own
    // discussion leans toward this being opt-in, not silently default) --
    // so `parse_cdxml`/`parse_cdxml_all` (no explicit options) must NOT
    // manufacture a chirality from it. `bond.order` still faithfully
    // records `BondOrder::Up`/`Down` either way (structural fidelity is
    // unconditional) -- only chirality PERCEPTION is gated.

    #[test]
    fn rdkit_9359_bold_default_does_not_perceive_chirality() {
        let (mol, _) = parse_cdxml(&rdkit_9359_cdxml(Some("Bold"))).unwrap();
        assert_eq!(
            mol.atom(RDKIT_9359_CENTER).chirality,
            Chirality::None,
            "Bold must NOT perceive chirality by default (opt-in required)"
        );
        assert_eq!(
            mol.bond(chematic_core::BondIdx(0)).order,
            BondOrder::Up,
            "the bond's own order must still faithfully record Up regardless \
             of whether chirality perception is opted in"
        );
    }

    #[test]
    fn rdkit_9359_hash_default_does_not_perceive_chirality() {
        let (mol, _) = parse_cdxml(&rdkit_9359_cdxml(Some("Hash"))).unwrap();
        assert_eq!(
            mol.atom(RDKIT_9359_CENTER).chirality,
            Chirality::None,
            "Hash must NOT perceive chirality by default (opt-in required)"
        );
        assert_eq!(mol.bond(chematic_core::BondIdx(0)).order, BondOrder::Down);
    }

    /// Directional wedges are NOT gated by the opt-in -- they're
    /// unambiguous stereo intent (a real tapered wedge), so the default
    /// (no options) already perceives chirality from them, matching the
    /// RDKit #9359 expectation (`center_atom(mol).GetChiralTag() !=
    /// CHI_UNSPECIFIED`) out of the box.
    #[test]
    fn rdkit_9359_wedge_begin_default_perceives_chirality() {
        let (mol, _) = parse_cdxml(&rdkit_9359_cdxml(Some("WedgeBegin"))).unwrap();
        assert_ne!(
            mol.atom(RDKIT_9359_CENTER).chirality,
            Chirality::None,
            "a real directional wedge must perceive chirality without any opt-in"
        );
        assert!(mol.stereo_neighbor_order(RDKIT_9359_CENTER).is_some());
    }

    // ── Opt-in behavior (`infer_nondirectional_stereo: true`) ─────────────

    #[test]
    fn rdkit_9359_bold_perceives_chirality_when_opted_in() {
        let (mol, _) = parse_infer(&rdkit_9359_cdxml(Some("Bold")));
        assert_ne!(
            mol.atom(RDKIT_9359_CENTER).chirality,
            Chirality::None,
            "RDKit #9359: Display=\"Bold\" must perceive a real chirality when \
             opted in, matching center_atom(mol).GetChiralTag() != CHI_UNSPECIFIED"
        );
        assert!(mol.stereo_neighbor_order(RDKIT_9359_CENTER).is_some());
    }

    #[test]
    fn rdkit_9359_hash_perceives_chirality_when_opted_in() {
        let (mol, _) = parse_infer(&rdkit_9359_cdxml(Some("Hash")));
        assert_ne!(
            mol.atom(RDKIT_9359_CENTER).chirality,
            Chirality::None,
            "RDKit #9359: Display=\"Hash\" (undirectional) must also perceive \
             chirality when opted in"
        );
    }

    /// The meaning check the reviewer required, not just non-None: Bold
    /// ("coming toward viewer") and Hash ("going away from viewer") are
    /// chemically OPPOSITE for the exact same 2D layout -- if perception
    /// were just guessing, there'd be no reason for these to differ let
    /// alone invert consistently.
    #[test]
    fn rdkit_9359_bold_and_hash_are_opposite_configurations() {
        let (mol_bold, _) = parse_infer(&rdkit_9359_cdxml(Some("Bold")));
        let (mol_hash, _) = parse_infer(&rdkit_9359_cdxml(Some("Hash")));
        let bold_chirality = mol_bold.atom(RDKIT_9359_CENTER).chirality;
        let hash_chirality = mol_hash.atom(RDKIT_9359_CENTER).chirality;
        assert_ne!(bold_chirality, Chirality::None);
        assert_ne!(hash_chirality, Chirality::None);
        assert_ne!(
            bold_chirality, hash_chirality,
            "Bold (toward viewer) and Hash (away from viewer) on the identical \
             2D layout must give opposite configurations, not the same one"
        );

        // Independent oracle cross-check: accurate CIP assignment must
        // agree that these two are opposite, not just chematic's own raw
        // Chirality enum by coincidence.
        let cip_bold = chematic_chem::assign_cip(&mol_bold).get(RDKIT_9359_CENTER);
        let cip_hash = chematic_chem::assign_cip(&mol_hash).get(RDKIT_9359_CENTER);
        assert!(
            cip_bold.is_some() && cip_hash.is_some(),
            "{cip_bold:?} {cip_hash:?}"
        );
        assert_ne!(
            cip_bold, cip_hash,
            "independent CIP oracle must also see Bold/Hash as opposite R/S"
        );
    }

    /// "Bold" and "WedgeBegin" both encode the same "coming toward viewer"
    /// convention -- for the identical layout they must produce the SAME
    /// chirality, not just both-non-None independently.
    #[test]
    fn rdkit_9359_bold_matches_wedge_begin_configuration() {
        let (mol_bold, _) = parse_infer(&rdkit_9359_cdxml(Some("Bold")));
        let (mol_wedge, _) = parse_cdxml(&rdkit_9359_cdxml(Some("WedgeBegin"))).unwrap();
        assert_eq!(
            mol_bold.atom(RDKIT_9359_CENTER).chirality,
            mol_wedge.atom(RDKIT_9359_CENTER).chirality,
        );
    }

    /// B/E-reversed variant, non-directional (required by review, and the
    /// review's own required *expected value* corrected from an earlier,
    /// wrong assumption): `Bold`/`Hash`/`Dash` have NO Begin/End
    /// convention in the CDXML spec at all -- unlike `WedgeBegin`/
    /// `WedgeEnd`, whose Begin atom is a real, drawn narrow end. Which
    /// atom a `<b>` element happens to list first as `B` for a `Bold`
    /// bond is an arbitrary artifact of how the file was written, not a
    /// stereo-relevant signal -- so swapping B/E for the SAME `Bold`
    /// drawing must PRESERVE the perceived configuration, not invert it.
    #[test]
    fn rdkit_9359_bold_be_reversed_preserves_configuration() {
        let (mol_forward, _) = parse_infer(&rdkit_9359_cdxml(Some("Bold")));
        let reversed = rdkit_9359_cdxml(Some("Bold")).replace(
            r#"<b id="101" B="1" E="2" Display="Bold"/>"#,
            r#"<b id="101" B="2" E="1" Display="Bold"/>"#,
        );
        let (mol_reversed, _) = parse_infer(&reversed);
        let forward = mol_forward.atom(RDKIT_9359_CENTER).chirality;
        let rev = mol_reversed.atom(RDKIT_9359_CENTER).chirality;
        assert_ne!(forward, Chirality::None);
        assert_eq!(
            forward, rev,
            "Bold has no Begin/End convention -- reversing B/E for the identical \
             drawing must NOT change the perceived configuration"
        );
    }

    /// Same B/E-reversal check for `Hash` (the other non-directional
    /// display), independently -- not inferred from the Bold case alone.
    #[test]
    fn rdkit_9359_hash_be_reversed_preserves_configuration() {
        let (mol_forward, _) = parse_infer(&rdkit_9359_cdxml(Some("Hash")));
        let reversed = rdkit_9359_cdxml(Some("Hash")).replace(
            r#"<b id="101" B="1" E="2" Display="Hash"/>"#,
            r#"<b id="101" B="2" E="1" Display="Hash"/>"#,
        );
        let (mol_reversed, _) = parse_infer(&reversed);
        let forward = mol_forward.atom(RDKIT_9359_CENTER).chirality;
        let rev = mol_reversed.atom(RDKIT_9359_CENTER).chirality;
        assert_ne!(forward, Chirality::None);
        assert_eq!(
            forward, rev,
            "Hash has no Begin/End convention -- reversing B/E for the identical \
             drawing must NOT change the perceived configuration"
        );
    }

    /// Negative control, directional: `WedgeBegin` DOES have a real
    /// Begin/End convention (the Begin atom is the drawn narrow end), so
    /// switching to `WedgeEnd` for the SAME B/E atom order must genuinely
    /// invert the configuration (the reference atom flips) -- proves the
    /// preserve-under-swap behavior above is specific to the
    /// non-directional set, not a blanket "wedges don't encode direction"
    /// regression. Uses the default parser (no opt-in needed, directional
    /// wedges are unaffected by the flag either way).
    #[test]
    fn rdkit_9359_wedge_begin_vs_wedge_end_inverts_configuration() {
        let (mol_begin, _) = parse_cdxml(&rdkit_9359_cdxml(Some("WedgeBegin"))).unwrap();
        let (mol_end, _) = parse_cdxml(&rdkit_9359_cdxml(Some("WedgeEnd"))).unwrap();
        let begin = mol_begin.atom(RDKIT_9359_CENTER).chirality;
        let end = mol_end.atom(RDKIT_9359_CENTER).chirality;
        assert_ne!(begin, Chirality::None);
        assert_ne!(end, Chirality::None);
        assert_ne!(
            begin, end,
            "WedgeBegin vs WedgeEnd (same B/E order) must invert -- the Begin \
             atom is a real, drawn reference, unlike Bold/Hash"
        );
    }

    /// Reflection (negate all Y coordinates) must invert the perceived
    /// configuration; rotation + translation must leave it unchanged.
    #[test]
    fn rdkit_9359_bold_reflection_inverts_rotation_translation_preserves() {
        let base = rdkit_9359_cdxml(Some("Bold"));
        let (mol_base, coords_base) = parse_infer(&base);
        let base_chirality = mol_base.atom(RDKIT_9359_CENTER).chirality;
        assert_ne!(base_chirality, Chirality::None);

        // Reflection: negate Y for every atom's `p` attribute.
        let reflect_re = regex_lite_reflect_y(&base);
        let (mol_reflected, _) = parse_infer(&reflect_re);
        assert_eq!(
            mol_reflected.atom(RDKIT_9359_CENTER).chirality,
            invert(base_chirality),
            "mirroring the 2D layout (negate Y) must invert the configuration"
        );

        // Rotation (37 degrees) + translation (+1000, -500): must be a
        // no-op for chirality -- recompute coords directly (not via the
        // CDXML string) since only chematic_perception's own math needs
        // exercising here, not the parser's coordinate formatting.
        let theta: f64 = 37.0_f64.to_radians();
        let (sin_t, cos_t) = theta.sin_cos();
        let rotated_translated: Vec<(f64, f64)> = coords_base
            .iter()
            .map(|&(x, y)| {
                let (rx, ry) = (x * cos_t - y * sin_t, x * sin_t + y * cos_t);
                (rx + 1000.0, ry - 500.0)
            })
            .collect();
        let mut mol_rt = mol_base.clone();
        chematic_perception::apply_local_parity_from_wedges(&mut mol_rt, &rotated_translated);
        assert_eq!(
            mol_rt.atom(RDKIT_9359_CENTER).chirality,
            base_chirality,
            "rotation + translation must not change the perceived configuration"
        );
    }

    fn invert(c: Chirality) -> Chirality {
        match c {
            Chirality::Clockwise => Chirality::CounterClockwise,
            Chirality::CounterClockwise => Chirality::Clockwise,
            Chirality::None => Chirality::None,
            sp @ Chirality::SquarePlanar(_) => sp,
        }
    }

    /// Negates the Y coordinate of every atom's `p="x y"` attribute in a
    /// CDXML fragment string -- test-only helper, not a general XML editor.
    fn regex_lite_reflect_y(cdxml: &str) -> String {
        let mut out = String::with_capacity(cdxml.len());
        let mut rest = cdxml;
        while let Some(start) = rest.find("p=\"") {
            out.push_str(&rest[..start + 3]);
            rest = &rest[start + 3..];
            let end = rest.find('"').expect("unterminated p attribute");
            let coord = &rest[..end];
            let mut parts = coord.split_whitespace();
            let x: f64 = parts.next().unwrap().parse().unwrap();
            let y: f64 = parts.next().unwrap().parse().unwrap();
            out.push_str(&format!("{x} {}", -y));
            rest = &rest[end..];
        }
        out.push_str(rest);
        out
    }

    // ── Ambiguity: non-directional wedge between two stereocenters ────────

    /// Both endpoints of a single "Bold" bond independently qualify as
    /// stereocenter candidates (degree 4 each) -- a non-directional wedge
    /// has no narrow/wide-end convention to say which one it actually
    /// describes, so BOTH must abstain (fail closed) rather than each
    /// independently claiming a configuration derived from contradictory
    /// height assignments of the same bond.
    #[test]
    fn bold_bond_between_two_stereocenters_abstains_both_sides() {
        let cdxml = r#"<CDXML><fragment>
<n id="1" p="0 0"/>
<n id="2" p="10 0"/>
<n id="3" p="0 10" Element="9"/>
<n id="4" p="-10 0" Element="17"/>
<n id="5" p="0 -10" Element="35"/>
<n id="6" p="20 10" Element="9"/>
<n id="7" p="20 -10" Element="17"/>
<n id="8" p="10 -20" Element="35"/>
<b B="1" E="2" Display="Bold"/>
<b B="1" E="3"/>
<b B="1" E="4"/>
<b B="1" E="5"/>
<b B="2" E="6"/>
<b B="2" E="7"/>
<b B="2" E="8"/>
</fragment></CDXML>"#;
        let (mol, _) = parse_infer(cdxml);
        let a1 = chematic_core::AtomIdx(0);
        let a2 = chematic_core::AtomIdx(1);
        assert_eq!(mol.neighbors(a1).count(), 4);
        assert_eq!(mol.neighbors(a2).count(), 4);
        assert_eq!(
            mol.atom(a1).chirality,
            Chirality::None,
            "ambiguous non-directional wedge between two stereocenters must abstain"
        );
        assert_eq!(mol.atom(a2).chirality, Chirality::None);
        assert_eq!(
            mol.bond(chematic_core::BondIdx(0)).order,
            BondOrder::Up,
            "the bond's own order still faithfully records Bold -> Up; only \
             chirality PERCEPTION abstains, structural fidelity is unconditional"
        );
    }

    /// Negative control: the exact same connectivity/geometry, but with a
    /// real directional "WedgeBegin" instead of "Bold" -- the Begin atom is
    /// unambiguously the reference by CDXML's own convention, so this is
    /// NOT ambiguous and must keep perceiving chirality normally (proves
    /// the downgrade in the test above is specific to non-directional
    /// displays, not a blanket "shared bond between two stereocenters"
    /// suppression that would also break the common, unambiguous case).
    #[test]
    fn wedge_begin_between_two_stereocenters_is_not_ambiguous() {
        let cdxml = r#"<CDXML><fragment>
<n id="1" p="0 0"/>
<n id="2" p="10 0"/>
<n id="3" p="0 10" Element="9"/>
<n id="4" p="-10 0" Element="17"/>
<n id="5" p="0 -10" Element="35"/>
<n id="6" p="20 10" Element="9"/>
<n id="7" p="20 -10" Element="17"/>
<n id="8" p="10 -20" Element="35"/>
<b B="1" E="2" Display="WedgeBegin"/>
<b B="1" E="3"/>
<b B="1" E="4"/>
<b B="1" E="5"/>
<b B="2" E="6"/>
<b B="2" E="7"/>
<b B="2" E="8"/>
</fragment></CDXML>"#;
        let (mol, _) = parse_cdxml(cdxml).unwrap();
        let a1 = chematic_core::AtomIdx(0);
        assert_ne!(
            mol.atom(a1).chirality,
            Chirality::None,
            "a directional WedgeBegin between two stereocenters is not ambiguous \
             and must still perceive chirality"
        );
        assert_eq!(
            mol.bond(chematic_core::BondIdx(0)).order,
            BondOrder::Up,
            "a directional wedge must never be downgraded by the ambiguity check"
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

    #[test]
    fn bounded_parser_rejects_input_line_and_atom_limits() {
        assert!(matches!(
            parse_cdxml_with_limits(
                ETHANOL_CDXML,
                &CdxmlParseLimits {
                    max_input_bytes: 8,
                    ..Default::default()
                }
            ),
            Err(CdxmlError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));
        assert!(matches!(
            parse_cdxml_with_limits(
                ETHANOL_CDXML,
                &CdxmlParseLimits {
                    max_lines: 1,
                    ..Default::default()
                }
            ),
            Err(CdxmlError::ResourceLimit {
                resource: "lines",
                ..
            })
        ));
        assert!(matches!(
            parse_cdxml_with_limits(
                ETHANOL_CDXML,
                &CdxmlParseLimits {
                    max_atoms: 1,
                    ..Default::default()
                }
            ),
            Err(CdxmlError::ResourceLimit {
                resource: "atoms",
                ..
            })
        ));
    }

    #[test]
    fn bounded_parser_rejects_bond_and_fragment_limits() {
        assert!(matches!(
            parse_cdxml_all_with_limits(
                ETHANOL_CDXML,
                &CdxmlParseLimits {
                    max_bonds: 1,
                    ..Default::default()
                }
            ),
            Err(CdxmlError::ResourceLimit {
                resource: "bonds",
                ..
            })
        ));
        assert!(matches!(
            parse_cdxml_all_with_limits(
                TWO_FRAGMENT_CDXML,
                &CdxmlParseLimits {
                    max_fragments: 1,
                    ..Default::default()
                }
            ),
            Err(CdxmlError::ResourceLimit {
                resource: "fragments",
                ..
            })
        ));
    }
}
