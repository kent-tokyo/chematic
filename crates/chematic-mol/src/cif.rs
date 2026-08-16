//! CIF (Crystallographic Information File) parser and writer.
//!
//! ## Scope
//!
//! Extracts **atomic positions** and converts fractional coordinates to
//! orthogonal Ångström coordinates.  Symmetry expansion (space group /
//! symmetry operations) is **not** performed — only atoms listed in the
//! `_atom_site_*` loop are returned (effectively P1 treatment).
//!
//! ## CIF structure
//!
//! ```text
//! data_NaCl
//! _cell_length_a   5.6402
//! _cell_length_b   5.6402
//! _cell_length_c   5.6402
//! _cell_angle_alpha  90.000
//! _cell_angle_beta   90.000
//! _cell_angle_gamma  90.000
//!
//! loop_
//! _atom_site_label
//! _atom_site_type_symbol
//! _atom_site_fract_x
//! _atom_site_fract_y
//! _atom_site_fract_z
//! Na1  Na  0.00000  0.00000  0.00000
//! Cl1  Cl  0.50000  0.50000  0.50000
//! ```
//!
//! ## Fractional → Orthogonal conversion (IUCr convention)
//!
//! ```text
//! X = a·fx + b·cos(γ)·fy + c·cos(β)·fz
//! Y =        b·sin(γ)·fy + c·(cos(α)−cos(β)cos(γ))/sin(γ)·fz
//! Z =                      (V/(a·b·sin(γ)))·fz
//! ```

use chematic_core::{Atom, AtomIdx, Element, Molecule, MoleculeBuilder};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Unit cell parameters (lengths in Å, angles in degrees).
#[derive(Debug, Clone, PartialEq)]
pub struct UnitCell {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
}

impl UnitCell {
    /// Volume in Å³.
    pub fn volume(&self) -> f64 {
        let (ca, cb, cg) = (
            self.alpha.to_radians().cos(),
            self.beta.to_radians().cos(),
            self.gamma.to_radians().cos(),
        );
        self.a * self.b * self.c * (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt()
    }

    /// Convert fractional coordinates to orthogonal Å.
    pub fn frac_to_cart(&self, fx: f64, fy: f64, fz: f64) -> (f64, f64, f64) {
        let ca = self.alpha.to_radians().cos();
        let cb = self.beta.to_radians().cos();
        let sg = self.gamma.to_radians().sin();
        let cg = self.gamma.to_radians().cos();
        let x = self.a * fx + self.b * cg * fy + self.c * cb * fz;
        let y = self.b * sg * fy + self.c * ((ca - cb * cg) / sg) * fz;
        let z = (self.volume() / (self.a * self.b * sg)) * fz;
        (x, y, z)
    }
}

impl Default for UnitCell {
    fn default() -> Self {
        Self {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        }
    }
}

/// Result of parsing a CIF file.
pub struct CifResult {
    /// Molecular topology (atoms; no bonds).
    pub mol: Molecule,
    /// Atomic coordinates in orthogonal Ångströms as `(x, y, z)` tuples.
    pub coords: Vec<(f64, f64, f64)>,
    /// Unit cell parameters, if present in the file.
    pub cell: Option<UnitCell>,
}

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Errors that can occur when parsing CIF files.
#[derive(Debug, Clone, PartialEq)]
pub enum CifError {
    /// No `_atom_site_*` loop with coordinates was found.
    NoAtomSiteLoop,
    /// The atom site loop lacked at least one coordinate column.
    MissingCoordinateColumns,
    /// An element symbol or label could not be resolved.
    UnknownElement(String),
    /// A coordinate value could not be parsed as a float.
    InvalidCoordinate(String),
    /// The unit cell has degenerate angles (e.g. γ = 0° or 180°) that make
    /// the fractional → Cartesian transformation undefined.
    InvalidCellParameters(String),
    /// The atom site loop uses fractional coordinates but the CIF contains no
    /// `_cell_length_*` / `_cell_angle_*` parameters to convert them.
    MissingCellParameters,
}

impl core::fmt::Display for CifError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoAtomSiteLoop => write!(f, "no _atom_site_* loop found in CIF"),
            Self::MissingCoordinateColumns => {
                write!(f, "atom_site loop missing fract_x/y/z columns")
            }
            Self::UnknownElement(s) => write!(f, "unknown element '{s}' in CIF"),
            Self::InvalidCoordinate(s) => write!(f, "invalid coordinate '{s}' in CIF"),
            Self::InvalidCellParameters(s) => write!(f, "invalid cell parameters: {s}"),
            Self::MissingCellParameters => write!(
                f,
                "fractional coordinates present but no _cell_length_*/_cell_angle_* parameters found"
            ),
        }
    }
}

impl std::error::Error for CifError {}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a CIF file and return atoms with orthogonal coordinates.
///
/// Only the first `data_` block is parsed.  Symmetry expansion is **not**
/// performed; returned atoms are exactly those listed in the `_atom_site_*`
/// loop.
/// Strip a CIF comment from one line, respecting single- and double-quoted
/// strings (a `#` inside quotes is not a comment delimiter).
///
/// `pub(crate)`: shared with [`crate::mmcif`], which reuses this crate's
/// existing STAR/CIF tokenizing layer rather than re-implementing it --
/// see that module's top-level docs.
pub(crate) fn strip_cif_comment(line: &str) -> &str {
    let mut in_single = false;
    let mut in_double = false;
    for (i, c) in line.char_indices() {
        match c {
            '\'' if !in_double => in_single = !in_single,
            '"' if !in_single => in_double = !in_double,
            '#' if !in_single && !in_double => return &line[..i],
            _ => {}
        }
    }
    line
}

/// Resolve an element from a raw `_atom_site_type_symbol` or
/// `_atom_site_label` token, stripping trailing digits and oxidation-state
/// signs (e.g. `"Na1"`, `"Cu2+"`, `"Fe3+"`, `"O2-"`).
///
/// Shared by [`parse_cif`]'s row loop and the `crystal`-feature adapter
/// (`parse_cif_periodic_structure`) so both resolve elements identically.
///
/// `pub(crate)`: also reused by [`crate::mmcif`] for `_atom_site.type_symbol`.
pub(crate) fn resolve_element(elem_raw: &str) -> Result<Element, CifError> {
    let elem_str = elem_raw.trim_end_matches(|c: char| c.is_ascii_digit() || c == '+' || c == '-');
    Element::from_symbol(elem_str).ok_or_else(|| CifError::UnknownElement(elem_str.to_string()))
}

/// Which of the six `_cell_length_*`/`_cell_angle_*` tags [`scan_cell`]
/// actually saw while scanning.
#[derive(Debug, Clone, Copy, Default)]
struct CellFieldsSeen {
    a: bool,
    b: bool,
    c: bool,
    alpha: bool,
    beta: bool,
    gamma: bool,
}

impl CellFieldsSeen {
    /// `true` only if every one of the six tags was present.
    ///
    /// Only consumed by the `crystal`-feature adapter today; allowed dead
    /// otherwise rather than gating the whole method behind `#[cfg]`.
    #[cfg_attr(not(feature = "crystal"), allow(dead_code))]
    fn all(&self) -> bool {
        self.a && self.b && self.c && self.alpha && self.beta && self.gamma
    }
}

/// Scan `tokens` for `_cell_length_*`/`_cell_angle_*` tags, building a
/// [`UnitCell`] from [`UnitCell::default`].
///
/// Returns `(cell, has_cell, seen)`. `has_cell` mirrors [`parse_cif`]'s
/// original behavior exactly (`true` iff `_cell_length_a` specifically was
/// seen -- a pre-existing quirk, not normalized here, so as not to change
/// `parse_cif`'s behavior on a CIF with some but not all six tags present).
/// `seen` is the honest per-tag signal: callers that need all six present
/// (e.g. the `crystal` adapter, which cannot build a `Lattice` from a
/// partially-defaulted cell) should check `seen.all()` instead of
/// `has_cell`.
fn scan_cell(tokens: &[String]) -> (UnitCell, bool, CellFieldsSeen) {
    let mut cell = UnitCell::default();
    let mut has_cell = false;
    let mut seen = CellFieldsSeen::default();
    let mut i = 0;
    while i + 1 < tokens.len() {
        match tokens[i].to_ascii_lowercase().as_str() {
            "_cell_length_a" => {
                cell.a = parse_esd(&tokens[i + 1]).unwrap_or(cell.a);
                has_cell = true;
                seen.a = true;
            }
            "_cell_length_b" => {
                cell.b = parse_esd(&tokens[i + 1]).unwrap_or(cell.b);
                seen.b = true;
            }
            "_cell_length_c" => {
                cell.c = parse_esd(&tokens[i + 1]).unwrap_or(cell.c);
                seen.c = true;
            }
            "_cell_angle_alpha" => {
                cell.alpha = parse_esd(&tokens[i + 1]).unwrap_or(cell.alpha);
                seen.alpha = true;
            }
            "_cell_angle_beta" => {
                cell.beta = parse_esd(&tokens[i + 1]).unwrap_or(cell.beta);
                seen.beta = true;
            }
            "_cell_angle_gamma" => {
                cell.gamma = parse_esd(&tokens[i + 1]).unwrap_or(cell.gamma);
                seen.gamma = true;
            }
            _ => {}
        }
        i += 1;
    }
    (cell, has_cell, seen)
}

/// Find the first `loop_` block whose header tags include an
/// `_atom_site_*` tag. Returns `(lowercased column headers, token index
/// where the data rows begin)`.
fn find_atom_site_loop(tokens: &[String]) -> Result<(Vec<String>, usize), CifError> {
    let mut i = 0;
    while i < tokens.len() {
        if tokens[i].as_str() == "loop_" {
            let mut j = i + 1;
            let mut headers: Vec<String> = Vec::new();
            while j < tokens.len() && tokens[j].starts_with('_') {
                headers.push(tokens[j].to_ascii_lowercase());
                j += 1;
            }
            if headers.iter().any(|h| h.starts_with("_atom_site_")) {
                if headers.is_empty() {
                    return Err(CifError::NoAtomSiteLoop);
                }
                return Ok((headers, j));
            }
            i = j;
        } else {
            i += 1;
        }
    }
    Err(CifError::NoAtomSiteLoop)
}

pub fn parse_cif(input: &str) -> Result<CifResult, CifError> {
    // Strip CIF comments (# to end of line), but not '#' inside quoted strings.
    let clean: String = input
        .lines()
        .map(strip_cif_comment)
        .collect::<Vec<_>>()
        .join("\n");

    let tokens = tokenize_cif(&clean);

    // --- Extract cell parameters ---
    let (cell, has_cell, _cell_fields_seen) = scan_cell(&tokens);

    // --- Find _atom_site loop ---
    let (col_headers, data_start) = find_atom_site_loop(&tokens)?;
    let ncols = col_headers.len();

    let col = |name: &str| -> Option<usize> { col_headers.iter().position(|h| h.as_str() == name) };

    let col_type = col("_atom_site_type_symbol");
    let col_label = col("_atom_site_label");
    // Prefer fractional; fall back to Cartesian.
    let use_cartesian = col("_atom_site_fract_x").is_none();
    let col_x = col("_atom_site_fract_x").or_else(|| col("_atom_site_cartn_x"));
    let col_y = col("_atom_site_fract_y").or_else(|| col("_atom_site_cartn_y"));
    let col_z = col("_atom_site_fract_z").or_else(|| col("_atom_site_cartn_z"));

    let (col_x, col_y, col_z) = match (col_x, col_y, col_z) {
        (Some(x), Some(y), Some(z)) => (x, y, z),
        _ => return Err(CifError::MissingCoordinateColumns),
    };

    // Validate that fractional→Cartesian conversion is well-defined.
    if !use_cartesian {
        if !has_cell {
            return Err(CifError::MissingCellParameters);
        }
        let sg = cell.gamma.to_radians().sin();
        if sg.abs() < 1e-10 {
            return Err(CifError::InvalidCellParameters(format!(
                "_cell_angle_gamma = {} makes sin(γ) ≈ 0, transformation undefined",
                cell.gamma
            )));
        }
        // Also guard against a non-physical cell volume (≤ 0 or NaN).
        let vol = cell.volume();
        if !vol.is_finite() || vol <= 0.0 {
            return Err(CifError::InvalidCellParameters(format!(
                "unit cell volume is non-positive ({vol:.6} Å³); check _cell_angle_* values"
            )));
        }
    }

    let mut builder = MoleculeBuilder::new();
    let mut coords: Vec<(f64, f64, f64)> = Vec::new();
    let data_tokens = &tokens[data_start..];

    let mut row = 0;
    while (row + 1) * ncols <= data_tokens.len() {
        let base = row * ncols;
        let tok = &data_tokens[base..base + ncols];

        // Stop at next loop_ or data_ block.
        if tok[0].as_str() == "loop_" || tok[0].starts_with("data_") {
            break;
        }

        // Resolve element.
        let elem_raw: &str = col_type
            .and_then(|c| tok.get(c))
            .map(|s| s.as_str())
            .or_else(|| col_label.and_then(|c| tok.get(c)).map(|s| s.as_str()))
            .unwrap_or("X");
        let elem = resolve_element(elem_raw)?;

        let parse_coord = |s: &str| -> Result<f64, CifError> {
            strip_esd(s)
                .parse::<f64>()
                .map_err(|_| CifError::InvalidCoordinate(s.to_string()))
        };
        let fx = parse_coord(&tok[col_x])?;
        let fy = parse_coord(&tok[col_y])?;
        let fz = parse_coord(&tok[col_z])?;

        let (x, y, z) = if use_cartesian {
            (fx, fy, fz)
        } else {
            cell.frac_to_cart(fx, fy, fz)
        };

        builder.add_atom(Atom::new(elem));
        coords.push((x, y, z));
        row += 1;
    }

    if coords.is_empty() {
        return Err(CifError::NoAtomSiteLoop);
    }

    Ok(CifResult {
        mol: builder.build(),
        coords,
        cell: if has_cell { Some(cell) } else { None },
    })
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Write a minimal CIF file (P1, no symmetry operators).
///
/// When `cell` is `None`, fractional coordinates are written as-is
/// (the coordinates are treated as already orthogonal with a=b=c=1).
pub fn write_cif(mol: &Molecule, coords: &[(f64, f64, f64)], cell: Option<&UnitCell>) -> String {
    let mut out = String::from("data_chematic\n\n");

    if let Some(c) = cell {
        out.push_str(&format!(
            "_cell_length_a   {:.4}\n\
             _cell_length_b   {:.4}\n\
             _cell_length_c   {:.4}\n\
             _cell_angle_alpha  {:.3}\n\
             _cell_angle_beta   {:.3}\n\
             _cell_angle_gamma  {:.3}\n\n\
             _symmetry_space_group_name_H-M  'P 1'\n\n",
            c.a, c.b, c.c, c.alpha, c.beta, c.gamma
        ));
    }

    out.push_str(
        "loop_\n\
         _atom_site_label\n\
         _atom_site_type_symbol\n\
         _atom_site_fract_x\n\
         _atom_site_fract_y\n\
         _atom_site_fract_z\n",
    );

    let fallback_cell = UnitCell::default();
    let cell_ref = cell.unwrap_or(&fallback_cell);

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let sym = mol.atom(idx).element.symbol();
        let (x, y, z) = coords.get(i).copied().unwrap_or((0.0, 0.0, 0.0));
        let (fx, fy, fz) = cart_to_frac(cell_ref, x, y, z);
        out.push_str(&format!(
            "{sym}{} {sym}  {fx:.5}  {fy:.5}  {fz:.5}\n",
            i + 1
        ));
    }
    out
}

fn cart_to_frac(cell: &UnitCell, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    // Exact inverse for orthogonal cells; approximate for triclinic.
    let sg = cell.gamma.to_radians().sin();
    let cg = cell.gamma.to_radians().cos();
    let cb = cell.beta.to_radians().cos();
    let ca = cell.alpha.to_radians().cos();
    let fz = z * cell.a * cell.b * sg / cell.volume();
    let fy = (y - cell.c * ((ca - cb * cg) / sg) * fz) / (cell.b * sg);
    let fx = (x - cell.b * cg * fy - cell.c * cb * fz) / cell.a;
    (fx, fy, fz)
}

// ---------------------------------------------------------------------------
// CIF tokenizer
// ---------------------------------------------------------------------------

/// `pub(crate)`: shared low-level STAR/CIF tokenizer, reused by
/// [`crate::mmcif`] rather than re-implemented (mmCIF is the same
/// underlying STAR syntax, just a different category-tag convention --
/// `_atom_site.Cartn_x` vs. this module's `_atom_site_fract_x`).
pub(crate) fn tokenize_cif(input: &str) -> Vec<String> {
    let mut tokens: Vec<String> = Vec::new();
    let chars: Vec<char> = input.chars().collect();
    let len = chars.len();
    let mut i = 0;

    while i < len {
        let ch = chars[i];
        // Skip whitespace.
        if ch.is_whitespace() {
            i += 1;
            continue;
        }

        // Semicolon text block (starts at column 0 on a new line).
        if ch == ';' && (i == 0 || chars[i - 1] == '\n') {
            i += 1;
            let start = i;
            while i < len {
                if chars[i] == '\n' && i + 1 < len && chars[i + 1] == ';' {
                    break;
                }
                i += 1;
            }
            // Skip the text content (we don't use it for atom sites).
            tokens.push(chars[start..i].iter().collect());
            i += 2; // skip '\n;'
            continue;
        }

        // Single-quoted string.
        if ch == '\'' {
            i += 1;
            let start = i;
            while i < len && chars[i] != '\'' {
                i += 1;
            }
            tokens.push(chars[start..i].iter().collect());
            if i < len {
                i += 1;
            }
            continue;
        }

        // Double-quoted string.
        if ch == '"' {
            i += 1;
            let start = i;
            while i < len && chars[i] != '"' {
                i += 1;
            }
            tokens.push(chars[start..i].iter().collect());
            if i < len {
                i += 1;
            }
            continue;
        }

        // Regular token.
        let start = i;
        while i < len && !chars[i].is_whitespace() {
            i += 1;
        }
        tokens.push(chars[start..i].iter().collect());
    }
    tokens
}

/// `pub(crate)`: shared with [`crate::mmcif`] (esd-suffix stripping is a
/// general CIF numeric-literal rule, not specific to small-molecule tags).
pub(crate) fn strip_esd(s: &str) -> &str {
    s.find('(').map_or(s, |pos| &s[..pos])
}

pub(crate) fn parse_esd(s: &str) -> Option<f64> {
    strip_esd(s).parse::<f64>().ok()
}

// ---------------------------------------------------------------------------
// `crystal` feature: PeriodicStructure adapter
// ---------------------------------------------------------------------------
//
// Bridges this module's existing tokenizer/cell-scanner/atom-site-loop
// parsing (above) to `chematic_crystal::PeriodicStructure`. Reuses
// `resolve_element`/`scan_cell`/`find_atom_site_loop`/`strip_cif_comment`/
// `tokenize_cif`/`strip_esd` rather than re-parsing CIF text -- see
// docs/rfcs/chematic_crystal_foundation.md's "CIF migration" section for
// the original sketch this implements. `parse_cif`/`write_cif`/`UnitCell`
// above are untouched by this section.

#[cfg(feature = "crystal")]
mod crystal_adapter {
    use super::{
        CellFieldsSeen, CifError, UnitCell, find_atom_site_loop, resolve_element, scan_cell,
        strip_cif_comment, strip_esd, tokenize_cif,
    };
    use crate::cif_symmetry::{self, CifSymmetryError, sites_within_tolerance};
    use chematic_core::Element;
    use chematic_crystal::{
        CartesianCoord, CrystalError, FractionalCoord, Lattice, Occupancy, PeriodicSite,
        PeriodicStructure, SiteSpecies,
    };

    /// How a CIF's declared symmetry relates to the sites
    /// [`parse_cif_periodic_structure`] (or
    /// [`parse_cif_periodic_structure_with_options`]) returns.
    ///
    /// **What `ExpandedExplicitOperations` does *not* claim**: it means
    /// only that every symmetry operation *literally written* in the CIF's
    /// own `_space_group_symop_operation_xyz` / `_symmetry_equiv_pos_as_xyz`
    /// loop (or dotted-tag alias) was applied to the asymmetric unit. The
    /// CIF specification technically requires that explicit list to be a
    /// *complete* set for the declared space group, but this implementation
    /// never cross-checks that against any space-group name, International
    /// Tables number, or external database -- it has no space-group
    /// database and does not generate operations from a name/number (see
    /// `crates/chematic-mol/src/cif_symmetry.rs`'s module docs and
    /// `docs/crystal_scope.md`). A file whose explicit list happens to be
    /// incomplete or wrong for its stated space group will still expand
    /// "successfully" against exactly the operations it wrote down.
    /// **"Expanded" here is a faithfulness claim about the CIF's own text,
    /// not a correctness/completeness claim against crystallographic
    /// convention.**
    #[derive(Debug, Clone, PartialEq)]
    pub enum CifSymmetryStatus {
        /// No symmetry beyond P1 was declared (or nothing was declared at
        /// all), and no explicit operation list was present: as far as
        /// this file states, the returned sites are already the complete
        /// cell content.
        P1,
        /// Every symmetry operation explicitly listed in the CIF's own
        /// operation-expression tag (loop or single standalone item) was
        /// successfully parsed and applied to the asymmetric unit -- see
        /// this type's own docs above for what this does and does not
        /// claim.
        ExpandedExplicitOperations {
            /// `_symmetry_space_group_name_H-M` / `_space_group_name_H-M_alt`
            /// value, if present.
            space_group_name: Option<String>,
            /// Number of operations in the (parsed, validated) explicit
            /// list that was applied.
            operation_count: usize,
            /// Number of sites in the CIF's original asymmetric unit,
            /// before expansion.
            asymmetric_site_count: usize,
            /// Number of sites in the expanded structure, after
            /// special-position deduplication (`<= operation_count *
            /// asymmetric_site_count`).
            expanded_site_count: usize,
        },
        /// The file declares symmetry beyond P1 that this adapter did
        /// **not** expand -- either because expansion was turned off
        /// (`CifPeriodicParseOptions::expand_explicit_symmetry = false`) or
        /// because there was no parseable explicit operation list to
        /// expand from (a space-group name/number alone, with no operation
        /// loop, is never enough -- see this crate's scope: no
        /// space-group-database lookup, no name/number-to-operations
        /// generation). The returned sites are only the raw asymmetric
        /// unit as literally listed, not a full unit cell.
        UnexpandedSymmetry {
            /// `_symmetry_space_group_name_H-M` / `_space_group_name_H-M_alt`
            /// value, if present (CIF-quoting already stripped by the
            /// tokenizer; not otherwise normalized).
            space_group_name: Option<String>,
            /// Row count of a `_space_group_symop_operation_xyz` /
            /// `_symmetry_equiv_pos_as_xyz` loop (or dotted-tag alias), if
            /// one was present (`0` if none was found -- e.g. a file
            /// naming its space group by number only).
            operation_count: usize,
        },
    }

    impl CifSymmetryStatus {
        /// `true` for [`Self::P1`] and [`Self::ExpandedExplicitOperations`]
        /// -- the two statuses under which the returned sites are a
        /// genuinely complete unit cell, safe to write back out as a
        /// (nominally P1) CIF via
        /// [`CifPeriodicResult::to_cif_checked`]/[`write_cif_periodic_structure`].
        /// `false` for [`Self::UnexpandedSymmetry`], where the sites are
        /// only an asymmetric unit.
        pub fn is_complete_cell(&self) -> bool {
            matches!(self, Self::P1 | Self::ExpandedExplicitOperations { .. })
        }
    }

    /// Options for [`parse_cif_periodic_structure_with_options`].
    #[derive(Debug, Clone, Copy)]
    pub struct CifPeriodicParseOptions {
        /// If `true` (the default via [`parse_cif_periodic_structure`]),
        /// a CIF with a parseable explicit symmetry-operation list is
        /// expanded into a full unit cell. If `false`, the operation list
        /// (if any) is only counted, never parsed/validated/applied -- the
        /// asymmetric unit is returned as-is with
        /// [`CifSymmetryStatus::UnexpandedSymmetry`], matching this
        /// adapter's pre-expansion behavior exactly.
        pub expand_explicit_symmetry: bool,
    }

    /// Result of [`parse_cif_periodic_structure`] /
    /// [`parse_cif_periodic_structure_with_options`].
    #[derive(Debug)]
    pub struct CifPeriodicResult {
        /// The parsed periodic structure (lattice + sites -- see
        /// [`CifSymmetryStatus`] for whether this is a full cell or only
        /// an asymmetric unit).
        pub structure: PeriodicStructure,
        /// Whether/how the CIF's declared symmetry was expanded.
        pub symmetry: CifSymmetryStatus,
    }

    impl CifPeriodicResult {
        /// Write as CIF text, refusing when `self.symmetry` is
        /// [`CifSymmetryStatus::UnexpandedSymmetry`] (only an asymmetric
        /// unit -- writing it as a nominal `P 1` cell would falsely
        /// declare it complete). Both [`CifSymmetryStatus::P1`] and
        /// [`CifSymmetryStatus::ExpandedExplicitOperations`] write freely:
        /// in both cases the sites genuinely are a complete cell, even
        /// though [`write_cif_periodic_structure`] always declares the
        /// output's space group as a nominal `P 1` (that is correct here
        /// -- the *geometry* is complete, which is the property that
        /// matters for round-tripping; the output not re-stating the
        /// original space group name/operations is a pre-existing,
        /// documented writer limitation, not new to this method).
        ///
        /// This is the single source of truth for the write-safety
        /// judgment -- `chematic-py`'s `PeriodicStructure.to_cif()`
        /// delegates here rather than re-implementing the check.
        pub fn to_cif_checked(&self) -> Result<String, CifPeriodicError> {
            match &self.symmetry {
                CifSymmetryStatus::P1 | CifSymmetryStatus::ExpandedExplicitOperations { .. } => {
                    Ok(write_cif_periodic_structure(&self.structure))
                }
                CifSymmetryStatus::UnexpandedSymmetry {
                    space_group_name,
                    operation_count,
                } => Err(CifPeriodicError::UnexpandedSymmetryWrite {
                    space_group_name: space_group_name.clone(),
                    operation_count: *operation_count,
                }),
            }
        }
    }

    /// Errors from [`parse_cif_periodic_structure`] /
    /// [`parse_cif_periodic_structure_with_options`] /
    /// [`CifPeriodicResult::to_cif_checked`].
    #[derive(Debug, Clone, PartialEq)]
    pub enum CifPeriodicError {
        /// A CIF-level parsing problem (same errors [`parse_cif`](super::parse_cif)
        /// can produce: no atom_site loop, unknown element, bad
        /// coordinate, etc.).
        Cif(CifError),
        /// A `chematic-crystal` validation failure -- non-finite/degenerate
        /// cell, empty species list, occupancy-sum exceeded, etc. Reused
        /// as-is from `chematic_crystal::CrystalError` rather than
        /// re-validated here, so there is exactly one place that decides
        /// what a valid `PeriodicStructure` is.
        Crystal(CrystalError),
        /// An `_atom_site_occupancy` value was present but was neither a
        /// CIF placeholder (`"."`/`"?"`, both treated as full occupancy)
        /// nor parseable as a number.
        InvalidOccupancy(String),
        /// A problem parsing, validating, or applying the CIF's explicit
        /// symmetry-operation list -- see [`CifSymmetryError`] for the
        /// specific failure modes (malformed expression, invalid rotation
        /// matrix, duplicate/missing/conflicting operations, a
        /// special-position collision between differently-composed sites,
        /// or an expansion too large to safely attempt).
        Symmetry(CifSymmetryError),
        /// [`CifPeriodicResult::to_cif_checked`] refused to write: the
        /// structure holds only the asymmetric unit of an unexpanded
        /// space group.
        UnexpandedSymmetryWrite {
            space_group_name: Option<String>,
            operation_count: usize,
        },
    }

    impl core::fmt::Display for CifPeriodicError {
        fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            match self {
                Self::Cif(e) => write!(f, "{e}"),
                Self::Crystal(e) => write!(f, "{e}"),
                Self::InvalidOccupancy(s) => write!(f, "invalid occupancy '{s}' in CIF"),
                Self::Symmetry(e) => write!(f, "{e}"),
                Self::UnexpandedSymmetryWrite {
                    space_group_name,
                    operation_count,
                } => write!(
                    f,
                    "cannot write CIF: this structure holds only the asymmetric unit of an \
                     unexpanded space group (space_group_name={space_group_name:?}, \
                     {operation_count} symmetry operation(s) declared in the source CIF) -- \
                     writing it back would falsely declare a complete P1 cell; expand the \
                     symmetry first (or re-parse with expand_explicit_symmetry=true) if you \
                     need a full-cell CIF"
                ),
            }
        }
    }

    impl std::error::Error for CifPeriodicError {}

    impl From<CifError> for CifPeriodicError {
        fn from(e: CifError) -> Self {
            Self::Cif(e)
        }
    }

    impl From<CrystalError> for CifPeriodicError {
        fn from(e: CrystalError) -> Self {
            Self::Crystal(e)
        }
    }

    impl From<CifSymmetryError> for CifPeriodicError {
        fn from(e: CifSymmetryError) -> Self {
            Self::Symmetry(e)
        }
    }

    /// Parse an `_atom_site_occupancy` token. `"."`/`"?"` (CIF's own
    /// "inapplicable"/"unknown" placeholders) map to full occupancy
    /// (`1.0`); anything else must parse as a real number -- deliberately
    /// not defaulting a malformed value to `1.0`, since that would
    /// silently misreport the source data as fully occupied.
    fn parse_occupancy_token(s: &str) -> Result<f64, CifPeriodicError> {
        if s == "." || s == "?" {
            return Ok(1.0);
        }
        strip_esd(s)
            .parse::<f64>()
            .map_err(|_| CifPeriodicError::InvalidOccupancy(s.to_string()))
    }

    /// Scan `tokens` for `_symmetry_space_group_name_H-M` (or its mmCIF
    /// alternate tag) and `_symmetry_Int_Tables_number` (or
    /// `_space_group_IT_number`) only -- the explicit operation list itself
    /// is scanned/parsed separately, by
    /// [`cif_symmetry::scan_operation_sources`]/[`cif_symmetry::resolve_symmetry_operations`].
    ///
    /// # Known gap
    ///
    /// ponytail: an *unquoted* space-group name (e.g. bare `P 1` with no
    /// surrounding `'`/`"`) truncates at the first whitespace under this
    /// crate's tokenizer (whitespace-delimited outside quotes), so it
    /// would read as `"P"` and get misclassified as non-P1. Real CIFs
    /// quote this value (and this module's own writer does too), so this
    /// is a narrow, known gap rather than a silent one -- widen if an
    /// unquoted fixture ever surfaces.
    fn scan_symmetry_name_and_number(tokens: &[String]) -> (Option<String>, Option<i64>) {
        let mut space_group_name: Option<String> = None;
        let mut it_number: Option<i64> = None;
        let mut i = 0;
        while i + 1 < tokens.len() {
            match tokens[i].to_ascii_lowercase().as_str() {
                "_symmetry_space_group_name_h-m" | "_space_group_name_h-m_alt" => {
                    space_group_name = Some(tokens[i + 1].clone());
                }
                "_symmetry_int_tables_number" | "_space_group_it_number" => {
                    it_number = tokens[i + 1].trim().parse::<i64>().ok();
                }
                _ => {}
            }
            i += 1;
        }
        (space_group_name, it_number)
    }

    /// Classify as [`CifSymmetryStatus::P1`] or
    /// [`CifSymmetryStatus::UnexpandedSymmetry`] from the declared
    /// name/number plus a raw operation-row count, **without** attempting
    /// to parse/validate the operations themselves. Used both when the
    /// caller opted out of expansion
    /// (`CifPeriodicParseOptions::expand_explicit_symmetry = false`) and
    /// when expansion was on but no explicit operation list was found at
    /// all (`operation_count == 0` in that case, so this reduces to the
    /// pure name/number heuristic).
    ///
    /// Any operation list at all (`operation_count >= 1`), regardless of
    /// its own content, is treated as non-P1 -- deliberately dropping this
    /// function's predecessor's old quirk of treating a lone listed
    /// operator as P1 unconditionally (i.e. even if that one operation
    /// were *not* the identity). That old quirk only mattered pre-expansion;
    /// now that `expand_explicit_symmetry = true` (the default) fully
    /// parses and expands any explicit list -- including a legitimate
    /// single-identity-operation list, which just expands trivially to
    /// itself -- this classification-only path exists purely for the
    /// opt-out case and should report "there is a declared list I did not
    /// expand" honestly, not re-derive expansion's own judgment about
    /// whether that list happens to be trivial.
    fn classify_p1_or_unexpanded(
        space_group_name: Option<String>,
        it_number: Option<i64>,
        operation_count: usize,
    ) -> CifSymmetryStatus {
        let is_p1_name = space_group_name
            .as_deref()
            .map(|s| {
                let norm: String = s.chars().filter(|c| !c.is_whitespace()).collect();
                norm.eq_ignore_ascii_case("p1")
            })
            .unwrap_or(true); // absent tag doesn't itself claim non-P1
        let is_p1_number = it_number.map(|n| n == 1).unwrap_or(true);

        if !is_p1_name || !is_p1_number || operation_count >= 1 {
            CifSymmetryStatus::UnexpandedSymmetry {
                space_group_name,
                operation_count,
            }
        } else {
            CifSymmetryStatus::P1
        }
    }

    /// One `_atom_site_*` row's data, before disorder grouping.
    struct RawAtomRow {
        label: Option<String>,
        element: Element,
        fractional: FractionalCoord,
        occupancy: f64,
    }

    /// Parse a CIF file into a [`PeriodicStructure`], reusing this module's
    /// existing tokenizer/cell-scanner/atom-site-loop-finder.
    ///
    /// Additionally reads `_atom_site_occupancy` (defaulting to `1.0` when
    /// absent) and groups rows that share a fractional position (within a
    /// small internal tolerance) into one [`PeriodicSite`] with multiple
    /// [`SiteSpecies`] -- the CIF convention for modeling positional
    /// disorder.
    ///
    /// Requires all six `_cell_length_*`/`_cell_angle_*` tags to build a
    /// [`Lattice`] (unlike [`parse_cif`](super::parse_cif), which silently
    /// defaults any tag it didn't see to `1.0`/`90.0` -- a
    /// `PeriodicStructure` cannot tolerate a partially-defaulted cell the
    /// way a plain `Molecule` + coordinate list can).
    ///
    /// Thin wrapper around [`parse_cif_periodic_structure_with_options`]
    /// with `expand_explicit_symmetry: true` -- the recommended default.
    /// See [`CifSymmetryStatus`] for what "expanded" does and does not
    /// claim.
    ///
    /// # Examples
    ///
    /// ```
    /// use chematic_mol::cif::parse_cif_periodic_structure;
    ///
    /// let cif = "data_NaCl\n\
    ///     _cell_length_a 5.6402\n_cell_length_b 5.6402\n_cell_length_c 5.6402\n\
    ///     _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
    ///     loop_\n_atom_site_label\n_atom_site_type_symbol\n\
    ///     _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
    ///     Na1 Na 0.0 0.0 0.0\nCl1 Cl 0.5 0.5 0.5\n";
    /// let result = parse_cif_periodic_structure(cif).unwrap();
    /// assert_eq!(result.structure.site_count(), 2);
    /// ```
    pub fn parse_cif_periodic_structure(
        input: &str,
    ) -> Result<CifPeriodicResult, CifPeriodicError> {
        parse_cif_periodic_structure_with_options(
            input,
            CifPeriodicParseOptions {
                expand_explicit_symmetry: true,
            },
        )
    }

    /// [`parse_cif_periodic_structure`] with explicit control over whether
    /// a parseable symmetry-operation list gets expanded (see
    /// [`CifPeriodicParseOptions`]).
    ///
    /// Only symmetry operations **literally written** in the CIF's own
    /// `_space_group_symop_operation_xyz` / `_symmetry_equiv_pos_as_xyz`
    /// tag (or dotted-tag alias) are ever applied -- this never generates
    /// operations from a space-group name or International Tables number,
    /// and never consults any space-group database (see
    /// `crates/chematic-mol/src/cif_symmetry.rs` and
    /// `docs/crystal_scope.md`).
    ///
    /// # Examples
    ///
    /// A CIF with an explicit inversion center expands to 2 sites:
    ///
    /// ```
    /// use chematic_mol::cif::{
    ///     CifPeriodicParseOptions, CifSymmetryStatus, parse_cif_periodic_structure_with_options,
    /// };
    ///
    /// let cif = "data_x\n\
    ///     _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
    ///     _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
    ///     loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'\n\
    ///     loop_\n_atom_site_label\n_atom_site_type_symbol\n\
    ///     _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
    ///     C1 C 0.1 0.2 0.3\n";
    /// let result = parse_cif_periodic_structure_with_options(
    ///     cif,
    ///     CifPeriodicParseOptions { expand_explicit_symmetry: true },
    /// ).unwrap();
    /// assert_eq!(result.structure.site_count(), 2);
    /// assert!(matches!(result.symmetry, CifSymmetryStatus::ExpandedExplicitOperations { .. }));
    /// ```
    pub fn parse_cif_periodic_structure_with_options(
        input: &str,
        options: CifPeriodicParseOptions,
    ) -> Result<CifPeriodicResult, CifPeriodicError> {
        let clean: String = input
            .lines()
            .map(strip_cif_comment)
            .collect::<Vec<_>>()
            .join("\n");
        let tokens = tokenize_cif(&clean);

        let (cell, _has_cell, seen): (UnitCell, bool, CellFieldsSeen) = scan_cell(&tokens);
        if !seen.all() {
            return Err(CifPeriodicError::Cif(CifError::MissingCellParameters));
        }
        let lattice =
            Lattice::from_parameters(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)?;

        let (space_group_name, it_number) = scan_symmetry_name_and_number(&tokens);

        let (col_headers, data_start) = find_atom_site_loop(&tokens)?;
        let ncols = col_headers.len();
        let col =
            |name: &str| -> Option<usize> { col_headers.iter().position(|h| h.as_str() == name) };
        let col_type = col("_atom_site_type_symbol");
        let col_label = col("_atom_site_label");
        let use_cartesian = col("_atom_site_fract_x").is_none();
        let col_x = col("_atom_site_fract_x").or_else(|| col("_atom_site_cartn_x"));
        let col_y = col("_atom_site_fract_y").or_else(|| col("_atom_site_cartn_y"));
        let col_z = col("_atom_site_fract_z").or_else(|| col("_atom_site_cartn_z"));
        let col_occ = col("_atom_site_occupancy");
        let (col_x, col_y, col_z) = match (col_x, col_y, col_z) {
            (Some(x), Some(y), Some(z)) => (x, y, z),
            _ => return Err(CifPeriodicError::Cif(CifError::MissingCoordinateColumns)),
        };

        let data_tokens = &tokens[data_start..];
        let mut rows: Vec<RawAtomRow> = Vec::new();
        let mut row_idx = 0;
        while (row_idx + 1) * ncols <= data_tokens.len() {
            let base = row_idx * ncols;
            let tok = &data_tokens[base..base + ncols];
            if tok[0].as_str() == "loop_" || tok[0].starts_with("data_") {
                break;
            }

            let label = col_label.and_then(|c| tok.get(c)).map(|s| s.to_string());
            let elem_raw: &str = col_type
                .and_then(|c| tok.get(c))
                .map(|s| s.as_str())
                .or_else(|| col_label.and_then(|c| tok.get(c)).map(|s| s.as_str()))
                .unwrap_or("X");
            let element = resolve_element(elem_raw)?;

            let parse_coord = |s: &str| -> Result<f64, CifError> {
                strip_esd(s)
                    .parse::<f64>()
                    .map_err(|_| CifError::InvalidCoordinate(s.to_string()))
            };
            let v0 = parse_coord(&tok[col_x])?;
            let v1 = parse_coord(&tok[col_y])?;
            let v2 = parse_coord(&tok[col_z])?;

            let fractional = if use_cartesian {
                lattice.cart_to_frac(CartesianCoord::new([v0, v1, v2]))
            } else {
                FractionalCoord::new([v0, v1, v2])
            };

            let occupancy = col_occ
                .and_then(|c| tok.get(c))
                .map(|s| parse_occupancy_token(s))
                .transpose()?
                .unwrap_or(1.0);

            rows.push(RawAtomRow {
                label,
                element,
                fractional,
                occupancy,
            });
            row_idx += 1;
        }

        if rows.is_empty() {
            return Err(CifPeriodicError::Cif(CifError::NoAtomSiteLoop));
        }

        // Group rows sharing a fractional position (within tolerance) into
        // one site with multiple species -- CIF's convention for
        // positional/substitutional disorder. First-seen position sets
        // site order and label; later rows at that position only add a
        // species.
        struct Grouped {
            fractional: FractionalCoord,
            label: Option<String>,
            species: Vec<SiteSpecies>,
        }
        let mut grouped: Vec<Grouped> = Vec::new();
        for raw in rows {
            // ponytail: O(n) scan per row (O(n^2) total) via
            // sites_within_tolerance's minimum_image call -- fine at the
            // tens-to-low-hundreds row counts real CIFs have; see
            // cif_symmetry::SITE_MERGE_TOLERANCE_ANGSTROM's docs for the
            // unified tolerance this and symmetry-expansion dedup both use.
            let existing = grouped
                .iter_mut()
                .find(|g| sites_within_tolerance(&lattice, g.fractional, raw.fractional));
            let species = SiteSpecies {
                element: raw.element,
                occupancy: Occupancy::new(raw.occupancy)?,
            };
            match existing {
                Some(g) => g.species.push(species),
                None => grouped.push(Grouped {
                    fractional: raw.fractional,
                    label: raw.label,
                    species: vec![species],
                }),
            }
        }

        let asymmetric_sites = grouped
            .into_iter()
            .map(|g| PeriodicSite::new(g.species, g.fractional, g.label))
            .collect::<Result<Vec<_>, CrystalError>>()?;

        if options.expand_explicit_symmetry {
            if let Some(operations) = cif_symmetry::resolve_symmetry_operations(&tokens)? {
                let expanded_sites =
                    cif_symmetry::expand_sites(&asymmetric_sites, &operations, &lattice)?;
                let structure = PeriodicStructure::new(lattice, expanded_sites)?;
                let symmetry = CifSymmetryStatus::ExpandedExplicitOperations {
                    space_group_name,
                    operation_count: operations.len(),
                    asymmetric_site_count: asymmetric_sites.len(),
                    expanded_site_count: structure.site_count(),
                };
                return Ok(CifPeriodicResult {
                    structure,
                    symmetry,
                });
            }
            // No parseable explicit operation list -- fall through to the
            // plain name/number classification (operation_count = 0).
            let structure = PeriodicStructure::new(lattice, asymmetric_sites)?;
            let symmetry = classify_p1_or_unexpanded(space_group_name, it_number, 0);
            return Ok(CifPeriodicResult {
                structure,
                symmetry,
            });
        }

        // expand_explicit_symmetry = false: count (but never parse/validate)
        // any explicit operation list, matching this adapter's pre-expansion
        // behavior exactly.
        let operation_count = cif_symmetry::scan_operation_sources(&tokens)
            .first()
            .map(Vec::len)
            .unwrap_or(0);
        let structure = PeriodicStructure::new(lattice, asymmetric_sites)?;
        let symmetry = classify_p1_or_unexpanded(space_group_name, it_number, operation_count);
        Ok(CifPeriodicResult {
            structure,
            symmetry,
        })
    }

    /// Write a minimal CIF file (P1, no symmetry operators) from a
    /// [`PeriodicStructure`].
    ///
    /// Mirrors [`write_cif`](super::write_cif)'s shape (same cell-parameter
    /// tags, same `'P 1'` declaration) but reads a structure's lattice and
    /// sites directly instead of a `Molecule` + coordinate list, and
    /// additionally emits `_atom_site_occupancy`. A site with more than one
    /// [`SiteSpecies`] (disorder) is written as one row per species, all
    /// sharing that site's fractional position; when a site has more than
    /// one species, each row's label gets a letter suffix (`A`, `B`, ...)
    /// appended so every emitted `_atom_site_label` is distinct, even if
    /// the site itself carries no label (in which case a
    /// `"{symbol}{site_number}"` label is generated).
    ///
    /// This writer never expands or claims symmetry beyond P1 -- it writes
    /// exactly the sites `structure` contains, which is honest regardless
    /// of how those sites were originally obtained (e.g. round-tripping a
    /// [`CifSymmetryStatus::P1`] structure from
    /// [`parse_cif_periodic_structure`] is lossless; round-tripping an
    /// [`CifSymmetryStatus::UnexpandedSymmetry`] one faithfully preserves
    /// only the asymmetric unit that was actually parsed, same as before).
    pub fn write_cif_periodic_structure(structure: &PeriodicStructure) -> String {
        let mut out = String::from("data_chematic\n\n");

        let [a, b, c] = structure.lattice().lengths();
        let [alpha, beta, gamma] = structure.lattice().angles_degrees();
        out.push_str(&format!(
            "_cell_length_a   {a:.4}\n\
             _cell_length_b   {b:.4}\n\
             _cell_length_c   {c:.4}\n\
             _cell_angle_alpha  {alpha:.3}\n\
             _cell_angle_beta   {beta:.3}\n\
             _cell_angle_gamma  {gamma:.3}\n\n\
             _symmetry_space_group_name_H-M  'P 1'\n\n"
        ));

        out.push_str(
            "loop_\n\
             _atom_site_label\n\
             _atom_site_type_symbol\n\
             _atom_site_fract_x\n\
             _atom_site_fract_y\n\
             _atom_site_fract_z\n\
             _atom_site_occupancy\n",
        );

        for (site_index, site) in structure.sites().iter().enumerate() {
            let [fx, fy, fz] = site.fractional.0;
            let multi_species = site.species.len() > 1;
            for (species_index, species) in site.species.iter().enumerate() {
                let sym = species.element.symbol();
                let base_label = site
                    .label
                    .clone()
                    .unwrap_or_else(|| format!("{sym}{}", site_index + 1));
                let label = if multi_species {
                    // ponytail: A-Z only (26 species at one site is not a
                    // realistic disorder case -- occupancies must sum to
                    // <=~1.0, so this ceiling is never reached in practice).
                    let suffix = (b'A' + (species_index as u8 % 26)) as char;
                    format!("{base_label}{suffix}")
                } else {
                    base_label
                };
                out.push_str(&format!(
                    "{label} {sym}  {fx:.5}  {fy:.5}  {fz:.5}  {:.4}\n",
                    species.occupancy.value()
                ));
            }
        }
        out
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        const NACL_P1_CIF: &str = "data_NaCl\n\
            _cell_length_a 5.6402\n_cell_length_b 5.6402\n_cell_length_c 5.6402\n\
            _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
            loop_\n_atom_site_label\n_atom_site_type_symbol\n\
            _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
            Na1 Na 0.0 0.0 0.0\nCl1 Cl 0.5 0.5 0.5\n";

        #[test]
        fn parses_site_count_and_lattice() {
            let r = parse_cif_periodic_structure(NACL_P1_CIF).unwrap();
            assert_eq!(r.structure.site_count(), 2);
            assert!((r.structure.lattice().lengths()[0] - 5.6402).abs() < 1e-4);
            assert_eq!(r.symmetry, CifSymmetryStatus::P1);
        }

        #[test]
        fn no_symmetry_tags_at_all_classifies_as_p1() {
            // NACL_P1_CIF has zero _symmetry_*/_space_group_* tags -- matches
            // the existing parser's "effectively P1" scope for such files.
            let r = parse_cif_periodic_structure(NACL_P1_CIF).unwrap();
            assert_eq!(r.symmetry, CifSymmetryStatus::P1);
        }

        /// Shared C2/c (space group No. 15) fixture: a realistic 8-operation
        /// symop list with one asymmetric-unit atom sitting on a special
        /// position. Used both to prove the default (expand=true) path
        /// genuinely expands+dedups, and (via `expand_explicit_symmetry:
        /// false`) that opting out still reports the pre-expansion
        /// asymmetric-unit-only behavior.
        const C2C_CIF: &str = "data_synthetic_c2c\n\
            _cell_length_a 10.0\n_cell_length_b 8.0\n_cell_length_c 12.0\n\
            _cell_angle_alpha 90\n_cell_angle_beta 105.0\n_cell_angle_gamma 90\n\
            _symmetry_space_group_name_H-M 'C 2/c'\n\
            loop_\n_space_group_symop_operation_xyz\n\
            'x, y, z'\n'-x, y, -z+1/2'\n'-x, -y, -z'\n'x, -y, z+1/2'\n\
            'x+1/2, y+1/2, z'\n'-x+1/2, y+1/2, -z+1/2'\n'-x+1/2, -y+1/2, -z'\n\
            'x+1/2, -y+1/2, z+1/2'\n\
            loop_\n_atom_site_label\n_atom_site_type_symbol\n\
            _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
            Ti1 Ti 0.0 0.25 0.25\n";

        #[test]
        fn single_identity_operation_cif_expands_trivially_by_default() {
            // A genuinely P1-equivalent CIF that explicitly lists its one
            // (identity) operation -- IUCr's own convention for how a P1
            // file states its symmetry explicitly. Under the default
            // (expand_explicit_symmetry=true), this is not P1-classified by
            // absence of an operation list; it is honestly reported as
            // ExpandedExplicitOperations, expanding trivially to itself
            // (never an error -- see CifSymmetryError::MissingIdentityOperation's
            // docs for why a *present* identity, even alone, is not an error).
            let cif = "data_x\n\
                _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 1'\n\
                loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0.0 0.0 0.0\n";
            let r = parse_cif_periodic_structure(cif).unwrap();
            assert_eq!(r.structure.site_count(), 1);
            assert_eq!(
                r.symmetry,
                CifSymmetryStatus::ExpandedExplicitOperations {
                    space_group_name: Some("P 1".to_string()),
                    operation_count: 1,
                    asymmetric_site_count: 1,
                    expanded_site_count: 1,
                }
            );
            assert!(r.symmetry.is_complete_cell());
        }

        #[test]
        fn single_identity_operation_cif_reports_unexpanded_when_expansion_is_opted_out() {
            let cif = "data_x\n\
                _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 1'\n\
                loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0.0 0.0 0.0\n";
            let r = parse_cif_periodic_structure_with_options(
                cif,
                CifPeriodicParseOptions {
                    expand_explicit_symmetry: false,
                },
            )
            .unwrap();
            assert_eq!(r.structure.site_count(), 1);
            assert_eq!(
                r.symmetry,
                CifSymmetryStatus::UnexpandedSymmetry {
                    space_group_name: Some("P 1".to_string()),
                    operation_count: 1,
                }
            );
            assert!(!r.symmetry.is_complete_cell());
        }

        #[test]
        fn declared_non_p1_space_group_expands_by_default_and_dedups_a_special_position() {
            let r = parse_cif_periodic_structure(C2C_CIF).unwrap();
            // Ti1 sits on a special position of C2/c (Wyckoff 4e): the
            // 8-operation general-position list collapses pairwise, so the
            // expanded cell has 4 distinct sites, not 8.
            assert_eq!(r.structure.site_count(), 4);
            match r.symmetry {
                CifSymmetryStatus::ExpandedExplicitOperations {
                    space_group_name,
                    operation_count,
                    asymmetric_site_count,
                    expanded_site_count,
                } => {
                    assert_eq!(space_group_name.as_deref(), Some("C 2/c"));
                    assert_eq!(operation_count, 8);
                    assert_eq!(asymmetric_site_count, 1);
                    assert_eq!(expanded_site_count, 4);
                }
                other => panic!("expected ExpandedExplicitOperations, got {other:?}"),
            }
        }

        #[test]
        fn declared_non_p1_space_group_is_not_silently_treated_as_p1_when_not_expanded() {
            let r = parse_cif_periodic_structure_with_options(
                C2C_CIF,
                CifPeriodicParseOptions {
                    expand_explicit_symmetry: false,
                },
            )
            .unwrap();
            assert_eq!(r.structure.site_count(), 1);
            match r.symmetry {
                CifSymmetryStatus::UnexpandedSymmetry {
                    space_group_name,
                    operation_count,
                } => {
                    assert_eq!(space_group_name.as_deref(), Some("C 2/c"));
                    assert_eq!(operation_count, 8);
                }
                other => panic!("expected UnexpandedSymmetry, got {other:?}"),
            }
        }

        #[test]
        fn it_number_alone_flags_non_p1_even_without_a_name_or_symop_loop() {
            let cif = "data_x\n\
                _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_Int_Tables_number 15\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                Ti1 Ti 0.0 0.25 0.25\n";
            let r = parse_cif_periodic_structure(cif).unwrap();
            match r.symmetry {
                CifSymmetryStatus::UnexpandedSymmetry {
                    space_group_name,
                    operation_count,
                } => {
                    assert_eq!(space_group_name, None);
                    assert_eq!(operation_count, 0);
                }
                other => panic!("expected UnexpandedSymmetry, got {other:?}"),
            }
        }

        #[test]
        fn disorder_rows_at_the_same_site_merge_into_multi_species_site() {
            // Realistic disorder loop shape: extra columns
            // (U_iso_or_equiv/adp_type/symmetry_multiplicity/calc_flag)
            // interleaved before occupancy, as commonly exported by COD/ICSD
            // -- exercises column-offset resolution, not just a 5-column toy.
            let cif = "data_synthetic_disordered_alloy\n\
                _cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n\
                _atom_site_label\n\
                _atom_site_type_symbol\n\
                _atom_site_fract_x\n\
                _atom_site_fract_y\n\
                _atom_site_fract_z\n\
                _atom_site_U_iso_or_equiv\n\
                _atom_site_adp_type\n\
                _atom_site_occupancy\n\
                _atom_site_symmetry_multiplicity\n\
                _atom_site_calc_flag\n\
                Fe1 Fe 0.0 0.0 0.0 0.012 Uiso 0.6 1 d\n\
                Ni1 Ni 0.0 0.0 0.0 0.012 Uiso 0.4 1 d\n\
                O1  O  0.5 0.5 0.5 0.020 Uiso 1.0 1 d\n";
            let r = parse_cif_periodic_structure(cif).unwrap();
            assert_eq!(
                r.structure.site_count(),
                2,
                "Fe1/Ni1 must merge into one site"
            );
            let disordered = &r.structure.sites()[0];
            assert_eq!(disordered.species.len(), 2);
            assert_eq!(disordered.label.as_deref(), Some("Fe1"));
            assert_eq!(disordered.species[0].element, Element::FE);
            assert_eq!(disordered.species[1].element, Element::NI);
            let occs: Vec<f64> = disordered
                .species
                .iter()
                .map(|s| s.occupancy.value())
                .collect();
            assert!((occs[0] - 0.6).abs() < 1e-9);
            assert!((occs[1] - 0.4).abs() < 1e-9);
            let ordered = &r.structure.sites()[1];
            assert_eq!(ordered.species.len(), 1);
            assert_eq!(ordered.species[0].element, Element::O);
            assert!((ordered.species[0].occupancy.value() - 1.0).abs() < 1e-9);
        }

        #[test]
        fn occupancy_sum_exceeded_is_rejected_by_crystals_own_validator() {
            let cif = "data_x\n\
                _cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n\
                Fe1 Fe 0.0 0.0 0.0 0.7\nNi1 Ni 0.0 0.0 0.0 0.5\n";
            let err = parse_cif_periodic_structure(cif).unwrap_err();
            // Raised by PeriodicSite::new itself (before PeriodicStructure
            // wraps per-site errors in InvalidSite{index,..}), proving this
            // adapter reuses chematic-crystal's own occupancy-sum validator
            // rather than re-checking the sum a second, different way.
            match err {
                CifPeriodicError::Crystal(CrystalError::OccupancySumExceeded { sum, .. }) => {
                    assert!((sum - 1.2).abs() < 1e-9);
                }
                other => panic!("expected Crystal(OccupancySumExceeded), got {other:?}"),
            }
        }

        #[test]
        fn occupancy_placeholder_dot_and_question_mark_default_to_full() {
            let cif = "data_x\n\
                _cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n\
                C1 C 0.0 0.0 0.0 .\nC2 C 0.5 0.5 0.5 ?\n";
            let r = parse_cif_periodic_structure(cif).unwrap();
            for site in r.structure.sites() {
                assert!((site.species[0].occupancy.value() - 1.0).abs() < 1e-9);
            }
        }

        #[test]
        fn malformed_occupancy_value_is_an_explicit_error_not_a_silent_default() {
            let cif = "data_x\n\
                _cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n\
                C1 C 0.0 0.0 0.0 0.5x\n";
            let err = parse_cif_periodic_structure(cif).unwrap_err();
            assert_eq!(err, CifPeriodicError::InvalidOccupancy("0.5x".to_string()));
        }

        #[test]
        fn missing_one_of_six_cell_tags_is_rejected_even_though_parse_cif_defaults_it() {
            // parse_cif would silently default the missing beta to 90.0;
            // the periodic adapter must not build a Lattice from a
            // partially-defaulted cell.
            let cif = "data_x\n\
                _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
                _cell_angle_alpha 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0.0 0.0 0.0\n";
            let err = parse_cif_periodic_structure(cif).unwrap_err();
            assert_eq!(err, CifPeriodicError::Cif(CifError::MissingCellParameters));
            // Sanity: parse_cif (the existing, unaffected API) does NOT
            // reject this file -- confirms the two functions genuinely
            // have different cell-completeness requirements.
            assert!(super::super::parse_cif(cif).is_ok());
        }

        #[test]
        fn cartesian_only_cif_converts_through_the_lattice() {
            let cif = "data_x\n\
                _cell_length_a 10.0\n_cell_length_b 10.0\n_cell_length_c 10.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_Cartn_x\n_atom_site_Cartn_y\n_atom_site_Cartn_z\n\
                C1 C 5.0 5.0 5.0\n";
            let r = parse_cif_periodic_structure(cif).unwrap();
            let f = r.structure.sites()[0].fractional.0;
            assert!((f[0] - 0.5).abs() < 1e-9);
            assert!((f[1] - 0.5).abs() < 1e-9);
            assert!((f[2] - 0.5).abs() < 1e-9);
        }

        #[test]
        fn write_then_parse_round_trips_as_p1_and_preserves_disorder() {
            let cif = "data_synthetic_disordered_alloy\n\
                _cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n\
                Fe1 Fe 0.0 0.0 0.0 0.6\nNi1 Ni 0.0 0.0 0.0 0.4\n\
                O1  O  0.5 0.5 0.5 1.0\n";
            let parsed = parse_cif_periodic_structure(cif).unwrap();
            let written = write_cif_periodic_structure(&parsed.structure);
            assert!(written.contains("_atom_site_occupancy"));
            let roundtripped = parse_cif_periodic_structure(&written).unwrap();
            assert_eq!(roundtripped.symmetry, CifSymmetryStatus::P1);
            assert_eq!(roundtripped.structure.site_count(), 2);
            let disordered = &roundtripped.structure.sites()[0];
            assert_eq!(disordered.species.len(), 2);
        }

        // -----------------------------------------------------------------
        // Explicit symmetry expansion -- end-to-end (full CIF text) fixtures.
        // Parser/expansion internals are exercised in isolation by
        // `crate::cif_symmetry::tests`; these confirm the whole pipeline
        // (tokenize -> resolve operations -> expand -> PeriodicStructure)
        // agrees, using real/citable space-group data where practical.
        // -----------------------------------------------------------------

        fn cif_with_ops(ops_loop: &str, atom_line: &str) -> String {
            format!(
                "data_x\n\
                 _cell_length_a 6.0\n_cell_length_b 8.0\n_cell_length_c 10.0\n\
                 _cell_angle_alpha 90\n_cell_angle_beta 100.0\n_cell_angle_gamma 90\n\
                 {ops_loop}\n\
                 loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                 _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                 {atom_line}\n"
            )
        }

        #[test]
        fn p1_single_identity_operation_expands_to_itself() {
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'",
                "C1 C 0.1 0.2 0.3",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            assert_eq!(r.structure.site_count(), 1);
            assert!(matches!(
                r.symmetry,
                CifSymmetryStatus::ExpandedExplicitOperations {
                    operation_count: 1,
                    asymmetric_site_count: 1,
                    expanded_site_count: 1,
                    ..
                }
            ));
        }

        #[test]
        fn p_minus_1_general_position_expands_to_two_sites_identity_first() {
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'",
                "C1 C 0.1 0.2 0.3",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            assert_eq!(r.structure.site_count(), 2);
            let sites = r.structure.sites();
            assert_eq!(sites[0].label.as_deref(), Some("C1"));
            let f0 = sites[0].fractional.0;
            assert!(
                (f0[0] - 0.1).abs() < 1e-9
                    && (f0[1] - 0.2).abs() < 1e-9
                    && (f0[2] - 0.3).abs() < 1e-9
            );
            assert_eq!(sites[1].label.as_deref(), Some("C1@sym2"));
            let f1 = sites[1].fractional.0;
            assert!(
                (f1[0] - 0.9).abs() < 1e-9
                    && (f1[1] - 0.8).abs() < 1e-9
                    && (f1[2] - 0.7).abs() < 1e-9
            );
        }

        #[test]
        fn p_minus_1_special_position_at_origin_is_not_duplicated() {
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'",
                "C1 C 0.0 0.0 0.0",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            assert_eq!(r.structure.site_count(), 1);
            assert_eq!(r.structure.sites()[0].label.as_deref(), Some("C1"));
        }

        #[test]
        fn cell_wrapping_translation_operation_wraps_into_unit_cell() {
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'x+1/2,y,z'",
                "C1 C 0.8 0.1 0.1",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            assert_eq!(r.structure.site_count(), 2);
            let f1 = r.structure.sites()[1].fractional.0;
            assert!(
                (f1[0] - 0.3).abs() < 1e-9,
                "expected wrap to 0.3, got {}",
                f1[0]
            );
        }

        #[test]
        fn non_trivial_rational_translation_and_variable_combination() {
            // "-y+x,-y,1/3+z": x' = x-y, y' = -y, z' = z+1/3.
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-y+x,-y,1/3+z'",
                "C1 C 0.5 0.6 0.7",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            assert_eq!(r.structure.site_count(), 2);
            let f1 = r.structure.sites()[1].fractional.0;
            assert!((f1[0] - (0.5 - 0.6f64).rem_euclid(1.0)).abs() < 1e-9);
            assert!((f1[1] - (-0.6f64).rem_euclid(1.0)).abs() < 1e-9);
            assert!((f1[2] - (0.7f64 + 1.0 / 3.0).rem_euclid(1.0)).abs() < 1e-9);
        }

        #[test]
        fn modern_dotted_operation_tag_is_supported() {
            let cif = cif_with_ops(
                "loop_\n_space_group_symop.operation_xyz\n'x,y,z'\n'-x,-y,-z'",
                "C1 C 0.1 0.2 0.3",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            assert_eq!(r.structure.site_count(), 2);
        }

        #[test]
        fn modern_dotted_equiv_pos_tag_is_supported() {
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv.pos_as_xyz\n'x,y,z'\n'-x,-y,-z'",
                "C1 C 0.1 0.2 0.3",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            assert_eq!(r.structure.site_count(), 2);
        }

        #[test]
        fn legacy_underscore_equiv_pos_tag_is_supported() {
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'",
                "C1 C 0.1 0.2 0.3",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            assert_eq!(r.structure.site_count(), 2);
        }

        #[test]
        fn trigonal_p3_three_operation_general_position_expands_to_three_sites() {
            // P3 (No. 143), standard general positions -- a real 3-operation
            // (3-fold rotation) case, distinct from the 2- and 4-operation
            // fixtures above.
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-y,x-y,z'\n'y-x,-x,z'",
                "C1 C 0.2 0.35 0.4",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            assert_eq!(r.structure.site_count(), 3);
            let labels: Vec<Option<String>> = r
                .structure
                .sites()
                .iter()
                .map(|s| s.label.clone())
                .collect();
            assert_eq!(
                labels,
                vec![
                    Some("C1".to_string()),
                    Some("C1@sym2".to_string()),
                    Some("C1@sym3".to_string()),
                ]
            );
        }

        #[test]
        fn fe_ni_disorder_survives_symmetry_expansion_intact() {
            // cif_with_ops's atom loop header has no occupancy column, so
            // this CIF is built directly instead of via that helper.
            let cif = "data_x\n\
                _cell_length_a 6.0\n_cell_length_b 8.0\n_cell_length_c 10.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 100.0\n_cell_angle_gamma 90\n\
                loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n\
                Fe1 Fe 0.1 0.2 0.3 0.6\nNi1 Ni 0.1 0.2 0.3 0.4\n";
            let r = parse_cif_periodic_structure(cif).unwrap();
            assert_eq!(
                r.structure.site_count(),
                2,
                "one disordered site, two images"
            );
            for site in r.structure.sites() {
                assert_eq!(site.species.len(), 2);
                let fe = site
                    .species
                    .iter()
                    .find(|s| s.element.symbol() == "Fe")
                    .unwrap();
                let ni = site
                    .species
                    .iter()
                    .find(|s| s.element.symbol() == "Ni")
                    .unwrap();
                assert!((fe.occupancy.value() - 0.6).abs() < 1e-9);
                assert!((ni.occupancy.value() - 0.4).abs() < 1e-9);
            }
        }

        #[test]
        fn duplicate_operation_is_a_typed_error() {
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'\n'-x,-y,-z'",
                "C1 C 0.1 0.2 0.3",
            );
            let err = parse_cif_periodic_structure(&cif).unwrap_err();
            assert!(matches!(
                err,
                CifPeriodicError::Symmetry(
                    cif_symmetry::CifSymmetryError::DuplicateSymmetryOperation { .. }
                )
            ));
        }

        #[test]
        fn missing_identity_is_a_typed_error() {
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv_pos_as_xyz\n'-x,-y,-z'",
                "C1 C 0.1 0.2 0.3",
            );
            let err = parse_cif_periodic_structure(&cif).unwrap_err();
            assert_eq!(
                err,
                CifPeriodicError::Symmetry(
                    cif_symmetry::CifSymmetryError::MissingIdentityOperation
                )
            );
        }

        #[test]
        fn declared_space_group_with_no_operation_list_at_all_stays_unexpanded() {
            // Only a name/number, no _space_group_symop*/_symmetry_equiv*
            // loop -- nothing to expand from, so this must stay
            // UnexpandedSymmetry even under the default expand=true.
            let cif = "data_x\n\
                _cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 5.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 21/c'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0.1 0.2 0.3\n";
            let r = parse_cif_periodic_structure(cif).unwrap();
            assert_eq!(r.structure.site_count(), 1);
            match r.symmetry {
                CifSymmetryStatus::UnexpandedSymmetry {
                    space_group_name,
                    operation_count,
                } => {
                    assert_eq!(space_group_name.as_deref(), Some("P 21/c"));
                    assert_eq!(operation_count, 0);
                }
                other => panic!("expected UnexpandedSymmetry, got {other:?}"),
            }
        }

        /// IUCr International Tables for Crystallography Vol. A, space
        /// group P2_1/c (No. 14), standard-setting general positions --
        /// public crystallographic convention, not copied from any tool's
        /// implementation. Self-verifying: composes every pair of the 4
        /// operations using test-local arithmetic (independent of the
        /// parser under test) and confirms the result matches one of the 4
        /// operations modulo integer translation, i.e. these 4 operators
        /// really do form a closed group before trusting them as a golden
        /// fixture.
        #[test]
        fn iucr_p21c_four_operation_golden_fixture() {
            type Mat = [[f64; 3]; 3];
            let identity: Mat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
            let screw: Mat = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]];
            let inversion: Mat = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];
            let glide: Mat = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]];
            let ops: [(Mat, [f64; 3]); 4] = [
                (identity, [0.0, 0.0, 0.0]),
                (screw, [0.0, 0.5, 0.5]),
                (inversion, [0.0, 0.0, 0.0]),
                (glide, [0.0, 0.5, 0.5]),
            ];

            fn apply(rot: Mat, trans: [f64; 3], p: [f64; 3]) -> [f64; 3] {
                let mut out = [0.0; 3];
                for i in 0..3 {
                    out[i] = trans[i] + (0..3).map(|k| rot[i][k] * p[k]).sum::<f64>();
                }
                out
            }
            fn compose(a: (Mat, [f64; 3]), b: (Mat, [f64; 3])) -> (Mat, [f64; 3]) {
                let (rot_a, t_a) = a;
                let (rot_b, t_b) = b;
                let mut rot = [[0.0; 3]; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        rot[i][j] = (0..3).map(|k| rot_a[i][k] * rot_b[k][j]).sum();
                    }
                }
                let t = apply(rot_a, t_a, t_b);
                (rot, t)
            }
            fn mats_eq(a: Mat, b: Mat) -> bool {
                (0..3).all(|i| (0..3).all(|j| (a[i][j] - b[i][j]).abs() < 1e-9))
            }
            fn trans_eq_mod1(a: [f64; 3], b: [f64; 3]) -> bool {
                (0..3).all(|i| {
                    let d = (a[i] - b[i]).rem_euclid(1.0);
                    d < 1e-9 || (1.0 - d) < 1e-9
                })
            }

            // Group closure: every pairwise composition matches one of the
            // 4 declared operations, modulo integer translation.
            for &a in &ops {
                for &b in &ops {
                    let (rot_c, t_c) = compose(a, b);
                    let matches_some_op = ops
                        .iter()
                        .any(|&(rot, t)| mats_eq(rot_c, rot) && trans_eq_mod1(t_c, t));
                    assert!(
                        matches_some_op,
                        "composition not closed within the declared set"
                    );
                }
            }

            // Now confirm the actual parser+expander agrees with this
            // independently-computed reference on a general position.
            let cif = cif_with_ops(
                "loop_\n_space_group_symop_id\n_space_group_symop_operation_xyz\n\
                 1 'x,y,z'\n2 '-x,y+1/2,-z+1/2'\n3 '-x,-y,-z'\n4 'x,-y+1/2,z+1/2'",
                "C1 C 0.2 0.3 0.4",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            assert_eq!(r.structure.site_count(), 4);
            assert_eq!(
                r.symmetry,
                CifSymmetryStatus::ExpandedExplicitOperations {
                    space_group_name: None,
                    operation_count: 4,
                    asymmetric_site_count: 1,
                    expanded_site_count: 4,
                }
            );
            let p0 = [0.2, 0.3, 0.4];
            let expected_labels = ["C1", "C1@sym2", "C1@sym3", "C1@sym4"];
            for (i, site) in r.structure.sites().iter().enumerate() {
                let (rot, t) = ops[i];
                let expected = apply(rot, t, p0).map(|c| c.rem_euclid(1.0));
                let actual = site.fractional.0;
                for k in 0..3 {
                    assert!(
                        (actual[k] - expected[k]).abs() < 1e-9,
                        "op {i}: expected {expected:?}, got {actual:?}"
                    );
                }
                assert_eq!(site.label.as_deref(), Some(expected_labels[i]));
            }
        }

        // -----------------------------------------------------------------
        // to_cif_checked
        // -----------------------------------------------------------------

        #[test]
        fn to_cif_checked_succeeds_on_expanded_structure() {
            let cif = cif_with_ops(
                "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'",
                "C1 C 0.1 0.2 0.3",
            );
            let r = parse_cif_periodic_structure(&cif).unwrap();
            let written = r.to_cif_checked().unwrap();
            assert!(written.contains("_symmetry_space_group_name_H-M"));
            let reparsed = parse_cif_periodic_structure(&written).unwrap();
            assert_eq!(reparsed.structure.site_count(), 2);
            assert_eq!(reparsed.symmetry, CifSymmetryStatus::P1);
        }

        #[test]
        fn to_cif_checked_refuses_unexpanded_structure() {
            let r = parse_cif_periodic_structure_with_options(
                C2C_CIF,
                CifPeriodicParseOptions {
                    expand_explicit_symmetry: false,
                },
            )
            .unwrap();
            let err = r.to_cif_checked().unwrap_err();
            assert!(matches!(
                err,
                CifPeriodicError::UnexpandedSymmetryWrite {
                    operation_count: 8,
                    ..
                }
            ));
        }

        #[test]
        fn to_cif_checked_succeeds_on_plain_p1() {
            let r = parse_cif_periodic_structure(NACL_P1_CIF).unwrap();
            assert!(r.to_cif_checked().is_ok());
        }
    }
}

#[cfg(feature = "crystal")]
pub use crate::cif_symmetry::CifSymmetryError;
#[cfg(feature = "crystal")]
pub use crystal_adapter::{
    CifPeriodicError, CifPeriodicParseOptions, CifPeriodicResult, CifSymmetryStatus,
    parse_cif_periodic_structure, parse_cif_periodic_structure_with_options,
    write_cif_periodic_structure,
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const NACL_CIF: &str = r#"data_NaCl
_cell_length_a   5.6402
_cell_length_b   5.6402
_cell_length_c   5.6402
_cell_angle_alpha  90.000
_cell_angle_beta   90.000
_cell_angle_gamma  90.000

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na1  Na  0.00000  0.00000  0.00000
Cl1  Cl  0.50000  0.50000  0.50000
"#;

    #[test]
    fn parse_nacl_atom_count() {
        let r = parse_cif(NACL_CIF).unwrap();
        assert_eq!(r.mol.atom_count(), 2);
        assert_eq!(r.coords.len(), 2);
    }

    #[test]
    fn parse_nacl_cell_params() {
        let r = parse_cif(NACL_CIF).unwrap();
        let cell = r.cell.as_ref().unwrap();
        assert!((cell.a - 5.6402).abs() < 1e-4);
        assert!((cell.alpha - 90.0).abs() < 1e-4);
    }

    #[test]
    fn parse_nacl_na_at_origin() {
        let r = parse_cif(NACL_CIF).unwrap();
        let (x, y, z) = r.coords[0];
        assert!(
            x.abs() < 1e-6 && y.abs() < 1e-6 && z.abs() < 1e-6,
            "Na at origin: got ({x}, {y}, {z})"
        );
    }

    #[test]
    fn parse_nacl_cl_at_body_center() {
        let r = parse_cif(NACL_CIF).unwrap();
        let (x, y, z) = r.coords[1];
        let half = 5.6402 / 2.0;
        assert!((x - half).abs() < 1e-3, "Cl x: got {x}, expected ~{half}");
        assert!((y - half).abs() < 1e-3);
        assert!((z - half).abs() < 1e-3);
    }

    #[test]
    fn unit_cell_volume_cubic() {
        let cell = UnitCell {
            a: 5.0,
            b: 5.0,
            c: 5.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        };
        assert!((cell.volume() - 125.0).abs() < 1e-6);
    }

    #[test]
    fn write_cif_roundtrip() {
        let r = parse_cif(NACL_CIF).unwrap();
        let out = write_cif(&r.mol, &r.coords, r.cell.as_ref());
        assert!(out.contains("data_chematic"));
        let r2 = parse_cif(&out).unwrap();
        assert_eq!(r2.mol.atom_count(), 2);
    }

    #[test]
    fn parse_esd_stripped() {
        let cif = "data_x\n_cell_length_a 5.64(2)\n_cell_length_b 5.64\n_cell_length_c 5.64\n\
                   _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                   loop_\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                   Na 0.0 0.0 0.0\n";
        let r = parse_cif(cif).unwrap();
        assert!((r.cell.unwrap().a - 5.64).abs() < 1e-3);
    }

    #[test]
    fn frac_to_cart_cubic_roundtrip() {
        let cell = UnitCell {
            a: 5.0,
            b: 5.0,
            c: 5.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        };
        let (x, y, z) = cell.frac_to_cart(0.3, 0.4, 0.5);
        let (fx, fy, fz) = cart_to_frac(&cell, x, y, z);
        assert!((fx - 0.3).abs() < 1e-10);
        assert!((fy - 0.4).abs() < 1e-10);
        assert!((fz - 0.5).abs() < 1e-10);
    }
}
