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

use chematic_core::{Atom, Element, Molecule, MoleculeBuilder};

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
    /// Validate lengths, angles and numerical invertibility of the cell.
    ///
    /// Lengths must be finite and positive; angles must lie strictly between
    /// 0 and 180 degrees. The dimensionless Gram determinant must exceed
    /// `1e-12`, and the volume must be finite and positive. Near-singular cells
    /// below this numerical threshold are rejected as well as singular cells.
    pub fn validate(&self) -> Result<(), CifError> {
        for (name, value) in [("a", self.a), ("b", self.b), ("c", self.c)] {
            if !value.is_finite() || value <= 0.0 {
                return Err(CifError::InvalidCellParameters(format!(
                    "length {name} must be finite and positive, got {value}"
                )));
            }
        }
        for (name, value) in [
            ("alpha", self.alpha),
            ("beta", self.beta),
            ("gamma", self.gamma),
        ] {
            if !value.is_finite() || value <= 0.0 || value >= 180.0 {
                return Err(CifError::InvalidCellParameters(format!(
                    "angle {name} must be finite and between 0 and 180 degrees, got {value}"
                )));
            }
        }
        let ca = self.alpha.to_radians().cos();
        let cb = self.beta.to_radians().cos();
        let cg = self.gamma.to_radians().cos();
        let determinant = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg;
        if determinant <= 1e-12 {
            return Err(CifError::InvalidCellParameters(format!(
                "cell is singular or numerically near-singular (Gram determinant {determinant})"
            )));
        }
        let volume = self.volume();
        if !volume.is_finite() || volume <= 0.0 {
            return Err(CifError::InvalidCellParameters(format!(
                "unit cell volume must be finite and positive, got {volume}"
            )));
        }
        Ok(())
    }

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

/// Resource limits applied while parsing a plain CIF document.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CifParseLimits {
    /// Maximum UTF-8 input size, in bytes.
    pub max_input_bytes: usize,
    /// Maximum physical line size, in bytes.
    pub max_line_bytes: usize,
    /// Maximum number of tokens produced by the CIF tokenizer.
    pub max_tokens: usize,
    /// Maximum number of atom-site rows returned.
    pub max_atoms: usize,
}

impl Default for CifParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_line_bytes: 1024 * 1024,
            max_tokens: 10_000_000,
            max_atoms: 1_000_000,
        }
    }
}

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Errors that can occur when parsing or writing CIF files.
#[derive(Debug, Clone, PartialEq)]
pub enum CifError {
    /// No `_atom_site_*` loop with coordinates was found.
    NoAtomSiteLoop,
    /// The atom site loop lacked at least one coordinate column.
    MissingCoordinateColumns,
    /// An element symbol or label could not be resolved.
    UnknownElement(String),
    /// A coordinate value is invalid, non-finite, or overflows during conversion.
    InvalidCoordinate(String),
    /// The unit cell has invalid lengths, angles, volume or numerical invertibility.
    InvalidCellParameters(String),
    /// The atom site loop uses fractional coordinates but the CIF contains no
    /// `_cell_length_*` / `_cell_angle_*` parameters to convert them.
    MissingCellParameters,
    /// The document exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl core::fmt::Display for CifError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoAtomSiteLoop => write!(f, "no _atom_site_* loop found in CIF"),
            Self::MissingCoordinateColumns => {
                write!(
                    f,
                    "atom_site loop requires a complete fract_x/y/z or Cartn_x/y/z set"
                )
            }
            Self::UnknownElement(s) => write!(f, "unknown element '{s}' in CIF"),
            Self::InvalidCoordinate(s) => write!(f, "invalid coordinate '{s}' in CIF"),
            Self::InvalidCellParameters(s) => write!(f, "invalid cell parameters: {s}"),
            Self::MissingCellParameters => write!(
                f,
                "fractional coordinates present but no _cell_length_*/_cell_angle_* parameters found"
            ),
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "CIF: {resource} has size {actual}, exceeding limit {limit}"
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

#[derive(Clone, Copy)]
struct AtomSiteColumns {
    width: usize,
    element: Option<usize>,
    label: Option<usize>,
    x: usize,
    y: usize,
    z: usize,
    #[cfg(feature = "crystal")]
    occupancy: Option<usize>,
    cartesian: bool,
}

impl AtomSiteColumns {
    fn from_headers(headers: &[String]) -> Result<Self, CifError> {
        let col = |name: &str| headers.iter().position(|header| header == name);
        // All three components must describe the same coordinate system.
        // Prefer a complete fractional set; a stray fractional column must
        // not override a complete Cartesian set or be mixed into it.
        let complete = |prefix: &str| {
            Some((
                col(&format!("_atom_site_{prefix}_x"))?,
                col(&format!("_atom_site_{prefix}_y"))?,
                col(&format!("_atom_site_{prefix}_z"))?,
            ))
        };
        let (x, y, z, cartesian) = if let Some((x, y, z)) = complete("fract") {
            (x, y, z, false)
        } else if let Some((x, y, z)) = complete("cartn") {
            (x, y, z, true)
        } else {
            return Err(CifError::MissingCoordinateColumns);
        };
        Ok(Self {
            width: headers.len(),
            element: col("_atom_site_type_symbol"),
            label: col("_atom_site_label"),
            x,
            y,
            z,
            #[cfg(feature = "crystal")]
            occupancy: col("_atom_site_occupancy"),
            cartesian,
        })
    }

    fn is_end(self, row: &[String]) -> bool {
        row[0] == "loop_" || row[0].starts_with("data_")
    }

    fn element(self, row: &[String]) -> Result<Element, CifError> {
        let raw = self
            .element
            .and_then(|column| row.get(column))
            .or_else(|| self.label.and_then(|column| row.get(column)))
            .map(String::as_str)
            .unwrap_or("X");
        resolve_element(raw)
    }
}

fn parse_cif_number(token: &str) -> Result<f64, CifError> {
    strip_esd(token)
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| CifError::InvalidCoordinate(token.to_owned()))
}

fn validate_fractional_cell(cell: &UnitCell, has_cell: bool) -> Result<(), CifError> {
    if !has_cell {
        return Err(CifError::MissingCellParameters);
    }
    cell.validate()
}

fn validate_plain_cif_input(input: &str, limits: &CifParseLimits) -> Result<Vec<String>, CifError> {
    if input.len() > limits.max_input_bytes {
        return Err(CifError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    if let Some(actual) = input.lines().map(str::len).max()
        && actual > limits.max_line_bytes
    {
        return Err(CifError::ResourceLimit {
            resource: "line bytes",
            actual,
            limit: limits.max_line_bytes,
        });
    }
    let clean = input
        .lines()
        .map(strip_cif_comment)
        .collect::<Vec<_>>()
        .join("\n");
    let tokens = tokenize_cif(&clean);
    if tokens.len() > limits.max_tokens {
        return Err(CifError::ResourceLimit {
            resource: "tokens",
            actual: tokens.len(),
            limit: limits.max_tokens,
        });
    }
    Ok(tokens)
}

struct PlainCifAtoms {
    mol: Molecule,
    coords: Vec<(f64, f64, f64)>,
}

fn parse_plain_cif_rows(
    tokens: &[String],
    start: usize,
    columns: AtomSiteColumns,
    cell: &UnitCell,
    limits: &CifParseLimits,
) -> Result<PlainCifAtoms, CifError> {
    let mut builder = MoleculeBuilder::new();
    let mut coords = Vec::new();
    for row in tokens[start..].chunks(columns.width) {
        if row.len() < columns.width || columns.is_end(row) {
            break;
        }
        if coords.len() >= limits.max_atoms {
            return Err(CifError::ResourceLimit {
                resource: "atom rows",
                actual: coords.len() + 1,
                limit: limits.max_atoms,
            });
        }
        let values = (
            parse_cif_number(&row[columns.x])?,
            parse_cif_number(&row[columns.y])?,
            parse_cif_number(&row[columns.z])?,
        );
        let coord = if columns.cartesian {
            values
        } else {
            cell.frac_to_cart(values.0, values.1, values.2)
        };
        if !coord.0.is_finite() || !coord.1.is_finite() || !coord.2.is_finite() {
            return Err(CifError::InvalidCoordinate(
                "fractional-to-Cartesian conversion produced a non-finite value".into(),
            ));
        }
        builder.add_atom(Atom::new(columns.element(row)?));
        coords.push(coord);
    }
    if coords.is_empty() {
        return Err(CifError::NoAtomSiteLoop);
    }
    Ok(PlainCifAtoms {
        mol: builder.build(),
        coords,
    })
}

pub fn parse_cif(input: &str) -> Result<CifResult, CifError> {
    parse_cif_with_limits(input, &CifParseLimits::default())
}

/// Parse a CIF file with explicit resource limits.
pub fn parse_cif_with_limits(input: &str, limits: &CifParseLimits) -> Result<CifResult, CifError> {
    let tokens = validate_plain_cif_input(input, limits)?;
    let (cell, has_cell, _cell_fields_seen) = scan_cell(&tokens);
    let (headers, data_start) = find_atom_site_loop(&tokens)?;
    let columns = AtomSiteColumns::from_headers(&headers)?;
    if !columns.cartesian {
        validate_fractional_cell(&cell, has_cell)?;
    }
    let parsed = parse_plain_cif_rows(&tokens, data_start, columns, &cell, limits)?;
    Ok(CifResult {
        mol: parsed.mol,
        coords: parsed.coords,
        cell: if has_cell { Some(cell) } else { None },
    })
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

mod bonds;
pub use bonds::{CifBondWriteError, write_cif_with_bonds_checked};
mod writer;
pub use writer::{write_cif, write_cif_checked};

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
        AtomSiteColumns, CellFieldsSeen, CifError, UnitCell, find_atom_site_loop, parse_cif_number,
        scan_cell, strip_cif_comment, strip_esd, tokenize_cif,
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

    struct GroupedAtomRows {
        fractional: FractionalCoord,
        label: Option<String>,
        species: Vec<SiteSpecies>,
    }

    fn parse_periodic_rows(
        tokens: &[String],
        start: usize,
        columns: AtomSiteColumns,
        lattice: &Lattice,
    ) -> Result<Vec<RawAtomRow>, CifPeriodicError> {
        let mut rows = Vec::new();
        for row in tokens[start..].chunks(columns.width) {
            if row.len() < columns.width || columns.is_end(row) {
                break;
            }
            let values = [
                parse_cif_number(&row[columns.x])?,
                parse_cif_number(&row[columns.y])?,
                parse_cif_number(&row[columns.z])?,
            ];
            let fractional = if columns.cartesian {
                lattice.cart_to_frac(CartesianCoord::new(values))
            } else {
                FractionalCoord::new(values)
            };
            let occupancy = columns
                .occupancy
                .and_then(|column| row.get(column))
                .map(|token| parse_occupancy_token(token))
                .transpose()?
                .unwrap_or(1.0);
            rows.push(RawAtomRow {
                label: columns.label.and_then(|column| row.get(column)).cloned(),
                element: columns.element(row)?,
                fractional,
                occupancy,
            });
        }
        if rows.is_empty() {
            return Err(CifPeriodicError::Cif(CifError::NoAtomSiteLoop));
        }
        Ok(rows)
    }

    fn group_periodic_rows(
        rows: Vec<RawAtomRow>,
        lattice: &Lattice,
    ) -> Result<Vec<PeriodicSite>, CifPeriodicError> {
        let mut grouped: Vec<GroupedAtomRows> = Vec::new();
        for raw in rows {
            let species = SiteSpecies {
                element: raw.element,
                occupancy: Occupancy::new(raw.occupancy)?,
            };
            if let Some(group) = grouped
                .iter_mut()
                .find(|group| sites_within_tolerance(lattice, group.fractional, raw.fractional))
            {
                group.species.push(species);
            } else {
                grouped.push(GroupedAtomRows {
                    fractional: raw.fractional,
                    label: raw.label,
                    species: vec![species],
                });
            }
        }
        grouped
            .into_iter()
            .map(|group| PeriodicSite::new(group.species, group.fractional, group.label))
            .collect::<Result<Vec<_>, _>>()
            .map_err(Into::into)
    }

    fn finish_periodic_result(
        tokens: &[String],
        lattice: Lattice,
        asymmetric_sites: Vec<PeriodicSite>,
        space_group_name: Option<String>,
        it_number: Option<i64>,
        expand: bool,
    ) -> Result<CifPeriodicResult, CifPeriodicError> {
        if expand && let Some(operations) = cif_symmetry::resolve_symmetry_operations(tokens)? {
            let expanded = cif_symmetry::expand_sites(&asymmetric_sites, &operations, &lattice)?;
            let structure = PeriodicStructure::new(lattice, expanded)?;
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
        let operation_count = if expand {
            0
        } else {
            cif_symmetry::scan_operation_sources(tokens)
                .first()
                .map(Vec::len)
                .unwrap_or(0)
        };
        let structure = PeriodicStructure::new(lattice, asymmetric_sites)?;
        let symmetry = classify_p1_or_unexpanded(space_group_name, it_number, operation_count);
        Ok(CifPeriodicResult {
            structure,
            symmetry,
        })
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

        let (headers, data_start) = find_atom_site_loop(&tokens)?;
        let columns = AtomSiteColumns::from_headers(&headers)?;
        let rows = parse_periodic_rows(&tokens, data_start, columns, &lattice)?;
        let asymmetric_sites = group_periodic_rows(rows, &lattice)?;
        finish_periodic_result(
            &tokens,
            lattice,
            asymmetric_sites,
            space_group_name,
            it_number,
            options.expand_explicit_symmetry,
        )
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
}

#[cfg(feature = "crystal")]
pub use crate::cif_symmetry::CifSymmetryError;
#[cfg(feature = "crystal")]
pub use crystal_adapter::{
    CifPeriodicError, CifPeriodicParseOptions, CifPeriodicResult, CifSymmetryStatus,
    parse_cif_periodic_structure, parse_cif_periodic_structure_with_options,
    write_cif_periodic_structure,
};

#[cfg(all(test, feature = "crystal"))]
#[path = "cif/tests.rs"]
mod crystal_tests;

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
        let (fx, fy, fz) = writer::cart_to_frac(&cell, x, y, z);
        assert!((fx - 0.3).abs() < 1e-10);
        assert!((fy - 0.4).abs() < 1e-10);
        assert!((fz - 0.5).abs() < 1e-10);
    }

    #[test]
    fn bounded_parser_rejects_oversized_input() {
        let limits = CifParseLimits {
            max_input_bytes: 8,
            ..Default::default()
        };
        assert!(matches!(
            parse_cif_with_limits(NACL_CIF, &limits),
            Err(CifError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));
    }

    #[test]
    fn bounded_parser_rejects_long_lines_before_tokenizing() {
        let cif = format!("data_x\n#{}\n{}", "x".repeat(32), NACL_CIF);
        let limits = CifParseLimits {
            max_line_bytes: 16,
            ..Default::default()
        };
        assert!(matches!(
            parse_cif_with_limits(&cif, &limits),
            Err(CifError::ResourceLimit {
                resource: "line bytes",
                ..
            })
        ));
    }

    #[test]
    fn bounded_parser_rejects_token_and_atom_budgets() {
        let token_limits = CifParseLimits {
            max_tokens: 2,
            ..Default::default()
        };
        assert!(matches!(
            parse_cif_with_limits(NACL_CIF, &token_limits),
            Err(CifError::ResourceLimit {
                resource: "tokens",
                ..
            })
        ));

        let atom_limits = CifParseLimits {
            max_atoms: 1,
            ..Default::default()
        };
        assert!(matches!(
            parse_cif_with_limits(NACL_CIF, &atom_limits),
            Err(CifError::ResourceLimit {
                resource: "atom rows",
                ..
            })
        ));
    }
}
