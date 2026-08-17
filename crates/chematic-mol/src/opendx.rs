//! OpenDX scalar-field format (`.dx`) reader and writer -- specifically the
//! regular-grid scalar-field subset that APBS/electrostatics tooling
//! actually produces, not the full general OpenDX/IBM Data Explorer format
//! family.
//!
//! ## Sources and confidence
//!
//! Implemented independently from APBS's own public format documentation
//! (<https://apbs.readthedocs.io/en/latest/formats/opendx.html>). No source
//! code, comments, or tables were copied from APBS, VMD, PyMOL, Chimera, or
//! any other tool; where real-world example files are cited below as
//! corroboration, they were read only as a behavioral oracle, never as a
//! source to copy from.
//!
//! **High confidence** -- quoted/paraphrased directly from the APBS format
//! page:
//! - `object 1 class gridpositions counts nx ny nz`
//! - `origin xmin ymin zmin`
//! - 3 `delta` lines: `delta hx 0.0 0.0` / `delta 0.0 hy 0.0` /
//!   `delta 0.0 0.0 hz` in APBS's own typical (axis-aligned) output -- but
//!   each line is syntactically just `delta` followed by an arbitrary
//!   3-vector, so this reader accepts a fully general (non-axis-aligned)
//!   vector on each `delta` line, not just the diagonal case APBS itself
//!   happens to emit.
//! - `object 2 class gridconnections counts nx ny nz` (redundant with
//!   object 1's counts in every real APBS file; this reader verifies the
//!   two agree rather than silently trusting one or the other).
//! - `object 3 class array type double rank 0 items n data follows`,
//!   followed by `n` data values.
//! - Data ordering: "the data values, ordered with the z-index increasing
//!   most quickly, followed by the y-index, and then the x-index" -- i.e.
//!   the *same* first-axis-outermost/third-axis-innermost order Gaussian
//!   Cube uses (see `crate::cube` and `crate::volumetric` module docs),
//!   which is why this crate's shared [`VolumetricGrid`] value ordering
//!   needs no transpose for either format.
//! - Trailing `attribute "dep" string "positions"` / `object "regular
//!   positions regular connections" class field` / `component ...` lines:
//!   the APBS doc page shows these following the data block. This module
//!   does not parse or preserve them (see "Trailing metadata" below).
//!
//! **Medium confidence** -- not shown on the APBS format page itself, but
//! corroborated by real APBS-generated `.dx` files quoted in independent
//! third-party sources (e.g. mailing-list/support-forum posts showing
//! actual APBS output beginning with lines like `# Data from APBS 0.3.2` /
//! `# POTENTIAL (kT/e)`): real files commonly have one or more leading
//! `#`-prefixed comment lines before `object 1`. This reader skips any
//! number of leading `#` lines defensively; the APBS spec page's own
//! grammar excerpt does not show or require them.
//!
//! ## Scope: APBS scalar-field subset only
//!
//! OpenDX/IBM Data Explorer is, in general, a much broader format (it can
//! describe arbitrary irregular meshes, vector fields, multiple composite
//! objects, etc.). This module implements *only* the narrow regular-grid,
//! rank-0 (scalar), `type double` subset described above, because that is
//! what APBS and the electrostatics tools that consume its output actually
//! produce and read. A `rank`/`type`/`class` combination other than
//! exactly `array type double rank 0` is rejected with a typed
//! [`OpenDxError::UnsupportedArrayDeclaration`], not silently
//! misinterpreted.
//!
//! ## No atom section, no unit tag
//!
//! OpenDX carries no atom list at all (the source molecule/PQR is a
//! separate file); [`VolumetricGrid::atoms`] is always empty here. The
//! format also carries no explicit unit tag for `origin`/`delta` -- by
//! convention (not something read from the file), APBS's own real-world
//! usage is Ångström, so every grid parsed by this module is tagged
//! [`crate::volumetric::GridUnits::Angstrom`]; this is a stated assumption
//! about the format's dominant real-world producer, not information
//! recovered from the file itself.
//!
//! ## Trailing metadata: out of scope
//!
//! The `attribute`/`object "field is ..."`/`component` footer lines after
//! the data block are neither parsed nor preserved -- [`write_opendx`]
//! always regenerates the standard boilerplate footer itself rather than
//! reconstructing whatever a source file had. This is a deliberate scope
//! decision (round-trip fidelity is claimed for `origin`/`axes`/`shape`/
//! `values`, not for byte-for-byte footer preservation), not an oversight.
//!
//! ## Out of scope
//!
//! VASP's CHGCAR/LOCPOT volumetric formats are a separate, later roadmap
//! step, as in `crate::cube`.

use crate::volumetric::{GridError, GridUnits, LineFeed, VolumetricGrid};

// ---------------------------------------------------------------------------
// Limits
// ---------------------------------------------------------------------------

/// Security/robustness limits enforced before any large allocation happens.
#[derive(Debug, Clone, Copy)]
pub struct OpenDxParseLimits {
    /// Cumulative byte budget across the whole input.
    pub max_input_bytes: usize,
    /// Cap on `shape[0] * shape[1] * shape[2]`, checked from the header
    /// *before* the (potentially huge) data block is read.
    pub max_grid_points: usize,
}

impl Default for OpenDxParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 1 << 30, // 1 GiB
            max_grid_points: 100_000_000,
        }
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors from [`parse_opendx`]/[`parse_opendx_with_limits`]/[`write_opendx`].
#[derive(Debug, Clone, PartialEq)]
pub enum OpenDxError {
    /// Input exceeded [`OpenDxParseLimits::max_input_bytes`].
    InputTooLarge { limit: usize },
    /// The input ended before a required section was fully read.
    UnexpectedEnd { context: &'static str },
    /// A required header line (`object 1 ...`/`origin`/`delta`/
    /// `object 2 ...`/`object 3 ...`) was missing or didn't match the
    /// expected fixed keyword shape.
    MalformedHeaderLine {
        line: usize,
        context: &'static str,
        detail: String,
    },
    /// `object 2 class gridconnections counts` didn't match `object 1`'s
    /// counts.
    GridConnectionsMismatch {
        gridpositions: [usize; 3],
        gridconnections: [usize; 3],
    },
    /// `object 3`'s `class`/`type`/`rank` were not exactly
    /// `array`/`double`/`0` -- outside this module's documented APBS
    /// scalar-field-only scope.
    UnsupportedArrayDeclaration { line: usize, detail: String },
    /// `object 3`'s declared `items` count did not equal
    /// `nx * ny * nz` from `object 1`/`object 2`.
    ItemsCountMismatch { declared: usize, expected: usize },
    /// A numeric token could not be parsed as a float at all.
    InvalidNumber {
        line: usize,
        context: &'static str,
        raw: String,
    },
    /// A numeric token parsed but was NaN/Infinite.
    NonFiniteValue {
        line: usize,
        context: &'static str,
        raw: String,
    },
    /// The data block had fewer whitespace-delimited numeric tokens than
    /// `items` declared.
    ValueCountMismatch { expected: usize, found: usize },
    /// The data block had at least one more *numeric* token than `items`
    /// declared (a non-numeric trailing token is treated as the start of
    /// the out-of-scope footer -- see module docs -- and is not an error).
    TrailingData { after_values: usize },
    /// [`write_opendx`] was given a grid tagged
    /// [`crate::volumetric::GridUnits::Bohr`]. OpenDX carries no unit tag
    /// of its own and every real-world DX consumer (APBS included)
    /// assumes Ångström (see module docs) -- writing Bohr-magnitude
    /// numbers into a DX file would silently reinterpret them as Ångström
    /// on the next read (a real, silent ~1.89x error), so the default
    /// writer refuses rather than doing that implicitly. Use
    /// [`write_opendx_lossy`] to opt into an explicit Bohr->Ångström
    /// conversion instead.
    NonAngstromUnits { units: GridUnits },
    /// [`write_opendx`]/[`write_opendx_lossy`] was given a grid with a
    /// non-empty `atoms` list. OpenDX has no atom section at all (see
    /// module docs) -- there is no lossy-but-acceptable way to write atom
    /// data into this format, so both writers refuse rather than silently
    /// dropping it. Clear `grid.atoms` first if that loss is intentional.
    AtomsNotSupported { count: usize },
    /// A [`VolumetricGrid`] structural-invariant check failed -- raised by
    /// [`write_opendx`], since a hand-built `VolumetricGrid` is not
    /// guaranteed to satisfy them.
    Grid(GridError),
}

impl std::fmt::Display for OpenDxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InputTooLarge { limit } => write!(f, "OpenDX input exceeds {limit}-byte limit"),
            Self::UnexpectedEnd { context } => {
                write!(f, "unexpected end of input while reading {context}")
            }
            Self::MalformedHeaderLine {
                line,
                context,
                detail,
            } => write!(f, "malformed {context} line at line {line}: '{detail}'"),
            Self::GridConnectionsMismatch {
                gridpositions,
                gridconnections,
            } => write!(
                f,
                "object 2 gridconnections counts {gridconnections:?} do not match object 1 gridpositions counts {gridpositions:?}"
            ),
            Self::UnsupportedArrayDeclaration { line, detail } => write!(
                f,
                "unsupported object 3 array declaration at line {line}: '{detail}' (only 'array type double rank 0' is supported)"
            ),
            Self::ItemsCountMismatch { declared, expected } => write!(
                f,
                "object 3 declares {declared} items but the grid shape implies {expected}"
            ),
            Self::InvalidNumber { line, context, raw } => {
                write!(f, "invalid {context} value '{raw}' at line {line}")
            }
            Self::NonFiniteValue { line, context, raw } => write!(
                f,
                "{context} value '{raw}' at line {line} is not finite (NaN/Infinite)"
            ),
            Self::ValueCountMismatch { expected, found } => write!(
                f,
                "OpenDX data block expects {expected} values, only found {found}"
            ),
            Self::TrailingData { after_values } => write!(
                f,
                "OpenDX data block has extra numeric values after the expected {after_values}"
            ),
            Self::NonAngstromUnits { units } => write!(
                f,
                "cannot write a {units:?}-unit grid to OpenDX (which has no unit tag and is universally read as Angstrom) without an explicit conversion -- use write_opendx_lossy"
            ),
            Self::AtomsNotSupported { count } => write!(
                f,
                "OpenDX has no atom section; grid carries {count} atom(s) that would be silently dropped"
            ),
            Self::Grid(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for OpenDxError {}

impl From<GridError> for OpenDxError {
    fn from(e: GridError) -> Self {
        OpenDxError::Grid(e)
    }
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

fn expect_tokens<'a>(
    line: &'a str,
    line_no: usize,
    context: &'static str,
    exact: &[&str],
) -> Result<Vec<&'a str>, OpenDxError> {
    let toks: Vec<&str> = line.split_whitespace().collect();
    if toks.len() < exact.len()
        || toks[..exact.len()]
            .iter()
            .zip(exact.iter())
            .any(|(a, b)| a != b)
    {
        return Err(OpenDxError::MalformedHeaderLine {
            line: line_no,
            context,
            detail: line.to_string(),
        });
    }
    Ok(toks)
}

fn parse_finite(raw: &str, line: usize, context: &'static str) -> Result<f64, OpenDxError> {
    let v: f64 = raw.parse().map_err(|_| OpenDxError::InvalidNumber {
        line,
        context,
        raw: raw.to_string(),
    })?;
    if !v.is_finite() {
        return Err(OpenDxError::NonFiniteValue {
            line,
            context,
            raw: raw.to_string(),
        });
    }
    Ok(v)
}

fn parse_usize_triplet(
    toks: &[&str],
    start: usize,
    line: usize,
    context: &'static str,
) -> Result<[usize; 3], OpenDxError> {
    let mut out = [0usize; 3];
    for (c, slot) in out.iter_mut().enumerate() {
        *slot = toks
            .get(start + c)
            .and_then(|t| t.parse::<usize>().ok())
            .ok_or_else(|| OpenDxError::MalformedHeaderLine {
                line,
                context,
                detail: toks.join(" "),
            })?;
    }
    Ok(out)
}

fn parse_opendx_from_feed<F>(
    next_line: F,
    limits: &OpenDxParseLimits,
) -> Result<VolumetricGrid, OpenDxError>
where
    F: FnMut() -> Result<Option<String>, OpenDxError>,
{
    let mut feed = LineFeed::new(next_line);

    // Skip any number of leading '#' comment lines (real-world APBS output
    // convention; see module docs -- medium confidence, not in the primary
    // grammar source itself).
    let object1_line = loop {
        let line = feed.line()?.ok_or(OpenDxError::UnexpectedEnd {
            context: "object 1 (gridpositions) line",
        })?;
        if line.trim_start().starts_with('#') {
            continue;
        }
        break line;
    };
    let line_no = feed.line_no;
    let toks = expect_tokens(
        &object1_line,
        line_no,
        "object 1 (gridpositions)",
        &["object", "1", "class", "gridpositions", "counts"],
    )?;
    let gridpositions = parse_usize_triplet(&toks, 5, line_no, "object 1 (gridpositions)")?;

    let expected_points = gridpositions[0]
        .checked_mul(gridpositions[1])
        .and_then(|v| v.checked_mul(gridpositions[2]))
        .ok_or(GridError::ShapeOverflow {
            shape: gridpositions,
        })?;
    if expected_points > limits.max_grid_points {
        return Err(GridError::GridTooLarge {
            points: expected_points,
            limit: limits.max_grid_points,
        }
        .into());
    }

    // origin
    let origin_line = feed.line()?.ok_or(OpenDxError::UnexpectedEnd {
        context: "origin line",
    })?;
    let line_no = feed.line_no;
    let toks = expect_tokens(&origin_line, line_no, "origin", &["origin"])?;
    let mut origin = [0.0f64; 3];
    for (c, slot) in origin.iter_mut().enumerate() {
        *slot = parse_finite(
            toks.get(1 + c).ok_or(OpenDxError::MalformedHeaderLine {
                line: line_no,
                context: "origin",
                detail: origin_line.clone(),
            })?,
            line_no,
            "origin",
        )?;
    }

    // 3 delta lines -- general 3-vectors, not required to be axis-aligned.
    let mut axes = [[0.0f64; 3]; 3];
    for row in axes.iter_mut() {
        let delta_line = feed.line()?.ok_or(OpenDxError::UnexpectedEnd {
            context: "delta line",
        })?;
        let line_no = feed.line_no;
        let toks = expect_tokens(&delta_line, line_no, "delta", &["delta"])?;
        for (c, slot) in row.iter_mut().enumerate() {
            *slot = parse_finite(
                toks.get(1 + c).ok_or(OpenDxError::MalformedHeaderLine {
                    line: line_no,
                    context: "delta",
                    detail: delta_line.clone(),
                })?,
                line_no,
                "delta",
            )?;
        }
    }

    // object 2: gridconnections, must match object 1's counts.
    let object2_line = feed.line()?.ok_or(OpenDxError::UnexpectedEnd {
        context: "object 2 (gridconnections) line",
    })?;
    let line_no = feed.line_no;
    let toks = expect_tokens(
        &object2_line,
        line_no,
        "object 2 (gridconnections)",
        &["object", "2", "class", "gridconnections", "counts"],
    )?;
    let gridconnections = parse_usize_triplet(&toks, 5, line_no, "object 2 (gridconnections)")?;
    if gridconnections != gridpositions {
        return Err(OpenDxError::GridConnectionsMismatch {
            gridpositions,
            gridconnections,
        });
    }

    // object 3: array type double rank 0 items N data follows
    let object3_line = feed.line()?.ok_or(OpenDxError::UnexpectedEnd {
        context: "object 3 (array) line",
    })?;
    let line_no = feed.line_no;
    let toks: Vec<&str> = object3_line.split_whitespace().collect();
    let expected_prefix = [
        "object", "3", "class", "array", "type", "double", "rank", "0",
    ];
    if toks.len() < expected_prefix.len()
        || toks[..expected_prefix.len()]
            .iter()
            .zip(expected_prefix.iter())
            .any(|(a, b)| a != b)
    {
        return Err(OpenDxError::UnsupportedArrayDeclaration {
            line: line_no,
            detail: object3_line,
        });
    }
    if toks.get(expected_prefix.len()) != Some(&"items")
        || toks.get(expected_prefix.len() + 2) != Some(&"data")
        || toks.get(expected_prefix.len() + 3) != Some(&"follows")
    {
        return Err(OpenDxError::MalformedHeaderLine {
            line: line_no,
            context: "object 3 (array) items/data follows",
            detail: object3_line,
        });
    }
    let declared_items: usize = toks
        .get(expected_prefix.len() + 1)
        .and_then(|t| t.parse().ok())
        .ok_or(OpenDxError::MalformedHeaderLine {
            line: line_no,
            context: "object 3 (array) items count",
            detail: object3_line.clone(),
        })?;
    if declared_items != expected_points {
        return Err(OpenDxError::ItemsCountMismatch {
            declared: declared_items,
            expected: expected_points,
        });
    }

    // Data block: exactly `expected_points` whitespace-delimited numeric
    // tokens, tolerant of any line-wrapping.
    let mut values = Vec::with_capacity(expected_points.min(1_000_000));
    for i in 0..expected_points {
        let tok = feed.token()?.ok_or(OpenDxError::ValueCountMismatch {
            expected: expected_points,
            found: i,
        })?;
        let line_no = feed.line_no;
        let v = parse_finite(&tok, line_no, "data value")?;
        values.push(v);
    }
    // Anything left that still parses as a number is too much data; a
    // non-numeric leftover token is the start of the (unsupported, see
    // module docs) footer and is not an error.
    if let Some(tok) = feed.token()?
        && tok.parse::<f64>().is_ok()
    {
        return Err(OpenDxError::TrailingData {
            after_values: expected_points,
        });
    }

    Ok(VolumetricGrid {
        origin,
        axes,
        shape: gridpositions,
        values,
        atoms: Vec::new(),
        units: GridUnits::Angstrom,
    })
}

/// Parse an OpenDX (APBS scalar-field subset) file with default limits
/// ([`OpenDxParseLimits::default`]).
pub fn parse_opendx(input: &str) -> Result<VolumetricGrid, OpenDxError> {
    parse_opendx_with_limits(input, &OpenDxParseLimits::default())
}

/// Parse an OpenDX (APBS scalar-field subset) file, enforcing `limits`.
pub fn parse_opendx_with_limits(
    input: &str,
    limits: &OpenDxParseLimits,
) -> Result<VolumetricGrid, OpenDxError> {
    if input.len() > limits.max_input_bytes {
        return Err(OpenDxError::InputTooLarge {
            limit: limits.max_input_bytes,
        });
    }
    let mut lines = input.lines();
    let next_line =
        move || -> Result<Option<String>, OpenDxError> { Ok(lines.next().map(|s| s.to_string())) };
    parse_opendx_from_feed(next_line, limits)
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// CODATA 2018 Bohr radius in Angstrom, matching
/// `qcschema.rs::BOHR_TO_ANGSTROM` (kept as a separate local constant
/// rather than importing across format modules -- see that module for the
/// same value's own derivation note).
const BOHR_TO_ANGSTROM: f64 = 0.529177210903;

/// Write a [`VolumetricGrid`] as an OpenDX (APBS scalar-field subset) file.
///
/// Fails closed rather than silently losing information: OpenDX has no
/// unit tag of its own (every real-world consumer, APBS included, assumes
/// Ångström -- see module docs) and no atom section at all, so this
/// refuses with a typed [`OpenDxError::NonAngstromUnits`] for a
/// [`crate::volumetric::GridUnits::Bohr`] grid (writing its numbers
/// as-is would have them silently reinterpreted as Ångström -- a real,
/// silent ~1.89x error -- the next time the file is read) and with
/// [`OpenDxError::AtomsNotSupported`] for a grid with any atoms (there is
/// no lossy-but-acceptable way to write those). Use
/// [`write_opendx_lossy`] to opt into an explicit Bohr->Ångström
/// conversion; atoms are never silently dropped by either writer -- clear
/// `grid.atoms` yourself first if that loss is intentional.
///
/// The trailing `attribute`/`object "field is ..."` footer is always the
/// standard boilerplate, not reconstructed from any source file.
pub fn write_opendx(grid: &VolumetricGrid) -> Result<String, OpenDxError> {
    grid.validate()?;
    if grid.units != GridUnits::Angstrom {
        return Err(OpenDxError::NonAngstromUnits { units: grid.units });
    }
    if !grid.atoms.is_empty() {
        return Err(OpenDxError::AtomsNotSupported {
            count: grid.atoms.len(),
        });
    }
    Ok(write_opendx_body(grid, grid.origin, grid.axes))
}

/// Like [`write_opendx`], but if `grid.units` is
/// [`crate::volumetric::GridUnits::Bohr`], explicitly converts `origin`
/// and `axes` (length quantities) to Ångström before writing, rather than
/// refusing. `values` (the scalar-field samples themselves -- e.g.
/// electron density or an electrostatic potential) are **never** rescaled
/// by this conversion: only `origin`/`axes` are lengths, and this module
/// has no way to know what physical unit the field values are in. Still
/// refuses with [`OpenDxError::AtomsNotSupported`] for a grid with any
/// atoms, same as [`write_opendx`] -- the unit conversion is the only
/// thing this function does that the default writer won't.
pub fn write_opendx_lossy(grid: &VolumetricGrid) -> Result<String, OpenDxError> {
    grid.validate()?;
    if !grid.atoms.is_empty() {
        return Err(OpenDxError::AtomsNotSupported {
            count: grid.atoms.len(),
        });
    }
    let (origin, axes) = if grid.units == GridUnits::Bohr {
        (
            grid.origin.map(|v| v * BOHR_TO_ANGSTROM),
            grid.axes.map(|row| row.map(|v| v * BOHR_TO_ANGSTROM)),
        )
    } else {
        (grid.origin, grid.axes)
    };
    Ok(write_opendx_body(grid, origin, axes))
}

/// Shared serialization body. `origin`/`axes` are passed separately
/// (rather than read from `grid`) so [`write_opendx_lossy`] can supply a
/// Bohr->Ångström-converted pair without cloning the rest of the grid
/// (`values` in particular can be large). Callers must already have
/// validated the grid and resolved the units/atoms preconditions
/// [`write_opendx`] and [`write_opendx_lossy`] each enforce their own way.
fn write_opendx_body(grid: &VolumetricGrid, origin: [f64; 3], axes: [[f64; 3]; 3]) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "object 1 class gridpositions counts {} {} {}\n",
        grid.shape[0], grid.shape[1], grid.shape[2]
    ));
    out.push_str(&format!(
        "origin {} {} {}\n",
        origin[0], origin[1], origin[2]
    ));
    for row in &axes {
        out.push_str(&format!("delta {} {} {}\n", row[0], row[1], row[2]));
    }
    out.push_str(&format!(
        "object 2 class gridconnections counts {} {} {}\n",
        grid.shape[0], grid.shape[1], grid.shape[2]
    ));
    out.push_str(&format!(
        "object 3 class array type double rank 0 items {} data follows\n",
        grid.values.len()
    ));
    for (idx, v) in grid.values.iter().enumerate() {
        out.push_str(&v.to_string());
        if (idx + 1) % 3 == 0 {
            out.push('\n');
        } else {
            out.push(' ');
        }
    }
    if grid.values.is_empty() || !grid.values.len().is_multiple_of(3) {
        out.push('\n');
    }
    out.push_str("attribute \"dep\" string \"positions\"\n");
    out.push_str("object \"regular positions regular connections\" class field\n");
    out.push_str("component \"positions\" value 1\n");
    out.push_str("component \"connections\" value 2\n");
    out.push_str("component \"data\" value 3\n");
    out
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn small_2x2x2() -> VolumetricGrid {
        VolumetricGrid {
            origin: [-1.0, -1.0, -1.0],
            axes: [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]],
            shape: [2, 2, 2],
            values: vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            atoms: Vec::new(),
            units: GridUnits::Angstrom,
        }
    }

    const HAND_WRITTEN_2X2X2_DX: &str = "# Data from APBS\n\
# POTENTIAL (kT/e)\n\
object 1 class gridpositions counts 2 2 2\n\
origin -1.0 -1.0 -1.0\n\
delta 0.5 0.0 0.0\n\
delta 0.0 0.5 0.0\n\
delta 0.0 0.0 0.5\n\
object 2 class gridconnections counts 2 2 2\n\
object 3 class array type double rank 0 items 8 data follows\n\
0.0 1.0 2.0\n\
3.0 4.0 5.0\n\
6.0 7.0\n\
attribute \"dep\" string \"positions\"\n\
object \"regular positions regular connections\" class field\n\
component \"positions\" value 1\n\
component \"connections\" value 2\n\
component \"data\" value 3\n";

    #[test]
    fn parse_hand_written_2x2x2_fixture_with_leading_comments() {
        let grid = parse_opendx(HAND_WRITTEN_2X2X2_DX).unwrap();
        assert_eq!(grid.shape, [2, 2, 2]);
        assert_eq!(grid.origin, [-1.0, -1.0, -1.0]);
        assert_eq!(grid.values, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
        assert!(grid.atoms.is_empty());
        assert_eq!(grid.units, GridUnits::Angstrom);
    }

    #[test]
    fn round_trip_2x2x2() {
        let grid = small_2x2x2();
        let text = write_opendx(&grid).unwrap();
        let parsed = parse_opendx(&text).unwrap();
        assert_eq!(parsed, grid);
    }

    #[test]
    fn round_trip_3x2x4_non_cubic_shape() {
        let n = 3 * 2 * 4;
        let grid = VolumetricGrid {
            origin: [0.0, 0.0, 0.0],
            axes: [[0.2, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.1]],
            shape: [3, 2, 4],
            values: (0..n).map(|i| i as f64 * 0.25).collect(),
            atoms: Vec::new(),
            units: GridUnits::Angstrom,
        };
        let text = write_opendx(&grid).unwrap();
        let parsed = parse_opendx(&text).unwrap();
        assert_eq!(parsed, grid);
    }

    #[test]
    fn round_trip_non_orthogonal_delta_vectors() {
        let mut grid = small_2x2x2();
        grid.axes = [[0.5, 0.1, 0.0], [0.05, 0.5, 0.1], [0.0, 0.05, 0.5]];
        let text = write_opendx(&grid).unwrap();
        let parsed = parse_opendx(&text).unwrap();
        assert_eq!(parsed.axes, grid.axes);
        assert_eq!(parsed, grid);
    }

    #[test]
    fn nan_data_value_is_typed_rejected() {
        let input = HAND_WRITTEN_2X2X2_DX.replace("0.0 1.0 2.0\n", "NaN 1.0 2.0\n");
        assert!(input.contains("NaN"));
        let err = parse_opendx(&input).unwrap_err();
        assert!(matches!(
            err,
            OpenDxError::NonFiniteValue {
                context: "data value",
                ..
            }
        ));
    }

    #[test]
    fn infinite_origin_is_typed_rejected() {
        let input = HAND_WRITTEN_2X2X2_DX.replace("origin -1.0 -1.0 -1.0", "origin inf -1.0 -1.0");
        let err = parse_opendx(&input).unwrap_err();
        assert!(matches!(
            err,
            OpenDxError::NonFiniteValue {
                context: "origin",
                ..
            }
        ));
    }

    #[test]
    fn gridconnections_mismatch_is_typed_error() {
        let input = HAND_WRITTEN_2X2X2_DX.replace(
            "object 2 class gridconnections counts 2 2 2",
            "object 2 class gridconnections counts 2 2 3",
        );
        let err = parse_opendx(&input).unwrap_err();
        assert_eq!(
            err,
            OpenDxError::GridConnectionsMismatch {
                gridpositions: [2, 2, 2],
                gridconnections: [2, 2, 3],
            }
        );
    }

    #[test]
    fn items_count_mismatch_is_typed_error() {
        let input = HAND_WRITTEN_2X2X2_DX.replace(
            "object 3 class array type double rank 0 items 8 data follows",
            "object 3 class array type double rank 0 items 999 data follows",
        );
        let err = parse_opendx(&input).unwrap_err();
        assert_eq!(
            err,
            OpenDxError::ItemsCountMismatch {
                declared: 999,
                expected: 8,
            }
        );
    }

    #[test]
    fn unsupported_rank_is_typed_rejected() {
        let input = HAND_WRITTEN_2X2X2_DX.replace(
            "object 3 class array type double rank 0 items 8 data follows",
            "object 3 class array type double rank 1 shape 3 items 8 data follows",
        );
        let err = parse_opendx(&input).unwrap_err();
        assert!(matches!(
            err,
            OpenDxError::UnsupportedArrayDeclaration { .. }
        ));
    }

    #[test]
    fn truncated_input_is_typed_error_not_panic() {
        let input = "object 1 class gridpositions counts 2 2 2\norigin -1.0 -1.0 -1.0\n";
        let err = parse_opendx(input).unwrap_err();
        assert!(matches!(err, OpenDxError::UnexpectedEnd { .. }));
    }

    #[test]
    fn too_few_values_is_typed_mismatch() {
        // Genuinely truncated (no footer following) -- if a footer were
        // present, its first non-numeric token would surface as an
        // `InvalidNumber` instead, which is a *different*, also-typed
        // failure mode this format's lack of an explicit data-block
        // terminator can produce; this test isolates the pure "ran out of
        // input entirely" case.
        let input = "object 1 class gridpositions counts 2 2 2\n\
origin -1.0 -1.0 -1.0\n\
delta 0.5 0.0 0.0\n\
delta 0.0 0.5 0.0\n\
delta 0.0 0.0 0.5\n\
object 2 class gridconnections counts 2 2 2\n\
object 3 class array type double rank 0 items 8 data follows\n\
0.0 1.0 2.0 3.0 4.0 5.0 6.0\n";
        let err = parse_opendx(input).unwrap_err();
        assert_eq!(
            err,
            OpenDxError::ValueCountMismatch {
                expected: 8,
                found: 7
            }
        );
    }

    #[test]
    fn too_many_values_is_typed_trailing_data_error() {
        let input = HAND_WRITTEN_2X2X2_DX.replace("6.0 7.0\n", "6.0 7.0 8.0\n");
        let err = parse_opendx(&input).unwrap_err();
        assert_eq!(err, OpenDxError::TrailingData { after_values: 8 });
    }

    #[test]
    fn pathological_header_overflow_is_typed_error() {
        let input = "object 1 class gridpositions counts 4000000000 4000000000 2\n";
        let err = parse_opendx(input).unwrap_err();
        assert!(matches!(
            err,
            OpenDxError::Grid(GridError::ShapeOverflow { .. })
        ));
    }

    #[test]
    fn pathological_header_within_usize_but_over_cap_is_typed_error() {
        let input = "object 1 class gridpositions counts 10000 10000 10000\n";
        let err = parse_opendx(input).unwrap_err();
        assert!(matches!(
            err,
            OpenDxError::Grid(GridError::GridTooLarge {
                points: 1_000_000_000_000,
                limit: 100_000_000
            })
        ));
    }

    #[test]
    fn write_rejects_nan_via_validate() {
        let mut grid = small_2x2x2();
        grid.values[0] = f64::NAN;
        let err = write_opendx(&grid).unwrap_err();
        assert!(matches!(
            err,
            OpenDxError::Grid(GridError::NonFiniteValue { .. })
        ));
    }

    #[test]
    fn footer_is_not_required_to_round_trip() {
        // write_opendx always emits its own standard footer; a hand-parsed
        // file whose footer differs (or is entirely absent) must still
        // parse, since the footer is documented out of scope.
        let no_footer = HAND_WRITTEN_2X2X2_DX
            .lines()
            .take_while(|l| !l.starts_with("attribute"))
            .collect::<Vec<_>>()
            .join("\n");
        let grid = parse_opendx(&no_footer).unwrap();
        assert_eq!(grid.shape, [2, 2, 2]);
    }

    #[test]
    fn write_rejects_bohr_units_rather_than_silently_reinterpreting() {
        // Writing Bohr-magnitude numbers into a format that's always read
        // back as Angstrom would be a real, silent ~1.89x error -- the
        // default writer must refuse, not convert or pass through.
        let mut grid = small_2x2x2();
        grid.units = GridUnits::Bohr;
        let err = write_opendx(&grid).unwrap_err();
        assert_eq!(
            err,
            OpenDxError::NonAngstromUnits {
                units: GridUnits::Bohr
            }
        );
    }

    #[test]
    fn write_rejects_nonempty_atoms() {
        let mut grid = small_2x2x2();
        grid.atoms.push(crate::volumetric::GridAtom {
            element: chematic_core::Element::from_symbol("C").unwrap(),
            charge: 6.0,
            position: [0.0, 0.0, 0.0],
        });
        let err = write_opendx(&grid).unwrap_err();
        assert_eq!(err, OpenDxError::AtomsNotSupported { count: 1 });
    }

    #[test]
    fn write_lossy_converts_bohr_origin_and_axes_but_not_values() {
        let mut grid = small_2x2x2();
        grid.units = GridUnits::Bohr;
        grid.origin = [1.0, 0.0, 0.0];
        grid.axes = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
        let text = write_opendx_lossy(&grid).unwrap();
        let parsed = parse_opendx(&text).unwrap();
        // Geometry is converted...
        assert!((parsed.origin[0] - 1.0 * BOHR_TO_ANGSTROM).abs() < 1e-9);
        assert!((parsed.axes[0][0] - 2.0 * BOHR_TO_ANGSTROM).abs() < 1e-9);
        // ...but the scalar-field samples themselves are untouched (not a
        // length quantity, so never rescaled by this conversion).
        assert_eq!(parsed.values, grid.values);
    }

    #[test]
    fn write_lossy_still_rejects_atoms() {
        let mut grid = small_2x2x2();
        grid.units = GridUnits::Bohr;
        grid.atoms.push(crate::volumetric::GridAtom {
            element: chematic_core::Element::from_symbol("C").unwrap(),
            charge: 6.0,
            position: [0.0, 0.0, 0.0],
        });
        let err = write_opendx_lossy(&grid).unwrap_err();
        assert_eq!(err, OpenDxError::AtomsNotSupported { count: 1 });
    }

    #[test]
    fn write_rejects_zero_dimension_shape() {
        let mut grid = small_2x2x2();
        grid.shape = [0, 2, 2];
        grid.values = Vec::new();
        let err = write_opendx(&grid).unwrap_err();
        assert_eq!(err, OpenDxError::Grid(GridError::ZeroDimension { axis: 0 }));
    }
}
