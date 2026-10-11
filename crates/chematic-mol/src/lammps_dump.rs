//! LAMMPS dump/trajectory file (`dump` command's default text style)
//! reader and writer.
//!
//! ## Provenance
//!
//! Implemented independently against LAMMPS's own official manual
//! (<https://docs.lammps.org/dump.html>, <https://docs.lammps.org/Howto_triclinic.html>).
//! Both pages were fetched and their content confirmed directly against
//! the live pages (each fetch mediated by a summarization step, not a raw
//! byte-for-byte scrape -- formulas quoted below were independently
//! cross-checked by re-fetching the specific sentence/formula in
//! isolation, but are presented here as *confirmed*, not as a guaranteed
//! character-exact transcription). No source code, comments, or tables
//! were copied from LAMMPS, VMD, OVITO, MDAnalysis, or any other tool.
//!
//! ## Structure
//!
//! Each frame is:
//! ```text
//! ITEM: TIMESTEP
//! <int>
//! ITEM: NUMBER OF ATOMS
//! <int>
//! ITEM: BOX BOUNDS <flags>
//! <xlo_or_bound> <xhi_or_bound> [<xy>]
//! <ylo_or_bound> <yhi_or_bound> [<xz>]
//! <zlo_or_bound> <zhi_or_bound> [<yz>]
//! ITEM: ATOMS <column-name-list>
//! <one row per atom>
//! ```
//! A trajectory file is simply N frames concatenated back to back, with
//! no separator beyond the next `ITEM: TIMESTEP`.
//!
//! ## Box bounds: orthogonal vs. triclinic
//!
//! The `ITEM: BOX BOUNDS` line's own trailing tokens say which shape
//! follows: dump.html's own example headers are
//! `ITEM: BOX BOUNDS xx yy zz` (orthogonal -- `xx`/`yy`/`zz` are 2-
//! character boundary-condition flag pairs like `pp`, e.g.
//! `ITEM: BOX BOUNDS pp pp pp`) versus
//! `ITEM: BOX BOUNDS xy xz yz xx yy zz` (triclinic -- the literal tokens
//! `xy xz yz` appear before the same 3 flag pairs, e.g.
//! `ITEM: BOX BOUNDS xy xz yz pp pp pp`), with the 3 data rows then
//! carrying `xlo_bound xhi_bound xy` / `ylo_bound yhi_bound xz` /
//! `zlo_bound zhi_bound yz` instead of plain `xlo xhi` / `ylo yhi` /
//! `zlo zhi`. This reader detects triclinic-ness from the header line's
//! `xy xz yz` tokens (not from counting the 3 data rows' own column
//! count) but cross-checks both signals -- a row with 3 columns after an
//! orthogonal header, or 2 columns after a triclinic header, is a typed
//! [`LammpsDumpError::InvalidBox`] (indicating a corrupt/hand-edited
//! file), not silently reinterpreted.
//!
//! ## Triclinic bounding-box <-> true-box conversion
//!
//! Confirmed against <https://docs.lammps.org/Howto_triclinic.html>,
//! which states the bounding-box extents are computed from the true box
//! as:
//! ```text
//! xlo_bound = xlo + MIN(0.0, xy, xz, xy+xz)
//! xhi_bound = xhi + MAX(0.0, xy, xz, xy+xz)
//! ylo_bound = ylo + MIN(0.0, yz)
//! yhi_bound = yhi + MAX(0.0, yz)
//! zlo_bound = zlo
//! zhi_bound = zhi
//! ```
//! and that this can be inverted (the page gives the `xlo` case
//! explicitly: `xlo = xlo_bound - MIN(0.0, xy, xz, xy+xz)`) to recover the
//! true box from a dump file's bounding-box values -- this is exactly the
//! formula this brief specified, independently confirmed rather than
//! taken on trust; see [`box_bounds_to_true`]/[`true_to_box_bounds`] and
//! this module's tests (both a bound->true->bound identity check over
//! several distinct tilt fixtures *and* a hand-computed-absolute-value
//! check, since the identity check alone cannot distinguish a correct
//! `MIN`/`MAX` formula from a subtly wrong one that happens to be its own
//! inverse).
//!
//! ## Scaled coordinates (`xs ys zs`)
//!
//! Per the same Howto_triclinic page, a restricted triclinic box's edge
//! vectors are **A** = `(xhi-xlo, 0, 0)`, **B** = `(xy, yhi-ylo, 0)`,
//! **C** = `(xz, yz, zhi-zlo)`, with origin `(xlo, ylo, zlo)`; dump.html
//! states "the actual unscaled (x,y,z) coordinate is `xs*A + ys*B +
//! zs*C`". Taken completely literally that sentence omits the box
//! origin, but it cannot be correct without it: dump.html separately
//! defines `xs`/`ys`/`zs` as "scaled to the box size so that each value
//! is 0.0 to 1.0", and `xs=ys=zs=0` must therefore map to the box's own
//! `(xlo,ylo,zlo)` corner, not to the coordinate origin `(0,0,0)`.
//! [`LammpsDumpFrame::cartesian_positions`] therefore computes
//! `(xlo,ylo,zlo) + xs*A + ys*B + zs*C`, expanding to:
//! ```text
//! x = xlo + xs*(xhi-xlo) + ys*xy + zs*xz
//! y = ylo + ys*(yhi-ylo) + zs*yz
//! z = zlo + zs*(zhi-zlo)
//! ```
//! (reducing to ordinary independent per-axis scaling when `tilt` is
//! `None`). **Medium confidence**: the `A`/`B`/`C` edge-vector
//! definitions and the `xs*A+ys*B+zs*C` sentence are both confirmed
//! against the live docs; the origin term is this module's own necessary
//! completion of that sentence, not a verbatim-quoted part of it --
//! flagged explicitly per this project's stop-and-report-rather-than-
//! silently-guess convention (see `crate::cube`'s module docs for the
//! same discipline applied to Cube's Bohr/Ångström ambiguity). Verified
//! against a hand-computed example in this module's tests.
//!
//! `xu yu zu` ("unwrapped") coordinates need no transform -- "unwrapped"
//! specifically means not passed through periodic-boundary wrapping (an
//! atom that crossed a periodic boundary many times can have an
//! unwrapped coordinate far outside the visible box), a materially
//! different physical quantity from a wrapped Cartesian position, not
//! merely an unscaled one. [`LammpsDumpFrame::cartesian_positions`]
//! deliberately does not resolve `xu`/`yu`/`zu` itself (returns `None` if
//! neither `x y z` nor `xs ys zs` is present, even when `xu yu zu` is) --
//! callers wanting them should use [`LammpsDumpFrame::column`]`("xu")`
//! etc. directly, since they need no transform, and conflating "current
//! position" with "unwrapped position" into one accessor would be
//! misleading.
//!
//! ## Arbitrary, self-declared columns -- no fixed schema
//!
//! The `ITEM: ATOMS` column list varies per `dump` command invocation.
//! [`LammpsDumpFrame`] stores it as parallel `column_names`/`rows`
//! (row-major: `rows[i][j]` is column `j`'s value for atom `i`), not a
//! `HashMap<String, Vec<f64>>` -- this project's standing determinism
//! rule -- which exactly preserves original column names and order.
//!
//! ## Streaming reader
//!
//! [`LammpsDumpReader`] follows this crate's `SdfFileReader` precedent
//! (`crate::sdf`) -- `impl Iterator<Item = Result<LammpsDumpFrame,
//! LammpsDumpError>>` over any `BufRead` source -- rather than
//! `crate::cube::CubeFileReader`'s one-shot `.read()`, since a dump
//! trajectory is inherently multi-record (like SDF) while a Cube file is
//! a single grid.
//!
//! ## Box type
//!
//! Reuses [`crate::lammps_data::LammpsBox`] -- both formats describe the
//! same physical concept (an axis-aligned box plus 3 tilt factors). Only
//! *this* module performs the bound<->true conversion: a LAMMPS **data**
//! file's `xlo xhi`/`xy xz yz` header lines already give the true box
//! directly (data files have no bounding-box concept at all -- that only
//! exists in dump output, for tools that expect an axis-aligned
//! rendering box around a tilted cell).

use crate::lammps_data::LammpsBox;
use crate::volumetric::LineFeed;

// ---------------------------------------------------------------------------
// Box-bounds <-> true-box conversion
// ---------------------------------------------------------------------------

/// `(min(0, xy, xz, xy+xz), max(0, xy, xz, xy+xz))` -- the x-axis
/// bounding-box shift term. See module docs for the citation.
fn x_shift(xy: f64, xz: f64) -> (f64, f64) {
    let candidates = [0.0, xy, xz, xy + xz];
    let min = candidates.iter().copied().fold(f64::INFINITY, f64::min);
    let max = candidates.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (min, max)
}

/// `(min(0, yz), max(0, yz))` -- the y-axis bounding-box shift term.
fn y_shift(yz: f64) -> (f64, f64) {
    (0.0f64.min(yz), 0.0f64.max(yz))
}

/// Convert a dump file's `ITEM: BOX BOUNDS` values into the true
/// simulation box. `tilt` is `None` for an orthogonal box (in which case
/// bound and true coincide). See module docs for the formula's citation.
pub fn box_bounds_to_true(
    bound_lo: [f64; 3],
    bound_hi: [f64; 3],
    tilt: Option<[f64; 3]>,
) -> LammpsBox {
    match tilt {
        None => LammpsBox {
            lo: bound_lo,
            hi: bound_hi,
            tilt: None,
        },
        Some([xy, xz, yz]) => {
            let (xmin, xmax) = x_shift(xy, xz);
            let (ymin, ymax) = y_shift(yz);
            LammpsBox {
                lo: [bound_lo[0] - xmin, bound_lo[1] - ymin, bound_lo[2]],
                hi: [bound_hi[0] - xmax, bound_hi[1] - ymax, bound_hi[2]],
                tilt: Some([xy, xz, yz]),
            }
        }
    }
}

/// Inverse of [`box_bounds_to_true`]: convert a true [`LammpsBox`] into
/// the `xlo_bound/xhi_bound/...` values a dump file's `ITEM: BOX BOUNDS`
/// section would show.
pub fn true_to_box_bounds(b: &LammpsBox) -> ([f64; 3], [f64; 3]) {
    match b.tilt {
        None => (b.lo, b.hi),
        Some([xy, xz, yz]) => {
            let (xmin, xmax) = x_shift(xy, xz);
            let (ymin, ymax) = y_shift(yz);
            (
                [b.lo[0] + xmin, b.lo[1] + ymin, b.lo[2]],
                [b.hi[0] + xmax, b.hi[1] + ymax, b.hi[2]],
            )
        }
    }
}

/// Real Cartesian position from box-scaled coordinates, each
/// conventionally in `[0, 1]`. See module docs for the derivation.
fn scaled_to_cartesian(b: &LammpsBox, xs: f64, ys: f64, zs: f64) -> [f64; 3] {
    match b.tilt {
        None => [
            b.lo[0] + xs * (b.hi[0] - b.lo[0]),
            b.lo[1] + ys * (b.hi[1] - b.lo[1]),
            b.lo[2] + zs * (b.hi[2] - b.lo[2]),
        ],
        Some([xy, xz, yz]) => [
            b.lo[0] + xs * (b.hi[0] - b.lo[0]) + ys * xy + zs * xz,
            b.lo[1] + ys * (b.hi[1] - b.lo[1]) + zs * yz,
            b.lo[2] + zs * (b.hi[2] - b.lo[2]),
        ],
    }
}

// ---------------------------------------------------------------------------
// LammpsDumpFrame
// ---------------------------------------------------------------------------

/// One frame of a LAMMPS dump/trajectory file.
#[derive(Debug, Clone, PartialEq)]
pub struct LammpsDumpFrame {
    pub timestep: i64,
    pub num_atoms: usize,
    pub box_bounds: LammpsBox,
    /// The 2-character boundary-condition flag pair for each axis (x, y,
    /// z) -- e.g. `"pp"`, `"ff"`, `"ss"`, `"mm"`, or a mixed pair --
    /// taken verbatim from the `ITEM: BOX BOUNDS ... <flags>` line.
    /// Stored, not interpreted: this module has no periodic-wrapping
    /// logic of its own.
    pub boundary_flags: [String; 3],
    pub column_names: Vec<String>,
    /// `rows[i][j]` = value of `column_names[j]` for atom `i`, in file
    /// order. Row-major, index-parallel with `column_names` -- not a
    /// `HashMap<String, Vec<f64>>` -- see module docs.
    pub rows: Vec<Vec<f64>>,
}

/// Resource limits for LAMMPS dump frame and trajectory parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LammpsDumpParseLimits {
    /// Maximum input size, in bytes. For a streaming reader this is the
    /// cumulative number of bytes consumed from the source.
    pub max_input_bytes: usize,
    /// Maximum physical line size, in bytes.
    pub max_line_bytes: usize,
    /// Maximum atoms in one frame.
    pub max_atoms_per_frame: usize,
    /// Maximum declared atom columns in one frame.
    pub max_columns: usize,
    /// Maximum frames yielded by a streaming reader.
    pub max_frames: usize,
}

impl Default for LammpsDumpParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_line_bytes: 1024 * 1024,
            max_atoms_per_frame: 1_000_000,
            max_columns: 256,
            max_frames: 100_000,
        }
    }
}

impl LammpsDumpFrame {
    /// Index of `name` within [`Self::column_names`], if present.
    pub fn column_index(&self, name: &str) -> Option<usize> {
        self.column_names.iter().position(|c| c == name)
    }

    /// The values of column `name` across every atom, in file order.
    ///
    /// Returns an owned `Vec<f64>`, not a borrowed slice: [`Self::rows`]
    /// is row-major, so a single column is not stored contiguously and
    /// cannot be borrowed without a copy. This is a deliberate deviation
    /// from a `&[f64]`-returning signature (which would require
    /// column-oriented storage, the very layout `rows`/`column_names`'s
    /// row-major shape was chosen to avoid per this module's determinism
    /// rule) -- documented here rather than silently changed without
    /// comment.
    pub fn column(&self, name: &str) -> Option<Vec<f64>> {
        let idx = self.column_index(name)?;
        Some(self.rows.iter().map(|r| r[idx]).collect())
    }

    /// Real Cartesian positions per atom, in file order, from whichever
    /// of the recognized coordinate-triple conventions is present:
    /// - `x y z`: already real Cartesian, passed through as-is.
    /// - `xs ys zs`: box-scaled, transformed through the box (including
    ///   the triclinic shear terms when [`LammpsBox::tilt`] is `Some` --
    ///   see module docs; **not** simply independent per-axis scaling).
    ///
    /// Returns `None` if neither triple is fully present -- including
    /// when only `xu yu zu` ("unwrapped") is present; see module docs for
    /// why those are deliberately not resolved by this method. Never
    /// drops or overwrites `column_names`/`rows`: this is a read-only
    /// convenience view over them.
    pub fn cartesian_positions(&self) -> Option<Vec<[f64; 3]>> {
        if let (Some(x), Some(y), Some(z)) = (self.column("x"), self.column("y"), self.column("z"))
        {
            return Some(
                x.into_iter()
                    .zip(y)
                    .zip(z)
                    .map(|((x, y), z)| [x, y, z])
                    .collect(),
            );
        }
        if let (Some(xs), Some(ys), Some(zs)) =
            (self.column("xs"), self.column("ys"), self.column("zs"))
        {
            return Some(
                xs.into_iter()
                    .zip(ys)
                    .zip(zs)
                    .map(|((xs, ys), zs)| scaled_to_cartesian(&self.box_bounds, xs, ys, zs))
                    .collect(),
            );
        }
        None
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors from [`parse_lammps_dump_frame`]/[`LammpsDumpReader`].
#[derive(Debug, Clone, PartialEq)]
pub enum LammpsDumpError {
    /// A header value line (`TIMESTEP`/`NUMBER OF ATOMS` value, or a
    /// `BOX BOUNDS` line with the wrong token shape) failed to parse.
    MalformedHeader { line: usize, detail: String },
    /// An `ITEM:` line did not match what was expected at this point in
    /// the frame.
    UnexpectedItem { expected: String, found: String },
    /// The box bounds were inconsistent (wrong row column count for the
    /// declared orthogonal/triclinic-ness, a non-finite value, or
    /// `hi <= lo`).
    InvalidBox { reason: String },
    /// An `ATOMS` data row's field count did not match the number of
    /// declared columns.
    ColumnCountMismatch {
        declared: usize,
        actual: usize,
        line: usize,
    },
    /// A field failed to parse as a number at all.
    InvalidNumber {
        column: String,
        row: usize,
        raw: String,
    },
    /// A field parsed but was NaN/Infinite.
    NonFiniteValue { column: String, row: usize },
    /// `ITEM: NUMBER OF ATOMS` declared more atoms than the file actually
    /// had before hitting the next `ITEM:` line.
    AtomCountMismatch { declared: usize, actual: usize },
    /// The input ended before a required section/field was fully read.
    TruncatedInput { context: &'static str },
    /// An IO error occurred while reading from a streaming
    /// [`LammpsDumpReader`] source.
    Io(String),
    /// The input exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl std::fmt::Display for LammpsDumpError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MalformedHeader { line, detail } => {
                write!(f, "malformed header line at line {line}: '{detail}'")
            }
            Self::UnexpectedItem { expected, found } => {
                write!(f, "expected '{expected}', found '{found}'")
            }
            Self::InvalidBox { reason } => write!(f, "invalid LAMMPS dump box: {reason}"),
            Self::ColumnCountMismatch {
                declared,
                actual,
                line,
            } => write!(
                f,
                "line {line}: expected {declared} column(s), found {actual}"
            ),
            Self::InvalidNumber { column, row, raw } => write!(
                f,
                "invalid value '{raw}' for column '{column}' at atom row {row}"
            ),
            Self::NonFiniteValue { column, row } => write!(
                f,
                "non-finite value for column '{column}' at atom row {row}"
            ),
            Self::AtomCountMismatch { declared, actual } => write!(
                f,
                "ITEM: NUMBER OF ATOMS declared {declared} but only {actual} row(s) were present"
            ),
            Self::TruncatedInput { context } => {
                write!(f, "unexpected end of input while reading {context}")
            }
            Self::Io(msg) => write!(f, "IO error: {msg}"),
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(f, "{resource} has size {actual}, exceeding limit {limit}"),
        }
    }
}

impl std::error::Error for LammpsDumpError {}

// ---------------------------------------------------------------------------
// Parser core
// ---------------------------------------------------------------------------

fn expect_item_exact(line: &str, expected: &str) -> Result<(), LammpsDumpError> {
    let want = format!("ITEM: {expected}");
    if line.trim_end() != want {
        return Err(LammpsDumpError::UnexpectedItem {
            expected: want,
            found: line.to_string(),
        });
    }
    Ok(())
}

fn strip_item_prefix<'a>(line: &'a str, expected: &str) -> Result<&'a str, LammpsDumpError> {
    let prefix = format!("ITEM: {expected}");
    match line.trim_end().strip_prefix(&prefix) {
        Some(rest) => Ok(rest.trim()),
        None => Err(LammpsDumpError::UnexpectedItem {
            expected: prefix,
            found: line.to_string(),
        }),
    }
}

fn required_line<F>(
    feed: &mut LineFeed<LammpsDumpError, F>,
    context: &'static str,
) -> Result<String, LammpsDumpError>
where
    F: FnMut() -> Result<Option<String>, LammpsDumpError>,
{
    feed.line()?
        .ok_or(LammpsDumpError::TruncatedInput { context })
}

fn read_timestep<F>(feed: &mut LineFeed<LammpsDumpError, F>) -> Result<Option<i64>, LammpsDumpError>
where
    F: FnMut() -> Result<Option<String>, LammpsDumpError>,
{
    let Some(line) = feed.line()? else {
        return Ok(None);
    };
    expect_item_exact(&line, "TIMESTEP")?;
    let value = required_line(feed, "TIMESTEP value")?;
    value
        .trim()
        .parse()
        .map(Some)
        .map_err(|_| LammpsDumpError::MalformedHeader {
            line: feed.line_no,
            detail: value,
        })
}

fn read_atom_count<F>(
    feed: &mut LineFeed<LammpsDumpError, F>,
    limits: &LammpsDumpParseLimits,
) -> Result<usize, LammpsDumpError>
where
    F: FnMut() -> Result<Option<String>, LammpsDumpError>,
{
    expect_item_exact(
        &required_line(feed, "ITEM: NUMBER OF ATOMS")?,
        "NUMBER OF ATOMS",
    )?;
    let value = required_line(feed, "NUMBER OF ATOMS value")?;
    let count = value
        .trim()
        .parse::<usize>()
        .map_err(|_| LammpsDumpError::MalformedHeader {
            line: feed.line_no,
            detail: value,
        })?;
    if count > limits.max_atoms_per_frame {
        return Err(LammpsDumpError::ResourceLimit {
            resource: "atoms per frame",
            actual: count,
            limit: limits.max_atoms_per_frame,
        });
    }
    Ok(count)
}

fn parse_box_header(line: &str, line_no: usize) -> Result<(bool, [String; 3]), LammpsDumpError> {
    let fields = strip_item_prefix(line, "BOX BOUNDS")?
        .split_whitespace()
        .collect::<Vec<_>>();
    if fields.len() == 6 && fields[..3] == ["xy", "xz", "yz"] {
        return Ok((true, [fields[3].into(), fields[4].into(), fields[5].into()]));
    }
    if fields.len() == 3 {
        return Ok((
            false,
            [fields[0].into(), fields[1].into(), fields[2].into()],
        ));
    }
    Err(LammpsDumpError::MalformedHeader {
        line: line_no,
        detail: line.into(),
    })
}

fn parse_box_row(line: &str, line_no: usize, triclinic: bool) -> Result<Vec<f64>, LammpsDumpError> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    let expected = if triclinic { 3 } else { 2 };
    if fields.len() != expected {
        return Err(LammpsDumpError::InvalidBox {
            reason: format!(
                "BOX BOUNDS row at line {line_no} has {} column(s), but the header's \
                 'xy xz yz' flag {} present, which requires {expected} -- \
                 inconsistent, possibly hand-edited file",
                fields.len(),
                if triclinic { "is" } else { "is not" }
            ),
        });
    }
    fields
        .into_iter()
        .map(|field| {
            let value = field
                .parse::<f64>()
                .map_err(|_| LammpsDumpError::MalformedHeader {
                    line: line_no,
                    detail: line.into(),
                })?;
            if value.is_finite() {
                Ok(value)
            } else {
                Err(LammpsDumpError::InvalidBox {
                    reason: format!("non-finite BOX BOUNDS value '{field}' at line {line_no}"),
                })
            }
        })
        .collect()
}

fn read_box<F>(
    feed: &mut LineFeed<LammpsDumpError, F>,
) -> Result<(LammpsBox, [String; 3]), LammpsDumpError>
where
    F: FnMut() -> Result<Option<String>, LammpsDumpError>,
{
    let header = required_line(feed, "ITEM: BOX BOUNDS")?;
    let (triclinic, boundary_flags) = parse_box_header(&header, feed.line_no)?;
    let mut rows = Vec::with_capacity(3);
    for _ in 0..3 {
        let line = required_line(feed, "BOX BOUNDS data row")?;
        rows.push(parse_box_row(&line, feed.line_no, triclinic)?);
    }
    let tilt = triclinic.then(|| [rows[0][2], rows[1][2], rows[2][2]]);
    let bounds = box_bounds_to_true(
        [rows[0][0], rows[1][0], rows[2][0]],
        [rows[0][1], rows[1][1], rows[2][1]],
        tilt,
    );
    bounds
        .validate()
        .map_err(|reason| LammpsDumpError::InvalidBox { reason })?;
    Ok((bounds, boundary_flags))
}

fn read_columns<F>(
    feed: &mut LineFeed<LammpsDumpError, F>,
    limits: &LammpsDumpParseLimits,
    frame_index: usize,
) -> Result<Vec<String>, LammpsDumpError>
where
    F: FnMut() -> Result<Option<String>, LammpsDumpError>,
{
    let line = required_line(feed, "ITEM: ATOMS")?;
    let columns = strip_item_prefix(&line, "ATOMS")?
        .split_whitespace()
        .map(str::to_string)
        .collect::<Vec<_>>();
    if columns.len() > limits.max_columns {
        return Err(LammpsDumpError::ResourceLimit {
            resource: "columns per frame",
            actual: columns.len(),
            limit: limits.max_columns,
        });
    }
    if frame_index >= limits.max_frames {
        return Err(LammpsDumpError::ResourceLimit {
            resource: "frames",
            actual: frame_index + 1,
            limit: limits.max_frames,
        });
    }
    Ok(columns)
}

fn parse_atom_row(
    line: &str,
    line_no: usize,
    row_index: usize,
    columns: &[String],
) -> Result<Vec<f64>, LammpsDumpError> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    if fields.len() != columns.len() {
        return Err(LammpsDumpError::ColumnCountMismatch {
            declared: columns.len(),
            actual: fields.len(),
            line: line_no,
        });
    }
    fields
        .into_iter()
        .enumerate()
        .map(|(column_index, field)| {
            let value = field
                .parse::<f64>()
                .map_err(|_| LammpsDumpError::InvalidNumber {
                    column: columns[column_index].clone(),
                    row: row_index,
                    raw: field.into(),
                })?;
            if value.is_finite() {
                Ok(value)
            } else {
                Err(LammpsDumpError::NonFiniteValue {
                    column: columns[column_index].clone(),
                    row: row_index,
                })
            }
        })
        .collect()
}

fn read_atom_rows<F>(
    feed: &mut LineFeed<LammpsDumpError, F>,
    count: usize,
    columns: &[String],
) -> Result<Vec<Vec<f64>>, LammpsDumpError>
where
    F: FnMut() -> Result<Option<String>, LammpsDumpError>,
{
    let mut rows = Vec::with_capacity(count);
    for row_index in 0..count {
        let line = required_line(feed, "ATOMS data row")?;
        if line.starts_with("ITEM:") {
            return Err(LammpsDumpError::AtomCountMismatch {
                declared: count,
                actual: row_index,
            });
        }
        rows.push(parse_atom_row(&line, feed.line_no, row_index, columns)?);
    }
    Ok(rows)
}

/// Read one frame from `feed`, or `Ok(None)` on a clean EOF *between*
/// frames (i.e. before any of this frame's lines have been read).
fn read_frame_from_feed<F>(
    feed: &mut LineFeed<LammpsDumpError, F>,
    limits: &LammpsDumpParseLimits,
    frame_index: usize,
) -> Result<Option<LammpsDumpFrame>, LammpsDumpError>
where
    F: FnMut() -> Result<Option<String>, LammpsDumpError>,
{
    let Some(timestep) = read_timestep(feed)? else {
        return Ok(None);
    };
    let num_atoms = read_atom_count(feed, limits)?;
    let (box_bounds, boundary_flags) = read_box(feed)?;
    let column_names = read_columns(feed, limits, frame_index)?;
    let rows = read_atom_rows(feed, num_atoms, &column_names)?;

    Ok(Some(LammpsDumpFrame {
        timestep,
        num_atoms,
        box_bounds,
        boundary_flags,
        column_names,
        rows,
    }))
}

// ---------------------------------------------------------------------------
// Whole-string entry point
// ---------------------------------------------------------------------------

/// Parse a single LAMMPS dump frame from `text`.
///
/// If `text` contains more than one frame (multiple `ITEM: TIMESTEP`
/// blocks concatenated, as a trajectory file would), only the first is
/// parsed; anything after it is silently ignored. Use
/// [`LammpsDumpReader`] to read every frame of a multi-frame trajectory.
pub fn parse_lammps_dump_frame(text: &str) -> Result<LammpsDumpFrame, LammpsDumpError> {
    parse_lammps_dump_frame_with_limits(text, &LammpsDumpParseLimits::default())
}

/// Parse a single LAMMPS dump frame with explicit resource limits.
pub fn parse_lammps_dump_frame_with_limits(
    text: &str,
    limits: &LammpsDumpParseLimits,
) -> Result<LammpsDumpFrame, LammpsDumpError> {
    if text.len() > limits.max_input_bytes {
        return Err(LammpsDumpError::ResourceLimit {
            resource: "input bytes",
            actual: text.len(),
            limit: limits.max_input_bytes,
        });
    }
    if let Some(line_bytes) = text.lines().map(str::len).max()
        && line_bytes > limits.max_line_bytes
    {
        return Err(LammpsDumpError::ResourceLimit {
            resource: "line bytes",
            actual: line_bytes,
            limit: limits.max_line_bytes,
        });
    }
    let mut lines = text.lines();
    let next_line =
        move || -> Result<Option<String>, LammpsDumpError> { Ok(lines.next().map(str::to_string)) };
    let mut feed = LineFeed::new(next_line);
    read_frame_from_feed(&mut feed, limits, 0)?.ok_or(LammpsDumpError::TruncatedInput {
        context: "frame (empty input)",
    })
}

// ---------------------------------------------------------------------------
// Streaming reader
// ---------------------------------------------------------------------------

/// Streaming-*input* reader over a LAMMPS dump/trajectory file, yielding
/// one [`LammpsDumpFrame`] at a time without loading the whole trajectory
/// into memory at once. Follows this crate's `SdfFileReader`
/// (`crate::sdf`) precedent -- `Iterator<Item = Result<LammpsDumpFrame,
/// LammpsDumpError>>` over any [`std::io::BufRead`] source -- since a
/// trajectory, like SDF, is inherently multi-record (unlike
/// `crate::cube::CubeFileReader`'s single-grid `.read()`).
pub struct LammpsDumpReader<R: std::io::BufRead> {
    reader: R,
    done: bool,
    limits: LammpsDumpParseLimits,
    bytes_read: usize,
    frames_read: usize,
}

impl<R: std::io::BufRead> LammpsDumpReader<R> {
    /// Wrap any `BufRead` source (e.g. `BufReader<File>`).
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            done: false,
            limits: LammpsDumpParseLimits::default(),
            bytes_read: 0,
            frames_read: 0,
        }
    }

    /// Wrap a source with explicit trajectory resource limits.
    pub fn with_limits(reader: R, limits: LammpsDumpParseLimits) -> Self {
        Self {
            reader,
            done: false,
            limits,
            bytes_read: 0,
            frames_read: 0,
        }
    }
}

impl<R: std::io::BufRead> Iterator for LammpsDumpReader<R> {
    type Item = Result<LammpsDumpFrame, LammpsDumpError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let reader = &mut self.reader;
        let limits = self.limits;
        let bytes_read = &mut self.bytes_read;
        let next_line = move || -> Result<Option<String>, LammpsDumpError> {
            let mut buf = String::new();
            match reader.read_line(&mut buf) {
                Ok(0) => Ok(None),
                Ok(n) => {
                    *bytes_read = bytes_read.saturating_add(n);
                    if *bytes_read > limits.max_input_bytes {
                        return Err(LammpsDumpError::ResourceLimit {
                            resource: "input bytes",
                            actual: *bytes_read,
                            limit: limits.max_input_bytes,
                        });
                    }
                    if buf.trim_end_matches(['\n', '\r']).len() > limits.max_line_bytes {
                        return Err(LammpsDumpError::ResourceLimit {
                            resource: "line bytes",
                            actual: buf.trim_end_matches(['\n', '\r']).len(),
                            limit: limits.max_line_bytes,
                        });
                    }
                    Ok(Some(buf.trim_end_matches(['\n', '\r']).to_string()))
                }
                Err(e) => Err(LammpsDumpError::Io(e.to_string())),
            }
        };
        let mut feed = LineFeed::new(next_line);
        match read_frame_from_feed(&mut feed, &limits, self.frames_read) {
            Ok(Some(frame)) => {
                self.frames_read += 1;
                Some(Ok(frame))
            }
            Ok(None) => {
                self.done = true;
                None
            }
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Write one [`LammpsDumpFrame`] as LAMMPS dump text.
pub fn write_lammps_dump_frame(frame: &LammpsDumpFrame) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    let _ = writeln!(out, "ITEM: TIMESTEP");
    let _ = writeln!(out, "{}", frame.timestep);
    let _ = writeln!(out, "ITEM: NUMBER OF ATOMS");
    let _ = writeln!(out, "{}", frame.num_atoms);

    let (bound_lo, bound_hi) = true_to_box_bounds(&frame.box_bounds);
    let [f0, f1, f2] = &frame.boundary_flags;
    match frame.box_bounds.tilt {
        Some([xy, xz, yz]) => {
            let _ = writeln!(out, "ITEM: BOX BOUNDS xy xz yz {f0} {f1} {f2}");
            let _ = writeln!(out, "{} {} {}", bound_lo[0], bound_hi[0], xy);
            let _ = writeln!(out, "{} {} {}", bound_lo[1], bound_hi[1], xz);
            let _ = writeln!(out, "{} {} {}", bound_lo[2], bound_hi[2], yz);
        }
        None => {
            let _ = writeln!(out, "ITEM: BOX BOUNDS {f0} {f1} {f2}");
            let _ = writeln!(out, "{} {}", bound_lo[0], bound_hi[0]);
            let _ = writeln!(out, "{} {}", bound_lo[1], bound_hi[1]);
            let _ = writeln!(out, "{} {}", bound_lo[2], bound_hi[2]);
        }
    }

    let _ = writeln!(out, "ITEM: ATOMS {}", frame.column_names.join(" "));
    for row in &frame.rows {
        let strs: Vec<String> = row.iter().map(f64::to_string).collect();
        let _ = writeln!(out, "{}", strs.join(" "));
    }
    out
}

/// Write multiple frames as one trajectory file (plain concatenation --
/// see module docs, a trajectory is just N frames back to back). Only the
/// *reader* streams; a caller assembling a `&[LammpsDumpFrame]` to pass
/// here already holds every frame in memory.
pub fn write_lammps_trajectory(frames: &[LammpsDumpFrame]) -> String {
    frames.iter().map(write_lammps_dump_frame).collect()
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn orthogonal_frame() -> LammpsDumpFrame {
        LammpsDumpFrame {
            timestep: 1000,
            num_atoms: 2,
            box_bounds: LammpsBox {
                lo: [0.0, 0.0, 0.0],
                hi: [10.0, 20.0, 30.0],
                tilt: None,
            },
            boundary_flags: ["pp".to_string(), "pp".to_string(), "pp".to_string()],
            column_names: vec![
                "id".to_string(),
                "type".to_string(),
                "x".to_string(),
                "y".to_string(),
                "z".to_string(),
            ],
            rows: vec![vec![1.0, 1.0, 1.0, 2.0, 3.0], vec![2.0, 1.0, 4.0, 5.0, 6.0]],
        }
    }

    fn triclinic_frame() -> LammpsDumpFrame {
        LammpsDumpFrame {
            timestep: 2000,
            num_atoms: 1,
            box_bounds: LammpsBox {
                lo: [0.0, 0.0, 0.0],
                hi: [10.0, 10.0, 10.0],
                tilt: Some([2.0, 1.0, 0.5]),
            },
            boundary_flags: ["pp".to_string(), "ff".to_string(), "ss".to_string()],
            column_names: vec![
                "id".to_string(),
                "xs".to_string(),
                "ys".to_string(),
                "zs".to_string(),
            ],
            rows: vec![vec![1.0, 0.5, 0.5, 0.5]],
        }
    }

    #[test]
    fn round_trip_orthogonal_frame() {
        let frame = orthogonal_frame();
        let text = write_lammps_dump_frame(&frame);
        let parsed = parse_lammps_dump_frame(&text).unwrap();
        assert_eq!(parsed, frame);
    }

    #[test]
    fn round_trip_triclinic_frame() {
        let frame = triclinic_frame();
        let text = write_lammps_dump_frame(&frame);
        let parsed = parse_lammps_dump_frame(&text).unwrap();
        assert_eq!(parsed, frame);
    }

    #[test]
    fn triclinic_bound_true_bound_identity_over_several_tilts() {
        // Cheap inverse-consistency check: bound -> true -> bound is the
        // identity for *any* tilt, including this formula's own sign
        // convention -- but see the hand-computed test below for why this
        // alone cannot prove the formula itself is correct (a
        // self-consistently-wrong MIN/MAX definition would also pass
        // this).
        let tilts = [[3.0, 4.0, 2.0], [-3.0, -4.0, -2.0], [3.0, 4.0, -2.0]];
        for tilt in tilts {
            let bound_lo = [1.0, 2.0, 3.0];
            let bound_hi = [11.0, 22.0, 33.0];
            let true_box = box_bounds_to_true(bound_lo, bound_hi, Some(tilt));
            let (round_lo, round_hi) = true_to_box_bounds(&true_box);
            for i in 0..3 {
                assert!(
                    (round_lo[i] - bound_lo[i]).abs() < 1e-9,
                    "lo[{i}] tilt={tilt:?}"
                );
                assert!(
                    (round_hi[i] - bound_hi[i]).abs() < 1e-9,
                    "hi[{i}] tilt={tilt:?}"
                );
            }
        }
    }

    #[test]
    fn triclinic_bound_conversion_matches_hand_computed_absolute_values() {
        // The identity test above cannot catch a formula that's
        // internally self-consistent but wrong (e.g. MAX(0,xy,xz)
        // instead of MAX(0,xy,xz,xy+xz)) -- both directions would use the
        // same wrong extremum and still round-trip. This test instead
        // asserts true_to_box_bounds against hand-computed absolute
        // values, picking tilts where xy+xz (not just xy or xz alone) is
        // the true extremum.
        //
        // true box: lo=[0,0,0] hi=[10,10,10], tilt=[xy=3, xz=4, yz=2]
        //   x candidates {0,3,4,7}: min=0, max=7 -> xlo_bound=0+0=0, xhi_bound=10+7=17
        //   y candidates {0,2}:     min=0, max=2 -> ylo_bound=0+0=0, yhi_bound=10+2=12
        //   z: no tilt term         -> zlo_bound=0, zhi_bound=10
        let b1 = LammpsBox {
            lo: [0.0, 0.0, 0.0],
            hi: [10.0, 10.0, 10.0],
            tilt: Some([3.0, 4.0, 2.0]),
        };
        let (lo1, hi1) = true_to_box_bounds(&b1);
        assert_eq!(lo1, [0.0, 0.0, 0.0]);
        assert_eq!(hi1, [17.0, 12.0, 10.0]);

        // tilt=[-3,-4,-2]: x candidates {0,-3,-4,-7}: min=-7,max=0 -> xlo_bound=-7, xhi_bound=10
        //                  y candidates {0,-2}: min=-2,max=0 -> ylo_bound=-2, yhi_bound=10
        let b2 = LammpsBox {
            lo: [0.0, 0.0, 0.0],
            hi: [10.0, 10.0, 10.0],
            tilt: Some([-3.0, -4.0, -2.0]),
        };
        let (lo2, hi2) = true_to_box_bounds(&b2);
        assert_eq!(lo2, [-7.0, -2.0, 0.0]);
        assert_eq!(hi2, [10.0, 10.0, 10.0]);

        // tilt=[3,4,-2]: x same as b1 (min=0,max=7) -> xlo_bound=0, xhi_bound=17
        //                y candidates {0,-2}: min=-2,max=0 -> ylo_bound=-2, yhi_bound=10
        let b3 = LammpsBox {
            lo: [0.0, 0.0, 0.0],
            hi: [10.0, 10.0, 10.0],
            tilt: Some([3.0, 4.0, -2.0]),
        };
        let (lo3, hi3) = true_to_box_bounds(&b3);
        assert_eq!(lo3, [0.0, -2.0, 0.0]);
        assert_eq!(hi3, [17.0, 10.0, 10.0]);

        // And the read-direction inverse recovers the original true box
        // from each hand-computed bound pair.
        for (b, lo, hi) in [(b1, lo1, hi1), (b2, lo2, hi2), (b3, lo3, hi3)] {
            let recovered = box_bounds_to_true(lo, hi, b.tilt);
            assert_eq!(recovered, b);
        }
    }

    #[test]
    fn multi_frame_streaming_round_trip() {
        use std::io::{BufReader, Cursor};
        let mut f1 = orthogonal_frame();
        f1.timestep = 0;
        let mut f2 = orthogonal_frame();
        f2.timestep = 100;
        let mut f3 = triclinic_frame();
        f3.timestep = 200;
        let frames = vec![f1, f2, f3];

        let text = write_lammps_trajectory(&frames);
        let cursor = Cursor::new(text.into_bytes());
        let read_back: Result<Vec<LammpsDumpFrame>, LammpsDumpError> =
            LammpsDumpReader::new(BufReader::new(cursor)).collect();
        assert_eq!(read_back.unwrap(), frames);
    }

    #[test]
    fn scaled_coordinates_orthogonal_hand_computed() {
        // box lo=[0,0,0] hi=[10,20,30], xs=0.5 ys=0.25 zs=0.1
        // x = 0 + 0.5*10 = 5.0
        // y = 0 + 0.25*20 = 5.0
        // z = 0 + 0.1*30 = 3.0
        let frame = LammpsDumpFrame {
            timestep: 0,
            num_atoms: 1,
            box_bounds: LammpsBox {
                lo: [0.0, 0.0, 0.0],
                hi: [10.0, 20.0, 30.0],
                tilt: None,
            },
            boundary_flags: ["pp".to_string(), "pp".to_string(), "pp".to_string()],
            column_names: vec![
                "id".to_string(),
                "xs".to_string(),
                "ys".to_string(),
                "zs".to_string(),
            ],
            rows: vec![vec![1.0, 0.5, 0.25, 0.1]],
        };
        let pos = frame.cartesian_positions().unwrap();
        assert!((pos[0][0] - 5.0).abs() < 1e-9);
        assert!((pos[0][1] - 5.0).abs() < 1e-9);
        assert!((pos[0][2] - 3.0).abs() < 1e-9);
    }

    #[test]
    fn scaled_coordinates_triclinic_hand_computed() {
        // box lo=[0,0,0] hi=[10,10,10], tilt xy=2 xz=1 yz=0.5, xs=ys=zs=0.5
        // x = 0 + 0.5*10 + 0.5*2 + 0.5*1 = 5 + 1 + 0.5 = 6.5
        // y = 0 + 0.5*10 + 0.5*0.5 = 5 + 0.25 = 5.25
        // z = 0 + 0.5*10 = 5.0
        let frame = triclinic_frame();
        let pos = frame.cartesian_positions().unwrap();
        assert!((pos[0][0] - 6.5).abs() < 1e-9);
        assert!((pos[0][1] - 5.25).abs() < 1e-9);
        assert!((pos[0][2] - 5.0).abs() < 1e-9);
    }

    #[test]
    fn cartesian_positions_passes_through_x_y_z() {
        let frame = orthogonal_frame();
        let pos = frame.cartesian_positions().unwrap();
        assert_eq!(pos, vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
    }

    #[test]
    fn cartesian_positions_none_for_unwrapped_only_columns() {
        let mut frame = orthogonal_frame();
        frame.column_names = vec![
            "id".to_string(),
            "xu".to_string(),
            "yu".to_string(),
            "zu".to_string(),
        ];
        frame.rows = vec![vec![1.0, 100.0, 200.0, 300.0], vec![2.0, 1.0, 2.0, 3.0]];
        assert_eq!(frame.cartesian_positions(), None);
        // But the raw column is still directly accessible.
        assert_eq!(frame.column("xu"), Some(vec![100.0, 1.0]));
    }

    #[test]
    fn custom_column_names_round_trip() {
        let mut frame = orthogonal_frame();
        frame.column_names = vec![
            "id".to_string(),
            "c_myCompute[1]".to_string(),
            "f_myFix".to_string(),
        ];
        frame.rows = vec![vec![1.0, 0.123, -4.5], vec![2.0, 0.456, 7.8]];
        let text = write_lammps_dump_frame(&frame);
        let parsed = parse_lammps_dump_frame(&text).unwrap();
        assert_eq!(parsed.column_names, frame.column_names);
        assert_eq!(parsed.rows, frame.rows);
    }

    #[test]
    fn zero_atom_frame_does_not_panic_or_error() {
        let mut frame = orthogonal_frame();
        frame.num_atoms = 0;
        frame.rows = Vec::new();
        let text = write_lammps_dump_frame(&frame);
        let parsed = parse_lammps_dump_frame(&text).unwrap();
        assert_eq!(parsed.num_atoms, 0);
        assert!(parsed.rows.is_empty());
    }

    #[test]
    fn unexpected_item_is_typed_error() {
        let text = "ITEM: NOT TIMESTEP\n0\n";
        let err = parse_lammps_dump_frame(text).unwrap_err();
        assert!(matches!(err, LammpsDumpError::UnexpectedItem { .. }));
    }

    #[test]
    fn column_count_mismatch_is_typed_error() {
        let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n\
                     0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id type x y z\n1 1 0.0 0.0\n";
        let err = parse_lammps_dump_frame(text).unwrap_err();
        assert_eq!(
            err,
            LammpsDumpError::ColumnCountMismatch {
                declared: 5,
                actual: 4,
                line: 10
            }
        );
    }

    #[test]
    fn non_finite_value_is_typed_error() {
        let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n\
                     0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id x\n1 NaN\n";
        let err = parse_lammps_dump_frame(text).unwrap_err();
        assert_eq!(
            err,
            LammpsDumpError::NonFiniteValue {
                column: "x".to_string(),
                row: 0
            }
        );
    }

    #[test]
    fn atom_count_mismatch_is_typed_error_not_corruption() {
        let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\nITEM: BOX BOUNDS pp pp pp\n\
                     0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id x\n1 0.0\n2 1.0\nITEM: TIMESTEP\n1\n";
        let err = parse_lammps_dump_frame(text).unwrap_err();
        assert_eq!(
            err,
            LammpsDumpError::AtomCountMismatch {
                declared: 3,
                actual: 2
            }
        );
    }

    #[test]
    fn box_bounds_column_count_disagrees_with_header_flag_is_typed_error() {
        // Header says orthogonal (no "xy xz yz"), but a data row has 3
        // columns -- inconsistent, must fail rather than silently
        // guessing which signal to trust.
        let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n\
                     0.0 10.0 1.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id x\n1 0.0\n";
        let err = parse_lammps_dump_frame(text).unwrap_err();
        assert!(matches!(err, LammpsDumpError::InvalidBox { .. }));
    }

    #[test]
    fn truncated_input_is_typed_error_not_panic() {
        let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n";
        let err = parse_lammps_dump_frame(text).unwrap_err();
        assert!(matches!(err, LammpsDumpError::TruncatedInput { .. }));
    }

    #[test]
    fn empty_input_is_truncated_input_not_panic() {
        let err = parse_lammps_dump_frame("").unwrap_err();
        assert!(matches!(err, LammpsDumpError::TruncatedInput { .. }));
    }

    #[test]
    fn bounded_frame_parser_rejects_input_and_line_limits() {
        let text = write_lammps_dump_frame(&orthogonal_frame());
        assert!(matches!(
            parse_lammps_dump_frame_with_limits(
                &text,
                &LammpsDumpParseLimits {
                    max_input_bytes: 8,
                    ..Default::default()
                }
            ),
            Err(LammpsDumpError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));

        let long_line = format!("ITEM: TIMESTEP{}\n", "x".repeat(32));
        assert!(matches!(
            parse_lammps_dump_frame_with_limits(
                &long_line,
                &LammpsDumpParseLimits {
                    max_line_bytes: 16,
                    ..Default::default()
                }
            ),
            Err(LammpsDumpError::ResourceLimit {
                resource: "line bytes",
                ..
            })
        ));
    }

    #[test]
    fn bounded_frame_parser_rejects_atom_and_column_limits() {
        let text = write_lammps_dump_frame(&orthogonal_frame());
        assert!(matches!(
            parse_lammps_dump_frame_with_limits(
                &text,
                &LammpsDumpParseLimits {
                    max_atoms_per_frame: 1,
                    ..Default::default()
                }
            ),
            Err(LammpsDumpError::ResourceLimit {
                resource: "atoms per frame",
                ..
            })
        ));
        assert!(matches!(
            parse_lammps_dump_frame_with_limits(
                &text,
                &LammpsDumpParseLimits {
                    max_columns: 2,
                    ..Default::default()
                }
            ),
            Err(LammpsDumpError::ResourceLimit {
                resource: "columns per frame",
                ..
            })
        ));
    }

    #[test]
    fn bounded_reader_rejects_a_frame_beyond_frame_budget() {
        use std::io::Cursor;

        let text = write_lammps_trajectory(&[orthogonal_frame(), orthogonal_frame()]);
        let mut reader = LammpsDumpReader::with_limits(
            Cursor::new(text),
            LammpsDumpParseLimits {
                max_frames: 1,
                ..Default::default()
            },
        );
        assert!(reader.next().unwrap().is_ok());
        assert!(matches!(
            reader.next().unwrap(),
            Err(LammpsDumpError::ResourceLimit {
                resource: "frames",
                ..
            })
        ));
    }
}
