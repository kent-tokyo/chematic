//! Gaussian Cube volumetric scalar-field format (`.cube`/`.cub`) reader and
//! writer.
//!
//! ## Sources and confidence
//!
//! Implemented independently from public documentation of the format. No
//! source code, comments, or tables were copied from any other tool (VMD,
//! Open Babel, ASE, Multiwfn, ...); where a real independent tool's
//! documented/observed behavior is cited below as corroboration, it was
//! read only as a behavioral oracle, never as a source to copy from.
//!
//! **High confidence** -- corroborated by multiple independent primary
//! sources that agree exactly:
//! - 2 leading comment/title lines, contents otherwise unconstrained
//!   (<https://paulbourke.net/dataformats/cube/>, a long-standing
//!   independent community reference; corroborated by
//!   <https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html>, a
//!   modern formal grammar written to be RFC-like, and by
//!   <https://gaussian.com/cubegen/>, Gaussian Inc.'s own `cubegen` tool
//!   documentation.)
//! - Line 3: `NAtoms X-Origin Y-Origin Z-Origin [NVal]` -- `NVal` is an
//!   optional 5th token (h5cube-spec; gaussian.com's own field list is
//!   `"NAtoms, X-Origin, Y-Origin, Z-Origin, NVal"`).
//! - 3 axis lines, each `N_axis Vx Vy Vz` (voxel count along that axis,
//!   then its full 3-vector step -- not just a scalar spacing, so
//!   non-orthogonal grids are representable).
//! - Atom lines: `AtomicNumber Charge X Y Z` -- gaussian.com states this
//!   field order explicitly (`"IA1, Chg1, X1, Y1, Z1"`); h5cube-spec
//!   clarifies `Charge` is the atom's **effective nuclear charge, which
//!   deviates from the atomic number when the calculation used an
//!   effective core potential (ECP)** -- i.e. it is *not* a partial/ionic
//!   charge, a common misreading of this field.
//! - Multiple values per voxel: a **negative** `NAtoms` requires a
//!   dataset-identifier line after the atom lines (`m` positive, followed
//!   by `m` arbitrary integer IDs -- typically MO indices) and *forbids*
//!   `NVal` on line 3 (h5cube-spec); gaussian.com corroborates negative
//!   `NAtoms` meaning "molecular orbital output ... an additional record
//!   follows the data for the final atom". A **positive** `NAtoms` with
//!   `NVal != 1` is the *other*, simpler multi-dataset form (no
//!   identifier line). This module does not support either: see
//!   [`CubeError::MultiDatasetUnsupported`].
//! - Voxel data block ordering: first axis outermost, third axis
//!   innermost (`"the x axis as the outer loop and the z axis as the
//!   inner loop"`, paulbourke.net; h5cube-spec's nested-loop pseudocode
//!   agrees). Values are whitespace/newline-separated; this reader does
//!   **not** hard-code "6 values per line" as a requirement (h5cube-spec
//!   describes 6-per-line as the conventional wrapping, not a mandated
//!   one) -- it reads tokens until the expected total count is reached,
//!   tolerant of any line-wrapping.
//!
//! **Medium confidence** -- a genuine, unresolved conflict between sources
//! on what a *negative voxel count on the first axis line* means, flagged
//! explicitly per this project's stop-and-report-rather-than-guess
//! discipline rather than silently picking one and asserting it:
//! - paulbourke.net and h5cube-spec both describe it as a **file-content**
//!   units flag: negative first-axis count => the file's origin/axis
//!   vectors (and, by this module's own extension -- see below -- atom
//!   positions) are in **Ångström**; positive => the format's native
//!   **Bohr** (atomic units).
//! - gaussian.com's `cubegen` documentation contains a sentence that reads
//!   the *opposite* direction ("If N1<0 the input cube coordinates are
//!   assumed to be in Bohr, otherwise ... Angstroms") -- but on fetching
//!   surrounding context directly, this sentence is in `cubegen`'s
//!   **custom-grid input-specification syntax** section (how a user types
//!   a grid request *to* `cubegen`), not a description of the *output*
//!   `.cube` file's own header convention. It is not in conflict with the
//!   file-format claim above once read in context; it is answering a
//!   different question.
//! - A real-world resolution, corroborated by a Jmol bug-tracker thread
//!   (<https://sourceforge.net/p/jmol/bugs/370/>) in which Jmol's
//!   maintainer quotes Gaussian's own documentation: **"All values in the
//!   cube file are in atomic units, regardless of the input units"** --
//!   i.e. a genuine Gaussian-authored `.cube` file is always Bohr on disk,
//!   and the sign is largely vestigial for such files. Consistent with
//!   this, ASE's cube reader (`ase.io.cube`, read only as a behavioral
//!   oracle, not copied) does not implement the sign-flag convention at
//!   all: it always takes `abs()` of the voxel count and always assumes
//!   Bohr.
//!
//!   This module implements the **defensive/compatibility reading**
//!   (negative first-axis count => Ångström), matching paulbourke.net and
//!   h5cube-spec, because real non-Gaussian tools in the wild are
//!   documented to have produced cube-like files using that convention
//!   (the cited Jmol bug report exists precisely because such a file was
//!   encountered) -- silently ignoring the sign risks silently
//!   misinterpreting those files' units by a factor of ~1.89
//!   (1 Bohr = 0.529177 Å). Rather than silently normalizing away this
//!   ambiguity, [`VolumetricGrid::units`] records exactly which unit was
//!   read (or, for [`write_cube`], which unit to write), so the file's
//!   own numbers are stored/emitted verbatim and never silently rescaled.
//!   No source states whether atom-line coordinates follow the same
//!   sign-flag convention as origin/axis vectors (only "the input cube
//!   coordinates" is discussed, ambiguous between "just the grid" and
//!   "the whole file") -- this module applies the flag to atom positions
//!   too, for internal single-unit-per-file consistency; **this specific
//!   extension is this module's own assumption, not a quoted source
//!   convention.**
//! - h5cube-spec additionally states NY and NZ (axes 2 and 3) "MUST always
//!   be positive" -- only the first axis line's sign is meaningful; this
//!   module rejects a non-positive count on axis 2 or 3 as a typed error
//!   rather than silently taking its absolute value.
//!
//! ## Multiple datasets and CHGCAR/LOCPOT: explicitly out of scope
//!
//! Multi-dataset Cube files (see above) are rejected with a typed error,
//! never silently truncated to the first dataset. VASP's CHGCAR/LOCPOT
//! volumetric formats are a separate, later roadmap step -- not
//! implemented here; both would build on the same [`VolumetricGrid`] type
//! once that step starts.

use crate::volumetric::{GridAtom, GridError, GridUnits, LineFeed, VolumetricGrid};
use chematic_core::Element;

// ---------------------------------------------------------------------------
// Limits
// ---------------------------------------------------------------------------

/// Security/robustness limits enforced before any chemistry or large
/// allocation happens. Shared by [`parse_cube_with_limits`] and
/// [`CubeFileReader`].
#[derive(Debug, Clone, Copy)]
pub struct CubeParseLimits {
    /// Cumulative byte budget across the whole input (checked incrementally
    /// for [`CubeFileReader`], upfront for [`parse_cube_with_limits`]).
    pub max_input_bytes: usize,
    pub max_atoms: usize,
    /// Cap on `shape[0] * shape[1] * shape[2]`, checked from the header
    /// *before* the (potentially huge) voxel data block is read.
    pub max_grid_points: usize,
}

impl Default for CubeParseLimits {
    fn default() -> Self {
        Self {
            // Real cube files for large grids are routinely multi-hundred-MB;
            // callers that need more (or that target a 32-bit `usize`
            // platform, where this constant is architecture-bounded) can
            // widen this via `parse_with_limits`/`CubeFileReader::with_limits`.
            max_input_bytes: 1 << 30, // 1 GiB
            max_atoms: 200_000,
            max_grid_points: 100_000_000, // ~800 MB of f64, far below a pathological 10^9-per-axis header
        }
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors from [`parse_cube`]/[`parse_cube_with_limits`]/[`CubeFileReader`]/
/// [`write_cube`].
#[derive(Debug, Clone, PartialEq)]
pub enum CubeError {
    /// Input exceeded [`CubeParseLimits::max_input_bytes`].
    InputTooLarge { limit: usize },
    /// An IO error occurred while reading from a streaming
    /// [`CubeFileReader`] source.
    Io(String),
    /// The input ended before a required section was fully read.
    UnexpectedEnd { context: &'static str },
    /// Line 3 (`NAtoms X0 Y0 Z0 [NVal]`) was missing, or did not have 4 or
    /// 5 whitespace-delimited tokens.
    InvalidCountLine { line: usize, detail: String },
    /// One of the 3 axis lines (`N Vx Vy Vz`) was missing or malformed.
    InvalidAxisLine {
        line: usize,
        axis: usize,
        detail: String,
    },
    /// Axis 2 or 3's voxel count was not strictly positive (only axis 1's
    /// sign carries the Bohr/Ångström convention -- see module docs).
    NonPositiveAxisCount { axis: usize, value: i64 },
    /// An atom line did not have exactly 5 whitespace-delimited tokens.
    InvalidAtomLine { line: usize, detail: String },
    /// An atom line's atomic number is not a recognized element (1-118).
    /// Note: real Cube files from counterpoise/ghost-atom calculations can
    /// legitimately use atomic number 0 for a dummy center with no basis
    /// functions -- this parser does not special-case that and rejects it
    /// the same as any other unrecognized value.
    UnknownAtomicNumber { line: usize, value: i64 },
    /// A numeric token could not be parsed as a float/int at all.
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
    /// Atom count exceeded [`CubeParseLimits::max_atoms`].
    TooManyAtoms { limit: usize },
    /// `NAtoms` was negative (dataset-identifier-list form), or line 3's
    /// optional `NVal` field was present and not `1` (positive-`NAtoms`
    /// form) -- both are real, documented multiple-values-per-voxel
    /// conventions. See module docs; not supported.
    MultiDatasetUnsupported {
        natoms_field: i64,
        nval: Option<i64>,
    },
    /// The voxel data block had fewer whitespace-delimited numeric tokens
    /// than `shape[0]*shape[1]*shape[2]` requires.
    ValueCountMismatch { expected: usize, found: usize },
    /// The voxel data block had at least one more numeric token than
    /// `shape[0]*shape[1]*shape[2]` requires. Cube has no defined trailing
    /// footer, so any leftover token here is malformed input.
    TrailingData { after_values: usize },
    /// [`write_cube`] found a `shape` dimension too large to represent as
    /// the signed `i64` an axis line's voxel-count field requires --
    /// raised instead of silently wrapping via an unchecked cast.
    DimensionOutOfRange { axis: usize, value: usize },
    /// A [`VolumetricGrid`] structural-invariant check failed (shape
    /// overflow/cap, value-count mismatch, or a non-finite field) --
    /// raised by [`write_cube`], since a hand-built `VolumetricGrid` is not
    /// guaranteed to satisfy these.
    Grid(GridError),
}

impl std::fmt::Display for CubeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InputTooLarge { limit } => write!(f, "cube input exceeds {limit}-byte limit"),
            Self::Io(msg) => write!(f, "IO error: {msg}"),
            Self::UnexpectedEnd { context } => {
                write!(f, "unexpected end of input while reading {context}")
            }
            Self::InvalidCountLine { line, detail } => {
                write!(
                    f,
                    "invalid atom-count/origin line at line {line}: '{detail}'"
                )
            }
            Self::InvalidAxisLine { line, axis, detail } => {
                write!(f, "invalid axis {axis} line at line {line}: '{detail}'")
            }
            Self::NonPositiveAxisCount { axis, value } => {
                write!(f, "axis {axis} voxel count {value} must be positive")
            }
            Self::InvalidAtomLine { line, detail } => {
                write!(f, "invalid atom line at line {line}: '{detail}'")
            }
            Self::UnknownAtomicNumber { line, value } => {
                write!(f, "unrecognized atomic number {value} at line {line}")
            }
            Self::InvalidNumber { line, context, raw } => {
                write!(f, "invalid {context} value '{raw}' at line {line}")
            }
            Self::NonFiniteValue { line, context, raw } => write!(
                f,
                "{context} value '{raw}' at line {line} is not finite (NaN/Infinite)"
            ),
            Self::TooManyAtoms { limit } => {
                write!(f, "cube atom count exceeds {limit}-atom limit")
            }
            Self::MultiDatasetUnsupported { natoms_field, nval } => write!(
                f,
                "multi-dataset cube file not supported (NAtoms={natoms_field}, NVal={nval:?})"
            ),
            Self::ValueCountMismatch { expected, found } => {
                write!(f, "cube grid expects {expected} values, only found {found}")
            }
            Self::TrailingData { after_values } => write!(
                f,
                "cube data block has extra values after the expected {after_values}"
            ),
            Self::DimensionOutOfRange { axis, value } => write!(
                f,
                "shape axis {axis}'s voxel count {value} exceeds i64::MAX and cannot be written as a Cube axis-line field"
            ),
            Self::Grid(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for CubeError {}

impl From<GridError> for CubeError {
    fn from(e: GridError) -> Self {
        CubeError::Grid(e)
    }
}

// ---------------------------------------------------------------------------
// Shared parse core
// ---------------------------------------------------------------------------

fn parse_finite(raw: &str, line: usize, context: &'static str) -> Result<f64, CubeError> {
    let v: f64 = raw.parse().map_err(|_| CubeError::InvalidNumber {
        line,
        context,
        raw: raw.to_string(),
    })?;
    if !v.is_finite() {
        return Err(CubeError::NonFiniteValue {
            line,
            context,
            raw: raw.to_string(),
        });
    }
    Ok(v)
}

fn parse_cube_from_feed<F>(
    next_line: F,
    limits: &CubeParseLimits,
) -> Result<VolumetricGrid, CubeError>
where
    F: FnMut() -> Result<Option<String>, CubeError>,
{
    let mut feed = LineFeed::new(next_line);

    // Lines 1-2: comment/title lines, contents unconstrained and not preserved.
    feed.line()?.ok_or(CubeError::UnexpectedEnd {
        context: "comment line 1",
    })?;
    feed.line()?.ok_or(CubeError::UnexpectedEnd {
        context: "comment line 2",
    })?;

    // Line 3: NAtoms X0 Y0 Z0 [NVal]
    let count_line = feed.line()?.ok_or(CubeError::UnexpectedEnd {
        context: "atom count / origin line",
    })?;
    let count_line_no = feed.line_no;
    let toks: Vec<&str> = count_line.split_whitespace().collect();
    if toks.len() != 4 && toks.len() != 5 {
        return Err(CubeError::InvalidCountLine {
            line: count_line_no,
            detail: count_line,
        });
    }
    let natoms: i64 = toks[0].parse().map_err(|_| CubeError::InvalidCountLine {
        line: count_line_no,
        detail: count_line.clone(),
    })?;
    let mut origin = [0.0f64; 3];
    for (c, slot) in origin.iter_mut().enumerate() {
        *slot = parse_finite(toks[1 + c], count_line_no, "origin")?;
    }
    let nval: Option<i64> = if toks.len() == 5 {
        Some(toks[4].parse().map_err(|_| CubeError::InvalidCountLine {
            line: count_line_no,
            detail: count_line.clone(),
        })?)
    } else {
        None
    };
    if natoms < 0 || matches!(nval, Some(n) if n != 1) {
        return Err(CubeError::MultiDatasetUnsupported {
            natoms_field: natoms,
            nval,
        });
    }
    let natoms_usize = usize::try_from(natoms).unwrap_or(usize::MAX);
    if natoms_usize > limits.max_atoms {
        return Err(CubeError::TooManyAtoms {
            limit: limits.max_atoms,
        });
    }

    // Lines 4-6: axis lines.
    let mut shape = [0usize; 3];
    let mut axes = [[0.0f64; 3]; 3];
    let mut units = GridUnits::Bohr;
    for (axis, shape_slot) in shape.iter_mut().enumerate() {
        let line = feed.line()?.ok_or(CubeError::UnexpectedEnd {
            context: "axis line",
        })?;
        let line_no = feed.line_no;
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() != 4 {
            return Err(CubeError::InvalidAxisLine {
                line: line_no,
                axis,
                detail: line,
            });
        }
        let n: i64 = toks[0].parse().map_err(|_| CubeError::InvalidAxisLine {
            line: line_no,
            axis,
            detail: line.clone(),
        })?;
        if axis == 0 {
            if n == 0 {
                return Err(CubeError::NonPositiveAxisCount { axis, value: n });
            }
            units = if n < 0 {
                GridUnits::Angstrom
            } else {
                GridUnits::Bohr
            };
        } else if n <= 0 {
            return Err(CubeError::NonPositiveAxisCount { axis, value: n });
        }
        *shape_slot = usize::try_from(n.unsigned_abs()).unwrap_or(usize::MAX);
        for (c, slot) in axes[axis].iter_mut().enumerate() {
            *slot = parse_finite(toks[1 + c], line_no, "axis vector")?;
        }
    }

    let expected_points = shape[0]
        .checked_mul(shape[1])
        .and_then(|v| v.checked_mul(shape[2]))
        .ok_or(GridError::ShapeOverflow { shape })?;
    if expected_points > limits.max_grid_points {
        return Err(GridError::GridTooLarge {
            points: expected_points,
            limit: limits.max_grid_points,
        }
        .into());
    }

    // `|NAtoms|` atom lines.
    let mut atoms = Vec::with_capacity(natoms_usize.min(1_000_000));
    for _ in 0..natoms_usize {
        let line = feed.line()?.ok_or(CubeError::UnexpectedEnd {
            context: "atom line",
        })?;
        let line_no = feed.line_no;
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() != 5 {
            return Err(CubeError::InvalidAtomLine {
                line: line_no,
                detail: line,
            });
        }
        let atomic_number: i64 = toks[0].parse().map_err(|_| CubeError::InvalidAtomLine {
            line: line_no,
            detail: line.clone(),
        })?;
        let element = u8::try_from(atomic_number)
            .ok()
            .and_then(Element::from_atomic_number)
            .ok_or(CubeError::UnknownAtomicNumber {
                line: line_no,
                value: atomic_number,
            })?;
        let charge = parse_finite(toks[1], line_no, "atom charge")?;
        let mut position = [0.0f64; 3];
        for (c, slot) in position.iter_mut().enumerate() {
            *slot = parse_finite(toks[2 + c], line_no, "atom position")?;
        }
        atoms.push(GridAtom {
            element,
            charge,
            position,
        });
    }

    // Voxel data block: exactly `expected_points` whitespace-delimited
    // tokens, tolerant of any line-wrapping (never assumes a fixed count
    // per line).
    let mut values = Vec::with_capacity(expected_points.min(1_000_000));
    for i in 0..expected_points {
        let tok = feed.token()?.ok_or(CubeError::ValueCountMismatch {
            expected: expected_points,
            found: i,
        })?;
        let line_no = feed.line_no;
        let v = parse_finite(&tok, line_no, "voxel value")?;
        values.push(v);
    }
    if feed.token()?.is_some() {
        return Err(CubeError::TrailingData {
            after_values: expected_points,
        });
    }

    Ok(VolumetricGrid {
        origin,
        axes,
        shape,
        values,
        atoms,
        units,
    })
}

// ---------------------------------------------------------------------------
// Whole-string entry points
// ---------------------------------------------------------------------------

/// Parse a Gaussian Cube file with default limits ([`CubeParseLimits::default`]).
pub fn parse_cube(input: &str) -> Result<VolumetricGrid, CubeError> {
    parse_cube_with_limits(input, &CubeParseLimits::default())
}

/// Parse a Gaussian Cube file, enforcing `limits`.
pub fn parse_cube_with_limits(
    input: &str,
    limits: &CubeParseLimits,
) -> Result<VolumetricGrid, CubeError> {
    if input.len() > limits.max_input_bytes {
        return Err(CubeError::InputTooLarge {
            limit: limits.max_input_bytes,
        });
    }
    let mut lines = input.lines();
    let next_line =
        move || -> Result<Option<String>, CubeError> { Ok(lines.next().map(|s| s.to_string())) };
    parse_cube_from_feed(next_line, limits)
}

// ---------------------------------------------------------------------------
// Streaming reader
// ---------------------------------------------------------------------------

/// Streaming-*input* Cube-file reader over any [`std::io::BufRead`] source.
///
/// Unlike [`parse_cube`]/[`parse_cube_with_limits`] (which require the
/// whole file in memory as one `&str`), `CubeFileReader` reads line-by-line
/// directly from the source rather than materializing the whole raw file
/// text as one in-memory `String` first, suitable for the multi-gigabyte
/// cube files real quantum-chemistry workflows produce for large grids.
/// Header validation (the grid-point cap, the atom cap) happens *before*
/// the voxel data block is read, so a pathological header is rejected
/// before any large allocation is attempted. Follows the same shape as
/// [`crate::sdf::SdfFileReader`].
///
/// **Only the input reading is streaming** -- [`Self::read`] still returns
/// one fully-materialized [`VolumetricGrid`], whose `values: Vec<f64>`
/// holds every voxel in memory at once (this is inherent to
/// `VolumetricGrid`'s shape, not a limitation specific to this reader).
/// What this type avoids is the *doubled* memory of first reading the
/// entire raw text into a `String` and only then parsing it, plus giving
/// callers a bounded-per-line read path instead of requiring
/// `std::fs::read_to_string` up front.
pub struct CubeFileReader<R: std::io::BufRead> {
    reader: R,
    limits: CubeParseLimits,
}

impl<R: std::io::BufRead> CubeFileReader<R> {
    /// Wrap any `BufRead` source (e.g. `BufReader<File>`) with default
    /// limits ([`CubeParseLimits::default`]).
    pub fn new(reader: R) -> Self {
        Self::with_limits(reader, CubeParseLimits::default())
    }

    /// Wrap any `BufRead` source, enforcing `limits`.
    pub fn with_limits(reader: R, limits: CubeParseLimits) -> Self {
        Self { reader, limits }
    }

    /// Parse the wrapped source as a single Cube grid.
    pub fn read(self) -> Result<VolumetricGrid, CubeError> {
        let limits = self.limits;
        let mut reader = self.reader;
        let mut bytes_read: usize = 0;
        let next_line = move || -> Result<Option<String>, CubeError> {
            let mut buf = String::new();
            match reader.read_line(&mut buf) {
                Ok(0) => Ok(None),
                Ok(n) => {
                    bytes_read = bytes_read.saturating_add(n);
                    if bytes_read > limits.max_input_bytes {
                        return Err(CubeError::InputTooLarge {
                            limit: limits.max_input_bytes,
                        });
                    }
                    Ok(Some(buf.trim_end_matches(['\n', '\r']).to_string()))
                }
                Err(e) => Err(CubeError::Io(e.to_string())),
            }
        };
        parse_cube_from_feed(next_line, &limits)
    }
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Write a [`VolumetricGrid`] as a Gaussian Cube file. Always emits a
/// single-dataset file (positive/plain `NAtoms`, no `NVal` field). The
/// first axis line's sign is set from `grid.units` (negative for
/// [`GridUnits::Angstrom`], positive for [`GridUnits::Bohr`]) -- see module
/// docs on this convention's confidence level; numbers are emitted exactly
/// as stored, never rescaled.
pub fn write_cube(grid: &VolumetricGrid) -> Result<String, CubeError> {
    // Checked before `validate()`: a shape dimension too large for the
    // signed i64 an axis-line field requires is a distinct failure mode
    // from "values doesn't match the shape's point count" (which, for a
    // dimension this large, would itself require an infeasible
    // allocation to construct a matching `values` -- this check must be
    // reachable without one).
    let mut signed_shape = [0i64; 3];
    for (axis, count) in grid.shape.iter().enumerate() {
        signed_shape[axis] = i64::try_from(*count).map_err(|_| CubeError::DimensionOutOfRange {
            axis,
            value: *count,
        })?;
    }
    grid.validate()?;

    let mut out = String::new();
    out.push_str("chematic Gaussian Cube file\n");
    out.push_str("Written by chematic-mol::cube::write_cube\n");
    out.push_str(&format!(
        "{} {} {} {}\n",
        grid.atoms.len(),
        grid.origin[0],
        grid.origin[1],
        grid.origin[2]
    ));
    for (axis, magnitude) in signed_shape.into_iter().enumerate() {
        let signed_count: i64 = if axis == 0 && grid.units == GridUnits::Angstrom {
            -magnitude
        } else {
            magnitude
        };
        out.push_str(&format!(
            "{signed_count} {} {} {}\n",
            grid.axes[axis][0], grid.axes[axis][1], grid.axes[axis][2]
        ));
    }
    for atom in &grid.atoms {
        out.push_str(&format!(
            "{} {} {} {} {}\n",
            atom.element.atomic_number(),
            atom.charge,
            atom.position[0],
            atom.position[1],
            atom.position[2]
        ));
    }
    for (idx, v) in grid.values.iter().enumerate() {
        out.push_str(&v.to_string());
        if (idx + 1) % 6 == 0 {
            out.push('\n');
        } else {
            out.push(' ');
        }
    }
    if grid.values.is_empty() || !grid.values.len().is_multiple_of(6) {
        out.push('\n');
    }
    Ok(out)
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::volumetric::GridAtom;

    fn small_2x2x2() -> VolumetricGrid {
        VolumetricGrid {
            origin: [0.0, 0.0, 0.0],
            axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            shape: [2, 2, 2],
            values: vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            atoms: vec![GridAtom {
                element: Element::from_symbol("H").unwrap(),
                charge: 1.0,
                position: [0.5, 0.5, 0.5],
            }],
            units: GridUnits::Bohr,
        }
    }

    const HAND_WRITTEN_2X2X2_CUBE: &str = "Water density\n\
Generated for chematic tests\n\
1    0.000000    0.000000    0.000000\n\
2    1.000000    0.000000    0.000000\n\
2    0.000000    1.000000    0.000000\n\
2    0.000000    0.000000    1.000000\n\
8    8.000000    0.500000    0.500000    0.500000\n\
0.0 1.0 2.0 3.0\n\
4.0 5.0 6.0 7.0\n";

    #[test]
    fn parse_hand_written_2x2x2_fixture() {
        let grid = parse_cube(HAND_WRITTEN_2X2X2_CUBE).unwrap();
        assert_eq!(grid.shape, [2, 2, 2]);
        assert_eq!(grid.units, GridUnits::Bohr);
        assert_eq!(grid.values, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
        assert_eq!(grid.atoms.len(), 1);
        assert_eq!(grid.atoms[0].element, Element::from_symbol("O").unwrap());
        assert_eq!(grid.atoms[0].charge, 8.0);
    }

    #[test]
    fn round_trip_2x2x2() {
        let grid = small_2x2x2();
        let text = write_cube(&grid).unwrap();
        let parsed = parse_cube(&text).unwrap();
        assert_eq!(parsed, grid);
    }

    #[test]
    fn round_trip_3x2x4_non_cubic_shape() {
        let n = 3 * 2 * 4;
        let grid = VolumetricGrid {
            origin: [-1.0, -0.5, 0.25],
            axes: [[0.2, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.1]],
            shape: [3, 2, 4],
            values: (0..n).map(|i| i as f64 * 0.5).collect(),
            atoms: vec![
                GridAtom {
                    element: Element::from_symbol("C").unwrap(),
                    charge: 6.0,
                    position: [0.0, 0.0, 0.0],
                },
                GridAtom {
                    element: Element::from_symbol("N").unwrap(),
                    charge: 7.0,
                    position: [1.0, 1.0, 1.0],
                },
            ],
            units: GridUnits::Bohr,
        };
        let text = write_cube(&grid).unwrap();
        let parsed = parse_cube(&text).unwrap();
        assert_eq!(parsed, grid);
    }

    #[test]
    fn round_trip_non_orthogonal_axes() {
        let mut grid = small_2x2x2();
        // A genuinely non-axis-aligned set of step vectors.
        grid.axes = [[1.0, 0.2, 0.0], [0.1, 1.0, 0.3], [0.0, 0.1, 1.0]];
        let text = write_cube(&grid).unwrap();
        let parsed = parse_cube(&text).unwrap();
        assert_eq!(parsed, grid);
        assert_eq!(parsed.axes, grid.axes);
    }

    #[test]
    fn round_trip_angstrom_units_negative_axis_sign() {
        let mut grid = small_2x2x2();
        grid.units = GridUnits::Angstrom;
        let text = write_cube(&grid).unwrap();
        // First axis line must carry the negative sign for Angstrom.
        let axis1_line = text.lines().nth(3).unwrap();
        assert!(axis1_line.trim_start().starts_with('-'));
        let parsed = parse_cube(&text).unwrap();
        assert_eq!(parsed, grid);
        assert_eq!(parsed.units, GridUnits::Angstrom);
    }

    #[test]
    fn negative_natoms_multi_dataset_is_typed_rejected() {
        let input = "title\ncomment\n-1    0.0 0.0 0.0\n1 1.0 0.0 0.0\n1 0.0 1.0 0.0\n1 0.0 0.0 1.0\n6 6.0 0.0 0.0 0.0\n1 5\n0.0\n";
        let err = parse_cube(input).unwrap_err();
        assert!(matches!(
            err,
            CubeError::MultiDatasetUnsupported {
                natoms_field: -1,
                ..
            }
        ));
    }

    #[test]
    fn positive_natoms_with_nval_multi_dataset_is_typed_rejected() {
        let input = "title\ncomment\n1    0.0 0.0 0.0 2\n1 1.0 0.0 0.0\n1 0.0 1.0 0.0\n1 0.0 0.0 1.0\n6 6.0 0.0 0.0 0.0\n0.0 0.0\n";
        let err = parse_cube(input).unwrap_err();
        assert!(matches!(
            err,
            CubeError::MultiDatasetUnsupported {
                natoms_field: 1,
                nval: Some(2),
            }
        ));
    }

    #[test]
    fn nan_voxel_value_is_typed_rejected() {
        // "NaN"/"nan" is accepted by f64::from_str, so this is a
        // NonFiniteValue (parses, but not finite), not an InvalidNumber
        // (fails to parse at all).
        let input = HAND_WRITTEN_2X2X2_CUBE.replace("0.0 1.0 2.0 3.0", "NaN 1.0 2.0 3.0");
        assert!(input.contains("NaN"));
        let err = parse_cube(&input).unwrap_err();
        assert!(matches!(
            err,
            CubeError::NonFiniteValue {
                context: "voxel value",
                ..
            }
        ));
    }

    #[test]
    fn infinite_origin_is_typed_rejected() {
        let input = "title\ncomment\n1    inf 0.0 0.0\n1 1.0 0.0 0.0\n1 0.0 1.0 0.0\n1 0.0 0.0 1.0\n6 6.0 0.0 0.0 0.0\n0.0\n";
        let err = parse_cube(input).unwrap_err();
        assert!(matches!(
            err,
            CubeError::NonFiniteValue {
                context: "origin",
                ..
            }
        ));
    }

    #[test]
    fn write_rejects_nan_via_validate() {
        let mut grid = small_2x2x2();
        grid.values[0] = f64::NAN;
        let err = write_cube(&grid).unwrap_err();
        assert!(matches!(
            err,
            CubeError::Grid(GridError::NonFiniteValue { .. })
        ));
    }

    #[test]
    fn write_rejects_zero_dimension_shape() {
        let mut grid = small_2x2x2();
        grid.shape = [2, 0, 2];
        grid.values = Vec::new();
        let err = write_cube(&grid).unwrap_err();
        assert_eq!(err, CubeError::Grid(GridError::ZeroDimension { axis: 1 }));
    }

    #[test]
    fn write_rejects_dimension_exceeding_i64_max() {
        let mut grid = small_2x2x2();
        // A shape[0] this large can't be written as a signed i64 axis-line
        // field -- must be a typed error, not a wrapped/truncated cast.
        // This check runs (and must be reachable) before the full
        // values.len()-matches-shape validate() check, since actually
        // allocating a matching `values` for a shape this large is
        // infeasible -- `values` here is deliberately left empty/
        // mismatched.
        let huge = i64::MAX as usize + 1;
        grid.shape = [huge, 1, 1];
        let err = write_cube(&grid).unwrap_err();
        assert_eq!(
            err,
            CubeError::DimensionOutOfRange {
                axis: 0,
                value: huge
            }
        );
    }

    #[test]
    fn truncated_input_is_typed_error_not_panic() {
        let input = "title\ncomment\n1 0.0 0.0 0.0\n1 1.0 0.0 0.0\n";
        let err = parse_cube(input).unwrap_err();
        assert!(matches!(err, CubeError::UnexpectedEnd { .. }));
    }

    #[test]
    fn too_few_values_is_typed_mismatch() {
        let input = "title\ncomment\n1 0.0 0.0 0.0\n2 1.0 0.0 0.0\n2 0.0 1.0 0.0\n2 0.0 0.0 1.0\n6 6.0 0.0 0.0 0.0\n0.0 1.0 2.0\n";
        let err = parse_cube(input).unwrap_err();
        assert_eq!(
            err,
            CubeError::ValueCountMismatch {
                expected: 8,
                found: 3
            }
        );
    }

    #[test]
    fn too_many_values_is_typed_trailing_data_error() {
        let input = "title\ncomment\n1 0.0 0.0 0.0\n2 1.0 0.0 0.0\n2 0.0 1.0 0.0\n2 0.0 0.0 1.0\n6 6.0 0.0 0.0 0.0\n0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0\n";
        let err = parse_cube(input).unwrap_err();
        assert_eq!(err, CubeError::TrailingData { after_values: 8 });
    }

    #[test]
    fn pathological_header_overflow_is_typed_error() {
        // Each axis count (4_000_000_000) comfortably fits the i64 the axis
        // line is parsed as, but 4e9 * 4e9 * 2 = 3.2e19 overflows a 64-bit
        // usize's checked_mul (usize::MAX is ~1.84e19) -- a genuine
        // arithmetic-overflow header, distinct from "fits usize but exceeds
        // the configured cap" below.
        let input = "title\ncomment\n1 0.0 0.0 0.0\n4000000000 1.0 0.0 0.0\n4000000000 0.0 1.0 0.0\n2 0.0 0.0 1.0\n6 6.0 0.0 0.0 0.0\n";
        let err = parse_cube(input).unwrap_err();
        assert!(matches!(
            err,
            CubeError::Grid(GridError::ShapeOverflow { .. })
        ));
    }

    #[test]
    fn pathological_header_within_usize_but_over_cap_is_typed_error() {
        // 10_000 * 10_000 * 10_000 = 10^12 fits usize but exceeds the default
        // 100M-point cap without overflowing arithmetic -- a distinct failure
        // mode from the overflow case above.
        let input = "title\ncomment\n1 0.0 0.0 0.0\n10000 1.0 0.0 0.0\n10000 0.0 1.0 0.0\n10000 0.0 0.0 1.0\n6 6.0 0.0 0.0 0.0\n";
        let err = parse_cube(input).unwrap_err();
        assert!(matches!(
            err,
            CubeError::Grid(GridError::GridTooLarge {
                points: 1_000_000_000_000,
                limit: 100_000_000
            })
        ));
    }

    #[test]
    fn too_many_atoms_is_typed_error() {
        let limits = CubeParseLimits {
            max_atoms: 2,
            ..CubeParseLimits::default()
        };
        let input = "title\ncomment\n3 0.0 0.0 0.0\n1 1.0 0.0 0.0\n1 0.0 1.0 0.0\n1 0.0 0.0 1.0\n";
        let err = parse_cube_with_limits(input, &limits).unwrap_err();
        assert_eq!(err, CubeError::TooManyAtoms { limit: 2 });
    }

    #[test]
    fn negative_ny_axis_count_is_typed_rejected() {
        let input = "title\ncomment\n1 0.0 0.0 0.0\n1 1.0 0.0 0.0\n-1 0.0 1.0 0.0\n1 0.0 0.0 1.0\n6 6.0 0.0 0.0 0.0\n0.0\n";
        let err = parse_cube(input).unwrap_err();
        assert_eq!(err, CubeError::NonPositiveAxisCount { axis: 1, value: -1 });
    }

    #[test]
    fn unknown_atomic_number_is_typed_rejected() {
        let input = "title\ncomment\n1 0.0 0.0 0.0\n1 1.0 0.0 0.0\n1 0.0 1.0 0.0\n1 0.0 0.0 1.0\n999 999.0 0.0 0.0 0.0\n0.0\n";
        let err = parse_cube(input).unwrap_err();
        assert_eq!(
            err,
            CubeError::UnknownAtomicNumber {
                line: 7,
                value: 999
            }
        );
    }

    #[test]
    fn streaming_reader_matches_whole_string_parse() {
        use std::io::{BufReader, Cursor};
        let grid_from_str = parse_cube(HAND_WRITTEN_2X2X2_CUBE).unwrap();
        let cursor = Cursor::new(HAND_WRITTEN_2X2X2_CUBE.as_bytes().to_vec());
        let grid_from_reader = CubeFileReader::new(BufReader::new(cursor)).read().unwrap();
        assert_eq!(grid_from_str, grid_from_reader);
    }

    #[test]
    fn streaming_reader_enforces_byte_limit() {
        use std::io::{BufReader, Cursor};
        let limits = CubeParseLimits {
            max_input_bytes: 10,
            ..CubeParseLimits::default()
        };
        let cursor = Cursor::new(HAND_WRITTEN_2X2X2_CUBE.as_bytes().to_vec());
        let err = CubeFileReader::with_limits(BufReader::new(cursor), limits)
            .read()
            .unwrap_err();
        assert!(matches!(err, CubeError::InputTooLarge { limit: 10 }));
    }

    #[test]
    fn malformed_count_line_field_count_is_typed_error() {
        let input = "title\ncomment\n1 0.0 0.0 0.0 1 1\n";
        let err = parse_cube(input).unwrap_err();
        assert!(matches!(err, CubeError::InvalidCountLine { .. }));
    }
}
