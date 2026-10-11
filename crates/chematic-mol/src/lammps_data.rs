//! LAMMPS data file (`read_data` command format) reader and writer.
//!
//! ## Provenance
//!
//! Implemented independently from LAMMPS's own official manual
//! (<https://docs.lammps.org/read_data.html>, fetched and confirmed
//! directly against the live page -- not from memory). No source code,
//! comments, or tables were copied from LAMMPS itself, VMD, OVITO,
//! MDAnalysis, or any other tool.
//!
//! ## Structure
//!
//! A data file has: an unconditionally ignored first comment line; a
//! header block of `<count> <label>` lines (e.g. `120 atoms`,
//! `4 atom types`) and box-bounds lines (`xlo xhi`, `ylo yhi`, `zlo zhi`,
//! plus an optional `xy xz yz` triclinic-tilt line); then a sequence of
//! named sections, each introduced by a bare section-name line, a
//! mandatory skipped line (per read_data.html: "the next line is always
//! skipped" after a section keyword -- conventionally blank, but this
//! reader discards whatever is there unconditionally, matching LAMMPS's
//! own actual behavior), then the section's data rows, ending at the next
//! section-name line or EOF. Section order in the file is not fixed.
//!
//! `#` trailing comments are permitted by read_data.html on header lines,
//! section-keyword lines, and individual data rows ("There must be at
//! least one blank between valid content and the comment") -- stripped
//! here for header lines, section-name lines, and the 4 typed sections'
//! rows below. **Not** stripped for opaque (`unparsed_sections`) rows,
//! which are preserved byte-for-byte including any such comment -- see
//! "Sections not semantically parsed" below.
//!
//! ## `atom_style` is not inferable from the file -- do not guess it
//!
//! read_data.html states plainly that the atom style must be defined
//! before reading a data file; it is genuinely not always recoverable
//! from the `Atoms` section's column count alone. For example,
//! `atom_style charge` rows (`atom-ID atom-type q x y z`) and
//! `atom_style molecular` rows (`atom-ID molecule-ID atom-type x y z`)
//! are both 6 fields -- same column count, different meaning, genuinely
//! ambiguous by shape alone. [`parse_lammps_data`] therefore takes
//! `atom_style: LammpsAtomStyle` as a required parameter rather than a
//! "guess from column count" fallback, which would silently mis-parse
//! that ambiguous case.
//!
//! ## Supported atom styles
//!
//! Exactly 4, matching read_data.html's documented column layouts
//! verbatim:
//! - [`LammpsAtomStyle::Atomic`]: `atom-ID atom-type x y z`
//! - [`LammpsAtomStyle::Charge`]: `atom-ID atom-type q x y z`
//! - [`LammpsAtomStyle::Molecular`]: `atom-ID molecule-ID atom-type x y z`
//! - [`LammpsAtomStyle::Full`]: `atom-ID molecule-ID atom-type q x y z`
//!
//! All 4 optionally carry 3 trailing image flags `ix iy iz`
//! (read_data.html: "atom lines (all lines or none of them) can
//! optionally list 3 trailing integer values (nx,ny,nz)"), detected per
//! row from whether it has exactly 3 more whitespace-delimited fields
//! than the style's base count.
//!
//! Any `atom_style` outside this set of 4 -- passed as
//! [`LammpsAtomStyle::Other`] -- is rejected with a typed
//! [`LammpsDataError::UnsupportedAtomStyle`]; there is no best-effort
//! fallback parse.
//!
//! ## Sections fully parsed (typed)
//!
//! `Masses` (`atom-type mass` pairs), `Atoms` (per `atom_style` above),
//! `Velocities` (`atom-ID vx vy vz`), `Bonds`
//! (`bond-ID bond-type atom1-ID atom2-ID`).
//!
//! ## Sections preserved but not semantically parsed
//!
//! `Angles`, `Dihedrals`, `Impropers`, `Pair Coeffs`, `Bond Coeffs`,
//! `Angle Coeffs`, `Dihedral Coeffs`, `Improper Coeffs`, and any other
//! section name not listed above: captured verbatim (name + raw row
//! text, byte-for-byte, including any inline `#` comment) into
//! [`LammpsData::unparsed_sections`], an ordered `Vec<(String, String)>`
//! (not a `HashMap` -- this project's standing determinism rule; order is
//! the sections' original relative order among themselves, and matters).
//! [`write_lammps_data`] always writes the 4 typed sections first (in a
//! fixed canonical order: Masses, Atoms, Velocities, Bonds), followed by
//! the opaque sections in their original relative order. The exact
//! *interleaving* of typed and opaque sections in the source file is
//! therefore not preserved -- only each section's own content is. Real
//! LAMMPS itself does not require or preserve a specific section order
//! either, and `parse(write(parse(text)))` is still a fixed point of
//! `parse(text)` under this scheme.
//!
//! ## Why opaque preservation is safe *right now* -- and a constraint on future work
//!
//! See the "Mutation safety" note on [`LammpsData`] below: this module is
//! read/write/round-trip only in v1 (no atom-removal or atom-renumbering
//! API), which is exactly why treating `Angles`/`Dihedrals`/etc. as
//! opaque byte blobs is safe rather than merely convenient.
//!
//! ## LAMMPS's type-label framework: unsupported, fails closed
//!
//! LAMMPS also supports an `Atom Type Labels` / `Bond Type Labels` /
//! `Angle Type Labels` / `Dihedral Type Labels` / `Improper Type Labels`
//! extension under which `Masses`/`Atoms`/`Bonds`/etc. rows can key off a
//! string label (e.g. a `Masses` row `C 12.011` instead of `1 12.011`)
//! rather than a numeric type/ID. This module's section-boundary
//! detection (see below) and typed-row parsing are numeric-ID-only and
//! cannot safely distinguish a label-keyed data row from the start of a
//! new section -- silently misreading such a file would drop or
//! mis-split real data rather than failing loudly. Rather than attempt
//! that, any section whose name ends in `"Type Labels"` causes the whole
//! parse to fail closed with [`LammpsDataError::TypeLabelsUnsupported`]
//! as soon as it's seen, since its presence signals the rest of the file
//! may use the same string-keyed convention this module cannot safely
//! interpret.
//!
//! ## Units: not tracked, by design
//!
//! LAMMPS data files do not declare their unit system in-file (`units` is
//! a separate *input-script* command, never present in the data file
//! itself). This module performs **no** unit conversion or inference:
//! every numeric value is stored exactly as read, unit-agnostic. This is
//! a deliberate, disclosed non-goal -- matching this project's own
//! Cube/OpenDX lesson about never silently assuming a unit system -- not
//! an oversight.
//!
//! ## Section-boundary detection
//!
//! After a section-name line and its mandatory skipped line, this reader
//! consumes data rows until either EOF or a line whose first
//! whitespace-delimited token (after stripping any `#` comment) fails to
//! parse as a number -- every documented LAMMPS section-name line begins
//! with a non-numeric word (`Masses`, `Atoms`, `Pair Coeffs`, ...) while
//! every documented data row begins with a numeric ID/type index, so this
//! is a safe, generic boundary heuristic *given* the type-label framework
//! is rejected up front (see above) -- without that guard, a label-keyed
//! row (`C 12.011`) would look exactly like a new section header to this
//! heuristic.
//!
//! ## Box: a new type, not `chematic_crystal::Lattice`
//!
//! See [`LammpsBox`]'s own doc comment.

// ---------------------------------------------------------------------------
// LammpsBox
// ---------------------------------------------------------------------------

/// The simulation box shape shared by LAMMPS data files (this module) and
/// dump files ([`crate::lammps_dump`]).
///
/// Deliberately **not** `chematic_crystal::Lattice`: a `Lattice` is a
/// crystallographic unit-cell matrix (row-vector fractional/Cartesian
/// convention, positive-length/non-degenerate-angle validation, a cached
/// inverse) built for periodic *crystal* structures. LAMMPS's lo/hi/tilt
/// convention is a different, simpler shape -- an axis-aligned box plus 3
/// shear/tilt components -- and forcing that mapping is unnecessary
/// coupling for no benefit. `chematic-mol` only depends on
/// `chematic-crystal` behind its optional `crystal` feature (for the CIF
/// adapter); this type must work without it, same reasoning as
/// [`crate::volumetric::VolumetricGrid`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LammpsBox {
    /// `[xlo, ylo, zlo]`.
    pub lo: [f64; 3],
    /// `[xhi, yhi, zhi]`.
    pub hi: [f64; 3],
    /// `(xy, xz, yz)` tilt factors. `None` for an orthogonal box.
    pub tilt: Option<[f64; 3]>,
}

impl LammpsBox {
    /// Validate structural invariants: every `hi[i] > lo[i]`, and every
    /// `lo`/`hi`/`tilt` component is finite. Returns a plain `String`
    /// reason rather than a format-specific error type, so both
    /// [`LammpsDataError::InvalidBox`] (this module) and
    /// `LammpsDumpError::InvalidBox` ([`crate::lammps_dump`]) can wrap it
    /// without this type depending on either error enum.
    pub fn validate(&self) -> Result<(), String> {
        for i in 0..3 {
            if !self.lo[i].is_finite() {
                return Err(format!("lo[{i}] = {} is not finite", self.lo[i]));
            }
            if !self.hi[i].is_finite() {
                return Err(format!("hi[{i}] = {} is not finite", self.hi[i]));
            }
            if self.hi[i] <= self.lo[i] {
                return Err(format!(
                    "axis {i}: hi ({}) must be greater than lo ({})",
                    self.hi[i], self.lo[i]
                ));
            }
        }
        if let Some(t) = self.tilt {
            for (i, v) in t.iter().enumerate() {
                if !v.is_finite() {
                    return Err(format!("tilt[{i}] = {v} is not finite"));
                }
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// LammpsAtomStyle
// ---------------------------------------------------------------------------

/// Which `atom_style` a data file's `Atoms` section rows follow. See
/// module docs for why this must be supplied by the caller rather than
/// inferred.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LammpsAtomStyle {
    /// `atom-ID atom-type x y z`.
    Atomic,
    /// `atom-ID atom-type q x y z`.
    Charge,
    /// `atom-ID molecule-ID atom-type x y z`.
    Molecular,
    /// `atom-ID molecule-ID atom-type q x y z`.
    Full,
    /// Any style name other than the 4 above (e.g. `"sphere"`, `"bond"`,
    /// `"granular"`, `"full/omp"`, ...). Carried as the raw style name so
    /// [`parse_lammps_data`] can fail closed with a typed
    /// [`LammpsDataError::UnsupportedAtomStyle`] rather than requiring
    /// every caller to pre-filter before calling in -- see module docs.
    Other(String),
}

impl LammpsAtomStyle {
    /// The literal LAMMPS keyword for this style, used for the
    /// `Atoms # <style>` comment [`write_lammps_data`] emits (matching
    /// real `write_data` output convention; purely cosmetic -- readers of
    /// this module strip and ignore it, see module docs on comments).
    fn keyword(&self) -> &str {
        match self {
            Self::Atomic => "atomic",
            Self::Charge => "charge",
            Self::Molecular => "molecular",
            Self::Full => "full",
            Self::Other(s) => s,
        }
    }

    /// Number of whitespace-delimited fields an `Atoms` row has for this
    /// style, *not* counting optional trailing image flags. `None` for
    /// [`Self::Other`] (out of scope; see module docs).
    fn base_field_count(&self) -> Option<usize> {
        match self {
            Self::Atomic => Some(5),    // id type x y z
            Self::Charge => Some(6),    // id type q x y z
            Self::Molecular => Some(6), // id mol-id type x y z
            Self::Full => Some(7),      // id mol-id type q x y z
            Self::Other(_) => None,
        }
    }
}

// ---------------------------------------------------------------------------
// Records
// ---------------------------------------------------------------------------

/// One `Masses` section row: `atom-type mass`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LammpsMass {
    pub atom_type: i64,
    pub mass: f64,
}

/// One `Atoms` section row. Which of `molecule_id`/`charge` are `Some`
/// depends on the file's [`LammpsAtomStyle`] (stored once on
/// [`LammpsData::atom_style`], not per-atom, since one file has exactly
/// one atom style).
#[derive(Debug, Clone, PartialEq)]
pub struct LammpsAtom {
    pub id: i64,
    /// `Some` for [`LammpsAtomStyle::Molecular`]/[`LammpsAtomStyle::Full`].
    pub molecule_id: Option<i64>,
    pub atom_type: i64,
    /// `Some` for [`LammpsAtomStyle::Charge`]/[`LammpsAtomStyle::Full`].
    pub charge: Option<f64>,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// Trailing `ix iy iz` periodic-image flags, if present on this row.
    pub image: Option<[i32; 3]>,
}

/// One `Velocities` section row: `atom-ID vx vy vz`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LammpsVelocity {
    pub atom_id: i64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}

/// One `Bonds` section row: `bond-ID bond-type atom1-ID atom2-ID`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LammpsBond {
    pub id: i64,
    pub bond_type: i64,
    pub atom1: i64,
    pub atom2: i64,
}

// ---------------------------------------------------------------------------
// LammpsData
// ---------------------------------------------------------------------------

/// A parsed LAMMPS data file. This is a **standalone document type**, not
/// integrated with [`chematic_core::Molecule`] -- LAMMPS bonds/angles/etc.
/// are raw atom-index topology, not chemically perceived bonds, and an MD
/// atom under some `atom_style`s is not necessarily even a chemical
/// element. Callers needing a `Molecule` should build one themselves from
/// [`LammpsData::atoms`]/[`LammpsData::bonds`] plus whatever
/// element-inference or bond-perception policy fits their use case; this
/// module does not attempt that.
///
/// ## Mutation safety (read this before adding an atom-removal/renumbering API)
///
/// This module is **read/write/round-trip only in v1** -- there is no
/// atom-removal or atom-renumbering operation. That is deliberate, not an
/// oversight, and it is exactly what makes opaque preservation of
/// `Angles`/`Dihedrals`/`Impropers`/`*Coeffs`/etc. (see
/// [`Self::unparsed_sections`]) safe: those sections reference atom/type
/// indices by raw integer, which this module never interprets. If a
/// future PR adds an API that removes or renumbers atoms, it **must**
/// also either remap every opaque section's atom-index references or
/// reject the operation outright -- otherwise renumbering would silently
/// corrupt `unparsed_sections` into referencing atoms that no longer
/// exist (or the wrong ones). Since no such mutation API exists yet, that
/// staleness trap is unreachable today; it becomes reachable the moment
/// one is added. The same reasoning applies to [`Self::counts`]: it is
/// preserved and written back verbatim as parsed, so a mutation API would
/// also need to keep it consistent with the mutated sections' actual row
/// counts.
#[derive(Debug, Clone, PartialEq)]
pub struct LammpsData {
    /// Header `<count> <label>` lines, in file order, exactly as declared
    /// (e.g. `[("atoms", 120), ("bonds", 60), ("atom types", 4), ...]`).
    /// Looked up via [`Self::count`] rather than a fixed field per
    /// topology kind, since LAMMPS data files declare an open-ended set
    /// of these. [`write_lammps_data`] writes this back verbatim -- it is
    /// never recomputed from e.g. `atoms.len()`, so a hand-mutated
    /// `LammpsData` whose `atoms` no longer matches its stored `"atoms"`
    /// count will silently write an inconsistent file (see "Mutation
    /// safety" above).
    pub counts: Vec<(String, i64)>,
    pub atom_style: LammpsAtomStyle,
    pub simulation_box: LammpsBox,
    pub masses: Vec<LammpsMass>,
    pub atoms: Vec<LammpsAtom>,
    pub velocities: Vec<LammpsVelocity>,
    pub bonds: Vec<LammpsBond>,
    /// `(section_name, raw_row_text)` for every section this module does
    /// not semantically parse, in the sections' original relative order.
    /// `raw_row_text` is every data row of that section joined by `\n`
    /// (no trailing newline), byte-for-byte as read -- `#` comments and
    /// all. See module docs.
    pub unparsed_sections: Vec<(String, String)>,
}

/// Resource limits applied while parsing a LAMMPS data file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LammpsDataParseLimits {
    /// Maximum UTF-8 input size, in bytes.
    pub max_input_bytes: usize,
    /// Maximum physical line size, in bytes.
    pub max_line_bytes: usize,
    /// Maximum header count entries.
    pub max_header_counts: usize,
    pub max_masses: usize,
    pub max_atoms: usize,
    pub max_velocities: usize,
    pub max_bonds: usize,
    /// Maximum raw bytes retained for one opaque section.
    pub max_opaque_section_bytes: usize,
    /// Maximum number of section headers.
    pub max_sections: usize,
}

impl Default for LammpsDataParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_line_bytes: 1024 * 1024,
            max_header_counts: 256,
            max_masses: 1_000_000,
            max_atoms: 1_000_000,
            max_velocities: 1_000_000,
            max_bonds: 2_000_000,
            max_opaque_section_bytes: 64 * 1024 * 1024,
            max_sections: 256,
        }
    }
}

impl LammpsData {
    /// Look up a header count by its exact label (e.g. `"atoms"`,
    /// `"atom types"`). `None` if the header never declared that label.
    pub fn count(&self, label: &str) -> Option<i64> {
        self.counts
            .iter()
            .find(|(k, _)| k == label)
            .map(|(_, v)| *v)
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors from [`parse_lammps_data`].
#[derive(Debug, Clone, PartialEq)]
pub enum LammpsDataError {
    /// The requested [`LammpsAtomStyle`] is [`LammpsAtomStyle::Other`] --
    /// outside this module's 4 supported styles. See module docs.
    UnsupportedAtomStyle { style: String },
    /// A header line (count or box-bounds) did not match any recognized
    /// shape.
    MalformedHeader { line: usize, detail: String },
    /// [`LammpsBox::validate`] failed -- missing bound line(s), or a
    /// non-finite/non-increasing bound.
    InvalidBox { reason: String },
    /// A numeric field parsed but was NaN/Infinite.
    NonFiniteValue {
        section: String,
        line: usize,
        raw: String,
    },
    /// The header's declared `"atoms"` count did not match the number of
    /// rows actually found in the `Atoms` section.
    AtomCountMismatch { declared: i64, actual: usize },
    /// Two `Atoms` rows declared the same `atom-ID`.
    DuplicateAtomId {
        id: i64,
        first_line: usize,
        line: usize,
    },
    /// A `Bonds` row referenced an `atom-ID` not present in `Atoms`.
    BondReferencesUnknownAtom { bond_id: i64, atom_id: i64 },
    /// The input ended before a required section/field was fully read.
    TruncatedInput { context: &'static str },
    /// A row in one of the 4 typed sections (`Masses`/`Atoms`/
    /// `Velocities`/`Bonds`) did not have the expected field count or
    /// shape for that section.
    MalformedRow {
        section: String,
        line: usize,
        detail: String,
    },
    /// A section name ending in `"Type Labels"` was found -- LAMMPS's
    /// type-label framework, which this module cannot safely parse. See
    /// module docs.
    TypeLabelsUnsupported { section: String },
    /// The input exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl std::fmt::Display for LammpsDataError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedAtomStyle { style } => {
                write!(f, "unsupported LAMMPS atom_style '{style}'")
            }
            Self::MalformedHeader { line, detail } => {
                write!(f, "malformed header line at line {line}: '{detail}'")
            }
            Self::InvalidBox { reason } => write!(f, "invalid LAMMPS box: {reason}"),
            Self::NonFiniteValue { section, line, raw } => write!(
                f,
                "non-finite value '{raw}' in section '{section}' at line {line}"
            ),
            Self::AtomCountMismatch { declared, actual } => write!(
                f,
                "header declares {declared} atoms but Atoms section has {actual} rows"
            ),
            Self::DuplicateAtomId {
                id,
                first_line,
                line,
            } => write!(
                f,
                "duplicate atom-ID {id} at line {line} (first seen at line {first_line})"
            ),
            Self::BondReferencesUnknownAtom { bond_id, atom_id } => {
                write!(f, "bond {bond_id} references unknown atom-ID {atom_id}")
            }
            Self::TruncatedInput { context } => {
                write!(f, "unexpected end of input while reading {context}")
            }
            Self::MalformedRow {
                section,
                line,
                detail,
            } => write!(
                f,
                "malformed row in section '{section}' at line {line}: '{detail}'"
            ),
            Self::TypeLabelsUnsupported { section } => write!(
                f,
                "section '{section}' uses LAMMPS's type-label framework, which is not supported"
            ),
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(f, "{resource} has size {actual}, exceeding limit {limit}"),
        }
    }
}

impl std::error::Error for LammpsDataError {}

// ---------------------------------------------------------------------------
// Parsing helpers
// ---------------------------------------------------------------------------

/// Strip a trailing `#` comment (read_data.html: permitted on any line,
/// with at least one blank before `#`) and trailing whitespace.
fn strip_comment(line: &str) -> &str {
    match line.find('#') {
        Some(i) => line[..i].trim_end(),
        None => line.trim_end(),
    }
}

/// Whether `tok` parses as a number at all (int or float) -- used only for
/// the section-boundary heuristic (see module docs), never to decide a
/// value's actual type.
fn looks_numeric(tok: &str) -> bool {
    tok.parse::<f64>().is_ok()
}

struct LineCursor<'a> {
    lines: Vec<&'a str>,
    idx: usize,
}

impl<'a> LineCursor<'a> {
    fn new(input: &'a str) -> Self {
        Self {
            lines: input.lines().collect(),
            idx: 0,
        }
    }

    /// 1-based line number of the line `peek()`/`next()` would return.
    fn line_no(&self) -> usize {
        self.idx + 1
    }

    fn peek(&self) -> Option<&'a str> {
        self.lines.get(self.idx).copied()
    }

    fn next(&mut self) -> Option<&'a str> {
        let l = self.peek();
        if l.is_some() {
            self.idx += 1;
        }
        l
    }
}

fn parse_finite(
    tok: &str,
    section: &str,
    line: usize,
    raw_line: &str,
) -> Result<f64, LammpsDataError> {
    let v: f64 = tok.parse().map_err(|_| LammpsDataError::MalformedRow {
        section: section.to_string(),
        line,
        detail: raw_line.to_string(),
    })?;
    if !v.is_finite() {
        return Err(LammpsDataError::NonFiniteValue {
            section: section.to_string(),
            line,
            raw: tok.to_string(),
        });
    }
    Ok(v)
}

fn parse_int(
    tok: &str,
    section: &str,
    line: usize,
    raw_line: &str,
) -> Result<i64, LammpsDataError> {
    tok.parse::<i64>()
        .map_err(|_| LammpsDataError::MalformedRow {
            section: section.to_string(),
            line,
            detail: raw_line.to_string(),
        })
}

// ---------------------------------------------------------------------------
// Header parsing
// ---------------------------------------------------------------------------

enum HeaderLine {
    Count { label: String, value: i64 },
    BoxPair { axis: usize, lo: f64, hi: f64 },
    Tilt { xy: f64, xz: f64, yz: f64 },
}

/// Classify one header-block line. `Ok(None)` means `raw` is not a header
/// line at all (blank, or the first section-name line) -- the caller
/// decides what to do next.
fn classify_header_line(raw: &str, line_no: usize) -> Result<Option<HeaderLine>, LammpsDataError> {
    let content = strip_comment(raw);
    let toks: Vec<&str> = content.split_whitespace().collect();
    if toks.is_empty() {
        // Blank (post-comment-stripping) line: not a header line, but not
        // a section header either -- caller skips it.
        return Ok(None);
    }

    // Tilt line: 3 leading floats + literal "xy xz yz" keywords.
    if toks.len() == 6 && toks[3] == "xy" && toks[4] == "xz" && toks[5] == "yz" {
        let xy = parse_finite(toks[0], "header", line_no, raw)?;
        let xz = parse_finite(toks[1], "header", line_no, raw)?;
        let yz = parse_finite(toks[2], "header", line_no, raw)?;
        return Ok(Some(HeaderLine::Tilt { xy, xz, yz }));
    }

    // Box-bound pair line: 2 leading floats + a known keyword pair.
    if toks.len() == 4 {
        let axis = match (toks[2], toks[3]) {
            ("xlo", "xhi") => Some(0),
            ("ylo", "yhi") => Some(1),
            ("zlo", "zhi") => Some(2),
            _ => None,
        };
        if let Some(axis) = axis {
            let lo = parse_finite(toks[0], "header", line_no, raw)?;
            let hi = parse_finite(toks[1], "header", line_no, raw)?;
            return Ok(Some(HeaderLine::BoxPair { axis, lo, hi }));
        }
    }

    // Count line: leading integer + a non-numeric label (>= 1 word).
    if let Ok(value) = toks[0].parse::<i64>()
        && toks.len() >= 2
        && !looks_numeric(toks[1])
    {
        let label = toks[1..].join(" ");
        return Ok(Some(HeaderLine::Count { label, value }));
    }

    if looks_numeric(toks[0]) {
        // Starts with a number but matches none of the known header
        // shapes -- genuinely malformed, not "end of header".
        return Err(LammpsDataError::MalformedHeader {
            line: line_no,
            detail: raw.to_string(),
        });
    }

    // Non-numeric leading token that isn't blank: this is the first
    // section-name line, not a header line.
    Ok(None)
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

struct ParsedLammpsHeader {
    counts: Vec<(String, i64)>,
    simulation_box: LammpsBox,
}

#[derive(Default)]
struct ParsedLammpsSections {
    masses: Vec<LammpsMass>,
    atoms: Vec<LammpsAtom>,
    atom_lines: Vec<usize>,
    velocities: Vec<LammpsVelocity>,
    bonds: Vec<LammpsBond>,
    unparsed: Vec<(String, String)>,
}

fn validate_lammps_input(
    text: &str,
    atom_style: &LammpsAtomStyle,
    limits: &LammpsDataParseLimits,
) -> Result<(), LammpsDataError> {
    if text.len() > limits.max_input_bytes {
        return Err(LammpsDataError::ResourceLimit {
            resource: "input bytes",
            actual: text.len(),
            limit: limits.max_input_bytes,
        });
    }
    if let Some(line_bytes) = text.lines().map(str::len).max()
        && line_bytes > limits.max_line_bytes
    {
        return Err(LammpsDataError::ResourceLimit {
            resource: "line bytes",
            actual: line_bytes,
            limit: limits.max_line_bytes,
        });
    }
    if let LammpsAtomStyle::Other(style) = atom_style {
        return Err(LammpsDataError::UnsupportedAtomStyle {
            style: style.clone(),
        });
    }
    Ok(())
}

fn read_lammps_header(
    cursor: &mut LineCursor<'_>,
    limits: &LammpsDataParseLimits,
) -> Result<ParsedLammpsHeader, LammpsDataError> {
    cursor.next();
    let mut counts = Vec::new();
    let mut lo = [0.0; 3];
    let mut hi = [0.0; 3];
    let mut have_axis = [false; 3];
    let mut tilt = None;
    loop {
        let Some(raw) = cursor.peek() else {
            return Err(LammpsDataError::TruncatedInput {
                context: "header block (no sections found)",
            });
        };
        match classify_header_line(raw, cursor.line_no())? {
            Some(HeaderLine::Count { label, value }) => {
                if counts.len() >= limits.max_header_counts {
                    return Err(LammpsDataError::ResourceLimit {
                        resource: "header counts",
                        actual: counts.len() + 1,
                        limit: limits.max_header_counts,
                    });
                }
                counts.push((label, value));
                cursor.next();
            }
            Some(HeaderLine::BoxPair { axis, lo: l, hi: h }) => {
                lo[axis] = l;
                hi[axis] = h;
                have_axis[axis] = true;
                cursor.next();
            }
            Some(HeaderLine::Tilt { xy, xz, yz }) => {
                tilt = Some([xy, xz, yz]);
                cursor.next();
            }
            None if strip_comment(raw).split_whitespace().next().is_none() => {
                cursor.next();
            }
            None => break,
        }
    }
    if !have_axis.iter().all(|axis| *axis) {
        return Err(LammpsDataError::InvalidBox {
            reason: "missing one or more of xlo/xhi, ylo/yhi, zlo/zhi bound lines".into(),
        });
    }
    let simulation_box = LammpsBox { lo, hi, tilt };
    simulation_box
        .validate()
        .map_err(|reason| LammpsDataError::InvalidBox { reason })?;
    Ok(ParsedLammpsHeader {
        counts,
        simulation_box,
    })
}

fn section_row_limit(name: &str, limits: &LammpsDataParseLimits) -> usize {
    match name {
        "Masses" => limits.max_masses,
        "Atoms" => limits.max_atoms,
        "Velocities" => limits.max_velocities,
        "Bonds" => limits.max_bonds,
        _ => usize::MAX,
    }
}

fn is_typed_lammps_section(name: &str) -> bool {
    matches!(name, "Masses" | "Atoms" | "Velocities" | "Bonds")
}

fn read_section_rows<'a>(
    cursor: &mut LineCursor<'a>,
    name: &str,
    limits: &LammpsDataParseLimits,
) -> Result<Vec<(usize, &'a str)>, LammpsDataError> {
    let mut rows = Vec::new();
    let row_limit = section_row_limit(name, limits);
    let mut opaque_bytes = 0usize;
    while let Some(raw) = cursor.peek() {
        let Some(first) = strip_comment(raw).split_whitespace().next() else {
            cursor.next();
            continue;
        };
        if !looks_numeric(first) {
            break;
        }
        if rows.len() >= row_limit {
            return Err(LammpsDataError::ResourceLimit {
                resource: "section rows",
                actual: rows.len() + 1,
                limit: row_limit,
            });
        }
        opaque_bytes = opaque_bytes.saturating_add(raw.len());
        if !is_typed_lammps_section(name) && opaque_bytes > limits.max_opaque_section_bytes {
            return Err(LammpsDataError::ResourceLimit {
                resource: "opaque section bytes",
                actual: opaque_bytes,
                limit: limits.max_opaque_section_bytes,
            });
        }
        rows.push((cursor.line_no(), raw));
        cursor.next();
    }
    Ok(rows)
}

fn malformed_row(section: &str, line: usize, raw: &str) -> LammpsDataError {
    LammpsDataError::MalformedRow {
        section: section.into(),
        line,
        detail: raw.into(),
    }
}

fn parse_mass_rows(
    name: &str,
    rows: Vec<(usize, &str)>,
) -> Result<Vec<LammpsMass>, LammpsDataError> {
    rows.into_iter()
        .map(|(line, raw)| {
            let fields = strip_comment(raw).split_whitespace().collect::<Vec<_>>();
            if fields.len() != 2 {
                return Err(malformed_row(name, line, raw));
            }
            Ok(LammpsMass {
                atom_type: parse_int(fields[0], name, line, raw)?,
                mass: parse_finite(fields[1], name, line, raw)?,
            })
        })
        .collect()
}

fn parse_image_flags(
    fields: &[&str],
    base: usize,
    section: &str,
    line: usize,
    raw: &str,
) -> Result<([i32; 3], usize), LammpsDataError> {
    let mut image = [0; 3];
    for (offset, value) in image.iter_mut().enumerate() {
        *value = fields[base + offset]
            .parse()
            .map_err(|_| malformed_row(section, line, raw))?;
    }
    Ok((image, base))
}

fn parse_atom_row(
    style: &LammpsAtomStyle,
    section: &str,
    line: usize,
    raw: &str,
) -> Result<LammpsAtom, LammpsDataError> {
    let fields = strip_comment(raw).split_whitespace().collect::<Vec<_>>();
    let base = style
        .base_field_count()
        .expect("unsupported atom styles are rejected before parsing");
    let image = if fields.len() == base {
        None
    } else if fields.len() == base + 3 {
        Some(parse_image_flags(&fields, base, section, line, raw)?.0)
    } else {
        return Err(malformed_row(section, line, raw));
    };
    let fields = &fields[..base];
    let id = parse_int(fields[0], section, line, raw)?;
    let (molecule_id, type_index) = match style {
        LammpsAtomStyle::Molecular | LammpsAtomStyle::Full => {
            (Some(parse_int(fields[1], section, line, raw)?), 2)
        }
        _ => (None, 1),
    };
    let atom_type = parse_int(fields[type_index], section, line, raw)?;
    let mut coordinate_index = type_index + 1;
    let charge = match style {
        LammpsAtomStyle::Charge | LammpsAtomStyle::Full => {
            let charge = parse_finite(fields[coordinate_index], section, line, raw)?;
            coordinate_index += 1;
            Some(charge)
        }
        _ => None,
    };
    Ok(LammpsAtom {
        id,
        molecule_id,
        atom_type,
        charge,
        x: parse_finite(fields[coordinate_index], section, line, raw)?,
        y: parse_finite(fields[coordinate_index + 1], section, line, raw)?,
        z: parse_finite(fields[coordinate_index + 2], section, line, raw)?,
        image,
    })
}

fn parse_atom_rows(
    style: &LammpsAtomStyle,
    name: &str,
    rows: Vec<(usize, &str)>,
) -> Result<(Vec<LammpsAtom>, Vec<usize>), LammpsDataError> {
    let mut atoms = Vec::with_capacity(rows.len());
    let mut lines = Vec::with_capacity(rows.len());
    for (line, raw) in rows {
        atoms.push(parse_atom_row(style, name, line, raw)?);
        lines.push(line);
    }
    Ok((atoms, lines))
}

fn parse_velocity_rows(
    name: &str,
    rows: Vec<(usize, &str)>,
) -> Result<Vec<LammpsVelocity>, LammpsDataError> {
    rows.into_iter()
        .map(|(line, raw)| {
            let fields = strip_comment(raw).split_whitespace().collect::<Vec<_>>();
            if fields.len() != 4 {
                return Err(malformed_row(name, line, raw));
            }
            Ok(LammpsVelocity {
                atom_id: parse_int(fields[0], name, line, raw)?,
                vx: parse_finite(fields[1], name, line, raw)?,
                vy: parse_finite(fields[2], name, line, raw)?,
                vz: parse_finite(fields[3], name, line, raw)?,
            })
        })
        .collect()
}

fn parse_bond_rows(
    name: &str,
    rows: Vec<(usize, &str)>,
) -> Result<Vec<LammpsBond>, LammpsDataError> {
    rows.into_iter()
        .map(|(line, raw)| {
            let fields = strip_comment(raw).split_whitespace().collect::<Vec<_>>();
            if fields.len() != 4 {
                return Err(malformed_row(name, line, raw));
            }
            Ok(LammpsBond {
                id: parse_int(fields[0], name, line, raw)?,
                bond_type: parse_int(fields[1], name, line, raw)?,
                atom1: parse_int(fields[2], name, line, raw)?,
                atom2: parse_int(fields[3], name, line, raw)?,
            })
        })
        .collect()
}

fn parse_lammps_section(
    parsed: &mut ParsedLammpsSections,
    style: &LammpsAtomStyle,
    name: String,
    rows: Vec<(usize, &str)>,
) -> Result<(), LammpsDataError> {
    match name.as_str() {
        "Masses" => parsed.masses.extend(parse_mass_rows(&name, rows)?),
        "Atoms" => {
            let (atoms, lines) = parse_atom_rows(style, &name, rows)?;
            parsed.atoms.extend(atoms);
            parsed.atom_lines.extend(lines);
        }
        "Velocities" => parsed.velocities.extend(parse_velocity_rows(&name, rows)?),
        "Bonds" => parsed.bonds.extend(parse_bond_rows(&name, rows)?),
        _ => parsed.unparsed.push((
            name,
            rows.iter()
                .map(|(_, row)| *row)
                .collect::<Vec<_>>()
                .join("\n"),
        )),
    }
    Ok(())
}

fn read_lammps_sections(
    cursor: &mut LineCursor<'_>,
    style: &LammpsAtomStyle,
    limits: &LammpsDataParseLimits,
) -> Result<ParsedLammpsSections, LammpsDataError> {
    let mut parsed = ParsedLammpsSections::default();
    let mut section_count = 0;
    while let Some(raw) = cursor.peek() {
        let name = strip_comment(raw)
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ");
        if name.is_empty() {
            cursor.next();
            continue;
        }
        cursor.next();
        section_count += 1;
        if section_count > limits.max_sections {
            return Err(LammpsDataError::ResourceLimit {
                resource: "sections",
                actual: section_count,
                limit: limits.max_sections,
            });
        }
        if name.ends_with("Type Labels") {
            return Err(LammpsDataError::TypeLabelsUnsupported { section: name });
        }
        if cursor.next().is_none() {
            return Err(LammpsDataError::TruncatedInput {
                context: "line after section header",
            });
        }
        let rows = read_section_rows(cursor, &name, limits)?;
        parse_lammps_section(&mut parsed, style, name, rows)?;
    }
    Ok(parsed)
}

fn validate_lammps_references(
    counts: &[(String, i64)],
    parsed: &ParsedLammpsSections,
) -> Result<(), LammpsDataError> {
    let Some(declared) = counts
        .iter()
        .find(|(key, _)| key == "atoms")
        .map(|(_, n)| *n)
    else {
        return Err(LammpsDataError::MalformedHeader {
            line: 0,
            detail: "missing required 'N atoms' header count line".into(),
        });
    };
    if declared != parsed.atoms.len() as i64 {
        return Err(LammpsDataError::AtomCountMismatch {
            declared,
            actual: parsed.atoms.len(),
        });
    }
    let mut by_id = parsed
        .atoms
        .iter()
        .zip(&parsed.atom_lines)
        .map(|(atom, &line)| (atom.id, line))
        .collect::<Vec<_>>();
    by_id.sort_unstable_by_key(|(id, _)| *id);
    if let Some(pair) = by_id.windows(2).find(|pair| pair[0].0 == pair[1].0) {
        return Err(LammpsDataError::DuplicateAtomId {
            id: pair[0].0,
            first_line: pair[0].1.min(pair[1].1),
            line: pair[0].1.max(pair[1].1),
        });
    }
    let mut ids = parsed.atoms.iter().map(|atom| atom.id).collect::<Vec<_>>();
    ids.sort_unstable();
    for bond in &parsed.bonds {
        for atom_id in [bond.atom1, bond.atom2] {
            if ids.binary_search(&atom_id).is_err() {
                return Err(LammpsDataError::BondReferencesUnknownAtom {
                    bond_id: bond.id,
                    atom_id,
                });
            }
        }
    }
    Ok(())
}

/// Parse a LAMMPS data file. `atom_style` must be supplied by the caller
/// (see module docs -- it cannot be recovered from the file alone).
pub fn parse_lammps_data(
    text: &str,
    atom_style: LammpsAtomStyle,
) -> Result<LammpsData, LammpsDataError> {
    parse_lammps_data_with_limits(text, atom_style, &LammpsDataParseLimits::default())
}

/// Parse a LAMMPS data file with explicit resource limits.
pub fn parse_lammps_data_with_limits(
    text: &str,
    atom_style: LammpsAtomStyle,
    limits: &LammpsDataParseLimits,
) -> Result<LammpsData, LammpsDataError> {
    validate_lammps_input(text, &atom_style, limits)?;
    let mut cur = LineCursor::new(text);
    let header = read_lammps_header(&mut cur, limits)?;
    let parsed = read_lammps_sections(&mut cur, &atom_style, limits)?;
    validate_lammps_references(&header.counts, &parsed)?;
    Ok(LammpsData {
        counts: header.counts,
        atom_style,
        simulation_box: header.simulation_box,
        masses: parsed.masses,
        atoms: parsed.atoms,
        velocities: parsed.velocities,
        bonds: parsed.bonds,
        unparsed_sections: parsed.unparsed,
    })
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Write a [`LammpsData`] back out as a LAMMPS data file. Writes
/// [`LammpsData::counts`] verbatim (never recomputed from the typed
/// sections' actual lengths -- see [`LammpsData::counts`]'s doc comment),
/// then the box-bounds lines, then the 4 typed sections in a fixed
/// canonical order (`Masses`, `Atoms`, `Velocities`, `Bonds` -- each
/// omitted entirely when empty), then [`LammpsData::unparsed_sections`]
/// in their original relative order. See module docs on why section
/// *order* is not treated as significant round-trip content.
pub fn write_lammps_data(data: &LammpsData) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    let _ = writeln!(out, "LAMMPS data file written by chematic");
    out.push('\n');

    for (label, value) in &data.counts {
        let _ = writeln!(out, "{value} {label}");
    }
    out.push('\n');

    let b = &data.simulation_box;
    let _ = writeln!(out, "{} {} xlo xhi", b.lo[0], b.hi[0]);
    let _ = writeln!(out, "{} {} ylo yhi", b.lo[1], b.hi[1]);
    let _ = writeln!(out, "{} {} zlo zhi", b.lo[2], b.hi[2]);
    if let Some([xy, xz, yz]) = b.tilt {
        let _ = writeln!(out, "{xy} {xz} {yz} xy xz yz");
    }
    out.push('\n');

    if !data.masses.is_empty() {
        out.push_str("Masses\n\n");
        for m in &data.masses {
            let _ = writeln!(out, "{} {}", m.atom_type, m.mass);
        }
        out.push('\n');
    }

    if !data.atoms.is_empty() {
        let _ = writeln!(out, "Atoms # {}", data.atom_style.keyword());
        out.push('\n');
        for a in &data.atoms {
            let mut line = a.id.to_string();
            if let Some(mol) = a.molecule_id {
                let _ = write!(line, " {mol}");
            }
            let _ = write!(line, " {}", a.atom_type);
            if let Some(q) = a.charge {
                let _ = write!(line, " {q}");
            }
            let _ = write!(line, " {} {} {}", a.x, a.y, a.z);
            if let Some([ix, iy, iz]) = a.image {
                let _ = write!(line, " {ix} {iy} {iz}");
            }
            out.push_str(&line);
            out.push('\n');
        }
        out.push('\n');
    }

    if !data.velocities.is_empty() {
        out.push_str("Velocities\n\n");
        for v in &data.velocities {
            let _ = writeln!(out, "{} {} {} {}", v.atom_id, v.vx, v.vy, v.vz);
        }
        out.push('\n');
    }

    if !data.bonds.is_empty() {
        out.push_str("Bonds\n\n");
        for bd in &data.bonds {
            let _ = writeln!(out, "{} {} {} {}", bd.id, bd.bond_type, bd.atom1, bd.atom2);
        }
        out.push('\n');
    }

    for (name, raw) in &data.unparsed_sections {
        out.push_str(name);
        out.push_str("\n\n");
        if !raw.is_empty() {
            out.push_str(raw);
            out.push('\n');
        }
        out.push('\n');
    }

    out
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
#[path = "lammps_data/tests.rs"]
mod tests;
