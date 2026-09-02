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
    if let LammpsAtomStyle::Other(style) = &atom_style {
        return Err(LammpsDataError::UnsupportedAtomStyle {
            style: style.clone(),
        });
    }

    let mut cur = LineCursor::new(text);
    // Line 1: unconditionally ignored comment line.
    cur.next();

    let mut counts: Vec<(String, i64)> = Vec::new();
    let mut lo = [0.0f64; 3];
    let mut hi = [0.0f64; 3];
    let mut have_axis = [false; 3];
    let mut tilt: Option<[f64; 3]> = None;

    // Header block: classify lines until we hit the first section-name
    // line (classify_header_line returns Ok(None) for both blanks and
    // that boundary; blanks are skipped, the boundary line is left
    // un-consumed for the section loop below).
    loop {
        let Some(raw) = cur.peek() else {
            return Err(LammpsDataError::TruncatedInput {
                context: "header block (no sections found)",
            });
        };
        let line_no = cur.line_no();
        match classify_header_line(raw, line_no)? {
            Some(HeaderLine::Count { label, value }) => {
                if counts.len() >= limits.max_header_counts {
                    return Err(LammpsDataError::ResourceLimit {
                        resource: "header counts",
                        actual: counts.len() + 1,
                        limit: limits.max_header_counts,
                    });
                }
                counts.push((label, value));
                cur.next();
            }
            Some(HeaderLine::BoxPair { axis, lo: l, hi: h }) => {
                lo[axis] = l;
                hi[axis] = h;
                have_axis[axis] = true;
                cur.next();
            }
            Some(HeaderLine::Tilt { xy, xz, yz }) => {
                tilt = Some([xy, xz, yz]);
                cur.next();
            }
            None => {
                if strip_comment(raw).split_whitespace().next().is_none() {
                    // Blank line -- skip and keep scanning the header.
                    cur.next();
                    continue;
                }
                // First section-name line -- stop, leave it for the
                // section loop.
                break;
            }
        }
    }

    if !have_axis.iter().all(|&a| a) {
        return Err(LammpsDataError::InvalidBox {
            reason: "missing one or more of xlo/xhi, ylo/yhi, zlo/zhi bound lines".to_string(),
        });
    }
    let simulation_box = LammpsBox { lo, hi, tilt };
    simulation_box
        .validate()
        .map_err(|reason| LammpsDataError::InvalidBox { reason })?;

    let mut masses: Vec<LammpsMass> = Vec::new();
    let mut atoms: Vec<LammpsAtom> = Vec::new();
    // Source line number for each `atoms[i]`, index-parallel -- used only
    // for the `DuplicateAtomId` diagnostic below (not retained on
    // `LammpsData` itself, which has no line-number fields).
    let mut atom_lines: Vec<usize> = Vec::new();
    let mut velocities: Vec<LammpsVelocity> = Vec::new();
    let mut bonds: Vec<LammpsBond> = Vec::new();
    let mut unparsed_sections: Vec<(String, String)> = Vec::new();
    let mut section_count = 0;

    // Body: sequence of sections.
    while let Some(raw) = cur.peek() {
        let name = strip_comment(raw)
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ");
        if name.is_empty() {
            // Blank line between sections.
            cur.next();
            continue;
        }
        cur.next(); // consume the section-name line
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

        // Mandatory "next line is always skipped" per read_data.html.
        if cur.next().is_none() {
            return Err(LammpsDataError::TruncatedInput {
                context: "line after section header",
            });
        }

        // Collect this section's data rows: everything up to (but not
        // including) the next line whose first token isn't numeric, or
        // EOF. Blank lines within a section are skipped, not counted as
        // rows or as section boundaries.
        let mut rows: Vec<(usize, &str)> = Vec::new();
        let row_limit = match name.as_str() {
            "Masses" => limits.max_masses,
            "Atoms" => limits.max_atoms,
            "Velocities" => limits.max_velocities,
            "Bonds" => limits.max_bonds,
            _ => usize::MAX,
        };
        let mut opaque_bytes = 0usize;
        while let Some(row_raw) = cur.peek() {
            let stripped = strip_comment(row_raw);
            let Some(first_tok) = stripped.split_whitespace().next() else {
                cur.next(); // blank line inside a section
                continue;
            };
            if !looks_numeric(first_tok) {
                break; // next section header
            }
            if rows.len() >= row_limit {
                return Err(LammpsDataError::ResourceLimit {
                    resource: "section rows",
                    actual: rows.len() + 1,
                    limit: row_limit,
                });
            }
            opaque_bytes = opaque_bytes.saturating_add(row_raw.len());
            if name != "Masses"
                && name != "Atoms"
                && name != "Velocities"
                && name != "Bonds"
                && opaque_bytes > limits.max_opaque_section_bytes
            {
                return Err(LammpsDataError::ResourceLimit {
                    resource: "opaque section bytes",
                    actual: opaque_bytes,
                    limit: limits.max_opaque_section_bytes,
                });
            }
            rows.push((cur.line_no(), row_raw));
            cur.next();
        }

        match name.as_str() {
            "Masses" => {
                for (line_no, row_raw) in rows {
                    let content = strip_comment(row_raw);
                    let toks: Vec<&str> = content.split_whitespace().collect();
                    if toks.len() != 2 {
                        return Err(LammpsDataError::MalformedRow {
                            section: name.clone(),
                            line: line_no,
                            detail: row_raw.to_string(),
                        });
                    }
                    let atom_type = parse_int(toks[0], &name, line_no, row_raw)?;
                    let mass = parse_finite(toks[1], &name, line_no, row_raw)?;
                    masses.push(LammpsMass { atom_type, mass });
                }
            }
            "Atoms" => {
                let base = atom_style.base_field_count().expect(
                    "atom_style is one of the 4 supported styles here -- Other was rejected above",
                );
                for (line_no, row_raw) in rows {
                    let content = strip_comment(row_raw);
                    let toks: Vec<&str> = content.split_whitespace().collect();
                    let (fields, image) = if toks.len() == base {
                        (&toks[..], None)
                    } else if toks.len() == base + 3 {
                        let mut img = [0i32; 3];
                        for (k, slot) in img.iter_mut().enumerate() {
                            *slot = toks[base + k].parse().map_err(|_| {
                                LammpsDataError::MalformedRow {
                                    section: name.clone(),
                                    line: line_no,
                                    detail: row_raw.to_string(),
                                }
                            })?;
                        }
                        (&toks[..base], Some(img))
                    } else {
                        return Err(LammpsDataError::MalformedRow {
                            section: name.clone(),
                            line: line_no,
                            detail: row_raw.to_string(),
                        });
                    };

                    let id = parse_int(fields[0], &name, line_no, row_raw)?;
                    let (molecule_id, type_idx) = match atom_style {
                        LammpsAtomStyle::Molecular | LammpsAtomStyle::Full => {
                            (Some(parse_int(fields[1], &name, line_no, row_raw)?), 2)
                        }
                        _ => (None, 1),
                    };
                    let atom_type = parse_int(fields[type_idx], &name, line_no, row_raw)?;
                    let mut next = type_idx + 1;
                    let charge = match atom_style {
                        LammpsAtomStyle::Charge | LammpsAtomStyle::Full => {
                            let q = parse_finite(fields[next], &name, line_no, row_raw)?;
                            next += 1;
                            Some(q)
                        }
                        _ => None,
                    };
                    let x = parse_finite(fields[next], &name, line_no, row_raw)?;
                    let y = parse_finite(fields[next + 1], &name, line_no, row_raw)?;
                    let z = parse_finite(fields[next + 2], &name, line_no, row_raw)?;

                    atoms.push(LammpsAtom {
                        id,
                        molecule_id,
                        atom_type,
                        charge,
                        x,
                        y,
                        z,
                        image,
                    });
                    atom_lines.push(line_no);
                }
            }
            "Velocities" => {
                for (line_no, row_raw) in rows {
                    let content = strip_comment(row_raw);
                    let toks: Vec<&str> = content.split_whitespace().collect();
                    if toks.len() != 4 {
                        return Err(LammpsDataError::MalformedRow {
                            section: name.clone(),
                            line: line_no,
                            detail: row_raw.to_string(),
                        });
                    }
                    velocities.push(LammpsVelocity {
                        atom_id: parse_int(toks[0], &name, line_no, row_raw)?,
                        vx: parse_finite(toks[1], &name, line_no, row_raw)?,
                        vy: parse_finite(toks[2], &name, line_no, row_raw)?,
                        vz: parse_finite(toks[3], &name, line_no, row_raw)?,
                    });
                }
            }
            "Bonds" => {
                for (line_no, row_raw) in rows {
                    let content = strip_comment(row_raw);
                    let toks: Vec<&str> = content.split_whitespace().collect();
                    if toks.len() != 4 {
                        return Err(LammpsDataError::MalformedRow {
                            section: name.clone(),
                            line: line_no,
                            detail: row_raw.to_string(),
                        });
                    }
                    bonds.push(LammpsBond {
                        id: parse_int(toks[0], &name, line_no, row_raw)?,
                        bond_type: parse_int(toks[1], &name, line_no, row_raw)?,
                        atom1: parse_int(toks[2], &name, line_no, row_raw)?,
                        atom2: parse_int(toks[3], &name, line_no, row_raw)?,
                    });
                }
            }
            _ => {
                let raw_text = rows.iter().map(|(_, r)| *r).collect::<Vec<_>>().join("\n");
                unparsed_sections.push((name.clone(), raw_text));
            }
        }
    }

    // Cross-checks.
    if let Some(declared) = counts.iter().find(|(k, _)| k == "atoms").map(|(_, v)| *v) {
        if declared != atoms.len() as i64 {
            return Err(LammpsDataError::AtomCountMismatch {
                declared,
                actual: atoms.len(),
            });
        }
    } else {
        return Err(LammpsDataError::MalformedHeader {
            line: 0,
            detail: "missing required 'N atoms' header count line".to_string(),
        });
    }

    // Duplicate atom-ID check: O(n log n) via a sorted (id, source line)
    // copy -- no hash container anywhere in this module (project
    // determinism rule). Uses `atom_lines` (the real source line each
    // `atoms[i]` was parsed from), not the atom's position in `atoms`, so
    // the reported line numbers are genuinely useful for locating the
    // offending rows in the original file.
    let mut by_id: Vec<(i64, usize)> = atoms
        .iter()
        .zip(&atom_lines)
        .map(|(a, &line)| (a.id, line))
        .collect();
    by_id.sort_unstable_by_key(|(id, _)| *id);
    for w in by_id.windows(2) {
        if w[0].0 == w[1].0 {
            // `sort_unstable_by_key` doesn't preserve original relative
            // order for equal ids, so recover chronological order via
            // min/max rather than assuming `w[0]` was read first.
            let (first_line, line) = (w[0].1.min(w[1].1), w[0].1.max(w[1].1));
            return Err(LammpsDataError::DuplicateAtomId {
                id: w[0].0,
                first_line,
                line,
            });
        }
    }

    let mut sorted_ids: Vec<i64> = atoms.iter().map(|a| a.id).collect();
    sorted_ids.sort_unstable();
    for bond in &bonds {
        if sorted_ids.binary_search(&bond.atom1).is_err() {
            return Err(LammpsDataError::BondReferencesUnknownAtom {
                bond_id: bond.id,
                atom_id: bond.atom1,
            });
        }
        if sorted_ids.binary_search(&bond.atom2).is_err() {
            return Err(LammpsDataError::BondReferencesUnknownAtom {
                bond_id: bond.id,
                atom_id: bond.atom2,
            });
        }
    }

    Ok(LammpsData {
        counts,
        atom_style,
        simulation_box,
        masses,
        atoms,
        velocities,
        bonds,
        unparsed_sections,
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
mod tests {
    use super::*;

    fn atomic_fixture() -> &'static str {
        "comment line, ignored\n\
\n\
3 atoms\n\
1 atom types\n\
2 bonds\n\
1 bond types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Masses\n\
\n\
1 12.011\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
2 1 1.5 0.0 0.0\n\
3 1 3.0 0.0 0.0\n\
\n\
Bonds\n\
\n\
1 1 1 2\n\
2 1 2 3\n\
"
    }

    #[test]
    fn parses_atomic_fixture() {
        let d = parse_lammps_data(atomic_fixture(), LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(d.atoms.len(), 3);
        assert_eq!(d.bonds.len(), 2);
        assert_eq!(d.masses.len(), 1);
        assert_eq!(d.atoms[0].molecule_id, None);
        assert_eq!(d.atoms[0].charge, None);
        assert_eq!(d.count("atoms"), Some(3));
        assert_eq!(d.count("bond types"), Some(1));
    }

    #[test]
    fn round_trip_atomic_style() {
        let d = parse_lammps_data(atomic_fixture(), LammpsAtomStyle::Atomic).unwrap();
        let written = write_lammps_data(&d);
        let d2 = parse_lammps_data(&written, LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(d, d2);
    }

    fn charge_fixture() -> String {
        "comment\n\
\n\
2 atoms\n\
1 atom types\n\
\n\
-5.0 5.0 xlo xhi\n\
-5.0 5.0 ylo yhi\n\
-5.0 5.0 zlo zhi\n\
\n\
Atoms # charge\n\
\n\
1 1 -0.5 0.0 0.0 0.0\n\
2 1 0.5 1.0 0.0 0.0\n\
"
        .to_string()
    }

    #[test]
    fn round_trip_charge_style() {
        let d = parse_lammps_data(&charge_fixture(), LammpsAtomStyle::Charge).unwrap();
        assert_eq!(d.atoms[0].charge, Some(-0.5));
        assert_eq!(d.atoms[0].molecule_id, None);
        let written = write_lammps_data(&d);
        let d2 = parse_lammps_data(&written, LammpsAtomStyle::Charge).unwrap();
        assert_eq!(d, d2);
    }

    fn molecular_fixture() -> String {
        "comment\n\
\n\
2 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # molecular\n\
\n\
1 1 1 0.0 0.0 0.0\n\
2 1 1 1.0 0.0 0.0\n\
"
        .to_string()
    }

    #[test]
    fn round_trip_molecular_style() {
        let d = parse_lammps_data(&molecular_fixture(), LammpsAtomStyle::Molecular).unwrap();
        assert_eq!(d.atoms[0].molecule_id, Some(1));
        assert_eq!(d.atoms[0].charge, None);
        let written = write_lammps_data(&d);
        let d2 = parse_lammps_data(&written, LammpsAtomStyle::Molecular).unwrap();
        assert_eq!(d, d2);
    }

    fn full_fixture_with_images() -> String {
        "comment\n\
\n\
2 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # full\n\
\n\
1 1 1 -0.834 0.0 0.0 0.0 1 0 0\n\
2 1 1 0.417 1.0 0.0 0.0 0 -1 2\n\
"
        .to_string()
    }

    #[test]
    fn round_trip_full_style_with_image_flags() {
        let d = parse_lammps_data(&full_fixture_with_images(), LammpsAtomStyle::Full).unwrap();
        assert_eq!(d.atoms[0].molecule_id, Some(1));
        assert_eq!(d.atoms[0].charge, Some(-0.834));
        assert_eq!(d.atoms[0].image, Some([1, 0, 0]));
        assert_eq!(d.atoms[1].image, Some([0, -1, 2]));
        let written = write_lammps_data(&d);
        let d2 = parse_lammps_data(&written, LammpsAtomStyle::Full).unwrap();
        assert_eq!(d, d2);
    }

    #[test]
    fn round_trip_triclinic_box() {
        let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
1.0 2.0 0.5 xy xz yz\n\
\n\
Atoms # atomic\n\
\n\
1 1 5.0 5.0 5.0\n\
";
        let d = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(d.simulation_box.tilt, Some([1.0, 2.0, 0.5]));
        let written = write_lammps_data(&d);
        let d2 = parse_lammps_data(&written, LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(d, d2);
    }

    #[test]
    fn opaque_sections_round_trip_as_exact_fixed_point() {
        let text = "comment\n\
\n\
2 atoms\n\
1 atom types\n\
1 angle types\n\
1 dihedral types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
2 1 1.0 0.0 0.0\n\
\n\
Angles\n\
\n\
1 1 1 2 1 # a comment on an opaque row\n\
\n\
Angle Coeffs\n\
\n\
1 100.0 109.5\n\
";
        let d1 = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(
            d1.unparsed_sections,
            vec![
                (
                    "Angles".to_string(),
                    "1 1 1 2 1 # a comment on an opaque row".to_string()
                ),
                ("Angle Coeffs".to_string(), "1 100.0 109.5".to_string()),
            ]
        );
        let written = write_lammps_data(&d1);
        let d2 = parse_lammps_data(&written, LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(d2.unparsed_sections, d1.unparsed_sections);
        // Fixed point: parsing the rewritten file again changes nothing further.
        let written2 = write_lammps_data(&d2);
        let d3 = parse_lammps_data(&written2, LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(d3, d2);
    }

    #[test]
    fn zero_row_opaque_section_is_a_fixed_point() {
        let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
\n\
Pair Coeffs\n\
\n\
";
        let d1 = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(
            d1.unparsed_sections,
            vec![("Pair Coeffs".to_string(), String::new())]
        );
        let written = write_lammps_data(&d1);
        let d2 = parse_lammps_data(&written, LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(d2.unparsed_sections, d1.unparsed_sections);
    }

    #[test]
    fn atom_count_mismatch_is_typed_error() {
        let text = "comment\n\
\n\
5 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
2 1 1.0 0.0 0.0\n\
";
        let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
        assert_eq!(
            err,
            LammpsDataError::AtomCountMismatch {
                declared: 5,
                actual: 2
            }
        );
    }

    #[test]
    fn malformed_header_is_typed_error() {
        // "not" doesn't parse as numeric, so this reads as "no header
        // lines at all, first section is named 'not a valid header line
        // at all @@@'" -- no box bounds were ever declared, so this
        // surfaces as InvalidBox (missing bound lines), not
        // TruncatedInput.
        let text = "comment\nnot a valid header line at all @@@\n";
        let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
        assert!(matches!(err, LammpsDataError::InvalidBox { .. }));

        // A genuinely malformed header line: a leading number that
        // matches none of the known header shapes (not a `<count>
        // <label>` count line, not a 2-float box-bound pair, not a
        // 3-float tilt line).
        let text2 = "comment\n5 3 3 3 3 3 3\n";
        let err2 = parse_lammps_data(text2, LammpsAtomStyle::Atomic).unwrap_err();
        assert!(matches!(err2, LammpsDataError::MalformedHeader { .. }));
    }

    #[test]
    fn non_finite_value_is_typed_error() {
        let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 NaN 0.0 0.0\n\
";
        let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
        assert!(matches!(
            err,
            LammpsDataError::NonFiniteValue { section, .. } if section == "Atoms"
        ));
    }

    #[test]
    fn unsupported_atom_style_is_typed_error() {
        let err = parse_lammps_data(
            atomic_fixture(),
            LammpsAtomStyle::Other("sphere".to_string()),
        )
        .unwrap_err();
        assert_eq!(
            err,
            LammpsDataError::UnsupportedAtomStyle {
                style: "sphere".to_string()
            }
        );
    }

    #[test]
    fn duplicate_atom_id_is_typed_error() {
        let text = "comment\n\
\n\
2 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
1 1 1.0 0.0 0.0\n\
";
        let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
        // Pin the actual source line numbers (line 12: first `1 1 0.0 0.0
        // 0.0` row; line 13: the duplicate `1 1 1.0 0.0 0.0` row) so a
        // regression back to reporting the atom's position instead of its
        // real source line is caught.
        assert_eq!(
            err,
            LammpsDataError::DuplicateAtomId {
                id: 1,
                first_line: 12,
                line: 13,
            }
        );
    }

    #[test]
    fn bond_referencing_unknown_atom_is_typed_error() {
        let text = "comment\n\
\n\
2 atoms\n\
1 atom types\n\
1 bond types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
2 1 1.0 0.0 0.0\n\
\n\
Bonds\n\
\n\
1 1 1 99\n\
";
        let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
        assert_eq!(
            err,
            LammpsDataError::BondReferencesUnknownAtom {
                bond_id: 1,
                atom_id: 99
            }
        );
    }

    #[test]
    fn type_labels_framework_is_rejected_not_misparsed() {
        let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atom Type Labels\n\
\n\
1 C\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
";
        let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
        assert_eq!(
            err,
            LammpsDataError::TypeLabelsUnsupported {
                section: "Atom Type Labels".to_string()
            }
        );
    }

    #[test]
    fn trailing_comments_are_stripped_on_typed_rows() {
        let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
0.0 10.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Masses\n\
\n\
1 12.011 # carbon\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0 # salt ion\n\
";
        let d = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap();
        assert_eq!(d.masses[0].mass, 12.011);
        assert_eq!(d.atoms[0].x, 0.0);
    }

    #[test]
    fn invalid_box_is_typed_error() {
        let text = "comment\n\
\n\
1 atoms\n\
1 atom types\n\
\n\
10.0 0.0 xlo xhi\n\
0.0 10.0 ylo yhi\n\
0.0 10.0 zlo zhi\n\
\n\
Atoms # atomic\n\
\n\
1 1 0.0 0.0 0.0\n\
";
        let err = parse_lammps_data(text, LammpsAtomStyle::Atomic).unwrap_err();
        assert!(matches!(err, LammpsDataError::InvalidBox { .. }));
    }

    #[test]
    fn bounded_parser_rejects_input_and_line_limits() {
        let text = atomic_fixture();
        assert!(matches!(
            parse_lammps_data_with_limits(
                text,
                LammpsAtomStyle::Atomic,
                &LammpsDataParseLimits {
                    max_input_bytes: 8,
                    ..Default::default()
                }
            ),
            Err(LammpsDataError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));

        let long_line = format!("{}\n{}", "x".repeat(32), text);
        assert!(matches!(
            parse_lammps_data_with_limits(
                &long_line,
                LammpsAtomStyle::Atomic,
                &LammpsDataParseLimits {
                    max_line_bytes: 16,
                    ..Default::default()
                }
            ),
            Err(LammpsDataError::ResourceLimit {
                resource: "line bytes",
                ..
            })
        ));
    }

    #[test]
    fn bounded_parser_rejects_typed_section_rows() {
        let err = parse_lammps_data_with_limits(
            atomic_fixture(),
            LammpsAtomStyle::Atomic,
            &LammpsDataParseLimits {
                max_atoms: 0,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(matches!(
            err,
            LammpsDataError::ResourceLimit {
                resource: "section rows",
                ..
            }
        ));
    }
}
