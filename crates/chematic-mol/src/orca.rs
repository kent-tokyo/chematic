//! ORCA quantum-chemistry input/output support.
//!
//! ORCA (<https://orcaforum.kofo.mpg.de>) is a widely used QM package. This
//! module implements:
//!
//! | Format | Read | Write | Description |
//! |--------|------|-------|-------------|
//! | Input (`.inp`)  | Y | Y | keywords, raw `%...end` blocks, `*` coordinate block |
//! | Output (`.out` / `.log`) | Y | -- | geometry/trajectory, energy, frequencies, termination |
//!
//! ## Sources and confidence
//!
//! The primary spec source is ORCA's own public manual
//! (<https://orca-manual.mpi-muelheim.mpg.de/>), specifically the "General
//! Structure of the Input File", "Input of Coordinates", "Geometry
//! Optimizations" and "Vibrational Frequencies" chapters. The following
//! markers are quoted directly from that manual and are **high confidence**:
//!
//! - `!` keyword line(s) (multiple allowed, concatenated) and `%name ...
//!   end` blocks (general input structure chapter).
//! - `* xyz <charge> <mult> ... *`, `* xyzfile <charge> <mult> <file>`,
//!   `* int <charge> <mult> ...` coordinate block forms (coordinate input
//!   chapter).
//! - `#` end-of-line comments (general input structure chapter).
//! - `FINAL SINGLE POINT ENERGY <value>` (single-point energy tutorials).
//! - The `***********************HURRAY********************` /
//!   `THE OPTIMIZATION HAS CONVERGED` banner, and the
//!   `The optimization did not converge but reached the maximum number of
//!   optimization cycles` warning for non-convergence (geometry
//!   optimization chapter). The manual explicitly warns that
//!   `ORCA TERMINATED NORMALLY` alone does **not** imply the optimization
//!   converged -- this module treats termination and optimization
//!   convergence as two independent, separately reported fields for
//!   exactly that reason.
//! - The `-----------------------` / `VIBRATIONAL FREQUENCIES` /
//!   `-----------------------` header, `N:  <value> cm**-1` entry format,
//!   and the `***imaginary mode***` suffix on negative modes (vibrational
//!   frequencies chapter).
//! - The `---------------------------------` / `CARTESIAN COORDINATES
//!   (ANGSTROEM)` / `---------------------------------` geometry-block
//!   header (confirmed against a real ORCA run's output, not just the
//!   manual prose).
//!
//! The following are corroborated by multiple independent secondary
//! sources (real output snippets quoted in forum threads / tutorials / a
//! real reference output file) rather than a verbatim manual quote, and
//! are **medium confidence** -- flagged explicitly per this project's
//! stop-and-report-rather-than-guess discipline:
//!
//! - The per-atom inline frozen-coordinate marker: a bare `$` appended
//!   directly (no space) to an x/y/z value inside a `* xyz` block, e.g.
//!   `0.000000$`, freezes that Cartesian component during optimization.
//!   Supported on read/write; if real-world ORCA files use a different
//!   spacing/placement convention, this will fail to round-trip that
//!   specific input (typed parse error, not silent corruption).
//! - `****ORCA TERMINATED NORMALLY****` as the exact success banner text
//!   (matched here via substring, not exact asterisk count, for
//!   robustness).
//! - `ORCA finished by error termination in <MODULE>` as the error-path
//!   banner (quoted verbatim in several real bug-report/forum threads;
//!   ORCA's manual does not appear to document the exact wording).
//! - The `Total Charge ... Charge .... <n>` / `Multiplicity ... Mult ....
//!   <n>` echo lines near the top of the output.
//! - `%maxcore <n>` (and similarly `%moinp "file"`) being a single-line
//!   "global variable" directive with **no** `end` at all, unlike block
//!   directives such as `%pal ... end` -- corroborated by several
//!   independent HPC-site ORCA guides and the official ORCA tutorials
//!   site, though not found spelled out as a general rule in the manual
//!   itself.
//!
//! Deferred / explicitly out of scope (see individual doc comments below
//! for the exact boundary): `%coords` block as an alternative to `*
//! xyz`/`* int`, internal-coordinate (Z-matrix) semantic parsing, ghost
//! atoms/dummy atoms/point charges/fragment tags/per-atom basis overrides
//! in the coordinate block (preserved verbatim as an opaque trailing
//! string where they don't stop clean parsing, but not semantically
//! interpreted), full normal-mode vectors, and any energy other than the
//! plain (non `(MM)`/`(QM/MM)`-suffixed) `FINAL SINGLE POINT ENERGY`.
//!
//! A `%name ... end` block's raw content is found by scanning for the
//! *last* `end`-terminated line before the next top-level construct
//! (`%`/`*`/`!`) or EOF, not the first -- with no semantic knowledge of
//! any block's internal grammar. This is what lets a nested sub-block
//! such as `%geom`'s `Constraints ... end` (which closes with its own
//! `end` before `%geom`'s real closing `end`) round-trip correctly as
//! ordinary raw content: see `input_block_with_nested_end_round_trips` in
//! the tests below. A block with no `end` anywhere in that window and
//! nothing on the block-name line itself is a genuinely unterminated
//! block and fails with a typed [`OrcaInputError::UnterminatedBlock`],
//! never silently or by panicking.
//!
//! No bond perception is performed anywhere in this module: ORCA input and
//! output carry no bond table, only geometry + charge/multiplicity, and
//! [`Molecule`] values produced here never have bonds added.

use chematic_core::{Atom, Element, Molecule, MoleculeBuilder};

/// Return type shared by [`OrcaCoords::to_molecule`]:
/// `(Molecule, coords, charge, multiplicity)` -- same shape as
/// `gaussian.rs`'s `GjfResult`.
type OrcaXyzMolecule = (Molecule, Vec<(f64, f64, f64)>, i32, u32);

/// Resource limits for ORCA input parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OrcaInputParseLimits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_keywords: usize,
    pub max_blocks: usize,
    pub max_block_bytes: usize,
    pub max_atoms: usize,
}

impl Default for OrcaInputParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 * 1024 * 1024,
            max_line_bytes: 1024 * 1024,
            max_lines: 1_000_000,
            max_keywords: 100_000,
            max_blocks: 256,
            max_block_bytes: 16 * 1024 * 1024,
            max_atoms: 1_000_000,
        }
    }
}

// ===========================================================================
// Input: errors
// ===========================================================================

/// Errors from [`parse_orca_input`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OrcaInputError {
    /// A `%` block line had no block name after the `%`.
    MalformedBlock { line: usize },
    /// A `%name ... end` block never found its closing `end`.
    UnterminatedBlock { name: String },
    /// A `* ...` coordinate header line was missing its charge/multiplicity
    /// or (for `xyzfile`/`gzmtfile`) its filename.
    MalformedCoordHeader { line: usize },
    /// A `* xyz`/`* int` coordinate block never found its closing `*`.
    UnterminatedCoordBlock,
    /// The coordinate-block type token (after `*`) was not one of
    /// `xyz`/`xyzfile`/`gzmtfile`/`int`.
    UnknownCoordType { kind: String, line: usize },
    /// An atom line inside a `* xyz` block had fewer than 4 whitespace
    /// tokens (element + x + y + z).
    InvalidAtomLine { line: usize, detail: String },
    /// An atom line's element symbol was not recognized.
    UnknownElement { symbol: String, line: usize },
    /// A coordinate value could not be parsed as a float.
    InvalidCoordinate { line: usize, value: String },
    /// A coordinate value parsed but is NaN or infinite -- rejected rather
    /// than silently accepted.
    NonFiniteCoordinate { line: usize, value: String },
    /// A non-blank top-level line matched none of `#`/`!`/`%`/`*`.
    UnexpectedLine { line: usize, content: String },
    /// The input exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl std::fmt::Display for OrcaInputError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MalformedBlock { line } => write!(f, "malformed '%' block at line {line}"),
            Self::UnterminatedBlock { name } => {
                write!(f, "block '%{name}' has no closing 'end'")
            }
            Self::MalformedCoordHeader { line } => {
                write!(f, "malformed '*' coordinate header at line {line}")
            }
            Self::UnterminatedCoordBlock => {
                write!(f, "coordinate block has no closing '*'")
            }
            Self::UnknownCoordType { kind, line } => {
                write!(f, "unknown coordinate block type '{kind}' at line {line}")
            }
            Self::InvalidAtomLine { line, detail } => {
                write!(f, "invalid atom line at line {line}: {detail}")
            }
            Self::UnknownElement { symbol, line } => {
                write!(f, "unknown element symbol '{symbol}' at line {line}")
            }
            Self::InvalidCoordinate { line, value } => {
                write!(f, "invalid coordinate value '{value}' at line {line}")
            }
            Self::NonFiniteCoordinate { line, value } => {
                write!(
                    f,
                    "coordinate value '{value}' at line {line} is not finite (NaN/Infinite)"
                )
            }
            Self::UnexpectedLine { line, content } => {
                write!(f, "unexpected top-level line {line}: '{content}'")
            }
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(f, "{resource} has size {actual}, exceeding limit {limit}"),
        }
    }
}

impl std::error::Error for OrcaInputError {}

// ===========================================================================
// Input: data model
// ===========================================================================

/// One atom line inside a `* xyz` coordinate block.
#[derive(Debug, Clone, PartialEq)]
pub struct OrcaAtom {
    pub element: Element,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// Per-axis (x, y, z) frozen-coordinate marker: a bare `$` appended
    /// directly to the number (e.g. `0.000000$`). See the module docs for
    /// this convention's confidence level.
    pub frozen: [bool; 3],
    /// Any whitespace-joined trailing tokens after x/y/z on the atom line
    /// (e.g. a per-atom `NewGTO "..." end` basis override, a fragment tag,
    /// a ghost-atom `:`/point-charge `Q` marker fused onto a neighboring
    /// token, etc.). Preserved verbatim for round-trip; never
    /// semantically interpreted.
    pub extra: Option<String>,
}

/// The `* ... *` coordinate specification of an ORCA input file.
#[derive(Debug, Clone, PartialEq)]
pub enum OrcaCoords {
    /// `* xyz <charge> <mult>` / atom lines / `*` -- fully parsed.
    Xyz {
        charge: i32,
        multiplicity: u32,
        atoms: Vec<OrcaAtom>,
    },
    /// `* xyzfile <charge> <mult> <filename>` -- geometry lives in an
    /// external file; no atoms are present in this input to parse.
    XyzFile {
        charge: i32,
        multiplicity: u32,
        filename: String,
    },
    /// `* gzmtfile <charge> <mult> <filename>`.
    GzmtFile {
        charge: i32,
        multiplicity: u32,
        filename: String,
    },
    /// `* int <charge> <mult>` / internal-coordinate (Z-matrix) lines / `*`.
    /// Z-matrix semantics are out of scope (see module docs); the
    /// coordinate lines are preserved verbatim.
    Internal {
        charge: i32,
        multiplicity: u32,
        raw: String,
    },
}

impl OrcaCoords {
    /// Convert an [`OrcaCoords::Xyz`] block into a chematic [`Molecule`]
    /// (atoms only -- no bonds, ORCA input carries no bond table) plus its
    /// Cartesian coordinates, charge, and multiplicity. Returns `None` for
    /// every other variant (no atom list to convert).
    pub fn to_molecule(&self) -> Option<OrcaXyzMolecule> {
        match self {
            OrcaCoords::Xyz {
                charge,
                multiplicity,
                atoms,
            } => {
                let mut builder = MoleculeBuilder::new();
                let mut coords = Vec::with_capacity(atoms.len());
                for a in atoms {
                    builder.add_atom(Atom::new(a.element));
                    coords.push((a.x, a.y, a.z));
                }
                Some((builder.build(), coords, *charge, *multiplicity))
            }
            _ => None,
        }
    }
}

/// A `%name ... end` block (or a `%name value` single-line directive with
/// no `end` at all, e.g. `%maxcore 3000` -- see [`Self::has_end`]),
/// preserved losslessly as opaque raw text keyed by its (lower-cased)
/// block name -- this parser does not attempt to understand the internal
/// structure of `%scf`, `%geom`, `%basis`, etc.
#[derive(Debug, Clone, PartialEq)]
pub struct OrcaBlock {
    /// Lower-cased block name (ORCA input is case-insensitive), e.g.
    /// `"scf"`, `"geom"`, `"pal"`, `"maxcore"`.
    pub name: String,
    /// Raw inner text, trimmed of leading/trailing whitespace but
    /// otherwise verbatim (internal line breaks and indentation
    /// preserved). For a `has_end: true` block this is everything between
    /// the block name and the closing `end` (which may itself contain
    /// other `end`-terminated lines belonging to a nested sub-block, such
    /// as `%geom`'s `Constraints ... end` -- this block's own `end` is
    /// found as the *last* `end`-terminated line before the next
    /// top-level construct, not the first, precisely so nested `end`s
    /// like that round-trip as ordinary raw content instead of being
    /// mistaken for the terminator). For a `has_end: false` block this is
    /// whatever followed the block name on its own line.
    pub raw: String,
    /// Whether this block is closed with an explicit `end` keyword.
    /// `false` for ORCA's single-line "global variable" directives that
    /// take no `end` at all (e.g. `%maxcore 3000`, `%moinp "file.gbw"`).
    pub has_end: bool,
}

/// A parsed ORCA input file.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct OrcaInput {
    /// Whole-line `#...` comments, verbatim (including the leading `#`),
    /// in original order. Inline trailing comments on `!`/block/atom lines
    /// are stripped during structured-field extraction and are not
    /// separately preserved.
    pub comments: Vec<String>,
    /// All keywords from every `!` line, concatenated in order (multiple
    /// `!` lines are allowed by ORCA and are flattened here).
    pub keywords: Vec<String>,
    /// All `%name ... end` blocks, in original order.
    pub blocks: Vec<OrcaBlock>,
    /// The `* ... *` coordinate specification, if present.
    pub coords: Option<OrcaCoords>,
}

// ===========================================================================
// Input: parser
// ===========================================================================

/// Parse an ORCA input file (`.inp`).
pub fn parse_orca_input(input: &str) -> Result<OrcaInput, OrcaInputError> {
    parse_orca_input_with_limits(input, &OrcaInputParseLimits::default())
}

/// Parse an ORCA input file with explicit resource limits.
pub fn parse_orca_input_with_limits(
    input: &str,
    limits: &OrcaInputParseLimits,
) -> Result<OrcaInput, OrcaInputError> {
    if input.len() > limits.max_input_bytes {
        return Err(OrcaInputError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let lines: Vec<&str> = input.lines().collect();
    if lines.len() > limits.max_lines {
        return Err(OrcaInputError::ResourceLimit {
            resource: "lines",
            actual: lines.len(),
            limit: limits.max_lines,
        });
    }
    if let Some(line_bytes) = lines.iter().map(|line| line.len()).max()
        && line_bytes > limits.max_line_bytes
    {
        return Err(OrcaInputError::ResourceLimit {
            resource: "line bytes",
            actual: line_bytes,
            limit: limits.max_line_bytes,
        });
    }
    let mut out = OrcaInput::default();
    let mut i = 0usize;
    while i < lines.len() {
        let trimmed = lines[i].trim();
        if trimmed.is_empty() {
            i += 1;
            continue;
        }
        if let Some(comment) = trimmed.strip_prefix('#') {
            out.comments.push(format!("#{comment}"));
            i += 1;
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('!') {
            let rest = rest.split('#').next().unwrap_or("");
            out.keywords
                .extend(rest.split_whitespace().map(|s| s.to_string()));
            if out.keywords.len() > limits.max_keywords {
                return Err(OrcaInputError::ResourceLimit {
                    resource: "keywords",
                    actual: out.keywords.len(),
                    limit: limits.max_keywords,
                });
            }
            i += 1;
            continue;
        }
        if trimmed.starts_with('%') {
            let (block, next_i) = parse_block(&lines, i)?;
            if out.blocks.len() >= limits.max_blocks {
                return Err(OrcaInputError::ResourceLimit {
                    resource: "blocks",
                    actual: out.blocks.len() + 1,
                    limit: limits.max_blocks,
                });
            }
            if block.raw.len() > limits.max_block_bytes {
                return Err(OrcaInputError::ResourceLimit {
                    resource: "block bytes",
                    actual: block.raw.len(),
                    limit: limits.max_block_bytes,
                });
            }
            out.blocks.push(block);
            i = next_i;
            continue;
        }
        if trimmed.starts_with('*') {
            let (coords, next_i) = parse_coords(&lines, i)?;
            if let OrcaCoords::Xyz { atoms, .. } = &coords
                && atoms.len() > limits.max_atoms
            {
                return Err(OrcaInputError::ResourceLimit {
                    resource: "atoms",
                    actual: atoms.len(),
                    limit: limits.max_atoms,
                });
            }
            out.coords = Some(coords);
            i = next_i;
            continue;
        }
        return Err(OrcaInputError::UnexpectedLine {
            line: i + 1,
            content: trimmed.to_string(),
        });
    }
    Ok(out)
}

/// Returns the byte offset (within `line`) of a trailing standalone `end`
/// token (case-insensitive), if the line's last whitespace-delimited token
/// is `end`. ORCA's block terminator keyword is always the final token on
/// its line, whether alone on its own line or trailing other content
/// (e.g. `%pal nprocs 4 end`).
fn line_ends_with_end(line: &str) -> Option<usize> {
    let trimmed_end = line.trim_end();
    let last_tok_start = trimmed_end
        .rfind(char::is_whitespace)
        .map(|p| p + 1)
        .unwrap_or(0);
    let last_tok = &trimmed_end[last_tok_start..];
    if !last_tok.is_empty() && last_tok.eq_ignore_ascii_case("end") {
        Some(last_tok_start)
    } else {
        None
    }
}

fn parse_block(lines: &[&str], start: usize) -> Result<(OrcaBlock, usize), OrcaInputError> {
    let trimmed = lines[start].trim_start();
    let after_pct = &trimmed[1..]; // drop leading '%'
    let name_end = after_pct
        .find(char::is_whitespace)
        .unwrap_or(after_pct.len());
    let name = after_pct[..name_end].to_string();
    if name.is_empty() {
        return Err(OrcaInputError::MalformedBlock { line: start + 1 });
    }
    let rest_of_first_line = after_pct[name_end..].trim_start();

    // Window: from `start` up to (not including) the next top-level
    // construct line (`%`/`*`/`!`) or EOF. Bounds how far we search for
    // this block's closing `end` so a later, unrelated block's `end` can
    // never be misattributed to this one.
    let mut window_end = lines.len();
    for (k, l) in lines.iter().enumerate().skip(start + 1) {
        let t = l.trim();
        if t.starts_with('%') || t.starts_with('*') || t.starts_with('!') {
            window_end = k;
            break;
        }
    }

    // `content[0]` is whatever followed the block name on its own line;
    // `content[k]` for `k >= 1` is `lines[start + k]`.
    let mut content: Vec<&str> = Vec::with_capacity(window_end - start);
    content.push(rest_of_first_line);
    content.extend_from_slice(&lines[start + 1..window_end]);

    // The block's true terminator is the LAST `end`-terminated line in
    // the window, not the first: some blocks contain a nested sub-block
    // (e.g. `%geom`'s `Constraints ... end`) that closes with its own
    // `end` before the block's real closing `end`, and there is no
    // purely syntactic way to tell those apart -- but the block's own
    // `end` is always the last one before the next top-level construct,
    // so taking the last match handles the nested case correctly without
    // needing any semantic knowledge of what "Constraints" means.
    match content
        .iter()
        .rposition(|l| line_ends_with_end(l).is_some())
    {
        Some(end_idx) => {
            let end_pos = line_ends_with_end(content[end_idx]).unwrap();
            let mut collected: Vec<String> =
                content[..end_idx].iter().map(|s| s.to_string()).collect();
            let before = content[end_idx][..end_pos].trim_end();
            if !before.is_empty() {
                collected.push(before.to_string());
            }
            let raw = collected.join("\n").trim().to_string();
            let next_i = if end_idx == 0 {
                start + 1
            } else {
                start + end_idx + 1
            };
            Ok((
                OrcaBlock {
                    name: name.to_lowercase(),
                    raw,
                    has_end: true,
                },
                next_i,
            ))
        }
        None if !rest_of_first_line.is_empty() => {
            // No `end` anywhere in the window, but the block name has a
            // value directly on its own line: one of ORCA's single-line
            // "global variable" directives with no `end` at all (e.g.
            // `%maxcore 3000`, `%moinp "file.gbw"`). Consumes only the
            // block-name line itself.
            Ok((
                OrcaBlock {
                    name: name.to_lowercase(),
                    raw: rest_of_first_line.to_string(),
                    has_end: false,
                },
                start + 1,
            ))
        }
        None => {
            // No `end` anywhere in the window, and nothing on the block's
            // own line either: a genuinely unterminated multi-line block.
            Err(OrcaInputError::UnterminatedBlock { name })
        }
    }
}

fn parse_charge_mult(
    charge_tok: Option<&str>,
    mult_tok: Option<&str>,
    line: usize,
) -> Result<(i32, u32), OrcaInputError> {
    let charge_tok = charge_tok.ok_or(OrcaInputError::MalformedCoordHeader { line })?;
    let mult_tok = mult_tok.ok_or(OrcaInputError::MalformedCoordHeader { line })?;
    let charge: i32 = charge_tok
        .parse()
        .map_err(|_| OrcaInputError::MalformedCoordHeader { line })?;
    let multiplicity: u32 = mult_tok
        .parse()
        .map_err(|_| OrcaInputError::MalformedCoordHeader { line })?;
    Ok((charge, multiplicity))
}

fn parse_coords(lines: &[&str], start: usize) -> Result<(OrcaCoords, usize), OrcaInputError> {
    let trimmed = lines[start].trim();
    let after_star = trimmed[1..].trim_start();
    let mut tokens = after_star.split_whitespace();
    let kind = tokens.next().unwrap_or("").to_ascii_lowercase();
    let (charge, multiplicity) = parse_charge_mult(tokens.next(), tokens.next(), start + 1)?;

    match kind.as_str() {
        "xyz" => {
            if tokens.next().is_some() {
                return Err(OrcaInputError::MalformedCoordHeader { line: start + 1 });
            }
            let mut atoms = Vec::new();
            let mut i = start + 1;
            loop {
                if i >= lines.len() {
                    return Err(OrcaInputError::UnterminatedCoordBlock);
                }
                let t = lines[i].trim();
                if t == "*" {
                    i += 1;
                    break;
                }
                if t.is_empty() {
                    i += 1;
                    continue;
                }
                atoms.push(parse_atom_line(t, i + 1)?);
                i += 1;
            }
            Ok((
                OrcaCoords::Xyz {
                    charge,
                    multiplicity,
                    atoms,
                },
                i,
            ))
        }
        "xyzfile" | "gzmtfile" => {
            let filename = tokens.collect::<Vec<_>>().join(" ");
            if filename.is_empty() {
                return Err(OrcaInputError::MalformedCoordHeader { line: start + 1 });
            }
            let coords = if kind == "xyzfile" {
                OrcaCoords::XyzFile {
                    charge,
                    multiplicity,
                    filename,
                }
            } else {
                OrcaCoords::GzmtFile {
                    charge,
                    multiplicity,
                    filename,
                }
            };
            Ok((coords, start + 1))
        }
        "int" => {
            let mut raw_lines = Vec::new();
            let mut i = start + 1;
            loop {
                if i >= lines.len() {
                    return Err(OrcaInputError::UnterminatedCoordBlock);
                }
                let t = lines[i].trim();
                if t == "*" {
                    i += 1;
                    break;
                }
                raw_lines.push(lines[i].to_string());
                i += 1;
            }
            Ok((
                OrcaCoords::Internal {
                    charge,
                    multiplicity,
                    raw: raw_lines.join("\n").trim().to_string(),
                },
                i,
            ))
        }
        other => Err(OrcaInputError::UnknownCoordType {
            kind: other.to_string(),
            line: start + 1,
        }),
    }
}

fn parse_atom_line(line: &str, line_no: usize) -> Result<OrcaAtom, OrcaInputError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 4 {
        return Err(OrcaInputError::InvalidAtomLine {
            line: line_no,
            detail: line.to_string(),
        });
    }
    let element =
        Element::from_symbol(tokens[0]).ok_or_else(|| OrcaInputError::UnknownElement {
            symbol: tokens[0].to_string(),
            line: line_no,
        })?;
    let mut coords = [0.0f64; 3];
    let mut frozen = [false; 3];
    for axis in 0..3 {
        let tok = tokens[axis + 1];
        let (numeric, is_frozen) = match tok.strip_suffix('$') {
            Some(n) => (n, true),
            None => (tok, false),
        };
        let v: f64 = numeric
            .parse()
            .map_err(|_| OrcaInputError::InvalidCoordinate {
                line: line_no,
                value: tok.to_string(),
            })?;
        if !v.is_finite() {
            return Err(OrcaInputError::NonFiniteCoordinate {
                line: line_no,
                value: tok.to_string(),
            });
        }
        coords[axis] = v;
        frozen[axis] = is_frozen;
    }
    let extra = if tokens.len() > 4 {
        Some(tokens[4..].join(" "))
    } else {
        None
    };
    Ok(OrcaAtom {
        element,
        x: coords[0],
        y: coords[1],
        z: coords[2],
        frozen,
        extra,
    })
}

// ===========================================================================
// Input: writer
// ===========================================================================

/// Write an [`OrcaInput`] back out as ORCA input-file text.
///
/// Ordering is deterministic: comments, then the keyword line, then blocks
/// in their stored order, then the coordinate block -- never dependent on
/// hash-map iteration order (no hash maps are used in this data model at
/// all).
pub fn write_orca_input(input: &OrcaInput) -> String {
    let mut out = String::new();
    for c in &input.comments {
        out.push_str(c);
        out.push('\n');
    }
    if !input.keywords.is_empty() {
        if !out.is_empty() {
            out.push('\n');
        }
        out.push_str("! ");
        out.push_str(&input.keywords.join(" "));
        out.push('\n');
    }
    for b in &input.blocks {
        out.push('\n');
        out.push('%');
        out.push_str(&b.name);
        if !b.has_end {
            // Single-line directive with no `end` at all (e.g. `%maxcore
            // 3000`).
            if !b.raw.is_empty() {
                out.push(' ');
                out.push_str(&b.raw);
            }
            out.push('\n');
        } else if b.raw.is_empty() {
            out.push_str(" end\n");
        } else {
            out.push('\n');
            out.push_str(&b.raw);
            out.push_str("\nend\n");
        }
    }
    if let Some(coords) = &input.coords {
        out.push('\n');
        write_coords(coords, &mut out);
    }
    out
}

fn write_coords(coords: &OrcaCoords, out: &mut String) {
    match coords {
        OrcaCoords::Xyz {
            charge,
            multiplicity,
            atoms,
        } => {
            out.push_str(&format!("* xyz {charge} {multiplicity}\n"));
            for a in atoms {
                out.push_str(a.element.symbol());
                for (v, frozen) in [a.x, a.y, a.z].iter().zip(a.frozen.iter()) {
                    out.push_str(&format!(" {v:.6}"));
                    if *frozen {
                        out.push('$');
                    }
                }
                if let Some(extra) = &a.extra {
                    out.push(' ');
                    out.push_str(extra);
                }
                out.push('\n');
            }
            out.push_str("*\n");
        }
        OrcaCoords::XyzFile {
            charge,
            multiplicity,
            filename,
        } => {
            out.push_str(&format!("* xyzfile {charge} {multiplicity} {filename}\n"));
        }
        OrcaCoords::GzmtFile {
            charge,
            multiplicity,
            filename,
        } => {
            out.push_str(&format!("* gzmtfile {charge} {multiplicity} {filename}\n"));
        }
        OrcaCoords::Internal {
            charge,
            multiplicity,
            raw,
        } => {
            out.push_str(&format!("* int {charge} {multiplicity}\n"));
            if !raw.is_empty() {
                out.push_str(raw);
                out.push('\n');
            }
            out.push_str("*\n");
        }
    }
}

// ===========================================================================
// Output: errors
// ===========================================================================

/// Errors from [`parse_orca_output`].
#[derive(Debug, Clone, PartialEq)]
pub enum OrcaOutputError {
    /// Input exceeded the byte-size limit.
    InputTooLarge { limit: usize },
    /// Input exceeded the line-count limit.
    TooManyLines { limit: usize },
    /// A `FINAL SINGLE POINT ENERGY` value parsed but is NaN/Infinite.
    NonFiniteEnergy(String),
    /// A Cartesian-coordinate value in a geometry block parsed but is
    /// NaN/Infinite.
    NonFiniteCoordinate { line: usize, value: String },
    /// A vibrational-frequency value parsed but is NaN/Infinite.
    NonFiniteFrequency(String),
    /// A physical line exceeded the configured size limit.
    LineTooLong { actual: usize, limit: usize },
    /// The output contained more geometry frames than configured.
    TooManyGeometryFrames { actual: usize, limit: usize },
    /// A geometry frame contained more atoms than configured.
    TooManyGeometryAtoms { actual: usize, limit: usize },
    /// The output contained more vibrational frequencies than configured.
    TooManyFrequencies { actual: usize, limit: usize },
}

impl std::fmt::Display for OrcaOutputError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InputTooLarge { limit } => {
                write!(f, "input exceeds the {limit}-byte limit")
            }
            Self::TooManyLines { limit } => {
                write!(f, "input exceeds the {limit}-line limit")
            }
            Self::NonFiniteEnergy(v) => {
                write!(
                    f,
                    "FINAL SINGLE POINT ENERGY value '{v}' is not finite (NaN/Infinite)"
                )
            }
            Self::NonFiniteCoordinate { line, value } => {
                write!(
                    f,
                    "coordinate value '{value}' at line {line} is not finite (NaN/Infinite)"
                )
            }
            Self::NonFiniteFrequency(v) => {
                write!(f, "frequency value '{v}' is not finite (NaN/Infinite)")
            }
            Self::LineTooLong { actual, limit } => {
                write!(
                    f,
                    "line has size {actual}, exceeding the {limit}-byte limit"
                )
            }
            Self::TooManyGeometryFrames { actual, limit } => {
                write!(
                    f,
                    "geometry frames have size {actual}, exceeding the {limit} limit"
                )
            }
            Self::TooManyGeometryAtoms { actual, limit } => {
                write!(
                    f,
                    "geometry atoms have size {actual}, exceeding the {limit} limit"
                )
            }
            Self::TooManyFrequencies { actual, limit } => {
                write!(
                    f,
                    "frequencies have size {actual}, exceeding the {limit} limit"
                )
            }
        }
    }
}

impl std::error::Error for OrcaOutputError {}

// ===========================================================================
// Output: data model
// ===========================================================================

/// One geometry snapshot parsed from a `CARTESIAN COORDINATES (ANGSTROEM)`
/// block in an ORCA output file.
pub struct GeometryFrame {
    /// Molecular topology (atoms only; no bonds -- ORCA output carries no
    /// bond table).
    pub mol: Molecule,
    /// Atomic coordinates in Angstroms `(x, y, z)`, same order as `mol`'s
    /// atoms.
    pub coords: Vec<(f64, f64, f64)>,
}

/// How an ORCA job's process terminated. Never inferred from the presence
/// of a final geometry or energy -- always taken from ORCA's own explicit
/// termination marker (or the absence of one).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OrcaTermination {
    /// `****ORCA TERMINATED NORMALLY****` was found. Per the ORCA manual,
    /// this does **not** by itself mean a requested geometry optimization
    /// converged -- check [`OrcaOptConvergence`] separately.
    Normal,
    /// An `ORCA finished by error termination in <module>` line was found;
    /// the string is that line, verbatim.
    Error(String),
    /// Neither a normal nor an error termination marker was found before
    /// the input ended -- the output is truncated (crashed job, killed
    /// job, or a log file copied while the run was still in progress).
    Incomplete,
}

/// Geometry-optimization convergence status, reported explicitly rather
/// than inferred from whether a final geometry block exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrcaOptConvergence {
    /// No `GEOMETRY OPTIMIZATION CYCLE` marker was found -- this wasn't an
    /// optimization job (e.g. a single point or frequency-only run).
    NotRequested,
    /// The `...HURRAY...` / `THE OPTIMIZATION HAS CONVERGED` banner was
    /// found.
    Converged,
    /// ORCA's explicit "did not converge but reached the maximum number of
    /// optimization cycles" warning was found.
    NotConverged,
    /// Optimization cycles were seen but neither the converged nor the
    /// not-converged marker appeared before the input ended (e.g. a
    /// truncated/crashed output mid-optimization).
    Unknown,
}

/// Data extracted from an ORCA output file.
pub struct OrcaOutput {
    /// Total molecular charge, echoed near the top of the output.
    pub charge: Option<i32>,
    /// Spin multiplicity (2S+1), echoed near the top of the output.
    pub multiplicity: Option<u32>,
    /// The last plain (non `(MM)`/`(QM/MM)`-suffixed) `FINAL SINGLE POINT
    /// ENERGY` value found, in Hartree.
    pub final_energy_hartree: Option<f64>,
    /// Every `CARTESIAN COORDINATES (ANGSTROEM)` block found, in file
    /// order -- this is the geometry trajectory when an optimization ran
    /// (best-effort: one frame per block ORCA printed, which in practice
    /// is once per optimization cycle plus a final one), or a single frame
    /// for a plain single-point/frequency job. Empty if none were found.
    pub trajectory: Vec<GeometryFrame>,
    /// Vibrational frequencies in cm⁻¹, in mode order (including the
    /// near-zero translational/rotational modes ORCA lists first).
    /// Negative values are imaginary modes. Empty if no frequency
    /// calculation ran (or none was found).
    pub frequencies_cm1: Vec<f64>,
    /// How the job's process terminated.
    pub termination: OrcaTermination,
    /// Geometry-optimization convergence status.
    pub optimization_convergence: OrcaOptConvergence,
}

impl OrcaOutput {
    /// The final geometry -- the last frame of [`Self::trajectory`], if
    /// any was found.
    pub fn final_geometry(&self) -> Option<&GeometryFrame> {
        self.trajectory.last()
    }
}

// ===========================================================================
// Output: parser
// ===========================================================================

/// Byte-size limit for [`parse_orca_output`] -- a malformed or huge log
/// file returns a typed error rather than hanging or exhausting memory.
pub const MAX_OUTPUT_BYTES: usize = 64 << 20; // 64 MiB

/// Line-count limit for [`parse_orca_output`], checked independently of
/// the byte limit (guards against a file that is mostly short/empty lines,
/// which the byte limit alone wouldn't bound as tightly).
pub const MAX_OUTPUT_LINES: usize = 1_000_000;

/// Resource limits for ORCA output parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OrcaOutputParseLimits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_geometry_frames: usize,
    pub max_geometry_atoms: usize,
    pub max_frequencies: usize,
}

impl Default for OrcaOutputParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: MAX_OUTPUT_BYTES,
            max_line_bytes: 16 * 1024 * 1024,
            max_lines: MAX_OUTPUT_LINES,
            max_geometry_frames: 100_000,
            max_geometry_atoms: 1_000_000,
            max_frequencies: 1_000_000,
        }
    }
}

fn last_token(line: &str) -> Option<&str> {
    line.split_whitespace().next_back()
}

/// Parse a `"FINAL SINGLE POINT ENERGY <value>"` line's tokens (already
/// whitespace-split). Only matches the plain 4-word-prefix + 1-value shape
/// (5 tokens total) -- deliberately excludes QM/MM-suffixed variants like
/// `FINAL SINGLE POINT ENERGY (QM/MM) <value>`, which are out of scope
/// (see module docs). Returns `Ok(None)` for a non-matching shape, and a
/// typed error for a NaN/Infinite value.
fn parse_final_energy_line(toks: &[&str]) -> Result<Option<f64>, OrcaOutputError> {
    if toks.len() != 5 {
        return Ok(None);
    }
    match toks[4].parse::<f64>() {
        Ok(v) if v.is_finite() => Ok(Some(v)),
        Ok(_) => Err(OrcaOutputError::NonFiniteEnergy(toks[4].to_string())),
        Err(_) => Ok(None),
    }
}

/// Parse an ORCA output file (`.out` / `.log`).
///
/// Never panics on malformed or truncated input; always returns a typed
/// `Result`. Truncated input (a job that crashed or was killed mid-run) is
/// reported via [`OrcaTermination::Incomplete`], not an error.
/// Unrecognized fragments without an explicit ORCA error marker follow the
/// same incomplete-result contract; resource-limit violations remain typed
/// errors.
pub fn parse_orca_output(input: &str) -> Result<OrcaOutput, OrcaOutputError> {
    parse_orca_output_with_limits(input, &OrcaOutputParseLimits::default())
}

fn validated_output_lines<'a>(
    input: &'a str,
    limits: &OrcaOutputParseLimits,
) -> Result<Vec<&'a str>, OrcaOutputError> {
    if input.len() > limits.max_input_bytes {
        return Err(OrcaOutputError::InputTooLarge {
            limit: limits.max_input_bytes,
        });
    }
    let lines: Vec<&str> = input.lines().collect();
    if lines.len() > limits.max_lines {
        return Err(OrcaOutputError::TooManyLines {
            limit: limits.max_lines,
        });
    }
    if let Some(line) = lines.iter().find(|line| line.len() > limits.max_line_bytes) {
        return Err(OrcaOutputError::LineTooLong {
            actual: line.len(),
            limit: limits.max_line_bytes,
        });
    }
    Ok(lines)
}

/// Parse an ORCA output file with explicit resource limits.
pub fn parse_orca_output_with_limits(
    input: &str,
    limits: &OrcaOutputParseLimits,
) -> Result<OrcaOutput, OrcaOutputError> {
    let lines = validated_output_lines(input, limits)?;

    let mut charge = None;
    let mut multiplicity = None;
    let mut final_energy_hartree = None;
    let mut trajectory = Vec::new();
    let mut frequencies_cm1 = Vec::new();
    let mut termination = OrcaTermination::Incomplete;
    let mut saw_opt_cycle = false;
    let mut convergence = OrcaOptConvergence::NotRequested;

    let mut i = 0usize;
    while i < lines.len() {
        let t = lines[i].trim();

        if charge.is_none() && t.starts_with("Total Charge") {
            charge = last_token(t).and_then(|s| s.parse::<i32>().ok());
        } else if multiplicity.is_none() && t.starts_with("Multiplicity") {
            multiplicity = last_token(t).and_then(|s| s.parse::<u32>().ok());
        } else if t.contains("GEOMETRY OPTIMIZATION CYCLE") {
            saw_opt_cycle = true;
        } else if t == "CARTESIAN COORDINATES (ANGSTROEM)" {
            let mut j = i + 1;
            if j < lines.len() && is_dashes(lines[j]) {
                j += 1;
            }
            let (frame, next_j) = parse_geometry_table(&lines, j, limits.max_geometry_atoms)?;
            if !frame.coords.is_empty() {
                if trajectory.len() >= limits.max_geometry_frames {
                    return Err(OrcaOutputError::TooManyGeometryFrames {
                        actual: trajectory.len() + 1,
                        limit: limits.max_geometry_frames,
                    });
                }
                trajectory.push(frame);
            }
            i = next_j;
            continue;
        } else if t.starts_with("FINAL SINGLE POINT ENERGY") {
            let toks: Vec<&str> = t.split_whitespace().collect();
            if let Some(v) = parse_final_energy_line(&toks)? {
                final_energy_hartree = Some(v);
            }
        } else if t == "VIBRATIONAL FREQUENCIES" {
            let (freqs, next_i) = parse_frequency_table(&lines, i + 1, limits.max_frequencies)?;
            frequencies_cm1 = freqs;
            i = next_i;
            continue;
        } else if t.contains("THE OPTIMIZATION HAS CONVERGED") {
            convergence = OrcaOptConvergence::Converged;
        } else if t.contains("The optimization did not converge") {
            convergence = OrcaOptConvergence::NotConverged;
        } else if t.contains("ORCA TERMINATED NORMALLY") {
            termination = OrcaTermination::Normal;
        } else if t.contains("ORCA finished by error termination") {
            termination = OrcaTermination::Error(t.to_string());
        }

        i += 1;
    }

    if saw_opt_cycle && convergence == OrcaOptConvergence::NotRequested {
        convergence = OrcaOptConvergence::Unknown;
    }

    Ok(OrcaOutput {
        charge,
        multiplicity,
        final_energy_hartree,
        trajectory,
        frequencies_cm1,
        termination,
        optimization_convergence: convergence,
    })
}

fn is_dashes(line: &str) -> bool {
    let t = line.trim();
    !t.is_empty() && t.chars().all(|c| c == '-')
}

fn parse_geometry_table(
    lines: &[&str],
    start: usize,
    max_atoms: usize,
) -> Result<(GeometryFrame, usize), OrcaOutputError> {
    let mut builder = MoleculeBuilder::new();
    let mut coords = Vec::new();
    let mut i = start;
    while i < lines.len() {
        let t = lines[i].trim();
        if t.is_empty() {
            break;
        }
        let toks: Vec<&str> = t.split_whitespace().collect();
        if toks.len() != 4 {
            break;
        }
        let Some(element) = Element::from_symbol(toks[0]) else {
            break;
        };
        let mut parsed = [0.0f64; 3];
        let mut row_is_numeric = true;
        for (axis, tok) in toks[1..4].iter().enumerate() {
            match tok.parse::<f64>() {
                Ok(v) if v.is_finite() => parsed[axis] = v,
                Ok(_) => {
                    return Err(OrcaOutputError::NonFiniteCoordinate {
                        line: i + 1,
                        value: tok.to_string(),
                    });
                }
                Err(_) => {
                    row_is_numeric = false;
                    break;
                }
            }
        }
        if !row_is_numeric {
            break;
        }
        if coords.len() >= max_atoms {
            return Err(OrcaOutputError::TooManyGeometryAtoms {
                actual: coords.len() + 1,
                limit: max_atoms,
            });
        }
        builder.add_atom(Atom::new(element));
        coords.push((parsed[0], parsed[1], parsed[2]));
        i += 1;
    }
    Ok((
        GeometryFrame {
            mol: builder.build(),
            coords,
        },
        i,
    ))
}

/// Parse one `"N:   <value> cm**-1"` frequency-table entry (optionally
/// followed by trailing text such as `***imaginary mode***`, which is
/// ignored). Returns `None` if the line doesn't match this shape at all,
/// `Some(Err(value))` if it matches but the value is NaN/Infinite.
fn parse_frequency_line(line: &str) -> Option<Result<f64, String>> {
    let (index_part, rest) = line.split_once(':')?;
    index_part.trim().parse::<u32>().ok()?;
    let mut toks = rest.split_whitespace();
    let value_tok = toks.next()?;
    let unit_tok = toks.next()?;
    if !unit_tok.eq_ignore_ascii_case("cm**-1") {
        return None;
    }
    match value_tok.parse::<f64>() {
        Ok(v) if v.is_finite() => Some(Ok(v)),
        Ok(_) => Some(Err(value_tok.to_string())),
        Err(_) => None,
    }
}

fn parse_frequency_table(
    lines: &[&str],
    start: usize,
    max_frequencies: usize,
) -> Result<(Vec<f64>, usize), OrcaOutputError> {
    let mut i = start;
    if i < lines.len() && is_dashes(lines[i]) {
        i += 1;
    }
    let mut freqs = Vec::new();
    // Bound how many non-matching lines (blank lines, an optional "Mode
    // freq (cm**-1)" sub-header and its own dashed underline) we skip
    // before giving up on ever finding an entry, so a spurious header
    // match can't silently swallow the rest of a large file.
    let mut skip_budget = 10;
    loop {
        if i >= lines.len() {
            break;
        }
        let t = lines[i].trim();
        if t.is_empty() {
            i += 1;
            if !freqs.is_empty() {
                break;
            }
            continue;
        }
        match parse_frequency_line(t) {
            Some(Ok(v)) => {
                if freqs.len() >= max_frequencies {
                    return Err(OrcaOutputError::TooManyFrequencies {
                        actual: freqs.len() + 1,
                        limit: max_frequencies,
                    });
                }
                freqs.push(v);
                i += 1;
            }
            Some(Err(bad)) => return Err(OrcaOutputError::NonFiniteFrequency(bad)),
            None => {
                if freqs.is_empty() && skip_budget > 0 {
                    skip_budget -= 1;
                    i += 1;
                    continue;
                }
                break;
            }
        }
    }
    Ok((freqs, i))
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
#[path = "orca/tests.rs"]
mod tests;
