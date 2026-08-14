//! POSCAR / CONTCAR (VASP structure file format) read/write.
//!
//! POSCAR is VASP's plain-text periodic structure format: a comment line, a
//! scale factor, three lattice-vector lines, a species-name line, a
//! per-species atom-count line, an optional `Selective dynamics` line, a
//! coordinate-mode line (`Direct`/`Cartesian`), then one line per atom.
//! CONTCAR (VASP's own *output* structure file, e.g. after a relaxation or
//! an MD run) is byte-for-byte the same format -- [`parse_poscar`] reads
//! both; [`parse_contcar`] is a discoverability-only alias.
//!
//! ```text
//! NaCl
//! 1.0
//!    5.6400000000000000    0.0000000000000000    0.0000000000000000
//!    0.0000000000000000    5.6400000000000000    0.0000000000000000
//!    0.0000000000000000    0.0000000000000000    5.6400000000000000
//! Na Cl
//! 1 1
//! Direct
//!   0.0 0.0 0.0
//!   0.5 0.5 0.5
//! ```
//!
//! # Scope decisions (read primary VASP docs / a reference implementation
//! before changing any of these)
//!
//! - **VASP 5 only.** An explicit species-name line before the atom-count
//!   line is required. A file that omits it (VASP 4's implicit
//!   POTCAR-derived species order) is rejected with
//!   [`PoscarError::Vasp4NotSupported`] rather than silently mis-parsed --
//!   detected the same way a real implementation
//!   ([pymatgen's `Poscar.from_str`](https://github.com/materialsproject/pymatgen/blob/v2024.1.27/pymatgen/io/vasp/inputs.py#L279))
//!   does: if every token on the species-name line parses as an integer,
//!   that line is actually the atom-count line (VASP 4 layout).
//! - **Species/count header is single-line only.** Real VASP output wraps
//!   these onto multiple lines only past ~20 distinct species groups in one
//!   cell -- a case this reader does not special-case (a wrapped header
//!   fails the atom-count line's own length check with a clear error,
//!   rather than being silently misread).
//! - **Scale factor**: both forms documented on the
//!   [VASP wiki POSCAR page](https://www.vasp.at/wiki/index.php/POSCAR) are
//!   supported on read: a single number (negative meaning "target cell
//!   volume in cubic Angstrom", the linear factor derived via
//!   `(|volume| / |det(raw matrix)|).cbrt()`), and three numbers ("individual
//!   scaling factors for the x-, y- and z-Cartesian components of the
//!   lattice vectors (and 'Cartesian' mode ion positions)", per that page --
//!   applied here as a per-column multiplier on both the raw lattice matrix
//!   and any Cartesian-mode ion coordinates, all three components required
//!   positive). [`write_poscar`] always emits the single-value form `1.0`
//!   with pre-scaled vectors -- the simplest form that is always exactly
//!   correct, and what most real POSCAR writers do.
//! - **Direct/Cartesian mode + Selective dynamics are matched the way VASP
//!   itself matches them**: only the *first character* is significant,
//!   case-insensitively (`C`/`c`/`K`/`k` -> Cartesian, anything else ->
//!   Direct; `S`/`s` -> Selective dynamics on) -- confirmed against the VASP
//!   wiki POSCAR page's own wording, not guessed. Per-atom selective-
//!   dynamics flags use the same first-character rule (`T`/`t` -> `true`);
//!   this one place does **not** fail closed on a malformed flag (e.g. `X`
//!   or an unexpected token) -- it silently reads as `false`, deliberately
//!   matching a real implementation's own leniency here (pymatgen's
//!   `tok.upper()[0] == "T"`) rather than rejecting files real tools accept.
//!   Every other field in this reader (coordinates, scale factor, atom
//!   counts, ...) does fail closed on a malformed value.
//! - **Write always emits `Direct` (fractional) coordinates.**
//!   [`crate::structure::PeriodicStructure`] stores fractional coordinates
//!   canonically; a Cartesian-mode POSCAR is converted to fractional on
//!   read (via [`Lattice::cart_to_frac`]) like every other consumer of this
//!   crate's coordinates, and is not round-tripped back out as Cartesian
//!   text.
//! - **No disorder / partial occupancy.** POSCAR/CONTCAR has no concept of
//!   either -- every atom is one fully-occupied species. [`parse_poscar`]
//!   always builds single-species, occupancy-1.0
//!   [`SiteSpecies`](crate::site::SiteSpecies)s; [`write_poscar`] rejects
//!   (via [`PoscarError::Unwritable`]) any site with more than one species
//!   or an occupancy other than 1.0 rather than silently dropping the
//!   disorder data.
//! - **No per-atom labels.** VASP has no per-atom name field beyond the
//!   species symbol (unlike a CIF's `Na1`/`Cl1` atom-site labels) --
//!   [`parse_poscar`] leaves every [`PeriodicSite::label`](crate::site::PeriodicSite::label)
//!   `None` rather than fabricating one.
//! - **Ion velocities**: read as a blank-line-separated block of exactly
//!   one `vx vy vz` line per atom, in "direct lattice vector / timestep"
//!   units, with **no** coordinate-mode line of its own -- this is what a
//!   real implementation
//!   ([pymatgen, same file, `if len(chunks) > 1`](https://github.com/materialsproject/pymatgen/blob/v2024.1.27/pymatgen/io/vasp/inputs.py#L445))
//!   assumes, and it fails closed (a typed parse error, not silent
//!   corruption) if a line doesn't parse as three finite numbers.
//! - **Predictor-corrector data** (`PredictorCorrector`, CONTCAR-only,
//!   MD-restart bookkeeping): VASP's own wiki explicitly does not give this
//!   section's numeric field layout ("cannot be entered by hand"), so this
//!   reader does not attempt to interpret it -- it stores the 3-line
//!   preamble and every remaining line **verbatim** (matching the section
//!   shape a real implementation
//!   ([pymatgen, same file, `if len(chunks) > 2`](https://github.com/materialsproject/pymatgen/blob/v2024.1.27/pymatgen/io/vasp/inputs.py#L451))
//!   uses: 3 preamble lines then `3 * atom_count` data lines) and replays
//!   it byte-for-byte on write. Because the data is opaque, [`write_poscar`]
//!   refuses (rather than silently reordering it wrongly) to write it
//!   together with a structure whose sites need reordering into contiguous
//!   species groups.
//! - **Lattice velocities** (an NPT/variable-cell MD-restart extension,
//!   detected by an initial `Lattice velocities and vectors` header
//!   appearing directly after the coordinate block with no blank line
//!   separator, again per pymatgen's parser) are out of scope: detected and
//!   rejected with [`PoscarError::UnsupportedSection`] rather than
//!   misread as an ion-velocities block.
//! - No VASP INCAR/KPOINTS/POTCAR parsing, no symmetry/space-group
//!   detection, no DFT -- pure structure I/O, matching this crate's own
//!   non-goals (see the crate-level docs).

use chematic_core::Element;

use crate::error::CrystalError;
use crate::lattice::{Lattice, cross3, dot3};
use crate::site::{CartesianCoord, FractionalCoord, Occupancy, PeriodicSite, SiteSpecies};
use crate::structure::PeriodicStructure;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// CONTCAR-only MD-restart "predictor-corrector" data trailing the ion
/// velocities. Opaque by design -- see the module docs for why this reader
/// does not interpret its numeric fields, only preserves them verbatim.
#[derive(Debug, Clone, PartialEq)]
pub struct PredictorCorrector {
    /// The section's first three lines verbatim: an initialization-state
    /// key, the MD timestep (`POTIM`), and the Nose-Hoover thermostat
    /// values for iterations `n`/`n+1`.
    pub preamble: [String; 3],
    /// Every remaining line, verbatim, in file order (real files hold `3 *
    /// atom_count` lines here, per a real implementation's parsing -- see
    /// module docs -- but this type does not assume or enforce that count).
    pub lines: Vec<String>,
}

/// A parsed POSCAR/CONTCAR document: the structural data
/// [`PeriodicStructure`] can represent, plus the format's other fields that
/// it cannot (comment, original species-group order, selective dynamics,
/// velocities, predictor-corrector).
#[derive(Debug, Clone, PartialEq)]
pub struct PoscarDocument {
    /// The parsed lattice + sites (see module docs: always fractional,
    /// always single-species/full-occupancy).
    pub structure: PeriodicStructure,
    /// The file's first line, verbatim (POSCAR's free-form comment/title).
    pub comment: String,
    /// The species-name line's symbols, in file order, one entry per
    /// contiguous atom group (may repeat an element if the file itself
    /// declared two separate groups of the same species). **Read-fidelity
    /// only** -- [`write_poscar`] does not consult this field; it re-derives
    /// species grouping from `structure.sites()`'s own order (first
    /// appearance = group order), per the format's own requirement that
    /// atoms of one species are contiguous. Kept here purely so a reader
    /// can inspect the file's original grouping without re-deriving it.
    pub species_order: Vec<Element>,
    /// Per-atom `[x, y, z]` selective-dynamics flags (`true` = allowed to
    /// move), same order as `structure.sites()`, if the file had a
    /// `Selective dynamics` line.
    pub selective_dynamics: Option<Vec<[bool; 3]>>,
    /// Per-atom `[vx, vy, vz]` ion velocities (direct lattice vector /
    /// timestep units), same order as `structure.sites()`, if present.
    pub velocities: Option<Vec<[f64; 3]>>,
    /// CONTCAR MD-restart predictor-corrector data, if present. Only
    /// possible when `velocities` is also `Some` (matches the file layout:
    /// this section always follows ion velocities).
    pub predictor_corrector: Option<PredictorCorrector>,
}

/// Errors from [`parse_poscar`]/[`parse_contcar`]/[`write_poscar`].
#[derive(Debug, Clone, PartialEq)]
pub enum PoscarError {
    /// The input has no lines at all.
    Empty,
    /// A required line is missing (input ended early).
    MissingLine {
        /// 1-indexed line number that should have been present.
        line: usize,
        /// What was expected there.
        detail: &'static str,
    },
    /// A line's content doesn't match the format's grammar at that
    /// position.
    InvalidLine {
        /// 1-indexed line number.
        line: usize,
        /// What was wrong.
        detail: String,
    },
    /// The species-name line's tokens all parsed as integers, meaning this
    /// is a VASP 4-style file (implicit POTCAR-derived species order) --
    /// out of scope, see module docs.
    Vasp4NotSupported,
    /// A species symbol is not a recognized element.
    UnknownElement {
        /// The unrecognized symbol, verbatim.
        symbol: String,
        /// 1-indexed line number.
        line: usize,
    },
    /// A numeric field parsed but is `NaN`/`Infinity`.
    NonFiniteValue {
        /// 1-indexed line number.
        line: usize,
        /// Which field.
        detail: String,
    },
    /// The atom-count line declared more atoms than coordinate lines were
    /// present for.
    AtomCountMismatch {
        /// Declared atom count (sum of the atom-count line's entries).
        declared: usize,
        /// Coordinate lines actually found before EOF/next section.
        found: usize,
    },
    /// Trailing content after the coordinate block did not match any
    /// section this reader understands (e.g. an NPT lattice-velocities
    /// block, or content following predictor-corrector data).
    UnsupportedSection {
        /// What was found and why it's unsupported.
        detail: String,
    },
    /// The parsed lattice or a site failed [`crate::structure`]'s own
    /// validation (e.g. a singular lattice matrix, an out-of-range
    /// component).
    Structure(CrystalError),
    /// [`write_poscar`] cannot represent `doc` as POSCAR text (disorder,
    /// partial occupancy, or opaque predictor-corrector data that would
    /// need reordering).
    Unwritable {
        /// What was wrong.
        detail: String,
    },
}

impl From<CrystalError> for PoscarError {
    fn from(e: CrystalError) -> Self {
        PoscarError::Structure(e)
    }
}

impl std::fmt::Display for PoscarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PoscarError::Empty => write!(f, "empty POSCAR/CONTCAR input"),
            PoscarError::MissingLine { line, detail } => {
                write!(f, "missing line {line}: expected {detail}")
            }
            PoscarError::InvalidLine { line, detail } => {
                write!(f, "invalid line {line}: {detail}")
            }
            PoscarError::Vasp4NotSupported => write!(
                f,
                "VASP 4-style POSCAR (species names implicit from POTCAR order, no species-name line) is not supported -- add an explicit species-name line (VASP 5 format)"
            ),
            PoscarError::UnknownElement { symbol, line } => {
                write!(f, "unknown element symbol '{symbol}' at line {line}")
            }
            PoscarError::NonFiniteValue { line, detail } => {
                write!(f, "non-finite {detail} at line {line}")
            }
            PoscarError::AtomCountMismatch { declared, found } => write!(
                f,
                "atom-count line declared {declared} atom(s), found {found} coordinate line(s)"
            ),
            PoscarError::UnsupportedSection { detail } => {
                write!(f, "unsupported section: {detail}")
            }
            PoscarError::Structure(e) => write!(f, "{e}"),
            PoscarError::Unwritable { detail } => write!(f, "cannot write POSCAR: {detail}"),
        }
    }
}

impl std::error::Error for PoscarError {}

// ---------------------------------------------------------------------------
// Read
// ---------------------------------------------------------------------------

/// Parse a POSCAR (or CONTCAR) file. See the module docs for exactly which
/// format variants are supported.
pub fn parse_poscar(text: &str) -> Result<PoscarDocument, PoscarError> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.is_empty() {
        return Err(PoscarError::Empty);
    }

    // Positional header parsing (no blank-line stripping): line 1 (the
    // comment) is legitimately free-form and may itself be blank, so
    // treating blank lines as section separators here (as this module does
    // for the *trailing* velocities/predictor-corrector sections) would
    // silently shift every following field by one line.
    let mut i = 0usize; // 0-indexed cursor into `lines`

    let comment = lines[i].to_string();
    i += 1;

    let scale_line = *lines.get(i).ok_or(PoscarError::MissingLine {
        line: i + 1,
        detail: "scale factor",
    })?;
    let scale_line_no = i + 1;
    i += 1;
    let scale_tokens: Vec<f64> = scale_line
        .split_whitespace()
        .map(|t| {
            t.parse::<f64>().map_err(|_| PoscarError::InvalidLine {
                line: scale_line_no,
                detail: format!("invalid scale factor token '{t}'"),
            })
        })
        .collect::<Result<_, _>>()?;
    if scale_tokens.is_empty() || scale_tokens.len() == 2 || scale_tokens.len() > 3 {
        return Err(PoscarError::InvalidLine {
            line: scale_line_no,
            detail: format!(
                "scale factor line must have 1 or 3 number(s), got {}",
                scale_tokens.len()
            ),
        });
    }
    if scale_tokens.iter().any(|s| !s.is_finite()) {
        return Err(PoscarError::NonFiniteValue {
            line: scale_line_no,
            detail: "scale factor".to_string(),
        });
    }
    if scale_tokens.len() == 3 && scale_tokens.iter().any(|s| *s <= 0.0) {
        return Err(PoscarError::InvalidLine {
            line: scale_line_no,
            detail: "a 3-component scale factor's entries must all be positive".to_string(),
        });
    }

    let mut raw_matrix = [[0.0f64; 3]; 3];
    for row in raw_matrix.iter_mut() {
        let line = *lines.get(i).ok_or(PoscarError::MissingLine {
            line: i + 1,
            detail: "lattice vector",
        })?;
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 3 {
            return Err(PoscarError::InvalidLine {
                line: i + 1,
                detail: format!("lattice vector needs 3 numbers, got '{line}'"),
            });
        }
        for (c, tok) in toks[..3].iter().enumerate() {
            let v: f64 = tok.parse().map_err(|_| PoscarError::InvalidLine {
                line: i + 1,
                detail: format!("cannot parse lattice component '{tok}'"),
            })?;
            if !v.is_finite() {
                return Err(PoscarError::NonFiniteValue {
                    line: i + 1,
                    detail: "lattice vector component".to_string(),
                });
            }
            row[c] = v;
        }
        i += 1;
    }

    // Per-column scale factor: uniform (all 3 = the single value, with the
    // negative/volume-target case pre-resolved to a linear factor) or the
    // 3-component per-axis form. Applying the *same* column-scale array to
    // both the lattice matrix and any Cartesian-mode ion coordinates below
    // unifies the single-value and 3-component cases into one code path,
    // and is what the VASP wiki's own wording for the 3-component form
    // describes (it names both "the lattice vectors" and "'Cartesian' mode
    // ion positions" as scaled by the same per-axis factors).
    let col_scale: [f64; 3] = if scale_tokens.len() == 3 {
        [scale_tokens[0], scale_tokens[1], scale_tokens[2]]
    } else {
        let s = scale_tokens[0];
        let factor = if s < 0.0 {
            let raw_det = dot3(raw_matrix[0], cross3(raw_matrix[1], raw_matrix[2]));
            (s.abs() / raw_det.abs()).cbrt()
        } else {
            s
        };
        [factor, factor, factor]
    };
    let mut matrix = raw_matrix;
    for row in matrix.iter_mut() {
        for c in 0..3 {
            row[c] *= col_scale[c];
        }
    }
    let lattice = Lattice::from_matrix(matrix)?;

    // Species-name line vs. VASP 4's implicit-order atom-count line: if
    // every token here parses as an integer, this is actually the
    // atom-count line (no species-name line present) -- see module docs.
    let header_line = *lines.get(i).ok_or(PoscarError::MissingLine {
        line: i + 1,
        detail: "species-name line",
    })?;
    let header_tokens: Vec<&str> = header_line.split_whitespace().collect();
    if header_tokens.is_empty() {
        return Err(PoscarError::InvalidLine {
            line: i + 1,
            detail: "species-name line is empty".to_string(),
        });
    }
    if header_tokens.iter().all(|t| t.parse::<u64>().is_ok()) {
        return Err(PoscarError::Vasp4NotSupported);
    }
    let species_line_no = i + 1;
    let species_order: Vec<Element> = header_tokens
        .iter()
        .map(|t| {
            Element::from_symbol(t).ok_or_else(|| PoscarError::UnknownElement {
                symbol: t.to_string(),
                line: species_line_no,
            })
        })
        .collect::<Result<_, _>>()?;
    i += 1;

    let counts_line = *lines.get(i).ok_or(PoscarError::MissingLine {
        line: i + 1,
        detail: "atom-count line",
    })?;
    let counts_line_no = i + 1;
    let counts: Vec<usize> = counts_line
        .split_whitespace()
        .map(|t| {
            t.parse::<usize>().map_err(|_| PoscarError::InvalidLine {
                line: counts_line_no,
                detail: format!("invalid atom count '{t}'"),
            })
        })
        .collect::<Result<_, _>>()?;
    i += 1;
    if counts.len() != species_order.len() {
        return Err(PoscarError::InvalidLine {
            line: counts_line_no,
            detail: format!(
                "species-name line has {} symbol(s) but atom-count line has {} count(s)",
                species_order.len(),
                counts.len()
            ),
        });
    }

    let n_sites: usize = counts.iter().sum();
    let mut atom_elements: Vec<Element> = Vec::with_capacity(n_sites);
    for (&count, &elem) in counts.iter().zip(species_order.iter()) {
        atom_elements.extend(std::iter::repeat_n(elem, count));
    }

    // Optional "Selective dynamics" line, then the coordinate-mode line.
    // Only the first character is significant, case-insensitively, exactly
    // matching VASP's own parsing rule (see module docs).
    let mut has_selective = false;
    let mut mode_line = *lines.get(i).ok_or(PoscarError::MissingLine {
        line: i + 1,
        detail: "coordinate-mode line",
    })?;
    if first_char_matches(mode_line, 's') {
        has_selective = true;
        i += 1;
        mode_line = *lines.get(i).ok_or(PoscarError::MissingLine {
            line: i + 1,
            detail: "coordinate-mode line",
        })?;
    }
    let cartesian = first_char_matches(mode_line, 'c') || first_char_matches(mode_line, 'k');
    i += 1;

    let mut fractional: Vec<[f64; 3]> = Vec::with_capacity(n_sites);
    let mut selective_flags: Vec<[bool; 3]> = Vec::with_capacity(n_sites);
    for _ in 0..n_sites {
        let Some(&line) = lines.get(i) else {
            return Err(PoscarError::AtomCountMismatch {
                declared: n_sites,
                found: fractional.len(),
            });
        };
        // A blank line here means the coordinate block ended early and what
        // follows is (or was meant to be) a blank-line-separated trailing
        // section (velocities/predictor-corrector) -- report the real
        // problem (too few coordinate lines) rather than the confusing
        // "coordinate line needs N fields, got ''" the empty-token check
        // below would otherwise produce.
        if line.trim().is_empty() {
            return Err(PoscarError::AtomCountMismatch {
                declared: n_sites,
                found: fractional.len(),
            });
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        let need = if has_selective { 6 } else { 3 };
        if toks.len() < need {
            return Err(PoscarError::InvalidLine {
                line: i + 1,
                detail: format!("atom coordinate line needs {need} field(s), got '{line}'"),
            });
        }
        let mut xyz = [0.0f64; 3];
        for (c, tok) in toks[..3].iter().enumerate() {
            let v: f64 = tok.parse().map_err(|_| PoscarError::InvalidLine {
                line: i + 1,
                detail: format!("cannot parse coordinate '{tok}'"),
            })?;
            if !v.is_finite() {
                return Err(PoscarError::NonFiniteValue {
                    line: i + 1,
                    detail: "atom coordinate".to_string(),
                });
            }
            xyz[c] = v;
        }
        let frac = if cartesian {
            let cart = [
                xyz[0] * col_scale[0],
                xyz[1] * col_scale[1],
                xyz[2] * col_scale[2],
            ];
            lattice.cart_to_frac(CartesianCoord::new(cart)).0
        } else {
            xyz
        };
        fractional.push(frac);
        if has_selective {
            let mut flags = [false; 3];
            for (c, tok) in toks[3..6].iter().enumerate() {
                flags[c] = first_char_matches(tok, 't');
            }
            selective_flags.push(flags);
        }
        i += 1;
    }

    let mut sites = Vec::with_capacity(n_sites);
    for (&elem, &frac) in atom_elements.iter().zip(fractional.iter()) {
        sites.push(PeriodicSite::new(
            vec![SiteSpecies::full(elem)],
            FractionalCoord::new(frac),
            None,
        )?);
    }
    let structure = PeriodicStructure::new(lattice, sites)?;

    // Trailing sections: blank-line-separated, unlike the positionally-
    // parsed header above -- this is where blank lines genuinely are
    // section separators.
    let (velocities, predictor_corrector) = parse_trailing_sections(&lines[i..], i, n_sites)?;

    Ok(PoscarDocument {
        structure,
        comment,
        species_order,
        selective_dynamics: has_selective.then_some(selective_flags),
        velocities,
        predictor_corrector,
    })
}

/// Parse a CONTCAR file. CONTCAR is byte-for-byte the same format as
/// POSCAR -- VASP's own convention is simply that CONTCAR is the *output*
/// structure file from a run (potentially with velocities/predictor-
/// corrector data appended). This is a thin, discoverability-only alias
/// for [`parse_poscar`].
pub fn parse_contcar(text: &str) -> Result<PoscarDocument, PoscarError> {
    parse_poscar(text)
}

/// `true` if `s`'s first character matches `expected`, case-insensitively
/// (VASP's own rule for `Direct`/`Cartesian`/`Selective dynamics` -- see
/// module docs). `false` for an empty/all-whitespace string.
fn first_char_matches(s: &str, expected: char) -> bool {
    s.trim_start()
        .chars()
        .next()
        .is_some_and(|c| c.eq_ignore_ascii_case(&expected))
}

/// Parse the optional ion-velocities and predictor-corrector sections
/// following the coordinate block. `rest` is every remaining line;
/// `start_line` is `rest[0]`'s 1-indexed line number (for error messages).
#[allow(clippy::type_complexity)]
fn parse_trailing_sections(
    rest: &[&str],
    start_line: usize,
    n_sites: usize,
) -> Result<(Option<Vec<[f64; 3]>>, Option<PredictorCorrector>), PoscarError> {
    let chunks = split_into_chunks(rest, start_line);
    let mut chunk_iter = chunks.into_iter();

    let Some(vel_chunk) = chunk_iter.next() else {
        return Ok((None, None));
    };

    // A "Lattice velocities and vectors" block (NPT/variable-cell MD
    // restart, out of scope -- see module docs) attaches directly after
    // the coordinate block with no blank-line separator, so it would land
    // right here as the first trailing chunk. Detect and reject it
    // explicitly rather than let it silently masquerade as ion velocities.
    if let Some(&(_, first_line)) = vel_chunk.first()
        && first_char_matches(first_line, 'l')
    {
        return Err(PoscarError::UnsupportedSection {
            detail: "lattice velocities / variable-cell MD restart data (an NPT CONTCAR extension) is not supported".to_string(),
        });
    }
    if vel_chunk.len() != n_sites {
        return Err(PoscarError::UnsupportedSection {
            detail: format!(
                "expected {n_sites} ion-velocity line(s) after the coordinate block, found {} non-blank line(s) -- unrecognized trailing section",
                vel_chunk.len()
            ),
        });
    }
    let mut velocities = Vec::with_capacity(n_sites);
    for (line_no, line) in &vel_chunk {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 3 {
            return Err(PoscarError::InvalidLine {
                line: *line_no,
                detail: format!("velocity line needs 3 numbers, got '{line}'"),
            });
        }
        let mut v = [0.0f64; 3];
        for (c, tok) in toks[..3].iter().enumerate() {
            let val: f64 = tok.parse().map_err(|_| PoscarError::InvalidLine {
                line: *line_no,
                detail: format!("cannot parse velocity component '{tok}'"),
            })?;
            if !val.is_finite() {
                return Err(PoscarError::NonFiniteValue {
                    line: *line_no,
                    detail: "velocity component".to_string(),
                });
            }
            v[c] = val;
        }
        velocities.push(v);
    }

    let predictor_corrector = match chunk_iter.next() {
        None => None,
        Some(pc_chunk) => {
            if pc_chunk.len() < 3 {
                return Err(PoscarError::InvalidLine {
                    line: pc_chunk.first().map_or(start_line, |(n, _)| *n),
                    detail: "predictor-corrector section needs at least 3 preamble lines"
                        .to_string(),
                });
            }
            let preamble = [
                pc_chunk[0].1.to_string(),
                pc_chunk[1].1.to_string(),
                pc_chunk[2].1.to_string(),
            ];
            let lines = pc_chunk[3..].iter().map(|(_, l)| l.to_string()).collect();
            Some(PredictorCorrector { preamble, lines })
        }
    };

    if let Some(extra) = chunk_iter.next() {
        return Err(PoscarError::UnsupportedSection {
            detail: format!(
                "unrecognized content at line {} after the predictor-corrector section",
                extra.first().map_or(start_line, |(n, _)| *n)
            ),
        });
    }

    Ok((Some(velocities), predictor_corrector))
}

/// Group `lines` (already offset so `lines[k]` is absolute line number
/// `start_line + k + 1`) into blank-line-separated chunks, pairing each
/// kept line with its 1-indexed line number.
fn split_into_chunks<'a>(lines: &[&'a str], start_line: usize) -> Vec<Vec<(usize, &'a str)>> {
    let mut chunks = Vec::new();
    let mut current: Vec<(usize, &str)> = Vec::new();
    for (offset, &line) in lines.iter().enumerate() {
        if line.trim().is_empty() {
            if !current.is_empty() {
                chunks.push(std::mem::take(&mut current));
            }
        } else {
            current.push((start_line + offset + 1, line));
        }
    }
    if !current.is_empty() {
        chunks.push(current);
    }
    chunks
}

// ---------------------------------------------------------------------------
// Write
// ---------------------------------------------------------------------------

/// Write a [`PoscarDocument`] as POSCAR text. Always emits scale factor
/// `1.0` with pre-scaled lattice vectors and `Direct` (fractional)
/// coordinates -- see the module docs for why. Species groups are written
/// in first-appearance order among `doc.structure.sites()` (stable: a
/// species' sites keep their relative order), **not** `doc.species_order`
/// (a read-fidelity record the writer does not consult) and not
/// alphabetical order.
///
/// Returns [`PoscarError::Unwritable`] if any site has more than one
/// species or an occupancy other than `1.0` (POSCAR cannot represent
/// disorder/partial occupancy), or if writing would require reordering
/// sites into contiguous species groups while `doc.predictor_corrector` is
/// present (its data is opaque verbatim text -- see module docs -- so it
/// cannot be safely permuted).
pub fn write_poscar(doc: &PoscarDocument) -> Result<String, PoscarError> {
    let sites = doc.structure.sites();
    let n = sites.len();

    let mut elements = Vec::with_capacity(n);
    for (idx, site) in sites.iter().enumerate() {
        if site.species.len() != 1 {
            return Err(PoscarError::Unwritable {
                detail: format!(
                    "site {idx} has {} species; POSCAR has no disorder/partial-occupancy representation (every site must be exactly one species)",
                    site.species.len()
                ),
            });
        }
        let sp = &site.species[0];
        if (sp.occupancy.value() - 1.0).abs() > Occupancy::SUM_TOLERANCE {
            return Err(PoscarError::Unwritable {
                detail: format!(
                    "site {idx} has occupancy {} (POSCAR requires full occupancy 1.0)",
                    sp.occupancy.value()
                ),
            });
        }
        elements.push(sp.element);
    }

    // Stable group order = first-appearance order of each species among
    // the sites (not alphabetical) -- required by the format itself
    // (atoms of one species must be contiguous).
    let mut group_order: Vec<Element> = Vec::new();
    for &e in &elements {
        if !group_order.contains(&e) {
            group_order.push(e);
        }
    }
    let group_index = |e: Element| {
        group_order
            .iter()
            .position(|&g| g == e)
            .expect("built from elements above")
    };

    let mut perm: Vec<usize> = (0..n).collect();
    perm.sort_by_key(|&idx| group_index(elements[idx])); // stable: preserves relative order within a group

    let is_identity = perm.iter().enumerate().all(|(i, &p)| i == p);
    if !is_identity && doc.predictor_corrector.is_some() {
        return Err(PoscarError::Unwritable {
            detail: "sites need reordering into contiguous species groups, but predictor_corrector holds opaque verbatim MD-restart data that cannot be safely permuted -- reorder sites (or clear predictor_corrector) before writing".to_string(),
        });
    }

    let mut counts = vec![0usize; group_order.len()];
    for &e in &elements {
        counts[group_index(e)] += 1;
    }

    if let Some(vel) = &doc.velocities
        && vel.len() != n
    {
        return Err(PoscarError::Unwritable {
            detail: format!(
                "velocities has {} entries but structure has {n} site(s)",
                vel.len()
            ),
        });
    }
    if let Some(sel) = &doc.selective_dynamics
        && sel.len() != n
    {
        return Err(PoscarError::Unwritable {
            detail: format!(
                "selective_dynamics has {} entries but structure has {n} site(s)",
                sel.len()
            ),
        });
    }
    if doc.predictor_corrector.is_some() && doc.velocities.is_none() {
        return Err(PoscarError::Unwritable {
            detail: "predictor_corrector is present but velocities is None -- VASP's CONTCAR layout always has ion velocities before predictor-corrector data".to_string(),
        });
    }

    let mut out = String::new();
    out.push_str(&doc.comment);
    out.push('\n');
    out.push_str("1.0\n");
    for row in doc.structure.lattice().matrix() {
        out.push_str(&format!("{} {} {}\n", row[0], row[1], row[2]));
    }
    out.push_str(
        &group_order
            .iter()
            .map(|e| e.symbol())
            .collect::<Vec<_>>()
            .join(" "),
    );
    out.push('\n');
    out.push_str(
        &counts
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(" "),
    );
    out.push('\n');
    if doc.selective_dynamics.is_some() {
        out.push_str("Selective dynamics\n");
    }
    out.push_str("Direct\n");

    for &orig_idx in &perm {
        let f = sites[orig_idx].fractional.0;
        out.push_str(&format!("{} {} {}", f[0], f[1], f[2]));
        if let Some(sel) = &doc.selective_dynamics {
            for flag in sel[orig_idx] {
                out.push_str(if flag { " T" } else { " F" });
            }
        }
        out.push('\n');
    }

    if let Some(vel) = &doc.velocities {
        out.push('\n');
        for &orig_idx in &perm {
            let v = vel[orig_idx];
            out.push_str(&format!("{} {} {}\n", v[0], v[1], v[2]));
        }
        if let Some(pc) = &doc.predictor_corrector {
            out.push('\n');
            for line in &pc.preamble {
                out.push_str(line);
                out.push('\n');
            }
            for line in &pc.lines {
                out.push_str(line);
                out.push('\n');
            }
        }
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const NACL_POSCAR: &str = "NaCl\n\
1.0\n\
   5.6400000000000000    0.0000000000000000    0.0000000000000000\n\
   0.0000000000000000    5.6400000000000000    0.0000000000000000\n\
   0.0000000000000000    0.0000000000000000    5.6400000000000000\n\
Na Cl\n\
1 1\n\
Direct\n\
  0.0 0.0 0.0\n\
  0.5 0.5 0.5\n";

    #[test]
    fn parses_basic_vasp5_direct() {
        let doc = parse_poscar(NACL_POSCAR).unwrap();
        assert_eq!(doc.comment, "NaCl");
        assert_eq!(doc.species_order, vec![Element::NA, Element::CL]);
        assert_eq!(doc.structure.site_count(), 2);
        assert_eq!(doc.structure.sites()[0].species[0].element, Element::NA);
        assert_eq!(doc.structure.sites()[1].species[0].element, Element::CL);
        assert_eq!(doc.structure.sites()[0].fractional.0, [0.0, 0.0, 0.0]);
        assert_eq!(doc.structure.sites()[1].fractional.0, [0.5, 0.5, 0.5]);
        assert!(doc.selective_dynamics.is_none());
        assert!(doc.velocities.is_none());
        assert!(doc.predictor_corrector.is_none());
        // No per-atom labels fabricated (see module docs).
        assert!(doc.structure.sites()[0].label.is_none());
    }

    #[test]
    fn empty_comment_line_is_preserved() {
        let text = NACL_POSCAR.replacen("NaCl\n", "\n", 1);
        let doc = parse_poscar(&text).unwrap();
        assert_eq!(doc.comment, "");
        // The rest of the header must still line up correctly.
        assert_eq!(doc.structure.site_count(), 2);
        assert_eq!(doc.structure.lattice().lengths()[0], 5.64);
    }

    #[test]
    fn parses_cartesian_mode() {
        let text = "cart\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C\n\
1\n\
Cartesian\n\
2.0 2.0 2.0\n";
        let doc = parse_poscar(text).unwrap();
        let f = doc.structure.sites()[0].fractional.0;
        assert!((f[0] - 0.5).abs() < 1e-12);
        assert!((f[1] - 0.5).abs() < 1e-12);
        assert!((f[2] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn cartesian_mode_case_insensitive_first_letter_only() {
        for mode in ["c", "C", "k", "K", "Kartesian", "cartesian"] {
            let text = format!(
                "t\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nC\n1\n{mode}\n2.0 2.0 2.0\n"
            );
            let doc = parse_poscar(&text).unwrap();
            let f = doc.structure.sites()[0].fractional.0;
            assert!(
                (f[0] - 0.5).abs() < 1e-9,
                "mode line '{mode}' not treated as Cartesian"
            );
        }
        for mode in ["d", "D", "Direct", "direct", "Fractional"] {
            let text = format!(
                "t\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nC\n1\n{mode}\n0.25 0.25 0.25\n"
            );
            let doc = parse_poscar(&text).unwrap();
            let f = doc.structure.sites()[0].fractional.0;
            assert_eq!(
                f,
                [0.25, 0.25, 0.25],
                "mode line '{mode}' not treated as Direct"
            );
        }
    }

    #[test]
    fn negative_scale_factor_is_target_volume() {
        // Raw cubic cell 2x2x2 = volume 8. Request volume 64 -> linear
        // factor 2 -> cell becomes 4x4x4.
        let text = "vol\n\
-64.0\n\
2.0 0.0 0.0\n\
0.0 2.0 0.0\n\
0.0 0.0 2.0\n\
C\n\
1\n\
Direct\n\
0.0 0.0 0.0\n";
        let doc = parse_poscar(text).unwrap();
        assert!((doc.structure.lattice().volume() - 64.0).abs() < 1e-9);
        for len in doc.structure.lattice().lengths() {
            assert!((len - 4.0).abs() < 1e-9);
        }
    }

    #[test]
    fn negative_scale_factor_with_cartesian_coords_scales_consistently() {
        // Same target-volume cell as above, but the ion position is given
        // in Cartesian mode using the *raw* (unscaled) vector's units --
        // must be scaled by the same derived linear factor as the lattice,
        // not left as a raw-magnitude mirror (a naive `scale * coord` with
        // scale still negative would flip the sign here).
        let text = "vol_cart\n\
-64.0\n\
2.0 0.0 0.0\n\
0.0 2.0 0.0\n\
0.0 0.0 2.0\n\
C\n\
1\n\
Cartesian\n\
1.0 1.0 1.0\n";
        let doc = parse_poscar(text).unwrap();
        let f = doc.structure.sites()[0].fractional.0;
        // Raw Cartesian (1,1,1) with raw cell edge 2 -> fractional (0.5,
        // 0.5, 0.5) regardless of the linear factor applied uniformly to
        // both lattice and coordinate -- and must NOT be negative/mirrored.
        for c in f {
            assert!(
                (c - 0.5).abs() < 1e-9,
                "fractional {f:?} should be (0.5,0.5,0.5)"
            );
        }
    }

    #[test]
    fn three_component_scale_factor_scales_columns_not_rows() {
        // Non-diagonal (triclinic-ish) raw matrix with 3 distinct factors:
        // a column-scale and a row-scale convention diverge here, and
        // produce different angles -- pinning down which one this reader
        // implements (column, per the VASP wiki wording).
        let text = "tri3\n\
2.0 3.0 5.0\n\
1.0 0.0 0.0\n\
0.3 1.0 0.0\n\
0.2 0.1 1.0\n\
Fe\n\
1\n\
Direct\n\
0.0 0.0 0.0\n";
        let doc = parse_poscar(text).unwrap();
        let m = doc.structure.lattice().matrix();
        // Column-scale convention: matrix[row][col] = raw[row][col] * scale[col].
        let expected = [
            [1.0 * 2.0, 0.0 * 3.0, 0.0 * 5.0],
            [0.3 * 2.0, 1.0 * 3.0, 0.0 * 5.0],
            [0.2 * 2.0, 0.1 * 3.0, 1.0 * 5.0],
        ];
        for r in 0..3 {
            for c in 0..3 {
                assert!(
                    (m[r][c] - expected[r][c]).abs() < 1e-9,
                    "m[{r}][{c}] = {}, expected {}",
                    m[r][c],
                    expected[r][c]
                );
            }
        }
        // Cross-check via lengths and angles too. The entry-by-entry matrix
        // assertion above already fully pins down column-vs-row scaling on
        // its own; lengths alone would not (a row-scale convention with
        // these particular factors happens to produce different lengths
        // too, just not by coincidence the *same* ones), so angles are the
        // more targeted secondary check -- row-scaling would change the
        // triangle each pair of (now differently-scaled) row vectors forms
        // to a different degree than column-scaling does.
        let lengths = doc.structure.lattice().lengths();
        let expected_lengths = [
            (expected[0][0].powi(2) + expected[0][1].powi(2) + expected[0][2].powi(2)).sqrt(),
            (expected[1][0].powi(2) + expected[1][1].powi(2) + expected[1][2].powi(2)).sqrt(),
            (expected[2][0].powi(2) + expected[2][1].powi(2) + expected[2][2].powi(2)).sqrt(),
        ];
        for i in 0..3 {
            assert!((lengths[i] - expected_lengths[i]).abs() < 1e-9);
        }
        let angles = doc.structure.lattice().angles_degrees();
        let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        let norm = |u: [f64; 3]| dot(u, u).sqrt();
        let angle_between = |u: [f64; 3], v: [f64; 3]| {
            (dot(u, v) / (norm(u) * norm(v)))
                .clamp(-1.0, 1.0)
                .acos()
                .to_degrees()
        };
        let expected_angles = [
            angle_between(expected[1], expected[2]),
            angle_between(expected[0], expected[2]),
            angle_between(expected[0], expected[1]),
        ];
        for i in 0..3 {
            assert!(
                (angles[i] - expected_angles[i]).abs() < 1e-6,
                "angle[{i}] = {}, expected {}",
                angles[i],
                expected_angles[i]
            );
        }
    }

    #[test]
    fn three_component_scale_rejects_non_positive() {
        let text = "bad3\n\
1.0 -2.0 3.0\n\
1.0 0.0 0.0\n\
0.0 1.0 0.0\n\
0.0 0.0 1.0\n\
C\n\
1\n\
Direct\n\
0.0 0.0 0.0\n";
        let err = parse_poscar(text).unwrap_err();
        assert!(matches!(err, PoscarError::InvalidLine { line: 2, .. }));
    }

    #[test]
    fn selective_dynamics_parsed() {
        let text = "sd\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C O\n\
1 1\n\
Selective dynamics\n\
Direct\n\
0.0 0.0 0.0 T F T\n\
0.5 0.5 0.5 F F F\n";
        let doc = parse_poscar(text).unwrap();
        let sd = doc.selective_dynamics.unwrap();
        assert_eq!(sd, vec![[true, false, true], [false, false, false]]);
    }

    #[test]
    fn velocities_parsed_and_no_mode_line() {
        let text = "v\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C\n\
2\n\
Direct\n\
0.0 0.0 0.0\n\
0.5 0.5 0.5\n\
\n\
0.1 0.2 0.3\n\
-0.1 0.0 0.05\n";
        let doc = parse_poscar(text).unwrap();
        let v = doc.velocities.unwrap();
        assert_eq!(v, vec![[0.1, 0.2, 0.3], [-0.1, 0.0, 0.05]]);
    }

    #[test]
    fn velocities_bad_line_fails_closed() {
        let text = "v\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C\n\
1\n\
Direct\n\
0.0 0.0 0.0\n\
\n\
not a velocity\n";
        let err = parse_poscar(text).unwrap_err();
        assert!(matches!(err, PoscarError::InvalidLine { .. }));
    }

    #[test]
    fn predictor_corrector_round_trips_verbatim() {
        let text = "pc\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C\n\
1\n\
Direct\n\
0.0 0.0 0.0\n\
\n\
0.1 0.2 0.3\n\
\n\
1\n\
0.5\n\
1.0 2.0\n\
0.01 0.02 0.03\n\
0.04 0.05 0.06\n\
0.07 0.08 0.09\n";
        let doc = parse_poscar(text).unwrap();
        let pc = doc.predictor_corrector.clone().unwrap();
        assert_eq!(
            pc.preamble,
            ["1".to_string(), "0.5".to_string(), "1.0 2.0".to_string()]
        );
        assert_eq!(
            pc.lines,
            vec![
                "0.01 0.02 0.03".to_string(),
                "0.04 0.05 0.06".to_string(),
                "0.07 0.08 0.09".to_string(),
            ]
        );

        let written = write_poscar(&doc).unwrap();
        let reparsed = parse_poscar(&written).unwrap();
        assert_eq!(reparsed.predictor_corrector, doc.predictor_corrector);
        assert_eq!(reparsed.velocities, doc.velocities);
    }

    #[test]
    fn lattice_velocities_block_rejected_as_unsupported() {
        let text = "lv\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C\n\
1\n\
Direct\n\
0.0 0.0 0.0\n\
Lattice velocities and vectors\n\
1\n\
0.0 0.0 0.0\n\
0.0 0.0 0.0\n\
0.0 0.0 0.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n";
        let err = parse_poscar(text).unwrap_err();
        assert!(matches!(err, PoscarError::UnsupportedSection { .. }));
    }

    #[test]
    fn vasp4_style_rejected_not_silently_misparsed() {
        // No species-name line: line 6 goes straight to atom counts.
        let text = "v4\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
1 1\n\
Direct\n\
0.0 0.0 0.0\n\
0.5 0.5 0.5\n";
        let err = parse_poscar(text).unwrap_err();
        assert_eq!(err, PoscarError::Vasp4NotSupported);
    }

    #[test]
    fn unknown_element_rejected() {
        let text = "bad\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nXx\n1\nDirect\n0.0 0.0 0.0\n";
        let err = parse_poscar(text).unwrap_err();
        assert!(matches!(err, PoscarError::UnknownElement { .. }));
    }

    #[test]
    fn atom_count_mismatch_detected() {
        let text = "short\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nC\n2\nDirect\n0.0 0.0 0.0\n";
        let err = parse_poscar(text).unwrap_err();
        assert_eq!(
            err,
            PoscarError::AtomCountMismatch {
                declared: 2,
                found: 1
            }
        );
    }

    #[test]
    fn atom_count_mismatch_detected_when_trailing_section_follows() {
        // Coordinate block is short by one, but (unlike the fixture above,
        // which simply ends at EOF) a properly blank-line-separated
        // trailing section follows -- must still report a clear
        // AtomCountMismatch, not consume the trailing section's line as
        // bogus coordinate data or surface a confusing "coordinate line
        // needs 3 fields, got ''" error from the blank separator itself.
        let text = "short-trailer\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nC\n2\nDirect\n0.0 0.0 0.0\n\n0.1 0.2 0.3\n";
        let err = parse_poscar(text).unwrap_err();
        assert_eq!(
            err,
            PoscarError::AtomCountMismatch {
                declared: 2,
                found: 1
            }
        );
    }

    #[test]
    fn contcar_alias_matches_poscar_parse() {
        assert_eq!(
            parse_contcar(NACL_POSCAR).unwrap(),
            parse_poscar(NACL_POSCAR).unwrap()
        );
    }

    // -- round-trip -----------------------------------------------------

    #[test]
    fn roundtrip_triclinic_fixture_text() {
        // A genuine triclinic POSCAR (not built via Lattice::from_parameters
        // -- this is literal file text, per the task's round-trip
        // requirement): a, b, c = 5, 6, 7 with distinct, non-90-degree
        // angles, two species, non-alphabetical group order (O before C)
        // to also exercise "preserve file order, don't alphabetize".
        let text = "triclinic fixture\n\
1.0\n\
5.0000000000000000 0.0000000000000000 0.0000000000000000\n\
1.5529142706151024 5.7996130342078438 0.0000000000000000\n\
-1.4790559332561876 1.4998727130474476 6.5850758980592793\n\
O C\n\
2 1\n\
Direct\n\
0.1 0.2 0.3\n\
0.6 0.1 0.9\n\
0.4 0.4 0.4\n";
        let doc = parse_poscar(text).unwrap();
        assert_eq!(doc.species_order, vec![Element::O, Element::C]);

        let written = write_poscar(&doc).unwrap();
        let reparsed = parse_poscar(&written).unwrap();

        assert_eq!(doc.structure, reparsed.structure);
        assert_eq!(doc.comment, reparsed.comment);
        assert_eq!(doc.selective_dynamics, reparsed.selective_dynamics);
        assert_eq!(doc.velocities, reparsed.velocities);
        assert_eq!(doc.predictor_corrector, reparsed.predictor_corrector);
        // Species header order is re-derived from site order (first
        // appearance), which for this fixture matches the original file's
        // own O-then-C group order -- but note write_poscar does not
        // consult doc.species_order to produce it (see docs).
        assert_eq!(reparsed.species_order, vec![Element::O, Element::C]);

        // Lattice vectors match to float precision (write emits exact
        // Display-formatted f64s, which round-trip losslessly).
        let (m0, m1) = (
            doc.structure.lattice().matrix(),
            reparsed.structure.lattice().matrix(),
        );
        for r in 0..3 {
            for c in 0..3 {
                assert!((m0[r][c] - m1[r][c]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn roundtrip_triclinic_cartesian_mode_fixture_text() {
        // Same triclinic lattice as `roundtrip_triclinic_fixture_text`, but
        // ion positions given in Cartesian mode -- this is the path
        // write_poscar's "always emit Direct" decision actually transforms
        // (cart_to_frac on read, never converted back), so it needs its own
        // coverage rather than relying on the Direct-mode fixture, and a
        // tolerance comparison rather than `==` since cart_to_frac's matrix
        // solve introduces ordinary floating-point rounding.
        let text = "triclinic cartesian fixture\n\
1.0\n\
5.0000000000000000 0.0000000000000000 0.0000000000000000\n\
1.5529142706151024 5.7996130342078438 0.0000000000000000\n\
-1.4790559332561876 1.4998727130474476 6.5850758980592793\n\
O C\n\
2 1\n\
Cartesian\n\
1.1 0.9 1.3\n\
2.4 3.1 4.2\n\
0.5 0.5 0.5\n";
        let doc = parse_poscar(text).unwrap();
        assert_eq!(doc.structure.site_count(), 3);

        let written = write_poscar(&doc).unwrap();
        let reparsed = parse_poscar(&written).unwrap();

        assert_eq!(doc.structure.site_count(), reparsed.structure.site_count());
        let (f0, f1) = (
            doc.structure.fractional_positions(),
            reparsed.structure.fractional_positions(),
        );
        for (a, b) in f0.iter().zip(f1.iter()) {
            for k in 0..3 {
                assert!(
                    (a.0[k] - b.0[k]).abs() < 1e-9,
                    "fractional coord mismatch: {a:?} vs {b:?}"
                );
            }
        }
        let (c0, c1) = (
            doc.structure.cartesian_positions(),
            reparsed.structure.cartesian_positions(),
        );
        for (a, b) in c0.iter().zip(c1.iter()) {
            for k in 0..3 {
                assert!(
                    (a.0[k] - b.0[k]).abs() < 1e-9,
                    "cartesian coord mismatch: {a:?} vs {b:?}"
                );
            }
        }
        let (m0, m1) = (
            doc.structure.lattice().matrix(),
            reparsed.structure.lattice().matrix(),
        );
        for r in 0..3 {
            for c in 0..3 {
                assert!((m0[r][c] - m1[r][c]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn roundtrip_with_selective_dynamics_and_velocities() {
        let text = "sdv\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
Fe Ni\n\
1 1\n\
Selective dynamics\n\
Direct\n\
0.0 0.0 0.0 T T F\n\
0.5 0.5 0.5 F F F\n\
\n\
0.01 0.0 0.0\n\
0.0 0.0 0.0\n";
        let doc = parse_poscar(text).unwrap();
        let written = write_poscar(&doc).unwrap();
        let reparsed = parse_poscar(&written).unwrap();
        assert_eq!(doc, reparsed);
    }

    #[test]
    fn write_species_grouped_by_first_appearance_not_alphabetical() {
        let lattice = Lattice::cubic(4.0).unwrap();
        let sites = vec![
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::ZN)],
                FractionalCoord::new([0.0, 0.0, 0.0]),
                None,
            )
            .unwrap(),
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::AL)],
                FractionalCoord::new([0.25, 0.25, 0.25]),
                None,
            )
            .unwrap(),
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::ZN)],
                FractionalCoord::new([0.5, 0.5, 0.5]),
                None,
            )
            .unwrap(),
        ];
        let structure = PeriodicStructure::new(lattice, sites).unwrap();
        let doc = PoscarDocument {
            structure,
            comment: "zn-al-zn".to_string(),
            species_order: vec![],
            selective_dynamics: None,
            velocities: None,
            predictor_corrector: None,
        };
        let written = write_poscar(&doc).unwrap();
        let lines: Vec<&str> = written.lines().collect();
        // Zn appears before Al in site order (first appearance), and
        // alphabetically Al < Zn -- confirms grouping is by first
        // appearance, not alphabetical.
        assert_eq!(lines[5], "Zn Al");
        assert_eq!(lines[6], "2 1");

        let reparsed = parse_poscar(&written).unwrap();
        assert_eq!(reparsed.structure.site_count(), 3);
        // Original site order was Zn, Al, Zn -- after grouping, Zn's two
        // sites (original fractional (0,0,0) and (0.5,0.5,0.5)) must stay
        // in that relative order, now contiguous, followed by Al.
        let frac: Vec<[f64; 3]> = reparsed
            .structure
            .sites()
            .iter()
            .map(|s| s.fractional.0)
            .collect();
        assert_eq!(frac[0], [0.0, 0.0, 0.0]);
        assert_eq!(frac[1], [0.5, 0.5, 0.5]);
        assert_eq!(frac[2], [0.25, 0.25, 0.25]);
    }

    #[test]
    fn write_rejects_disordered_site() {
        let lattice = Lattice::cubic(4.0).unwrap();
        let sites = vec![
            PeriodicSite::new(
                vec![
                    SiteSpecies {
                        element: Element::FE,
                        occupancy: Occupancy::new(0.6).unwrap(),
                    },
                    SiteSpecies {
                        element: Element::NI,
                        occupancy: Occupancy::new(0.4).unwrap(),
                    },
                ],
                FractionalCoord::new([0.0, 0.0, 0.0]),
                None,
            )
            .unwrap(),
        ];
        let structure = PeriodicStructure::new(lattice, sites).unwrap();
        let doc = PoscarDocument {
            structure,
            comment: "disordered".to_string(),
            species_order: vec![],
            selective_dynamics: None,
            velocities: None,
            predictor_corrector: None,
        };
        let err = write_poscar(&doc).unwrap_err();
        assert!(matches!(err, PoscarError::Unwritable { .. }));
    }

    #[test]
    fn write_rejects_partial_occupancy() {
        let lattice = Lattice::cubic(4.0).unwrap();
        let sites = vec![
            PeriodicSite::new(
                vec![SiteSpecies {
                    element: Element::FE,
                    occupancy: Occupancy::new(0.5).unwrap(),
                }],
                FractionalCoord::new([0.0, 0.0, 0.0]),
                None,
            )
            .unwrap(),
        ];
        let structure = PeriodicStructure::new(lattice, sites).unwrap();
        let doc = PoscarDocument {
            structure,
            comment: "vacancy".to_string(),
            species_order: vec![],
            selective_dynamics: None,
            velocities: None,
            predictor_corrector: None,
        };
        let err = write_poscar(&doc).unwrap_err();
        assert!(matches!(err, PoscarError::Unwritable { .. }));
    }

    #[test]
    fn write_rejects_predictor_corrector_when_reorder_needed() {
        // Site order Al, Zn, Al is NOT contiguous per species -- writing
        // requires reordering, which is unsafe with opaque
        // predictor_corrector data present.
        let lattice = Lattice::cubic(4.0).unwrap();
        let sites = vec![
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::AL)],
                FractionalCoord::new([0.0, 0.0, 0.0]),
                None,
            )
            .unwrap(),
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::ZN)],
                FractionalCoord::new([0.25, 0.25, 0.25]),
                None,
            )
            .unwrap(),
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::AL)],
                FractionalCoord::new([0.5, 0.5, 0.5]),
                None,
            )
            .unwrap(),
        ];
        let structure = PeriodicStructure::new(lattice, sites).unwrap();
        let doc = PoscarDocument {
            structure,
            comment: "needs-reorder".to_string(),
            species_order: vec![],
            selective_dynamics: None,
            velocities: Some(vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            predictor_corrector: Some(PredictorCorrector {
                preamble: ["1".to_string(), "0.5".to_string(), "1.0 2.0".to_string()],
                lines: vec!["0.0 0.0 0.0".to_string()],
            }),
        };
        let err = write_poscar(&doc).unwrap_err();
        assert!(matches!(err, PoscarError::Unwritable { .. }));
    }

    #[test]
    fn write_rejects_predictor_corrector_without_velocities() {
        let lattice = Lattice::cubic(4.0).unwrap();
        let sites = vec![
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::AL)],
                FractionalCoord::new([0.0, 0.0, 0.0]),
                None,
            )
            .unwrap(),
        ];
        let structure = PeriodicStructure::new(lattice, sites).unwrap();
        let doc = PoscarDocument {
            structure,
            comment: "pc-no-vel".to_string(),
            species_order: vec![],
            selective_dynamics: None,
            velocities: None,
            predictor_corrector: Some(PredictorCorrector {
                preamble: ["1".to_string(), "0.5".to_string(), "1.0 2.0".to_string()],
                lines: vec!["0.0 0.0 0.0".to_string()],
            }),
        };
        let err = write_poscar(&doc).unwrap_err();
        assert!(matches!(err, PoscarError::Unwritable { .. }));
    }

    #[test]
    fn empty_input_rejected() {
        assert_eq!(parse_poscar("").unwrap_err(), PoscarError::Empty);
    }
}
