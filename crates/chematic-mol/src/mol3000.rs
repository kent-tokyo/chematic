//! MOL V3000 (Extended Ctab) parser.
//!
//! Reference: MDL/Dassault Systèmes CTfile Formats specification, V3000 section.
//!
//! MOL V3000 (Extended Ctab) parser and writer.
//!
//! Reference: MDL/Dassault Systèmes CTfile Formats specification, V3000 section.

use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Coords3D, Element, Molecule, MoleculeBuilder, Point3,
    StereoGroup, StereoGroupKind,
};
use chematic_perception::{
    apply_ez_directions_from_2d_ex, apply_local_parity_from_wedges_with_diagnostics,
};

use crate::error::MolParseError;
use crate::mol2000::{
    CoordinateDimension, GeometryRank, MolMetadata, MolReadReport, Stereo3DDiagnostic,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Prefix that marks every V3000 continuation line.
const V30_PREFIX: &str = "M  V30 ";

/// Return a `V3000ParseError` for the given 1-based line number and message.
#[inline]
fn v3k_err(line: usize, msg: impl Into<String>) -> MolParseError {
    MolParseError::V3000ParseError {
        line,
        msg: msg.into(),
    }
}

/// True when the first two tokens match the given block marker (e.g. `END CTAB`).
fn is_marker(tokens: &[&str], kw1: &str, kw2: &str) -> bool {
    tokens.len() >= 2 && tokens[0] == kw1 && tokens[1] == kw2
}

// ---------------------------------------------------------------------------
// Line-continuation pre-pass
// ---------------------------------------------------------------------------

/// A logical V30 line with its 1-based source line number (the first physical
/// line that makes up this logical line).
struct LogicalLine {
    /// 1-based line number of the first physical line.
    line_num: usize,
    /// Payload after stripping the `M  V30 ` prefix and joining continuations.
    payload: String,
}

/// Pre-process the raw input lines into logical V30 lines.
///
/// Lines that do **not** start with `M  V30 ` are silently skipped; we only
/// care about V30 content lines here (the header and counts line are read
/// separately before calling this function).
///
/// When a V30 payload ends with a hyphen `-` the hyphen is stripped and the
/// payload of the next `M  V30 ` line is appended (with a single space
/// separator so that tokens do not merge).
fn collect_v30_lines(lines: &[(usize, &str)]) -> Vec<LogicalLine> {
    let mut result: Vec<LogicalLine> = Vec::new();

    let mut iter = lines.iter().peekable();
    while let Some(&(lineno, raw)) = iter.next() {
        if let Some(payload) = raw.strip_prefix(V30_PREFIX) {
            let mut text = payload.to_string();
            let first_line = lineno;

            // Handle line continuations: a payload ending with `-` means the
            // next `M  V30 ` line is a continuation.
            while text.ends_with('-') {
                // Drop the trailing hyphen.
                text.pop();
                // Consume the next V30 line if available.
                match iter.next() {
                    Some(&(_, cont_raw)) => {
                        if let Some(cont_payload) = cont_raw.strip_prefix(V30_PREFIX) {
                            // Append with a separating space so that tokens
                            // that were split across lines re-join cleanly.
                            text.push(' ');
                            text.push_str(cont_payload);
                        }
                        // If the continuation line is not a V30 line, stop.
                    }
                    None => break,
                }
            }

            result.push(LogicalLine {
                line_num: first_line,
                payload: text,
            });
        }
        // Non-V30 lines (header lines, blank lines, M  END) are ignored here.
    }

    result
}

// ---------------------------------------------------------------------------
// Key=value parsing helper
// ---------------------------------------------------------------------------

/// Parse optional `KEY=VALUE` pairs from the tail of a V30 atom or bond line.
///
/// Tokens are space-separated. Any token containing `=` is treated as a
/// key-value pair; tokens without `=` that appear after the required positional
/// fields are silently ignored.
fn parse_kv(tokens: &[&str], key: &str) -> Option<String> {
    for tok in tokens {
        if let Some(rest) = tok.strip_prefix(key)
            && let Some(val) = rest.strip_prefix('=')
        {
            return Some(val.to_string());
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse a MOL V3000 (Extended Ctab) block.
///
/// `input` should be the full text of a V3000 MOL block, starting from the
/// molecule name (header line 1) through `M  END`.
///
/// Returns `(Molecule, MolMetadata)` on success.
/// Parse a MOL V3000 string, running stereo perception and returning every
/// rejected wedge/hash center as a structured `StereoDiagnostic`.
///
/// `coords[i]` is the `(x, y)` position for atom `i` extracted from the V30 atom block.
/// Z-coordinates are not captured (V3000 stores 3D; we retain only the 2D projection).
/// This is the one parsing core for V3000 MOL text --
/// [`parse_mol_v3000_with_coords`]/[`parse_mol_v3000`] are thin wrappers that
/// discard `stereo_diagnostics`. Only bond-line `CFG` (wedge direction) is
/// decoded; atom-line `CFG` (parity) is out of scope, matching RDKit's own
/// primary `MolFromMolBlock` path (bond-CFG/wedge-based).
pub fn read_mol_v3000_with_diagnostics(input: &str) -> Result<MolReadReport, MolParseError> {
    // Collect all physical lines with 1-based numbering.
    let all_lines: Vec<(usize, &str)> =
        input.lines().enumerate().map(|(i, l)| (i + 1, l)).collect();

    if all_lines.len() < 4 {
        return Err(MolParseError::UnexpectedEnd);
    }

    // -- Header lines 1–3 ----------------------------------------------------

    let name = all_lines[0].1.to_string();
    // Line 2 (program/date) is mostly discarded, except for the
    // dimensional-code field (columns 20..22, "2D"/"3D") -- same column as
    // V2000's header (see `crate::mol2000::parse_dimension_code`); V3000
    // keeps the identical 3-line header before its own counts line.
    let line2_raw = all_lines[1].1;
    let comment = all_lines[2].1.to_string();

    let metadata = MolMetadata { name, comment };

    // -- Counts line (line 4) — verify V3000 tag ----------------------------

    let (counts_lineno, counts_line) = all_lines[3];
    if !counts_line.contains("V3000") {
        return Err(MolParseError::InvalidCountLine {
            line: counts_lineno,
            detail: "missing V3000 version tag".to_string(),
        });
    }

    // -- Pre-process V30 logical lines --------------------------------------

    let v30_lines = collect_v30_lines(&all_lines);

    // -- State machine over logical V30 lines --------------------------------

    // Expected sequence:
    //   BEGIN CTAB
    //   COUNTS <na> <nb> ...
    //   BEGIN ATOM
    //     <atom lines>
    //   END ATOM
    //   BEGIN BOND          (optional — may be absent if nb == 0)
    //     <bond lines>
    //   END BOND
    //   END CTAB

    let mut builder = MoleculeBuilder::new();

    // Map from 1-based V3000 atom index to the 0-based builder index.
    // V3000 atom indices are not required to be contiguous (though they
    // usually are), so we track them explicitly.
    let mut atom_idx_map: Vec<(u32, AtomIdx)> = Vec::new(); // (v3k_idx, builder_idx)

    // Collect (x, y) coordinates in the order atoms are added to builder.
    let mut coords: Vec<(f64, f64)> = Vec::new();
    // Real z coordinates, parallel to `coords` -- previously discarded
    // entirely (root cause of the 3D-coordinate-loss bug this PR fixes).
    let mut raw_z: Vec<f64> = Vec::new();

    enum State {
        BeforeCtab,
        InCtab,
        InAtomBlock,
        AfterAtomBlock,
        InBondBlock,
        AfterBondBlock,
        InCollection,
        Done,
    }

    let mut state = State::BeforeCtab;
    let mut expected_atoms: usize = 0;
    let mut stereo_groups: Vec<StereoGroup> = Vec::new();
    // Double bonds whose `CFG=2` marks explicitly unspecified E/Z (the same
    // token V3000 also uses for a single bond's "either" wedge) -- confirmed
    // against a live RDKit 2026.03.3 oracle (B0 diagnosis): RDKit's own
    // V3000 writer emits `CFG=2` on a double-bond line for `STEREOANY`.
    // Threaded into `apply_ez_directions_from_2d_ex` below, same as V2000's
    // stereo-code-3 handling in `mol2000.rs`.
    let mut explicitly_unspecified_ez: std::collections::HashSet<BondIdx> =
        std::collections::HashSet::new();

    for LogicalLine { line_num, payload } in &v30_lines {
        let lnum = *line_num;
        let tokens: Vec<&str> = payload.split_whitespace().collect();
        if tokens.is_empty() {
            continue;
        }

        match state {
            State::BeforeCtab => {
                if is_marker(&tokens, "BEGIN", "CTAB") {
                    state = State::InCtab;
                }
            }

            State::InCtab => {
                if tokens[0] == "COUNTS" {
                    if tokens.len() < 3 {
                        return Err(v3k_err(lnum, "COUNTS line has fewer than 2 values"));
                    }
                    expected_atoms = tokens[1].parse::<usize>().map_err(|_| {
                        v3k_err(
                            lnum,
                            format!("cannot parse atom count from '{}'", tokens[1]),
                        )
                    })?;
                    // Parse bond count to surface malformed COUNTS lines, even
                    // though the value is not used downstream.
                    tokens[2].parse::<usize>().map_err(|_| {
                        v3k_err(
                            lnum,
                            format!("cannot parse bond count from '{}'", tokens[2]),
                        )
                    })?;
                } else if is_marker(&tokens, "BEGIN", "ATOM") {
                    state = State::InAtomBlock;
                } else if is_marker(&tokens, "END", "CTAB") {
                    state = State::Done;
                }
            }

            State::InAtomBlock => {
                if is_marker(&tokens, "END", "ATOM") {
                    state = State::AfterAtomBlock;
                    if builder.atom_count() != expected_atoms {
                        return Err(v3k_err(
                            lnum,
                            format!(
                                "expected {} atoms, found {}",
                                expected_atoms,
                                builder.atom_count()
                            ),
                        ));
                    }
                    continue;
                }

                // Atom line: <idx> <symbol> <x> <y> <z> <aamap> [KEY=VAL ...]
                if tokens.len() < 6 {
                    return Err(MolParseError::InvalidAtomLine {
                        line: lnum,
                        detail: format!(
                            "V3000 atom line needs at least 6 fields, got {}",
                            tokens.len()
                        ),
                    });
                }

                let v3k_idx =
                    tokens[0]
                        .parse::<u32>()
                        .map_err(|_| MolParseError::InvalidAtomLine {
                            line: lnum,
                            detail: format!("cannot parse atom index from '{}'", tokens[0]),
                        })?;

                // Strip bracket notation e.g. "[OH]" → "OH", "[R]" → "R".
                let sym = tokens[1].trim_start_matches('[').trim_end_matches(']');

                let element =
                    Element::from_symbol(sym).ok_or_else(|| MolParseError::UnknownElement {
                        symbol: sym.to_string(),
                        line: lnum,
                    })?;

                // Parse x, y coordinates (tokens[2] and tokens[3]) --
                // unchanged, still lenient (a malformed x/y silently
                // defaults to 0.0, matching pre-existing behavior).
                let x: f64 = tokens[2].parse().unwrap_or(0.0);
                let y: f64 = tokens[3].parse().unwrap_or(0.0);

                // z coordinate (tokens[4]): previously silently discarded
                // entirely -- root cause of the 3D-coordinate-loss bug this
                // PR fixes. Unlike x/y above, a garbled or non-finite z is a
                // typed error rather than a silent 0.0 default: nothing
                // downstream reads z today, so a malformed value can only be
                // file corruption. The token always exists syntactically
                // (the `tokens.len() < 6` check above already guarantees
                // it), so unlike V2000's "line too short" leniency, there is
                // no "missing field" case here to default instead of error.
                let z: f64 = tokens[4]
                    .parse()
                    .map_err(|_| MolParseError::InvalidAtomLine {
                        line: lnum,
                        detail: format!("cannot parse z coordinate from '{}'", tokens[4]),
                    })?;
                if !z.is_finite() {
                    return Err(MolParseError::InvalidAtomLine {
                        line: lnum,
                        detail: format!(
                            "z coordinate is not finite (NaN/Infinite): '{}'",
                            tokens[4]
                        ),
                    });
                }

                // Atom-map number (positional field 6, 0 = no mapping).
                let aamap_raw = tokens[5].parse::<u16>().unwrap_or(0);
                let atom_map = if aamap_raw == 0 {
                    None
                } else {
                    Some(aamap_raw)
                };

                let kv_tokens = tokens.get(6..).unwrap_or(&[]);

                let charge: i8 = parse_kv(kv_tokens, "CHG")
                    .and_then(|v| v.parse::<i8>().ok())
                    .unwrap_or(0);

                let isotope: Option<u16> =
                    parse_kv(kv_tokens, "MASS").and_then(|v| v.parse::<u16>().ok());

                // HCOUNT: -1 means unspecified; treat as None.
                let hydrogen_count: Option<u8> = parse_kv(kv_tokens, "HCOUNT").and_then(|v| {
                    let n: i32 = v.parse().ok()?;
                    if n < 0 { None } else { Some(n as u8) }
                });

                let mut atom = Atom::new(element);
                atom.charge = charge;
                atom.isotope = isotope;
                atom.hydrogen_count = hydrogen_count;
                atom.atom_map = atom_map;

                let builder_idx = builder.add_atom(atom);
                atom_idx_map.push((v3k_idx, builder_idx));
                coords.push((x, y));
                raw_z.push(z);
            }

            State::AfterAtomBlock => {
                if is_marker(&tokens, "BEGIN", "BOND") {
                    state = State::InBondBlock;
                } else if is_marker(&tokens, "END", "CTAB") {
                    state = State::Done;
                }
            }

            State::InBondBlock => {
                if is_marker(&tokens, "END", "BOND") {
                    state = State::AfterBondBlock;
                    continue;
                }

                // Bond line: <idx> <type> <atom1> <atom2> [KEY=VAL ...]
                if tokens.len() < 4 {
                    return Err(MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: format!(
                            "V3000 bond line needs at least 4 fields, got {}",
                            tokens.len()
                        ),
                    });
                }

                let btype_raw =
                    tokens[1]
                        .parse::<u8>()
                        .map_err(|_| MolParseError::InvalidBondLine {
                            line: lnum,
                            detail: format!("cannot parse bond type from '{}'", tokens[1]),
                        })?;

                let a1_v3k =
                    tokens[2]
                        .parse::<u32>()
                        .map_err(|_| MolParseError::InvalidBondLine {
                            line: lnum,
                            detail: format!("cannot parse atom1 index from '{}'", tokens[2]),
                        })?;

                let a2_v3k =
                    tokens[3]
                        .parse::<u32>()
                        .map_err(|_| MolParseError::InvalidBondLine {
                            line: lnum,
                            detail: format!("cannot parse atom2 index from '{}'", tokens[3]),
                        })?;

                let a1 = resolve_atom_idx(a1_v3k, &atom_idx_map).ok_or_else(|| {
                    MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: format!("atom index {a1_v3k} not found in atom block"),
                    }
                })?;

                let a2 = resolve_atom_idx(a2_v3k, &atom_idx_map).ok_or_else(|| {
                    MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: format!("atom index {a2_v3k} not found in atom block"),
                    }
                })?;

                let kv_tokens = tokens.get(4..).unwrap_or(&[]);

                let order = match btype_raw {
                    0 => BondOrder::Zero,
                    // Bond CFG (wedge direction): 1=Up, 3=Down, 2=Either --
                    // only meaningful for a single bond, mirroring V2000's
                    // own gating. CFG=2 ("either"/unspecified) is left as
                    // Single -- a defined wedge/hash needs a definite
                    // direction, not "unknown" (same policy as V2000's own
                    // code-4 handling).
                    1 => match parse_kv(kv_tokens, "CFG").as_deref() {
                        Some("1") => BondOrder::Up,
                        Some("3") => BondOrder::Down,
                        _ => BondOrder::Single,
                    },
                    2 => BondOrder::Double,
                    3 => BondOrder::Triple,
                    4 => BondOrder::Aromatic,
                    5 => BondOrder::QuerySingleOrDouble,
                    6 => BondOrder::QuerySingleOrAromatic,
                    7 => BondOrder::QueryDoubleOrAromatic,
                    8 => BondOrder::QueryAny,
                    // 9 = dative/coordinate bond. This is how RDKit itself
                    // writes `Bond::BondType::DATIVE` in V3000 (it cannot
                    // represent a dative bond in V2000 at all, so it upgrades
                    // automatically) -- `a1`/`a2` are already in the same
                    // donor/acceptor order the file encodes, matching
                    // `BondOrder::Dative`'s documented atom1(donor) ->
                    // atom2(acceptor) convention. Previously fell through to
                    // `Single`, silently discarding coordination bonds from
                    // any RDKit- (or other-tool-) generated V3000 molfile.
                    9 => BondOrder::Dative,
                    _ => BondOrder::Single,
                };

                let bidx = builder.add_bond(a1, a2, order).map_err(|e| {
                    MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: e.to_string(),
                    }
                })?;

                if btype_raw == 2 && parse_kv(kv_tokens, "CFG").as_deref() == Some("2") {
                    explicitly_unspecified_ez.insert(bidx);
                }
            }

            State::AfterBondBlock => {
                if is_marker(&tokens, "BEGIN", "COLLECTION") {
                    state = State::InCollection;
                } else if is_marker(&tokens, "END", "CTAB") {
                    state = State::Done;
                }
            }

            State::InCollection => {
                if is_marker(&tokens, "END", "COLLECTION") {
                    state = State::AfterBondBlock;
                } else if let Some(group) = parse_stereo_group_line(payload, &atom_idx_map) {
                    stereo_groups.push(group);
                }
            }

            State::Done => {}
        }
    }

    // Validate that we reached a terminal state with the expected atom block.
    match state {
        State::Done | State::AfterBondBlock => {}
        State::InAtomBlock => {
            return Err(MolParseError::V3000ParseError {
                line: 0,
                msg: "missing M  V30 END ATOM".to_string(),
            });
        }
        State::InBondBlock => {
            return Err(MolParseError::V3000ParseError {
                line: 0,
                msg: "missing M  V30 END BOND".to_string(),
            });
        }
        _ => {
            return Err(MolParseError::UnexpectedEnd);
        }
    }

    let mut mol = builder.build();
    if !stereo_groups.is_empty() {
        mol.set_stereo_groups(stereo_groups);
    }
    // Tetrahedral parity first, then E/Z direction (side channel only) --
    // same ordering rationale as `mol2000.rs`.
    let stereo_diagnostics = apply_local_parity_from_wedges_with_diagnostics(&mut mol, &coords);
    let ez_diagnostics =
        apply_ez_directions_from_2d_ex(&mut mol, &coords, &explicitly_unspecified_ez);

    // Same 3D bookkeeping as `mol2000.rs::read_mol_with_diagnostics` -- see
    // that function's comments for the full rationale (kept in one place,
    // reused via `pub(crate)` helpers, not re-derived here).
    let coordinate_dimension = crate::mol2000::parse_dimension_code(line2_raw);
    let points: Vec<Point3> = coords
        .iter()
        .zip(raw_z.iter())
        .map(|(&(x, y), &z)| Point3::new(x, y, z))
        .collect();
    let geometry_rank = crate::mol2000::classify_geometry_rank(&points);
    let conformer = match geometry_rank {
        GeometryRank::Coplanar | GeometryRank::ThreeD => Some(Coords3D { points }),
        GeometryRank::FlatZero | GeometryRank::Indeterminate => None,
    };

    let mut stereo3d_diagnostics = Vec::new();
    match (coordinate_dimension, geometry_rank) {
        (CoordinateDimension::TwoD, GeometryRank::Coplanar | GeometryRank::ThreeD) => {
            stereo3d_diagnostics.push(Stereo3DDiagnostic::DeclaredTwoDButNonzeroZ {
                observed: geometry_rank,
            });
        }
        (CoordinateDimension::ThreeD, GeometryRank::FlatZero | GeometryRank::Coplanar) => {
            stereo3d_diagnostics.push(Stereo3DDiagnostic::DeclaredThreeDButFlat {
                observed: geometry_rank,
            });
        }
        _ => {}
    }
    // Same ordering rationale as `mol2000.rs::read_mol_with_diagnostics`:
    // square-planar reperception first (may set `Chirality::SquarePlanar`),
    // then `wedge_vs_3d_conflicts` (whose `is_tetrahedral()` gate must see
    // the final chirality).
    let mut square_planar_diagnostics = Vec::new();
    if let Some(ref conf) = conformer {
        square_planar_diagnostics = crate::mol2000::perceive_square_planar_from_3d(&mut mol, conf);
        stereo3d_diagnostics.extend(crate::mol2000::wedge_vs_3d_conflicts(&mol, conf));
    }

    Ok(MolReadReport {
        mol,
        metadata,
        coords,
        stereo_diagnostics,
        ez_diagnostics,
        conformer,
        coordinate_dimension,
        geometry_rank,
        stereo3d_diagnostics,
        square_planar_diagnostics,
    })
}

/// Parse a MOL V3000 string into a `(Molecule, MolMetadata, Vec<(f64, f64)>)` triple.
///
/// Thin wrapper around [`read_mol_v3000_with_diagnostics`] that discards
/// `stereo_diagnostics` -- signature and behavior unchanged from before
/// stereo perception was wired in.
#[allow(clippy::type_complexity)]
pub fn parse_mol_v3000_with_coords(
    input: &str,
) -> Result<(Molecule, MolMetadata, Vec<(f64, f64)>), MolParseError> {
    read_mol_v3000_with_diagnostics(input).map(|r| (r.mol, r.metadata, r.coords))
}

/// Parse a MOL V3000 string into a `(Molecule, MolMetadata)` pair.
///
/// Coordinates from the atom block are discarded. Use `parse_mol_v3000_with_coords`
/// to retain 2D coordinates.
pub fn parse_mol_v3000(input: &str) -> Result<(Molecule, MolMetadata), MolParseError> {
    parse_mol_v3000_with_coords(input).map(|(mol, meta, _coords)| (mol, meta))
}

/// Look up the builder `AtomIdx` for a given V3000 1-based index.
fn resolve_atom_idx(v3k_idx: u32, map: &[(u32, AtomIdx)]) -> Option<AtomIdx> {
    map.iter().find(|&&(k, _)| k == v3k_idx).map(|&(_, v)| v)
}

/// Parse a V3000 COLLECTION line into a [`StereoGroup`].
///
/// Expected formats:
/// - `MDLV30/STEABS ATOMS=(3 1 2 3)` → `StereoGroupKind::Absolute`
/// - `MDLV30/STEOR1 ATOMS=(1 4)`     → `StereoGroupKind::Or(1)`
/// - `MDLV30/STEAND2 ATOMS=(1 6)`    → `StereoGroupKind::And(2)`
fn parse_stereo_group_line(payload: &str, atom_idx_map: &[(u32, AtomIdx)]) -> Option<StereoGroup> {
    // First token is the group kind key.
    let first_tok = payload.split_whitespace().next()?;

    let kind = if first_tok == "MDLV30/STEABS" {
        StereoGroupKind::Absolute
    } else if let Some(n_str) = first_tok.strip_prefix("MDLV30/STEOR") {
        let n: u32 = n_str.parse().ok()?;
        StereoGroupKind::Or(n)
    } else {
        let n_str = first_tok.strip_prefix("MDLV30/STEAND")?;
        let n: u32 = n_str.parse().ok()?;
        StereoGroupKind::And(n)
    };

    // Extract the ATOMS=(...) value from the remainder of the payload.
    let atoms_start = payload.find("ATOMS=(")?;
    let after_paren = &payload[atoms_start + "ATOMS=(".len()..];
    let close = after_paren.find(')')?;
    let inner = &after_paren[..close];

    // First number is the count; the rest are 1-based V3000 atom indices.
    let mut nums = inner.split_whitespace();
    let _count: usize = nums.next()?.parse().ok()?;
    // Collect and deduplicate atom indices.
    // A malformed V3000 file may list the same atom twice in one group.
    // Deduplication here mirrors the defensive check in StereoGroup::new()
    // and makes the intent explicit at the parser level (RDKit PR #9258).
    let mut seen = std::collections::HashSet::new();
    let atom_indices: Vec<AtomIdx> = nums
        .filter_map(|s| {
            let v3k: u32 = s.parse().ok()?;
            resolve_atom_idx(v3k, atom_idx_map)
        })
        .filter(|idx| seen.insert(*idx))
        .collect();

    if atom_indices.is_empty() {
        return None;
    }

    Some(StereoGroup::new(kind, atom_indices))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{AtomIdx, BondOrder, Element};

    // -----------------------------------------------------------------------
    // Test fixtures
    // -----------------------------------------------------------------------

    const METHANE_V3K: &str = "\
methane
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";

    const ETHANOL_V3K: &str = "\
ethanol
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 C 1.5 0.0 0.0 0
M  V30 3 O 3.0 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
";

    // -----------------------------------------------------------------------
    // Test 1: methane — 1 atom, 0 bonds
    // -----------------------------------------------------------------------
    #[test]
    fn test_methane_counts() {
        let (mol, _) = parse_mol_v3000(METHANE_V3K).expect("parse methane");
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.bond_count(), 0);
    }

    // -----------------------------------------------------------------------
    // Test 2: ethanol — 3 atoms, 2 bonds
    // -----------------------------------------------------------------------
    #[test]
    fn test_ethanol_counts() {
        let (mol, _) = parse_mol_v3000(ETHANOL_V3K).expect("parse ethanol");
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
    }

    // -----------------------------------------------------------------------
    // Test 3: ethanol — bond between atoms 0 and 1 is Single
    // -----------------------------------------------------------------------
    #[test]
    fn test_ethanol_bond_0_1_single() {
        let (mol, _) = parse_mol_v3000(ETHANOL_V3K).expect("parse ethanol");
        let (_, bond) = mol
            .bond_between(AtomIdx(0), AtomIdx(1))
            .expect("bond 0-1 exists");
        assert_eq!(bond.order, BondOrder::Single);
    }

    // -----------------------------------------------------------------------
    // Test 4: ethanol — bond between atoms 1 and 2 is Single
    // -----------------------------------------------------------------------
    #[test]
    fn test_ethanol_bond_1_2_single() {
        let (mol, _) = parse_mol_v3000(ETHANOL_V3K).expect("parse ethanol");
        let (_, bond) = mol
            .bond_between(AtomIdx(1), AtomIdx(2))
            .expect("bond 1-2 exists");
        assert_eq!(bond.order, BondOrder::Single);
    }

    // -----------------------------------------------------------------------
    // Test 5: CHG=1 → atom.charge == 1
    // -----------------------------------------------------------------------
    #[test]
    fn test_positive_charge() {
        let mol_str = "\
charged_pos
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.0 0.0 0.0 0 CHG=1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
        let (mol, _) = parse_mol_v3000(mol_str).expect("parse charged_pos");
        assert_eq!(mol.atom(AtomIdx(0)).charge, 1);
    }

    // -----------------------------------------------------------------------
    // Test 6: CHG=-1 → atom.charge == -1
    // -----------------------------------------------------------------------
    #[test]
    fn test_negative_charge() {
        let mol_str = "\
charged_neg
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 0.0 0.0 0.0 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
        let (mol, _) = parse_mol_v3000(mol_str).expect("parse charged_neg");
        assert_eq!(mol.atom(AtomIdx(0)).charge, -1);
    }

    // -----------------------------------------------------------------------
    // Test 7: MASS=13 → atom.isotope == Some(13)
    // -----------------------------------------------------------------------
    #[test]
    fn test_isotope() {
        let mol_str = "\
isotope
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0 MASS=13
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
        let (mol, _) = parse_mol_v3000(mol_str).expect("parse isotope");
        assert_eq!(mol.atom(AtomIdx(0)).isotope, Some(13));
    }

    // -----------------------------------------------------------------------
    // Test 8: bond type=4 → BondOrder::Aromatic
    // -----------------------------------------------------------------------
    #[test]
    fn test_aromatic_bond() {
        let mol_str = "\
aromatic
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 C 1.5 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
";
        let (mol, _) = parse_mol_v3000(mol_str).expect("parse aromatic");
        let (_, bond) = mol
            .bond_between(AtomIdx(0), AtomIdx(1))
            .expect("bond exists");
        assert_eq!(bond.order, BondOrder::Aromatic);
    }

    // -----------------------------------------------------------------------
    // Test 9: bond type=2 → BondOrder::Double
    // -----------------------------------------------------------------------
    #[test]
    fn test_double_bond() {
        let mol_str = "\
double_bond
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 O 1.2 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
";
        let (mol, _) = parse_mol_v3000(mol_str).expect("parse double_bond");
        let (_, bond) = mol
            .bond_between(AtomIdx(0), AtomIdx(1))
            .expect("bond exists");
        assert_eq!(bond.order, BondOrder::Double);
    }

    // -----------------------------------------------------------------------
    // Test 10: metadata — name and comment parsed correctly
    // -----------------------------------------------------------------------
    #[test]
    fn test_metadata() {
        let mol_str = "\
my_molecule
  some_prog
my comment line
  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
        let (_, meta) = parse_mol_v3000(mol_str).expect("parse metadata");
        assert_eq!(meta.name, "my_molecule");
        assert_eq!(meta.comment, "my comment line");
    }

    // -----------------------------------------------------------------------
    // Test 11: line continuation (atom line split with trailing `-`)
    // -----------------------------------------------------------------------
    #[test]
    fn test_line_continuation() {
        // The atom line for atom 1 is split after the z-coordinate.
        // The HCOUNT keyword appears on the continuation line.
        let mol_str = "\
continuation
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0 MASS=12 -
M  V30 HCOUNT=3
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
        let (mol, _) = parse_mol_v3000(mol_str).expect("parse continuation");
        // atom 0 should have isotope 12 and hydrogen_count 3 (set via continuation)
        assert_eq!(mol.atom(AtomIdx(0)).element, Element::C);
        assert_eq!(mol.atom(AtomIdx(0)).isotope, Some(12));
        assert_eq!(mol.atom(AtomIdx(0)).hydrogen_count, Some(3));
    }

    // -----------------------------------------------------------------------
    // Test 12: missing END ATOM → returns Err
    // -----------------------------------------------------------------------
    #[test]
    fn test_missing_end_atom_is_error() {
        let mol_str = "\
bad
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  END
";
        let result = parse_mol_v3000(mol_str);
        assert!(
            matches!(result, Err(MolParseError::V3000ParseError { .. })),
            "expected V3000ParseError but got a different result"
        );
    }

    // -----------------------------------------------------------------------
    // Test 13: elements parsed correctly for ethanol (C, C, O)
    // -----------------------------------------------------------------------
    #[test]
    fn test_ethanol_elements() {
        let (mol, _) = parse_mol_v3000(ETHANOL_V3K).expect("parse ethanol");
        let atoms: Vec<_> = mol.atoms().collect();
        assert_eq!(atoms[0].1.element, Element::C);
        assert_eq!(atoms[1].1.element, Element::C);
        assert_eq!(atoms[2].1.element, Element::O);
    }

    // -----------------------------------------------------------------------
    // Test 14: triple bond (type=3 → BondOrder::Triple)
    // -----------------------------------------------------------------------
    #[test]
    fn test_triple_bond() {
        let mol_str = "\
triple_bond
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 N 1.2 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 3 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
";
        let (mol, _) = parse_mol_v3000(mol_str).expect("parse triple_bond");
        let (_, bond) = mol
            .bond_between(AtomIdx(0), AtomIdx(1))
            .expect("bond exists");
        assert_eq!(bond.order, BondOrder::Triple);
    }

    #[test]
    fn test_query_bond_types_preserved() {
        let mol_str = "\
query_bonds
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 C 1.0 0.0 0.0 0
M  V30 3 C 2.0 0.0 0.0 0
M  V30 4 C 3.0 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 5 1 2
M  V30 2 8 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
";
        let (mol, meta) = parse_mol_v3000(mol_str).expect("parse query bonds");
        let bonds: Vec<_> = mol.bonds().collect();
        assert_eq!(bonds[0].1.order, BondOrder::QuerySingleOrDouble);
        assert_eq!(bonds[1].1.order, BondOrder::QueryAny);

        let written = write_mol_v3000(&mol, &meta, &[]);
        assert!(written.contains("M  V30 1 5 1 2"), "{written}");
        assert!(written.contains("M  V30 2 8 3 4"), "{written}");
    }

    // -----------------------------------------------------------------------
    // Test: MDL bond type 9 (dative/coordinate) round-trips, not silently
    // corrupted to Single -- regression test for the platinum coordination-
    // chemistry benchmark (validation/platinum/FEASIBILITY.md). The molblock
    // below is not hand-crafted: it is exactly what RDKit 2026.03.3 writes
    // for `Chem.MolFromSmiles` + `AddBond(n, pt, Chem.BondType.DATIVE)` +
    // `Chem.MolToMolBlock(mol, forceV3000=True)` -- RDKit auto-upgrades to
    // V3000 for any dative bond (it cannot express one in V2000 either).
    // Before this fix, bond type 9 fell through the reader's `_ =>
    // BondOrder::Single` catch-all with no error, no warning: a
    // structurally valid file silently produced a different molecule.
    // -----------------------------------------------------------------------
    #[test]
    fn test_dative_bond_type_9_round_trips() {
        // NOT built with a `"\` line-continuation: the real molblock's second
        // header line is blank-then-indented (RDKit writes no molecule name),
        // and a leading `\`-continuation would silently swallow that blank
        // line along with its indentation (rustc even warns "multiple lines
        // skipped by escaped newline" if written that way) -- explicit `\n`
        // keeps this byte-for-byte identical to RDKit's own output.
        let mol_str = "\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 N 0.000000 0.000000 0.000000 0\nM  V30 2 Pt 1.299038 0.750000 0.000000 0 VAL=2\nM  V30 3 Cl 2.598076 1.500000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 9 1 2\nM  V30 2 1 2 3\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
        let (mol, meta) = parse_mol_v3000(mol_str).expect("parse RDKit dative molblock");
        let bonds: Vec<_> = mol.bonds().collect();
        assert_eq!(
            bonds[0].1.order,
            BondOrder::Dative,
            "bond type 9 must read as Dative, not silently fall through to Single"
        );
        // atom1/atom2 preserve the file's donor/acceptor order (N -> Pt).
        assert_eq!(mol.atom(bonds[0].1.atom1).element, Element::N);
        assert_eq!(mol.atom(bonds[0].1.atom2).element, Element::PT);
        assert_eq!(bonds[1].1.order, BondOrder::Single);

        let written = write_mol_v3000(&mol, &meta, &[]);
        assert!(
            written.contains("M  V30 1 9 1 2"),
            "Dative must write back out as bond type 9: {written}"
        );

        // Full round trip: re-parsing the freshly-written molblock still
        // gives a Dative bond with the same donor/acceptor order.
        let (rt_mol, _) = parse_mol_v3000(&written).expect("re-parse written molblock");
        let rt_bonds: Vec<_> = rt_mol.bonds().collect();
        assert_eq!(rt_bonds[0].1.order, BondOrder::Dative);
        assert_eq!(rt_mol.atom(rt_bonds[0].1.atom1).element, Element::N);
        assert_eq!(rt_mol.atom(rt_bonds[0].1.atom2).element, Element::PT);
    }

    // -----------------------------------------------------------------------
    // Test 15: atom-map number stored when nonzero
    // -----------------------------------------------------------------------
    #[test]
    fn test_atom_map() {
        let mol_str = "\
atommapped
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 3
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
        let (mol, _) = parse_mol_v3000(mol_str).expect("parse atommapped");
        assert_eq!(mol.atom(AtomIdx(0)).atom_map, Some(3));
    }
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Serialise `mol` to a MOL V3000 (Extended Ctab) string.
///
/// `coords[i]` is the `(x, y)` position for atom `i`.  Atoms beyond
/// `coords.len()` receive `(0.0, 0.0, 0.0)`.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo.** This is a
/// 2D-only writer with no z channel; MOL/CTfile has no other field for a
/// non-tetrahedral stereo tag either (see
/// `docs/rfcs/square_planar_mol_io_rfc.md`). A square-planar-tagged atom is
/// written with no indication anything was dropped. If `mol` may carry
/// `Chirality::SquarePlanar`, use
/// [`write_mol_v3000_with_conformer_checked`] instead, which fails closed
/// with a typed error rather than silently discarding the tag.
pub fn write_mol_v3000(mol: &Molecule, metadata: &MolMetadata, coords: &[(f64, f64)]) -> String {
    let natoms = mol.atom_count();
    let nbonds = mol.bond_count();

    let mut out = String::new();

    // Header (same 3-line format as V2000)
    out.push_str(&metadata.name);
    out.push('\n');
    out.push_str("  chematic\n");
    out.push_str(&metadata.comment);
    out.push('\n');

    // Counts line with V3000 tag (no atom/bond counts — they go in M  V30 COUNTS)
    out.push_str("  0  0  0  0  0  0  0  0  0  0999 V3000\n");

    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!("M  V30 COUNTS {natoms} {nbonds} 0 0 0\n"));

    // Atom block
    out.push_str("M  V30 BEGIN ATOM\n");
    for (idx, atom) in mol.atoms() {
        let (x, y) = coords.get(idx.0 as usize).copied().unwrap_or((0.0, 0.0));
        let sym = atom.element.symbol();
        let atom_map = atom.atom_map.unwrap_or(0);
        let i = idx.0 + 1; // 1-based

        let mut line = format!("M  V30 {i} {sym} {x:.4} {y:.4} 0.0000 {atom_map}");
        if atom.charge != 0 {
            line.push_str(&format!(" CHG={}", atom.charge));
        }
        if let Some(iso) = atom.isotope {
            line.push_str(&format!(" MASS={iso}"));
        }
        if let Some(h) = atom.hydrogen_count {
            line.push_str(&format!(" HCOUNT={h}"));
        }
        out.push_str(&line);
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");

    // Bond block
    out.push_str("M  V30 BEGIN BOND\n");
    for (bidx, bond) in mol.bonds() {
        let a1 = bond.atom1.0 + 1;
        let a2 = bond.atom2.0 + 1;
        let order = match bond.order {
            BondOrder::Zero => 0,
            BondOrder::Single | BondOrder::Up | BondOrder::Down => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            BondOrder::QuerySingleOrDouble => 5,
            BondOrder::QuerySingleOrAromatic => 6,
            BondOrder::QueryDoubleOrAromatic => 7,
            BondOrder::QueryAny => 8,
            // Matches RDKit's own V3000 convention for `Bond::BondType::DATIVE`
            // (see the reader's matching case above); `atom1`/`atom2` are
            // already in donor/acceptor order. V2000's writer (mol2000.rs)
            // still collapses `Dative` to a plain single bond, since RDKit
            // itself cannot express a dative bond in V2000 either -- full
            // dative round-tripping through chematic's own MOL writer
            // requires V3000, same as RDKit.
            BondOrder::Dative => 9,
            BondOrder::Quadruple => 4,
        };
        let i = bidx.0 + 1;
        // V3000 bond CFG: 1=Up, 3=Down (NOT V2000's stereo-field codes 1/6 --
        // `CFG=6` is not a valid V3000 value).
        let stereo = match bond.order {
            BondOrder::Up => " CFG=1",
            BondOrder::Down => " CFG=3",
            _ => "",
        };
        out.push_str(&format!("M  V30 {i} {order} {a1} {a2}{stereo}\n"));
    }
    out.push_str("M  V30 END BOND\n");

    // Optional COLLECTION block for enhanced stereo groups.
    let groups = mol.stereo_groups();
    if !groups.is_empty() {
        out.push_str("M  V30 BEGIN COLLECTION\n");
        for group in groups {
            let key = match &group.kind {
                StereoGroupKind::Absolute => "MDLV30/STEABS".to_string(),
                StereoGroupKind::Or(n) => format!("MDLV30/STEOR{n}"),
                StereoGroupKind::And(n) => format!("MDLV30/STEAND{n}"),
            };
            let n = group.atom_indices.len();
            let idxs: Vec<String> = group
                .atom_indices
                .iter()
                .map(|ai| (ai.0 + 1).to_string()) // 0-based → 1-based
                .collect();
            out.push_str(&format!("M  V30 {key} ATOMS=({n} {})\n", idxs.join(" ")));
        }
        out.push_str("M  V30 END COLLECTION\n");
    }

    out.push_str("M  V30 END CTAB\n");
    out.push_str("M  END\n");

    out
}

/// Serialize `mol` to MOL V3000 (Extended Ctab) format using `conformer`'s
/// real 3D coordinates, stamping the header's line-2 dimensional code as
/// `3D` -- the V3000 counterpart of
/// [`crate::mol2000::write_mol_with_conformer`]. As with that function, no
/// wedge/hash `CFG` is ever emitted on a bond line here: a real 3D geometry
/// makes a 2D wedge symbol redundant at best, and round-tripping this output
/// back through [`read_mol_v3000_with_diagnostics`] must not manufacture a
/// fresh [`crate::mol2000::Stereo3DDiagnostic::WedgeVs3DParityConflict`] on
/// its own output. Enhanced stereo groups (`COLLECTION`/`STEABS`/`STEOR`/
/// `STEAND`) are unaffected -- they label which atoms form a stereo group,
/// not a direction, and remain meaningful for a 3D record.
///
/// **Does not validate `Chirality::SquarePlanar` stereo against
/// `conformer`** -- it writes whatever coordinates it is given, trusting the
/// caller. If `conformer`'s geometry doesn't actually match a declared
/// square-planar tag, this will silently write a self-inconsistent file.
/// Use [`write_mol_v3000_with_conformer_checked`] instead to fail closed on
/// that mismatch (or on a missing/flat conformer) rather than trusting it.
pub fn write_mol_v3000_with_conformer(
    mol: &Molecule,
    metadata: &MolMetadata,
    conformer: &Coords3D,
) -> String {
    let natoms = mol.atom_count();
    let nbonds = mol.bond_count();

    let mut out = String::new();

    out.push_str(&metadata.name);
    out.push('\n');
    // Same column convention as `mol2000::write_mol_with_conformer`.
    out.push_str("  chematic          3D\n");
    out.push_str(&metadata.comment);
    out.push('\n');

    out.push_str("  0  0  0  0  0  0  0  0  0  0999 V3000\n");

    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!("M  V30 COUNTS {natoms} {nbonds} 0 0 0\n"));

    out.push_str("M  V30 BEGIN ATOM\n");
    for (idx, atom) in mol.atoms() {
        let p = conformer
            .points
            .get(idx.0 as usize)
            .copied()
            .unwrap_or(Point3::zero());
        let sym = atom.element.symbol();
        let atom_map = atom.atom_map.unwrap_or(0);
        let i = idx.0 + 1;

        let mut line = format!(
            "M  V30 {i} {sym} {:.4} {:.4} {:.4} {atom_map}",
            p.x, p.y, p.z
        );
        if atom.charge != 0 {
            line.push_str(&format!(" CHG={}", atom.charge));
        }
        if let Some(iso) = atom.isotope {
            line.push_str(&format!(" MASS={iso}"));
        }
        if let Some(h) = atom.hydrogen_count {
            line.push_str(&format!(" HCOUNT={h}"));
        }
        out.push_str(&line);
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");

    out.push_str("M  V30 BEGIN BOND\n");
    for (bidx, bond) in mol.bonds() {
        let a1 = bond.atom1.0 + 1;
        let a2 = bond.atom2.0 + 1;
        let order = match bond.order {
            BondOrder::Zero => 0,
            BondOrder::Single | BondOrder::Up | BondOrder::Down => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            BondOrder::QuerySingleOrDouble => 5,
            BondOrder::QuerySingleOrAromatic => 6,
            BondOrder::QueryDoubleOrAromatic => 7,
            BondOrder::QueryAny => 8,
            // Matches RDKit's own V3000 convention for `Bond::BondType::DATIVE`
            // (see the reader's matching case above); `atom1`/`atom2` are
            // already in donor/acceptor order. V2000's writer (mol2000.rs)
            // still collapses `Dative` to a plain single bond, since RDKit
            // itself cannot express a dative bond in V2000 either -- full
            // dative round-tripping through chematic's own MOL writer
            // requires V3000, same as RDKit.
            BondOrder::Dative => 9,
            BondOrder::Quadruple => 4,
        };
        let i = bidx.0 + 1;
        // No wedge/hash CFG here -- see doc comment above.
        out.push_str(&format!("M  V30 {i} {order} {a1} {a2}\n"));
    }
    out.push_str("M  V30 END BOND\n");

    let groups = mol.stereo_groups();
    if !groups.is_empty() {
        out.push_str("M  V30 BEGIN COLLECTION\n");
        for group in groups {
            let key = match &group.kind {
                StereoGroupKind::Absolute => "MDLV30/STEABS".to_string(),
                StereoGroupKind::Or(n) => format!("MDLV30/STEOR{n}"),
                StereoGroupKind::And(n) => format!("MDLV30/STEAND{n}"),
            };
            let n = group.atom_indices.len();
            let idxs: Vec<String> = group
                .atom_indices
                .iter()
                .map(|ai| (ai.0 + 1).to_string())
                .collect();
            out.push_str(&format!("M  V30 {key} ATOMS=({n} {})\n", idxs.join(" ")));
        }
        out.push_str("M  V30 END COLLECTION\n");
    }

    out.push_str("M  V30 END CTAB\n");
    out.push_str("M  END\n");

    out
}

/// [`write_mol_v3000_with_conformer`], but fails closed with a typed
/// [`crate::mol2000::MolStereoWriteError`] instead of silently writing
/// coordinates that don't actually match a molecule's declared
/// square-planar stereo (or don't exist at all) -- the V3000 counterpart of
/// [`crate::mol2000::write_mol_with_conformer_checked`]. See
/// [`crate::mol2000::validate_square_planar_for_write`].
pub fn write_mol_v3000_with_conformer_checked(
    mol: &Molecule,
    metadata: &MolMetadata,
    conformer: &Coords3D,
) -> Result<String, crate::mol2000::MolStereoWriteError> {
    crate::mol2000::validate_square_planar_for_write(
        mol,
        Some(conformer),
        crate::mol2000::MolFormat::V3000,
    )?;
    Ok(write_mol_v3000_with_conformer(mol, metadata, conformer))
}

#[cfg(test)]
mod write_tests {
    use super::*;
    use crate::mol3000::parse_mol_v3000;
    use chematic_core::Element;
    use chematic_core::{Atom, MoleculeBuilder};

    fn ethanol() -> Molecule {
        use chematic_core::BondOrder;
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        let c2 = b.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        let o = b.add_atom(Atom::new(Element::from_symbol("O").unwrap()));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, o, BondOrder::Single).unwrap();
        b.build()
    }

    #[test]
    fn write_v3000_roundtrip_ethanol() {
        let mol = ethanol();
        let meta = MolMetadata {
            name: "ethanol".into(),
            comment: String::new(),
        };
        let v3k = write_mol_v3000(&mol, &meta, &[]);
        let (mol2, meta2) = parse_mol_v3000(&v3k).expect("round-trip parse");
        assert_eq!(mol.atom_count(), mol2.atom_count());
        assert_eq!(mol.bond_count(), mol2.bond_count());
        assert_eq!(meta2.name, "ethanol");
    }

    #[test]
    fn write_v3000_contains_v3000_tag() {
        let mol = ethanol();
        let meta = MolMetadata::default();
        let v3k = write_mol_v3000(&mol, &meta, &[]);
        assert!(v3k.contains("V3000"), "output should contain V3000 tag");
        assert!(
            v3k.contains("M  V30 BEGIN CTAB"),
            "should contain CTAB block"
        );
    }

    #[test]
    fn write_v3000_stereo_group_roundtrip() {
        use chematic_core::{AtomIdx, BondOrder, StereoGroup, StereoGroupKind};

        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        let c2 = b.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        let c3 = b.add_atom(Atom::new(Element::from_symbol("N").unwrap()));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        let mut mol = b.build();

        // Absolute group on atom 0; OR group on atoms 1 and 2.
        mol.set_stereo_groups(vec![
            StereoGroup::new(StereoGroupKind::Absolute, vec![c1]),
            StereoGroup::new(StereoGroupKind::Or(1), vec![c2, c3]),
        ]);

        let meta = MolMetadata {
            name: "stereo_test".into(),
            comment: String::new(),
        };
        let v3k = write_mol_v3000(&mol, &meta, &[]);

        // Verify COLLECTION block is present.
        assert!(
            v3k.contains("BEGIN COLLECTION"),
            "should have COLLECTION block"
        );
        assert!(v3k.contains("MDLV30/STEABS"), "should have STEABS entry");
        assert!(v3k.contains("MDLV30/STEOR1"), "should have STEOR1 entry");

        // Round-trip: parse back and verify stereo groups are preserved.
        let (mol2, _) = parse_mol_v3000(&v3k).expect("round-trip parse");
        assert_eq!(
            mol2.stereo_groups().len(),
            2,
            "should have 2 stereo groups after round-trip"
        );

        let abs_group = mol2
            .stereo_groups()
            .iter()
            .find(|g| g.kind == StereoGroupKind::Absolute)
            .expect("Absolute group should exist");
        assert_eq!(abs_group.atom_indices, vec![AtomIdx(0)]);

        let or_group = mol2
            .stereo_groups()
            .iter()
            .find(|g| g.kind == StereoGroupKind::Or(1))
            .expect("Or(1) group should exist");
        assert_eq!(or_group.atom_indices, vec![AtomIdx(1), AtomIdx(2)]);
    }

    #[test]
    fn test_stereo_group_duplicate_atoms_deduplicated() {
        // RDKit PR #9258: duplicate atom indices in a StereoGroup must be silently
        // removed. In RDKit this caused heap-use-after-free; in chematic the
        // duplicate indices produce incorrect stereo group membership.
        //
        // Craft a minimal V3000 MOL block with STEABS ATOMS=(3 1 1 2)
        // — atom 1 appears twice, which should be collapsed to [AtomIdx(0)].
        let v3k = "\n\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\n\
M  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0 CFG=2\n\
M  V30 2 C 1 0 0 0 CFG=2\n\
M  V30 3 C 2 0 0 0\n\
M  V30 END ATOM\nM  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 2 1 2 3\n\
M  V30 END BOND\nM  V30 BEGIN COLLECTION\n\
M  V30 MDLV30/STEABS ATOMS=(3 1 1 2)\n\
M  V30 END COLLECTION\nM  V30 END CTAB\nM  END\n";

        let (mol, _) = parse_mol_v3000(v3k).expect("should parse despite duplicate atom");
        assert_eq!(mol.stereo_groups().len(), 1);
        let group = &mol.stereo_groups()[0];
        assert_eq!(
            group.kind,
            StereoGroupKind::Absolute,
            "group kind must be Absolute"
        );
        // Atom 1 should appear only ONCE (deduplicated), plus atom 2.
        assert_eq!(
            group.atom_indices.len(),
            2,
            "duplicate atom index must be removed: got {:?}",
            group.atom_indices
        );
        assert!(
            group.atom_indices.contains(&AtomIdx(0)),
            "atom 0 must be in group"
        );
        assert!(
            group.atom_indices.contains(&AtomIdx(1)),
            "atom 1 must be in group"
        );
    }
}
