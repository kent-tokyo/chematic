//! MOL V3000 (Extended Ctab) parser.
//!
//! Reference: MDL/Dassault Systèmes CTfile Formats specification, V3000 section.
//!
//! MOL V3000 (Extended Ctab) parser and writer.
//!
//! Reference: MDL/Dassault Systèmes CTfile Formats specification, V3000 section.

use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Coords3D, Element, Molecule, MoleculeBuilder, Point3,
    RGroupLabel, StereoGroup, StereoGroupKind,
};
use chematic_perception::{
    apply_ctab_local_parity_from_wedges_with_diagnostics, apply_ez_directions_from_2d_ex,
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

fn parse_rgroup_number(tokens: &[&str]) -> Result<Option<u16>, String> {
    let joined = tokens.join(" ");
    let Some(rest) = joined.split("RGROUPS=(").nth(1) else {
        return Ok(None);
    };
    let Some(body) = rest.split(')').next() else {
        return Err("RGROUPS property is missing ')'".to_string());
    };
    let values = body.split_whitespace().collect::<Vec<_>>();
    let count = values
        .first()
        .and_then(|value| value.parse::<usize>().ok())
        .ok_or_else(|| "RGROUPS property lacks a valid count".to_string())?;
    if count != 1 || values.len() != 2 {
        return Err("chematic supports exactly one R-group number per atom".to_string());
    }
    let number = values[1]
        .parse::<u16>()
        .ok()
        .filter(|number| (1..=9999).contains(number))
        .ok_or_else(|| "R-group number must be in 1..=9999".to_string())?;
    Ok(Some(number))
}

fn opaque_atom_properties(tokens: &[&str]) -> String {
    let mut opaque = Vec::new();
    let mut skipping_rgroups = false;
    for token in tokens {
        if skipping_rgroups {
            if token.ends_with(')') {
                skipping_rgroups = false;
            }
            continue;
        }
        if token.starts_with("RGROUPS=(") {
            skipping_rgroups = !token.ends_with(')');
            continue;
        }
        if token.starts_with("CHG=") || token.starts_with("MASS=") || token.starts_with("HCOUNT=") {
            continue;
        }
        opaque.push(*token);
    }
    opaque.join(" ")
}

/// Parse one V3000 atom record and append it to the parser's parallel state.
fn parse_v3000_atom_line(
    tokens: &[&str],
    line: usize,
    builder: &mut MoleculeBuilder,
    atom_idx_map: &mut Vec<(u32, AtomIdx)>,
    coords: &mut Vec<(f64, f64)>,
    raw_z: &mut Vec<f64>,
    atom_properties: &mut Vec<(u32, String)>,
) -> Result<(), MolParseError> {
    // Atom line: <idx> <symbol> <x> <y> <z> <aamap> [KEY=VAL ...]
    if tokens.len() < 6 {
        return Err(MolParseError::InvalidAtomLine {
            line,
            detail: format!(
                "V3000 atom line needs at least 6 fields, got {}",
                tokens.len()
            ),
        });
    }

    let v3k_idx = tokens[0]
        .parse::<u32>()
        .map_err(|_| MolParseError::InvalidAtomLine {
            line,
            detail: format!("cannot parse atom index from '{}'", tokens[0]),
        })?;

    // Strip bracket notation e.g. "[OH]" -> "OH".
    let sym = tokens[1].trim_start_matches('[').trim_end_matches(']');
    let kv_tokens = tokens.get(6..).unwrap_or(&[]);
    let rgroup_number = parse_rgroup_number(kv_tokens)
        .map_err(|detail| MolParseError::InvalidAtomLine { line, detail })?;
    // V3000 writers, including Indigo, may spell deuterium as `D` rather
    // than hydrogen with `MASS=2`. Internally retain the standard element
    // representation plus its isotope so the record is usable everywhere.
    let deuterium = sym == "D";
    let element_symbol = if deuterium { "H" } else { sym };

    // Coordinates are part of the serialized graph: malformed or non-finite
    // values must not silently turn into a different structure.
    let x: f64 = tokens[2]
        .parse()
        .map_err(|_| MolParseError::InvalidAtomLine {
            line,
            detail: format!("cannot parse x coordinate from '{}'", tokens[2]),
        })?;
    let y: f64 = tokens[3]
        .parse()
        .map_err(|_| MolParseError::InvalidAtomLine {
            line,
            detail: format!("cannot parse y coordinate from '{}'", tokens[3]),
        })?;
    if !x.is_finite() || !y.is_finite() {
        return Err(MolParseError::InvalidAtomLine {
            line,
            detail: "x/y coordinate is not finite (NaN/Infinite)".to_string(),
        });
    }

    // Unlike x/y, z is retained as 3D input and therefore must be a finite number.
    let z: f64 = tokens[4]
        .parse()
        .map_err(|_| MolParseError::InvalidAtomLine {
            line,
            detail: format!("cannot parse z coordinate from '{}'", tokens[4]),
        })?;
    if !z.is_finite() {
        return Err(MolParseError::InvalidAtomLine {
            line,
            detail: format!("z coordinate is not finite (NaN/Infinite): '{}'", tokens[4]),
        });
    }

    let aamap_raw = tokens[5].parse::<u16>().unwrap_or(0);
    let atom_map = (aamap_raw != 0).then_some(aamap_raw);
    let charge: i8 = parse_kv(kv_tokens, "CHG")
        .and_then(|v| v.parse::<i8>().ok())
        .unwrap_or(0);
    let isotope = parse_kv(kv_tokens, "MASS")
        .and_then(|v| v.parse::<u16>().ok())
        .or_else(|| deuterium.then_some(2));
    let hydrogen_count = parse_kv(kv_tokens, "HCOUNT").and_then(|v| {
        let n: i32 = v.parse().ok()?;
        if n < 0 { None } else { Some(n as u8) }
    });

    let (mut atom, r_group) = if sym == "*" {
        (Atom::wildcard(), None)
    } else if sym == "R" {
        (Atom::wildcard(), Some(RGroupLabel::unnumbered()))
    } else if sym == "R#" {
        match rgroup_number {
            Some(number) => (
                Atom::wildcard(),
                Some(RGroupLabel::numbered(number).expect("validated non-zero number")),
            ),
            None => (Atom::wildcard(), Some(RGroupLabel::unnumbered())),
        }
    } else if let Some(number) = sym.strip_prefix('R')
        && !number.is_empty()
        && number.bytes().all(|byte| byte.is_ascii_digit())
    {
        let label = number
            .parse::<u16>()
            .ok()
            .and_then(RGroupLabel::numbered)
            .ok_or_else(|| MolParseError::UnknownElement {
                symbol: sym.to_string(),
                line,
            })?;
        (Atom::wildcard(), Some(label))
    } else if let Some(number) = rgroup_number {
        (
            Atom::wildcard(),
            Some(RGroupLabel::numbered(number).expect("validated non-zero number")),
        )
    } else {
        let element =
            Element::from_symbol(element_symbol).ok_or_else(|| MolParseError::UnknownElement {
                symbol: sym.to_string(),
                line,
            })?;
        (Atom::new(element), None)
    };
    atom.charge = charge;
    atom.isotope = isotope;
    atom.hydrogen_count = hydrogen_count;
    atom.atom_map = atom_map;

    let builder_idx = builder.add_atom(atom);
    if let Some(label) = r_group {
        builder.set_r_group(builder_idx, label);
    }
    atom_idx_map.push((v3k_idx, builder_idx));
    coords.push((x, y));
    raw_z.push(z);
    // Keep attributes that do not have a lossless core representation opaque.
    // In particular, V3000 query constraints such as SUBST/RBCNT/UNSAT must
    // survive a V3000 -> V3000 round trip.  Conversion to a format that has no
    // matching query contract is rejected by the caller rather than silently
    // treating such constraints as an ordinary molecule.
    let opaque_properties = opaque_atom_properties(kv_tokens);
    if !opaque_properties.is_empty() {
        // Writers emit atoms in builder order. Retain attributes under that
        // output id so a valid non-contiguous source id cannot make them
        // silently disappear on the next write.
        atom_properties.push((builder_idx.0 + 1, opaque_properties));
    }
    Ok(())
}

/// Parse one V3000 bond record and append it to the molecule builder.
fn parse_v3000_bond_line(
    tokens: &[&str],
    line: usize,
    builder: &mut MoleculeBuilder,
    atom_idx_map: &[(u32, AtomIdx)],
    explicitly_unspecified_ez: &mut std::collections::HashSet<BondIdx>,
    bond_properties: &mut Vec<(u32, String)>,
) -> Result<(), MolParseError> {
    // Bond line: <idx> <type> <atom1> <atom2> [KEY=VAL ...]
    if tokens.len() < 4 {
        return Err(MolParseError::InvalidBondLine {
            line,
            detail: format!(
                "V3000 bond line needs at least 4 fields, got {}",
                tokens.len()
            ),
        });
    }

    let btype_raw = tokens[1]
        .parse::<u8>()
        .map_err(|_| MolParseError::InvalidBondLine {
            line,
            detail: format!("cannot parse bond type from '{}'", tokens[1]),
        })?;
    let a1_v3k = tokens[2]
        .parse::<u32>()
        .map_err(|_| MolParseError::InvalidBondLine {
            line,
            detail: format!("cannot parse atom1 index from '{}'", tokens[2]),
        })?;
    let a2_v3k = tokens[3]
        .parse::<u32>()
        .map_err(|_| MolParseError::InvalidBondLine {
            line,
            detail: format!("cannot parse atom2 index from '{}'", tokens[3]),
        })?;
    let a1 =
        resolve_atom_idx(a1_v3k, atom_idx_map).ok_or_else(|| MolParseError::InvalidBondLine {
            line,
            detail: format!("atom index {a1_v3k} not found in atom block"),
        })?;
    let a2 =
        resolve_atom_idx(a2_v3k, atom_idx_map).ok_or_else(|| MolParseError::InvalidBondLine {
            line,
            detail: format!("atom index {a2_v3k} not found in atom block"),
        })?;

    let kv_tokens = tokens.get(4..).unwrap_or(&[]);
    let order = match btype_raw {
        0 => BondOrder::Zero,
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
        9 => BondOrder::Dative,
        _ => BondOrder::Single,
    };

    let bidx = builder
        .add_bond(a1, a2, order)
        .map_err(|e| MolParseError::InvalidBondLine {
            line,
            detail: e.to_string(),
        })?;
    if btype_raw == 2 && parse_kv(kv_tokens, "CFG").as_deref() == Some("2") {
        explicitly_unspecified_ez.insert(bidx);
    }
    // Preserve richer V3000 bond attributes without pretending to interpret
    // them. CFG is represented by BondOrder and is regenerated by the
    // writer. ENDPTS is the one supported reference-bearing attribute: its
    // atom ids must follow the normalized output order when source ids are
    // non-contiguous. ATTACH and future keys remain opaque.
    let opaque_properties = kv_tokens
        .iter()
        .filter(|token| !token.starts_with("CFG="))
        .copied()
        .collect::<Vec<_>>()
        .join(" ");
    let opaque_properties = normalize_v3000_endpoints(&opaque_properties, atom_idx_map);
    if !opaque_properties.is_empty() {
        // As with atom metadata, normalize a source bond id to the writer's
        // sequential id. V3000 permits non-contiguous ids, whereas the
        // writer always emits `BondIdx + 1`.
        bond_properties.push((bidx.0 + 1, opaque_properties));
    }
    Ok(())
}

/// Convert an `ENDPTS=(n atom-id...)` property from source V3000 ids to the
/// sequential ids emitted by this writer. Leave malformed or unknown-shaped
/// values unchanged: they remain opaque metadata rather than becoming a
/// guessed semantic edit.
fn normalize_v3000_endpoints(properties: &str, atom_idx_map: &[(u32, AtomIdx)]) -> String {
    let Some(start) = properties.find("ENDPTS=(") else {
        return properties.to_owned();
    };
    let value_start = start + "ENDPTS=(".len();
    let Some(close_offset) = properties[value_start..].find(')') else {
        return properties.to_owned();
    };
    let value_end = value_start + close_offset;
    let values: Vec<&str> = properties[value_start..value_end]
        .split_whitespace()
        .collect();
    let Some((declared, atom_ids)) = values.split_first() else {
        return properties.to_owned();
    };
    let Ok(declared) = declared.parse::<usize>() else {
        return properties.to_owned();
    };
    if atom_ids.len() != declared {
        return properties.to_owned();
    }
    let Some(writer_ids) = atom_ids
        .iter()
        .map(|id| {
            id.parse::<u32>()
                .ok()
                .and_then(|source_id| resolve_atom_idx(source_id, atom_idx_map))
                .map(|idx| (idx.0 + 1).to_string())
        })
        .collect::<Option<Vec<_>>>()
    else {
        return properties.to_owned();
    };
    let normalized = format!("ENDPTS=({declared} {})", writer_ids.join(" "));
    format!(
        "{}{}{}",
        &properties[..start],
        normalized,
        &properties[value_end + 1..]
    )
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
    let terminal = all_lines
        .iter()
        .rev()
        .find(|(_, line)| !line.trim().is_empty())
        .ok_or(MolParseError::UnexpectedEnd)?;
    if terminal.1 != "M  END" {
        return Err(v3k_err(
            terminal.0,
            "V3000 MOL block must terminate with an exact 'M  END' line",
        ));
    }
    for &(line_num, line) in all_lines.iter().skip(4) {
        if line.trim().is_empty() || line.starts_with(V30_PREFIX) {
            continue;
        }
        if line == "M  END" && line_num == terminal.0 {
            continue;
        }
        return Err(v3k_err(
            line_num,
            "expected an 'M  V30 ' record or the terminal 'M  END' line",
        ));
    }

    // -- Header lines 1–3 ----------------------------------------------------

    let name = all_lines[0].1.to_string();
    // Line 2 (program/date) is mostly discarded, except for the
    // dimensional-code field (columns 20..22, "2D"/"3D") -- same column as
    // V2000's header (see `crate::mol2000::parse_dimension_code`); V3000
    // keeps the identical 3-line header before its own counts line.
    let line2_raw = all_lines[1].1;
    let comment = all_lines[2].1.to_string();

    let mut metadata = MolMetadata {
        name,
        comment,
        v3000_sgroups: Vec::new(),
        v3000_bond_properties: Vec::new(),
        v3000_atom_properties: Vec::new(),
    };

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
        InSgroup,
        Done,
    }

    let mut state = State::BeforeCtab;
    let mut expected_atoms: usize = 0;
    let mut expected_bonds: usize = 0;
    let mut parsed_bonds: usize = 0;
    let mut expected_sgroups: Option<usize> = None;
    let mut stereo_groups: Vec<StereoGroup> = Vec::new();
    let mut sgroups: Vec<String> = Vec::new();
    // Known SGROUP kinds have a typed atom-reference view. Keep their IDs
    // unique while parsing so an otherwise lossless round trip cannot retain
    // two incompatible definitions for the same group. Unknown kinds remain
    // opaque because this crate does not claim their semantics.
    let mut typed_sgroup_ids: std::collections::HashSet<u32> = std::collections::HashSet::new();
    let mut bond_properties: Vec<(u32, String)> = Vec::new();
    let mut atom_properties: Vec<(u32, String)> = Vec::new();
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
                    if tokens.len() < 4 {
                        return Err(v3k_err(lnum, "COUNTS line has fewer than 3 values"));
                    }
                    expected_atoms = tokens[1].parse::<usize>().map_err(|_| {
                        v3k_err(
                            lnum,
                            format!("cannot parse atom count from '{}'", tokens[1]),
                        )
                    })?;
                    expected_bonds = tokens[2].parse::<usize>().map_err(|_| {
                        v3k_err(
                            lnum,
                            format!("cannot parse bond count from '{}'", tokens[2]),
                        )
                    })?;
                    expected_sgroups = Some(tokens[3].parse::<usize>().map_err(|_| {
                        v3k_err(
                            lnum,
                            format!("cannot parse SGROUP count from '{}'", tokens[3]),
                        )
                    })?);
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

                parse_v3000_atom_line(
                    &tokens,
                    lnum,
                    &mut builder,
                    &mut atom_idx_map,
                    &mut coords,
                    &mut raw_z,
                    &mut atom_properties,
                )?;
            }

            State::AfterAtomBlock => {
                if is_marker(&tokens, "BEGIN", "BOND") {
                    state = State::InBondBlock;
                } else if is_marker(&tokens, "END", "CTAB") {
                    if expected_bonds != 0 {
                        return Err(v3k_err(
                            lnum,
                            format!(
                                "COUNTS declares {expected_bonds} bonds but the BOND block is absent"
                            ),
                        ));
                    }
                    state = State::Done;
                } else {
                    return Err(v3k_err(
                        lnum,
                        "expected BEGIN BOND or END CTAB after END ATOM",
                    ));
                }
            }

            State::InBondBlock => {
                if is_marker(&tokens, "END", "BOND") {
                    if parsed_bonds != expected_bonds {
                        return Err(v3k_err(
                            lnum,
                            format!(
                                "COUNTS declares {expected_bonds} bonds but found {}",
                                parsed_bonds
                            ),
                        ));
                    }
                    state = State::AfterBondBlock;
                    continue;
                }

                parse_v3000_bond_line(
                    &tokens,
                    lnum,
                    &mut builder,
                    &atom_idx_map,
                    &mut explicitly_unspecified_ez,
                    &mut bond_properties,
                )?;
                parsed_bonds += 1;
            }

            State::AfterBondBlock => {
                if is_marker(&tokens, "BEGIN", "COLLECTION") {
                    state = State::InCollection;
                } else if is_marker(&tokens, "BEGIN", "SGROUP") {
                    state = State::InSgroup;
                } else if is_marker(&tokens, "END", "CTAB") {
                    state = State::Done;
                }
            }

            State::InCollection => {
                if is_marker(&tokens, "END", "COLLECTION") {
                    state = State::AfterBondBlock;
                } else {
                    match parse_stereo_group_line(payload, &atom_idx_map) {
                        Ok(Some(group)) => stereo_groups.push(group),
                        Ok(None) => {}
                        Err(detail) => return Err(v3k_err(lnum, detail)),
                    }
                }
            }

            State::InSgroup => {
                if is_marker(&tokens, "END", "SGROUP") {
                    state = State::AfterBondBlock;
                } else {
                    validate_known_sgroup_references(payload, &atom_idx_map, &mut typed_sgroup_ids)
                        .map_err(|detail| v3k_err(lnum, detail))?;
                    // Preserve the logical payload without pretending to
                    // understand polymer/query semantics. The writer restores
                    // the standard V30 prefix and emits this block before
                    // COLLECTION, matching the interoperability ordering used
                    // by current reference writers.
                    sgroups.push(payload.to_string());
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

    if let Some(expected_sgroups) = expected_sgroups
        && sgroups.len() != expected_sgroups
    {
        return Err(MolParseError::V3000ParseError {
            line: 0,
            msg: format!(
                "COUNTS declares {expected_sgroups} SGROUP records but found {}",
                sgroups.len()
            ),
        });
    }

    let mut mol = builder.build();
    crate::mol2000::flag_aromatic_bond_atoms(&mut mol);
    metadata.v3000_sgroups = sgroups;
    metadata.v3000_bond_properties = bond_properties;
    metadata.v3000_atom_properties = atom_properties;
    if !stereo_groups.is_empty() {
        mol.set_stereo_groups(stereo_groups);
    }
    // Tetrahedral parity first, then E/Z direction (side channel only) --
    // same ordering rationale as `mol2000.rs`.
    let stereo_diagnostics =
        apply_ctab_local_parity_from_wedges_with_diagnostics(&mut mol, &coords);
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

/// Validate references for the SGROUP kinds with a declared typed view.
///
/// Other SGROUP records deliberately remain opaque: accepting a new vendor
/// extension must not turn it into a guessed polymer or query interpretation.
/// Known records, however, must not retain dangling atom ids or duplicate
/// group ids, because either condition makes the typed projection ambiguous.
fn validate_known_sgroup_references(
    payload: &str,
    atom_idx_map: &[(u32, AtomIdx)],
    seen_ids: &mut std::collections::HashSet<u32>,
) -> Result<(), String> {
    let Some(kind) = payload.split_whitespace().nth(1) else {
        return Ok(());
    };
    if !matches!(kind, "SUP" | "GEN" | "COP" | "DAT" | "EXT") {
        return Ok(());
    }

    let group = parse_v3000_sgroup_line(payload).map_err(|error| error.to_string())?;
    if !seen_ids.insert(group.id) {
        return Err(format!("duplicate typed SGROUP id {}", group.id));
    }
    for atom_id in group.atom_ids {
        if resolve_atom_idx(atom_id, atom_idx_map).is_none() {
            return Err(format!("SGROUP ATOMS references unknown atom id {atom_id}"));
        }
    }
    Ok(())
}

/// Look up the builder `AtomIdx` for a given V3000 1-based index.
fn resolve_atom_idx(v3k_idx: u32, map: &[(u32, AtomIdx)]) -> Option<AtomIdx> {
    map.iter().find(|&&(k, _)| k == v3k_idx).map(|&(_, v)| v)
}

/// The semantic type token of a V3000 SGROUP record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum V3000SGroupKind {
    Sup,
    Gen,
    Cop,
    Dat,
    Ext,
    Other(String),
}

/// A validated, typed view of one V3000 SGROUP logical line.
///
/// The original logical line remains in [`MolMetadata::v3000_sgroups`] for
/// lossless writing. Atom references stay as the file's 1-based V3000 IDs;
/// this view does not claim polymer expansion semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct V3000SGroup {
    pub id: u32,
    pub kind: V3000SGroupKind,
    pub parent_id: Option<u32>,
    pub atom_ids: Vec<u32>,
    pub attributes: Vec<(String, String)>,
}

/// Parse and validate one V3000 SGROUP logical line.
pub fn parse_v3000_sgroup_line(line: &str) -> Result<V3000SGroup, MolParseError> {
    // V3000 list values contain spaces, e.g. `ATOMS=(2 1 2)` and
    // `BRKXYZ=(4 1 2 3 4)`, so whitespace is a separator only outside
    // balanced parentheses.
    let mut tokens = Vec::new();
    let mut start = None;
    let mut depth = 0usize;
    for (offset, ch) in line.char_indices() {
        if start.is_none() && !ch.is_whitespace() {
            start = Some(offset);
        }
        match ch {
            '(' => depth += 1,
            ')' => depth = depth.saturating_sub(1),
            _ => {}
        }
        if ch.is_whitespace()
            && depth == 0
            && let Some(begin) = start.take()
        {
            tokens.push(&line[begin..offset]);
        }
    }
    if let Some(begin) = start {
        tokens.push(&line[begin..]);
    }
    if tokens.len() < 3 {
        return Err(v3k_err(0, "SGROUP line needs id, type, and parent"));
    }
    let id = tokens[0]
        .parse::<u32>()
        .map_err(|_| v3k_err(0, "SGROUP id is not a non-negative integer"))?;
    if id == 0 {
        return Err(v3k_err(0, "SGROUP id must be non-zero"));
    }
    let kind = match tokens[1] {
        "SUP" => V3000SGroupKind::Sup,
        "GEN" => V3000SGroupKind::Gen,
        "COP" => V3000SGroupKind::Cop,
        "DAT" => V3000SGroupKind::Dat,
        "EXT" => V3000SGroupKind::Ext,
        other => V3000SGroupKind::Other(other.to_owned()),
    };
    let parent_raw = tokens[2]
        .parse::<u32>()
        .map_err(|_| v3k_err(0, "SGROUP parent is not a non-negative integer"))?;
    let atoms_start = line
        .find("ATOMS=(")
        .ok_or_else(|| v3k_err(0, "SGROUP line is missing ATOMS=(...)"))?;
    let atoms_value_start = atoms_start + "ATOMS=(".len();
    let close = line[atoms_value_start..]
        .find(')')
        .map(|offset| atoms_value_start + offset)
        .ok_or_else(|| v3k_err(0, "SGROUP ATOMS value is not closed"))?;
    let atoms_inner = &line[atoms_value_start..close];
    let mut atom_values = atoms_inner.split_whitespace();
    let declared = atom_values
        .next()
        .ok_or_else(|| v3k_err(0, "SGROUP ATOMS value is missing its count"))?
        .parse::<usize>()
        .map_err(|_| v3k_err(0, "SGROUP ATOMS count is not an integer"))?;
    let atom_ids: Result<Vec<u32>, MolParseError> = atom_values
        .map(|value| {
            let id = value
                .parse::<u32>()
                .map_err(|_| v3k_err(0, "SGROUP ATOMS contains a non-integer atom id"))?;
            if id == 0 {
                return Err(v3k_err(0, "SGROUP ATOMS contains atom id zero"));
            }
            Ok(id)
        })
        .collect();
    let atom_ids = atom_ids?;
    if atom_ids.len() != declared {
        return Err(v3k_err(
            0,
            format!(
                "SGROUP ATOMS declares {declared} ids but contains {}",
                atom_ids.len()
            ),
        ));
    }
    let attributes = tokens[3..]
        .iter()
        .filter_map(|token| token.split_once('='))
        .filter(|(key, _)| *key != "ATOMS")
        .map(|(key, value)| (key.to_owned(), value.to_owned()))
        .collect();
    Ok(V3000SGroup {
        id,
        kind,
        parent_id: (parent_raw != 0).then_some(parent_raw),
        atom_ids,
        attributes,
    })
}

/// Parse a V3000 COLLECTION line into a [`StereoGroup`].
///
/// Expected formats:
/// - `MDLV30/STEABS ATOMS=(3 1 2 3)` → `StereoGroupKind::Absolute`
/// - `MDLV30/STEREL1 ATOMS=(1 4)`    → `StereoGroupKind::Or(1)`
/// - `MDLV30/STERAC2 ATOMS=(1 6)`    → `StereoGroupKind::And(2)`
///
/// `STEOR` and `STEAND` are also accepted as legacy input spellings, but the
/// writer uses the MDL tokens emitted by current RDKit.
fn parse_stereo_group_line(
    payload: &str,
    atom_idx_map: &[(u32, AtomIdx)],
) -> Result<Option<StereoGroup>, String> {
    // First token is the group kind key.
    let Some(first_tok) = payload.split_whitespace().next() else {
        return Ok(None);
    };

    let kind = if first_tok == "MDLV30/STEABS" {
        StereoGroupKind::Absolute
    } else if let Some(n_str) = first_tok
        .strip_prefix("MDLV30/STEREL")
        .or_else(|| first_tok.strip_prefix("MDLV30/STEOR"))
    {
        let n: u32 = n_str
            .parse()
            .map_err(|_| "COLLECTION OR group id is not an integer".to_string())?;
        StereoGroupKind::Or(n)
    } else if let Some(n_str) = first_tok
        .strip_prefix("MDLV30/STERAC")
        .or_else(|| first_tok.strip_prefix("MDLV30/STEAND"))
    {
        let n: u32 = n_str
            .parse()
            .map_err(|_| "COLLECTION AND group id is not an integer".to_string())?;
        StereoGroupKind::And(n)
    } else {
        return Ok(None);
    };

    // Extract the ATOMS=(...) value from the remainder of the payload.
    let atoms_start = payload
        .find("ATOMS=(")
        .ok_or_else(|| "COLLECTION line is missing ATOMS=(...)".to_string())?;
    let after_paren = &payload[atoms_start + "ATOMS=(".len()..];
    let close = after_paren
        .find(')')
        .ok_or_else(|| "COLLECTION ATOMS value is not closed".to_string())?;
    let inner = &after_paren[..close];

    // First number is the count; the rest are 1-based V3000 atom indices.
    let values: Vec<&str> = inner.split_whitespace().collect();
    let Some((count, atom_ids)) = values.split_first() else {
        return Err("COLLECTION ATOMS value is empty".to_string());
    };
    let declared: usize = count
        .parse()
        .map_err(|_| "COLLECTION ATOMS count is not an integer".to_string())?;
    if atom_ids.len() != declared {
        return Err(format!(
            "COLLECTION ATOMS declares {declared} ids but contains {}",
            atom_ids.len()
        ));
    }
    // Collect and deduplicate atom indices.
    // A malformed V3000 file may list the same atom twice in one group.
    // Deduplication here mirrors the defensive check in StereoGroup::new()
    // and makes the intent explicit at the parser level (RDKit PR #9258).
    let mut seen = std::collections::HashSet::new();
    let atom_indices: Vec<AtomIdx> = atom_ids
        .iter()
        .map(|s| {
            let v3k: u32 = s
                .parse()
                .map_err(|_| "COLLECTION ATOMS contains a non-integer atom id".to_string())?;
            resolve_atom_idx(v3k, atom_idx_map)
                .ok_or_else(|| format!("COLLECTION ATOMS references unknown atom id {v3k}"))
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .filter(|idx| seen.insert(*idx))
        .collect();

    if atom_indices.is_empty() {
        return Err("COLLECTION ATOMS has no atom ids".to_string());
    }

    Ok(Some(StereoGroup::new(kind, atom_indices)))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "mol3000/tests.rs"]
mod tests;

mod writer;
pub use writer::*;
