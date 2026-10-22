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
use chematic_perception::apply_ctab_local_parity_from_wedges_with_diagnostics;

use crate::error::MolParseError;
use crate::mol2000::{
    CoordinateDimension, GeometryRank, MolMetadata, MolReadReport, Stereo3DDiagnostic,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Prefix that marks every V3000 continuation line.
const V30_PREFIX: &str = "M  V30 ";

mod coordinates;
pub use coordinates::parse_mol_v3000_with_coords_3d;
mod source;
pub use source::{
    V3000SourceLimits, roundtrip_mol_v3000_source, roundtrip_mol_v3000_source_with_limits,
};

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
struct LogicalLine<'a> {
    /// 1-based line number of the first physical line.
    line_num: usize,
    /// Payload after stripping the `M  V30 ` prefix and joining continuations.
    payload: std::borrow::Cow<'a, str>,
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
fn collect_v30_lines<'a>(lines: &[(usize, &'a str)]) -> Vec<LogicalLine<'a>> {
    let mut result = Vec::new();

    let mut iter = lines.iter().peekable();
    while let Some(&(lineno, raw)) = iter.next() {
        if let Some(payload) = raw.strip_prefix(V30_PREFIX) {
            let mut text = std::borrow::Cow::Borrowed(payload);
            let first_line = lineno;

            // Handle line continuations: a payload ending with `-` means the
            // next `M  V30 ` line is a continuation.
            while text.ends_with('-') {
                // Drop the trailing hyphen.
                text.to_mut().pop();
                // Consume the next V30 line if available.
                match iter.next() {
                    Some(&(_, cont_raw)) => {
                        if let Some(cont_payload) = cont_raw.strip_prefix(V30_PREFIX) {
                            // Append with a separating space so that tokens
                            // that were split across lines re-join cleanly.
                            text.to_mut().push(' ');
                            text.to_mut().push_str(cont_payload);
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

/// Source ids are often exactly 1..=len in parser order. Keep that proof
/// while appending; sparse, permuted, zero and duplicate ids retain the
/// historical first-match lookup instead of guessing an index.
struct V3000IndexMap<T> {
    entries: Vec<(u32, T)>,
    dense: bool,
}

impl<T: Copy> V3000IndexMap<T> {
    fn new() -> Self {
        Self {
            entries: Vec::new(),
            dense: true,
        }
    }

    fn push(&mut self, entry: (u32, T)) {
        self.dense &= entry.0 as usize == self.entries.len() + 1;
        self.entries.push(entry);
    }

    fn resolve(&self, source: u32) -> Option<T> {
        if self.dense {
            source
                .checked_sub(1)
                .and_then(|index| self.entries.get(index as usize))
                .map(|&(_, value)| value)
        } else {
            self.entries
                .iter()
                .find(|&&(id, _)| id == source)
                .map(|&(_, value)| value)
        }
    }
}

impl<T> std::ops::Deref for V3000IndexMap<T> {
    type Target = [(u32, T)];

    fn deref(&self) -> &Self::Target {
        &self.entries
    }
}

// ---------------------------------------------------------------------------
// Key=value parsing helper
// ---------------------------------------------------------------------------

/// Parse optional `KEY=VALUE` pairs from the tail of a V30 atom or bond line.
///
/// Tokens are space-separated. Any token containing `=` is treated as a
/// key-value pair; tokens without `=` that appear after the required positional
/// fields are silently ignored.
fn parse_kv<'a>(tokens: &[&'a str], key: &str) -> Option<&'a str> {
    for tok in tokens {
        if let Some(rest) = tok.strip_prefix(key)
            && let Some(val) = rest.strip_prefix('=')
        {
            return Some(val);
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

fn parse_v3000_coords(tokens: &[&str], line: usize) -> Result<(f64, f64, f64), MolParseError> {
    let parse = |index: usize, axis: &str| {
        tokens[index]
            .parse::<f64>()
            .map_err(|_| MolParseError::InvalidAtomLine {
                line,
                detail: format!("cannot parse {axis} coordinate from '{}'", tokens[index]),
            })
            .and_then(|value| {
                value
                    .is_finite()
                    .then_some(value)
                    .ok_or_else(|| MolParseError::InvalidAtomLine {
                        line,
                        detail: format!(
                            "{axis} coordinate is not finite (NaN/Infinite): '{}'",
                            tokens[index]
                        ),
                    })
            })
    };
    Ok((
        parse(2, "x")?,
        parse(3, "y")?,
        // Common finite positive zeros need no general float parsing.
        // Keep every other spelling, including negative zero, on that path.
        match tokens[4] {
            "0" | "0.0000" | "0.000000" => 0.0,
            _ => parse(4, "z")?,
        },
    ))
}

fn v3000_atom_identity(
    symbol: &str,
    rgroup_number: Option<u16>,
    line: usize,
) -> Result<(Atom, Option<RGroupLabel>), MolParseError> {
    if symbol == "*" {
        return Ok((Atom::wildcard(), None));
    }
    if symbol == "R" {
        return Ok((Atom::wildcard(), Some(RGroupLabel::unnumbered())));
    }
    if symbol == "R#" {
        return Ok((
            Atom::wildcard(),
            Some(
                rgroup_number.map_or_else(RGroupLabel::unnumbered, |number| {
                    RGroupLabel::numbered(number).expect("validated non-zero number")
                }),
            ),
        ));
    }
    if let Some(number) = symbol.strip_prefix('R')
        && !number.is_empty()
        && number.bytes().all(|byte| byte.is_ascii_digit())
    {
        let label = number
            .parse::<u16>()
            .ok()
            .and_then(RGroupLabel::numbered)
            .ok_or_else(|| MolParseError::UnknownElement {
                symbol: symbol.to_string(),
                line,
            })?;
        return Ok((Atom::wildcard(), Some(label)));
    }
    if let Some(number) = rgroup_number {
        return Ok((
            Atom::wildcard(),
            Some(RGroupLabel::numbered(number).expect("validated non-zero number")),
        ));
    }
    let element_symbol = if symbol == "D" { "H" } else { symbol };
    let element =
        Element::from_symbol(element_symbol).ok_or_else(|| MolParseError::UnknownElement {
            symbol: symbol.to_string(),
            line,
        })?;
    Ok((Atom::new(element), None))
}

/// Parse one V3000 atom record and append it to the parser's parallel state.
fn parse_v3000_atom_line(
    tokens: &[&str],
    line: usize,
    builder: &mut MoleculeBuilder,
    atom_idx_map: &mut V3000IndexMap<AtomIdx>,
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
    let deuterium = sym == "D";
    let (x, y, z) = parse_v3000_coords(tokens, line)?;

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

    let (mut atom, r_group) = v3000_atom_identity(sym, rgroup_number, line)?;
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
    atom_idx_map: &V3000IndexMap<AtomIdx>,
    bond_idx_map: &mut V3000IndexMap<BondIdx>,
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

    let v3k_idx = tokens[0]
        .parse::<u32>()
        .map_err(|_| MolParseError::InvalidBondLine {
            line,
            detail: format!("cannot parse bond index from '{}'", tokens[0]),
        })?;
    if bond_idx_map.resolve(v3k_idx).is_some() {
        return Err(MolParseError::InvalidBondLine {
            line,
            detail: format!("duplicate V3000 bond index {v3k_idx}"),
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
    let a1 = atom_idx_map
        .resolve(a1_v3k)
        .ok_or_else(|| MolParseError::InvalidBondLine {
            line,
            detail: format!("atom index {a1_v3k} not found in atom block"),
        })?;
    let a2 = atom_idx_map
        .resolve(a2_v3k)
        .ok_or_else(|| MolParseError::InvalidBondLine {
            line,
            detail: format!("atom index {a2_v3k} not found in atom block"),
        })?;

    let kv_tokens = tokens.get(4..).unwrap_or(&[]);
    let order = match btype_raw {
        0 => BondOrder::Zero,
        1 => match parse_kv(kv_tokens, "CFG") {
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
    bond_idx_map.push((v3k_idx, bidx));
    if kv_tokens.is_empty() {
        return Ok(());
    }
    if btype_raw == 2 && parse_kv(kv_tokens, "CFG") == Some("2") {
        explicitly_unspecified_ez.insert(bidx);
    }
    // Preserve richer V3000 bond attributes without pretending to interpret
    // them. Wedge CFG is regenerated from BondOrder; double-bond CFG=2
    // has no core tag, so retain it explicitly for lossless rewriting. ENDPTS is the one supported reference-bearing attribute: its
    // atom ids must follow the normalized output order when source ids are
    // non-contiguous. ATTACH and future keys remain opaque.
    let opaque_properties = kv_tokens
        .iter()
        .filter(|token| !token.starts_with("CFG=") || (btype_raw == 2 && token.trim() == "CFG=2"))
        .copied()
        .collect::<Vec<_>>()
        .join(" ");
    let opaque_properties = normalize_v3000_endpoints(opaque_properties, atom_idx_map);
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
fn normalize_v3000_endpoints(properties: String, atom_idx_map: &[(u32, AtomIdx)]) -> String {
    let Some(start) = properties.find("ENDPTS=(") else {
        return properties;
    };
    let value_start = start + "ENDPTS=(".len();
    let Some(close_offset) = properties[value_start..].find(')') else {
        return properties;
    };
    let value_end = value_start + close_offset;
    let values: Vec<&str> = properties[value_start..value_end]
        .split_whitespace()
        .collect();
    let Some((declared, atom_ids)) = values.split_first() else {
        return properties;
    };
    let Ok(declared) = declared.parse::<usize>() else {
        return properties;
    };
    if atom_ids.len() != declared {
        return properties;
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
        return properties;
    };
    let normalized = format!("ENDPTS=({declared} {})", writer_ids.join(" "));
    format!(
        "{}{}{}",
        &properties[..start],
        normalized,
        &properties[value_end + 1..]
    )
}

#[derive(Clone, Copy)]
enum V3000State {
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

struct V3000Parser {
    state: V3000State,
    builder: MoleculeBuilder,
    atom_idx_map: V3000IndexMap<AtomIdx>,
    bond_idx_map: V3000IndexMap<BondIdx>,
    coords: Vec<(f64, f64)>,
    raw_z: Vec<f64>,
    allocation_limit: usize,
    expected_atoms: usize,
    expected_bonds: usize,
    parsed_bonds: usize,
    expected_sgroups: Option<usize>,
    stereo_groups: Vec<ParsedStereoGroup>,
    sgroups: Vec<String>,
    typed_sgroup_ids: std::collections::HashSet<u32>,
    bond_properties: Vec<(u32, String)>,
    atom_properties: Vec<(u32, String)>,
    unspecified_ez: std::collections::HashSet<BondIdx>,
}

impl V3000Parser {
    fn new() -> Self {
        Self {
            state: V3000State::BeforeCtab,
            builder: MoleculeBuilder::new(),
            atom_idx_map: V3000IndexMap::new(),
            bond_idx_map: V3000IndexMap::new(),
            coords: Vec::new(),
            raw_z: Vec::new(),
            allocation_limit: 0,
            expected_atoms: 0,
            expected_bonds: 0,
            parsed_bonds: 0,
            expected_sgroups: None,
            stereo_groups: Vec::new(),
            sgroups: Vec::new(),
            typed_sgroup_ids: std::collections::HashSet::new(),
            bond_properties: Vec::new(),
            atom_properties: Vec::new(),
            unspecified_ez: std::collections::HashSet::new(),
        }
    }

    fn consume(&mut self, line: &LogicalLine, tokens: &[&str]) -> Result<(), MolParseError> {
        if tokens.is_empty() {
            return Ok(());
        }
        match self.state {
            V3000State::BeforeCtab => self.before_ctab(tokens),
            V3000State::InCtab => self.in_ctab(tokens, line.line_num),
            V3000State::InAtomBlock => self.in_atom_block(tokens, line.line_num),
            V3000State::AfterAtomBlock => self.after_atom_block(tokens, line.line_num),
            V3000State::InBondBlock => self.in_bond_block(tokens, line.line_num),
            V3000State::AfterBondBlock => self.after_bond_block(tokens),
            V3000State::InCollection => self.in_collection(tokens, line),
            V3000State::InSgroup => self.in_sgroup(tokens, line),
            V3000State::Done => Ok(()),
        }
    }

    fn before_ctab(&mut self, tokens: &[&str]) -> Result<(), MolParseError> {
        if is_marker(tokens, "BEGIN", "CTAB") {
            self.state = V3000State::InCtab;
        }
        Ok(())
    }

    fn in_ctab(&mut self, tokens: &[&str], line: usize) -> Result<(), MolParseError> {
        if tokens[0] == "COUNTS" {
            if tokens.len() < 4 {
                return Err(v3k_err(line, "COUNTS line has fewer than 3 values"));
            }
            self.expected_atoms = parse_v3000_count(tokens[1], line, "atom")?;
            self.expected_bonds = parse_v3000_count(tokens[2], line, "bond")?;
            self.expected_sgroups = Some(parse_v3000_count(tokens[3], line, "SGROUP")?);
        } else if is_marker(tokens, "BEGIN", "ATOM") {
            // Allocate once when entering ATOM, after the final COUNTS row.
            // A declared count is only an allocation hint: bound it by the
            // actual logical row count before trusting potentially huge input.
            let atoms = v3000_storage_capacity(self.expected_atoms, self.allocation_limit);
            let bonds = v3000_storage_capacity(self.expected_bonds, self.allocation_limit);
            self.builder = MoleculeBuilder::with_capacity(atoms, bonds);
            self.atom_idx_map.entries.reserve_exact(atoms);
            self.bond_idx_map.entries.reserve_exact(bonds);
            self.coords.reserve_exact(atoms);
            self.raw_z.reserve_exact(atoms);
            self.state = V3000State::InAtomBlock;
        } else if is_marker(tokens, "END", "CTAB") {
            self.state = V3000State::Done;
        }
        Ok(())
    }

    fn in_atom_block(&mut self, tokens: &[&str], line: usize) -> Result<(), MolParseError> {
        if is_marker(tokens, "END", "ATOM") {
            if self.builder.atom_count() != self.expected_atoms {
                return Err(v3k_err(
                    line,
                    format!(
                        "expected {} atoms, found {}",
                        self.expected_atoms,
                        self.builder.atom_count()
                    ),
                ));
            }
            self.state = V3000State::AfterAtomBlock;
            return Ok(());
        }
        parse_v3000_atom_line(
            tokens,
            line,
            &mut self.builder,
            &mut self.atom_idx_map,
            &mut self.coords,
            &mut self.raw_z,
            &mut self.atom_properties,
        )
    }

    fn after_atom_block(&mut self, tokens: &[&str], line: usize) -> Result<(), MolParseError> {
        if is_marker(tokens, "BEGIN", "BOND") {
            self.state = V3000State::InBondBlock;
        } else if is_marker(tokens, "END", "CTAB") {
            if self.expected_bonds != 0 {
                return Err(v3k_err(
                    line,
                    format!(
                        "COUNTS declares {} bonds but the BOND block is absent",
                        self.expected_bonds
                    ),
                ));
            }
            self.state = V3000State::Done;
        } else {
            return Err(v3k_err(
                line,
                "expected BEGIN BOND or END CTAB after END ATOM",
            ));
        }
        Ok(())
    }

    fn in_bond_block(&mut self, tokens: &[&str], line: usize) -> Result<(), MolParseError> {
        if is_marker(tokens, "END", "BOND") {
            if self.parsed_bonds != self.expected_bonds {
                return Err(v3k_err(
                    line,
                    format!(
                        "COUNTS declares {} bonds but found {}",
                        self.expected_bonds, self.parsed_bonds
                    ),
                ));
            }
            self.state = V3000State::AfterBondBlock;
            return Ok(());
        }
        parse_v3000_bond_line(
            tokens,
            line,
            &mut self.builder,
            &self.atom_idx_map,
            &mut self.bond_idx_map,
            &mut self.unspecified_ez,
            &mut self.bond_properties,
        )?;
        self.parsed_bonds += 1;
        Ok(())
    }

    fn after_bond_block(&mut self, tokens: &[&str]) -> Result<(), MolParseError> {
        if is_marker(tokens, "BEGIN", "COLLECTION") {
            self.state = V3000State::InCollection;
        } else if is_marker(tokens, "BEGIN", "SGROUP") {
            self.state = V3000State::InSgroup;
        } else if is_marker(tokens, "END", "CTAB") {
            self.state = V3000State::Done;
        }
        Ok(())
    }

    fn in_collection(&mut self, tokens: &[&str], line: &LogicalLine) -> Result<(), MolParseError> {
        if is_marker(tokens, "END", "COLLECTION") {
            self.state = V3000State::AfterBondBlock;
        } else if let Some(group) =
            parse_stereo_group_line(&line.payload, &self.atom_idx_map, &self.bond_idx_map)
                .map_err(|detail| v3k_err(line.line_num, detail))?
        {
            if let Some(existing) = self
                .stereo_groups
                .iter_mut()
                .find(|existing| existing.group.kind == group.group.kind)
            {
                existing.group.atom_indices.extend(group.group.atom_indices);
                existing
                    .group
                    .atom_indices
                    .sort_unstable_by_key(|idx| idx.0);
                existing.group.atom_indices.dedup();
                existing.bond_indices.extend(group.bond_indices);
                existing.bond_indices.sort_unstable_by_key(|idx| idx.0);
                existing.bond_indices.dedup();
            } else {
                self.stereo_groups.push(group);
            }
        }
        Ok(())
    }

    fn in_sgroup(&mut self, tokens: &[&str], line: &LogicalLine) -> Result<(), MolParseError> {
        if is_marker(tokens, "END", "SGROUP") {
            self.state = V3000State::AfterBondBlock;
        } else {
            validate_known_sgroup_references(
                &line.payload,
                &self.atom_idx_map,
                &mut self.typed_sgroup_ids,
            )
            .map_err(|detail| v3k_err(line.line_num, detail))?;
            self.sgroups.push(line.payload.to_string());
        }
        Ok(())
    }
}

// Preserve normal Vec growth capacities where the real input budget allows.
// Storage shape after parsing can affect subsequent prepared operations.
fn v3000_storage_capacity(count: usize, limit: usize) -> usize {
    let bounded = count.min(limit);
    if bounded == 0 {
        0
    } else {
        bounded
            .checked_next_power_of_two()
            .unwrap_or(limit)
            .max(4)
            .min(limit)
    }
}

fn parse_v3000_count(raw: &str, line: usize, kind: &str) -> Result<usize, MolParseError> {
    raw.parse::<usize>()
        .map_err(|_| v3k_err(line, format!("cannot parse {kind} count from '{raw}'")))
}

fn validate_v3000_physical_lines(input: &str) -> Result<Vec<(usize, &str)>, MolParseError> {
    let lines = input
        .lines()
        .enumerate()
        .map(|(index, line)| (index + 1, line))
        .collect::<Vec<_>>();
    if lines.len() < 4 {
        return Err(MolParseError::UnexpectedEnd);
    }
    let terminal = lines
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
    for &(line_number, line) in lines.iter().skip(4) {
        if line.trim().is_empty()
            || line.starts_with(V30_PREFIX)
            || (line == "M  END" && line_number == terminal.0)
        {
            continue;
        }
        return Err(v3k_err(
            line_number,
            "expected an 'M  V30 ' record or the terminal 'M  END' line",
        ));
    }
    Ok(lines)
}

fn v3000_metadata<'a>(
    lines: &'a [(usize, &'a str)],
) -> Result<(MolMetadata, &'a str), MolParseError> {
    let (counts_line, counts) = lines[3];
    if !counts.contains("V3000") {
        return Err(MolParseError::InvalidCountLine {
            line: counts_line,
            detail: "missing V3000 version tag".to_string(),
        });
    }
    Ok((
        MolMetadata {
            name: lines[0].1.to_string(),
            comment: lines[2].1.to_string(),
            v3000_sgroups: Vec::new(),
            v3000_bond_properties: Vec::new(),
            v3000_atom_properties: Vec::new(),
        },
        lines[1].1,
    ))
}

fn validate_v3000_terminal(parser: &V3000Parser) -> Result<(), MolParseError> {
    match parser.state {
        V3000State::Done | V3000State::AfterBondBlock => {}
        V3000State::InAtomBlock => return Err(v3k_err(0, "missing M  V30 END ATOM")),
        V3000State::InBondBlock => return Err(v3k_err(0, "missing M  V30 END BOND")),
        _ => return Err(MolParseError::UnexpectedEnd),
    }
    if let Some(expected) = parser.expected_sgroups
        && parser.sgroups.len() != expected
    {
        return Err(v3k_err(
            0,
            format!(
                "COUNTS declares {expected} SGROUP records but found {}",
                parser.sgroups.len()
            ),
        ));
    }
    Ok(())
}

fn v3000_geometry(
    mol: &mut Molecule,
    coords: &[(f64, f64)],
    raw_z: &[f64],
    dimension: CoordinateDimension,
) -> (
    Option<Coords3D>,
    GeometryRank,
    Vec<Stereo3DDiagnostic>,
    Vec<crate::mol2000::SquarePlanarPerceptionDiagnostic>,
) {
    let points = coords
        .iter()
        .zip(raw_z)
        .map(|(&(x, y), &z)| Point3::new(x, y, z))
        .collect::<Vec<_>>();
    let rank = crate::mol2000::classify_geometry_rank(&points);
    let conformer = matches!(rank, GeometryRank::Coplanar | GeometryRank::ThreeD)
        .then_some(Coords3D { points });
    let mut diagnostics = match (dimension, rank) {
        (CoordinateDimension::TwoD, GeometryRank::Coplanar | GeometryRank::ThreeD) => {
            vec![Stereo3DDiagnostic::DeclaredTwoDButNonzeroZ { observed: rank }]
        }
        (CoordinateDimension::ThreeD, GeometryRank::FlatZero | GeometryRank::Coplanar) => {
            vec![Stereo3DDiagnostic::DeclaredThreeDButFlat { observed: rank }]
        }
        _ => Vec::new(),
    };
    let square_planar = conformer.as_ref().map_or_else(Vec::new, |coords| {
        let perceived = crate::mol2000::perceive_square_planar_from_3d(mol, coords);
        let conflicts = crate::mol2000::wedge_vs_3d_conflicts(mol, coords);
        if rank == GeometryRank::ThreeD && dimension != CoordinateDimension::TwoD {
            crate::mol2000::correct_conflicting_wedges_from_3d(mol, &conflicts);
            crate::mol2000::perceive_tetrahedral_from_3d(mol, coords);
            crate::mol2000::correct_ring_alkene_directions_from_3d(mol, coords);
        }
        diagnostics.extend(conflicts);
        perceived
    });
    (conformer, rank, diagnostics, square_planar)
}

/// Bond order spent at a V3000 atom's total-valence field. A coordinate
/// bond spends valence at its acceptor, but not its donor.
fn v3000_valence_bond_order_sum(mol: &Molecule, atom: AtomIdx) -> u8 {
    mol.neighbors(atom)
        .map(|(_, bond_idx)| {
            let bond = mol.bond(bond_idx);
            if bond.order == BondOrder::Dative && bond.atom1 == atom {
                0
            } else {
                bond.order.order_int()
            }
        })
        .fold(0u8, u8::saturating_add)
}

/// Apply supported CTfile total-valence overrides after bonds are known.
/// VAL=0 leaves implicit hydrogen inference alone; -1 and the legacy zero
/// sentinel 15 suppress hydrogens. Retain the source attribute for rewriting.
fn apply_v3000_valence_overrides(mol: &mut Molecule, properties: &[(u32, String)]) {
    for (atom_id, properties) in properties {
        let value = properties
            .split_whitespace()
            .find_map(|token| token.strip_prefix("VAL="))
            .and_then(|value| value.parse::<i32>().ok());
        let Some(value) = value.filter(|&value| value != 0) else {
            continue;
        };
        let valence = if value == -1 || value == 15 {
            0
        } else {
            let Ok(value) = u8::try_from(value) else {
                continue;
            };
            value
        };
        let atom = AtomIdx(atom_id - 1);
        // Keep the established V2000 boundary for aromatic/query valences.
        if mol
            .neighbors(atom)
            .any(|(_, bond)| mol.bond(bond).order == BondOrder::Aromatic)
        {
            continue;
        }
        let drawn = v3000_valence_bond_order_sum(mol, atom);
        mol.set_hydrogen_count(atom, Some(valence.saturating_sub(drawn)));
    }
}

/// The core molecule stores explicit hydrogen count but not radical spin.
/// For supported neutral main-group atoms, RAD constrains the hydrogens that
/// would otherwise be added by ordinary valence inference. Keep the original
/// RAD attribute in metadata for V3000 reserialization.
fn apply_v3000_radical_hydrogens(mol: &mut Molecule, properties: &[(u32, String)]) {
    for (atom_id, attributes) in properties {
        let tokens: Vec<_> = attributes.split_whitespace().collect();
        let radical_deficit = match parse_kv(&tokens, "RAD") {
            Some("2") => 1u8, // doublet
            Some("3") => 2u8, // triplet
            _ => continue,
        };
        if parse_kv(&tokens, "VAL").is_some() {
            continue;
        }
        let atom = AtomIdx(atom_id - 1);
        let identity = mol.atom(atom);
        if identity.charge != 0
            || identity.hydrogen_count.is_some()
            || !matches!(identity.element, Element::C | Element::N | Element::O)
        {
            continue;
        }
        let Some(&normal_valence) = identity.element.normal_valences().first() else {
            continue;
        };
        let aromatic_bonds = mol
            .neighbors(atom)
            .filter(|&(_, bond)| mol.bond(bond).order == BondOrder::Aromatic)
            .count() as u8;
        let drawn = v3000_valence_bond_order_sum(mol, atom).saturating_add(aromatic_bonds / 2);
        mol.set_hydrogen_count(
            atom,
            Some(normal_valence.saturating_sub(drawn.saturating_add(radical_deficit))),
        );
    }
}

fn finish_v3000_report(
    parser: V3000Parser,
    mut metadata: MolMetadata,
    line2: &str,
) -> Result<(MolReadReport, Vec<f64>), MolParseError> {
    validate_v3000_terminal(&parser)?;
    let mut mol = parser.builder.build();
    crate::mol2000::flag_aromatic_bond_atoms(&mut mol);
    apply_v3000_valence_overrides(&mut mol, &parser.atom_properties);
    apply_v3000_radical_hydrogens(&mut mol, &parser.atom_properties);
    metadata.v3000_sgroups = parser.sgroups;
    metadata.v3000_bond_properties = parser.bond_properties;
    metadata.v3000_atom_properties = parser.atom_properties;
    if !parser.stereo_groups.is_empty() {
        mol.set_stereo_groups_with_bonds(
            parser
                .stereo_groups
                .into_iter()
                .map(|group| (group.group, group.bond_indices))
                .collect(),
        );
    }
    let stereo_diagnostics =
        apply_ctab_local_parity_from_wedges_with_diagnostics(&mut mol, &parser.coords);
    let dimension = crate::mol2000::parse_dimension_code(line2);
    let ez_diagnostics = crate::mol2000::apply_ctab_ez_directions(
        &mut mol,
        &parser.coords,
        &parser.raw_z,
        dimension,
        &parser.unspecified_ez,
    );
    let (conformer, geometry_rank, stereo3d_diagnostics, square_planar_diagnostics) =
        v3000_geometry(&mut mol, &parser.coords, &parser.raw_z, dimension);
    let report = MolReadReport {
        mol,
        metadata,
        coords: parser.coords,
        stereo_diagnostics,
        ez_diagnostics,
        conformer,
        coordinate_dimension: dimension,
        geometry_rank,
        stereo3d_diagnostics,
        square_planar_diagnostics,
    };
    Ok((report, parser.raw_z))
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
/// This field is the 2D projection; the optional conformer is selected for
/// perception. Use [`parse_mol_v3000_with_coords_3d`] to retain every source z
/// value independently of dimension headers and perception tolerances.
/// This is the one parsing core for V3000 MOL text --
/// [`parse_mol_v3000_with_coords`]/[`parse_mol_v3000`] are thin wrappers that
/// discard `stereo_diagnostics`. Only bond-line `CFG` (wedge direction) is
/// decoded; atom-line `CFG` (parity) is out of scope, matching RDKit's own
/// primary `MolFromMolBlock` path (bond-CFG/wedge-based).
pub fn read_mol_v3000_with_diagnostics(input: &str) -> Result<MolReadReport, MolParseError> {
    read_v3000_with_source_z(input).map(|(report, _)| report)
}

fn read_v3000_with_source_z(input: &str) -> Result<(MolReadReport, Vec<f64>), MolParseError> {
    let lines = validate_v3000_physical_lines(input)?;
    let (metadata, line2) = v3000_metadata(&lines)?;
    let logical_lines = collect_v30_lines(&lines);
    let mut parser = V3000Parser::new();
    parser.allocation_limit = logical_lines.len();
    // Logical payloads live until the entire parse ends, including owned
    // continuation text. Reuse only the token vector's storage; clear every
    // row before borrowing that row's exact Unicode-whitespace tokens.
    let mut tokens = Vec::new();
    for line in &logical_lines {
        tokens.clear();
        tokens.extend(line.payload.split_whitespace());
        parser.consume(line, &tokens)?;
    }
    finish_v3000_report(parser, metadata, line2)
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
/// - `MDLV30/STEBABS BONDS=(1 2)`     → absolute bond stereo
/// - `MDLV30/STEBREL1 BONDS=(1 4)`    → `StereoGroupKind::Or(1)` bond stereo
/// - `MDLV30/STEBRAC2 BONDS=(1 6)`    → `StereoGroupKind::And(2)` bond stereo
///
/// `STEOR` and `STEAND` are also accepted as legacy input spellings, but the
/// writer uses the MDL tokens emitted by current RDKit.
struct ParsedStereoGroup {
    group: StereoGroup,
    bond_indices: Vec<BondIdx>,
}

fn parse_stereo_group_kind(token: &str) -> Result<Option<(StereoGroupKind, &'static str)>, String> {
    let parsed = if token == "MDLV30/STEABS" {
        (StereoGroupKind::Absolute, "ATOMS")
    } else if token == "MDLV30/STEBABS" {
        (StereoGroupKind::Absolute, "BONDS")
    } else if let Some(id) = token
        .strip_prefix("MDLV30/STEREL")
        .or_else(|| token.strip_prefix("MDLV30/STEOR"))
    {
        (
            StereoGroupKind::Or(
                id.parse()
                    .map_err(|_| "COLLECTION OR group id is not an integer".to_string())?,
            ),
            "ATOMS",
        )
    } else if let Some(id) = token.strip_prefix("MDLV30/STEBREL") {
        (
            StereoGroupKind::Or(
                id.parse()
                    .map_err(|_| "COLLECTION OR group id is not an integer".to_string())?,
            ),
            "BONDS",
        )
    } else if let Some(id) = token
        .strip_prefix("MDLV30/STERAC")
        .or_else(|| token.strip_prefix("MDLV30/STEAND"))
    {
        (
            StereoGroupKind::And(
                id.parse()
                    .map_err(|_| "COLLECTION AND group id is not an integer".to_string())?,
            ),
            "ATOMS",
        )
    } else if let Some(id) = token.strip_prefix("MDLV30/STEBRAC") {
        (
            StereoGroupKind::And(
                id.parse()
                    .map_err(|_| "COLLECTION AND group id is not an integer".to_string())?,
            ),
            "BONDS",
        )
    } else {
        return Ok(None);
    };
    Ok(Some(parsed))
}

fn parse_collection_member_ids(payload: &str, member_label: &str) -> Result<Vec<u32>, String> {
    let marker = format!("{member_label}=(");
    let members_start = payload
        .find(&marker)
        .ok_or_else(|| format!("COLLECTION line is missing {member_label}=(...)"))?;
    let after_paren = &payload[members_start + marker.len()..];
    let close = after_paren
        .find(')')
        .ok_or_else(|| format!("COLLECTION {member_label} value is not closed"))?;
    let values = after_paren[..close].split_whitespace().collect::<Vec<_>>();
    let Some((count, member_ids)) = values.split_first() else {
        return Err(format!("COLLECTION {member_label} value is empty"));
    };
    let declared = count
        .parse::<usize>()
        .map_err(|_| format!("COLLECTION {member_label} count is not an integer"))?;
    if member_ids.len() != declared {
        return Err(format!(
            "COLLECTION {member_label} declares {declared} ids but contains {}",
            member_ids.len()
        ));
    }
    let ids = member_ids
        .iter()
        .map(|id| {
            id.parse::<u32>()
                .map_err(|_| format!("COLLECTION {member_label} contains a non-integer member id"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    if ids.is_empty() {
        let member_kind = if member_label == "ATOMS" {
            "atom"
        } else {
            "bond"
        };
        return Err(format!(
            "COLLECTION {member_label} has no {member_kind} ids"
        ));
    }
    Ok(ids)
}

fn parse_stereo_group_line(
    payload: &str,
    atom_idx_map: &[(u32, AtomIdx)],
    bond_idx_map: &[(u32, BondIdx)],
) -> Result<Option<ParsedStereoGroup>, String> {
    // First token is the group kind key.
    let Some(first_tok) = payload.split_whitespace().next() else {
        return Ok(None);
    };

    let Some((kind, member_label)) = parse_stereo_group_kind(first_tok)? else {
        return Ok(None);
    };
    let ids = parse_collection_member_ids(payload, member_label)?;

    let (atom_indices, bond_indices) = if member_label == "ATOMS" {
        let atoms = ids
            .into_iter()
            .map(|v3k| {
                resolve_atom_idx(v3k, atom_idx_map)
                    .ok_or_else(|| format!("COLLECTION ATOMS references unknown atom id {v3k}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        (atoms, Vec::new())
    } else {
        let bonds = ids
            .into_iter()
            .map(|v3k| {
                bond_idx_map
                    .iter()
                    .find(|&&(source, _)| source == v3k)
                    .map(|&(_, idx)| idx)
                    .ok_or_else(|| format!("COLLECTION BONDS references unknown bond id {v3k}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        (Vec::new(), bonds)
    };

    let mut seen = std::collections::HashSet::new();
    let bond_indices = bond_indices
        .into_iter()
        .filter(|idx| seen.insert(*idx))
        .collect();
    Ok(Some(ParsedStereoGroup {
        group: StereoGroup::new(kind, atom_indices),
        bond_indices,
    }))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "mol3000/tests.rs"]
mod tests;

mod writer;
pub use writer::*;

#[cfg(test)]
#[path = "mol3000/token_buffer_tests.rs"]
mod token_buffer_tests;

#[cfg(test)]
#[path = "mol3000/zero_z_tests.rs"]
mod zero_z_tests;

#[cfg(test)]
#[path = "mol3000/count_capacity_tests.rs"]
mod count_capacity_tests;

#[cfg(test)]
#[path = "mol3000/attribute_tail_tests.rs"]
mod attribute_tail_tests;
