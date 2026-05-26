//! MOL V3000 (Extended Ctab) parser.
//!
//! Reference: MDL/Dassault Systèmes CTfile Formats specification, V3000 section.
//!
//! Only parsing is implemented; writing V3000 is not required.

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

use crate::error::MolParseError;
use crate::mol2000::MolMetadata;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Prefix that marks every V3000 continuation line.
const V30_PREFIX: &str = "M  V30 ";

/// Return a `V3000ParseError` for the given 1-based line number and message.
#[inline]
fn v3k_err(line: usize, msg: impl Into<String>) -> MolParseError {
    MolParseError::V3000ParseError { line, msg: msg.into() }
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

            result.push(LogicalLine { line_num: first_line, payload: text });
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
        if let Some(rest) = tok.strip_prefix(key) {
            if let Some(val) = rest.strip_prefix('=') {
                return Some(val.to_string());
            }
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
pub fn parse_mol_v3000(input: &str) -> Result<(Molecule, MolMetadata), MolParseError> {
    // Collect all physical lines with 1-based numbering.
    let all_lines: Vec<(usize, &str)> = input
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l))
        .collect();

    if all_lines.len() < 4 {
        return Err(MolParseError::UnexpectedEnd);
    }

    // -- Header lines 1–3 ----------------------------------------------------

    let name = all_lines[0].1.to_string();
    // Line 2 (program/date) is read but discarded.
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

    enum State {
        BeforeCtab,
        InCtab,
        InAtomBlock,
        AfterAtomBlock,
        InBondBlock,
        AfterBondBlock,
        Done,
    }

    let mut state = State::BeforeCtab;
    let mut expected_atoms: usize = 0;
    let mut _expected_bonds: usize = 0;

    for LogicalLine { line_num, payload } in &v30_lines {
        let lnum = *line_num;
        let tokens: Vec<&str> = payload.split_whitespace().collect();
        if tokens.is_empty() {
            continue;
        }

        match state {
            State::BeforeCtab => {
                // Expect: BEGIN CTAB
                if tokens.len() >= 2 && tokens[0] == "BEGIN" && tokens[1] == "CTAB" {
                    state = State::InCtab;
                }
                // Ignore any lines before BEGIN CTAB.
            }

            State::InCtab => {
                // Expect: COUNTS <na> <nb> ... or BEGIN ATOM
                if tokens[0] == "COUNTS" {
                    if tokens.len() < 3 {
                        return Err(v3k_err(lnum, "COUNTS line has fewer than 2 values"));
                    }
                    expected_atoms = tokens[1].parse::<usize>().map_err(|_| {
                        v3k_err(lnum, format!("cannot parse atom count from '{}'", tokens[1]))
                    })?;
                    _expected_bonds = tokens[2].parse::<usize>().map_err(|_| {
                        v3k_err(lnum, format!("cannot parse bond count from '{}'", tokens[2]))
                    })?;
                } else if tokens.len() >= 2 && tokens[0] == "BEGIN" && tokens[1] == "ATOM" {
                    state = State::InAtomBlock;
                } else if tokens.len() >= 2 && tokens[0] == "END" && tokens[1] == "CTAB" {
                    state = State::Done;
                }
            }

            State::InAtomBlock => {
                if tokens.len() >= 2 && tokens[0] == "END" && tokens[1] == "ATOM" {
                    state = State::AfterAtomBlock;
                    // Validate atom count.
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
                        detail: format!("V3000 atom line needs at least 6 fields, got {}", tokens.len()),
                    });
                }

                let v3k_idx = tokens[0].parse::<u32>().map_err(|_| {
                    MolParseError::InvalidAtomLine {
                        line: lnum,
                        detail: format!("cannot parse atom index from '{}'", tokens[0]),
                    }
                })?;

                // Strip bracket notation e.g. "[OH]" → "OH", "[R]" → "R".
                let raw_sym = tokens[1];
                let sym = raw_sym.trim_start_matches('[').trim_end_matches(']');

                let element = Element::from_symbol(sym).ok_or_else(|| MolParseError::UnknownElement {
                    symbol: sym.to_string(),
                    line: lnum,
                })?;

                // Atom-map number (positional field 6, 0 = no mapping).
                let aamap_raw = tokens[5].parse::<u16>().unwrap_or(0);
                let atom_map = if aamap_raw == 0 { None } else { Some(aamap_raw) };

                // Parse optional key=value pairs from tokens[6..].
                let kv_tokens = if tokens.len() > 6 { &tokens[6..] } else { &[] as &[&str] };

                let charge: i8 = parse_kv(kv_tokens, "CHG")
                    .and_then(|v| v.parse::<i8>().ok())
                    .unwrap_or(0);

                let isotope: Option<u16> = parse_kv(kv_tokens, "MASS")
                    .and_then(|v| v.parse::<u16>().ok());

                // HCOUNT: -1 means unspecified; treat as None.
                let hydrogen_count: Option<u8> = parse_kv(kv_tokens, "HCOUNT")
                    .and_then(|v| {
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
            }

            State::AfterAtomBlock => {
                if tokens.len() >= 2 && tokens[0] == "BEGIN" && tokens[1] == "BOND" {
                    state = State::InBondBlock;
                } else if tokens.len() >= 2 && tokens[0] == "END" && tokens[1] == "CTAB" {
                    state = State::Done;
                }
            }

            State::InBondBlock => {
                if tokens.len() >= 2 && tokens[0] == "END" && tokens[1] == "BOND" {
                    state = State::AfterBondBlock;
                    continue;
                }

                // Bond line: <idx> <type> <atom1> <atom2> [KEY=VAL ...]
                if tokens.len() < 4 {
                    return Err(MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: format!("V3000 bond line needs at least 4 fields, got {}", tokens.len()),
                    });
                }

                let btype_raw = tokens[1].parse::<u8>().map_err(|_| {
                    MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: format!("cannot parse bond type from '{}'", tokens[1]),
                    }
                })?;

                let a1_v3k = tokens[2].parse::<u32>().map_err(|_| {
                    MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: format!("cannot parse atom1 index from '{}'", tokens[2]),
                    }
                })?;

                let a2_v3k = tokens[3].parse::<u32>().map_err(|_| {
                    MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: format!("cannot parse atom2 index from '{}'", tokens[3]),
                    }
                })?;

                // Resolve V3000 1-based atom indices → builder AtomIdx.
                let a1 = resolve_atom_idx(a1_v3k, &atom_idx_map).ok_or_else(|| {
                    MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: format!("atom index {} not found in atom block", a1_v3k),
                    }
                })?;

                let a2 = resolve_atom_idx(a2_v3k, &atom_idx_map).ok_or_else(|| {
                    MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: format!("atom index {} not found in atom block", a2_v3k),
                    }
                })?;

                let order = match btype_raw {
                    1 => BondOrder::Single,
                    2 => BondOrder::Double,
                    3 => BondOrder::Triple,
                    4 => BondOrder::Aromatic,
                    // Unknown/query bond types fall back to Single.
                    _ => BondOrder::Single,
                };

                builder.add_bond(a1, a2, order).map_err(|e| {
                    MolParseError::InvalidBondLine {
                        line: lnum,
                        detail: format!("{e}"),
                    }
                })?;
            }

            State::AfterBondBlock => {
                if tokens.len() >= 2 && tokens[0] == "END" && tokens[1] == "CTAB" {
                    state = State::Done;
                }
            }

            State::Done => {
                // Ignore everything after END CTAB.
            }
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

    Ok((builder.build(), metadata))
}

/// Look up the builder `AtomIdx` for a given V3000 1-based index.
fn resolve_atom_idx(v3k_idx: u32, map: &[(u32, AtomIdx)]) -> Option<AtomIdx> {
    map.iter().find(|&&(k, _)| k == v3k_idx).map(|&(_, v)| v)
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
        let (_, bond) = mol.bond_between(AtomIdx(0), AtomIdx(1)).expect("bond 0-1 exists");
        assert_eq!(bond.order, BondOrder::Single);
    }

    // -----------------------------------------------------------------------
    // Test 4: ethanol — bond between atoms 1 and 2 is Single
    // -----------------------------------------------------------------------
    #[test]
    fn test_ethanol_bond_1_2_single() {
        let (mol, _) = parse_mol_v3000(ETHANOL_V3K).expect("parse ethanol");
        let (_, bond) = mol.bond_between(AtomIdx(1), AtomIdx(2)).expect("bond 1-2 exists");
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
        let (_, bond) = mol.bond_between(AtomIdx(0), AtomIdx(1)).expect("bond exists");
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
        let (_, bond) = mol.bond_between(AtomIdx(0), AtomIdx(1)).expect("bond exists");
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
        let (_, bond) = mol.bond_between(AtomIdx(0), AtomIdx(1)).expect("bond exists");
        assert_eq!(bond.order, BondOrder::Triple);
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
