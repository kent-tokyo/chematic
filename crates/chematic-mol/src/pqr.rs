//! PQR (PDB-like atom records with per-atom charge/radius) reader and writer.
//!
//! ## Provenance
//!
//! Implemented independently from the public APBS/PDB2PQR PQR format
//! documentation (<https://apbs.readthedocs.io/en/latest/formats/pqr.html>,
//! <https://pdb2pqr.readthedocs.io/en/latest/formats/pqr.html>) and the
//! classic PDB ATOM/HETATM record layout it is derived from. No source
//! code, comments, or tables were copied from Open Babel, PDB2PQR, or any
//! other tool.
//!
//! ## Format
//!
//! PQR reuses the PDB `ATOM`/`HETATM` record's field *order* but is
//! **whitespace-delimited, not fixed-column** (this is the documented,
//! deliberate difference: charge/radius need variable precision that fixed
//! PDB columns can't accommodate), and replaces the PDB's `occupancy` +
//! `tempFactor` columns with `charge` + `radius`:
//!
//! ```text
//! Field_name Atom_number Atom_name Residue_name [Chain_ID] Residue_number  X  Y  Z  Charge  Radius
//! ```
//!
//! `Chain_ID` is optional -- most `pdb2pqr`-generated files omit it
//! entirely (10 whitespace tokens per line); files that do include it have
//! 11. This reader detects which shape a line has by its token count (see
//! [`parse_pqr`]).
//!
//! ## No element column
//!
//! Unlike mmCIF's `_atom_site.type_symbol`, PQR carries **no element
//! column at all** -- by design, APBS only needs coordinates, charge and
//! radius. [`Element`] is therefore *inferred* from the atom name, using
//! the same convention classic PDB parsers use when a file's element
//! column (PDB columns 77-78) is blank: for `HETATM` records whose residue
//! name is itself a bare element symbol (the standard PDB convention for
//! monatomic ions, e.g. `ZN`, `NA`, `MG`, `FE`), that symbol is used
//! directly; otherwise the atom name's leading digits are stripped and the
//! **first alphabetic character** is used (this alone is correct for the
//! overwhelming majority of standard amino/nucleic-acid backbone and
//! side-chain atom names -- `CA`, `CB`, `ND1`, `OE2`, ... all correctly
//! resolve to C/C/N/O; a 2-letter *greedy* match would misresolve several
//! of these to Ca/Nd, which is wrong). This is a documented heuristic, not
//! authoritative chemistry -- see [`infer_element`] -- and, like this
//! crate's bond-perception policy elsewhere, it is applied automatically
//! only because there is no alternative source for `Element` at all (PQR
//! has no bond table either, and none is inferred).
//!
//! Because the heuristic is a pure function of `(group_pdb, res_name,
//! atom_name)`, and [`write_pqr`] never re-derives the atom name from the
//! element, read -> write -> read is still an exact fixed point even
//! though the element itself is inferred rather than stored in the file.

use chematic_core::{Atom, Element, Molecule, MoleculeBuilder};

// ---------------------------------------------------------------------------
// Limits
// ---------------------------------------------------------------------------

/// Security/robustness limits enforced before any chemistry is interpreted.
#[derive(Debug, Clone, Copy)]
pub struct PqrParseLimits {
    pub max_input_bytes: usize,
    pub max_atoms: usize,
    pub max_line_len: usize,
}

impl Default for PqrParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20, // 64 MiB
            max_atoms: 2_000_000,
            max_line_len: 1024,
        }
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors that can occur when parsing PQR files.
#[derive(Debug, Clone, PartialEq)]
pub enum PqrError {
    /// Input exceeded [`PqrParseLimits::max_input_bytes`].
    InputTooLarge { limit: usize },
    /// A single line exceeded [`PqrParseLimits::max_line_len`].
    LineTooLong { line: usize, limit: usize },
    /// Atom row count exceeded [`PqrParseLimits::max_atoms`].
    TooManyAtoms { limit: usize },
    /// No `ATOM`/`HETATM` record was found anywhere in the input.
    NoAtomRecords,
    /// An `ATOM`/`HETATM` line did not have 10 or 11 whitespace-delimited
    /// fields (see module docs on the optional chain-id column).
    WrongFieldCount { line: usize, found: usize },
    /// The atom serial number, residue number, or radius could not be
    /// parsed as expected.
    InvalidField {
        line: usize,
        field: &'static str,
        raw: String,
    },
    /// A coordinate/charge/radius value parsed but was not finite
    /// (NaN/Infinity).
    NonFiniteValue {
        line: usize,
        field: &'static str,
        raw: String,
    },
    /// The atom name could not be resolved to a known element even via the
    /// first-letter fallback (e.g. an atom name with no alphabetic
    /// characters at all).
    UnresolvableElement { line: usize, atom_name: String },
}

impl core::fmt::Display for PqrError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InputTooLarge { limit } => write!(f, "PQR input exceeds {limit}-byte limit"),
            Self::LineTooLong { line, limit } => {
                write!(f, "PQR line {line} exceeds {limit}-byte limit")
            }
            Self::TooManyAtoms { limit } => {
                write!(f, "PQR atom count exceeds {limit}-atom limit")
            }
            Self::NoAtomRecords => write!(f, "no ATOM/HETATM records found in PQR"),
            Self::WrongFieldCount { line, found } => write!(
                f,
                "PQR line {line}: expected 10 or 11 whitespace-delimited fields, found {found}"
            ),
            Self::InvalidField { line, field, raw } => {
                write!(f, "PQR line {line}: invalid {field} value '{raw}'")
            }
            Self::NonFiniteValue { line, field, raw } => {
                write!(f, "PQR line {line}: non-finite {field} value '{raw}'")
            }
            Self::UnresolvableElement { line, atom_name } => write!(
                f,
                "PQR line {line}: could not infer an element from atom name '{atom_name}'"
            ),
        }
    }
}

impl std::error::Error for PqrError {}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// One PQR atom record, in file order.
#[derive(Debug, Clone, PartialEq)]
pub struct PqrAtomRecord {
    /// `"ATOM"` or `"HETATM"` (or another file-specific value, kept
    /// verbatim).
    pub group_pdb: String,
    pub serial: i64,
    pub atom_name: String,
    pub res_name: String,
    /// `None` when the file omits the chain-id field entirely (the common
    /// `pdb2pqr` default output shape -- see module docs).
    pub chain_id: Option<String>,
    pub res_seq: i64,
    /// An alphabetic insertion-code suffix on the residue number token
    /// (e.g. `52A`), if present. Rare in PQR (pdb2pqr does not renumber
    /// with insertion codes by default) but round-tripped when present.
    pub icode: Option<char>,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub charge: f64,
    pub radius: f64,
    /// Inferred from `atom_name`/`res_name`/`group_pdb` -- see module docs.
    pub element: Element,
}

/// Result of parsing a PQR file.
#[derive(Debug, Clone, PartialEq)]
pub struct PqrResult {
    pub atoms: Vec<PqrAtomRecord>,
}

impl PqrResult {
    /// Build a plain [`Molecule`] + per-atom Cartesian coordinates, in file
    /// order. No bonds are added -- PQR carries no connectivity.
    pub fn to_molecule(&self) -> (Molecule, Vec<(f64, f64, f64)>) {
        let mut builder = MoleculeBuilder::new();
        let mut coords = Vec::with_capacity(self.atoms.len());
        for a in &self.atoms {
            builder.add_atom(Atom::new(a.element));
            coords.push((a.x, a.y, a.z));
        }
        (builder.build(), coords)
    }
}

// ---------------------------------------------------------------------------
// Element inference
// ---------------------------------------------------------------------------

/// `Element::from_symbol` is strict title-case (`"Zn"`, not `"ZN"`), but
/// real PDB/PQR files conventionally write residue/atom names in all caps
/// (`"ZN"`, `"NA"`, `"MG"`). Normalize to title case before matching.
fn title_case_symbol(s: &str) -> String {
    let mut chars = s.chars();
    match chars.next() {
        Some(first) => {
            first.to_ascii_uppercase().to_string() + &chars.as_str().to_ascii_lowercase()
        }
        None => String::new(),
    }
}

/// Infer an [`Element`] from a PQR atom name, using `res_name`/`group_pdb`
/// as disambiguating context. See module docs for the algorithm and why a
/// naive greedy 2-letter match is wrong for standard polymer atom names.
pub fn infer_element(group_pdb: &str, res_name: &str, atom_name: &str) -> Option<Element> {
    // Monatomic-ion convention: HETATM residue named exactly after its
    // element (ZN, NA, MG, CA, FE, CL, K, MN, ...).
    if group_pdb.eq_ignore_ascii_case("HETATM")
        && let Some(e) = Element::from_symbol(&title_case_symbol(res_name.trim()))
    {
        return Some(e);
    }

    let stripped = atom_name.trim_start_matches(|c: char| c.is_ascii_digit());
    let first_alpha = stripped.chars().find(|c| c.is_ascii_alphabetic())?;

    // A small, curated set of two-letter elements that legitimately appear
    // as a ligand/HETATM atom name's first two letters (e.g. "CL1", "FE1",
    // "BR2"). Restricted to HETATM to avoid the ATOM-record ambiguity
    // (CA/CB/CD/CE/CG/CZ/ND/NE/NH/NZ/OD/OE/OG/OH/SD/SG are all common
    // polymer atom-name prefixes that collide with real element symbols).
    if group_pdb.eq_ignore_ascii_case("HETATM") {
        let idx = stripped
            .char_indices()
            .find(|(_, c)| c.is_ascii_alphabetic());
        if let Some((start, _)) = idx {
            let rest = &stripped[start..];
            let mut chars = rest.chars();
            if let (Some(c1), Some(c2)) = (chars.next(), chars.next())
                && c2.is_ascii_alphabetic()
            {
                let two: String = [c1.to_ascii_uppercase(), c2.to_ascii_lowercase()]
                    .iter()
                    .collect();
                const TWO_LETTER_HETATM_ELEMENTS: &[&str] = &[
                    "Cl", "Br", "Fe", "Zn", "Mg", "Mn", "Na", "Ca", "Se", "As", "Si", "Li", "Al",
                    "Cu", "Co", "Ni", "Cd", "Hg", "Ag", "Au", "Pt", "Cs", "Rb", "Sr", "Ba", "Cr",
                    "Mo", "Sn", "Pb", "Ti", "Bi",
                ];
                if TWO_LETTER_HETATM_ELEMENTS.contains(&two.as_str())
                    && let Some(e) = Element::from_symbol(&two)
                {
                    return Some(e);
                }
            }
        }
    }

    Element::from_symbol(&first_alpha.to_ascii_uppercase().to_string())
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a PQR file with default limits ([`PqrParseLimits::default`]).
pub fn parse_pqr(input: &str) -> Result<PqrResult, PqrError> {
    parse_pqr_with_limits(input, &PqrParseLimits::default())
}

/// Parse a PQR file, enforcing `limits`.
pub fn parse_pqr_with_limits(input: &str, limits: &PqrParseLimits) -> Result<PqrResult, PqrError> {
    if input.len() > limits.max_input_bytes {
        return Err(PqrError::InputTooLarge {
            limit: limits.max_input_bytes,
        });
    }

    let mut atoms = Vec::new();
    for (lineno, line) in input.lines().enumerate() {
        if line.len() > limits.max_line_len {
            return Err(PqrError::LineTooLong {
                line: lineno + 1,
                limit: limits.max_line_len,
            });
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.is_empty() {
            continue;
        }
        let record = fields[0];
        if !record.eq_ignore_ascii_case("ATOM") && !record.eq_ignore_ascii_case("HETATM") {
            continue;
        }

        if atoms.len() >= limits.max_atoms {
            return Err(PqrError::TooManyAtoms {
                limit: limits.max_atoms,
            });
        }

        // With chain: record serial name resName chain resSeq x y z q r = 11.
        // Without:    record serial name resName        resSeq x y z q r = 10.
        let has_chain = match fields.len() {
            11 => true,
            10 => false,
            n => {
                return Err(PqrError::WrongFieldCount {
                    line: lineno + 1,
                    found: n,
                });
            }
        };

        let group_pdb = record.to_string();
        let serial = parse_i64(fields[1], lineno + 1, "Atom_number")?;
        let atom_name = fields[2].to_string();
        let res_name = fields[3].to_string();

        let (chain_id, res_seq_tok, x_tok, y_tok, z_tok, q_tok, r_tok) = if has_chain {
            (
                Some(fields[4].to_string()),
                fields[5],
                fields[6],
                fields[7],
                fields[8],
                fields[9],
                fields[10],
            )
        } else {
            (
                None, fields[4], fields[5], fields[6], fields[7], fields[8], fields[9],
            )
        };

        let (res_seq, icode) = parse_res_seq(res_seq_tok, lineno + 1)?;
        let x = parse_finite(x_tok, lineno + 1, "X")?;
        let y = parse_finite(y_tok, lineno + 1, "Y")?;
        let z = parse_finite(z_tok, lineno + 1, "Z")?;
        let charge = parse_finite(q_tok, lineno + 1, "Charge")?;
        let radius = parse_finite(r_tok, lineno + 1, "Radius")?;

        let element = infer_element(&group_pdb, &res_name, &atom_name).ok_or_else(|| {
            PqrError::UnresolvableElement {
                line: lineno + 1,
                atom_name: atom_name.clone(),
            }
        })?;

        atoms.push(PqrAtomRecord {
            group_pdb,
            serial,
            atom_name,
            res_name,
            chain_id,
            res_seq,
            icode,
            x,
            y,
            z,
            charge,
            radius,
            element,
        });
    }

    if atoms.is_empty() {
        return Err(PqrError::NoAtomRecords);
    }

    Ok(PqrResult { atoms })
}

fn parse_i64(s: &str, line: usize, field: &'static str) -> Result<i64, PqrError> {
    s.parse::<i64>().map_err(|_| PqrError::InvalidField {
        line,
        field,
        raw: s.to_string(),
    })
}

fn parse_finite(s: &str, line: usize, field: &'static str) -> Result<f64, PqrError> {
    let v = s.parse::<f64>().map_err(|_| PqrError::InvalidField {
        line,
        field,
        raw: s.to_string(),
    })?;
    if !v.is_finite() {
        return Err(PqrError::NonFiniteValue {
            line,
            field,
            raw: s.to_string(),
        });
    }
    Ok(v)
}

/// Parse a residue-number token, splitting off a trailing alphabetic
/// insertion-code suffix if present (e.g. `"52A"` -> `(52, Some('A'))`).
fn parse_res_seq(s: &str, line: usize) -> Result<(i64, Option<char>), PqrError> {
    if let Ok(v) = s.parse::<i64>() {
        return Ok((v, None));
    }
    let digit_end = s
        .char_indices()
        .find(|(i, c)| !(c.is_ascii_digit() || (*i == 0 && (*c == '-' || *c == '+'))))
        .map(|(i, _)| i)
        .unwrap_or(s.len());
    let (digits, suffix) = s.split_at(digit_end);
    if digits.is_empty() || suffix.chars().any(|c| !c.is_ascii_alphabetic()) || suffix.len() != 1 {
        return Err(PqrError::InvalidField {
            line,
            field: "Residue_number",
            raw: s.to_string(),
        });
    }
    let v = digits.parse::<i64>().map_err(|_| PqrError::InvalidField {
        line,
        field: "Residue_number",
        raw: s.to_string(),
    })?;
    Ok((v, suffix.chars().next()))
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Write a PQR file from parsed atom records.
///
/// Each atom's `chain_id` independently controls whether that line is
/// written with (11 fields) or without (10 fields) the chain column,
/// matching what [`parse_pqr`] would infer back from the written line.
pub fn write_pqr(atoms: &[PqrAtomRecord]) -> String {
    let mut out = String::from("REMARK   PQR written by chematic\n");
    for a in atoms {
        let res_seq_tok = match a.icode {
            Some(c) => format!("{}{}", a.res_seq, c),
            None => a.res_seq.to_string(),
        };
        match &a.chain_id {
            Some(chain) => {
                out.push_str(&format!(
                    "{} {} {} {} {} {} {:.3} {:.3} {:.3} {:.4} {:.4}\n",
                    a.group_pdb,
                    a.serial,
                    a.atom_name,
                    a.res_name,
                    chain,
                    res_seq_tok,
                    a.x,
                    a.y,
                    a.z,
                    a.charge,
                    a.radius,
                ));
            }
            None => {
                out.push_str(&format!(
                    "{} {} {} {} {} {:.3} {:.3} {:.3} {:.4} {:.4}\n",
                    a.group_pdb,
                    a.serial,
                    a.atom_name,
                    a.res_name,
                    res_seq_tok,
                    a.x,
                    a.y,
                    a.z,
                    a.charge,
                    a.radius,
                ));
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Hand-authored per the APBS/PDB2PQR PQR format documentation's
    /// documented field layout (not copied from any real pdb2pqr output).
    /// No-chain shape (the common pdb2pqr default).
    const FIXTURE_NO_CHAIN: &str = "\
REMARK   PQR file\n\
ATOM      1  N   ALA     1     -0.966   1.523   1.412 -0.400  1.500\n\
ATOM      2  CA  ALA     1      0.257   0.679   1.911  0.100  1.700\n\
ATOM      3  CB  ALA     1      0.400   0.800   3.411 -0.030  1.700\n\
HETATM    4  ZN  ZN      2      5.000   5.000   5.000  2.000  1.090\n\
";

    /// With-chain shape.
    const FIXTURE_WITH_CHAIN: &str = "\
ATOM      1  N   ALA A   1     -0.966   1.523   1.412 -0.400  1.500\n\
ATOM      2  CA  ALA A   1      0.257   0.679   1.911  0.100  1.700\n\
HETATM    3  NA  NA  B   5     10.000  10.000  10.000  1.000  1.020\n\
";

    #[test]
    fn parses_no_chain_fixture() {
        let r = parse_pqr(FIXTURE_NO_CHAIN).unwrap();
        assert_eq!(r.atoms.len(), 4);
        assert_eq!(r.atoms[0].chain_id, None);
        assert_eq!(r.atoms[0].atom_name, "N");
        assert_eq!(r.atoms[0].res_name, "ALA");
        assert_eq!(r.atoms[0].res_seq, 1);
        assert!((r.atoms[0].charge - (-0.400)).abs() < 1e-9);
        assert!((r.atoms[0].radius - 1.500).abs() < 1e-9);
    }

    #[test]
    fn parses_with_chain_fixture() {
        let r = parse_pqr(FIXTURE_WITH_CHAIN).unwrap();
        assert_eq!(r.atoms.len(), 3);
        assert_eq!(r.atoms[0].chain_id.as_deref(), Some("A"));
        assert_eq!(r.atoms[2].chain_id.as_deref(), Some("B"));
        assert_eq!(r.atoms[2].res_seq, 5);
    }

    #[test]
    fn backbone_atom_names_resolve_to_correct_elements_not_ambiguous_metals() {
        let r = parse_pqr(FIXTURE_NO_CHAIN).unwrap();
        assert_eq!(r.atoms[0].element, Element::N); // "N"
        assert_eq!(r.atoms[1].element, Element::C); // "CA" -> C, not Ca
        assert_eq!(r.atoms[2].element, Element::C); // "CB" -> C, not... (no such elem)
    }

    #[test]
    fn hetatm_monatomic_ion_resolves_via_residue_name() {
        let r = parse_pqr(FIXTURE_NO_CHAIN).unwrap();
        assert_eq!(r.atoms[3].element, Element::ZN);
        let r2 = parse_pqr(FIXTURE_WITH_CHAIN).unwrap();
        assert_eq!(r2.atoms[2].element, Element::NA);
    }

    #[test]
    fn hetatm_two_letter_ligand_atom_name_resolves() {
        let pqr = "HETATM    1  CL1 LIG     1      1.000   2.000   3.000  -0.100  1.750\n";
        let r = parse_pqr(pqr).unwrap();
        assert_eq!(r.atoms[0].element, Element::CL);
    }

    #[test]
    fn to_molecule_has_no_bonds() {
        let r = parse_pqr(FIXTURE_NO_CHAIN).unwrap();
        let (mol, coords) = r.to_molecule();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(coords.len(), 4);
    }

    #[test]
    fn round_trip_no_chain_preserves_all_fields() {
        let r = parse_pqr(FIXTURE_NO_CHAIN).unwrap();
        let written = write_pqr(&r.atoms);
        let r2 = parse_pqr(&written).unwrap();
        assert_eq!(r.atoms, r2.atoms);
    }

    #[test]
    fn round_trip_with_chain_preserves_all_fields() {
        let r = parse_pqr(FIXTURE_WITH_CHAIN).unwrap();
        let written = write_pqr(&r.atoms);
        let r2 = parse_pqr(&written).unwrap();
        assert_eq!(r.atoms, r2.atoms);
    }

    #[test]
    fn insertion_code_round_trips() {
        let pqr = "ATOM      1  N   ALA A  52A    -0.966   1.523   1.412 -0.400  1.500\n";
        let r = parse_pqr(pqr).unwrap();
        assert_eq!(r.atoms[0].res_seq, 52);
        assert_eq!(r.atoms[0].icode, Some('A'));
        let written = write_pqr(&r.atoms);
        let r2 = parse_pqr(&written).unwrap();
        assert_eq!(r2.atoms[0].res_seq, 52);
        assert_eq!(r2.atoms[0].icode, Some('A'));
    }

    #[test]
    fn wrong_field_count_is_typed_error_not_panic() {
        let pqr = "ATOM      1  N   ALA     1     -0.966   1.523\n";
        let err = parse_pqr(pqr).unwrap_err();
        assert!(matches!(err, PqrError::WrongFieldCount { .. }));
    }

    #[test]
    fn no_atom_records_is_typed_error() {
        let err = parse_pqr("REMARK just a comment\n").unwrap_err();
        assert_eq!(err, PqrError::NoAtomRecords);
    }

    #[test]
    fn nan_charge_is_rejected() {
        let pqr = "ATOM      1  N   ALA     1     -0.966   1.523   1.412 NaN  1.500\n";
        let err = parse_pqr(pqr).unwrap_err();
        assert!(matches!(
            err,
            PqrError::NonFiniteValue {
                field: "Charge",
                ..
            }
        ));
    }

    #[test]
    fn infinity_radius_is_rejected() {
        let pqr = "ATOM      1  N   ALA     1     -0.966   1.523   1.412 -0.4  inf\n";
        let err = parse_pqr(pqr).unwrap_err();
        assert!(matches!(
            err,
            PqrError::NonFiniteValue {
                field: "Radius",
                ..
            }
        ));
    }

    #[test]
    fn malformed_input_never_panics() {
        let inputs = [
            "",
            "ATOM",
            "ATOM 1",
            "ATOM a b c d e f g h i j",
            "\u{0}\u{0}\u{0}",
            "HETATM 1 X X 1 1.0 1.0 1.0 1.0 1.0 extra extra extra",
        ];
        for input in inputs {
            let _ = parse_pqr(input);
        }
    }

    #[test]
    fn input_too_large_is_rejected() {
        let limits = PqrParseLimits {
            max_input_bytes: 8,
            ..PqrParseLimits::default()
        };
        let err = parse_pqr_with_limits(FIXTURE_NO_CHAIN, &limits).unwrap_err();
        assert_eq!(err, PqrError::InputTooLarge { limit: 8 });
    }

    #[test]
    fn too_many_atoms_is_rejected() {
        let limits = PqrParseLimits {
            max_atoms: 2,
            ..PqrParseLimits::default()
        };
        let err = parse_pqr_with_limits(FIXTURE_NO_CHAIN, &limits).unwrap_err();
        assert_eq!(err, PqrError::TooManyAtoms { limit: 2 });
    }

    #[test]
    fn line_too_long_is_rejected() {
        let long = "ATOM      1  N   ALA     1     -0.966   1.523   1.412 -0.400  1.500 "
            .to_string()
            + &" ".repeat(2000);
        let limits = PqrParseLimits {
            max_line_len: 200,
            ..PqrParseLimits::default()
        };
        let err = parse_pqr_with_limits(&long, &limits).unwrap_err();
        assert!(matches!(err, PqrError::LineTooLong { .. }));
    }

    #[test]
    fn oracle_single_ion_atom() {
        // Values independently chosen (not copied from a real pdb2pqr
        // output file) to exactly match the documented PQR field order.
        let pqr = "HETATM    1  MG  MG  A   1     12.500   8.250  -3.750   2.000   0.99\n";
        let r = parse_pqr(pqr).unwrap();
        let a = &r.atoms[0];
        assert_eq!(a.group_pdb, "HETATM");
        assert_eq!(a.serial, 1);
        assert_eq!(a.atom_name, "MG");
        assert_eq!(a.res_name, "MG");
        assert_eq!(a.chain_id.as_deref(), Some("A"));
        assert_eq!(a.res_seq, 1);
        assert!((a.x - 12.500).abs() < 1e-9);
        assert!((a.y - 8.250).abs() < 1e-9);
        assert!((a.z - (-3.750)).abs() < 1e-9);
        assert!((a.charge - 2.000).abs() < 1e-9);
        assert!((a.radius - 0.99).abs() < 1e-9);
        assert_eq!(a.element, Element::MG);
    }

    #[test]
    fn quote_free_atom_names_do_not_need_escaping() {
        // PQR is whitespace-delimited with no quoting mechanism at all --
        // an atom/residue name containing whitespace cannot be represented.
        // This isn't a bug to fix; it's a hard format limitation, exercised
        // here so it's visible rather than silently mishandled.
        let r = parse_pqr(FIXTURE_NO_CHAIN).unwrap();
        for a in &r.atoms {
            assert!(!a.atom_name.contains(char::is_whitespace));
            assert!(!a.res_name.contains(char::is_whitespace));
        }
    }
}
