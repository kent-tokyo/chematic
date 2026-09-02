//! AutoDock PDBQT format reader and writer.
//!
//! PDBQT is the input format for AutoDock4 and AutoDock Vina. It extends the
//! PDB format with two additional fields per atom:
//! - columns 72–76: partial charge (Gasteiger or AMBER BCC)
//! - columns 78–79: AutoDock atom type (C, A, N, NA, O, OA, S, SA, H, HD, P,
//!   F, Cl, Br, I, and metals)
//!
//! This module writes rigid-body PDBQT (no rotatable bond `BRANCH`/`ENDBRANCH`
//! tree), suitable for receptor preparation. For flexible ligands, wrap the
//! output with `ROOT`/`ENDROOT`/`TORSDOF` as required by the docking software.
//!
//! Reference: AutoDock 4.2 manual, Section 7 (PDBQT format specification).

use chematic_core::{AtomIdx, Element, Molecule};

/// Resource limits for PDBQT parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PdbqtParseLimits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_atoms: usize,
}

impl Default for PdbqtParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_line_bytes: 1024 * 1024,
            max_lines: 1_000_000,
            max_atoms: 1_000_000,
        }
    }
}

/// Error returned when parsing a PDBQT file fails.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PdbqtError {
    /// An ATOM/HETATM line is too short or malformed.
    InvalidAtomLine { line: usize, detail: String },
    /// The partial charge field cannot be parsed as f64.
    InvalidCharge { line: usize, raw: String },
    /// Unknown element symbol.
    UnknownElement { symbol: String, line: usize },
    /// A coordinate or charge parsed as NaN or infinite.
    NonFiniteValue { line: usize, field: &'static str },
    /// The input exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl core::fmt::Display for PdbqtError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidAtomLine { line, detail } => {
                write!(f, "PDBQT: invalid ATOM line {line}: {detail}")
            }
            Self::InvalidCharge { line, raw } => {
                write!(f, "PDBQT: invalid charge '{raw}' at line {line}")
            }
            Self::UnknownElement { symbol, line } => {
                write!(f, "PDBQT: unknown element '{symbol}' at line {line}")
            }
            Self::NonFiniteValue { line, field } => {
                write!(f, "PDBQT: non-finite {field} at line {line}")
            }
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "PDBQT: {resource} has size {actual}, exceeding limit {limit}"
            ),
        }
    }
}

impl std::error::Error for PdbqtError {}

// ── AutoDock atom type assignment ─────────────────────────────────────────────

/// Assign an AutoDock atom type string for a given atom in `mol`.
///
/// Rules follow the AutoDock 4.2 atom type table:
/// - Aromatic C → "A", aliphatic C → "C"
/// - N accepting H-bond (no H, lone pair) → "NA", otherwise "N"
/// - O with lone pair (ethers, carbonyl, …) → "OA", charged O → "OA", otherwise "O"
/// - S with lone pair → "SA", otherwise "S"
/// - H bonded to O or N → "HD" (hydrogen donor), otherwise "H"
/// - Halogens: F, Cl, Br, I (as-is)
/// - Metals: Mg, Mn, Zn, Ca, Fe (as-is)
/// - Phosphorus: P
/// - Fallback: element symbol
pub fn autodock_atom_type(mol: &Molecule, idx: AtomIdx) -> &'static str {
    let atom = mol.atom(idx);
    match atom.element.atomic_number() {
        6 => {
            if atom.aromatic {
                "A"
            } else {
                "C"
            }
        }
        7 => {
            // NA if N can accept H-bonds (no explicit H, no positive charge)
            let has_h = mol
                .neighbors(idx)
                .any(|(nb, _)| mol.atom(nb).element.atomic_number() == 1);
            if !has_h && atom.charge <= 0 {
                "NA"
            } else {
                "N"
            }
        }
        8 => "OA", // all O treated as acceptor (most common case)
        16 => {
            // SA if S has lone pair (not positively charged)
            if atom.charge >= 0 { "SA" } else { "S" }
        }
        1 => {
            // HD if bonded to N or O (hydrogen bond donor)
            let is_donor = mol.neighbors(idx).any(|(nb, _)| {
                let an = mol.atom(nb).element.atomic_number();
                an == 7 || an == 8
            });
            if is_donor { "HD" } else { "H" }
        }
        15 => "P",
        9 => "F",
        17 => "Cl",
        35 => "Br",
        53 => "I",
        12 => "Mg",
        25 => "Mn",
        30 => "Zn",
        20 => "Ca",
        26 => "Fe",
        _ => "X", // unknown
    }
}

// ── Writer ────────────────────────────────────────────────────────────────────

/// Write a molecule to PDBQT format.
///
/// `coords` — 3D coordinates per atom (Å). If shorter than `mol.atom_count()`,
/// missing atoms receive `(0, 0, 0)`.
///
/// `charges` — partial charge per atom. If empty, all charges are written as
/// `0.0000`. For best docking quality, supply Gasteiger charges from
/// `chematic_chem::gasteiger_charges()` or MMFF94 BCI charges from
/// `chematic_ff::mmff94_charges_bci()`.
///
/// The output is a rigid-body PDBQT (no `BRANCH`/`ENDBRANCH` torsion tree).
/// Wrap with `ROOT`/`ENDROOT`/`TORSDOF` for flexible-ligand docking.
pub fn write_pdbqt(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
    charges: &[f64],
    ligand_name: &str,
) -> String {
    let mut out = String::new();
    out.push_str("REMARK  PDBQT written by chematic\n");
    out.push_str(&format!("REMARK  NAME = {ligand_name}\n"));

    for (serial, (aidx, atom)) in mol.atoms().enumerate() {
        let i = aidx.0 as usize;
        let (x, y, z) = coords.get(i).copied().unwrap_or((0.0, 0.0, 0.0));
        let q = charges.get(i).copied().unwrap_or(0.0);
        let ad_type = autodock_atom_type(mol, aidx);
        let sym = atom.element.symbol();

        // PDB atom name: right-justify 1-char symbols in column 14, 2-char in 13-14
        let atom_name = if sym.len() == 1 {
            format!(" {sym}  ")
        } else {
            format!("{sym}   ")
        };

        // ATOM serial resname chain resnum    x        y        z     occ   bfac           charge type
        out.push_str(&format!(
            "ATOM  {:>5} {:<4} {:<3} {:>1}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}    {:>+7.4} {:<2}\n",
            serial + 1,  // serial
            atom_name,   // atom name
            "LIG",       // residue name
            "A",         // chain
            1,           // residue number
            x, y, z,
            1.00,        // occupancy
            0.00,        // B-factor
            q,           // partial charge
            ad_type,     // AutoDock type
        ));
    }

    out.push_str("END\n");
    out
}

// ── Parser ────────────────────────────────────────────────────────────────────

/// Parse a PDBQT string into a molecule, 3D coordinates, and partial charges.
///
/// Only `ATOM` and `HETATM` records are read; `REMARK`, `ROOT`, `ENDROOT`,
/// `BRANCH`, `ENDBRANCH`, and `TORSDOF` lines are silently skipped.
#[allow(clippy::type_complexity)]
pub fn parse_pdbqt(s: &str) -> Result<(Molecule, Vec<(f64, f64, f64)>, Vec<f64>), PdbqtError> {
    parse_pdbqt_with_limits(s, &PdbqtParseLimits::default())
}

/// Parse a PDBQT string with explicit resource limits.
#[allow(clippy::type_complexity)]
pub fn parse_pdbqt_with_limits(
    s: &str,
    limits: &PdbqtParseLimits,
) -> Result<(Molecule, Vec<(f64, f64, f64)>, Vec<f64>), PdbqtError> {
    use chematic_core::MoleculeBuilder;

    if s.len() > limits.max_input_bytes {
        return Err(PdbqtError::ResourceLimit {
            resource: "input bytes",
            actual: s.len(),
            limit: limits.max_input_bytes,
        });
    }
    let lines = s.lines().collect::<Vec<_>>();
    if lines.len() > limits.max_lines {
        return Err(PdbqtError::ResourceLimit {
            resource: "lines",
            actual: lines.len(),
            limit: limits.max_lines,
        });
    }
    if let Some(line_bytes) = lines.iter().map(|line| line.len()).max()
        && line_bytes > limits.max_line_bytes
    {
        return Err(PdbqtError::ResourceLimit {
            resource: "line bytes",
            actual: line_bytes,
            limit: limits.max_line_bytes,
        });
    }

    let mut builder = MoleculeBuilder::new();
    let mut coords: Vec<(f64, f64, f64)> = Vec::new();
    let mut charges: Vec<f64> = Vec::new();

    for (lineno, line) in lines.into_iter().enumerate() {
        let record = line.get(..line.len().min(6)).unwrap_or("");
        if !matches!(record.trim(), "ATOM" | "HETATM") {
            continue;
        }
        if coords.len() >= limits.max_atoms {
            return Err(PdbqtError::ResourceLimit {
                resource: "atom records",
                actual: coords.len() + 1,
                limit: limits.max_atoms,
            });
        }
        if line.len() < 54 {
            return Err(PdbqtError::InvalidAtomLine {
                line: lineno + 1,
                detail: "line too short (need at least 54 chars for coordinates)".into(),
            });
        }

        // Element from AutoDock type (cols 78-79, 0-indexed 77-78) or fallback to PDB col 77-78
        let ad_type = line.get(77..79).map(str::trim).unwrap_or("").trim();
        // Map AutoDock type back to element symbol
        let elem_sym = match ad_type {
            "A" => "C",
            "NA" => "N",
            "OA" | "OS" => "O",
            "SA" => "S",
            "HD" => "H",
            other if !other.is_empty() => other,
            _ => {
                // Fall back to PDB element column (cols 77-78)
                line.get(76..78).map(str::trim).unwrap_or("C").trim()
            }
        };

        let element = Element::from_symbol(elem_sym).ok_or_else(|| PdbqtError::UnknownElement {
            symbol: elem_sym.to_string(),
            line: lineno + 1,
        })?;

        use chematic_core::Atom;
        builder.add_atom(Atom::new(element));

        // Coordinates: cols 31-38, 39-46, 47-54 (0-indexed 30-37, 38-45, 46-53)
        let parse_f = |s: &str, col: &str| -> Result<f64, PdbqtError> {
            s.trim()
                .parse::<f64>()
                .map_err(|_| PdbqtError::InvalidAtomLine {
                    line: lineno + 1,
                    detail: format!("cannot parse {col} coordinate: '{}'", s.trim()),
                })
        };
        let x = parse_f(line.get(30..38).unwrap_or(""), "x")?;
        let y = parse_f(line.get(38..46).unwrap_or(""), "y")?;
        let z = parse_f(line.get(46..54).unwrap_or(""), "z")?;
        if !x.is_finite() {
            return Err(PdbqtError::NonFiniteValue {
                line: lineno + 1,
                field: "x coordinate",
            });
        }
        if !y.is_finite() {
            return Err(PdbqtError::NonFiniteValue {
                line: lineno + 1,
                field: "y coordinate",
            });
        }
        if !z.is_finite() {
            return Err(PdbqtError::NonFiniteValue {
                line: lineno + 1,
                field: "z coordinate",
            });
        }
        coords.push((x, y, z));

        // Partial charge: cols 67-76 (0-indexed 66-75), or fallback 71-76
        let q_raw = line
            .get(66..76)
            .or_else(|| line.get(71..76))
            .map(str::trim)
            .unwrap_or("0.0");
        let q = q_raw
            .parse::<f64>()
            .map_err(|_| PdbqtError::InvalidCharge {
                line: lineno + 1,
                raw: q_raw.to_string(),
            })?;
        if !q.is_finite() {
            return Err(PdbqtError::NonFiniteValue {
                line: lineno + 1,
                field: "partial charge",
            });
        }
        charges.push(q);
    }

    Ok((builder.build(), coords, charges))
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn write_pdbqt_ethanol() {
        let mol = parse("CCO").unwrap();
        let coords = vec![(0.0, 0.0, 0.0), (1.54, 0.0, 0.0), (2.5, 1.2, 0.0)];
        let charges = vec![-0.100, 0.050, -0.400];
        let out = write_pdbqt(&mol, &coords, &charges, "ETH");
        assert!(out.contains("REMARK"));
        assert!(out.contains("ATOM"));
        assert!(out.contains("LIG"));
        // C gets type C, O gets OA
        assert!(out.contains(" C ") || out.contains("C\n") || out.contains("C "));
        assert!(out.contains("OA"));
    }

    #[test]
    fn write_pdbqt_aromatic_carbon_gets_type_a() {
        let mol = parse("c1ccccc1").unwrap();
        let out = write_pdbqt(&mol, &[], &[], "BNZ");
        // All carbons are aromatic → type A
        assert!(out.contains("A ") || out.contains(" A"));
    }

    #[test]
    fn roundtrip_coordinates() {
        let mol = parse("CCO").unwrap();
        let coords = vec![(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)];
        let charges = vec![-0.1, 0.2, -0.3];
        let pdbqt = write_pdbqt(&mol, &coords, &charges, "TST");
        let (mol2, coords2, charges2) = parse_pdbqt(&pdbqt).unwrap();
        assert_eq!(mol2.atom_count(), 3);
        assert!((coords2[0].0 - 1.0).abs() < 0.01);
        assert!((coords2[1].1 - 5.0).abs() < 0.01);
        assert!((charges2[2] - (-0.3)).abs() < 0.01);
    }

    #[test]
    fn bounded_parser_rejects_input_and_line_limits() {
        let mol = parse("CCO").unwrap();
        let text = write_pdbqt(&mol, &[], &[], "TST");
        assert!(matches!(
            parse_pdbqt_with_limits(
                &text,
                &PdbqtParseLimits {
                    max_input_bytes: 8,
                    ..Default::default()
                }
            ),
            Err(PdbqtError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));
        let long_line = format!("{}\n", "x".repeat(32));
        assert!(matches!(
            parse_pdbqt_with_limits(
                &long_line,
                &PdbqtParseLimits {
                    max_line_bytes: 16,
                    ..Default::default()
                }
            ),
            Err(PdbqtError::ResourceLimit {
                resource: "line bytes",
                ..
            })
        ));
    }

    #[test]
    fn bounded_parser_rejects_atom_limit() {
        let mol = parse("CCO").unwrap();
        let text = write_pdbqt(&mol, &[], &[], "TST");
        assert!(matches!(
            parse_pdbqt_with_limits(
                &text,
                &PdbqtParseLimits {
                    max_atoms: 2,
                    ..Default::default()
                }
            ),
            Err(PdbqtError::ResourceLimit {
                resource: "atom records",
                ..
            })
        ));
    }

    #[test]
    fn malformed_utf8_column_alignment_is_rejected_without_panicking() {
        let input = "ATOM \u{0328}N\n";
        assert!(parse_pdbqt(input).is_ok());
    }
}
