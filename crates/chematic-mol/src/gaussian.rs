//! Gaussian 16 / Gaussian 09 file format support.
//!
//! ## Supported formats
//!
//! | Format | Read | Write | Description |
//! |--------|------|-------|-------------|
//! | `.gjf` / `.com` | ✅ | ✅ | Gaussian input file |
//! | `.log` / `.out` | ✅ | — | Gaussian output file (geometry + energy) |
//!
//! ## GJF structure
//!
//! ```text
//! %chk=ethanol.chk
//! # B3LYP/6-31G* opt
//!
//! Ethanol geometry optimization
//!
//! 0 1
//! C   0.000000  0.000000  0.000000
//! C   1.531000  0.000000  0.000000
//! O   2.058000  1.198000  0.000000
//! H   ...
//! ```
//!
//! Coordinates are returned as `Vec<(f64, f64, f64)>` (Ångströms, same
//! convention as `chematic-mol`'s other parsers).

use chematic_core::{Atom, AtomIdx, Element, Molecule, MoleculeBuilder};

/// Return type shared by the GJF parser and its coordinate helper:
/// `(Molecule, coords, charge, multiplicity)`.
type GjfResult = (Molecule, Vec<(f64, f64, f64)>, i32, u32);

/// Resource limits for Gaussian input and log parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GaussianParseLimits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_sections: usize,
    pub max_atoms: usize,
}

impl Default for GaussianParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 * 1024 * 1024,
            max_line_bytes: 1024 * 1024,
            max_lines: 1_000_000,
            max_sections: 100_000,
            max_atoms: 1_000_000,
        }
    }
}

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Errors that can occur when parsing Gaussian files.
#[derive(Debug, Clone, PartialEq)]
pub enum GaussianError {
    /// The charge/multiplicity line was missing or malformed.
    MissingChargeMultiplicity,
    /// No atom coordinates were found in the file.
    NoAtoms,
    /// An element symbol or atomic number was not recognised.
    UnknownElement(String),
    /// A coordinate value could not be parsed as a float.
    InvalidCoordinate(String),
    /// The log file contained no `Standard orientation:` block.
    NoStandardOrientation,
    /// A parsed coordinate was NaN or infinite.
    NonFiniteCoordinate(String),
    /// The input exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl core::fmt::Display for GaussianError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MissingChargeMultiplicity => {
                write!(f, "charge/multiplicity line not found in GJF")
            }
            Self::NoAtoms => write!(f, "no atom coordinates found in file"),
            Self::UnknownElement(s) => write!(f, "unknown element '{s}'"),
            Self::InvalidCoordinate(s) => write!(f, "invalid coordinate value '{s}'"),
            Self::NoStandardOrientation => {
                write!(f, "no 'Standard orientation' block found in Gaussian log")
            }
            Self::NonFiniteCoordinate(s) => write!(f, "non-finite coordinate '{s}'"),
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(f, "{resource} has size {actual}, exceeding limit {limit}"),
        }
    }
}

impl std::error::Error for GaussianError {}

// ---------------------------------------------------------------------------
// Result type for LOG files
// ---------------------------------------------------------------------------

/// Data extracted from a Gaussian output (.log / .out) file.
pub struct GaussianLogResult {
    /// Molecular topology (atoms only; no bonds).
    pub mol: Molecule,
    /// Atomic coordinates in Ångströms `(x, y, z)`.
    pub coords: Vec<(f64, f64, f64)>,
    /// Final SCF energy in hartree, if present.
    pub scf_energy: Option<f64>,
}

// ---------------------------------------------------------------------------
// GJF parser
// ---------------------------------------------------------------------------

/// Parse a Gaussian input file (`.gjf` / `.com`).
///
/// Returns `(Molecule, coords, charge, multiplicity)`.
/// `coords` is a `Vec<(f64, f64, f64)>` in Ångströms.
/// The `Molecule` has no bonds; use `determine_bonds` to infer connectivity.
pub fn parse_gjf(input: &str) -> Result<GjfResult, GaussianError> {
    parse_gjf_with_limits(input, &GaussianParseLimits::default())
}

/// Parse a Gaussian input file with explicit resource limits.
pub fn parse_gjf_with_limits(
    input: &str,
    limits: &GaussianParseLimits,
) -> Result<GjfResult, GaussianError> {
    check_input_limits(input, limits)?;
    // Split into sections at blank lines.
    let sections: Vec<Vec<&str>> = {
        let mut secs: Vec<Vec<&str>> = Vec::new();
        let mut current: Vec<&str> = Vec::new();
        for line in input.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                if !current.is_empty() {
                    secs.push(current.clone());
                    current.clear();
                }
            } else {
                current.push(trimmed);
            }
        }
        if !current.is_empty() {
            secs.push(current);
        }
        secs
    };

    // GJF structure: [%directives +] #route [blank] title [blank] charge/mult coords...
    // Find the route section (a section containing a line that starts with '#'),
    // then the charge/multiplicity section is two sections later.
    // This avoids falsely matching a title line that happens to be "0 1".
    let route_idx = sections
        .iter()
        .position(|sec| sec.iter().any(|l| l.starts_with('#')))
        .ok_or(GaussianError::MissingChargeMultiplicity)?;
    let cm_idx = route_idx + 2;
    if cm_idx >= sections.len() {
        return Err(GaussianError::MissingChargeMultiplicity);
    }
    // Verify the section looks like a charge/multiplicity line.
    {
        let parts: Vec<&str> = sections[cm_idx][0].split_whitespace().collect();
        if parts.len() < 2 || parts[0].parse::<i32>().is_err() || parts[1].parse::<u32>().is_err() {
            return Err(GaussianError::MissingChargeMultiplicity);
        }
    }

    let parts: Vec<&str> = sections[cm_idx][0].split_whitespace().collect();
    let charge: i32 = parts[0].parse().unwrap();
    let multiplicity: u32 = parts[1].parse().unwrap();

    // Coordinates follow the charge line in the same or next section.
    let coord_lines: &[&str] = if sections[cm_idx].len() > 1 {
        &sections[cm_idx][1..]
    } else if let Some(next) = sections.get(cm_idx + 1) {
        next.as_slice()
    } else {
        return Err(GaussianError::NoAtoms);
    };

    if sections.len() > limits.max_sections {
        return Err(GaussianError::ResourceLimit {
            resource: "sections",
            actual: sections.len(),
            limit: limits.max_sections,
        });
    }
    parse_atom_coords(coord_lines, charge, multiplicity, limits.max_atoms)
}

fn parse_atom_coords(
    lines: &[&str],
    charge: i32,
    multiplicity: u32,
    max_atoms: usize,
) -> Result<GjfResult, GaussianError> {
    let mut builder = MoleculeBuilder::new();
    let mut coords: Vec<(f64, f64, f64)> = Vec::new();

    for line in lines {
        let line = line.trim();
        if line.is_empty() || line.starts_with('!') {
            continue;
        }
        if line.starts_with("Variables:") || line.starts_with("Constants:") {
            break;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let raw_sym = parts[0].trim_end_matches(|c: char| c.is_ascii_digit());
        let elem = if raw_sym.is_empty() {
            // Bare atomic number (e.g. "6" for carbon — valid Gaussian input).
            let atomic_num: u8 = parts[0]
                .parse()
                .map_err(|_| GaussianError::UnknownElement(parts[0].to_string()))?;
            Element::from_atomic_number(atomic_num)
                .ok_or_else(|| GaussianError::UnknownElement(parts[0].to_string()))?
        } else {
            Element::from_symbol(raw_sym)
                .ok_or_else(|| GaussianError::UnknownElement(raw_sym.to_string()))?
        };

        let x: f64 = parts[1]
            .parse()
            .map_err(|_| GaussianError::InvalidCoordinate(parts[1].to_string()))?;
        let y: f64 = parts[2]
            .parse()
            .map_err(|_| GaussianError::InvalidCoordinate(parts[2].to_string()))?;
        let z: f64 = parts[3]
            .parse()
            .map_err(|_| GaussianError::InvalidCoordinate(parts[3].to_string()))?;
        if !x.is_finite() {
            return Err(GaussianError::NonFiniteCoordinate(parts[1].to_string()));
        }
        if !y.is_finite() {
            return Err(GaussianError::NonFiniteCoordinate(parts[2].to_string()));
        }
        if !z.is_finite() {
            return Err(GaussianError::NonFiniteCoordinate(parts[3].to_string()));
        }
        if coords.len() >= max_atoms {
            return Err(GaussianError::ResourceLimit {
                resource: "atoms",
                actual: coords.len() + 1,
                limit: max_atoms,
            });
        }

        builder.add_atom(Atom::new(elem));
        coords.push((x, y, z));
    }

    if coords.is_empty() {
        return Err(GaussianError::NoAtoms);
    }

    Ok((builder.build(), coords, charge, multiplicity))
}

// ---------------------------------------------------------------------------
// GJF writer
// ---------------------------------------------------------------------------

/// Write a Gaussian input file (`.gjf`) string.
///
/// - `method`: route section keywords, e.g. `"B3LYP/6-31G* opt"`.
/// - `title`: job title comment (uses `"chematic"` if empty).
pub fn write_gjf(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
    charge: i32,
    multiplicity: u32,
    method: &str,
    title: &str,
) -> String {
    let method = if method.is_empty() {
        "B3LYP/6-31G* opt"
    } else {
        method
    };
    let title = if title.is_empty() { "chematic" } else { title };

    let mut out = format!("# {method}\n\n{title}\n\n{charge} {multiplicity}\n");

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let sym = mol.atom(idx).element.symbol();
        let (x, y, z) = coords.get(i).copied().unwrap_or((0.0, 0.0, 0.0));
        out.push_str(&format!("{sym:<3}  {:12.6}  {:12.6}  {:12.6}\n", x, y, z));
    }
    out.push('\n');
    out
}

// ---------------------------------------------------------------------------
// Gaussian LOG parser
// ---------------------------------------------------------------------------

/// Parse a Gaussian output file (`.log` / `.out`).
///
/// Extracts:
/// - The **last** `Standard orientation:` geometry block.
/// - The **last** `SCF Done:` energy in hartree (if present).
///
/// Returns [`GaussianLogResult`] with no bond information.
pub fn parse_gaussian_log(input: &str) -> Result<GaussianLogResult, GaussianError> {
    parse_gaussian_log_with_limits(input, &GaussianParseLimits::default())
}

/// Parse a Gaussian output log with explicit resource limits.
pub fn parse_gaussian_log_with_limits(
    input: &str,
    limits: &GaussianParseLimits,
) -> Result<GaussianLogResult, GaussianError> {
    check_input_limits(input, limits)?;
    let lines: Vec<&str> = input.lines().collect();

    // Locate last Standard orientation: block.
    let last_orient_start = lines
        .iter()
        .rposition(|l| l.contains("Standard orientation:"))
        .ok_or(GaussianError::NoStandardOrientation)?;

    let mut builder = MoleculeBuilder::new();
    let mut coords: Vec<(f64, f64, f64)> = Vec::new();
    let mut dashes_seen = 0usize;
    let mut in_table = false;

    for line in &lines[last_orient_start + 1..] {
        let trimmed = line.trim();
        if trimmed.starts_with("---") {
            dashes_seen += 1;
            if dashes_seen == 2 {
                in_table = true; // after 2nd dashes line
            } else if in_table {
                break; // closing dashes = end of table
            }
            continue;
        }
        if !in_table {
            continue;
        }
        // Gaussian 09/16: Center_Num Atomic_Num Atomic_Type X Y Z  (6 cols)
        // Gaussian 03:    Center_Num Atomic_Num X Y Z             (5 cols)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        let (an_col, x_col) = match parts.len() {
            n if n >= 6 => (1, 3), // standard 6-column format
            5 => (1, 2),           // Gaussian 03 5-column format (no Atomic_Type)
            _ => continue,
        };
        let atomic_num: u8 = parts[an_col]
            .parse()
            .map_err(|_| GaussianError::UnknownElement(parts[an_col].to_string()))?;
        let elem = Element::from_atomic_number(atomic_num)
            .ok_or_else(|| GaussianError::UnknownElement(parts[an_col].to_string()))?;
        let x: f64 = parts[x_col]
            .parse()
            .map_err(|_| GaussianError::InvalidCoordinate(parts[x_col].to_string()))?;
        let y: f64 = parts[x_col + 1]
            .parse()
            .map_err(|_| GaussianError::InvalidCoordinate(parts[x_col + 1].to_string()))?;
        let z: f64 = parts[x_col + 2]
            .parse()
            .map_err(|_| GaussianError::InvalidCoordinate(parts[x_col + 2].to_string()))?;
        if !x.is_finite() || !y.is_finite() || !z.is_finite() {
            return Err(GaussianError::NonFiniteCoordinate(
                parts[x_col..x_col + 3].join(" "),
            ));
        }
        if coords.len() >= limits.max_atoms {
            return Err(GaussianError::ResourceLimit {
                resource: "atoms",
                actual: coords.len() + 1,
                limit: limits.max_atoms,
            });
        }

        builder.add_atom(Atom::new(elem));
        coords.push((x, y, z));
    }

    if coords.is_empty() {
        return Err(GaussianError::NoAtoms);
    }

    // Extract last SCF Done energy.
    let scf_energy = lines
        .iter()
        .rfind(|l| l.contains("SCF Done:"))
        .and_then(|l| {
            let after = l[l.find('=')? + 1..].trim();
            after.split_whitespace().next()?.parse::<f64>().ok()
        });

    Ok(GaussianLogResult {
        mol: builder.build(),
        coords,
        scf_energy,
    })
}

fn check_input_limits(input: &str, limits: &GaussianParseLimits) -> Result<(), GaussianError> {
    if input.len() > limits.max_input_bytes {
        return Err(GaussianError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let lines = input.lines().count();
    if lines > limits.max_lines {
        return Err(GaussianError::ResourceLimit {
            resource: "lines",
            actual: lines,
            limit: limits.max_lines,
        });
    }
    if let Some(line_bytes) = input.lines().map(str::len).max()
        && line_bytes > limits.max_line_bytes
    {
        return Err(GaussianError::ResourceLimit {
            resource: "line bytes",
            actual: line_bytes,
            limit: limits.max_line_bytes,
        });
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const ETHANOL_GJF: &str = r#"# B3LYP/6-31G* opt

Ethanol

0 1
C   0.000000   0.000000   0.000000
C   1.531000   0.000000   0.000000
O   2.058000   1.198000   0.000000
H  -0.390000   1.020000   0.000000
H  -0.390000  -0.510000   0.884000
H  -0.390000  -0.510000  -0.884000
H   1.921000  -1.020000   0.000000
H   1.921000   0.510000  -0.884000
H   2.028000   1.604000   0.886000

"#;

    const LOG_SNIPPET: &str = r#"
 Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0        0.000000    0.000000    0.000000
      2          6           0        1.531000    0.000000    0.000000
      3          8           0        2.058000    1.198000    0.000000
 ---------------------------------------------------------------------
 SCF Done:  E(RB3LYP) =  -154.0987765432  A.U. after   12 cycles
"#;

    #[test]
    fn parse_gjf_atom_count() {
        let (mol, coords, charge, mult) = parse_gjf(ETHANOL_GJF).unwrap();
        assert_eq!(mol.atom_count(), 9);
        assert_eq!(coords.len(), 9);
        assert_eq!(charge, 0);
        assert_eq!(mult, 1);
    }

    #[test]
    fn parse_gjf_first_coord() {
        let (_, coords, _, _) = parse_gjf(ETHANOL_GJF).unwrap();
        let (x, y, z) = coords[0];
        assert!(x.abs() < 1e-6 && y.abs() < 1e-6 && z.abs() < 1e-6);
    }

    #[test]
    fn write_gjf_roundtrip() {
        let (mol, coords, charge, mult) = parse_gjf(ETHANOL_GJF).unwrap();
        let out = write_gjf(&mol, &coords, charge, mult, "B3LYP/6-31G* opt", "Ethanol");
        assert!(out.contains("# B3LYP/6-31G* opt"));
        assert!(out.contains("0 1"));
        let (mol2, coords2, _, _) = parse_gjf(&out).unwrap();
        assert_eq!(mol2.atom_count(), 9);
        assert_eq!(coords2.len(), 9);
    }

    #[test]
    fn parse_gjf_charged_molecule() {
        let gjf =
            "# HF/STO-3G\n\nwater cation\n\n1 2\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n";
        let (mol, _, charge, mult) = parse_gjf(gjf).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(charge, 1);
        assert_eq!(mult, 2);
    }

    #[test]
    fn parse_log_geometry() {
        let r = parse_gaussian_log(LOG_SNIPPET).unwrap();
        assert_eq!(r.mol.atom_count(), 3);
        assert_eq!(r.coords.len(), 3);
        let (x, y, z) = r.coords[0];
        assert!(x.abs() < 1e-6 && y.abs() < 1e-6 && z.abs() < 1e-6);
    }

    #[test]
    fn parse_log_energy() {
        let r = parse_gaussian_log(LOG_SNIPPET).unwrap();
        let e = r.scf_energy.expect("energy should be present");
        assert!((e - (-154.0987765432)).abs() < 1e-8);
    }

    #[test]
    fn parse_log_last_orientation() {
        // When multiple Standard orientation blocks exist, last one wins.
        let input = format!(
            "{}\n Standard orientation:\n ---\n ---\n      1          6           0        9.0    9.0    9.0\n ---\n{}",
            LOG_SNIPPET, ""
        );
        let r = parse_gaussian_log(&input).unwrap();
        // Last block has atom at (9, 9, 9)
        assert!((r.coords[0].0 - 9.0).abs() < 1e-6);
    }

    #[test]
    fn bounded_gaussian_parsers_reject_input_and_line_limits() {
        let limits = GaussianParseLimits {
            max_input_bytes: 8,
            ..Default::default()
        };
        assert!(matches!(
            parse_gjf_with_limits(ETHANOL_GJF, &limits),
            Err(GaussianError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));
        assert!(matches!(
            parse_gaussian_log_with_limits(LOG_SNIPPET, &limits),
            Err(GaussianError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));

        let long_line = format!("{}\n", "x".repeat(32));
        assert!(matches!(
            parse_gjf_with_limits(
                &long_line,
                &GaussianParseLimits {
                    max_line_bytes: 16,
                    ..Default::default()
                }
            ),
            Err(GaussianError::ResourceLimit {
                resource: "line bytes",
                ..
            })
        ));
    }

    #[test]
    fn bounded_gaussian_parsers_reject_atom_limit() {
        assert!(matches!(
            parse_gjf_with_limits(
                ETHANOL_GJF,
                &GaussianParseLimits {
                    max_atoms: 2,
                    ..Default::default()
                }
            ),
            Err(GaussianError::ResourceLimit {
                resource: "atoms",
                ..
            })
        ));
        assert!(matches!(
            parse_gaussian_log_with_limits(
                LOG_SNIPPET,
                &GaussianParseLimits {
                    max_atoms: 2,
                    ..Default::default()
                }
            ),
            Err(GaussianError::ResourceLimit {
                resource: "atoms",
                ..
            })
        ));
    }
}
