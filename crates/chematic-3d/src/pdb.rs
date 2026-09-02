//! PDB file format parser and writer.
//!
//! Supports `ATOM` and `HETATM` record types only.  Other record types are
//! silently ignored.  Atoms are read using the standard PDB fixed-column layout
//! (columns are 1-indexed in the PDB specification; below we use 0-indexed
//! Rust slice indices).
//!
//! Column mapping (0-indexed):
//!   [0..6]   record type ("ATOM  " or "HETATM")
//!   [6..11]  serial (right-justified integer)
//!   [12..16] atom name
//!   [17..20] residue name
//!   [21]     chain ID
//!   [22..26] residue sequence number
//!   [30..38] x (8.3f)
//!   [38..46] y (8.3f)
//!   [46..54] z (8.3f)
//!   [54..60] occupancy (6.2f)
//!   [60..66] temp factor (6.2f)
//!   [76..78] element symbol

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

use crate::coords::{Coords3D, Point3};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Data from a single ATOM or HETATM record in a PDB file.
#[derive(Debug, Clone)]
pub struct PdbAtom {
    /// Atom serial number.
    pub serial: u32,
    /// Atom name field (e.g., "CA", " N ", "OXT").
    pub name: String,
    /// Residue name (e.g., "ALA").
    pub res_name: String,
    /// Chain identifier.
    pub chain_id: char,
    /// Residue sequence number.
    pub res_seq: i32,
    /// Cartesian X coordinate (angstroms).
    pub x: f64,
    /// Cartesian Y coordinate (angstroms).
    pub y: f64,
    /// Cartesian Z coordinate (angstroms).
    pub z: f64,
    /// Occupancy factor.
    pub occupancy: f64,
    /// Temperature factor (B-factor).
    pub temp_factor: f64,
    /// Element symbol (right-justified, 2 chars in PDB; stored trimmed here).
    pub element: String,
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse all `ATOM` and `HETATM` records from a PDB-format string.
///
/// Lines that are not `ATOM` / `HETATM` records (or whose record type field
/// is not recognised) are silently skipped.
///
/// Atoms are capped at PDB_MAX_ATOMS (10,000) to prevent O(n²) bond inference DoS.
pub fn parse_pdb_atoms(input: &str) -> Vec<PdbAtom> {
    parse_pdb_atoms_capped(input, 10_000)
}

/// Resource limits for parsing untrusted PDB input.
#[derive(Debug, Clone, Copy)]
pub struct PdbParseLimits {
    /// Maximum input size in bytes.
    pub max_input_bytes: usize,
    /// Maximum physical line size in bytes.
    pub max_line_bytes: usize,
    /// Maximum ATOM/HETATM records.
    pub max_atoms: usize,
    /// Maximum MODEL records.
    pub max_models: usize,
}

impl Default for PdbParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 << 20,
            max_line_bytes: 1 << 20,
            max_atoms: 10_000,
            max_models: 1_000,
        }
    }
}

/// Error returned when a configured PDB resource limit is exceeded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbResourceLimitError {
    pub resource: &'static str,
    pub actual: usize,
    pub limit: usize,
}

impl std::fmt::Display for PdbResourceLimitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "PDB {} exceeds limit {} (got {})",
            self.resource, self.limit, self.actual
        )
    }
}

impl std::error::Error for PdbResourceLimitError {}

/// Parse PDB atoms while enforcing input, line, atom, and model limits.
pub fn parse_pdb_atoms_with_limits(
    input: &str,
    limits: &PdbParseLimits,
) -> Result<Vec<PdbAtom>, PdbResourceLimitError> {
    if input.len() > limits.max_input_bytes {
        return Err(PdbResourceLimitError {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let mut atom_count = 0usize;
    let mut model_count = 0usize;
    for line in input.lines() {
        if line.len() > limits.max_line_bytes {
            return Err(PdbResourceLimitError {
                resource: "line bytes",
                actual: line.len(),
                limit: limits.max_line_bytes,
            });
        }
        let record = line.get(0..6).unwrap_or("").trim_end();
        if record == "MODEL" {
            model_count = model_count.saturating_add(1);
            if model_count > limits.max_models {
                return Err(PdbResourceLimitError {
                    resource: "models",
                    actual: model_count,
                    limit: limits.max_models,
                });
            }
        } else if record == "ATOM" || record == "HETATM" {
            atom_count = atom_count.saturating_add(1);
            if atom_count > limits.max_atoms {
                return Err(PdbResourceLimitError {
                    resource: "atoms",
                    actual: atom_count,
                    limit: limits.max_atoms,
                });
            }
        }
    }
    Ok(parse_pdb_atoms_capped(input, limits.max_atoms))
}

fn parse_pdb_atoms_capped(input: &str, max_atoms: usize) -> Vec<PdbAtom> {
    let mut atoms = Vec::new();

    for line in input.lines() {
        if atoms.len() >= max_atoms {
            break;
        }

        let record = line.get(0..6).unwrap_or("").trim_end();
        if record != "ATOM" && record != "HETATM" {
            continue;
        }

        // Extract a fixed-width substring, tolerating short lines.
        let field = |start: usize, end: usize| -> &str {
            let end = end.min(line.len());
            line.get(start..end).unwrap_or("")
        };

        atoms.push(PdbAtom {
            serial: field(6, 11).trim().parse::<u32>().unwrap_or(0),
            name: field(12, 16).to_string(),
            res_name: field(17, 20).trim().to_string(),
            chain_id: field(21, 22).chars().next().unwrap_or(' '),
            res_seq: field(22, 26).trim().parse::<i32>().unwrap_or(0),
            x: field(30, 38).trim().parse::<f64>().unwrap_or(0.0),
            y: field(38, 46).trim().parse::<f64>().unwrap_or(0.0),
            z: field(46, 54).trim().parse::<f64>().unwrap_or(0.0),
            occupancy: field(54, 60).trim().parse::<f64>().unwrap_or(1.0),
            temp_factor: field(60, 66).trim().parse::<f64>().unwrap_or(0.0),
            element: field(76, 78).trim().to_string(),
        });
    }

    atoms
}

// ---------------------------------------------------------------------------
// PDB to Molecule
// ---------------------------------------------------------------------------

/// Build a [`Molecule`] and [`Coords3D`] from a slice of [`PdbAtom`] records.
///
/// Connectivity is inferred by distance: two atoms are bonded when their
/// distance is less than `1.3 × (r_a + r_b)` where `r_a` and `r_b` are
/// their respective covalent radii.
pub fn pdb_to_molecule(atoms: &[PdbAtom]) -> (Molecule, Coords3D) {
    let mut builder = MoleculeBuilder::new();
    let mut points: Vec<Point3> = Vec::with_capacity(atoms.len());
    let mut elements: Vec<Element> = Vec::with_capacity(atoms.len());

    for pdb_atom in atoms {
        // Use the element field when populated, else fall back to the first
        // character of the atom name.
        let sym = if !pdb_atom.element.is_empty() {
            pdb_atom.element.clone()
        } else {
            pdb_atom
                .name
                .trim()
                .chars()
                .next()
                .unwrap_or('C')
                .to_string()
        };
        let element = Element::from_symbol(&sym).unwrap_or(Element::C);
        elements.push(element);
        builder.add_atom(Atom::new(element));
        points.push(Point3::new(pdb_atom.x, pdb_atom.y, pdb_atom.z));
    }

    // Distance-based bond inference: bonded when d < 1.3 × (r_a + r_b).
    let n = points.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let threshold =
                1.3 * (elements[i].covalent_radius() as f64 + elements[j].covalent_radius() as f64);
            if points[i].distance(&points[j]) < threshold {
                // Ignore duplicate-bond errors (shouldn't occur for distinct pairs).
                let _ = builder.add_bond(AtomIdx(i as u32), AtomIdx(j as u32), BondOrder::Single);
            }
        }
    }

    (builder.build(), Coords3D { points })
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Write a molecule and its 3D coordinates as a PDB string (HETATM records).
///
/// All atoms are written as `HETATM` records in a single residue ("LIG", chain A,
/// residue 1).  No `CONECT` records are written.
///
/// Column layout (0-indexed):
///   [0..6]   "HETATM"
///   [6..11]  serial (5 chars, right-justified)
///   [11]     space
///   [12..16] atom name (4 chars, e.g. " C  ")
///   [16]     altLoc (space)
///   [17..20] resName (3 chars)
///   [20]     space
///   [21]     chainID
///   [22..26] resSeq (4 chars, right-justified)
///   [26]     iCode (space)
///   [27..30] 3 spaces
///   [30..38] x (8.3f)
///   [38..46] y (8.3f)
///   [46..54] z (8.3f)
///   [54..60] occupancy (6.2f)
///   [60..66] tempFactor (6.2f)
///   [66..76] 10 spaces
///   [76..78] element (2 chars, right-justified)
pub fn write_pdb(mol: &Molecule, coords: &Coords3D) -> String {
    let mut out = String::new();

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        let p = coords.get(idx);

        let serial = (i + 1) as u32;
        let elem_sym = atom.element.symbol();

        // PDB atom name field (4 chars): single-char elements left-padded with space.
        let atom_name = if elem_sym.len() == 1 {
            format!(" {:<3}", elem_sym)
        } else {
            format!("{:<4}", elem_sym)
        };

        // Build the line by assembling fixed-width fields directly.
        // Total width must reach at least col 78 for the element field.
        let line = format!(
            //  0-5    6-10   11  12-15         16   17-19  20  21   22-25  26  27-29   30-37     38-45     46-53     54-59     60-65       66-75         76-77
            "{:<6}{:>5} {:<4} {:>3} {:1}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}\n",
            "HETATM", // 0-5
            serial,   // 6-10
            // space at 11 (literal in format)
            atom_name, // 12-15
            // altLoc at 16 (space, literal in format)
            "LIG", // 17-19
            // space at 20 (literal in format)
            'A',   // 21
            1_i32, // 22-25
            // iCode + 3 spaces (27-29) handled by 4-space literal in format
            p.x,      // 30-37
            p.y,      // 38-45
            p.z,      // 46-53
            1.00_f64, // 54-59
            0.00_f64, // 60-65
            elem_sym, // 76-77
        );
        out.push_str(&line);
    }

    out.push_str("END\n");
    out
}

#[cfg(test)]
mod tests {
    use super::{PdbParseLimits, PdbResourceLimitError, parse_pdb_atoms_with_limits};

    const ATOM: &str =
        "ATOM      1  CA  ALA A   1      10.000  11.000  12.000  1.00 20.00           C  ";

    #[test]
    fn bounded_parser_rejects_input_lines_atoms_and_models() {
        let limits = PdbParseLimits {
            max_input_bytes: 2,
            ..PdbParseLimits::default()
        };
        assert!(matches!(
            parse_pdb_atoms_with_limits(ATOM, &limits),
            Err(PdbResourceLimitError {
                resource: "input bytes",
                ..
            })
        ));

        let limits = PdbParseLimits {
            max_line_bytes: 2,
            ..PdbParseLimits::default()
        };
        assert!(matches!(
            parse_pdb_atoms_with_limits(ATOM, &limits),
            Err(PdbResourceLimitError {
                resource: "line bytes",
                ..
            })
        ));

        let limits = PdbParseLimits {
            max_atoms: 0,
            ..PdbParseLimits::default()
        };
        assert!(matches!(
            parse_pdb_atoms_with_limits(ATOM, &limits),
            Err(PdbResourceLimitError {
                resource: "atoms",
                ..
            })
        ));

        let limits = PdbParseLimits {
            max_models: 0,
            ..PdbParseLimits::default()
        };
        assert!(matches!(
            parse_pdb_atoms_with_limits("MODEL        1\n", &limits),
            Err(PdbResourceLimitError {
                resource: "models",
                ..
            })
        ));
    }
}
