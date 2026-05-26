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
// Covalent radii (angstroms) for distance-based bonding.
// ---------------------------------------------------------------------------

fn covalent_radius(element: &Element) -> f64 {
    match element.atomic_number() {
        1  => 0.31, // H
        6  => 0.77, // C
        7  => 0.75, // N
        8  => 0.73, // O
        9  => 0.71, // F
        15 => 1.07, // P
        16 => 1.03, // S
        17 => 0.99, // Cl
        35 => 1.14, // Br
        53 => 1.33, // I
        _  => 0.77, // default
    }
}

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
pub fn parse_pdb_atoms(input: &str) -> Vec<PdbAtom> {
    let mut atoms = Vec::new();

    for line in input.lines() {
        let record = &line.get(0..6).unwrap_or("").trim_end();
        if *record != "ATOM" && *record != "HETATM" {
            continue;
        }

        // Helper: extract a fixed-width substring, padding short lines with spaces.
        let field = |start: usize, end: usize| -> &str {
            let end = end.min(line.len());
            if start >= end {
                ""
            } else {
                &line[start..end]
            }
        };

        let serial = field(6, 11).trim().parse::<u32>().unwrap_or(0);
        let name = field(12, 16).to_string();
        let res_name = field(17, 20).trim().to_string();
        let chain_id = field(21, 22).chars().next().unwrap_or(' ');
        let res_seq = field(22, 26).trim().parse::<i32>().unwrap_or(0);
        let x = field(30, 38).trim().parse::<f64>().unwrap_or(0.0);
        let y = field(38, 46).trim().parse::<f64>().unwrap_or(0.0);
        let z = field(46, 54).trim().parse::<f64>().unwrap_or(0.0);
        let occupancy = field(54, 60).trim().parse::<f64>().unwrap_or(1.0);
        let temp_factor = field(60, 66).trim().parse::<f64>().unwrap_or(0.0);
        let element = field(76, 78).trim().to_string();

        atoms.push(PdbAtom {
            serial,
            name,
            res_name,
            chain_id,
            res_seq,
            x,
            y,
            z,
            occupancy,
            temp_factor,
            element,
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
        // Try element field first, then fall back to first character of atom name.
        let sym = if !pdb_atom.element.is_empty() {
            pdb_atom.element.clone()
        } else {
            pdb_atom.name.trim().chars().next().unwrap_or('C').to_string()
        };

        let element = Element::from_symbol(&sym).unwrap_or(Element::C);
        elements.push(element);
        builder.add_atom(Atom::new(element));
        points.push(Point3::new(pdb_atom.x, pdb_atom.y, pdb_atom.z));
    }

    // Distance-based bond inference.
    let n = points.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let r_i = covalent_radius(&elements[i]);
            let r_j = covalent_radius(&elements[j]);
            let threshold = 1.3 * (r_i + r_j);
            let dist = points[i].distance(&points[j]);
            if dist < threshold {
                // Attempt to add bond; ignore duplicate-bond errors (shouldn't occur).
                let _ = builder.add_bond(AtomIdx(i as u32), AtomIdx(j as u32), BondOrder::Single);
            }
        }
    }

    let mol = builder.build();
    let coords = Coords3D { points };
    (mol, coords)
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
            "HETATM",  // 0-5
            serial,    // 6-10
            // space at 11 (literal in format)
            atom_name, // 12-15
            // altLoc at 16 (space, literal in format)
            "LIG",     // 17-19
            // space at 20 (literal in format)
            'A',       // 21
            1_i32,     // 22-25
            // iCode + 3 spaces (27-29) handled by 4-space literal in format
            p.x,       // 30-37
            p.y,       // 38-45
            p.z,       // 46-53
            1.00_f64,  // 54-59
            0.00_f64,  // 60-65
            elem_sym,  // 76-77
        );
        out.push_str(&line);
    }

    out.push_str("END\n");
    out
}
