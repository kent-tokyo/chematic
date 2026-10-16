//! Immutable Tripos MOL2 source records with a strict bounded projection.
//!
//! Supports the three mandatory sections of the Tripos MOL2 format:
//! `@<TRIPOS>MOLECULE`, `@<TRIPOS>ATOM`, and `@<TRIPOS>BOND`.
//!
//! Reference: Tripos MOL2 format specification (SYBYL 7.x).

use std::collections::{HashMap, HashSet};

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

/// Source atom metadata, aligned with the record's molecule and coordinates.
/// Partial charge is deliberately separate from [`Atom::charge`].
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct Mol2SourceAtomRecord {
    pub source_id: u32,
    pub name: String,
    pub atom_type: String,
    pub substructure_id: Option<u32>,
    pub substructure_name: Option<String>,
    /// The numeric source column, if present. Consult `charge_type()` before
    /// treating it as an assigned partial charge (e.g. `NO_CHARGES`).
    pub partial_charge: Option<f64>,
}

/// Source bond metadata, aligned with the record's molecule bonds.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct Mol2SourceBondRecord {
    pub source_id: u32,
    pub atom_ids: [u32; 2],
    /// Retains distinctions such as `am` versus `1` and `du` versus `nc`.
    pub bond_type: String,
}

/// Immutable source-preserving MOL2 record with a bounded graph projection.
///
/// [`write_mol2_source_record`] returns the original source, including names, typing,
/// optional columns, comments, and opaque sections. This is source retention,
/// not a claim that every retained construct is chemically interpreted.
/// No edit API is exposed, so opaque references cannot become stale.
///
/// The graph does not infer formal charge from partial charge. It recognizes
/// explicit `UNITY_ATOM_ATTR` charge attributes, falling back to +1 for `N.4`
/// when no explicit charge is present. Other atom-type rules and
/// non-tetrahedral coordinate-derived stereochemistry is not interpreted.
#[derive(Clone)]
pub struct Mol2SourceRecord {
    source: String,
    molecule: Molecule,
    coordinates: Vec<(f64, f64, f64)>,
    atoms: Vec<Mol2SourceAtomRecord>,
    bonds: Vec<Mol2SourceBondRecord>,
    opaque_sections: Vec<String>,
}

impl Mol2SourceRecord {
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }

    pub fn coordinates(&self) -> &[(f64, f64, f64)] {
        &self.coordinates
    }

    pub fn atoms(&self) -> &[Mol2SourceAtomRecord] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[Mol2SourceBondRecord] {
        &self.bonds
    }

    /// Sections containing source data beyond the graph projection.
    /// `UNITY_ATOM_ATTR` is listed only if it contains attributes other than charge.
    pub fn opaque_sections(&self) -> &[String] {
        &self.opaque_sections
    }

    pub fn source(&self) -> &str {
        &self.source
    }

    fn header_field(&self, index: usize) -> Option<&str> {
        self.source
            .lines()
            .skip_while(|line| !line.trim().eq_ignore_ascii_case("@<TRIPOS>MOLECULE"))
            .skip(1)
            .take_while(|line| tripos_section(line).is_none())
            .nth(index)
            .map(str::trim)
    }

    pub fn name(&self) -> Option<&str> {
        self.header_field(0)
    }

    pub fn molecule_type(&self) -> Option<&str> {
        self.header_field(2)
    }

    pub fn charge_type(&self) -> Option<&str> {
        self.header_field(3)
    }

    /// Shared inspection schema for bindings. Opaque retention and the limited
    /// graph projection are explicit; the source remains the write authority.
    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({
            "schema": "chematic.mol2-record.v1",
            "name": self.name(),
            "molecule_type": self.molecule_type(),
            "charge_type": self.charge_type(),
            "atoms": self.atoms,
            "bonds": self.bonds,
            "coordinates": self.coordinates,
            "formal_charges": self.molecule.atoms().map(|(_, atom)| atom.charge).collect::<Vec<_>>(),
            "opaque_sections": self.opaque_sections,
            "graph_projection": {
                "formal_charge_policy": "UNITY_charge_then_N.4",
                "coordinate_stereo": "noncoplanar_tetrahedral_only",
                "opaque_sections": "not_interpreted"
            },
            "source": self.source,
        })
    }
}

use crate::mol2_tripos::{Mol2Error, Mol2ParseLimits};

fn tripos_section(line: &str) -> Option<&str> {
    let line = line.trim();
    line.get(..9)
        .filter(|prefix| prefix.eq_ignore_ascii_case("@<TRIPOS>"))
        .map(|_| &line[9..])
}

/// Extract all lines belonging to a named `@<TRIPOS>SECTION` block.
///
/// Returns the (1-based line number, content) pairs starting immediately
/// after the section header, stopping at the next `@<TRIPOS>` header or EOF.
fn section_lines<'a>(lines: &'a [(usize, &'a str)], name: &str) -> Vec<(usize, &'a str)> {
    let header = format!("@<TRIPOS>{name}");
    let mut in_section = false;
    let mut result = Vec::new();
    for &(lineno, line) in lines {
        let trimmed = line.trim();
        if trimmed.eq_ignore_ascii_case(&header) {
            in_section = true;
            continue;
        }
        if in_section {
            if tripos_section(trimmed).is_some() {
                break;
            }
            // Skip blank lines and comment lines (starting with #).
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                result.push((lineno, line));
            }
        }
    }
    result
}

/// Strip a Tripos atom type suffix (e.g. `C.ar` → `C`, `N.am` → `N`, `N.3` → `N`).
fn strip_atom_type(sym: &str) -> &str {
    sym.split('.').next().unwrap_or(sym)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse a Tripos MOL2 string into a `(Molecule, 3D-coordinates)` pair.
///
/// Only the `MOLECULE`, `ATOM`, and `BOND` sections are required; all other
/// sections (SUBSTRUCTURE, SET, …) are ignored, except `UNITY_ATOM_ATTR`
/// formal-charge attributes, which are validated and projected by source ID.
///
/// 3D coordinates are returned as `Vec<(f64, f64, f64)>` aligned with the
/// molecule's atom indices.
/// Partial charges are validated but not rounded into formal charges. `N.4`
/// supplies +1 when no explicit UNITY charge exists; general formal-charge
/// perception is not performed. Use
/// [`parse_mol2_source_record`] to retain charges, atom typing, and source metadata.
#[allow(clippy::type_complexity)]
pub fn parse_mol2_source(s: &str) -> Result<(Molecule, Vec<(f64, f64, f64)>), Mol2Error> {
    parse_mol2_source_with_limits(s, &Mol2ParseLimits::default())
}

/// Parse a Tripos MOL2 string with explicit resource limits.
#[allow(clippy::type_complexity)]
pub fn parse_mol2_source_with_limits(
    s: &str,
    limits: &Mol2ParseLimits,
) -> Result<(Molecule, Vec<(f64, f64, f64)>), Mol2Error> {
    let record = parse_mol2_impl(s, limits, false)?;
    Ok((record.molecule, record.coordinates))
}

/// Parse one MOL2 record while retaining the source and typed atom/bond metadata.
/// Multi-record input and duplicate core sections are rejected rather than
/// silently taking the first molecule. Extra sections remain explicitly opaque.
pub fn parse_mol2_source_record(s: &str) -> Result<Mol2SourceRecord, Mol2Error> {
    parse_mol2_source_record_with_limits(s, &Mol2ParseLimits::default())
}

/// Source-preserving parse with the same resource limits as the graph reader.
pub fn parse_mol2_source_record_with_limits(
    s: &str,
    limits: &Mol2ParseLimits,
) -> Result<Mol2SourceRecord, Mol2Error> {
    parse_mol2_impl(s, limits, true)
}

/// Return the original record text exactly. Opaque sections are retained, not
/// interpreted. For regenerated graph-only output, use [`write_mol2_source`].
pub fn write_mol2_source_record(record: &Mol2SourceRecord) -> String {
    record.source.clone()
}

fn validate_source_sections(
    all_lines: &[(usize, &str)],
    retain_source: bool,
    unity_has_opaque: bool,
) -> Result<Vec<String>, Mol2Error> {
    let mut core_sections = HashSet::new();
    let mut opaque_sections = Vec::new();
    for &(line, text) in all_lines {
        if let Some(section) = tripos_section(text) {
            let section = section.to_ascii_uppercase();
            if matches!(section.as_str(), "MOLECULE" | "ATOM" | "BOND") {
                if !core_sections.insert(section.clone()) {
                    return Err(Mol2Error::InvalidAtomLine {
                        line,
                        detail: format!("multiple {section} sections/records are unsupported"),
                    });
                }
            } else if retain_source && (section != "UNITY_ATOM_ATTR" || unity_has_opaque) {
                opaque_sections.push(section);
            }
        }
    }
    if !core_sections.contains("MOLECULE") {
        return Err(Mol2Error::MissingSection("MOLECULE".into()));
    }

    Ok(opaque_sections)
}

fn parse_mol2_impl(
    s: &str,
    limits: &Mol2ParseLimits,
    retain_source: bool,
) -> Result<Mol2SourceRecord, Mol2Error> {
    if s.len() > limits.max_input_bytes {
        return Err(Mol2Error::ResourceLimit {
            resource: "input bytes",
            actual: s.len(),
            limit: limits.max_input_bytes,
        });
    }
    let lines = s.lines().collect::<Vec<_>>();
    if lines.len() > limits.max_lines {
        return Err(Mol2Error::ResourceLimit {
            resource: "lines",
            actual: lines.len(),
            limit: limits.max_lines,
        });
    }
    if let Some(line_bytes) = lines.iter().map(|line| line.len()).max()
        && line_bytes > limits.max_line_bytes
    {
        return Err(Mol2Error::ResourceLimit {
            resource: "line bytes",
            actual: line_bytes,
            limit: limits.max_line_bytes,
        });
    }
    let section_count = lines
        .iter()
        .filter(|line| tripos_section(line).is_some())
        .count();
    if section_count > limits.max_sections {
        return Err(Mol2Error::ResourceLimit {
            resource: "sections",
            actual: section_count,
            limit: limits.max_sections,
        });
    }
    let all_lines: Vec<(usize, &str)> = lines
        .into_iter()
        .enumerate()
        .map(|(i, l)| (i + 1, l))
        .collect();
    let (formal_charges, unity_attributes) =
        crate::mol2_tripos::parse_unity_formal_charges(&all_lines)?;
    let unity_has_opaque = unity_attributes.iter().any(|group| {
        group
            .attributes
            .iter()
            .any(|(name, _)| !name.eq_ignore_ascii_case("charge"))
    });

    let opaque_sections = validate_source_sections(&all_lines, retain_source, unity_has_opaque)?;

    // -- ATOM section ----------------------------------------------------------
    let atom_lines = section_lines(&all_lines, "ATOM");
    if atom_lines.is_empty() {
        return Err(Mol2Error::MissingSection("ATOM".into()));
    }

    let header_index = all_lines
        .iter()
        .position(|(_, text)| text.trim().eq_ignore_ascii_case("@<TRIPOS>MOLECULE"))
        .expect("MOLECULE section checked above");
    let counts_line = header_index + 2;
    let counts: Vec<&str> = all_lines
        .get(counts_line)
        .map(|(_, text)| text.split_whitespace().collect())
        .unwrap_or_default();
    let declared_atoms = counts.first().and_then(|value| value.parse::<usize>().ok());
    let declared_bonds = counts.get(1).and_then(|value| value.parse::<usize>().ok());
    if declared_atoms != Some(atom_lines.len()) {
        return Err(Mol2Error::InvalidAtomLine {
            line: counts_line + 1,
            detail: format!(
                "invalid or mismatched atom count; found {} atom records",
                atom_lines.len()
            ),
        });
    }

    let mut builder = MoleculeBuilder::new();
    let mut coords: Vec<(f64, f64, f64)> = Vec::new();
    // Map from MOL2 1-based atom_id → builder AtomIdx.
    let mut atom_id_map: HashMap<u32, AtomIdx> = HashMap::new();
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();

    if atom_lines.len() > limits.max_atoms {
        return Err(Mol2Error::ResourceLimit {
            resource: "atom records",
            actual: atom_lines.len(),
            limit: limits.max_atoms,
        });
    }

    for (lineno, line) in &atom_lines {
        // MOL2 ATOM line format (space-separated):
        // atom_id  atom_name  x  y  z  atom_type  [subst_id  subst_name  charge]
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 6 {
            return Err(Mol2Error::InvalidAtomLine {
                line: *lineno,
                detail: format!("expected at least 6 fields, got {}", parts.len()),
            });
        }

        let atom_id: u32 = parts[0].parse().map_err(|_| Mol2Error::InvalidAtomLine {
            line: *lineno,
            detail: format!("cannot parse atom_id from '{}'", parts[0]),
        })?;
        if atom_id == 0 || atom_id_map.contains_key(&atom_id) {
            return Err(Mol2Error::InvalidAtomLine {
                line: *lineno,
                detail: format!("atom_id {atom_id} must be positive and unique"),
            });
        }

        let x: f64 = parts[2].parse().map_err(|_| Mol2Error::InvalidAtomLine {
            line: *lineno,
            detail: format!("cannot parse x from '{}'", parts[2]),
        })?;
        let y: f64 = parts[3].parse().map_err(|_| Mol2Error::InvalidAtomLine {
            line: *lineno,
            detail: format!("cannot parse y from '{}'", parts[3]),
        })?;
        let z: f64 = parts[4].parse().map_err(|_| Mol2Error::InvalidAtomLine {
            line: *lineno,
            detail: format!("cannot parse z from '{}'", parts[4]),
        })?;
        if !x.is_finite() || !y.is_finite() || !z.is_finite() {
            return Err(Mol2Error::InvalidAtomLine {
                line: *lineno,
                detail: "coordinates must be finite".into(),
            });
        }

        // atom_type may be "C.ar", "N.3", "O.2", etc. Strip the suffix.
        let atom_type_raw = parts[5];
        let sym_raw = strip_atom_type(atom_type_raw);
        // Capitalise first letter for element lookup.
        let sym = {
            let mut c = sym_raw.chars();
            match c.next() {
                None => String::new(),
                Some(first) => first.to_uppercase().collect::<String>() + c.as_str(),
            }
        };

        // Optional partial charge (field 9, 0-indexed as parts[8]).
        let partial_charge = parts
            .get(8)
            .map(|raw| {
                raw.parse::<f64>()
                    .ok()
                    .filter(|charge| charge.is_finite())
                    .ok_or_else(|| Mol2Error::InvalidAtomLine {
                        line: *lineno,
                        detail: format!("partial charge '{raw}' must be a finite number"),
                    })
            })
            .transpose()?;
        let substructure_id = parts
            .get(6)
            .map(|raw| {
                raw.parse::<u32>().map_err(|_| Mol2Error::InvalidAtomLine {
                    line: *lineno,
                    detail: format!("invalid substructure id '{raw}'"),
                })
            })
            .transpose()?;

        let element = Element::from_symbol(&sym).ok_or_else(|| Mol2Error::UnknownElement {
            symbol: sym.clone(),
            line: *lineno,
        })?;

        let mut atom = Atom::new(element);
        // An electrostatic partial charge is not a formal charge. The N.4
        // atom type explicitly identifies positively charged nitrogen.
        atom.charge = formal_charges
            .get(&atom_id)
            .copied()
            .unwrap_or_else(|| i8::from(atom_type_raw == "N.4"));

        let builder_idx = builder.add_atom(atom);
        atom_id_map.insert(atom_id, builder_idx);
        coords.push((x, y, z));
        if retain_source {
            atoms.push(Mol2SourceAtomRecord {
                source_id: atom_id,
                name: parts[1].to_string(),
                atom_type: atom_type_raw.to_string(),
                substructure_id,
                substructure_name: parts.get(7).map(|value| value.to_string()),
                partial_charge,
            });
        }
    }

    // -- BOND section ----------------------------------------------------------
    let bond_lines = section_lines(&all_lines, "BOND");
    // BOND section is optional when the header declares zero bonds.
    if declared_bonds != Some(bond_lines.len()) {
        return Err(Mol2Error::InvalidBondLine {
            line: counts_line + 1,
            detail: format!(
                "invalid or mismatched bond count; found {} bond records",
                bond_lines.len()
            ),
        });
    }

    if bond_lines.len() > limits.max_bonds {
        return Err(Mol2Error::ResourceLimit {
            resource: "bond records",
            actual: bond_lines.len(),
            limit: limits.max_bonds,
        });
    }
    let mut bond_ids = HashSet::new();
    for (lineno, line) in &bond_lines {
        // bond_id  origin_atom_id  target_atom_id  bond_type
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: format!("expected at least 4 fields, got {}", parts.len()),
            });
        }
        let bond_id: u32 = parts[0].parse().map_err(|_| Mol2Error::InvalidBondLine {
            line: *lineno,
            detail: format!("invalid bond_id '{}'", parts[0]),
        })?;
        if bond_id == 0 || !bond_ids.insert(bond_id) {
            return Err(Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: format!("bond_id {bond_id} must be positive and unique"),
            });
        }

        let a1_id: u32 = parts[1].parse().map_err(|_| Mol2Error::InvalidBondLine {
            line: *lineno,
            detail: format!("cannot parse origin atom_id from '{}'", parts[1]),
        })?;
        let a2_id: u32 = parts[2].parse().map_err(|_| Mol2Error::InvalidBondLine {
            line: *lineno,
            detail: format!("cannot parse target atom_id from '{}'", parts[2]),
        })?;

        let a1 = atom_id_map
            .get(&a1_id)
            .copied()
            .ok_or_else(|| Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: format!("atom_id {a1_id} not found"),
            })?;
        let a2 = atom_id_map
            .get(&a2_id)
            .copied()
            .ok_or_else(|| Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: format!("atom_id {a2_id} not found"),
            })?;

        let bond_type = parts[3];
        let order = match bond_type {
            "1" | "am" => BondOrder::Single,
            "2" => BondOrder::Double,
            "3" => BondOrder::Triple,
            "ar" => BondOrder::Aromatic,
            "un" => BondOrder::QueryAny,
            "du" | "nc" => BondOrder::Zero,
            _ => {
                return Err(Mol2Error::InvalidBondLine {
                    line: *lineno,
                    detail: format!("unsupported bond type '{bond_type}'"),
                });
            }
        };

        if a1 == a2 {
            return Err(Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: "self bonds are unsupported".into(),
            });
        }
        builder
            .add_bond(a1, a2, order)
            .map_err(|error| Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: error.to_string(),
            })?;
        if retain_source {
            bonds.push(Mol2SourceBondRecord {
                source_id: bond_id,
                atom_ids: [a1_id, a2_id],
                bond_type: bond_type.to_string(),
            });
        }
    }

    let mut molecule = builder.build();
    crate::mol2_tripos::perceive_mol2_stereo_from_coordinates(&mut molecule, &coords);
    Ok(Mol2SourceRecord {
        source: if retain_source {
            s.to_string()
        } else {
            String::new()
        },
        molecule,
        coordinates: coords,
        atoms,
        bonds,
        opaque_sections,
    })
}

/// Write a molecule and its 3D coordinates to Tripos MOL2 format.
/// This graph-only writer does not compute partial charges or preserve source
/// names/types. It emits `NO_CHARGES`. Use [`write_mol2_source_record`] for an exact
/// source-preserving round trip, including the original charge method/values.
/// Nonzero formal charges are emitted as `UNITY_ATOM_ATTR` charge attributes;
/// `NO_CHARGES` refers to the partial-charge column.
///
/// `coords` must have one entry per atom in `mol`.  Missing entries use
/// `(0.0, 0.0, 0.0)`.
pub fn write_mol2_source(mol: &Molecule, coords: &[(f64, f64, f64)]) -> String {
    let mut out = String::new();

    // MOLECULE section.
    out.push_str("@<TRIPOS>MOLECULE\n");
    out.push_str("chematic\n");
    out.push_str(&format!(
        "{} {} 0 0 0\n",
        mol.atom_count(),
        mol.bond_count()
    ));
    out.push_str("SMALL\n");
    out.push_str("NO_CHARGES\n\n");

    // ATOM section.
    out.push_str("@<TRIPOS>ATOM\n");
    for (idx, atom) in mol.atoms() {
        let i = idx.0 + 1; // 1-based
        let sym = atom.element.symbol();
        let (x, y, z) = coords
            .get(idx.0 as usize)
            .copied()
            .unwrap_or((0.0, 0.0, 0.0));
        // Use symbol as atom type (simplified; no hybridisation-based suffix).
        out.push_str(&format!(
            "{i:>6} {sym:<4} {x:>10.4} {y:>10.4} {z:>10.4} {sym:<8} 1  LIG  0.0000\n"
        ));
    }

    // BOND section.
    out.push_str("@<TRIPOS>BOND\n");
    for (bidx, bond) in mol.bonds() {
        let bi = bidx.0 + 1;
        let a1 = bond.atom1.0 + 1;
        let a2 = bond.atom2.0 + 1;
        let btype = match bond.order {
            BondOrder::Zero => "nc",
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => "1",
            BondOrder::Double => "2",
            BondOrder::Triple => "3",
            BondOrder::Aromatic => "ar",
            BondOrder::Quadruple => "4",
            BondOrder::QueryAny
            | BondOrder::QuerySingleOrDouble
            | BondOrder::QuerySingleOrAromatic
            | BondOrder::QueryDoubleOrAromatic => "un",
        };
        out.push_str(&format!("{bi:>6} {a1:>6} {a2:>6} {btype}\n"));
    }

    if mol.atoms().any(|(_, atom)| atom.charge != 0) {
        out.push_str("@<TRIPOS>UNITY_ATOM_ATTR\n");
        for (idx, atom) in mol.atoms().filter(|(_, atom)| atom.charge != 0) {
            out.push_str(&format!("{} 1\ncharge {}\n", idx.0 + 1, atom.charge));
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

    const ETHANOL_MOL2: &str = "\
@<TRIPOS>MOLECULE
ethanol
 3 2 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
      2 C2          1.5000    0.0000    0.0000 C.3     1  LIG1        0.0000
      3 O1          3.0000    0.0000    0.0000 O.3     1  LIG1       -0.3940
@<TRIPOS>BOND
     1     1     2    1
     2     2     3    1
";

    #[test]
    fn test_parse_mol2_atom_count() {
        let (mol, coords) = parse_mol2_source(ETHANOL_MOL2).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        assert_eq!(coords.len(), 3);
    }

    #[test]
    fn test_parse_mol2_coords() {
        let (_, coords) = parse_mol2_source(ETHANOL_MOL2).unwrap();
        assert!((coords[0].0 - 0.0).abs() < 1e-4);
        assert!((coords[1].0 - 1.5).abs() < 1e-4);
        assert!((coords[2].0 - 3.0).abs() < 1e-4);
    }

    #[test]
    fn test_parse_mol2_elements() {
        let (mol, _) = parse_mol2_source(ETHANOL_MOL2).unwrap();
        assert_eq!(mol.atom(AtomIdx(0)).element.symbol(), "C");
        assert_eq!(mol.atom(AtomIdx(1)).element.symbol(), "C");
        assert_eq!(mol.atom(AtomIdx(2)).element.symbol(), "O");
    }

    #[test]
    fn partial_charges_are_not_formal_charges() {
        for charge in ["-0.7000", "0.7000", "-1.4000", "1.4000"] {
            let input = ETHANOL_MOL2.replace("-0.3940", charge);
            let (mol, _) = parse_mol2_source(&input).unwrap();
            assert_eq!(mol.atom(AtomIdx(2)).charge, 0, "partial charge {charge}");
        }
    }

    #[test]
    fn amide_bonds_are_single_not_aromatic() {
        let input = ETHANOL_MOL2.replace("2     2     3    1", "2     2     3    am");
        let (mol, _) = parse_mol2_source(&input).unwrap();
        assert_eq!(mol.bonds().nth(1).unwrap().1.order, BondOrder::Single);
    }

    #[test]
    fn unknown_bond_type_does_not_silently_become_single() {
        let input = ETHANOL_MOL2.replace("2     2     3    1", "2     2     3    vendor");
        assert!(matches!(
            parse_mol2_source(&input),
            Err(Mol2Error::InvalidBondLine { .. })
        ));
    }

    #[test]
    fn graph_writer_does_not_claim_computed_partial_charges() {
        let (mol, coords) = parse_mol2_source(ETHANOL_MOL2).unwrap();
        let output = write_mol2_source(&mol, &coords);
        assert!(output.contains("\nNO_CHARGES\n"));
        assert!(!output.contains("GASTEIGER"));
    }

    #[test]
    fn test_write_mol2_roundtrip() {
        let (mol, coords) = parse_mol2_source(ETHANOL_MOL2).unwrap();
        let written = write_mol2_source(&mol, &coords);
        let (mol2, coords2) = parse_mol2_source(&written).unwrap();
        assert_eq!(mol.atom_count(), mol2.atom_count());
        assert_eq!(mol.bond_count(), mol2.bond_count());
        // Coordinates should survive round-trip.
        for (a, b) in coords.iter().zip(coords2.iter()) {
            assert!((a.0 - b.0).abs() < 0.01);
            assert!((a.1 - b.1).abs() < 0.01);
            assert!((a.2 - b.2).abs() < 0.01);
        }
    }

    #[test]
    fn shared_mol2_roundtrip_contract_matches() {
        let document: serde_json::Value = serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/cross_binding_contract.json"
        )))
        .expect("contract JSON");
        let contract = &document["mol2_contract"];
        let (mol, coords) = parse_mol2_source(contract["input"].as_str().unwrap()).unwrap();
        let serialized = write_mol2_source(&mol, &coords);
        let (roundtripped, _) = parse_mol2_source(&serialized).unwrap();
        assert_eq!(
            roundtripped.atom_count(),
            contract["expected"]["atom_count"]
        );
        assert_eq!(
            roundtripped.bond_count(),
            contract["expected"]["bond_count"]
        );
    }

    #[test]
    fn test_aromatic_bond_type() {
        let benzene_mol2 = "\
@<TRIPOS>MOLECULE
benzene
 6 6 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C1  0.0000  1.4000  0.0000 C.ar  1  BNZ  0.0000
      2 C2  1.2124  0.7000  0.0000 C.ar  1  BNZ  0.0000
      3 C3  1.2124 -0.7000  0.0000 C.ar  1  BNZ  0.0000
      4 C4  0.0000 -1.4000  0.0000 C.ar  1  BNZ  0.0000
      5 C5 -1.2124 -0.7000  0.0000 C.ar  1  BNZ  0.0000
      6 C6 -1.2124  0.7000  0.0000 C.ar  1  BNZ  0.0000
@<TRIPOS>BOND
     1     1     2   ar
     2     2     3   ar
     3     3     4   ar
     4     4     5   ar
     5     5     6   ar
     6     6     1   ar
";
        let (mol, _) = parse_mol2_source(benzene_mol2).unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        // All bonds should be aromatic.
        for (_, bond) in mol.bonds() {
            assert_eq!(bond.order, BondOrder::Aromatic);
        }
    }

    #[test]
    fn test_missing_atom_section() {
        let bad = "@<TRIPOS>MOLECULE\nbad\n";
        assert!(parse_mol2_source(bad).is_err());
    }

    #[test]
    fn bounded_parser_rejects_input_and_line_limits() {
        assert!(matches!(
            parse_mol2_source_with_limits(
                ETHANOL_MOL2,
                &Mol2ParseLimits {
                    max_input_bytes: 8,
                    ..Default::default()
                }
            ),
            Err(Mol2Error::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));
        let long_line = format!("{}\n", "x".repeat(32));
        assert!(matches!(
            parse_mol2_source_with_limits(
                &long_line,
                &Mol2ParseLimits {
                    max_line_bytes: 16,
                    ..Default::default()
                }
            ),
            Err(Mol2Error::ResourceLimit {
                resource: "line bytes",
                ..
            })
        ));
    }

    #[test]
    fn bounded_parser_rejects_atom_and_bond_limits() {
        assert!(matches!(
            parse_mol2_source_with_limits(
                ETHANOL_MOL2,
                &Mol2ParseLimits {
                    max_atoms: 2,
                    ..Default::default()
                }
            ),
            Err(Mol2Error::ResourceLimit {
                resource: "atom records",
                ..
            })
        ));
        assert!(matches!(
            parse_mol2_source_with_limits(
                ETHANOL_MOL2,
                &Mol2ParseLimits {
                    max_bonds: 1,
                    ..Default::default()
                }
            ),
            Err(Mol2Error::ResourceLimit {
                resource: "bond records",
                ..
            })
        ));
    }
}
