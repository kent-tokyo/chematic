//! Tripos MOL2 format reader and writer.
//!
//! Supports the three mandatory sections of the Tripos MOL2 format:
//! `@<TRIPOS>MOLECULE`, `@<TRIPOS>ATOM`, and `@<TRIPOS>BOND`.
//!
//! Reference: Tripos MOL2 format specification (SYBYL 7.x).

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

/// One `@<TRIPOS>ATOM` row whose interchange fields are not representable by
/// [`Molecule`]. The vector is index-aligned with [`Mol2Record::molecule`].
#[derive(Debug, Clone, PartialEq)]
pub struct Mol2AtomRecord {
    pub atom_id: u32,
    pub atom_name: String,
    pub atom_type: String,
    pub subst_id: Option<i32>,
    pub subst_name: Option<String>,
    pub partial_charge: Option<f64>,
    pub status_bits: Vec<String>,
}

/// One `@<TRIPOS>BOND` row whose source identifier and status bits are not
/// represented by [`Molecule`]. The vector is index-aligned with its bonds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mol2BondRecord {
    pub bond_id: u32,
    pub atom1_id: u32,
    pub atom2_id: u32,
    pub bond_type: String,
    pub status_bits: Vec<String>,
}

/// An untyped Tripos section retained verbatim across a checked MOL2
/// parse/write round trip.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mol2OpaqueSection {
    pub name: String,
    pub lines: Vec<String>,
}

/// One `UNITY_ATOM_ATTR` group. Unknown attributes remain explicit instead of
/// being dropped while the recognized `charge` attribute feeds `Atom.charge`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mol2UnityAtomAttributes {
    pub atom_id: u32,
    pub attributes: Vec<(String, String)>,
}

/// Loss-aware Tripos MOL2 record.
///
/// The ordinary [`parse_mol2`] API remains available for graph-only callers.
/// Interchange code should use this type so atom types, partial charges,
/// residue labels, status bits, and extension sections are not silently lost.
#[derive(Clone)]
pub struct Mol2Record {
    pub molecule: Molecule,
    pub coords: Vec<(f64, f64, f64)>,
    pub name: String,
    pub molecule_type: String,
    pub charge_type: String,
    pub molecule_tail: Vec<String>,
    pub atoms: Vec<Mol2AtomRecord>,
    pub bonds: Vec<Mol2BondRecord>,
    pub duplicate_bonds: Vec<Mol2BondRecord>,
    pub unity_atom_attributes: Vec<Mol2UnityAtomAttributes>,
    pub opaque_sections: Vec<Mol2OpaqueSection>,
}

/// Resource limits for Tripos MOL2 parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mol2ParseLimits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_sections: usize,
    pub max_atoms: usize,
    pub max_bonds: usize,
}

impl Default for Mol2ParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_line_bytes: 1024 * 1024,
            max_lines: 1_000_000,
            max_sections: 256,
            max_atoms: 1_000_000,
            max_bonds: 2_000_000,
        }
    }
}

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Error returned when parsing a Tripos MOL2 file fails.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Mol2Error {
    /// Required section not found.
    MissingSection(String),
    /// This one-record API received a multi-molecule MOL2 stream.
    MultipleMolecules { count: usize },
    /// An atom line is malformed.
    InvalidAtomLine { line: usize, detail: String },
    /// A bond line is malformed.
    InvalidBondLine { line: usize, detail: String },
    /// Unknown element symbol.
    UnknownElement { symbol: String, line: usize },
    /// The input exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl core::fmt::Display for Mol2Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MissingSection(s) => write!(f, "MOL2: missing section @<TRIPOS>{s}"),
            Self::MultipleMolecules { count } => write!(
                f,
                "MOL2: this API accepts one molecule, but the input contains {count} MOLECULE sections"
            ),
            Self::InvalidAtomLine { line, detail } => {
                write!(f, "MOL2: invalid atom line {line}: {detail}")
            }
            Self::InvalidBondLine { line, detail } => {
                write!(f, "MOL2: invalid bond line {line}: {detail}")
            }
            Self::UnknownElement { symbol, line } => {
                write!(f, "MOL2: unknown element '{symbol}' at line {line}")
            }
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "MOL2: {resource} has size {actual}, exceeding limit {limit}"
            ),
        }
    }
}

impl std::error::Error for Mol2Error {}

// ---------------------------------------------------------------------------
// Parser helpers
// ---------------------------------------------------------------------------

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
            if trimmed.starts_with("@<TRIPOS>") {
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

fn parse_unity_formal_charges(
    lines: &[(usize, &str)],
) -> Result<
    (
        std::collections::HashMap<u32, i8>,
        Vec<Mol2UnityAtomAttributes>,
    ),
    Mol2Error,
> {
    let rows = section_lines(lines, "UNITY_ATOM_ATTR");
    let mut charges = std::collections::HashMap::new();
    let mut groups = Vec::new();
    let mut cursor = 0usize;
    while cursor < rows.len() {
        let (line_no, header) = rows[cursor];
        let fields = header.split_whitespace().collect::<Vec<_>>();
        if fields.len() != 2 {
            return Err(Mol2Error::InvalidAtomLine {
                line: line_no,
                detail: "UNITY_ATOM_ATTR header must be '<atom_id> <attribute_count>'".into(),
            });
        }
        let atom_id = fields[0]
            .parse::<u32>()
            .map_err(|_| Mol2Error::InvalidAtomLine {
                line: line_no,
                detail: format!("cannot parse UNITY atom_id from '{}'", fields[0]),
            })?;
        let count = fields[1]
            .parse::<usize>()
            .map_err(|_| Mol2Error::InvalidAtomLine {
                line: line_no,
                detail: format!("cannot parse UNITY attribute count from '{}'", fields[1]),
            })?;
        cursor += 1;
        if cursor.saturating_add(count) > rows.len() {
            return Err(Mol2Error::InvalidAtomLine {
                line: line_no,
                detail: "UNITY_ATOM_ATTR ends before all declared attributes".into(),
            });
        }
        let mut attributes = Vec::with_capacity(count);
        for &(attribute_line, attribute) in &rows[cursor..cursor + count] {
            let mut parts = attribute.split_whitespace();
            let Some(name) = parts.next() else { continue };
            let Some(value) = parts.next() else { continue };
            attributes.push((name.to_string(), value.to_string()));
            if name.eq_ignore_ascii_case("charge") {
                let charge = value
                    .parse::<i8>()
                    .map_err(|_| Mol2Error::InvalidAtomLine {
                        line: attribute_line,
                        detail: format!("cannot parse formal charge from '{value}'"),
                    })?;
                charges.insert(atom_id, charge);
            }
        }
        groups.push(Mol2UnityAtomAttributes {
            atom_id,
            attributes,
        });
        cursor += count;
    }
    Ok((charges, groups))
}

fn opaque_sections(lines: &[(usize, &str)]) -> Vec<Mol2OpaqueSection> {
    let mut sections = Vec::new();
    let mut current: Option<Mol2OpaqueSection> = None;
    for &(_, line) in lines {
        let trimmed = line.trim();
        if let Some(name) = trimmed.strip_prefix("@<TRIPOS>") {
            if let Some(section) = current.take() {
                sections.push(section);
            }
            current =
                (!matches!(name, "MOLECULE" | "ATOM" | "BOND" | "UNITY_ATOM_ATTR")).then(|| {
                    Mol2OpaqueSection {
                        name: name.to_string(),
                        lines: Vec::new(),
                    }
                });
        } else if let Some(section) = current.as_mut() {
            section.lines.push(line.to_string());
        }
    }
    if let Some(section) = current {
        sections.push(section);
    }
    sections
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse a Tripos MOL2 string into a `(Molecule, 3D-coordinates)` pair.
///
/// Only the `MOLECULE`, `ATOM`, and `BOND` sections are required; all other
/// sections (SUBSTRUCTURE, SET, …) are silently ignored.
///
/// 3D coordinates are returned as `Vec<(f64, f64, f64)>` aligned with the
/// molecule's atom indices.
#[allow(clippy::type_complexity)]
pub fn parse_mol2(s: &str) -> Result<(Molecule, Vec<(f64, f64, f64)>), Mol2Error> {
    parse_mol2_record(s).map(|record| (record.molecule, record.coords))
}

/// Parse a Tripos MOL2 string with explicit resource limits.
#[allow(clippy::type_complexity)]
pub fn parse_mol2_with_limits(
    s: &str,
    limits: &Mol2ParseLimits,
) -> Result<(Molecule, Vec<(f64, f64, f64)>), Mol2Error> {
    parse_mol2_record_with_limits(s, limits).map(|record| (record.molecule, record.coords))
}

/// Parse a MOL2 string without discarding interchange-only fields.
pub fn parse_mol2_record(s: &str) -> Result<Mol2Record, Mol2Error> {
    parse_mol2_record_with_limits(s, &Mol2ParseLimits::default())
}

/// Parse a loss-aware MOL2 record with explicit resource limits.
pub fn parse_mol2_record_with_limits(
    s: &str,
    limits: &Mol2ParseLimits,
) -> Result<Mol2Record, Mol2Error> {
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
        .filter(|line| line.trim_start().starts_with("@<TRIPOS>"))
        .count();
    if section_count > limits.max_sections {
        return Err(Mol2Error::ResourceLimit {
            resource: "sections",
            actual: section_count,
            limit: limits.max_sections,
        });
    }
    let molecule_count = lines
        .iter()
        .filter(|line| line.trim().eq_ignore_ascii_case("@<TRIPOS>MOLECULE"))
        .count();
    if molecule_count > 1 {
        return Err(Mol2Error::MultipleMolecules {
            count: molecule_count,
        });
    }
    let all_lines: Vec<(usize, &str)> = lines
        .into_iter()
        .enumerate()
        .map(|(i, l)| (i + 1, l))
        .collect();

    let molecule_lines = section_lines(&all_lines, "MOLECULE");
    if molecule_lines.is_empty() {
        return Err(Mol2Error::MissingSection("MOLECULE".into()));
    }
    let name = molecule_lines
        .first()
        .map(|(_, line)| line.trim().to_string())
        .unwrap_or_default();
    let molecule_type = molecule_lines
        .get(2)
        .map(|(_, line)| line.trim().to_string())
        .unwrap_or_else(|| "SMALL".to_string());
    let charge_type = molecule_lines
        .get(3)
        .map(|(_, line)| line.trim().to_string())
        .unwrap_or_else(|| "NO_CHARGES".to_string());
    let molecule_tail = molecule_lines
        .iter()
        .skip(4)
        .map(|(_, line)| (*line).to_string())
        .collect::<Vec<_>>();
    let (unity_formal_charges, unity_atom_attributes) = parse_unity_formal_charges(&all_lines)?;

    // -- ATOM section ----------------------------------------------------------
    let atom_lines = section_lines(&all_lines, "ATOM");
    if atom_lines.is_empty() {
        return Err(Mol2Error::MissingSection("ATOM".into()));
    }

    let mut builder = MoleculeBuilder::new();
    let mut coords: Vec<(f64, f64, f64)> = Vec::new();
    let mut atom_records = Vec::new();
    // Map from MOL2 1-based atom_id → builder AtomIdx.
    let mut atom_id_map: Vec<(u32, AtomIdx)> = Vec::new();

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

        let partial_charge = parts
            .get(8)
            .map(|value| {
                value
                    .parse::<f64>()
                    .map_err(|_| Mol2Error::InvalidAtomLine {
                        line: *lineno,
                        detail: format!("cannot parse partial charge from '{value}'"),
                    })
            })
            .transpose()?;

        let element = Element::from_symbol(&sym).ok_or_else(|| Mol2Error::UnknownElement {
            symbol: sym.clone(),
            line: *lineno,
        })?;

        let mut atom = Atom::new(element);
        atom.charge = unity_formal_charges.get(&atom_id).copied().unwrap_or(0);
        atom.aromatic = atom_type_raw
            .split_once('.')
            .is_some_and(|(_, suffix)| suffix.eq_ignore_ascii_case("ar"));

        let builder_idx = builder.add_atom(atom);
        atom_id_map.push((atom_id, builder_idx));
        coords.push((x, y, z));
        atom_records.push(Mol2AtomRecord {
            atom_id,
            atom_name: parts[1].to_string(),
            atom_type: atom_type_raw.to_string(),
            subst_id: parts.get(6).and_then(|value| value.parse().ok()),
            subst_name: parts.get(7).map(|value| (*value).to_string()),
            partial_charge,
            status_bits: parts
                .iter()
                .skip(9)
                .map(|value| (*value).to_string())
                .collect(),
        });
    }

    // -- BOND section ----------------------------------------------------------
    let bond_lines = section_lines(&all_lines, "BOND");
    // BOND section is optional (0-atom molecules are valid).

    if bond_lines.len() > limits.max_bonds {
        return Err(Mol2Error::ResourceLimit {
            resource: "bond records",
            actual: bond_lines.len(),
            limit: limits.max_bonds,
        });
    }
    let mut bond_records = Vec::new();
    let mut duplicate_bonds = Vec::new();
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
            detail: format!("cannot parse bond_id from '{}'", parts[0]),
        })?;
        let a1_id: u32 = parts[1].parse().map_err(|_| Mol2Error::InvalidBondLine {
            line: *lineno,
            detail: format!("cannot parse origin atom_id from '{}'", parts[1]),
        })?;
        let a2_id: u32 = parts[2].parse().map_err(|_| Mol2Error::InvalidBondLine {
            line: *lineno,
            detail: format!("cannot parse target atom_id from '{}'", parts[2]),
        })?;

        let a1 = atom_id_map
            .iter()
            .find(|&&(k, _)| k == a1_id)
            .map(|&(_, v)| v)
            .ok_or_else(|| Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: format!("atom_id {a1_id} not found"),
            })?;
        let a2 = atom_id_map
            .iter()
            .find(|&&(k, _)| k == a2_id)
            .map(|&(_, v)| v)
            .ok_or_else(|| Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: format!("atom_id {a2_id} not found"),
            })?;

        let bond_type = parts[3];
        let order = match bond_type {
            "1" | "1.5" => BondOrder::Single,
            "2" => BondOrder::Double,
            "3" => BondOrder::Triple,
            "4" => BondOrder::Quadruple,
            "ar" | "am" => BondOrder::Aromatic,
            "un" => BondOrder::QueryAny,
            "du" | "nc" => BondOrder::Zero,
            _ => BondOrder::Single,
        };

        let source = Mol2BondRecord {
            bond_id,
            atom1_id: a1_id,
            atom2_id: a2_id,
            bond_type: bond_type.to_string(),
            status_bits: parts
                .iter()
                .skip(4)
                .map(|value| (*value).to_string())
                .collect(),
        };
        // Keep the graph API's historical tolerance while retaining repeated
        // rows for a loss-aware same-format rewrite.
        if builder.add_bond(a1, a2, order).is_ok() {
            bond_records.push(source);
        } else {
            duplicate_bonds.push(source);
        }
    }

    Ok(Mol2Record {
        molecule: builder.build(),
        coords,
        name,
        molecule_type,
        charge_type,
        molecule_tail,
        atoms: atom_records,
        bonds: bond_records,
        duplicate_bonds,
        unity_atom_attributes,
        opaque_sections: opaque_sections(&all_lines),
    })
}

/// Write a molecule and its 3D coordinates to Tripos MOL2 format.
///
/// `coords` must have one entry per atom in `mol`.  Missing entries use
/// `(0.0, 0.0, 0.0)`.
pub fn write_mol2(mol: &Molecule, coords: &[(f64, f64, f64)]) -> String {
    let record = Mol2Record {
        molecule: mol.clone(),
        coords: coords.to_vec(),
        name: "chematic".into(),
        molecule_type: "SMALL".into(),
        charge_type: "NO_CHARGES".into(),
        molecule_tail: Vec::new(),
        atoms: mol
            .atoms()
            .map(|(idx, atom)| Mol2AtomRecord {
                atom_id: idx.0 + 1,
                atom_name: format!("{}{}", atom.element.symbol(), idx.0 + 1),
                atom_type: atom.element.symbol().to_string(),
                subst_id: Some(1),
                subst_name: Some("LIG".into()),
                partial_charge: Some(0.0),
                status_bits: Vec::new(),
            })
            .collect(),
        bonds: mol
            .bonds()
            .map(|(idx, bond)| Mol2BondRecord {
                bond_id: idx.0 + 1,
                atom1_id: bond.atom1.0 + 1,
                atom2_id: bond.atom2.0 + 1,
                bond_type: bond_type_token(bond.order).to_string(),
                status_bits: Vec::new(),
            })
            .collect(),
        duplicate_bonds: Vec::new(),
        unity_atom_attributes: Vec::new(),
        opaque_sections: Vec::new(),
    };
    write_mol2_record(&record)
}

fn bond_type_token(order: BondOrder) -> &'static str {
    match order {
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
    }
}

/// Serialize a loss-aware MOL2 record, retaining atom types, partial charges,
/// residue fields, status bits, formal-charge attributes, and opaque sections.
pub fn write_mol2_record(record: &Mol2Record) -> String {
    let mol = &record.molecule;
    let coords = &record.coords;
    let mut out = String::new();

    // MOLECULE section.
    out.push_str("@<TRIPOS>MOLECULE\n");
    out.push_str(&record.name);
    out.push('\n');
    let substructures = record
        .opaque_sections
        .iter()
        .find(|section| section.name.eq_ignore_ascii_case("SUBSTRUCTURE"))
        .map_or(0, |section| {
            section
                .lines
                .iter()
                .filter(|line| !line.trim().is_empty())
                .count()
        });
    out.push_str(&format!(
        "{} {} {substructures} 0 0\n",
        mol.atom_count(),
        mol.bond_count() + record.duplicate_bonds.len()
    ));
    out.push_str(&record.molecule_type);
    out.push('\n');
    out.push_str(&record.charge_type);
    out.push('\n');
    for line in &record.molecule_tail {
        out.push_str(line);
        out.push('\n');
    }
    out.push('\n');

    // ATOM section.
    out.push_str("@<TRIPOS>ATOM\n");
    for (idx, atom) in mol.atoms() {
        let fallback_id = idx.0 + 1;
        let source = record.atoms.get(idx.0 as usize);
        let atom_id = source.map_or(fallback_id, |row| row.atom_id);
        let atom_name = source
            .map(|row| row.atom_name.as_str())
            .unwrap_or_else(|| atom.element.symbol());
        let atom_type = source
            .map(|row| row.atom_type.as_str())
            .unwrap_or_else(|| atom.element.symbol());
        let (x, y, z) = coords
            .get(idx.0 as usize)
            .copied()
            .unwrap_or((0.0, 0.0, 0.0));
        let subst_id = source.and_then(|row| row.subst_id).unwrap_or(1);
        let subst_name = source
            .and_then(|row| row.subst_name.as_deref())
            .unwrap_or("LIG");
        let partial_charge = source.and_then(|row| row.partial_charge).unwrap_or(0.0);
        let status = source
            .filter(|row| !row.status_bits.is_empty())
            .map(|row| format!(" {}", row.status_bits.join(" ")))
            .unwrap_or_default();
        out.push_str(&format!(
            "{atom_id:>6} {atom_name:<8} {x:>10.4} {y:>10.4} {z:>10.4} {atom_type:<8} {subst_id:>4}  {subst_name:<8} {partial_charge:>10.4}{status}\n"
        ));
    }

    let mut unity_groups = record.unity_atom_attributes.clone();
    for group in &mut unity_groups {
        let Some((idx, atom)) = mol.atoms().find(|(idx, _)| {
            record
                .atoms
                .get(idx.0 as usize)
                .is_some_and(|row| row.atom_id == group.atom_id)
        }) else {
            continue;
        };
        let _ = idx;
        if let Some((_, value)) = group
            .attributes
            .iter_mut()
            .find(|(name, _)| name.eq_ignore_ascii_case("charge"))
        {
            *value = atom.charge.to_string();
        } else if atom.charge != 0 {
            group
                .attributes
                .push(("charge".into(), atom.charge.to_string()));
        }
    }
    for (idx, atom) in mol.atoms().filter(|(_, atom)| atom.charge != 0) {
        let atom_id = record
            .atoms
            .get(idx.0 as usize)
            .map_or(idx.0 + 1, |row| row.atom_id);
        if !unity_groups.iter().any(|group| group.atom_id == atom_id) {
            unity_groups.push(Mol2UnityAtomAttributes {
                atom_id,
                attributes: vec![("charge".into(), atom.charge.to_string())],
            });
        }
    }
    if !unity_groups.is_empty() {
        out.push_str("@<TRIPOS>UNITY_ATOM_ATTR\n");
        for group in &unity_groups {
            out.push_str(&format!("{} {}\n", group.atom_id, group.attributes.len()));
            for (name, value) in &group.attributes {
                out.push_str(&format!("{name} {value}\n"));
            }
        }
    }

    // BOND section.
    out.push_str("@<TRIPOS>BOND\n");
    for (bidx, bond) in mol.bonds() {
        let source = record.bonds.get(bidx.0 as usize);
        let bi = source.map_or(bidx.0 + 1, |row| row.bond_id);
        let a1 = record
            .atoms
            .get(bond.atom1.0 as usize)
            .map_or(bond.atom1.0 + 1, |row| row.atom_id);
        let a2 = record
            .atoms
            .get(bond.atom2.0 as usize)
            .map_or(bond.atom2.0 + 1, |row| row.atom_id);
        let btype = source
            .map(|row| row.bond_type.as_str())
            .unwrap_or_else(|| bond_type_token(bond.order));
        let status = source
            .filter(|row| !row.status_bits.is_empty())
            .map(|row| format!(" {}", row.status_bits.join(" ")))
            .unwrap_or_default();
        out.push_str(&format!("{bi:>6} {a1:>6} {a2:>6} {btype}{status}\n"));
    }
    for bond in &record.duplicate_bonds {
        out.push_str(&format!(
            "{:>6} {:>6} {:>6} {}{}\n",
            bond.bond_id,
            bond.atom1_id,
            bond.atom2_id,
            bond.bond_type,
            if bond.status_bits.is_empty() {
                String::new()
            } else {
                format!(" {}", bond.status_bits.join(" "))
            }
        ));
    }

    for section in &record.opaque_sections {
        out.push_str("@<TRIPOS>");
        out.push_str(&section.name);
        out.push('\n');
        for line in &section.lines {
            out.push_str(line);
            out.push('\n');
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
        let (mol, coords) = parse_mol2(ETHANOL_MOL2).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        assert_eq!(coords.len(), 3);
    }

    #[test]
    fn test_parse_mol2_coords() {
        let (_, coords) = parse_mol2(ETHANOL_MOL2).unwrap();
        assert!((coords[0].0 - 0.0).abs() < 1e-4);
        assert!((coords[1].0 - 1.5).abs() < 1e-4);
        assert!((coords[2].0 - 3.0).abs() < 1e-4);
    }

    #[test]
    fn test_parse_mol2_elements() {
        let (mol, _) = parse_mol2(ETHANOL_MOL2).unwrap();
        assert_eq!(mol.atom(AtomIdx(0)).element.symbol(), "C");
        assert_eq!(mol.atom(AtomIdx(1)).element.symbol(), "C");
        assert_eq!(mol.atom(AtomIdx(2)).element.symbol(), "O");
    }

    #[test]
    fn test_write_mol2_roundtrip() {
        let (mol, coords) = parse_mol2(ETHANOL_MOL2).unwrap();
        let written = write_mol2(&mol, &coords);
        let (mol2, coords2) = parse_mol2(&written).unwrap();
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
    fn partial_charge_is_not_guessed_as_formal_charge() {
        let input = ETHANOL_MOL2.replace("-0.3940", "-0.7500");
        let record = parse_mol2_record(&input).unwrap();
        assert_eq!(record.molecule.atom(AtomIdx(2)).charge, 0);
        assert_eq!(record.atoms[2].partial_charge, Some(-0.75));
    }

    #[test]
    fn record_roundtrip_preserves_interchange_fields() {
        let input = ETHANOL_MOL2.replace(
            "@<TRIPOS>BOND",
            "@<TRIPOS>SUBSTRUCTURE\n     1 LIG1        1 RESIDUE\n@<TRIPOS>BOND",
        );
        let record = parse_mol2_record(&input).unwrap();
        assert_eq!(record.name, "ethanol");
        assert_eq!(record.molecule_type, "SMALL");
        assert_eq!(record.charge_type, "GASTEIGER");
        assert_eq!(record.atoms[0].atom_name, "C1");
        assert_eq!(record.atoms[0].atom_type, "C.3");
        assert_eq!(record.atoms[2].partial_charge, Some(-0.394));
        assert_eq!(record.atoms[0].subst_name.as_deref(), Some("LIG1"));
        assert_eq!(record.opaque_sections[0].name, "SUBSTRUCTURE");

        let rewritten = write_mol2_record(&record);
        let reparsed = parse_mol2_record(&rewritten).unwrap();
        assert_eq!(reparsed.name, record.name);
        assert_eq!(reparsed.molecule_type, record.molecule_type);
        assert_eq!(reparsed.charge_type, record.charge_type);
        assert_eq!(reparsed.atoms, record.atoms);
        assert_eq!(reparsed.bonds, record.bonds);
        assert_eq!(reparsed.opaque_sections, record.opaque_sections);
    }

    #[test]
    fn unity_formal_charge_is_distinct_from_partial_charge() {
        let input = "@<TRIPOS>MOLECULE\nammonium\n 1 0 0 0 0\nSMALL\nGASTEIGER\n\n@<TRIPOS>ATOM\n      1 N 0.0 0.0 0.0 N.4 1 UNL1 0.2500\n@<TRIPOS>UNITY_ATOM_ATTR\n1 2\ncharge 1\ncolor red\n@<TRIPOS>BOND\n";
        let record = parse_mol2_record(input).unwrap();
        assert_eq!(record.molecule.atom(AtomIdx(0)).charge, 1);
        assert_eq!(record.atoms[0].partial_charge, Some(0.25));
        assert_eq!(
            record.unity_atom_attributes[0].attributes[1],
            ("color".into(), "red".into())
        );

        let reparsed = parse_mol2_record(&write_mol2_record(&record)).unwrap();
        assert_eq!(reparsed.molecule.atom(AtomIdx(0)).charge, 1);
        assert_eq!(reparsed.atoms[0].partial_charge, Some(0.25));
        assert_eq!(reparsed.unity_atom_attributes, record.unity_atom_attributes);
    }

    #[test]
    fn duplicate_bond_rows_are_retained_for_same_format_rewrite() {
        let input = "@<TRIPOS>MOLECULE\nduplicate\n 2 2 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 LIG 0\n2 C2 1 0 0 C.3 1 LIG 0\n@<TRIPOS>BOND\n1 1 2 1\n2 1 2 1 DUP\n";
        let record = parse_mol2_record(input).unwrap();
        assert_eq!(record.molecule.bond_count(), 1);
        assert_eq!(record.duplicate_bonds.len(), 1);
        let reparsed = parse_mol2_record(&write_mol2_record(&record)).unwrap();
        assert_eq!(reparsed.molecule.bond_count(), 1);
        assert_eq!(reparsed.duplicate_bonds, record.duplicate_bonds);
    }

    #[test]
    fn aromatic_atom_type_sets_aromatic_atom_state() {
        let record = parse_mol2_record(
            "@<TRIPOS>MOLECULE\nbenzene-fragment\n 2 1 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.ar 1 BNZ 0\n2 C2 1 0 0 C.ar 1 BNZ 0\n@<TRIPOS>BOND\n1 1 2 ar\n",
        )
        .unwrap();
        assert!(record.molecule.atom(AtomIdx(0)).aromatic);
        assert!(record.molecule.atom(AtomIdx(1)).aromatic);
    }

    #[test]
    fn shared_mol2_roundtrip_contract_matches() {
        let document: serde_json::Value = serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/cross_binding_contract.json"
        )))
        .expect("contract JSON");
        let contract = &document["mol2_contract"];
        let (mol, coords) = parse_mol2(contract["input"].as_str().unwrap()).unwrap();
        let serialized = write_mol2(&mol, &coords);
        let (roundtripped, _) = parse_mol2(&serialized).unwrap();
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
        let (mol, _) = parse_mol2(benzene_mol2).unwrap();
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
        assert!(parse_mol2(bad).is_err());
    }

    #[test]
    fn multi_molecule_stream_is_a_typed_refusal() {
        let input = format!("{ETHANOL_MOL2}{ETHANOL_MOL2}");
        assert!(matches!(
            parse_mol2_record(&input),
            Err(Mol2Error::MultipleMolecules { count: 2 })
        ));
    }

    #[test]
    fn bounded_parser_rejects_input_and_line_limits() {
        assert!(matches!(
            parse_mol2_with_limits(
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
            parse_mol2_with_limits(
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
            parse_mol2_with_limits(
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
            parse_mol2_with_limits(
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
