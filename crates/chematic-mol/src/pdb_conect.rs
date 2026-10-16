//! Explicit PDB CONECT graph projection without distance-based bond inference.
use std::collections::BTreeMap;

use chematic_3d::{PdbParseLimits, PdbStrictError, parse_pdb_atoms_strict};
use chematic_core::{
    Atom, AtomIdx, BondOrder, Coords3D, Element, Molecule, MoleculeBuilder, Point3,
};

#[derive(Debug, Clone, Copy, Default)]
pub struct PdbConectOptions {
    /// Accept target slots after the four modern covalent fields. Legacy PDB
    /// uses those columns for other contacts, so this requires explicit consent.
    pub allow_extended_targets: bool,
}
#[derive(Debug, Clone, Copy)]
pub struct PdbConectLimits {
    pub atoms: PdbParseLimits,
    pub max_bonds: usize,
    pub max_targets: usize,
}
impl Default for PdbConectLimits {
    fn default() -> Self {
        Self {
            atoms: PdbParseLimits::default(),
            max_bonds: 30_000,
            max_targets: 180_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum PdbConectError {
    Atoms(PdbStrictError),
    InvalidRecord {
        line: usize,
        field: &'static str,
    },
    UnsupportedRecord {
        line: usize,
        record: &'static str,
    },
    DuplicateSerial {
        line: usize,
        serial: u32,
    },
    UnknownSerial {
        line: usize,
        serial: u32,
    },
    Multiplicity {
        source: u32,
        target: u32,
        forward: u8,
        reverse: u8,
    },
    Limit {
        resource: &'static str,
        limit: usize,
    },
}
impl core::fmt::Display for PdbConectError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Atoms(e) => e.fmt(f),
            Self::InvalidRecord { line, field } => write!(f, "invalid PDB {field} on line {line}"),
            Self::UnsupportedRecord { line, record } => {
                write!(f, "unsupported PDB {record} on line {line}")
            }
            Self::DuplicateSerial { line, serial } => {
                write!(f, "duplicate PDB serial {serial} on line {line}")
            }
            Self::UnknownSerial { line, serial } => {
                write!(f, "unknown PDB serial {serial} on line {line}")
            }
            Self::Multiplicity {
                source,
                target,
                forward,
                reverse,
            } => write!(
                f,
                "invalid PDB CONECT multiplicity for {source}-{target}: {forward}/{reverse}"
            ),
            Self::Limit { resource, limit } => {
                write!(f, "PDB CONECT {resource} exceeds limit {limit}")
            }
        }
    }
}
impl std::error::Error for PdbConectError {}

pub struct PdbConectFrame {
    pub molecule: Molecule,
    /// Cartesian source coordinates, in atom-record order (including H).
    pub coords: Coords3D,
    /// PDB serials in exactly the same order as molecule atoms and coordinates.
    pub serials: Vec<u32>,
}

fn invalid(line: usize, field: &'static str) -> PdbConectError {
    PdbConectError::InvalidRecord { line, field }
}
fn serial(text: &str, line: usize) -> Result<u32, PdbConectError> {
    text.trim()
        .parse()
        .ok()
        .filter(|&n| n > 0)
        .ok_or_else(|| invalid(line, "serial"))
}
fn element(text: &str) -> Option<Element> {
    Element::from_symbol(text).or_else(|| {
        let [a, b] = text.as_bytes() else {
            return None;
        };
        if !a.is_ascii_uppercase() || !b.is_ascii_uppercase() {
            return None;
        }
        let normalized = [*a, b.to_ascii_lowercase()];
        std::str::from_utf8(&normalized)
            .ok()
            .and_then(Element::from_symbol)
    })
}
fn charge(line: &str, number: usize) -> Result<i8, PdbConectError> {
    let text = line
        .get(78..80)
        .or_else(|| line.get(78..))
        .unwrap_or("")
        .trim();
    match text.as_bytes() {
        [] => Ok(0),
        [n @ b'1'..=b'9', sign @ (b'+' | b'-')] => {
            Ok((*n - b'0') as i8 * if *sign == b'+' { 1 } else { -1 })
        }
        _ => Err(invalid(number, "formal charge")),
    }
}

/// Read the unmodelled, single-location PDB CONECT profile with default limits.
/// Repeated targets mean bond multiplicity 1..=3. Reciprocal declarations must
/// agree; they are not added together. Missing edges remain absent. Atom names
/// never determine elements. No proximity, aromaticity or stereo perception is
/// performed; implicit H follows the existing native molecule valence model.
pub fn parse_pdb_conect(input: &str) -> Result<PdbConectFrame, PdbConectError> {
    parse_pdb_conect_with_options(input, &PdbConectOptions::default())
}
pub fn parse_pdb_conect_with_options(
    input: &str,
    options: &PdbConectOptions,
) -> Result<PdbConectFrame, PdbConectError> {
    parse_pdb_conect_with_options_and_limits(input, options, &PdbConectLimits::default())
}
pub fn parse_pdb_conect_with_options_and_limits(
    input: &str,
    options: &PdbConectOptions,
    limits: &PdbConectLimits,
) -> Result<PdbConectFrame, PdbConectError> {
    let atoms = parse_pdb_atoms_strict(input, &limits.atoms).map_err(PdbConectError::Atoms)?;
    let mut builder = MoleculeBuilder::with_capacity(atoms.len(), atoms.len());
    let mut indices = BTreeMap::new();
    let mut serials = Vec::new();
    let mut atom_index = 0;
    let mut ended = false;
    for (i, line) in input.lines().enumerate() {
        let record = line.get(..6).unwrap_or(line).trim_end();
        match record {
            "MODEL" | "ENDMDL" => {
                return Err(PdbConectError::UnsupportedRecord {
                    line: i + 1,
                    record: "MODEL/ENDMDL",
                });
            }
            "END" => ended = true,
            "ATOM" | "HETATM" => {
                if ended || !line.is_ascii() || line.as_bytes().get(16).is_some_and(|&b| b != b' ')
                {
                    return Err(invalid(
                        i + 1,
                        "atom profile (after END, non-ASCII or alternate location)",
                    ));
                }
                let parsed = atoms
                    .get(atom_index)
                    .ok_or_else(|| invalid(i + 1, "atom width"))?;
                let elem = element(line.get(76..78).unwrap_or("").trim())
                    .ok_or_else(|| invalid(i + 1, "element"))?;
                let id = serial(line.get(6..11).unwrap_or(""), i + 1)?;
                if indices.contains_key(&id) {
                    return Err(PdbConectError::DuplicateSerial {
                        line: i + 1,
                        serial: id,
                    });
                }
                let mut atom = Atom::new(elem);
                atom.charge = charge(line, i + 1)?;
                indices.insert(id, builder.add_atom(atom));
                serials.push(parsed.serial);
                atom_index += 1;
            }
            "CONECT" if ended => return Err(invalid(i + 1, "CONECT after END")),
            _ => {}
        }
    }
    let bonds = read_bonds(input, options, limits, &indices)?;
    for ((a, b), (forward, reverse, line)) in bonds {
        if forward > 0 && reverse > 0 && forward != reverse {
            return Err(PdbConectError::Multiplicity {
                source: a,
                target: b,
                forward,
                reverse,
            });
        }
        let order = match forward.max(reverse) {
            1 => BondOrder::Single,
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            _ => {
                return Err(PdbConectError::Multiplicity {
                    source: a,
                    target: b,
                    forward,
                    reverse,
                });
            }
        };
        builder
            .add_bond(indices[&a], indices[&b], order)
            .map_err(|_| invalid(line, "bond"))?;
    }
    Ok(PdbConectFrame {
        molecule: builder.build(),
        coords: Coords3D {
            points: atoms.iter().map(|a| Point3::new(a.x, a.y, a.z)).collect(),
        },
        serials,
    })
}

type Bonds = BTreeMap<(u32, u32), (u8, u8, usize)>;
fn read_bonds(
    input: &str,
    options: &PdbConectOptions,
    limits: &PdbConectLimits,
    indices: &BTreeMap<u32, AtomIdx>,
) -> Result<Bonds, PdbConectError> {
    let mut bonds = BTreeMap::new();
    let mut targets = 0;
    for (i, line) in input.lines().enumerate() {
        if !line.starts_with("CONECT") {
            continue;
        }
        if !line.is_ascii() {
            return Err(invalid(i + 1, "CONECT ASCII"));
        }
        let source = serial(
            line.get(6..11)
                .ok_or_else(|| invalid(i + 1, "CONECT source"))?,
            i + 1,
        )?;
        if !indices.contains_key(&source) {
            return Err(PdbConectError::UnknownSerial {
                line: i + 1,
                serial: source,
            });
        }
        for (slot, chunk) in line.as_bytes()[11..].chunks(5).enumerate() {
            if chunk.iter().all(u8::is_ascii_whitespace) {
                continue;
            }
            if chunk.len() != 5 {
                return Err(invalid(i + 1, "CONECT target width"));
            }
            if slot >= 4 && !options.allow_extended_targets {
                return Err(PdbConectError::UnsupportedRecord {
                    line: i + 1,
                    record: "extended CONECT target slots",
                });
            }
            if targets == limits.max_targets {
                return Err(PdbConectError::Limit {
                    resource: "targets",
                    limit: limits.max_targets,
                });
            }
            targets += 1;
            let target = serial(
                std::str::from_utf8(chunk).map_err(|_| invalid(i + 1, "CONECT target"))?,
                i + 1,
            )?;
            if !indices.contains_key(&target) {
                return Err(PdbConectError::UnknownSerial {
                    line: i + 1,
                    serial: target,
                });
            }
            if source == target {
                return Err(invalid(i + 1, "CONECT self bond"));
            }
            let pair = (source.min(target), source.max(target));
            if !bonds.contains_key(&pair) && bonds.len() == limits.max_bonds {
                return Err(PdbConectError::Limit {
                    resource: "bonds",
                    limit: limits.max_bonds,
                });
            }
            let entry = bonds.entry(pair).or_insert((0u8, 0u8, i + 1));
            let count = if source < target {
                &mut entry.0
            } else {
                &mut entry.1
            };
            *count += 1;
            if *count > 3 {
                return Err(PdbConectError::Multiplicity {
                    source: pair.0,
                    target: pair.1,
                    forward: entry.0,
                    reverse: entry.1,
                });
            }
        }
    }
    Ok(bonds)
}
