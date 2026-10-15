//! Checked single-model atom_site to PDB coordinate-record projection.
//! Independently implemented from wwPDB format 3.3; no topology inference.
use crate::{MmcifAtomRecord, MmcifResidueDocument, MmcifResidueIdentity};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PdbNamespace {
    Author,
    Label,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifPdbError {
    EmptyInput,
    InvalidDocument,
    UnsupportedModels,
    MetadataLossRequiresConsent,
    InvalidField { atom: usize, field: &'static str },
    DuplicateSerial { atom: usize, serial: i64 },
    ResidueCollision { first: usize, second: usize },
    AtomCollision { first: usize, second: usize },
}
impl core::fmt::Display for MmcifPdbError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "PDB coordinate output requires at least one atom"),
            Self::InvalidDocument => write!(
                f,
                "mmCIF residue document has inconsistent atom/residue references"
            ),
            Self::UnsupportedModels => {
                write!(f, "checked PDB coordinate output supports model 1 only")
            }
            Self::MetadataLossRequiresConsent => write!(
                f,
                "PDB cannot preserve complete mmCIF identity namespaces; set allow_metadata_loss=true and retain the source sidecar"
            ),
            Self::InvalidField { atom, field } => write!(
                f,
                "mmCIF atom {atom}: {field} is missing, invalid or cannot fit its PDB field"
            ),
            Self::DuplicateSerial { atom, serial } => {
                write!(f, "mmCIF atom {atom}: duplicate PDB serial {serial}")
            }
            Self::ResidueCollision { first, second } => write!(
                f,
                "mmCIF residues {first} and {second} collapse to the same PDB residue identifier"
            ),
            Self::AtomCollision { first, second } => write!(
                f,
                "mmCIF atoms {first} and {second} collapse to the same PDB atom name/alternate location"
            ),
        }
    }
}
impl std::error::Error for MmcifPdbError {}

#[derive(Debug, Clone, PartialEq)]
pub struct MmcifPdbReport {
    pub pdb: String,
    /// Stable codes for information absent or rounded in the PDB block.
    pub losses: Vec<&'static str>,
    pub max_coordinate_rounding_error: f64,
    pub max_occupancy_rounding_error: f64,
    pub max_b_iso_rounding_error: f64,
}

#[derive(Clone)]
struct Position {
    chain: String,
    seq: i64,
    icode: char,
    name: String,
}

fn text(
    value: Option<&str>,
    width: usize,
    atom: usize,
    field: &'static str,
) -> Result<String, MmcifPdbError> {
    value
        .filter(|s| !s.is_empty() && s.len() <= width && s.bytes().all(|b| b.is_ascii_graphic()))
        .map(str::to_owned)
        .ok_or(MmcifPdbError::InvalidField { atom, field })
}
fn character(value: Option<&str>, atom: usize, field: &'static str) -> Result<char, MmcifPdbError> {
    match value {
        None => Ok(' '),
        Some(s) => text(Some(s), 1, atom, field).and_then(|s| {
            if s.as_bytes()[0].is_ascii_alphanumeric() {
                Ok(s.as_bytes()[0] as char)
            } else {
                Err(MmcifPdbError::InvalidField { atom, field })
            }
        }),
    }
}
fn position(
    id: &MmcifResidueIdentity,
    namespace: PdbNamespace,
    atom: usize,
) -> Result<Position, MmcifPdbError> {
    let (chain, name, seq) = match namespace {
        PdbNamespace::Author => (
            id.auth_asym_id.as_deref(),
            id.auth_comp_id.as_deref(),
            id.auth_seq_id
                .as_ref()
                .and_then(|s| s.parse::<i64>().ok().filter(|n| n.to_string() == *s)),
        ),
        PdbNamespace::Label => (
            id.label_asym_id.as_deref(),
            id.label_comp_id.as_deref(),
            id.label_seq_id,
        ),
    };
    let seq = seq
        .filter(|v| (-999..=9999).contains(v))
        .ok_or(MmcifPdbError::InvalidField {
            atom,
            field: "residue number",
        })?;
    let chain = text(chain, 1, atom, "chain identifier")?;
    if !chain.as_bytes()[0].is_ascii_alphanumeric() {
        return Err(MmcifPdbError::InvalidField {
            atom,
            field: "chain identifier",
        });
    }
    Ok(Position {
        chain,
        name: text(name, 3, atom, "residue name")?,
        seq,
        icode: character(id.insertion_code.as_deref(), atom, "insertion code")?,
    })
}
fn number(
    v: f64,
    width: usize,
    precision: usize,
    atom: usize,
    field: &'static str,
) -> Result<(String, f64), MmcifPdbError> {
    let invalid = || MmcifPdbError::InvalidField { atom, field };
    if !v.is_finite() {
        return Err(invalid());
    }
    let s = format!("{v:>width$.precision$}");
    if s.len() != width {
        return Err(invalid());
    }
    let rounded: f64 = s.trim().parse().map_err(|_| invalid())?;
    Ok((s, (rounded - v).abs()))
}
fn row(
    a: &MmcifAtomRecord,
    name: &str,
    alt: char,
    p: &Position,
    i: usize,
) -> Result<(String, [f64; 3]), MmcifPdbError> {
    if !matches!(a.group_pdb.as_str(), "ATOM" | "HETATM") {
        return Err(MmcifPdbError::InvalidField {
            atom: i,
            field: "record kind",
        });
    }
    if !(1..=99999).contains(&a.serial) {
        return Err(MmcifPdbError::InvalidField {
            atom: i,
            field: "serial",
        });
    }
    let name = if name.len() < 4
        && a.element.symbol().len() == 1
        && !name.as_bytes()[0].is_ascii_digit()
    {
        format!(" {name:<3}")
    } else {
        format!("{name:<4}")
    };
    let (x, dx) = number(a.x, 8, 3, i, "x")?;
    let (y, dy) = number(a.y, 8, 3, i, "y")?;
    let (z, dz) = number(a.z, 8, 3, i, "z")?;
    if !(0.0..=1.0).contains(&a.occupancy) {
        return Err(MmcifPdbError::InvalidField {
            atom: i,
            field: "occupancy",
        });
    }
    let (occupancy, do_) = number(a.occupancy, 6, 2, i, "occupancy")?;
    let (b, db) = number(a.b_iso, 6, 2, i, "B factor")?;
    let q = a.formal_charge.unwrap_or(0);
    if !(-9..=9).contains(&q) {
        return Err(MmcifPdbError::InvalidField {
            atom: i,
            field: "formal charge",
        });
    }
    let charge = if q == 0 {
        "  ".into()
    } else {
        format!("{}{}", q.abs(), if q > 0 { '+' } else { '-' })
    };
    let out = format!(
        "{:<6}{:>5} {}{}{:>3} {}{:>4}{}   {}{}{}{}{}          {:>2}{}\n",
        a.group_pdb,
        a.serial,
        name,
        alt,
        p.name,
        p.chain,
        p.seq,
        p.icode,
        x,
        y,
        z,
        occupancy,
        b,
        a.element.symbol().to_ascii_uppercase(),
        charge
    );
    Ok((out, [dx.max(dy).max(dz), do_, db]))
}

fn validate_document(doc: &MmcifResidueDocument) -> Result<(), MmcifPdbError> {
    let n = doc.parsed.atoms.len();
    if n == 0 {
        return Err(MmcifPdbError::EmptyInput);
    }
    if n > 99999 {
        return Err(MmcifPdbError::InvalidField {
            atom: 99999,
            field: "atom count",
        });
    }
    if doc.atom_ids.len() != n || doc.atom_to_residue.len() != n {
        return Err(MmcifPdbError::InvalidDocument);
    }
    let mut seen = vec![false; n];
    for (r, residue) in doc.residues.iter().enumerate() {
        if residue.atom_indices.is_empty() {
            return Err(MmcifPdbError::InvalidDocument);
        }
        if residue.identity.model_num != 1 {
            return Err(MmcifPdbError::UnsupportedModels);
        }
        for &a in &residue.atom_indices {
            if a >= n
                || seen[a]
                || doc.atom_to_residue[a] != r
                || doc.parsed.atoms[a].model_num != 1
            {
                return Err(MmcifPdbError::InvalidDocument);
            }
            seen[a] = true;
        }
    }
    if seen.iter().any(|s| !*s) {
        return Err(MmcifPdbError::InvalidDocument);
    }
    Ok(())
}

/// Project one model's atom_site records to fixed-width PDB coordinate lines.
/// Names and identities never truncate or renumber, even after loss consent.
/// Complete original namespaces live in the input document; retain it as the
/// sidecar. The block contains ATOM/HETATM and END only (no topology/TER).
pub fn write_mmcif_pdb_checked(
    doc: &MmcifResidueDocument,
    namespace: PdbNamespace,
    allow_metadata_loss: bool,
) -> Result<MmcifPdbReport, MmcifPdbError> {
    validate_document(doc)?;
    let mut positions = Vec::with_capacity(doc.residues.len());
    let mut projected = HashMap::new();
    for (r, residue) in doc.residues.iter().enumerate() {
        let p = position(&residue.identity, namespace, residue.atom_indices[0])?;
        let key = (p.chain.clone(), p.seq, p.icode);
        if let Some(first) = projected.insert(key, r) {
            return Err(MmcifPdbError::ResidueCollision { first, second: r });
        }
        positions.push(p);
    }
    let mut serials = HashSet::new();
    let mut atom_names = HashMap::new();
    let mut report = MmcifPdbReport {
        pdb: String::new(),
        losses: vec!["identity_namespaces"],
        max_coordinate_rounding_error: 0.0,
        max_occupancy_rounding_error: 0.0,
        max_b_iso_rounding_error: 0.0,
    };
    for (i, a) in doc.parsed.atoms.iter().enumerate() {
        if !serials.insert(a.serial) {
            return Err(MmcifPdbError::DuplicateSerial {
                atom: i,
                serial: a.serial,
            });
        }
        let id = &doc.atom_ids[i];
        let selected = match namespace {
            PdbNamespace::Author => id.auth_atom_id.as_deref(),
            PdbNamespace::Label => id.label_atom_id.as_deref(),
        };
        let name = text(selected, 4, i, "atom name")?;
        let alt = character(id.alt_id.as_deref(), i, "alternate location")?;
        if let Some(first) = atom_names.insert((doc.atom_to_residue[i], name.clone(), alt), i) {
            return Err(MmcifPdbError::AtomCollision { first, second: i });
        }
        let (line, errors) = row(a, &name, alt, &positions[doc.atom_to_residue[i]], i)?;
        report.pdb.push_str(&line);
        report.max_coordinate_rounding_error = report.max_coordinate_rounding_error.max(errors[0]);
        report.max_occupancy_rounding_error = report.max_occupancy_rounding_error.max(errors[1]);
        report.max_b_iso_rounding_error = report.max_b_iso_rounding_error.max(errors[2]);
    }
    if doc.residues.iter().any(|r| r.identity.entity_id.is_some()) {
        report.losses.push("entity_identifiers");
    }
    if doc.parsed.cell.is_some() || doc.parsed.space_group.is_some() {
        report.losses.push("cell_and_symmetry");
    }
    if !doc.parsed.unhandled_columns.is_empty() {
        report.losses.push("unhandled_atom_site_columns");
    }
    if report.max_coordinate_rounding_error > 0.0
        || report.max_occupancy_rounding_error > 0.0
        || report.max_b_iso_rounding_error > 0.0
    {
        report.losses.push("numeric_rounding");
    }
    if doc.parsed.atoms.iter().any(|a| a.formal_charge == Some(0)) {
        report.losses.push("optional_formal_charge_presence");
    }
    if !allow_metadata_loss {
        return Err(MmcifPdbError::MetadataLossRequiresConsent);
    }
    report.pdb.push_str("END\n");
    Ok(report)
}
