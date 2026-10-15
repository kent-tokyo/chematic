//! General residue grouping from atom_site identity, including non-polymers.
//! Independent implementation from the wwPDB dictionary; no bond perception.
use super::*;
use std::collections::{HashMap, HashSet};

/// Complete stored residue namespace. Both author and label identifiers matter.
/// Alternate component names at a position form separate residue variants.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MmcifResidueIdentity {
    pub model_num: i32,
    pub label_asym_id: Option<String>,
    pub auth_asym_id: Option<String>,
    pub entity_id: Option<String>,
    pub label_seq_id: Option<i64>,
    pub auth_seq_id: Option<String>,
    pub insertion_code: Option<String>,
    pub label_comp_id: Option<String>,
    pub auth_comp_id: Option<String>,
}

/// One residue variant and zero-based source row indices, in first-seen order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifResidue {
    pub identity: MmcifResidueIdentity,
    pub atom_indices: Vec<usize>,
}

/// Full atom-name aliases and alternate-location code, without char truncation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifAtomIdentity {
    pub label_atom_id: Option<String>,
    pub auth_atom_id: Option<String>,
    pub alt_id: Option<String>,
}

/// Residue view over every atom_site record (including water/ligands/glycans).
/// `parsed` retains the existing numeric/Cartesian record representation;
/// residue/atom identities here are authoritative for the original aliases.
#[derive(Debug, Clone, PartialEq)]
pub struct MmcifResidueDocument {
    pub parsed: MmcifResult,
    pub residues: Vec<MmcifResidue>,
    pub atom_ids: Vec<MmcifAtomIdentity>,
    pub atom_to_residue: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum MmcifResidueError {
    Parse(MmcifError),
    MissingIdentity { atom: usize, field: &'static str },
    DuplicateAtom { atom: usize, residue: usize },
    InvalidField { atom: usize, field: &'static str },
    IncompleteRow { atom: usize },
}
impl core::fmt::Display for MmcifResidueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Parse(e) => e.fmt(f),
            Self::MissingIdentity { atom, field } => write!(
                f,
                "mmCIF atom {atom} lacks {field}; residue identity cannot be inferred"
            ),
            Self::DuplicateAtom { atom, residue } => write!(
                f,
                "mmCIF atom {atom} repeats an atom name/alternate location in residue {residue}"
            ),
            Self::InvalidField { atom, field } => {
                write!(f, "mmCIF atom {atom} has an invalid {field}")
            }
            Self::IncompleteRow { atom } => write!(f, "mmCIF atom_site row {atom} is incomplete"),
        }
    }
}
impl std::error::Error for MmcifResidueError {}
impl From<MmcifError> for MmcifResidueError {
    fn from(value: MmcifError) -> Self {
        Self::Parse(value)
    }
}

fn value(row: &[String], columns: &MmcifColumns, col: Option<usize>) -> Option<String> {
    columns
        .get(row, col)
        .filter(|s| !is_cif_placeholder(s))
        .map(str::to_owned)
}

fn identities(
    row: &[String],
    c: &MmcifColumns,
    atom: usize,
) -> Result<(MmcifResidueIdentity, MmcifAtomIdentity), MmcifResidueError> {
    let model = optional_mmcif_int(c.get(row, c.model), "pdbx_PDB_model_num")?.unwrap_or(1);
    let model_num =
        i32::try_from(model)
            .ok()
            .filter(|m| *m > 0)
            .ok_or(MmcifResidueError::InvalidField {
                atom,
                field: "model number",
            })?;
    if let Some(charge) = optional_mmcif_int(c.get(row, c.charge), "pdbx_formal_charge")?
        && i32::try_from(charge).is_err()
    {
        return Err(MmcifResidueError::InvalidField {
            atom,
            field: "formal charge",
        });
    }
    let residue = MmcifResidueIdentity {
        model_num,
        label_asym_id: value(row, c, c.label_asym_id),
        auth_asym_id: value(row, c, c.auth_asym_id),
        entity_id: value(row, c, c.entity_id),
        label_seq_id: optional_mmcif_int(c.get(row, c.label_seq_id), "label_seq_id")?,
        auth_seq_id: value(row, c, c.auth_seq_id),
        insertion_code: value(row, c, c.ins_code),
        label_comp_id: value(row, c, c.label_comp_id),
        auth_comp_id: value(row, c, c.auth_comp_id),
    };
    let atom_id = MmcifAtomIdentity {
        label_atom_id: value(row, c, c.label_atom_id),
        auth_atom_id: value(row, c, c.auth_atom_id),
        alt_id: value(row, c, c.alt_id),
    };
    for (missing, field) in [
        (
            residue.label_asym_id.is_none() && residue.auth_asym_id.is_none(),
            "chain identifier",
        ),
        (
            residue.label_seq_id.is_none() && residue.auth_seq_id.is_none(),
            "residue sequence identifier",
        ),
        (
            residue.label_comp_id.is_none() && residue.auth_comp_id.is_none(),
            "component identifier",
        ),
        (
            atom_id.label_atom_id.is_none() && atom_id.auth_atom_id.is_none(),
            "atom identifier",
        ),
    ] {
        if missing {
            return Err(MmcifResidueError::MissingIdentity { atom, field });
        }
    }
    Ok((residue, atom_id))
}

/// Parse atom_site into residues with default mmCIF resource limits.
/// Missing label_seq_id never means residue zero; auth_seq_id distinguishes
/// non-polymer residues. Missing both sequence identifiers is an explicit error.
pub fn parse_mmcif_residues(input: &str) -> Result<MmcifResidueDocument, MmcifResidueError> {
    parse_mmcif_residues_with_limits(input, &MmcifParseLimits::default())
}

/// Parse with caller-specified byte, line and atom limits. Identity matches
/// require every stored namespace field to match; partial aliases are not guessed.
/// Like the existing numeric atom reader, noninteger auth_seq_id is unsupported.
pub fn parse_mmcif_residues_with_limits(
    input: &str,
    limits: &MmcifParseLimits,
) -> Result<MmcifResidueDocument, MmcifResidueError> {
    let tokens = mmcif_tokens(input, limits)?;
    let (headers, start) =
        find_loop_by_category(&tokens, "atom_site").ok_or(MmcifError::NoAtomSiteLoop)?;
    let c = MmcifColumns::from_headers(&headers)?;
    let (cell, space_group) = scan_cell_and_symmetry(&tokens);
    let known = c.indices();
    let unhandled_columns = headers
        .iter()
        .enumerate()
        .filter(|(i, _)| !known.contains(i))
        .map(|(_, h)| format!("_atom_site.{h}"))
        .collect();
    let mut doc = MmcifResidueDocument {
        parsed: MmcifResult {
            atoms: Vec::new(),
            cell,
            space_group,
            unhandled_columns,
        },
        residues: Vec::new(),
        atom_ids: Vec::new(),
        atom_to_residue: Vec::new(),
    };
    let mut lookup = HashMap::new();
    let mut seen = HashSet::new();
    for row in tokens[start..].chunks(c.width) {
        if row[0] == "loop_" || row[0].starts_with("data_") || row[0].starts_with('_') {
            break;
        }
        let atom = doc.parsed.atoms.len();
        if row.len() != c.width {
            return Err(MmcifResidueError::IncompleteRow { atom });
        }
        if atom >= limits.max_atoms {
            return Err(MmcifError::TooManyAtoms {
                limit: limits.max_atoms,
            }
            .into());
        }
        let (residue, atom_id) = identities(row, &c, atom)?;
        let record = parse_mmcif_atom(row, atom, &c)?;
        let next = doc.residues.len();
        let index = *lookup.entry(residue.clone()).or_insert_with(|| {
            doc.residues.push(MmcifResidue {
                identity: residue,
                atom_indices: Vec::new(),
            });
            next
        });
        // Label atom names are authoritative when supplied; author names are
        // only the fallback. Different alt IDs stay separate atoms in a residue.
        let name = atom_id
            .label_atom_id
            .as_ref()
            .or(atom_id.auth_atom_id.as_ref())
            .cloned();
        if !seen.insert((index, name, atom_id.alt_id.clone())) {
            return Err(MmcifResidueError::DuplicateAtom {
                atom,
                residue: index,
            });
        }
        doc.residues[index].atom_indices.push(atom);
        doc.atom_to_residue.push(index);
        doc.atom_ids.push(atom_id);
        doc.parsed.atoms.push(record);
    }
    if doc.parsed.atoms.is_empty() {
        return Err(MmcifError::NoAtomSiteLoop.into());
    }
    Ok(doc)
}
