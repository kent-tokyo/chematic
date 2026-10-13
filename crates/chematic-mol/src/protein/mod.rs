//! Protein-focused projections for PDB and mmCIF structures.
//!
//! This module deliberately does not call the result a complete biological
//! structure.  [`ProteinProjection`] keeps amino-acid residues only; waters,
//! ligands, ions, and nucleic acids remain available from the underlying
//! format readers but are excluded from this view.

use std::collections::{HashMap, HashSet};

use chematic_core::Element;

use crate::{MmcifError, MmcifParseLimits, parse_mmcif_with_limits};

/// One atom in a protein-focused structural projection.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct ProteinAtom {
    pub serial: i64,
    pub name: String,
    pub element: Element,
    pub position: [f64; 3],
    pub occupancy: f64,
    pub b_iso: f64,
    pub alt_loc: Option<char>,
}

/// One amino-acid residue in source order.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct ProteinResidue {
    pub name: String,
    pub sequence_number: i64,
    pub insertion_code: Option<char>,
    pub canonical_one_letter_code: char,
    pub modified: bool,
    pub atoms: Vec<ProteinAtom>,
}

/// One chain in one coordinate model.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct ProteinChain {
    pub model_number: i32,
    pub id: String,
    pub residues: Vec<ProteinResidue>,
}

impl ProteinChain {
    /// Canonical one-letter amino-acid sequence for this chain.
    pub fn sequence(&self) -> String {
        self.residues
            .iter()
            .map(|residue| residue.canonical_one_letter_code)
            .collect()
    }

    pub fn atom_count(&self) -> usize {
        self.residues
            .iter()
            .map(|residue| residue.atoms.len())
            .sum()
    }
}

/// Amino-acid-only view of a PDB or mmCIF structure.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ProteinProjection {
    pub chains: Vec<ProteinChain>,
}

impl ProteinProjection {
    pub fn model_count(&self) -> usize {
        self.chains
            .iter()
            .map(|chain| chain.model_number)
            .collect::<HashSet<_>>()
            .len()
    }

    pub fn chain_count(&self) -> usize {
        self.chains.len()
    }

    pub fn residue_count(&self) -> usize {
        self.chains.iter().map(|chain| chain.residues.len()).sum()
    }

    pub fn atom_count(&self) -> usize {
        self.chains.iter().map(ProteinChain::atom_count).sum()
    }
}

/// Failure while building a protein projection.
#[non_exhaustive]
#[derive(Debug)]
pub enum ProteinProjectionError {
    Pdb(chematic_3d::PdbStrictError),
    Mmcif(MmcifError),
    UnknownPdbElement { serial: u32, symbol: String },
}

impl core::fmt::Display for ProteinProjectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Pdb(error) => write!(f, "PDB protein projection failed: {error}"),
            Self::Mmcif(error) => write!(f, "mmCIF protein projection failed: {error}"),
            Self::UnknownPdbElement { serial, symbol } => write!(
                f,
                "PDB protein projection could not resolve element '{symbol}' for atom {serial}"
            ),
        }
    }
}

impl std::error::Error for ProteinProjectionError {}

#[derive(Debug, Clone, Copy)]
struct ResidueClass {
    canonical: char,
    modified: bool,
}

/// Return the canonical one-letter code and modification status for a residue.
///
/// The table covers the 20 standard amino acids, selenocysteine,
/// pyrrolysine, `UNK`, and common PDB modifications/protonation variants. It
/// is intentionally bounded; an unrecognized residue is not guessed to be a
/// protein residue.
fn classify_residue(name: &str) -> Option<ResidueClass> {
    let (canonical, modified) = match name.trim().to_ascii_uppercase().as_str() {
        "ALA" => ('A', false),
        "ARG" => ('R', false),
        "ASN" => ('N', false),
        "ASP" => ('D', false),
        "CYS" => ('C', false),
        "GLN" => ('Q', false),
        "GLU" => ('E', false),
        "GLY" => ('G', false),
        "HIS" => ('H', false),
        "ILE" => ('I', false),
        "LEU" => ('L', false),
        "LYS" => ('K', false),
        "MET" => ('M', false),
        "PHE" => ('F', false),
        "PRO" => ('P', false),
        "SER" => ('S', false),
        "THR" => ('T', false),
        "TRP" => ('W', false),
        "TYR" => ('Y', false),
        "VAL" => ('V', false),
        "SEC" => ('U', false),
        "PYL" => ('O', false),
        "UNK" => ('X', false),
        "MSE" => ('M', true),
        "HYP" => ('P', true),
        "SEP" => ('S', true),
        "TPO" => ('T', true),
        "PTR" => ('Y', true),
        "CSO" | "CME" | "CSD" | "CYM" | "CYX" => ('C', true),
        "KCX" | "MLY" | "LYN" => ('K', true),
        "HSD" | "HSE" | "HSP" => ('H', true),
        "ASH" => ('D', true),
        "GLH" => ('E', true),
        _ => return None,
    };
    Some(ResidueClass {
        canonical,
        modified,
    })
}

#[derive(Debug)]
struct ProjectionBuilder {
    projection: ProteinProjection,
    chain_indices: HashMap<ChainKey, usize>,
    residue_indices: Vec<HashMap<ResidueKey, usize>>,
}

type ChainKey = (i32, String);
type ResidueKey = (i64, Option<char>, String);

impl ProjectionBuilder {
    fn new() -> Self {
        Self {
            projection: ProteinProjection::default(),
            chain_indices: HashMap::new(),
            residue_indices: Vec::new(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn push(
        &mut self,
        model_number: i32,
        chain_id: &str,
        residue_name: &str,
        sequence_number: i64,
        insertion_code: Option<char>,
        atom: ProteinAtom,
    ) {
        let Some(class) = classify_residue(residue_name) else {
            return;
        };
        let chain_key = (model_number, chain_id.to_owned());
        let chain_index = *self.chain_indices.entry(chain_key).or_insert_with(|| {
            let index = self.projection.chains.len();
            self.projection.chains.push(ProteinChain {
                model_number,
                id: chain_id.to_owned(),
                residues: Vec::new(),
            });
            self.residue_indices.push(HashMap::new());
            index
        });
        let residue_key = (
            sequence_number,
            insertion_code,
            residue_name.trim().to_ascii_uppercase(),
        );
        let residue_index = *self.residue_indices[chain_index]
            .entry(residue_key)
            .or_insert_with(|| {
                let index = self.projection.chains[chain_index].residues.len();
                self.projection.chains[chain_index]
                    .residues
                    .push(ProteinResidue {
                        name: residue_name.trim().to_ascii_uppercase(),
                        sequence_number,
                        insertion_code,
                        canonical_one_letter_code: class.canonical,
                        modified: class.modified,
                        atoms: Vec::new(),
                    });
                index
            });
        self.projection.chains[chain_index].residues[residue_index]
            .atoms
            .push(atom);
    }

    fn finish(self) -> ProteinProjection {
        self.projection
    }
}

fn nonblank_char(value: &str) -> Option<char> {
    value.chars().next().filter(|ch| !ch.is_whitespace())
}

fn title_case_symbol(value: &str) -> String {
    let mut chars = value.trim().chars();
    match chars.next() {
        Some(first) => {
            first.to_ascii_uppercase().to_string() + &chars.as_str().to_ascii_lowercase()
        }
        None => String::new(),
    }
}

fn pdb_element(
    source: &chematic_3d::PdbAtom,
    record: &str,
) -> Result<Element, ProteinProjectionError> {
    if !source.element.trim().is_empty() {
        let symbol = title_case_symbol(&source.element);
        return Element::from_symbol(&symbol).ok_or(ProteinProjectionError::UnknownPdbElement {
            serial: source.serial,
            symbol,
        });
    }
    crate::infer_element(record, &source.res_name, &source.name).ok_or_else(|| {
        ProteinProjectionError::UnknownPdbElement {
            serial: source.serial,
            symbol: source.name.trim().to_owned(),
        }
    })
}

/// Parse PDB text into an amino-acid-only hierarchy.
pub fn project_protein_from_pdb(input: &str) -> Result<ProteinProjection, ProteinProjectionError> {
    project_protein_from_pdb_with_limits(input, &Default::default())
}

/// Parse PDB text with explicit resource limits.
pub fn project_protein_from_pdb_with_limits(
    input: &str,
    limits: &chematic_3d::PdbParseLimits,
) -> Result<ProteinProjection, ProteinProjectionError> {
    let parsed =
        chematic_3d::parse_pdb_atoms_strict(input, limits).map_err(ProteinProjectionError::Pdb)?;
    let mut parsed_atoms = parsed.into_iter();
    let mut builder = ProjectionBuilder::new();
    let mut model_number = 1_i32;

    for line in input.lines() {
        let record = line.get(0..6).unwrap_or("").trim();
        if record == "MODEL" {
            model_number = line
                .get(10..14)
                .unwrap_or("")
                .trim()
                .parse()
                .unwrap_or(model_number);
            continue;
        }
        if record != "ATOM" && record != "HETATM" {
            continue;
        }
        let Some(source) = parsed_atoms.next() else {
            break;
        };
        let element = pdb_element(&source, record)?;
        builder.push(
            model_number,
            &source.chain_id.to_string(),
            &source.res_name,
            i64::from(source.res_seq),
            nonblank_char(line.get(26..27).unwrap_or("")),
            ProteinAtom {
                serial: i64::from(source.serial),
                name: source.name.trim().to_owned(),
                element,
                position: [source.x, source.y, source.z],
                occupancy: source.occupancy,
                b_iso: source.temp_factor,
                alt_loc: nonblank_char(line.get(16..17).unwrap_or("")),
            },
        );
    }
    Ok(builder.finish())
}

/// Parse mmCIF text into the same amino-acid-only hierarchy as
/// [`project_protein_from_pdb`].
pub fn project_protein_from_mmcif(
    input: &str,
) -> Result<ProteinProjection, ProteinProjectionError> {
    project_protein_from_mmcif_with_limits(input, &MmcifParseLimits::default())
}

/// Parse mmCIF text with explicit resource limits.
pub fn project_protein_from_mmcif_with_limits(
    input: &str,
    limits: &MmcifParseLimits,
) -> Result<ProteinProjection, ProteinProjectionError> {
    let parsed = parse_mmcif_with_limits(input, limits).map_err(ProteinProjectionError::Mmcif)?;
    let mut builder = ProjectionBuilder::new();
    for source in parsed.atoms {
        builder.push(
            source.model_num,
            &source.chain_id,
            &source.res_name,
            source.res_seq,
            source.icode,
            ProteinAtom {
                serial: source.serial,
                name: source.atom_name,
                element: source.element,
                position: [source.x, source.y, source.z],
                occupancy: source.occupancy,
                b_iso: source.b_iso,
                alt_loc: source.alt_loc,
            },
        );
    }
    Ok(builder.finish())
}

#[cfg(test)]
mod tests;
