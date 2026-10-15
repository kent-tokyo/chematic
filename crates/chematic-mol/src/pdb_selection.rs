//! Identity-preserving single-model PDB coordinate selection (Open Babel #2481).
//! Independently implemented from wwPDB 3.3; no residue or bond perception.
use std::collections::HashMap;
mod reader;
mod writer;

#[derive(Debug, Clone, Copy)]
pub struct PdbSelectionLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
    pub max_records: usize,
}
impl Default for PdbSelectionLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 1024,
            max_atoms: 99_999,
            max_records: 1_000_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbSelectionError {
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for PdbSelectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "PDB selection line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for PdbSelectionError {}
fn error(line: usize, detail: impl Into<String>) -> PdbSelectionError {
    PdbSelectionError {
        line,
        detail: detail.into(),
    }
}
fn reserve<T>(v: &mut Vec<T>, n: usize) -> Result<(), PdbSelectionError> {
    v.try_reserve(n)
        .map_err(|e| error(0, format!("selection allocation failed: {e}")))
}

#[derive(Debug, Clone, PartialEq)]
pub struct PdbSelectedAtom {
    pub source_index: usize,
    pub serial: u32,
    pub record_type: String,
    pub atom_name: String,
    pub alt_loc: char,
    pub residue_name: String,
    pub chain_id: char,
    pub residue_number: i32,
    pub insertion_code: char,
    pub position_angstrom: [f64; 3],
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbSelectionLoss {
    pub line: usize,
    pub record_name: String,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PdbSelectionReport {
    pub pdb: String,
    pub source_text: String,
    pub source_atom_count: usize,
    pub selected_source_indices: Vec<usize>,
    pub old_to_new: Vec<Option<usize>>,
    pub atoms: Vec<PdbSelectedAtom>,
    pub losses: Vec<PdbSelectionLoss>,
    pub retargeted_ter_lines: Vec<usize>,
    pub filtered_conect_lines: Vec<usize>,
    pub model_num: Option<u32>,
}
enum Kind {
    Atom(usize),
    Anisou(usize),
    Ter(u32),
    Conect(Vec<u32>),
    Model,
    EndModel,
    End,
}
struct Row<'a> {
    line: usize,
    text: &'a str,
    kind: Kind,
}
struct Document<'a> {
    rows: Vec<Row<'a>>,
    atom_lines: Vec<&'a str>,
    atoms: Vec<chematic_3d::pdb::PdbAtom>,
    by_serial: HashMap<u32, usize>,
    losses: Vec<PdbSelectionLoss>,
    model_num: Option<u32>,
}

/// Select zero-based source ATOM/HETATM rows. None keeps all, [] keeps none.
/// Output stays in source order regardless of index-list order; duplicates fail.
/// Atom/ANISOU line contents and serials are retained without inference. Unknown
/// records require explicit metadata-loss consent and remain in source_text.
pub fn select_pdb_atoms(
    text: &str,
    indices: Option<&[usize]>,
    allow_metadata_loss: bool,
) -> Result<PdbSelectionReport, PdbSelectionError> {
    select_pdb_atoms_with_limits(
        text,
        indices,
        allow_metadata_loss,
        &PdbSelectionLimits::default(),
    )
}
pub fn select_pdb_atoms_with_limits(
    text: &str,
    indices: Option<&[usize]>,
    allow_metadata_loss: bool,
    limits: &PdbSelectionLimits,
) -> Result<PdbSelectionReport, PdbSelectionError> {
    if indices.is_some_and(|v| {
        v.len()
            > limits
                .max_atoms
                .min(PdbSelectionLimits::default().max_atoms)
    }) {
        return Err(error(0, "selection count exceeds atom limit"));
    }
    let doc = reader::read(text, allow_metadata_loss, limits)?;
    let mut selected = Vec::new();
    reserve(&mut selected, doc.atoms.len())?;
    selected.resize(doc.atoms.len(), indices.is_none());
    if let Some(indices) = indices {
        for &i in indices {
            let entry = selected
                .get_mut(i)
                .ok_or_else(|| error(0, format!("atom index {i} outside source atom range")))?;
            if *entry {
                return Err(error(0, "duplicate selected atom index"));
            }
            *entry = true;
        }
    }
    writer::select(text, &doc, &selected)
}
