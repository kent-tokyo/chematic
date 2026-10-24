//! DL-POLY CONFIG records, implemented from the public DL_POLY manual.
//! Numbers retain source units: CONFIG alone cannot distinguish internal units
//! from DPD units. Names remain opaque; no chemistry or force field is inferred.

use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

mod reader;
pub use reader::{parse_dlpoly_config, parse_dlpoly_config_with_limits};

#[derive(Debug, Clone, Copy)]
pub struct DlpolyConfigLimits {
    pub max_input_bytes: usize,
    pub max_atoms: usize,
}
impl Default for DlpolyConfigLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_atoms: 100_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DlpolyError {
    /// One-based source line, or zero for record validation/global limits.
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for DlpolyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "DL-POLY CONFIG line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for DlpolyError {}
fn error(line: usize, detail: impl Into<String>) -> DlpolyError {
    DlpolyError {
        line,
        detail: detail.into(),
    }
}

/// All vector components use the file's units, without conversion or wrapping.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DlpolyAtom {
    pub name: String,
    /// Explicit ID, or None for the Classic name-only variant.
    pub index: Option<usize>,
    /// Optional atom-header suffix words; never interpreted as masses/charges.
    pub header_extra: Vec<String>,
    pub position: [f64; 3],
    pub velocity: Option<[f64; 3]>,
    pub force: Option<[f64; 3]>,
}

/// Editable single CONFIG/REVCON record. FIELD/CONTROL context is external.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DlpolyConfig {
    pub title: String,
    /// 0: position; 1: position/velocity; 2: position/velocity/force.
    pub levcfg: u8,
    /// Supported profile: 0, 1, 2, 3. Other boundary geometries are refused.
    pub imcon: u8,
    pub declared_atom_count: Option<usize>,
    /// Rows are the a, b, c Cartesian vectors, retaining orientation/handedness.
    pub cell_vectors: Option<[[f64; 3]; 3]>,
    pub atoms: Vec<DlpolyAtom>,
}

impl DlpolyConfig {
    pub fn validate(&self) -> Result<(), DlpolyError> {
        self.validate_with_limits(&DlpolyConfigLimits::default())
    }
    pub(crate) fn validate_with_limits(
        &self,
        limits: &DlpolyConfigLimits,
    ) -> Result<(), DlpolyError> {
        if !self.title.is_ascii()
            || self.title.len() > 72
            || self.title.chars().any(char::is_control)
        {
            return Err(error(
                0,
                "title must be ASCII, at most 72 bytes, without controls",
            ));
        }
        if self.levcfg > 2 || self.imcon > 3 {
            return Err(error(0, "supported levcfg is 0..2 and imcon is 0..3"));
        }
        if self.atoms.is_empty()
            || self.atoms.len() > limits.max_atoms
            || self
                .declared_atom_count
                .is_some_and(|n| n != self.atoms.len())
        {
            return Err(error(
                0,
                "atom count is empty, exceeds limit, or disagrees with declaration",
            ));
        }
        match (self.imcon, &self.cell_vectors) {
            (0, None) => {}
            (1..=3, Some(cell)) => validate_cell(cell, self.imcon)?,
            _ => {
                return Err(error(
                    0,
                    "cell vectors must be present exactly when imcon is nonzero",
                ));
            }
        }
        let indexed = self.atoms[0].index.is_some();
        let mut seen = BTreeSet::new();
        for atom in &self.atoms {
            validate_atom(atom, self.levcfg)?;
            if atom.index.is_some() != indexed {
                return Err(error(
                    0,
                    "mixed indexed and name-only atom headers are unsupported",
                ));
            }
            if let Some(id) = atom.index
                && (id == 0 || id > self.atoms.len() || !seen.insert(id))
            {
                return Err(error(
                    0,
                    "IDs must be unique and cover 1..N; source order is retained",
                ));
            }
        }
        Ok(())
    }
}

fn validate_atom(atom: &DlpolyAtom, level: u8) -> Result<(), DlpolyError> {
    if atom.name.is_empty()
        || atom.name.len() > 8
        || !atom.name.bytes().all(|b| b.is_ascii_graphic())
        || atom.name.starts_with(['#', '!'])
    {
        return Err(error(
            0,
            "atom names need 1..8 printable ASCII non-whitespace characters, excluding comment prefixes",
        ));
    }
    if atom.velocity.is_some() != (level >= 1) || atom.force.is_some() != (level == 2) {
        return Err(error(0, "velocity/force presence disagrees with levcfg"));
    }
    if (!atom.header_extra.is_empty() && atom.index.is_none())
        || atom
            .header_extra
            .iter()
            .any(|s| s.is_empty() || s.len() > 40 || !s.bytes().all(|b| b.is_ascii_graphic()))
    {
        return Err(error(
            0,
            "opaque atom-header suffix needs an index and nonempty ASCII words of at most 40 bytes",
        ));
    }
    for vector in std::iter::once(&atom.position)
        .chain(atom.velocity.as_ref())
        .chain(atom.force.as_ref())
    {
        finite_vector(vector)?;
    }
    Ok(())
}

fn finite_vector(v: &[f64; 3]) -> Result<(), DlpolyError> {
    if v.iter().any(|x| !x.is_finite()) {
        return Err(error(0, "non-finite vector component"));
    }
    Ok(())
}

fn validate_cell(cell: &[[f64; 3]; 3], imcon: u8) -> Result<(), DlpolyError> {
    let mut normalized = *cell;
    for row in &mut normalized {
        finite_vector(row)?;
        let scale = row.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return Err(error(0, "zero-length cell vector"));
        }
        for v in row {
            *v /= scale;
        }
    }
    let [a, b, c] = normalized;
    let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);
    if det.abs() <= 1e-12 {
        return Err(error(0, "singular or near-singular cell"));
    }
    if imcon <= 2 {
        for (i, row) in cell.iter().enumerate() {
            if row[i] <= 0.0 || row.iter().enumerate().any(|(j, v)| i != j && *v != 0.0) {
                return Err(error(
                    0,
                    "imcon 1/2 profile requires positive axis-aligned cell vectors",
                ));
            }
        }
        if imcon == 1 && (cell[0][0] != cell[1][1] || cell[0][0] != cell[2][2]) {
            return Err(error(0, "imcon 1 requires equal cubic side lengths"));
        }
    }
    Ok(())
}

/// Strict typed JSON: unknown fields and invalid current data are refused.
pub fn dlpoly_config_from_json(text: &str) -> Result<DlpolyConfig, DlpolyError> {
    if text.len() > DlpolyConfigLimits::default().max_input_bytes {
        return Err(error(0, "JSON exceeds byte limit"));
    }
    let record: DlpolyConfig = serde_json::from_str(text).map_err(|e| error(0, e.to_string()))?;
    record.validate()?;
    Ok(record)
}

/// Write current data, keeping every declared vector and optional ID/count.
/// Refuses numbers whose shortest round-trip text cannot fit a 72-byte record.
pub fn write_dlpoly_config(record: &DlpolyConfig) -> Result<String, DlpolyError> {
    record.validate()?;
    let mut out = format!("{}\n{} {}", record.title, record.levcfg, record.imcon);
    if let Some(n) = record.declared_atom_count {
        out.push_str(&format!(" {n}"));
    }
    out.push('\n');
    if let Some(cell) = &record.cell_vectors {
        for v in cell {
            write_vector(&mut out, v)?;
        }
    }
    for atom in &record.atoms {
        write_atom_header(&mut out, atom)?;
        out.push('\n');
        for v in std::iter::once(&atom.position)
            .chain(atom.velocity.as_ref())
            .chain(atom.force.as_ref())
        {
            write_vector(&mut out, v)?;
        }
        if out.len() > DlpolyConfigLimits::default().max_input_bytes {
            return Err(error(0, "output exceeds byte limit"));
        }
    }
    Ok(out)
}

fn write_atom_header(out: &mut String, atom: &DlpolyAtom) -> Result<(), DlpolyError> {
    let start = out.len();
    out.push_str(&atom.name);
    if let Some(id) = atom.index {
        out.push_str(&format!(" {id}"));
    }
    for word in &atom.header_extra {
        out.push(' ');
        out.push_str(word);
    }
    if out.len() - start > 72 {
        return Err(error(0, "atom header exceeds 72-byte record limit"));
    }
    Ok(())
}

fn write_vector(out: &mut String, v: &[f64; 3]) -> Result<(), DlpolyError> {
    let line = v
        .map(|n| {
            let plain = n.to_string();
            let scientific = format!("{n:e}");
            if scientific.len() < plain.len() {
                scientific
            } else {
                plain
            }
        })
        .join(" ");
    if line.len() > 72 {
        return Err(error(
            0,
            "numeric vector cannot fit the 72-byte CONFIG record without precision loss",
        ));
    }
    out.push_str(&line);
    out.push('\n');
    Ok(())
}
