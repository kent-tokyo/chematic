//! Editable formatted DL_POLY Classic HISTORY trajectories.
//! Source units stay external (FIELD/CONTROL); labels, charges and vectors are
//! never projected into inferred elements, formal charges or wrapped positions.
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

mod reader;
mod writer;
pub use reader::{parse_dlpoly_history, parse_dlpoly_history_with_limits};
pub use writer::write_dlpoly_history;

#[derive(Debug, Clone, Copy)]
pub struct DlpolyHistoryLimits {
    pub max_input_bytes: usize,
    pub max_atoms: usize,
    pub max_frames: usize,
    /// Sum of atom records across all frames.
    pub max_atom_records: usize,
}
impl Default for DlpolyHistoryLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_atoms: 100_000,
            max_frames: 10_000,
            max_atom_records: 1_000_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DlpolyHistoryError {
    /// One-based source line, or zero for record validation/global limits.
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for DlpolyHistoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "DL-POLY HISTORY line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for DlpolyHistoryError {}
fn error(line: usize, detail: impl Into<String>) -> DlpolyHistoryError {
    DlpolyHistoryError {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DlpolyHistoryAtom {
    pub name: String,
    pub index: usize,
    pub mass: f64,
    /// Source atomic charge, not molecular formal charge.
    pub charge: f64,
    pub position: [f64; 3],
    pub velocity: Option<[f64; 3]>,
    pub force: Option<[f64; 3]>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DlpolyHistoryFrame {
    pub step: u64,
    /// Integration timestep in source units. Not elapsed time.
    pub timestep: f64,
    /// Rows are the source a,b,c vectors; handedness/orientation are retained.
    pub cell_vectors: Option<[[f64; 3]; 3]>,
    pub atoms: Vec<DlpolyHistoryAtom>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DlpolyHistory {
    pub title: String,
    /// 0 positions, 1 positions/velocities, 2 positions/velocities/forces.
    pub keytrj: u8,
    /// Bounded profile: 0 nonperiodic, 1 cubic, 2 orthorhombic, 3 parallelepiped.
    pub imcon: u8,
    pub atom_count: usize,
    pub frames: Vec<DlpolyHistoryFrame>,
}
impl DlpolyHistory {
    pub fn validate(&self) -> Result<(), DlpolyHistoryError> {
        self.validate_with_limits(&DlpolyHistoryLimits::default())
    }
    pub(crate) fn validate_with_limits(
        &self,
        limits: &DlpolyHistoryLimits,
    ) -> Result<(), DlpolyHistoryError> {
        if self.title.len() > 80
            || !self.title.is_ascii()
            || self.title.chars().any(char::is_control)
        {
            return Err(error(
                0,
                "title must be printable ASCII of at most 80 bytes",
            ));
        }
        if self.keytrj > 2
            || self.imcon > 3
            || self.atom_count == 0
            || self.atom_count > limits.max_atoms
        {
            return Err(error(
                0,
                "expected keytrj 0..2, imcon 0..3 and bounded nonzero atom count",
            ));
        }
        if self.frames.is_empty()
            || self.frames.len() > limits.max_frames
            || self
                .frames
                .len()
                .checked_mul(self.atom_count)
                .is_none_or(|n| n > limits.max_atom_records)
        {
            return Err(error(0, "frame/total atom-record count outside limits"));
        }
        for frame in &self.frames {
            validate_frame(frame, self)?;
        }
        Ok(())
    }
}
fn validate_frame(f: &DlpolyHistoryFrame, g: &DlpolyHistory) -> Result<(), DlpolyHistoryError> {
    if f.step > 9_999_999_999
        || f.atoms.len() != g.atom_count
        || !f.timestep.is_finite()
        || f.timestep <= 0.0
    {
        return Err(error(0, "invalid frame step, timestep or atom count"));
    }
    real(f.timestep)?;
    match (g.imcon, &f.cell_vectors) {
        (0, None) => {}
        (1..=3, Some(cell)) => validate_cell(cell, g.imcon)?,
        _ => return Err(error(0, "cell presence disagrees with imcon")),
    }
    let mut ids = BTreeSet::new();
    for a in &f.atoms {
        if a.index == 0
            || a.index > g.atom_count
            || !ids.insert(a.index)
            || a.name.is_empty()
            || a.name.len() > 8
            || !a.name.bytes().all(|b| b.is_ascii_graphic())
        {
            return Err(error(
                0,
                "names need 1..8 ASCII nonspace bytes; IDs must uniquely cover 1..N",
            ));
        }
        if a.mass < 0.0
            || a.velocity.is_some() != (g.keytrj >= 1)
            || a.force.is_some() != (g.keytrj == 2)
        {
            return Err(error(
                0,
                "negative mass or velocity/force presence disagrees with keytrj",
            ));
        }
        real(a.mass)?;
        real(a.charge)?;
        for v in std::iter::once(&a.position)
            .chain(a.velocity.as_ref())
            .chain(a.force.as_ref())
        {
            vector(v)?;
        }
    }
    Ok(())
}
fn vector(v: &[f64; 3]) -> Result<(), DlpolyHistoryError> {
    for &n in v {
        real(n)?;
    }
    Ok(())
}
fn validate_cell(cell: &[[f64; 3]; 3], imcon: u8) -> Result<(), DlpolyHistoryError> {
    let mut normalized = *cell;
    for row in &mut normalized {
        vector(row)?;
        let scale = row.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return Err(error(0, "zero cell vector"));
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
                return Err(error(0, "imcon 1/2 requires positive axis-aligned vectors"));
            }
        }
        if imcon == 1 && (cell[0][0] != cell[1][1] || cell[0][0] != cell[2][2]) {
            return Err(error(0, "imcon 1 requires a cube"));
        }
    }
    Ok(())
}
/// Exact, explicit-decimal real text fitting a twelve-column Fortran input
/// field. Never silently rounds edited values to a lower output precision.
fn real(n: f64) -> Result<String, DlpolyHistoryError> {
    if !n.is_finite() {
        return Err(error(0, "nonfinite scalar/vector value"));
    }
    let mut choices = Vec::new();
    for s in [n.to_string(), format!("{n:e}")] {
        let text = if let Some((m, e)) = s.split_once('e') {
            let m = if m.contains('.') {
                m.to_owned()
            } else {
                format!("{m}.0")
            };
            let e: i32 = e.parse().map_err(|_| error(0, "invalid exponent"))?;
            format!("{m}E{e:+03}")
        } else if s.contains('.') {
            s
        } else {
            format!("{s}.0")
        };
        if text.len() <= 12 && text.parse::<f64>().ok() == Some(n) {
            choices.push(text);
        }
    }
    choices
        .into_iter()
        .min_by_key(String::len)
        .map(|s| format!("{s:>12}"))
        .ok_or_else(|| error(0, "number cannot fit exactly in a twelve-column real field"))
}
pub fn dlpoly_history_from_json(text: &str) -> Result<DlpolyHistory, DlpolyHistoryError> {
    if text.len() > 64 << 20 {
        return Err(error(0, "JSON byte limit exceeded"));
    }
    let g: DlpolyHistory = serde_json::from_str(text).map_err(|e| error(0, e.to_string()))?;
    g.validate()?;
    Ok(g)
}
