//! Molden geometry/vibration profile, from the public Molden format description.
//! Explicit source units and mode IDs are retained. Electronic/trajectory
//! sections are refused, never read past and silently discarded.
use chematic_core::Element;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

mod reader;
mod writer;
pub use reader::{parse_molden, parse_molden_with_limits};
pub use writer::write_molden;

#[derive(Debug, Clone, Copy)]
pub struct MoldenLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
    pub max_modes: usize,
    /// Total displacement triplets across all modes, not a per-mode limit.
    pub max_displacements: usize,
}
impl Default for MoldenLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 65_536,
            max_atoms: 100_000,
            max_modes: 100_000,
            max_displacements: 1_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MoldenError {
    /// One-based source line, or zero for record validation/global limits.
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for MoldenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Molden line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for MoldenError {}
fn error(line: usize, detail: impl Into<String>) -> MoldenError {
    MoldenError {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MoldenLengthUnit {
    Angstrom,
    Bohr,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MoldenAtom {
    pub symbol: String,
    pub index: usize,
    pub atomic_number: u8,
    /// In the explicit enclosing [Atoms] unit; never automatically converted.
    pub position: [f64; 3],
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MoldenAtoms {
    pub unit: MoldenLengthUnit,
    pub atoms: Vec<MoldenAtom>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MoldenVibrationAtom {
    pub symbol: String,
    pub position_bohr: [f64; 3],
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MoldenNormalMode {
    /// One-based index into the frequency list, independent of block order.
    pub index: usize,
    pub displacements_bohr: Vec<[f64; 3]>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MoldenIntensity {
    pub ir: f64,
    pub raman: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MoldenRecord {
    pub atoms: Option<MoldenAtoms>,
    /// Source frequency values, including zero and negative values; no scaling.
    pub frequencies: Option<Vec<f64>>,
    /// [FR-COORD] has its own Bohr unit even when [Atoms] uses angstroms.
    pub fr_coords: Option<Vec<MoldenVibrationAtom>>,
    pub normal_modes: Option<Vec<MoldenNormalMode>>,
    /// Source [INT] values; absent Raman data stays absent, never zero-filled.
    pub intensities: Option<Vec<MoldenIntensity>>,
}

impl MoldenRecord {
    pub fn validate(&self) -> Result<(), MoldenError> {
        self.validate_with_limits(&MoldenLimits::default())
    }
    pub(crate) fn validate_with_limits(&self, limits: &MoldenLimits) -> Result<(), MoldenError> {
        if self.atoms.is_none() && self.fr_coords.is_none() {
            return Err(error(0, "missing geometry"));
        }
        if let Some(atoms) = &self.atoms {
            validate_atoms(atoms, limits)?;
        }
        if let Some(coords) = &self.fr_coords {
            if coords.is_empty() || coords.len() > limits.max_atoms {
                return Err(error(0, "FR-COORD atom count outside limits"));
            }
            for atom in coords {
                canonical_symbol(&atom.symbol)?;
                finite_vector(&atom.position_bohr)?;
            }
        }
        if let (Some(atoms), Some(fr)) = (&self.atoms, &self.fr_coords)
            && (atoms.atoms.len() != fr.len()
                || atoms
                    .atoms
                    .iter()
                    .zip(fr)
                    .any(|(a, b)| a.symbol != b.symbol))
        {
            return Err(error(
                0,
                "Atoms and FR-COORD must have the same source-order element sequence",
            ));
        }
        validate_vibrations(self, limits)
    }
}

fn canonical_symbol(symbol: &str) -> Result<Element, MoldenError> {
    Element::from_symbol(symbol)
        .filter(|e| e.symbol() == symbol)
        .ok_or_else(|| {
            error(
                0,
                "expected canonical element symbol; dummy/ghost atoms are unsupported",
            )
        })
}
fn finite_vector(v: &[f64; 3]) -> Result<(), MoldenError> {
    if v.iter().any(|n| !n.is_finite()) {
        return Err(error(0, "non-finite coordinate/displacement"));
    }
    Ok(())
}
fn validate_atoms(atoms: &MoldenAtoms, limits: &MoldenLimits) -> Result<(), MoldenError> {
    if atoms.atoms.is_empty() || atoms.atoms.len() > limits.max_atoms {
        return Err(error(0, "Atoms count outside limits"));
    }
    let mut ids = BTreeSet::new();
    for atom in &atoms.atoms {
        if canonical_symbol(&atom.symbol)?.atomic_number() != atom.atomic_number {
            return Err(error(0, "element symbol disagrees with atomic number"));
        }
        if atom.index == 0 || atom.index > atoms.atoms.len() || !ids.insert(atom.index) {
            return Err(error(
                0,
                "atom IDs must uniquely cover 1..N; source order is retained",
            ));
        }
        finite_vector(&atom.position)?;
    }
    Ok(())
}
fn validate_vibrations(record: &MoldenRecord, limits: &MoldenLimits) -> Result<(), MoldenError> {
    let any = record.frequencies.is_some()
        || record.normal_modes.is_some()
        || record.intensities.is_some();
    if !any {
        return Ok(());
    }
    let frequencies = record
        .frequencies
        .as_ref()
        .ok_or_else(|| error(0, "vibrational data requires FREQ"))?;
    let coords = record
        .fr_coords
        .as_ref()
        .ok_or_else(|| error(0, "vibrational data requires FR-COORD"))?;
    if frequencies.is_empty()
        || frequencies.len() > limits.max_modes
        || frequencies.iter().any(|n| !n.is_finite())
    {
        return Err(error(
            0,
            "frequency count outside limits or non-finite frequency",
        ));
    }
    if let Some(modes) = &record.normal_modes {
        validate_modes(modes, frequencies.len(), coords.len(), limits)?;
    }
    if let Some(intensities) = &record.intensities
        && (intensities.len() != frequencies.len()
            || intensities
                .iter()
                .any(|i| !i.ir.is_finite() || i.raman.is_some_and(|v| !v.is_finite())))
    {
        return Err(error(
            0,
            "INT count must match FREQ and all intensities must be finite",
        ));
    }
    Ok(())
}
fn validate_modes(
    modes: &[MoldenNormalMode],
    count: usize,
    atoms: usize,
    limits: &MoldenLimits,
) -> Result<(), MoldenError> {
    if modes.len() != count
        || atoms
            .checked_mul(count)
            .is_none_or(|n| n > limits.max_displacements)
    {
        return Err(error(
            0,
            "mode count must match FREQ and total displacements must fit limit",
        ));
    }
    let mut ids = BTreeSet::new();
    for mode in modes {
        if mode.index == 0
            || mode.index > count
            || !ids.insert(mode.index)
            || mode.displacements_bohr.len() != atoms
        {
            return Err(error(
                0,
                "mode IDs must cover 1..frequency count; every mode needs one triplet per FR-COORD atom",
            ));
        }
        for v in &mode.displacements_bohr {
            finite_vector(v)?;
        }
    }
    Ok(())
}

pub fn molden_record_from_json(json: &str) -> Result<MoldenRecord, MoldenError> {
    if json.len() > MoldenLimits::default().max_input_bytes {
        return Err(error(0, "JSON exceeds byte limit"));
    }
    let record: MoldenRecord = serde_json::from_str(json).map_err(|e| error(0, e.to_string()))?;
    record.validate()?;
    Ok(record)
}
