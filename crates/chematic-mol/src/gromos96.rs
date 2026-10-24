//! Checked GROMACS-compatible GROMOS96 blocks. Atom names are labels, never
//! elements. Length units are caller-owned context: the file has no unit tag.
use serde::{Deserialize, Serialize};
mod fields;
mod frames;
mod text;
use fields::{check_line, error};
pub use text::{
    parse_gromos96, parse_gromos96_with_limits, write_gromos96, write_gromos96_with_limits,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gromos96Error {
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for Gromos96Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "GROMOS96 line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for Gromos96Error {}
#[derive(Debug, Clone, Copy)]
pub struct Gromos96ParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
    pub max_frames: usize,
}
impl Default for Gromos96ParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 8192,
            max_atoms: 1_000_000,
            max_frames: 1024,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum Gromos96LengthUnit {
    #[default]
    Unspecified,
    Nm,
    Angstrom,
}
impl core::str::FromStr for Gromos96LengthUnit {
    type Err = Gromos96Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "unspecified" => Ok(Self::Unspecified),
            "nm" => Ok(Self::Nm),
            "angstrom" => Ok(Self::Angstrom),
            _ => Err(error(0, "expected length unit unspecified, nm or angstrom")),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Gromos96Block {
    pub kind: String,
    /// Blank/comment layout preceding this block, retained in original order.
    pub leading: Vec<String>,
    /// Authoritative current payload, including comments, padding and decimals.
    pub lines: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Gromos96Document {
    /// External interpretation only; write/read requires the same caller context.
    pub length_unit: Gromos96LengthUnit,
    pub blocks: Vec<Gromos96Block>,
    pub trailing: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Gromos96Identity {
    pub residue_number: String,
    pub residue_name: String,
    pub atom_name: String,
    pub atom_number: String,
}
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Gromos96Geometry {
    pub identities: Option<Vec<Gromos96Identity>>,
    pub positions_angstrom: Vec<[f64; 3]>,
    pub box_vectors_angstrom: Option<[[f64; 3]; 3]>,
}
impl Gromos96Document {
    pub fn validate(&self) -> Result<(), Gromos96Error> {
        self.validate_with_limits(&Default::default())
    }
    pub fn validate_with_limits(&self, limits: &Gromos96ParseLimits) -> Result<(), Gromos96Error> {
        frames::collect(self, limits).map(|_| ())
    }
    pub fn geometry_angstrom(&self) -> Result<Vec<Gromos96Geometry>, Gromos96Error> {
        self.geometry_angstrom_with_limits(&Default::default())
    }
    pub fn geometry_angstrom_with_limits(
        &self,
        limits: &Gromos96ParseLimits,
    ) -> Result<Vec<Gromos96Geometry>, Gromos96Error> {
        let frames = frames::collect(self, limits)?;
        let factor = match self.length_unit {
            Gromos96LengthUnit::Nm => 10.,
            Gromos96LengthUnit::Angstrom => 1.,
            Gromos96LengthUnit::Unspecified => {
                return Err(error(
                    0,
                    "length unit is unspecified; caller must declare nm or angstrom",
                ));
            }
        };
        frames
            .into_iter()
            .map(|frame| {
                let mut positions = frame.positions;
                for row in &mut positions {
                    for v in row {
                        *v = scaled(*v, factor)?;
                    }
                }
                let mut cell = frame.cell;
                if let Some(cell) = &mut cell {
                    for row in cell {
                        for v in row {
                            *v = scaled(*v, factor)?;
                        }
                    }
                }
                Ok(Gromos96Geometry {
                    identities: frame.identities,
                    positions_angstrom: positions,
                    box_vectors_angstrom: cell,
                })
            })
            .collect()
    }
}
fn scaled(value: f64, factor: f64) -> Result<f64, Gromos96Error> {
    let result = value * factor;
    if !result.is_finite() || (value != 0. && result == 0.) {
        return Err(error(
            0,
            "Angstrom projection overflows or underflows to zero",
        ));
    }
    Ok(result)
}
pub fn parse_gromos96_json(
    text: &str,
    limits: &Gromos96ParseLimits,
) -> Result<Gromos96Document, Gromos96Error> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "JSON byte limit exceeded"));
    }
    let value: crate::cjson_record::UniqueValue =
        serde_json::from_str(text).map_err(|e| error(0, e.to_string()))?;
    let document: Gromos96Document =
        serde_json::from_value(value.0).map_err(|e| error(0, e.to_string()))?;
    document.validate_with_limits(limits)?;
    Ok(document)
}
