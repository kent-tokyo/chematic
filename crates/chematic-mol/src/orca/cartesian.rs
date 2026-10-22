//! Unit-aware Cartesian input projection; no coordinate conversion or job execution.
use super::{
    OrcaAtom, OrcaCoords, OrcaInputError, OrcaInputParseLimits, parse_orca_input_with_limits,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrcaCoordinateUnits {
    Angstrom,
    Bohr,
}
impl OrcaCoordinateUnits {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Angstrom => "angstrom",
            Self::Bohr => "bohr",
        }
    }
}

/// Refusal from the bounded parser or the explicit Cartesian/unit profile.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OrcaCartesianError {
    Input(OrcaInputError),
    ConflictingUnits,
    UnsupportedUnits(String),
    CoordsBlock,
    NoEmbeddedCartesian,
}
impl core::fmt::Display for OrcaCartesianError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input(error) => error.fmt(f),
            Self::ConflictingUnits => write!(
                f,
                "ORCA Cartesian input has conflicting Angs/Bohrs declarations"
            ),
            Self::UnsupportedUnits(word) => {
                write!(f, "unsupported ORCA coordinate unit keyword: {word}")
            }
            Self::CoordsBlock => write!(
                f,
                "ORCA Cartesian projection does not interpret %coords blocks"
            ),
            Self::NoEmbeddedCartesian => write!(
                f,
                "ORCA Cartesian projection requires one nonempty embedded xyz block"
            ),
        }
    }
}
impl std::error::Error for OrcaCartesianError {}
impl From<OrcaInputError> for OrcaCartesianError {
    fn from(error: OrcaInputError) -> Self {
        Self::Input(error)
    }
}

/// Input numeric values with explicit units and ordered state/freeze/extra fields.
/// Extra atom tokens are opaque. This is not a complete calculation/job model.
#[derive(Debug, Clone, PartialEq)]
pub struct OrcaCartesianGeometry {
    pub atoms: Vec<OrcaAtom>,
    pub charge: i32,
    pub multiplicity: u32,
    pub coordinate_units: OrcaCoordinateUnits,
    /// False when no Angs/Bohrs keyword was present (default: Angstrom).
    pub units_explicit: bool,
}
impl OrcaCartesianGeometry {
    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({
            "schema": "chematic.orca-cartesian.v1",
            "coordinate_units": self.coordinate_units.as_str(),
            "units_explicit": self.units_explicit,
            "charge": self.charge, "multiplicity": self.multiplicity,
            "elements": self.atoms.iter().map(|a| a.element.symbol()).collect::<Vec<_>>(),
            "coords": self.atoms.iter().map(|a| [a.x,a.y,a.z]).collect::<Vec<_>>(),
            "frozen": self.atoms.iter().map(|a| a.frozen).collect::<Vec<_>>(),
            "atom_extras": self.atoms.iter().map(|a| &a.extra).collect::<Vec<_>>(),
        })
    }
}

fn coordinate_units(
    keywords: &[String],
) -> Result<(OrcaCoordinateUnits, bool), OrcaCartesianError> {
    let mut declared = None;
    for word in keywords {
        let lower = word.to_ascii_lowercase();
        let units = match lower.as_str() {
            "angs" => OrcaCoordinateUnits::Angstrom,
            "bohrs" => OrcaCoordinateUnits::Bohr,
            _ if lower.starts_with("ang")
                || lower.starts_with("bohr")
                || lower.starts_with("units") =>
            {
                return Err(OrcaCartesianError::UnsupportedUnits(word.clone()));
            }
            _ => continue,
        };
        if declared.is_some_and(|old| old != units) {
            return Err(OrcaCartesianError::ConflictingUnits);
        }
        declared = Some(units);
    }
    Ok((
        declared.unwrap_or(OrcaCoordinateUnits::Angstrom),
        declared.is_some(),
    ))
}

pub fn parse_orca_cartesian(input: &str) -> Result<OrcaCartesianGeometry, OrcaCartesianError> {
    parse_orca_cartesian_with_limits(input, &OrcaInputParseLimits::default())
}

/// Read one explicit nonempty Cartesian geometry and its numeric coordinate unit.
/// Only exact Angs/Bohrs keywords are interpreted, case insensitively; absent
/// units default to Angstrom. Conflicting declarations and unit-like unsupported
/// tokens are refused. Any %coords block is refused because it may override
/// geometry/units. Internal/external geometries are not resolved or selected.
/// Values are not rescaled. Atom extras remain opaque; non-coordinate job
/// settings are omitted from this projection. No job execution is validated.
pub fn parse_orca_cartesian_with_limits(
    input: &str,
    limits: &OrcaInputParseLimits,
) -> Result<OrcaCartesianGeometry, OrcaCartesianError> {
    let parsed = parse_orca_input_with_limits(input, limits)?;
    if parsed.blocks.iter().any(|block| block.name == "coords") {
        return Err(OrcaCartesianError::CoordsBlock);
    }
    let (coordinate_units, units_explicit) = coordinate_units(&parsed.keywords)?;
    let Some(OrcaCoords::Xyz {
        charge,
        multiplicity,
        atoms,
    }) = parsed.coords
    else {
        return Err(OrcaCartesianError::NoEmbeddedCartesian);
    };
    Ok(OrcaCartesianGeometry {
        atoms,
        charge,
        multiplicity,
        coordinate_units,
        units_explicit,
    })
}
