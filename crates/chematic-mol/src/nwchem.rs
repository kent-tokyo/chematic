//! NWChem Cartesian GEOMETRY blocks, independently implemented from the manual.
//! Tags and nuclear charges remain source values, never graph/formal charges.
use chematic_core::Element;
use serde::{Deserialize, Serialize};

mod reader;
mod writer;
pub use reader::{parse_nwchem_geometry, parse_nwchem_geometry_with_limits};
pub use writer::write_nwchem_geometry;

#[derive(Debug, Clone, Copy)]
pub struct NwchemLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_centers: usize,
}
impl Default for NwchemLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 << 20,
            max_line_len: 1023,
            max_centers: 100_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NwchemError {
    /// One-based source line, or zero for record/global validation.
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for NwchemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "NWChem geometry line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for NwchemError {}
fn error(line: usize, detail: impl Into<String>) -> NwchemError {
    NwchemError {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NwchemNucleus {
    Point,
    Finite,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NwchemPrint {
    Normal,
    Xyz,
    None,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NwchemCenter {
    /// Case-sensitive basis/ECP/symmetry tag, including dummy-center tags.
    pub tag: String,
    /// In the declared geometry unit; absence means NWChem's angstrom default.
    pub position: [f64; 3],
    /// Atomic units, independent of the coordinate unit.
    pub velocity_au: Option<[f64; 3]>,
    /// Nuclear/center charge, NOT a partial charge or a molecular formal charge.
    pub nuclear_charge: Option<f64>,
    pub mass_amu: Option<f64>,
    pub nucleus: Option<NwchemNucleus>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NwchemGeometry {
    pub name: Option<String>,
    /// Source unit alias (lowercase), or None for an omitted declaration.
    pub units: Option<String>,
    pub angstrom_to_au: Option<f64>,
    pub center: Option<bool>,
    pub autosym: Option<bool>,
    pub symmetry_tolerance: Option<f64>,
    pub autoz: Option<bool>,
    pub print: Option<NwchemPrint>,
    pub bqbq: bool,
    pub nucleus: Option<NwchemNucleus>,
    pub centers: Vec<NwchemCenter>,
    /// Ordered comment text after '#'. Canonical output places comments before
    /// the block; inline/source comment positions are not represented.
    pub comments: Vec<String>,
}
const UNITS: &[&str] = &[
    "angstroms",
    "an",
    "au",
    "atomic",
    "bohr",
    "nm",
    "nanometers",
    "pm",
    "picometers",
];
const OPTIONS: &[&str] = &[
    "units",
    "angstrom_to_au",
    "ang2au",
    "center",
    "nocenter",
    "autosym",
    "noautosym",
    "autoz",
    "noautoz",
    "print",
    "noprint",
    "bqbq",
    "nuc",
    "nucl",
    "nucleus",
];
const UNSUPPORTED: &[&str] = &[
    "adjust", "symmetry", "load", "system", "zmatrix", "zmt", "zmat", "zcoord",
];

impl NwchemGeometry {
    pub fn validate(&self) -> Result<(), NwchemError> {
        self.validate_with_limits(&NwchemLimits::default())
    }
    pub(crate) fn validate_with_limits(&self, limits: &NwchemLimits) -> Result<(), NwchemError> {
        if self.centers.is_empty() || self.centers.len() > limits.max_centers {
            return Err(error(0, "center count outside limits"));
        }
        if let Some(name) = &self.name {
            token(name, 64)?;
            if OPTIONS.contains(&name.to_ascii_lowercase().as_str())
                || UNSUPPORTED.contains(&name.to_ascii_lowercase().as_str())
            {
                return Err(error(0, "geometry name collides with an option"));
            }
        }
        if self
            .units
            .as_ref()
            .is_some_and(|u| !UNITS.contains(&u.as_str()))
        {
            return Err(error(0, "unsupported geometry unit"));
        }
        positive(self.angstrom_to_au, "conversion factor")?;
        positive(self.symmetry_tolerance, "symmetry tolerance")?;
        if self.symmetry_tolerance.is_some() && self.autosym != Some(true) {
            return Err(error(0, "symmetry tolerance requires explicit autosym"));
        }
        for c in &self.centers {
            token(&c.tag, 16)?;
            validate_tag(&c.tag)?;
            if c.position
                .iter()
                .chain(c.velocity_au.iter().flatten())
                .chain(c.nuclear_charge.iter())
                .any(|v| !v.is_finite())
            {
                return Err(error(0, "nonfinite position, velocity or nuclear charge"));
            }
            positive(c.mass_amu, "mass")?;
        }
        for c in &self.comments {
            if !c.is_ascii()
                || c.bytes().any(|b| b.is_ascii_control())
                || c.len() + 1 > limits.max_line_len
            {
                return Err(error(
                    0,
                    "comment exceeds line limit or contains control/non-ASCII characters",
                ));
            }
        }
        Ok(())
    }
}
fn token(s: &str, max: usize) -> Result<(), NwchemError> {
    if s.is_empty()
        || s.len() > max
        || !s
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_.-".contains(&b))
    {
        return Err(error(0, "expected bounded unquoted ASCII name/tag"));
    }
    Ok(())
}
fn validate_tag(tag: &str) -> Result<(), NwchemError> {
    let lower = tag.to_ascii_lowercase();
    if UNSUPPORTED.contains(&lower.as_str())
        || ["include", "eof", "end", "geometry"].contains(&lower.as_str())
    {
        return Err(error(
            0,
            "reserved directive cannot be used as a center tag",
        ));
    }
    let valid = lower.starts_with('x')
        || lower.starts_with("bq")
        || (1..=118).any(|z| {
            Element::from_atomic_number(z)
                .is_some_and(|e| lower.starts_with(&e.symbol().to_ascii_lowercase()))
        });
    if !valid {
        return Err(error(0, "tag has no supported element/dummy-center prefix"));
    }
    Ok(())
}
fn positive(n: Option<f64>, label: &str) -> Result<(), NwchemError> {
    if n.is_some_and(|v| !v.is_finite() || v <= 0.0) {
        return Err(error(0, format!("{label} must be positive and finite")));
    }
    Ok(())
}
pub fn nwchem_geometry_from_json(json: &str) -> Result<NwchemGeometry, NwchemError> {
    if json.len() > 64 << 20 {
        return Err(error(0, "geometry JSON exceeds byte limit"));
    }
    let record: NwchemGeometry = serde_json::from_str(json).map_err(|e| error(0, e.to_string()))?;
    record.validate()?;
    Ok(record)
}
