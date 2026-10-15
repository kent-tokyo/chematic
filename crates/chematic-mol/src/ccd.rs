//! Checked single-component CCD reader, independently implemented from the
//! wwPDB chem_comp_atom/chem_comp_bond dictionary (Open Babel #1516).
//! This is a molecular-template reader, separate from the atom_site adapter.
use std::collections::{BTreeMap, HashMap, HashSet};

use chematic_core::{Atom, BondOrder, Element, Molecule, MoleculeBuilder};
use serde::Serialize;

mod fields;
mod reader;
mod tokens;
pub use reader::{parse_ccd_component, parse_ccd_component_with_limits};

/// Bounds for one CCD data block. The whole dictionary is not a single component.
#[derive(Debug, Clone, Copy)]
pub struct CcdParseLimits {
    pub max_input_bytes: usize,
    pub max_tokens: usize,
    pub max_atoms: usize,
    pub max_bonds: usize,
}
impl Default for CcdParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 8 << 20,
            max_tokens: 500_000,
            max_atoms: 10_000,
            max_bonds: 40_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CcdError {
    ResourceLimit {
        resource: &'static str,
        limit: usize,
    },
    Invalid(String),
    StereoLoss,
    UnknownCharge {
        atom_id: String,
    },
}
impl std::fmt::Display for CcdError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ResourceLimit { resource, limit } => {
                write!(f, "CCD {resource} exceeds limit {limit}")
            }
            Self::Invalid(message) => write!(f, "invalid CCD: {message}"),
            Self::UnknownCharge { atom_id } => write!(
                f,
                "CCD formal charge is unknown at atom {atom_id}; graph projection refuses charge imputation"
            ),
            Self::StereoLoss => write!(
                f,
                "CCD stereochemical annotations require allow_stereo_loss=true for an unassigned graph projection"
            ),
        }
    }
}
impl std::error::Error for CcdError {}
fn invalid(message: impl Into<String>) -> CcdError {
    CcdError::Invalid(message.into())
}

type Row = BTreeMap<String, String>;

/// One explicit template atom; no implicit hydrogens are inferred.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CcdAtom {
    pub atom_id: String,
    pub element: String,
    pub isotope: Option<u16>,
    pub charge: Option<i8>,
    pub aromatic: bool,
    pub stereo_config: Option<String>,
    pub model_coordinates: Option<[f64; 3]>,
    pub ideal_coordinates: Option<[f64; 3]>,
    /// Unmodeled atom fields, keyed by their complete lower-case CIF tags.
    pub extra_fields: BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CcdBond {
    pub atom_id_1: String,
    pub atom_id_2: String,
    pub value_order: String,
    pub aromatic: bool,
    pub stereo_config: Option<String>,
    pub extra_fields: BTreeMap<String, String>,
}

/// Validated template in source row order. Non-atom/bond categories are listed,
/// not interpreted. Retains both coordinate families and stereo annotations.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CcdComponent {
    pub component_id: String,
    pub atoms: Vec<CcdAtom>,
    pub bonds: Vec<CcdBond>,
    pub unhandled_tags: Vec<String>,
}
impl CcdComponent {
    /// Project the explicit atom/bond template. R/S and E/Z are not assigned
    /// to graph parity; refusing is the default policy of language bindings.
    /// With consent they remain on this source record, but not the graph.
    pub fn to_molecule(&self, allow_stereo_loss: bool) -> Result<Molecule, CcdError> {
        if !allow_stereo_loss && self.has_unassigned_stereo() {
            return Err(CcdError::StereoLoss);
        }
        let aromatic_ids: HashSet<_> = self
            .bonds
            .iter()
            .filter(|b| b.aromatic)
            .flat_map(|b| [b.atom_id_1.as_str(), b.atom_id_2.as_str()])
            .collect();
        let mut builder = MoleculeBuilder::new();
        let mut ids = HashMap::new();
        for source in &self.atoms {
            let element =
                Element::from_symbol(&source.element).ok_or_else(|| invalid("unknown element"))?;
            let mut atom = Atom::new(element);
            atom.hydrogen_count = Some(0);
            atom.charge = source.charge.ok_or_else(|| CcdError::UnknownCharge {
                atom_id: source.atom_id.clone(),
            })?;
            atom.isotope = source.isotope;
            atom.aromatic = source.aromatic || aromatic_ids.contains(source.atom_id.as_str());
            let idx = builder.add_atom(atom);
            if ids.insert(source.atom_id.as_str(), idx).is_some() {
                return Err(invalid("duplicate atom ID"));
            }
        }
        for source in &self.bonds {
            let a = *ids
                .get(source.atom_id_1.as_str())
                .ok_or_else(|| invalid("unknown bond endpoint"))?;
            let b = *ids
                .get(source.atom_id_2.as_str())
                .ok_or_else(|| invalid("unknown bond endpoint"))?;
            let order = if source.aromatic {
                BondOrder::Aromatic
            } else {
                bond_order(&source.value_order)?
            };
            builder
                .add_bond(a, b, order)
                .map_err(|e| invalid(e.to_string()))?;
        }
        Ok(builder.build())
    }
    pub fn has_unassigned_stereo(&self) -> bool {
        self.atoms
            .iter()
            .filter_map(|a| a.stereo_config.as_deref())
            .any(|s| s != "N")
            || self
                .bonds
                .iter()
                .filter_map(|b| b.stereo_config.as_deref())
                .any(|s| s != "N")
    }
}
fn bond_order(value: &str) -> Result<BondOrder, CcdError> {
    match value {
        "SING" => Ok(BondOrder::Single),
        "DOUB" => Ok(BondOrder::Double),
        "TRIP" => Ok(BondOrder::Triple),
        "QUAD" => Ok(BondOrder::Quadruple),
        "AROM" => Ok(BondOrder::Aromatic),
        _ => Err(invalid(format!("unsupported bond value_order {value}"))),
    }
}
