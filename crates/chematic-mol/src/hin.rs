//! Checked, single-molecule HyperChem HIN records.
//!
//! The supported profile is `mol`, reciprocal `atom` connectivity and `endmol`,
//! with blank lines and `;` comments. Atom names, type/flag tokens, partial
//! charges and Cartesian coordinates (Å) are retained independently of the
//! graph. Bond tokens `s/d/t/a` mean single/double/triple/aromatic; `a` never
//! becomes the MDL query bond "single or double" (Open Babel issue #3000).
//!
//! Residue/system directives, multiple molecules and unknown bond tokens reject
//! rather than disappearing. The graph projection does not interpret force-field
//! types/flags, assign formal charges from partial charges, infer stereochemistry
//! from coordinates or add hydrogen atoms. Formal charges default to zero;
//! implicit-H queries retain the core's ordinary valence rules. Use the immutable
//! source record for exact HIN retention, not the projected graph alone.

mod json;
mod reader;

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

/// Resource bounds applied before constructing a source record.
#[derive(Debug, Clone, Copy)]
pub struct HinParseLimits {
    pub max_bytes: usize,
    pub max_atoms: usize,
    pub max_bonds: usize,
}
impl Default for HinParseLimits {
    fn default() -> Self {
        Self {
            max_bytes: 1 << 20,
            max_atoms: 4096,
            max_bonds: 16384,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HinError {
    ResourceLimit,
    InvalidRecord { line: usize, reason: &'static str },
}
impl std::fmt::Display for HinError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ResourceLimit => f.write_str("HIN resource limit exceeded"),
            Self::InvalidRecord { line, reason } => {
                write!(f, "invalid HIN at line {line}: {reason}")
            }
        }
    }
}
impl std::error::Error for HinError {}

/// Literal source atom; source IDs need not be consecutive or sorted.
#[derive(Debug, Clone, PartialEq)]
pub struct HinAtom {
    pub source_id: u32,
    pub name: String,
    pub element: Element,
    pub atom_type: String,
    pub flags: String,
    pub partial_charge: f64,
    pub coordinates: [f64; 3],
}

/// Bond endpoints are zero-based indices in source atom-row order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HinBond {
    pub atom1: usize,
    pub atom2: usize,
    pub order: BondOrder,
}

/// Validated immutable source plus a typed, bounded molecular view.
#[derive(Debug, Clone)]
pub struct HinRecord {
    source: String,
    molecule_id: u32,
    title: String,
    atoms: Vec<HinAtom>,
    bonds: Vec<HinBond>,
}
impl HinRecord {
    pub fn source(&self) -> &str {
        &self.source
    }
    pub fn molecule_id(&self) -> u32 {
        self.molecule_id
    }
    pub fn title(&self) -> &str {
        &self.title
    }
    pub fn atoms(&self) -> &[HinAtom] {
        &self.atoms
    }
    pub fn bonds(&self) -> &[HinBond] {
        &self.bonds
    }

    /// Literal graph projection in atom-row order; see module-level boundaries.
    pub fn to_molecule(&self) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let mut aromatic = vec![false; self.atoms.len()];
        for bond in &self.bonds {
            if bond.order == BondOrder::Aromatic {
                aromatic[bond.atom1] = true;
                aromatic[bond.atom2] = true;
            }
        }
        for (entry, aromatic) in self.atoms.iter().zip(aromatic) {
            let mut atom = Atom::new(entry.element);
            atom.aromatic = aromatic;
            builder.add_atom(atom);
        }
        for bond in &self.bonds {
            builder
                .add_bond(
                    AtomIdx(bond.atom1 as u32),
                    AtomIdx(bond.atom2 as u32),
                    bond.order,
                )
                .expect("private HIN record has validated distinct bond endpoints");
        }
        builder.build()
    }
}

/// Parse the bounded single-molecule profile with default resource limits.
pub fn parse_hin_record(source: &str) -> Result<HinRecord, HinError> {
    parse_hin_record_with_limits(source, HinParseLimits::default())
}

/// Parse with caller-selected bounds; malformed or unsupported content refuses.
pub fn parse_hin_record_with_limits(
    source: &str,
    limits: HinParseLimits,
) -> Result<HinRecord, HinError> {
    reader::parse(source, limits)
}

/// Exact original UTF-8 text, including comments, line endings and numeric spelling.
pub fn write_hin_record(record: &HinRecord) -> String {
    record.source.clone()
}
