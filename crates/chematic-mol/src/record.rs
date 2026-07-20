//! [`MoleculeRecord`] — a format-agnostic single-molecule record shared by
//! streaming readers/writers for record-oriented formats (SMILES table
//! files, TDT, MRV) that are not MOL/SDF's own Ctab-based `SdfRecord`.
//!
//! Properties are an insertion-ordered `Vec`, not a `HashMap` — several
//! writers in this module family must reproduce the exact column/tag order
//! a caller configured, and `HashMap` iteration order is not deterministic
//! across runs, which would make round-trip output non-reproducible.

use chematic_core::Molecule;

/// A single molecule plus its file-format-level metadata: name, arbitrary
/// key/value properties (in file order), and optional 2D/3D coordinates.
///
/// Not every record-oriented format populates every field — e.g. a SMILES
/// table record never has `coordinates_2d`/`coordinates_3d` (SMILES text
/// carries no geometry), while a TDT record may have either or both.
#[derive(Clone)]
pub struct MoleculeRecord {
    /// The parsed molecule (heavy atoms only, no explicit H).
    pub mol: Molecule,
    /// The record's name/title, or an empty string if the format/record has none.
    pub name: String,
    /// Arbitrary key/value properties, in the order they appeared in (or
    /// should be written to) the file.
    pub properties: Vec<(String, String)>,
    /// 2D atom coordinates in file order, one entry per atom, if present.
    pub coordinates_2d: Option<Vec<[f64; 2]>>,
    /// 3D atom coordinates in file order, one entry per atom, if present.
    pub coordinates_3d: Option<Vec<[f64; 3]>>,
}

impl MoleculeRecord {
    /// A record with no name, no properties, and no coordinates.
    pub fn new(mol: Molecule) -> Self {
        Self {
            mol,
            name: String::new(),
            properties: Vec::new(),
            coordinates_2d: None,
            coordinates_3d: None,
        }
    }

    /// Look up a property's value by key (first match, in insertion order).
    pub fn get_property(&self, key: &str) -> Option<&str> {
        self.properties
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
}
