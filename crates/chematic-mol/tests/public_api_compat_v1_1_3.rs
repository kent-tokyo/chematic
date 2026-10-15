//! Compile-time source-compatibility guard for public error enums from v1.1.3.
//!
//! These exhaustive matches intentionally list the complete v1.1.3 surfaces.
//! Adding a variant in a patch release makes this test fail to compile, which
//! catches the same E0004 break experienced by downstream Rust users.

use chematic_mol::{CifError, CmlError};

fn classify_cml(error: CmlError) -> &'static str {
    match error {
        CmlError::MissingMolecule => "missing_molecule",
        CmlError::EmptyMolecule => "empty_molecule",
        CmlError::MalformedXml(_) => "malformed_xml",
        CmlError::UnknownElement(_) => "unknown_element",
        CmlError::UnknownAtomRef(_) => "unknown_atom_ref",
        CmlError::InvalidAtomRefs2(_) => "invalid_atom_refs2",
        CmlError::InvalidBondOrder(_) => "invalid_bond_order",
        CmlError::NonFiniteCoordinate(_) => "non_finite_coordinate",
        CmlError::ResourceLimit { .. } => "resource_limit",
    }
}

fn classify_cif(error: CifError) -> &'static str {
    match error {
        CifError::NoAtomSiteLoop => "no_atom_site_loop",
        CifError::MissingCoordinateColumns => "missing_coordinate_columns",
        CifError::UnknownElement(_) => "unknown_element",
        CifError::InvalidCoordinate(_) => "invalid_coordinate",
        CifError::InvalidCellParameters(_) => "invalid_cell_parameters",
        CifError::MissingCellParameters => "missing_cell_parameters",
        CifError::ResourceLimit { .. } => "resource_limit",
    }
}

#[test]
fn v1_1_3_error_enum_matches_remain_exhaustive() {
    assert_eq!(
        classify_cml(CmlError::MalformedXml("x".into())),
        "malformed_xml"
    );
    assert_eq!(
        classify_cif(CifError::InvalidCoordinate("x".into())),
        "invalid_coordinate"
    );
}
