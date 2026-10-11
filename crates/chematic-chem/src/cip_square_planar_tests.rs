use super::*;
use chematic_smiles::parse;

#[test]
fn accurate_mode_error_preserves_comparison_diagnostics() {
    let error = CipModeError::Accurate(chematic_cip::CipCompareError::InvalidDigraph(
        "broken fixture".into(),
    ));
    assert_eq!(
        error.to_string(),
        "accurate CIP engine error: invalid digraph: broken fixture"
    );
}

/// Regression for the square-planar-stereo PR's required CIP safety fix:
/// `tetrahedral_stereo_neighbors`/`assign_tetrahedral` gated on
/// `chirality == Chirality::None` (equality), not an exhaustive match --
/// adding `Chirality::SquarePlanar` did NOT force a compile error there, so
/// an `@SP1`-tagged 4-neighbor Pt center would have silently fallen through
/// into the tetrahedral R/S algorithm and produced a bogus CIP code. Must
/// return no code at all, not a wrong one.
#[test]
fn square_planar_center_never_gets_a_bogus_tetrahedral_cip_code() {
    let mol = parse("N->[Pt@SP1](<-N)(Cl)Cl").unwrap();
    let pt = (0..mol.atom_count())
        .map(|i| chematic_core::AtomIdx(i as u32))
        .find(|&i| mol.atom(i).element == chematic_core::Element::PT)
        .expect("has a Pt atom");
    assert!(
        matches!(mol.atom(pt).chirality, Chirality::SquarePlanar(_)),
        "fixture must actually carry a SquarePlanar tag"
    );
    assert_eq!(
        tetrahedral_stereo_neighbors(&mol, pt),
        None,
        "a square-planar center must never resolve a tetrahedral CIP code"
    );
    assert_eq!(
        assign_cip(&mol).get(pt),
        None,
        "assign_cip must not assign any code to a square-planar center"
    );
}
