use chematic_3d::{Coords3D, Point3, constraints::build_constraints};
use chematic_core::AtomIdx;
use chematic_smiles::parse;

#[test]
fn bond_constraints_use_element_and_bond_order_targets() {
    for (smiles, target) in [("C#C", 1.20), ("CF", 1.35), ("CBr", 1.94), ("N[H]", 1.01)] {
        let mol = parse(smiles).unwrap();
        let constraints = build_constraints(&mol);
        assert_eq!(constraints.bonds.len(), 1, "{smiles}");
        assert!(
            (constraints.bonds[0].target_distance - target).abs() < 1e-12,
            "{smiles}"
        );
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(1), Point3::new(target, 0.0, 0.0));
        assert_eq!(constraints.violated_count(&coords), 0, "{smiles}");
        coords.set(AtomIdx(1), Point3::new(target + 0.2, 0.0, 0.0));
        assert_eq!(constraints.violated_count(&coords), 1, "{smiles}");
        assert!(
            (constraints.max_violation(&coords) - 0.15).abs() < 1e-12,
            "{smiles}"
        );
    }
}

#[test]
fn maximum_violation_accounts_for_angle_constraints() {
    let mol = parse("CCC").unwrap();
    let constraints = build_constraints(&mol);
    let mut coords = Coords3D::new_zeroed(3);
    coords.set(AtomIdx(0), Point3::new(1.54, 0.0, 0.0));
    coords.set(AtomIdx(2), Point3::new(0.0, 1.54, 0.0));
    assert_eq!(constraints.angles.len(), 1);
    assert_eq!(constraints.violated_count(&coords), 1);
    assert!((constraints.max_violation(&coords) - 14.47_f64.to_radians()).abs() < 1e-12);
}
