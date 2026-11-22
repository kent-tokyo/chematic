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
