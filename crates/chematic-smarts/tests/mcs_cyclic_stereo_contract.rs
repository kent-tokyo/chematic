use chematic_smarts::{McsConfig, McsObjective, McsOutcome, find_mcs_with_objective_checked};

#[test]
fn cyclic_chiral_mcs_matches_reference_atom_and_bond_counts() {
    // Independently generated with RDKit 2025.09.3, MatchChiralTag enabled.
    let cases = [
        (
            "F[C@]12CC[C@](Cl)(C1)CC2",
            "F[C@@]12CC[C@@](Cl)(C1)CC2",
            9,
            10,
        ),
        ("C1[C@H]2CC[C@H]1C2", "C1[C@@H]2CC[C@H]1C2", 6, 6),
    ];
    let config = McsConfig {
        match_chiral_tag: true,
        ..Default::default()
    };
    for (left, right, atoms, bonds) in cases {
        let a = chematic_smiles::parse(left).unwrap();
        let b = chematic_smiles::parse(right).unwrap();
        let result = find_mcs_with_objective_checked(&[&a, &b], &config, McsObjective::BondCount);
        assert!(
            matches!(result, McsOutcome::Exhaustive(_)),
            "{left}, {right}"
        );
        let query = result.into_query();
        assert_eq!(
            (query.atom_count(), query.bonds.len()),
            (atoms, bonds),
            "{left}, {right}"
        );
    }
}
