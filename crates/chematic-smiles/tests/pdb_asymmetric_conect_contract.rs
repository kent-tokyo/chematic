use chematic_core::BondOrder;
use chematic_smiles::rdkit_mol_from_pdb_block;

#[test]
fn asymmetric_bidirectional_conect_follows_pdb_record_order() {
    // RDKit 2025.09.3 sanitize=false, removeHs=false, proximityBonding=false.
    let atoms = "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\nHETATM    2  C2  LIG A   1       1.300   0.000   0.000  1.00  0.00           C\n";
    for (records, expected) in [
        (
            "CONECT    1    2    2\nCONECT    2    1    1    1\n",
            BondOrder::Double,
        ),
        (
            "CONECT    2    1    1    1\nCONECT    1    2    2\n",
            BondOrder::Triple,
        ),
        (
            "CONECT    1    2    2    2    2\nCONECT    2    1    1    1    1\n",
            BondOrder::Quadruple,
        ),
    ] {
        let parsed =
            rdkit_mol_from_pdb_block(&format!("{atoms}{records}END\n"), false, false, 0, false)
                .unwrap()
                .unwrap();
        assert_eq!(
            parsed.molecule.bonds().next().unwrap().1.order,
            expected,
            "{records}"
        );
    }
}

#[test]
fn reciprocal_conect_multiplicity_does_not_override_first_direction() {
    // RDKit 2025.09.3, sanitize=false, removeHs=false, proximityBonding=false.
    let atoms = "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\nHETATM    2  C2  LIG A   1       1.300   0.000   0.000  1.00  0.00           C\n";
    let orders = [
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Quadruple,
    ];
    for (first_count, expected) in (1..=4).zip(orders) {
        for second_count in 1..=4 {
            for (source, destination) in [(1, 2), (2, 1)] {
                let record = |from, to, count: usize| {
                    format!("CONECT{from:>5}{}\n", format!("{to:>5}").repeat(count))
                };
                let records = format!(
                    "{}{}",
                    record(source, destination, first_count),
                    record(destination, source, second_count)
                );
                let parsed = rdkit_mol_from_pdb_block(
                    &format!("{atoms}{records}END\n"),
                    false,
                    false,
                    0,
                    false,
                )
                .unwrap()
                .unwrap();
                assert_eq!(
                    parsed.molecule.bonds().next().unwrap().1.order,
                    expected,
                    "{records}"
                );
            }
        }
    }
}
