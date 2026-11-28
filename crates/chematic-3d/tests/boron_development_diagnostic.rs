use chematic_3d::{generate_and_minimize_dreiding, get_bond_angle_deg};

#[test]
fn overlapping_boron_substituents_are_separated_before_minimization() {
    for smiles in [
        "CC(C)(C)c1sc(C(C)(C)C)c(-c2c(F)c(F)c(F)c(F)c2F)c1B(c1c(F)c(F)c(F)c(F)c1F)c1c(F)c(F)c(F)c(F)c1F",
        "Cc1cc(C)c(B(c2ccc(-n3c(-c4ccccc4)nc4c5ccccc5c5ccccc5c43)cc2)c2c(C)cc(C)cc2C)c(C)c1",
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let boron = mol
            .atoms()
            .find(|(_, atom)| atom.element.atomic_number() == 5)
            .unwrap()
            .0;
        let neighbors: Vec<_> = mol.neighbors(boron).map(|(idx, _)| idx).collect();
        let angles = |coords: &chematic_3d::Coords3D| {
            let mut values = Vec::new();
            for i in 0..neighbors.len() {
                for j in i + 1..neighbors.len() {
                    values.push(get_bond_angle_deg(
                        coords,
                        neighbors[i],
                        boron,
                        neighbors[j],
                    ));
                }
            }
            values.sort_by(f64::total_cmp);
            values
        };
        let result = angles(&generate_and_minimize_dreiding(&mol));
        assert!(
            result.iter().all(|angle| (110.0..=130.0).contains(angle)),
            "{result:?}"
        );
    }
}
