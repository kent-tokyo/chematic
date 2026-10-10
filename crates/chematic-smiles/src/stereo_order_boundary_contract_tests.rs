use crate::{canonical_smiles, parse, rdkit_canonical_smiles};
use chematic_core::Chirality;

#[test]
fn tetrahedral_output_preserves_each_recorded_hydrogen_slot() {
    for source in [
        "F[C@H](Cl)Br",
        "[C@@H](F)(Cl)Br",
        "N[C@H](C)O",
        "F[C@H]1CCCC1Cl",
        "[2H][C@](F)(Cl)Br",
    ] {
        let mol = parse(source).unwrap();
        let expected = rdkit_canonical_smiles(&mol).unwrap();
        let expected_native = canonical_smiles(&mol);
        for (center, atom) in mol.atoms().filter(|(_, a)| a.chirality.is_tetrahedral()) {
            let original = mol.stereo_neighbor_order(center).unwrap();
            assert_eq!(original.len(), 4);
            for a in 0..4 {
                for b in 0..4 {
                    if a == b {
                        continue;
                    }
                    for c in 0..4 {
                        if c == a || c == b {
                            continue;
                        }
                        let d = (0..4).find(|&d| d != a && d != b && d != c).unwrap();
                        let indices = [a, b, c, d];
                        let inversions = (0..4)
                            .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
                            .filter(|&(i, j)| indices[i] > indices[j])
                            .count();
                        let mut changed = mol.clone();
                        changed.set_stereo_neighbor_order(
                            center,
                            indices.iter().map(|&i| original[i]).collect(),
                        );
                        let clockwise =
                            (atom.chirality == Chirality::Clockwise) != (inversions % 2 == 1);
                        changed.set_chirality(
                            center,
                            if clockwise {
                                Chirality::Clockwise
                            } else {
                                Chirality::CounterClockwise
                            },
                        );
                        assert_eq!(
                            rdkit_canonical_smiles(&changed).unwrap(),
                            expected,
                            "{source}: {indices:?}"
                        );
                        assert_eq!(
                            canonical_smiles(&changed),
                            expected_native,
                            "{source}: {indices:?}"
                        );
                    }
                }
            }
        }
    }
}
