//! Dump the canonical atom sets of the six #337 macrocycle families.

use chematic_perception::find_symmetrized_sssr;

const FIXTURES: &[(&str, &str)] = &[
    (
        "0009",
        "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
    ),
    (
        "0023",
        "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)-c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
    ),
    (
        "0028",
        "C1=C\\c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCNc2cc[n+](c3ccccc23)Cc2ccc/1cc2",
    ),
    (
        "0029",
        "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)CCc2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
    ),
    (
        "0030",
        "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)Cc2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
    ),
    (
        "0034",
        "c1ccc2c(c1)c1cc[n+]2Cc2ccc3c(c2)Cc2cc(ccc2-3)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
    ),
];

fn main() {
    for &(name, smiles) in FIXTURES {
        let mol = chematic_smiles::parse(smiles).expect(name);
        let kekule =
            chematic_core::apply_kekule(&mol, &chematic_core::kekulize(&mol).expect("kekulize"));
        let doubles: Vec<_> = kekule
            .bonds()
            .filter_map(|(_, bond)| {
                (bond.order == chematic_core::BondOrder::Double)
                    .then_some((bond.atom1.0, bond.atom2.0))
            })
            .collect();
        let mut rings: Vec<Vec<u32>> = find_symmetrized_sssr(&mol)
            .rings()
            .iter()
            .filter(|ring| ring.len() >= 20)
            .map(|ring| {
                let mut atoms: Vec<u32> = ring.iter().map(|atom| atom.0).collect();
                atoms.sort_unstable();
                atoms
            })
            .collect();
        rings.sort();
        let mut six_rings: Vec<Vec<u32>> = find_symmetrized_sssr(&mol)
            .rings()
            .iter()
            .filter(|ring| ring.len() == 6)
            .map(|ring| {
                let mut atoms: Vec<u32> = ring.iter().map(|atom| atom.0).collect();
                atoms.sort_unstable();
                atoms
            })
            .collect();
        six_rings.sort();
        let six_paths: Vec<Vec<u32>> = find_symmetrized_sssr(&mol)
            .rings()
            .iter()
            .filter(|ring| ring.len() == 6)
            .map(|ring| ring.iter().map(|atom| atom.0).collect())
            .collect();
        println!(
            "{name} macrocycles={} six_rings={} six_paths={} doubles={doubles:?}",
            serde_json::to_string(&rings).unwrap(),
            serde_json::to_string(&six_rings).unwrap(),
            serde_json::to_string(&six_paths).unwrap()
        );
    }
}
