// Reproduction harness for a clean baseline checkout. Append this module to
// crates/chematic-3d/src/rdkit_embed/bounds.rs, then run:
// cargo test -p chematic-3d --lib stereo_relation_must_not_change --offline
// The baseline compatibility predicate fails this assertion. The additive
// public query is tested separately; that predicate remains a parity detail.
#[cfg(test)]
mod issue9719_order_probe {
    use super::*;

    #[test]
    fn stereo_relation_must_not_change_when_the_path_is_reversed() {
        let mol = chematic_smiles::parse("F/C(Cl)=C(Br)/I").unwrap();
        let view = chematic_smiles::rdkit_mol_view(&mol).unwrap();
        let (b, bond) = view.bonds.iter().enumerate()
            .find(|(_, bond)| bond.stereo > STEREO_ANY).unwrap();
        for a in view.neighbors(bond.begin).filter(|&a| a != bond.end) {
            for d in view.neighbors(bond.end).filter(|&d| d != bond.begin) {
                assert_eq!(atom_stereo(&view, b, a, d), atom_stereo(&view, b, d, a));
            }
        }
    }
}
