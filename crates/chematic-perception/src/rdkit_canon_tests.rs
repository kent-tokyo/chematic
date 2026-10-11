#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> Molecule {
        chematic_smiles::parse(s).expect("parses")
    }

    #[test]
    fn ranks_match_rdkit_on_small_molecules() {
        // Chem.CanonicalRankAtoms(Chem.MolFromSmiles(s)), RDKit 2026.03.6.
        let cases: [(&str, &[u32]); 3] = [
            ("O=[N+]([O-])c1ccccc1", &[0, 8, 1, 7, 5, 3, 2, 4, 6]),
            ("CC(=O)O", &[0, 3, 1, 2]),
            ("N[C@@H](C)C(=O)O", &[1, 5, 0, 4, 2, 3]),
        ];
        for (s, want) in cases {
            assert_eq!(rdkit_canonical_atom_ranks(&parse(s)), want, "{s}");
        }
    }

    #[test]
    fn kekule_is_a_valid_structure() {
        let mol = parse("c1ccc2ccccc2c1");
        let valid = chematic_core::kekulize(&mol).expect("kekulizes");
        let kek = rdkit_canonical_kekule(&mol, &valid).expect("rdkit walk");
        let doubles = kek.values().filter(|o| **o == BondOrder::Double).count();
        assert_eq!(doubles, 5);
        let mut deg = vec![0; mol.atom_count()];
        for (b, o) in &kek {
            if *o == BondOrder::Double {
                let bond = mol.bond(*b);
                deg[bond.atom1.0 as usize] += 1;
                deg[bond.atom2.0 as usize] += 1;
            }
        }
        assert!(deg.iter().all(|&d| d == 1));
    }
}
