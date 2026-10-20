use chematic_chem::standardize::normalize_groups;
use chematic_smiles::parse;

#[test]
fn normalize_groups_preserves_mapped_nitro_oxygen_identity() {
    for input in ["C[N+](=O)[O-:28]", "C[N+](=[O:28])[O-]"] {
        let mol = parse(input).unwrap();
        let normalized = normalize_groups(&mol);
        for (idx, atom) in mol.atoms() {
            assert_eq!(normalized.atom(idx).charge, atom.charge, "{input}");
            assert_eq!(normalized.atom(idx).atom_map, atom.atom_map, "{input}");
        }
        for (idx, bond) in mol.bonds() {
            assert_eq!(normalized.bond(idx).order, bond.order, "{input}");
        }
    }
}

#[test]
fn normalize_groups_keeps_sulfites_and_sulfinates_neutral() {
    for input in ["COS(=O)OC", "CS(=O)OC", "O=S1OCCO1", "O=S1O[C@H](C)CO1"] {
        let mol = parse(input).unwrap();
        let normalized = normalize_groups(&mol);
        assert_eq!(
            chematic_smiles::canonical_smiles(&normalized),
            chematic_smiles::canonical_smiles(&mol),
            "{input}"
        );
        assert_eq!(
            chematic_smiles::canonical_smiles(&normalize_groups(&normalized)),
            chematic_smiles::canonical_smiles(&normalized),
            "idempotence: {input}"
        );
    }
}
