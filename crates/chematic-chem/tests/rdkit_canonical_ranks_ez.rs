//! `rdkit_canonical_atom_ranks` with double-bond stereo against RDKit
//! 2026.03.6 `Chem.CanonicalRankAtoms(Chem.AddHs(mol), breakTies=True)`.
//! A molecule whose E and Z arms are otherwise identical ranks them apart
//! only through RDKit's `STEREOE`/`STEREOZ` bond invariant.

fn ranks(smiles: &str) -> Vec<u32> {
    let mol = chematic_smiles::parse(smiles).expect("parses");
    chematic_chem::rdkit_canonical_atom_ranks(&chematic_chem::add_hydrogens(&mol))
}

#[test]
fn e_and_z_arms_rank_apart_as_in_rdkit() {
    assert_eq!(
        ranks("C(/C=C/Cl)(/C=C\\Cl)O"),
        [13, 12, 11, 7, 10, 9, 6, 8, 5, 4, 3, 2, 1, 0]
    );
    assert_eq!(
        ranks("OC(/C=C/C)/C=C\\C"),
        [
            12, 19, 16, 15, 18, 14, 13, 17, 0, 11, 4, 3, 8, 9, 10, 2, 1, 5, 6, 7
        ]
    );
    assert_eq!(
        ranks("F/C=C/C(/C=C\\F)C"),
        [9, 12, 13, 15, 11, 10, 8, 14, 2, 3, 7, 1, 0, 4, 5, 6]
    );
}

#[test]
fn without_stereo_bonds_the_perception_ranking_is_unchanged() {
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse("CC(O)C(=O)N").unwrap());
    assert_eq!(
        chematic_chem::rdkit_canonical_atom_ranks(&mol),
        chematic_perception::rdkit_canonical_atom_ranks(&mol)
    );
}
