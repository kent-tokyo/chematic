//! CDK #1353 demonstrates that aromaticity depends on the chosen model and
//! sulfur-oxide representation. These assertions document CheMatic's two
//! explicit policies; neither count is a model-independent chemical verdict.

use chematic_perception::{apply_aromaticity, apply_aromaticity_rdkit_parity_experimental};

fn aromatic_atom_count(mol: &chematic_core::Molecule) -> usize {
    mol.atoms().filter(|(_, atom)| atom.aromatic).count()
}

#[test]
fn thiophene_oxide_forms_have_explicit_model_specific_results() {
    for (smiles, native_count, rdkit_count) in [
        ("O=S1C=CC=C1", 0, 0),
        ("[O-][S+]1C=CC=C1", 0, 5),
        ("S1C=CC=C1", 5, 5),
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let native = apply_aromaticity(&mol);
        let rdkit = apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
        assert_eq!(aromatic_atom_count(&native), native_count, "{smiles}");
        assert_eq!(aromatic_atom_count(&rdkit), rdkit_count, "{smiles}");
    }
}
