use chematic_chem::standardize::{CleanupError, cleanup, cleanup_checked, cleanup_unchecked};

#[test]
fn cleanup_errors_have_stable_codes_and_descriptions() {
    for (error, code, description) in [
        (
            CleanupError::AmbiguousMetalHydrogen,
            "cleanup_metal_hydrogen",
            "Cleanup cannot account for the metal hydrogen change",
        ),
        (
            CleanupError::UnstableMetalDisconnection,
            "cleanup_unstable_metal",
            "Cleanup would charge the detached metal on a subsequent call",
        ),
        (
            CleanupError::UnstableMetalStereo,
            "cleanup_metal_stereo",
            "Cleanup cannot reliably represent a tetrahedral center with five ligands including a metal",
        ),
    ] {
        assert_eq!(error.code(), code);
        assert_eq!(error.to_string(), description);
    }
}

#[test]
fn cleanup_alias_and_checked_cleanup_agree_on_ordinary_molecules() {
    for smiles in ["CCO", "CC(=O)O", "[Na+].[Cl-]", "C[N+](=O)[O-]"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let unchecked = cleanup_unchecked(&mol);
        let alias = cleanup(&mol);
        let checked = cleanup_checked(&mol).unwrap();
        let canonical = |mol| chematic_smiles::rdkit_smiles(mol, &Default::default()).unwrap();
        assert_eq!(canonical(&unchecked), canonical(&alias));
        assert_eq!(canonical(&unchecked), canonical(&checked));
    }
}
