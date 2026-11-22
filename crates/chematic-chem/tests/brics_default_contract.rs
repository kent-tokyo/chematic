use chematic_chem::{BricsConfig, brics_fragments, brics_fragments_with_config};
use chematic_smiles::{canonical_smiles, parse};

#[test]
fn default_brics_configuration_retains_all_fragments() {
    let mol = parse("CCCC").unwrap();
    let names = |molecules: Vec<chematic_core::Molecule>| {
        let mut smiles: Vec<_> = molecules.iter().map(canonical_smiles).collect();
        smiles.sort();
        smiles
    };
    let all = names(brics_fragments(&mol));
    assert_eq!(all.len(), 2);
    assert_eq!(
        names(brics_fragments_with_config(&mol, &BricsConfig::default())),
        all
    );
    assert!(
        brics_fragments_with_config(
            &mol,
            &BricsConfig {
                min_fragment_size: 3
            }
        )
        .is_empty()
    );
}
