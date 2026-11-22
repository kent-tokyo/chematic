use chematic_smarts::rdkit_parity_match::rdkit_relevant_ring_count_presence;
use chematic_smarts::{RdkitParityConfig, RdkitParityError};

#[test]
fn presence_distinguishes_acyclic_ring_and_fused_atoms() {
    for (smiles, mask) in [("CCC", 1), ("C1CCCCC1", 2), ("CC1CCC2CCCCC2C1", 7)] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let presence =
            rdkit_relevant_ring_count_presence(&mol, &RdkitParityConfig::default()).unwrap();
        assert_eq!(presence, [mask, 0, 0, 0], "{smiles}");
    }
}

#[test]
fn presence_respects_ring_search_budget_after_a_cached_success() {
    let mol = chematic_smiles::parse("C1CCCCC1").unwrap();
    let mut config = RdkitParityConfig::default();
    assert_eq!(
        rdkit_relevant_ring_count_presence(&mol, &config).unwrap(),
        [2, 0, 0, 0]
    );
    config.ring_model_budget.max_candidates = 0;
    assert!(matches!(
        rdkit_relevant_ring_count_presence(&mol, &config),
        Err(RdkitParityError::RingModelBudgetExceeded { cap: 0, .. })
    ));
}
