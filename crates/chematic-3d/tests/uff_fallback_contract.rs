use chematic_3d::{
    dg::generate_coords,
    minimize::{ForceFieldBridgeError, ForceFieldPolicy, MinimizeConfig, minimize_with_policy},
};

#[test]
fn unsupported_mmff_selenium_uses_uff_and_reports_reason() {
    let mol = chematic_smiles::parse("C[Se]C").unwrap();
    let result = minimize_with_policy(
        &mol,
        generate_coords(&mol),
        ForceFieldPolicy::Mmff94WithUffFallback,
        &MinimizeConfig::default(),
    )
    .unwrap();
    assert_eq!(
        result.requested_force_field,
        ForceFieldPolicy::Mmff94WithUffFallback
    );
    assert_eq!(result.actual_force_field_used, ForceFieldPolicy::UffOnly);
    assert!(matches!(
        result.fallback_reason,
        Some(ForceFieldBridgeError::UnsupportedAtomType(_))
    ));
    assert!(result.coords.is_finite());
    assert!(result.energy_after.total().is_finite());
}
