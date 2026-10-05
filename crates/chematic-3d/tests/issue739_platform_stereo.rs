//! #739: the two A6 rows whose stereo-safe MMFF94 run failed with
//! `FinalStereoViolation` on Linux but succeeded on macOS. The force field
//! now takes its transcendental functions from the `libm` crate (as wasm32
//! already did), so the trajectory no longer depends on the host C library;
//! with the A6 benchmark configuration both rows succeed.

use chematic_3d::{PipelineV2Config, RingTorsionApplicationPolicy, embed_pipeline_v2};

fn a6_config() -> PipelineV2Config {
    let mut config = PipelineV2Config::stereo_safe(
        chematic_3d::minimize::ForceFieldPolicy::Mmff94BondAngleStrict,
    );
    config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
    config.fail_on_unevaluable_stereo = false;
    config.embed.random_seed = 20_260_801;
    config.embed.max_attempts = 8;
    config.embed.use_exp_torsions = true;
    config.embed.use_small_ring_torsions = true;
    config.embed.use_macrocycle_torsions = true;
    config.embed.use_macrocycle_14_bounds = true;
    config.include_legacy_torsion_heuristic = false;
    config.force_field_max_iterations = 200;
    config.gate_mmff94_torsion_oop = false;
    config.gate_mmff94_stretch_bend = false;
    config.total_timeout_ms = Some(20_000);
    config
}

#[test]
fn penam_rows_keep_their_stereo() {
    for (row, smiles) in [
        (53, "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O"),
        (
            246,
            "CC1(C)[C@H](C(=O)O)N2C(=O)[C@@H](Cc3cn(-c4ccc(F)cc4)nn3)[C@H]2S1(=O)=O",
        ),
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let result = embed_pipeline_v2(&mol, &a6_config());
        assert!(
            result.is_ok(),
            "row {row}: {:?}",
            result.err().map(|e| e.cause)
        );
    }
}
