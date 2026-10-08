//! #739: the two A6 rows whose stereo-safe MMFF94 run failed with
//! `FinalStereoViolation` on Linux but succeeded on macOS. The force field
//! now takes its transcendental functions from the `libm` crate (as wasm32
//! already did), so the trajectory no longer depends on the host C library;
//! with the A6 benchmark configuration both rows succeed.
//!
//! The libm change only moved the failure to other seeds: from about a
//! quarter of embeddings, MMFF94 relaxes these penam bridgeheads through
//! inversion, on every platform. A stereo-safe run that ends in a final
//! stereo violation now embeds again from other seeds; the seeds below all
//! failed before that.

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

#[test]
fn penams_keep_their_stereo_from_seeds_that_used_to_invert_them() {
    // Seeds whose run ended in `FinalStereoViolation` with the published
    // v1.0.37 wheel (libm math), at the A6 configuration.
    for (name, smiles, seeds) in [
        (
            "penicillin_core",
            "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O",
            &[10_u64, 11, 13, 14, 25][..],
        ),
        (
            "chembl_tier_b_0181",
            "CC1(C)[C@H](C(=O)O)N2C(=O)[C@@H](Cc3cn(-c4ccc(F)cc4)nn3)[C@H]2S1(=O)=O",
            &[2, 7, 8][..],
        ),
        (
            "chembl_tier_b_0175",
            "CC1(C)[C@H](C(=O)O)N2C(=O)[C@@H](Cc3cn(CC4CCCN4)nn3)[C@H]2S1(=O)=O",
            &[1, 2, 4, 5][..],
        ),
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        for &seed in seeds {
            let mut config = a6_config();
            config.embed.random_seed = seed;
            let result = embed_pipeline_v2(&mol, &config)
                .unwrap_or_else(|e| panic!("{name} seed {seed}: {:?}", e.cause));
            assert!(
                result.final_stereo.is_fully_satisfied(),
                "{name} seed {seed}"
            );
            assert!(result.final_validation.sound, "{name} seed {seed}");
            assert_eq!(
                result.final_validation.gross_clash_count, 0,
                "{name} seed {seed}"
            );
        }
    }
}

#[test]
fn verify_only_still_reports_the_violation_without_reseeding() {
    // The reseed is part of the stereo-safe (`RepairAndVerify`) contract
    // only; `VerifyOnly` keeps its typed failure at the seed it was given.
    let mol = chematic_smiles::parse("CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O").unwrap();
    let mut config = a6_config();
    config.stereo_policy = chematic_3d::StereoPolicy::VerifyOnly;
    let mut failures = 0;
    for seed in [10_u64, 11, 13, 14, 25] {
        config.embed.random_seed = seed;
        if let Err(e) = embed_pipeline_v2(&mol, &config) {
            assert!(
                matches!(
                    e.cause,
                    chematic_3d::PipelineV2FailureCause::FinalStereoViolation
                ),
                "seed {seed}: {:?}",
                e.cause
            );
            failures += 1;
        }
    }
    assert!(
        failures > 0,
        "VerifyOnly must still surface violations at these seeds"
    );
}
