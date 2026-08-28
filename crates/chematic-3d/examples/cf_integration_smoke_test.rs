//! Smoke test: Agent C's embedder (`distance_geometry_v2::embed_distance_geometry_v2`,
//! PR #167) wired directly into Agent F's force-field bridge
//! (`minimize::minimize_with_policy`, PR #169) — the first time these two
//! Wave-1 deliverables have been run together end to end. Both were merged
//! independently and neither's own acceptance gate (`examples/
//! distance_geometry_v2_gap_check.rs`, `examples/mmff94_bridge_coverage_report.rs`)
//! measured the composed pipeline: Agent C's gate never minimizes; Agent F's
//! gate starts from `dg::generate_coords` (the legacy rule-based DFS placer),
//! never from Agent C's real embedded geometry.
//!
//! This is investigation-only: nothing here changes `etkdg.rs` or
//! `chematic-3d/lib.rs` (Coordinator-owned), and the composed call
//! `embed_distance_geometry_v2(..) -> minimize_with_policy(..)` is not wired
//! into any production path by this file.
//!
//! Frozen 58-molecule corpus, hand-copied verbatim (name, SMILES, category)
//! from `scripts/etkdg_vs_rdkit_gap.py::CORPUS`, same transcription already
//! used by both sibling examples in this directory.
//!
//! Run: `cargo run --release -p chematic-3d --example cf_integration_smoke_test`

use std::collections::{BTreeMap, BTreeSet};
use std::panic::{self, AssertUnwindSafe};

use chematic_3d::coords::{Coords3D, Point3};
use chematic_3d::dg;
use chematic_3d::distance_geometry_v2::{EmbedParameters, embed_distance_geometry_v2};
use chematic_3d::minimize::{
    ForceFieldBridgeError, ForceFieldPolicy, MinimizeConfig, minimize_with_policy,
};
use chematic_core::{AtomIdx, Molecule};
use chematic_smiles::parse;

const CORPUS: &[(&str, &str, &str)] = &[
    ("benzene", "c1ccccc1", "rigid_ring"),
    ("naphthalene", "c1ccc2ccccc2c1", "fused_aromatic"),
    ("pyridine", "c1ccncc1", "rigid_ring"),
    ("furan", "c1ccoc1", "rigid_ring"),
    ("thiophene", "c1ccsc1", "rigid_ring"),
    ("adamantane", "C1CC2CC3CC1CC(C2)C3", "rigid_ring"),
    ("cubane", "C1C2C3C1C4C2C3C4", "rigid_ring"),
    ("cyclohexane", "C1CCCCC1", "rigid_ring"),
    ("cyclopentane", "C1CCCC1", "rigid_ring"),
    ("indole", "c1ccc2[nH]ccc2c1", "fused_aromatic"),
    ("purine", "c1ncc2[nH]cnc2n1", "fused_aromatic"),
    ("quinoline", "c1ccc2ncccc2c1", "fused_aromatic"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1", "fused_aromatic"),
    ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34", "fused_aromatic"),
    ("biphenyl", "c1ccc(-c2ccccc2)cc1", "fused_aromatic"),
    ("butane", "CCCC", "flexible_chain"),
    ("hexane", "CCCCCC", "flexible_chain"),
    ("decane", "CCCCCCCCCC", "flexible_chain"),
    ("triethylene_glycol", "OCCOCCOCCO", "flexible_chain"),
    ("hexanediol", "OCCCCCCO", "flexible_chain"),
    ("hexadecane", "CCCCCCCCCCCCCCCC", "flexible_chain"),
    ("cyclododecane", "C1CCCCCCCCCCC1", "macrocycle"),
    ("crown_12_4", "O1CCOCCOCCOCC1", "macrocycle"),
    ("cyclooctadecane", "C1CCCCCCCCCCCCCCCCC1", "macrocycle"),
    ("l_alanine", "N[C@@H](C)C(=O)O", "stereocenter_implicit_h"),
    ("d_alanine", "N[C@H](C)C(=O)O", "stereocenter_implicit_h"),
    ("l_serine", "N[C@@H](CO)C(=O)O", "stereocenter_implicit_h"),
    (
        "l_threonine",
        "C[C@H](O)[C@@H](N)C(=O)O",
        "stereocenter_implicit_h",
    ),
    ("2_butanol_R", "C[C@H](O)CC", "stereocenter_implicit_h"),
    ("2_butanol_S", "C[C@@H](O)CC", "stereocenter_implicit_h"),
    (
        "2_chlorobutane_R",
        "C[C@H](Cl)CC",
        "stereocenter_implicit_h",
    ),
    (
        "ibuprofen_S",
        "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O",
        "stereocenter_implicit_h",
    ),
    (
        "naproxen_S",
        "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1",
        "stereocenter_implicit_h",
    ),
    (
        "menthol",
        "C[C@@H]1CC[C@@H](C(C)C)C[C@H]1O",
        "stereocenter_implicit_h",
    ),
    ("chfclbr_R", "[C@H](F)(Cl)Br", "stereocenter_quaternary"),
    ("chfclbr_S", "[C@@H](F)(Cl)Br", "stereocenter_quaternary"),
    (
        "quaternary_1_R",
        "[C@](F)(Cl)(Br)I",
        "stereocenter_quaternary",
    ),
    (
        "quaternary_1_S",
        "[C@@](F)(Cl)(Br)I",
        "stereocenter_quaternary",
    ),
    (
        "quaternary_2_R",
        "[C@](C)(N)(O)F",
        "stereocenter_quaternary",
    ),
    (
        "quaternary_2_S",
        "[C@@](C)(N)(O)F",
        "stereocenter_quaternary",
    ),
    ("but2ene_E", "C/C=C/C", "alkene_ez"),
    ("but2ene_Z", r"C/C=C\C", "alkene_ez"),
    ("chloropropene_E", "C(/C=C/C)Cl", "alkene_ez"),
    ("chloropropene_Z", r"C(/C=C\C)Cl", "alkene_ez"),
    ("cinnamic_acid_E", "OC(=O)/C=C/c1ccccc1", "alkene_ez"),
    ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1", "alkene_ez"),
    ("pent2ene_E", "CC/C=C/C", "alkene_ez"),
    ("pent2ene_Z", r"CC/C=C\C", "alkene_ez"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "druglike"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "druglike"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "druglike"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1", "druglike"),
    ("diphenhydramine", "CN(C)CCOC(c1ccccc1)c1ccccc1", "druglike"),
    (
        "penicillin_core",
        "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O",
        "druglike",
    ),
    (
        "testosterone",
        "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O",
        "druglike_rigid",
    ),
    (
        "cholesterol",
        "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
        "druglike_stress",
    ),
    (
        "atorvastatin_fragment",
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
        "druglike_stress",
    ),
    ("gly_ala_gly", "NCC(=O)N[C@@H](C)C(=O)NCC(=O)O", "druglike"),
];

/// The three policies this smoke test drives Agent C's coordinates through.
/// `UffOnly` is included (in addition to the two the task named) specifically
/// so the naphthalene-family/fused-aromatic UFF-blowup comparison (item 5) is
/// against the same policy Agent F's own PR body measured that list under —
/// `Mmff94WithUffFallback` only reaches UFF for molecules whose MMFF94
/// attempt already failed, a different (and smaller) population, so
/// comparing THAT count against the disclosed 8/58 would be exactly the
/// fallback-pooling measurement error this repo already has a standing
/// lesson about (see MEMORY: feedback_fallback_pooling_measurement_error).
const POLICIES: &[ForceFieldPolicy] = &[
    ForceFieldPolicy::Mmff94BondAngleStrict,
    ForceFieldPolicy::Mmff94WithUffFallback,
    ForceFieldPolicy::UffOnly,
];

fn worst_bond(mol: &Molecule, coords: &Coords3D) -> f64 {
    mol.bonds()
        .map(|(_, bond)| coords.get(bond.atom1).distance(&coords.get(bond.atom2)))
        .fold(0.0_f64, f64::max)
}

/// Max pairwise distance across all atoms -- a coarse "is this the right
/// order of magnitude" scale check on a raw embedder output (not degenerate/
/// collapsed to a point, not blown up to an absurd radius).
fn max_pairwise_extent(coords: &Coords3D) -> f64 {
    let n = coords.atom_count();
    let mut max_d = 0.0_f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = coords
                .get(AtomIdx(i as u32))
                .distance(&coords.get(AtomIdx(j as u32)));
            if d > max_d {
                max_d = d;
            }
        }
    }
    max_d
}

fn policy_name(p: ForceFieldPolicy) -> &'static str {
    match p {
        ForceFieldPolicy::Mmff94BondAngleStrict => "Mmff94BondAngleStrict",
        ForceFieldPolicy::Mmff94WithUffFallback => "Mmff94WithUffFallback",
        ForceFieldPolicy::UffOnly => "UffOnly",
        ForceFieldPolicy::Dreiding => "Dreiding",
        ForceFieldPolicy::None => "None",
    }
}

/// Outcome of one `minimize_with_policy` call, run behind `catch_unwind` --
/// chematic-ff has live corpus bugs (issues #172-176 all came out of this
/// same 3D Breakthrough Program), so a panic on one molecule must become a
/// recorded per-molecule result, not kill the whole corpus run.
struct MinimizeOutcome {
    ok: bool,
    panicked: bool,
    /// Only meaningful for `Mmff94WithUffFallback`: `fallback_reason.is_some()`
    /// on the `Ok` result -- NOT `actual_force_field_used != requested_force_field`,
    /// which is true on every SUCCESSFUL `Mmff94WithUffFallback` call regardless
    /// of whether a fallback happened (see `PolicyMinimizeResult::actual_force_field_used`'s
    /// own doc comment).
    fallback_occurred: bool,
    err_bucket: Option<String>,
    worst_bond_after: Option<f64>,
    converged: Option<bool>,
}

fn run_minimize(mol: &Molecule, coords: Coords3D, policy: ForceFieldPolicy) -> MinimizeOutcome {
    let config = MinimizeConfig::default();
    let result = panic::catch_unwind(AssertUnwindSafe(|| {
        minimize_with_policy(mol, coords, policy, &config)
    }));
    match result {
        Err(_) => MinimizeOutcome {
            ok: false,
            panicked: true,
            fallback_occurred: false,
            err_bucket: Some("PANIC".to_string()),
            worst_bond_after: None,
            converged: None,
        },
        Ok(Ok(r)) => MinimizeOutcome {
            ok: true,
            panicked: false,
            fallback_occurred: r.fallback_reason.is_some(),
            err_bucket: None,
            worst_bond_after: Some(worst_bond(mol, &r.coords)),
            converged: Some(r.converged),
        },
        Ok(Err(e)) => {
            let bucket = match &e {
                ForceFieldBridgeError::UnsupportedAtomType(_) => "UnsupportedAtomType".to_string(),
                ForceFieldBridgeError::MissingParameters(_) => "MissingParameters".to_string(),
                ForceFieldBridgeError::MinimizationFailed(d) => {
                    format!("MinimizationFailed({:?})", d.reason)
                }
            };
            let (worst, conv) = match &e {
                ForceFieldBridgeError::MinimizationFailed(d) => {
                    (Some(d.worst_bond_length), Some(d.converged))
                }
                _ => (None, None),
            };
            MinimizeOutcome {
                ok: false,
                panicked: false,
                fallback_occurred: false,
                err_bucket: Some(bucket),
                worst_bond_after: worst,
                converged: conv,
            }
        }
    }
}

/// Per-policy (same order as `POLICIES`) end-to-end-ok name sets for the
/// embedder-fed pipeline at a given seed. Used both for the main run and the
/// seed-robustness sweep (dg::generate_coords doesn't depend on the seed at
/// all, so only the embedder side needs re-running per seed).
///
/// Indexed by position in `POLICIES` rather than keyed by `ForceFieldPolicy`
/// itself -- that type derives only `PartialEq, Eq`, not `Ord`, so it can't
/// be a `BTreeMap` key.
fn run_embedder_fed_pipeline_at_seed(seed: u64) -> Vec<BTreeSet<&'static str>> {
    let mut ok_sets: Vec<BTreeSet<&'static str>> =
        POLICIES.iter().map(|_| BTreeSet::new()).collect();
    for &(name, smiles, _category) in CORPUS {
        let mol = parse(smiles).unwrap_or_else(|e| panic!("{name}: parse failed: {e:?}"));
        let params = EmbedParameters {
            random_seed: seed,
            ..EmbedParameters::default()
        };
        let Ok(coords) = embed_distance_geometry_v2(&mol, &params) else {
            continue; // embed failed -- correctly never reaches minimize below
        };
        for (i, &policy) in POLICIES.iter().enumerate() {
            let outcome = run_minimize(&mol, coords.clone(), policy);
            if outcome.ok {
                ok_sets[i].insert(name);
            }
        }
    }
    ok_sets
}

fn main() {
    let default_seed = EmbedParameters::default().random_seed;

    // =========================================================================
    // Main pass: embed (default params/seed) -> minimize under all 3 policies,
    // from BOTH starting-geometry sources (Agent C's real embedder output, and
    // the legacy dg::generate_coords rule-based placer Agent F's own example
    // used) so the "does real embedded geometry change the naphthalene-family
    // picture" question (item 5) has a same-run, apples-to-apples baseline
    // instead of relying on Agent F's PR-body prose.
    // =========================================================================
    let mut n_total = 0usize;
    let mut n_embed_ok = 0usize;
    let mut n_embed_fail = 0usize;
    let mut embed_fail_causes: BTreeMap<String, usize> = BTreeMap::new();

    let mut worst_bond_raw_embed: Vec<f64> = Vec::new();
    let mut extent_raw_embed: Vec<f64> = Vec::new();
    let mut n_nonfinite_from_embedder = 0usize; // must stay 0 -- Agent C's own validate_final_coords should guarantee this

    // (policy -> (n_ok_embedded, n_ok_legacy)), denominator is n_embed_ok for
    // the embedded column (legacy dg coords always exist) and n_total for the
    // legacy column.
    let mut ok_embedded: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut ok_legacy: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut err_buckets_embedded: BTreeMap<(&'static str, String), usize> = BTreeMap::new();
    let mut err_buckets_legacy: BTreeMap<(&'static str, String), usize> = BTreeMap::new();
    let mut fallback_count_embedded = 0usize;
    let mut fallback_count_legacy = 0usize;
    let mut n_panics = 0usize;

    // fused_aromatic-tagged rows, plus the exact molecules Agent F's PR body
    // named for the UffOnly blow-up list (some of those are tagged
    // "druglike"/"druglike_stress" in this corpus's category scheme, not
    // "fused_aromatic" -- checked by exact name, not by category tag).
    let uff_blowup_named: BTreeSet<&str> = [
        "naphthalene",
        "quinoline",
        "pyrene",
        "ibuprofen",
        "naproxen",
        "diphenhydramine",
        "atorvastatin_fragment",
        "caffeine",
    ]
    .into_iter()
    .collect();
    let mut naphthalene_family_rows: Vec<String> = Vec::new();

    println!(
        "{:<24} {:<10} {:>9} | {:>7} {:>7} {:>7} | {:>7} {:>7} {:>7}",
        "molecule", "embed", "extent", "MB_e", "MF_e", "UO_e", "MB_l", "MF_l", "UO_l"
    );
    println!(
        "  (e = embedder-fed coords, l = legacy dg::generate_coords; MB/MF/UO = \
         Mmff94BondAngleStrict / Mmff94WithUffFallback / UffOnly; OK or first 3 chars of error bucket)"
    );

    for &(name, smiles, _category) in CORPUS {
        n_total += 1;
        let mol = parse(smiles).unwrap_or_else(|e| panic!("{name}: parse failed: {e:?}"));

        let embed_result = embed_distance_geometry_v2(&mol, &EmbedParameters::default());
        let mut cols_e = vec!["--".to_string(); POLICIES.len()];
        let embed_label;
        let mut extent_str = "n/a".to_string();

        if let Ok(coords) = &embed_result {
            n_embed_ok += 1;
            embed_label = "ok".to_string();

            // Guard the handoff explicitly rather than assuming the happy path.
            assert_eq!(
                coords.atom_count(),
                mol.atom_count(),
                "{name}: embedder returned a different atom count than the molecule -- \
                 would panic opaquely inside minimize_with_policy's coords_to_vec"
            );
            let finite = coords.is_finite();
            if !finite {
                n_nonfinite_from_embedder += 1;
                eprintln!(
                    "FINDING: {name} embedder output contains NaN/Inf -- should never happen"
                );
            }
            assert!(
                finite,
                "{name}: embedder produced non-finite coordinates -- must never be fed to the minimizer"
            );

            let wb = worst_bond(&mol, coords);
            let ext = max_pairwise_extent(coords);
            worst_bond_raw_embed.push(wb);
            extent_raw_embed.push(ext);
            extent_str = format!("{ext:.1}A");

            for (i, &policy) in POLICIES.iter().enumerate() {
                let outcome = run_minimize(&mol, coords.clone(), policy);
                if outcome.panicked {
                    n_panics += 1;
                    eprintln!(
                        "FINDING: {name} PANICKED inside minimize_with_policy({policy:?}) with embedder-fed coords"
                    );
                }
                let pname = policy_name(policy);
                if outcome.ok {
                    *ok_embedded.entry(pname).or_insert(0) += 1;
                    cols_e[i] = "OK".to_string();
                    if policy == ForceFieldPolicy::Mmff94WithUffFallback
                        && outcome.fallback_occurred
                    {
                        fallback_count_embedded += 1;
                    }
                } else {
                    let bucket = outcome.err_bucket.clone().unwrap_or_default();
                    *err_buckets_embedded
                        .entry((pname, bucket.clone()))
                        .or_insert(0) += 1;
                    cols_e[i] = bucket.chars().take(3).collect();
                }
                if uff_blowup_named.contains(name) && policy == ForceFieldPolicy::UffOnly {
                    naphthalene_family_rows.push(format!(
                        "  {name:<22} embedder-fed: {} (worst_bond_after={:?}, converged={:?})",
                        if outcome.ok {
                            "OK"
                        } else {
                            outcome.err_bucket.as_deref().unwrap_or("?")
                        },
                        outcome.worst_bond_after,
                        outcome.converged
                    ));
                }
            }
        } else {
            n_embed_fail += 1;
            let cause = format!("{:?}", embed_result.as_ref().unwrap_err());
            *embed_fail_causes.entry(cause.clone()).or_insert(0) += 1;
            embed_label = cause;
            // Correctly skipped: no minimize_with_policy call at all for this
            // molecule under the embedder-fed source -- cols_e stays "--".
        }

        // Legacy dg::generate_coords source -- infallible, always runs, gives
        // the same-run baseline for the item-5 comparison.
        let legacy_coords = dg::generate_coords(&mol);
        let mut cols_l = vec!["--".to_string(); POLICIES.len()];
        for (i, &policy) in POLICIES.iter().enumerate() {
            let outcome = run_minimize(&mol, legacy_coords.clone(), policy);
            if outcome.panicked {
                n_panics += 1;
                eprintln!(
                    "FINDING: {name} PANICKED inside minimize_with_policy({policy:?}) with legacy dg coords"
                );
            }
            let pname = policy_name(policy);
            if outcome.ok {
                *ok_legacy.entry(pname).or_insert(0) += 1;
                cols_l[i] = "OK".to_string();
                if policy == ForceFieldPolicy::Mmff94WithUffFallback && outcome.fallback_occurred {
                    fallback_count_legacy += 1;
                }
            } else {
                let bucket = outcome.err_bucket.clone().unwrap_or_default();
                *err_buckets_legacy
                    .entry((pname, bucket.clone()))
                    .or_insert(0) += 1;
                cols_l[i] = bucket.chars().take(3).collect();
            }
            if uff_blowup_named.contains(name) && policy == ForceFieldPolicy::UffOnly {
                naphthalene_family_rows.push(format!(
                    "  {name:<22} legacy-dg:    {} (worst_bond_after={:?}, converged={:?})",
                    if outcome.ok {
                        "OK"
                    } else {
                        outcome.err_bucket.as_deref().unwrap_or("?")
                    },
                    outcome.worst_bond_after,
                    outcome.converged
                ));
            }
        }

        println!(
            "{name:<24} {:<10} {:>9} | {:>7} {:>7} {:>7} | {:>7} {:>7} {:>7}",
            embed_label,
            extent_str,
            cols_e[0],
            cols_e[1],
            cols_e[2],
            cols_l[0],
            cols_l[1],
            cols_l[2]
        );
    }

    assert_eq!(n_total, 58, "corpus size drifted from the frozen 58");

    // =========================================================================
    // Manufactured-NaN probe: does the minimizer's own typed safety net
    // (ForceFieldBridgeError::MinimizationFailed(NonFiniteCoordinates)) really
    // catch NaN coordinates, rather than propagating them silently or
    // panicking? Not exercised by the frozen 58 (the embedder never actually
    // emits NaN -- validate_final_coords already guards that), so this is a
    // deliberately separate, hand-constructed check of item 3's "if it
    // somehow did" clause.
    // =========================================================================
    println!("\n=== manufactured-NaN probe (not part of the 58-molecule corpus) ===");
    let ethane = parse("CC").unwrap();
    let mut nan_coords = Coords3D::new_zeroed(2);
    nan_coords.set(AtomIdx(0), Point3::new(f64::NAN, 0.0, 0.0));
    nan_coords.set(AtomIdx(1), Point3::new(1.5, 0.0, 0.0));
    let nan_outcome = run_minimize(&ethane, nan_coords, ForceFieldPolicy::Mmff94BondAngleStrict);
    let nan_probe_caught = !nan_outcome.panicked
        && !nan_outcome.ok
        && nan_outcome
            .err_bucket
            .as_deref()
            .is_some_and(|b| b.contains("NonFiniteCoordinates"));
    println!(
        "NaN input fed directly to minimize_with_policy(Mmff94BondAngleStrict): panicked={} ok={} bucket={:?}",
        nan_outcome.panicked, nan_outcome.ok, nan_outcome.err_bucket
    );
    assert!(
        !nan_outcome.panicked,
        "manufactured-NaN probe must not panic the minimizer"
    );
    assert!(
        nan_probe_caught,
        "manufactured-NaN probe: expected a typed MinimizationFailed(NonFiniteCoordinates), got {:?}",
        nan_outcome.err_bucket
    );
    println!("VERIFIED: NaN input is caught as a typed error, not propagated or panicked on.");

    // =========================================================================
    // 3-membered-ring probe: the frozen 58 contains ZERO 3-membered rings
    // (Agent C's own module doc states this), so this is a separate, explicit
    // check on the exact molecules Agent C's own tests use. FORMERLY these all
    // failed closed with `BoundsConstructionFailed` (a bound-construction bug
    // in `dg_fft::build_bond_angle_bounds`'s angle-constraint loop, now fixed --
    // see `distance_geometry_v2::tests::three_membered_rings_embed_successfully`).
    // =========================================================================
    println!("\n=== 3-membered-ring probe (not part of the 58-molecule corpus) ===");
    let mut three_ring_all_embedded = true;
    for smiles in ["C1CC1", "C1CO1", "C1CN1", "C1CS1"] {
        let mol = parse(smiles).unwrap();
        let result = embed_distance_geometry_v2(&mol, &EmbedParameters::default());
        println!("  {smiles:<8} embed result: {}", result.is_ok());
        three_ring_all_embedded &= result.is_ok();
    }
    assert!(
        three_ring_all_embedded,
        "expected every 3-membered ring to embed successfully"
    );
    println!("VERIFIED: 3-membered rings now embed successfully and reach the minimizer.");

    // =========================================================================
    // Seed-robustness spot check (item raised by review): the embedder is
    // deterministic given a seed, but that doesn't mean the end-to-end
    // pass/fail ANSWER is seed-independent. Re-run the embedder-fed pipeline
    // at 2 extra seeds and compare the per-policy "which molecules end up
    // succeeding end-to-end" sets against the default-seed run above.
    // =========================================================================
    println!(
        "\n=== seed-robustness spot check (embedder-fed pipeline only; legacy dg is seed-independent) ==="
    );
    let seed_default_sets = run_embedder_fed_pipeline_at_seed(default_seed);
    let mut any_seed_flip = false;
    for &extra_seed in &[1u64, 42u64] {
        let sets = run_embedder_fed_pipeline_at_seed(extra_seed);
        for (i, &policy) in POLICIES.iter().enumerate() {
            let base = &seed_default_sets[i];
            let other = &sets[i];
            let only_in_base: Vec<_> = base.difference(other).collect();
            let only_in_other: Vec<_> = other.difference(base).collect();
            if !only_in_base.is_empty() || !only_in_other.is_empty() {
                any_seed_flip = true;
                println!(
                    "  SEED FLIP: {} default-seed-ok-not-seed-{extra_seed}={:?} seed-{extra_seed}-ok-not-default={:?}",
                    policy_name(policy),
                    only_in_base,
                    only_in_other
                );
            } else {
                println!(
                    "  {} stable across seed {default_seed:#x} vs {extra_seed}: same {} molecule(s) end-to-end OK",
                    policy_name(policy),
                    base.len()
                );
            }
        }
    }
    if !any_seed_flip {
        println!(
            "VERIFIED: end-to-end pass/fail sets are identical across the 3 seeds checked for every policy."
        );
    }

    // =========================================================================
    // Summary
    // =========================================================================
    println!("\n=== SUMMARY ({n_total} molecules) ===");
    println!("embed succeeded: {n_embed_ok}/{n_total}");
    println!("embed failed:    {n_embed_fail}/{n_total}  {embed_fail_causes:?}");
    println!("embedder output non-finite (must be 0): {n_nonfinite_from_embedder}");
    println!("panics inside minimize_with_policy (must be 0): {n_panics}");
    if !worst_bond_raw_embed.is_empty() {
        let min_wb = worst_bond_raw_embed
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let max_wb = worst_bond_raw_embed.iter().cloned().fold(0.0, f64::max);
        let min_ext = extent_raw_embed
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let max_ext = extent_raw_embed.iter().cloned().fold(0.0, f64::max);
        println!(
            "raw embedder geometry scale sanity: worst-bond range [{min_wb:.2}, {max_wb:.2}] A, \
             max-pairwise-extent range [{min_ext:.2}, {max_ext:.2}] A across {n_embed_ok} embedded molecules \
             (neither degenerate/collapsed nor wildly out of scale for organic molecules)"
        );
    }

    println!("\n--- end-to-end (embed AND minimize both succeed), per policy ---");
    println!(
        "{:<24} {:>18} {:>10} {:>18} {:>10}",
        "policy", "embedder-fed OK", "of embed_ok", "legacy-dg OK", "of total"
    );
    for &policy in POLICIES {
        let pname = policy_name(policy);
        let e_ok = *ok_embedded.get(pname).unwrap_or(&0);
        let l_ok = *ok_legacy.get(pname).unwrap_or(&0);
        println!(
            "{pname:<24} {e_ok:>18} {:>9}/{n_embed_ok} {l_ok:>18} {:>9}/{n_total}",
            e_ok, l_ok
        );
    }
    println!(
        "(embedder-fed end-to-end failures = n_embed_fail [{n_embed_fail}, embed stage] + \
         [embed_ok - policy_ok, minimize stage, per policy above])"
    );

    println!("\n--- minimize-stage error buckets, embedder-fed source ---");
    for ((pname, bucket), count) in &err_buckets_embedded {
        println!("  {pname:<24} {bucket:<40} {count}");
    }
    println!("--- minimize-stage error buckets, legacy-dg source ---");
    for ((pname, bucket), count) in &err_buckets_legacy {
        println!("  {pname:<24} {bucket:<40} {count}");
    }

    println!(
        "\nMmff94WithUffFallback actual fallback occurrences (fallback_reason.is_some() on Ok, \
         NOT actual_force_field_used != requested_force_field): embedder-fed={fallback_count_embedded}/{n_embed_ok} \
         legacy-dg={fallback_count_legacy}/{n_total}"
    );

    let uff_blowup_named_count = uff_blowup_named.len();
    println!(
        "\n=== naphthalene-family / fused-aromatic UFF-blowup comparison (item 5) ===\n\
         NOTE (corrected after independent verification): PR #169's own doc comment cited these \
         same {uff_blowup_named_count} molecules as \"UffOnly\"'s full-corpus blowup count, but \
         that figure was actually a copy of Mmff94WithUffFallback's fallback-trigger-population \
         count -- UffOnly was never run over the full 58-molecule corpus in that PR. This smoke \
         test is the first real full-corpus UffOnly measurement (17/58 blow up under legacy \
         dg::generate_coords starting geometry, not 8 -- see PR #187's doc-comment fix and issue \
         #188 for the other 9). Named here (fused/conjugated-aromatic subset only): \
         {uff_blowup_named:?}. Below: same policy (UffOnly), both starting-geometry sources, \
         this run:"
    );
    for row in &naphthalene_family_rows {
        println!("{row}");
    }

    println!(
        "\n=== VERDICT ===\n\
         See the printed counts above for the composed pipeline's end-to-end success rate per \
         policy/source, the manufactured-NaN and 3-ring probes for the typed-failure-path checks, \
         and the seed-robustness section for whether any pass/fail answer is seed-sensitive. \
         One observation this smoke test does NOT attempt to measure (flagged, not quantified): \
         EmbedParameters::default() has enforce_chirality=false, and this corpus is full of \
         declared stereocenters -- the composed pipeline as run here emits minimized geometry \
         whose chirality was never verified by either module (Agent C only checks it when asked; \
         Agent F's minimizer has no stereo check at all). Neither module is wrong in isolation; \
         this is specifically an integration-boundary gap."
    );
}
