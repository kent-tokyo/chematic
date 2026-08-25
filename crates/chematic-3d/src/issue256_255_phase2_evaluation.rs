//! Issue #256/#255 Phase 2: differential evaluation of `generate_coords_connectivity_ordered`
//! (now [`crate::dg_connectivity_ordered::generate_coords_connectivity_ordered`], `pub`)
//! against legacy `generate_coords` (RFC `docs/rfcs/dg_connectivity_ordered_placement_rfc.md`
//! §5's Phase 2, spec given directly by the user 2026-08-25).
//!
//! Originally written while the engine was still `pub(crate)` (it had a
//! known, unfixed new-island ring-entry-direction regression at the time and
//! hadn't cleared Phase 3's go/no-go criteria), hence an ignored test rather
//! than an example crate. The engine has since cleared those criteria (fixed
//! by the new-island fix this same harness was used to measure) and moved to
//! its own public module, `crate::dg_connectivity_ordered` -- kept as an
//! internal ignored test regardless, since that's still the fastest way to
//! re-run this specific differential comparison. Run manually with:
//! `cargo test -p chematic-3d issue256_255_phase2_differential_evaluation --release -- --ignored --nocapture`
//!
//! Compares, per molecule, four states: legacy raw, new-engine raw, legacy →
//! `ForceFieldPolicy::UffOnly`, new-engine → `ForceFieldPolicy::UffOnly` (same
//! post-processing applied to both engines' output, so any difference in the
//! minimized result traces back to the starting geometry, not to a different
//! force-field path). `generate_coords` itself is not modified anywhere in
//! this crate by this module -- purely a read-only measurement.
//!
//! Population: #277's own 17-molecule regression set
//! (`chembl_tier_b_0003/_0008/_0010/_0021/_0025/_0026/_0027/_0035/_0065/_0066/
//! _0095/_0143/_0146/_0147/_0153/_0156/_0157`, SMILES pulled directly from
//! `validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json`, the same
//! source the RFC itself cites) + the RFC's own 8 known-broken topologies +
//! anthracene (Phase 1's extra coverage) + spiro/biphenyl/terphenyl/meta-linked
//! positive controls + a ring-with-tail-substituent case + a pure chain + a
//! disconnected-component case. Three of the 17 chembl molecules
//! (_0143/_0156/_0157) already carry declared `@`/`@@` stereocenters, so no
//! separate "stereo-bearing control" fixture was added -- see the module doc
//! on `dg.rs` for why both engines are stereo-blind by design either way.
//!
//! Required metrics per state (as specified): worst bond-length ratio
//! (actual/ideal for whichever bond deviates most), 15%/50% bond-length
//! violation rate, non-bonded gross clash count, min pairwise distance,
//! finite/NaN, stereo violation count, `sound` (finite + atom count unchanged +
//! worst bond length within a sane absolute bound). Also: UFF success rate (for
//! the two minimized states), determinism (new engine run twice per molecule,
//! exact-equality expected -- it's a deterministic, non-RNG algorithm), wall
//! time, whether #277's currently-sound population stays sound, and whether
//! any previously-sound molecule newly breaks.
//!
//! Separately: an atom-order permutation-robustness check on 5 representative
//! molecules, each run under 2 non-identity relabelings. Per
//! `feedback_permutation_invariance_test_template` (this project's own prior
//! lesson): the correct invariant here is that AGGREGATE QUALITY METRICS stay
//! consistent across relabelings, NOT that coordinates are literally equal --
//! `seed_ring_system_index`'s lowest-`AtomIdx` tie-break means the new
//! engine's seed choice (and everything built from it) can legitimately differ
//! under a relabeling.
//!
use std::time::Instant;

use chematic_core::{AtomIdx, Molecule, MoleculeBuilder};
use chematic_smiles::parse;

use crate::Coords3D;
use crate::dg::generate_coords;
use crate::dg_connectivity_ordered::generate_coords_connectivity_ordered;
use crate::dg_fft::ideal_bond_length;
use crate::minimize::{
    ForceFieldPolicy, MAX_SANE_BOND_LENGTH, MinimizeConfig, minimize_with_policy,
};
use crate::stereo_constraints::verify_stereo;

fn population() -> Vec<(&'static str, &'static str)> {
    let mut v: Vec<(&'static str, &'static str)> = vec![
        // #277's 17-molecule regression population (SMILES from
        // validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json).
        (
            "chembl_tier_b_0003",
            "Cc1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2)cc1",
        ),
        (
            "chembl_tier_b_0008",
            "CS(=O)(=O)c1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2Cl)cc1",
        ),
        (
            "chembl_tier_b_0010",
            "CSc1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2Cl)cc1",
        ),
        (
            "chembl_tier_b_0021",
            "COc1ccc(C(=O)c2ccc(-n3ncc(=O)[nH]c3=O)cc2)cc1",
        ),
        (
            "chembl_tier_b_0025",
            "Nc1cc[n+](Cc2ccc(CCc3ccc(C[n+]4ccc(N)c5ccccc54)cc3)cc2)c2ccccc12",
        ),
        (
            "chembl_tier_b_0026",
            r"Nc1cc[n+](Cc2ccc(/C=C\c3ccc(C[n+]4ccc(N)c5ccccc54)cc3)cc2)c2ccccc12",
        ),
        (
            "chembl_tier_b_0027",
            "c1ccc(C[n+]2ccc(NCCCCCCCCCCNc3cc[n+](Cc4ccccc4)c4ccccc34)c3ccccc32)cc1",
        ),
        (
            "chembl_tier_b_0035",
            "Nc1cc[n+](Cc2ccc(Cc3ccc(C[n+]4ccc(N)c5ccccc54)cc3)cc2)c2ccccc12",
        ),
        ("chembl_tier_b_0065", "CN1CCCCC1CN1CCN(Cc2ccncc2)CC1"),
        ("chembl_tier_b_0066", "CN1CCCCC1CN1CCN(Cc2cccnc2)CC1"),
        ("chembl_tier_b_0095", "O=C(C1CCCCN1)N1CCN(Cc2cccnc2)CC1"),
        (
            "chembl_tier_b_0143",
            "Cc1nccn1CCCCc1ccc(CC(=O)N[C@H](CO)Cc2ccc(OCCC3CCCCC3)c(CCCN)c2)cc1",
        ),
        ("chembl_tier_b_0146", "O=C(C1CCCCN1)N1CCN(Cc2ccncc2)CC1"),
        ("chembl_tier_b_0147", "O=C(C1CCCCN1)N1CCN(Cc2ccccn2)CC1"),
        (
            "chembl_tier_b_0153",
            "CCC(CC)C(=O)OCC1CN(Cc2cc(OC)c(OC)c(OC)c2)CCN1Cc1cc(OC)c(OC)c(OC)c1",
        ),
        (
            "chembl_tier_b_0156",
            "Cc1nccn1CCCCc1ccc(CC(=O)N[C@H](CO)Cc2ccc(OCCC3CCCCC3)c(CCCCN)c2)cc1",
        ),
        (
            "chembl_tier_b_0157",
            "Cc1nccn1CCCCc1ccc(CC(=O)N[C@@H](CO)Cc2ccc(OCCC3CCCCC3)c(CCCN)c2)cc1",
        ),
        // RFC's own 8 known-broken topologies.
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("quinoline", "c1ccc2ncccc2c1"),
        ("phenanthrene", "c1ccc2c(c1)ccc1ccccc12"),
        ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34"),
        ("diphenylmethane_bridge1", "c1ccccc1Cc1ccccc1"),
        ("bibenzyl_bridge2", "c1ccccc1CCc1ccccc1"),
        ("diphenylpropane_bridge3", "c1ccccc1CCCc1ccccc1"),
        ("diphenylbutane_bridge4", "c1ccccc1CCCCc1ccccc1"),
        // Extra coverage / positive controls.
        ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
        ("spiro_5_5_undecane", "C1CCC2(CC1)CCCCC2"),
        ("biphenyl", "c1ccc(cc1)-c1ccccc1"),
        ("terphenyl", "c1ccc(cc1)-c1ccc(cc1)-c1ccccc1"),
        ("3_phenylpyridine_meta", "c1ccc(cc1)-c1cccnc1"),
        ("ibuprofen_ring_with_tail", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("pentane_pure_chain", "CCCCC"),
        ("disconnected_ethanol_benzene", "CCO.c1ccccc1"),
    ];
    v.sort_by_key(|&(name, _)| name);
    v
}

#[derive(Clone, Copy, Debug, Default)]
struct GeometryMetrics {
    all_finite: bool,
    atom_count_unchanged: bool,
    worst_bond_length: f64,
    /// actual/ideal for whichever bond has the largest |actual/ideal - 1|.
    /// >1.0 means stretched, <1.0 means compressed.
    worst_bond_ratio: f64,
    bond_violation_rate_15pct: f64,
    bond_violation_rate_50pct: f64,
    gross_clash_count: usize,
    min_pairwise_distance: f64,
    stereo_violations: usize,
    /// Declared stereo elements `verify_stereo` couldn't evaluate at all
    /// (typically an implicit-H ring-fused center -- see issue #291's
    /// `phantom_neighbor_position` history). If this equals `stereo_declared`
    /// for a molecule, `stereo_violations == 0` means "nothing was checked,"
    /// not "stereo is correct" -- both must be read together.
    stereo_unevaluable: usize,
    stereo_declared: usize,
    sound: bool,
}

const NONBONDED_CLASH_THRESHOLD_ANGSTROM: f64 = 1.2;

fn compute_metrics(mol: &Molecule, coords: &Coords3D) -> GeometryMetrics {
    let n = mol.atom_count();
    let all_finite = coords.is_finite();
    let atom_count_unchanged = coords.atom_count() == n;

    let mut worst_bond_length = 0.0_f64;
    let mut worst_bond_ratio = 1.0_f64;
    let mut worst_dev = 0.0_f64;
    let mut total_bonds = 0usize;
    let mut violations_15 = 0usize;
    let mut violations_50 = 0usize;
    for (_, bond) in mol.bonds() {
        let actual = coords.get(bond.atom1).distance(&coords.get(bond.atom2));
        if actual > worst_bond_length {
            worst_bond_length = actual;
        }
        let ideal = ideal_bond_length(mol, bond.atom1, bond.atom2);
        if ideal > 0.0 && ideal.is_finite() && actual.is_finite() {
            total_bonds += 1;
            let dev = (actual - ideal).abs() / ideal;
            if dev > 0.15 {
                violations_15 += 1;
            }
            if dev > 0.50 {
                violations_50 += 1;
            }
            if dev > worst_dev {
                worst_dev = dev;
                worst_bond_ratio = actual / ideal;
            }
        }
    }
    let bond_violation_rate_15pct = if total_bonds == 0 {
        0.0
    } else {
        violations_15 as f64 / total_bonds as f64
    };
    let bond_violation_rate_50pct = if total_bonds == 0 {
        0.0
    } else {
        violations_50 as f64 / total_bonds as f64
    };

    let mut gross_clash_count = 0usize;
    let mut min_pairwise_distance = f64::MAX;
    for i in 0..n {
        for j in (i + 1)..n {
            let a = AtomIdx(i as u32);
            let b = AtomIdx(j as u32);
            let d = coords.get(a).distance(&coords.get(b));
            if d < min_pairwise_distance {
                min_pairwise_distance = d;
            }
            if mol.bond_between(a, b).is_none() && d < NONBONDED_CLASH_THRESHOLD_ANGSTROM {
                gross_clash_count += 1;
            }
        }
    }
    if n < 2 {
        min_pairwise_distance = f64::INFINITY;
    }

    let stereo = verify_stereo(mol, coords);
    let stereo_violations = stereo.n_violations();
    let stereo_unevaluable = stereo.n_unevaluable();
    let stereo_declared = stereo.n_declared();
    let sound = all_finite && atom_count_unchanged && worst_bond_length <= MAX_SANE_BOND_LENGTH;

    GeometryMetrics {
        all_finite,
        atom_count_unchanged,
        worst_bond_length,
        worst_bond_ratio,
        bond_violation_rate_15pct,
        bond_violation_rate_50pct,
        gross_clash_count,
        min_pairwise_distance,
        stereo_violations,
        stereo_unevaluable,
        stereo_declared,
        sound,
    }
}

/// Build a molecule with atoms relabeled through `perm` (`perm[new_idx] =
/// old_idx`), preserving every bond -- ring perception (`find_sssr`) is
/// recomputed fresh from the rebuilt graph topology inside both engines, so
/// this genuinely exercises a different `AtomIdx` order end to end, not just
/// a coordinate relabeling after the fact.
fn permute_molecule(mol: &Molecule, perm: &[usize]) -> Molecule {
    let mut old_to_new = vec![0u32; perm.len()];
    for (new_idx, &old_idx) in perm.iter().enumerate() {
        old_to_new[old_idx] = new_idx as u32;
    }
    let mut b = MoleculeBuilder::new();
    for &old_idx in perm {
        b.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
    }
    for (_, bond) in mol.bonds() {
        let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
        let c = AtomIdx(old_to_new[bond.atom2.0 as usize]);
        let _ = b.add_bond(a, c, bond.order);
    }
    b.build()
}

/// Reversal and a fixed-stride riffle -- both non-identity, both easy to hand
/// -verify, no new `rand`-crate dependency needed (this project's own
/// ladder-first convention: don't add a dependency for what a few lines can
/// do).
fn permutations_for(n: usize) -> Vec<Vec<usize>> {
    let reversed: Vec<usize> = (0..n).rev().collect();
    let mut riffle = Vec::with_capacity(n);
    let half = n.div_ceil(2);
    for i in 0..half {
        riffle.push(i);
        if half + i < n {
            riffle.push(half + i);
        }
    }
    vec![reversed, riffle]
}

#[test]
#[ignore = "manual Phase 2 differential evaluation -- cargo test -p chematic-3d issue256_255_phase2_differential_evaluation --release -- --ignored --nocapture"]
fn issue256_255_phase2_differential_evaluation() {
    let ffconfig = MinimizeConfig::default();

    println!(
        "name\tengine\tstate\tfinite\tatoms_ok\tworst_bond_len\tworst_bond_ratio\tviol15\tviol50\tclash\tmin_pair\tstereo_viol/unevaluable/declared\tsound\tuff_ok"
    );

    struct Row {
        name: &'static str,
        engine: &'static str,
        state: &'static str,
        m: Option<GeometryMetrics>,
        uff_ok: Option<bool>,
    }
    let mut rows: Vec<Row> = Vec::new();

    let mut wall_old = std::time::Duration::ZERO;
    let mut wall_new = std::time::Duration::ZERO;
    let mut determinism_mismatches: Vec<&'static str> = Vec::new();

    for (name, smiles) in population() {
        let mol = match parse(smiles) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("SKIP {name}: parse failed: {e:?}");
                continue;
            }
        };

        let t0 = Instant::now();
        let old_raw = generate_coords(&mol);
        wall_old += t0.elapsed();

        let t1 = Instant::now();
        let new_raw = generate_coords_connectivity_ordered(&mol);
        wall_new += t1.elapsed();

        let new_raw_again = generate_coords_connectivity_ordered(&mol);
        let identical = (0..mol.atom_count()).all(|i| {
            let a = new_raw.get(AtomIdx(i as u32));
            let b = new_raw_again.get(AtomIdx(i as u32));
            a.x == b.x && a.y == b.y && a.z == b.z
        });
        if !identical {
            determinism_mismatches.push(name);
        }

        let old_metrics = compute_metrics(&mol, &old_raw);
        let new_metrics = compute_metrics(&mol, &new_raw);

        let old_min = minimize_with_policy(&mol, old_raw, ForceFieldPolicy::UffOnly, &ffconfig);
        let new_min = minimize_with_policy(&mol, new_raw, ForceFieldPolicy::UffOnly, &ffconfig);

        let old_min_metrics = old_min
            .as_ref()
            .ok()
            .map(|r| compute_metrics(&mol, &r.coords));
        let new_min_metrics = new_min
            .as_ref()
            .ok()
            .map(|r| compute_metrics(&mol, &r.coords));

        for (engine, state, m, uff_ok) in [
            ("legacy", "raw", Some(old_metrics), None),
            ("new", "raw", Some(new_metrics), None),
            ("legacy", "uff_only", old_min_metrics, Some(old_min.is_ok())),
            ("new", "uff_only", new_min_metrics, Some(new_min.is_ok())),
        ] {
            if let Some(mm) = m {
                println!(
                    "{name}\t{engine}\t{state}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{:.4}\t{}/{}/{}\t{}\t{}",
                    mm.all_finite,
                    mm.atom_count_unchanged,
                    mm.worst_bond_length,
                    mm.worst_bond_ratio,
                    mm.bond_violation_rate_15pct,
                    mm.bond_violation_rate_50pct,
                    mm.gross_clash_count,
                    mm.min_pairwise_distance,
                    mm.stereo_violations,
                    mm.stereo_unevaluable,
                    mm.stereo_declared,
                    mm.sound,
                    uff_ok.map(|b| b.to_string()).unwrap_or_default(),
                );
            } else {
                println!(
                    "{name}\t{engine}\t{state}\t-\t-\t-\t-\t-\t-\t-\t-\t-/-/-\tfalse\t{}",
                    uff_ok.map(|b| b.to_string()).unwrap_or_default(),
                );
            }
            rows.push(Row {
                name,
                engine,
                state,
                m,
                uff_ok,
            });
        }
    }

    // ---- Summary ----
    eprintln!("\n=== Summary ===");
    eprintln!(
        "Raw-geometry wall time: legacy={:.3}s new={:.3}s",
        wall_old.as_secs_f64(),
        wall_new.as_secs_f64()
    );
    eprintln!(
        "Determinism (new engine, run twice per molecule): {}/{} identical{}",
        population().len() - determinism_mismatches.len(),
        population().len(),
        if determinism_mismatches.is_empty() {
            String::new()
        } else {
            format!(", mismatches: {determinism_mismatches:?}")
        }
    );

    for state in ["raw", "uff_only"] {
        for engine in ["legacy", "new"] {
            let sound_count = rows
                .iter()
                .filter(|r| r.state == state && r.engine == engine)
                .filter(|r| r.m.map(|m| m.sound).unwrap_or(false))
                .count();
            let total = rows
                .iter()
                .filter(|r| r.state == state && r.engine == engine)
                .count();
            eprintln!("{engine}/{state}: sound {sound_count}/{total}");
        }
    }

    // Fine-grained aggregation -- NOT just sound/uff_ok. `sound` is a coarse
    // pass/fail gate; the RFC's own worry (and the reason this section
    // exists) is that a strictly-better raw geometry can still settle into a
    // WORSE post-minimization local optimum. Only measuring sound/uff_ok
    // would silently miss that -- both engines can pass the coarse gate while
    // differing materially on every metric below.
    eprintln!("\n--- Per-engine/state aggregates (successfully-computed rows only) ---");
    for state in ["raw", "uff_only"] {
        for engine in ["legacy", "new"] {
            let metrics: Vec<GeometryMetrics> = rows
                .iter()
                .filter(|r| r.state == state && r.engine == engine)
                .filter_map(|r| r.m)
                .collect();
            let n = metrics.len();
            if n == 0 {
                eprintln!("{engine}/{state}: no successful rows");
                continue;
            }
            let clash_sum: usize = metrics.iter().map(|m| m.gross_clash_count).sum();
            let mean_viol15 = metrics
                .iter()
                .map(|m| m.bond_violation_rate_15pct)
                .sum::<f64>()
                / n as f64;
            let mean_viol50 = metrics
                .iter()
                .map(|m| m.bond_violation_rate_50pct)
                .sum::<f64>()
                / n as f64;
            let worst_min_pair = metrics
                .iter()
                .map(|m| m.min_pairwise_distance)
                .fold(f64::MAX, f64::min);
            let worst_bond_len = metrics
                .iter()
                .map(|m| m.worst_bond_length)
                .fold(0.0_f64, f64::max);
            let max_worst_bond_ratio_dev = metrics
                .iter()
                .map(|m| (m.worst_bond_ratio - 1.0).abs())
                .fold(0.0_f64, f64::max);
            eprintln!(
                "{engine}/{state}: n={n} clash_sum={clash_sum} mean_viol15={mean_viol15:.4} \
                 mean_viol50={mean_viol50:.4} worst_min_pair={worst_min_pair:.4} \
                 worst_bond_len={worst_bond_len:.4} max_bond_ratio_dev={max_worst_bond_ratio_dev:.4}"
            );
        }
    }

    for engine in ["legacy", "new"] {
        let ok_count = rows
            .iter()
            .filter(|r| r.state == "uff_only" && r.engine == engine)
            .filter(|r| r.uff_ok == Some(true))
            .count();
        let total = rows
            .iter()
            .filter(|r| r.state == "uff_only" && r.engine == engine)
            .count();
        eprintln!("{engine}/uff_only: UFF Ok {ok_count}/{total}");
    }

    // Which molecules regressed (legacy sound -> new unsound) or improved
    // (legacy unsound -> new sound), per state.
    for state in ["raw", "uff_only"] {
        let mut regressed = Vec::new();
        let mut improved = Vec::new();
        for name in population().iter().map(|&(n, _)| n) {
            let legacy_sound = rows
                .iter()
                .find(|r| r.name == name && r.engine == "legacy" && r.state == state)
                .and_then(|r| r.m)
                .map(|m| m.sound)
                .unwrap_or(false);
            let new_sound = rows
                .iter()
                .find(|r| r.name == name && r.engine == "new" && r.state == state)
                .and_then(|r| r.m)
                .map(|m| m.sound)
                .unwrap_or(false);
            if legacy_sound && !new_sound {
                regressed.push(name);
            }
            if !legacy_sound && new_sound {
                improved.push(name);
            }
        }
        eprintln!(
            "{state}: regressed(sound->unsound)={regressed:?} improved(unsound->sound)={improved:?}"
        );
    }

    // ---- Atom-order permutation robustness ----
    eprintln!("\n=== Atom-order permutation robustness ===");
    let perm_targets = [
        "naphthalene",
        "phenanthrene",
        "chembl_tier_b_0003",
        "ibuprofen_ring_with_tail",
        "biphenyl",
    ];
    for target in perm_targets {
        let Some(&(_, smiles)) = population().iter().find(|&&(n, _)| n == target) else {
            eprintln!("{target}: not found in population, skipping");
            continue;
        };
        let mol = parse(smiles).expect("perm target parses");
        let n = mol.atom_count();
        let base_new = compute_metrics(&mol, &generate_coords_connectivity_ordered(&mol));
        eprintln!(
            "{target} (n={n}): base sound={} worst_bond_ratio={:.4} clash={}",
            base_new.sound, base_new.worst_bond_ratio, base_new.gross_clash_count
        );
        for perm in permutations_for(n) {
            let permuted = permute_molecule(&mol, &perm);
            let m = compute_metrics(&permuted, &generate_coords_connectivity_ordered(&permuted));
            // The one real gate: soundness must survive relabeling -- this
            // CAN fail (unlike a symmetric-metric-vs-symmetric-metric
            // comparison, which is always trivially true). Everything else
            // is reported as a raw delta, not smoothed into a synthetic
            // pass/fail verdict this harness has no principled threshold for.
            let gate = if m.sound { "SOUND" } else { "**UNSOUND**" };
            let clash_delta = m.gross_clash_count as i64 - base_new.gross_clash_count as i64;
            let bond_ratio_dev_delta =
                (m.worst_bond_ratio - 1.0).abs() - (base_new.worst_bond_ratio - 1.0).abs();
            eprintln!(
                "  perm: {gate}  worst_bond_ratio={:.4} (base {:.4})  clash={} (delta {:+})  bond_ratio_dev_delta={bond_ratio_dev_delta:+.4}",
                m.worst_bond_ratio, base_new.worst_bond_ratio, m.gross_clash_count, clash_delta,
            );
        }
    }
}
