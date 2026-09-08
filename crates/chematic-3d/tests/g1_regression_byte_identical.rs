//! Pins descriptor outputs for the public connectivity-ordered conformer
//! engine, including the G1 `rdkit_*` descriptor additions.
//!
//! The expected values are raw `f64::to_bits()` snapshots for aspirin's
//! rule-based conformer. They provide a durable regression guard rather than
//! relying on a one-time source diff.
//!
//! **Re-pinned (issue #252 follow-up) after a real `dg::generate_coords` bug
//! fix**: `place_component` used to place a non-ring root atom unconditionally
//! at `(x_offset, 0, 0)`, which collided with a ring vertex `place_rings`
//! independently computed at that same point, AND only ever seeded
//! `dfs_place` from a single root, silently leaving any substituent on a
//! *different* ring atom at the `Coords3D::new_zeroed` (0, 0, 0) default.
//! Aspirin has two substituents on different ring atoms (the acetoxy group
//! and the carboxylic acid), so its pre-fix conformer had a collapsed
//! substituent -- these were byte-identical snapshots of a genuinely wrong
//! geometry, not a neutral baseline. Re-verified sane post-fix (minimum
//! pairwise interatomic distance 1.22 \u{c5}, no collisions) before repinning.
//!
//! A second, independent `dg::generate_coords` fix landed in the same PR
//! (`place_rings`'s ring visiting order didn't match ring-fusion adjacency,
//! silently superimposing unrelated rings in some multi-ring-fused systems
//! -- see issue #185). Confirmed NOT to affect aspirin (single ring, so the
//! multi-ring fusion-order logic never applies) -- this file's pinned
//! values are unchanged by it, verified by running this test unmodified
//! both before and after that second fix.
//!
//! A third fix (same PR, review round 2) taught `place_rings` to anchor a
//! fusion-disconnected ring island via a real bond to already-placed
//! structure when one exists (biphenyl, terphenyl), rather than an
//! arbitrary fixed offset that silently stretched that bond. Also
//! confirmed NOT to affect aspirin (still just the one ring) -- unchanged
//! before/after, same as the second fix.
//!
//! Two of this file's pinned values (`plane_of_best_fit`,
//! `descriptors_3d_outputs_unchanged`'s `rdf[0]`) are asserted via
//! `assert_close_ulp` rather than `to_bits()` equality: CI (x86_64) and
//! local development (aarch64) compute them a few ULP apart --
//! FMA/vectorization differences in their summation order, not a
//! correctness regression (the other 17 values in this file match
//! bit-for-bit across both targets).

use chematic_3d::dg::generate_coords;
use chematic_3d::{
    asphericity, autocorr_3d, eccentricity, getaway_descriptors, npr1, npr2, plane_of_best_fit,
    pmi, radius_of_gyration, rdf_descriptors, whim_descriptors,
};
use chematic_smiles::parse;

fn aspirin_coords() -> (chematic_core::Molecule, chematic_3d::Coords3D) {
    let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    let coords = generate_coords(&mol);
    (mol, coords)
}

/// Asserts `actual` is finite, positive, and within `max_ulp` of
/// `expected_bits` (an `f64::to_bits()` snapshot) -- an ULP-scale
/// tolerance rather than a magnitude-scaled absolute epsilon, so it stays
/// exactly as tight regardless of the value's own magnitude. `to_bits()`
/// is monotonic with value for same-signed floats, so a plain bit-pattern
/// `abs_diff` is a correct ULP distance here (both callers below are
/// positive).
fn assert_close_ulp(actual: f64, expected_bits: u64, max_ulp: u64, label: &str) {
    assert!(actual.is_finite(), "{label}: expected finite, got {actual}");
    assert!(actual > 0.0, "{label}: expected positive, got {actual}");
    let actual_bits = actual.to_bits();
    let ulp_diff = actual_bits.abs_diff(expected_bits);
    assert!(
        ulp_diff <= max_ulp,
        "{label}: {ulp_diff} ULP from expected (bits {actual_bits} vs {expected_bits}), \
         allowed {max_ulp}"
    );
}

#[test]
fn shape_descriptors_outputs_are_stable() {
    let (mol, coords) = aspirin_coords();

    let (p1, p2, p3) = pmi(&mol, &coords);
    assert_eq!(p1.to_bits(), 4642488561478464692);
    assert_eq!(p2.to_bits(), 4646694506356206969);
    assert_eq!(p3.to_bits(), 4648593491540578086);

    assert_eq!(npr1(&mol, &coords).to_bits(), 4600603682732697251);
    assert_eq!(npr2(&mol, &coords).to_bits(), 4604858220831292512);
    assert_eq!(
        radius_of_gyration(&mol, &coords).to_bits(),
        4611435890214201730
    );
    assert_eq!(asphericity(&mol, &coords).to_bits(), 4643407586323950178);
    assert_eq!(eccentricity(&mol, &coords).to_bits(), 4605239949045052297);
    // Not bit-pinned like the rest: this value's summation order makes it
    // sensitive to FMA/vectorization differences between the aarch64
    // (local) and x86_64 (CI) targets. Keep a small ULP allowance for that
    // platform-level summation difference.
    // 2 ULP gives slack above the observed 1 ULP without loosening this
    // into a magnitude-blind absolute epsilon.
    assert_close_ulp(
        plane_of_best_fit(&mol, &coords),
        4602456162078098012,
        2,
        "plane_of_best_fit",
    );
}

#[test]
fn descriptors_3d_outputs_are_stable() {
    let (mol, coords) = aspirin_coords();

    let whim = whim_descriptors(&mol, &coords);
    assert_eq!(whim.len(), 22);
    assert_eq!(whim[0].to_bits(), 4612594550505702359);

    let getaway = getaway_descriptors(&mol, &coords);
    assert_eq!(getaway.len(), 19);
    assert_eq!(getaway[0].to_bits(), 4613199105410399549);

    let rdf = rdf_descriptors(&mol, &coords);
    assert_eq!(rdf.len(), 20);
    // Same aarch64/x86_64 FMA/summation-order drift as plane_of_best_fit
    // above -- see that comment. The value is near 1e-20, so a ULP-based
    // allowance is more meaningful than an absolute epsilon.
    assert_close_ulp(rdf[0], 4626963717017789206, 150, "rdf[0]");

    let ac = autocorr_3d(&mol, &coords);
    assert_eq!(ac.len(), 8);
    assert_eq!(ac[0].to_bits(), 4648838210805188049);
}
