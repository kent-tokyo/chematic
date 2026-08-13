//! Periodicity integration tests: verify `chematic_crystal::minimum_image`
//! against an independent, from-scratch brute-force oracle (no shared code
//! with the production bounded-search implementation, other than the
//! `Lattice::frac_to_cart` primitive both necessarily use to turn a
//! fractional candidate into a Cartesian distance -- the property under
//! test is the *image search strategy*, not that primitive).

use chematic_crystal::{FractionalCoord, Lattice, minimum_image};

/// Independent brute-force minimum-image oracle: enumerate every integer
/// image shift in `-half..=half` per axis directly, keep the closest.
/// Deliberately hardcoded rather than calling any of `chematic_crystal`'s
/// own search-bound derivation.
fn brute_force_minimum_image(
    lattice: &Lattice,
    from: FractionalCoord,
    to: FractionalCoord,
    half: i32,
) -> (f64, [i32; 3]) {
    // (to - from) + n, not "(to + n) - from" -- see the comment on the
    // analogous oracle in tests/neighbor.rs for why the order matters at
    // exact-distance boundaries.
    let raw = [
        to.0[0] - from.0[0],
        to.0[1] - from.0[1],
        to.0[2] - from.0[2],
    ];
    let mut best_dist = f64::INFINITY;
    let mut best_image = [0i32; 3];
    for i in -half..=half {
        for j in -half..=half {
            for k in -half..=half {
                let frac = [
                    raw[0] + f64::from(i),
                    raw[1] + f64::from(j),
                    raw[2] + f64::from(k),
                ];
                let cart = lattice.frac_to_cart(FractionalCoord::new(frac));
                let dist = (cart.0[0].powi(2) + cart.0[1].powi(2) + cart.0[2].powi(2)).sqrt();
                if dist < best_dist {
                    best_dist = dist;
                    best_image = [i, j, k];
                }
            }
        }
    }
    (best_dist, best_image)
}

fn assert_matches_oracle(lattice: &Lattice, from: FractionalCoord, to: FractionalCoord, half: i32) {
    let (oracle_dist, oracle_image) = brute_force_minimum_image(lattice, from, to, half);
    let got = minimum_image(lattice, from, to);
    assert!(
        (got.distance - oracle_dist).abs() < 1e-7,
        "distance mismatch: impl={} oracle={} (lattice={:?} from={:?} to={:?})",
        got.distance,
        oracle_dist,
        lattice.matrix(),
        from.0,
        to.0
    );
    // Multiple images can tie at the same minimal distance (e.g. a face
    // midpoint in a cubic cell); only assert the image matches when the
    // oracle's own best isn't tied with a neighboring candidate.
    let raw = [
        to.0[0] - from.0[0],
        to.0[1] - from.0[1],
        to.0[2] - from.0[2],
    ];
    let tie_free = {
        let mut count = 0;
        for i in -half..=half {
            for j in -half..=half {
                for k in -half..=half {
                    let frac = [
                        raw[0] + f64::from(i),
                        raw[1] + f64::from(j),
                        raw[2] + f64::from(k),
                    ];
                    let cart = lattice.frac_to_cart(FractionalCoord::new(frac));
                    let dist = (cart.0[0].powi(2) + cart.0[1].powi(2) + cart.0[2].powi(2)).sqrt();
                    if (dist - oracle_dist).abs() < 1e-9 {
                        count += 1;
                    }
                }
            }
        }
        count == 1
    };
    if tie_free {
        assert_eq!(
            got.image,
            oracle_image,
            "image mismatch (untied): lattice={:?} from={:?} to={:?}",
            lattice.matrix(),
            from.0,
            to.0
        );
    }
}

// -- fixed-shape oracle comparisons --------------------------------------

#[test]
fn cubic_minimum_image_matches_brute_force() {
    let l = Lattice::cubic(5.0).unwrap();
    let pairs = [
        ([0.05, 0.5, 0.5], [0.95, 0.5, 0.5]),
        ([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]),
        ([0.1, 0.2, 0.3], [0.9, 0.8, 0.7]),
        ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
    ];
    for (from, to) in pairs {
        assert_matches_oracle(&l, FractionalCoord::new(from), FractionalCoord::new(to), 4);
    }
}

#[test]
fn orthorhombic_minimum_image_matches_brute_force() {
    let l = Lattice::orthorhombic(4.0, 9.0, 15.0).unwrap();
    let pairs = [
        ([0.02, 0.02, 0.02], [0.98, 0.98, 0.98]),
        ([0.3, 0.6, 0.1], [0.7, 0.1, 0.95]),
    ];
    for (from, to) in pairs {
        assert_matches_oracle(&l, FractionalCoord::new(from), FractionalCoord::new(to), 4);
    }
}

#[test]
fn triclinic_minimum_image_matches_brute_force() {
    let l = Lattice::from_parameters(5.0, 6.0, 7.0, 75.0, 100.0, 60.0).unwrap();
    let pairs = [
        ([0.02, 0.02, 0.02], [0.98, 0.98, 0.98]),
        ([0.1, 0.9, 0.4], [0.85, 0.05, 0.6]),
        ([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]),
    ];
    for (from, to) in pairs {
        assert_matches_oracle(&l, FractionalCoord::new(from), FractionalCoord::new(to), 4);
    }
}

// -- pinned regression fixture: naive round() disagrees with the exact
// answer, exact answer independently confirmed against the brute-force
// oracle above (half=8) --------------------------------------------------

#[test]
fn skewed_triclinic_naive_round_disagrees_with_exact_minimum_image() {
    // Found by randomized search (fixed seed, see
    // randomized_triclinic_matches_brute_force_oracle below for the
    // generator shape) and pinned here as a permanent regression fixture:
    // naive `delta_frac -= delta_frac.round()` reports distance ~8.92
    // Angstrom for this pair, but the true nearest periodic image (image
    // [1, 0, 1]) is ~1.73 Angstrom away -- confirmed independently against
    // a hardcoded brute-force oracle over -8..=8 per axis (see this
    // module's search history; reproduced directly below via
    // brute_force_minimum_image with half=8).
    let l = Lattice::from_parameters(
        7.987030736455244,
        7.905173387101117,
        6.17331345211384,
        40.620511882721246,
        63.965671339359986,
        29.95313760584825,
    )
    .unwrap();
    let from = FractionalCoord::new([0.79075334224289, 0.8185506937359817, 0.9768055390549911]);
    let to = FractionalCoord::new([0.2902896769377683, 0.3139302234302791, 0.2814707286846695]);

    // Naive round()-only "minimum image" (the rejected approach).
    let raw = [
        to.0[0] - from.0[0],
        to.0[1] - from.0[1],
        to.0[2] - from.0[2],
    ];
    let naive_frac = [
        raw[0] - raw[0].round(),
        raw[1] - raw[1].round(),
        raw[2] - raw[2].round(),
    ];
    let naive_cart = l.frac_to_cart(FractionalCoord::new(naive_frac));
    let naive_dist =
        (naive_cart.0[0].powi(2) + naive_cart.0[1].powi(2) + naive_cart.0[2].powi(2)).sqrt();
    assert!(
        (naive_dist - 8.923828466561345).abs() < 1e-6,
        "naive_dist={naive_dist}"
    );

    let got = minimum_image(&l, from, to);
    assert!(
        (got.distance - 1.7281501267279693).abs() < 1e-6,
        "distance={}",
        got.distance
    );
    assert_eq!(got.image, [1, 0, 1]);

    // The whole point: naive and exact must disagree by a wide margin, not
    // a rounding-noise-sized amount.
    assert!(naive_dist - got.distance > 5.0);

    // Independent oracle confirmation, same fixture.
    let (oracle_dist, oracle_image) = brute_force_minimum_image(&l, from, to, 8);
    assert!((oracle_dist - got.distance).abs() < 1e-7);
    assert_eq!(oracle_image, got.image);
}

#[test]
fn near_singular_triclinic_requires_large_image_shift() {
    // The randomized-search draw that motivated widening the property
    // test's oracle half-width (see the comment there): a cell with
    // condition_indicator ~0.02 (near, but above, the RFC's 1e-3 rejection
    // floor) whose true nearest periodic image has magnitude 11 on two
    // axes. Independently reconfirmed with a hardcoded brute force up to
    // half=30 (converges to the same answer from half=15 onward) --
    // pinned here as a permanent regression fixture for "the bounded
    // search box must actually grow to cover cells this skewed, not just
    // a few units either way".
    let matrix = [
        [4.059173556199132, 0.0, 0.0],
        [-3.500362048473072, 3.554383756790896, 0.0],
        [1.602131464217347, 7.317385649094703, 0.21424309247456755],
    ];
    let l = Lattice::from_matrix(matrix).unwrap();
    assert!(l.condition_indicator() < 0.03);
    let from = FractionalCoord::new([0.6289642968793616, 0.7796459005054238, 0.1371601305857918]);
    let to = FractionalCoord::new([0.7690665050911141, 0.2146789376436573, 0.09206061526606202]);

    let got = minimum_image(&l, from, to);
    assert!(
        (got.distance - 1.253297300461551).abs() < 1e-6,
        "distance={}",
        got.distance
    );
    assert_eq!(got.image, [11, 11, -5]);

    let (oracle_dist, oracle_image) = brute_force_minimum_image(&l, from, to, 15);
    assert!((oracle_dist - got.distance).abs() < 1e-7);
    assert_eq!(oracle_image, got.image);
}

// -- translation invariance / wrapping -----------------------------------

#[test]
fn fractional_integer_translation_invariance() {
    let l = Lattice::from_parameters(6.0, 7.0, 8.0, 70.0, 100.0, 110.0).unwrap();
    let from = FractionalCoord::new([0.15, 0.25, 0.35]);
    let to = FractionalCoord::new([0.65, 0.85, 0.05]);
    let base = minimum_image(&l, from, to).distance;
    for shift in [[3, -2, 1], [-1, -1, -1], [7, 0, -4]] {
        let shifted = minimum_image(&l, from, to.translated(shift)).distance;
        assert!((base - shifted).abs() < 1e-9, "shift={shift:?}");
    }
}

#[test]
fn wrapped_fractional_coordinates_land_in_unit_range() {
    let coords = [
        FractionalCoord::new([1.7, -0.3, 4.999999]),
        FractionalCoord::new([-3.2, 2.4, -0.0001]),
        FractionalCoord::new([0.0, 1.0, 0.5]),
    ];
    for c in coords {
        let w = c.wrapped();
        for component in w.0 {
            assert!(
                (0.0..1.0).contains(&component),
                "component {component} out of [0,1) for input {c:?}"
            );
        }
    }
}

#[test]
fn lattice_origin_translation_invariance() {
    let l = Lattice::from_parameters(6.0, 7.0, 8.0, 70.0, 100.0, 110.0).unwrap();
    let from = FractionalCoord::new([0.15, 0.25, 0.35]);
    let to = FractionalCoord::new([0.65, 0.85, 0.05]);
    let base = minimum_image(&l, from, to).distance;
    let origin_shift = [1.37, -2.81, 0.44];
    let from2 = FractionalCoord::new([
        from.0[0] + origin_shift[0],
        from.0[1] + origin_shift[1],
        from.0[2] + origin_shift[2],
    ]);
    let to2 = FractionalCoord::new([
        to.0[0] + origin_shift[0],
        to.0[1] + origin_shift[1],
        to.0[2] + origin_shift[2],
    ]);
    let shifted = minimum_image(&l, from2, to2).distance;
    assert!((base - shifted).abs() < 1e-9);
}

// -- property-based / randomized (fixed seed, no external RNG dependency) --

/// Minimal SplitMix64 PRNG -- avoids adding a `rand` dependency (not used
/// anywhere else in this workspace either) for a handful of fixed-seed
/// property tests. Not cryptographic; fine for test-fixture generation.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn f64_in(&mut self, lo: f64, hi: f64) -> f64 {
        let unit = (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
        lo + unit * (hi - lo)
    }
}

#[test]
fn randomized_triclinic_matches_brute_force_oracle() {
    // Fixed seed -> fully reproducible. Deliberately the same generator
    // shape used to originally discover the pinned regression fixture
    // above, so this test is the reason to trust that fixture's "exact"
    // value wasn't cherry-picked from a broken implementation: every one of
    // these 500 random triclinic cells is independently checked against
    // the brute-force oracle too.
    //
    // condition_indicator floor is 0.15, not the RFC's validation floor of
    // 1e-3: a *fixed*-half-width brute-force oracle (needed here so the
    // oracle stays independent of `chematic_crystal`'s own reciprocal-
    // vector-derived search-box logic) is only trustworthy when it's wide
    // enough to contain the true minimum. An earlier version of this test
    // used condition_indicator >= 0.02 and half=6, and hit a real draw
    // (condition_indicator ~0.0204) whose true nearest image was [11, 11,
    // -5] -- half=6 (and even half=10) missed it entirely, while
    // `chematic_crystal::minimum_image` found it correctly (independently
    // reconfirmed with a half=30 brute force). That was the oracle being
    // too narrow, not a production bug -- but it means a *fixed*-half
    // oracle can't safely cover cells that close to singular. The
    // near-threshold regime is covered separately by the pinned regression
    // fixture above, whose "exact" distance was independently reconfirmed
    // up to half=8 (that fixture's true image only had magnitude 1, so
    // half=8 was already generous there).
    let mut rng = Rng::new(0x00C0_FFEE_1234_5678);
    let mut checked = 0;
    for _ in 0..500 {
        let a = rng.f64_in(2.0, 8.0);
        let b = rng.f64_in(2.0, 8.0);
        let c = rng.f64_in(2.0, 8.0);
        let alpha = rng.f64_in(25.0, 155.0);
        let beta = rng.f64_in(25.0, 155.0);
        let gamma = rng.f64_in(25.0, 155.0);
        let Ok(l) = Lattice::from_parameters(a, b, c, alpha, beta, gamma) else {
            continue;
        };
        if l.condition_indicator() < 0.15 {
            continue; // skip near-degenerate draws -- a fixed-half brute-force oracle isn't trustworthy there, see comment above
        }
        let from = FractionalCoord::new([
            rng.f64_in(0.0, 1.0),
            rng.f64_in(0.0, 1.0),
            rng.f64_in(0.0, 1.0),
        ]);
        let to = FractionalCoord::new([
            rng.f64_in(0.0, 1.0),
            rng.f64_in(0.0, 1.0),
            rng.f64_in(0.0, 1.0),
        ]);
        assert_matches_oracle(&l, from, to, 8);
        checked += 1;
    }
    assert!(
        checked > 300,
        "too many draws skipped as near-degenerate: {checked}/500"
    );
}
