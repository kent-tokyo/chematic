//! Neighbor-enumeration integration tests, including comparison against an
//! independent brute-force oracle (hardcoded fixed-range triple loop per
//! pair, no shared code with `chematic_crystal::neighbor`'s search-box
//! logic).

use chematic_core::Element;
use chematic_crystal::{FractionalCoord, Lattice, PeriodicSite, PeriodicStructure, SiteSpecies};
use std::collections::HashSet;

fn brute_force_neighbors_within(
    structure: &PeriodicStructure,
    cutoff: f64,
    half: i32,
) -> HashSet<(usize, usize, [i32; 3])> {
    let lattice = structure.lattice();
    let sites = structure.sites();
    let n = sites.len();
    let mut found = HashSet::new();
    for i in 0..n {
        for j in 0..n {
            // Computed as (to - from) + n, matching the natural additive
            // decomposition the RFC's derivation itself uses (base + n) --
            // not "(to + n) - from", which for i == j (from == to exactly)
            // can lose the `x - x == 0.0` exactness IEEE754 guarantees and
            // land a boundary case (e.g. an image whose true distance is
            // exactly the cutoff) on the other side of `<=` purely from
            // summation-order rounding, unrelated to the search-box logic
            // this oracle exists to check.
            let raw = [
                sites[j].fractional.0[0] - sites[i].fractional.0[0],
                sites[j].fractional.0[1] - sites[i].fractional.0[1],
                sites[j].fractional.0[2] - sites[i].fractional.0[2],
            ];
            for a in -half..=half {
                for b in -half..=half {
                    for c in -half..=half {
                        if i == j && a == 0 && b == 0 && c == 0 {
                            continue;
                        }
                        let frac = [
                            raw[0] + f64::from(a),
                            raw[1] + f64::from(b),
                            raw[2] + f64::from(c),
                        ];
                        let cart = lattice.frac_to_cart(FractionalCoord::new(frac));
                        let dist =
                            (cart.0[0].powi(2) + cart.0[1].powi(2) + cart.0[2].powi(2)).sqrt();
                        if dist <= cutoff {
                            found.insert((i, j, [a, b, c]));
                        }
                    }
                }
            }
        }
    }
    found
}

fn assert_matches_oracle(structure: &PeriodicStructure, cutoff: f64, half: i32) {
    let got: HashSet<_> = structure
        .neighbors_within(cutoff)
        .unwrap()
        .into_iter()
        .map(|nb| (nb.center_index, nb.neighbor_index, nb.image))
        .collect();
    let oracle = brute_force_neighbors_within(structure, cutoff, half);
    assert_eq!(got, oracle, "cutoff={cutoff} half={half}");
}

fn two_site_cubic(a: f64) -> PeriodicStructure {
    PeriodicStructure::new(
        Lattice::cubic(a).unwrap(),
        vec![
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::NA)],
                FractionalCoord::new([0.0, 0.0, 0.0]),
                Some("A".to_string()),
            )
            .unwrap(),
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::CL)],
                FractionalCoord::new([0.5, 0.5, 0.5]),
                Some("B".to_string()),
            )
            .unwrap(),
        ],
    )
    .unwrap()
}

#[test]
fn cubic_matches_brute_force_oracle() {
    let s = two_site_cubic(4.0);
    for cutoff in [1.0, 3.5, 4.0, 6.0] {
        assert_matches_oracle(&s, cutoff, 4);
    }
}

#[test]
fn triclinic_matches_brute_force_oracle() {
    let lattice = Lattice::from_parameters(5.0, 6.0, 7.0, 75.0, 100.0, 60.0).unwrap();
    let s = PeriodicStructure::new(
        lattice,
        vec![
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::SI)],
                FractionalCoord::new([0.1, 0.2, 0.3]),
                None,
            )
            .unwrap(),
            PeriodicSite::new(
                vec![SiteSpecies::full(Element::O)],
                FractionalCoord::new([0.6, 0.4, 0.9]),
                None,
            )
            .unwrap(),
        ],
    )
    .unwrap();
    for cutoff in [2.0, 4.0, 6.0] {
        assert_matches_oracle(&s, cutoff, 4);
    }
}

#[test]
fn triclinic_sweep_matches_brute_force_oracle() {
    // Exercise different skew directions and fractional site distributions.
    // The oracle intentionally remains the fixed-range fractional triple
    // loop above, so this is independent of the cell-list bin construction.
    let cases = [
        (
            (4.5, 5.5, 6.5, 68.0, 82.0, 57.0),
            [[0.03, 0.17, 0.29], [0.41, 0.63, 0.88], [0.77, 0.22, 0.54]],
        ),
        (
            (7.0, 4.8, 5.2, 103.0, 71.0, 119.0),
            [[0.11, 0.91, 0.36], [0.58, 0.08, 0.79], [0.86, 0.47, 0.14]],
        ),
        (
            (6.2, 6.7, 4.4, 89.0, 96.0, 74.0),
            [[0.24, 0.05, 0.95], [0.69, 0.34, 0.12], [0.92, 0.81, 0.56]],
        ),
    ];
    for ((a, b, c, alpha, beta, gamma), coordinates) in cases {
        let lattice = Lattice::from_parameters(a, b, c, alpha, beta, gamma).unwrap();
        let sites = coordinates
            .into_iter()
            .enumerate()
            .map(|(index, fractional)| {
                PeriodicSite::new(
                    vec![SiteSpecies::full(if index == 1 {
                        Element::O
                    } else {
                        Element::SI
                    })],
                    FractionalCoord::new(fractional),
                    None,
                )
                .unwrap()
            })
            .collect();
        let structure = PeriodicStructure::new(lattice, sites).unwrap();
        for cutoff in [1.75, 3.25, 5.5] {
            assert_matches_oracle(&structure, cutoff, 5);
        }
    }
}

#[test]
fn seeded_triclinic_sweep_matches_brute_force_oracle() {
    // This is deliberately a small deterministic pseudo-random sweep rather
    // than a property-test dependency.  The seed and the generated cases are
    // part of the local proof: changing the implementation cannot silently
    // change the fixture population or make the oracle share cell-list logic.
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let next = |state: &mut u64| {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((*state >> 11) as f64) / ((1_u64 << 53) as f64)
    };

    for case_index in 0..12 {
        let a = 4.5 + 3.0 * next(&mut state);
        let b = 4.5 + 3.0 * next(&mut state);
        let c = 4.5 + 3.0 * next(&mut state);
        let alpha = 68.0 + 24.0 * next(&mut state);
        let beta = 68.0 + 24.0 * next(&mut state);
        let gamma = 68.0 + 24.0 * next(&mut state);
        let lattice =
            Lattice::from_parameters(a, b, c, alpha, beta, gamma).unwrap_or_else(|error| {
                panic!("seeded case {case_index} has invalid lattice: {error}")
            });
        let sites = (0..4)
            .map(|site_index| {
                PeriodicSite::new(
                    vec![SiteSpecies::full(if site_index % 2 == 0 {
                        Element::SI
                    } else {
                        Element::O
                    })],
                    FractionalCoord::new([next(&mut state), next(&mut state), next(&mut state)]),
                    Some(format!("seed-{case_index}-{site_index}")),
                )
                .unwrap()
            })
            .collect();
        let structure = PeriodicStructure::new(lattice, sites).unwrap();
        for cutoff in [1.5, 2.75, 4.25, 5.75] {
            assert_matches_oracle(&structure, cutoff, 6);
        }
    }
}

#[test]
fn site_order_change_does_not_change_the_physical_neighbor_set() {
    // Permute the two sites' insertion order; the *set* of physically
    // distinct (unordered pair + image, up to relabeling) neighbor facts
    // must be identical, only the index labels swap.
    let lattice = Lattice::cubic(4.0).unwrap();
    let site_a = PeriodicSite::new(
        vec![SiteSpecies::full(Element::NA)],
        FractionalCoord::new([0.0, 0.0, 0.0]),
        None,
    )
    .unwrap();
    let site_b = PeriodicSite::new(
        vec![SiteSpecies::full(Element::CL)],
        FractionalCoord::new([0.5, 0.5, 0.5]),
        None,
    )
    .unwrap();

    let s1 = PeriodicStructure::new(lattice.clone(), vec![site_a.clone(), site_b.clone()]).unwrap();
    let s2 = PeriodicStructure::new(lattice, vec![site_b, site_a]).unwrap();

    let n1 = s1.neighbors_within(6.0).unwrap();
    let n2 = s2.neighbors_within(6.0).unwrap();
    assert_eq!(n1.len(), n2.len());

    // Relabel s2's indices back to s1's labeling (0<->1 swap) and compare
    // as sets.
    let relabeled: HashSet<_> = n2
        .iter()
        .map(|nb| (1 - nb.center_index, 1 - nb.neighbor_index, nb.image))
        .collect();
    let original: HashSet<_> = n1
        .iter()
        .map(|nb| (nb.center_index, nb.neighbor_index, nb.image))
        .collect();
    assert_eq!(relabeled, original);
}

#[test]
fn no_duplicates_across_a_wider_cutoff_sweep() {
    let s = two_site_cubic(4.0);
    for cutoff in [1.0, 2.0, 3.0, 5.0, 8.0] {
        let neighbors = s.neighbors_within(cutoff).unwrap();
        let mut seen = HashSet::new();
        for nb in &neighbors {
            assert!(seen.insert((nb.center_index, nb.neighbor_index, nb.image)));
        }
    }
}
