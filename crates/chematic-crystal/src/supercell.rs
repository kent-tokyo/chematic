//! Diagonal supercell generation.
//!
//! `v0.1` supports only diagonal `[nx, ny, nz]` expansion, each `>= 1`.
//! Arbitrary 3x3 integer supercell transforms are out of scope (see
//! `docs/crystal_scope.md`).

use crate::error::CrystalError;
use crate::lattice::{Lattice, scale3};
use crate::site::{FractionalCoord, PeriodicSite};
use crate::structure::PeriodicStructure;

/// Build a diagonal supercell: lattice vectors scaled by `[nx, ny, nz]`,
/// sites replicated across every translated cell and re-normalized into
/// the new (larger) cell's fractional frame.
///
/// Site order is deterministic: outer loop over `i in 0..nx`, then `j in
/// 0..ny`, then `k in 0..nz`, then the original sites in their existing
/// order -- so site `((i*ny + j)*nz + k) * original_site_count +
/// original_index` in the result corresponds to translated-cell `(i, j,
/// k)`'s copy of `original_index`. Species, occupancy, and labels are
/// copied unchanged. `structure` itself is not modified.
///
/// # Errors
///
/// [`CrystalError::NonPositiveSupercellMultiplier`] if any of `nx, ny, nz`
/// is `0`.
pub fn make_supercell(
    structure: &PeriodicStructure,
    mult: [u32; 3],
) -> Result<PeriodicStructure, CrystalError> {
    for (axis, &m) in mult.iter().enumerate() {
        if m < 1 {
            return Err(CrystalError::NonPositiveSupercellMultiplier { axis, value: m });
        }
    }
    let [nx, ny, nz] = mult;

    let old_matrix = structure.lattice().matrix();
    let new_matrix = [
        scale3(old_matrix[0], f64::from(nx)),
        scale3(old_matrix[1], f64::from(ny)),
        scale3(old_matrix[2], f64::from(nz)),
    ];
    // Scaling each row independently by a positive factor can't introduce
    // a singularity or change the (scale-invariant) condition indicator --
    // see supercell_never_fails_lattice_validation_when_original_is_valid
    // below -- but the Result is still propagated rather than unwrapped,
    // matching this crate's pattern of never asserting away a fallible
    // constructor's own check.
    let new_lattice = Lattice::from_matrix(new_matrix)?;

    let orig_sites = structure.sites();
    let mut new_sites = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                for site in orig_sites {
                    let f = site.fractional.0;
                    let new_frac = FractionalCoord::new([
                        (f[0] + f64::from(i)) / f64::from(nx),
                        (f[1] + f64::from(j)) / f64::from(ny),
                        (f[2] + f64::from(k)) / f64::from(nz),
                    ]);
                    new_sites.push(PeriodicSite {
                        species: site.species.clone(),
                        fractional: new_frac,
                        label: site.label.clone(),
                    });
                }
            }
        }
    }

    PeriodicStructure::new(new_lattice, new_sites)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::periodic::minimum_image;
    use crate::site::SiteSpecies;
    use chematic_core::Element;

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
    fn identity_supercell_is_equivalent() {
        let s = two_site_cubic(4.0);
        let sc = make_supercell(&s, [1, 1, 1]).unwrap();
        assert_eq!(sc.site_count(), s.site_count());
        assert_eq!(sc.lattice().matrix(), s.lattice().matrix());
        for (a, b) in sc.sites().iter().zip(s.sites()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn supercell_2_1_1_doubles_a_axis_and_site_count() {
        let s = two_site_cubic(4.0);
        let sc = make_supercell(&s, [2, 1, 1]).unwrap();
        assert_eq!(sc.site_count(), 4);
        let m = sc.lattice().matrix();
        assert!((m[0][0] - 8.0).abs() < 1e-9);
        assert!((m[1][1] - 4.0).abs() < 1e-9);
        assert!((m[2][2] - 4.0).abs() < 1e-9);
    }

    #[test]
    fn supercell_2_2_2_site_count_is_original_times_multiplier_product() {
        let s = two_site_cubic(4.0);
        let sc = make_supercell(&s, [2, 2, 2]).unwrap();
        assert_eq!(sc.site_count(), s.site_count() * 8);
    }

    #[test]
    fn supercell_volume_scales_by_multiplier_product() {
        let s = two_site_cubic(4.0);
        for mult in [[1, 1, 1], [2, 1, 1], [2, 2, 2], [3, 1, 2]] {
            let sc = make_supercell(&s, mult).unwrap();
            let expected = s.lattice().volume() * f64::from(mult[0] * mult[1] * mult[2]);
            assert!(
                (sc.lattice().volume() - expected).abs() < 1e-6,
                "mult={mult:?} got={} expected={expected}",
                sc.lattice().volume()
            );
        }
    }

    #[test]
    fn supercell_preserves_species_and_occupancy() {
        let s = two_site_cubic(4.0);
        let sc = make_supercell(&s, [2, 2, 2]).unwrap();
        let orig_species: Vec<_> = s.sites().iter().map(|site| &site.species).collect();
        for site in sc.sites() {
            assert!(orig_species.contains(&&site.species));
        }
    }

    #[test]
    fn supercell_site_order_is_deterministic_across_calls() {
        let s = two_site_cubic(4.0);
        let a = make_supercell(&s, [2, 1, 2]).unwrap();
        let b = make_supercell(&s, [2, 1, 2]).unwrap();
        assert_eq!(a.sites(), b.sites());
    }

    #[test]
    fn supercell_preserves_nearest_neighbor_geometry() {
        // The nearest-neighbor distance from a site to its periodic images
        // must be identical before and after supercell expansion (the
        // physical structure hasn't changed, only how many periodic copies
        // are stored explicitly).
        let s = two_site_cubic(4.0);
        let sc = make_supercell(&s, [2, 2, 2]).unwrap();
        let orig_dist = minimum_image(
            s.lattice(),
            s.sites()[0].fractional,
            s.sites()[1].fractional,
        )
        .distance;
        let sc_dist = minimum_image(
            sc.lattice(),
            sc.sites()[0].fractional,
            sc.sites()[1].fractional,
        )
        .distance;
        assert!(
            (orig_dist - sc_dist).abs() < 1e-9,
            "orig={orig_dist} supercell={sc_dist}"
        );
    }

    #[test]
    fn original_structure_is_not_modified() {
        let s = two_site_cubic(4.0);
        let before = s.clone();
        let _ = make_supercell(&s, [3, 2, 1]).unwrap();
        assert_eq!(s, before);
    }

    #[test]
    fn rejects_zero_multiplier() {
        let s = two_site_cubic(4.0);
        let err = make_supercell(&s, [0, 1, 1]).unwrap_err();
        assert!(matches!(
            err,
            CrystalError::NonPositiveSupercellMultiplier { axis: 0, value: 0 }
        ));
    }

    #[test]
    fn supercell_never_fails_lattice_validation_when_original_is_valid() {
        // Anisotropic diagonal scaling preserves the condition indicator
        // (volume and bounding-box both scale by the same nx*ny*nz
        // factor), so a supercell of any validated lattice must itself
        // validate, for any positive integer multiplier triple.
        let l = Lattice::from_parameters(5.0, 6.0, 7.0, 75.0, 100.0, 60.0).unwrap();
        let s = PeriodicStructure::new(
            l,
            vec![
                PeriodicSite::new(
                    vec![SiteSpecies::full(Element::C)],
                    FractionalCoord::new([0.1, 0.2, 0.3]),
                    None,
                )
                .unwrap(),
            ],
        )
        .unwrap();
        for mult in [[1, 1, 1], [5, 1, 1], [1, 5, 1], [3, 4, 2]] {
            assert!(make_supercell(&s, mult).is_ok(), "mult={mult:?}");
        }
    }
}
