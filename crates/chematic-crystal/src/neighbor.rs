//! Cutoff-radius periodic neighbor enumeration.
//!
//! Reuses `crate::periodic::axis_bound`/`crate::periodic::padded_axis_range`
//! -- the same reciprocal-vector-derived search-box machinery
//! [`crate::periodic::minimum_image`] uses, with the cutoff radius itself as
//! the bound (no naive-round bootstrap candidate needed here, since the
//! bound is already given rather than something to be discovered).

use crate::error::CrystalError;
use crate::lattice::norm3;
use crate::periodic::{axis_bound, padded_axis_range};
use crate::site::FractionalCoord;
use crate::structure::PeriodicStructure;

/// Safety cap on the number of candidate periodic images examined per site
/// pair. The search-box width is identical for every pair in a structure
/// (it depends only on the cutoff and the lattice's reciprocal vectors, not
/// on any pair's own offset), so this is checked once per
/// [`neighbors_within`] call rather than per pair. Exists to turn "cutoff
/// accidentally far larger than the cell" (a likely unit-conversion
/// mistake -- e.g. an Angstrom cutoff passed as picometers) into a prompt
/// error instead of a search that silently examines billions of image
/// candidates.
pub const MAX_NEIGHBOR_IMAGE_CANDIDATES: u64 = 5_000_000;

/// One periodic neighbor relationship: `neighbor_index`'s image at `image`
/// is within cutoff of `center_index`'s zero image.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeriodicNeighbor {
    pub center_index: usize,
    pub neighbor_index: usize,
    /// Which periodic image of the neighbor site this entry describes.
    pub image: [i32; 3],
    /// Cartesian displacement from the center site to this neighbor image.
    pub displacement: [f64; 3],
    /// Euclidean length of `displacement`, in Angstrom.
    pub distance: f64,
}

/// Enumerate every periodic neighbor pair within `cutoff` (Angstrom,
/// inclusive: `distance <= cutoff`).
///
/// This is a **full** neighbor list: for two distinct sites `i != j`, both
/// `(center=i, neighbor=j, image=n)` and `(center=j, neighbor=i,
/// image=-n)` can appear as independent entries (each is a different
/// center's neighbor shell, not a duplicate of the other). `(same site,
/// zero image)` is always excluded; `(same site, nonzero image)` is kept
/// when within cutoff (physically valid for cells smaller than the
/// cutoff). Output is sorted by `(center_index, neighbor_index, image[0],
/// image[1], image[2])` for determinism.
///
/// # Errors
///
/// [`CrystalError::InvalidCutoff`] if `cutoff` is not finite and positive;
/// [`CrystalError::NeighborSearchTooLarge`] if the implied per-pair search
/// box exceeds [`MAX_NEIGHBOR_IMAGE_CANDIDATES`].
pub fn neighbors_within(
    structure: &PeriodicStructure,
    cutoff: f64,
) -> Result<Vec<PeriodicNeighbor>, CrystalError> {
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(CrystalError::InvalidCutoff { value: cutoff });
    }

    let lattice = structure.lattice();
    let recip_norms = lattice.reciprocal_row_norms();

    // Per-axis search width is the same for every pair (only the box's
    // center shifts per pair, via `raw`, not its size) -- so the guard
    // only needs computing once. Stay in f64 until the very end so a
    // pathological cutoff can't overflow an integer type first.
    let mut axis_width_f = [0.0f64; 3];
    for (axis, width) in axis_width_f.iter_mut().enumerate() {
        // 2*half is axis_bound's raw span; +3 matches padded_axis_range's
        // +-1 padding on each side plus one for inclusive-ceiling slack.
        *width = 2.0 * cutoff * recip_norms[axis] + 3.0;
    }
    let candidate_count_f = axis_width_f[0] * axis_width_f[1] * axis_width_f[2];
    if !candidate_count_f.is_finite() || candidate_count_f > MAX_NEIGHBOR_IMAGE_CANDIDATES as f64 {
        return Err(CrystalError::NeighborSearchTooLarge {
            candidate_count: if candidate_count_f.is_finite() {
                candidate_count_f as u64
            } else {
                u64::MAX
            },
            limit: MAX_NEIGHBOR_IMAGE_CANDIDATES,
        });
    }

    let sites = structure.sites();
    let n = sites.len();
    let mut results = Vec::new();

    for i in 0..n {
        for j in 0..n {
            let raw = [
                sites[j].fractional.0[0] - sites[i].fractional.0[0],
                sites[j].fractional.0[1] - sites[i].fractional.0[1],
                sites[j].fractional.0[2] - sites[i].fractional.0[2],
            ];
            let mut ranges = [(0i32, 0i32); 3];
            for axis in 0..3 {
                let (lo, hi) = axis_bound(raw[axis], cutoff, recip_norms[axis]);
                ranges[axis] = padded_axis_range(lo, hi);
            }

            for n0 in ranges[0].0..=ranges[0].1 {
                for n1 in ranges[1].0..=ranges[1].1 {
                    for n2 in ranges[2].0..=ranges[2].1 {
                        if i == j && n0 == 0 && n1 == 0 && n2 == 0 {
                            continue; // exclude only the trivial self-pair
                        }
                        let frac = [
                            raw[0] + f64::from(n0),
                            raw[1] + f64::from(n1),
                            raw[2] + f64::from(n2),
                        ];
                        let cart = lattice.frac_to_cart(FractionalCoord::new(frac));
                        let dist = norm3(cart.0);
                        if dist <= cutoff {
                            results.push(PeriodicNeighbor {
                                center_index: i,
                                neighbor_index: j,
                                image: [n0, n1, n2],
                                displacement: cart.0,
                                distance: dist,
                            });
                        }
                    }
                }
            }
        }
    }

    results.sort_by_key(|nb| {
        (
            nb.center_index,
            nb.neighbor_index,
            nb.image[0],
            nb.image[1],
            nb.image[2],
        )
    });

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::Lattice;
    use crate::site::{PeriodicSite, SiteSpecies};
    use chematic_core::Element;
    use std::collections::HashSet;

    fn single_site_cubic(a: f64) -> PeriodicStructure {
        PeriodicStructure::new(
            Lattice::cubic(a).unwrap(),
            vec![
                PeriodicSite::new(
                    vec![SiteSpecies::full(Element::AR)],
                    FractionalCoord::new([0.0, 0.0, 0.0]),
                    None,
                )
                .unwrap(),
            ],
        )
        .unwrap()
    }

    #[test]
    fn simple_cubic_one_site_finds_six_face_neighbors() {
        let s = single_site_cubic(3.0);
        let neighbors = neighbors_within(&s, 3.01).unwrap();
        // 6 face-adjacent periodic self-images at distance exactly 3.0.
        assert_eq!(neighbors.len(), 6);
        for nb in &neighbors {
            assert!((nb.distance - 3.0).abs() < 1e-9);
            assert_eq!(nb.center_index, 0);
            assert_eq!(nb.neighbor_index, 0);
            assert_ne!(nb.image, [0, 0, 0]);
        }
    }

    #[test]
    fn periodic_self_image_neighbor_not_excluded() {
        let s = single_site_cubic(2.0);
        let neighbors = neighbors_within(&s, 2.5).unwrap();
        assert!(
            neighbors
                .iter()
                .all(|nb| nb.center_index == 0 && nb.neighbor_index == 0)
        );
        assert!(!neighbors.is_empty());
    }

    #[test]
    fn zero_image_self_pair_always_excluded() {
        let s = single_site_cubic(2.0);
        // Even with a huge cutoff, (0,0,[0,0,0]) must never appear.
        let neighbors = neighbors_within(&s, 20.0).unwrap();
        assert!(
            !neighbors
                .iter()
                .any(|nb| nb.center_index == 0 && nb.neighbor_index == 0 && nb.image == [0, 0, 0])
        );
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
    fn two_site_structure_finds_cross_neighbors() {
        let s = two_site_cubic(4.0);
        // Site 0 <-> site 1 body-diagonal distance = sqrt(3)*2 ~ 3.464.
        let neighbors = neighbors_within(&s, 3.5).unwrap();
        let cross: Vec<_> = neighbors
            .iter()
            .filter(|nb| nb.center_index == 0 && nb.neighbor_index == 1)
            .collect();
        assert!(!cross.is_empty());
        for nb in &cross {
            assert!((nb.distance - 3.4641016151).abs() < 1e-6);
        }
    }

    #[test]
    fn cutoff_boundary_is_inclusive() {
        let s = single_site_cubic(3.0);
        let exact = neighbors_within(&s, 3.0).unwrap();
        assert_eq!(exact.len(), 6, "distance==cutoff must be included");

        let below = neighbors_within(&s, 2.999).unwrap();
        assert_eq!(below.len(), 0, "distance>cutoff must be excluded");

        let above = neighbors_within(&s, 3.001).unwrap();
        assert_eq!(above.len(), 6);
    }

    #[test]
    fn no_duplicate_entries() {
        let s = two_site_cubic(4.0);
        let neighbors = neighbors_within(&s, 6.0).unwrap();
        let mut seen = HashSet::new();
        for nb in &neighbors {
            let key = (nb.center_index, nb.neighbor_index, nb.image);
            assert!(seen.insert(key), "duplicate entry: {key:?}");
        }
    }

    #[test]
    fn deterministic_ordering() {
        let s = two_site_cubic(4.0);
        let a = neighbors_within(&s, 6.0).unwrap();
        let b = neighbors_within(&s, 6.0).unwrap();
        assert_eq!(a, b);
        // sorted by (center, neighbor, image) ascending
        for w in a.windows(2) {
            let ka = (
                w[0].center_index,
                w[0].neighbor_index,
                w[0].image[0],
                w[0].image[1],
                w[0].image[2],
            );
            let kb = (
                w[1].center_index,
                w[1].neighbor_index,
                w[1].image[0],
                w[1].image[1],
                w[1].image[2],
            );
            assert!(ka <= kb);
        }
    }

    #[test]
    fn triclinic_structure_neighbors_are_finite_and_within_cutoff() {
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
        let neighbors = neighbors_within(&s, 4.0).unwrap();
        assert!(!neighbors.is_empty());
        for nb in &neighbors {
            assert!(nb.distance <= 4.0 + 1e-9);
            assert!(nb.distance.is_finite());
        }
    }

    #[test]
    fn invalid_cutoff_rejected() {
        let s = single_site_cubic(3.0);
        assert!(matches!(
            neighbors_within(&s, 0.0),
            Err(CrystalError::InvalidCutoff { .. })
        ));
        assert!(matches!(
            neighbors_within(&s, -1.0),
            Err(CrystalError::InvalidCutoff { .. })
        ));
        assert!(matches!(
            neighbors_within(&s, f64::NAN),
            Err(CrystalError::InvalidCutoff { .. })
        ));
        assert!(matches!(
            neighbors_within(&s, f64::INFINITY),
            Err(CrystalError::InvalidCutoff { .. })
        ));
    }

    #[test]
    fn absurd_cutoff_rejected_as_too_large() {
        let s = single_site_cubic(1.0);
        // A cutoff ~1e8x the cell edge implies an astronomically large
        // search box.
        let err = neighbors_within(&s, 1e8).unwrap_err();
        assert!(matches!(err, CrystalError::NeighborSearchTooLarge { .. }));
    }
}
