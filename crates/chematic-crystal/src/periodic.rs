//! Exact periodic-boundary displacement and minimum-image distance.
//!
//! # Why not `delta_frac -= delta_frac.round()`
//!
//! That one-line reduction is only guaranteed to find the *true* nearest
//! periodic image when the cell is close to orthogonal. For a sufficiently
//! skewed triclinic cell, "nearest per fractional axis" and "nearest in
//! Cartesian space" diverge -- rounding each fractional component
//! independently can miss an image that is farther in fractional space but
//! closer in Cartesian space. This module instead derives a finite,
//! provably sufficient search box from the lattice's own reciprocal
//! vectors (see `axis_bound` below) and brute-force-checks every candidate
//! inside it, which is exact for any lattice [`Lattice::from_matrix`]
//! accepts (validated non-singular). See
//! `docs/rfcs/chematic_crystal_foundation.md` for the full derivation and
//! `tests/periodicity.rs` for verification against an independent
//! brute-force oracle, including a pinned triclinic regression fixture
//! where naive `round()` disagrees with the exact result.

use crate::lattice::{Lattice, norm3};
use crate::site::FractionalCoord;

/// A minimum-image periodic displacement between two fractional positions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeriodicDisplacement {
    /// Cartesian displacement vector (Angstrom), `to`'s chosen image minus
    /// `from`.
    pub cartesian: [f64; 3],
    /// The same displacement in fractional coordinates.
    pub fractional: [f64; 3],
    /// Which periodic image of `to` was chosen (integer lattice-vector
    /// shift added to `to`'s fractional position).
    pub image: [i32; 3],
    /// Euclidean length of `cartesian`, in Angstrom.
    pub distance: f64,
}

/// Real-valued (not yet floored/ceiled/padded) bound on the integer
/// image-shift component along one axis: any image achieving Cartesian
/// distance `<= bound` from a point at fractional offset `base_component`
/// (on this axis) must have its integer shift `m` on this axis inside
/// `[lo, hi]`.
///
/// Derivation (full version in the RFC): for `r = (base + m) . M` and
/// reciprocal row `b_j` (`a_i . b_j = delta_ij`), `(base + m)_j = r . b_j`,
/// so by Cauchy-Schwarz `|base_j + m_j| <= |r| * |b_j| <= bound * |b_j|`
/// whenever `|r| <= bound`.
///
/// Shared by [`minimum_image`] (`bound` = the naive round()-based
/// candidate's own distance) and [`crate::neighbor`] (`bound` = the
/// cutoff radius) -- same bounding-box shape, different `bound` value and
/// different post-processing (this module keeps only the minimum;
/// `neighbor` keeps every candidate inside the box that also passes the
/// cutoff after an exact distance check).
pub(crate) fn axis_bound(base_component: f64, bound: f64, reciprocal_row_norm: f64) -> (f64, f64) {
    let half = bound * reciprocal_row_norm;
    (-half - base_component, half - base_component)
}

/// Convert a real-valued `axis_bound` result into an inclusive `i32` range,
/// padded by one extra integer on each side as a defensive margin against
/// floating-point error landing a true boundary case just outside the
/// unpadded `floor`/`ceil` (cheap: the search box stays small either way
/// for any lattice that passed [`Lattice`]'s own validation).
pub(crate) fn padded_axis_range(lo: f64, hi: f64) -> (i32, i32) {
    (lo.floor() as i32 - 1, hi.ceil() as i32 + 1)
}

/// The exact minimum-image displacement from fractional position `from` to
/// fractional position `to` under `lattice`'s periodic boundary conditions.
///
/// Neither `from` nor `to` need to be pre-wrapped into `[0, 1)` -- only
/// their difference matters. Ties (multiple images at the same minimal
/// distance, e.g. a perfectly cubic cell's face-centered midpoint) are
/// broken deterministically by iteration order: the first candidate
/// encountered while enumerating image-shift components in ascending order
/// is kept. Since `image = image0 + m` and `m` is enumerated ascending per
/// axis (`image0` a fixed per-call offset), this is equivalent to keeping
/// the lexicographically smallest tied `image`.
///
/// # Examples
///
/// ```
/// use chematic_crystal::{FractionalCoord, Lattice, minimum_image};
///
/// let lattice = Lattice::cubic(10.0)?;
/// let from = FractionalCoord::new([0.05, 0.5, 0.5]);
/// let to = FractionalCoord::new([0.95, 0.5, 0.5]);
/// // Direct fractional difference is 0.9 (9.0 Angstrom); the periodic
/// // image wraps around to 0.1 (1.0 Angstrom).
/// let displacement = minimum_image(&lattice, from, to);
/// assert!((displacement.distance - 1.0).abs() < 1e-9);
/// assert_eq!(displacement.image, [-1, 0, 0]);
/// # Ok::<(), chematic_crystal::CrystalError>(())
/// ```
pub fn minimum_image(
    lattice: &Lattice,
    from: FractionalCoord,
    to: FractionalCoord,
) -> PeriodicDisplacement {
    let raw = [
        to.0[0] - from.0[0],
        to.0[1] - from.0[1],
        to.0[2] - from.0[2],
    ];
    // n0 is the naive round()-based image shift; base = raw + n0 is that
    // candidate's fractional displacement, already reduced into [-0.5, 0.5]
    // per axis.
    let n0 = [raw[0].round(), raw[1].round(), raw[2].round()];
    let base = [raw[0] - n0[0], raw[1] - n0[1], raw[2] - n0[2]];
    let image0 = [-(n0[0] as i32), -(n0[1] as i32), -(n0[2] as i32)];

    let cart0 = lattice.frac_to_cart(FractionalCoord::new(base));
    let dist0 = norm3(cart0.0);

    let recip_norms = lattice.reciprocal_row_norms();
    let mut ranges = [(0i32, 0i32); 3];
    for axis in 0..3 {
        let (lo, hi) = axis_bound(base[axis], dist0, recip_norms[axis]);
        ranges[axis] = padded_axis_range(lo, hi);
    }

    // Not seeded with the naive candidate: ties must be broken by ascending
    // search order (first-found wins), not by which candidate happened to
    // be computed first. The naive candidate (m = [0, 0, 0]) is itself
    // always inside `ranges` (its own distance is the `bound` the ranges
    // were derived from), so it's still considered -- just not favored.
    let mut best: Option<PeriodicDisplacement> = None;

    for m0 in ranges[0].0..=ranges[0].1 {
        for m1 in ranges[1].0..=ranges[1].1 {
            for m2 in ranges[2].0..=ranges[2].1 {
                let frac = [
                    base[0] + f64::from(m0),
                    base[1] + f64::from(m1),
                    base[2] + f64::from(m2),
                ];
                let cart = lattice.frac_to_cart(FractionalCoord::new(frac));
                let dist = norm3(cart.0);
                if best.as_ref().is_none_or(|current| dist < current.distance) {
                    best = Some(PeriodicDisplacement {
                        cartesian: cart.0,
                        fractional: frac,
                        image: [image0[0] + m0, image0[1] + m1, image0[2] + m2],
                        distance: dist,
                    });
                }
            }
        }
    }

    best.expect("search ranges always contain at least one candidate (m = [0, 0, 0])")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::Lattice;

    #[test]
    fn cubic_zero_displacement_is_zero_distance_zero_image() {
        let l = Lattice::cubic(5.0).unwrap();
        let f = FractionalCoord::new([0.3, 0.3, 0.3]);
        let d = minimum_image(&l, f, f);
        assert_eq!(d.image, [0, 0, 0]);
        assert!(d.distance < 1e-12);
    }

    #[test]
    fn cubic_wraps_to_nearest_image() {
        let l = Lattice::cubic(10.0).unwrap();
        let from = FractionalCoord::new([0.05, 0.5, 0.5]);
        let to = FractionalCoord::new([0.95, 0.5, 0.5]);
        // Direct distance would be 0.9*10=9; the periodic image is 0.1*10=1.
        let d = minimum_image(&l, from, to);
        assert!((d.distance - 1.0).abs() < 1e-9, "distance={}", d.distance);
        assert_eq!(d.image[0], -1);
    }

    #[test]
    fn translated_fractional_produces_same_minimum_image_distance() {
        // Integer translation invariance: shifting `to` by a whole lattice
        // vector must not change the minimum-image distance.
        let l = Lattice::from_parameters(6.0, 7.0, 8.0, 85.0, 95.0, 100.0).unwrap();
        let from = FractionalCoord::new([0.1, 0.2, 0.3]);
        let to = FractionalCoord::new([0.6, 0.4, 0.9]);
        let base_dist = minimum_image(&l, from, to).distance;
        for shift in [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-2, 3, -1], [5, -5, 5]] {
            let to_shifted = to.translated(shift);
            let shifted_dist = minimum_image(&l, from, to_shifted).distance;
            assert!(
                (base_dist - shifted_dist).abs() < 1e-9,
                "shift {shift:?}: {base_dist} vs {shifted_dist}"
            );
        }
    }

    #[test]
    fn origin_translation_of_both_points_does_not_change_distance() {
        let l = Lattice::from_parameters(6.0, 7.0, 8.0, 85.0, 95.0, 100.0).unwrap();
        let from = FractionalCoord::new([0.1, 0.2, 0.3]);
        let to = FractionalCoord::new([0.6, 0.4, 0.9]);
        let base_dist = minimum_image(&l, from, to).distance;
        let shift = [0.37, -1.21, 2.05];
        let from2 = FractionalCoord::new([
            from.0[0] + shift[0],
            from.0[1] + shift[1],
            from.0[2] + shift[2],
        ]);
        let to2 =
            FractionalCoord::new([to.0[0] + shift[0], to.0[1] + shift[1], to.0[2] + shift[2]]);
        let shifted_dist = minimum_image(&l, from2, to2).distance;
        assert!((base_dist - shifted_dist).abs() < 1e-9);
    }

    /// Regression for a tie-break bug: `best` used to be seeded with the
    /// naive `round()`-based candidate before the search loop, so a
    /// later candidate at the *same* distance never replaced it even
    /// when it was lexicographically smaller. Here `from = [0.5, 0, 0]`,
    /// `to = [0, 0, 0]` in a cubic cell: `image = [0, 0, 0]` and
    /// `image = [1, 0, 0]` are equidistant (2.0 Angstrom), and
    /// `round(-0.5) == -1` makes `[1, 0, 0]` the naive candidate. The
    /// correct, ascending-search-order tie-break picks `[0, 0, 0]`.
    #[test]
    fn cubic_half_cell_tie_uses_lexicographically_smallest_image() {
        let lattice = Lattice::cubic(4.0).unwrap();
        let from = FractionalCoord::new([0.5, 0.0, 0.0]);
        let to = FractionalCoord::new([0.0, 0.0, 0.0]);

        let result = minimum_image(&lattice, from, to);

        assert!((result.distance - 2.0).abs() < 1e-12);
        assert_eq!(result.image, [0, 0, 0]);
    }

    /// Same tie, `from`/`to` swapped. The invariant is "lexicographically
    /// smallest image wins," not "mirror the forward case": swapping
    /// `from`/`to` shifts which two images are tied (`image` is defined
    /// relative to `to`), so here the tied pair is `[-1, 0, 0]` /
    /// `[0, 0, 0]` and the smaller one, `[-1, 0, 0]`, wins -- still
    /// confirming the distance is symmetric either way.
    #[test]
    fn cubic_half_cell_tie_reverse_direction_also_resolves_deterministically() {
        let lattice = Lattice::cubic(4.0).unwrap();
        let from = FractionalCoord::new([0.0, 0.0, 0.0]);
        let to = FractionalCoord::new([0.5, 0.0, 0.0]);

        let result = minimum_image(&lattice, from, to);

        assert!((result.distance - 2.0).abs() < 1e-12);
        assert_eq!(result.image, [-1, 0, 0]);
    }
}
