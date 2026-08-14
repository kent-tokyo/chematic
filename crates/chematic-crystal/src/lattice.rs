//! [`Lattice`]: a validated 3x3 lattice matrix + cached inverse.
//!
//! # Matrix convention
//!
//! `matrix()` rows are the lattice vectors **a, b, c** (`matrix[0]` = a,
//! `matrix[1]` = b, `matrix[2]` = c). A Cartesian point is obtained from a
//! fractional one by row-vector x matrix:
//!
//! ```text
//! cartesian = fractional . matrix
//! cartesian_k = sum_j fractional_j * matrix[j][k]
//! ```
//!
//! This matches `chematic_ewald::BoxVectors`'s row convention (see
//! `docs/rfcs/chematic_crystal_foundation.md`), without `chematic-crystal`
//! depending on `chematic-ewald`.
//!
//! Reciprocal vectors use the crystallographic (no `2*pi`) convention,
//! `a_i . b_j = delta_ij`; under the row convention, reciprocal row `b_j` is
//! **column `j`** of `inverse`, not row `j` -- asserted directly by this
//! module's tests rather than trusted by inspection.

use crate::error::CrystalError;
use crate::site::{CartesianCoord, FractionalCoord};
use crate::validation::{require_finite, require_finite3};

// ---------------------------------------------------------------------------
// Small vector helpers (crate-internal; periodic.rs/neighbor.rs reuse these)
// ---------------------------------------------------------------------------

pub(crate) fn dot3(u: [f64; 3], v: [f64; 3]) -> f64 {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}

pub(crate) fn cross3(u: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

pub(crate) fn norm3(u: [f64; 3]) -> f64 {
    dot3(u, u).sqrt()
}

pub(crate) fn scale3(u: [f64; 3], s: f64) -> [f64; 3] {
    [u[0] * s, u[1] * s, u[2] * s]
}

/// Signed volume `a . (b x c)` of the parallelepiped with rows `m[0..3]`.
/// Equal to `det(m)`.
fn signed_volume(m: &[[f64; 3]; 3]) -> f64 {
    dot3(m[0], cross3(m[1], m[2]))
}

// ---------------------------------------------------------------------------
// Lattice
// ---------------------------------------------------------------------------

/// A validated 3x3 lattice matrix (rows = lattice vectors a, b, c) with its
/// cached inverse.
///
/// All constructors reject `NaN`/`Infinity`, non-positive crystallographic
/// lengths, angles that don't define a non-degenerate parallelepiped, and
/// matrices that are exactly or numerically singular (see
/// [`Self::MIN_CONDITION_INDICATOR`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Lattice {
    matrix: [[f64; 3]; 3],
    inverse: [[f64; 3]; 3],
}

impl Lattice {
    /// Minimum accepted crystallographic length, in Angstrom.
    ///
    /// Below this, a length is essentially guaranteed to be a unit-
    /// conversion or parsing error rather than a physically meaningful
    /// lattice constant -- real crystal lattice parameters are on the order
    /// of a few Angstrom at the very smallest; `1e-6` is far below any
    /// physical crystal while staying safely above float noise for any
    /// input that was ever a real length.
    pub const MIN_LENGTH: f64 = 1e-6;

    /// Minimum accepted dimensionless condition indicator, `|det(M)| /
    /// (|a| |b| |c|)`.
    ///
    /// This ratio is the fraction of the bounding box volume that the
    /// actual parallelepiped occupies: `1.0` for a fully orthogonal cell,
    /// `-> 0` as the three lattice vectors approach coplanarity (a singular
    /// matrix). `1e-3` is a conservative floor well below any physically
    /// reasonable crystallographic cell (even a very acute/obtuse but real
    /// cell rarely drops under ~0.05) while still catching near-degenerate
    /// matrices before the inverse -- needed by every fractional<->Cartesian
    /// conversion and the minimum-image algorithm -- becomes numerically
    /// unstable.
    pub const MIN_CONDITION_INDICATOR: f64 = 1e-3;

    /// Tolerance (relative to `c^2`) for the `cz^2 >= 0` consistency check
    /// in [`Self::from_parameters`]: floating-point evaluation of the IUCr
    /// formula can land a hair below zero for angles that are, up to
    /// rounding, exactly degenerate. Genuine angle incompatibilities miss
    /// by far more than this.
    const ANGLE_CONSISTENCY_TOLERANCE: f64 = 1e-9;

    /// Build a lattice directly from a 3x3 matrix (rows = a, b, c).
    pub fn from_matrix(matrix: [[f64; 3]; 3]) -> Result<Self, CrystalError> {
        require_finite3(matrix[0], "matrix[0]")?;
        require_finite3(matrix[1], "matrix[1]")?;
        require_finite3(matrix[2], "matrix[2]")?;

        let lengths = [norm3(matrix[0]), norm3(matrix[1]), norm3(matrix[2])];
        for (len, axis) in lengths.iter().zip(["a", "b", "c"]) {
            if *len <= Self::MIN_LENGTH {
                return Err(CrystalError::NonPositiveLength { axis, value: *len });
            }
        }

        let signed_vol = signed_volume(&matrix);
        if !signed_vol.is_finite() {
            return Err(CrystalError::NonFiniteVolume);
        }
        let volume = signed_vol.abs();
        if volume == 0.0 {
            return Err(CrystalError::SingularMatrix);
        }

        let bounding_box = lengths[0] * lengths[1] * lengths[2];
        let condition = volume / bounding_box;
        if condition < Self::MIN_CONDITION_INDICATOR {
            return Err(CrystalError::NearSingularMatrix {
                condition,
                threshold: Self::MIN_CONDITION_INDICATOR,
            });
        }

        // Reciprocal rows b1, b2, b3 (a_i . b_j = delta_ij), computed from
        // the *signed* volume so the algebraic identity M . M^-1 = I holds
        // exactly (not just up to sign) for left-handed input matrices too.
        let b1 = scale3(cross3(matrix[1], matrix[2]), 1.0 / signed_vol);
        let b2 = scale3(cross3(matrix[2], matrix[0]), 1.0 / signed_vol);
        let b3 = scale3(cross3(matrix[0], matrix[1]), 1.0 / signed_vol);
        // inverse's *columns* are b1, b2, b3 (see module doc: reciprocal
        // row j = inverse column j).
        let inverse = [
            [b1[0], b2[0], b3[0]],
            [b1[1], b2[1], b3[1]],
            [b1[2], b2[2], b3[2]],
        ];
        for row in inverse {
            require_finite3(row, "inverse")?;
        }

        Ok(Self { matrix, inverse })
    }

    /// Build a lattice from cell parameters: lengths `a, b, c` in Angstrom,
    /// angles `alpha, beta, gamma` in degrees (`alpha` = angle between b
    /// and c, `beta` = angle between a and c, `gamma` = angle between a and
    /// b -- standard crystallographic convention).
    ///
    /// Uses the same IUCr placement `chematic_mol::cif::UnitCell` derives
    /// its Cartesian conversion from (a along x, b in the xy-plane, c
    /// completing the frame), so a CIF-parsed cell and a `from_parameters`
    /// cell built from the same six numbers agree on Cartesian output.
    ///
    /// # Examples
    ///
    /// ```
    /// use chematic_crystal::Lattice;
    ///
    /// // Quartz-like triclinic-ish cell (illustrative numbers, not a
    /// // real quartz refinement).
    /// let lattice = Lattice::from_parameters(4.9, 4.9, 5.4, 90.0, 90.0, 120.0)?;
    /// assert!(lattice.volume() > 0.0);
    /// let [alpha, beta, gamma] = lattice.angles_degrees();
    /// assert!((gamma - 120.0).abs() < 1e-6);
    /// assert!((alpha - 90.0).abs() < 1e-6 && (beta - 90.0).abs() < 1e-6);
    /// # Ok::<(), chematic_crystal::CrystalError>(())
    /// ```
    pub fn from_parameters(
        a: f64,
        b: f64,
        c: f64,
        alpha: f64,
        beta: f64,
        gamma: f64,
    ) -> Result<Self, CrystalError> {
        for (len, axis) in [a, b, c].into_iter().zip(["a", "b", "c"]) {
            require_finite(len, axis)?;
            if len <= Self::MIN_LENGTH {
                return Err(CrystalError::NonPositiveLength { axis, value: len });
            }
        }
        for (angle, name) in [alpha, beta, gamma]
            .into_iter()
            .zip(["alpha", "beta", "gamma"])
        {
            require_finite(angle, name)?;
            if !(angle > 0.0 && angle < 180.0) {
                return Err(CrystalError::InvalidAngle {
                    angle: name,
                    value: angle,
                });
            }
        }

        let (ca, cb, cg, sg) = (
            alpha.to_radians().cos(),
            beta.to_radians().cos(),
            gamma.to_radians().cos(),
            gamma.to_radians().sin(),
        );

        let a_vec = [a, 0.0, 0.0];
        let b_vec = [b * cg, b * sg, 0.0];
        let cx = c * cb;
        let cy = c * (ca - cb * cg) / sg;
        let cz_sq = c * c - cx * cx - cy * cy;
        if cz_sq < -Self::ANGLE_CONSISTENCY_TOLERANCE * c * c {
            return Err(CrystalError::IncompatibleAngles { alpha, beta, gamma });
        }
        let cz = cz_sq.max(0.0).sqrt();
        let c_vec = [cx, cy, cz];

        Self::from_matrix([a_vec, b_vec, c_vec])
    }

    /// A cubic cell with edge length `a` (all angles 90 degrees).
    pub fn cubic(a: f64) -> Result<Self, CrystalError> {
        require_finite(a, "a")?;
        if a <= Self::MIN_LENGTH {
            return Err(CrystalError::NonPositiveLength {
                axis: "a",
                value: a,
            });
        }
        Self::from_matrix([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]])
    }

    /// An orthorhombic cell with edge lengths `a, b, c` (all angles 90
    /// degrees).
    pub fn orthorhombic(a: f64, b: f64, c: f64) -> Result<Self, CrystalError> {
        for (len, axis) in [a, b, c].into_iter().zip(["a", "b", "c"]) {
            require_finite(len, axis)?;
            if len <= Self::MIN_LENGTH {
                return Err(CrystalError::NonPositiveLength { axis, value: len });
            }
        }
        Self::from_matrix([[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]])
    }

    /// The lattice matrix (rows = a, b, c).
    #[inline]
    pub fn matrix(&self) -> [[f64; 3]; 3] {
        self.matrix
    }

    /// The cached matrix inverse (`matrix() . inverse_matrix() == I`).
    #[inline]
    pub fn inverse_matrix(&self) -> [[f64; 3]; 3] {
        self.inverse
    }

    /// Cell volume in cubic Angstrom (always `>= 0`; the signed determinant
    /// is not exposed since it depends on lattice-vector handedness, which
    /// is not a meaningful distinction for volume).
    pub fn volume(&self) -> f64 {
        signed_volume(&self.matrix).abs()
    }

    /// Dimensionless condition indicator `|det(M)| / (|a| |b| |c|)` -- see
    /// [`Self::MIN_CONDITION_INDICATOR`].
    pub fn condition_indicator(&self) -> f64 {
        let lengths = self.lengths();
        self.volume() / (lengths[0] * lengths[1] * lengths[2])
    }

    /// Lattice vector lengths `[|a|, |b|, |c|]`, in Angstrom.
    pub fn lengths(&self) -> [f64; 3] {
        [
            norm3(self.matrix[0]),
            norm3(self.matrix[1]),
            norm3(self.matrix[2]),
        ]
    }

    /// Lattice angles `[alpha, beta, gamma]` in degrees (alpha = angle(b,
    /// c), beta = angle(a, c), gamma = angle(a, b)).
    pub fn angles_degrees(&self) -> [f64; 3] {
        let angle = |u: [f64; 3], v: [f64; 3]| -> f64 {
            let cos_theta = (dot3(u, v) / (norm3(u) * norm3(v))).clamp(-1.0, 1.0);
            cos_theta.acos().to_degrees()
        };
        [
            angle(self.matrix[1], self.matrix[2]), // alpha: b,c
            angle(self.matrix[0], self.matrix[2]), // beta: a,c
            angle(self.matrix[0], self.matrix[1]), // gamma: a,b
        ]
    }

    /// Fractional -> Cartesian: `cartesian = fractional . matrix`.
    ///
    /// # Examples
    ///
    /// ```
    /// use chematic_crystal::{FractionalCoord, Lattice};
    ///
    /// let lattice = Lattice::cubic(4.0)?;
    /// let cart = lattice.frac_to_cart(FractionalCoord::new([0.5, 0.0, 0.0]));
    /// assert!((cart.0[0] - 2.0).abs() < 1e-12);
    /// # Ok::<(), chematic_crystal::CrystalError>(())
    /// ```
    pub fn frac_to_cart(&self, frac: FractionalCoord) -> CartesianCoord {
        let f = frac.0;
        let m = self.matrix;
        CartesianCoord([
            f[0] * m[0][0] + f[1] * m[1][0] + f[2] * m[2][0],
            f[0] * m[0][1] + f[1] * m[1][1] + f[2] * m[2][1],
            f[0] * m[0][2] + f[1] * m[1][2] + f[2] * m[2][2],
        ])
    }

    /// Cartesian -> Fractional: `fractional = cartesian . inverse`.
    pub fn cart_to_frac(&self, cart: CartesianCoord) -> FractionalCoord {
        let r = cart.0;
        let inv = self.inverse;
        FractionalCoord([
            r[0] * inv[0][0] + r[1] * inv[1][0] + r[2] * inv[2][0],
            r[0] * inv[0][1] + r[1] * inv[1][1] + r[2] * inv[2][1],
            r[0] * inv[0][2] + r[1] * inv[1][2] + r[2] * inv[2][2],
        ])
    }

    /// Reciprocal lattice matrix (rows = reciprocal vectors b1, b2, b3,
    /// crystallographic convention `a_i . b_j = delta_ij`, no `2*pi`
    /// factor). Equal to the transpose of [`Self::inverse_matrix`].
    pub fn reciprocal_matrix(&self) -> [[f64; 3]; 3] {
        let inv = self.inverse;
        [
            [inv[0][0], inv[1][0], inv[2][0]],
            [inv[0][1], inv[1][1], inv[2][1]],
            [inv[0][2], inv[1][2], inv[2][2]],
        ]
    }

    /// Row norms `[|b1|, |b2|, |b3|]` of the reciprocal matrix. Used by the
    /// exact bounded-search minimum-image algorithm ([`crate::periodic`])
    /// and cutoff neighbor search ([`crate::neighbor`]) to derive a finite,
    /// provably sufficient search box -- see
    /// `docs/rfcs/chematic_crystal_foundation.md`.
    pub(crate) fn reciprocal_row_norms(&self) -> [f64; 3] {
        let recip = self.reciprocal_matrix();
        [norm3(recip[0]), norm3(recip[1]), norm3(recip[2])]
    }
}

/// Serializes/deserializes only `matrix` -- `inverse` is cached derived
/// state, not persisted (persisting it would let a hand-edited JSON file
/// carry a `matrix`/`inverse` pair that no longer agree). Deserializing
/// re-derives `inverse` and re-runs every constructor validation via
/// [`Lattice::from_matrix`], so a deserialized `Lattice` has exactly the
/// same guarantees as one built through the normal constructors.
#[cfg(feature = "serde")]
mod serde_impl {
    use super::Lattice;
    use crate::error::CrystalError;
    use serde::de::Error as _;
    use serde::ser::SerializeStruct;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    impl Serialize for Lattice {
        fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            let mut state = serializer.serialize_struct("Lattice", 1)?;
            state.serialize_field("matrix", &self.matrix)?;
            state.end()
        }
    }

    #[derive(Deserialize)]
    struct LatticeData {
        matrix: [[f64; 3]; 3],
    }

    impl<'de> Deserialize<'de> for Lattice {
        fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            let data = LatticeData::deserialize(deserializer)?;
            Lattice::from_matrix(data.matrix)
                .map_err(|e: CrystalError| D::Error::custom(e.to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn assert_close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() < tol, "{a} vs {b} (tol {tol})");
    }

    // -- volumes --------------------------------------------------------

    #[test]
    fn cubic_volume() {
        let l = Lattice::cubic(3.0).unwrap();
        assert_close(l.volume(), 27.0, EPS);
    }

    #[test]
    fn orthorhombic_volume() {
        let l = Lattice::orthorhombic(2.0, 3.0, 4.0).unwrap();
        assert_close(l.volume(), 24.0, EPS);
    }

    #[test]
    fn triclinic_volume() {
        // Known cell, cross-checked against the standard IUCr formula
        // V = abc*sqrt(1 - cos^2a - cos^2b - cos^2g + 2 cos(a)cos(b)cos(g)).
        let (a, b, c) = (5.0, 6.0, 7.0);
        let (alpha, beta, gamma) = (80.0, 95.0, 110.0);
        let l = Lattice::from_parameters(a, b, c, alpha, beta, gamma).unwrap();
        let (ca, cb, cg) = (
            alpha.to_radians().cos(),
            beta.to_radians().cos(),
            gamma.to_radians().cos(),
        );
        let expected = a * b * c * (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt();
        assert_close(l.volume(), expected, 1e-6);
    }

    #[test]
    fn from_parameters_known_value_cubic() {
        let l = Lattice::from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0).unwrap();
        assert_close(l.volume(), 64.0, EPS);
        let angles = l.angles_degrees();
        for a in angles {
            assert_close(a, 90.0, 1e-6);
        }
        let lengths = l.lengths();
        for len in lengths {
            assert_close(len, 4.0, 1e-9);
        }
    }

    // -- round-trips ------------------------------------------------------

    #[test]
    fn frac_cart_frac_round_trip_triclinic() {
        let l = Lattice::from_parameters(5.0, 6.0, 7.0, 75.0, 100.0, 60.0).unwrap();
        let f0 = FractionalCoord::new([0.13, 0.71, 0.42]);
        let cart = l.frac_to_cart(f0);
        let f1 = l.cart_to_frac(cart);
        for i in 0..3 {
            assert_close(f0.0[i], f1.0[i], 1e-9);
        }
    }

    #[test]
    fn cart_frac_cart_round_trip_triclinic() {
        let l = Lattice::from_parameters(5.0, 6.0, 7.0, 75.0, 100.0, 60.0).unwrap();
        let c0 = CartesianCoord::new([1.23, -4.56, 2.34]);
        let frac = l.cart_to_frac(c0);
        let c1 = l.frac_to_cart(frac);
        for i in 0..3 {
            assert_close(c0.0[i], c1.0[i], 1e-9);
        }
    }

    // -- reciprocal relation ------------------------------------------------

    #[test]
    // Cross-product of two independent 3x3 matrices' rows, keyed off both
    // numeric indices for the `i == j` Kronecker-delta comparison -- not a
    // single-array walk `enumerate()` would simplify.
    #[allow(clippy::needless_range_loop)]
    fn reciprocal_relation_a_dot_b_is_kronecker_delta() {
        let l = Lattice::from_parameters(5.0, 6.0, 7.0, 80.0, 95.0, 110.0).unwrap();
        let m = l.matrix();
        let recip = l.reciprocal_matrix();
        for i in 0..3 {
            for j in 0..3 {
                let d = dot3(m[i], recip[j]);
                let expected = if i == j { 1.0 } else { 0.0 };
                assert_close(d, expected, 1e-9);
            }
        }
    }

    #[test]
    // See the allow above: `recip[i][j]` vs `inv[j][i]` needs both indices.
    #[allow(clippy::needless_range_loop)]
    fn reciprocal_matrix_is_inverse_transpose() {
        let l = Lattice::from_parameters(5.0, 6.0, 7.0, 80.0, 95.0, 110.0).unwrap();
        let inv = l.inverse_matrix();
        let recip = l.reciprocal_matrix();
        for i in 0..3 {
            for j in 0..3 {
                assert_close(recip[i][j], inv[j][i], 1e-12);
            }
        }
    }

    // -- rejections -------------------------------------------------------

    #[test]
    fn singular_matrix_rejected() {
        // Three coplanar (in fact collinear-ish) rows: c = a + b.
        let m = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        assert_eq!(Lattice::from_matrix(m), Err(CrystalError::SingularMatrix));
    }

    #[test]
    fn near_singular_matrix_rejected() {
        // c almost coplanar with a,b (tiny out-of-plane component).
        let m = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 1e-8]];
        let err = Lattice::from_matrix(m).unwrap_err();
        assert!(matches!(err, CrystalError::NearSingularMatrix { .. }));
    }

    #[test]
    fn nan_rejected() {
        let m = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, f64::NAN]];
        assert!(matches!(
            Lattice::from_matrix(m),
            Err(CrystalError::NonFinite { .. })
        ));
    }

    #[test]
    fn infinity_rejected() {
        let m = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, f64::INFINITY]];
        assert!(matches!(
            Lattice::from_matrix(m),
            Err(CrystalError::NonFinite { .. })
        ));
    }

    #[test]
    fn zero_length_vector_rejected() {
        let m = [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(matches!(
            Lattice::from_matrix(m),
            Err(CrystalError::NonPositiveLength { .. })
        ));
    }

    #[test]
    fn from_parameters_rejects_non_positive_length() {
        assert!(matches!(
            Lattice::from_parameters(0.0, 4.0, 4.0, 90.0, 90.0, 90.0),
            Err(CrystalError::NonPositiveLength { .. })
        ));
        assert!(matches!(
            Lattice::from_parameters(-1.0, 4.0, 4.0, 90.0, 90.0, 90.0),
            Err(CrystalError::NonPositiveLength { .. })
        ));
    }

    #[test]
    fn from_parameters_rejects_out_of_range_angle() {
        assert!(matches!(
            Lattice::from_parameters(4.0, 4.0, 4.0, 0.0, 90.0, 90.0),
            Err(CrystalError::InvalidAngle { .. })
        ));
        assert!(matches!(
            Lattice::from_parameters(4.0, 4.0, 4.0, 180.0, 90.0, 90.0),
            Err(CrystalError::InvalidAngle { .. })
        ));
        assert!(matches!(
            Lattice::from_parameters(4.0, 4.0, 4.0, 190.0, 90.0, 90.0),
            Err(CrystalError::InvalidAngle { .. })
        ));
    }

    #[test]
    fn from_parameters_rejects_incompatible_angles() {
        // Three angles each individually valid (< 180) but jointly
        // impossible for a real 3D cell.
        let err = Lattice::from_parameters(4.0, 4.0, 4.0, 179.0, 179.0, 179.0).unwrap_err();
        assert!(matches!(err, CrystalError::IncompatibleAngles { .. }));
    }

    #[test]
    fn cubic_rejects_non_positive() {
        assert!(matches!(
            Lattice::cubic(0.0),
            Err(CrystalError::NonPositiveLength { .. })
        ));
        assert!(matches!(
            Lattice::cubic(-2.0),
            Err(CrystalError::NonPositiveLength { .. })
        ));
        assert!(matches!(
            Lattice::cubic(f64::NAN),
            Err(CrystalError::NonFinite { .. })
        ));
    }
}
