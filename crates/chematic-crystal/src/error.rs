//! Error type for `chematic-crystal`.

use std::fmt;

/// Errors returned by `chematic-crystal`'s validated constructors.
///
/// Every variant that carries a numeric field is worded so the message is
/// self-contained (no separate lookup table needed to interpret it).
#[derive(Debug, Clone, PartialEq)]
pub enum CrystalError {
    /// A field that must be finite (no `NaN`, no `+-Infinity`) was not.
    NonFinite {
        /// Name of the offending field/quantity (e.g. `"matrix[1][2]"`).
        field: &'static str,
    },
    /// A crystallographic length (`a`, `b`, or `c`) was not a positive
    /// finite number.
    NonPositiveLength {
        /// Which lattice vector (`"a"`, `"b"`, or `"c"`).
        axis: &'static str,
        /// The rejected value.
        value: f64,
    },
    /// A crystallographic angle (alpha, beta, gamma) was outside the open
    /// interval `(0, 180)` degrees, i.e. does not define a non-degenerate
    /// parallelepiped.
    InvalidAngle {
        /// Which angle (`"alpha"`, `"beta"`, or `"gamma"`).
        angle: &'static str,
        /// The rejected value, in degrees.
        value: f64,
    },
    /// The 3x3 lattice matrix is exactly singular (zero volume): the three
    /// lattice vectors are linearly dependent.
    SingularMatrix,
    /// The lattice matrix is not exactly singular but is numerically too
    /// close to it for a stable inverse -- see
    /// [`Lattice::MIN_CONDITION_INDICATOR`](crate::lattice::Lattice::MIN_CONDITION_INDICATOR)
    /// for the threshold and its rationale.
    NearSingularMatrix {
        /// The computed `|det(M)| / (|a| |b| |c|)` condition indicator.
        condition: f64,
        /// The rejection threshold it fell below.
        threshold: f64,
    },
    /// The lattice volume computed to a non-finite value (can only happen
    /// from non-finite inputs that individually passed other checks, e.g.
    /// via extreme cancellation -- guarded defensively).
    NonFiniteVolume,
    /// An [`Occupancy`](crate::site::Occupancy) value was negative.
    NegativeOccupancy {
        /// The rejected value.
        value: f64,
    },
    /// An [`Occupancy`](crate::site::Occupancy) value was non-finite.
    NonFiniteOccupancy,
    /// The sum of occupancies at one [`PeriodicSite`](crate::site::PeriodicSite)
    /// exceeded `1.0` by more than the configured tolerance.
    OccupancySumExceeded {
        /// The computed sum.
        sum: f64,
        /// The tolerance above `1.0` that was allowed.
        tolerance: f64,
    },
    /// A [`PeriodicSite`](crate::site::PeriodicSite) was constructed with an
    /// empty species list.
    EmptySpeciesList,
    /// Lattice angles that individually satisfy `InvalidAngle`'s per-angle
    /// range but do not jointly define a real (non-degenerate) 3D
    /// parallelepiped -- e.g. `alpha=170, beta=170, gamma=170` fails the
    /// triangle-like consistency constraint between the three angles.
    IncompatibleAngles {
        /// alpha, in degrees.
        alpha: f64,
        /// beta, in degrees.
        beta: f64,
        /// gamma, in degrees.
        gamma: f64,
    },
    /// A [`PeriodicSite`](crate::site::PeriodicSite) at the given index
    /// failed its own validation; wraps the underlying reason.
    InvalidSite {
        /// Index of the offending site (in structure order).
        index: usize,
        /// The underlying validation failure.
        source: Box<CrystalError>,
    },
    /// A diagonal supercell multiplier was not `>= 1`.
    NonPositiveSupercellMultiplier {
        /// Which axis (`0`, `1`, or `2`).
        axis: usize,
        /// The rejected value.
        value: u32,
    },
    /// A requested cutoff-radius neighbor search was rejected: the exact
    /// bounded image search implies more candidate periodic images than
    /// [`neighbor::MAX_NEIGHBOR_IMAGE_CANDIDATES`](crate::neighbor::MAX_NEIGHBOR_IMAGE_CANDIDATES)
    /// -- almost always a cutoff far larger than the cell (a likely unit or
    /// input error) rather than a legitimate search.
    NeighborSearchTooLarge {
        /// The computed candidate-image count that triggered the guard.
        candidate_count: u64,
        /// The configured limit.
        limit: u64,
    },
    /// A cutoff radius was not finite and positive.
    InvalidCutoff {
        /// The rejected value.
        value: f64,
    },
}

impl fmt::Display for CrystalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CrystalError::NonFinite { field } => {
                write!(f, "{field} must be finite (got NaN or Infinity)")
            }
            CrystalError::NonPositiveLength { axis, value } => {
                write!(
                    f,
                    "lattice length {axis} must be a positive finite number, got {value}"
                )
            }
            CrystalError::InvalidAngle { angle, value } => {
                write!(
                    f,
                    "lattice angle {angle} must be in the open interval (0, 180) degrees, got {value}"
                )
            }
            CrystalError::SingularMatrix => {
                write!(
                    f,
                    "lattice matrix is singular (zero volume): lattice vectors are linearly dependent"
                )
            }
            CrystalError::NearSingularMatrix {
                condition,
                threshold,
            } => {
                write!(
                    f,
                    "lattice matrix is near-singular: condition indicator {condition} is below the minimum {threshold}"
                )
            }
            CrystalError::NonFiniteVolume => {
                write!(f, "lattice volume computed to a non-finite value")
            }
            CrystalError::NegativeOccupancy { value } => {
                write!(f, "occupancy must be >= 0, got {value}")
            }
            CrystalError::NonFiniteOccupancy => {
                write!(f, "occupancy must be finite (got NaN or Infinity)")
            }
            CrystalError::OccupancySumExceeded { sum, tolerance } => {
                write!(
                    f,
                    "species occupancies sum to {sum}, which exceeds 1.0 + tolerance ({tolerance})"
                )
            }
            CrystalError::EmptySpeciesList => {
                write!(f, "species list must not be empty")
            }
            CrystalError::IncompatibleAngles { alpha, beta, gamma } => {
                write!(
                    f,
                    "angles alpha={alpha}, beta={beta}, gamma={gamma} do not define a valid (non-degenerate) parallelepiped cell"
                )
            }
            CrystalError::InvalidSite { index, source } => {
                write!(f, "site {index}: {source}")
            }
            CrystalError::NonPositiveSupercellMultiplier { axis, value } => {
                write!(
                    f,
                    "supercell multiplier for axis {axis} must be >= 1, got {value}"
                )
            }
            CrystalError::NeighborSearchTooLarge {
                candidate_count,
                limit,
            } => {
                write!(
                    f,
                    "neighbor search would examine {candidate_count} candidate periodic images, exceeding the limit of {limit} -- check that the cutoff is not far larger than the cell"
                )
            }
            CrystalError::InvalidCutoff { value } => {
                write!(f, "cutoff must be a finite, positive number, got {value}")
            }
        }
    }
}

impl std::error::Error for CrystalError {}
