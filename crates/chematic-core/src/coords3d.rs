//! Minimal, algorithm-free 3D coordinate storage shared across crates.
//!
//! [`Coords3D`] holds one [`Point3`] per atom, indexed by atom insertion
//! order (`AtomIdx.0` as `usize`). It carries no embedding, force-field, or
//! alignment logic -- see the 3D Breakthrough Program's master plan
//! (`docs/rfcs/3d_breakthrough_master_plan.md`, decision 1a) for why this type
//! lives in `chematic-core` rather than `chematic-mol` importing
//! `chematic_3d::coords::Coords3D` directly: `chematic-3d` pulls in
//! `chematic-ff`/`chematic-chem`/`chematic-fp`/`chematic-smarts` transitively,
//! an unwanted dependency footprint for `chematic-mol` and every other
//! `chematic-core` consumer that never touches 3D generation. Every crate
//! already depends on `chematic-core`, so this requires zero new dependency
//! edges.
//!
//! **Field layout and method names deliberately mirror
//! `chematic_3d::coords::{Point3, Coords3D}` as closely as possible** (same
//! `x`/`y`/`z` fields, same `points: Vec<Point3>` layout, same
//! `new_zeroed`/`get`/`set`/`atom_count` names and panics-on-out-of-range
//! indexing convention) so that a later Coordinator-authored bridge PR can
//! make `chematic_3d::coords` re-export this type
//! (`pub use chematic_core::{Coords3D, Point3};`) with minimal changes to
//! existing `chematic-3d` call sites. `chematic-3d`'s own `Coords3D` is not
//! touched by this PR (read-only reference); the two types are reconciled at
//! Wave 1->2 integration time, not here.
//!
//! This module adds two things `chematic_3d::coords` doesn't have:
//! - `is_finite()` on both types. **Not currently called by this PR's own
//!   reader code** -- `chematic-mol`'s V2000/V3000 z-coordinate parsers
//!   validate NaN/Inf with a direct `f64::is_finite()` check on the raw
//!   parsed value, *before* a `Point3`/`Coords3D` is ever constructed (see
//!   `mol2000.rs`/`mol3000.rs`), so these methods have zero non-test callers
//!   today. Provided for future consumers (e.g. Wave 2 work that may
//!   construct a `Coords3D` from elsewhere and want to validate it after the
//!   fact) rather than because this PR's own code needs them.
//! - `Default`/`PartialEq` derives on `Coords3D` (`chematic_3d::coords::Coords3D`
//!   derives only `Debug, Clone`) -- needed for this crate's own tests and
//!   generally harmless additions for a plain data holder.

use crate::molecule::AtomIdx;

/// A 3D point in Cartesian space (Angstrom).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    /// Create a new point.
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The origin (0, 0, 0).
    #[inline]
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// `true` when all three components are finite (not NaN, not Infinite).
    #[inline]
    pub fn is_finite(&self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Euclidean distance to another point.
    #[inline]
    pub fn distance(&self, other: &Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Component-wise addition.
    #[inline]
    pub fn add(&self, other: &Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    /// Component-wise subtraction.
    #[inline]
    pub fn sub(&self, other: &Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    /// Scalar multiplication.
    #[inline]
    pub fn scale(&self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }

    /// Dot product.
    #[inline]
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product.
    #[inline]
    pub fn cross(&self, other: &Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Euclidean norm (length).
    #[inline]
    pub fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Normalize to a unit vector.
    ///
    /// # Panics
    /// Panics if the vector has zero length.
    pub fn normalize(&self) -> Self {
        let n = self.norm();
        assert!(n > 0.0, "cannot normalize a zero-length vector");
        self.scale(1.0 / n)
    }

    /// Try to normalize to a unit vector, returning `None` if the vector has
    /// zero length.
    pub fn try_normalize(&self) -> Option<Self> {
        let n = self.norm();
        if n > 0.0 {
            Some(self.scale(1.0 / n))
        } else {
            None
        }
    }
}

/// 3D coordinates for all heavy atoms in a molecule.
///
/// Indexed by atom insertion order (`AtomIdx.0` as `usize`). A dumb data
/// holder only -- no embedding, minimization, or alignment logic (those stay
/// in `chematic-3d`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Coords3D {
    pub points: Vec<Point3>,
}

impl Coords3D {
    /// Create zeroed coordinates for `n` atoms.
    pub fn new_zeroed(n: usize) -> Self {
        Self {
            points: vec![Point3::zero(); n],
        }
    }

    /// Get the coordinate of atom `idx`.
    ///
    /// # Panics
    /// Panics if `idx` is out of range (matching `chematic_3d::coords::Coords3D`'s
    /// existing indexing convention).
    pub fn get(&self, idx: AtomIdx) -> Point3 {
        self.points[idx.0 as usize]
    }

    /// Set the coordinate of atom `idx`.
    ///
    /// # Panics
    /// Panics if `idx` is out of range.
    pub fn set(&mut self, idx: AtomIdx, p: Point3) {
        self.points[idx.0 as usize] = p;
    }

    /// Number of atom coordinate slots.
    pub fn atom_count(&self) -> usize {
        self.points.len()
    }

    /// `true` when every point has all-finite components.
    pub fn is_finite(&self) -> bool {
        self.points.iter().all(Point3::is_finite)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_zeroed_has_n_points_at_origin() {
        let c = Coords3D::new_zeroed(3);
        assert_eq!(c.atom_count(), 3);
        for i in 0..3 {
            assert_eq!(c.get(AtomIdx(i)), Point3::zero());
        }
    }

    #[test]
    fn get_set_roundtrip() {
        let mut c = Coords3D::new_zeroed(2);
        c.set(AtomIdx(1), Point3::new(1.0, 2.0, 3.0));
        assert_eq!(c.get(AtomIdx(1)), Point3::new(1.0, 2.0, 3.0));
        assert_eq!(c.get(AtomIdx(0)), Point3::zero());
    }

    #[test]
    fn is_finite_true_for_normal_coords() {
        let mut c = Coords3D::new_zeroed(1);
        c.set(AtomIdx(0), Point3::new(1.0, -2.5, 0.0));
        assert!(c.is_finite());
    }

    #[test]
    fn is_finite_false_for_nan() {
        let mut c = Coords3D::new_zeroed(1);
        c.set(AtomIdx(0), Point3::new(f64::NAN, 0.0, 0.0));
        assert!(!c.is_finite());
    }

    #[test]
    fn is_finite_false_for_infinite() {
        let mut c = Coords3D::new_zeroed(1);
        c.set(AtomIdx(0), Point3::new(0.0, f64::INFINITY, 0.0));
        assert!(!c.is_finite());
    }

    #[test]
    fn point3_basic_vector_ops() {
        let a = Point3::new(1.0, 0.0, 0.0);
        let b = Point3::new(0.0, 1.0, 0.0);
        assert_eq!(a.cross(&b), Point3::new(0.0, 0.0, 1.0));
        assert_eq!(a.dot(&b), 0.0);
        assert_eq!(a.distance(&b), std::f64::consts::SQRT_2);
        assert_eq!(a.add(&b), Point3::new(1.0, 1.0, 0.0));
    }
}
