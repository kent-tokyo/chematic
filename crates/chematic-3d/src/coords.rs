//! 3D coordinate types for molecular structures.

use chematic_core::AtomIdx;

/// A 3D point in Cartesian space (angstroms).
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

    /// Try to normalize to a unit vector, returning None if the vector has zero length.
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
/// Indexed by atom insertion order (`AtomIdx.0` as `usize`).
#[derive(Debug, Clone)]
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
    pub fn get(&self, idx: AtomIdx) -> Point3 {
        self.points[idx.0 as usize]
    }

    /// Set the coordinate of atom `idx`.
    pub fn set(&mut self, idx: AtomIdx, p: Point3) {
        self.points[idx.0 as usize] = p;
    }

    /// Number of atom coordinate slots.
    pub fn atom_count(&self) -> usize {
        self.points.len()
    }
}
