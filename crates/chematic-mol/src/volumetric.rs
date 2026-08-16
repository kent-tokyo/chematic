//! Shared volumetric scalar-field grid type for [`crate::cube`] (Gaussian
//! Cube) and [`crate::opendx`] (the APBS/electrostatics subset of OpenDX).
//!
//! This type was designed independently for chematic; it is not copied from
//! any other tool's data model. It intentionally does **not** reuse
//! `chematic_crystal::Lattice`: a lattice's 3x3 matrix describes a unit
//! cell's edge vectors (with positive-length/non-degenerate-angle
//! validation and a cached inverse), while [`VolumetricGrid::axes`]
//! describes per-voxel step vectors on an open (non-periodic) grid -- a
//! different concept that happens to share the "3x3 matrix + origin" shape.
//! `chematic-mol` also does not depend on `chematic-crystal` outside its
//! optional `crystal` feature, and this type must work without it.
//!
//! ## Value ordering
//!
//! `values` is a flat `Vec<f64>` in **row-major, third-axis-fastest**
//! order: `index = (i * shape[1] + j) * shape[2] + k` for grid coordinates
//! `(i, j, k)` (see [`VolumetricGrid::checked_index`]). This is not an arbitrary
//! third convention -- it is the *native* storage order of both formats
//! this module serves:
//! - Gaussian Cube: "the x axis as the outer loop and the z axis as the
//!   inner loop" (<https://paulbourke.net/dataformats/cube/>).
//! - OpenDX/APBS: "the data values, ordered with the z-index increasing
//!   most quickly, followed by the y-index, and then the x-index"
//!   (<https://apbs.readthedocs.io/en/latest/formats/opendx.html>).
//!
//! Choosing this order means neither `cube.rs` nor `opendx.rs` needs to
//! transpose data relative to how its source file stores it.
//!
//! ## Units
//!
//! Gaussian Cube has a real ambiguity over whether a file's `origin`/`axes`
//! (and, by this module's extension, atom positions) are in Bohr or
//! Ångström -- see `cube.rs`'s module docs for the full citation trail.
//! Rather than silently normalizing to one unit (and possibly guessing
//! wrong for an ambiguous file) or silently assuming Bohr always, this type
//! carries an explicit [`GridUnits`] tag: numbers are stored exactly as
//! read, never converted, and [`GridUnits`] records which unit they are in.
//!
//! ## Multiple datasets per voxel
//!
//! Gaussian Cube has a real, documented convention for storing more than
//! one value per voxel (e.g. several molecular-orbital cubes in one file):
//! either a negative `NAtoms` header field followed by a dataset-identifier
//! line, or a positive `NAtoms` with an `NVal != 1` trailing field on the
//! same line. `VolumetricGrid` is single-valued-per-voxel only; `cube.rs`
//! detects both forms and rejects them with a typed
//! [`crate::cube::CubeError::MultiDatasetUnsupported`] rather than silently
//! reading only the first dataset or misinterpreting the data layout.
//!
//! ## Out of scope
//!
//! VASP's CHGCAR/LOCPOT volumetric formats are explicitly out of scope for
//! this module (a later roadmap step) -- they would naturally build on this
//! same `VolumetricGrid` type once a need arises, but are not implemented
//! here.

use chematic_core::{Atom, Element, Molecule, MoleculeBuilder};

// ---------------------------------------------------------------------------
// GridUnits
// ---------------------------------------------------------------------------

/// The length unit `VolumetricGrid::origin`, `VolumetricGrid::axes`, and
/// every [`GridAtom::position`] are expressed in. See module docs for why
/// this is an explicit tag rather than a silent normalization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GridUnits {
    /// Atomic units (a0). Gaussian Cube's native/default unit.
    Bohr,
    /// Ångström.
    Angstrom,
}

// ---------------------------------------------------------------------------
// GridAtom
// ---------------------------------------------------------------------------

/// One atom carried alongside a [`VolumetricGrid`] (Cube's per-file atom
/// list; always empty for OpenDX, which has no atom section at all -- see
/// `opendx.rs` module docs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridAtom {
    pub element: Element,
    /// Cube's "charge" atom-line field: the effective **nuclear** charge,
    /// not a partial/formal/ionic charge. Equal to
    /// `element.atomic_number() as f64` for a standard all-electron atom;
    /// differs when the source calculation used an effective core
    /// potential (ECP/pseudopotential) that replaces some core electrons.
    /// Always 0.0 for OpenDX (no atom section to carry a value at all).
    pub charge: f64,
    /// Same unit as the parent [`VolumetricGrid::units`].
    pub position: [f64; 3],
}

// ---------------------------------------------------------------------------
// GridError — structural invariants shared by both formats
// ---------------------------------------------------------------------------

/// Structural-invariant errors shared by [`crate::cube::CubeError`] and
/// [`crate::opendx::OpenDxError`] (each wraps this in its own `Grid`
/// variant). Parsers reject a pathological header via these *before*
/// attempting any large allocation or read; writers call
/// [`VolumetricGrid::validate`] (which returns these) before serializing,
/// since a hand-built `VolumetricGrid` is not guaranteed to satisfy them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GridError {
    /// `shape[0] * shape[1] * shape[2]` overflows `usize`.
    ShapeOverflow { shape: [usize; 3] },
    /// The grid-point count is within `usize` range but exceeds the
    /// caller's configured limit.
    GridTooLarge { points: usize, limit: usize },
    /// A `shape` dimension is `0` -- a degenerate grid that no format this
    /// module writes can represent losslessly (Cube's axis-line voxel
    /// count must be strictly positive; a `0`-length axis leaves nothing
    /// for this crate's own parser to read back).
    ZeroDimension { axis: usize },
    /// `values.len()` does not equal `shape[0] * shape[1] * shape[2]`.
    ValueCountMismatch { expected: usize, found: usize },
    /// An `origin` component is NaN/Infinite.
    NonFiniteOrigin { component: usize, value: f64 },
    /// An `axes` component is NaN/Infinite.
    NonFiniteAxisVector {
        axis: usize,
        component: usize,
        value: f64,
    },
    /// A `values` entry is NaN/Infinite.
    NonFiniteValue { index: usize, value: f64 },
    /// A [`GridAtom::position`] component is NaN/Infinite.
    NonFiniteAtomPosition {
        atom_index: usize,
        component: usize,
        value: f64,
    },
    /// A [`GridAtom::charge`] is NaN/Infinite.
    NonFiniteAtomCharge { atom_index: usize, value: f64 },
}

impl std::fmt::Display for GridError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ShapeOverflow { shape } => {
                write!(f, "grid shape {shape:?} overflows usize when multiplied")
            }
            Self::GridTooLarge { points, limit } => {
                write!(
                    f,
                    "grid has {points} points, exceeding the {limit}-point limit"
                )
            }
            Self::ZeroDimension { axis } => {
                write!(f, "shape axis {axis} has a zero-length dimension")
            }
            Self::ValueCountMismatch { expected, found } => write!(
                f,
                "grid shape implies {expected} values but {found} were found"
            ),
            Self::NonFiniteOrigin { component, value } => {
                write!(f, "origin component {component} is not finite: {value}")
            }
            Self::NonFiniteAxisVector {
                axis,
                component,
                value,
            } => write!(
                f,
                "axis {axis} component {component} is not finite: {value}"
            ),
            Self::NonFiniteValue { index, value } => {
                write!(f, "grid value at flat index {index} is not finite: {value}")
            }
            Self::NonFiniteAtomPosition {
                atom_index,
                component,
                value,
            } => write!(
                f,
                "atom {atom_index} position component {component} is not finite: {value}"
            ),
            Self::NonFiniteAtomCharge { atom_index, value } => {
                write!(f, "atom {atom_index} charge is not finite: {value}")
            }
        }
    }
}

impl std::error::Error for GridError {}

// ---------------------------------------------------------------------------
// VolumetricGrid
// ---------------------------------------------------------------------------

/// A scalar field sampled on a regular (possibly non-orthogonal) 3D grid,
/// shared by [`crate::cube`] and [`crate::opendx`]. See module docs for the
/// value-ordering, unit, and multi-dataset design decisions.
#[derive(Debug, Clone, PartialEq)]
pub struct VolumetricGrid {
    /// Grid origin (position of voxel `(0, 0, 0)`), in `units`.
    pub origin: [f64; 3],
    /// Per-axis step vectors, in `units`. General 3-vectors, not
    /// necessarily axis-aligned -- both Cube and OpenDX permit
    /// non-orthogonal grids.
    pub axes: [[f64; 3]; 3],
    /// Voxel counts along each axis, `[nx, ny, nz]`.
    pub shape: [usize; 3],
    /// Flat scalar-field samples; see module docs for the index ordering.
    pub values: Vec<f64>,
    /// Atoms carried alongside the grid (Cube only; always empty for
    /// OpenDX).
    pub atoms: Vec<GridAtom>,
    /// The length unit `origin`, `axes`, and every atom position are
    /// expressed in.
    pub units: GridUnits,
}

impl VolumetricGrid {
    /// `shape[0] * shape[1] * shape[2]`, computed with checked arithmetic.
    pub fn point_count(&self) -> Result<usize, GridError> {
        let [nx, ny, nz] = self.shape;
        nx.checked_mul(ny)
            .and_then(|v| v.checked_mul(nz))
            .ok_or(GridError::ShapeOverflow { shape: self.shape })
    }

    /// Flat index into `values` for grid coordinates `(i, j, k)`, or `None`
    /// if the coordinates are out of bounds *or* the index arithmetic
    /// itself would overflow `usize`. Every `VolumetricGrid` field is
    /// `pub`, so a caller can construct a grid whose `shape` product
    /// overflows `usize` even though any individual `(i, j, k)` looks
    /// in-bounds -- this is checked arithmetic throughout specifically so
    /// that case can never panic. See module docs for why the third axis
    /// (`k`) varies fastest.
    pub fn checked_index(&self, i: usize, j: usize, k: usize) -> Option<usize> {
        if i >= self.shape[0] || j >= self.shape[1] || k >= self.shape[2] {
            return None;
        }
        i.checked_mul(self.shape[1])?
            .checked_add(j)?
            .checked_mul(self.shape[2])?
            .checked_add(k)
    }

    /// The value at grid coordinates `(i, j, k)`, or `None` if out of
    /// bounds (or, per [`Self::checked_index`], if the index arithmetic
    /// would overflow).
    pub fn get(&self, i: usize, j: usize, k: usize) -> Option<f64> {
        let idx = self.checked_index(i, j, k)?;
        self.values.get(idx).copied()
    }

    /// Validate structural invariants: the shape doesn't overflow,
    /// `values.len()` matches the shape's point count, and every numeric
    /// field (origin, axes, values, atom positions/charges) is finite.
    /// Called by both `write_cube` and `write_opendx` before serializing --
    /// only a grid produced by this module's own parsers is guaranteed to
    /// satisfy these already.
    pub fn validate(&self) -> Result<(), GridError> {
        for (axis, dim) in self.shape.iter().enumerate() {
            if *dim == 0 {
                return Err(GridError::ZeroDimension { axis });
            }
        }
        let expected = self.point_count()?;
        if self.values.len() != expected {
            return Err(GridError::ValueCountMismatch {
                expected,
                found: self.values.len(),
            });
        }
        for (component, value) in self.origin.iter().enumerate() {
            if !value.is_finite() {
                return Err(GridError::NonFiniteOrigin {
                    component,
                    value: *value,
                });
            }
        }
        for (axis, row) in self.axes.iter().enumerate() {
            for (component, value) in row.iter().enumerate() {
                if !value.is_finite() {
                    return Err(GridError::NonFiniteAxisVector {
                        axis,
                        component,
                        value: *value,
                    });
                }
            }
        }
        for (index, value) in self.values.iter().enumerate() {
            if !value.is_finite() {
                return Err(GridError::NonFiniteValue {
                    index,
                    value: *value,
                });
            }
        }
        for (atom_index, atom) in self.atoms.iter().enumerate() {
            for (component, value) in atom.position.iter().enumerate() {
                if !value.is_finite() {
                    return Err(GridError::NonFiniteAtomPosition {
                        atom_index,
                        component,
                        value: *value,
                    });
                }
            }
            if !atom.charge.is_finite() {
                return Err(GridError::NonFiniteAtomCharge {
                    atom_index,
                    value: atom.charge,
                });
            }
        }
        Ok(())
    }

    /// Build a plain [`Molecule`] + per-atom Cartesian coordinates (in
    /// `self.units`) from `atoms`, in file order. No bonds are added --
    /// neither Cube nor OpenDX carries a bond table. Follows the same
    /// shape as `PqrResult::to_molecule`/`OrcaCoords::to_molecule`
    /// elsewhere in this crate.
    pub fn to_molecule(&self) -> (Molecule, Vec<(f64, f64, f64)>) {
        let mut builder = MoleculeBuilder::new();
        let mut coords = Vec::with_capacity(self.atoms.len());
        for a in &self.atoms {
            builder.add_atom(Atom::new(a.element));
            coords.push((a.position[0], a.position[1], a.position[2]));
        }
        (builder.build(), coords)
    }
}

// ---------------------------------------------------------------------------
// LineFeed — shared line/token source for cube.rs and opendx.rs
// ---------------------------------------------------------------------------

/// Minimal line-then-token feed shared by the Cube and OpenDX parsers:
/// fixed positional header lines are read via [`Self::line`], then the
/// remaining input is flattened into whitespace-delimited tokens via
/// [`Self::token`] for the voxel/data block, which both formats allow to
/// wrap across an unspecified number of lines per line (never hard-coded
/// to e.g. exactly 6 values per line). Generic over the caller's error
/// type `E` so it carries no format-specific knowledge; both a whole-string
/// `&str` source (no IO errors possible) and a streaming `BufRead` source
/// (IO errors mapped to `E` by the caller) plug into the same `next_line`
/// closure shape, so the actual grammar logic is written once and shared
/// by both entry points.
pub(crate) struct LineFeed<E, F: FnMut() -> Result<Option<String>, E>> {
    next_line: F,
    pub(crate) line_no: usize,
    tokens: std::vec::IntoIter<String>,
}

impl<E, F: FnMut() -> Result<Option<String>, E>> LineFeed<E, F> {
    pub(crate) fn new(next_line: F) -> Self {
        Self {
            next_line,
            line_no: 0,
            tokens: Vec::new().into_iter(),
        }
    }

    /// Pull the next raw line, for positional header parsing. Discards any
    /// unconsumed tokens left over from a partially-read line (callers
    /// should not mix `line()` and `token()` mid-line).
    pub(crate) fn line(&mut self) -> Result<Option<String>, E> {
        self.tokens = Vec::new().into_iter();
        let l = (self.next_line)()?;
        if l.is_some() {
            self.line_no += 1;
        }
        Ok(l)
    }

    /// Pull the next whitespace-delimited token, spanning line boundaries
    /// as needed, for the voxel/data block.
    pub(crate) fn token(&mut self) -> Result<Option<String>, E> {
        loop {
            if let Some(t) = self.tokens.next() {
                return Ok(Some(t));
            }
            match (self.next_line)()? {
                None => return Ok(None),
                Some(line) => {
                    self.line_no += 1;
                    self.tokens = line
                        .split_whitespace()
                        .map(str::to_string)
                        .collect::<Vec<_>>()
                        .into_iter();
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny_grid() -> VolumetricGrid {
        VolumetricGrid {
            origin: [0.0, 0.0, 0.0],
            axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            shape: [2, 2, 2],
            values: vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            atoms: Vec::new(),
            units: GridUnits::Angstrom,
        }
    }

    #[test]
    fn index_matches_third_axis_fastest_ordering() {
        let g = tiny_grid();
        // (i, j, k) -> (i*2 + j)*2 + k
        assert_eq!(g.checked_index(0, 0, 0), Some(0));
        assert_eq!(g.checked_index(0, 0, 1), Some(1));
        assert_eq!(g.checked_index(0, 1, 0), Some(2));
        assert_eq!(g.checked_index(1, 0, 0), Some(4));
        assert_eq!(g.get(1, 0, 0), Some(4.0));
        assert_eq!(g.get(2, 0, 0), None);
    }

    #[test]
    fn get_on_overflowing_shape_returns_none_not_panic() {
        // Every field is `pub`, so a caller can hand-build a grid whose
        // shape product overflows usize even though (i, j, k) itself looks
        // in-bounds -- checked_index/get must return None, never panic via
        // the internal multiply.
        let mut g = tiny_grid();
        g.shape = [usize::MAX, usize::MAX, 2];
        g.values = Vec::new();
        assert_eq!(g.checked_index(1, 1, 1), None);
        assert_eq!(g.get(1, 1, 1), None);
    }

    #[test]
    fn point_count_overflow_is_typed_error() {
        let mut g = tiny_grid();
        g.shape = [usize::MAX, usize::MAX, 2];
        assert_eq!(
            g.point_count(),
            Err(GridError::ShapeOverflow { shape: g.shape })
        );
    }

    #[test]
    fn validate_rejects_zero_dimension() {
        let mut g = tiny_grid();
        g.shape = [0, 2, 2];
        g.values = Vec::new();
        assert_eq!(g.validate(), Err(GridError::ZeroDimension { axis: 0 }));
    }

    #[test]
    fn validate_rejects_value_count_mismatch() {
        let mut g = tiny_grid();
        g.values.pop();
        assert_eq!(
            g.validate(),
            Err(GridError::ValueCountMismatch {
                expected: 8,
                found: 7
            })
        );
    }

    #[test]
    fn validate_rejects_nan_value() {
        // NaN != NaN, so this can't be an `assert_eq!` against a literal
        // `GridError::NonFiniteValue { value: f64::NAN, .. }` -- match on
        // shape instead and check `is_nan()` separately.
        let mut g = tiny_grid();
        g.values[3] = f64::NAN;
        match g.validate() {
            Err(GridError::NonFiniteValue { index, value }) => {
                assert_eq!(index, 3);
                assert!(value.is_nan());
            }
            other => panic!("expected NonFiniteValue, got {other:?}"),
        }
    }

    #[test]
    fn validate_rejects_infinite_origin() {
        let mut g = tiny_grid();
        g.origin[1] = f64::INFINITY;
        assert_eq!(
            g.validate(),
            Err(GridError::NonFiniteOrigin {
                component: 1,
                value: f64::INFINITY
            })
        );
    }

    #[test]
    fn validate_passes_for_well_formed_grid() {
        assert_eq!(tiny_grid().validate(), Ok(()));
    }
}
