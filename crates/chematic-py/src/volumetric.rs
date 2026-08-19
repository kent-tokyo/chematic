//! Python bindings for `chematic-mol`'s shared `VolumetricGrid` type
//! (Gaussian Cube `.cube`/`.cub` and OpenDX `.dx` scalar-field grids).
//!
//! One pyclass, [`PyVolumetricGrid`] (`VolumetricGrid` in Python), shared by
//! both formats -- they read/write the same underlying Rust type
//! (`chematic_mol::VolumetricGrid`), which has real behavior worth exposing
//! as methods (`get`/`checked_index`/`point_count`/`to_molecule`) and a
//! large array-shaped `values` field, unlike the "read a molecule + some
//! metadata" formats in `formats.rs` that return a plain dict. This follows
//! `crystal.rs`'s precedent (a dedicated pyclass for a Rust type with real
//! methods) rather than `formats.rs`'s `parse_cif`-style plain-function/dict
//! precedent.
//!
//! `values` (a potentially large flat `f64` array) and `axes` (a `(3, 3)`
//! matrix) are returned as numpy arrays, matching `crystal.rs`'s
//! `Lattice.matrix`/`PeriodicStructure.cartesian_positions` convention for
//! comparable data -- a fresh copy per call, not zero-copy, same tradeoff
//! `crystal.rs` already made.
//!
//! `CubeFileReader` (streaming *input* reading) is intentionally not bound:
//! it still materializes one full `VolumetricGrid` in memory either way (a
//! Cube file is a single grid, not a multi-record stream -- see that type's
//! own Rust docs), so the only thing streaming saves a Python caller is
//! avoiding the doubled peak memory of first reading the whole file into a
//! `str`. Python callers who need that already have `open(path).read()`
//! either way; true zero-copy streaming for a single-grid format would
//! require a much larger binding surface for a marginal benefit.

use std::sync::Arc;

use ndarray::{Array1, Array2, Array3};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArray3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::Mol;

fn cube_err(e: chematic_mol::CubeError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn opendx_err(e: chematic_mol::OpenDxError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn grid_err(e: chematic_mol::GridError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn units_to_str(u: chematic_mol::GridUnits) -> &'static str {
    match u {
        chematic_mol::GridUnits::Bohr => "bohr",
        chematic_mol::GridUnits::Angstrom => "angstrom",
    }
}

fn units_from_str(s: &str) -> PyResult<chematic_mol::GridUnits> {
    match s {
        "bohr" => Ok(chematic_mol::GridUnits::Bohr),
        "angstrom" => Ok(chematic_mol::GridUnits::Angstrom),
        other => Err(PyValueError::new_err(format!(
            "unknown units {other:?}, expected 'bohr' or 'angstrom'"
        ))),
    }
}

fn element_from_symbol(sym: &str) -> PyResult<chematic_core::Element> {
    chematic_core::Element::from_symbol(sym)
        .ok_or_else(|| PyValueError::new_err(format!("unknown element symbol: {sym:?}")))
}

fn axes_to_pyarray<'py>(py: Python<'py>, axes: [[f64; 3]; 3]) -> Bound<'py, PyArray2<f64>> {
    let flat: Vec<f64> = axes.into_iter().flatten().collect();
    Array2::from_shape_vec((3, 3), flat)
        .expect("3x3 fixed-size matrix always reshapes cleanly")
        .into_pyarray(py)
}

/// A scalar field sampled on a regular (possibly non-orthogonal) 3D grid,
/// shared by Gaussian Cube and OpenDX. `values` is a flat array in
/// row-major, third-axis-fastest order: `index = (i * shape[1] + j) *
/// shape[2] + k` for grid coordinates `(i, j, k)` -- see :meth:`get`/
/// :meth:`checked_index`.
#[pyclass(name = "VolumetricGrid")]
pub struct PyVolumetricGrid {
    inner: chematic_mol::VolumetricGrid,
}

#[pymethods]
impl PyVolumetricGrid {
    /// Construct directly. Routes through `VolumetricGrid::validate`, so an
    /// inconsistent shape/value-count, non-finite number, or unknown
    /// element symbol raises `ValueError` here rather than constructing an
    /// invalid grid.
    #[new]
    #[pyo3(signature = (origin, axes, shape, values, atoms=Vec::new(), units="angstrom"))]
    fn new(
        origin: (f64, f64, f64),
        axes: [[f64; 3]; 3],
        shape: (usize, usize, usize),
        values: Vec<f64>,
        atoms: Vec<(String, f64, (f64, f64, f64))>,
        units: &str,
    ) -> PyResult<Self> {
        let units = units_from_str(units)?;
        let atoms = atoms
            .into_iter()
            .map(|(sym, charge, pos)| -> PyResult<chematic_mol::GridAtom> {
                Ok(chematic_mol::GridAtom {
                    element: element_from_symbol(&sym)?,
                    charge,
                    position: [pos.0, pos.1, pos.2],
                })
            })
            .collect::<PyResult<_>>()?;
        let inner = chematic_mol::VolumetricGrid {
            origin: [origin.0, origin.1, origin.2],
            axes,
            shape: [shape.0, shape.1, shape.2],
            values,
            atoms,
            units,
        };
        inner.validate().map_err(grid_err)?;
        Ok(PyVolumetricGrid { inner })
    }

    /// Parse a Gaussian Cube file (text, not a path).
    ///
    /// Args:
    ///     text: `.cube`/`.cub` file contents.
    ///     max_input_bytes: byte-size limit (default 1 GiB).
    ///     max_atoms: atom-count limit (default 200,000).
    ///     max_grid_points: cap on ``shape[0] * shape[1] * shape[2]``
    ///         (default 100,000,000), checked before the voxel data block
    ///         is read.
    ///
    /// Raises:
    ///     ValueError: on parse failure (see `chematic_mol::CubeError`,
    ///     including a rejected multi-dataset file -- unsupported, see
    ///     module docs).
    #[staticmethod]
    #[pyo3(signature = (text, max_input_bytes=None, max_atoms=None, max_grid_points=None))]
    fn from_cube(
        text: &str,
        max_input_bytes: Option<usize>,
        max_atoms: Option<usize>,
        max_grid_points: Option<usize>,
    ) -> PyResult<Self> {
        let mut limits = chematic_mol::CubeParseLimits::default();
        if let Some(v) = max_input_bytes {
            limits.max_input_bytes = v;
        }
        if let Some(v) = max_atoms {
            limits.max_atoms = v;
        }
        if let Some(v) = max_grid_points {
            limits.max_grid_points = v;
        }
        let inner = chematic_mol::parse_cube_with_limits(text, &limits).map_err(cube_err)?;
        Ok(PyVolumetricGrid { inner })
    }

    /// Parse an OpenDX (APBS scalar-field subset) file (text, not a path).
    ///
    /// Args:
    ///     text: `.dx` file contents.
    ///     max_input_bytes: byte-size limit (default 1 GiB).
    ///     max_grid_points: cap on ``shape[0] * shape[1] * shape[2]``
    ///         (default 100,000,000).
    ///
    /// The returned grid always has ``units="angstrom"`` (OpenDX has no
    /// in-file unit tag; Ångström is APBS's own real-world convention --
    /// see `chematic_mol::opendx`'s module docs) and ``atoms=[]`` (OpenDX
    /// has no atom section at all).
    ///
    /// Raises:
    ///     ValueError: on parse failure.
    #[staticmethod]
    #[pyo3(signature = (text, max_input_bytes=None, max_grid_points=None))]
    fn from_opendx(
        text: &str,
        max_input_bytes: Option<usize>,
        max_grid_points: Option<usize>,
    ) -> PyResult<Self> {
        let mut limits = chematic_mol::OpenDxParseLimits::default();
        if let Some(v) = max_input_bytes {
            limits.max_input_bytes = v;
        }
        if let Some(v) = max_grid_points {
            limits.max_grid_points = v;
        }
        let inner = chematic_mol::parse_opendx_with_limits(text, &limits).map_err(opendx_err)?;
        Ok(PyVolumetricGrid { inner })
    }

    /// Write as a Gaussian Cube file. Always a single-dataset file.
    ///
    /// Raises:
    ///     ValueError: on an invalid grid (see `VolumetricGrid.validate`).
    fn to_cube(&self) -> PyResult<String> {
        chematic_mol::write_cube(&self.inner).map_err(cube_err)
    }

    /// Write as an OpenDX file. **Fails closed**: raises `ValueError` for a
    /// ``units="bohr"`` grid (OpenDX has no unit tag of its own; every
    /// real-world consumer assumes Ångström, so writing Bohr numbers as-is
    /// would have them silently reinterpreted ~1.89x wrong on next read --
    /// see `chematic_mol::opendx::write_opendx`'s Rust docs) and for a grid
    /// with any `atoms` (OpenDX has no atom section). Use
    /// :meth:`to_opendx_lossy` to explicitly opt into a Bohr->Ångström
    /// conversion instead of refusing.
    fn to_opendx(&self) -> PyResult<String> {
        chematic_mol::write_opendx(&self.inner).map_err(opendx_err)
    }

    /// Like :meth:`to_opendx`, but if ``units="bohr"``, explicitly converts
    /// `origin`/`axes` to Ångström instead of refusing (`values` -- the
    /// scalar-field samples themselves -- are never rescaled: only
    /// `origin`/`axes` are lengths). Still raises `ValueError` for a grid
    /// with any `atoms`, same as :meth:`to_opendx`.
    fn to_opendx_lossy(&self) -> PyResult<String> {
        chematic_mol::write_opendx_lossy(&self.inner).map_err(opendx_err)
    }

    #[getter]
    fn origin(&self) -> (f64, f64, f64) {
        let [x, y, z] = self.inner.origin;
        (x, y, z)
    }

    /// Per-axis step vectors, as a `(3, 3)` numpy array (rows are the 3
    /// axis vectors).
    #[getter]
    fn axes<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        axes_to_pyarray(py, self.inner.axes)
    }

    #[getter]
    fn shape(&self) -> (usize, usize, usize) {
        let [nx, ny, nz] = self.inner.shape;
        (nx, ny, nz)
    }

    /// Flat scalar-field samples as a 1-D numpy array of length
    /// ``shape[0] * shape[1] * shape[2]``. See the class docs for the index
    /// ordering. See :attr:`values_3d` for the same data reshaped to
    /// ``self.shape``.
    #[getter]
    fn values<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        Array1::from_vec(self.inner.values.clone()).into_pyarray(py)
    }

    /// :attr:`values` reshaped to a 3-D numpy array of shape ``self.shape``
    /// (``(nx, ny, nz)``), so ``values_3d[i, j, k] == get(i, j, k)``. A
    /// plain copy of :attr:`values`, not a zero-copy view.
    ///
    /// Raises:
    ///     ValueError: if ``shape``'s point count doesn't match
    ///     ``len(values)`` (only reachable for a hand-built, invalid grid;
    ///     see [`Self::new`]).
    #[getter]
    fn values_3d<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let [nx, ny, nz] = self.inner.shape;
        let arr = Array3::from_shape_vec((nx, ny, nz), self.inner.values.clone())
            .map_err(|e| PyValueError::new_err(format!("values does not match shape: {e}")))?;
        Ok(arr.into_pyarray(py))
    }

    #[getter]
    fn units(&self) -> &'static str {
        units_to_str(self.inner.units)
    }

    /// Atoms carried alongside the grid (Cube only; always empty for
    /// OpenDX), as ``[(element_symbol, charge, (x, y, z)), ...]``. `charge`
    /// is the effective nuclear charge (not a partial/formal charge -- see
    /// `chematic_mol::GridAtom::charge`'s Rust docs); position is in
    /// `self.units`.
    #[getter]
    fn atoms(&self) -> Vec<(String, f64, (f64, f64, f64))> {
        self.inner
            .atoms
            .iter()
            .map(|a| {
                (
                    a.element.symbol().to_string(),
                    a.charge,
                    (a.position[0], a.position[1], a.position[2]),
                )
            })
            .collect()
    }

    /// `shape[0] * shape[1] * shape[2]`.
    ///
    /// Raises:
    ///     ValueError: if that product overflows (only reachable for a
    ///     hand-built grid with a pathological `shape`).
    fn point_count(&self) -> PyResult<usize> {
        self.inner.point_count().map_err(grid_err)
    }

    /// Flat index into `values` for grid coordinates `(i, j, k)`, or `None`
    /// if out of bounds.
    fn checked_index(&self, i: usize, j: usize, k: usize) -> Option<usize> {
        self.inner.checked_index(i, j, k)
    }

    /// The value at grid coordinates `(i, j, k)`, or `None` if out of
    /// bounds.
    fn get(&self, i: usize, j: usize, k: usize) -> Option<f64> {
        self.inner.get(i, j, k)
    }

    /// Build a plain `Mol` + per-atom Cartesian coordinates (in
    /// `self.units`) from `atoms`, in file order. No bonds (neither Cube
    /// nor OpenDX carries a bond table).
    fn to_molecule(&self) -> (Mol, Vec<Vec<f64>>) {
        let (mol, coords) = self.inner.to_molecule();
        let py_coords = coords.iter().map(|&(x, y, z)| vec![x, y, z]).collect();
        (
            Mol {
                inner: Arc::new(mol),
                props: Default::default(),
            },
            py_coords,
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "VolumetricGrid(shape={:?}, units={:?}, atoms={})",
            self.inner.shape,
            units_to_str(self.inner.units),
            self.inner.atoms.len()
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyVolumetricGrid>()?;
    Ok(())
}
