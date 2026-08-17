//! Python bindings for `chematic-crystal`'s `Lattice`/`PeriodicStructure` —
//! periodic (crystal) structure representation, CIF/POSCAR read+write.
//!
//! # Design notes (see the PR body for the full write-up)
//!
//! - **Immutable wrappers.** `Lattice`/`PeriodicStructure`/`Site` mirror the
//!   Rust side's "immutable by convention" design: no setters, every
//!   transform (`wrap_into_cell`, `make_supercell`) returns a *new*
//!   `PeriodicStructure` rather than mutating in place.
//! - **Disorder is never collapsed.** `Site.species` is always the full
//!   `list[(element_symbol, occupancy)]` — never "the first" or "the
//!   highest-occupancy" element.
//! - **Every constructor routes through the Rust side's own fallible
//!   constructors** (`Lattice::from_matrix`/`from_parameters`/`cubic`/
//!   `orthorhombic`, `Occupancy::new`, `PeriodicSite::new`,
//!   `PeriodicStructure::new`) — there is no parallel validation logic in
//!   this module, so a Python caller can never construct an object that
//!   violates the Rust side's own invariants (occupancy sums, coordinate
//!   finiteness, degenerate lattices, ...).
//! - **Error mapping.** Every `CrystalError`/`CifPeriodicError`/`PoscarError`
//!   variant is, in Python-exception terms, a "this input value is invalid"
//!   error — none of them are I/O, missing-key, or index-range errors — so
//!   all three map uniformly to `ValueError`, matching this crate's existing
//!   convention of using builtin exceptions directly (see `rwmol.rs`/
//!   `io.rs`) rather than inventing a new exception hierarchy for a single
//!   error family (the one precedent for a *custom* exception,
//!   `pipeline_v2.rs`'s `PipelineV2Error`, exists because that call site
//!   needed rich structured diagnostics attached to the error itself; these
//!   error enums are already self-describing `Display` strings, so no
//!   custom exception adds anything here).
//! - **CIF symmetry status.** `PeriodicStructure.from_cif(text,
//!   expand_symmetry=True)` (the default) expands every symmetry operation
//!   literally written in the CIF's own operation-list tag into a full
//!   unit cell — never generated from a space-group name/number, never
//!   cross-checked against any space-group database (see
//!   `chematic_mol::CifSymmetryStatus`'s Rust docs for the full caveat on
//!   what "expanded" does and doesn't claim). `expand_symmetry=False`
//!   restores the pre-expansion behavior (asymmetric unit only).
//!   `PeriodicStructure.symmetry_status` is `None` for any structure not
//!   read via `from_cif` (direct construction, `from_poscar`, ...); for a
//!   CIF-sourced structure it is always `Some(CifSymmetryStatus(...))`,
//!   whose `is_complete_cell` property is `True` for `is_p1` or
//!   `is_expanded`, `False` for a genuinely unexpanded (or
//!   expansion-opted-out) asymmetric unit. `to_cif()` refuses (raises
//!   `ValueError`, via `chematic_mol::CifPeriodicResult::to_cif_checked` —
//!   the single Rust-side source of truth for this judgment, not
//!   re-implemented here) to re-emit a structure whose
//!   `is_complete_cell=False` — `write_cif_periodic_structure` always
//!   writes a literal `P 1` tag, so writing back an asymmetric-unit-only
//!   structure would falsely declare it complete (chematic's standing
//!   "never silently treat undeclared/unexpanded symmetry as P1" rule, PR
//!   #323). `symmetry_status` (and this `to_cif` guard) survive
//!   `wrap_into_cell()` (site count/order unchanged) and `make_supercell()`
//!   (site count changes, but the underlying "only the asymmetric unit was
//!   ever read" problem for an `UnexpandedSymmetry` structure does not go
//!   away, so the guard must not either — an already-complete-cell
//!   structure, `is_p1` or `is_expanded`, correctly stays writable after a
//!   supercell multiply too); they do not survive being passed through the
//!   direct `PeriodicStructure()` constructor (a fresh structure with no
//!   CIF provenance).
//! - **POSCAR extras out of scope as individual attributes.** `selective_dynamics`/
//!   `velocities`/`predictor_corrector` (and the file `comment`) are not
//!   exposed as separate Python-visible fields in this first version —
//!   but a bare `from_poscar(text)` immediately followed by `to_poscar()`
//!   (no intervening transform) does round-trip them faithfully, since
//!   they are kept internally and only cleared by a transform that could
//!   invalidate their per-site correspondence (`make_supercell`; not
//!   `wrap_into_cell`, which preserves site count and order). Verified by
//!   `test_crystal.py::test_poscar_roundtrip_preserves_extras`.
//! - **No `__eq__`, no custom serialization support**, matching this
//!   crate's existing baseline (only `__repr__`/`__str__` exist anywhere
//!   in `chematic-py` today).

use std::collections::BTreeMap;

use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use chematic_core::Element;
use chematic_crystal::poscar::{self, PoscarDocument, PoscarError, PredictorCorrector};
use chematic_crystal::{
    CartesianCoord, CrystalError, FractionalCoord, Lattice, Occupancy, PeriodicNeighbor,
    PeriodicSite, PeriodicStructure, SiteSpecies,
};
use chematic_mol::{
    CifPeriodicError, CifPeriodicParseOptions, CifPeriodicResult, CifSymmetryStatus,
    parse_cif_periodic_structure_with_options,
};

// ---------------------------------------------------------------------------
// Error mapping — see module docs for the "why ValueError uniformly" call.
// ---------------------------------------------------------------------------

fn crystal_err(e: CrystalError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn cif_err(e: CifPeriodicError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn poscar_err(e: PoscarError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

// ---------------------------------------------------------------------------
// Small numpy helpers (mirrors bulk.rs's IntoPyArray pattern — a fresh copy
// per call, not zero-copy; simple and safe, matching this crate's existing
// choice everywhere else it returns a numpy array).
// ---------------------------------------------------------------------------

fn matrix_to_pyarray<'py>(py: Python<'py>, m: [[f64; 3]; 3]) -> Bound<'py, PyArray2<f64>> {
    let flat: Vec<f64> = m.into_iter().flatten().collect();
    Array2::from_shape_vec((3, 3), flat)
        .expect("3x3 fixed-size matrix always reshapes cleanly")
        .into_pyarray(py)
}

fn points_to_pyarray<'py>(py: Python<'py>, rows: Vec<[f64; 3]>) -> Bound<'py, PyArray2<f64>> {
    let n = rows.len();
    let flat: Vec<f64> = rows.into_iter().flatten().collect();
    Array2::from_shape_vec((n, 3), flat)
        .expect("row width fixed at 3 always reshapes cleanly")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// CifSymmetryStatus
// ---------------------------------------------------------------------------

/// How a CIF's declared symmetry relates to a `PeriodicStructure.sites`.
///
/// - `is_p1=True`: no symmetry beyond P1 was declared (or nothing was
///   declared at all) -- the returned sites are already the complete cell.
/// - `is_expanded=True`: every symmetry operation *literally written* in
///   the CIF's own operation-list tag was applied. This is a faithfulness
///   claim about the CIF's own text, not a claim that the list is complete
///   or correct for the named/numbered space group -- see
///   `chematic_mol::CifSymmetryStatus`'s Rust docs for the full caveat.
/// - `is_complete_cell` (`is_p1 or is_expanded`): `True` iff `sites` is a
///   genuinely complete unit cell (safe to round-trip through `to_cif()`).
/// - `asymmetric_site_count`/`expanded_site_count`: only set (not `None`)
///   when `is_expanded` is `True`.
#[pyclass(name = "CifSymmetryStatus", from_py_object)]
#[derive(Clone)]
pub struct PyCifSymmetryStatus {
    #[pyo3(get)]
    is_p1: bool,
    #[pyo3(get)]
    is_expanded: bool,
    #[pyo3(get)]
    is_complete_cell: bool,
    #[pyo3(get)]
    space_group_name: Option<String>,
    #[pyo3(get)]
    operation_count: usize,
    #[pyo3(get)]
    asymmetric_site_count: Option<usize>,
    #[pyo3(get)]
    expanded_site_count: Option<usize>,
}

#[pymethods]
impl PyCifSymmetryStatus {
    fn __repr__(&self) -> String {
        format!(
            "CifSymmetryStatus(is_p1={}, is_expanded={}, space_group_name={:?}, \
             operation_count={}, asymmetric_site_count={:?}, expanded_site_count={:?})",
            self.is_p1,
            self.is_expanded,
            self.space_group_name,
            self.operation_count,
            self.asymmetric_site_count,
            self.expanded_site_count
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<CifSymmetryStatus> for PyCifSymmetryStatus {
    fn from(status: CifSymmetryStatus) -> Self {
        let is_complete_cell = status.is_complete_cell();
        match status {
            CifSymmetryStatus::P1 => PyCifSymmetryStatus {
                is_p1: true,
                is_expanded: false,
                is_complete_cell,
                space_group_name: None,
                operation_count: 0,
                asymmetric_site_count: None,
                expanded_site_count: None,
            },
            CifSymmetryStatus::ExpandedExplicitOperations {
                space_group_name,
                operation_count,
                asymmetric_site_count,
                expanded_site_count,
            } => PyCifSymmetryStatus {
                is_p1: false,
                is_expanded: true,
                is_complete_cell,
                space_group_name,
                operation_count,
                asymmetric_site_count: Some(asymmetric_site_count),
                expanded_site_count: Some(expanded_site_count),
            },
            CifSymmetryStatus::UnexpandedSymmetry {
                space_group_name,
                operation_count,
            } => PyCifSymmetryStatus {
                is_p1: false,
                is_expanded: false,
                is_complete_cell,
                space_group_name,
                operation_count,
                asymmetric_site_count: None,
                expanded_site_count: None,
            },
        }
    }
}

// ---------------------------------------------------------------------------
// PeriodicNeighbor
// ---------------------------------------------------------------------------

/// One periodic neighbor relationship — see the Python-facing
/// `PeriodicStructure.neighbors` method.
#[pyclass(name = "PeriodicNeighbor", from_py_object)]
#[derive(Clone)]
pub struct PyPeriodicNeighbor {
    #[pyo3(get)]
    center_index: usize,
    #[pyo3(get)]
    neighbor_index: usize,
    #[pyo3(get)]
    image: (i32, i32, i32),
    #[pyo3(get)]
    displacement: (f64, f64, f64),
    #[pyo3(get)]
    distance: f64,
}

#[pymethods]
impl PyPeriodicNeighbor {
    fn __repr__(&self) -> String {
        format!(
            "PeriodicNeighbor(center_index={}, neighbor_index={}, image={:?}, distance={:.4})",
            self.center_index, self.neighbor_index, self.image, self.distance
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<PeriodicNeighbor> for PyPeriodicNeighbor {
    fn from(n: PeriodicNeighbor) -> Self {
        PyPeriodicNeighbor {
            center_index: n.center_index,
            neighbor_index: n.neighbor_index,
            image: (n.image[0], n.image[1], n.image[2]),
            displacement: (n.displacement[0], n.displacement[1], n.displacement[2]),
            distance: n.distance,
        }
    }
}

// ---------------------------------------------------------------------------
// Site
// ---------------------------------------------------------------------------

/// A periodic site: one or more `(element_symbol, occupancy)` species (more
/// than one models positional/substitutional disorder), a fractional
/// position, and an optional label.
#[pyclass(name = "Site", from_py_object)]
#[derive(Clone)]
pub struct PySite {
    inner: PeriodicSite,
}

#[pymethods]
impl PySite {
    /// `species` is a list of `(element_symbol, occupancy)` pairs — more
    /// than one entry models a disordered site. Routes through
    /// `Element::from_symbol`, `Occupancy::new`, and `PeriodicSite::new`,
    /// so an unknown symbol, a non-finite/negative occupancy, an
    /// occupancy sum over `1.0 + tolerance`, an empty species list, or a
    /// non-finite fractional position all raise `ValueError` here rather
    /// than constructing an invalid `Site`.
    #[new]
    #[pyo3(signature = (species, fractional, label=None))]
    fn new(
        species: Vec<(String, f64)>,
        fractional: (f64, f64, f64),
        label: Option<String>,
    ) -> PyResult<Self> {
        let species = species
            .into_iter()
            .map(|(symbol, occupancy)| -> PyResult<SiteSpecies> {
                let element = Element::from_symbol(&symbol).ok_or_else(|| {
                    PyValueError::new_err(format!("unknown element symbol: {symbol:?}"))
                })?;
                let occupancy = Occupancy::new(occupancy).map_err(crystal_err)?;
                Ok(SiteSpecies { element, occupancy })
            })
            .collect::<PyResult<Vec<_>>>()?;
        let frac = FractionalCoord::new([fractional.0, fractional.1, fractional.2]);
        let inner = PeriodicSite::new(species, frac, label).map_err(crystal_err)?;
        Ok(PySite { inner })
    }

    /// `[(element_symbol, occupancy), ...]` — never collapsed to a single
    /// entry, even for a fully-ordered (single-species) site.
    #[getter]
    fn species(&self) -> Vec<(String, f64)> {
        self.inner
            .species
            .iter()
            .map(|s| (s.element.symbol().to_string(), s.occupancy.value()))
            .collect()
    }

    #[getter]
    fn fractional(&self) -> (f64, f64, f64) {
        let [x, y, z] = self.inner.fractional.0;
        (x, y, z)
    }

    #[getter]
    fn label(&self) -> Option<String> {
        self.inner.label.clone()
    }

    fn __repr__(&self) -> String {
        let species = self
            .species()
            .iter()
            .map(|(symbol, occ)| format!("{symbol}:{occ:.4}"))
            .collect::<Vec<_>>()
            .join(", ");
        let (x, y, z) = self.fractional();
        format!(
            "Site([{species}], fractional=({x:.5}, {y:.5}, {z:.5}), label={:?})",
            self.label()
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

// ---------------------------------------------------------------------------
// Lattice
// ---------------------------------------------------------------------------

/// A validated 3x3 lattice matrix (rows = lattice vectors a, b, c). See
/// `chematic_crystal::Lattice`'s own docs for the row/matrix convention.
#[pyclass(name = "Lattice", from_py_object)]
#[derive(Clone)]
pub struct PyLattice {
    inner: Lattice,
}

#[pymethods]
impl PyLattice {
    /// Build directly from a 3x3 matrix (rows = a, b, c).
    #[staticmethod]
    fn from_matrix(matrix: [[f64; 3]; 3]) -> PyResult<Self> {
        Ok(PyLattice {
            inner: Lattice::from_matrix(matrix).map_err(crystal_err)?,
        })
    }

    /// Build from crystallographic parameters (lengths in Angstrom, angles
    /// in degrees).
    #[staticmethod]
    fn from_parameters(
        a: f64,
        b: f64,
        c: f64,
        alpha: f64,
        beta: f64,
        gamma: f64,
    ) -> PyResult<Self> {
        Ok(PyLattice {
            inner: Lattice::from_parameters(a, b, c, alpha, beta, gamma).map_err(crystal_err)?,
        })
    }

    /// A cubic cell with edge length `a` (all angles 90 degrees).
    #[staticmethod]
    fn cubic(a: f64) -> PyResult<Self> {
        Ok(PyLattice {
            inner: Lattice::cubic(a).map_err(crystal_err)?,
        })
    }

    /// An orthorhombic cell with edge lengths `a, b, c` (all angles 90
    /// degrees).
    #[staticmethod]
    fn orthorhombic(a: f64, b: f64, c: f64) -> PyResult<Self> {
        Ok(PyLattice {
            inner: Lattice::orthorhombic(a, b, c).map_err(crystal_err)?,
        })
    }

    /// The lattice matrix (rows = a, b, c), as a `(3, 3)` numpy array.
    #[getter]
    fn matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        matrix_to_pyarray(py, self.inner.matrix())
    }

    /// The cached matrix inverse, as a `(3, 3)` numpy array.
    #[getter]
    fn inverse_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        matrix_to_pyarray(py, self.inner.inverse_matrix())
    }

    /// Reciprocal-lattice matrix (crystallographic, no `2*pi`), as a `(3, 3)`
    /// numpy array.
    #[getter]
    fn reciprocal_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        matrix_to_pyarray(py, self.inner.reciprocal_matrix())
    }

    #[getter]
    fn volume(&self) -> f64 {
        self.inner.volume()
    }

    /// `(a, b, c)` lengths in Angstrom.
    #[getter]
    fn lengths(&self) -> (f64, f64, f64) {
        let [a, b, c] = self.inner.lengths();
        (a, b, c)
    }

    /// `(alpha, beta, gamma)` angles in degrees.
    #[getter]
    fn angles_degrees(&self) -> (f64, f64, f64) {
        let [a, b, c] = self.inner.angles_degrees();
        (a, b, c)
    }

    /// Convert one fractional point to Cartesian (Angstrom).
    fn frac_to_cart(&self, point: (f64, f64, f64)) -> (f64, f64, f64) {
        let c = self
            .inner
            .frac_to_cart(FractionalCoord::new([point.0, point.1, point.2]));
        (c.0[0], c.0[1], c.0[2])
    }

    /// Convert one Cartesian point (Angstrom) to fractional.
    fn cart_to_frac(&self, point: (f64, f64, f64)) -> (f64, f64, f64) {
        let f = self
            .inner
            .cart_to_frac(CartesianCoord::new([point.0, point.1, point.2]));
        (f.0[0], f.0[1], f.0[2])
    }

    fn __repr__(&self) -> String {
        let (a, b, c) = self.lengths();
        let (alpha, beta, gamma) = self.angles_degrees();
        format!(
            "Lattice(a={a:.4}, b={b:.4}, c={c:.4}, alpha={alpha:.2}, beta={beta:.2}, gamma={gamma:.2})"
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

// ---------------------------------------------------------------------------
// PeriodicStructure
// ---------------------------------------------------------------------------

/// POSCAR-only fields this binding doesn't expose as individual Python
/// attributes (see module docs). Kept only long enough to make a bare
/// `from_poscar()` -> `to_poscar()` round trip faithful.
#[derive(Clone)]
struct PoscarExtra {
    comment: String,
    selective_dynamics: Option<Vec<[bool; 3]>>,
    velocities: Option<Vec<[f64; 3]>>,
    predictor_corrector: Option<PredictorCorrector>,
}

/// A periodic structure: a [`PyLattice`] plus an ordered list of
/// [`PySite`]s. Immutable by convention (matching the Rust side) —
/// `wrap_into_cell`/`make_supercell` return a new `PeriodicStructure`
/// rather than mutating in place.
#[pyclass(name = "PeriodicStructure")]
pub struct PyPeriodicStructure {
    inner: PeriodicStructure,
    /// The raw Rust status, not the Python-facing summary struct — kept
    /// this way so `to_cif()` can delegate to
    /// `chematic_mol::CifPeriodicResult::to_cif_checked` (the single source
    /// of truth for the write-safety judgment) instead of re-implementing
    /// it against a lossy Python-side copy; `symmetry_status` derives
    /// `PyCifSymmetryStatus` from this on demand.
    cif_symmetry_status: Option<CifSymmetryStatus>,
    poscar_extra: Option<PoscarExtra>,
}

impl PyPeriodicStructure {
    fn bare(inner: PeriodicStructure) -> Self {
        PyPeriodicStructure {
            inner,
            cif_symmetry_status: None,
            poscar_extra: None,
        }
    }
}

#[pymethods]
impl PyPeriodicStructure {
    /// Construct directly from a lattice and a list of sites. Routes
    /// through `PeriodicStructure::new`, so an invalid site (caught
    /// already by `Site.__new__`) can't reach here, but this also re-runs
    /// that same validation.
    #[new]
    fn new(lattice: PyLattice, sites: Vec<PySite>) -> PyResult<Self> {
        let sites: Vec<PeriodicSite> = sites.into_iter().map(|s| s.inner).collect();
        let inner = PeriodicStructure::new(lattice.inner, sites).map_err(crystal_err)?;
        Ok(PyPeriodicStructure::bare(inner))
    }

    /// Parse a CIF file (text, not a path). `expand_symmetry` (default
    /// `True`) expands every symmetry operation literally written in the
    /// CIF's own operation-list tag into a full unit cell — see
    /// `symmetry_status` on the returned structure for whether/how that
    /// happened (never generated from a space-group name/number, and never
    /// cross-checked against any space-group database — see
    /// `CifSymmetryStatus`'s docs). `expand_symmetry=False` restores the
    /// pre-expansion behavior: only the asymmetric unit as literally
    /// listed is returned, and `symmetry_status.is_expanded` is always
    /// `False`.
    #[staticmethod]
    #[pyo3(signature = (text, expand_symmetry=true))]
    fn from_cif(text: &str, expand_symmetry: bool) -> PyResult<Self> {
        let options = CifPeriodicParseOptions {
            expand_explicit_symmetry: expand_symmetry,
        };
        let result = parse_cif_periodic_structure_with_options(text, options).map_err(cif_err)?;
        Ok(PyPeriodicStructure {
            inner: result.structure,
            cif_symmetry_status: Some(result.symmetry),
            poscar_extra: None,
        })
    }

    /// Parse a POSCAR/CONTCAR file (text, not a path). Both formats are
    /// byte-identical grammars (CONTCAR is just VASP's *output*
    /// convention for the same layout), so there is no separate
    /// `from_contcar`.
    #[staticmethod]
    fn from_poscar(text: &str) -> PyResult<Self> {
        let doc = poscar::parse_poscar(text).map_err(poscar_err)?;
        Ok(PyPeriodicStructure {
            inner: doc.structure,
            cif_symmetry_status: None,
            poscar_extra: Some(PoscarExtra {
                comment: doc.comment,
                selective_dynamics: doc.selective_dynamics,
                velocities: doc.velocities,
                predictor_corrector: doc.predictor_corrector,
            }),
        })
    }

    #[getter]
    fn lattice(&self) -> PyLattice {
        PyLattice {
            inner: self.inner.lattice().clone(),
        }
    }

    #[getter]
    fn sites(&self) -> Vec<PySite> {
        self.inner
            .sites()
            .iter()
            .cloned()
            .map(|inner| PySite { inner })
            .collect()
    }

    fn site_count(&self) -> usize {
        self.inner.site_count()
    }

    /// Cartesian position of every site, in site order, as an `(N, 3)`
    /// numpy array.
    fn cartesian_positions<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let rows = self
            .inner
            .cartesian_positions()
            .into_iter()
            .map(|c| c.0)
            .collect();
        points_to_pyarray(py, rows)
    }

    /// Fractional position of every site, in site order, as an `(N, 3)`
    /// numpy array.
    fn fractional_positions<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let rows = self
            .inner
            .fractional_positions()
            .into_iter()
            .map(|c| c.0)
            .collect();
        points_to_pyarray(py, rows)
    }

    /// Every periodic neighbor pair within `cutoff` Angstrom (inclusive).
    /// Raises `ValueError` for a non-finite/non-positive cutoff, or for a
    /// cutoff so large relative to the cell that the exact search would
    /// examine an unreasonable number of periodic images (fails closed
    /// rather than hanging — see `chematic_crystal::neighbor::MAX_NEIGHBOR_IMAGE_CANDIDATES`).
    fn neighbors(&self, cutoff: f64) -> PyResult<Vec<PyPeriodicNeighbor>> {
        let neighbors = self.inner.neighbors_within(cutoff).map_err(crystal_err)?;
        Ok(neighbors
            .into_iter()
            .map(PyPeriodicNeighbor::from)
            .collect())
    }

    /// Build a diagonal `(nx, ny, nz)` supercell as a new `PeriodicStructure`.
    fn make_supercell(&self, mult: (u32, u32, u32)) -> PyResult<Self> {
        let inner = self
            .inner
            .make_supercell([mult.0, mult.1, mult.2])
            .map_err(crystal_err)?;
        Ok(PyPeriodicStructure {
            inner,
            // The asymmetric-unit-vs-full-cell problem an `UnexpandedSymmetry`
            // status flags doesn't go away under a supercell expansion, so
            // the status (and the `to_cif` guard it drives) must survive
            // too. Conversely, an already-`ExpandedExplicitOperations` or
            // `P1` structure is still a genuinely complete cell after a
            // supercell multiply (that transform is defined in terms of
            // whole unit cells) -- `is_complete_cell` correctly stays
            // `True` for those.
            cif_symmetry_status: self.cif_symmetry_status.clone(),
            // Site count changes, so `poscar_extra`'s per-site vectors
            // (selective_dynamics/velocities) would no longer correspond
            // positionally -- must not carry those over.
            poscar_extra: None,
        })
    }

    /// A new `PeriodicStructure` with every site's fractional position
    /// reduced into `[0, 1)`. Site count and order are unchanged, so
    /// `symmetry_status` and the POSCAR round-trip extras both survive
    /// this transform.
    fn wrap_into_cell(&self) -> Self {
        PyPeriodicStructure {
            inner: self.inner.wrapped(),
            cif_symmetry_status: self.cif_symmetry_status.clone(),
            poscar_extra: self.poscar_extra.clone(),
        }
    }

    /// `None` unless this structure was produced by `from_cif`, in which
    /// case it reports how (if at all) the file's declared symmetry was
    /// expanded (see `CifSymmetryStatus`'s docs).
    #[getter]
    fn symmetry_status(&self) -> Option<PyCifSymmetryStatus> {
        self.cif_symmetry_status
            .clone()
            .map(PyCifSymmetryStatus::from)
    }

    /// Hill-order (C, H, then alphabetical) formula string summed over
    /// every site's occupancy-weighted species — the raw, *unreduced* cell
    /// content (e.g. a 2x2x2 supercell of NaCl reports `"Cl8Na8"`, not
    /// `"ClNa"`). A disordered site contributes its species' occupancies
    /// as fractional counts (e.g. `"Fe0.6Ni0.4"`).
    #[getter]
    fn formula(&self) -> String {
        format_crystal_formula(&self.inner)
    }

    /// Write as CIF text. Raises `ValueError` if this structure's
    /// `symmetry_status.is_complete_cell` is `False` — see
    /// `CifSymmetryStatus`'s docs for why (writing would falsely declare a
    /// complete cell for what is only an asymmetric unit). Delegates to
    /// `chematic_mol::CifPeriodicResult::to_cif_checked` (the single
    /// Rust-side source of truth for this judgment) rather than
    /// re-implementing the check here; a structure with no CIF provenance
    /// (`symmetry_status is None`) is always writable, same as before.
    fn to_cif(&self) -> PyResult<String> {
        let symmetry = self
            .cif_symmetry_status
            .clone()
            .unwrap_or(CifSymmetryStatus::P1);
        let result = CifPeriodicResult {
            structure: self.inner.clone(),
            symmetry,
        };
        result.to_cif_checked().map_err(cif_err)
    }

    /// Write as POSCAR text (VASP 5 format: explicit species-name line,
    /// scale factor `1.0`, pre-scaled lattice vectors, `Direct`
    /// (fractional) coordinates). Raises `ValueError` if any site has more
    /// than one species or occupancy other than `1.0` -- POSCAR has no
    /// disorder/partial-occupancy representation.
    fn to_poscar(&self) -> PyResult<String> {
        let (comment, selective_dynamics, velocities, predictor_corrector) =
            match &self.poscar_extra {
                Some(extra) => (
                    extra.comment.clone(),
                    extra.selective_dynamics.clone(),
                    extra.velocities.clone(),
                    extra.predictor_corrector.clone(),
                ),
                None => ("chematic".to_string(), None, None, None),
            };
        let doc = PoscarDocument {
            structure: self.inner.clone(),
            comment,
            // Read-fidelity-only field; write_poscar re-derives species
            // grouping from `structure.sites()` and never consults this.
            species_order: Vec::new(),
            selective_dynamics,
            velocities,
            predictor_corrector,
        };
        poscar::write_poscar(&doc).map_err(poscar_err)
    }

    fn __repr__(&self) -> String {
        format!(
            "PeriodicStructure(sites={}, formula='{}')",
            self.inner.site_count(),
            self.formula()
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

/// Hill order (C first, H second, then alphabetical), matching
/// `chematic_core::Molecule::formula`'s own convention -- see the
/// `formula` getter's docs for the occupancy-weighted / unreduced
/// semantics specific to a periodic structure.
fn format_crystal_formula(structure: &PeriodicStructure) -> String {
    let mut counts: BTreeMap<&'static str, f64> = BTreeMap::new();
    for site in structure.sites() {
        for species in &site.species {
            *counts.entry(species.element.symbol()).or_insert(0.0) += species.occupancy.value();
        }
    }
    let mut out = String::new();
    if let Some(c) = counts.remove("C") {
        push_formula_count("C", c, &mut out);
    }
    if let Some(h) = counts.remove("H") {
        push_formula_count("H", h, &mut out);
    }
    for (symbol, count) in &counts {
        push_formula_count(symbol, *count, &mut out);
    }
    out
}

/// Integer counts (within float tolerance) print as `Sym` (count `1`
/// omitted) or `SymN`; non-integer counts (disordered/partial-occupancy
/// sites) always print the number, trimmed of trailing zeros.
fn push_formula_count(symbol: &str, count: f64, out: &mut String) {
    out.push_str(symbol);
    if (count - count.round()).abs() < 1e-6 {
        let n = count.round() as i64;
        if n != 1 {
            out.push_str(&n.to_string());
        }
    } else {
        let s = format!("{count:.4}");
        out.push_str(s.trim_end_matches('0').trim_end_matches('.'));
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyLattice>()?;
    m.add_class::<PySite>()?;
    m.add_class::<PyPeriodicStructure>()?;
    m.add_class::<PyPeriodicNeighbor>()?;
    m.add_class::<PyCifSymmetryStatus>()?;
    Ok(())
}
