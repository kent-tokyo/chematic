//! Python bindings for `chematic-mol`'s LAMMPS data (`read_data`/`write_data`
//! format) and dump/trajectory (`dump` command's default text style)
//! support.
//!
//! Both are **standalone document types**, not integrated with
//! `chematic_core::Molecule` -- LAMMPS bonds/atom-types are raw index
//! topology, not chemically perceived bonds, and an MD atom under some
//! `atom_style`s is not necessarily even a chemical element (see
//! `chematic_mol::lammps_data`'s module docs). Callers needing a `Mol`
//! should build one themselves from the returned atom/bond data.
//!
//! `LammpsData` (a single document: header + typed sections) is bound as a
//! plain function + dict pair (`parse_lammps_data`/`write_lammps_data`),
//! matching `formats.rs`'s `parse_cif` precedent -- its only "behavior"
//! beyond a bag of parsed sections is `LammpsData::count(label)`, a lookup
//! Python's own `dict.get` on the returned `"counts"` dict already covers.
//!
//! `LammpsDumpFrame` (one frame of a trajectory) is instead a small
//! pyclass, [`PyLammpsDumpFrame`], because it has real computed behavior
//! worth exposing as methods -- `column()`/`cartesian_positions()`, the
//! latter doing a real (triclinic-aware) coordinate transform, not just a
//! field lookup -- following `crystal.rs`'s precedent for "Rust type with
//! real methods -> dedicated pyclass" rather than `formats.rs`'s
//! plain-dict precedent.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2};

use crate::formats::{dict_get, dict_get_opt, dict_get_required};

fn atom_style_from_str(s: &str) -> chematic_mol::LammpsAtomStyle {
    match s {
        "atomic" => chematic_mol::LammpsAtomStyle::Atomic,
        "charge" => chematic_mol::LammpsAtomStyle::Charge,
        "molecular" => chematic_mol::LammpsAtomStyle::Molecular,
        "full" => chematic_mol::LammpsAtomStyle::Full,
        // Not pre-validated here -- `parse_lammps_data` itself rejects any
        // `Other` style with a typed `UnsupportedAtomStyle` error; letting
        // that be the single source of truth for "which styles are
        // supported" avoids two places that could drift apart.
        other => chematic_mol::LammpsAtomStyle::Other(other.to_string()),
    }
}

fn atom_style_to_string(s: &chematic_mol::LammpsAtomStyle) -> String {
    match s {
        chematic_mol::LammpsAtomStyle::Atomic => "atomic".to_string(),
        chematic_mol::LammpsAtomStyle::Charge => "charge".to_string(),
        chematic_mol::LammpsAtomStyle::Molecular => "molecular".to_string(),
        chematic_mol::LammpsAtomStyle::Full => "full".to_string(),
        chematic_mol::LammpsAtomStyle::Other(s) => s.clone(),
    }
}

fn lammps_box_to_pydict<'py>(
    py: Python<'py>,
    b: &chematic_mol::LammpsBox,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("lo", b.lo)?;
    d.set_item("hi", b.hi)?;
    d.set_item("tilt", b.tilt)?;
    Ok(d)
}

fn pydict_to_lammps_box(d: &Bound<PyDict>) -> PyResult<chematic_mol::LammpsBox> {
    let b = chematic_mol::LammpsBox {
        lo: dict_get_required(d, "lo")?,
        hi: dict_get_required(d, "hi")?,
        tilt: dict_get_opt(d, "tilt")?,
    };
    b.validate().map_err(PyValueError::new_err)?;
    Ok(b)
}

// ---------------------------------------------------------------------------
// LammpsData (plain function + dict)
// ---------------------------------------------------------------------------

fn lammps_data_to_pydict<'py>(
    py: Python<'py>,
    data: &chematic_mol::LammpsData,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);

    let counts = PyDict::new(py);
    for (label, count) in &data.counts {
        counts.set_item(label, count)?;
    }
    d.set_item("counts", counts)?;
    d.set_item("atom_style", atom_style_to_string(&data.atom_style))?;
    d.set_item("box", lammps_box_to_pydict(py, &data.simulation_box)?)?;

    let masses: Vec<Bound<PyDict>> = data
        .masses
        .iter()
        .map(|m| {
            let md = PyDict::new(py);
            md.set_item("atom_type", m.atom_type)?;
            md.set_item("mass", m.mass)?;
            Ok::<_, PyErr>(md)
        })
        .collect::<PyResult<_>>()?;
    d.set_item("masses", masses)?;

    let atoms: Vec<Bound<PyDict>> = data
        .atoms
        .iter()
        .map(|a| {
            let ad = PyDict::new(py);
            ad.set_item("id", a.id)?;
            ad.set_item("molecule_id", a.molecule_id)?;
            ad.set_item("atom_type", a.atom_type)?;
            ad.set_item("charge", a.charge)?;
            ad.set_item("x", a.x)?;
            ad.set_item("y", a.y)?;
            ad.set_item("z", a.z)?;
            ad.set_item("image", a.image)?;
            Ok::<_, PyErr>(ad)
        })
        .collect::<PyResult<_>>()?;
    d.set_item("atoms", atoms)?;

    let velocities: Vec<Bound<PyDict>> = data
        .velocities
        .iter()
        .map(|v| {
            let vd = PyDict::new(py);
            vd.set_item("atom_id", v.atom_id)?;
            vd.set_item("vx", v.vx)?;
            vd.set_item("vy", v.vy)?;
            vd.set_item("vz", v.vz)?;
            Ok::<_, PyErr>(vd)
        })
        .collect::<PyResult<_>>()?;
    d.set_item("velocities", velocities)?;

    let bonds: Vec<Bound<PyDict>> = data
        .bonds
        .iter()
        .map(|b| {
            let bd = PyDict::new(py);
            bd.set_item("id", b.id)?;
            bd.set_item("bond_type", b.bond_type)?;
            bd.set_item("atom1", b.atom1)?;
            bd.set_item("atom2", b.atom2)?;
            Ok::<_, PyErr>(bd)
        })
        .collect::<PyResult<_>>()?;
    d.set_item("bonds", bonds)?;

    d.set_item("unparsed_sections", data.unparsed_sections.clone())?;
    Ok(d)
}

fn pydict_to_lammps_data(d: &Bound<PyDict>) -> PyResult<chematic_mol::LammpsData> {
    let counts: Vec<(String, i64)> = match d.get_item("counts")? {
        Some(v) if !v.is_none() => {
            let cd = v
                .cast::<PyDict>()
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            cd.iter()
                .map(|(k, v)| Ok::<_, PyErr>((k.extract::<String>()?, v.extract::<i64>()?)))
                .collect::<PyResult<_>>()?
        }
        _ => Vec::new(),
    };
    let atom_style_str: String = dict_get(d, "atom_style", "atomic".to_string())?;
    let box_dict: Bound<PyDict> = dict_get_required(d, "box")?;

    let masses_list: Vec<Bound<PyDict>> = dict_get(d, "masses", Vec::new())?;
    let masses = masses_list
        .iter()
        .map(|md| {
            Ok::<_, PyErr>(chematic_mol::LammpsMass {
                atom_type: dict_get_required(md, "atom_type")?,
                mass: dict_get_required(md, "mass")?,
            })
        })
        .collect::<PyResult<_>>()?;

    let atoms_list: Vec<Bound<PyDict>> = dict_get(d, "atoms", Vec::new())?;
    let atoms = atoms_list
        .iter()
        .map(|ad| {
            Ok::<_, PyErr>(chematic_mol::LammpsAtom {
                id: dict_get_required(ad, "id")?,
                molecule_id: dict_get_opt(ad, "molecule_id")?,
                atom_type: dict_get_required(ad, "atom_type")?,
                charge: dict_get_opt(ad, "charge")?,
                x: dict_get_required(ad, "x")?,
                y: dict_get_required(ad, "y")?,
                z: dict_get_required(ad, "z")?,
                image: dict_get_opt(ad, "image")?,
            })
        })
        .collect::<PyResult<_>>()?;

    let velocities_list: Vec<Bound<PyDict>> = dict_get(d, "velocities", Vec::new())?;
    let velocities = velocities_list
        .iter()
        .map(|vd| {
            Ok::<_, PyErr>(chematic_mol::LammpsVelocity {
                atom_id: dict_get_required(vd, "atom_id")?,
                vx: dict_get_required(vd, "vx")?,
                vy: dict_get_required(vd, "vy")?,
                vz: dict_get_required(vd, "vz")?,
            })
        })
        .collect::<PyResult<_>>()?;

    let bonds_list: Vec<Bound<PyDict>> = dict_get(d, "bonds", Vec::new())?;
    let bonds = bonds_list
        .iter()
        .map(|bd| {
            Ok::<_, PyErr>(chematic_mol::LammpsBond {
                id: dict_get_required(bd, "id")?,
                bond_type: dict_get_required(bd, "bond_type")?,
                atom1: dict_get_required(bd, "atom1")?,
                atom2: dict_get_required(bd, "atom2")?,
            })
        })
        .collect::<PyResult<_>>()?;

    let unparsed_sections: Vec<(String, String)> = dict_get(d, "unparsed_sections", Vec::new())?;

    Ok(chematic_mol::LammpsData {
        counts,
        atom_style: atom_style_from_str(&atom_style_str),
        simulation_box: pydict_to_lammps_box(&box_dict)?,
        masses,
        atoms,
        velocities,
        bonds,
        unparsed_sections,
    })
}

/// Parse a LAMMPS data file (`read_data` command format).
///
/// Args:
///     text: data file contents.
///     atom_style: one of ``"atomic"``, ``"charge"``, ``"molecular"``,
///         ``"full"`` -- **required**, since it cannot be recovered from
///         the file alone (e.g. `atom_style charge` and `atom_style
///         molecular` rows are both 6 fields, genuinely ambiguous by shape
///         alone -- see `chematic_mol::lammps_data`'s module docs). Any
///         other string is rejected with ``ValueError`` (matching
///         `LammpsAtomStyle::Other`'s Rust-side rejection).
///
/// Returns:
///     dict: ``{"counts": dict[str, int], "atom_style": str, "box": dict,
///     "masses": list[dict], "atoms": list[dict], "velocities": list[dict],
///     "bonds": list[dict], "unparsed_sections": list[tuple[str, str]]}``.
///
///     ``"box"`` is ``{"lo": (x, y, z), "hi": (x, y, z), "tilt":
///     (xy, xz, yz) | None}``.
///
///     ``"unparsed_sections"`` carries every section this module doesn't
///     semantically parse (``Angles``, ``Dihedrals``, `` *Coeffs``, ...)
///     verbatim, as ``(section_name, raw_row_text)`` pairs -- never
///     silently dropped.
///
/// Raises:
///     ValueError: on parse failure (see `chematic_mol::LammpsDataError`,
///     including LAMMPS's string-keyed "Type Labels" extension, which is
///     unsupported and fails closed).
///
/// Example::
///
///     data = chematic.parse_lammps_data(open("in.data").read(), "full")
///     print(data["counts"]["atoms"], len(data["atoms"]))
#[pyfunction]
fn parse_lammps_data<'py>(
    py: Python<'py>,
    text: &str,
    atom_style: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let style = atom_style_from_str(atom_style);
    let data = chematic_mol::parse_lammps_data(text, style)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    lammps_data_to_pydict(py, &data)
}

/// Write a LAMMPS data file from a dict, same shape as
/// :func:`parse_lammps_data`'s return value.
///
/// Raises:
///     ValueError: on a missing required key or an invalid box (see
///     `LammpsBox.validate`).
#[pyfunction]
fn write_lammps_data(data: &Bound<PyDict>) -> PyResult<String> {
    let data = pydict_to_lammps_data(data)?;
    Ok(chematic_mol::write_lammps_data(&data))
}

/// Convert a dump file's `ITEM: BOX BOUNDS` values into the true
/// simulation box (identity for an orthogonal box; undoes the triclinic
/// bounding-box shift otherwise -- see
/// `chematic_mol::lammps_dump`'s module docs).
#[pyfunction]
fn box_bounds_to_true<'py>(
    py: Python<'py>,
    bound_lo: (f64, f64, f64),
    bound_hi: (f64, f64, f64),
    tilt: Option<(f64, f64, f64)>,
) -> PyResult<Bound<'py, PyDict>> {
    let b = chematic_mol::box_bounds_to_true(
        [bound_lo.0, bound_lo.1, bound_lo.2],
        [bound_hi.0, bound_hi.1, bound_hi.2],
        tilt.map(|(xy, xz, yz)| [xy, xz, yz]),
    );
    lammps_box_to_pydict(py, &b)
}

/// Inverse of :func:`box_bounds_to_true`: convert a true box dict (``{"lo":
/// ..., "hi": ..., "tilt": ...}``) into the `xlo_bound`/`xhi_bound`/...
/// values a dump file's `ITEM: BOX BOUNDS` section would show.
///
/// Returns:
///     tuple: ``((xlo_bound, ylo_bound, zlo_bound), (xhi_bound, yhi_bound,
///     zhi_bound))``.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn true_to_box_bounds(box_dict: &Bound<PyDict>) -> PyResult<((f64, f64, f64), (f64, f64, f64))> {
    let b = pydict_to_lammps_box(box_dict)?;
    let (lo, hi) = chematic_mol::true_to_box_bounds(&b);
    Ok(((lo[0], lo[1], lo[2]), (hi[0], hi[1], hi[2])))
}

// ---------------------------------------------------------------------------
// LammpsDumpFrame (pyclass -- see module docs for why)
// ---------------------------------------------------------------------------

/// One frame of a LAMMPS dump/trajectory file.
#[pyclass(name = "LammpsDumpFrame", from_py_object)]
#[derive(Clone)]
pub struct PyLammpsDumpFrame {
    inner: chematic_mol::LammpsDumpFrame,
}

#[pymethods]
impl PyLammpsDumpFrame {
    /// `num_atoms` defaults to `len(rows)` when omitted.
    #[new]
    #[pyo3(signature = (timestep, box_bounds, column_names, rows, boundary_flags=("pp".to_string(), "pp".to_string(), "pp".to_string()), num_atoms=None))]
    fn new(
        timestep: i64,
        box_bounds: &Bound<PyDict>,
        column_names: Vec<String>,
        rows: Vec<Vec<f64>>,
        boundary_flags: (String, String, String),
        num_atoms: Option<usize>,
    ) -> PyResult<Self> {
        let box_bounds = pydict_to_lammps_box(box_bounds)?;
        let num_atoms = num_atoms.unwrap_or(rows.len());
        Ok(PyLammpsDumpFrame {
            inner: chematic_mol::LammpsDumpFrame {
                timestep,
                num_atoms,
                box_bounds,
                boundary_flags: [boundary_flags.0, boundary_flags.1, boundary_flags.2],
                column_names,
                rows,
            },
        })
    }

    #[getter]
    fn timestep(&self) -> i64 {
        self.inner.timestep
    }

    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms
    }

    #[getter]
    fn box_bounds<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        lammps_box_to_pydict(py, &self.inner.box_bounds)
    }

    /// The 2-character boundary-condition flag pair for each axis (x, y,
    /// z), e.g. `("pp", "pp", "pp")` -- taken verbatim from the file, not
    /// interpreted.
    #[getter]
    fn boundary_flags(&self) -> (String, String, String) {
        let [x, y, z] = &self.inner.boundary_flags;
        (x.clone(), y.clone(), z.clone())
    }

    #[getter]
    fn column_names(&self) -> Vec<String> {
        self.inner.column_names.clone()
    }

    /// `(num_atoms, num_columns)` numpy array; `rows[i][j]` is
    /// `column_names[j]` for atom `i`, in file order.
    #[getter]
    fn rows<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let n = self.inner.rows.len();
        let m = self.inner.column_names.len();
        let flat: Vec<f64> = self.inner.rows.iter().flatten().copied().collect();
        Array2::from_shape_vec((n, m), flat)
            .expect("LammpsDumpFrame.rows is always rectangular (enforced at parse)")
            .into_pyarray(py)
    }

    /// Index of `name` within `column_names`, if present.
    fn column_index(&self, name: &str) -> Option<usize> {
        self.inner.column_index(name)
    }

    /// The values of column `name` across every atom, in file order, or
    /// `None` if `name` isn't a declared column.
    fn column(&self, name: &str) -> Option<Vec<f64>> {
        self.inner.column(name)
    }

    /// Real Cartesian positions per atom, as an `(N, 3)` numpy array, from
    /// whichever of `x y z` (passed through) or `xs ys zs` (box-scaled,
    /// including triclinic shear terms) is present. `None` if neither
    /// triple is fully present -- including when only `xu yu zu`
    /// ("unwrapped") is present; use :meth:`column` for those directly (see
    /// `chematic_mol::LammpsDumpFrame::cartesian_positions`'s Rust docs for
    /// why "unwrapped" isn't resolved here).
    fn cartesian_positions<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        let pts = self.inner.cartesian_positions()?;
        let n = pts.len();
        let flat: Vec<f64> = pts.into_iter().flatten().collect();
        Some(
            Array2::from_shape_vec((n, 3), flat)
                .expect("row width fixed at 3 always reshapes cleanly")
                .into_pyarray(py),
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "LammpsDumpFrame(timestep={}, num_atoms={}, columns={:?})",
            self.inner.timestep, self.inner.num_atoms, self.inner.column_names
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

/// Parse a single LAMMPS dump frame (the first frame, if `text` holds a
/// multi-frame trajectory -- use :func:`parse_lammps_dump_all` for the
/// rest).
///
/// Raises:
///     ValueError: on parse failure (see `chematic_mol::LammpsDumpError`).
#[pyfunction]
fn parse_lammps_dump_frame(text: &str) -> PyResult<PyLammpsDumpFrame> {
    let inner = chematic_mol::parse_lammps_dump_frame(text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(PyLammpsDumpFrame { inner })
}

/// Parse every frame of a LAMMPS dump/trajectory file.
///
/// Materializes the whole trajectory as a list rather than exposing
/// `chematic_mol::LammpsDumpReader`'s true streaming iteration to Python: a
/// dump trajectory is inherently multi-frame (unlike a single Gaussian Cube
/// grid, where materializing is essentially free either way -- see
/// :func:`VolumetricGrid.from_cube`'s docs), so this is a genuine, disclosed
/// scope decision, not a silently dropped capability: most Python callers
/// materialize a full frame list anyway (``list(reader)``), and real
/// streaming (a Python iterator yielding one frame at a time) can be added
/// later if a use case needs bounded memory for a huge trajectory.
///
/// Raises:
///     ValueError: on the first parse failure (stops there).
///
/// Example::
///
///     frames = chematic.parse_lammps_dump_all(open("dump.lammpstrj").read())
///     positions = frames[-1].cartesian_positions()
#[pyfunction]
fn parse_lammps_dump_all(text: &str) -> PyResult<Vec<PyLammpsDumpFrame>> {
    let reader = std::io::BufReader::new(std::io::Cursor::new(text.as_bytes()));
    chematic_mol::LammpsDumpReader::new(reader)
        .map(|r| {
            r.map(|inner| PyLammpsDumpFrame { inner })
                .map_err(|e| PyValueError::new_err(e.to_string()))
        })
        .collect()
}

#[pyfunction]
fn write_lammps_dump_frame(frame: &PyLammpsDumpFrame) -> String {
    chematic_mol::write_lammps_dump_frame(&frame.inner)
}

/// Write multiple frames as one trajectory file (plain concatenation -- a
/// LAMMPS trajectory is just N frames back to back, with no separator
/// beyond the next frame's own `ITEM: TIMESTEP`).
#[pyfunction]
fn write_lammps_trajectory(frames: Vec<PyRef<PyLammpsDumpFrame>>) -> String {
    let inner: Vec<chematic_mol::LammpsDumpFrame> =
        frames.iter().map(|f| f.inner.clone()).collect();
    chematic_mol::write_lammps_trajectory(&inner)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyLammpsDumpFrame>()?;
    m.add_function(wrap_pyfunction!(parse_lammps_data, m)?)?;
    m.add_function(wrap_pyfunction!(write_lammps_data, m)?)?;
    m.add_function(wrap_pyfunction!(box_bounds_to_true, m)?)?;
    m.add_function(wrap_pyfunction!(true_to_box_bounds, m)?)?;
    m.add_function(wrap_pyfunction!(parse_lammps_dump_frame, m)?)?;
    m.add_function(wrap_pyfunction!(parse_lammps_dump_all, m)?)?;
    m.add_function(wrap_pyfunction!(write_lammps_dump_frame, m)?)?;
    m.add_function(wrap_pyfunction!(write_lammps_trajectory, m)?)?;
    Ok(())
}
