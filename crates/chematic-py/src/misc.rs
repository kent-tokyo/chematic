//! Miscellaneous bindings that don't fit another domain file (SMARTS matching, colors, abbreviations, depiction helpers).

use crate::Mol;
use crate::formats::flat_to_coords3d;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

/// Parse `smarts`, memoizing successful parses process-wide.
///
/// Python callers typically pass the same pattern string for every molecule
/// in a loop; re-parsing it each time dominated cheap existence checks.
/// Parse errors are never cached, so error values/messages are unchanged.
pub(crate) fn cached_smarts(
    smarts: &str,
) -> Result<Arc<chematic_smarts::QueryMolecule>, chematic_smarts::SmartsError> {
    use rustc_hash::FxHashMap;
    use std::sync::{Mutex, OnceLock};
    const CAPACITY: usize = 4096;
    // Keys are pattern strings chosen by the caller; Fx hashing keeps the
    // per-call lookup cheap (the cache is bounded, so no DoS concern).
    static CACHE: OnceLock<Mutex<FxHashMap<String, Arc<chematic_smarts::QueryMolecule>>>> =
        OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(FxHashMap::default()));
    if let Some(q) = cache.lock().ok().and_then(|c| c.get(smarts).cloned()) {
        return Ok(q);
    }
    let query = Arc::new(chematic_smarts::parse_smarts(smarts)?);
    if let Ok(mut c) = cache.lock() {
        if c.len() >= CAPACITY {
            c.clear();
        }
        c.insert(smarts.to_owned(), Arc::clone(&query));
    }
    Ok(query)
}

/// A SMARTS query compiled once and reusable across molecules.
///
/// Use this for a query that is applied repeatedly; it avoids the process-wide
/// string-cache lookup performed by the convenience string APIs.
#[pyclass(name = "SmartsQuery", frozen)]
pub(crate) struct PySmartsQuery {
    source: String,
    query: Arc<chematic_smarts::QueryMolecule>,
}

/// Molecules retained in Rust for repeated batch operations.
#[pyclass(name = "MoleculeBatch", frozen)]
pub(crate) struct PyMoleculeBatch {
    molecules: Vec<Arc<chematic_core::Molecule>>,
}

#[pymethods]
impl PyMoleculeBatch {
    #[new]
    fn new(molecules: Vec<PyRef<'_, Mol>>) -> Self {
        Self {
            molecules: molecules
                .into_iter()
                .map(|mol| Arc::clone(&mol.inner))
                .collect(),
        }
    }

    fn __len__(&self) -> usize {
        self.molecules.len()
    }
}

#[pymethods]
impl PySmartsQuery {
    #[new]
    fn new(smarts: &str) -> PyResult<Self> {
        let query = chematic_smarts::parse_smarts(smarts)
            .map_err(|e| crate::errors::malformed("smarts", e.to_string()))?;
        Ok(Self {
            source: smarts.to_owned(),
            query: Arc::new(query),
        })
    }

    #[getter]
    fn source(&self) -> &str {
        &self.source
    }

    /// Return whether this query has at least one match in ``mol``.
    fn matches(&self, mol: &Mol) -> bool {
        let config = chematic_smarts::MatchConfig {
            max_matches: Some(1),
            uniquify: false,
            ..chematic_smarts::MatchConfig::default()
        };
        chematic_smarts::has_match_perceived(&self.query, &mol.inner, &config)
    }

    /// Match this query against many molecules with one Python/Rust call.
    fn matches_many(&self, molecules: Vec<PyRef<'_, Mol>>) -> Vec<bool> {
        let config = chematic_smarts::MatchConfig {
            max_matches: Some(1),
            uniquify: false,
            ..chematic_smarts::MatchConfig::default()
        };
        molecules
            .into_iter()
            .map(|mol| chematic_smarts::has_match_perceived(&self.query, &mol.inner, &config))
            .collect()
    }

    /// Match a Rust-retained molecule batch without re-extracting Python objects.
    fn matches_batch(&self, batch: &PyMoleculeBatch) -> Vec<bool> {
        let config = chematic_smarts::MatchConfig {
            max_matches: Some(1),
            uniquify: false,
            ..chematic_smarts::MatchConfig::default()
        };
        batch
            .molecules
            .iter()
            .map(|mol| chematic_smarts::has_match_perceived(&self.query, mol, &config))
            .collect()
    }

    /// Return sorted target-atom indices for every unique match.
    fn find_matches(&self, mol: &Mol) -> Vec<Vec<usize>> {
        chematic_smarts::find_match_atom_sets_perceived(
            &self.query,
            &mol.inner,
            &chematic_smarts::MatchConfig::default(),
        )
    }

    fn __repr__(&self) -> String {
        format!("SmartsQuery({:?})", self.source)
    }
}

/// Compile a SMARTS query for repeated matching.
#[pyfunction]
fn compile_smarts(smarts: &str) -> PyResult<PySmartsQuery> {
    PySmartsQuery::new(smarts)
}

/// Test whether a SMARTS pattern matches a molecule.
///
///     if chematic.smarts_match("[OH]", mol):
///         print("has hydroxyl")
#[pyfunction]
fn smarts_match(smarts: &str, mol: &Mol) -> PyResult<bool> {
    let query =
        cached_smarts(smarts).map_err(|e| crate::errors::malformed("smarts", e.to_string()))?;
    // Stop at the first embedding instead of enumerating every match — an
    // existence check doesn't need the full match set or the dedup pass.
    let config = chematic_smarts::MatchConfig {
        max_matches: Some(1),
        uniquify: false,
        ..chematic_smarts::MatchConfig::default()
    };
    Ok(chematic_smarts::has_match_perceived(
        &query, &mol.inner, &config,
    ))
}

/// Return all substructure matches of a SMARTS pattern in a molecule.
///
/// Each match is a list of atom indices (in query-atom order).
/// Returns an empty list when there are no matches.
///
///     matches = chematic.smarts_find("[OH]", mol)
///     # → [[3], [7], ...]   (one list per match; each element is a mol atom index)
#[pyfunction]
fn smarts_find(smarts: &str, mol: &Mol) -> PyResult<Vec<Vec<usize>>> {
    let query =
        cached_smarts(smarts).map_err(|e| crate::errors::malformed("smarts", e.to_string()))?;
    let n = query.atom_count();
    Ok(chematic_smarts::find_matches_perceived(
        &query,
        &mol.inner,
        &chematic_smarts::MatchConfig::default(),
    )
    .into_iter()
    .map(|map| {
        (0..n)
            .filter_map(|qi| map.get(&qi).map(|a| a.0 as usize))
            .collect()
    })
    .collect())
}

/// Render a molecule SVG with atoms coloured by a weight vector.
///
/// ``mol``: :class:`Mol` to render.
/// ``weights``: list of floats, one per heavy atom.  Positive → blue, negative → red, zero → white.
///
///     weights = mol.logp_per_atom()
///     svg = chematic.similarity_map_svg(mol, weights)
#[pyfunction]
fn similarity_map_svg(mol: &Mol, weights: Vec<f64>) -> String {
    chematic_mol::stereo_depiction::with_stereo_depiction(&mol.inner, |m, layout| {
        chematic_depict::render_svg_opts(
            m,
            layout,
            &chematic_depict::similarity_map_options(
                &mol.inner,
                &weights,
                &chematic_depict::RenderOptions::default(),
            ),
        )
    })
}

/// Return all known chemical abbreviations as a dict ``{symbol: SMILES}``.
///
/// Symbols include ``"Boc"``, ``"Cbz"``, ``"Ts"``, ``"Ph"``, ``"OMe"``, …
///
///     abbrevs = chematic.abbreviations()
///     print(abbrevs.get("Ph"))  # "c1ccccc1"
#[pyfunction]
fn abbreviations<'py>(py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let d = pyo3::types::PyDict::new(py);
    for (sym, smi) in chematic_chem::abbreviations() {
        d.set_item(sym, smi)?;
    }
    Ok(d)
}

/// Expand a chemical abbreviation to a :class:`Mol`.
///
/// Returns ``None`` if the symbol is unknown.
///
///     mol = chematic.expand_abbreviation("Ph")  # phenyl → Mol
///     if mol:
///         print(mol.smiles)  # c1ccccc1
#[pyfunction]
fn expand_abbreviation(symbol: &str) -> Option<Mol> {
    chematic_chem::expand_abbreviation(symbol).map(|mol| Mol {
        inner: Arc::new(mol),
        props: Default::default(),
    })
}

/// Translate all atom coordinates so the centroid is at the origin.
///
/// ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
/// Returns new centered coordinates.
///
///     centered = chematic.center_on_origin(mol.generate_3d())
#[pyfunction]
fn center_on_origin(coords: Vec<[f64; 3]>) -> Vec<Vec<f64>> {
    let c3d = flat_to_coords3d(&coords);
    let out = chematic_3d::center_on_origin(&c3d);
    out.points.iter().map(|p| vec![p.x, p.y, p.z]).collect()
}

/// Apply a 4×4 affine transformation matrix to 3D coordinates.
///
/// ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
/// ``matrix``: 4×4 homogeneous transformation matrix (row-major).
/// Returns new transformed coordinates.
///
///     import numpy as np
///     R = np.eye(4); R[:3, 3] = [1, 0, 0]  # translation by 1 Å in x
///     new_coords = chematic.transform_conformer(coords, R.tolist())
#[pyfunction]
fn transform_conformer(coords: Vec<[f64; 3]>, matrix: Vec<Vec<f64>>) -> PyResult<Vec<Vec<f64>>> {
    if matrix.len() != 4 || matrix.iter().any(|row| row.len() != 4) {
        return Err(PyValueError::new_err("matrix must be 4×4"));
    }
    let mat: [[f64; 4]; 4] = [
        [matrix[0][0], matrix[0][1], matrix[0][2], matrix[0][3]],
        [matrix[1][0], matrix[1][1], matrix[1][2], matrix[1][3]],
        [matrix[2][0], matrix[2][1], matrix[2][2], matrix[2][3]],
        [matrix[3][0], matrix[3][1], matrix[3][2], matrix[3][3]],
    ];
    let c3d = flat_to_coords3d(&coords);
    let out = chematic_3d::transform_conformer(&c3d, &mat);
    Ok(out.points.iter().map(|p| vec![p.x, p.y, p.z]).collect())
}

/// Look up a built-in named SMARTS pattern by name.
///
/// Returns the SMARTS string for well-known pharmacophore and functional group
/// patterns, or ``None`` if the name is unknown.
///
/// Available names (partial list):
///   ``"donor"``, ``"donor_strict"``, ``"acceptor"``, ``"acceptor_strict"``,
///   ``"aromatic"``, ``"aromatic_ring"``, ``"hydrophobic"``,
///   ``"positive"``, ``"negative"``.
///
///     if smarts := chematic.named_pattern("donor"):
///         hits = chematic.smarts_find(smarts, mol)
#[pyfunction]
fn named_pattern(name: &str) -> Option<&'static str> {
    chematic_smarts::named_pattern(name)
}

/// CSS color string for an element by atomic number.
///
/// Returns the CPK/standard coloring used by chematic's SVG renderer.
/// Useful for custom visualization code.
///
///     print(chematic.atom_color(8))   # "#FF0000" (oxygen = red)
///     print(chematic.atom_color(6))   # "#808080" (carbon = grey)
///     print(chematic.atom_color(7))   # "#0000FF" (nitrogen = blue)
#[pyfunction]
fn atom_color(atomic_num: u8) -> &'static str {
    chematic_depict::atom_color(atomic_num)
}

/// RGB color triple for an element by atomic number.
///
/// Returns the same color as :func:`atom_color` as a ``(R, G, B)`` tuple (0–255).
///
///     r, g, b = chematic.atom_color_rgb(8)   # (255, 0, 0) for oxygen
#[pyfunction]
fn atom_color_rgb(atomic_num: u8) -> (u8, u8, u8) {
    let [r, g, b] = chematic_depict::atom_color_rgb(atomic_num);
    (r, g, b)
}

/// Render a list of molecules as a grid SVG.
///
///     svg = chematic.depict_grid([mol1, mol2, mol3], cols=3)
#[pyfunction]
fn depict_grid(mols: Vec<Mol>, cols: usize) -> String {
    let (copies, layouts): (Vec<_>, Vec<chematic_depict::Layout>) = mols
        .iter()
        .map(|m| chematic_mol::stereo_depiction::depiction_with_stereo(&m.inner))
        .unzip();
    let refs: Vec<&chematic_core::Molecule> = mols
        .iter()
        .zip(&copies)
        .map(|(m, copy)| copy.as_ref().unwrap_or(m.inner.as_ref()))
        .collect();
    chematic_depict::depict_svg_grid_with_layouts(&refs, &layouts, cols)
}

/// Look up an element's atomic number by symbol (e.g. ``"O"`` → 8).
///
/// Raises ``ValueError`` for an unrecognized symbol. Used by
/// ``rdkit_compat.RWMol.AddAtom`` to accept element symbols.
///
///     chematic.element_atomic_number("O")  # 8
#[pyfunction]
fn element_atomic_number(symbol: &str) -> PyResult<u8> {
    chematic_core::Element::from_symbol(symbol)
        .map(|e| e.atomic_number())
        .ok_or_else(|| PyValueError::new_err(format!("unknown element symbol: {symbol:?}")))
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------

/// RDKit 2026.03.1's ``Chem.DetectChemistryProblems(
/// Chem.MolFromSmiles(smiles, sanitize=False))``: ``[(type, atom_indices)]``
/// in RDKit's order, ``type`` one of ``"AtomValenceException"``,
/// ``"AtomKekulizeException"``, ``"KekulizeException"``; atom indices in the
/// SMILES atom order. Raises ``ValueError`` for unparsable SMILES::
///
///     chematic.rdkit_detect_chemistry_problems("CC(C)(C)(C)C")
///     # [('AtomValenceException', [1])]
#[pyfunction]
fn rdkit_detect_chemistry_problems(smiles: &str) -> PyResult<Vec<(String, Vec<usize>)>> {
    let mol = chematic_smiles::parse_template(smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let problems = chematic_smiles::rdkit_detect_chemistry_problems(&mol)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(problems.into_iter().map(|p| (p.kind, p.atoms)).collect())
}

/// RDKit 2026.03.1's ``rdChemReactions.ReactionToSmarts(
/// rdChemReactions.ReactionFromSmarts(reaction_smarts))``: every template
/// re-written as RDKit writes query molecules (``[OH]`` -> ``[O&H1]``,
/// chirality relative to the output order, ``(A.B)`` groups in
/// parentheses). Raises ``ValueError`` for unsupported syntax (directional
/// or dative bonds, chirality classes)::
///
///     chematic.rdkit_reaction_to_smarts("[C:1](=[O:2])[OH]>>[C:1](=[O:2])N")
///     # '[C:1](=[O:2])[O&H1]>>[C:1](=[O:2])N'
#[pyfunction]
fn rdkit_reaction_to_smarts(reaction_smarts: &str) -> PyResult<String> {
    chematic_smiles::rdkit_reaction_to_smarts(reaction_smarts)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// RDKit 2026.03.1's ``Chem.MolToSmarts(Chem.MolFromSmarts(smarts))``.
#[pyfunction]
fn rdkit_smarts_to_smarts(smarts: &str) -> PyResult<String> {
    chematic_smiles::rdkit_smarts_to_smarts(smarts)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMoleculeBatch>()?;
    m.add_class::<PySmartsQuery>()?;
    m.add_function(wrap_pyfunction!(compile_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(rdkit_reaction_to_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(rdkit_smarts_to_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(rdkit_detect_chemistry_problems, m)?)?;
    m.add_function(wrap_pyfunction!(smarts_match, m)?)?;
    m.add_function(wrap_pyfunction!(smarts_find, m)?)?;
    m.add_function(wrap_pyfunction!(similarity_map_svg, m)?)?;
    m.add_function(wrap_pyfunction!(abbreviations, m)?)?;
    m.add_function(wrap_pyfunction!(expand_abbreviation, m)?)?;
    m.add_function(wrap_pyfunction!(center_on_origin, m)?)?;
    m.add_function(wrap_pyfunction!(transform_conformer, m)?)?;
    m.add_function(wrap_pyfunction!(named_pattern, m)?)?;
    m.add_function(wrap_pyfunction!(atom_color, m)?)?;
    m.add_function(wrap_pyfunction!(atom_color_rgb, m)?)?;
    m.add_function(wrap_pyfunction!(depict_grid, m)?)?;
    m.add_function(wrap_pyfunction!(element_atomic_number, m)?)?;
    Ok(())
}
