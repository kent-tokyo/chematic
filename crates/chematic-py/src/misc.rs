//! Miscellaneous bindings that don't fit another domain file (SMARTS matching, colors, abbreviations, depiction helpers).

use crate::Mol;
use crate::formats::flat_to_coords3d;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

/// Test whether a SMARTS pattern matches a molecule.
///
///     if chematic.smarts_match("[OH]", mol):
///         print("has hydroxyl")
#[pyfunction]
fn smarts_match(smarts: &str, mol: &Mol) -> PyResult<bool> {
    let query =
        chematic_smarts::parse_smarts(smarts).map_err(|e| PyValueError::new_err(e.to_string()))?;
    // Stop at the first embedding instead of enumerating every match — an
    // existence check doesn't need the full match set or the dedup pass.
    let config = chematic_smarts::MatchConfig {
        max_matches: Some(1),
        uniquify: false,
        ..chematic_smarts::MatchConfig::default()
    };
    Ok(!chematic_smarts::find_matches_with_config(&query, &mol.inner, &config).is_empty())
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
        chematic_smarts::parse_smarts(smarts).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let n = query.atom_count();
    Ok(chematic_smarts::find_matches(&query, &mol.inner)
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
    chematic_depict::similarity_map_svg(&mol.inner, &weights)
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
    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    chematic_depict::depict_svg_grid(&refs, cols)
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

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
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
