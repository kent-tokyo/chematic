use ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rayon::prelude::*;
use std::sync::Arc;

use crate::Mol;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse_parallel(smiles: &[String]) -> Vec<Option<chematic_core::Molecule>> {
    smiles
        .par_iter()
        .map(|s| chematic_smiles::parse(s).ok())
        .collect()
}

fn bitvec2048_to_bits(fp: &chematic_fp::bitvec::BitVec2048) -> Vec<u8> {
    (0..2048usize).map(|i| u8::from(fp.get(i))).collect()
}

// ---------------------------------------------------------------------------
// bulk.parse — parallel SMILES parsing
// ---------------------------------------------------------------------------

/// Parse a list of SMILES strings in parallel.
///
/// Returns a list of ``Mol`` objects (or ``None`` for invalid SMILES).
/// Invalid entries are silently skipped — check for ``None`` if needed.
///
///     mols = chematic.bulk.parse(["CCO", "c1ccccc1", "INVALID"])
///     valid = [m for m in mols if m is not None]
#[pyfunction]
pub fn parse(smiles: Vec<String>) -> Vec<Option<Mol>> {
    parse_parallel(&smiles)
        .into_iter()
        .map(|opt| {
            opt.map(|mol| Mol {
                inner: Arc::new(mol),
            })
        })
        .collect()
}

// ---------------------------------------------------------------------------
// bulk.ecfp4 — batch ECFP4 fingerprints as numpy (N, 2048) uint8
// ---------------------------------------------------------------------------

/// Compute ECFP4 fingerprints for a list of SMILES in parallel.
///
/// Returns a numpy array of shape ``(N, 2048)`` with ``dtype=uint8`` (0/1).
/// Invalid SMILES are silently skipped; ``N`` equals the number of valid molecules.
///
///     fps = chematic.bulk.ecfp4(smiles_list)  # shape (N, 2048)
///
///     # Use directly with scikit-learn:
///     from sklearn.decomposition import PCA
///     pca = PCA(n_components=50)
///     coords = pca.fit_transform(fps.astype(float))
#[pyfunction]
pub fn ecfp4<'py>(py: Python<'py>, smiles: Vec<String>) -> Bound<'py, PyArray2<u8>> {
    let bits: Vec<Vec<u8>> = smiles
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| bitvec2048_to_bits(&chematic_fp::ecfp4(&mol)))
        .collect();

    let n = bits.len();
    if n == 0 {
        return Array2::<u8>::zeros((0, 2048)).into_pyarray(py);
    }
    let flat: Vec<u8> = bits.into_iter().flatten().collect();
    Array2::from_shape_vec((n, 2048), flat)
        .expect("shape mismatch")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// bulk.maccs — batch MACCS 166-bit as numpy (N, 166) uint8
// ---------------------------------------------------------------------------

/// Compute MACCS 166-bit fingerprints for a list of SMILES in parallel.
///
/// Returns a numpy array of shape ``(N, 166)`` with ``dtype=uint8`` (0/1).
#[pyfunction]
pub fn maccs<'py>(py: Python<'py>, smiles: Vec<String>) -> Bound<'py, PyArray2<u8>> {
    let bits: Vec<Vec<u8>> = smiles
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| {
            let fp = chematic_fp::maccs(&mol);
            (0..166usize).map(|i| u8::from(fp.get(i))).collect()
        })
        .collect();

    let n = bits.len();
    if n == 0 {
        return Array2::<u8>::zeros((0, 166)).into_pyarray(py);
    }
    let flat: Vec<u8> = bits.into_iter().flatten().collect();
    Array2::from_shape_vec((n, 166), flat)
        .expect("shape mismatch")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// bulk.descriptors — parallel descriptor computation
// ---------------------------------------------------------------------------

/// Compute 55+ descriptors for a list of SMILES in parallel.
///
/// Returns a list of dicts. Each dict matches the output of ``mol.descriptors()``.
/// Invalid SMILES are silently skipped.
///
///     descs = chematic.bulk.descriptors(smiles_list)
///
///     # Build a Pandas DataFrame in one line:
///     import pandas as pd
///     df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))
#[pyfunction]
pub fn descriptors<'py>(py: Python<'py>, smiles: Vec<String>) -> PyResult<Vec<Bound<'py, PyDict>>> {
    // Phase 1: parallel computation (no GIL)
    struct Desc {
        mw: f64,
        exact_mass: f64,
        tpsa: f64,
        logp: f64,
        mr: f64,
        hbd: usize,
        hba: usize,
        rb: usize,
        hac: usize,
        rc: usize,
        arc: usize,
        nh: usize,
        nsc: usize,
        nsp: usize,
        nbh: usize,
        fsp3: f64,
        qed: f64,
        sa: f64,
        fc: i32,
        asa: f64,
        bertz: f64,
        wi: f64,
        k1: f64,
        k2: f64,
        k3: f64,
        c0: f64,
        c1: f64,
        c2: f64,
        c3: f64,
        c4: f64,
        c0v: f64,
        c1v: f64,
        c2v: f64,
        c3v: f64,
        c4v: f64,
        n_ah: usize,
        n_alh: usize,
        n_sr: usize,
        n_ar: usize,
        n_usc: usize,
        sum_e: f64,
        max_e: f64,
        min_e: f64,
        lip: bool,
        veb: bool,
        egan: bool,
        ghose: bool,
        reos: bool,
        pains: bool,
        bbb: f64,
        bbp: bool,
        caco: f64,
        herg: f64,
        cyp: f64,
        pka_acid: Option<f64>,
        pka_base: Option<f64>,
    }

    let descs: Vec<Desc> = smiles
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| {
            let m = &mol;
            // Pre-compute shared data: ring_bundle (1 SSSR), logp+MR (1 Crippen pass),
            // kappa (1 heavy_indices), chi (1 heavy_indices), estate (1 BFS pass), pKa (1 scan).
            let rb_data = chematic_chem::ring_bundle(m);
            let (pka_a, pka_b) = chematic_chem::pka_both(m);
            let (logp_val, mr_val) = chematic_chem::logp_and_mr(m);
            let (k1, k2, k3) = chematic_chem::kappa_all(m);
            let (c0, c1, c2, c3, c4, c0v, c1v, c2v, c3v, c4v) = chematic_chem::chi_all(m);
            let (sum_e, max_e, min_e) = chematic_chem::estate_all(m);
            Desc {
                mw: chematic_chem::molecular_weight(m),
                exact_mass: chematic_chem::exact_mass(m),
                tpsa: chematic_chem::tpsa(m),
                logp: logp_val,
                mr: mr_val,
                hbd: chematic_chem::hbd_count(m),
                hba: rb_data.hba_count,
                rb: rb_data.rotatable_bond_count,
                hac: chematic_chem::heavy_atom_count(m),
                rc: rb_data.ring_count,
                arc: rb_data.aromatic_ring_count,
                nh: chematic_chem::num_heteroatoms(m),
                nsc: chematic_chem::num_stereocenters(m),
                nsp: rb_data.num_spiro_atoms,
                nbh: rb_data.num_bridgehead_atoms,
                fsp3: chematic_chem::fsp3(m),
                qed: chematic_chem::qed(m),
                sa: chematic_chem::sa_score(m),
                fc: chematic_chem::formal_charge_sum(m),
                asa: chematic_chem::labute_asa(m),
                bertz: chematic_chem::bertz_ct(m),
                wi: chematic_chem::wiener_index(m),
                k1,
                k2,
                k3,
                c0,
                c1,
                c2,
                c3,
                c4,
                c0v,
                c1v,
                c2v,
                c3v,
                c4v,
                n_ah: rb_data.num_aromatic_heterocycles,
                n_alh: rb_data.num_aliphatic_heterocycles,
                n_sr: rb_data.num_saturated_rings,
                n_ar: rb_data.num_aliphatic_rings,
                n_usc: chematic_chem::num_unspecified_stereocenters(m),
                sum_e,
                max_e,
                min_e,
                lip: chematic_chem::lipinski_passes(m),
                veb: chematic_chem::veber_passes(m),
                egan: chematic_chem::egan_passes(m),
                ghose: chematic_chem::ghose_passes(m),
                reos: chematic_chem::reos_passes(m),
                pains: chematic_chem::pains_passes(m),
                bbb: chematic_chem::bbb_score(m),
                bbp: chematic_chem::bbb_passes(m),
                caco: chematic_chem::caco2_permeability(m),
                herg: chematic_chem::herg_risk_score(m),
                cyp: chematic_chem::cyp3a4_inhibition_risk(m),
                pka_acid: pka_a,
                pka_base: pka_b,
            }
        })
        .collect();

    // Phase 2: convert to Python dicts (sequential, GIL held)
    descs
        .into_iter()
        .map(|d| {
            let dict = PyDict::new(py);
            dict.set_item("mw", d.mw)?;
            dict.set_item("exact_mass", d.exact_mass)?;
            dict.set_item("tpsa", d.tpsa)?;
            dict.set_item("logp", d.logp)?;
            dict.set_item("molar_refractivity", d.mr)?;
            dict.set_item("hbd", d.hbd)?;
            dict.set_item("hba", d.hba)?;
            dict.set_item("rotatable_bonds", d.rb)?;
            dict.set_item("heavy_atoms", d.hac)?;
            dict.set_item("ring_count", d.rc)?;
            dict.set_item("aromatic_ring_count", d.arc)?;
            dict.set_item("num_heteroatoms", d.nh)?;
            dict.set_item("num_stereocenters", d.nsc)?;
            dict.set_item("num_spiro_atoms", d.nsp)?;
            dict.set_item("num_bridgehead_atoms", d.nbh)?;
            dict.set_item("fsp3", d.fsp3)?;
            dict.set_item("qed", d.qed)?;
            dict.set_item("sa_score", d.sa)?;
            dict.set_item("formal_charge", d.fc)?;
            dict.set_item("labute_asa", d.asa)?;
            dict.set_item("bertz_ct", d.bertz)?;
            dict.set_item("wiener_index", d.wi)?;
            dict.set_item("kappa1", d.k1)?;
            dict.set_item("kappa2", d.k2)?;
            dict.set_item("kappa3", d.k3)?;
            dict.set_item("chi0", d.c0)?;
            dict.set_item("chi1", d.c1)?;
            dict.set_item("chi2", d.c2)?;
            dict.set_item("chi3", d.c3)?;
            dict.set_item("chi4", d.c4)?;
            dict.set_item("chi0v", d.c0v)?;
            dict.set_item("chi1v", d.c1v)?;
            dict.set_item("chi2v", d.c2v)?;
            dict.set_item("chi3v", d.c3v)?;
            dict.set_item("chi4v", d.c4v)?;
            dict.set_item("num_aromatic_heterocycles", d.n_ah)?;
            dict.set_item("num_aliphatic_heterocycles", d.n_alh)?;
            dict.set_item("num_saturated_rings", d.n_sr)?;
            dict.set_item("num_aliphatic_rings", d.n_ar)?;
            dict.set_item("num_unspecified_stereocenters", d.n_usc)?;
            dict.set_item("sum_estate", d.sum_e)?;
            dict.set_item("max_estate", d.max_e)?;
            dict.set_item("min_estate", d.min_e)?;
            dict.set_item("lipinski_passes", d.lip)?;
            dict.set_item("veber_passes", d.veb)?;
            dict.set_item("egan_passes", d.egan)?;
            dict.set_item("ghose_passes", d.ghose)?;
            dict.set_item("reos_passes", d.reos)?;
            dict.set_item("pains_passes", d.pains)?;
            dict.set_item("bbb_score", d.bbb)?;
            dict.set_item("bbb_passes", d.bbp)?;
            dict.set_item("caco2", d.caco)?;
            dict.set_item("herg_risk", d.herg)?;
            dict.set_item("cyp3a4_risk", d.cyp)?;
            match d.pka_acid {
                Some(v) => dict.set_item("pka_acid", v)?,
                None => dict.set_item("pka_acid", py.None())?,
            }
            match d.pka_base {
                Some(v) => dict.set_item("pka_base", v)?,
                None => dict.set_item("pka_base", py.None())?,
            }
            Ok(dict)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// bulk.tanimoto — pairwise Tanimoto similarity matrix
// ---------------------------------------------------------------------------

/// Compute an M×N Tanimoto similarity matrix (ECFP4) in parallel.
///
/// Args:
///     smiles_a: list of M SMILES strings (query molecules)
///     smiles_b: list of N SMILES strings (database molecules)
///
/// Returns:
///     numpy array of shape ``(M, N)`` with ``dtype=float32``.
///     Invalid SMILES are silently skipped.
///
///     # Self-similarity matrix:
///     sim = chematic.bulk.tanimoto(smiles, smiles)
///     # → (N, N) matrix; diagonal ≈ 1.0
///
///     # Query vs database:
///     sim = chematic.bulk.tanimoto(query_smiles, db_smiles)
///     # → (len(query_smiles), len(db_smiles))
#[pyfunction]
pub fn tanimoto<'py>(
    py: Python<'py>,
    smiles_a: Vec<String>,
    smiles_b: Vec<String>,
) -> Bound<'py, PyArray2<f32>> {
    let fps_a: Vec<chematic_fp::bitvec::BitVec2048> = smiles_a
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| chematic_fp::ecfp4(&mol))
        .collect();

    let fps_b: Vec<chematic_fp::bitvec::BitVec2048> = smiles_b
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| chematic_fp::ecfp4(&mol))
        .collect();

    let m = fps_a.len();
    let n = fps_b.len();

    if m == 0 || n == 0 {
        return Array2::<f32>::zeros((m, n)).into_pyarray(py);
    }

    // Re-use the existing Rust bulk Tanimoto matrix (row-major flat Vec)
    let flat = chematic_fp::bulk::tanimoto_matrix(&fps_a, &fps_b);
    Array2::from_shape_vec((m, n), flat)
        .expect("shape mismatch")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// bulk.tanimoto_search — 1:N similarity search
// ---------------------------------------------------------------------------

/// Compute Tanimoto similarity between one query molecule and N database molecules.
///
/// Returns a numpy array of shape ``(N,)`` with ``dtype=float32``.
///
///     scores = chematic.bulk.tanimoto_search("c1ccccc1", db_smiles)
///     top10_idx = np.argsort(scores)[::-1][:10]
#[pyfunction]
pub fn tanimoto_search<'py>(
    py: Python<'py>,
    query: String,
    smiles: Vec<String>,
) -> PyResult<Bound<'py, PyArray1<f32>>> {
    let query_mol = chematic_smiles::parse(&query)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    let query_fp = chematic_fp::ecfp4(&query_mol);

    let db_fps: Vec<chematic_fp::bitvec::BitVec2048> = smiles
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| chematic_fp::ecfp4(&mol))
        .collect();

    if db_fps.is_empty() {
        return Ok(Array1::<f32>::zeros(0).into_pyarray(py));
    }

    let scores = chematic_fp::bulk::tanimoto_slice(&query_fp, &db_fps);
    Ok(Array1::from_vec(scores).into_pyarray(py))
}

// ---------------------------------------------------------------------------
// bulk.map4 — batch MAP4 fingerprints as numpy (N, 1024) uint32
// ---------------------------------------------------------------------------

/// Compute MAP4 fingerprints for a list of SMILES in parallel.
///
/// Returns a numpy array of shape ``(N, 1024)`` with ``dtype=uint32``.
/// Use :func:`chematic.tanimoto_map4` (not :func:`chematic.tanimoto`) for similarity.
/// Invalid SMILES are silently skipped.
///
///     fps = chematic.bulk.map4(smiles_list)  # shape (N, 1024), dtype uint32
#[pyfunction]
pub fn map4<'py>(py: Python<'py>, smiles: Vec<String>) -> Bound<'py, PyArray2<u32>> {
    let fps: Vec<Vec<u32>> = smiles
        .par_iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|mol| chematic_fp::map4_default(&mol))
        .collect();

    let n = fps.len();
    if n == 0 {
        return Array2::<u32>::zeros((0, 1024)).into_pyarray(py);
    }
    let flat: Vec<u32> = fps.into_iter().flatten().collect();
    Array2::from_shape_vec((n, 1024), flat)
        .expect("shape mismatch")
        .into_pyarray(py)
}

// ---------------------------------------------------------------------------
// bulk.substructure_search — parallel SMARTS screen
// ---------------------------------------------------------------------------

/// Screen a list of SMILES against a SMARTS pattern in parallel.
///
/// Returns a ``list[bool]`` of length N (same order as ``smiles``).
/// ``True`` if the molecule matches; ``False`` if it does not match or if the
/// SMILES is invalid.  The SMARTS is compiled once and shared across threads.
///
///     hits = chematic.bulk.substructure_search("[nH]", smiles_list)
///     actives = [s for s, h in zip(smiles_list, hits) if h]
#[pyfunction]
#[pyo3(signature = (smarts, smiles))]
pub fn substructure_search(smarts: &str, smiles: Vec<String>) -> PyResult<Vec<bool>> {
    let query = chematic_smarts::parse_smarts(smarts)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let results: Vec<bool> = smiles
        .par_iter()
        .map(|smi| {
            chematic_smiles::parse(smi)
                .map(|mol| !chematic_smarts::find_matches(&query, &mol).is_empty())
                .unwrap_or(false)
        })
        .collect();

    Ok(results)
}

// ---------------------------------------------------------------------------
// Register the bulk submodule
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_function(wrap_pyfunction!(ecfp4, m)?)?;
    m.add_function(wrap_pyfunction!(maccs, m)?)?;
    m.add_function(wrap_pyfunction!(map4, m)?)?;
    m.add_function(wrap_pyfunction!(descriptors, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_search, m)?)?;
    m.add_function(wrap_pyfunction!(substructure_search, m)?)?;
    Ok(())
}
