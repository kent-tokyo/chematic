//! Similarity/fingerprint comparison bindings (Tanimoto/Dice/Tversky variants, clustering, alignment).

use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Tanimoto similarity between two fingerprint byte arrays.
///
/// Works with any equal-length ``bytes`` objects (ECFP4, ECFP6, MACCS, …)::
///
///     sim = chematic.tanimoto(mol1.ecfp4(), mol2.ecfp4())
#[pyfunction]
fn tanimoto(a: &[u8], b: &[u8]) -> PyResult<f64> {
    if a.len() != b.len() {
        return Err(PyValueError::new_err(format!(
            "fingerprints must be the same length ({} vs {})",
            a.len(),
            b.len()
        )));
    }
    let and_bits: u32 = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x & y).count_ones())
        .sum();
    let or_bits: u32 = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x | y).count_ones())
        .sum();
    if or_bits == 0 {
        return Ok(0.0);
    }
    Ok(and_bits as f64 / or_bits as f64)
}

/// Screen a SMILES library by 3D shape similarity to a query molecule.
///
/// Returns ``[(index, similarity), ...]`` sorted by decreasing similarity.
/// 3D coordinates are auto-generated via distance geometry for each molecule.
///
///     hits = chematic.shape_screen(query, smiles_list)
///     for idx, sim in hits[:10]:
///         print(f"{smiles_list[idx]}  sim={sim:.3f}")
#[pyfunction]
fn shape_screen(query: &Mol, smiles_list: Vec<String>) -> Vec<(usize, f64)> {
    let mols: Vec<chematic_core::Molecule> = smiles_list
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();
    let refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    chematic_3d::shape_screen(&query.inner, &refs)
}

/// Estimate MAP4 Tanimoto similarity between two MAP4 fingerprints.
///
/// ``a`` and ``b`` must be lists of 1024 integers as returned by :meth:`Mol.map4`.
/// Returns a value in [0, 1].
///
///     sim = chematic.tanimoto_map4(mol1.map4(), mol2.map4())
#[pyfunction]
fn tanimoto_map4(a: Vec<u32>, b: Vec<u32>) -> f64 {
    chematic_fp::tanimoto_map4(&a, &b)
}

/// Tanimoto-like similarity between two Spectrophores fingerprints.
///
/// Uses the USR formula ``S = 1 / (1 + mean|a − b|)``, returning values in (0, 1].
/// Both vectors must have the same length (typically 48).
///
///     coords1 = mol1.generate_3d()
///     coords2 = mol2.generate_3d()
///     fp1 = mol1.spectrophores(coords1)
///     fp2 = mol2.spectrophores(coords2)
///     sim = chematic.tanimoto_spectrophores(fp1, fp2)
#[pyfunction]
fn tanimoto_spectrophores(a: Vec<f64>, b: Vec<f64>) -> f64 {
    chematic_3d::tanimoto_spectrophores(&a, &b)
}

/// Butina clustering — group molecules by ECFP4 Tanimoto similarity.
///
/// Returns a list of clusters; each cluster is a list of SMILES indices (centroid first).
/// Clusters are sorted by size (largest first).
///
/// Args:
///     smiles: list of SMILES strings.
///     cutoff: Tanimoto similarity threshold (default 0.65 — molecules ≥ cutoff → same cluster).
///
///     clusters = chematic.butina_cluster(smiles, 0.65)
///     for c in clusters[:5]:
///         print(f"cluster centroid: {smiles[c[0]]}, size: {len(c)}")
#[pyfunction]
#[pyo3(signature = (smiles, cutoff = 0.65))]
fn butina_cluster(smiles: Vec<String>, cutoff: f64) -> Vec<Vec<usize>> {
    let mols: Vec<chematic_core::Molecule> = smiles
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();
    chematic_chem::butina_cluster(&mols, cutoff, chematic_fp::tanimoto_ecfp4)
}

/// MaxMin diversity picking — select `n` maximally diverse molecules.
///
/// Returns a list of indices into the ``smiles`` list, in selection order.
/// Uses ECFP4 Tanimoto distance.
///
///     picks = chematic.maxmin_picks(smiles, 100)
///     diverse_set = [smiles[i] for i in picks]
#[pyfunction]
fn maxmin_picks(smiles: Vec<String>, n: usize) -> Vec<usize> {
    let mols: Vec<chematic_core::Molecule> = smiles
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();
    chematic_chem::maxmin_picks(&mols, n, chematic_fp::tanimoto_ecfp4)
}

/// Cosine similarity between two ERG feature vectors.
///
/// Both ``a`` and ``b`` must have length 315 (from :meth:`Mol.erg_vec`).
/// Returns a value in [0, 1].
///
///     sim = chematic.cosine_erg_vec(mol1.erg_vec(), mol2.erg_vec())
#[pyfunction]
fn cosine_erg_vec(a: Vec<f64>, b: Vec<f64>) -> PyResult<f64> {
    const LEN: usize = chematic_fp::ERG_VEC_LEN;
    if a.len() != LEN || b.len() != LEN {
        return Err(PyValueError::new_err(format!(
            "erg_vec must have length {LEN}, got {} and {}",
            a.len(),
            b.len()
        )));
    }
    let a_arr: &[f64; LEN] = a.as_slice().try_into().unwrap();
    let b_arr: &[f64; LEN] = b.as_slice().try_into().unwrap();
    Ok(chematic_fp::cosine_erg_vec(a_arr, b_arr))
}

/// Tanimoto similarity between two ERG feature vectors.
///
/// Both ``a`` and ``b`` must have length 315 (from :meth:`Mol.erg_vec`).
/// Returns a value in [0, 1].
///
///     sim = chematic.tanimoto_erg_vec(mol1.erg_vec(), mol2.erg_vec())
#[pyfunction]
fn tanimoto_erg_vec(a: Vec<f64>, b: Vec<f64>) -> PyResult<f64> {
    const LEN: usize = chematic_fp::ERG_VEC_LEN;
    if a.len() != LEN || b.len() != LEN {
        return Err(PyValueError::new_err(format!(
            "erg_vec must have length {LEN}, got {} and {}",
            a.len(),
            b.len()
        )));
    }
    let a_arr: &[f64; LEN] = a.as_slice().try_into().unwrap();
    let b_arr: &[f64; LEN] = b.as_slice().try_into().unwrap();
    Ok(chematic_fp::tanimoto_erg_vec(a_arr, b_arr))
}

/// Find the top-K most similar molecules using a selectable fingerprint type.
///
/// Like :func:`top_k_similar` but lets you choose the fingerprint used for
/// Tanimoto scoring. Supported ``fp`` values:
///
///   - ``"ecfp4"`` (default) — ECFP4, 2048-bit
///   - ``"ecfp6"`` — ECFP6, 2048-bit
///   - ``"ecfp4_chiral"`` — ECFP4 with chirality
///   - ``"fcfp4"`` — FCFP4 feature-based
///   - ``"maccs"`` — 166-bit MACCS keys
///   - ``"topo_path"`` — topological path FP
///
///     results = chematic.top_k_similar_fp("c1ccccc1", smiles_list, k=5, fp="maccs")
///     for idx, score in results:
///         print(smiles_list[idx], score)
#[pyfunction]
#[pyo3(signature = (query, smiles, k=10, fp=None))]
fn top_k_similar_fp(
    query: &str,
    smiles: Vec<String>,
    k: usize,
    fp: Option<&str>,
) -> PyResult<Vec<(usize, f64)>> {
    use chematic_fp::search::FpType;
    let fp_type = match fp.unwrap_or("ecfp4") {
        "ecfp6" => FpType::Ecfp6,
        "ecfp4_chiral" => FpType::Ecfp4Chiral,
        "fcfp4" => FpType::Fcfp4,
        "maccs" => FpType::Maccs,
        "topo_path" => FpType::TopoPath,
        _ => FpType::Ecfp4,
    };
    let query_mol =
        chematic_smiles::parse(query).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let db: Vec<chematic_core::Molecule> = smiles
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();
    Ok(chematic_fp::search::nearest_neighbors(
        &query_mol, &db, k, fp_type,
    ))
}

/// Find the top-K most similar molecules to a query by ECFP4 Tanimoto.
///
/// More memory-efficient than computing the full similarity matrix for large libraries.
///
/// Returns a list of ``(index, score)`` tuples sorted by descending similarity.
/// Invalid SMILES are silently skipped; returned indices refer to the original list.
///
///     hits = chematic.top_k_similar("c1ccccc1", smiles_library, k=10)
///     for idx, score in hits:
///         print(f"{smiles_library[idx]}: {score:.3f}")
#[pyfunction]
fn top_k_similar(query: &str, smiles: Vec<String>, k: usize) -> PyResult<Vec<(usize, f32)>> {
    let query_mol =
        chematic_smiles::parse(query).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let query_fp = chematic_fp::ecfp4(&query_mol);
    let db_fps: Vec<chematic_fp::bitvec::BitVec2048> = smiles
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .map(|m| chematic_fp::ecfp4(&m))
        .collect();
    Ok(chematic_fp::top_k_similar(&query_fp, &db_fps, k))
}

/// Dice similarity between two fingerprint byte arrays.
///
/// Dice = 2 * |A∩B| / (|A| + |B|). Works with the same byte fingerprints as
/// :func:`tanimoto` (ECFP4, ECFP6, MACCS, ERG, pharmacophore, …).
/// Returns a value in [0, 1].
///
///     sim = chematic.dice_similarity(mol1.ecfp4(), mol2.ecfp4())
#[pyfunction]
fn dice_similarity(a: &[u8], b: &[u8]) -> f64 {
    if a.len() != b.len() {
        return 0.0;
    }
    let and_bits: u32 = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x & y).count_ones())
        .sum();
    let a_bits: u32 = a.iter().map(|x| x.count_ones()).sum();
    let b_bits: u32 = b.iter().map(|x| x.count_ones()).sum();
    if a_bits + b_bits == 0 {
        0.0
    } else {
        2.0 * and_bits as f64 / (a_bits + b_bits) as f64
    }
}

/// Tversky similarity between two fingerprint byte arrays.
///
/// Tversky(α, β) = |A∩B| / (α|A\B| + β|B\A| + |A∩B|).
///
/// - α=β=0.5 → Dice similarity
/// - α=β=1.0 → Tanimoto similarity
/// - α=0, β=1 → recall-oriented (sub-structure search bias)
///
///     sim = chematic.tversky_similarity(query.ecfp4(), target.ecfp4(), 0.0, 1.0)
#[pyfunction]
fn tversky_similarity(a: &[u8], b: &[u8], alpha: f64, beta: f64) -> f64 {
    if a.len() != b.len() {
        return 0.0;
    }
    let and_bits: f64 = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x & y).count_ones() as f64)
        .sum();
    let a_not_b: f64 = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x & !y).count_ones() as f64)
        .sum();
    let b_not_a: f64 = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (!x & y).count_ones() as f64)
        .sum();
    let denom = alpha * a_not_b + beta * b_not_a + and_bits;
    if denom == 0.0 { 0.0 } else { and_bits / denom }
}

/// Estimate MHFP Tanimoto similarity between two MHFP fingerprints.
///
/// ``a`` and ``b`` must be lists of 128 u64 values as returned by :meth:`Mol.mhfp`.
/// Uses position-wise matching (not bitwise AND/OR).
///
///     sim = chematic.tanimoto_mhfp(mol1.mhfp(), mol2.mhfp())
#[pyfunction]
fn tanimoto_mhfp(a: Vec<u64>, b: Vec<u64>) -> f64 {
    if a.len() != b.len() || a.is_empty() {
        return 0.0;
    }
    let matches = a.iter().zip(b.iter()).filter(|(x, y)| x == y).count();
    matches as f64 / a.len() as f64
}

/// Align ``probe`` coordinates onto ``reference`` using the Kabsch algorithm.
///
/// Both ``probe`` and ``reference`` must be lists of ``[x,y,z]`` lists with
/// the same number of atoms (atom-to-atom correspondence is assumed).
///
/// Returns ``(aligned_coords, rmsd)`` where ``aligned_coords`` is a new
/// ``[[x,y,z], ...]`` list with the probe optimally superposed on reference.
///
///     coords = mol.generate_3d()
///     aligned, rmsd = chematic.align_coords(coords, ref_coords)
#[pyfunction]
fn align_coords(probe: Vec<[f64; 3]>, reference: Vec<[f64; 3]>) -> (Vec<Vec<f64>>, f64) {
    let result = chematic_3d::align_coords(&reference, &probe);
    let aligned = chematic_3d::apply_alignment(&probe, &result);
    let py_coords: Vec<Vec<f64>> = aligned.iter().map(|c| vec![c[0], c[1], c[2]]).collect();
    (py_coords, result.rmsd)
}

/// Compute RMSD between two sets of paired 3D coordinates **without** alignment.
///
/// ``coords_a`` and ``coords_b`` must have the same number of atoms.
/// Returns RMSD in the same units as the input (typically Å).
///
///     rmsd = chematic.rmsd(mol.generate_3d(), ref_coords)
#[pyfunction]
fn rmsd(coords_a: Vec<[f64; 3]>, coords_b: Vec<[f64; 3]>) -> f64 {
    chematic_3d::rmsd_no_align(&coords_a, &coords_b)
}

/// Tanimoto similarity between two molecules using ERG fingerprints.
///
/// Convenience alternative to ``chematic.tanimoto_erg_vec(m1.erg_vec(), m2.erg_vec())``.
/// Both ERG fingerprints are computed internally.
///
///     sim = chematic.tanimoto_erg(mol1, mol2)
#[pyfunction]
fn tanimoto_erg(mol1: &Mol, mol2: &Mol) -> f64 {
    chematic_fp::tanimoto_erg(&mol1.inner, &mol2.inner)
}

/// Compute an M×N Tanimoto similarity matrix from two lists of fingerprint byte arrays.
///
/// Returns a list of M rows, each row containing N Tanimoto scores:
/// ``result[i][j] = Tanimoto(fps_a[i], fps_b[j])``.
/// All fingerprints must have the same byte length (e.g., all from :meth:`Mol.ecfp4`).
///
///     matrix = chematic.tanimoto_matrix(
///         [m.ecfp4() for m in queries],
///         [m.ecfp4() for m in library],
///     )
///     # matrix[i][j] = similarity of query i against library compound j
#[pyfunction]
fn tanimoto_matrix(fps_a: Vec<Vec<u8>>, fps_b: Vec<Vec<u8>>) -> Vec<Vec<f32>> {
    let db_counts: Vec<u32> = fps_b
        .iter()
        .map(|fp| fp.iter().map(|b| b.count_ones()).sum())
        .collect();
    fps_a
        .iter()
        .map(|qa| {
            let qa_count: u32 = qa.iter().map(|b| b.count_ones()).sum();
            fps_b
                .iter()
                .zip(db_counts.iter())
                .map(|(qb, &db_cnt)| {
                    if qa.len() != qb.len() {
                        return 0.0;
                    }
                    let and: u32 = qa
                        .iter()
                        .zip(qb.iter())
                        .map(|(a, b)| (a & b).count_ones())
                        .sum();
                    let or = qa_count + db_cnt - and;
                    if or == 0 { 0.0 } else { and as f32 / or as f32 }
                })
                .collect()
        })
        .collect()
}

/// Compute Tanimoto similarity of one fingerprint against a list of fingerprints.
///
/// All byte arrays must be the same length (e.g., all from :meth:`Mol.ecfp4`).
/// More efficient than repeated :func:`tanimoto` calls for virtual screening.
///
///     db_fps = [mol.ecfp4() for mol in library]
///     scores = chematic.tanimoto_slice(query.ecfp4(), db_fps)
///     top = sorted(enumerate(scores), key=lambda x: -x[1])[:10]
#[pyfunction]
fn tanimoto_slice(query: &[u8], db: Vec<Vec<u8>>) -> Vec<f32> {
    let qa: u32 = query.iter().map(|b| b.count_ones()).sum();
    db.iter()
        .map(|fp| {
            if fp.len() != query.len() {
                return 0.0;
            }
            let and: u32 = query
                .iter()
                .zip(fp.iter())
                .map(|(a, b)| (a & b).count_ones())
                .sum();
            let db_a: u32 = fp.iter().map(|b| b.count_ones()).sum();
            let or = qa + db_a - and;
            if or == 0 { 0.0 } else { and as f32 / or as f32 }
        })
        .collect()
}

/// Tanimoto similarity between two 3D pharmacophore fingerprints.
///
/// Both ``a`` and ``b`` must be byte arrays from :meth:`Mol.pharmacophore_fp_3d`.
/// Returns a value in [0, 1].
///
///     fp1 = mol1.pharmacophore_fp_3d(coords1)
///     fp2 = mol2.pharmacophore_fp_3d(coords2)
///     sim = chematic.tanimoto_pharmacophore_3d(fp1, fp2)
#[pyfunction]
fn tanimoto_pharmacophore_3d(a: &[u8], b: &[u8]) -> PyResult<f64> {
    if a.len() != b.len() {
        return Err(PyValueError::new_err(format!(
            "fingerprints must be the same length ({} vs {})",
            a.len(),
            b.len()
        )));
    }
    let and_bits: u32 = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x & y).count_ones())
        .sum();
    let or_bits: u32 = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x | y).count_ones())
        .sum();
    if or_bits == 0 {
        return Ok(0.0);
    }
    Ok(and_bits as f64 / or_bits as f64)
}

/// `(pairs, aligned_coords2, rmsd, score)` — see [`o3a_align`].
type O3AAlignReturn = (Vec<(usize, usize)>, Vec<Vec<f64>>, f64, f64);

/// Find an atom correspondence between two molecules (O3A-style) and
/// superpose ``mol2`` onto ``mol1`` using it — unlike :func:`align_coords`,
/// this does not require the atom pairs to already be known.
///
/// Returns ``(pairs, aligned_coords2, rmsd, score)``:
///
/// - ``pairs``: list of ``(mol1_atom_idx, mol2_atom_idx)`` tuples
/// - ``aligned_coords2``: all of ``coords2`` superposed onto ``mol1``'s frame
/// - ``rmsd``: fit RMSD (Å), over the paired atoms only
/// - ``score``: Gaussian overlap score of the paired atoms after alignment
///   (higher means a tighter, more extensive overlap)
///
///     pairs, aligned, rmsd, score = chematic.o3a_align(mol1, coords1, mol2, coords2)
#[pyfunction]
fn o3a_align(
    mol1: &Mol,
    coords1: Vec<[f64; 3]>,
    mol2: &Mol,
    coords2: Vec<[f64; 3]>,
) -> PyResult<O3AAlignReturn> {
    let result = chematic_3d::o3a_align(&mol1.inner, &coords1, &mol2.inner, &coords2)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let aligned = chematic_3d::apply_alignment(&coords2, &result.alignment);
    let py_coords: Vec<Vec<f64>> = aligned.iter().map(|c| vec![c[0], c[1], c[2]]).collect();
    Ok((result.pairs, py_coords, result.alignment.rmsd, result.score))
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(tanimoto, m)?)?;
    m.add_function(wrap_pyfunction!(shape_screen, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_map4, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_spectrophores, m)?)?;
    m.add_function(wrap_pyfunction!(butina_cluster, m)?)?;
    m.add_function(wrap_pyfunction!(maxmin_picks, m)?)?;
    m.add_function(wrap_pyfunction!(cosine_erg_vec, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_erg_vec, m)?)?;
    m.add_function(wrap_pyfunction!(top_k_similar_fp, m)?)?;
    m.add_function(wrap_pyfunction!(top_k_similar, m)?)?;
    m.add_function(wrap_pyfunction!(dice_similarity, m)?)?;
    m.add_function(wrap_pyfunction!(tversky_similarity, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_mhfp, m)?)?;
    m.add_function(wrap_pyfunction!(align_coords, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_erg, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_slice, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_pharmacophore_3d, m)?)?;
    m.add_function(wrap_pyfunction!(o3a_align, m)?)?;
    Ok(())
}
