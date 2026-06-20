use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

// ---------------------------------------------------------------------------
// SimilarityIndex — MinHash LSH index for fast similarity search
// ---------------------------------------------------------------------------

/// A locality-sensitive hashing (LSH) index for fast approximate similarity search.
///
/// Uses MinHash fingerprints (MHFP) with band decomposition to find molecules
/// with Tanimoto similarity ≥ threshold in sub-linear time.
///
/// **RDKit does not have this built-in.** This is a unique chematic feature.
///
/// Example::
///
///     import chematic
///
///     # Build index from a library of SMILES
///     idx = chematic.SimilarityIndex()
///     smiles_db = ["c1ccccc1", "CCO", "CC(=O)O", "c1cccnc1", "c1ccoc1"]
///     for smi in smiles_db:
///         idx.add(smi)
///
///     # Search: find molecules with Tanimoto ≥ 0.5 to the query
///     hits = idx.search("c1ccccc1N", threshold=0.5)
///     # → [(index_in_db, tanimoto_score), ...]
///
///     # Or build in bulk and search:
///     idx2 = chematic.SimilarityIndex.from_smiles(smiles_db)
///     hits = idx2.search("c1ccccc1N", threshold=0.3, k=10)
#[pyclass(name = "SimilarityIndex")]
pub struct PySimilarityIndex {
    inner: chematic_fp::MhfpLshIndex,
    smiles: Vec<String>,
}

#[pymethods]
impl PySimilarityIndex {
    /// Create a new empty index.
    ///
    /// ``num_hashes`` controls the fingerprint length (128 by default).
    /// Higher values improve recall at the cost of memory.
    #[new]
    #[pyo3(signature = (num_hashes = 128))]
    fn new(num_hashes: usize) -> Self {
        PySimilarityIndex {
            inner: chematic_fp::MhfpLshIndex::new(num_hashes),
            smiles: Vec::new(),
        }
    }

    /// Create an index pre-loaded with a list of SMILES strings.
    ///
    ///     idx = chematic.SimilarityIndex.from_smiles(["CCO", "c1ccccc1", ...])
    #[staticmethod]
    fn from_smiles(smiles: Vec<String>) -> PyResult<Self> {
        let mut idx = PySimilarityIndex::new(128);
        for smi in &smiles {
            idx.add(smi.clone())?;
        }
        Ok(idx)
    }

    /// Add one molecule (by SMILES) to the index. Returns its 0-based index.
    ///
    ///     idx = chematic.SimilarityIndex()
    ///     i = idx.add("c1ccccc1")   # → 0
    ///     j = idx.add("CCO")         # → 1
    fn add(&mut self, smiles: String) -> PyResult<usize> {
        let mol =
            chematic_smiles::parse(&smiles).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let fp = chematic_fp::mhfp_128(&mol);
        let idx = self.inner.add(fp);
        self.smiles.push(smiles);
        Ok(idx)
    }

    /// Search for molecules similar to ``query`` (SMILES).
    ///
    /// Args:
    ///     query: SMILES string of the query molecule.
    ///     threshold: Minimum Tanimoto similarity (0.0–1.0). Default: 0.7.
    ///     k: Maximum number of results to return. Default: all matches.
    ///
    /// Returns:
    ///     List of ``(index, similarity)`` pairs, sorted by descending similarity.
    ///
    ///     hits = idx.search("c1ccccc1N", threshold=0.5)
    ///     top_smiles = [smiles_db[i] for i, _ in hits]
    #[pyo3(signature = (query, threshold = 0.7, k = None))]
    fn search(&self, query: &str, threshold: f64, k: Option<usize>) -> PyResult<Vec<(usize, f64)>> {
        let mol =
            chematic_smiles::parse(query).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let fp = chematic_fp::mhfp_128(&mol);
        let mut results = self.inner.query(&fp, threshold);
        if let Some(max_k) = k {
            results.truncate(max_k);
        }
        Ok(results)
    }

    /// Return the SMILES of a molecule by its index in the database.
    fn get_smiles(&self, index: usize) -> PyResult<&str> {
        self.smiles
            .get(index)
            .map(|s| s.as_str())
            .ok_or_else(|| PyValueError::new_err(format!("index {index} out of range")))
    }

    /// Number of molecules in the index.
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!("SimilarityIndex(n={})", self.inner.len())
    }
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySimilarityIndex>()?;
    Ok(())
}
