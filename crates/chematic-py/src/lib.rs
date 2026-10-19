use pyo3::prelude::*;
use std::collections::HashMap;
use std::ops::{Deref, DerefMut};
use std::sync::Arc;

/// `(fingerprint bytes, {bit: [(atom_idx, radius), ...]})` for RDKit bitInfo.
type EcfpBitInfo = (
    Vec<u8>,
    std::collections::HashMap<usize, Vec<(usize, usize)>>,
);

/// `(fingerprint bytes, {raw_id: count}, {raw_id: [(atom_idx, radius), ...]}, {folded_bit:
/// [(atom_idx, radius), ...]})` for the RDKit-bit-exact Morgan/ECFP detail methods
/// (`rdkit_ecfp4_detail`, `rdkit_ecfp_config_detail`).
type RdkitMorganDetail = (
    Vec<u8>,
    std::collections::HashMap<u32, u32>,
    std::collections::HashMap<u32, Vec<(u32, u32)>>,
    std::collections::HashMap<usize, Vec<(u32, u32)>>,
);

/// `(rmsd, 4x4 transform rows, [(probe atom, reference atom)])`.
type RdkitBestAlignment = (f64, [[f64; 4]; 4], Vec<(usize, usize)>);

/// Register a list of module-level Python functions without repeating the
/// same error-propagation plumbing in every binding domain.
macro_rules! register_pyfunctions {
    ($module:expr, $($function:ident),+ $(,)?) => {{
        $(
            $module.add_function(pyo3::wrap_pyfunction!($function, $module)?)?;
        )+
        pyo3::PyResult::Ok(())
    }};
}

mod cml_named;
mod crystal;
mod descriptor_distances;
mod ensemble_v2;
mod errors;
mod feature_patterns;
mod filesystem;
mod fingerprint_similarity;
mod fold_density;
mod force_field;
mod formats;
mod lammps;
mod misc;
mod mol_methods;
mod nucleic_acid;
mod pickle;
mod pipeline_v2;
mod reactions;
mod reports;
mod rwmol;
mod similarity;
mod sparse_similarity;
mod tautomer_pick;
mod volumetric;

mod bulk;
mod index;
mod io;

// ---------------------------------------------------------------------------
// Mol — the main Python-facing molecule wrapper
// ---------------------------------------------------------------------------

/// A parsed molecule. Create with `chematic.from_smiles()` or `chematic.from_mol_block()`.
#[pyclass(name = "Mol", from_py_object)]
#[derive(Clone)]
struct Mol {
    inner: Arc<chematic_core::Molecule>,
    props: PropertyStore,
    num_rings_cache: CachedRingCount,
}

/// The graph behind `Mol.inner` is immutable. Keep its numeric ring count
/// beside the Python wrapper, avoiding another pointer chase on hot reads.
/// New graph-producing constructors start cold; `Mol::clone` shares the same
/// graph, so it can safely copy the numeric value. Property edits do not edit it.
struct CachedRingCount(std::sync::atomic::AtomicUsize);

impl CachedRingCount {
    fn new(value: usize) -> Self {
        Self(std::sync::atomic::AtomicUsize::new(value))
    }

    #[inline(always)]
    fn get_or_compute(&self, compute: impl FnOnce() -> usize) -> usize {
        use std::sync::atomic::Ordering;
        let value = self.0.load(Ordering::Relaxed);
        if value != usize::MAX {
            value
        } else {
            let value = compute();
            debug_assert_ne!(value, usize::MAX);
            // No pointer or other data is published; concurrent cold readers
            // may compute the same deterministic scalar, as in the core cache.
            self.0.store(value, Ordering::Relaxed);
            value
        }
    }
}

impl Default for CachedRingCount {
    fn default() -> Self {
        Self::new(usize::MAX)
    }
}

impl Clone for CachedRingCount {
    fn clone(&self) -> Self {
        Self::new(self.0.load(std::sync::atomic::Ordering::Relaxed))
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum PropertyType {
    String,
    Int,
    Unsigned,
    Double,
    Bool,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct PropertyStore {
    values: HashMap<String, String>,
    types: HashMap<String, PropertyType>,
}

impl PropertyStore {
    fn insert_typed(&mut self, key: String, value: String, property_type: PropertyType) {
        self.values.insert(key.clone(), value);
        if property_type == PropertyType::String {
            self.types.remove(&key);
        } else {
            self.types.insert(key, property_type);
        }
    }

    fn property_type(&self, key: &str) -> PropertyType {
        self.types.get(key).copied().unwrap_or(PropertyType::String)
    }

    fn clear(&mut self, key: &str) {
        self.values.remove(key);
        self.types.remove(key);
    }
}

impl From<HashMap<String, String>> for PropertyStore {
    fn from(values: HashMap<String, String>) -> Self {
        Self {
            values,
            types: HashMap::new(),
        }
    }
}

impl FromIterator<(String, String)> for PropertyStore {
    fn from_iter<T: IntoIterator<Item = (String, String)>>(iter: T) -> Self {
        Self::from(HashMap::from_iter(iter))
    }
}

impl<'a> IntoIterator for &'a PropertyStore {
    type Item = (&'a String, &'a String);
    type IntoIter = std::collections::hash_map::Iter<'a, String, String>;

    fn into_iter(self) -> Self::IntoIter {
        self.values.iter()
    }
}

impl Deref for PropertyStore {
    type Target = HashMap<String, String>;

    fn deref(&self) -> &Self::Target {
        &self.values
    }
}

impl DerefMut for PropertyStore {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.values
    }
}

impl Mol {
    fn bare(mol: chematic_core::Molecule) -> Self {
        Mol {
            num_rings_cache: Default::default(),
            inner: Arc::new(mol),
            props: Default::default(),
        }
    }

    /// Wrap a molecule parsed from SMILES after preparing the same two
    /// perception results RDKit materializes during `MolFromSmiles`.
    ///
    /// Keeping this at the Python boundary preserves the deliberately lean
    /// Rust SMILES parser while making Python parse/descriptor comparisons
    /// symmetric: aromaticity-dependent scalar descriptors and
    /// `CalcNumRings` no longer pay sanitization work after parsing that
    /// RDKit already charged to parsing.
    fn from_smiles_parsed(mol: chematic_core::Molecule) -> Self {
        // Keep the parser's resource contract for unusually large inputs.
        // RDKit and COSMolKit both materialize sanitized aromatic/ring state
        // while parsing, so the Python compatibility boundary does the same
        // within a fixed graph-size budget.  Ring-closure count alone is not
        // a useful complexity proxy: compact fullerenes routinely contain
        // 32 closures yet stay well inside the bounded relevant-cycle and
        // aromaticity implementations.  Every derived cache is still
        // invalidated by graph mutation.
        const MAX_EAGER_ATOMS: usize = 128;
        const MAX_EAGER_BONDS: usize = 192;
        let mut prepared_rings = usize::MAX;
        if mol.atom_count() <= MAX_EAGER_ATOMS && mol.bond_count() <= MAX_EAGER_BONDS {
            let model_prepared = chematic_smiles::prepare_rdkit_scalar_model(&mol);
            if !model_prepared && !chematic_perception::rdkit_parity_view_is_identity(&mol) {
                // Computing only the identity verdict leaves Kekule/fused
                // systems to pay the full view construction on their first
                // descriptor call.  Materialize the shared view now so
                // parse and first-use comparisons account for the same work
                // boundary as the comparator parsers.
                let _ = chematic_perception::apply_aromaticity_rdkit_parity_shared(&mol);
            }
            prepared_rings = chematic_chem::rdkit_num_rings(&mol);
        }
        Mol {
            inner: Arc::new(mol),
            props: Default::default(),
            num_rings_cache: CachedRingCount::new(prepared_rings),
        }
    }
}

/// Pure-Rust cheminformatics — SMILES, fingerprints, 70+ descriptors, pKa, ADMET.
///
/// Quick start::
///
///     import chematic
///
///     mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
///     print(mol.mw)        # 180.16
///     print(mol.logp)      # 1.31
///     print(mol.admet())   # {"bbb": False, "caco2": ..., "herg_risk": ..., ...}
///
///     fp = mol.ecfp4()     # bytes (2048-bit fingerprint)
///     print(chematic.tanimoto(fp, chematic.from_smiles("c1ccccc1").ecfp4()))
#[pymodule]
fn chematic(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Mol>()?;
    m.add_class::<force_field::PyMolecularForceField>()?;
    m.add_class::<rwmol::RWMol>()?;
    m.add_class::<pipeline_v2::PyPipelineV2Config>()?;
    m.add(
        "PipelineV2Error",
        m.py().get_type::<pipeline_v2::PipelineV2Error>(),
    )?;
    m.add(
        "ChematicInputError",
        m.py().get_type::<errors::ChematicInputError>(),
    )?;
    m.add(
        "StereoLossWarning",
        m.py().get_type::<mol_methods::StereoLossWarning>(),
    )?;
    m.add_class::<ensemble_v2::PyEnsembleV2Config>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    // bulk submodule
    let bulk_mod = PyModule::new(m.py(), "bulk")?;
    bulk::register(&bulk_mod)?;
    m.add_submodule(&bulk_mod)?;

    // io functions (iter_sdf, iter_sdf_str)
    io::register(m)?;
    pickle::register(m)?;

    // SimilarityIndex class
    index::register(m)?;

    formats::register(m)?;
    reactions::register(m)?;
    similarity::register(m)?;
    fold_density::register(m)?;
    reports::register(m)?;
    misc::register(m)?;
    nucleic_acid::register(m)?;
    crystal::register(m)?;
    volumetric::register(m)?;
    lammps::register(m)?;

    Ok(())
}

#[cfg(test)]
mod python_ring_cache_tests {
    use super::CachedRingCount;

    #[test]
    fn cold_values_are_computed_and_clones_keep_the_same_numeric_value() {
        let cache = CachedRingCount::default();
        assert_eq!(cache.get_or_compute(|| 6), 6);
        assert_eq!(cache.get_or_compute(|| panic!("cached")), 6);
        assert_eq!(cache.clone().get_or_compute(|| panic!("cloned")), 6);
        assert_eq!(CachedRingCount::default().get_or_compute(|| 0), 0);
    }

    #[test]
    fn concurrent_cold_readers_publish_only_the_numeric_value() {
        let cache = CachedRingCount::default();
        std::thread::scope(|scope| {
            for _ in 0..8 {
                let cache = &cache;
                scope.spawn(move || assert_eq!(cache.get_or_compute(|| 6), 6));
            }
        });
        assert_eq!(cache.get_or_compute(|| panic!("cached")), 6);
    }
}
