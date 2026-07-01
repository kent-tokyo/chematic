use pyo3::prelude::*;
use std::sync::Arc;

/// `(fingerprint bytes, {bit: [(atom_idx, radius), ...]})` for RDKit bitInfo.
type EcfpBitInfo = (
    Vec<u8>,
    std::collections::HashMap<usize, Vec<(usize, usize)>>,
);

mod formats;
mod misc;
mod mol_methods;
mod reactions;
mod reports;
mod similarity;

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
    props: std::collections::HashMap<String, String>,
}

impl Mol {
    fn bare(mol: chematic_core::Molecule) -> Self {
        Mol {
            inner: Arc::new(mol),
            props: Default::default(),
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
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    // bulk submodule
    let bulk_mod = PyModule::new(m.py(), "bulk")?;
    bulk::register(&bulk_mod)?;
    m.add_submodule(&bulk_mod)?;

    // io functions (iter_sdf, iter_sdf_str)
    io::register(m)?;

    // SimilarityIndex class
    index::register(m)?;

    formats::register(m)?;
    reactions::register(m)?;
    similarity::register(m)?;
    reports::register(m)?;
    misc::register(m)?;

    Ok(())
}
