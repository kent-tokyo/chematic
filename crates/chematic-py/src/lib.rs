use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyDict;
use numpy::{IntoPyArray, PyArray1};
use ndarray::Array1;
use std::sync::Arc;

mod bulk;
mod io;
mod index;

// ---------------------------------------------------------------------------
// Mol — the main Python-facing molecule wrapper
// ---------------------------------------------------------------------------

/// A parsed molecule. Create with `chematic.from_smiles()` or `chematic.from_mol_block()`.
#[pyclass(name = "Mol", from_py_object)]
#[derive(Clone)]
struct Mol {
    inner: Arc<chematic_core::Molecule>,
}

#[pymethods]
impl Mol {
    // -----------------------------------------------------------------------
    // Identity / structure
    // -----------------------------------------------------------------------

    /// Canonical SMILES string.
    #[getter]
    fn smiles(&self) -> String {
        chematic_smiles::canonical_smiles(&self.inner)
    }

    /// Molecular formula in Hill notation (C first, H second, then alphabetical).
    #[getter]
    fn formula(&self) -> String {
        chematic_chem::calc_mol_formula(&self.inner)
    }

    /// Number of heavy atoms (explicit atoms; does not count implicit H).
    #[getter]
    fn heavy_atoms(&self) -> usize {
        self.inner.atom_count()
    }

    /// Non-standard InChI string (pure-Rust approximation, not IUPAC-compliant).
    ///
    /// For standard IUPAC InChI, use :attr:`standard_inchi` (requires the
    /// ``native-inchi`` feature at build time).
    #[getter]
    fn inchi(&self) -> String {
        chematic_inchi::inchi(&self.inner)
    }

    /// Non-standard InChIKey (pure-Rust approximation, not IUPAC-compliant).
    ///
    /// For standard IUPAC InChIKey, use :attr:`standard_inchikey`.
    #[getter]
    fn inchikey(&self) -> String {
        let s = chematic_inchi::inchi(&self.inner);
        chematic_inchi::inchi_key(&s)
    }

    /// Standard IUPAC InChI string via the vendored InChI C library (v1.07.5).
    ///
    /// Returns a :exc:`RuntimeError` if generation fails or if the ``native-inchi``
    /// Cargo feature was not enabled at build time.
    #[cfg(feature = "native-inchi")]
    #[getter]
    fn standard_inchi(&self) -> pyo3::PyResult<String> {
        chematic_inchi::standard_inchi(&self.inner)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    #[cfg(not(feature = "native-inchi"))]
    #[getter]
    fn standard_inchi(&self) -> pyo3::PyResult<String> {
        Err(pyo3::exceptions::PyRuntimeError::new_err(
            "standard_inchi requires the native-inchi feature \
             (rebuild chematic from source with --features native-inchi)",
        ))
    }

    /// Standard IUPAC InChIKey (27 characters) via the vendored InChI C library.
    ///
    /// Returns a :exc:`RuntimeError` if generation fails or if the ``native-inchi``
    /// Cargo feature was not enabled at build time.
    #[cfg(feature = "native-inchi")]
    #[getter]
    fn standard_inchikey(&self) -> pyo3::PyResult<String> {
        let inchi = chematic_inchi::standard_inchi(&self.inner)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        chematic_inchi::standard_inchi_key(&inchi)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    #[cfg(not(feature = "native-inchi"))]
    #[getter]
    fn standard_inchikey(&self) -> pyo3::PyResult<String> {
        Err(pyo3::exceptions::PyRuntimeError::new_err(
            "standard_inchikey requires the native-inchi feature \
             (rebuild chematic from source with --features native-inchi)",
        ))
    }

    /// IUPAC systematic name. Returns an empty string for unsupported structures.
    #[getter]
    fn iupac_name(&self) -> String {
        chematic_iupac::name(&self.inner).unwrap_or_default()
    }

    /// Serialize to Tripos MOL2 format string.
    ///
    /// Example::
    ///
    ///     mol = chematic.from_smiles("CCO")
    ///     with open("ethanol.mol2", "w") as f:
    ///         f.write(mol.to_mol2())
    fn to_mol2(&self) -> String {
        chematic_mol::write_mol2(&self.inner, &[])
    }

    // -----------------------------------------------------------------------
    // Core physicochemical descriptors
    // -----------------------------------------------------------------------

    /// Average molecular weight (Da).
    #[getter]
    fn mw(&self) -> f64 {
        chematic_chem::molecular_weight(&self.inner)
    }

    /// Monoisotopic (exact) mass.
    #[getter]
    fn exact_mass(&self) -> f64 {
        chematic_chem::exact_mass(&self.inner)
    }

    /// Crippen–Wildman LogP.
    #[getter]
    fn logp(&self) -> f64 {
        chematic_chem::logp_crippen(&self.inner)
    }

    /// Topological polar surface area (Å²).
    #[getter]
    fn tpsa(&self) -> f64 {
        chematic_chem::tpsa(&self.inner)
    }

    /// Quantitative Estimate of Drug-likeness [0, 1].
    #[getter]
    fn qed(&self) -> f64 {
        chematic_chem::qed(&self.inner)
    }

    /// Hydrogen bond donors.
    #[getter]
    fn hbd(&self) -> usize {
        chematic_chem::hbd_count(&self.inner)
    }

    /// Hydrogen bond acceptors.
    #[getter]
    fn hba(&self) -> usize {
        chematic_chem::hba_count(&self.inner)
    }

    /// Number of rotatable bonds.
    #[getter]
    fn rotatable_bonds(&self) -> usize {
        chematic_chem::rotatable_bond_count(&self.inner)
    }

    /// Fraction of sp3 carbons (Fsp3).
    #[getter]
    fn fsp3(&self) -> f64 {
        chematic_chem::fsp3(&self.inner)
    }

    /// Synthetic Accessibility Score [1–10]; lower = easier to synthesize.
    #[getter]
    fn sa_score(&self) -> f64 {
        chematic_chem::sa_score(&self.inner)
    }

    /// Wildman–Crippen molar refractivity.
    #[getter]
    fn molar_refractivity(&self) -> f64 {
        chematic_chem::molar_refractivity(&self.inner)
    }

    /// Sum of formal charges.
    #[getter]
    fn formal_charge(&self) -> i32 {
        chematic_chem::formal_charge_sum(&self.inner)
    }

    // -----------------------------------------------------------------------
    // Ring / stereo counts
    // -----------------------------------------------------------------------

    /// Total number of rings (SSSR count).
    #[getter]
    fn ring_count(&self) -> usize {
        chematic_chem::ring_count(&self.inner)
    }

    /// Number of aromatic rings.
    #[getter]
    fn aromatic_ring_count(&self) -> usize {
        chematic_chem::aromatic_ring_count(&self.inner)
    }

    /// Number of assigned stereocenters (R/S).
    #[getter]
    fn num_stereocenters(&self) -> usize {
        chematic_chem::num_stereocenters(&self.inner)
    }

    // -----------------------------------------------------------------------
    // Drug-likeness rules
    // -----------------------------------------------------------------------

    /// True if Lipinski's Rule of Five passes.
    #[getter]
    fn lipinski_passes(&self) -> bool {
        chematic_chem::lipinski_passes(&self.inner)
    }

    /// True if Veber's oral bioavailability criteria pass.
    #[getter]
    fn veber_passes(&self) -> bool {
        chematic_chem::veber_passes(&self.inner)
    }

    /// True if no PAINS structural alerts are present.
    #[getter]
    fn pains_passes(&self) -> bool {
        chematic_chem::pains_passes(&self.inner)
    }

    /// True if Ghose drug-likeness criteria pass (MW 160–480, LogP −0.4–5.6, MR 40–130, atoms 20–70).
    #[getter]
    fn ghose_passes(&self) -> bool {
        chematic_chem::ghose_passes(&self.inner)
    }

    /// True if Egan absorption criteria pass (TPSA ≤ 131.6, LogP ≤ 5.88).
    #[getter]
    fn egan_passes(&self) -> bool {
        chematic_chem::egan_passes(&self.inner)
    }

    /// True if REOS drug-likeness filter passes.
    #[getter]
    fn reos_passes(&self) -> bool {
        chematic_chem::reos_passes(&self.inner)
    }

    /// True if no Brenk structural alerts are present.
    #[getter]
    fn brenk_passes(&self) -> bool {
        chematic_chem::brenk_passes(&self.inner)
    }

    // -----------------------------------------------------------------------
    // pKa and ADMET
    // -----------------------------------------------------------------------

    /// pKa prediction.
    ///
    /// Returns a dict with keys ``most_acidic`` and ``most_basic``
    /// (float or ``None`` when no such site is found).
    fn pka<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        let m = &self.inner;
        match chematic_chem::pka_acid(m) {
            Some(v) => d.set_item("most_acidic", v)?,
            None => d.set_item("most_acidic", py.None())?,
        }
        match chematic_chem::pka_base(m) {
            Some(v) => d.set_item("most_basic", v)?,
            None => d.set_item("most_basic", py.None())?,
        }
        Ok(d)
    }

    /// ADMET profile.
    ///
    /// Returns a dict with keys:
    /// ``bbb`` (bool), ``bbb_score`` (float),
    /// ``caco2`` (float, logPCaco2),
    /// ``herg_risk`` (float, 0–1),
    /// ``cyp3a4_risk`` (float, 0–1),
    /// ``ames_risk`` (float, 0–1),
    /// ``ppb`` (float, plasma protein binding %),
    /// ``clearance`` (str: ``"Low"`` / ``"Medium"`` / ``"High"``),
    /// ``gi_absorbed`` (bool), ``bbb_penetrant`` (bool).
    fn admet<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        let m = &self.inner;
        d.set_item("bbb", chematic_chem::bbb_passes(m))?;
        d.set_item("bbb_score", chematic_chem::bbb_score(m))?;
        d.set_item("caco2", chematic_chem::caco2_permeability(m))?;
        d.set_item("herg_risk", chematic_chem::herg_risk_score(m))?;
        d.set_item("cyp3a4_risk", chematic_chem::cyp3a4_inhibition_risk(m))?;
        d.set_item("ames_risk", chematic_chem::ames_risk_score(m))?;
        d.set_item("ppb", chematic_chem::ppb_percent(m))?;
        d.set_item("clearance", match chematic_chem::clearance_class(m) {
            chematic_chem::ClearanceClass::Low    => "Low",
            chematic_chem::ClearanceClass::Medium => "Medium",
            chematic_chem::ClearanceClass::High   => "High",
        })?;
        let egg = chematic_chem::boiled_egg(m);
        d.set_item("gi_absorbed", egg.gi_absorbed)?;
        d.set_item("bbb_penetrant", egg.bbb_penetrant)?;
        Ok(d)
    }

    /// Predict GI absorption and BBB penetration using the BOILED-Egg method
    /// (Daina & Zoete 2016).
    ///
    /// Returns a dict with keys:
    /// ``gi_absorbed`` (bool), ``bbb_penetrant`` (bool),
    /// ``logp`` (float), ``tpsa`` (float).
    ///
    /// Example::
    ///
    ///     egg = mol.boiled_egg()
    ///     print(egg["gi_absorbed"])    # True / False
    ///     print(egg["bbb_penetrant"])  # True / False
    fn boiled_egg<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let e = chematic_chem::boiled_egg(&self.inner);
        let d = PyDict::new(py);
        d.set_item("gi_absorbed", e.gi_absorbed)?;
        d.set_item("bbb_penetrant", e.bbb_penetrant)?;
        d.set_item("logp", e.logp)?;
        d.set_item("tpsa", e.tpsa)?;
        Ok(d)
    }

    // -----------------------------------------------------------------------
    // All descriptors in one call
    // -----------------------------------------------------------------------

    /// Return all scalar descriptors as a dict (70+ keys).
    ///
    /// Useful for building a Pandas DataFrame row::
    ///
    ///     import pandas as pd
    ///     import chematic
    ///     smiles = ["CCO", "c1ccccc1", "CC(=O)O"]
    ///     df = pd.DataFrame([chematic.from_smiles(s).descriptors() for s in smiles])
    fn descriptors<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        let m = &self.inner;
        d.set_item("mw", chematic_chem::molecular_weight(m))?;
        d.set_item("exact_mass", chematic_chem::exact_mass(m))?;
        d.set_item("tpsa", chematic_chem::tpsa(m))?;
        d.set_item("logp", chematic_chem::logp_crippen(m))?;
        d.set_item("molar_refractivity", chematic_chem::molar_refractivity(m))?;
        d.set_item("hbd", chematic_chem::hbd_count(m))?;
        d.set_item("hba", chematic_chem::hba_count(m))?;
        d.set_item("rotatable_bonds", chematic_chem::rotatable_bond_count(m))?;
        d.set_item("heavy_atoms", chematic_chem::heavy_atom_count(m))?;
        d.set_item("ring_count", chematic_chem::ring_count(m))?;
        d.set_item("aromatic_ring_count", chematic_chem::aromatic_ring_count(m))?;
        d.set_item("num_heteroatoms", chematic_chem::num_heteroatoms(m))?;
        d.set_item("num_stereocenters", chematic_chem::num_stereocenters(m))?;
        d.set_item("num_spiro_atoms", chematic_chem::num_spiro_atoms(m))?;
        d.set_item("num_bridgehead_atoms", chematic_chem::num_bridgehead_atoms(m))?;
        d.set_item("fsp3", chematic_chem::fsp3(m))?;
        d.set_item("qed", chematic_chem::qed(m))?;
        d.set_item("sa_score", chematic_chem::sa_score(m))?;
        d.set_item("formal_charge", chematic_chem::formal_charge_sum(m))?;
        d.set_item("labute_asa", chematic_chem::labute_asa(m))?;
        d.set_item("bertz_ct", chematic_chem::bertz_ct(m))?;
        d.set_item("wiener_index", chematic_chem::wiener_index(m))?;
        d.set_item("kappa1", chematic_chem::kappa1(m))?;
        d.set_item("kappa2", chematic_chem::kappa2(m))?;
        d.set_item("kappa3", chematic_chem::kappa3(m))?;
        d.set_item("chi0", chematic_chem::chi0(m))?;
        d.set_item("chi1", chematic_chem::chi1(m))?;
        d.set_item("chi2", chematic_chem::chi2(m))?;
        d.set_item("chi3", chematic_chem::chi3(m))?;
        d.set_item("chi4", chematic_chem::chi4(m))?;
        d.set_item("chi0v", chematic_chem::chi0v(m))?;
        d.set_item("chi1v", chematic_chem::chi1v(m))?;
        d.set_item("chi2v", chematic_chem::chi2v(m))?;
        d.set_item("chi3v", chematic_chem::chi3v(m))?;
        d.set_item("chi4v", chematic_chem::chi4v(m))?;
        d.set_item("num_aromatic_heterocycles", chematic_chem::num_aromatic_heterocycles(m))?;
        d.set_item("num_aliphatic_heterocycles", chematic_chem::num_aliphatic_heterocycles(m))?;
        d.set_item("num_saturated_rings", chematic_chem::num_saturated_rings(m))?;
        d.set_item("num_aliphatic_rings", chematic_chem::num_aliphatic_rings(m))?;
        d.set_item("num_unspecified_stereocenters", chematic_chem::num_unspecified_stereocenters(m))?;
        d.set_item("sum_estate", chematic_chem::sum_estate(m))?;
        d.set_item("max_estate", chematic_chem::max_estate(m))?;
        d.set_item("min_estate", chematic_chem::min_estate(m))?;
        d.set_item("lipinski_passes", chematic_chem::lipinski_passes(m))?;
        d.set_item("veber_passes", chematic_chem::veber_passes(m))?;
        d.set_item("egan_passes", chematic_chem::egan_passes(m))?;
        d.set_item("ghose_passes", chematic_chem::ghose_passes(m))?;
        d.set_item("reos_passes", chematic_chem::reos_passes(m))?;
        d.set_item("pains_passes", chematic_chem::pains_passes(m))?;
        d.set_item("bbb_score", chematic_chem::bbb_score(m))?;
        d.set_item("bbb_passes", chematic_chem::bbb_passes(m))?;
        d.set_item("caco2", chematic_chem::caco2_permeability(m))?;
        d.set_item("herg_risk", chematic_chem::herg_risk_score(m))?;
        d.set_item("cyp3a4_risk", chematic_chem::cyp3a4_inhibition_risk(m))?;
        let egg = chematic_chem::boiled_egg(m);
        d.set_item("gi_absorbed", egg.gi_absorbed)?;
        d.set_item("bbb_penetrant", egg.bbb_penetrant)?;
        match chematic_chem::pka_acid(m) {
            Some(v) => d.set_item("pka_acid", v)?,
            None => d.set_item("pka_acid", py.None())?,
        }
        match chematic_chem::pka_base(m) {
            Some(v) => d.set_item("pka_base", v)?,
            None => d.set_item("pka_base", py.None())?,
        }
        Ok(d)
    }

    // -----------------------------------------------------------------------
    // Fingerprints
    // -----------------------------------------------------------------------

    /// ECFP4 fingerprint as bytes (256 bytes = 2048 bits, LSB-first).
    fn ecfp4(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::ecfp4(&self.inner))
    }

    /// ECFP6 fingerprint as bytes (256 bytes = 2048 bits, LSB-first).
    fn ecfp6(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::ecfp6(&self.inner))
    }

    /// FCFP4 (functional-class ECFP4) fingerprint as bytes (256 bytes = 2048 bits, LSB-first).
    fn fcfp4(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::fcfp4(&self.inner))
    }

    /// Atom-pair fingerprint as bytes (256 bytes = 2048 bits, LSB-first).
    fn atom_pair_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::atom_pair_fp(&self.inner))
    }

    /// Topological torsion fingerprint as bytes (256 bytes = 2048 bits, LSB-first).
    fn torsion_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::torsion_fp(&self.inner))
    }

    /// ECFP4 with chirality fingerprint as bytes.
    fn ecfp4_chiral(&self) -> Vec<u8> {
        let config = chematic_fp::EcfpConfig {
            radius: 2,
            nbits: 2048,
            use_chirality: true,
            use_double_fold: false,
        };
        bitvec2048_to_bytes(&chematic_fp::ecfp(&self.inner, &config))
    }

    /// MACCS 166-bit keys as bytes (21 bytes, LSB-first).
    fn maccs(&self) -> Vec<u8> {
        let fp = chematic_fp::maccs(&self.inner);
        (0..21usize)
            .map(|byte_idx| {
                let mut byte = 0u8;
                for bit in 0..8usize {
                    let global_bit = byte_idx * 8 + bit;
                    if global_bit < 166 && fp.get(global_bit) {
                        byte |= 1 << bit;
                    }
                }
                byte
            })
            .collect()
    }

    // -----------------------------------------------------------------------
    // numpy fingerprint outputs
    // -----------------------------------------------------------------------

    /// ECFP4 fingerprint as a numpy array of shape ``(2048,)`` with ``dtype=uint8`` (0/1).
    ///
    /// Useful for direct use with scikit-learn or PyTorch::
    ///
    ///     fp = mol.ecfp4_numpy()   # shape (2048,), dtype uint8
    ///     X = np.stack([m.ecfp4_numpy() for m in mols])  # (N, 2048)
    fn ecfp4_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u8>> {
        let fp = chematic_fp::ecfp4(&self.inner);
        let bits: Vec<u8> = (0..2048usize).map(|i| u8::from(fp.get(i))).collect();
        Array1::from_vec(bits).into_pyarray(py)
    }

    /// MACCS 166-bit fingerprint as a numpy array of shape ``(166,)`` with ``dtype=uint8`` (0/1).
    fn maccs_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u8>> {
        let fp = chematic_fp::maccs(&self.inner);
        let bits: Vec<u8> = (0..166usize).map(|i| u8::from(fp.get(i))).collect();
        Array1::from_vec(bits).into_pyarray(py)
    }

    // -----------------------------------------------------------------------
    // Visualization
    // -----------------------------------------------------------------------

    /// 2D SVG depiction as a string.
    ///
    /// In Jupyter notebooks this renders automatically::
    ///
    ///     from IPython.display import SVG
    ///     SVG(mol.svg())
    fn svg(&self) -> String {
        chematic_depict::depict_svg(&self.inner)
    }

    /// 2D SVG depiction with highlighted atoms.
    ///
    /// ``atom_indices``: zero-based atom indices to highlight.
    /// ``color``: CSS color string (default ``"#FFFF00"`` yellow).
    ///
    ///     svg = mol.svg_highlighted([0, 1, 2])
    ///     svg = mol.svg_highlighted([0], color="#FF0000")
    #[pyo3(signature = (atom_indices, color = "#FFFF00"))]
    fn svg_highlighted(&self, atom_indices: Vec<usize>, color: &str) -> String {
        use chematic_core::AtomIdx;
        let mut opts = chematic_depict::RenderOptions {
            highlight_color: color.to_string(),
            ..Default::default()
        };
        for i in atom_indices {
            opts.highlight_atoms.insert(AtomIdx(i as u32));
        }
        chematic_depict::depict_svg_opts(&self.inner, &opts)
    }

    // -----------------------------------------------------------------------
    // ADMET / property predictions
    // -----------------------------------------------------------------------

    /// ESOL estimated aqueous solubility (log mol/L).
    ///
    /// Based on Delaney (2004). Negative values = less soluble.
    /// Typical range: −6 (insoluble) to 0 (freely soluble).
    #[getter]
    fn esol(&self) -> f64 {
        chematic_chem::esol_solubility(&self.inner)
    }

    // -----------------------------------------------------------------------
    // Transformations
    // -----------------------------------------------------------------------

    /// Return the standardized molecule (largest fragment, charges neutralized,
    /// tautomer canonicalized, isotopes/stereo preserved by default).
    fn standardize(&self) -> Mol {
        let opts = chematic_chem::StandardizeOptions::default();
        Mol {
            inner: Arc::new(chematic_chem::standardize(&self.inner, &opts)),
        }
    }

    /// Return the Murcko scaffold as a new Mol.
    fn scaffold(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::murcko_scaffold(&self.inner)),
        }
    }

    /// Return the canonical tautomer as a new Mol.
    fn canonical_tautomer(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::canonical_tautomer(&self.inner)),
        }
    }

    /// Return all tautomers as a list of Mol objects.
    fn enumerate_tautomers(&self) -> Vec<Mol> {
        chematic_chem::enumerate_tautomers(&self.inner)
            .into_iter()
            .map(|m| Mol { inner: Arc::new(m) })
            .collect()
    }

    /// Return all stereoisomers as a list of Mol objects.
    ///
    /// Returns an empty list when the molecule has more than 6 unspecified
    /// stereocenters (combinatorial explosion guard). Use
    /// ``mol.num_stereocenters`` to distinguish this from having no centers.
    fn enumerate_stereoisomers(&self) -> Vec<Mol> {
        chematic_chem::enumerate_stereoisomers(&self.inner)
            .into_iter()
            .map(|m| Mol { inner: Arc::new(m) })
            .collect()
    }

    /// Return a copy with all implicit hydrogens made explicit.
    fn add_hydrogens(&self) -> Mol {
        Mol { inner: Arc::new(chematic_chem::add_hydrogens(&self.inner)) }
    }

    /// Return a copy with all explicit hydrogen atoms removed.
    fn remove_hydrogens(&self) -> Mol {
        Mol { inner: Arc::new(chematic_chem::remove_hydrogens(&self.inner)) }
    }

    /// Return a copy with all stereochemistry assignments removed.
    fn remove_stereo(&self) -> Mol {
        Mol { inner: Arc::new(chematic_chem::remove_stereo(&self.inner)) }
    }

    /// Return a copy with all isotope labels removed.
    fn remove_isotopes(&self) -> Mol {
        Mol { inner: Arc::new(chematic_chem::remove_isotopes(&self.inner)) }
    }

    /// Return the largest covalently connected fragment.
    fn largest_fragment(&self) -> Mol {
        Mol { inner: Arc::new(chematic_chem::largest_fragment(&self.inner)) }
    }

    /// Return a charge-neutralized copy.
    fn neutralize(&self) -> Mol {
        Mol { inner: Arc::new(chematic_chem::neutralize_charges(&self.inner)) }
    }

    /// Return the generic Murcko scaffold (all atoms replaced with carbons, all bonds single).
    fn generic_scaffold(&self) -> Mol {
        Mol { inner: Arc::new(chematic_chem::generic_murcko_scaffold(&self.inner)) }
    }

    /// Fragment the molecule using BRICS rules. Returns a list of fragment Mol objects.
    ///
    /// When no BRICS-breakable bonds are found, returns a list containing the
    /// original molecule (not an empty list).
    fn brics_fragments(&self) -> Vec<Mol> {
        chematic_chem::brics_fragments(&self.inner)
            .into_iter()
            .map(|m| Mol { inner: Arc::new(m) })
            .collect()
    }

    // -----------------------------------------------------------------------
    // B5: Layered fingerprint
    // -----------------------------------------------------------------------

    /// Layered fingerprint as bytes (256 bytes = 2048 bits).
    ///
    /// 7-layer structural fingerprint encoding atom type, aromaticity, ring
    /// membership, charge, and connectivity. RDKit-compatible layer design.
    fn layered_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::layered_fp(&self.inner))
    }

    /// Layered fingerprint as a numpy array of shape ``(2048,)`` with ``dtype=uint8`` (0/1).
    fn layered_fp_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u8>> {
        let fp = chematic_fp::layered_fp(&self.inner);
        let bits: Vec<u8> = (0..2048usize).map(|i| u8::from(fp.get(i))).collect();
        Array1::from_vec(bits).into_pyarray(py)
    }

    // -----------------------------------------------------------------------
    // B8: 3D SASA
    // -----------------------------------------------------------------------

    /// Solvent-Accessible Surface Area (Å²) via Shrake-Rupley algorithm.
    ///
    /// Uses distance-geometry 3D coordinates (no prior minimization).
    /// probe_radius = 1.4 Å (water), sphere_points = 100.
    fn sasa(&self) -> f64 {
        chematic_3d::sasa_from_dg(&self.inner).unwrap_or(0.0)
    }

    /// Per-atom SASA as a list of float values (Å²).
    fn sasa_per_atom(&self) -> Vec<f64> {
        chematic_3d::sasa_per_atom_from_dg(&self.inner).unwrap_or_default()
    }

    // -----------------------------------------------------------------------
    // Task 1: 3D shape similarity (USR)
    // -----------------------------------------------------------------------

    /// 12 USR shape descriptors computed from auto-generated 3D coordinates.
    ///
    /// Based on Ballester & Richards (2007) Ultrafast Shape Recognition.
    /// Returns a list of 12 floats (3 statistical moments × 4 reference points).
    fn usr_descriptors(&self) -> Vec<f64> {
        chematic_3d::usr_from_dg(&self.inner).to_vec()
    }

    /// USR 3D shape similarity to another molecule (0.0–1.0).
    ///
    /// 1.0 = identical shape, approaches 0.0 for very different shapes.
    fn usr_similarity(&self, other: &Mol) -> f64 {
        let a = chematic_3d::usr_from_dg(&self.inner);
        let b = chematic_3d::usr_from_dg(&other.inner);
        chematic_3d::usr_similarity(&a, &b)
    }

    // -----------------------------------------------------------------------
    // Task 2: Additional ADMET — Ames, PPB, clearance
    // -----------------------------------------------------------------------

    /// Ames mutagenicity risk score (0.0–1.0).
    ///
    /// Based on SMARTS structural alerts (Kazius 2005, simplified).
    /// Score > 0 indicates at least one alert pattern was matched.
    fn ames_risk(&self) -> f64 {
        chematic_chem::ames_risk_score(&self.inner)
    }

    /// Returns `True` if no Ames structural alerts are present.
    fn ames_passes(&self) -> bool {
        chematic_chem::ames_passes(&self.inner)
    }

    /// Predicted plasma protein binding (%).
    ///
    /// Logistic model based on LogP (Arnott 2012 heuristic).
    /// Interpretation: > 90% = highly bound, < 20% = low binding.
    fn ppb(&self) -> f64 {
        chematic_chem::ppb_percent(&self.inner)
    }

    /// Predicted hepatic clearance class: ``"Low"``, ``"Medium"``, or ``"High"``.
    fn clearance(&self) -> &'static str {
        match chematic_chem::clearance_class(&self.inner) {
            chematic_chem::ClearanceClass::Low    => "Low",
            chematic_chem::ClearanceClass::Medium => "Medium",
            chematic_chem::ClearanceClass::High   => "High",
        }
    }

    // -----------------------------------------------------------------------
    // C1: Atropisomer detection
    // -----------------------------------------------------------------------

    /// Detect atropisomeric axes in the molecule.
    ///
    /// Returns a list of ``(bond_idx, kind)`` tuples where ``kind`` is
    /// ``"Biaryl"`` or ``"Allene"``.
    fn atropisomers(&self) -> Vec<(usize, String)> {
        chematic_chem::detect_atropisomers(&self.inner)
            .into_iter()
            .map(|(bidx, kind)| {
                let label = match kind {
                    chematic_chem::AtropisomerType::Biaryl      => "Biaryl",
                    chematic_chem::AtropisomerType::Allene      => "Allene",
                    chematic_chem::AtropisomerType::Constrained => "Constrained",
                };
                (bidx.0 as usize, label.to_string())
            })
            .collect()
    }

    // -----------------------------------------------------------------------
    // Dunder methods
    // -----------------------------------------------------------------------

    fn __repr__(&self) -> String {
        format!("Mol('{}')", chematic_smiles::canonical_smiles(&self.inner))
    }

    fn __str__(&self) -> String {
        chematic_smiles::canonical_smiles(&self.inner)
    }
}

// ---------------------------------------------------------------------------
// Helper: convert a 2048-bit fingerprint to a 256-byte Vec<u8>
// ---------------------------------------------------------------------------

fn bitvec2048_to_bytes(fp: &chematic_fp::bitvec::BitVec2048) -> Vec<u8> {
    (0..256usize)
        .map(|byte_idx| {
            let mut byte = 0u8;
            for bit in 0..8usize {
                if fp.get(byte_idx * 8 + bit) {
                    byte |= 1 << bit;
                }
            }
            byte
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Module-level functions
// ---------------------------------------------------------------------------

/// Parse a SMILES string and return a Mol.
///
/// Raises ``ValueError`` on invalid SMILES.
#[pyfunction]
fn from_smiles(smiles: &str) -> PyResult<Mol> {
    chematic_smiles::parse(smiles)
        .map(|mol| Mol { inner: Arc::new(mol) })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a MOL/SDF block and return a Mol.
///
/// Raises ``ValueError`` on parse failure.
#[pyfunction]
fn from_mol_block(block: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol(block)
        .map(|(mol, _meta)| Mol { inner: Arc::new(mol) })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Tripos MOL2 string into a ``Mol`` object.
///
/// Example::
///
///     with open("ligand.mol2") as f:
///         mol = chematic.from_mol2(f.read())
///     print(mol.mw)
///
/// Raises ``ValueError`` on parse failure.
#[pyfunction]
fn from_mol2(mol2_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol2(mol2_str)
        .map(|(mol, _coords)| Mol { inner: Arc::new(mol) })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Return True if the SMILES can be parsed without error.
#[pyfunction]
fn is_valid_smiles(smiles: &str) -> bool {
    chematic_smiles::parse(smiles).is_ok()
}

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
    let and_bits: u32 = a.iter().zip(b.iter()).map(|(x, y)| (x & y).count_ones()).sum();
    let or_bits: u32 = a.iter().zip(b.iter()).map(|(x, y)| (x | y).count_ones()).sum();
    if or_bits == 0 {
        return Ok(0.0);
    }
    Ok(and_bits as f64 / or_bits as f64)
}

/// Test whether a SMARTS pattern matches a molecule.
///
///     if chematic.smarts_match("[OH]", mol):
///         print("has hydroxyl")
#[pyfunction]
fn smarts_match(smarts: &str, mol: &Mol) -> PyResult<bool> {
    let query = chematic_smarts::parse_smarts(smarts)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(!chematic_smarts::find_matches(&query, &mol.inner).is_empty())
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
    let query = chematic_smarts::parse_smarts(smarts)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
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

/// Parse an InChI string and return a Mol.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol = chematic.from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
#[pyfunction]
fn from_inchi(inchi: &str) -> PyResult<Mol> {
    chematic_inchi::parse_inchi(inchi)
        .map(|mol| Mol { inner: Arc::new(mol) })
        .map_err(|e| PyValueError::new_err(e.to_string()))
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

/// Test whether a reaction SMARTS pattern matches a reaction SMILES.
///
/// ``smarts``: reaction SMARTS (e.g. ``"[OH:1]>>[O-:1]"``).
/// ``reaction_smiles``: reaction SMILES in ``"R>>P"`` or ``"R>A>P"`` format.
///
///     ok = chematic.reaction_smarts_match("[OH]>>[O-]", "CCO>>CC[O-]")
#[pyfunction]
fn reaction_smarts_match(smarts: &str, reaction_smiles: &str) -> PyResult<bool> {
    let query = chematic_rxn::parse_reaction_query(smarts)
        .map_err(|e| PyValueError::new_err(format!("invalid reaction SMARTS: {e}")))?;
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(format!("invalid reaction SMILES: {e}")))?;
    Ok(chematic_rxn::has_reaction_substructure_match(&rxn, &query))
}

/// Render a list of molecules as a grid SVG.
///
///     svg = chematic.depict_grid([mol1, mol2, mol3], cols=3)
#[pyfunction]
fn depict_grid(mols: Vec<Mol>, cols: usize) -> String {
    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    chematic_depict::depict_svg_grid(&refs, cols)
}

/// Apply a SMIRKS reaction template to a list of reactant molecules.
///
/// Returns a list of product sets; each set is a list of Mol.
/// Raises ``ValueError`` on SMIRKS parse failure or reactant count mismatch.
///
///     products = chematic.run_smirks("[OH:1]>>[O-:1]", [mol])
///     # → [[product_mol], ...]
#[pyfunction]
fn run_smirks(smirks: &str, reactants: Vec<Mol>) -> PyResult<Vec<Vec<Mol>>> {
    let refs: Vec<&chematic_core::Molecule> = reactants.iter().map(|m| m.inner.as_ref()).collect();
    chematic_rxn::run_reactants(smirks, &refs)
        .map(|sets| {
            sets.into_iter()
                .map(|set| {
                    set.into_iter()
                        .map(|m| Mol { inner: Arc::new(m) })
                        .collect()
                })
                .collect()
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Find the Maximum Common Substructure (MCS) of a list of molecules.
///
/// Returns the MCS as a Mol, or ``None`` when there is no common substructure.
///
///     mcs = chematic.find_mcs([mol1, mol2])
///     if mcs: print(mcs.smiles)
///
///     # Ring-aware scaffold extraction
///     scaffold = chematic.find_mcs(mols, ring_matches_ring_only=True, complete_rings_only=True)
#[pyfunction]
#[pyo3(signature = (mols, ring_matches_ring_only=false, complete_rings_only=false))]
fn find_mcs(mols: Vec<Mol>, ring_matches_ring_only: bool, complete_rings_only: bool) -> Option<Mol> {
    use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
    use chematic_smarts::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, McsConfig};

    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    let config = McsConfig { ring_matches_ring_only, complete_rings_only, ..McsConfig::default() };
    let qmol = chematic_smarts::find_mcs_with_config(&refs, &config);

    if qmol.atoms.is_empty() {
        return None;
    }

    fn extract_atomic_num(q: &AtomQuery) -> Option<u8> {
        match q {
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)) => Some(*n),
            AtomQuery::And(lhs, rhs) => {
                extract_atomic_num(lhs).or_else(|| extract_atomic_num(rhs))
            }
            _ => None,
        }
    }

    let mut builder = MoleculeBuilder::new();
    for qa in &qmol.atoms {
        let elem = extract_atomic_num(&qa.query)
            .and_then(Element::from_atomic_number)
            .unwrap_or(Element::C);
        builder.add_atom(Atom::new(elem));
    }
    for (atom_idx, neighbors) in qmol.adj.iter().enumerate() {
        for (bond_idx, neighbor_idx) in neighbors {
            if atom_idx < *neighbor_idx {
                let order = match &qmol.bonds[*bond_idx].query {
                    BondQuery::Primitive(BondPrimitive::Double) => BondOrder::Double,
                    BondQuery::Primitive(BondPrimitive::Triple) => BondOrder::Triple,
                    BondQuery::Primitive(BondPrimitive::Aromatic) => BondOrder::Aromatic,
                    _ => BondOrder::Single,
                };
                let _ = builder.add_bond(
                    AtomIdx(atom_idx as u32),
                    AtomIdx(*neighbor_idx as u32),
                    order,
                );
            }
        }
    }
    Some(Mol { inner: Arc::new(builder.build()) })
}

// ---------------------------------------------------------------------------
// Module definition
// ---------------------------------------------------------------------------

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
    m.add_function(wrap_pyfunction!(from_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_block, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol2, m)?)?;
    m.add_function(wrap_pyfunction!(from_inchi, m)?)?;
    m.add_function(wrap_pyfunction!(is_valid_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto, m)?)?;
    m.add_function(wrap_pyfunction!(smarts_match, m)?)?;
    m.add_function(wrap_pyfunction!(smarts_find, m)?)?;
    m.add_function(wrap_pyfunction!(depict_grid, m)?)?;
    m.add_function(wrap_pyfunction!(run_smirks, m)?)?;
    m.add_function(wrap_pyfunction!(find_mcs, m)?)?;
    m.add_function(wrap_pyfunction!(reaction_smarts_match, m)?)?;
    m.add_function(wrap_pyfunction!(shape_screen, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    // bulk submodule
    let bulk_mod = PyModule::new(m.py(), "bulk")?;
    bulk::register(&bulk_mod)?;
    m.add_submodule(&bulk_mod)?;

    // io functions (iter_sdf, iter_sdf_str)
    io::register(m)?;

    // SimilarityIndex class
    index::register(m)?;

    Ok(())
}
