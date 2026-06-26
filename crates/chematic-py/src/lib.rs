use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::sync::Arc;

mod bulk;
mod index;
mod io;

// ---------------------------------------------------------------------------
// Mol — the main Python-facing molecule wrapper
// ---------------------------------------------------------------------------

/// A rendered HTML report returned by `chematic.report()` and `chematic.compare()`.
///
/// In Jupyter, writing ``report`` in a cell renders the HTML automatically.
/// Use ``report.save("path.html")`` to write to disk, or ``str(report)`` to get the HTML string.
///
///     report = chematic.report(mols, names=["aspirin", "ibuprofen"])
///     report.save("report.html")   # write to disk
///     display(report)              # Jupyter: renders inline
#[pyclass]
struct Report {
    html: String,
}

#[pymethods]
impl Report {
    fn _repr_html_(&self) -> &str {
        &self.html
    }

    fn save(&self, path: &str) -> PyResult<()> {
        std::fs::write(path, &self.html)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    fn __str__(&self) -> &str {
        &self.html
    }

    fn __repr__(&self) -> String {
        format!("Report({} bytes)", self.html.len())
    }
}

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

    /// Generate a random (non-canonical) SMILES string for this molecule.
    ///
    /// ``seed`` controls the atom traversal order deterministically.
    /// Different seeds produce syntactically different but chemically identical SMILES.
    /// Useful for ML data augmentation (SMILES enumeration).
    ///
    ///     smi = mol.random_smiles(42)   # e.g. "OCC" for ethanol
    fn random_smiles(&self, seed: u64) -> String {
        chematic_smiles::random_smiles(&self.inner, seed)
    }

    /// Generate ``n`` unique random SMILES strings from this molecule.
    ///
    /// Uses sequential seeds starting from ``seed``. Returns up to ``n`` distinct
    /// SMILES (may be fewer if the molecule has limited symmetry variation).
    ///
    ///     variants = mol.random_smiles_n(n=5, seed=42)
    #[pyo3(signature = (n, seed=0))]
    fn random_smiles_n(&self, n: usize, seed: u64) -> Vec<String> {
        chematic_smiles::random_smiles_vect(&self.inner, n, seed)
    }

    /// Return ``True`` if atoms ``a`` and ``b`` are symmetry-equivalent.
    ///
    /// Two atoms are equivalent when they have the same Morgan canonical rank
    /// (i.e., they are interchangeable by a molecular symmetry operation).
    ///
    ///     assert mol.are_atoms_equivalent(0, 1)  # in benzene, all C are equivalent
    fn are_atoms_equivalent(&self, a: usize, b: usize) -> bool {
        chematic_smiles::are_atoms_equivalent(
            &self.inner,
            chematic_core::AtomIdx(a as u32),
            chematic_core::AtomIdx(b as u32),
        )
    }

    /// Total MMFF94 force field energy for the given 3D coordinates.
    ///
    /// Returns energy in kcal/mol. Returns ``0.0`` if MMFF94 typing fails
    /// (e.g., for molecules with unsupported elements).
    ///
    /// Complements :meth:`mmff94_energy_breakdown` (which returns per-component energies).
    ///
    ///     coords = mol.generate_3d()
    ///     e = mol.mmff94_total_energy(coords)  # kcal/mol
    fn mmff94_total_energy(&self, coords: Vec<[f64; 3]>) -> f64 {
        chematic_ff::mmff94_total_energy(&self.inner, &coords).unwrap_or(0.0)
    }

    /// Per-atom MMFF94 force field type names.
    ///
    /// Returns one string per heavy atom describing the MMFF94 atom type
    /// (e.g. ``"C_sp3"``, ``"C=O"``, ``"N_amide"``).
    /// Returns an empty list if MMFF94 typing is not supported for this molecule.
    ///
    /// Equivalent to inspecting RDKit's ``AllChem.MMFFGetMoleculeProperties(mol).GetMMFFAtomType(i)``.
    ///
    ///     types = mol.mmff94_atom_types()
    fn mmff94_atom_types(&self) -> Vec<String> {
        chematic_ff::assign_mmff94_types(&self.inner)
            .map(|types| types.iter().map(|t| format!("{t}")).collect())
            .unwrap_or_default()
    }

    /// IUPAC systematic name with CIP stereochemistry prefix.
    ///
    /// Like :attr:`iupac_name` but prepends ``(R)-``/``(S)-`` (or multi-centre
    /// descriptors like ``(1R,2S)-``) when stereocentres are present.
    /// Returns an empty string for structures outside the IUPAC naming scope.
    ///
    ///     name = mol.iupac_name_stereo()   # "(R)-butan-2-ol"
    fn iupac_name_stereo(&self) -> String {
        chematic_chem::iupac_name_stereo(&self.inner).unwrap_or_default()
    }

    /// Number of RECAP-breakable bonds (C–N, C–O, C–S single bonds).
    ///
    /// RECAP (Retrosynthetic Combinatorial Analysis Procedure) breaks bonds
    /// representing common synthetic transformations. Useful for SAR analysis.
    #[getter]
    fn recap_breakable_bond_count(&self) -> usize {
        chematic_chem::recap_breakable_bond_count(&self.inner)
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

    /// Serialize to MDL MOL V2000 format (without 3D coordinates).
    ///
    /// Equivalent to RDKit's ``Chem.MolToMolBlock(mol)``.
    /// Use :meth:`to_mol2` for Tripos format, :meth:`to_pdb` for PDB with 3D coords.
    ///
    ///     block = mol.to_mol_block()
    ///     with open("molecule.mol", "w") as f:
    ///         f.write(block)
    fn to_mol_block(&self) -> String {
        chematic_mol::write_mol(&self.inner, &chematic_mol::MolMetadata::default())
    }

    /// Serialize this molecule to MDL MOL V2000 format with 2D layout coordinates.
    ///
    /// Preserves the 2D layout obtained from :func:`from_mol_block_with_coords`.
    /// Each element of ``coords`` is an ``[x, y]`` pair in Å.
    ///
    ///     mol, name, coords_2d = chematic.from_mol_block_with_coords(block)
    ///     new_block = mol.to_mol_block_2d(coords_2d, name=name)
    #[pyo3(signature = (coords, name = None))]
    fn to_mol_block_2d(&self, coords: Vec<[f64; 2]>, name: Option<&str>) -> String {
        let metadata = chematic_mol::MolMetadata {
            name: name.unwrap_or("").to_string(),
            comment: String::new(),
        };
        let coords_2d: Vec<(f64, f64)> = coords.iter().map(|c| (c[0], c[1])).collect();
        chematic_mol::write_mol_with_coords(&self.inner, &metadata, &coords_2d)
    }

    /// Serialize this molecule to MDL MOL V3000 format with 2D layout coordinates.
    ///
    /// V3000 supports >999 atoms and extended atom/bond features.
    /// Accepts the same ``[[x, y], ...]`` coordinate format as :meth:`to_mol_block_2d`.
    /// Pass an empty list to write zero coordinates.
    ///
    /// Equivalent to RDKit ``Chem.MolToV3KMolBlock(mol)``.
    ///
    ///     block = mol.to_mol_v3000(coords_2d, name="my_mol")
    #[pyo3(signature = (coords, name = None))]
    fn to_mol_v3000(&self, coords: Vec<[f64; 2]>, name: Option<&str>) -> String {
        let metadata = chematic_mol::MolMetadata {
            name: name.unwrap_or("").to_string(),
            comment: String::new(),
        };
        let coords_2d: Vec<(f64, f64)> = coords.iter().map(|c| (c[0], c[1])).collect();
        chematic_mol::write_mol_v3000(&self.inner, &metadata, &coords_2d)
    }

    /// Serialize this molecule to Chemical Markup Language (CML) XML.
    ///
    /// ``coords``: optional list of ``[x, y]`` pairs (Å) — one per heavy atom.
    /// If omitted or ``None``, no coordinate attributes are written.
    ///
    /// Equivalent to RDKit ``Chem.MolToCMLBlock(mol)``.
    ///
    ///     cml = mol.to_cml()
    ///     cml_with_layout = mol.to_cml(coords_2d)
    #[pyo3(signature = (coords = None))]
    fn to_cml(&self, coords: Option<Vec<[f64; 2]>>) -> String {
        let coords_2d: Option<Vec<(f64, f64)>> =
            coords.map(|c| c.iter().map(|xy| (xy[0], xy[1])).collect());
        chematic_mol::write_cml(&self.inner, coords_2d.as_deref())
    }

    /// Write this molecule to AutoDock PDBQT format.
    ///
    /// Args:
    ///     coords: list of ``(x, y, z)`` tuples (Å) — one per heavy atom.
    ///             Pass ``[]`` to write zero coordinates.
    ///     charges: list of partial charges (float) — one per heavy atom.
    ///              Pass ``[]`` to write zeros. Use
    ///              ``chematic_ff.gasteiger_charges()`` or MMFF94 BCI for
    ///              best docking accuracy.
    ///     name: ligand name embedded in the REMARK header (e.g. ``"LIG"``).
    ///
    /// Returns:
    ///     str: PDBQT-format string (rigid body, no torsion tree).
    ///
    /// Example::
    ///
    ///     mol = chematic.from_smiles("CCO")
    ///     pdbqt = mol.to_pdbqt([], [], "ETH")
    ///     with open("ethanol.pdbqt", "w") as f:
    ///         f.write(pdbqt)
    fn to_pdbqt(&self, coords: Vec<(f64, f64, f64)>, charges: Vec<f64>, name: &str) -> String {
        chematic_mol::write_pdbqt(&self.inner, &coords, &charges, name)
    }

    /// Minimise geometry using the Universal Force Field (UFF, Rappé 1992).
    ///
    /// UFF handles all elements including metals, making it suitable for
    /// metal-ligand complexes where MMFF94 is not parameterised.
    ///
    /// Args:
    ///     coords: list of ``[x, y, z]`` lists (Å) — initial 3D coordinates.
    ///     max_iter: maximum steepest-descent iterations (default 500).
    ///
    /// Returns:
    ///     dict with keys:
    ///     ``coords`` (list of [x,y,z]), ``energy`` (float, kcal/mol),
    ///     ``iterations`` (int), ``converged`` (bool).
    ///
    /// Example::
    ///
    ///     mol = chematic.from_smiles("CCO")
    ///     result = mol.minimize_uff([[0,0,0],[1.54,0,0],[2.5,1.2,0]])
    ///     print(result["energy"], result["converged"])
    fn minimize_uff<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
        max_iter: Option<usize>,
    ) -> PyResult<pyo3::Bound<'py, pyo3::types::PyDict>> {
        let types = chematic_ff::assign_uff_types(&self.inner);
        let result =
            chematic_ff::minimize_uff(&self.inner, &types, coords, max_iter.unwrap_or(500));
        let d = pyo3::types::PyDict::new(py);
        let py_coords: Vec<Vec<f64>> = result.coords.iter().map(|c| c.to_vec()).collect();
        d.set_item("coords", py_coords)?;
        d.set_item("energy", result.energy)?;
        d.set_item("iterations", result.iterations)?;
        d.set_item("converged", result.converged)?;
        Ok(d)
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

    /// DrugScore — composite drug-likeness score [0, 1].
    ///
    /// Product of four sigmoidal factors:
    /// - **cLogP**: optimal 0–5
    /// - **logS** (ESOL): optimal > −5 log mol/L
    /// - **MW**: optimal < 500 Da
    /// - **Toxicity**: 0.5× per PAINS alert, 0.75× per Brenk alert
    ///
    /// Analogous to OCL DrugScore.  Use alongside :attr:`qed` for a second
    /// opinion on drug-likeness.
    #[getter]
    fn drug_score(&self) -> f64 {
        chematic_chem::drug_score(&self.inner)
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

    /// Number of distinct connected ring systems.
    ///
    /// Two SSSR rings form the same system when they share at least one atom
    /// (fused, bridged, or spiro).  Differs from :attr:`ring_count` which
    /// counts SSSR rings individually.
    ///
    ///     naphthalene = chematic.from_smiles("c1ccc2ccccc2c1")
    ///     naphthalene.ring_system_count  # 1  (two fused rings = one system)
    ///     biphenyl = chematic.from_smiles("c1ccc(-c2ccccc2)cc1")
    ///     biphenyl.ring_system_count     # 2  (two independent ring systems)
    #[getter]
    fn ring_system_count(&self) -> usize {
        chematic_chem::ring_system_count(&self.inner)
    }

    /// Lipinski (1997) HBA count — total number of N and O heavy atoms.
    ///
    /// This is the original Rule-of-Five HBA definition: count all N and O atoms
    /// regardless of hybridisation.  For the more accurate Ertl (2000) definition
    /// use :attr:`hba`.
    ///
    ///     caffeine = chematic.from_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
    ///     caffeine.hba_count_lipinski  # 5  (2 O + 3 N)
    #[getter]
    fn hba_count_lipinski(&self) -> usize {
        chematic_chem::hba_count_lipinski(&self.inner)
    }

    /// Fraction of heavy atoms that are rotatable bonds (0.0–1.0).
    ///
    /// `fraction_rotatable_bonds = rotatable_bond_count / heavy_atom_count`.
    /// Returns ``0.0`` for acyclic molecules with no rotatable bonds, or when
    /// the molecule has no heavy atoms.
    ///
    ///     benzene = chematic.from_smiles("c1ccccc1")
    ///     benzene.fraction_rotatable_bonds  # 0.0
    #[getter]
    fn fraction_rotatable_bonds(&self) -> f64 {
        chematic_chem::fraction_rotatable_bonds(&self.inner)
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

    /// Med-Chem Friendly (MCF) composite filter.
    ///
    /// Returns ``True`` when all of the following hold:
    /// no PAINS alerts, no Brenk alerts, Lipinski Ro5, and Veber oral bioavailability.
    /// Mirrors the "MCF" concept in the `medchemfilters` Python library.
    ///
    ///     caffeine = chematic.from_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
    ///     caffeine.mcf_passes  # True
    #[getter]
    fn mcf_passes(&self) -> bool {
        chematic_chem::mcf_passes(&self.inner)
    }

    /// True if Rule of Three criteria pass (fragment-based drug discovery, Congreve 2003).
    ///
    /// Passes when MW ≤ 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3, RotBonds ≤ 3.
    #[getter]
    fn ro3_passes(&self) -> bool {
        chematic_chem::ro3_passes(&self.inner)
    }

    /// True if lead-like criteria pass (Oprea 2001).
    ///
    /// Passes when MW ≤ 450, LogP −3.5–4.5, RotBonds ≤ 10, RingCount 1–4.
    #[getter]
    fn lead_like_passes(&self) -> bool {
        chematic_chem::lead_like_passes(&self.inner)
    }

    /// True if compound is NOT in the Pfizer 3/75 high-metabolic-liability zone.
    ///
    /// The danger zone is ``LogP > 3 AND TPSA < 75``; compounds there have higher CYP3A4
    /// metabolic clearance risk (Leeson & Springthorpe 2007).
    #[getter]
    fn pfizer_3_75_passes(&self) -> bool {
        chematic_chem::pfizer_3_75_passes(&self.inner)
    }

    /// CNS Multi-Parameter Optimisation (MPO) score (Wager 2010), range 0–6.
    ///
    /// Combines desirability functions for cLogP, cLogD, MW, TPSA, HBD, and pKa.
    /// Higher scores indicate better CNS drug-like properties.
    #[getter]
    fn cns_mpo_score(&self) -> f64 {
        chematic_chem::cns_mpo_score(&self.inner)
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
        d.set_item(
            "clearance",
            match chematic_chem::clearance_class(m) {
                chematic_chem::ClearanceClass::Low => "Low",
                chematic_chem::ClearanceClass::Medium => "Medium",
                chematic_chem::ClearanceClass::High => "High",
            },
        )?;
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
        // Pre-compute all ring-derived values with a single find_sssr call.
        let rb = chematic_chem::ring_bundle(m);
        let mw = chematic_chem::molecular_weight(m);
        // Share the 117-pattern Crippen SMARTS pass between logp and molar_refractivity.
        let (logp, mr) = chematic_chem::logp_and_mr(m);
        let tpsa = chematic_chem::tpsa(m);
        d.set_item("mw", mw)?;
        d.set_item("exact_mass", chematic_chem::exact_mass(m))?;
        d.set_item("tpsa", tpsa)?;
        d.set_item("logp", logp)?;
        d.set_item("molar_refractivity", mr)?;
        let hbd = chematic_chem::hbd_count(m);
        d.set_item("hbd", hbd)?;
        d.set_item("hba", rb.hba_count)?;
        d.set_item("rotatable_bonds", rb.rotatable_bond_count)?;
        d.set_item("heavy_atoms", chematic_chem::heavy_atom_count(m))?;
        d.set_item("ring_count", rb.ring_count)?;
        d.set_item("ring_system_count", rb.ring_system_count)?;
        d.set_item("aromatic_ring_count", rb.aromatic_ring_count)?;
        d.set_item("hba_lipinski", chematic_chem::hba_count_lipinski(m))?;
        d.set_item("fraction_rotatable_bonds", rb.fraction_rotatable_bonds)?;
        d.set_item("num_heteroatoms", chematic_chem::num_heteroatoms(m))?;
        d.set_item("num_stereocenters", chematic_chem::num_stereocenters(m))?;
        d.set_item("num_spiro_atoms", rb.num_spiro_atoms)?;
        d.set_item("num_bridgehead_atoms", rb.num_bridgehead_atoms)?;
        d.set_item("fsp3", chematic_chem::fsp3(m))?;
        d.set_item("qed", chematic_chem::qed_with_bundle(m, &rb))?;
        d.set_item("sa_score", chematic_chem::sa_score_with_bundle(m, &rb))?;
        d.set_item("formal_charge", chematic_chem::formal_charge_sum(m))?;
        d.set_item("labute_asa", chematic_chem::labute_asa(m))?;
        d.set_item("bertz_ct", chematic_chem::bertz_ct(m))?;
        d.set_item("wiener_index", chematic_chem::wiener_index(m))?;
        d.set_item("schultz_mti", chematic_chem::schultz_mti(m))?;
        d.set_item("gutman_mti", chematic_chem::gutman_mti(m))?;
        d.set_item("vabc", chematic_chem::vabc(m))?;
        d.set_item("gravitational_index", chematic_chem::gravitational_index(m))?;
        d.set_item("zagreb_m2", chematic_chem::zagreb_index_m2(m))?;
        d.set_item("kappa1", chematic_chem::kappa1(m))?;
        d.set_item("kappa2", chematic_chem::kappa2(m))?;
        d.set_item("kappa3", chematic_chem::kappa3(m))?;
        // Compute all 10 chi indices in a single heavy_indices pass.
        let (c0, c1, c2, c3, c4, c0v, c1v, c2v, c3v, c4v) = chematic_chem::chi_all(m);
        d.set_item("chi0", c0)?;
        d.set_item("chi1", c1)?;
        d.set_item("chi2", c2)?;
        d.set_item("chi3", c3)?;
        d.set_item("chi4", c4)?;
        d.set_item("chi0v", c0v)?;
        d.set_item("chi1v", c1v)?;
        d.set_item("chi2v", c2v)?;
        d.set_item("chi3v", c3v)?;
        d.set_item("chi4v", c4v)?;
        d.set_item("num_aromatic_heterocycles", rb.num_aromatic_heterocycles)?;
        d.set_item("num_aliphatic_heterocycles", rb.num_aliphatic_heterocycles)?;
        d.set_item("num_saturated_rings", rb.num_saturated_rings)?;
        d.set_item("num_aliphatic_rings", rb.num_aliphatic_rings)?;
        d.set_item(
            "num_unspecified_stereocenters",
            chematic_chem::num_unspecified_stereocenters(m),
        )?;
        let (sum_e, max_e, min_e) = chematic_chem::estate_all(m);
        d.set_item("sum_estate", sum_e)?;
        d.set_item("max_estate", max_e)?;
        d.set_item("min_estate", min_e)?;
        // Inline only filters that use rotatable_bond_count or hba_count (→ find_sssr).
        d.set_item(
            "lipinski_passes",
            mw <= 500.0 && hbd <= 5 && rb.hba_count <= 10 && logp <= 5.0,
        )?;
        d.set_item(
            "veber_passes",
            rb.rotatable_bond_count <= 10 && tpsa <= 140.0,
        )?;
        d.set_item("egan_passes", tpsa <= 131.6 && logp <= 5.88)?;
        d.set_item("ghose_passes", chematic_chem::ghose_passes(m))?;
        d.set_item("reos_passes", chematic_chem::reos_passes(m))?;
        d.set_item("pains_passes", chematic_chem::pains_passes(m))?;
        d.set_item("ro3_passes", chematic_chem::ro3_passes(m))?;
        d.set_item("lead_like_passes", chematic_chem::lead_like_passes(m))?;
        d.set_item("pfizer_3_75_passes", chematic_chem::pfizer_3_75_passes(m))?;
        // Reuse pre-computed logp/tpsa/mw/hbd; only pka_base is a new computation.
        let pka_b = chematic_chem::pka_base(m).unwrap_or(0.0);
        d.set_item(
            "cns_mpo_score",
            chematic_chem::cns_mpo_from_parts(m, logp, tpsa, mw, hbd, pka_b),
        )?;
        d.set_item("mcf_passes", chematic_chem::mcf_passes(m))?;
        d.set_item("bbb_score", chematic_chem::bbb_score_from_parts(tpsa, logp))?;
        d.set_item("bbb_passes", tpsa < 90.0 && mw < 400.0 && hbd <= 3)?;
        d.set_item("caco2", chematic_chem::caco2_precomputed(tpsa, logp))?;
        d.set_item(
            "herg_risk",
            chematic_chem::herg_risk_precomputed(m, logp, mw),
        )?;
        d.set_item(
            "cyp3a4_risk",
            chematic_chem::cyp3a4_precomputed(mw, logp, rb.num_aromatic_heterocycles, rb.hba_count),
        )?;
        let egg = chematic_chem::boiled_egg_from(logp, tpsa);
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
        // VSA descriptors (47 values: SlogP_VSA1-12, SMR_VSA1-10, PEOE_VSA1-14, EState_VSA1-11)
        for (i, v) in chematic_chem::slogp_vsa(m).into_iter().enumerate() {
            d.set_item(format!("SlogP_VSA{}", i + 1), v)?;
        }
        for (i, v) in chematic_chem::smr_vsa(m).into_iter().enumerate() {
            d.set_item(format!("SMR_VSA{}", i + 1), v)?;
        }
        for (i, v) in chematic_chem::peoe_vsa(m).into_iter().enumerate() {
            d.set_item(format!("PEOE_VSA{}", i + 1), v)?;
        }
        for (i, v) in chematic_chem::estate_vsa(m).into_iter().enumerate() {
            d.set_item(format!("EState_VSA{}", i + 1), v)?;
        }
        d.set_item("num_valence_electrons", chematic_chem::num_valence_electrons(m))?;
        d.set_item("hall_kier_alpha", chematic_chem::hall_kier_alpha(m))?;
        // BCUT2D descriptors (8 eigenvalue-based values)
        let bc = chematic_chem::bcut2d(m);
        d.set_item("bcut2d_chghi", bc.chghi)?;
        d.set_item("bcut2d_chglo", bc.chglo)?;
        d.set_item("bcut2d_logphi", bc.logphi)?;
        d.set_item("bcut2d_logplo", bc.logplo)?;
        d.set_item("bcut2d_mrhi", bc.mrhi)?;
        d.set_item("bcut2d_mrlo", bc.mrlo)?;
        d.set_item("bcut2d_mwhi", bc.mwhi)?;
        d.set_item("bcut2d_mwlo", bc.mwlo)?;
        // Carbon type breakdown (hybridisation × degree)
        let ct = chematic_chem::carbon_types(m);
        d.set_item("c1sp1", ct.c1sp1)?;
        d.set_item("c2sp1", ct.c2sp1)?;
        d.set_item("c1sp2", ct.c1sp2)?;
        d.set_item("c2sp2", ct.c2sp2)?;
        d.set_item("c3sp2", ct.c3sp2)?;
        d.set_item("c1sp3", ct.c1sp3)?;
        d.set_item("c2sp3", ct.c2sp3)?;
        d.set_item("c3sp3", ct.c3sp3)?;
        // Connectivity index + bond type counts
        d.set_item("balaban_j", chematic_chem::balaban_j(m))?;
        d.set_item("num_amide_bonds", chematic_chem::num_amide_bonds(m))?;
        d.set_item("num_ester_bonds", chematic_chem::num_ester_bonds(m))?;
        // Element-wise heavy atom counts
        d.set_item("num_carbons", chematic_chem::num_carbons(m))?;
        d.set_item("num_nitrogens", chematic_chem::num_nitrogens(m))?;
        d.set_item("num_oxygens", chematic_chem::num_oxygens(m))?;
        d.set_item("num_sulfurs", chematic_chem::num_sulfurs(m))?;
        d.set_item("num_phosphorus", chematic_chem::num_phosphorus(m))?;
        d.set_item("num_fluorines", chematic_chem::num_fluorines(m))?;
        d.set_item("num_chlorines", chematic_chem::num_chlorines(m))?;
        d.set_item("num_bromines", chematic_chem::num_bromines(m))?;
        d.set_item("num_iodines", chematic_chem::num_iodines(m))?;
        // ADME / solubility / alternative LogP (not in RDKit standard)
        d.set_item("esol", chematic_chem::esol_solubility(m))?;
        d.set_item("logd_7_4", chematic_chem::logd_simple(m, 7.4))?;
        d.set_item("xlogp3", chematic_chem::xlogp3(m))?;
        d.set_item("drug_score", chematic_chem::drug_score(m))?;
        // MQN (Molecular Quantum Numbers): 42 integer descriptors (Ertl 2010)
        for (i, &v) in chematic_chem::mqn(m).iter().enumerate() {
            d.set_item(format!("MQN{}", i + 1), u32::from(v))?;
        }
        Ok(d)
    }

    /// Human-readable summary of key properties — suitable for LLM prompts or MCP responses.
    ///
    /// ```python
    /// mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    /// print(mol.describe())
    /// ```
    fn describe(&self) -> String {
        let m = &self.inner;
        let rb = chematic_chem::ring_bundle(m);
        let mw = chematic_chem::molecular_weight(m);
        let (logp, _) = chematic_chem::logp_and_mr(m);
        let tpsa = chematic_chem::tpsa(m);
        let hbd = chematic_chem::hbd_count(m);
        let hba = rb.hba_count;
        let hba_lipinski = chematic_chem::hba_count_lipinski(m);
        let rot = rb.rotatable_bond_count;
        let arom = rb.aromatic_ring_count;
        let qed = chematic_chem::qed_with_bundle(m, &rb);
        let pains_ok = chematic_chem::pains_passes(m);
        let brenk_ok = chematic_chem::brenk_passes(m);

        // Lipinski rule-of-5 (original paper uses N+O count for HBA)
        let lipinski_violations: u8 = [(mw > 500.0) as u8, (hbd > 5) as u8, (hba_lipinski > 10) as u8, (logp > 5.0) as u8].iter().sum();

        // Lipopolicity characterisation
        let lipophilicity = if logp < 0.0 {
            "hydrophilic"
        } else if logp < 2.0 {
            "mildly lipophilic"
        } else if logp < 4.0 {
            "moderately lipophilic"
        } else {
            "highly lipophilic"
        };

        // Oral absorption proxy (Veber)
        let oral_absorption = if tpsa <= 140.0 && rot <= 10 {
            "Likely orally bioavailable (passes Veber criteria)"
        } else {
            "Oral bioavailability may be limited (fails Veber criteria)"
        };

        // Drug-likeness description
        let drug_like = if lipinski_violations == 0 {
            "no Lipinski rule-of-5 violations"
        } else if lipinski_violations == 1 {
            "1 Lipinski violation (borderline drug-like)"
        } else {
            "multiple Lipinski violations (likely not orally drug-like)"
        };

        let mut lines = Vec::new();
        lines.push(format!(
            "Molecular weight {mw:.1} Da, formula {}.",
            chematic_chem::calc_mol_formula(m)
        ));
        lines.push(format!(
            "LogP {logp:.2} ({lipophilicity}), TPSA {tpsa:.1} Å²."
        ));
        lines.push(format!(
            "HBD {hbd}, HBA {hba}, {rot} rotatable bond(s), {arom} aromatic ring(s)."
        ));
        lines.push(format!("Drug-likeness: {drug_like}. {oral_absorption}."));
        lines.push(format!("QED {qed:.2} (0 = non-drug-like, 1 = ideal)."));

        let mut alerts = Vec::new();
        if !pains_ok { alerts.push("PAINS alert"); }
        if !brenk_ok { alerts.push("Brenk alert"); }
        if alerts.is_empty() {
            lines.push("No structural alerts (PAINS / Brenk clean).".to_string());
        } else {
            lines.push(format!("Structural alerts: {}.", alerts.join(", ")));
        }

        lines.join("\n")
    }

    /// Return a structured Markdown analysis suitable for LLM prompts, Jupyter, or standalone reports.
    ///
    /// Sections: Structure · Physical Properties · Drug-likeness · ADMET Predictions.
    ///
    /// ```python
    /// md = mol.review()
    /// print(md)                            # human-readable Markdown
    /// prompt = f"Evaluate:\n\n{md}"        # LLM input
    /// open("review.md","w").write(md)      # save to file
    ///
    /// from IPython.display import Markdown, display
    /// display(Markdown(mol.review()))      # Jupyter rendered
    /// ```
    fn review(&self) -> String {
        use chematic_chem::{
            admet_profile, brenk_passes, calc_mol_formula, hba_count_lipinski, hbd_count,
            logp_and_mr, molecular_weight, pains_passes, qed_with_bundle, ring_bundle, tpsa,
        };

        let m = &self.inner;
        let mw = molecular_weight(m);
        let (logp, _) = logp_and_mr(m);
        let tpsa_val = tpsa(m);
        let hbd = hbd_count(m);
        let hba_lip = hba_count_lipinski(m);
        let rb = ring_bundle(m);
        let qed = qed_with_bundle(m, &rb);
        let pains_ok = pains_passes(m);
        let brenk_ok = brenk_passes(m);
        let admet = admet_profile(m);
        let smiles = chematic_smiles::canonical_smiles(m);
        let formula = calc_mol_formula(m);
        let n_heavy = m.atom_count();
        let ver = env!("CARGO_PKG_VERSION");

        // Lipinski violations
        let violations: u8 = [mw > 500.0, hbd > 5, hba_lip > 10, logp > 5.0]
            .iter()
            .filter(|&&b| b)
            .count() as u8;
        let lip_result = match violations {
            0 => "✓ Pass (0 violations)",
            1 => "⚠ Borderline (1 violation)",
            _ => "✗ Fail (multiple violations)",
        };
        let veber_ok = tpsa_val <= 140.0 && rb.rotatable_bond_count <= 10;
        let veber_result = if veber_ok { "✓ Likely (TPSA ≤ 140, rot ≤ 10)" } else { "✗ Unlikely" };

        let lipophilicity = if logp < 0.0 { "hydrophilic" }
            else if logp < 2.0 { "mildly lipophilic" }
            else if logp < 4.0 { "moderately lipophilic" }
            else { "highly lipophilic" };

        // ADMET labels
        let bbb = if admet.bbb_passes { "✓ CNS penetrant" } else { "✗ CNS non-penetrant" };
        let caco2 = if admet.caco2 > -5.5 { "High (well absorbed)" }
            else if admet.caco2 > -7.0 { "Moderate" }
            else { "Low (poor absorption)" };
        let herg = if admet.herg_risk < 0.3 { "Low" }
            else if admet.herg_risk < 0.6 { "Moderate" }
            else { "High ⚠" };
        let cyp = if admet.cyp3a4_risk < 0.3 { "Low" }
            else if admet.cyp3a4_risk < 0.6 { "Moderate" }
            else { "High ⚠" };

        format!(
            "# Molecular Review\n\n\
             ## Structure\n\
             - Formula: {formula}\n\
             - SMILES: `{smiles}`\n\
             - Heavy atoms: {n_heavy}\n\n\
             ## Physical Properties\n\
             | Property | Value |\n\
             |---|---|\n\
             | MW | {mw:.1} Da |\n\
             | LogP | {logp:.2} ({lipophilicity}) |\n\
             | TPSA | {tpsa_val:.1} Å² |\n\
             | HBD / HBA | {hbd} / {hba} |\n\
             | Rotatable bonds | {rot} |\n\
             | Aromatic rings | {arom} |\n\
             | QED | {qed:.2} |\n\n\
             ## Drug-likeness\n\
             | Filter | Result |\n\
             |---|---|\n\
             | Lipinski RO5 | {lip_result} |\n\
             | Oral bioavailability (Veber) | {veber_result} |\n\
             | PAINS | {pains_s} |\n\
             | Brenk | {brenk_s} |\n\n\
             ## ADMET Predictions\n\
             | Property | Prediction |\n\
             |---|---|\n\
             | BBB penetration | {bbb} |\n\
             | Caco-2 permeability | {caco2} |\n\
             | hERG cardiac risk | {herg} |\n\
             | CYP3A4 inhibition | {cyp} |\n\n\
             ---\n\
             *Generated by chematic v{ver}*",
            hba = rb.hba_count,
            rot = rb.rotatable_bond_count,
            arom = rb.aromatic_ring_count,
            pains_s = if pains_ok { "✓ Clean" } else { "✗ Alert detected" },
            brenk_s = if brenk_ok { "✓ Clean" } else { "⚠ Alert detected" },
        )
    }

    /// Compare this molecule to another and return element-level and descriptor-level differences.
    ///
    /// The diff is directional: `mol1.diff(mol2)` reports changes going *from* mol1 *to* mol2.
    /// Positive delta values mean mol2 has more; negative values mean mol1 has more.
    ///
    /// ```python
    /// aspirin = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    /// ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")
    /// d = aspirin.diff(ibuprofen)
    /// print(d["summary"])
    /// # "+C7, -O2. ΔLogP +2.75, ΔTPSA -26.3 Å², ΔMW +66.1 Da."
    /// ```
    fn diff<'py>(&self, py: Python<'py>, other: &Mol) -> PyResult<Bound<'py, PyDict>> {
        use std::collections::BTreeMap;
        let mol1 = self.inner.as_ref();
        let mol2 = other.inner.as_ref();

        // MCS size as common-core reference
        let config = chematic_smarts::McsConfig::default();
        let qmol = chematic_smarts::find_mcs_with_config(&[mol1, mol2], &config);
        let common_atoms = qmol.atoms.len();

        // Count heavy-atom elements in each molecule
        let mut counts1: BTreeMap<&'static str, i32> = BTreeMap::new();
        for i in 0..mol1.atom_count() {
            let sym = mol1.atom(chematic_core::AtomIdx(i as u32)).element.symbol();
            *counts1.entry(sym).or_insert(0) += 1;
        }
        let mut counts2: BTreeMap<&'static str, i32> = BTreeMap::new();
        for i in 0..mol2.atom_count() {
            let sym = mol2.atom(chematic_core::AtomIdx(i as u32)).element.symbol();
            *counts2.entry(sym).or_insert(0) += 1;
        }

        // Element-level delta (mol2 - mol1); zero-delta elements are omitted
        let all_elems: std::collections::BTreeSet<_> =
            counts1.keys().chain(counts2.keys()).copied().collect();
        let mut delta_elements: BTreeMap<&'static str, i32> = BTreeMap::new();
        for elem in &all_elems {
            let d = counts2.get(elem).copied().unwrap_or(0)
                - counts1.get(elem).copied().unwrap_or(0);
            if d != 0 {
                delta_elements.insert(elem, d);
            }
        }

        // Descriptor deltas
        let rb1 = chematic_chem::ring_bundle(mol1);
        let rb2 = chematic_chem::ring_bundle(mol2);
        let (logp1, _) = chematic_chem::logp_and_mr(mol1);
        let (logp2, _) = chematic_chem::logp_and_mr(mol2);
        let tpsa1 = chematic_chem::tpsa(mol1);
        let tpsa2 = chematic_chem::tpsa(mol2);
        let mw1 = chematic_chem::molecular_weight(mol1);
        let mw2 = chematic_chem::molecular_weight(mol2);
        let hbd1 = chematic_chem::hbd_count(mol1) as i32;
        let hbd2 = chematic_chem::hbd_count(mol2) as i32;

        // Human-readable element summary
        let elem_parts: Vec<String> = delta_elements.iter().map(|(e, d)| {
            if *d > 0 { format!("+{d}{e}") } else { format!("{d}{e}") }
        }).collect();
        let elem_str = if elem_parts.is_empty() {
            "Same elemental composition".to_string()
        } else {
            elem_parts.join(", ")
        };
        let summary = format!(
            "{}. \u{0394}LogP {:+.2}, \u{0394}TPSA {:+.1} \u{00c5}\u{00b2}, \u{0394}MW {:+.1} Da.",
            elem_str,
            logp2 - logp1,
            tpsa2 - tpsa1,
            mw2 - mw1,
        );

        let d = PyDict::new(py);
        d.set_item("common_atoms", common_atoms)?;
        d.set_item("delta_mw", mw2 - mw1)?;
        d.set_item("delta_logp", logp2 - logp1)?;
        d.set_item("delta_tpsa", tpsa2 - tpsa1)?;
        d.set_item("delta_hbd", hbd2 - hbd1)?;
        d.set_item("delta_hba", rb2.hba_count as i32 - rb1.hba_count as i32)?;
        let elem_dict = PyDict::new(py);
        for (elem, delta) in &delta_elements {
            elem_dict.set_item(elem.to_string(), delta)?;
        }
        d.set_item("delta_elements", elem_dict)?;
        d.set_item("summary", summary)?;
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

    /// Export the 2D structure as an EPS string.
    ///
    /// Generates a self-contained PostScript EPS document.
    /// Suitable for publication-quality vector graphics in LaTeX, Illustrator, etc.
    ///
    ///     eps_str = mol.to_eps()
    ///     open("molecule.eps", "w").write(eps_str)
    fn to_eps(&self) -> String {
        chematic_depict::depict_eps(&self.inner)
    }

    /// Export as ChemicalJSON (.cjson) string.
    ///
    /// ``coords``: optional ``[[x,y,z], ...]`` list (Å), one per heavy atom.
    /// Pass an empty list (default) to export topology only (no coordinates).
    ///
    ///     coords = mol.generate_3d()
    ///     cjson_str = mol.to_cjson(coords)
    ///     open("mol.cjson", "w").write(cjson_str)
    ///
    ///     # Topology-only (no 3D):
    ///     cjson_str = mol.to_cjson()
    #[pyo3(signature = (coords = vec![]))]
    fn to_cjson(&self, coords: Vec<Vec<f64>>) -> String {
        let c3d: Vec<(f64, f64, f64)> = coords
            .iter()
            .map(|c| {
                (
                    c.first().copied().unwrap_or(0.0),
                    c.get(1).copied().unwrap_or(0.0),
                    c.get(2).copied().unwrap_or(0.0),
                )
            })
            .collect();
        chematic_mol::write_cjson(&self.inner, &c3d)
    }

    /// Export the 2D structure as a PDF document (bytes).
    ///
    ///     pdf_bytes = mol.to_pdf()
    ///     open("molecule.pdf", "wb").write(pdf_bytes)
    fn to_pdf(&self) -> Vec<u8> {
        chematic_depict::depict_pdf(&self.inner)
    }

    /// Jupyter Notebook / JupyterLab の自動描画フック。
    ///
    /// セルに ``mol`` と書くだけで 2D 構造が表示される。手動で
    /// ``IPython.display.SVG(mol.svg())`` と書く必要はない。
    fn _repr_svg_(&self) -> String {
        chematic_depict::depict_svg(&self.inner)
    }

    /// Return ``True`` if this molecule matches the given SMARTS pattern.
    ///
    /// Equivalent to ``chematic.smarts_match(smarts, mol)`` but as a method::
    ///
    ///     if mol.has_substructure("[OH]"):
    ///         print("has hydroxyl")
    ///
    /// Raises ``ValueError`` for invalid SMARTS.
    fn has_substructure(&self, smarts: &str) -> PyResult<bool> {
        let query = chematic_smarts::parse_smarts(smarts)
            .map_err(|e| PyValueError::new_err(format!("invalid SMARTS '{smarts}': {e}")))?;
        Ok(!chematic_smarts::find_matches(&query, &self.inner).is_empty())
    }

    /// Return atom-index lists for all SMARTS matches in this molecule.
    ///
    /// Equivalent to ``chematic.smarts_find(smarts, mol)`` but as a method::
    ///
    ///     for match_atoms in mol.find_matches("[CX3](=O)[OH]"):
    ///         print("carboxyl atoms:", match_atoms)
    ///
    /// Returns an empty list when there are no matches.
    /// Raises ``ValueError`` for invalid SMARTS.
    fn find_matches(&self, smarts: &str) -> PyResult<Vec<Vec<usize>>> {
        let query = chematic_smarts::parse_smarts(smarts)
            .map_err(|e| PyValueError::new_err(format!("invalid SMARTS '{smarts}': {e}")))?;
        Ok(chematic_smarts::find_matches(&query, &self.inner)
            .into_iter()
            .map(|m| {
                let mut v: Vec<usize> = m.values().map(|idx| idx.0 as usize).collect();
                v.sort_unstable();
                v
            })
            .collect())
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

    /// 2D SVG depiction with atoms coloured by LogP contribution.
    ///
    /// Positive (lipophilic) atoms → blue, negative (hydrophilic) atoms → red,
    /// zero-contribution atoms → white.  Uses :meth:`logp_per_atom` as weights.
    fn logp_map_svg(&self) -> String {
        let weights = chematic_chem::logp_crippen_per_atom(&self.inner);
        chematic_depict::similarity_map_svg(&self.inner, &weights)
    }

    /// 2D SVG depiction with atoms coloured by TPSA contribution.
    ///
    /// Atoms contributing to TPSA (N, O, S, P) are coloured blue; zero-contribution
    /// atoms (C, halogens, …) remain white.  Uses :meth:`tpsa_per_atom` as weights.
    fn tpsa_map_svg(&self) -> String {
        let weights = chematic_chem::tpsa_per_atom(&self.inner);
        chematic_depict::similarity_map_svg(&self.inner, &weights)
    }

    /// 2D SVG depiction with atoms coloured by custom weights.
    ///
    /// ``weights``: list of floats, one per heavy atom (length = :attr:`heavy_atoms`).
    /// Positive → blue, negative → red, zero → white.
    fn similarity_map_svg(&self, weights: Vec<f64>) -> String {
        chematic_depict::similarity_map_svg(&self.inner, &weights)
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

    /// Predict aqueous solubility (logS, log mol/L) using an ECFP4-based MLP.
    ///
    /// When trained weights are installed (`MLP_SOLUBILITY_TRAINED = true`),
    /// runs a neural-network forward pass on the molecule's ECFP4 fingerprint.
    /// Until then, transparently falls back to the Delaney ESOL linear regression.
    ///
    /// To install trained weights, run `scripts/train_solubility_mlp.py` and
    /// follow the instructions printed at the end.
    ///
    ///     mol = chematic.from_smiles("c1ccccc1")
    ///     logs = mol.ml_solubility   # same as .esol until weights are trained
    #[getter]
    fn ml_solubility(&self) -> f64 {
        chematic_chem::mlp_solubility(&self.inner)
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

    /// Return the canonical SMILES under the given canonicalization mode.
    ///
    /// ``mode`` is one of:
    ///
    /// - ``"normal"`` — standard canonical SMILES (stereo, charges, isotopes preserved)
    /// - ``"nostereo"`` — stereo information stripped
    /// - ``"backbone"`` — charges, isotopes, and stereo stripped (element + topology only)
    /// - ``"tautomer"`` — canonical tautomer, then canonical SMILES
    /// - ``"nostereo_tautomer"`` — tautomer normalization + stereo removal
    ///
    /// Analogous to OCL's IDcode five-mode canonicalization, but outputs SMILES.
    ///
    ///     mol = chematic.from_smiles("[C@@H](N)(C)C(=O)O")  # L-alanine
    ///     mol.canonical_smiles_mode("nostereo")  # → "CC(N)C(=O)O"
    ///     mol.canonical_smiles_mode("backbone")  # → "CC(N)C(=O)O"
    fn canonical_smiles_mode(&self, mode: &str) -> PyResult<String> {
        let m = match mode {
            "normal" => chematic_chem::CanonicalMode::Normal,
            "nostereo" => chematic_chem::CanonicalMode::NoStereo,
            "backbone" => chematic_chem::CanonicalMode::Backbone,
            "tautomer" => chematic_chem::CanonicalMode::Tautomer,
            "nostereo_tautomer" => chematic_chem::CanonicalMode::NoStereoTautomer,
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown canonical mode '{other}'; expected one of: \
                     normal, nostereo, backbone, tautomer, nostereo_tautomer"
                )));
            }
        };
        Ok(chematic_chem::canonical_smiles_mode(&self.inner, m))
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
        Mol {
            inner: Arc::new(chematic_chem::add_hydrogens(&self.inner)),
        }
    }

    /// Return a copy with all explicit hydrogen atoms removed.
    fn remove_hydrogens(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::remove_hydrogens(&self.inner)),
        }
    }

    /// Return a copy with all stereochemistry assignments removed.
    fn remove_stereo(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::remove_stereo(&self.inner)),
        }
    }

    /// Return a copy with all isotope labels removed.
    fn remove_isotopes(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::remove_isotopes(&self.inner)),
        }
    }

    /// Return the largest covalently connected fragment.
    fn largest_fragment(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::largest_fragment(&self.inner)),
        }
    }

    /// Split this molecule into its connected components (fragments).
    ///
    /// Returns a list of :class:`Mol` objects, one per connected component.
    /// A fully connected molecule returns a single-element list.
    /// Useful after :func:`run_reactants` to obtain clean individual products.
    ///
    /// Equivalent to RDKit's ``Chem.GetMolFrags(mol, asMols=True)``.
    ///
    ///     mol = chematic.from_smiles("CC.[NH3]")  # disconnected salt
    ///     parts = mol.connected_components()
    ///     # [Mol("CC"), Mol("N")]
    fn connected_components(&self) -> Vec<Mol> {
        self.inner
            .fragments()
            .into_iter()
            .map(|m| Mol { inner: Arc::new(m) })
            .collect()
    }

    /// Return ``True`` if this molecule and ``other`` represent the same chemical structure.
    ///
    /// Uses canonical SMILES comparison — reliable after v0.4.11 (Morgan bond-order fix, #14).
    /// Equivalent to :func:`chematic.are_identical`.
    ///
    ///     m1 = chematic.from_smiles("CC(=O)O")
    ///     m2 = chematic.from_smiles("OC(C)=O")
    ///     assert m1.is_same_as(m2)   # True — same acetic acid
    fn is_same_as(&self, other: &Mol) -> bool {
        chematic_chem::are_identical(&self.inner, &other.inner)
    }

    /// Return a charge-neutralized copy.
    fn neutralize(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::neutralize_charges(&self.inner)),
        }
    }

    /// Return the generic Murcko scaffold (all atoms replaced with carbons, all bonds single).
    fn generic_scaffold(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::generic_murcko_scaffold(&self.inner)),
        }
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

    /// Template-based retrosynthetic disconnection.
    ///
    /// Applies a library of 60 reverse-SMIRKS templates to the molecule and
    /// returns a ranked list of one-step precursor sets.  Each result dict
    /// contains the template name, reaction class, precursor SMILES, and SA
    /// scores (1 = easy, 10 = hard to synthesise).
    ///
    /// Parameters
    /// ----------
    /// max_results : int, optional
    ///     Maximum number of disconnections to return (0 = unlimited).
    ///     Default: 20.
    /// reaction_class : str, optional
    ///     Filter to a single reaction class.  Valid values:
    ///     ``"AmideBond"``, ``"Ester"``, ``"Ether"``, ``"CNBond"``,
    ///     ``"CCBond"``, ``"CSBond"``, ``"Other"``.  Default: all classes.
    ///
    /// Returns
    /// -------
    /// list[dict]
    ///     Each dict has keys:
    ///     - ``template``       : str   — template name (e.g. ``"amide_secondary"``)
    ///     - ``reaction_class`` : str   — reaction class (e.g. ``"AmideBond"``)
    ///     - ``precursors``     : list[str] — canonical SMILES of precursors
    ///     - ``sa_scores``      : list[float] — SA score per precursor
    ///     - ``max_sa_score``   : float — max SA score across precursors
    #[pyo3(signature = (max_results=20, reaction_class=None))]
    fn retro_disconnect<'py>(
        &self,
        py: Python<'py>,
        max_results: usize,
        reaction_class: Option<&str>,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        use chematic_rxn::retro::{DEFAULT_TEMPLATES, RetroClass, retro_disconnect};

        let filter_class: Option<RetroClass> = match reaction_class {
            None => None,
            Some("AmideBond") => Some(RetroClass::AmideBond),
            Some("Ester") => Some(RetroClass::Ester),
            Some("Ether") => Some(RetroClass::Ether),
            Some("CNBond") => Some(RetroClass::CNBond),
            Some("CCBond") => Some(RetroClass::CCBond),
            Some("CSBond") => Some(RetroClass::CSBond),
            Some("Other") => Some(RetroClass::Other),
            Some(other) => {
                return Err(PyValueError::new_err(format!(
                    "unknown reaction_class '{other}'; valid: AmideBond, Ester, Ether, CNBond, CCBond, CSBond, Other"
                )));
            }
        };

        // Filter and collect owned templates.
        let owned: Vec<chematic_rxn::retro::RetroTemplate> = DEFAULT_TEMPLATES
            .iter()
            .filter(|t| filter_class.map(|c| c == t.reaction_class).unwrap_or(true))
            .map(|t| chematic_rxn::retro::RetroTemplate {
                name: t.name,
                smirks: t.smirks,
                reaction_class: t.reaction_class,
            })
            .collect();

        let results = retro_disconnect(&self.inner, &owned, max_results);

        results
            .into_iter()
            .map(|r| {
                let d = PyDict::new(py);
                d.set_item("template", &r.template_name)?;
                d.set_item("reaction_class", r.reaction_class.as_str())?;
                d.set_item("precursors", &r.precursor_smiles)?;

                // Compute SA scores for each precursor using chematic-chem.
                let sa_scores: Vec<f64> =
                    r.precursors.iter().map(chematic_chem::sa_score).collect();
                let max_sa = sa_scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                d.set_item("sa_scores", sa_scores)?;
                d.set_item(
                    "max_sa_score",
                    if max_sa.is_finite() { max_sa } else { 0.0 },
                )?;

                Ok(d)
            })
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

    /// Layered fingerprint decomposed into 7 individual layers.
    ///
    /// Each layer is a 2048-bit (256-byte) fingerprint encoding progressively
    /// more structural detail:
    ///
    /// - Layer 0: raw atom types (element, H count, charge)
    /// - Layer 1: + bond orders
    /// - Layer 2: + aromaticity
    /// - Layer 3: + ring membership
    /// - Layer 4: + is-ring-bond
    /// - Layer 5: + stereochemistry
    /// - Layer 6: all features combined
    ///
    /// Each element is a ``bytes`` object compatible with :func:`tanimoto`,
    /// :func:`dice_similarity`, etc.
    /// Equivalent to RDKit's ``Chem.LayeredFingerprint(mol, layerFlags=0x7F)``.
    ///
    ///     layers = mol.layered_fp_layers()
    ///     sim = chematic.tanimoto(layers[3], other.layered_fp_layers()[3])
    fn layered_fp_layers(&self) -> Vec<Vec<u8>> {
        let layers = chematic_fp::layered_fp_by_layer(&self.inner);
        layers.iter().map(bitvec2048_to_bytes).collect()
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

    /// Blood-brain barrier penetration score (0 = low, 1 = high penetration).
    ///
    /// See :meth:`admet` for the full ADMET profile dict.
    #[getter]
    fn bbb_score(&self) -> f64 {
        chematic_chem::bbb_score(&self.inner)
    }

    /// Predicted Caco-2 intestinal permeability (nm/s).
    ///
    /// Higher values indicate better oral absorption.
    #[getter]
    fn caco2(&self) -> f64 {
        chematic_chem::caco2_permeability(&self.inner)
    }

    /// hERG cardiac toxicity risk score (0 = low, 1 = high risk).
    ///
    /// hERG channel blockade can cause cardiac arrhythmias (QT prolongation).
    #[getter]
    fn herg_risk(&self) -> f64 {
        chematic_chem::herg_risk_score(&self.inner)
    }

    /// CYP3A4 inhibition risk score (0 = low, 1 = high).
    ///
    /// CYP3A4 is the primary drug-metabolising enzyme; inhibition causes drug-drug interactions.
    #[getter]
    fn cyp3a4_risk(&self) -> f64 {
        chematic_chem::cyp3a4_inhibition_risk(&self.inner)
    }

    /// Names of Ames mutagenicity SMARTS alerts matched by this molecule.
    ///
    /// An empty list means the molecule passes all Ames filters (equivalent to
    /// :attr:`ames_passes` returning ``True``).
    /// Complements :meth:`pains_alerts` and :meth:`brenk_alerts`.
    ///
    ///     alerts = mol.ames_alerts()   # ["aromatic_amine", ...]
    fn ames_alerts(&self) -> Vec<String> {
        chematic_chem::ames_alerts(&self.inner)
            .iter()
            .map(|s| s.to_string())
            .collect()
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
            chematic_chem::ClearanceClass::Low => "Low",
            chematic_chem::ClearanceClass::Medium => "Medium",
            chematic_chem::ClearanceClass::High => "High",
        }
    }

    /// Predicted hepatic clearance score (raw float 0.0–1.0).
    ///
    /// Lower = slower clearance; higher = faster.
    /// Complements :attr:`clearance_class` (returns discretised ``"Low"``/``"Medium"``/``"High"``).
    /// Useful for ML pipelines that need a continuous target variable.
    #[getter]
    fn clearance_score(&self) -> f64 {
        chematic_chem::clearance_score(&self.inner)
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
                    chematic_chem::AtropisomerType::Biaryl => "Biaryl",
                    chematic_chem::AtropisomerType::Allene => "Allene",
                    chematic_chem::AtropisomerType::Constrained => "Constrained",
                };
                (bidx.0 as usize, label.to_string())
            })
            .collect()
    }

    // -----------------------------------------------------------------------
    // Fingerprints (extended)
    // -----------------------------------------------------------------------

    /// MAP4 fingerprint (MinHashed Atom-Pair, Minervini 2020) — 1024 u32 hash values.
    ///
    /// Use :func:`tanimoto_map4` for similarity, not the bitwise :func:`tanimoto`.
    fn map4(&self) -> Vec<u32> {
        chematic_fp::map4_default(&self.inner)
    }

    /// MAP4 fingerprint as a numpy array of shape ``(1024,)`` with ``dtype=uint32``.
    fn map4_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        let fp = chematic_fp::map4_default(&self.inner);
        Array1::from_vec(fp).into_pyarray(py)
    }

    /// Extended Reduced Graph (ERG) fingerprint as bytes (256 bytes = 2048 bits).
    fn erg(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::erg(&self.inner).bits)
    }

    // -----------------------------------------------------------------------
    // Descriptors (extended)
    // -----------------------------------------------------------------------

    /// LogD at a given pH — accounts for ionization of basic/acidic groups.
    ///
    /// Default pH is 7.4 (physiological). More predictive than LogP for ADMET.
    #[pyo3(signature = (ph = 7.4))]
    fn logd(&self, ph: f64) -> f64 {
        chematic_chem::logd_simple(&self.inner, ph)
    }

    /// LogD profile — list of ``(pH, LogD)`` pairs from pH 0 to 14 (28 steps).
    fn logd_profile(&self) -> Vec<(f64, f64)> {
        chematic_chem::logd_profile(&self.inner, 0.0, 14.0, 28)
    }

    /// Molecular Quantum Numbers (MQN) — 42-element topological descriptor vector.
    ///
    /// Encodes atom counts, bond counts, ring counts, and degree statistics.
    /// Reference: Ertl et al., *J. Chem. Inf. Model.* 2009.
    fn mqn(&self) -> Vec<u8> {
        chematic_chem::mqn(&self.inner)
    }

    /// Per-atom Crippen LogP contributions — one float per heavy atom, in atom order.
    fn logp_per_atom(&self) -> Vec<f64> {
        chematic_chem::logp_crippen_per_atom(&self.inner)
    }

    /// Per-atom TPSA contributions (Ertl 2000) — one float per heavy atom, in atom order.
    ///
    /// Only N, O, S, P atoms can have non-zero contributions. Sums to :attr:`tpsa`.
    fn tpsa_per_atom(&self) -> Vec<f64> {
        chematic_chem::tpsa_per_atom(&self.inner)
    }

    /// Per-atom hybridization state: 1 = sp, 2 = sp2, 3 = sp3, 0 = other.
    ///
    /// Aromatic atoms → 2, triple-bond atoms → 1, double-bond atoms → 2, otherwise 3.
    /// Useful for scaffold modification (PromptSMILES-style) and atom featurization.
    ///
    ///     mol = chematic.from_smiles("CC=O")
    ///     mol.hybridization_per_atom()  # [3, 2, 2] (CH3=sp3, C=sp2, O=sp2)
    fn hybridization_per_atom(&self) -> Vec<u8> {
        chematic_chem::hybridization_per_atom(&self.inner)
    }

    /// Per-atom formal charge — one ``int`` per heavy atom.
    ///
    /// All values are 0 for neutral molecules. Non-zero for charged atoms ([NH4+] etc.).
    ///
    ///     mol = chematic.from_smiles("[NH4+]")
    ///     mol.formal_charge_per_atom()  # [1]   (N has +1)
    fn formal_charge_per_atom(&self) -> Vec<i8> {
        chematic_chem::formal_charge_per_atom(&self.inner)
    }

    /// Per-atom implicit hydrogen count — one ``int`` per heavy atom.
    ///
    ///     mol = chematic.from_smiles("CC")
    ///     mol.implicit_hcount_per_atom()  # [3, 3]
    fn implicit_hcount_per_atom(&self) -> Vec<u8> {
        chematic_chem::implicit_hcount_per_atom(&self.inner)
    }

    /// Isotopic distribution — list of ``(mass, relative_intensity)`` pairs.
    ///
    /// The highest-intensity peak is normalised to 1.0.
    fn isotope_distribution(&self) -> Vec<(f64, f64)> {
        chematic_chem::isotope_distribution(&self.inner, 0.001)
    }

    // -----------------------------------------------------------------------
    // Chemical analysis (extended)
    // -----------------------------------------------------------------------

    /// Identify functional groups (Ertl 2017 algorithm).
    ///
    /// Returns a list of dicts, each with:
    ///   ``atom_indices`` (list of int), ``atom_types`` (str, e.g. ``"N,O"``).
    fn functional_groups<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        chematic_chem::identify_functional_groups(&self.inner)
            .into_iter()
            .map(|fg| {
                let d = PyDict::new(py);
                d.set_item("atom_indices", fg.atom_indices)?;
                d.set_item("atom_types", fg.atom_types)?;
                Ok(d)
            })
            .collect()
    }

    /// Schuffenhauer scaffold parents — list of SMILES (outermost scaffold first).
    ///
    /// Each entry is a simpler scaffold obtained by removing one ring. Useful for
    /// scaffold-hopping and SAR analysis.
    fn scaffold_network(&self) -> Vec<String> {
        chematic_chem::schuffenhauer_parents(&self.inner)
            .into_iter()
            .map(|m| chematic_smiles::canonical_smiles(&m))
            .collect()
    }

    // -----------------------------------------------------------------------
    // 3D coordinate generation
    // -----------------------------------------------------------------------

    /// Generate 3D coordinates using distance geometry + DREIDING minimization.
    ///
    /// Returns a list of ``[x, y, z]`` lists (Å), one per heavy atom.
    /// Use the returned coords with :meth:`whim`, :meth:`getaway`,
    /// :meth:`mmff94_energy_breakdown`, :meth:`to_pdb`, etc.
    fn generate_3d(&self) -> Vec<Vec<f64>> {
        let coords = chematic_3d::generate_and_minimize_dreiding(&self.inner);
        coords.points.iter().map(|p| vec![p.x, p.y, p.z]).collect()
    }

    /// Generate multiple conformers with RMSD-based pruning.
    ///
    /// Returns a list of coordinate arrays — each is a ``[[x,y,z], ...]`` list.
    ///
    /// Args:
    ///     n: Number of conformers to attempt.
    ///     rmsd_threshold: Minimum RMSD (Å) between conformers (default 0.5).
    /// Generate a conformer ensemble using ETKDG + force-field minimization + RMSD pruning.
    ///
    /// Args:
    ///     n: Number of conformers to attempt.
    ///     rmsd_threshold: Minimum Kabsch-aligned RMSD (Å) between retained conformers (default 0.5).
    ///     force_field: ``"dreiding"`` (fast, default) or ``"mmff94"`` (higher accuracy).
    ///     noise_sigma_deg: Gaussian torsion noise σ in degrees (default 30.0).
    #[pyo3(signature = (n, rmsd_threshold = 0.5, force_field = "dreiding", noise_sigma_deg = 30.0))]
    fn conformer_ensemble(
        &self,
        n: usize,
        rmsd_threshold: f64,
        force_field: &str,
        noise_sigma_deg: f64,
    ) -> Vec<Vec<Vec<f64>>> {
        let smiles = chematic_smiles::canonical_smiles(&self.inner);
        let mol = match chematic_smiles::parse(&smiles) {
            Ok(m) => m,
            Err(_) => return Vec::new(),
        };
        let ff = if force_field.eq_ignore_ascii_case("mmff94") {
            chematic_3d::ConformerForceField::Mmff94
        } else {
            chematic_3d::ConformerForceField::Dreiding
        };
        let config = chematic_3d::ConformerConfig {
            count: n,
            rmsd_threshold,
            force_field: ff,
            noise_sigma_deg,
        };
        match chematic_3d::generate_conformer_ensemble_with_config(mol, &config) {
            Ok(ensemble) => (0..ensemble.conformer_count())
                .filter_map(|i| ensemble.get_conformer(i))
                .map(|c3d| c3d.points.iter().map(|p| vec![p.x, p.y, p.z]).collect())
                .collect(),
            Err(_) => Vec::new(),
        }
    }

    // -----------------------------------------------------------------------
    // 3D descriptors
    // -----------------------------------------------------------------------

    /// WHIM 3D descriptors (Todeschini & Gramatica 1997).
    ///
    /// ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
    /// Returns a flat list of floats (shape/symmetry descriptors).
    fn whim(&self, coords: Vec<[f64; 3]>) -> Vec<f64> {
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::whim_descriptors(&self.inner, &c3d)
    }

    /// GETAWAY 3D descriptors (Consonni et al. 2002).
    ///
    /// ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
    /// Returns a flat list of floats.
    fn getaway(&self, coords: Vec<[f64; 3]>) -> Vec<f64> {
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::getaway_descriptors(&self.inner, &c3d)
    }

    /// 3D autocorrelation descriptors.
    ///
    /// ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
    fn autocorr_3d(&self, coords: Vec<[f64; 3]>) -> Vec<f64> {
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::autocorr_3d(&self.inner, &c3d)
    }

    /// Spectrophores 3D fingerprint — 48-element vector.
    ///
    /// Encodes the electrostatic, lipophilic, aromatic, and H-bond character
    /// of the molecule's 3D surface into a fixed-size numerical vector suitable
    /// for 3D QSAR, shape-based screening, and virtual screening.
    ///
    /// Requires 3D coordinates (one ``[x, y, z]`` per heavy atom, in Å).
    /// Use :meth:`generate_3d` to obtain coordinates.
    ///
    /// Returns a list of 48 floats organised as four blocks of 12 probe values:
    /// electrostatic, lipophilic, aromatic, H-bond (in that order).
    ///
    /// Reference: Silicos-it Spectrophores (patent expired 2024).
    ///
    ///     coords = mol.generate_3d()
    ///     fp = mol.spectrophores(coords)          # len == 48
    ///     fp_z = mol.spectrophores(coords, normalize="zscore")
    ///     sim = chematic.tanimoto_spectrophores(fp1, fp2)
    #[pyo3(signature = (coords, normalize = "none"))]
    fn spectrophores(&self, coords: Vec<[f64; 3]>, normalize: &str) -> Vec<f64> {
        let c3d = flat_to_coords3d(&coords);
        let norm = match normalize.to_lowercase().as_str() {
            "zscore" | "z-score" => chematic_3d::SpectrophoresNorm::ZScore,
            "l2" => chematic_3d::SpectrophoresNorm::L2,
            _ => chematic_3d::SpectrophoresNorm::None,
        };
        let config = chematic_3d::SpectrophoresConfig {
            normalize: norm,
            ..Default::default()
        };
        chematic_3d::spectrophores(&self.inner, &c3d, &config)
    }

    // -----------------------------------------------------------------------
    // 3D file I/O
    // -----------------------------------------------------------------------

    /// Write this molecule's 3D structure to PDB format.
    ///
    /// ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
    fn to_pdb(&self, coords: Vec<[f64; 3]>) -> String {
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::write_pdb(&self.inner, &c3d)
    }

    /// Write this molecule's 3D structure to XYZ format.
    ///
    /// ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
    #[pyo3(signature = (coords, comment = ""))]
    fn to_xyz(&self, coords: Vec<[f64; 3]>, comment: &str) -> String {
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::write_xyz(&self.inner, &c3d, comment)
    }

    // -----------------------------------------------------------------------
    // Force field analysis
    // -----------------------------------------------------------------------

    /// MMFF94 energy breakdown for given 3D coordinates.
    ///
    /// Returns a dict with keys: ``bond``, ``angle``, ``stretch_bend``,
    /// ``torsion``, ``oop``, ``vdw``, ``electrostatic``, ``total`` (kcal/mol).
    ///
    /// ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
    /// Raises ``ValueError`` for atoms not parameterised by MMFF94.
    fn mmff94_energy_breakdown<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let b = chematic_ff::mmff94_energy_breakdown(&self.inner, &coords)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let d = PyDict::new(py);
        d.set_item("bond", b.bond)?;
        d.set_item("angle", b.angle)?;
        d.set_item("stretch_bend", b.stretch_bend)?;
        d.set_item("torsion", b.torsion)?;
        d.set_item("oop", b.oop)?;
        d.set_item("vdw", b.vdw)?;
        d.set_item("electrostatic", b.electrostatic)?;
        d.set_item("total", b.total)?;
        Ok(d)
    }

    /// Scan a torsion dihedral (atoms i–j–k–l) over 360° in ``steps`` increments.
    ///
    /// Returns a list of ``(angle_deg, energy_kcal)`` pairs. The molecule is not modified.
    ///
    /// ``coords``: ``[[x,y,z], ...]`` list (Å), one per heavy atom.
    /// ``steps``: number of scan points (default 36 = 10° per step).
    #[pyo3(signature = (coords, atom_i, atom_j, atom_k, atom_l, steps = 36))]
    fn mmff94_torsion_scan(
        &self,
        coords: Vec<[f64; 3]>,
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        atom_l: usize,
        steps: usize,
    ) -> PyResult<Vec<(f64, f64)>> {
        chematic_ff::mmff94_torsion_scan(
            &self.inner,
            &coords,
            atom_i,
            atom_j,
            atom_k,
            atom_l,
            steps,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    // -----------------------------------------------------------------------
    // Sprint 2: additional descriptors & fingerprints
    // -----------------------------------------------------------------------

    /// XLogP3 — alternative logP (Wang et al. 2000).
    #[getter]
    fn xlogp3(&self) -> f64 {
        chematic_chem::xlogp3(&self.inner)
    }

    /// Per-atom XLogP3 contributions — one float per heavy atom.
    fn xlogp3_per_atom(&self) -> Vec<f64> {
        chematic_chem::xlogp3_per_atom(&self.inner)
    }

    /// 2D autocorrelation descriptors — flat list of floats (mordred/Dragon compatible).
    fn autocorr_2d(&self) -> Vec<f64> {
        chematic_chem::autocorr_2d(&self.inner)
    }

    /// Hall-Kier alpha — correction term for kappa shape indices.
    #[getter]
    fn hall_kier_alpha(&self) -> f64 {
        chematic_chem::hall_kier_alpha(&self.inner)
    }

    /// PEOE VSA — bins of van der Waals surface area by partial charge.
    fn peoe_vsa(&self) -> Vec<f64> {
        chematic_chem::peoe_vsa(&self.inner)
    }

    /// SLogP VSA — bins of van der Waals surface area by SLogP contribution.
    fn slogp_vsa(&self) -> Vec<f64> {
        chematic_chem::slogp_vsa(&self.inner)
    }

    /// SMR VSA — bins of van der Waals surface area by molar refractivity.
    fn smr_vsa(&self) -> Vec<f64> {
        chematic_chem::smr_vsa(&self.inner)
    }

    /// EState VSA — bins of van der Waals surface area by E-state index.
    fn estate_vsa(&self) -> Vec<f64> {
        chematic_chem::estate_vsa(&self.inner)
    }

    /// USRCAT — 42-element topological shape/pharmacophore descriptor.
    fn usrcat(&self) -> Vec<f64> {
        chematic_chem::usrcat(&self.inner).to_vec()
    }

    /// Write this molecule as a 3D SDF record with partial charges.
    ///
    /// ``coords``: ``[[x,y,z], ...]`` (Å), one per heavy atom.
    /// ``charges``: partial charges, one per heavy atom.
    ///
    /// The charges are written as a ``> <PARTIAL_CHARGES>`` SD data field.
    fn to_sdf_with_charges(&self, coords: Vec<[f64; 3]>, charges: Vec<f64>) -> String {
        use chematic_core::BondOrder;
        let mol = &self.inner;
        let mut out = String::new();
        out.push_str("\n  chematic\n\n");
        let natoms = mol.atom_count();
        let bonds: Vec<_> = mol.bonds().collect();
        let nbonds = bonds.len();
        out.push_str(&format!(
            "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
            natoms, nbonds
        ));
        for (idx, atom) in mol.atoms() {
            let i = idx.0 as usize;
            let (x, y, z) = if i < coords.len() {
                (coords[i][0], coords[i][1], coords[i][2])
            } else {
                (0.0, 0.0, 0.0)
            };
            let sym = atom.element.symbol();
            out.push_str(&format!(
                "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
                x, y, z, sym
            ));
        }
        for (_, bond) in &bonds {
            let a1 = bond.atom1.0 + 1;
            let a2 = bond.atom2.0 + 1;
            let btype = match bond.order {
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Aromatic => 4,
                _ => 1,
            };
            out.push_str(&format!("{:>3}{:>3}{:>3}  0\n", a1, a2, btype));
        }
        out.push_str("M  END\n");
        if !charges.is_empty() {
            out.push_str("> <PARTIAL_CHARGES>\n");
            let vals: Vec<String> = charges.iter().map(|q| format!("{q:.4}")).collect();
            out.push_str(&vals.join(" "));
            out.push_str("\n\n");
        }
        out.push_str("$$$$\n");
        out
    }

    /// Pharmacophore 2D fingerprint as bytes (256 bytes = 2048 bits).
    ///
    /// Compatible with :func:`tanimoto`. Encodes HBD, HBA, hydrophobic,
    /// aromatic, and positive/negative pharmacophore features.
    fn pharmacophore_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::pharmacophore_fp_2d(&self.inner))
    }

    /// MHFP (MinHash Fingerprint) — 128 u64 hash values.
    ///
    /// Use :func:`tanimoto_mhfp` for similarity (position-wise, not bitwise).
    fn mhfp(&self) -> Vec<u64> {
        chematic_fp::mhfp_128(&self.inner).hashes
    }

    /// MinHash fingerprint with custom parameters.
    ///
    /// Returns a list of ``num_hashes`` unsigned 64-bit integers.
    /// Use :func:`tanimoto_mhfp` for similarity comparison.
    ///
    /// Args:
    ///     radius: circular subgraph radius (default 2, equivalent to ECFP4 radius)
    ///     num_hashes: fingerprint length in hash slots (default 128)
    ///     seed: hash seed for reproducibility (default 0)
    ///
    ///     fp = mol.mhfp_config(radius=3, num_hashes=256)
    ///     sim = chematic.tanimoto_mhfp(fp, other.mhfp_config(radius=3, num_hashes=256))
    #[pyo3(signature = (radius = 2, num_hashes = 128, seed = 0))]
    fn mhfp_config(&self, radius: u32, num_hashes: usize, seed: u64) -> Vec<u64> {
        let config = chematic_fp::MhfpConfig {
            radius,
            num_hashes,
            seed,
        };
        chematic_fp::mhfp_with_config(&self.inner, &config).hashes
    }

    // -----------------------------------------------------------------------
    // Sprint 3: charges, structure analysis, 3D shape, conformer tools
    // -----------------------------------------------------------------------

    /// Gasteiger–Marsili partial charges — one float per heavy atom.
    ///
    /// These are the standard open-source partial charges for docking prep.
    /// Use together with :meth:`to_pdbqt` to complete the docking pipeline::
    ///
    ///     coords = mol.generate_3d()
    ///     charges = mol.gasteiger_charges()
    ///     pdbqt = mol.to_pdbqt(coords, charges, "LIG")
    fn gasteiger_charges(&self) -> Vec<f64> {
        chematic_chem::gasteiger_charges(&self.inner)
    }

    /// CIP stereochemistry assignments — list of ``{"atom_idx": int, "descriptor": str}`` dicts.
    ///
    /// ``descriptor`` is ``"R"``, ``"S"``, ``"E"``, or ``"Z"``.
    /// Only assigned stereocenters / double bonds are returned.
    fn cip_stereo<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        use chematic_core::CipCode;
        chematic_chem::assign_cip(&self.inner)
            .assignments
            .iter()
            .map(|(idx, code)| {
                let d = PyDict::new(py);
                d.set_item("atom_idx", idx.0 as usize)?;
                let label = match code {
                    CipCode::R => "R",
                    CipCode::S => "S",
                    CipCode::E => "E",
                    CipCode::Z => "Z",
                };
                d.set_item("descriptor", label)?;
                Ok(d)
            })
            .collect()
    }

    /// Names of PAINS structural alerts matched by this molecule.
    ///
    /// Returns an empty list when no alerts match (i.e. same as :attr:`pains_passes` == True).
    fn pains_alerts(&self) -> Vec<String> {
        chematic_chem::pains_matches(&self.inner)
            .iter()
            .map(|s| s.to_string())
            .collect()
    }

    /// Names of Brenk structural alerts (instability / toxicity) matched by this molecule.
    fn brenk_alerts(&self) -> Vec<String> {
        chematic_chem::brenk_matches(&self.inner)
            .iter()
            .map(|s| s.to_string())
            .collect()
    }

    /// PAINS alerts with matched heavy-atom indices for substructure highlighting.
    ///
    /// Returns ``list[tuple[str, list[int]]]``: each entry is
    /// ``(alert_name, [atom_idx, ...])`` where the atom indices refer to
    /// heavy atoms in this molecule.  Use the indices with the ``highlight``
    /// parameter of :meth:`svg` or pass directly to
    /// ``chematic.depict.render_svg_highlighted``.
    ///
    ///     for name, atoms in mol.pains_alerts_detailed():
    ///         print(f"{name}: atoms {atoms}")
    fn pains_alerts_detailed(&self) -> Vec<(String, Vec<usize>)> {
        chematic_chem::pains_matches_detailed(&self.inner)
            .into_iter()
            .map(|(name, idxs)| {
                (
                    name.to_string(),
                    idxs.into_iter().map(|a| a.0 as usize).collect(),
                )
            })
            .collect()
    }

    /// Brenk structural alerts with matched heavy-atom indices for highlighting.
    ///
    /// Returns ``list[tuple[str, list[int]]]`` — same format as
    /// :meth:`pains_alerts_detailed`.
    ///
    ///     for name, atoms in mol.brenk_alerts_detailed():
    ///         print(f"{name}: atoms {atoms}")
    fn brenk_alerts_detailed(&self) -> Vec<(String, Vec<usize>)> {
        chematic_chem::brenk_matches_detailed(&self.inner)
            .into_iter()
            .map(|(name, idxs)| {
                (
                    name.to_string(),
                    idxs.into_iter().map(|a| a.0 as usize).collect(),
                )
            })
            .collect()
    }

    /// SVG depiction with structural alert atoms highlighted in red.
    ///
    /// Combines PAINS and Brenk alerts.  Atoms belonging to flagged
    /// substructures are coloured red; all others use the standard
    /// CPK colour scheme.
    ///
    ///     svg_str = mol.svg_with_alerts()
    ///     svg_str = mol.svg_with_alerts(width=600, height=400)
    #[pyo3(signature = (width = None, height = None))]
    fn svg_with_alerts(&self, width: Option<u32>, height: Option<u32>) -> String {
        use chematic_core::AtomIdx;

        let pains = chematic_chem::pains_matches_detailed(&self.inner);
        let brenk = chematic_chem::brenk_matches_detailed(&self.inner);

        let mut opts = chematic_depict::RenderOptions {
            width,
            height,
            highlight_color: "#FF4444".to_string(),
            ..Default::default()
        };
        for atom in pains.into_iter().chain(brenk).flat_map(|(_, atoms)| atoms) {
            opts.highlight_atoms.insert(AtomIdx(atom.0));
        }
        chematic_depict::depict_svg_opts(&self.inner, &opts)
    }

    /// Named functional groups detected in this molecule — list of group names.
    ///
    /// Names include: hydroxyl, carbonyl, carboxyl, aldehyde, ketone, amine,
    /// amide, ester, nitrile, sulfide, sulfonyl, phosphate, and more.
    fn named_functional_groups(&self) -> Vec<String> {
        chematic_chem::detect_named_functional_groups(&self.inner)
            .into_iter()
            .map(|g| g.name.to_string())
            .collect()
    }

    // ---- 3D shape descriptors ----

    /// Principal Moments of Inertia — ``[PMI1, PMI2, PMI3]`` (ascending order).
    ///
    /// ``coords``: ``[[x,y,z], ...]`` (Å), one per heavy atom.
    /// Used to classify shape as rod / disc / sphere in PMI plots.
    fn pmi(&self, coords: Vec<[f64; 3]>) -> [f64; 3] {
        let c3d = flat_to_coords3d(&coords);
        let (p1, p2, p3) = chematic_3d::pmi(&self.inner, &c3d);
        [p1, p2, p3]
    }

    /// Normalised Principal Moments — ``[NPR1, NPR2]``.
    ///
    /// NPR1 = PMI1/PMI3, NPR2 = PMI2/PMI3. Values in [0, 1].
    /// Rod: NPR1≈0, NPR2≈0.5; Sphere: NPR1≈NPR2≈1; Disc: NPR1≈0.5, NPR2≈1.
    fn npr(&self, coords: Vec<[f64; 3]>) -> [f64; 2] {
        let c3d = flat_to_coords3d(&coords);
        [
            chematic_3d::npr1(&self.inner, &c3d),
            chematic_3d::npr2(&self.inner, &c3d),
        ]
    }

    /// Asphericity — deviation from a perfect sphere. Range [0, 1].
    fn asphericity(&self, coords: Vec<[f64; 3]>) -> f64 {
        chematic_3d::asphericity(&self.inner, &flat_to_coords3d(&coords))
    }

    /// Eccentricity — elongation measure. Range [0, 1].
    fn eccentricity(&self, coords: Vec<[f64; 3]>) -> f64 {
        chematic_3d::eccentricity(&self.inner, &flat_to_coords3d(&coords))
    }

    /// Radius of gyration (Å).
    fn radius_of_gyration(&self, coords: Vec<[f64; 3]>) -> f64 {
        chematic_3d::radius_of_gyration(&self.inner, &flat_to_coords3d(&coords))
    }

    /// Plane of Best Fit (PBF) — deviation of atoms from the least-squares plane (Å).
    fn plane_of_best_fit(&self, coords: Vec<[f64; 3]>) -> f64 {
        chematic_3d::plane_of_best_fit(&self.inner, &flat_to_coords3d(&coords))
    }

    // ---- Conformer editing tools ----

    /// Generate 3D coordinates using the ETKDG algorithm (higher quality than basic DG).
    ///
    /// Returns ``[[x,y,z], ...]`` (Å), one per heavy atom.
    fn generate_3d_etkdg(&self) -> Vec<Vec<f64>> {
        let c3d = chematic_3d::generate_coords_etkdg(&self.inner);
        c3d.points.iter().map(|p| vec![p.x, p.y, p.z]).collect()
    }

    /// Measure dihedral angle i–j–k–l in degrees.
    ///
    /// Returns ``None`` when three atoms are collinear (undefined dihedral).
    fn get_dihedral(
        &self,
        coords: Vec<[f64; 3]>,
        i: usize,
        j: usize,
        k: usize,
        l: usize,
    ) -> Option<f64> {
        use chematic_core::AtomIdx;
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::get_dihedral_deg(
            &c3d,
            AtomIdx(i as u32),
            AtomIdx(j as u32),
            AtomIdx(k as u32),
            AtomIdx(l as u32),
        )
    }

    /// Set dihedral angle i–j–k–l to ``angle_deg`` degrees.
    ///
    /// Returns new ``[[x,y,z], ...]`` coordinates with atoms on the l-side
    /// of the j–k bond rotated to the target angle.
    fn set_dihedral(
        &self,
        coords: Vec<[f64; 3]>,
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        angle_deg: f64,
    ) -> Vec<Vec<f64>> {
        use chematic_core::AtomIdx;
        let c3d = flat_to_coords3d(&coords);
        let new_c3d = chematic_3d::set_dihedral(
            &c3d,
            &self.inner,
            AtomIdx(i as u32),
            AtomIdx(j as u32),
            AtomIdx(k as u32),
            AtomIdx(l as u32),
            angle_deg.to_radians(),
        );
        new_c3d.points.iter().map(|p| vec![p.x, p.y, p.z]).collect()
    }

    /// Measure bond length between atoms i and j (Å).
    fn get_bond_length(&self, coords: Vec<[f64; 3]>, i: usize, j: usize) -> f64 {
        use chematic_core::AtomIdx;
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::get_bond_length(&c3d, AtomIdx(i as u32), AtomIdx(j as u32))
    }

    /// Measure bond angle i–j–k in degrees (j is the central atom).
    fn get_bond_angle(&self, coords: Vec<[f64; 3]>, i: usize, j: usize, k: usize) -> f64 {
        use chematic_core::AtomIdx;
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::get_bond_angle_deg(
            &c3d,
            AtomIdx(i as u32),
            AtomIdx(j as u32),
            AtomIdx(k as u32),
        )
    }

    // -----------------------------------------------------------------------
    // Sprint 4: additional fingerprints
    // -----------------------------------------------------------------------

    /// RDKit-style Daylight path fingerprint as bytes (256 bytes = 2048 bits).
    ///
    /// Equivalent to RDKit's ``RDKFingerprint()``. Compatible with :func:`tanimoto`.
    fn path_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::rdkit_path_fp(&self.inner))
    }

    /// Topological path fingerprint as bytes (256 bytes = 2048 bits).
    fn topo_path_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::topo_path(
            &self.inner,
            &chematic_fp::TopoPathConfig::default(),
        ))
    }

    /// 3D pharmacophore fingerprint as bytes (256 bytes = 2048 bits).
    ///
    /// Requires 3D coordinates. Complements the 2D version from :meth:`pharmacophore_fp`.
    ///
    /// ``coords``: ``[[x,y,z], ...]`` (Å), one per heavy atom.
    fn pharmacophore_fp_3d(&self, coords: Vec<[f64; 3]>) -> Vec<u8> {
        let c3d = flat_to_coords3d(&coords);
        bitvec2048_to_bytes(&chematic_3d::pharmacophore_fp_3d(&self.inner, &c3d))
    }

    /// Reaction fingerprint — ``(reactant_fp, product_fp, combined_fp)`` as bytes.
    ///
    /// Each component is 256 bytes (2048 bits). The combined FP captures the full
    /// transformation and is suitable for reaction similarity search.
    ///
    /// Args:
    ///     reaction_smiles: Reaction SMILES string ``"R>>P"`` or ``"R>A>P"``.
    ///
    ///     rfp = chematic.from_smiles("CCO").reaction_fp("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
    fn reaction_fp<'py>(
        &self,
        py: Python<'py>,
        reaction_smiles: &str,
    ) -> PyResult<Bound<'py, PyDict>> {
        let rxn = chematic_rxn::parse_reaction(reaction_smiles)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let rfp = chematic_fp::reaction_fp(&rxn);
        let d = PyDict::new(py);
        d.set_item("reactant_fp", bitvec2048_to_bytes(&rfp.reactant_fp))?;
        d.set_item("product_fp", bitvec2048_to_bytes(&rfp.product_fp))?;
        d.set_item("combined_fp", bitvec2048_to_bytes(&rfp.combined_fp))?;
        Ok(d)
    }

    // -----------------------------------------------------------------------
    // Sprint 5: standardization steps, ADMET, SASA, FP, topology
    // -----------------------------------------------------------------------

    /// Normalize functional groups (nitro, azide, diazo, sulfoxide, etc.).
    ///
    /// Converts charge-separated forms like ``[N+](=O)[O-]`` → ``N(=O)=O``.
    /// Returns a new :class:`Mol` with normalized groups.
    fn normalize_groups(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::normalize_groups(&self.inner)),
        }
    }

    /// Keep the largest organic fragment; remove inorganic counterions.
    ///
    /// Useful after salt removal when the largest fragment is the drug molecule.
    fn prefer_organic(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::prefer_organic(&self.inner)),
        }
    }

    /// Re-apply ionization rules based on pKa.
    ///
    /// Transfers protons to maximize negative charge on the strongest acids.
    fn reionize(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::reionize(&self.inner)),
        }
    }

    /// Remove all formal charges (set every atom to neutral).
    fn uncharge(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::uncharge(&self.inner)),
        }
    }

    /// Predicted hepatic clearance class: ``"Low"``, ``"Medium"``, or ``"High"``.
    #[getter]
    fn clearance_class(&self) -> String {
        match chematic_chem::clearance_class(&self.inner) {
            chematic_chem::ClearanceClass::Low => "Low".to_string(),
            chematic_chem::ClearanceClass::Medium => "Medium".to_string(),
            chematic_chem::ClearanceClass::High => "High".to_string(),
        }
    }

    /// SASA (Å²) with a custom probe radius (default 1.4 Å = water probe).
    ///
    /// ``coords``: ``[[x,y,z], ...]`` (Å), one per heavy atom.
    #[pyo3(signature = (coords, probe_radius = 1.4))]
    fn sasa_with_probe(&self, coords: Vec<[f64; 3]>, probe_radius: f64) -> f64 {
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::calc_mol_sasa_with_probe(&self.inner, &c3d, probe_radius)
    }

    /// Per-atom Solvent-Accessible Surface Area (Å²) from explicit 3D coordinates.
    ///
    /// Uses the Shrake-Rupley algorithm (probe 1.4 Å, 100 sphere points/atom).
    /// Unlike :meth:`sasa_per_atom` (which generates DG coords internally), this method
    /// accepts coordinates from :meth:`generate_3d`, :meth:`minimize_mmff94`, etc.
    ///
    ///     coords = mol.generate_3d()
    ///     per_atom = mol.sasa_per_atom_3d(coords)
    ///     for i, sa in enumerate(per_atom):
    ///         print(f"atom {i}: {sa:.2f} Å²")
    fn sasa_per_atom_3d(&self, coords: Vec<[f64; 3]>) -> Vec<f64> {
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::sasa_per_atom(&self.inner, &c3d, 1.4, 100)
    }

    /// Per-element SASA breakdown — dict mapping element symbol to SASA (Å²).
    ///
    /// ``coords``: ``[[x,y,z], ...]`` (Å), one per heavy atom.
    fn sasa_per_element<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
    ) -> PyResult<Bound<'py, PyDict>> {
        use chematic_core::Element;
        let c3d = flat_to_coords3d(&coords);
        let per_elem = chematic_3d::sasa_per_element(&self.inner, &c3d);
        let d = PyDict::new(py);
        for (idx, val) in per_elem.by_element.iter().enumerate() {
            if *val > 0.0 {
                let sym = Element::from_atomic_number(idx as u8)
                    .map(|e| e.symbol().to_string())
                    .unwrap_or_else(|| idx.to_string());
                d.set_item(sym, val)?;
            }
        }
        Ok(d)
    }

    /// FCFP6 (Functional-Class ECFP6) fingerprint as bytes (256 bytes = 2048 bits).
    fn fcfp6(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::fcfp6(&self.inner))
    }

    /// Pattern fingerprint as bytes (256 bytes = 2048 bits).
    ///
    /// SMARTS-pattern-based structural fingerprint. Compatible with :func:`tanimoto`.
    fn pattern_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::pattern_fp(&self.inner))
    }

    /// Randić connectivity index — measures molecular branching.
    #[getter]
    fn randic_index(&self) -> f64 {
        chematic_chem::randic_index(&self.inner)
    }

    // -----------------------------------------------------------------------
    // Sprint 6: MMFF94 charges, topology, SASA stats
    // -----------------------------------------------------------------------

    /// MMFF94 BCI partial charges — one float per heavy atom.
    ///
    /// More accurate than :meth:`gasteiger_charges` for organic drug-like molecules.
    /// Use with :meth:`to_pdbqt` for best docking accuracy.
    fn mmff94_charges(&self) -> Vec<f64> {
        chematic_chem::mmff94_charges_bci(&self.inner)
    }

    /// Balaban J topological complexity index.
    ///
    /// High J → more complex/branched; related to molecular uniqueness.
    #[getter]
    fn balaban_j(&self) -> f64 {
        chematic_chem::balaban_j(&self.inner)
    }

    /// Information-theoretic connectivity index (IPC).
    #[getter]
    fn ipc(&self) -> f64 {
        chematic_chem::ipc(&self.inner)
    }

    /// Zagreb M1 index — sum of squared vertex degrees.
    #[getter]
    fn zagreb_m1(&self) -> u32 {
        chematic_chem::zagreb_index_m1(&self.inner)
    }

    /// Second Zagreb index M2 — Σ(deg(a) × deg(b)) over all heavy-atom bonds.
    ///
    /// Complements :attr:`zagreb_m1` (Σ deg(v)²); both quantify molecular branching.
    ///
    ///     ethane = chematic.from_smiles("CC")
    ///     ethane.zagreb_m2  # 1 (one bond between degree-1 atoms)
    #[getter]
    fn zagreb_m2(&self) -> u32 {
        chematic_chem::zagreb_index_m2(&self.inner)
    }

    /// Per-atom Labute approximate surface area contributions.
    fn labute_asa_per_atom(&self) -> Vec<f64> {
        chematic_chem::labute_asa_per_atom(&self.inner)
    }

    /// Per-atom molar refractivity contributions.
    ///
    /// Returns one float per heavy atom. The sum equals :attr:`molar_refractivity`.
    /// Useful for fragment-based drug design and QSAR feature generation.
    ///
    ///     mr = mol.mr_per_atom()   # [2.5, 1.1, ...]
    ///     assert abs(sum(mr) - mol.molar_refractivity) < 0.01
    fn mr_per_atom(&self) -> Vec<f64> {
        chematic_chem::mr_per_atom(&self.inner)
    }

    /// MMFF94 partial charges incorporating 3D polarization effects.
    ///
    /// Requires 3D coordinates (from :meth:`generate_3d`).
    /// Complements :meth:`mmff94_charges` (2D topology only).
    /// Returns one float per heavy atom.
    ///
    ///     coords = mol.generate_3d()
    ///     charges = mol.mmff94_charges_3d(coords)
    ///
    /// Raises ``ValueError`` if MMFF94 typing fails or coords length mismatches.
    fn mmff94_charges_3d(&self, coords: Vec<[f64; 3]>) -> PyResult<Vec<f64>> {
        let tuple_coords: Vec<(f64, f64, f64)> =
            coords.iter().map(|[x, y, z]| (*x, *y, *z)).collect();
        chematic_ff::mmff94_charges_3d(&self.inner, &tuple_coords)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// SASA statistics descriptor.
    ///
    /// Returns a dict with keys: ``total`` (Å²), ``mean`` (Å²/atom),
    /// ``std_dev``, ``per_atom`` (list of per-atom values).
    ///
    /// More informative than :meth:`sasa` (total only).
    ///
    /// ``coords``: ``[[x,y,z], ...]`` (Å), one per heavy atom.
    fn sasa_descriptor<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let c3d = flat_to_coords3d(&coords);
        let sd = chematic_3d::sasa_descriptor(&self.inner, &c3d);
        let d = PyDict::new(py);
        d.set_item("total", sd.total)?;
        d.set_item("mean", sd.mean)?;
        d.set_item("std_dev", sd.std_dev)?;
        d.set_item("per_atom", sd.per_atom)?;
        Ok(d)
    }

    // -----------------------------------------------------------------------
    // Sprint 8: BRICS bonds, full pKa, topology distance matrix
    // -----------------------------------------------------------------------

    /// Identify BRICS-breakable bonds — list of ``(atom_i, atom_j)`` index pairs.
    ///
    /// Complements :meth:`brics_fragments` which returns the fragment molecules.
    /// Use the bond indices to identify SAR hotspots without fragmenting.
    ///
    /// Equivalent to RDKit's ``BRICS.FindBRICSBonds(mol)``.
    fn brics_bonds(&self) -> Vec<(usize, usize)> {
        chematic_chem::brics_bonds(&self.inner)
            .into_iter()
            .map(|(a, b)| (a.0 as usize, b.0 as usize))
            .collect()
    }

    /// Predict all ionizable sites and their pKa values.
    ///
    /// Returns a list of dicts, each with keys:
    ///   ``atom_idx`` (int), ``pka`` (float), ``site_type`` (``"Acid"`` or ``"Base"``),
    ///   ``group_name`` (str, e.g. ``"carboxylic_acid"``).
    ///
    /// Complements :attr:`pka_acid` and :attr:`pka_base` which return only the
    /// strongest single site.
    fn predict_pka<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        chematic_chem::predict_pka(&self.inner)
            .into_iter()
            .map(|site| {
                let d = PyDict::new(py);
                d.set_item("atom_idx", site.atom_idx.0 as usize)?;
                d.set_item("pka", site.pka)?;
                let site_type = match site.site_type {
                    chematic_chem::PkaSiteType::Acid => "Acid",
                    chematic_chem::PkaSiteType::Base => "Base",
                };
                d.set_item("site_type", site_type)?;
                d.set_item("group_name", site.group_name)?;
                Ok(d)
            })
            .collect()
    }

    /// Topological (graph) distance matrix — N×N shortest-path distances.
    ///
    /// Returns a list of N lists of N integers (atom graph distances).
    /// Equivalent to RDKit's ``rdmolops.GetDistanceMatrix(mol)``.
    ///
    /// Note: O(N²) memory — use with caution on large molecules (>200 atoms).
    fn topological_distance_matrix(&self) -> Vec<Vec<u32>> {
        chematic_chem::topological_distance_matrix(&self.inner)
    }

    // -----------------------------------------------------------------------
    // Sprint 9: estate_indices, minimize_mmff94, num_unspecified_stereocenters, whim_getaway
    // -----------------------------------------------------------------------

    /// Per-atom E-state electrotopological indices (Kier & Hall).
    ///
    /// Returns one float per heavy atom. Equivalent to RDKit's
    /// ``EState.EState.EStateIndices(mol)``.
    ///
    ///     idx = mol.estate_indices()   # [1.5, 0.2, ...]
    fn estate_indices(&self) -> Vec<f64> {
        chematic_chem::estate_indices(&self.inner)
    }

    /// Minimize 3D coordinates with the MMFF94 force field.
    ///
    /// Complements :meth:`minimize_uff`. Returns minimized coords as
    /// ``[[x,y,z], ...]``.
    ///
    ///     coords = mol.generate_3d()
    ///     minimized = mol.minimize_mmff94(coords)
    fn minimize_mmff94(&self, coords: Vec<[f64; 3]>) -> Vec<Vec<f64>> {
        let c3d = flat_to_coords3d(&coords);
        let out = chematic_3d::minimize_mmff94(&self.inner, c3d);
        out.points.iter().map(|p| vec![p.x, p.y, p.z]).collect()
    }

    /// Minimize 3D coordinates with the DREIDING force field.
    ///
    /// Returns minimized coordinates as ``[[x, y, z], ...]``.
    /// Complements :meth:`minimize_mmff94` and :meth:`minimize_uff`.
    ///
    ///     coords = mol.generate_3d()
    ///     minimized = mol.minimize_dreiding(coords)
    fn minimize_dreiding(&self, coords: Vec<[f64; 3]>) -> Vec<Vec<f64>> {
        let c3d = flat_to_coords3d(&coords);
        let out = chematic_3d::minimize_dreiding(&self.inner, c3d);
        out.points.iter().map(|p| vec![p.x, p.y, p.z]).collect()
    }

    /// Number of stereocenters with unspecified (unknown) configuration.
    ///
    /// Equivalent to RDKit's ``rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)``.
    #[getter]
    fn num_unspecified_stereocenters(&self) -> usize {
        chematic_chem::num_unspecified_stereocenters(&self.inner)
    }

    /// Combined WHIM + GETAWAY 3D descriptor vector.
    ///
    /// Equivalent to calling :meth:`whim` and :meth:`getaway` and concatenating.
    /// Useful for single-call 3D featurisation pipelines (mordred compatible).
    ///
    ///     coords = mol.generate_3d()
    ///     vec = mol.whim_getaway(coords)   # len 41
    fn whim_getaway(&self, coords: Vec<[f64; 3]>) -> Vec<f64> {
        let c3d = flat_to_coords3d(&coords);
        chematic_3d::whim_getaway_combined(&self.inner, &c3d)
    }

    // -----------------------------------------------------------------------
    // Sprint 10: element/bond counts, ring topology, ERG vec, canonical utils
    // -----------------------------------------------------------------------

    /// Number of fluorine atoms.
    #[getter]
    fn num_fluorines(&self) -> usize {
        chematic_chem::num_fluorines(&self.inner)
    }

    /// Number of chlorine atoms.
    #[getter]
    fn num_chlorines(&self) -> usize {
        chematic_chem::num_chlorines(&self.inner)
    }

    /// Number of bromine atoms.
    #[getter]
    fn num_bromines(&self) -> usize {
        chematic_chem::num_bromines(&self.inner)
    }

    /// Number of iodine atoms.
    #[getter]
    fn num_iodines(&self) -> usize {
        chematic_chem::num_iodines(&self.inner)
    }

    /// Number of phosphorus atoms.
    #[getter]
    fn num_phosphorus(&self) -> usize {
        chematic_chem::num_phosphorus(&self.inner)
    }

    /// Number of heteroatoms (non-C, non-H heavy atoms).
    ///
    /// Equivalent to RDKit's ``rdMolDescriptors.CalcNumHeteroatoms(mol)``.
    #[getter]
    fn num_heteroatoms(&self) -> usize {
        chematic_chem::num_heteroatoms(&self.inner)
    }

    /// Number of carbon atoms (heavy, not including implicit H).
    #[getter]
    fn num_carbons(&self) -> usize {
        chematic_chem::num_carbons(&self.inner)
    }

    /// Number of nitrogen atoms.
    #[getter]
    fn num_nitrogens(&self) -> usize {
        chematic_chem::num_nitrogens(&self.inner)
    }

    /// Number of oxygen atoms.
    #[getter]
    fn num_oxygens(&self) -> usize {
        chematic_chem::num_oxygens(&self.inner)
    }

    /// Number of sulfur atoms.
    #[getter]
    fn num_sulfurs(&self) -> usize {
        chematic_chem::num_sulfurs(&self.inner)
    }

    /// Total implicit + explicit hydrogen count.
    #[getter]
    fn num_hydrogens(&self) -> usize {
        chematic_chem::num_hydrogens(&self.inner)
    }

    /// Number of amide bonds (–C(=O)–N–).
    ///
    /// Equivalent to RDKit's ``rdMolDescriptors.CalcNumAmideBonds(mol)``.
    #[getter]
    fn num_amide_bonds(&self) -> usize {
        chematic_chem::num_amide_bonds(&self.inner)
    }

    /// Number of ester bonds (–C(=O)–O–).
    #[getter]
    fn num_ester_bonds(&self) -> usize {
        chematic_chem::num_ester_bonds(&self.inner)
    }

    /// Number of spiro atoms (atoms shared between two rings via a single atom).
    ///
    /// Equivalent to RDKit's ``rdMolDescriptors.CalcNumSpiroAtoms(mol)``.
    #[getter]
    fn num_spiro_atoms(&self) -> usize {
        chematic_chem::num_spiro_atoms(&self.inner)
    }

    /// Number of bridgehead atoms (atoms at the junction of two or more ring systems).
    ///
    /// Equivalent to RDKit's ``rdMolDescriptors.CalcNumBridgeheadAtoms(mol)``.
    #[getter]
    fn num_bridgehead_atoms(&self) -> usize {
        chematic_chem::num_bridgehead_atoms(&self.inner)
    }

    /// Number of aromatic heterocyclic rings.
    #[getter]
    fn num_aromatic_heterocycles(&self) -> usize {
        chematic_chem::num_aromatic_heterocycles(&self.inner)
    }

    /// Number of aliphatic (non-aromatic) heterocyclic rings.
    #[getter]
    fn num_aliphatic_heterocycles(&self) -> usize {
        chematic_chem::num_aliphatic_heterocycles(&self.inner)
    }

    /// Number of saturated heterocyclic rings.
    #[getter]
    fn num_saturated_heterocycles(&self) -> usize {
        chematic_chem::num_saturated_heterocycles(&self.inner)
    }

    /// Number of aliphatic (non-aromatic) carbocyclic rings.
    #[getter]
    fn num_aliphatic_rings(&self) -> usize {
        chematic_chem::num_aliphatic_rings(&self.inner)
    }

    /// ERG (Extended Reduced Graph) continuous feature vector.
    ///
    /// Returns a ``list[float]`` of length 315. Use with
    /// :func:`cosine_erg_vec` or :func:`tanimoto_erg_vec` for similarity.
    ///
    ///     v1 = mol1.erg_vec()
    ///     v2 = mol2.erg_vec()
    ///     sim = chematic.cosine_erg_vec(v1, v2)
    fn erg_vec(&self) -> Vec<f64> {
        chematic_fp::erg_vec(&self.inner).to_vec()
    }

    /// Count-based Morgan (ECFP) fingerprint — hash-count map.
    ///
    /// Returns a dict mapping each substructure environment hash to how many
    /// times it appears at atoms within ``radius`` bonds.
    ///
    /// Unlike :meth:`ecfp4` (bit vector), this preserves multiplicity, making
    /// it more informative for ML regression targets and Tversky similarity.
    /// Equivalent to RDKit's ``GetMorganFingerprint(mol, radius)``.
    ///
    ///     counts = mol.morgan_fp_counts(radius=2)  # {hash: count, ...}
    #[pyo3(signature = (radius = 2))]
    fn morgan_fp_counts<'py>(
        &self,
        py: Python<'py>,
        radius: u32,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        let counts = chematic_fp::morgan_fp_counts(&self.inner, radius);
        let d = pyo3::types::PyDict::new(py);
        for (k, v) in counts {
            d.set_item(k, v)?;
        }
        Ok(d)
    }

    /// Count of each pharmacophore feature type present in this molecule.
    ///
    /// Returns a list of 6 integers in this order:
    /// ``[donor, acceptor, aromatic, hydrophobic, positive, negative]``
    ///
    ///     counts = mol.pharmacophore_feature_counts()
    ///     donors, acceptors = counts[0], counts[1]
    fn pharmacophore_feature_counts(&self) -> Vec<usize> {
        chematic_fp::pharmacophore_feature_counts(&self.inner).to_vec()
    }

    /// MMFF94 partial charges using the atom-type BCI model.
    ///
    /// A more accurate alternative to :meth:`mmff94_charges` (element-pair BCI).
    /// Returns one float per heavy atom.
    ///
    ///     charges = mol.mmff94_charges_typed()
    fn mmff94_charges_typed(&self) -> Vec<f64> {
        chematic_chem::mmff94_charges_typed(&self.inner)
    }

    /// Morgan canonical ranks for each heavy atom.
    ///
    /// Returns one ``int`` per heavy atom (unique ranks based on extended connectivity).
    /// Equivalent to RDKit's ``rdmolfiles.CanonicalRankAtoms(mol)``.
    ///
    ///     ranks = mol.morgan_ranks()
    fn morgan_ranks(&self) -> Vec<u64> {
        chematic_smiles::morgan_ranks(&self.inner)
    }

    /// Canonical atom order — permutation mapping original→canonical index.
    ///
    /// The returned list has length N (heavy atoms); ``order[i]`` is the
    /// canonical position of atom ``i``.
    fn canonical_atom_order(&self) -> Vec<usize> {
        chematic_smiles::canonical_atom_order(&self.inner)
    }

    /// Equivalent atom classes — atoms in the same class are graph-equivalent.
    ///
    /// Returns one class ID per heavy atom (0-indexed). Atoms with the same
    /// class ID are symmetry-equivalent (same canonical rank).
    fn equivalent_atom_classes(&self) -> Vec<usize> {
        chematic_smiles::equivalent_atom_classes(&self.inner)
    }

    // -----------------------------------------------------------------------
    // Sprint 12: num_saturated_rings, zwitterion, remove_salts, invert_stereocenter
    // -----------------------------------------------------------------------

    /// Number of saturated (fully sp³) rings.
    ///
    /// Equivalent to RDKit's ``rdMolDescriptors.CalcNumSaturatedRings(mol)``.
    #[getter]
    fn num_saturated_rings(&self) -> usize {
        chematic_chem::num_saturated_rings(&self.inner)
    }

    /// Return ``True`` if the molecule has simultaneous positive and negative formal charges.
    ///
    ///     assert chematic.from_smiles("[NH3+]CC([O-])=O").has_zwitterion()
    fn has_zwitterion(&self) -> bool {
        chematic_chem::has_zwitterion(&self.inner)
    }

    /// Normalize a zwitterion to neutral form via proton transfer.
    ///
    /// Equivalent to RDKit ``MolStandardize.Standardizer`` zwitterion step.
    ///
    ///     neutral = mol.normalize_zwitterion()
    fn normalize_zwitterion(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::normalize_zwitterion(&self.inner)),
        }
    }

    /// Remove salt fragments using the built-in salt catalog.
    ///
    /// Unlike :meth:`largest_fragment` (which keeps the heaviest fragment),
    /// this uses a catalog of known counterions/salts to remove them specifically.
    ///
    ///     desalted = chematic.from_smiles("CC(=O)[O-].[Na+]").remove_salts()
    fn remove_salts(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::standardize::remove_salts(&self.inner)),
        }
    }

    /// Invert the stereochemistry at atom ``atom_idx`` (flip wedge/dash bonds).
    ///
    /// Returns a new molecule with Up↔Down bonds at the specified center swapped,
    /// generating the enantiomer at that center. Returns an unchanged copy if the
    /// atom has no explicit stereo annotation.
    ///
    ///     inverted = mol.invert_stereocenter(0)
    fn invert_stereocenter(&self, atom_idx: usize) -> Mol {
        let idx = chematic_core::AtomIdx(atom_idx as u32);
        Mol {
            inner: Arc::new(chematic_chem::invert_stereocenter(&self.inner, idx)),
        }
    }

    // -----------------------------------------------------------------------
    // Sprint 11: topo descriptors, ring perception, stereo validation, pharmacophore
    // -----------------------------------------------------------------------

    /// Hall-Kier κ₁ shape index — encodes molecular size relative to a linear chain.
    #[getter]
    fn kappa1(&self) -> f64 {
        chematic_chem::kappa1(&self.inner)
    }

    /// Hall-Kier κ₂ shape index — encodes branching degree.
    #[getter]
    fn kappa2(&self) -> f64 {
        chematic_chem::kappa2(&self.inner)
    }

    /// Hall-Kier κ₃ shape index — encodes centrality of branching.
    #[getter]
    fn kappa3(&self) -> f64 {
        chematic_chem::kappa3(&self.inner)
    }

    /// Wiener index — sum of all topological distances between heavy atom pairs.
    #[getter]
    fn wiener_index(&self) -> f64 {
        chematic_chem::wiener_index(&self.inner)
    }

    /// Padmakar-Ivan (PI) topological index (Khadikar et al. 2001).
    ///
    /// For each bond e = (u, v): PI += n_u(e) + n_v(e), where n_u(e) is the
    /// number of heavy atoms strictly closer to u than to v.
    /// Reference: ethane = 2, propane = 6, butane = 12, benzene = 36.
    #[getter]
    fn padmakar_ivan_index(&self) -> u64 {
        chematic_chem::padmakar_ivan_index(&self.inner)
    }

    /// Schultz molecular topological index (MTI).
    ///
    /// MTI = Σ_{i<j} (δᵢ + δⱼ) × dᵢⱼ, where δᵢ is heavy-atom degree.
    /// Ref: Schultz, *J. Chem. Inf. Comput. Sci.* **1989**, 29, 227.
    #[getter]
    fn schultz_mti(&self) -> u64 {
        chematic_chem::schultz_mti(&self.inner)
    }

    /// Gutman molecular topological index (MTI*).
    ///
    /// MTI* = Σ_{i<j} δᵢ × δⱼ × dᵢⱼ, where δᵢ is heavy-atom degree.
    /// Ref: Gutman, *J. Serb. Chem. Soc.* **1994**, 59, 619.
    #[getter]
    fn gutman_mti(&self) -> u64 {
        chematic_chem::gutman_mti(&self.inner)
    }

    /// VABC van der Waals volume approximation (Å³).
    ///
    /// Estimated from Bondi radii with spherical-cap overlap corrections for bonds.
    /// Does not require 3D coordinates.
    /// Ref: Zhao et al., *J. Org. Chem.* **2003**, 68, 7368.
    #[getter]
    fn vabc(&self) -> f64 {
        chematic_chem::vabc(&self.inner)
    }

    /// Gravitational topological index.
    ///
    /// G = Σ_{i<j} mᵢ × mⱼ / dᵢⱼ², where mᵢ is average atomic mass.
    #[getter]
    fn gravitational_index(&self) -> f64 {
        chematic_chem::gravitational_index(&self.inner)
    }

    /// Bertz complexity index — information-theoretic graph complexity.
    #[getter]
    fn bertz_ct(&self) -> f64 {
        chematic_chem::bertz_ct(&self.inner)
    }

    /// Zero-order path connectivity index χ⁰ (Kier & Hall).
    #[getter]
    fn chi0(&self) -> f64 {
        chematic_chem::chi0(&self.inner)
    }

    /// First-order path connectivity index χ¹.
    #[getter]
    fn chi1(&self) -> f64 {
        chematic_chem::chi1(&self.inner)
    }

    /// Second-order path connectivity index χ².
    #[getter]
    fn chi2(&self) -> f64 {
        chematic_chem::chi2(&self.inner)
    }

    /// Third-order path connectivity index χ³.
    #[getter]
    fn chi3(&self) -> f64 {
        chematic_chem::chi3(&self.inner)
    }

    /// Fourth-order path connectivity index χ⁴.
    #[getter]
    fn chi4(&self) -> f64 {
        chematic_chem::chi4(&self.inner)
    }

    /// Zero-order valence connectivity index χ⁰ᵥ.
    #[getter]
    fn chi0v(&self) -> f64 {
        chematic_chem::chi0v(&self.inner)
    }

    /// First-order valence connectivity index χ¹ᵥ.
    #[getter]
    fn chi1v(&self) -> f64 {
        chematic_chem::chi1v(&self.inner)
    }

    /// Second-order valence connectivity index χ²ᵥ.
    #[getter]
    fn chi2v(&self) -> f64 {
        chematic_chem::chi2v(&self.inner)
    }

    /// Third-order valence connectivity index χ³ᵥ.
    #[getter]
    fn chi3v(&self) -> f64 {
        chematic_chem::chi3v(&self.inner)
    }

    /// Fourth-order valence connectivity index χ⁴ᵥ.
    #[getter]
    fn chi4v(&self) -> f64 {
        chematic_chem::chi4v(&self.inner)
    }

    /// SSSR ring membership per atom.
    ///
    /// Returns a list of N lists (one per heavy atom). Each inner list contains
    /// the 0-based SSSR ring indices to which that atom belongs.
    /// Acyclic atoms have an empty inner list.
    ///
    ///     membership = mol.ring_membership()
    ///     ring_idxs = membership[atom_i]  # list of ring indices
    fn ring_membership(&self) -> Vec<Vec<usize>> {
        chematic_perception::ring_membership(&self.inner)
    }

    /// Ring sizes of all SSSR rings that contain atom ``atom_idx``.
    ///
    /// Returns an empty list for acyclic atoms.
    ///
    ///     sizes = mol.ring_sizes_for_atom(0)  # e.g. [6] for benzene C
    fn ring_sizes_for_atom(&self, atom_idx: usize) -> Vec<usize> {
        chematic_perception::ring_sizes_for_atom(&self.inner, atom_idx)
    }

    /// Return ``True`` if the molecule contains a fused ring system.
    ///
    /// Two rings are fused when they share at least one bond (≥ 2 adjacent atoms).
    /// Spiro systems (sharing exactly one atom) return ``False``.
    ///
    /// Equivalent to checking whether indene, naphthalene, etc. are present.
    fn is_fused_ring_system(&self) -> bool {
        chematic_perception::is_fused_ring_system(&self.inner)
    }

    /// Classify the ring systems (families) in this molecule by topology.
    ///
    /// Returns a list of dicts, one per connected ring system, each with:
    ///
    /// - ``kind``: ``"simple"`` | ``"fused"`` | ``"spiro"`` | ``"bridged"``
    /// - ``atom_indices``: list of heavy-atom indices belonging to this ring system
    /// - ``ring_count``: number of SSSR rings in this family
    ///
    /// Two SSSR rings belong to the same family if they share at least one atom.
    ///
    ///     for fam in mol.ring_families():
    ///         print(fam['kind'], fam['ring_count'])
    fn ring_families<'py>(&self, py: Python<'py>) -> Vec<Bound<'py, PyDict>> {
        use chematic_perception::{RingSystemKind, find_ring_families, find_sssr};
        let sssr = find_sssr(&self.inner);
        find_ring_families(&self.inner, &sssr)
            .into_iter()
            .map(|fam| {
                let d = PyDict::new(py);
                let kind_str = match fam.kind {
                    RingSystemKind::Simple => "simple",
                    RingSystemKind::Fused => "fused",
                    RingSystemKind::Spiro => "spiro",
                    RingSystemKind::Bridged => "bridged",
                };
                d.set_item("kind", kind_str).unwrap();
                d.set_item(
                    "atom_indices",
                    fam.atoms.iter().map(|a| a.0 as usize).collect::<Vec<_>>(),
                )
                .unwrap();
                d.set_item("ring_count", fam.ring_indices.len()).unwrap();
                d
            })
            .collect()
    }

    /// Validate stereochemistry and return any detected errors.
    ///
    /// Each error is a dict with keys ``atom_idx`` (int) and ``kind`` (str).
    /// Possible ``kind`` values:
    ///   - ``"ImpossibleCenter"`` — fewer than 4 distinct neighbours
    ///   - ``"ConflictingWedges"`` — both Up and Down bonds from same center
    ///   - ``"RedundantStereo"`` — all neighbours have identical rank
    ///
    /// An empty list means the stereo annotations are chemically consistent.
    ///
    ///     errors = mol.validate_stereo()
    fn validate_stereo<'py>(&self, py: Python<'py>) -> Vec<Bound<'py, PyDict>> {
        chematic_perception::validate_stereo(&self.inner)
            .into_iter()
            .map(|e| {
                let d = PyDict::new(py);
                d.set_item("atom_idx", e.atom_idx).unwrap();
                d.set_item("kind", format!("{:?}", e.kind)).unwrap();
                d
            })
            .collect()
    }

    /// Summarise stereocenters — how many are specified vs unspecified.
    ///
    /// Returns a dict with keys ``specified``, ``unspecified``, ``total_centers``.
    ///
    ///     sc = mol.stereo_completeness()
    ///     print(sc["specified"], "/", sc["total_centers"])
    fn stereo_completeness<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let sc = chematic_perception::stereo_completeness(&self.inner);
        let d = PyDict::new(py);
        d.set_item("specified", sc.specified)?;
        d.set_item("unspecified", sc.unspecified)?;
        d.set_item("total_centers", sc.total_centers)?;
        Ok(d)
    }

    /// Perceive stereochemistry (R/S and E/Z) from 3D coordinates.
    ///
    /// Returns a list of dicts, one per assigned stereocentre:
    ///
    /// - ``atom_idx``: heavy-atom index of the chiral centre or E/Z double-bond atom
    /// - ``code``: ``"R"``, ``"S"``, ``"E"``, or ``"Z"``
    ///
    /// An empty list means no stereocentres could be assigned.
    ///
    /// .. note::
    ///     Only assigns R/S for atoms with **four heavy-atom** neighbours (no
    ///     implicit H). Chiral centres with an implicit H (e.g. amino acids)
    ///     are not currently assigned by this function.
    ///
    /// Equivalent to RDKit ``Chem.AssignStereochemistryFrom3D(mol)``.
    ///
    ///     coords = mol.generate_3d()
    ///     for a in mol.stereo_from_coords(coords):
    ///         print(a['atom_idx'], a['code'])  # e.g. 0 R
    fn stereo_from_coords<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
    ) -> Vec<Bound<'py, PyDict>> {
        use chematic_core::CipCode;
        let c3d = flat_to_coords3d(&coords);
        let assignment = chematic_3d::assign_stereo_from_3d(&self.inner, &c3d);
        assignment
            .assignments
            .iter()
            .map(|(idx, code)| {
                let d = PyDict::new(py);
                d.set_item("atom_idx", idx.0 as usize).unwrap();
                let code_str = match code {
                    CipCode::R => "R",
                    CipCode::S => "S",
                    CipCode::E => "E",
                    CipCode::Z => "Z",
                };
                d.set_item("code", code_str).unwrap();
                d
            })
            .collect()
    }

    /// Perceive stereochemistry (R/S and E/Z) from 2D layout coordinates.
    ///
    /// Returns a list of dicts, one per assigned stereocentre:
    ///
    /// - ``atom_idx``: heavy-atom index of the chiral centre or E/Z bond atom
    /// - ``code``: ``"R"``, ``"S"``, ``"E"``, or ``"Z"``
    ///
    /// Coordinates are typically obtained from :func:`from_mol_block_with_coords`.
    /// For 3D-coordinate-based assignment use :meth:`stereo_from_coords`.
    ///
    ///     mol, name, coords_2d = chematic.from_mol_block_with_coords(block)
    ///     for a in mol.stereo_from_2d_coords(coords_2d):
    ///         print(a['atom_idx'], a['code'])
    fn stereo_from_2d_coords<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 2]>,
    ) -> Vec<Bound<'py, PyDict>> {
        use chematic_core::CipCode;
        let coords_2d: Vec<(f64, f64)> = coords.iter().map(|c| (c[0], c[1])).collect();
        let assignment = chematic_perception::assign_stereo_from_2d(&self.inner, &coords_2d);
        assignment
            .assignments
            .iter()
            .map(|(idx, code)| {
                let d = PyDict::new(py);
                d.set_item("atom_idx", idx.0 as usize).unwrap();
                let code_str = match code {
                    CipCode::R => "R",
                    CipCode::S => "S",
                    CipCode::E => "E",
                    CipCode::Z => "Z",
                };
                d.set_item("code", code_str).unwrap();
                d
            })
            .collect()
    }

    /// Detect pharmacophore features in the molecule.
    ///
    /// Returns a list of dicts, each with:
    ///   - ``"type"`` (str): ``"Donor"``, ``"Acceptor"``, ``"Aromatic"``,
    ///     ``"Hydrophobic"``, ``"Positive"``, or ``"Negative"``
    ///   - ``"atom_idx"`` (int): primary atom index
    ///   - ``"neighbor_indices"`` (list[int]): secondary atoms (e.g., all ring atoms
    ///     for aromatic features)
    ///
    ///     feats = mol.pharmacophore_features()
    ///     donors = [f for f in feats if f["type"] == "Donor"]
    fn pharmacophore_features<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        chematic_perception::detect_features(&self.inner)
            .into_iter()
            .map(|f| {
                let d = PyDict::new(py);
                d.set_item("type", format!("{:?}", f.ftype))?;
                d.set_item("atom_idx", f.atom.0 as usize)?;
                let nb: Vec<usize> = f.neighbors.iter().map(|a| a.0 as usize).collect();
                d.set_item("neighbor_indices", nb)?;
                Ok(d)
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
// Helpers
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

fn flat_to_coords3d(coords: &[[f64; 3]]) -> chematic_3d::Coords3D {
    chematic_3d::Coords3D {
        points: coords
            .iter()
            .map(|c| chematic_3d::Point3::new(c[0], c[1], c[2]))
            .collect(),
    }
}

// ---------------------------------------------------------------------------
// Helper: convert MoleculeReport to Python dict
// ---------------------------------------------------------------------------

fn mol_report_to_dict<'py>(
    py: Python<'py>,
    report: &chematic_chem::MoleculeReport,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("canonical_smiles", &report.canonical_smiles)?;
    d.set_item("formula", &report.formula)?;
    match &report.murcko_scaffold_smiles {
        Some(s) => d.set_item("murcko_scaffold", s)?,
        None => d.set_item("murcko_scaffold", py.None())?,
    }
    let desc = PyDict::new(py);
    desc.set_item("mw", report.descriptors.molecular_weight)?;
    desc.set_item("exact_mass", report.descriptors.exact_mass)?;
    desc.set_item("tpsa", report.descriptors.tpsa)?;
    desc.set_item("logp", report.descriptors.logp)?;
    desc.set_item("molar_refractivity", report.descriptors.molar_refractivity)?;
    desc.set_item("hbd", report.descriptors.hbd)?;
    desc.set_item("hba", report.descriptors.hba)?;
    desc.set_item("rotatable_bonds", report.descriptors.rotatable_bonds)?;
    desc.set_item("heavy_atoms", report.descriptors.heavy_atom_count)?;
    desc.set_item("ring_count", report.descriptors.ring_count)?;
    desc.set_item("num_heteroatoms", report.descriptors.num_heteroatoms)?;
    desc.set_item("num_stereocenters", report.descriptors.num_stereocenters)?;
    desc.set_item("fsp3", report.descriptors.fsp3)?;
    desc.set_item("qed", report.descriptors.qed)?;
    desc.set_item("sa_score", report.descriptors.sa_score)?;
    desc.set_item("formal_charge", report.descriptors.formal_charge_sum)?;
    desc.set_item("labute_asa", report.descriptors.labute_asa)?;
    desc.set_item("bertz_ct", report.descriptors.bertz_ct)?;
    desc.set_item("wiener_index", report.descriptors.wiener_index)?;
    d.set_item("descriptors", desc)?;
    let filters = PyDict::new(py);
    filters.set_item("lipinski_passes", report.filters.lipinski_passes)?;
    filters.set_item("veber_passes", report.filters.veber_passes)?;
    filters.set_item("egan_passes", report.filters.egan_passes)?;
    filters.set_item("ghose_passes", report.filters.ghose_passes)?;
    filters.set_item("reos_passes", report.filters.reos_passes)?;
    filters.set_item("pains_passes", report.filters.pains_passes)?;
    let alerts: Vec<&str> = report
        .filters
        .pains_alerts
        .iter()
        .map(|s| s.as_str())
        .collect();
    filters.set_item("pains_alerts", alerts)?;
    d.set_item("filters", filters)?;
    let fgs: Vec<Bound<'py, PyDict>> = report
        .functional_groups
        .iter()
        .map(|fg| {
            let fd = PyDict::new(py);
            fd.set_item("name", &fg.name)?;
            fd.set_item("atom_indices", &fg.atom_indices)?;
            Ok(fd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("functional_groups", fgs)?;
    let ngs: Vec<Bound<'py, PyDict>> = report
        .named_groups
        .iter()
        .map(|ng| {
            let nd = PyDict::new(py);
            nd.set_item("name", &ng.name)?;
            nd.set_item("atom_indices", &ng.atom_indices)?;
            Ok(nd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("named_groups", ngs)?;
    Ok(d)
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
        .map(|mol| Mol {
            inner: Arc::new(mol),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a CXSMILES string and return the molecule with CX metadata.
///
/// Returns a 2-tuple ``(mol, cx)`` where ``cx`` is a dict with:
///
/// - ``atom_labels``: list of atom label strings (or ``None`` per atom)
/// - ``atom_props``: list of ``{"atom_idx", "key", "value"}`` dicts
/// - ``atom_radicals``: list of radical class integers (or ``None`` per atom)
///
/// Raises ``ValueError`` on parse failure. CXSMILES without a CX extension
/// block behaves like :func:`from_smiles` (all CX fields are empty).
///
///     mol, cx = chematic.from_cxsmiles("CC |$R1;R2$|")
///     print(cx['atom_labels'])   # ['R1', 'R2']
#[pyfunction]
fn from_cxsmiles<'py>(py: Python<'py>, s: &str) -> PyResult<(Mol, Bound<'py, PyDict>)> {
    let cx =
        chematic_smiles::parse_cxsmiles(s).map_err(|e| PyValueError::new_err(e.to_string()))?;

    let d = PyDict::new(py);
    let labels: Vec<Option<&str>> = cx.atom_labels.iter().map(|l| l.as_deref()).collect();
    d.set_item("atom_labels", labels)?;
    let props: Vec<Bound<'py, PyDict>> = cx
        .atom_props
        .iter()
        .map(|p| {
            let pd = PyDict::new(py);
            pd.set_item("atom_idx", p.atom.0 as usize).unwrap();
            pd.set_item("key", &p.key).unwrap();
            pd.set_item("value", &p.value).unwrap();
            pd
        })
        .collect();
    d.set_item("atom_props", props)?;
    let radicals: Vec<Option<u8>> = cx.atom_radicals.clone();
    d.set_item("atom_radicals", radicals)?;

    Ok((
        Mol {
            inner: Arc::new(cx.mol),
        },
        d,
    ))
}

/// Parse a MOL/SDF block and return a Mol.
///
/// Parse a condensed molecular formula (e.g., ``"CH3OH"``, ``"C6H12O6"``) into a :class:`Mol`.
///
/// Returns ``None`` if the formula is unknown or ambiguous.
/// Unlike SMILES parsing, condensed formulas may not encode connectivity uniquely;
/// this function uses a built-in formula→structure dictionary.
///
/// Equivalent to chempy's condensed formula support.
///
///     mol = chematic.from_condensed("CH3OH")  # methanol
///     if mol:
///         print(mol.smiles)  # CO
///
///     mol = chematic.from_condensed("C6H12O6")  # glucose
#[pyfunction]
fn from_condensed(formula: &str) -> Option<Mol> {
    chematic_chem::parse_condensed(formula).ok().map(|mol| Mol {
        inner: Arc::new(mol),
    })
}

/// Raises ``ValueError`` on parse failure.
#[pyfunction]
fn from_mol_block(block: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol(block)
        .map(|(mol, _meta)| Mol {
            inner: Arc::new(mol),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a MDL MOL V2000 block and return the molecule with its 2D layout coordinates.
///
/// Returns a 3-tuple ``(mol, name, coords_2d)`` where:
///
/// - ``mol``: :class:`Mol` object
/// - ``name``: molecule name from the MOL header (may be empty)
/// - ``coords_2d``: list of ``[x, y]`` pairs (one per heavy atom, Å)
///
/// Raises ``ValueError`` on parse failure.
///
/// Use :func:`from_mol_block` if you only need the molecule graph.
/// Use this function when you want to preserve the 2D layout for display or
/// round-trip back to MOL format via :meth:`Mol.to_mol_block_2d`.
///
///     mol, name, coords_2d = chematic.from_mol_block_with_coords(block)
///     new_block = mol.to_mol_block_2d(coords_2d, name=name)
#[pyfunction]
fn from_mol_block_with_coords(block: &str) -> PyResult<(Mol, String, Vec<Vec<f64>>)> {
    chematic_mol::parse_mol_with_coords(block)
        .map(|(mol, meta, coords)| {
            let py_coords: Vec<Vec<f64>> = coords.iter().map(|(x, y)| vec![*x, *y]).collect();
            (
                Mol {
                    inner: Arc::new(mol),
                },
                meta.name,
                py_coords,
            )
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a multi-record SDF string and return all molecules with their 2D layout coordinates.
///
/// Returns a list of 3-tuples ``(mol, name, coords_2d)`` — one per SDF record.
/// Invalid records are silently skipped (same behaviour as :func:`iter_sdf`).
///
/// This is the batch equivalent of :func:`from_mol_block_with_coords`.
///
///     with open("library.sdf") as f:
///         records = chematic.parse_sdf_with_coords(f.read())
///     for mol, name, coords_2d in records:
///         new_block = mol.to_mol_block_2d(coords_2d, name=name)
#[pyfunction]
fn parse_sdf_with_coords(text: &str) -> Vec<(Mol, String, Vec<Vec<f64>>)> {
    // Split SDF by $$$$ delimiter and parse each block with parse_mol_with_coords.
    // This avoids the Rust parse_sdf_with_coords leading-blank-line stripping issue.
    let mut results = Vec::new();
    let mut remaining = text;
    loop {
        let (block, rest) = match remaining.find("$$$$") {
            Some(pos) => {
                let after = &remaining[pos + 4..];
                let after = after
                    .strip_prefix("\r\n")
                    .or_else(|| after.strip_prefix('\n'))
                    .unwrap_or(after);
                (&remaining[..pos], after)
            }
            None => (remaining, ""),
        };
        if !block.trim().is_empty() {
            if let Ok((mol, meta, coords)) = chematic_mol::parse_mol_with_coords(block) {
                let py_coords: Vec<Vec<f64>> = coords.iter().map(|(x, y)| vec![*x, *y]).collect();
                results.push((
                    Mol {
                        inner: Arc::new(mol),
                    },
                    meta.name,
                    py_coords,
                ));
            }
        }
        if remaining.find("$$$$").is_none() || rest.is_empty() {
            break;
        }
        remaining = rest;
    }
    results
}

/// Parse a Chemical Markup Language (CML) string into a ``Mol`` object.
///
/// Raises ``ValueError`` on parse failure.
///
///     with open("molecule.cml") as f:
///         mol = chematic.from_cml(f.read())
#[pyfunction]
fn from_cml(cml_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_cml(cml_str)
        .map(|(mol, _coords)| Mol {
            inner: Arc::new(mol),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a ChemicalJSON (.cjson) string.
///
/// Returns ``(mol, coords)`` where ``coords`` is a list of ``[x, y, z]``
/// coordinate triples (Å), one per heavy atom.  ``coords`` is empty when
/// the file has no ``atoms.coords.3d`` field.
///
/// ChemicalJSON is the native format of Avogadro 2 and the MolSSI
/// Open Chemistry toolkit.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol, coords = chematic.from_cjson(open("mol.cjson").read())
///     print(mol.smiles)
///     # Round-trip:
///     open("out.cjson", "w").write(mol.to_cjson(coords))
#[pyfunction]
fn from_cjson(cjson_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let (mol, coords) = chematic_mol::parse_cjson(cjson_str)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let py_coords = coords.iter().map(|&(x, y, z)| vec![x, y, z]).collect();
    Ok((Mol { inner: Arc::new(mol) }, py_coords))
}

/// Parse a ChemDraw XML (CDXML) string into a ``Mol`` object.
///
/// Raises ``ValueError`` on parse failure.
///
///     with open("molecule.cdxml") as f:
///         mol = chematic.from_cdxml(f.read())
#[pyfunction]
fn from_cdxml(cdxml_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_cdxml(cdxml_str)
        .map(|(mol, _coords)| Mol {
            inner: Arc::new(mol),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse an MDL MOL V3000 (``V3000``) block into a ``Mol`` object.
///
/// Raises ``ValueError`` on parse failure.
///
///     with open("ligand_v3000.mol") as f:
///         mol = chematic.from_mol_v3000(f.read())
#[pyfunction]
fn from_mol_v3000(block: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol_v3000(block)
        .map(|(mol, _meta)| Mol {
            inner: Arc::new(mol),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a MDL MOL V3000 block and return the molecule with its 2D layout coordinates.
///
/// Returns a 3-tuple ``(mol, name, coords_2d)`` identical to
/// :func:`from_mol_block_with_coords` but for V3000 input.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol, name, coords_2d = chematic.from_mol_v3000_with_coords(block)
///     new_block = mol.to_mol_v3000(coords_2d, name=name)
#[pyfunction]
fn from_mol_v3000_with_coords(block: &str) -> PyResult<(Mol, String, Vec<Vec<f64>>)> {
    chematic_mol::parse_mol_v3000_with_coords(block)
        .map(|(mol, meta, coords)| {
            let py_coords: Vec<Vec<f64>> = coords.iter().map(|(x, y)| vec![*x, *y]).collect();
            (
                Mol {
                    inner: Arc::new(mol),
                },
                meta.name,
                py_coords,
            )
        })
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
        .map(|(mol, _coords)| Mol {
            inner: Arc::new(mol),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse an AutoDock PDBQT string and return a :class:`Mol`.
///
/// Only the molecular graph (elements and bonds) is extracted; 3D coordinates
/// and partial charges are discarded.  To retain them, use the lower-level
/// :func:`chematic_mol.parse_pdbqt` Rust API directly.
///
/// Raises:
///     ValueError: on parse failure.
///
/// Example::
///
///     with open("ligand.pdbqt") as f:
///         mol = chematic.from_pdbqt(f.read())
///     print(mol.mw)
#[pyfunction]
fn from_pdbqt(pdbqt_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_pdbqt(pdbqt_str)
        .map(|(mol, _coords, _charges)| Mol {
            inner: Arc::new(mol),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Gaussian input file (`.gjf` / `.com`) and return a :class:`Mol`.
///
/// Raises:
///     ValueError: on parse failure.
///
/// Example::
///
///     mol = chematic.from_gjf(open("mol.gjf").read())
#[pyfunction]
fn from_gjf(gjf_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_gjf(gjf_str)
        .map(|(mol, _coords, _charge, _mult)| Mol {
            inner: Arc::new(mol),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Gaussian output file (`.log` / `.out`) and return a dict with
/// ``mol``, ``coords`` and ``scf_energy`` fields.
///
/// Returns:
///     dict: ``{"mol": Mol, "coords": list[list[float]], "scf_energy": float | None}``
///
/// Raises:
///     ValueError: when no `Standard orientation:` block is found.
#[pyfunction]
fn parse_gaussian_log<'py>(
    py: Python<'py>,
    log_str: &str,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let result = chematic_mol::parse_gaussian_log(log_str)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let mol = Mol {
        inner: Arc::new(result.mol),
    };
    let coords: Vec<Vec<f64>> = result
        .coords
        .iter()
        .map(|&(x, y, z)| vec![x, y, z])
        .collect();
    let d = pyo3::types::PyDict::new(py);
    d.set_item("mol", mol)?;
    d.set_item("coords", coords)?;
    d.set_item("scf_energy", result.scf_energy)?;
    Ok(d)
}

/// Generate a Gaussian input file (`.gjf`) string from a molecule.
///
/// Args:
///     mol: The molecule to write.
///     coords: Atomic coordinates as ``[[x, y, z], ...]`` in Ångströms.
///     charge: Formal charge (default 0).
///     multiplicity: Spin multiplicity (default 1).
///     method: Route section keywords (default ``"B3LYP/6-31G* opt"``).
///     title: Job title comment (default ``"chematic"``).
///
/// Returns:
///     str: GJF file contents.
#[pyfunction]
#[pyo3(signature = (mol, coords, charge=0, multiplicity=1, method="B3LYP/6-31G* opt", title="chematic"))]
fn write_gjf(
    mol: &Mol,
    coords: Vec<[f64; 3]>,
    charge: i32,
    multiplicity: u32,
    method: &str,
    title: &str,
) -> String {
    let c: Vec<(f64, f64, f64)> = coords.into_iter().map(|[x, y, z]| (x, y, z)).collect();
    chematic_mol::write_gjf(&mol.inner, &c, charge, multiplicity, method, title)
}

/// Parse a CIF (Crystallographic Information File) string and return a dict.
///
/// Returns:
///     dict: ``{"mol": Mol, "coords": list[list[float]], "cell": dict | None}``
///     where ``cell`` has keys ``a, b, c, alpha, beta, gamma``.
///
/// Raises:
///     ValueError: on parse failure.
///
/// Example::
///
///     result = chematic.parse_cif(open("structure.cif").read())
///     mol = result["mol"]
///     coords = result["coords"]
#[pyfunction]
fn parse_cif<'py>(py: Python<'py>, cif_str: &str) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let result =
        chematic_mol::parse_cif(cif_str).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let mol = Mol {
        inner: Arc::new(result.mol),
    };
    let coords: Vec<Vec<f64>> = result
        .coords
        .iter()
        .map(|&(x, y, z)| vec![x, y, z])
        .collect();
    let d = pyo3::types::PyDict::new(py);
    d.set_item("mol", mol)?;
    d.set_item("coords", coords)?;
    if let Some(cell) = result.cell {
        let cd = pyo3::types::PyDict::new(py);
        cd.set_item("a", cell.a)?;
        cd.set_item("b", cell.b)?;
        cd.set_item("c", cell.c)?;
        cd.set_item("alpha", cell.alpha)?;
        cd.set_item("beta", cell.beta)?;
        cd.set_item("gamma", cell.gamma)?;
        d.set_item("cell", cd)?;
    } else {
        d.set_item("cell", py.None())?;
    }
    Ok(d)
}

/// Return True if the SMILES can be parsed without error.
#[pyfunction]
fn is_valid_smiles(smiles: &str) -> bool {
    chematic_smiles::parse(smiles).is_ok()
}

/// Return ``True`` if ``smarts`` is a valid SMARTS pattern, ``False`` otherwise.
///
/// Mirrors :func:`is_valid_smiles` for SMARTS pattern validation.
/// Useful for validating user-supplied SMARTS before calling
/// :func:`smarts_match` or :func:`smarts_find`.
///
///     chematic.is_valid_smarts("c1ccccc1")  # True
///     chematic.is_valid_smarts("[invalid")  # False
///     chematic.is_valid_smarts("[#6]-[#7]") # True
#[pyfunction]
fn is_valid_smarts(smarts: &str) -> bool {
    chematic_smarts::parse_smarts(smarts).is_ok()
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

/// Test whether a SMARTS pattern matches a molecule.
///
///     if chematic.smarts_match("[OH]", mol):
///         print("has hydroxyl")
#[pyfunction]
fn smarts_match(smarts: &str, mol: &Mol) -> PyResult<bool> {
    let query =
        chematic_smarts::parse_smarts(smarts).map_err(|e| PyValueError::new_err(e.to_string()))?;
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

/// Parse an InChI string and return a Mol.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol = chematic.from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
#[pyfunction]
fn from_inchi(inchi: &str) -> PyResult<Mol> {
    chematic_inchi::parse_inchi(inchi)
        .map(|mol| Mol {
            inner: Arc::new(mol),
        })
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

/// Parse a PDB string and return ``(Mol, coords)`` where coords is a list of ``[x,y,z]``.
///
/// Bond information is inferred from inter-atom distances; only ATOM/HETATM records
/// are used.  Returns ``ValueError`` when no atoms are found.
///
///     mol, coords = chematic.from_pdb(open("ligand.pdb").read())
#[pyfunction]
fn from_pdb(pdb_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let atoms = chematic_3d::parse_pdb_atoms(pdb_str);
    if atoms.is_empty() {
        return Err(PyValueError::new_err(
            "no ATOM/HETATM records found in PDB input",
        ));
    }
    let (mol, c3d) = chematic_3d::pdb_to_molecule(&atoms);
    let coords = c3d.points.iter().map(|p| vec![p.x, p.y, p.z]).collect();
    Ok((
        Mol {
            inner: Arc::new(mol),
        },
        coords,
    ))
}

/// Parse an XYZ string and return ``(Mol, coords)`` where coords is a list of ``[x,y,z]``.
///
/// Bond information is inferred from inter-atom distances.
///
///     mol, coords = chematic.from_xyz(open("molecule.xyz").read())
#[pyfunction]
fn from_xyz(xyz_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let (mol, c3d) =
        chematic_3d::parse_xyz(xyz_str).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let coords = c3d.points.iter().map(|p| vec![p.x, p.y, p.z]).collect();
    Ok((
        Mol {
            inner: Arc::new(mol),
        },
        coords,
    ))
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

/// Generate a complete single-molecule analysis report.
///
/// Returns a dict with keys:
///   ``canonical_smiles``, ``formula``, ``murcko_scaffold``,
///   ``descriptors`` (dict), ``filters`` (dict),
///   ``functional_groups`` (list of dicts), ``named_groups`` (list of dicts).
///
///     report = chematic.molecule_report("CC(=O)Oc1ccccc1C(=O)O")
///     print(report["descriptors"]["mw"])
///     print(report["filters"]["lipinski_passes"])
#[pyfunction]
fn molecule_report<'py>(smiles: &str, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let report =
        chematic_chem::molecule_report(smiles).map_err(|e| PyValueError::new_err(e.to_string()))?;
    mol_report_to_dict(py, &report)
}

/// Screen a list of SMILES and return a batch report with diversity analysis.
///
/// Returns a dict with keys:
///   ``records`` (list of per-molecule dicts with ``input_index``, ``smiles``, ``report``, ``error``),
///   ``maxmin_picks`` (list of indices — most diverse subset),
///   ``butina_clusters`` (list of cluster index lists).
///
///     result = chematic.screen_smiles(smiles_list)
///     for rec in result["records"]:
///         if rec["error"] is None:
///             print(rec["report"]["descriptors"]["mw"])
#[pyfunction]
fn screen_smiles<'py>(smiles: Vec<String>, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let refs: Vec<&str> = smiles.iter().map(|s| s.as_str()).collect();
    let report = chematic_chem::screen_smiles(&refs);
    let d = PyDict::new(py);
    let records: Vec<Bound<'py, PyDict>> = report
        .records
        .iter()
        .map(|rec| {
            let r = PyDict::new(py);
            r.set_item("input_index", rec.input_index)?;
            r.set_item("smiles", &rec.input_smiles)?;
            match &rec.report {
                Some(mr) => {
                    r.set_item("report", mol_report_to_dict(py, mr)?)?;
                    r.set_item("error", py.None())?;
                }
                None => {
                    r.set_item("report", py.None())?;
                    r.set_item("error", rec.error.as_deref().unwrap_or("unknown error"))?;
                }
            }
            Ok(r)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("records", records)?;
    d.set_item("maxmin_picks", report.maxmin_picks)?;
    d.set_item("butina_clusters", report.butina_clusters)?;
    Ok(d)
}

/// Compare two or more SMILES strings and return pairwise similarity + descriptor deltas.
///
/// Returns a dict with keys:
///   ``reports`` (list of molecule report dicts),
///   ``pairwise`` (list of pairwise similarity dicts),
///   ``descriptor_deltas`` (list of delta dicts),
///   ``mcs_smiles`` (str or None).
///
///     result = chematic.compare_molecules(["c1ccccc1", "Cc1ccccc1"])
///     print(result["pairwise"][0]["ecfp4_tanimoto"])
#[pyfunction]
fn compare_molecules<'py>(smiles: Vec<String>, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let refs: Vec<&str> = smiles.iter().map(|s| s.as_str()).collect();
    let cmp = chematic_chem::compare_molecules(&refs)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let d = PyDict::new(py);
    let reports: Vec<Bound<'py, PyDict>> = cmp
        .reports
        .iter()
        .map(|r| mol_report_to_dict(py, r))
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("reports", reports)?;
    let pairwise: Vec<Bound<'py, PyDict>> = cmp
        .pairwise
        .iter()
        .map(|p| {
            let pd = PyDict::new(py);
            pd.set_item("left_index", p.left_index)?;
            pd.set_item("right_index", p.right_index)?;
            pd.set_item("ecfp4_tanimoto", p.similarities.ecfp4_tanimoto)?;
            pd.set_item("maccs_tanimoto", p.similarities.maccs_tanimoto)?;
            pd.set_item("atom_pair_tanimoto", p.similarities.atom_pair_tanimoto)?;
            Ok(pd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("pairwise", pairwise)?;
    let deltas: Vec<Bound<'py, PyDict>> = cmp
        .descriptor_deltas
        .iter()
        .map(|delta| {
            let dd = PyDict::new(py);
            dd.set_item("left_index", delta.left_index)?;
            dd.set_item("right_index", delta.right_index)?;
            dd.set_item("mw", delta.molecular_weight)?;
            dd.set_item("logp", delta.logp)?;
            dd.set_item("tpsa", delta.tpsa)?;
            dd.set_item("hbd", delta.hbd)?;
            dd.set_item("hba", delta.hba)?;
            dd.set_item("rotatable_bonds", delta.rotatable_bonds)?;
            dd.set_item("qed", delta.qed)?;
            dd.set_item("sa_score", delta.sa_score)?;
            Ok(dd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("descriptor_deltas", deltas)?;
    match &cmp.mcs_smiles {
        Some(s) => d.set_item("mcs_smiles", s)?,
        None => d.set_item("mcs_smiles", py.None())?,
    }
    Ok(d)
}

/// Find all Matched Molecular Pairs (MMP) in a list of SMILES strings.
///
/// Returns a list of dicts, each with keys:
///   ``mol_a``, ``mol_b`` (canonical SMILES), ``core`` (shared scaffold),
///   ``fragment_a``, ``fragment_b`` (substituent SMILES containing ``[*]``).
///
/// Uses BRICS single-bond cuts. Pairs are deduplicated.
///
///     pairs = chematic.find_mmp(["c1ccccc1", "Cc1ccccc1", "Nc1ccccc1"])
///     for p in pairs:
///         print(f"{p['fragment_a']} → {p['fragment_b']} on {p['core']}")
#[pyfunction]
fn find_mmp<'py>(smiles: Vec<String>, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
    let mols: Vec<chematic_core::Molecule> = smiles
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();
    let refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    chematic_chem::find_mmp(&refs)
        .into_iter()
        .map(|pair| {
            let d = PyDict::new(py);
            d.set_item("mol_a", &pair.mol_a)?;
            d.set_item("mol_b", &pair.mol_b)?;
            d.set_item("core", &pair.core)?;
            d.set_item("fragment_a", &pair.fragment_a)?;
            d.set_item("fragment_b", &pair.fragment_b)?;
            Ok(d)
        })
        .collect()
}

/// R-group decomposition — split molecules into a scaffold core and variable R-groups.
///
/// ``scaffold_smarts``: SMARTS pattern defining the common scaffold.
/// ``mols``: list of :class:`Mol` objects to decompose.
///
/// Returns a list of dicts (one per input molecule), or ``None`` when the scaffold
/// does not match a particular molecule.  Each dict contains:
///   ``mol_idx``   (int)  — index in the input list,
///   ``core``      (str)  — scaffold SMILES with ``[*]`` at attachment points,
///   ``R1``, ``R2``, … (str) — SMILES for each R-group (``[*]`` marks attachment).
///
///     mols = [chematic.from_smiles(s) for s in ["CCc1ccccc1", "CCCc1ccccc1"]]
///     results = chematic.rgroup_decompose("c1ccccc1", mols)
///     # [{"mol_idx": 0, "core": "...", "R1": "[*]CC"}, ...]
#[pyfunction]
fn rgroup_decompose<'py>(
    scaffold_smarts: &str,
    mols: Vec<Mol>,
    py: Python<'py>,
) -> PyResult<Vec<Option<Bound<'py, PyDict>>>> {
    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    let results = chematic_chem::rgroup_decompose(scaffold_smarts, &refs)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    results
        .into_iter()
        .map(|opt| match opt {
            None => Ok(None),
            Some(r) => {
                let d = PyDict::new(py);
                d.set_item("mol_idx", r.mol_idx)?;
                d.set_item("core", &r.core_smiles)?;
                for (k, v) in &r.r_groups {
                    d.set_item(format!("R{k}"), v)?;
                }
                Ok(Some(d))
            }
        })
        .collect()
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

/// Parse a Hill-notation molecular formula string into an element count dictionary.
///
/// Returns a ``dict[str, int]`` mapping element symbol → atom count.
/// Mirrors the API of PyPI libraries **chemparse** and **chemformula**.
///
/// Supported syntax:
///   - Simple formulas: ``"H2O"``, ``"C6H12O6"``
///   - Parentheses with multipliers: ``"Ca(OH)2"`` → ``{"Ca":1,"O":2,"H":2}``
///   - SMILES-style brackets: ``"[NH4]+"`` → ``{"N":1,"H":4}``
///   - Trailing charge signs are ignored: ``"NH4+"`` → same as ``"NH4"``
///
/// Raises:
///     ValueError: on empty formula or unbalanced parentheses.
///
///     chematic.parse_formula("C6H12O6")  # {"C": 6, "H": 12, "O": 6}
///     chematic.parse_formula("Ca(OH)2")  # {"Ca": 1, "O": 2, "H": 2}
///     chematic.parse_formula("[NH4]+")   # {"N": 1, "H": 4}
#[pyfunction]
fn parse_formula<'py>(formula: &str, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let counts =
        chematic_chem::parse_formula(formula).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let d = PyDict::new(py);
    let mut sorted: Vec<(String, u32)> = counts.into_iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));
    for (elem, cnt) in sorted {
        d.set_item(elem, cnt)?;
    }
    Ok(d)
}

/// Detect activity cliffs in a set of molecules with known activity values.
///
/// An activity cliff is a structurally similar pair with a large activity difference —
/// a classic signal of SAR sensitivity. Common in MolScore and mol-eval type analyses.
///
/// ``mols``: list of :class:`Mol` objects.
/// ``activities``: list of floats (one per mol), e.g. pIC50 values.
/// ``sim_threshold``: minimum ECFP4 Tanimoto similarity to consider a pair (default 0.65).
/// ``cliff_delta``: minimum ``|activity_i − activity_j|`` to be a cliff (default 2.0).
///
/// Returns a list of dicts sorted by similarity descending, each containing:
///   ``mol_a_idx`` (int), ``mol_b_idx`` (int), ``similarity`` (float), ``activity_delta`` (float).
///
///     mols = [chematic.from_smiles(s) for s in ["c1ccccc1", "Cc1ccccc1"]]
///     cliffs = chematic.activity_cliffs(mols, [5.0, 8.5], sim_threshold=0.0, cliff_delta=2.0)
///     # [{"mol_a_idx": 0, "mol_b_idx": 1, "similarity": 0.xx, "activity_delta": 3.5}]
#[pyfunction]
#[pyo3(signature = (mols, activities, sim_threshold = 0.65, cliff_delta = 2.0))]
fn activity_cliffs<'py>(
    mols: Vec<Mol>,
    activities: Vec<f64>,
    sim_threshold: f32,
    cliff_delta: f64,
    py: Python<'py>,
) -> PyResult<Vec<Bound<'py, PyDict>>> {
    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    let cliffs = chematic_chem::activity_cliffs(&refs, &activities, sim_threshold, cliff_delta);
    cliffs
        .into_iter()
        .map(|c| {
            let d = PyDict::new(py);
            d.set_item("mol_a_idx", c.mol_a_idx)?;
            d.set_item("mol_b_idx", c.mol_b_idx)?;
            d.set_item("similarity", c.similarity)?;
            d.set_item("activity_delta", c.activity_delta)?;
            Ok(d)
        })
        .collect()
}

/// Identify the reaction center: bonds broken/formed and atoms changed.
///
/// Returns a dict with keys:
///   ``broken_bonds`` (list of ``[i, j]`` atom index pairs),
///   ``formed_bonds`` (list of ``[i, j]`` atom index pairs),
///   ``changed_atoms`` (list of atom indices).
///
/// Atom indices use reactant-side numbering.
///
///     rc = chematic.find_reaction_center("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
///     print("broken:", rc["broken_bonds"])
///     print("formed:", rc["formed_bonds"])
#[pyfunction]
fn find_reaction_center<'py>(
    reaction_smiles: &str,
    py: Python<'py>,
) -> PyResult<Bound<'py, PyDict>> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let center = chematic_rxn::find_reaction_center(&rxn);
    let d = PyDict::new(py);
    let broken: Vec<[usize; 2]> = center
        .broken_bonds
        .iter()
        .map(|(a, b)| [a.0 as usize, b.0 as usize])
        .collect();
    let formed: Vec<[usize; 2]> = center
        .formed_bonds
        .iter()
        .map(|(a, b)| [a.0 as usize, b.0 as usize])
        .collect();
    let changed: Vec<usize> = center.changed_atoms.iter().map(|a| a.0 as usize).collect();
    d.set_item("broken_bonds", broken)?;
    d.set_item("formed_bonds", formed)?;
    d.set_item("changed_atoms", changed)?;
    Ok(d)
}

/// Compute atom economy of a reaction (green chemistry metric).
/// E-factor (Environmental Factor) — waste-to-product mass ratio.
///
/// E-factor = waste_mass / product_mass.  Lower is greener.
/// Fine chemicals typically E=5–50; pharmaceuticals E=25–100.
///
///     ef = chematic.e_factor(waste_kg=90.0, product_kg=10.0)  # → 9.0
/// Fast structural hash for deduplication — one int per molecule.
///
/// Molecules with the same canonical graph return identical hashes.
/// Use with :func:`are_identical` to confirm true equivalence (no hash collisions).
///
/// Equivalent to RDKit's ``rdMolHash.MolHash()``.
///
///     seen = set()
///     unique = [m for m in mols if (h := chematic.mol_hash(m)) not in seen and not seen.add(h)]
#[pyfunction]
fn mol_hash(mol: &Mol) -> u64 {
    chematic_chem::mol_hash(&mol.inner)
}

/// Check whether two molecules are graph-isomorphic (exact structural identity).
///
/// More reliable than comparing SMILES strings (which depend on canonicalization).
/// Equivalent to RDKit's ``Chem.MolToInchiKey(m1) == Chem.MolToInchiKey(m2)``.
///
///     assert chematic.are_identical(
///         chematic.from_smiles("c1ccccc1"),
///         chematic.from_smiles("C1=CC=CC=C1"),  # kekulé form
///     )
#[pyfunction]
fn are_identical(mol1: &Mol, mol2: &Mol) -> bool {
    chematic_chem::are_identical(&mol1.inner, &mol2.inner)
}

/// Normalize and re-serialize a reaction SMILES.
///
/// Parses the reaction SMILES and writes it back in canonical form.
/// Useful for standardizing reaction data before storing or comparing.
///
///     canon = chematic.write_reaction("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
#[pyfunction]
fn write_reaction(reaction_smiles: &str) -> PyResult<String> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_rxn::write_reaction(&rxn))
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
    })
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

/// E-factor (Environmental Factor) — waste-to-product mass ratio.
///
/// E-factor = waste_mass / product_mass.  Lower is greener.
/// Fine chemicals typically E=5–50; pharmaceuticals E=25–100.
///
///     ef = chematic.e_factor(waste_kg=90.0, product_kg=10.0)  # → 9.0
#[pyfunction]
fn e_factor(waste_mass: f64, product_mass: f64) -> f64 {
    chematic_rxn::e_factor(waste_mass, product_mass)
}

/// Process Mass Intensity (PMI) — total mass used per unit product mass.
///
/// PMI = (sum of all input masses) / product_mass. Lower is greener.
///
///     pmi = chematic.pmi_rxn([solvent_kg, reagent1_kg, reagent2_kg], product_kg)
#[pyfunction]
fn pmi_rxn(all_masses: Vec<f64>, product_mass: f64) -> f64 {
    chematic_rxn::pmi_rxn(&all_masses, product_mass)
}

/// Reaction Mass Efficiency (RME) — fraction of reactant mass in the product.
///
/// RME = product_mass / sum(reactant_masses). Range [0, 1].
///
///     rme = chematic.reaction_mass_efficiency([reactant1_g, reactant2_g], product_g)
#[pyfunction]
fn reaction_mass_efficiency(reactant_masses: Vec<f64>, product_mass: f64) -> f64 {
    chematic_rxn::reaction_mass_efficiency(&reactant_masses, product_mass)
}

///
/// atom_economy = MW(desired products) / MW(all reactants) × 100.
/// Parse a MDL RXN V2000 file and return the canonical reaction SMILES.
///
/// Raises ``ValueError`` on parse failure.
///
///     rxn_smiles = chematic.from_rxn_file(text)
///     ae = chematic.atom_economy(rxn_smiles)
#[pyfunction]
fn from_rxn_file(text: &str) -> PyResult<String> {
    let rxn =
        chematic_mol::parse_rxn_file(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_rxn::write_reaction(&rxn))
}

/// Convert a reaction SMILES string to MDL RXN V2000 format.
///
/// Raises ``ValueError`` on invalid reaction SMILES.
///
///     block = chematic.to_rxn_file("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
#[pyfunction]
fn to_rxn_file(reaction_smiles: &str) -> PyResult<String> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_mol::write_rxn_file(&rxn))
}

/// A value of 100% means all atoms in reactants appear in the product.
///
///     ae = chematic.atom_economy("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
#[pyfunction]
fn atom_economy(reaction_smiles: &str) -> PyResult<f64> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_rxn::atom_economy(&rxn))
}

/// Check whether a reaction SMILES is atom-balanced.
///
/// Returns a dict with keys:
///   ``balanced`` (bool), ``diff`` (list of str describing imbalances).
///
///     result = chematic.balance_check("C+O>>CO")
///     print(result["balanced"], result["diff"])
#[pyfunction]
fn balance_check<'py>(reaction_smiles: &str, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let result = chematic_rxn::balance_check(&rxn);
    let d = PyDict::new(py);
    d.set_item("balanced", result.balanced)?;
    d.set_item("diff", result.diff())?;
    Ok(d)
}

/// Enumerate a combinatorial library from a SMIRKS template and fragment sets.
///
/// Args:
///     smirks: Reaction SMIRKS template (e.g. ``"[C:1]Cl.[N:2]>>[C:1][N:2]"``).
///     fragment_sets: List of SMILES lists — one list per reactant slot.
///                    All combinations across sets are generated.
///     max_size: Maximum library size (default 1_000_000).
///
/// Returns a list of product SMILES strings.
///
///     products = chematic.enumerate_library(
///         "[C:1]C(=O)Cl.[N:2]>>[C:1]C(=O)[N:2]",
///         [["c1ccccc1", "CC"], ["N", "CN"]],
///     )
#[pyfunction]
#[pyo3(signature = (smirks, fragment_sets, max_size = 1_000_000))]
fn enumerate_library(
    smirks: &str,
    fragment_sets: Vec<Vec<String>>,
    max_size: usize,
) -> PyResult<Vec<String>> {
    let parsed_sets: Vec<Vec<chematic_core::Molecule>> = fragment_sets
        .iter()
        .map(|set| {
            set.iter()
                .filter_map(|s| chematic_smiles::parse(s).ok())
                .collect()
        })
        .collect();
    let config = chematic_rxn::LibraryConfig {
        skip_failures: true,
        max_size: Some(max_size),
    };
    chematic_rxn::enumerate_library(smirks, parsed_sets, &config)
        .map(|mols| mols.iter().map(chematic_smiles::canonical_smiles).collect())
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Enumerate a 2-fragment combinatorial library (scaffold × building block).
///
/// Convenience alternative to ``enumerate_library(smirks, [scaffolds, building_blocks])``.
/// The most common combinatorial chemistry pattern: one scaffold set reacted with one
/// building-block set to produce all pairwise products.
///
///     products = chematic.enumerate_library_2way(
///         "[C:1]C(=O)Cl.[N:2]>>[C:1]C(=O)[N:2]",
///         scaffolds=["c1ccccc1C(=O)Cl", "CC(=O)Cl"],
///         building_blocks=["N", "CN"],
///     )
#[pyfunction]
#[pyo3(signature = (smirks, scaffolds, building_blocks, max_size = 1_000_000))]
fn enumerate_library_2way(
    smirks: &str,
    scaffolds: Vec<String>,
    building_blocks: Vec<String>,
    max_size: usize,
) -> PyResult<Vec<String>> {
    let parse_smiles_set = |set: Vec<String>| -> Vec<chematic_core::Molecule> {
        set.iter()
            .filter_map(|s| chematic_smiles::parse(s).ok())
            .collect()
    };
    let config = chematic_rxn::LibraryConfig {
        skip_failures: true,
        max_size: Some(max_size),
    };
    chematic_rxn::enumerate_library_2way(
        smirks,
        parse_smiles_set(scaffolds),
        parse_smiles_set(building_blocks),
        &config,
    )
    .map(|mols| mols.iter().map(chematic_smiles::canonical_smiles).collect())
    .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Enumerate a 3-fragment combinatorial library (scaffold × R1 × R2).
///
/// Convenience alternative to ``enumerate_library(smirks, [scaffolds, r1_set, r2_set])``.
/// Covers the common scaffold-decoration pattern with two variable positions.
///
///     products = chematic.enumerate_library_3way(
///         "[C:1]C(=O)Cl.[N:2]>>[C:1]C(=O)[N:2]",
///         scaffolds=["CC(=O)Cl"],
///         r1_set=["N", "CN"],
///         r2_set=["c1ccccc1", "CC"],
///     )
#[pyfunction]
#[pyo3(signature = (smirks, scaffolds, r1_set, r2_set, max_size = 1_000_000))]
fn enumerate_library_3way(
    smirks: &str,
    scaffolds: Vec<String>,
    r1_set: Vec<String>,
    r2_set: Vec<String>,
    max_size: usize,
) -> PyResult<Vec<String>> {
    let parse_smiles_set = |set: Vec<String>| -> Vec<chematic_core::Molecule> {
        set.iter()
            .filter_map(|s| chematic_smiles::parse(s).ok())
            .collect()
    };
    let config = chematic_rxn::LibraryConfig {
        skip_failures: true,
        max_size: Some(max_size),
    };
    chematic_rxn::enumerate_library_3way(
        smirks,
        parse_smiles_set(scaffolds),
        parse_smiles_set(r1_set),
        parse_smiles_set(r2_set),
        &config,
    )
    .map(|mols| mols.iter().map(chematic_smiles::canonical_smiles).collect())
    .map_err(|e| PyValueError::new_err(e.to_string()))
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

/// Find top-K nearest neighbors from precomputed fingerprint byte arrays.
///
/// More efficient than :func:`top_k_similar_fp` when the same ``db_fps`` list is
/// reused across multiple queries (fingerprints computed only once).
///
///     db_fps = [mol.ecfp4() for mol in library]   # compute once
///     for query in queries:
///         hits = chematic.nearest_neighbors_from_fp(query.ecfp4(), db_fps, k=10)
///         for idx, score in hits:
///             print(library_smiles[idx], score)
#[pyfunction]
#[pyo3(signature = (query_fp, db_fps, k = 10))]
fn nearest_neighbors_from_fp(query_fp: &[u8], db_fps: Vec<Vec<u8>>, k: usize) -> Vec<(usize, f64)> {
    let qa: u32 = query_fp.iter().map(|b| b.count_ones()).sum();
    let mut scores: Vec<(usize, f64)> = db_fps
        .iter()
        .enumerate()
        .filter_map(|(i, fp)| {
            if fp.len() != query_fp.len() {
                return None;
            }
            let and: u32 = query_fp
                .iter()
                .zip(fp.iter())
                .map(|(a, b)| (a & b).count_ones())
                .sum();
            let db_cnt: u32 = fp.iter().map(|b| b.count_ones()).sum();
            let or = qa + db_cnt - and;
            if or == 0 {
                return None;
            }
            Some((i, and as f64 / or as f64))
        })
        .collect();
    scores.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    scores.truncate(k);
    scores
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

/// Compute Tanimoto similarity between two reaction SMILES using reaction fingerprints.
///
/// Parses both reaction SMILES, computes reaction fingerprints, and returns
/// the Tanimoto coefficient.
///
///     sim = chematic.tanimoto_reaction_fp("CC>>CO", "c1ccccc1>>c1ccccc1N")
///
/// Raises ``ValueError`` on invalid reaction SMILES.
#[pyfunction]
fn tanimoto_reaction_fp(rxn1: &str, rxn2: &str) -> PyResult<f64> {
    let r1 =
        chematic_rxn::parse_reaction(rxn1).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let r2 =
        chematic_rxn::parse_reaction(rxn2).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_fp::tanimoto_reaction_fp(&r1, &r2))
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

/// Print a version and accuracy summary — useful for debugging and reporting.
///
/// ```python
/// chematic.doctor()
/// # chematic v0.4.21
/// # Python 3.13  |  darwin arm64
/// # ...
/// ```
#[pyfunction]
fn doctor(py: Python<'_>) {
    let ver = env!("CARGO_PKG_VERSION");
    let vi = py.version_info();
    let py_ver = format!("{}.{}.{}", vi.major, vi.minor, vi.patch);

    let platform = std::env::consts::OS;
    let arch = std::env::consts::ARCH;

    println!("chematic v{ver}");
    println!("Python {py_ver}  |  {platform} {arch}");
    println!();
    println!("Descriptor accuracy (benchmark 2026-06, v0.4.22 vs RDKit 2026.03.3):");
    println!("  MW / HBA / HBD / ARC  100%   (4,999-mol ChEMBL subset)");
    println!("  TPSA                  98.1%  (4,999-mol, ±0.1 Å²)");
    println!("  LogP (Crippen)        96.5%  (4,999-mol, ±0.01) · ~99% (175-mol, ±0.3)");
    println!("  ECFP4 throughput      3.6 µs/mol  (5–14× faster than RDKit)");
    println!("  WASM bundle           504 KB gzip");
    println!();
    println!("Feature stability:");
    println!("  Stable      SMILES · MW/HBA/HBD/TPSA/LogP · ECFP4/MACCS · SDF/MOL · SMARTS");
    println!("  Stable      Tanimoto · PAINS/Brenk · 2D SVG · QED");
    println!("  Experimental   3D conformer (not RDKit ETKDGv3 equivalent)");
    println!("  Rule-based     pKa · ADMET (screening use only, not clinical)");
    println!("  Partial        IUPAC name generation · pure-Rust InChI (approx.)");
    println!();
    println!("Docs:       https://kent-tokyo.github.io/chematic/");
    println!("Validation: https://github.com/kent-tokyo/chematic/tree/main/validation/");
    println!("Benchmarks: https://github.com/kent-tokyo/chematic/tree/main/benchmarks/");
}


/// Parse a ``.smi`` file (tab/space-separated SMILES + name) into (Mol, name) pairs.
///
/// Each line is ``SMILES[<tab>name]``. Lines with invalid SMILES are silently skipped.
/// Comment lines starting with ``#`` and blank lines are ignored.
/// Equivalent to RDKit's ``Chem.SmilesMolSupplier``.
///
///     records = chematic.parse_smi_file(open("library.smi").read())
///     for mol, name in records:
///         print(name, mol.mw)
#[pyfunction]
fn parse_smi_file(content: &str) -> Vec<(Mol, String)> {
    chematic_smiles::parse_smi_file(content)
        .into_iter()
        .filter_map(|r| r.ok())
        .map(|(mol, name)| {
            (
                Mol {
                    inner: Arc::new(mol),
                },
                name,
            )
        })
        .collect()
}

/// Write (Mol, name) pairs to ``.smi`` format.
///
/// Output format: ``SMILES<TAB>name<NEWLINE>`` per record (name omitted if empty).
/// Equivalent to RDKit's ``Chem.SmilesWriter``.
///
///     text = chematic.write_smi_file([(mol1, "cpd1"), (mol2, "cpd2")])
///     with open("output.smi", "w") as f:
///         f.write(text)
#[pyfunction]
fn write_smi_file(records: Vec<(Mol, String)>) -> String {
    let mut out = String::new();
    for (mol, name) in &records {
        let smiles = chematic_smiles::canonical_smiles(&mol.inner);
        if name.is_empty() {
            out.push_str(&smiles);
        } else {
            out.push_str(&smiles);
            out.push('\t');
            out.push_str(name);
        }
        out.push('\n');
    }
    out
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

/// Check whether a reaction SMILES matches a reaction SMARTS pattern.
///
/// Returns ``True`` if the reaction matches the query pattern, ``False`` otherwise.
/// Equivalent to ``chematic.reaction_smarts_match`` but returns a simple bool
/// via SMARTS-based pattern matching rather than substructure query.
///
///     matched = chematic.query_reaction("CC>>CO", "[C:1]>>[C:1]O")
///
/// Raises ``ValueError`` on invalid reaction SMILES or SMARTS.
#[pyfunction]
fn query_reaction(reaction_smiles: &str, smarts: &str) -> PyResult<bool> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let result = chematic_rxn::query_reaction(&rxn, smarts)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(result.is_complete_match)
}

/// Query a list of reaction SMILES against a single SMARTS pattern.
///
/// Returns a dict:
///   - ``total`` (int): total reactions processed
///   - ``matching`` (int): reactions that matched
///   - ``match_pct`` (float): match percentage (0–100)
///   - ``matches`` (list[(int, bool)]): per-reaction results as (original_index, matched)
///
/// Invalid SMILES are silently skipped (their indices will not appear in ``matches``).
///
/// Raises ``ValueError`` on invalid SMARTS.
///
///     rxns = ["CC>>CO", "CCCC>>CCCCO", "c1ccccc1>>c1ccccc1N"]
///     r = chematic.batch_query_reactions(rxns, "[C:1]>>[C:1]O")
///     print(r["matching"], "/", r["total"])
#[pyfunction]
fn batch_query_reactions<'py>(
    reactions: Vec<String>,
    smarts: &str,
    py: Python<'py>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut valid: Vec<chematic_rxn::Reaction> = Vec::new();
    let mut original_indices: Vec<usize> = Vec::new();
    for (i, s) in reactions.iter().enumerate() {
        if let Ok(rxn) = chematic_rxn::parse_reaction(s) {
            valid.push(rxn);
            original_indices.push(i);
        }
    }
    let result = chematic_rxn::batch_query_reactions(&valid, smarts)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    // Remap indices from filtered slice back to original input indices
    let matches: Vec<(usize, bool)> = result
        .matches
        .iter()
        .map(|&(idx, matched)| (original_indices[idx], matched))
        .collect();
    let d = PyDict::new(py);
    d.set_item("total", result.total_reactions)?;
    d.set_item("matching", result.matching_reactions)?;
    d.set_item("match_pct", result.match_percentage)?;
    d.set_item("matches", matches)?;
    Ok(d)
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

/// Render a reaction SMILES as an SVG diagram.
///
/// Returns an SVG string showing reactants → products with an arrow.
/// Equivalent to RDKit's ``Draw.ReactionToImage(rxn)``.
///
///     svg = chematic.reaction_svg("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
///     with open("reaction.svg", "w") as f:
///         f.write(svg)
///
/// Raises ``ValueError`` on invalid reaction SMILES.
#[pyfunction]
fn reaction_svg(reaction_smiles: &str) -> PyResult<String> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_depict::depict_reaction_svg(&rxn))
}

/// Compute scaffold network statistics across a molecule library.
///
/// Returns a dict with three parallel lists:
///   - ``scaffolds``: canonical SMILES of each unique scaffold
///   - ``counts``: how many input molecules contain each scaffold
///   - ``parents``: index of the parent (simpler) scaffold, or ``None`` for root
///
/// Invalid SMILES are silently skipped.
///
///     result = chematic.scaffold_network_counts(smiles_list)
///     for smi, n in zip(result["scaffolds"], result["counts"]):
///         print(smi, n)
#[pyfunction]
fn scaffold_network_counts<'py>(
    smiles: Vec<String>,
    py: Python<'py>,
) -> PyResult<Bound<'py, PyDict>> {
    let mols: Vec<chematic_core::Molecule> = smiles
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();
    let net = chematic_chem::scaffold_network_with_counts(&mols);
    let scaffolds: Vec<String> = net
        .scaffolds
        .iter()
        .map(chematic_smiles::canonical_smiles)
        .collect();
    let d = PyDict::new(py);
    d.set_item("scaffolds", scaffolds)?;
    d.set_item("counts", net.counts)?;
    d.set_item("parents", net.parents)?;
    Ok(d)
}

/// Render a list of molecules as a grid SVG.
///
///     svg = chematic.depict_grid([mol1, mol2, mol3], cols=3)
#[pyfunction]
fn depict_grid(mols: Vec<Mol>, cols: usize) -> String {
    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    chematic_depict::depict_svg_grid(&refs, cols)
}

/// Generate a self-contained HTML report for a list of molecules.
///
/// Returns an HTML string. If ``output`` is given, also writes the file.
/// Cards are sorted by QED descending (most drug-like first).
///
/// ```python
/// mols = [chematic.from_smiles(s) for s in ["CC(=O)Oc1ccccc1C(=O)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"]]
/// html = chematic.report(mols, names=["aspirin", "caffeine"], output="report.html")
/// ```
#[pyfunction]
#[pyo3(signature = (mols, names=None, title="chematic report", output=None))]
fn report(
    mols: Vec<Mol>,
    names: Option<Vec<Option<String>>>,
    title: &str,
    output: Option<&str>,
) -> PyResult<Report> {
    use chematic_chem::{
        brenk_passes, hbd_count, logp_and_mr, molecular_weight, pains_passes,
        qed_with_bundle, ring_bundle, tpsa,
    };

    // Build (qed, card_html) pairs so we can sort by QED
    let mut cards: Vec<(f64, String)> = mols
        .iter()
        .enumerate()
        .map(|(i, mol)| {
            let m = mol.inner.as_ref();
            let mw = molecular_weight(m);
            let (logp, _) = logp_and_mr(m);
            let tpsa_val = tpsa(m);
            let hbd = hbd_count(m);
            let rb = ring_bundle(m);
            let qed = qed_with_bundle(m, &rb);
            let lip = mw <= 500.0 && hbd <= 5 && rb.hba_count <= 10 && logp <= 5.0;
            let pains_ok = pains_passes(m);
            let brenk_ok = brenk_passes(m);

            let label = names
                .as_ref()
                .and_then(|ns| ns.get(i))
                .and_then(|n| n.as_deref())
                .unwrap_or("");

            let svg = chematic_depict::depict_svg(m);

            let lip_badge = if lip {
                r#"<span class="badge pass">Lipinski ✓</span>"#
            } else {
                r#"<span class="badge fail">Lipinski ✗</span>"#
            };
            let pains_badge = if pains_ok {
                r#"<span class="badge pass">PAINS ✓</span>"#
            } else {
                r#"<span class="badge fail">PAINS ✗</span>"#
            };
            let brenk_badge = if brenk_ok {
                r#"<span class="badge pass">Brenk ✓</span>"#
            } else {
                r#"<span class="badge warn">Brenk ⚠</span>"#
            };

            let card = format!(
                r#"<div class="card"><div class="svg">{svg}</div><div class="name">{label}</div><div class="desc">MW: {mw:.1} Da &nbsp;|&nbsp; LogP: {logp:.2} &nbsp;|&nbsp; TPSA: {tpsa_val:.1} Å²<br>HBD: {hbd} &nbsp;|&nbsp; HBA: {hba} &nbsp;|&nbsp; QED: {qed:.2}</div><div class="badges">{lip_badge}{pains_badge}{brenk_badge}</div></div>"#,
                hba = rb.hba_count,
            );
            (qed, card)
        })
        .collect();

    cards.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let cards_html: String = cards.into_iter().map(|(_, c)| c).collect();
    let n = mols.len();
    let ver = env!("CARGO_PKG_VERSION");

    let html = format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
body{{font-family:system-ui,sans-serif;background:#f8f9fa;padding:24px;margin:0}}
h1{{font-size:1.4rem;color:#333;margin-bottom:4px}}
.meta{{font-size:.85rem;color:#666;margin-bottom:20px}}
.grid{{display:flex;flex-wrap:wrap;gap:16px}}
.card{{background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:12px;width:220px;box-shadow:0 1px 3px rgba(0,0,0,.06)}}
.svg{{width:100%;height:160px;overflow:hidden}}
.svg svg{{width:100%;height:100%}}
.name{{font-weight:600;font-size:.9rem;margin:6px 0 4px;color:#333;min-height:1.1em}}
.desc{{font-size:.78rem;color:#555;line-height:1.8;margin:4px 0}}
.badges{{display:flex;flex-wrap:wrap;gap:4px;margin-top:6px}}
.badge{{font-size:.7rem;padding:2px 7px;border-radius:10px;font-weight:500}}
.pass{{background:#d1e7dd;color:#0a3622}}
.fail{{background:#f8d7da;color:#58151c}}
.warn{{background:#fff3cd;color:#664d03}}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="meta">{n} molecule{plural} &middot; generated by chematic v{ver}</p>
<div class="grid">{cards_html}</div>
</body>
</html>"#,
        plural = if n == 1 { "" } else { "s" },
    );

    let rep = Report { html };
    if let Some(path) = output {
        rep.save(path)?;
    }
    Ok(rep)
}

/// Compare two molecules side-by-side and return a self-contained HTML report.
///
/// Returns a ``Report`` object. In Jupyter, writing ``report`` renders it inline.
///
/// ```python
/// aspirin   = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
/// ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")
/// report = chematic.compare(aspirin, ibuprofen, names=("Aspirin", "Ibuprofen"))
/// report.save("compare.html")
/// ```
#[pyfunction]
#[pyo3(signature = (mol1, mol2, names=None, title=None))]
fn compare(
    mol1: &Mol,
    mol2: &Mol,
    names: Option<(String, String)>,
    title: Option<&str>,
) -> Report {
    use chematic_chem::{
        brenk_passes, hbd_count, logp_and_mr, molecular_weight, pains_passes,
        qed_with_bundle, ring_bundle, tpsa,
    };

    let m1 = mol1.inner.as_ref();
    let m2 = mol2.inner.as_ref();

    let (name1, name2) = names
        .unwrap_or_else(|| ("Molecule A".into(), "Molecule B".into()));

    let heading = title
        .map(|t| t.to_string())
        .unwrap_or_else(|| format!("{name1} vs {name2}"));

    let svg1 = chematic_depict::depict_svg(m1);
    let svg2 = chematic_depict::depict_svg(m2);

    let mw1 = molecular_weight(m1);
    let mw2 = molecular_weight(m2);
    let (logp1, _) = logp_and_mr(m1);
    let (logp2, _) = logp_and_mr(m2);
    let tpsa1 = tpsa(m1);
    let tpsa2 = tpsa(m2);
    let hbd1 = hbd_count(m1);
    let hbd2 = hbd_count(m2);
    let rb1 = ring_bundle(m1);
    let rb2 = ring_bundle(m2);
    let qed1 = qed_with_bundle(m1, &rb1);
    let qed2 = qed_with_bundle(m2, &rb2);
    let lip1 = mw1 <= 500.0 && hbd1 <= 5 && rb1.hba_count <= 10 && logp1 <= 5.0;
    let lip2 = mw2 <= 500.0 && hbd2 <= 5 && rb2.hba_count <= 10 && logp2 <= 5.0;
    let pains1 = pains_passes(m1);
    let pains2 = pains_passes(m2);
    let brenk1 = brenk_passes(m1);
    let brenk2 = brenk_passes(m2);

    // MCS common atoms (reuse diff logic)
    let config = chematic_smarts::McsConfig::default();
    let qmol = chematic_smarts::find_mcs_with_config(&[m1, m2], &config);
    let common = qmol.atoms.len();

    // Delta summary (same logic as Mol::diff)
    let elem_parts: Vec<String> = {
        use std::collections::BTreeMap;
        let mut c1: BTreeMap<&'static str, i32> = BTreeMap::new();
        let mut c2: BTreeMap<&'static str, i32> = BTreeMap::new();
        for i in 0..m1.atom_count() {
            *c1.entry(m1.atom(chematic_core::AtomIdx(i as u32)).element.symbol()).or_insert(0) += 1;
        }
        for i in 0..m2.atom_count() {
            *c2.entry(m2.atom(chematic_core::AtomIdx(i as u32)).element.symbol()).or_insert(0) += 1;
        }
        let all: std::collections::BTreeSet<_> = c1.keys().chain(c2.keys()).copied().collect();
        all.iter()
            .filter_map(|e| {
                let d = c2.get(e).copied().unwrap_or(0) - c1.get(e).copied().unwrap_or(0);
                if d != 0 { Some(if d > 0 { format!("+{d}{e}") } else { format!("{d}{e}") }) } else { None }
            })
            .collect()
    };
    let elem_str = if elem_parts.is_empty() { "Same elemental composition".into() } else { elem_parts.join(", ") };
    let summary = format!(
        "{}. \u{0394}LogP {:+.2}, \u{0394}TPSA {:+.1} \u{00c5}\u{00b2}, \u{0394}MW {:+.1} Da.",
        elem_str, logp2 - logp1, tpsa2 - tpsa1, mw2 - mw1,
    );

    fn delta_class(d: f64) -> &'static str {
        if d > 0.0 { "pos" } else if d < 0.0 { "neg" } else { "" }
    }
    fn flag(v: bool, ok: &str, fail: &str) -> String {
        if v { format!(r#"<span class="pass">{ok}</span>"#) } else { format!(r#"<span class="fail">{fail}</span>"#) }
    }
    fn warn_flag(v: bool) -> String {
        if v { r#"<span class="pass">✓</span>"#.into() } else { r#"<span class="warn">⚠</span>"#.into() }
    }

    let ver = env!("CARGO_PKG_VERSION");
    let html = format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{heading}</title>
<style>
body{{font-family:system-ui,sans-serif;background:#f8f9fa;padding:24px;margin:0}}
h1{{font-size:1.4rem;color:#333;margin-bottom:4px}}
.meta{{font-size:.85rem;color:#666;margin-bottom:6px}}
.summary{{font-size:.85rem;color:#444;background:#fff;border:1px solid #dee2e6;border-radius:6px;padding:8px 12px;margin-bottom:20px;display:inline-block}}
table{{border-collapse:collapse;background:#fff;border-radius:8px;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.06)}}
th,td{{padding:10px 16px;text-align:left;border-bottom:1px solid #f0f0f0;font-size:.88rem}}
th{{background:#f8f9fa;font-weight:600;color:#555}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}
td.delta{{text-align:right;font-size:.8rem;font-weight:600}}
.pos{{color:#0a6640}}
.neg{{color:#8b1c1c}}
.pass{{color:#0a6640;font-weight:600}}
.fail{{color:#8b1c1c;font-weight:600}}
.warn{{color:#7d5a00;font-weight:600}}
.svg-cell svg{{width:180px;height:140px}}
</style>
</head>
<body>
<h1>{heading}</h1>
<p class="meta">Common scaffold: {common} atoms &middot; generated by chematic v{ver}</p>
<p class="summary">{summary}</p>
<table>
<tr><th>Property</th><th>{name1}</th><th>{name2}</th><th>Delta</th></tr>
<tr><td>Structure</td>
    <td class="svg-cell">{svg1}</td>
    <td class="svg-cell">{svg2}</td>
    <td></td></tr>
<tr><td>MW (Da)</td>
    <td class="num">{mw1:.1}</td><td class="num">{mw2:.1}</td>
    <td class="delta {dc_mw}">{dmw:+.1}</td></tr>
<tr><td>LogP</td>
    <td class="num">{logp1:.2}</td><td class="num">{logp2:.2}</td>
    <td class="delta {dc_lp}">{dlp:+.2}</td></tr>
<tr><td>TPSA (Å²)</td>
    <td class="num">{tpsa1:.1}</td><td class="num">{tpsa2:.1}</td>
    <td class="delta {dc_tp}">{dtp:+.1}</td></tr>
<tr><td>HBD</td>
    <td class="num">{hbd1}</td><td class="num">{hbd2}</td>
    <td class="delta {dc_hbd}">{dhbd:+}</td></tr>
<tr><td>HBA</td>
    <td class="num">{hba1}</td><td class="num">{hba2}</td>
    <td class="delta {dc_hba}">{dhba:+}</td></tr>
<tr><td>QED</td>
    <td class="num">{qed1:.2}</td><td class="num">{qed2:.2}</td>
    <td class="delta {dc_qed}">{dqed:+.2}</td></tr>
<tr><td>Lipinski</td>
    <td>{lip1_s}</td><td>{lip2_s}</td><td></td></tr>
<tr><td>PAINS</td>
    <td>{pains1_s}</td><td>{pains2_s}</td><td></td></tr>
<tr><td>Brenk</td>
    <td>{brenk1_s}</td><td>{brenk2_s}</td><td></td></tr>
</table>
</body>
</html>"#,
        dc_mw  = delta_class(mw2 - mw1),   dmw  = mw2 - mw1,
        dc_lp  = delta_class(logp2-logp1),  dlp  = logp2 - logp1,
        dc_tp  = delta_class(tpsa2-tpsa1),  dtp  = tpsa2 - tpsa1,
        dc_hbd = delta_class((hbd2 as f64)-(hbd1 as f64)), dhbd = hbd2 as i32 - hbd1 as i32,
        dc_hba = delta_class((rb2.hba_count as f64)-(rb1.hba_count as f64)),
        dhba   = rb2.hba_count as i32 - rb1.hba_count as i32,
        hba1   = rb1.hba_count,
        hba2   = rb2.hba_count,
        dc_qed = delta_class(qed2-qed1),    dqed = qed2 - qed1,
        lip1_s  = flag(lip1, "✓ Lipinski", "✗ Lipinski"),
        lip2_s  = flag(lip2, "✓ Lipinski", "✗ Lipinski"),
        pains1_s = flag(pains1, "✓ PAINS", "✗ PAINS"),
        pains2_s = flag(pains2, "✓ PAINS", "✗ PAINS"),
        brenk1_s = warn_flag(brenk1),
        brenk2_s = warn_flag(brenk2),
    );

    Report { html }
}

/// Apply a SMIRKS reaction template to a list of reactant molecules.
///
/// Returns a list of product sets; each set is a list of Mol.
/// Raises ``ValueError`` on SMIRKS parse failure or reactant count mismatch.
///
/// **Stereochemistry**: when the reactant template contains ``@``/``@@`` stereo
/// descriptors, only reactant molecules whose chiral centres match the template
/// configuration are accepted (parity-aware comparison, SMILES write-order
/// independent). Templates without stereo descriptors match both enantiomers.
///
///     products = chematic.run_smirks("[OH:1]>>[O-:1]", [mol])
///     # → [[product_mol], ...]
///
///     # Stereo-selective: only L-amino acids match this template
///     l_products = chematic.run_smirks("[N:1][C@@H:2](C)C(=O)O>>[N:1].[C@@H:2](C)C(=O)O", [mol])
#[pyfunction]
fn run_smirks(smirks: &str, reactants: Vec<Mol>) -> PyResult<Vec<Vec<Mol>>> {
    for mol in &reactants {
        if mol.inner.atom_count() > 300 {
            return Err(PyValueError::new_err(
                "reactant too large for run_smirks (max 300 heavy atoms)",
            ));
        }
    }
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

/// Apply a SMIRKS reaction template to a list of reactant molecules (strict mode).
///
/// Like :func:`run_smirks` but **does not carry substituents** into products.
/// Only atoms explicitly mapped in the product template are included.
/// Stereo filtering behaviour is identical to :func:`run_smirks`.
///
///     products = chematic.run_smirks_strict("[N:1][C:2]>>[N:1].[C:2]", [mol])
///     # → only the mapped N and C atoms; no R-groups attached
#[pyfunction]
fn run_smirks_strict(smirks: &str, reactants: Vec<Mol>) -> PyResult<Vec<Vec<Mol>>> {
    for mol in &reactants {
        if mol.inner.atom_count() > 300 {
            return Err(PyValueError::new_err(
                "reactant too large for run_smirks_strict (max 300 heavy atoms)",
            ));
        }
    }
    let refs: Vec<&chematic_core::Molecule> = reactants.iter().map(|m| m.inner.as_ref()).collect();
    chematic_rxn::run_reactants_strict(smirks, &refs)
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
fn find_mcs(
    mols: Vec<Mol>,
    ring_matches_ring_only: bool,
    complete_rings_only: bool,
) -> Option<Mol> {
    use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
    use chematic_smarts::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, McsConfig};

    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    let config = McsConfig {
        ring_matches_ring_only,
        complete_rings_only,
        ..McsConfig::default()
    };
    let qmol = chematic_smarts::find_mcs_with_config(&refs, &config);

    if qmol.atoms.is_empty() {
        return None;
    }

    fn extract_atomic_num(q: &AtomQuery) -> Option<u8> {
        match q {
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)) => Some(*n),
            AtomQuery::And(lhs, rhs) => extract_atomic_num(lhs).or_else(|| extract_atomic_num(rhs)),
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
    Some(Mol {
        inner: Arc::new(builder.build()),
    })
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
    m.add_function(wrap_pyfunction!(from_cxsmiles, m)?)?;
    m.add_function(wrap_pyfunction!(named_pattern, m)?)?;
    m.add_function(wrap_pyfunction!(doctor, m)?)?;
    m.add_function(wrap_pyfunction!(parse_smi_file, m)?)?;
    m.add_function(wrap_pyfunction!(write_smi_file, m)?)?;
    m.add_function(wrap_pyfunction!(atom_color, m)?)?;
    m.add_function(wrap_pyfunction!(atom_color_rgb, m)?)?;
    m.add_function(wrap_pyfunction!(from_cml, m)?)?;
    m.add_function(wrap_pyfunction!(from_cjson, m)?)?;
    m.add_function(wrap_pyfunction!(from_cdxml, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_v3000, m)?)?;
    m.add_function(wrap_pyfunction!(from_condensed, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_block, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_block_with_coords, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_v3000_with_coords, m)?)?;
    m.add_function(wrap_pyfunction!(parse_sdf_with_coords, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol2, m)?)?;
    m.add_function(wrap_pyfunction!(from_pdbqt, m)?)?;
    m.add_function(wrap_pyfunction!(from_gjf, m)?)?;
    m.add_function(wrap_pyfunction!(parse_gaussian_log, m)?)?;
    m.add_function(wrap_pyfunction!(write_gjf, m)?)?;
    m.add_function(wrap_pyfunction!(parse_cif, m)?)?;
    m.add_function(wrap_pyfunction!(from_inchi, m)?)?;
    m.add_function(wrap_pyfunction!(is_valid_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(is_valid_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_erg, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(nearest_neighbors_from_fp, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_slice, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_reaction_fp, m)?)?;
    m.add_function(wrap_pyfunction!(smarts_match, m)?)?;
    m.add_function(wrap_pyfunction!(smarts_find, m)?)?;
    m.add_function(wrap_pyfunction!(depict_grid, m)?)?;
    m.add_class::<Report>()?;
    m.add_function(wrap_pyfunction!(report, m)?)?;
    m.add_function(wrap_pyfunction!(compare, m)?)?;
    m.add_function(wrap_pyfunction!(reaction_svg, m)?)?;
    m.add_function(wrap_pyfunction!(scaffold_network_counts, m)?)?;
    m.add_function(wrap_pyfunction!(run_smirks, m)?)?;
    m.add_function(wrap_pyfunction!(run_smirks_strict, m)?)?;
    m.add_function(wrap_pyfunction!(find_mcs, m)?)?;
    m.add_function(wrap_pyfunction!(reaction_smarts_match, m)?)?;
    m.add_function(wrap_pyfunction!(query_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(batch_query_reactions, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_pharmacophore_3d, m)?)?;
    m.add_function(wrap_pyfunction!(shape_screen, m)?)?;
    m.add_function(wrap_pyfunction!(from_pdb, m)?)?;
    m.add_function(wrap_pyfunction!(from_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_map4, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_mhfp, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_spectrophores, m)?)?;
    m.add_function(wrap_pyfunction!(butina_cluster, m)?)?;
    m.add_function(wrap_pyfunction!(maxmin_picks, m)?)?;
    m.add_function(wrap_pyfunction!(molecule_report, m)?)?;
    m.add_function(wrap_pyfunction!(screen_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(compare_molecules, m)?)?;
    m.add_function(wrap_pyfunction!(find_mmp, m)?)?;
    m.add_function(wrap_pyfunction!(rgroup_decompose, m)?)?;
    m.add_function(wrap_pyfunction!(similarity_map_svg, m)?)?;
    m.add_function(wrap_pyfunction!(activity_cliffs, m)?)?;
    m.add_function(wrap_pyfunction!(parse_formula, m)?)?;
    m.add_function(wrap_pyfunction!(find_reaction_center, m)?)?;
    m.add_function(wrap_pyfunction!(mol_hash, m)?)?;
    m.add_function(wrap_pyfunction!(are_identical, m)?)?;
    m.add_function(wrap_pyfunction!(write_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(abbreviations, m)?)?;
    m.add_function(wrap_pyfunction!(expand_abbreviation, m)?)?;
    m.add_function(wrap_pyfunction!(top_k_similar, m)?)?;
    m.add_function(wrap_pyfunction!(top_k_similar_fp, m)?)?;
    m.add_function(wrap_pyfunction!(center_on_origin, m)?)?;
    m.add_function(wrap_pyfunction!(transform_conformer, m)?)?;
    m.add_function(wrap_pyfunction!(e_factor, m)?)?;
    m.add_function(wrap_pyfunction!(pmi_rxn, m)?)?;
    m.add_function(wrap_pyfunction!(reaction_mass_efficiency, m)?)?;
    m.add_function(wrap_pyfunction!(from_rxn_file, m)?)?;
    m.add_function(wrap_pyfunction!(to_rxn_file, m)?)?;
    m.add_function(wrap_pyfunction!(atom_economy, m)?)?;
    m.add_function(wrap_pyfunction!(balance_check, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_library, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_library_2way, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_library_3way, m)?)?;
    m.add_function(wrap_pyfunction!(dice_similarity, m)?)?;
    m.add_function(wrap_pyfunction!(tversky_similarity, m)?)?;
    m.add_function(wrap_pyfunction!(cosine_erg_vec, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_erg_vec, m)?)?;
    m.add_function(wrap_pyfunction!(align_coords, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, m)?)?;
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
