//! The `#[pymethods] impl Mol` block — all instance methods/properties on `Mol`.
//!
//! Kept as a single `impl` block (not further split) because PyO3 0.29 requires
//! the `multiple-pymethods` Cargo feature (untested in this crate) to spread
//! `#[pymethods]` for one `#[pyclass]` across more than one block.

use crate::EcfpBitInfo;
use crate::Mol;
use crate::RdkitMorganDetail;
use crate::formats::{bitvec2048_to_bytes, flat_to_coords3d};
use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::sync::Arc;

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

    /// Number of heavy atoms (non-hydrogen atoms, including isotopic H like [3H]).
    #[getter]
    fn heavy_atoms(&self) -> usize {
        chematic_chem::heavy_atom_count(&self.inner)
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

    /// Serialize this molecule to ChemAxon Marvin (MRV) XML.
    ///
    /// ``kekulize=True`` (RDKit's own default) writes alternating
    /// single/double bonds for aromatic rings instead of an aromatic bond
    /// order token; re-reading that output back does not restore the
    /// aromatic flag (a different, chemically-equivalent representation --
    /// re-perceive aromaticity if the flag itself must survive a round
    /// trip). Only single/double/triple/aromatic bonds are supported;
    /// any other bond order raises ``ValueError``.
    ///
    /// Equivalent to RDKit ``Chem.MolToMrvBlock(mol)``.
    ///
    ///     mrv = mol.to_mrv_block()
    #[pyo3(signature = (kekulize = true, include_stereo = true))]
    fn to_mrv_block(&self, kekulize: bool, include_stereo: bool) -> PyResult<String> {
        let record = chematic_mol::MoleculeRecord::new((*self.inner).clone());
        let options = chematic_mol::MrvWriteOptions {
            precision: 4,
            kekulize,
            include_stereo,
        };
        chematic_mol::write_mrv(&record, &options).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Serialize this molecule to MRV XML, preserving 2D/3D layout coordinates.
    ///
    /// Pass an empty list for either ``coords_2d``/``coords_3d`` to omit
    /// that dimensionality. Use with :func:`from_mrv_block_with_coords` for
    /// a coordinate-preserving round trip.
    ///
    ///     mol, coords_2d, coords_3d = chematic.from_mrv_block_with_coords(block)
    ///     new_block = mol.to_mrv_block_with_coords(coords_2d, coords_3d)
    #[pyo3(signature = (coords_2d, coords_3d, kekulize = true, include_stereo = true))]
    fn to_mrv_block_with_coords(
        &self,
        coords_2d: Vec<[f64; 2]>,
        coords_3d: Vec<[f64; 3]>,
        kekulize: bool,
        include_stereo: bool,
    ) -> PyResult<String> {
        let mut record = chematic_mol::MoleculeRecord::new((*self.inner).clone());
        if !coords_2d.is_empty() {
            record.coordinates_2d = Some(coords_2d);
        }
        if !coords_3d.is_empty() {
            record.coordinates_3d = Some(coords_3d);
        }
        let options = chematic_mol::MrvWriteOptions {
            precision: 4,
            kekulize,
            include_stereo,
        };
        chematic_mol::write_mrv(&record, &options).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Serialize this molecule to ChemDraw XML (CDXML).
    ///
    /// ``coords``: optional list of ``[x, y]`` pairs — one per heavy atom,
    /// in ChemDraw's Y-down convention. If omitted, atoms are written at
    /// ``(0, 0)``.
    ///
    /// Targets self-round-trip correctness (``chematic.from_cdxml`` can read
    /// what this writes), not full ChemDraw-application compatibility.
    ///
    ///     cdxml = mol.to_cdxml()
    ///     cdxml_with_layout = mol.to_cdxml(coords_2d)
    #[pyo3(signature = (coords = None))]
    fn to_cdxml(&self, coords: Option<Vec<[f64; 2]>>) -> String {
        let coords_2d: Vec<(f64, f64)> = coords
            .unwrap_or_default()
            .iter()
            .map(|xy| (xy[0], xy[1]))
            .collect();
        chematic_mol::write_cdxml(&self.inner, &coords_2d)
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
    ///     ``iterations`` (int), ``converged`` (bool), ``sound`` (bool —
    ///     all-finite coordinates and no bond stretched past a sane
    ///     covalent-bond length; independent of ``converged``, since
    ///     steepest descent often reports ``converged=False`` on geometries
    ///     that are perfectly fine but simply haven't hit the tight
    ///     RMS-gradient threshold yet. Check this, not just ``converged``,
    ///     before trusting a result).
    ///
    /// Example::
    ///
    ///     mol = chematic.from_smiles("CCO")
    ///     result = mol.minimize_uff([[0,0,0],[1.54,0,0],[2.5,1.2,0]])
    ///     print(result["energy"], result["converged"], result["sound"])
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
        d.set_item("sound", result.sound)?;
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

    /// Apply aromaticity perception and return a new Mol.
    ///
    /// Args:
    ///     mode: ``"huckel"`` (default) or ``"rdkit"``. The ``"rdkit"`` mode
    ///           additionally recognises Se and Te as lone-pair donors in aromatic rings.
    ///
    /// Returns:
    ///     A new :class:`Mol` with aromatic flags set.
    ///
    ///     mol2 = mol.apply_aromaticity(mode="rdkit")
    #[pyo3(signature = (mode = "huckel"))]
    fn apply_aromaticity(&self, mode: &str) -> PyResult<Mol> {
        let algo = match mode {
            "huckel" | "hückel" | "" => chematic_perception::AromaticityAlgorithm::Huckel,
            "rdkit" | "rdkit_like" | "rdkit-like" => {
                chematic_perception::AromaticityAlgorithm::RdkitLike
            }
            other => {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "Unknown aromaticity mode {other:?}. Choose 'huckel' or 'rdkit'."
                )));
            }
        };
        let new_mol = chematic_perception::apply_aromaticity_ex(&self.inner, algo);
        Ok(Mol {
            inner: Arc::new(new_mol),
            props: Default::default(),
        })
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
        d.set_item("num_saturated_heterocycles", rb.num_saturated_heterocycles)?;
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
        d.set_item(
            "num_valence_electrons",
            chematic_chem::num_valence_electrons(m),
        )?;
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
        let lipinski_violations: u8 = [
            (mw > 500.0) as u8,
            (hbd > 5) as u8,
            (hba_lipinski > 10) as u8,
            (logp > 5.0) as u8,
        ]
        .iter()
        .sum();

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
        if !pains_ok {
            alerts.push("PAINS alert");
        }
        if !brenk_ok {
            alerts.push("Brenk alert");
        }
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
        let veber_result = if veber_ok {
            "✓ Likely (TPSA ≤ 140, rot ≤ 10)"
        } else {
            "✗ Unlikely"
        };

        let lipophilicity = if logp < 0.0 {
            "hydrophilic"
        } else if logp < 2.0 {
            "mildly lipophilic"
        } else if logp < 4.0 {
            "moderately lipophilic"
        } else {
            "highly lipophilic"
        };

        // ADMET labels
        let bbb = if admet.bbb_passes {
            "✓ CNS penetrant"
        } else {
            "✗ CNS non-penetrant"
        };
        let caco2 = if admet.caco2 > -5.5 {
            "High (well absorbed)"
        } else if admet.caco2 > -7.0 {
            "Moderate"
        } else {
            "Low (poor absorption)"
        };
        let herg = if admet.herg_risk < 0.3 {
            "Low"
        } else if admet.herg_risk < 0.6 {
            "Moderate"
        } else {
            "High ⚠"
        };
        let cyp = if admet.cyp3a4_risk < 0.3 {
            "Low"
        } else if admet.cyp3a4_risk < 0.6 {
            "Moderate"
        } else {
            "High ⚠"
        };

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
            pains_s = if pains_ok {
                "✓ Clean"
            } else {
                "✗ Alert detected"
            },
            brenk_s = if brenk_ok {
                "✓ Clean"
            } else {
                "⚠ Alert detected"
            },
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
            let d =
                counts2.get(elem).copied().unwrap_or(0) - counts1.get(elem).copied().unwrap_or(0);
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
        let elem_parts: Vec<String> = delta_elements
            .iter()
            .map(|(e, d)| {
                if *d > 0 {
                    format!("+{d}{e}")
                } else {
                    format!("{d}{e}")
                }
            })
            .collect();
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

    /// ECFP fingerprint with bitInfo for RDKit-style ``bitInfo`` maps.
    ///
    /// Returns ``(bytes, {bit: [(atom_idx, radius), ...]})``. The fingerprint
    /// bits are identical to ``ecfp4()`` (radius 2) / ``ecfp6()`` (radius 3),
    /// so each recorded ``(atom_idx, radius)`` is the true origin of its bit.
    fn ecfp_bitinfo(&self, radius: u32) -> EcfpBitInfo {
        let cfg = chematic_fp::EcfpConfig {
            radius,
            ..Default::default()
        };
        let (fp, info) = chematic_fp::ecfp_with_bitinfo(&self.inner, &cfg);
        let bytes = bitvec2048_to_bytes(&fp);
        let dict = info
            .into_iter()
            .map(|(bit, v)| {
                (
                    bit,
                    v.into_iter()
                        .map(|(a, r)| (a as usize, r as usize))
                        .collect(),
                )
            })
            .collect();
        (bytes, dict)
    }

    /// FCFP4 (functional-class ECFP4) fingerprint as bytes (256 bytes = 2048 bits, LSB-first).
    fn fcfp4(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::fcfp4(&self.inner))
    }

    /// FCFP fingerprint with bitInfo for RDKit-style ``useFeatures=True, bitInfo`` maps.
    ///
    /// Returns ``(bytes, {bit: [(atom_idx, radius), ...]})``, mirroring
    /// [`Mol::ecfp_bitinfo`] but using pharmacophore feature-class atom
    /// invariants instead of plain atomic properties.
    fn fcfp_bitinfo(&self, radius: u32) -> EcfpBitInfo {
        let cfg = chematic_fp::EcfpConfig {
            radius,
            ..Default::default()
        };
        let (fp, info) = chematic_fp::fcfp_with_bitinfo(&self.inner, &cfg);
        let bytes = bitvec2048_to_bytes(&fp);
        let dict = info
            .into_iter()
            .map(|(bit, v)| {
                (
                    bit,
                    v.into_iter()
                        .map(|(a, r)| (a as usize, r as usize))
                        .collect(),
                )
            })
            .collect();
        (bytes, dict)
    }

    /// Atom-pair fingerprint as bytes (256 bytes = 2048 bits, LSB-first).
    fn atom_pair_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::atom_pair_fp(&self.inner))
    }

    /// Topological torsion fingerprint as bytes (256 bytes = 2048 bits, LSB-first).
    fn torsion_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::torsion_fp(&self.inner))
    }

    /// RDKit-*compatible* (not yet fully bit-exact) Topological Torsion
    /// fingerprint as bytes (256 bytes = 2048 bits).
    ///
    /// A from-scratch Rust port of RDKit's
    /// ``rdkit.Chem.AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=2048)``.
    /// Measured 87.2% bit-exact on a 1000-molecule general corpus sample against a live
    /// RDKit oracle. Two known, documented residual sources remain (see
    /// ``chematic_fp::rdkit_torsion``'s own module doc comment for detail): RDKit's
    /// hybridization-gated pi-electron count for hypervalent atoms (e.g. P/S) is not yet
    /// replicated, and asymmetrically-substituted 3-membered rings can need more closure
    /// entries than this implementation generates. A separate, opt-in function from
    /// :meth:`torsion_fp` (chematic's own native, similarity-preserving scheme, unchanged);
    /// neither affects the other's output. Does not support chirality (``includeChirality=True``
    /// is not implemented).
    fn rdkit_torsion_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::rdkit_torsion_fp(&self.inner))
    }

    /// RDKit-*compatible* (not yet fully bit-exact) Atom Pair fingerprint as
    /// bytes (256 bytes = 2048 bits).
    ///
    /// A from-scratch Rust port of RDKit's
    /// ``rdkit.Chem.AllChem.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=2048)``.
    /// Measured 87.2% bit-exact on a 1000-molecule general corpus sample against a live
    /// RDKit oracle -- every mismatching molecule in that sample contains a hypervalent
    /// S or P atom, confirming the *only* residual is the one already documented for
    /// :meth:`rdkit_torsion_fp` (RDKit's hybridization-gated pi-electron count for
    /// hypervalent atoms is not yet replicated; see ``chematic_fp::rdkit_torsion``'s
    /// module doc comment). A separate, opt-in function from :meth:`atom_pair_fp`
    /// (chematic's own native scheme, unchanged); neither affects the other's output.
    /// Does not support chirality (``includeChirality=True`` is not implemented).
    fn rdkit_atom_pair_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::rdkit_atom_pair_fp(&self.inner))
    }

    /// RDKit-compatible "Layered fingerprint" (``rdkit.Chem.LayeredFingerprint``)
    /// as bytes (256 bytes = 2048 bits).
    ///
    /// A from-scratch Rust port of RDKit's own branched-subgraph, 6-layer
    /// fingerprint (``layerFlags=0xFFFFFFFF, minPath=1, maxPath=7,
    /// fpSize=2048``); upstream itself documents this fingerprint as
    /// experimental. A separate, opt-in function from :meth:`layered_fp`
    /// (chematic's own pre-existing, non-bit-exact scheme); neither affects
    /// the other's output. See ``chematic_fp::rdkit_layered``'s module doc
    /// comment for the full algorithm and current bit-exactness figures.
    fn rdkit_layered_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::rdkit_layered_fp(&self.inner))
    }

    /// RDKit-*compatible* (not yet fully bit-exact) Pattern fingerprint as
    /// bytes (256 bytes = 2048 bits).
    ///
    /// A from-scratch Rust port of RDKit's
    /// ``rdkit.Chem.PatternFingerprint(mol, fpSize=2048)``. Unlike
    /// :meth:`rdkit_torsion_fp`/:meth:`rdkit_atom_pair_fp`, this uses SMARTS
    /// substructure matching against 13 fixed patterns rather than a path/pair
    /// enumeration. Measured against a live RDKit oracle on three corpora with
    /// different chemical distributions: 100% bit-exact on a general and a
    /// ChEMBL sample, 99.6% on an NCI sample -- the residual traces entirely to
    /// chematic's own aromaticity-perception model disagreeing with RDKit's on
    /// specific ring systems (see ``chematic_fp::rdkit_pattern``'s module doc
    /// comment for detail), not a defect in this port's own logic. A separate,
    /// opt-in function from :meth:`pattern_fp` (chematic's own native scheme,
    /// unchanged); neither affects the other's output. Does not support
    /// ``tautomericFingerprint=True``.
    fn rdkit_pattern_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::rdkit_pattern_fp(&self.inner))
    }

    /// RDKit-compatible "RDKit fingerprint" (``rdkit.Chem.RDKFingerprint``) as
    /// bytes (256 bytes = 2048 bits).
    ///
    /// A from-scratch Rust port of RDKit's default branched-subgraph
    /// fingerprint (``minPath=1, maxPath=7, fpSize=2048, numBitsPerFeature=2``),
    /// including RDKit's own deliberately weakened Mersenne Twister variant
    /// (with a boost-library quirk that silently changes its effective seeding
    /// constant from what RDKit's own source suggests) used for the second bit
    /// per feature. 100% bit-exact against a live RDKit oracle on
    /// ``descriptor_census_corpus.smi`` and ``chembl_accuracy_corpus_4999.smi``,
    /// 99.44% on ``nci_first_5k_smiles_only.smi`` (every residual is a known,
    /// pre-existing ``chematic-perception`` aromaticity-model gap on fused
    /// polyheteroaromatic dyes and metal-coordination complexes, not a defect
    /// in this port). See
    /// ``chematic_fp::rdkit_rdk``'s module doc comment for the full algorithm.
    /// A separate function from :meth:`path_fp` (chematic's own pre-existing,
    /// non-bit-exact linear-path approximation, itself Rust-side named
    /// ``rdkit_path_fp``); neither affects the other's output.
    fn rdkit_rdk_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::rdkit_rdk_fp(&self.inner))
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

    /// RDKit-bit-exact ECFP4 (radius=2, 2048 bits, ``useChirality=False``,
    /// ``useBondTypes=True``, RDKit's default atom invariant) as bytes (256 bytes =
    /// 2048 bits, LSB-first) -- bit-for-bit identical to
    /// ``rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol)``
    /// for every input this preprocessing handles. **Not** the same bits as
    /// :meth:`ecfp4` (that path uses chematic's own FNV-1a hash and is not RDKit-bit-
    /// compatible by design -- the two are never silently interchanged).
    ///
    /// Raises ``ValueError`` if RDKit-parity aromaticity preprocessing fails (a real,
    /// currently-unfixed structural class -- certain bridgehead-nitrogen fused
    /// heterocycles that cannot be kekulized under either engine). Never silently
    /// falls back to :meth:`ecfp4`'s Hückel-based engine on such input -- the two
    /// engines are not bit-compatible, so a silent substitution would look successful
    /// while actually returning the wrong hash. See ``docs/rfcs/ecfp4_bitexact_api_rfc.md``.
    fn rdkit_ecfp4(&self) -> PyResult<Vec<u8>> {
        let result = chematic_fp::rdkit_morgan_ecfp4_experimental(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(bitvec2048_to_bytes(&result.fingerprint))
    }

    /// Same fingerprint as :meth:`rdkit_ecfp4`, plus the raw (unfolded) data behind it.
    ///
    /// Returns ``(fingerprint, sparse_counts, raw_bit_info, folded_bit_info)``:
    ///
    /// - ``fingerprint``: ``bytes`` (256 bytes = 2048 bits, LSB-first) -- identical to
    ///   :meth:`rdkit_ecfp4`.
    /// - ``sparse_counts``: ``{raw_id: count}`` -- RDKit's ``GetSparseCountFingerprint``
    ///   shape (unfolded 32-bit identifiers).
    /// - ``raw_bit_info``: ``{raw_id: [(atom_idx, radius), ...]}`` -- RDKit's
    ///   ``AdditionalOutput.GetBitInfoMap()`` on the unfolded fingerprint.
    /// - ``folded_bit_info``: ``{bit: [(atom_idx, radius), ...]}`` -- the same, folded
    ///   to the 2048-bit fingerprint.
    ///
    /// Raises ``ValueError`` on the same preprocessing failures as :meth:`rdkit_ecfp4`.
    fn rdkit_ecfp4_detail(&self) -> PyResult<RdkitMorganDetail> {
        let result = chematic_fp::rdkit_morgan_ecfp4_experimental(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok((
            bitvec2048_to_bytes(&result.fingerprint),
            result.sparse_counts.into_iter().collect(),
            result.raw_bit_info.into_iter().collect(),
            result.folded_bit_info.into_iter().collect(),
        ))
    }

    /// RDKit-bit-exact Morgan/ECFP fingerprint at a caller-chosen radius/bit-width, as
    /// bytes (``nbits / 8`` bytes, LSB-first).
    ///
    /// ``radius`` must be one of 0, 1, 2 (:meth:`rdkit_ecfp4`'s ECFP4), or 3.
    /// ``nbits`` must be one of 128, 256, 512, 1024, or 2048. Each of these 20
    /// combinations is independently re-verified against a live RDKit oracle (not
    /// assumed to generalize from radius=2/2048 bits alone) -- see
    /// ``validation/ecfp4_rdkit_stable_api_fixtures.json``. Any other value raises
    /// ``ValueError`` rather than being silently coerced to the nearest supported one.
    ///
    /// Raises ``ValueError`` on the same preprocessing failures as :meth:`rdkit_ecfp4`
    /// (regardless of ``radius``/``nbits`` -- the failure happens before folding).
    /// ``include_chirality`` enables RDKit-compatible tetrahedral chirality.
    /// E/Z bond stereo is not included by this API yet.
    #[pyo3(signature = (radius = 2, nbits = 2048, include_chirality = false))]
    fn rdkit_ecfp_config(
        &self,
        radius: u32,
        nbits: usize,
        include_chirality: bool,
    ) -> PyResult<Vec<u8>> {
        let config = python_rdkit_morgan_config(radius, nbits, include_chirality)?;
        let result = chematic_fp::rdkit_morgan_fingerprint(&self.inner, &config)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(bitvecn_to_bytes(&result.fingerprint))
    }

    /// Same fingerprint as :meth:`rdkit_ecfp_config`, plus the raw (unfolded) data --
    /// see :meth:`rdkit_ecfp4_detail` for the return shape (identical, generalized to
    /// this method's ``radius``/``nbits``).
    #[pyo3(signature = (radius = 2, nbits = 2048, include_chirality = false))]
    fn rdkit_ecfp_config_detail(
        &self,
        radius: u32,
        nbits: usize,
        include_chirality: bool,
    ) -> PyResult<RdkitMorganDetail> {
        let config = python_rdkit_morgan_config(radius, nbits, include_chirality)?;
        let result = chematic_fp::rdkit_morgan_fingerprint(&self.inner, &config)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok((
            bitvecn_to_bytes(&result.fingerprint),
            result.sparse_counts.into_iter().collect(),
            result.raw_bit_info.into_iter().collect(),
            result.folded_bit_info.into_iter().collect(),
        ))
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

    /// Structured 2D depiction data (atoms + bonds with layout coordinates).
    ///
    /// Use this instead of ``svg()`` when you want to drive your own
    /// renderer (e.g. matplotlib, a custom canvas) rather than parse SVG.
    ///
    /// Returns:
    ///     dict with keys:
    ///     ``atoms`` (list of dicts: ``idx``, ``element`` (symbol string),
    ///     ``x``, ``y``, ``label`` (``None`` when suppressed), ``color``
    ///     (CSS hex string), ``charge``) and ``bonds`` (list of dicts:
    ///     ``idx``, ``atom1``, ``atom2``, ``kind`` — one of ``"Single"``,
    ///     ``"Double"``, ``"Triple"``, ``"Aromatic"``, ``"Up"``, ``"Down"``).
    ///
    /// Example::
    ///
    ///     mol = chematic.from_smiles("CCO")
    ///     data = mol.depict_data()
    ///     for atom in data["atoms"]:
    ///         print(atom["element"], atom["x"], atom["y"])
    fn depict_data<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let data = chematic_depict::compute_depict_data(&self.inner);

        let atoms = data
            .atoms
            .iter()
            .map(|a| {
                let ad = PyDict::new(py);
                ad.set_item("idx", a.idx.0)?;
                ad.set_item("element", a.element.symbol())?;
                ad.set_item("x", a.pos.x)?;
                ad.set_item("y", a.pos.y)?;
                ad.set_item("label", a.label.as_deref())?;
                ad.set_item("color", &a.color)?;
                ad.set_item("charge", a.charge)?;
                Ok::<_, PyErr>(ad)
            })
            .collect::<PyResult<Vec<_>>>()?;

        let bond_kind = |k: &chematic_depict::DepictBondKind| match k {
            chematic_depict::DepictBondKind::Single => "Single",
            chematic_depict::DepictBondKind::Double => "Double",
            chematic_depict::DepictBondKind::Triple => "Triple",
            chematic_depict::DepictBondKind::Aromatic => "Aromatic",
            chematic_depict::DepictBondKind::Up => "Up",
            chematic_depict::DepictBondKind::Down => "Down",
        };

        let bonds = data
            .bonds
            .iter()
            .map(|b| {
                let bd = PyDict::new(py);
                bd.set_item("idx", b.idx.0)?;
                bd.set_item("atom1", b.atom1.0)?;
                bd.set_item("atom2", b.atom2.0)?;
                bd.set_item("kind", bond_kind(&b.kind))?;
                Ok::<_, PyErr>(bd)
            })
            .collect::<PyResult<Vec<_>>>()?;

        let d = PyDict::new(py);
        d.set_item("atoms", atoms)?;
        d.set_item("bonds", bonds)?;
        Ok(d)
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

    /// Serialize the molecule to a MolJSON string (pretty-printed).
    ///
    /// MolJSON is a JSON-based molecular representation designed for LLM
    /// (large language model) compatibility.  Unlike SMILES, it makes atoms,
    /// bonds, and connectivity explicit without requiring domain-specific
    /// parsing rules.
    ///
    ///     json_str = mol.to_moljson()
    ///     mol2 = chematic.from_moljson(json_str)
    fn to_moljson(&self) -> String {
        chematic_mol::write_moljson(&self.inner)
    }

    /// Return the molecule in a text format suited for LLM prompts.
    ///
    /// Based on arXiv 2026 "Rethinking Molecular Text Representations for LLMs":
    /// CML and MolJSON outperform SMILES on structural reasoning tasks.
    ///
    /// Args:
    ///     format: one of ``"canonical_smiles"``, ``"smiles"``, ``"inchi"``,
    ///             ``"inchikey"``, ``"moljson"``, ``"cml"``, ``"markdown"``
    ///
    ///     # Task-aware: use chematic.best_representation(task) to pick format
    ///     fmt = chematic.best_representation("structural_reasoning")  # → "moljson"
    ///     text = mol.to_llm_text(fmt)
    ///
    ///     # Direct format selection:
    ///     mol.to_llm_text("moljson")
    ///     mol.to_llm_text("inchi")
    ///     mol.to_llm_text("markdown")   # multi-field summary
    fn to_llm_text(&self, format: &str) -> PyResult<String> {
        match format {
            "canonical_smiles" | "smiles" => Ok(chematic_smiles::canonical_smiles(&self.inner)),
            "inchi" => Ok(chematic_inchi::inchi(&self.inner)),
            "inchikey" => {
                let i = chematic_inchi::inchi(&self.inner);
                Ok(chematic_inchi::inchi_key(&i))
            }
            "moljson" => Ok(chematic_mol::write_moljson(&self.inner)),
            "cml" => Ok(chematic_mol::write_cml(&self.inner, None)),
            "markdown" => {
                let smi = chematic_smiles::canonical_smiles(&self.inner);
                let inchi = chematic_inchi::inchi(&self.inner);
                let mw = chematic_chem::molecular_weight(&self.inner);
                let hac = chematic_chem::heavy_atom_count(&self.inner);
                Ok(format!(
                    "**SMILES**: {smi}\n**InChI**: {inchi}\n**MW**: {mw:.2} Da\n**Heavy atoms**: {hac}\n**MolJSON**:\n{}",
                    chematic_mol::write_moljson(&self.inner)
                ))
            }
            _ => Err(PyValueError::new_err(format!(
                "Unknown format '{format}'. Choose: canonical_smiles, inchi, inchikey, moljson, cml, markdown"
            ))),
        }
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
        // Stop at the first embedding instead of enumerating every match — an
        // existence check doesn't need the full match set or the dedup pass.
        let config = chematic_smarts::MatchConfig {
            max_matches: Some(1),
            uniquify: false,
            ..chematic_smarts::MatchConfig::default()
        };
        Ok(!chematic_smarts::find_matches_with_config(&query, &self.inner, &config).is_empty())
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
            props: Default::default(),
        }
    }

    /// Return the fragment parent, retaining the largest chemically relevant
    /// fragment and its transformation semantics.
    fn fragment_parent(&self) -> Mol {
        Mol::bare(chematic_chem::fragment_parent(&self.inner).0)
    }

    /// Return the charge parent (fragment selection followed by charge
    /// normalization), without removing isotopes or stereochemistry.
    fn charge_parent(&self) -> Mol {
        Mol::bare(chematic_chem::charge_parent(&self.inner).0)
    }

    /// Return the isotope parent while preserving stereochemistry.
    fn isotope_parent(&self) -> Mol {
        Mol::bare(chematic_chem::isotope_parent(&self.inner).0)
    }

    /// Return the stereo parent with stereochemical annotations removed.
    fn stereo_parent(&self) -> Mol {
        Mol::bare(chematic_chem::stereo_parent(&self.inner).0)
    }

    /// Return the Murcko scaffold as a new Mol.
    fn scaffold(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::murcko_scaffold(&self.inner)),
            props: Default::default(),
        }
    }

    /// Return the canonical tautomer as a new Mol.
    fn canonical_tautomer(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::canonical_tautomer(&self.inner)),
            props: Default::default(),
        }
    }

    /// Return the tautomer parent and its computation status.
    #[pyo3(signature = (max_transforms=16, max_tautomers=32, timeout_ms=None))]
    fn tautomer_parent(
        &self,
        max_transforms: usize,
        max_tautomers: usize,
        timeout_ms: Option<u64>,
    ) -> (Mol, String) {
        let mut limits = chematic_chem::TautomerLimits::default();
        limits.max_transforms = max_transforms;
        limits.max_tautomers = max_tautomers;
        limits.timeout_ms = timeout_ms;
        let result = chematic_chem::tautomer_parent(&self.inner, &limits);
        (Mol::bare(result.molecule), format!("{:?}", result.status))
    }

    /// Return the composed super parent and its computation status.
    #[pyo3(signature = (max_transforms=16, max_tautomers=32, timeout_ms=None))]
    fn super_parent(
        &self,
        max_transforms: usize,
        max_tautomers: usize,
        timeout_ms: Option<u64>,
    ) -> (Mol, String) {
        let mut limits = chematic_chem::TautomerLimits::default();
        limits.max_transforms = max_transforms;
        limits.max_tautomers = max_tautomers;
        limits.timeout_ms = timeout_ms;
        limits.timeout_ms = timeout_ms;
        let result = chematic_chem::super_parent(&self.inner, &limits);
        (Mol::bare(result.molecule), format!("{:?}", result.status))
    }

    /// Return the composed Parent result with every intermediate stage.
    ///
    /// The returned dictionary contains ``smiles``, ``status``, and a
    /// ``stages`` list with the five ordered Parent transformations. This is
    /// the binding-level counterpart of Rust's ``ParentAudit::Composed`` and
    /// keeps provenance inspectable without exposing internal Rust enums.
    #[pyo3(signature = (max_transforms=16, max_tautomers=32, timeout_ms=None))]
    fn super_parent_report<'py>(
        &self,
        py: Python<'py>,
        max_transforms: usize,
        max_tautomers: usize,
        timeout_ms: Option<u64>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let mut limits = chematic_chem::TautomerLimits::default();
        limits.max_transforms = max_transforms;
        limits.max_tautomers = max_tautomers;
        limits.timeout_ms = timeout_ms;

        let (fragment, _) = chematic_chem::fragment_parent(&self.inner);
        let (charge, _) = chematic_chem::charge_parent(&fragment);
        let (isotope, _) = chematic_chem::isotope_parent(&charge);
        let (stereo, _) = chematic_chem::stereo_parent(&isotope);
        let result = chematic_chem::super_parent(&self.inner, &limits);
        let stages = PyList::empty(py);
        for (name, molecule) in [
            ("fragment", &fragment),
            ("charge", &charge),
            ("isotope", &isotope),
            ("stereo", &stereo),
            ("tautomer", &result.molecule),
        ] {
            let stage = PyDict::new(py);
            stage.set_item("name", name)?;
            stage.set_item("smiles", chematic_smiles::canonical_smiles(molecule))?;
            stages.append(stage)?;
        }
        let report = PyDict::new(py);
        report.set_item(
            "smiles",
            chematic_smiles::canonical_smiles(&result.molecule),
        )?;
        report.set_item("status", format!("{:?}", result.status))?;
        report.set_item("stages", stages)?;
        Ok(report)
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
            .map(Mol::bare)
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
            .map(Mol::bare)
            .collect()
    }

    /// Return a copy with all implicit hydrogens made explicit.
    fn add_hydrogens(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::add_hydrogens(&self.inner)),
            props: Default::default(),
        }
    }

    /// Return a copy with all explicit hydrogen atoms removed.
    fn remove_hydrogens(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::remove_hydrogens(&self.inner)),
            props: Default::default(),
        }
    }

    /// Return a copy with all stereochemistry assignments removed.
    fn remove_stereo(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::remove_stereo(&self.inner)),
            props: Default::default(),
        }
    }

    /// Return a copy with all isotope labels removed.
    fn remove_isotopes(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::remove_isotopes(&self.inner)),
            props: Default::default(),
        }
    }

    /// Return the largest covalently connected fragment.
    fn largest_fragment(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::largest_fragment(&self.inner)),
            props: Default::default(),
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
        self.inner.fragments().into_iter().map(Mol::bare).collect()
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
            props: Default::default(),
        }
    }

    /// Return the generic Murcko scaffold (all atoms replaced with carbons, all bonds single).
    fn generic_scaffold(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::generic_murcko_scaffold(&self.inner)),
            props: Default::default(),
        }
    }

    /// Fragment the molecule using BRICS rules. Returns a list of fragment Mol objects.
    ///
    /// When no BRICS-breakable bonds are found, returns a list containing the
    /// original molecule (not an empty list).
    fn brics_fragments(&self) -> Vec<Mol> {
        chematic_chem::brics_fragments(&self.inner)
            .into_iter()
            .map(Mol::bare)
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

    /// Hyper-Dimensional Fingerprint (HDF) as a numpy float32 array of shape ``(dim,)``.
    ///
    /// HDF encodes a molecule as a unit-norm float32 vector using Hyperdimensional
    /// Computing (HDC).  Unlike hash-based fingerprints (ECFP), there is no hash
    /// collision: every distinct atom environment maps to a unique HD vector.
    /// Cosine similarity between two normalized HDF vectors gives molecular similarity.
    ///
    /// Based on: "Hyper-Dimensional Fingerprints as Molecular Representations" (arXiv 2026).
    ///
    ///     fp = mol.hdf()              # shape (1024,), dtype float32, unit norm
    ///     fp = mol.hdf(dim=512)       # smaller vector
    ///     fp = mol.hdf(dim=2048, radius=3)
    ///
    ///     sim = float(np.dot(mol1.hdf(), mol2.hdf()))  # cosine similarity
    #[pyo3(signature = (dim = 1024, radius = 2, seed = 42))]
    fn hdf<'py>(
        &self,
        py: Python<'py>,
        dim: usize,
        radius: usize,
        seed: u64,
    ) -> Bound<'py, PyArray1<f32>> {
        let config = chematic_fp::HdfConfig { dim, radius, seed };
        let fp = chematic_fp::hdf(&self.inner, &config);
        Array1::from_vec(fp.0).into_pyarray(py)
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

    /// Per-atom structural data for RDKit-compat wrappers.
    ///
    /// Returns one tuple per heavy atom:
    /// ``(symbol, atomic_num, formal_charge, is_aromatic, implicit_h, heavy_degree, is_in_ring)``
    #[getter]
    fn atom_table(&self) -> Vec<(String, u8, i8, bool, u8, usize, bool)> {
        use chematic_core::AtomIdx;
        use chematic_perception::find_sssr;
        use std::collections::HashSet;
        let mol = &self.inner;
        let ring_atoms: HashSet<AtomIdx> = find_sssr(mol)
            .rings()
            .iter()
            .flat_map(|r| r.iter().copied())
            .collect();
        (0..mol.atom_count())
            .map(|i| {
                let idx = AtomIdx(i as u32);
                let atom = mol.atom(idx);
                (
                    atom.element.symbol().to_string(),
                    atom.element.atomic_number(),
                    atom.charge,
                    atom.aromatic,
                    chematic_core::implicit_hcount(mol, idx),
                    mol.degree(idx),
                    ring_atoms.contains(&idx),
                )
            })
            .collect()
    }

    /// Per-bond structural data for RDKit-compat wrappers.
    ///
    /// Returns one tuple per bond:
    /// ``(atom1_idx, atom2_idx, bond_type_str, is_aromatic)``
    #[getter]
    fn bond_table(&self) -> Vec<(usize, usize, &'static str, bool)> {
        use chematic_core::BondOrder;
        let mol = &self.inner;
        mol.bonds()
            .map(|(_, b)| {
                let (type_str, is_aromatic) = match b.order {
                    BondOrder::Single | BondOrder::Up | BondOrder::Down => ("SINGLE", false),
                    BondOrder::Double => ("DOUBLE", false),
                    BondOrder::Triple => ("TRIPLE", false),
                    BondOrder::Aromatic => ("AROMATIC", true),
                    _ => ("OTHER", false),
                };
                (
                    b.atom1.0 as usize,
                    b.atom2.0 as usize,
                    type_str,
                    is_aromatic,
                )
            })
            .collect()
    }

    /// SSSR rings as lists of atom indices. Used by the Python RingInfo wrapper.
    #[getter]
    fn sssr_atom_rings(&self) -> Vec<Vec<usize>> {
        use chematic_perception::find_sssr;
        find_sssr(&self.inner)
            .rings()
            .iter()
            .map(|r| r.iter().map(|idx| idx.0 as usize).collect())
            .collect()
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

    /// Generate a conformer ensemble using ETKDG + force-field minimization + RMSD pruning.
    ///
    /// .. deprecated::
    ///     Prefer :meth:`conformer_ensemble_v2`, which is deterministic
    ///     (explicit seed), energy-ranks kept conformers, reports full
    ///     per-attempt provenance (kept / duplicate-pruned / failed), and
    ///     avoids a known soundness defect in this method's underlying
    ///     MMFF94 path (silently zero energy/gradient for atom-type pairs
    ///     its tables don't cover -- see PR #369). Also, on any internal
    ///     failure this method silently returns an empty list rather than
    ///     raising. This method is unchanged and not scheduled for removal
    ///     in this release.
    ///
    /// Returns a list of coordinate arrays — each is a ``[[x,y,z], ...]`` list.
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
        // Generate directly on self.inner's own atom order -- do NOT re-parse
        // from canonical SMILES here. canonical_smiles() routinely reorders
        // atoms (any branch or ring), and generating on that reparsed molecule
        // while returning coordinates "as-is" desyncs them from the Mol the
        // caller already holds (atom_table, cip_stereo(), bond_table, ...
        // all stay indexed by self.inner's original order). Matches the
        // existing, correct pattern in generate_3d()/generate_3d_etkdg() below.
        // See issue #172.
        let mol = (*self.inner).clone();
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

    /// Run the opt-in v2 embedding pipeline (torsion-knowledge-aware distance
    /// geometry + stereo verification/repair + policy-gated force field).
    ///
    /// Applied directly to this ``Mol``'s own atom order -- never
    /// canonicalizes/reparses first (see issue #172; the fix in
    /// :meth:`conformer_ensemble` is the precedent this follows), so the
    /// returned ``coords`` and every atom/bond index in the result dict
    /// correspond 1:1 to this ``Mol``'s existing atom/bond tables.
    ///
    /// ``config``: a :class:`PipelineV2Config`.
    ///
    /// Returns a dict with keys: ``coords``, ``embed_stats``,
    /// ``bound_adjustment_report``, ``torsion_knowledge_report``,
    /// ``ring_torsion_evidence``, ``torsion_optimization_report``,
    /// ``stereo_before``, ``stereo_repair``, ``stereo_after_repair``,
    /// ``force_field``, ``final_stereo``, ``final_validation``,
    /// ``elapsed_ms_by_stage``.
    ///
    /// Raises :class:`PipelineV2Error` (a ``ValueError`` subclass) on
    /// failure. ``error.diagnostics`` carries the same per-stage partial
    /// evidence a Rust caller sees on ``PipelineV2Failure`` --
    /// ``diagnostics['last_known_coords']`` is diagnostic only, never a
    /// usable result (see ``diagnostics['coords_are_diagnostic_only']``).
    fn embed_pipeline_v2<'py>(
        &self,
        py: Python<'py>,
        config: &crate::pipeline_v2::PyPipelineV2Config,
    ) -> PyResult<Bound<'py, PyDict>> {
        crate::pipeline_v2::run_embed_pipeline_v2(py, &self.inner, config)
    }

    /// Generate a multi-conformer ensemble by calling the v2 embedding
    /// pipeline (:meth:`embed_pipeline_v2`) ``config.count`` times, once
    /// per deterministically derived seed, then selecting kept
    /// representatives by ascending energy within each force-field group.
    ///
    /// Unlike :meth:`embed_pipeline_v2`, a call here does **not** raise
    /// just because no conformer was kept -- every per-attempt outcome,
    /// including an ensemble where every attempt failed, is a normal,
    /// fully-diagnosable result. This only raises ``ValueError`` for a
    /// ``config`` that could never succeed regardless of the molecule
    /// (currently: an invalid ``rmsd_threshold``). Always check
    /// ``len(result["conformers"])`` rather than relying on "no exception."
    ///
    /// Applied directly to this ``Mol``'s own atom order -- never
    /// canonicalizes/reparses first (see issue #172).
    ///
    /// ``config``: an :class:`EnsembleV2Config`.
    ///
    /// Returns a dict with keys:
    ///     conformers: kept conformers only, as ``[[[x,y,z], ...], ...]``,
    ///         ordered group-by-group (by force field actually used),
    ///         ascending energy within each group. Never a single
    ///         cross-group energy sort -- MMFF94 and UFF energies are not
    ///         on a comparable scale (see ``mixed_force_field`` below).
    ///     conformer_provenance: one dict per entry of ``conformers``, same
    ///         order: ``attempt_index``, ``seed``, ``energy``,
    ///         ``actual_force_field_used``.
    ///     attempts: every attempt, success or failure, in order. Each has
    ///         ``attempt_index``, ``seed``, ``outcome`` (``"success"`` or
    ///         ``"failure"``), and one of ``success``/``failure`` populated.
    ///         A successful attempt's dict has ``energy``,
    ///         ``actual_force_field_used``, ``fallback_reason``, and
    ///         ``disposition`` (``{"kind": "kept", "conformer_index"}`` or
    ///         ``{"kind": "pruned_as_duplicate", "representative_attempt_index",
    ///         "rmsd", "symmetric"}``). A failed attempt's dict has the same
    ///         shape :meth:`embed_pipeline_v2` raises on failure.
    ///     mixed_force_field: ``True`` iff kept conformers span more than
    ///         one force field actually used.
    ///     termination: ``"completed"`` or ``"timed_out"``
    ///         (``ensemble_timeout_ms`` exhausted before all attempts ran).
    ///     requested_count: ``config.count`` at call time.
    fn conformer_ensemble_v2<'py>(
        &self,
        py: Python<'py>,
        config: &crate::ensemble_v2::PyEnsembleV2Config,
    ) -> PyResult<Bound<'py, PyDict>> {
        crate::ensemble_v2::run_embed_ensemble_v2(py, &self.inner, config)
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
    /// ``descriptor`` is ``"R"``, ``"S"``, ``"E"``, ``"Z"``, ``"r"``, or ``"s"``.
    /// Only assigned stereocenters / double bonds are returned.
    ///
    /// ``mode`` selects the CIP engine:
    ///
    /// - ``"legacy"`` (default) — the original fast engine, ~96.3% agreement with
    ///   RDKit's modern ``rdCIPLabeler`` oracle. Unchanged from every prior release.
    /// - ``"accurate"`` — a hierarchical-digraph engine for tetrahedral R/S
    ///   (~99.6% oracle-agreement on the representation-stable subset; see
    ///   ``docs/rfcs/cip_accurate_rfc.md``), merged with legacy's E/Z and allene answers
    ///   (the accurate engine doesn't compute either). Atoms it explicitly can't
    ///   resolve (a genuine tie, or exceeding its computation budget) are omitted
    ///   here and reported instead via :meth:`cip_stereo_unresolved` — never a
    ///   silently-guessed label.
    #[pyo3(signature = (mode = "legacy"))]
    fn cip_stereo<'py>(&self, py: Python<'py>, mode: &str) -> PyResult<Vec<Bound<'py, PyDict>>> {
        use chematic_core::CipCode;
        let assignments = match mode {
            "legacy" => chematic_chem::assign_cip(&self.inner).assignments,
            "accurate" => {
                chematic_chem::assign_cip_with_mode(&self.inner, chematic_chem::CipMode::Accurate)
                    .map_err(|e| PyValueError::new_err(e.to_string()))?
                    .assignments
            }
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown CIP mode '{other}' -- expected 'legacy' or 'accurate'"
                )));
            }
        };
        assignments
            .iter()
            .map(|(idx, code)| {
                let d = PyDict::new(py);
                d.set_item("atom_idx", idx.0 as usize)?;
                let label = match code {
                    CipCode::R => "R",
                    CipCode::S => "S",
                    CipCode::E => "E",
                    CipCode::Z => "Z",
                    CipCode::LowerR => "r",
                    CipCode::LowerS => "s",
                };
                d.set_item("descriptor", label)?;
                Ok(d)
            })
            .collect()
    }

    /// Atoms :meth:`cip_stereo`\ (``mode="accurate"``) could not resolve a tetrahedral
    /// R/S for — list of ``{"atom_idx": int, "reason": str}`` dicts, ``reason`` is
    /// ``"tied"`` (a genuine CIP-rule tie, not a missing rule) or ``"budget_exceeded"``.
    /// Always empty for ``mode="legacy"`` (that engine never reports "I don't know").
    fn cip_stereo_unresolved<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let result =
            chematic_chem::assign_cip_with_mode(&self.inner, chematic_chem::CipMode::Accurate)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        result
            .unresolved
            .iter()
            .map(|(idx, reason)| {
                let d = PyDict::new(py);
                d.set_item("atom_idx", idx.0 as usize)?;
                let label = match reason {
                    chematic_chem::CipUnresolvedReason::Tied => "tied",
                    chematic_chem::CipUnresolvedReason::BudgetExceeded => "budget_exceeded",
                };
                d.set_item("reason", label)?;
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

    /// Avalon-style structural fingerprint as bytes (256 bytes = 2048 bits).
    ///
    /// A broad mix of atom, bond, ring, and path features, loosely modelled
    /// on RDKit's Avalon fingerprint (``rdkit.Avalon.pyAvalonTools.GetAvalonFP``).
    /// Bit positions are not RDKit-identical (see :mod:`chematic.rdkit_compat`
    /// notes on Morgan fingerprints for the same caveat).
    fn avalon_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::avalon_fp(&self.inner))
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
            props: Default::default(),
        }
    }

    /// Keep the largest organic fragment; remove inorganic counterions.
    ///
    /// Useful after salt removal when the largest fragment is the drug molecule.
    fn prefer_organic(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::prefer_organic(&self.inner)),
            props: Default::default(),
        }
    }

    /// Re-apply ionization rules based on pKa.
    ///
    /// Transfers protons to maximize negative charge on the strongest acids.
    fn reionize(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::reionize(&self.inner)),
            props: Default::default(),
        }
    }

    /// Remove all formal charges (set every atom to neutral).
    fn uncharge(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::uncharge(&self.inner)),
            props: Default::default(),
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
            props: Default::default(),
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
            props: Default::default(),
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
            props: Default::default(),
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
                    CipCode::LowerR => "r",
                    CipCode::LowerS => "s",
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
                    CipCode::LowerR => "r",
                    CipCode::LowerS => "s",
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
    // SD properties (RDKit-compatible)
    // -----------------------------------------------------------------------

    /// Get an SD property by name. Raises ``KeyError`` if not present.
    #[pyo3(name = "GetProp")]
    fn get_prop(&self, key: &str) -> PyResult<String> {
        self.props
            .get(key)
            .cloned()
            .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(format!("'{key}' not found")))
    }

    /// Set an SD property.
    #[pyo3(name = "SetProp")]
    fn set_prop(&mut self, key: String, value: String) {
        self.props.insert(key, value);
    }

    /// Return ``True`` if the property exists.
    #[pyo3(name = "HasProp")]
    fn has_prop(&self, key: &str) -> bool {
        self.props.contains_key(key)
    }

    /// Return all SD properties as a Python dict.
    #[pyo3(name = "GetPropsAsDict")]
    fn get_props_as_dict<'py>(&self, py: Python<'py>) -> Bound<'py, PyDict> {
        let d = PyDict::new(py);
        for (k, v) in &self.props {
            d.set_item(k, v).ok();
        }
        d
    }

    /// Return a list of all property names.
    #[pyo3(name = "GetPropNames")]
    fn get_prop_names(&self) -> Vec<String> {
        self.props.keys().cloned().collect()
    }

    /// Remove a property by name (no-op if not present).
    #[pyo3(name = "ClearProp")]
    fn clear_prop(&mut self, key: &str) {
        self.props.remove(key);
    }

    /// Set an integer property (stored as its string representation).
    #[pyo3(name = "SetIntProp")]
    fn set_int_prop(&mut self, key: String, val: i64) {
        self.props.insert(key, val.to_string());
    }

    /// Set a float property (stored as its string representation).
    #[pyo3(name = "SetDoubleProp")]
    fn set_double_prop(&mut self, key: String, val: f64) {
        self.props.insert(key, val.to_string());
    }

    /// Set a boolean property (stored as ``"1"`` / ``"0"``).
    #[pyo3(name = "SetBoolProp")]
    fn set_bool_prop(&mut self, key: String, val: bool) {
        self.props
            .insert(key, if val { "1" } else { "0" }.to_string());
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

/// Maps the plain `radius`/`nbits` integers `rdkit_ecfp_config`/`rdkit_ecfp_config_detail`
/// accept from Python into chematic-fp's closed `RdkitMorganConfig` enums. An
/// unsupported value raises `ValueError` explicitly rather than being coerced to the
/// nearest supported one -- there is no guessed conversion for an option this API
/// doesn't (yet) support.
fn python_rdkit_morgan_config(
    radius: u32,
    nbits: usize,
    include_chirality: bool,
) -> PyResult<chematic_fp::RdkitMorganConfig> {
    let radius = match radius {
        0 => chematic_fp::RdkitMorganRadius::R0,
        1 => chematic_fp::RdkitMorganRadius::R1,
        2 => chematic_fp::RdkitMorganRadius::R2,
        3 => chematic_fp::RdkitMorganRadius::R3,
        other => {
            return Err(PyValueError::new_err(format!(
                "unsupported radius {other} -- rdkit_ecfp_config only supports 0, 1, 2, or 3 \
                 (each independently verified against a live RDKit oracle)"
            )));
        }
    };
    let fp_size = match nbits {
        128 => chematic_fp::RdkitMorganFpSize::B128,
        256 => chematic_fp::RdkitMorganFpSize::B256,
        512 => chematic_fp::RdkitMorganFpSize::B512,
        1024 => chematic_fp::RdkitMorganFpSize::B1024,
        2048 => chematic_fp::RdkitMorganFpSize::B2048,
        other => {
            return Err(PyValueError::new_err(format!(
                "unsupported nbits {other} -- rdkit_ecfp_config only supports 128, 256, 512, \
                 1024, or 2048"
            )));
        }
    };
    Ok(chematic_fp::RdkitMorganConfig {
        radius,
        fp_size,
        include_chirality,
    })
}

/// Bit-pack a variable-width `BitVecN` into `bit_width()/8` bytes, LSB-first --
/// generalizes `formats::bitvec2048_to_bytes` to `rdkit_ecfp_config`'s caller-chosen
/// widths (always a multiple of 8: 128/256/512/1024/2048).
fn bitvecn_to_bytes(fp: &chematic_fp::bitvec::BitVecN) -> Vec<u8> {
    let byte_count = fp.bit_width() / 8;
    (0..byte_count)
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
