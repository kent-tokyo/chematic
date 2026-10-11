use super::*;

#[pymethods]
impl Mol {
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
    /// RDKit's ``rdMolDescriptors.CalcMQNs``: atom, polarity, charge,
    /// topology, bond and ring-size counts (Nguyen et al. 2009 order),
    /// returned as a list of ints.
    fn mqn(&self) -> Vec<u32> {
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
    /// RDKit's model (lone pairs count as orbitals): amide/enamine N and aryl
    /// ether O are sp2; s, sp3d and sp3d2 report 0. Useful for scaffold
    /// modification (PromptSMILES-style) and atom featurization.
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

    /// Symmetrized SSSR rings as atom-index lists for parity diagnostics.
    fn symmetrized_sssr_atom_rings(&self) -> Vec<Vec<usize>> {
        chematic_perception::find_symmetrized_sssr(&self.inner)
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
    /// ``ignore_interfrag_interactions=True`` leaves out van der Waals and
    /// electrostatic pairs between disconnected fragments, as RDKit's
    /// ``MMFFGetMoleculeForceField`` does by default; chematic keeps them by
    /// default.
    /// Raises ``ValueError`` for atoms not parameterised by MMFF94.
    #[pyo3(signature = (coords, ignore_interfrag_interactions = false))]
    fn mmff94_energy_breakdown<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let b = chematic_ff::mmff94_energy_breakdown_with_options(
            &self.inner,
            &coords,
            chematic_ff::Mmff94Options {
                ignore_interfragment_interactions: ignore_interfrag_interactions,
            },
        )
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

    /// Experimental analytic MMFF94 gradient for given 3D coordinates.
    ///
    /// This is the gradient of the same prepared objective used by
    /// :meth:`mmff94_energy_breakdown`. It is exposed for numerical validation;
    /// the default production minimizer continues to use finite differences
    /// until the full A6 gradient and neighbor-list gates are complete.
    fn mmff94_bounded_analytic_gradient(&self, coords: Vec<[f64; 3]>) -> PyResult<Vec<[f64; 3]>> {
        if coords.len() != self.inner.atom_count() {
            return Err(PyValueError::new_err(format!(
                "coordinate count mismatch: expected {}, got {}",
                self.inner.atom_count(),
                coords.len()
            )));
        }
        let model = chematic_ff::Mmff94EnergyModel::new(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(model.bounded_analytic_gradient(&coords))
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
}
