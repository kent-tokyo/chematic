use super::*;

#[pymethods]
impl Mol {
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

    /// Checked legacy MMFF94 BCI charges.
    ///
    /// Unlike :meth:`mmff94_charges`, this raises ``ValueError`` when the
    /// legacy element-pair table cannot represent an atom.
    fn mmff94_charges_checked(&self) -> PyResult<Vec<f64>> {
        chematic_chem::try_mmff94_charges_bci(&self.inner)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Balaban J topological complexity index.
    ///
    /// High J → more complex/branched; related to molecular uniqueness.
    #[getter]
    fn balaban_j(&self) -> f64 {
        chematic_chem::balaban_j(&self.inner)
    }

    /// Checked Balaban J calculation.
    fn balaban_j_checked(&self) -> PyResult<f64> {
        chematic_chem::try_balaban_j(&self.inner)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Information-theoretic connectivity index (IPC).
    #[getter]
    fn ipc(&self) -> f64 {
        chematic_chem::ipc(&self.inner)
    }

    /// Checked IPC calculation.
    fn ipc_checked(&self) -> PyResult<f64> {
        chematic_chem::try_ipc(&self.inner)
            .map_err(|error| PyValueError::new_err(error.to_string()))
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

    /// MMFF94 partial charges from the numeric PBCI+CHG implementation used
    /// internally by the force-field energy model.
    ///
    /// This is intentionally separate from ``mmff94_charges_typed``: the
    /// latter is the public atom-typed BCI descriptor, while this method
    /// exposes the force-field's RDKit-parity path for diagnostics.
    fn mmff94_charges_numeric(&self) -> PyResult<Vec<f64>> {
        chematic_ff::mmff94_charges_numeric(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Return the numeric MMFF94 atom type id used by the force-field model.
    ///
    /// This is a parity diagnostic corresponding to RDKit's
    /// ``GetMMFFAtomType``. It intentionally exposes numeric ids separately
    /// from :meth:`mmff94_atom_types`, whose strings belong to the older
    /// semantic typing API.
    fn mmff94_numeric_atom_types(&self) -> PyResult<Vec<u8>> {
        chematic_ff::assign_mmff94_numeric_types(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Return the per-atom aromatic flags used by the MMFF94-specific view.
    ///
    /// This diagnostic makes representation-dependent aromaticity decisions
    /// auditable without conflating them with chematic's general aromaticity.
    fn mmff94_numeric_aromatic_flags(&self) -> PyResult<Vec<bool>> {
        let rings = chematic_perception::find_symmetrized_sssr(&self.inner)
            .rings()
            .to_vec();
        let view = chematic_ff::compute_mmff94_aromatic_view(&self.inner, &rings)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok((0..view.atom_count())
            .map(|i| view.atom(chematic_core::AtomIdx(i as u32)).aromatic)
            .collect())
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

    /// Labute approximate surface area (RDKit ``CalcLabuteASA``), the same
    /// value as ``descriptors()["labute_asa"]`` without the other descriptors.
    #[getter]
    fn labute_asa(&self) -> f64 {
        chematic_chem::labute_asa(&self.inner)
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
}
