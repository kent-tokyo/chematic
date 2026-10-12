use super::*;

use chematic_3d::rdkit_embed::{RdkitEmbedError, RdkitEmbedOptions};

type PyRdkitMmffConformer = (i32, f64, Vec<[f64; 3]>);
type PyRdkitMmffConformerDiagnostics = (i32, f64, u32, Option<f64>, Vec<[f64; 3]>);

#[allow(clippy::too_many_arguments)]
fn rdkit_embed_options(
    random_seed: i32,
    max_iterations: u32,
    use_exp_torsion_angle_prefs: bool,
    use_basic_knowledge: bool,
    et_version: u32,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
    use_macrocycle_14_config: bool,
    enforce_chirality: bool,
    ignore_smoothing_failures: bool,
) -> PyResult<RdkitEmbedOptions> {
    if !(1..=2).contains(&et_version) {
        return Err(PyValueError::new_err(
            "Only version 1 and 2 of the experimental torsion-angle preferences (ETversion) supported",
        ));
    }
    Ok(RdkitEmbedOptions {
        random_seed,
        max_iterations,
        use_exp_torsion_angle_prefs,
        use_basic_knowledge,
        et_version,
        use_small_ring_torsions,
        use_macrocycle_torsions,
        use_macrocycle14config: use_macrocycle_14_config,
        enforce_chirality,
        force_trans_amides: true,
        ignore_smoothing_failures,
    })
}

fn map_rdkit_embed_error(error: RdkitEmbedError) -> PyErr {
    match error {
        RdkitEmbedError::Failed => pyo3::exceptions::PyRuntimeError::new_err(error.to_string()),
        _ => PyValueError::new_err(error.to_string()),
    }
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

    /// Canonical SMILES exactly as RDKit 2026.03.1 writes it.
    ///
    /// For a molecule read with :func:`from_smiles` this is the string
    /// ``Chem.MolToSmiles(Chem.MolFromSmiles(s))`` returns (isomeric,
    /// canonical, default parameters) — a port of RDKit's parse,
    /// sanitization, legacy stereo perception and canonical SMILES writer.
    /// It is separate from :attr:`smiles` (chematic's own canonical SMILES),
    /// which it does not change.
    ///
    /// Raises ``ValueError`` instead of returning a string when the molecule
    /// uses a feature the port does not model (message starting
    /// ``"RDKit-compatible SMILES: unsupported input"``, e.g. a molecule not
    /// read from SMILES) or when RDKit's
    /// sanitization would reject it (``"RDKit-compatible SMILES:
    /// sanitization failed"``; ``Chem.MolFromSmiles`` returns ``None``).
    ///
    ///     chematic.from_smiles("OC(=O)[C@@H]1CCCN1").rdkit_smiles
    ///     # 'O=C(O)[C@@H]1CCCN1'
    #[getter]
    fn rdkit_smiles(&self) -> PyResult<String> {
        chematic_smiles::rdkit_canonical_smiles(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// SMILES exactly as RDKit 2026.03.1's ``Chem.MolToSmiles`` writes it
    /// with the given options (``doRandom`` is always false):
    /// ``isomeric`` = ``isomericSmiles``, ``kekule`` = ``kekuleSmiles``,
    /// ``canonical``, ``all_bonds_explicit`` = ``allBondsExplicit``,
    /// ``all_hs_explicit`` = ``allHsExplicit`` and ``rooted_at_atom`` =
    /// ``rootedAtAtom`` (``None``: -1). With the defaults it equals
    /// :attr:`rdkit_smiles`. Raises ``ValueError`` like :attr:`rdkit_smiles`.
    ///
    ///     chematic.from_smiles("c1ccccc1O").rdkit_smiles_with(kekule=True)
    ///     # 'OC1=CC=CC=C1'
    #[pyo3(signature = (*, isomeric = true, kekule = false, canonical = true,
                        all_bonds_explicit = false, all_hs_explicit = false,
                        rooted_at_atom = None))]
    fn rdkit_smiles_with(
        &self,
        isomeric: bool,
        kekule: bool,
        canonical: bool,
        all_bonds_explicit: bool,
        all_hs_explicit: bool,
        rooted_at_atom: Option<usize>,
    ) -> PyResult<String> {
        let params = chematic_smiles::RdkitSmilesParams {
            isomeric,
            kekule,
            canonical,
            all_bonds_explicit,
            all_hs_explicit,
            rooted_at_atom,
        };
        chematic_smiles::rdkit_smiles(&self.inner, &params)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit 2026.09.1's ``AllChem.EmbedMolecule(mol, randomSeed=...)``
    /// (ETKDGv3 defaults), bit for bit: the 3D coordinates of every atom,
    /// one ``[x, y, z]`` per atom in this molecule's atom order.
    ///
    /// Call it on the explicit-hydrogen molecule, as RDKit is called on
    /// ``Chem.AddHs(m)``::
    ///
    ///     mh = chematic.from_smiles("CCO").add_hydrogens()
    ///     xyz = mh.rdkit_embed(random_seed=42)
    ///     # == Chem.AddHs(Chem.MolFromSmiles("CCO")) embedded with
    ///     #    AllChem.EmbedMolecule(mh, randomSeed=42)
    ///
    /// The keyword arguments are ``EmbedMolecule``'s (``maxAttempts`` =
    /// ``max_iterations``, 0: RDKit's default of 10 x atoms). Raises
    /// ``RuntimeError`` when RDKit's embedding fails (``EmbedMolecule``
    /// returns -1) and ``ValueError`` for inputs the port does not model
    /// (e.g. ``random_seed=-1``) or that make RDKit raise. Molecules with
    /// several fragments are embedded fragment by fragment, as RDKit does.
    /// chematic's own 3D methods are unchanged.
    #[pyo3(signature = (*, random_seed = 42, max_iterations = 0,
                        use_exp_torsion_angle_prefs = true, use_basic_knowledge = true,
                        et_version = 2, use_small_ring_torsions = false,
                        use_macrocycle_torsions = true, use_macrocycle_14_config = true,
                        enforce_chirality = true, ignore_smoothing_failures = false))]
    #[allow(clippy::too_many_arguments)]
    fn rdkit_embed(
        &self,
        py: Python<'_>,
        random_seed: i32,
        max_iterations: u32,
        use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool,
        et_version: u32,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle_14_config: bool,
        enforce_chirality: bool,
        ignore_smoothing_failures: bool,
    ) -> PyResult<Vec<[f64; 3]>> {
        use chematic_3d::rdkit_embed::rdkit_embed_molecule;
        let opts = rdkit_embed_options(
            random_seed,
            max_iterations,
            use_exp_torsion_angle_prefs,
            use_basic_knowledge,
            et_version,
            use_small_ring_torsions,
            use_macrocycle_torsions,
            use_macrocycle_14_config,
            enforce_chirality,
            ignore_smoothing_failures,
        )?;
        let mol = self.inner.clone();
        let res = py.detach(move || rdkit_embed_molecule(&mol, &opts));
        res.map_err(map_rdkit_embed_error)
    }

    /// Seeded metric-matrix coordinates before DG/ETKDG minimization,
    /// matching RDKit 2026.09.1's ``onlyInitialEmbedding=True`` boundary.
    #[pyo3(signature = (random_seed = 42, max_iterations = 0))]
    fn rdkit_embed_initial(
        &self,
        py: Python<'_>,
        random_seed: i32,
        max_iterations: u32,
    ) -> PyResult<Vec<[f64; 3]>> {
        use chematic_3d::rdkit_embed::rdkit_embed_initial_molecule;
        let opts = rdkit_embed_options(
            random_seed,
            max_iterations,
            true,
            true,
            2,
            false,
            true,
            true,
            true,
            false,
        )?;
        let mol = self.inner.clone();
        py.detach(move || rdkit_embed_initial_molecule(&mol, &opts))
            .map_err(map_rdkit_embed_error)
    }

    /// RDKit 2026.09.1's ``AllChem.EmbedMultipleConfs`` with ETKDGv3
    /// defaults. Bounds and torsion preferences are shared across attempts.
    /// Failed attempts are omitted, as in RDKit's returned conformer-id list.
    #[pyo3(signature = (num_confs, *, random_seed = 42, max_iterations = 0,
                        use_exp_torsion_angle_prefs = true, use_basic_knowledge = true,
                        et_version = 2, use_small_ring_torsions = false,
                        use_macrocycle_torsions = true, use_macrocycle_14_config = true,
                        enforce_chirality = true, ignore_smoothing_failures = false))]
    #[allow(clippy::too_many_arguments)]
    fn rdkit_embed_multiple_confs(
        &self,
        py: Python<'_>,
        num_confs: usize,
        random_seed: i32,
        max_iterations: u32,
        use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool,
        et_version: u32,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle_14_config: bool,
        enforce_chirality: bool,
        ignore_smoothing_failures: bool,
    ) -> PyResult<Vec<Vec<[f64; 3]>>> {
        use chematic_3d::rdkit_embed::rdkit_embed_multiple_confs;
        let opts = rdkit_embed_options(
            random_seed,
            max_iterations,
            use_exp_torsion_angle_prefs,
            use_basic_knowledge,
            et_version,
            use_small_ring_torsions,
            use_macrocycle_torsions,
            use_macrocycle_14_config,
            enforce_chirality,
            ignore_smoothing_failures,
        )?;
        let mol = self.inner.clone();
        py.detach(move || rdkit_embed_multiple_confs(&mol, num_confs, &opts))
            .map_err(map_rdkit_embed_error)
    }

    /// Generate a deterministic ETKDG ensemble and optimize all conformers
    /// with RDKit-compatible MMFF94. Returns ``(status, energy, coords)`` per
    /// retained conformer in candidate order.
    #[pyo3(signature = (num_confs, *, random_seed = 42, max_iterations = 0,
                        max_mmff_iters = 200, retry_mmff_iters = 0,
                        retry_candidates = 0, retry_only_if_none_converged = false,
                        use_exp_torsion_angle_prefs = true,
                        use_basic_knowledge = true, et_version = 2,
                        use_small_ring_torsions = false, use_macrocycle_torsions = true,
                        use_macrocycle_14_config = true, enforce_chirality = true,
                        ignore_smoothing_failures = false))]
    #[allow(clippy::too_many_arguments)]
    fn rdkit_embed_multiple_confs_mmff(
        &self,
        py: Python<'_>,
        num_confs: usize,
        random_seed: i32,
        max_iterations: u32,
        max_mmff_iters: u32,
        retry_mmff_iters: u32,
        retry_candidates: usize,
        retry_only_if_none_converged: bool,
        use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool,
        et_version: u32,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle_14_config: bool,
        enforce_chirality: bool,
        ignore_smoothing_failures: bool,
    ) -> PyResult<Vec<PyRdkitMmffConformer>> {
        use chematic_3d::rdkit_embed::{
            rdkit_embed_multiple_confs_mmff, rdkit_embed_multiple_confs_mmff_adaptive,
            rdkit_embed_multiple_confs_mmff_rescue,
        };
        let opts = rdkit_embed_options(
            random_seed,
            max_iterations,
            use_exp_torsion_angle_prefs,
            use_basic_knowledge,
            et_version,
            use_small_ring_torsions,
            use_macrocycle_torsions,
            use_macrocycle_14_config,
            enforce_chirality,
            ignore_smoothing_failures,
        )?;
        let mol = self.inner.clone();
        py.detach(move || {
            if retry_mmff_iters == 0 || retry_candidates == 0 {
                rdkit_embed_multiple_confs_mmff(&mol, num_confs, &opts, max_mmff_iters)
            } else if retry_only_if_none_converged {
                rdkit_embed_multiple_confs_mmff_rescue(
                    &mol,
                    num_confs,
                    &opts,
                    max_mmff_iters,
                    retry_mmff_iters,
                    retry_candidates,
                )
            } else {
                rdkit_embed_multiple_confs_mmff_adaptive(
                    &mol,
                    num_confs,
                    &opts,
                    max_mmff_iters,
                    retry_mmff_iters,
                    retry_candidates,
                )
            }
        })
        .map(|conformers| {
            conformers
                .into_iter()
                .map(|conformer| (conformer.status, conformer.energy, conformer.coords))
                .collect()
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Generate and optimize an ETKDG ensemble while returning BFGS
    /// diagnostics as ``(status, energy, iterations, gradient_test, coords)``.
    #[pyo3(signature = (num_confs, *, random_seed = 42, max_iterations = 0,
                        max_mmff_iters = 200,
                        use_exp_torsion_angle_prefs = true,
                        use_basic_knowledge = true, et_version = 2,
                        use_small_ring_torsions = false, use_macrocycle_torsions = true,
                        use_macrocycle_14_config = true, enforce_chirality = true,
                        ignore_smoothing_failures = false))]
    #[allow(clippy::too_many_arguments)]
    fn rdkit_embed_multiple_confs_mmff_diagnostics(
        &self,
        py: Python<'_>,
        num_confs: usize,
        random_seed: i32,
        max_iterations: u32,
        max_mmff_iters: u32,
        use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool,
        et_version: u32,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle_14_config: bool,
        enforce_chirality: bool,
        ignore_smoothing_failures: bool,
    ) -> PyResult<Vec<PyRdkitMmffConformerDiagnostics>> {
        use chematic_3d::rdkit_embed::rdkit_embed_multiple_confs_mmff_diagnostics;
        let opts = rdkit_embed_options(
            random_seed,
            max_iterations,
            use_exp_torsion_angle_prefs,
            use_basic_knowledge,
            et_version,
            use_small_ring_torsions,
            use_macrocycle_torsions,
            use_macrocycle_14_config,
            enforce_chirality,
            ignore_smoothing_failures,
        )?;
        let mol = self.inner.clone();
        py.detach(move || {
            rdkit_embed_multiple_confs_mmff_diagnostics(&mol, num_confs, &opts, max_mmff_iters)
        })
        .map(|conformers| {
            conformers
                .into_iter()
                .map(|result| {
                    let conformer = result.conformer;
                    (
                        conformer.status,
                        conformer.energy,
                        result.iterations,
                        result.final_gradient_test,
                        conformer.coords,
                    )
                })
                .collect()
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Generate and MMFF94-optimize an ETKDG ensemble, retrying only
    /// high-gradient candidates when no initial conformer converged.
    #[pyo3(signature = (num_confs, *, random_seed = 42, max_iterations = 0,
                        max_mmff_iters = 200, retry_mmff_iters = 100,
                        retry_candidates = 1, minimum_gradient_test = 0.05,
                        use_exp_torsion_angle_prefs = true,
                        use_basic_knowledge = true, et_version = 2,
                        use_small_ring_torsions = false, use_macrocycle_torsions = true,
                        use_macrocycle_14_config = true, enforce_chirality = true,
                        ignore_smoothing_failures = false))]
    #[allow(clippy::too_many_arguments)]
    fn rdkit_embed_multiple_confs_mmff_gradient_rescue(
        &self,
        py: Python<'_>,
        num_confs: usize,
        random_seed: i32,
        max_iterations: u32,
        max_mmff_iters: u32,
        retry_mmff_iters: u32,
        retry_candidates: usize,
        minimum_gradient_test: f64,
        use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool,
        et_version: u32,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle_14_config: bool,
        enforce_chirality: bool,
        ignore_smoothing_failures: bool,
    ) -> PyResult<Vec<PyRdkitMmffConformer>> {
        use chematic_3d::rdkit_embed::rdkit_embed_multiple_confs_mmff_gradient_rescue;
        if !minimum_gradient_test.is_finite() || minimum_gradient_test < 0.0 {
            return Err(PyValueError::new_err(
                "minimum_gradient_test must be finite and non-negative",
            ));
        }
        let opts = rdkit_embed_options(
            random_seed,
            max_iterations,
            use_exp_torsion_angle_prefs,
            use_basic_knowledge,
            et_version,
            use_small_ring_torsions,
            use_macrocycle_torsions,
            use_macrocycle_14_config,
            enforce_chirality,
            ignore_smoothing_failures,
        )?;
        let mol = self.inner.clone();
        py.detach(move || {
            rdkit_embed_multiple_confs_mmff_gradient_rescue(
                &mol,
                num_confs,
                &opts,
                max_mmff_iters,
                retry_mmff_iters,
                retry_candidates,
                minimum_gradient_test,
            )
        })
        .map(|conformers| {
            conformers
                .into_iter()
                .map(|conformer| (conformer.status, conformer.energy, conformer.coords))
                .collect()
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// RDKit 2026.09.1's ``rdDistGeom.GetMoleculeBoundsMatrix(mol, ...)``
    /// (``scaleVDW=False``), bit for bit: the distance-geometry bounds
    /// matrix as a list of rows, upper bounds above the diagonal and lower
    /// bounds below it. Call it on the explicit-hydrogen molecule, as RDKit
    /// is called on ``Chem.AddHs(m)``::
    ///
    ///     bm = chematic.from_smiles("CCO").add_hydrogens().rdkit_bounds_matrix()
    ///
    /// Raises ``ValueError`` for inputs the port does not model.
    #[pyo3(signature = (*, set15bounds = true, do_triangle_smoothing = true,
                        use_macrocycle_14_config = false))]
    fn rdkit_bounds_matrix(
        &self,
        set15bounds: bool,
        do_triangle_smoothing: bool,
        use_macrocycle_14_config: bool,
    ) -> PyResult<Vec<Vec<f64>>> {
        chematic_3d::rdkit_embed::rdkit_bounds_matrix(
            &self.inner,
            set15bounds,
            do_triangle_smoothing,
            use_macrocycle_14_config,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// 2D coordinates exactly as RDKit 2026.03.1's default depiction gives
    /// them: for a molecule read with :func:`from_smiles` this is
    /// ``[[p.x, p.y] for p in m.GetConformer().GetPositions()]`` after
    /// ``m = Chem.MolFromSmiles(s); rdDepictor.Compute2DCoords(m)`` (RDKit's
    /// own depictor with its defaults, not CoordGen), bit for bit.
    ///
    /// One ``[x, y]`` per atom of RDKit's molecule, in RDKit's atom order —
    /// this molecule's atom order unless ``MolFromSmiles`` removes hydrogen
    /// atoms written in brackets. chematic's own 2D layout (used by
    /// :meth:`depict_data`, the SVG and MOL writers) is separate and
    /// unchanged. Raises ``ValueError`` like :attr:`rdkit_smiles`.
    ///
    ///     chematic.from_smiles("CCO").rdkit_2d_coords()[1]
    ///     # [0.0, 0.5000000000000001]
    fn rdkit_2d_coords(&self) -> PyResult<Vec<[f64; 2]>> {
        chematic_smiles::rdkit_2d_coords(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// The MOL block RDKit 2026.03.1 writes for this molecule with its
    /// default 2D depiction: ``Chem.MolToMolBlock(m)`` after
    /// ``m = Chem.MolFromSmiles(s); rdDepictor.Compute2DCoords(m)`` for the
    /// SMILES ``s`` the molecule was read from (V2000, or V3000 where RDKit
    /// switches to it: dative bonds, more than 999 atoms or bonds;
    /// coordinates of :meth:`rdkit_2d_coords`, RDKit's kekulization and
    /// wedge bonds).
    ///
    /// Raises ``ValueError`` like :attr:`rdkit_smiles`.
    ///
    ///     print(chematic.from_smiles("C[C@H](O)F").rdkit_mol_block_2d())
    fn rdkit_mol_block_2d(&self) -> PyResult<String> {
        chematic_smiles::rdkit_mol_block_2d(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Canonical SMILES together with the atom output order.
    ///
    /// Returns ``(smiles, order)`` where ``smiles`` equals :attr:`smiles` and
    /// ``order[k]`` is the index in this molecule of the ``k``-th atom written
    /// in the string — also that atom's index after
    /// ``chematic.from_smiles(smiles)``. This is the SMILES DFS visit order
    /// (RDKit's ``_smilesAtomOutputOrder``), not a canonical rank.
    ///
    ///     smi, order = chematic.from_smiles("OCC").smiles_with_atom_order()
    ///     # ("C(C)O", [1, 2, 0])
    fn smiles_with_atom_order(&self) -> (String, Vec<usize>) {
        let (smiles, order) = chematic_smiles::canonical_smiles_with_atom_order(&self.inner);
        (smiles, order.into_iter().map(|a| a.0 as usize).collect())
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
    ///
    /// ``ignore_interfrag_interactions`` as in :meth:`mmff94_energy_breakdown`.
    #[pyo3(signature = (coords, ignore_interfrag_interactions = false))]
    fn mmff94_total_energy(
        &self,
        coords: Vec<[f64; 3]>,
        ignore_interfrag_interactions: bool,
    ) -> f64 {
        chematic_ff::mmff94_energy_breakdown_with_options(
            &self.inner,
            &coords,
            chematic_ff::Mmff94Options {
                ignore_interfragment_interactions: ignore_interfrag_interactions,
            },
        )
        .map(|b| b.total)
        .unwrap_or(0.0)
    }

    /// RDKit-compatible MMFF94 energy: ``MMFFGetMoleculeForceField(m,
    /// MMFFGetMoleculeProperties(m), nonBondedThresh,
    /// ignoreInterfragInteractions).CalcEnergy()`` for ``coords`` (all atoms,
    /// in the order of :meth:`add_hydrogens`, which is RDKit's ``AddHs``
    /// order). ``variant`` is ``"MMFF94"`` or ``"MMFF94s"``. Raises
    /// ``ValueError`` when MMFF94 typing fails.
    #[pyo3(signature = (coords, variant = "MMFF94", non_bonded_thresh = 100.0, ignore_interfrag_interactions = true))]
    fn rdkit_mmff_energy(
        &self,
        coords: Vec<[f64; 3]>,
        variant: &str,
        non_bonded_thresh: f64,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<f64> {
        let ff = rdkit_mmff_field(
            &self.inner,
            &coords,
            mmff_variant(variant)?,
            non_bonded_thresh,
            ignore_interfrag_interactions,
            [true; 7],
        )?;
        Ok(ff.energy(coords.as_flattened()))
    }

    /// RDKit-compatible ``MMFFOptimizeMolecule(m, maxIters, nonBondedThresh,
    /// ignoreInterfragInteractions)`` from ``coords``: returns
    /// ``(status, energy, coords)`` with RDKit's status (0 converged, 1 more
    /// iterations needed) and final energy. ``variant`` is ``"MMFF94"`` or
    /// ``"MMFF94s"`` (``"MMFF94S"`` is accepted as ``"MMFF94s"``; note RDKit
    /// itself silently treats any spelling other than ``"MMFF94s"`` as
    /// MMFF94).
    #[pyo3(signature = (coords, max_iters = 200, variant = "MMFF94", non_bonded_thresh = 100.0, ignore_interfrag_interactions = true))]
    fn rdkit_mmff_optimize(
        &self,
        coords: Vec<[f64; 3]>,
        max_iters: u32,
        variant: &str,
        non_bonded_thresh: f64,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<(i32, f64, Vec<[f64; 3]>)> {
        let ff = rdkit_mmff_field(
            &self.inner,
            &coords,
            mmff_variant(variant)?,
            non_bonded_thresh,
            ignore_interfrag_interactions,
            [true; 7],
        )?;
        let mut pos = coords.as_flattened().to_vec();
        let (status, energy) = ff.optimize(&mut pos, max_iters);
        if status < 0 {
            return Err(PyValueError::new_err("bad direction in linearSearch"));
        }
        Ok((
            status,
            energy,
            pos.chunks(3).map(|c| [c[0], c[1], c[2]]).collect(),
        ))
    }

    /// Add implicit hydrogens at RDKit-compatible 3D positions, refine the
    /// complete molecule with MMFF94, and return the refined coordinates for
    /// the original atoms. This avoids the cost of embedding every hydrogen
    /// while still ranking conformers with all X-H and non-bonded terms.
    ///
    /// Returns ``(heavy_coords, all_atom_coords, energy, converged,
    /// iterations, max_residual_force)``. Pair ``all_atom_coords`` with
    /// ``self.add_hydrogens()``. ``hydrogen_variant=0`` uses RDKit's default
    /// placement; positive values deterministically sample terminal
    /// heteroatom-H rotamers. Raises ``ValueError`` when hydrogen placement or
    /// MMFF94 typing fails.
    #[pyo3(signature = (coords, max_iters = 200, hydrogen_variant = 0))]
    fn refine_mmff94_with_explicit_hydrogens(
        &self,
        py: Python<'_>,
        coords: Vec<[f64; 3]>,
        max_iters: usize,
        hydrogen_variant: u64,
    ) -> PyResult<PyExplicitHydrogenMmffResult> {
        let mol = self.inner.clone();
        let coords = flat_to_coords3d(&coords);
        let result = py.detach(move || {
            chematic_3d::refine_mmff94_with_explicit_hydrogens_variant(
                &mol,
                &coords,
                max_iters,
                hydrogen_variant,
            )
        });
        result
            .map(|result| {
                (
                    result
                        .heavy_coords
                        .points
                        .into_iter()
                        .map(|point| [point.x, point.y, point.z])
                        .collect(),
                    result
                        .all_atom_coords
                        .points
                        .into_iter()
                        .map(|point| [point.x, point.y, point.z])
                        .collect(),
                    result.energy,
                    result.converged,
                    result.iterations,
                    result.max_residual_force,
                )
            })
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// RDKit-compatible UFF energy: ``UFFGetMoleculeForceField(m,
    /// vdwThresh, ignoreInterfragInteractions=...).CalcEnergy()`` for
    /// ``coords`` (all atoms, in :meth:`add_hydrogens` order, which is
    /// RDKit's ``AddHs`` order).
    #[pyo3(signature = (coords, vdw_thresh = 10.0, ignore_interfrag_interactions = true))]
    fn rdkit_uff_energy(
        &self,
        coords: Vec<[f64; 3]>,
        vdw_thresh: f64,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<f64> {
        let ff = rdkit_uff_field(
            &self.inner,
            &coords,
            vdw_thresh,
            ignore_interfrag_interactions,
        )?;
        Ok(ff.energy(coords.as_flattened()))
    }

    /// RDKit-compatible UFF gradient (``CalcGrad()``, unscaled) for
    /// ``coords``, flattened.
    #[pyo3(signature = (coords, vdw_thresh = 10.0, ignore_interfrag_interactions = true))]
    fn rdkit_uff_gradient(
        &self,
        coords: Vec<[f64; 3]>,
        vdw_thresh: f64,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<Vec<f64>> {
        let ff = rdkit_uff_field(
            &self.inner,
            &coords,
            vdw_thresh,
            ignore_interfrag_interactions,
        )?;
        let pos = coords.as_flattened();
        let mut grad = vec![0.0; pos.len()];
        ff.gradient(pos, &mut grad);
        Ok(grad)
    }

    /// RDKit-compatible ``UFFOptimizeMolecule(m, maxIters, vdwThresh,
    /// ignoreInterfragInteractions=...)`` from ``coords``: returns
    /// ``(status, energy, coords)`` with RDKit's status (0 converged, 1 more
    /// iterations needed) and final energy.
    #[pyo3(signature = (coords, max_iters = 200, vdw_thresh = 10.0, ignore_interfrag_interactions = true))]
    fn rdkit_uff_optimize(
        &self,
        coords: Vec<[f64; 3]>,
        max_iters: u32,
        vdw_thresh: f64,
        ignore_interfrag_interactions: bool,
    ) -> PyResult<(i32, f64, Vec<[f64; 3]>)> {
        let ff = rdkit_uff_field(
            &self.inner,
            &coords,
            vdw_thresh,
            ignore_interfrag_interactions,
        )?;
        let mut pos = coords.as_flattened().to_vec();
        let (status, energy) = ff.optimize(&mut pos, max_iters);
        if status < 0 {
            return Err(PyValueError::new_err("bad direction in linearSearch"));
        }
        Ok((
            status,
            energy,
            pos.chunks(3).map(|c| [c[0], c[1], c[2]]).collect(),
        ))
    }

    /// Diagnostic: the RDKit UFF terms ``(kind, atoms, params)`` in RDKit's
    /// contribution order.
    #[pyo3(signature = (coords, vdw_thresh = 10.0))]
    fn _rdkit_uff_terms(&self, coords: Vec<[f64; 3]>, vdw_thresh: f64) -> PyResult<Vec<PyUffTerm>> {
        let ff = rdkit_uff_field(&self.inner, &coords, vdw_thresh, true)?;
        Ok(ff
            .debug_terms()
            .into_iter()
            .map(|(k, a, p)| (k.to_string(), a, p))
            .collect())
    }

    /// RDKit ``UFFHasAllMoleculeParams``: whether every atom (of this
    /// explicit-hydrogen molecule) has a UFF atom type with parameters.
    fn rdkit_uff_has_all_params(&self) -> bool {
        chematic_ff::rdkit_uff::rdkit_uff_has_all_params(&self.inner)
    }

    /// RDKit's UFF atom labels (``UFF::Tools::getAtomLabel``), one per atom.
    fn rdkit_uff_atom_labels(&self) -> Vec<String> {
        chematic_ff::rdkit_uff::rdkit_uff_atom_labels(&self.inner)
    }

    /// Diagnostic: RDKit-compatible MMFF94 energy and unscaled gradient with
    /// only the term classes in ``terms`` (``bond``, ``angle``,
    /// ``stretch_bend``, ``oop``, ``torsion``, ``vdw``, ``ele``).
    #[pyo3(signature = (coords, terms = None, non_bonded_thresh = 100.0, variant = "MMFF94"))]
    fn _rdkit_mmff_terms(
        &self,
        coords: Vec<[f64; 3]>,
        terms: Option<Vec<String>>,
        non_bonded_thresh: f64,
        variant: &str,
    ) -> PyResult<(f64, Vec<f64>)> {
        let mut mask = [true; 7];
        if let Some(t) = terms {
            for (k, name) in chematic_ff::rdkit_mmff::RDKIT_MMFF_TERMS.iter().enumerate() {
                mask[k] = t.iter().any(|x| x == name);
            }
        }
        let ff = rdkit_mmff_field(
            &self.inner,
            &coords,
            mmff_variant(variant)?,
            non_bonded_thresh,
            true,
            mask,
        )?;
        let pos = coords.as_flattened();
        // Match RDKit's public call order used by the comparison surface:
        // CalcEnergy() followed by CalcGrad(). Keeping the energy walk first
        // also leaves the term and coordinate data warm for the gradient.
        let energy = ff.energy(pos);
        let mut grad = vec![0.0; pos.len()];
        ff.gradient(pos, &mut grad);
        Ok((energy, grad))
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

    /// Serialize to MDL MOL V2000 format with 2D coordinates (Å).
    ///
    /// Equivalent to RDKit's ``Chem.MolToMolBlock(mol)``: every molecule is
    /// written on a 2D layout, with stereo drawn as wedges and E/Z geometry.
    /// Use :meth:`to_mol2` for Tripos format, :meth:`to_pdb` for PDB with 3D coords.
    ///
    ///     block = mol.to_mol_block()
    ///     with open("molecule.mol", "w") as f:
    ///         f.write(block)
    ///
    /// A block that loses stereo (a centre left unwedged, an E/Z bond written
    /// "either") is still returned, with a :class:`StereoLossWarning` naming
    /// what is lost; ``strict=True`` raises ``ValueError`` instead.
    /// :meth:`to_mol_block_with_report` returns the block with the list and
    /// does not warn.
    #[pyo3(signature = (strict = false))]
    fn to_mol_block(&self, py: Python<'_>, strict: bool) -> PyResult<String> {
        // Laid out like RDKit's MolToMolBlock (and the WASM binding): the
        // stereo depiction, or the plain layout without stereo (this wrote
        // every atom at the origin for a molecule without stereo).
        // Undeclared E/Z double bonds are written "either", so the layout
        // never reads back as a configuration the molecule does not have.
        let (block, loss) =
            chematic_mol::write_laid_out_mol(&self.inner, &chematic_mol::MolMetadata::default());
        strict_mol_block(py, block, &loss, strict)
    }

    /// :meth:`to_mol_block` plus what the block does not carry of the
    /// molecule's stereo: ``(block, report)`` where ``report`` has
    /// ``centres`` (atom indices written without a wedge), ``double_bonds``
    /// (bond indices of declared E/Z written "either"),
    /// ``non_tetrahedral_centres`` (square-planar atoms) and
    /// ``stereo_groups_dropped``. All empty / ``False``: nothing lost.
    ///
    ///     block, report = mol.to_mol_block_with_report()
    ///     if report["centres"]:
    ///         ...
    fn to_mol_block_with_report<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<(String, Bound<'py, PyDict>)> {
        // Laid out like RDKit's MolToMolBlock (and the WASM binding): the
        // stereo depiction, or the plain layout without stereo (this wrote
        // every atom at the origin for a molecule without stereo).
        // Undeclared E/Z double bonds are written "either", so the layout
        // never reads back as a configuration the molecule does not have.
        let (block, loss) =
            chematic_mol::write_laid_out_mol(&self.inner, &chematic_mol::MolMetadata::default());
        Ok((block, stereo_loss_dict(py, &loss)?))
    }

    /// Serialize this molecule to MDL MOL V2000 format with 2D layout coordinates.
    ///
    /// Preserves the 2D layout obtained from :func:`from_mol_block_with_coords`.
    /// Each element of ``coords`` is an ``[x, y]`` pair in Å.
    ///
    ///     mol, name, coords_2d = chematic.from_mol_block_with_coords(block)
    ///     new_block = mol.to_mol_block_2d(coords_2d, name=name)
    ///
    /// A block that loses stereo warns with :class:`StereoLossWarning`;
    /// ``strict=True`` raises ``ValueError`` instead (see :meth:`to_mol_block`).
    #[pyo3(signature = (coords, name = None, strict = false))]
    fn to_mol_block_2d(
        &self,
        py: Python<'_>,
        coords: Vec<[f64; 2]>,
        name: Option<&str>,
        strict: bool,
    ) -> PyResult<String> {
        let metadata = chematic_mol::MolMetadata {
            name: name.unwrap_or("").to_string(),
            comment: String::new(),
            ..Default::default()
        };
        let coords_2d: Vec<(f64, f64)> = coords.iter().map(|c| (c[0], c[1])).collect();
        let (block, loss) =
            chematic_mol::write_mol_with_stereo_report(&self.inner, &metadata, &coords_2d);
        strict_mol_block(py, block, &loss, strict)
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
    ///
    /// A block that loses stereo warns with :class:`StereoLossWarning`;
    /// ``strict=True`` raises ``ValueError`` instead (see :meth:`to_mol_block`).
    #[pyo3(signature = (coords, name = None, strict = false))]
    fn to_mol_v3000(
        &self,
        py: Python<'_>,
        coords: Vec<[f64; 2]>,
        name: Option<&str>,
        strict: bool,
    ) -> PyResult<String> {
        let metadata = chematic_mol::MolMetadata {
            name: name.unwrap_or("").to_string(),
            comment: String::new(),
            ..Default::default()
        };
        let coords_2d: Vec<(f64, f64)> = coords.iter().map(|c| (c[0], c[1])).collect();
        let (block, loss) =
            chematic_mol::write_mol_v3000_with_stereo_report(&self.inner, &metadata, &coords_2d);
        strict_mol_block(py, block, &loss, strict)
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
    ///     before trusting a result), ``worst_bond_length`` (float, Å), and
    ///     ``rejected_unsound_step`` (bool — the line search rejected at
    ///     least one energy-decreasing step for violating the bond-length
    ///     soundness bound).
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
        d.set_item("worst_bond_length", result.worst_bond_length)?;
        d.set_item("rejected_unsound_step", result.rejected_unsound_step)?;
        Ok(d)
    }
}
