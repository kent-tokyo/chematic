use super::*;

#[pymethods]
impl Mol {
    // -----------------------------------------------------------------------
    // Transformations
    // -----------------------------------------------------------------------

    /// Return the standardized molecule with an optional largest-fragment profile.
    #[pyo3(signature = (largest_fragment_only = false))]
    fn standardize(&self, largest_fragment_only: bool) -> Mol {
        let opts = chematic_chem::StandardizeOptions {
            largest_fragment_only,
            ..Default::default()
        };
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

    /// Enumerate local N/O/S-C-N/O/S cation resonance (RDKit #9458).
    ///
    /// Returns ``forms`` (Mol objects), ``status``, ``scope``, and
    /// ``transforms_attempted``. Indexed symmetry-equivalent placements are
    /// retained. This is not general resonance or tautomer enumeration.
    /// Stereo/query inputs and invalid limits raise ValueError; stereo is
    /// never silently removed. Atom indices, isotopes and H placement survive.
    #[pyo3(signature = (max_forms=32, max_transforms=1024, max_atoms=256))]
    fn local_cation_resonance<'py>(
        &self,
        py: Python<'py>,
        max_forms: usize,
        max_transforms: usize,
        max_atoms: usize,
    ) -> PyResult<Bound<'py, PyDict>> {
        use chematic_chem::LocalResonanceStatus;
        let result = chematic_chem::enumerate_local_cation_resonance(
            &self.inner,
            &chematic_chem::LocalResonanceLimits {
                max_forms,
                max_transforms,
                max_atoms,
            },
        )
        .map_err(|err| PyValueError::new_err(err.to_string()))?;
        let status = match result.status {
            LocalResonanceStatus::CompleteWithinScope => "complete_within_scope",
            LocalResonanceStatus::FormLimitReached => "form_limit_reached",
            LocalResonanceStatus::TransformLimitReached => "transform_limit_reached",
        };
        let forms = PyList::empty(py);
        for form in result.forms {
            forms.append(Mol::bare(form))?;
        }
        let report = PyDict::new(py);
        report.set_item("forms", forms)?;
        report.set_item("status", status)?;
        report.set_item("scope", "local_nos_cation_shift_v1")?;
        report.set_item("transforms_attempted", result.transforms_attempted)?;
        Ok(report)
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

    /// RDKit-compatible stereoisomer enumeration: the isomers
    /// ``rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers`` yields
    /// with ``StereoEnumerationOptions(onlyUnassigned=True, unique=True,
    /// maxIsomers=max_isomers, tryEmbedding=False)``, as Mols parsed from
    /// their RDKit canonical SMILES, sorted by that SMILES.
    ///
    /// Only ``only_unassigned=True`` and ``unique=True`` are supported.
    /// With more flip combinations than ``max_isomers`` RDKit yields a random
    /// sample seeded deterministically from the molecule (``rand=None``);
    /// that sample is reproduced (CPython ``hash`` + ``random.Random``).
    /// Raises ``ValueError`` when the molecule is outside the
    /// RDKit-compatible model (e.g. enhanced stereo groups), or for
    /// ``max_isomers=0`` with more than 16 flips.
    #[pyo3(signature = (only_unassigned = true, unique = true, max_isomers = 1024))]
    fn rdkit_stereoisomers(
        &self,
        only_unassigned: bool,
        unique: bool,
        max_isomers: usize,
    ) -> PyResult<Vec<Mol>> {
        if !only_unassigned || !unique {
            return Err(PyValueError::new_err(
                "only only_unassigned=True and unique=True are supported",
            ));
        }
        let smiles = chematic_smiles::rdkit_stereoisomer_smiles(&self.inner, max_isomers)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        smiles
            .iter()
            .map(|s| {
                chematic_smiles::parse(s)
                    .map(Mol::bare)
                    .map_err(|e| PyValueError::new_err(format!("{s}: {e}")))
            })
            .collect()
    }

    /// RDKit-compatible ``Chem.MolToSmarts(m, isomericSmiles=isomeric,
    /// rootedAtAtom=rooted_at_atom)`` (input atom order, every atom
    /// bracketed, every bond explicit).
    #[pyo3(signature = (isomeric = true, rooted_at_atom = None))]
    fn rdkit_smarts(&self, isomeric: bool, rooted_at_atom: Option<usize>) -> PyResult<String> {
        chematic_smiles::rdkit_smarts(&self.inner, isomeric, rooted_at_atom)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit-compatible ``Chem.MolToPDBBlock(m)``; ``coords`` (one
    /// ``[x, y, z]`` per atom) stands for a conformer, otherwise zero
    /// coordinates are written as RDKit does.
    #[pyo3(signature = (coords = None))]
    fn rdkit_pdb_block(&self, coords: Option<Vec<[f64; 3]>>) -> PyResult<String> {
        chematic_smiles::rdkit_pdb_block(&self.inner, coords.as_deref())
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// ``rdMolAlign.GetAlignmentTransform`` / ``AlignMol`` (RDKit 2026.03.1)
    /// for two conformers of this molecule (coordinates in
    /// ``Chem.MolFromSmiles`` atom order): ``(rmsd, transform)`` with the
    /// 4x4 transform (rows) taking ``probe_coords`` onto ``ref_coords``.
    /// ``atom_map`` is a list of ``(probe atom, reference atom)`` pairs;
    /// without it the first substructure match is used.
    #[pyo3(signature = (probe_coords, ref_coords, atom_map = None, weights = None, reflect = false, max_iterations = 50))]
    fn rdkit_align(
        &self,
        probe_coords: Vec<[f64; 3]>,
        ref_coords: Vec<[f64; 3]>,
        atom_map: Option<Vec<(usize, usize)>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> PyResult<(f64, [[f64; 4]; 4])> {
        let res = chematic_smiles::rdkit_align_mol(
            &self.inner,
            &probe_coords,
            &ref_coords,
            atom_map.as_deref(),
            weights.as_deref(),
            reflect,
            max_iterations,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok((res.rmsd, res.transform.rows()))
    }

    /// ``rdMolAlign.GetBestRMS(prb, ref, maxMatches=max_matches,
    /// symmetrizeConjugatedTerminalGroups=..., weights=weights)`` (RDKit
    /// 2026.03.1) for two conformers of this molecule.
    #[pyo3(signature = (probe_coords, ref_coords, max_matches = 1_000_000, symmetrize_conjugated_terminal_groups = true, weights = None))]
    fn rdkit_best_rms(
        &self,
        probe_coords: Vec<[f64; 3]>,
        ref_coords: Vec<[f64; 3]>,
        max_matches: usize,
        symmetrize_conjugated_terminal_groups: bool,
        weights: Option<Vec<f64>>,
    ) -> PyResult<f64> {
        chematic_smiles::rdkit_best_rms(
            &self.inner,
            &probe_coords,
            &ref_coords,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            weights.as_deref(),
        )
        .map(|res| res.rmsd)
        .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// ``rdMolAlign.GetAllConformerBestRMS`` (RDKit 2026.09.1): best-fit
    /// RMSDs for every lower-triangle pair of ``conformers``. The returned
    /// order is ``(1,0), (2,0), (2,1), ...`` and atom mappings are computed
    /// once for the whole batch.
    #[pyo3(signature = (conformers, max_matches = 1_000_000, symmetrize_conjugated_terminal_groups = true, weights = None))]
    fn rdkit_all_conformer_best_rms(
        &self,
        conformers: Vec<Vec<[f64; 3]>>,
        max_matches: usize,
        symmetrize_conjugated_terminal_groups: bool,
        weights: Option<Vec<f64>>,
    ) -> PyResult<Vec<f64>> {
        chematic_smiles::rdkit_all_conformer_best_rms(
            &self.inner,
            &conformers,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            weights.as_deref(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit 2026.09.1 ``rdMolAlign.GetAllConformerBestRMSToRef``: best-fit
    /// RMSDs from every conformer of this probe molecule to every conformer
    /// of ``reference``. Results are reference-major:
    /// ``(ref0, probe0), (ref0, probe1), ...``. Atom mappings are computed
    /// once for the whole Cartesian product.
    #[pyo3(signature = (reference, probe_conformers, reference_conformers, max_matches = 1_000_000, symmetrize_conjugated_terminal_groups = true, weights = None))]
    fn rdkit_all_conformer_best_rms_to_reference(
        &self,
        reference: &Mol,
        probe_conformers: Vec<Vec<[f64; 3]>>,
        reference_conformers: Vec<Vec<[f64; 3]>>,
        max_matches: usize,
        symmetrize_conjugated_terminal_groups: bool,
        weights: Option<Vec<f64>>,
    ) -> PyResult<Vec<f64>> {
        chematic_smiles::rdkit_all_conformer_best_rms_to_ref(
            &self.inner,
            &reference.inner,
            &probe_conformers,
            &reference_conformers,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            weights.as_deref(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// ``rdMolAlign.GetBestAlignmentTransform`` (RDKit 2026.03.1):
    /// ``(rmsd, transform, atom_map)`` of the best fit over all matches,
    /// ``atom_map`` as ``(probe atom, reference atom)`` pairs.
    #[pyo3(signature = (probe_coords, ref_coords, max_matches = 1_000_000, symmetrize_conjugated_terminal_groups = true, weights = None))]
    fn rdkit_best_alignment(
        &self,
        probe_coords: Vec<[f64; 3]>,
        ref_coords: Vec<[f64; 3]>,
        max_matches: usize,
        symmetrize_conjugated_terminal_groups: bool,
        weights: Option<Vec<f64>>,
    ) -> PyResult<RdkitBestAlignment> {
        let res = chematic_smiles::rdkit_best_rms(
            &self.inner,
            &probe_coords,
            &ref_coords,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            weights.as_deref(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok((res.rmsd, res.transform.rows(), res.atom_map))
    }

    /// ``rdMolAlign.CalcRMS(prb, ref, maxMatches=max_matches,
    /// symmetrizeConjugatedTerminalGroups=..., weights=weights)`` (RDKit
    /// 2026.03.1): the smallest RMSD over all matches, without aligning.
    #[pyo3(signature = (probe_coords, ref_coords, max_matches = 1_000_000, symmetrize_conjugated_terminal_groups = true, weights = None))]
    fn rdkit_calc_rms(
        &self,
        probe_coords: Vec<[f64; 3]>,
        ref_coords: Vec<[f64; 3]>,
        max_matches: usize,
        symmetrize_conjugated_terminal_groups: bool,
        weights: Option<Vec<f64>>,
    ) -> PyResult<f64> {
        chematic_smiles::rdkit_calc_rms(
            &self.inner,
            &probe_coords,
            &ref_coords,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            weights.as_deref(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit 2026.03.1's ``Chem.MolToCXSmiles(m, params)`` (all CXSMILES
    /// fields) for ``m = Chem.MolFromSmiles(s)``: the SMILES with dative
    /// bonds written as single bonds plus the extension (radicals ``^n:``,
    /// coordinate bonds ``C:``, ring double bonds of unknown geometry
    /// ``ctu:``). Molecules with enhanced stereo groups are refused
    /// (``ValueError``)::
    ///
    ///     chematic.from_smiles("C[CH2]").rdkit_cx_smiles()  # '[CH2]C |^1:0|'
    #[pyo3(signature = (*, isomeric = true, kekule = false, canonical = true,
                        all_bonds_explicit = false, all_hs_explicit = false,
                        rooted_at_atom = None))]
    fn rdkit_cx_smiles(
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
        chematic_smiles::rdkit_cx_smiles(&self.inner, &params)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit 2026.03.1's ``Chem.MolToRandomSmilesVect(m, n, randomSeed,
    /// isomericSmiles, kekuleSmiles, allBondsExplicit, allHsExplicit)`` for
    /// ``m = Chem.MolFromSmiles(s)``, bit for bit (RDKit's
    /// ``boost::minstd_rand`` generator). ``random_seed > 0`` reseeds the
    /// generator; ``0`` continues a process-wide generator that starts as
    /// RDKit's (seed 42)::
    ///
    ///     chematic.from_smiles("CCO").rdkit_random_smiles(5, 42)
    #[pyo3(signature = (n, random_seed = 0, *, isomeric = true, kekule = false,
                        all_bonds_explicit = false, all_hs_explicit = false))]
    fn rdkit_random_smiles(
        &self,
        n: usize,
        random_seed: u32,
        isomeric: bool,
        kekule: bool,
        all_bonds_explicit: bool,
        all_hs_explicit: bool,
    ) -> PyResult<Vec<String>> {
        let params = chematic_smiles::RdkitSmilesParams {
            isomeric,
            kekule,
            canonical: false,
            all_bonds_explicit,
            all_hs_explicit,
            rooted_at_atom: None,
        };
        chematic_smiles::rdkit_random_smiles(&self.inner, n, random_seed, &params)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit 2026.03.1's ``Chem.MolFragmentToSmiles(m, atomsToUse,
    /// bondsToUse, ...)`` for ``m = Chem.MolFromSmiles(s)``. Atom and bond
    /// indices are RDKit's (see :meth:`rdkit_num_atoms`); ``bonds=None``
    /// uses every bond between two of the atoms. Raises ``ValueError``
    /// where RDKit raises (e.g. "neither end atom traversed") and for
    /// ``kekule=True`` on fragments with aromatic atoms::
    ///
    ///     chematic.from_smiles("OC(=O)c1ccccc1").rdkit_fragment_smiles([0, 1, 2])  # 'O=CO'
    #[pyo3(signature = (atoms, bonds = None, isomeric = true, kekule = false,
                        canonical = true, all_bonds_explicit = false,
                        all_hs_explicit = false, rooted_at_atom = None))]
    #[allow(clippy::too_many_arguments)]
    fn rdkit_fragment_smiles(
        &self,
        atoms: Vec<usize>,
        bonds: Option<Vec<usize>>,
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
        chematic_smiles::rdkit_fragment_smiles(&self.inner, &atoms, bonds.as_deref(), &params)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit 2026.03.1's ``Chem.GetDistanceMatrix(m, useBO, useAtomWts)``
    /// for ``m = Chem.MolFromSmiles(s)`` (``Chem.AddHs(m)`` for a molecule
    /// from :meth:`add_hydrogens`): rows in RDKit's atom order, ``1e8``
    /// between fragments.
    /// Inputs/prepared models above 1024 atoms raise ValueError before dense
    /// allocation; the input limit is checked before hydrogen removal.
    #[pyo3(signature = (use_bo = false, use_atom_wts = false))]
    fn rdkit_distance_matrix(&self, use_bo: bool, use_atom_wts: bool) -> PyResult<Vec<Vec<f64>>> {
        chematic_smiles::rdkit_distance_matrix(&self.inner, use_bo, use_atom_wts)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit 2026.03.1's ``Chem.Get3DDistanceMatrix(m, useAtomWts=...)``
    /// with the conformer ``coords`` (one ``[x, y, z]`` per atom of RDKit's
    /// molecule; for an :meth:`add_hydrogens` molecule its atom order, which
    /// is ``AddHs``'s).
    /// Inputs/prepared models above 1024 atoms raise ValueError before model
    /// preparation or dense allocation, including explicitly added H atoms.
    #[pyo3(signature = (coords, use_atom_wts = false))]
    fn rdkit_distance_matrix_3d(
        &self,
        coords: Vec<[f64; 3]>,
        use_atom_wts: bool,
    ) -> PyResult<Vec<Vec<f64>>> {
        chematic_smiles::rdkit_distance_matrix_3d(&self.inner, &coords, use_atom_wts)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit 2026.03.1's ``Chem.DetectChemistryProblems(Chem.MolFromSmiles(s))``
    /// on RDKit's sanitized molecule: ``[(type, atom_indices)]`` (see
    /// :func:`chematic.rdkit_detect_chemistry_problems` for unsanitized input).
    fn rdkit_chemistry_problems(&self) -> PyResult<Vec<(String, Vec<usize>)>> {
        let problems = chematic_smiles::rdkit_detect_chemistry_problems_sanitized(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(problems.into_iter().map(|p| (p.kind, p.atoms)).collect())
    }

    /// ``Chem.MolFromSmiles(s).GetNumAtoms()``: the number of atoms in
    /// RDKit's molecule (after its hydrogen removal).
    fn rdkit_num_atoms(&self) -> PyResult<usize> {
        chematic_smiles::rdkit_num_atoms(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit-compatible ``Chem.MolToCXSmarts(m)``.
    fn rdkit_cx_smarts(&self) -> PyResult<String> {
        chematic_smiles::rdkit_cx_smarts(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit-compatible Bemis-Murcko scaffold:
    /// ``Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))``.
    fn rdkit_murcko_scaffold(&self) -> PyResult<String> {
        chematic_smiles::rdkit_murcko_scaffold(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit's canonical tautomer as RDKit canonical SMILES:
    /// ``Chem.MolToSmiles(rdMolStandardize.TautomerEnumerator().Canonicalize(m))``
    /// for ``m = Chem.MolFromSmiles(smiles)`` (RDKit 2026.03.1, default
    /// parameters). Unlike :meth:`canonical_tautomer` (chematic's own
    /// algorithm) this is a port of RDKit's ``TautomerEnumerator``.
    ///
    /// Raises ``ValueError`` for molecules the RDKit port does not model and
    /// for molecules RDKit's sanitization rejects.
    ///
    ///     chematic.from_smiles("Oc1ccccn1").rdkit_canonical_tautomer()
    ///     # 'O=c1cccc[nH]1'
    fn rdkit_canonical_tautomer(&self) -> PyResult<String> {
        chematic_smiles::rdkit_canonical_tautomer(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// The tautomers RDKit enumerates, as sorted RDKit canonical SMILES:
    /// ``sorted(Chem.MolToSmiles(t) for t in
    /// rdMolStandardize.TautomerEnumerator().Enumerate(m))`` (RDKit
    /// 2026.03.1, default parameters: at most 1000 tautomers and 1000
    /// transforms).
    ///
    ///     chematic.from_smiles("Oc1ccccn1").rdkit_tautomers()
    ///     # ['O=C1CC=CC=N1', 'O=c1cccc[nH]1', 'Oc1ccccn1']
    fn rdkit_tautomers(&self) -> PyResult<Vec<String>> {
        let mut res = chematic_smiles::rdkit_enumerate_tautomers(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?
            .smiles;
        res.sort();
        Ok(res)
    }

    /// The status of RDKit's tautomer enumeration
    /// (``TautomerEnumeratorResult.status``): ``"Completed"``,
    /// ``"MaxTautomersReached"``, ``"MaxTransformsReached"`` or
    /// ``"Canceled"``.
    fn rdkit_tautomer_status(&self) -> PyResult<String> {
        let res = chematic_smiles::rdkit_enumerate_tautomers(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(format!("{:?}", res.status))
    }

    /// RDKit's tautomer score of this molecule:
    /// ``rdMolStandardize.TautomerEnumerator.ScoreTautomer(m)`` (aromatic
    /// rings, functional-group SMARTS terms and the H-on-P/S/Se/Te penalty).
    ///
    ///     chematic.from_smiles("O=c1cccc[nH]1").rdkit_tautomer_score()
    ///     # 102
    fn rdkit_tautomer_score(&self) -> PyResult<i32> {
        chematic_smiles::rdkit_tautomer_score(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit's ``rdMolHash.MolHash(m, function, useCXSmiles)`` (RDKit
    /// 2026.03.1) for ``m = Chem.MolFromSmiles(s)``, ``s`` the SMILES this
    /// molecule was read from. ``function`` names a member of
    /// ``rdMolHash.HashFunction`` (case-insensitive): ``"AnonymousGraph"``,
    /// ``"ElementGraph"``, ``"CanonicalSmiles"``, ``"MurckoScaffold"``,
    /// ``"ExtendedMurcko"``, ``"MolFormula"``, ``"AtomBondCounts"``,
    /// ``"DegreeVector"``, ``"Mesomer"``, ``"HetAtomTautomer"``,
    /// ``"HetAtomProtomer"``, ``"RedoxPair"``, ``"Regioisomer"``,
    /// ``"NetCharge"``, ``"SmallWorldIndexBR"``, ``"SmallWorldIndexBRL"``,
    /// ``"ArthorSubstructureOrder"``, ``"HetAtomTautomerv2"`` or
    /// ``"HetAtomProtomerv2"``.
    ///
    /// Raises ``ValueError`` for an unknown function name, for molecules the
    /// RDKit port does not model and for molecules RDKit's sanitization
    /// rejects.
    ///
    ///     chematic.from_smiles("Cc1ccccc1CC(=O)O").rdkit_mol_hash("ExtendedMurcko")
    ///     # '*c1ccccc1*'
    #[pyo3(signature = (function, use_cx_smiles=false))]
    fn rdkit_mol_hash(&self, function: &str, use_cx_smiles: bool) -> PyResult<String> {
        let f = chematic_smiles::RdkitHashFunction::from_name(function).ok_or_else(|| {
            let names: Vec<&str> = chematic_smiles::RdkitHashFunction::ALL
                .iter()
                .map(|f| f.name())
                .collect();
            PyValueError::new_err(format!(
                "unknown hash function {function:?}; expected one of {}",
                names.join(", ")
            ))
        })?;
        chematic_smiles::rdkit_mol_hash(&self.inner, f, use_cx_smiles)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// The RDKit canonical SMILES (``Chem.MolToSmiles(isomer)``) of the
    /// isomers :meth:`rdkit_stereoisomers` returns, sorted. Unlike
    /// re-reading them (``MolFromSmiles`` drops e.g. the chirality RDKit
    /// gives an aromatic ``[s+]([O-])``), these are RDKit's strings as
    /// written.
    #[pyo3(signature = (max_isomers = 1024))]
    fn rdkit_stereoisomer_smiles(&self, max_isomers: usize) -> PyResult<Vec<String>> {
        chematic_smiles::rdkit_stereoisomer_smiles(&self.inner, max_isomers)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit-compatible ``EnumerateStereoisomers.GetStereoisomerCount(m)``
    /// (default options): ``2 ** (number of flippable centres/bonds)``.
    fn rdkit_stereoisomer_count(&self) -> PyResult<u128> {
        chematic_smiles::rdkit_stereoisomer_count(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// RDKit-compatible ``Chem.FindMolChiralCenters(m, force=True,
    /// includeUnassigned=include_unassigned)`` (legacy stereo perception,
    /// RDKit 2026.03's default): ``[(atom_index, "R" | "S" | "?")]``.
    #[pyo3(signature = (include_unassigned = true))]
    fn rdkit_chiral_centers(&self, include_unassigned: bool) -> PyResult<Vec<(usize, String)>> {
        chematic_smiles::rdkit_chiral_centers(&self.inner, include_unassigned)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Return a copy with all implicit hydrogens made explicit.
    fn add_hydrogens(&self) -> Mol {
        Mol {
            inner: Arc::new(chematic_chem::add_hydrogens(&self.inner)),
            props: Default::default(),
        }
    }

    /// RDKit 2026.03.1's ``Chem.AddHs(Chem.MolFromSmiles(s))`` for a
    /// molecule read from ``s``: hydrogen atoms written in the SMILES
    /// (``[H]``) are first removed and stereo tags RDKit's reader drops are
    /// dropped (RDKit's molecule, in RDKit's atom order), then every
    /// hydrogen is added after the heavy atoms::
    ///
    ///     chematic.from_smiles("[H]C([H])([H])[H]").rdkit_add_hydrogens().rdkit_smiles
    ///     # '[H]C([H])([H])[H]'
    fn rdkit_add_hydrogens(&self) -> Mol {
        // Hydrogens on dummy atoms (`[*H]`), which chematic's own
        // `add_hydrogens` leaves implicit: the port's AddHs.
        if self
            .inner
            .atoms()
            .any(|(_, a)| a.wildcard && a.hydrogen_count.is_some_and(|h| h > 0))
            && let Some(m) = chematic_smiles::rdkit_added_hs_molecule(&self.inner)
        {
            return Mol {
                inner: Arc::new(m),
                props: Default::default(),
            };
        }
        let base = chematic_smiles::rdkit_parsed_molecule(&self.inner);
        let base = base.as_ref().unwrap_or(&self.inner);
        Mol {
            inner: Arc::new(chematic_chem::add_hydrogens(base)),
            props: Default::default(),
        }
    }

    /// RDKit 2026.03.1's ``Chem.MolFromSmiles(s)`` for a molecule read from
    /// ``s``: RDKit's atoms in RDKit's order (``[H]`` atoms removed), its
    /// aromaticity and hydrogen counts, and only the stereo tags RDKit's
    /// reader keeps. The molecule itself where the port cannot model it::
    ///
    ///     chematic.from_smiles("[H][C@](F)(Cl)Br").rdkit_parsed().smiles
    ///     # 'F[C@@H](Cl)Br'
    fn rdkit_parsed(&self) -> Mol {
        let inner = match chematic_smiles::rdkit_parsed_molecule(&self.inner) {
            Some(m) => Arc::new(m),
            None => Arc::clone(&self.inner),
        };
        Mol {
            inner,
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

    /// :meth:`connected_components` with the original atom indices.
    ///
    /// Returns a list of ``(fragment, source)`` pairs, where ``source[i]`` is
    /// the index in this molecule of fragment atom ``i``. Fragments are
    /// ordered by their lowest original atom index and keep the original
    /// relative atom order.
    ///
    ///     parts = chematic.from_smiles("O.CC").connected_components_with_atom_indices()
    ///     # [(Mol("O"), [0]), (Mol("CC"), [1, 2])]
    fn connected_components_with_atom_indices(&self) -> Vec<(Mol, Vec<usize>)> {
        self.inner
            .fragments_with_source_atoms()
            .into_iter()
            .map(|(fragment, source)| {
                (
                    Mol::bare(fragment),
                    source.into_iter().map(|a| a.0 as usize).collect(),
                )
            })
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
}
