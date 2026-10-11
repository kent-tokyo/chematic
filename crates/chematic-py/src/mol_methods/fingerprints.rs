use super::*;

#[pymethods]
impl Mol {
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

    /// RDKit's hashed atom-pair count fingerprint
    /// (``rdMolDescriptors.GetHashedAtomPairFingerprint(m, nBits=nbits)``)
    /// as a sorted list of ``(bucket, count)`` nonzero elements.
    #[pyo3(signature = (nbits = 2048))]
    fn rdkit_atom_pair_counts(&self, nbits: u32) -> PyResult<Vec<(u32, u32)>> {
        if nbits == 0 {
            return Err(PyValueError::new_err("nbits must be positive"));
        }
        Ok(chematic_fp::rdkit_atom_pair_counts(&self.inner, nbits))
    }

    /// RDKit's hashed topological-torsion count fingerprint
    /// (``rdMolDescriptors.GetHashedTopologicalTorsionFingerprint(m, nBits=nbits)``)
    /// as a sorted list of ``(bucket, count)`` nonzero elements.
    #[pyo3(signature = (nbits = 2048))]
    fn rdkit_torsion_counts(&self, nbits: u32) -> PyResult<Vec<(u32, u32)>> {
        if nbits == 0 {
            return Err(PyValueError::new_err("nbits must be positive"));
        }
        Ok(chematic_fp::rdkit_torsion_counts(&self.inner, nbits))
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
        // Same bits and errors as `rdkit_morgan_ecfp4_experimental`, without
        // building the provenance maps this method never returns.
        let fingerprint = chematic_fp::rdkit_morgan_ecfp4_bitvec(&self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(bitvec2048_to_bytes(&fingerprint))
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
    /// ``include_chirality`` is RDKit's ``includeChirality``: tetrahedral
    /// centres and E/Z double bonds as RDKit's (legacy) stereo perception
    /// sees them. ``count_simulation`` is RDKit's ``countSimulation`` with
    /// the default ``countBounds`` ``[1, 2, 4, 8]``.
    #[pyo3(signature = (radius = 2, nbits = 2048, include_chirality = false, count_simulation = false))]
    fn rdkit_ecfp_config(
        &self,
        radius: u32,
        nbits: usize,
        include_chirality: bool,
        count_simulation: bool,
    ) -> PyResult<Vec<u8>> {
        let config = python_rdkit_morgan_config(radius, nbits, include_chirality)?;
        if count_simulation {
            let fp = chematic_fp::rdkit_morgan_count_simulation(&self.inner, &config)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            return Ok(bitvecn_to_bytes(&fp));
        }
        let fingerprint = chematic_fp::rdkit_morgan_bitvec(&self.inner, &config)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(bitvecn_to_bytes(&fingerprint))
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

    /// RDKit 2026.03.1's ``rdFingerprintGenerator.GetMorganGenerator(
    /// radius=radius, includeChirality=include_chirality)
    /// .GetSparseCountFingerprint(m).GetNonzeroElements()``: unfolded
    /// 32-bit Morgan identifiers and their counts (radius 0-3).
    #[pyo3(signature = (radius = 2, include_chirality = false))]
    fn rdkit_morgan_sparse_counts(
        &self,
        radius: u32,
        include_chirality: bool,
    ) -> PyResult<std::collections::BTreeMap<u32, u32>> {
        let config = python_rdkit_morgan_config(radius, 2048, include_chirality)?;
        let result = chematic_fp::rdkit_morgan_fingerprint(&self.inner, &config)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(result.sparse_counts.into_iter().collect())
    }

    /// RDKit 2026.03.1's ``rdFingerprintGenerator.GetAtomPairGenerator()
    /// .GetSparseCountFingerprint(m).GetNonzeroElements()``: unhashed
    /// atom-pair codes and their counts.
    fn rdkit_atom_pair_sparse_counts(&self) -> std::collections::BTreeMap<u32, u32> {
        chematic_fp::rdkit_atom_pair_sparse_counts(&self.inner)
            .into_iter()
            .collect()
    }

    /// RDKit 2026.03.1's ``rdFingerprintGenerator.GetTopologicalTorsionGenerator()
    /// .GetSparseCountFingerprint(m).GetNonzeroElements()``: unhashed
    /// torsion codes and their counts.
    fn rdkit_torsion_sparse_counts(&self) -> std::collections::BTreeMap<u64, u32> {
        chematic_fp::rdkit_torsion_sparse_counts(&self.inner)
            .into_iter()
            .collect()
    }

    /// RDKit 2026.03.1's legacy
    /// ``rdMolDescriptors.GetTopologicalTorsionFingerprint(m).GetNonzeroElements()``
    /// (unhashed 64-bit torsion codes).
    fn rdkit_legacy_torsion_counts(&self) -> std::collections::BTreeMap<i64, u32> {
        chematic_fp::rdkit_legacy_torsion_counts(&self.inner)
            .into_iter()
            .collect()
    }

    /// RDKit's Morgan ``bitInfo`` for ``n_bits``-bit fingerprints of radius
    /// ``radius`` (``rdFingerprintGenerator.GetMorganGenerator(radius=radius,
    /// fpSize=n_bits, includeChirality=include_chirality)`` with
    /// ``AdditionalOutput.GetBitInfoMap()``): ``{bit: [(atom, radius), ...]}``
    /// with each list sorted.
    #[pyo3(signature = (radius = 2, n_bits = 2048, include_chirality = false))]
    fn rdkit_morgan_bit_info(
        &self,
        radius: u32,
        n_bits: usize,
        include_chirality: bool,
    ) -> PyResult<std::collections::BTreeMap<usize, Vec<(u32, u32)>>> {
        let config = python_rdkit_morgan_config(radius, n_bits, include_chirality)?;
        let result = chematic_fp::rdkit_morgan_fingerprint(&self.inner, &config)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(result
            .folded_bit_info
            .into_iter()
            .map(|(bit, mut env)| {
                env.sort_unstable();
                (bit, env)
            })
            .collect())
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
        chematic_mol::stereo_depiction::with_stereo_depiction(
            &self.inner,
            chematic_depict::render_svg,
        )
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
        let data =
            chematic_mol::stereo_depiction::with_stereo_depiction(&self.inner, |m, layout| {
                let coords: Vec<(f64, f64)> = layout.coords.iter().map(|p| (p.x, p.y)).collect();
                chematic_depict::depict_data_with_coords(m, &coords)
            });

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
        chematic_mol::stereo_depiction::with_stereo_depiction(
            &self.inner,
            chematic_depict::render_eps,
        )
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
        chematic_depict::svg_to_pdf(&chematic_mol::stereo_depiction::with_stereo_depiction(
            &self.inner,
            chematic_depict::render_svg,
        ))
    }

    /// Jupyter Notebook / JupyterLab の自動描画フック。
    ///
    /// セルに ``mol`` と書くだけで 2D 構造が表示される。手動で
    /// ``IPython.display.SVG(mol.svg())`` と書く必要はない。
    fn _repr_svg_(&self) -> String {
        chematic_mol::stereo_depiction::with_stereo_depiction(
            &self.inner,
            chematic_depict::render_svg,
        )
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
        let query = crate::misc::cached_smarts(smarts).map_err(|e| {
            crate::errors::malformed("smarts", format!("invalid SMARTS '{smarts}': {e}"))
        })?;
        // Stop at the first embedding instead of enumerating every match — an
        // existence check doesn't need the full match set or the dedup pass.
        let config = chematic_smarts::MatchConfig {
            max_matches: Some(1),
            uniquify: false,
            ..chematic_smarts::MatchConfig::default()
        };
        Ok(chematic_smarts::has_match_perceived(
            &query,
            &self.inner,
            &config,
        ))
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
        let query = crate::misc::cached_smarts(smarts).map_err(|e| {
            crate::errors::malformed("smarts", format!("invalid SMARTS '{smarts}': {e}"))
        })?;
        Ok(chematic_smarts::find_match_atom_sets_perceived(
            &query,
            &self.inner,
            &chematic_smarts::MatchConfig::default(),
        ))
    }

    /// Opt-in, RDKit 2026.03.6-style SMARTS match sets with typed refusals.
    ///
    /// The native :meth:`find_matches` contract is unchanged. In particular,
    /// ambiguous ring-count systems and bounded searches are never presented
    /// as an empty match set. This mode is not general RDKit SMARTS parity.
    ///
    /// The default profile counts ``[R<n>]`` rings in RDKit 2026.03's
    /// symmetrized SSSR, using the exact port of RDKit's atom-order-dependent
    /// ring list where that list depends on atom order; it refuses as
    /// ``ring_model_ambiguous`` only where that port declines.
    /// ``profile="2026.09.1"`` counts ``[R<n>]`` rings in RDKit 2026.09.1's
    /// ring list (all relevant cycles) instead.
    #[pyo3(signature = (smarts, profile = "2026.03.6"))]
    fn find_matches_rdkit_parity<'py>(
        &self,
        smarts: &str,
        profile: &str,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, PyDict>> {
        use chematic_smarts::{RdkitParityError, RdkitRingCountModel};

        let ring_count_model = match profile {
            "2026.03.6" => RdkitRingCountModel::SymmetrizedSssr,
            "2026.09.1" => RdkitRingCountModel::RelevantCycles,
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown profile '{other}': expected '2026.03.6' or '2026.09.1'"
                )));
            }
        };

        let query = crate::misc::cached_smarts(smarts).map_err(|e| {
            crate::errors::malformed("smarts", format!("invalid SMARTS '{smarts}': {e}"))
        })?;
        let result = PyDict::new(py);
        let config = chematic_smarts::RdkitParityConfig {
            use_rdkit_parity_aromaticity: true,
            ring_count_model,
            ..chematic_smarts::RdkitParityConfig::default()
        };
        // RDKit matches its hydrogen-suppressed molecule (`MolFromSmiles`
        // removes most hydrogen graph atoms), in its atom numbering.
        let suppressed = chematic_smiles::rdkit_hydrogen_suppressed(&self.inner);
        let target = suppressed.as_ref().unwrap_or(&self.inner);
        match chematic_smarts::find_match_atom_sets_rdkit_parity(&query, target, &config) {
            Ok((atom_sets, false)) => {
                result.set_item("status", "ok")?;
                result.set_item("reason", py.None())?;
                result.set_item("matches", atom_sets)?;
            }
            Ok((_, true)) => {
                result.set_item("status", "typed_refusal")?;
                result.set_item("reason", "matcher_budget_exceeded")?;
                result.set_item("matches", py.None())?;
            }
            Err(error) => {
                let (status, reason) = match &error {
                    RdkitParityError::RingModelAmbiguous { .. } => {
                        ("typed_unsupported", "ring_model_ambiguous")
                    }
                    RdkitParityError::RingModelBudgetExceeded { .. } => {
                        ("typed_refusal", "ring_model_budget_exceeded")
                    }
                    RdkitParityError::Aromaticity(_) => ("typed_refusal", "aromaticity_failure"),
                };
                result.set_item("status", status)?;
                result.set_item("reason", reason)?;
                result.set_item("matches", py.None())?;
            }
        }
        Ok(result)
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
        chematic_mol::stereo_depiction::with_stereo_depiction(&self.inner, |m, layout| {
            chematic_depict::render_svg_opts(m, layout, &opts)
        })
    }

    /// 2D SVG depiction with atoms coloured by LogP contribution.
    ///
    /// Positive (lipophilic) atoms → blue, negative (hydrophilic) atoms → red,
    /// zero-contribution atoms → white.  Uses :meth:`logp_per_atom` as weights.
    fn logp_map_svg(&self) -> String {
        let weights = chematic_chem::logp_crippen_per_atom(&self.inner);
        chematic_mol::stereo_depiction::with_stereo_depiction(&self.inner, |m, layout| {
            chematic_depict::render_svg_opts(
                m,
                layout,
                &chematic_depict::similarity_map_options(
                    &self.inner,
                    &weights,
                    &chematic_depict::RenderOptions::default(),
                ),
            )
        })
    }

    /// 2D SVG depiction with atoms coloured by TPSA contribution.
    ///
    /// Atoms contributing to TPSA (N, O, S, P) are coloured blue; zero-contribution
    /// atoms (C, halogens, …) remain white.  Uses :meth:`tpsa_per_atom` as weights.
    fn tpsa_map_svg(&self) -> String {
        let weights = chematic_chem::tpsa_per_atom(&self.inner);
        chematic_mol::stereo_depiction::with_stereo_depiction(&self.inner, |m, layout| {
            chematic_depict::render_svg_opts(
                m,
                layout,
                &chematic_depict::similarity_map_options(
                    &self.inner,
                    &weights,
                    &chematic_depict::RenderOptions::default(),
                ),
            )
        })
    }

    /// 2D SVG depiction with atoms coloured by custom weights.
    ///
    /// ``weights``: list of floats, one per heavy atom (length = :attr:`heavy_atoms`).
    /// Positive → blue, negative → red, zero → white.
    fn similarity_map_svg(&self, weights: Vec<f64>) -> String {
        chematic_mol::stereo_depiction::with_stereo_depiction(&self.inner, |m, layout| {
            chematic_depict::render_svg_opts(
                m,
                layout,
                &chematic_depict::similarity_map_options(
                    &self.inner,
                    &weights,
                    &chematic_depict::RenderOptions::default(),
                ),
            )
        })
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
}
