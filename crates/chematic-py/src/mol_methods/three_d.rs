use super::*;

fn checked_transform_coords(
    mol: &chematic_core::Molecule,
    coords: &[[f64; 3]],
    atoms: &[usize],
) -> PyResult<chematic_3d::Coords3D> {
    if coords.len() != mol.atom_count() {
        return Err(PyValueError::new_err(format!(
            "expected {} coordinates, got {}",
            mol.atom_count(),
            coords.len()
        )));
    }
    if let Some(&atom) = atoms.iter().find(|&&atom| atom >= mol.atom_count()) {
        return Err(PyValueError::new_err(format!(
            "atom index {atom} is out of range for {} atoms",
            mol.atom_count()
        )));
    }
    if coords.iter().flatten().any(|value| !value.is_finite()) {
        return Err(PyValueError::new_err("coordinates must all be finite"));
    }
    Ok(flat_to_coords3d(coords))
}

fn transform_coords_to_python(coords: &chematic_3d::Coords3D) -> Vec<Vec<f64>> {
    coords.points.iter().map(|p| vec![p.x, p.y, p.z]).collect()
}

#[pymethods]
impl Mol {
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
    /// An ``"E"``/``"Z"`` entry's ``atom_idx`` is the double bond's first atom;
    /// such entries also carry ``"bond_idx"`` and ``"bond_atoms"`` (the
    /// double bond's index and its two atoms).
    /// Only assigned stereocenters / double bonds are returned.
    ///
    /// ``mode`` selects the CIP engine:
    ///
    /// - ``"legacy"`` (default) — the original fast engine, ~96.3% agreement with
    ///   RDKit's modern ``rdCIPLabeler`` oracle. Unchanged from every prior release.
    /// - ``"accurate"`` — a hierarchical-digraph engine for tetrahedral R/S
    ///   (~99.6% oracle-agreement on the representation-stable subset; see
    ///   ``docs/rfcs/cip_accurate_rfc.md``); E/Z from the same engine's substituent
    ///   ranking, allenes from legacy. Atoms it explicitly can't
    ///   resolve (a genuine tie, or exceeding its computation budget) are omitted
    ///   here and reported instead via :meth:`cip_stereo_unresolved` — never a
    ///   silently-guessed label. A phosphorus on an unsaturated ring
    ///   (cyclophosphazene) gets RDKit's CIPLabeler label, which flips between
    ///   the ring's Kekulé spellings in both libraries; its entry carries
    ///   ``"kekule_dependent": True``. A ring centre with three single bonds
    ///   and a lone pair (bridgehead amine, cyclic phosphine or sulfonium)
    ///   gets RDKit's label for the parsed SMILES spelling (RDKit inverts its
    ///   ``@``/``@@`` when one ring-closure digit is on the centre); its entry
    ///   carries ``"spelling_dependent": True``.
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
        // E/Z entries are keyed by the double bond's first atom; attach the
        // bond itself so callers need not guess which bond is meant (#634).
        let mut ez_bonds = chematic_chem::assign_ez_bonds_with_mode(
            &self.inner,
            if mode == "accurate" {
                chematic_chem::CipMode::Accurate
            } else {
                chematic_chem::CipMode::LegacyFast
            },
        );
        assignments
            .iter()
            .map(|(idx, code)| {
                let d = PyDict::new(py);
                d.set_item("atom_idx", idx.0 as usize)?;
                if matches!(code, CipCode::E | CipCode::Z)
                    && let Some(pos) = ez_bonds
                        .iter()
                        .position(|&(bidx, c)| c == *code && self.inner.bond(bidx).atom1 == *idx)
                {
                    let (bidx, _) = ez_bonds.remove(pos);
                    let bond = self.inner.bond(bidx);
                    d.set_item("bond_idx", bidx.0 as usize)?;
                    d.set_item("bond_atoms", (bond.atom1.0 as usize, bond.atom2.0 as usize))?;
                }
                let label = match code {
                    CipCode::R => "R",
                    CipCode::S => "S",
                    CipCode::E => "E",
                    CipCode::Z => "Z",
                    CipCode::LowerR => "r",
                    CipCode::LowerS => "s",
                };
                d.set_item("descriptor", label)?;
                if mode == "accurate"
                    && matches!(
                        code,
                        CipCode::R | CipCode::S | CipCode::LowerR | CipCode::LowerS
                    )
                    && chematic_chem::cip_label_depends_on_kekule_spelling(&self.inner, *idx)
                {
                    d.set_item("kekule_dependent", true)?;
                }
                if mode == "accurate"
                    && matches!(
                        code,
                        CipCode::R | CipCode::S | CipCode::LowerR | CipCode::LowerS
                    )
                    && chematic_chem::cip_label_depends_on_smiles_spelling(&self.inner, *idx)
                {
                    d.set_item("spelling_dependent", true)?;
                }
                Ok(d)
            })
            .collect()
    }

    /// Atoms :meth:`cip_stereo`\ (``mode="accurate"``) could not resolve a tetrahedral
    /// R/S for — list of ``{"atom_idx": int, "reason": str}`` dicts, ``reason`` is
    /// ``"tied"`` (a genuine CIP-rule tie, not a missing rule),
    /// ``"budget_exceeded"``, ``"oracle_unstable"``, or ``"lone_pair_center"`` (a
    /// stereo-tagged three-ligand centre on an aromatic ring, or one with a
    /// hydrogen ligand; saturated ring centres such as bridgehead amines are
    /// labelled, see :meth:`cip_stereo`).
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
                    chematic_chem::CipUnresolvedReason::OracleUnstable => "oracle_unstable",
                    chematic_chem::CipUnresolvedReason::LonePairCenter => "lone_pair_center",
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
        chematic_mol::stereo_depiction::with_stereo_depiction(&self.inner, |m, layout| {
            chematic_depict::render_svg_opts(m, layout, &opts)
        })
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
    ) -> PyResult<Option<f64>> {
        use chematic_core::AtomIdx;
        let c3d = checked_transform_coords(&self.inner, &coords, &[i, j, k, l])?;
        match chematic_3d::try_get_dihedral_deg(
            &c3d,
            &self.inner,
            AtomIdx(i as u32),
            AtomIdx(j as u32),
            AtomIdx(k as u32),
            AtomIdx(l as u32),
        ) {
            Ok(value) => Ok(Some(value)),
            Err(chematic_3d::MolTransformError::UndefinedDihedral) => Ok(None),
            Err(error) => Err(PyValueError::new_err(error.to_string())),
        }
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
    ) -> PyResult<Vec<Vec<f64>>> {
        use chematic_core::AtomIdx;
        let c3d = checked_transform_coords(&self.inner, &coords, &[i, j, k, l])?;
        let new_c3d = chematic_3d::try_set_dihedral(
            &c3d,
            &self.inner,
            AtomIdx(i as u32),
            AtomIdx(j as u32),
            AtomIdx(k as u32),
            AtomIdx(l as u32),
            angle_deg.to_radians(),
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(transform_coords_to_python(&new_c3d))
    }

    /// Measure bond length between atoms i and j (Å).
    fn get_bond_length(&self, coords: Vec<[f64; 3]>, i: usize, j: usize) -> PyResult<f64> {
        use chematic_core::AtomIdx;
        let c3d = checked_transform_coords(&self.inner, &coords, &[i, j])?;
        chematic_3d::try_get_bond_length(&c3d, &self.inner, AtomIdx(i as u32), AtomIdx(j as u32))
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Measure bond angle i–j–k in degrees (j is the central atom).
    fn get_bond_angle(&self, coords: Vec<[f64; 3]>, i: usize, j: usize, k: usize) -> PyResult<f64> {
        use chematic_core::AtomIdx;
        let c3d = checked_transform_coords(&self.inner, &coords, &[i, j, k])?;
        chematic_3d::try_get_bond_angle_deg(
            &c3d,
            &self.inner,
            AtomIdx(i as u32),
            AtomIdx(j as u32),
            AtomIdx(k as u32),
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Set bond length i-j to ``length`` Å.
    ///
    /// Atom i stays fixed; atom j and its side of the molecular graph move as
    /// one rigid fragment. Ring bonds are rejected, matching RDKit's
    /// ``rdMolTransforms.SetBondLength`` contract.
    fn set_bond_length(
        &self,
        coords: Vec<[f64; 3]>,
        i: usize,
        j: usize,
        length: f64,
    ) -> PyResult<Vec<Vec<f64>>> {
        use chematic_core::AtomIdx;
        let c3d = checked_transform_coords(&self.inner, &coords, &[i, j])?;
        let changed = chematic_3d::set_bond_length(
            &c3d,
            &self.inner,
            AtomIdx(i as u32),
            AtomIdx(j as u32),
            length,
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(transform_coords_to_python(&changed))
    }

    /// Set bond angle i-j-k to ``angle_deg`` degrees.
    ///
    /// Atoms i and j stay fixed; atom k and its side of the molecular graph
    /// rotate as one rigid fragment. Two ring bonds are rejected, matching
    /// RDKit's ``rdMolTransforms.SetAngleDeg`` contract.
    fn set_bond_angle(
        &self,
        coords: Vec<[f64; 3]>,
        i: usize,
        j: usize,
        k: usize,
        angle_deg: f64,
    ) -> PyResult<Vec<Vec<f64>>> {
        use chematic_core::AtomIdx;
        let c3d = checked_transform_coords(&self.inner, &coords, &[i, j, k])?;
        let changed = chematic_3d::set_bond_angle(
            &c3d,
            &self.inner,
            AtomIdx(i as u32),
            AtomIdx(j as u32),
            AtomIdx(k as u32),
            angle_deg.to_radians(),
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(transform_coords_to_python(&changed))
    }

    /// Compute the RDKit-style principal-axis transform for a conformer.
    ///
    /// The result is a row-major 4x4 homogeneous matrix. Eigenvector signs
    /// are solver-dependent, so equivalent frames may differ by a
    /// 180-degree proper rotation across platforms.
    #[pyo3(signature = (coords, center = None, normalize_covar = false, ignore_hs = true))]
    fn compute_canonical_transform(
        &self,
        coords: Vec<[f64; 3]>,
        center: Option<[f64; 3]>,
        normalize_covar: bool,
        ignore_hs: bool,
    ) -> PyResult<Vec<Vec<f64>>> {
        let c3d = checked_transform_coords(&self.inner, &coords, &[])?;
        let transform = chematic_3d::compute_canonical_transform(
            &c3d,
            &self.inner,
            center,
            normalize_covar,
            ignore_hs,
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(transform.matrix.into_iter().map(Vec::from).collect())
    }

    /// Return a conformer centered and aligned to its principal axes.
    #[pyo3(signature = (coords, center = None, normalize_covar = false, ignore_hs = true))]
    fn canonicalize_conformer(
        &self,
        coords: Vec<[f64; 3]>,
        center: Option<[f64; 3]>,
        normalize_covar: bool,
        ignore_hs: bool,
    ) -> PyResult<Vec<Vec<f64>>> {
        let c3d = checked_transform_coords(&self.inner, &coords, &[])?;
        let canonical = chematic_3d::canonicalize_conformer(
            &c3d,
            &self.inner,
            center,
            normalize_covar,
            ignore_hs,
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(transform_coords_to_python(&canonical))
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
    /// notes on Morgan fingerprints for the same caveat); :meth:`rdkit_avalon_fp`
    /// is the RDKit-identical fingerprint.
    fn avalon_fp(&self) -> Vec<u8> {
        bitvec2048_to_bytes(&chematic_fp::avalon_fp(&self.inner))
    }

    /// RDKit's Avalon fingerprint as bytes (``n_bits / 8`` bytes, bit ``i`` at
    /// byte ``i // 8``, mask ``1 << (i % 8)``).
    ///
    /// Bit-identical to ``rdkit.Avalon.pyAvalonTools.GetAvalonFP(mol,
    /// nBits=n_bits)`` (default ``bitFlags``, ``isQuery=False``) for the RDKit
    /// molecule ``Chem.MolFromSmiles`` builds from the SMILES this molecule was
    /// read from: a port of the Avalon toolkit's fingerprint code (version
    /// 2.0.5-pre.3, the one RDKit 2026.03 builds) applied to the connection
    /// table RDKit's ``MolToMolBlock`` writes.
    ///
    /// Raises ``ValueError`` for molecules RDKit would reject or that the RDKit
    /// model does not cover.
    #[pyo3(signature = (n_bits = 2048))]
    fn rdkit_avalon_fp(&self, n_bits: usize) -> PyResult<Vec<u8>> {
        chematic_fp::rdkit_avalon_fp(&self.inner, n_bits)
            .map_err(|e| PyValueError::new_err(e.to_string()))
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

    /// Reaction fingerprint — ``(reactant_fp, product_fp, formed_fp, broken_fp,
    /// combined_fp)`` as bytes.
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
        d.set_item("formed_fp", bitvec2048_to_bytes(&rfp.formed_fp))?;
        d.set_item("broken_fp", bitvec2048_to_bytes(&rfp.broken_fp))?;
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
}
