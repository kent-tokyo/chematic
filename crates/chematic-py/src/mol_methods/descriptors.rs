use super::*;

#[pymethods]
impl Mol {
    // -----------------------------------------------------------------------
    // Core physicochemical descriptors
    // -----------------------------------------------------------------------

    /// Average molecular weight (Da).
    #[getter]
    fn mw(&self) -> f64 {
        chematic_chem::molecular_weight(&self.inner)
    }

    /// RDKit-compatible average molecular weight (Da), including explicit
    /// isotope labels and RDKit's mass-number fallback for unknown nuclides.
    #[getter]
    fn rdkit_mw(&self) -> f64 {
        chematic_chem::rdkit_molecular_weight(&self.inner)
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

    /// RDKit-default TPSA (Å²): N and O contributions only, matching
    /// ``rdMolDescriptors.CalcTPSA(mol)`` / ``Descriptors.TPSA(mol)``
    /// (``includeSandP=False``). :attr:`tpsa` also includes S and P.
    #[getter]
    fn rdkit_tpsa(&self) -> f64 {
        chematic_chem::rdkit_tpsa(&self.inner)
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

    /// RDKit-compatible hydrogen-bond acceptor count.
    #[getter]
    fn rdkit_hba(&self) -> usize {
        chematic_chem::rdkit_hba_count(&self.inner)
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

    /// Number of rings as RDKit's ``rdMolDescriptors.CalcNumRings`` counts
    /// them (its symmetrized SSSR); larger than :attr:`ring_count` for cages
    /// and bridged systems.
    #[getter]
    #[inline(always)]
    fn num_rings(&self) -> usize {
        self.num_rings_cache
            .get_or_compute(|| chematic_chem::rdkit_num_rings(&self.inner))
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

    /// RDKit-compatible aromatic ring count using the opt-in aromaticity model.
    #[getter]
    fn rdkit_aromatic_ring_count(&self) -> usize {
        chematic_chem::rdkit_aromatic_ring_count(&self.inner)
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
            num_rings_cache: Default::default(),
            inner: Arc::new(new_mol),
            props: Default::default(),
        })
    }

    /// Number of potential tetrahedral stereocenters (specified or unspecified).
    #[getter]
    fn num_stereocenters(&self) -> usize {
        chematic_chem::num_stereocenters(&self.inner)
    }

    /// Atom indices of potential tetrahedral stereocenters.
    #[getter]
    fn potential_stereocenter_indices(&self) -> Vec<usize> {
        chematic_chem::potential_stereocenter_indices(&self.inner)
            .into_iter()
            .map(|idx| idx.0 as usize)
            .collect()
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
        molecule_descriptors_dict(py, &self.inner)
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

        let lipophilicity = lipophilicity_label(logp);
        let bbb = if admet.bbb_passes {
            "✓ CNS penetrant"
        } else {
            "✗ CNS non-penetrant"
        };
        let caco2 = permeability_label(admet.caco2);
        let herg = risk_label(admet.herg_risk);
        let cyp = risk_label(admet.cyp3a4_risk);

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
}

fn lipophilicity_label(logp: f64) -> &'static str {
    if logp < 0.0 {
        "hydrophilic"
    } else if logp < 2.0 {
        "mildly lipophilic"
    } else if logp < 4.0 {
        "moderately lipophilic"
    } else {
        "highly lipophilic"
    }
}

fn permeability_label(value: f64) -> &'static str {
    if value > -5.5 {
        "High (well absorbed)"
    } else if value > -7.0 {
        "Moderate"
    } else {
        "Low (poor absorption)"
    }
}

fn risk_label(value: f64) -> &'static str {
    if value < 0.3 {
        "Low"
    } else if value < 0.6 {
        "Moderate"
    } else {
        "High ⚠"
    }
}

struct MoleculeDescriptorParts<'a> {
    molecule: &'a chematic_core::Molecule,
    rings: chematic_chem::RingBundle,
    mw: f64,
    logp: f64,
    mr: f64,
    tpsa: f64,
    hbd: usize,
}

impl<'a> MoleculeDescriptorParts<'a> {
    fn new(molecule: &'a chematic_core::Molecule) -> Self {
        let rings = chematic_chem::ring_bundle(molecule);
        let mw = chematic_chem::molecular_weight(molecule);
        let (logp, mr) = chematic_chem::logp_and_mr(molecule);
        let tpsa = chematic_chem::tpsa(molecule);
        let hbd = chematic_chem::hbd_count(molecule);
        Self {
            molecule,
            rings,
            mw,
            logp,
            mr,
            tpsa,
            hbd,
        }
    }
}

fn molecule_descriptors_dict<'py>(
    py: Python<'py>,
    molecule: &chematic_core::Molecule,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let parts = MoleculeDescriptorParts::new(molecule);
    set_primary_descriptors(&d, &parts)?;
    set_topology_descriptors(&d, &parts)?;
    set_filter_descriptors(&d, &parts)?;
    set_vsa_descriptors(&d, molecule)?;
    set_extended_descriptors(&d, molecule)?;
    Ok(d)
}

fn set_primary_descriptors(
    d: &Bound<'_, PyDict>,
    parts: &MoleculeDescriptorParts<'_>,
) -> PyResult<()> {
    let m = parts.molecule;
    let rings = &parts.rings;
    d.set_item("mw", parts.mw)?;
    d.set_item("exact_mass", chematic_chem::exact_mass(m))?;
    d.set_item("tpsa", parts.tpsa)?;
    d.set_item("logp", parts.logp)?;
    d.set_item("molar_refractivity", parts.mr)?;
    d.set_item("hbd", parts.hbd)?;
    d.set_item("hba", rings.hba_count)?;
    d.set_item("rotatable_bonds", rings.rotatable_bond_count)?;
    d.set_item("heavy_atoms", chematic_chem::heavy_atom_count(m))?;
    d.set_item("ring_count", rings.ring_count)?;
    d.set_item("ring_system_count", rings.ring_system_count)?;
    d.set_item("aromatic_ring_count", rings.aromatic_ring_count)?;
    d.set_item("hba_lipinski", chematic_chem::hba_count_lipinski(m))?;
    d.set_item("fraction_rotatable_bonds", rings.fraction_rotatable_bonds)?;
    d.set_item("num_heteroatoms", chematic_chem::num_heteroatoms(m))?;
    d.set_item("num_stereocenters", chematic_chem::num_stereocenters(m))?;
    d.set_item("num_spiro_atoms", rings.num_spiro_atoms)?;
    d.set_item("num_bridgehead_atoms", rings.num_bridgehead_atoms)?;
    d.set_item("fsp3", chematic_chem::fsp3(m))?;
    d.set_item("qed", chematic_chem::qed_with_bundle(m, rings))?;
    d.set_item("sa_score", chematic_chem::sa_score_with_bundle(m, rings))?;
    d.set_item("formal_charge", chematic_chem::formal_charge_sum(m))?;
    d.set_item("labute_asa", chematic_chem::labute_asa(m))?;
    d.set_item("bertz_ct", chematic_chem::bertz_ct(m))?;
    Ok(())
}

fn set_topology_descriptors(
    d: &Bound<'_, PyDict>,
    parts: &MoleculeDescriptorParts<'_>,
) -> PyResult<()> {
    let m = parts.molecule;
    let rings = &parts.rings;
    d.set_item("wiener_index", chematic_chem::wiener_index(m))?;
    d.set_item("schultz_mti", chematic_chem::schultz_mti(m))?;
    d.set_item("gutman_mti", chematic_chem::gutman_mti(m))?;
    d.set_item("vabc", chematic_chem::vabc(m))?;
    d.set_item("gravitational_index", chematic_chem::gravitational_index(m))?;
    d.set_item("zagreb_m2", chematic_chem::zagreb_index_m2(m))?;
    d.set_item("kappa1", chematic_chem::kappa1(m))?;
    d.set_item("kappa2", chematic_chem::kappa2(m))?;
    d.set_item("kappa3", chematic_chem::kappa3(m))?;
    let chi = chematic_chem::chi_all(m);
    for (name, value) in [
        ("chi0", chi.0),
        ("chi1", chi.1),
        ("chi2", chi.2),
        ("chi3", chi.3),
        ("chi4", chi.4),
        ("chi0v", chi.5),
        ("chi1v", chi.6),
        ("chi2v", chi.7),
        ("chi3v", chi.8),
        ("chi4v", chi.9),
    ] {
        d.set_item(name, value)?;
    }
    d.set_item("num_aromatic_heterocycles", rings.num_aromatic_heterocycles)?;
    d.set_item(
        "num_aliphatic_heterocycles",
        rings.num_aliphatic_heterocycles,
    )?;
    d.set_item(
        "num_saturated_heterocycles",
        rings.num_saturated_heterocycles,
    )?;
    d.set_item("num_saturated_rings", rings.num_saturated_rings)?;
    d.set_item("num_aliphatic_rings", rings.num_aliphatic_rings)?;
    d.set_item(
        "num_unspecified_stereocenters",
        chematic_chem::num_unspecified_stereocenters(m),
    )?;
    let estate = chematic_chem::estate_all(m);
    d.set_item("sum_estate", estate.0)?;
    d.set_item("max_estate", estate.1)?;
    d.set_item("min_estate", estate.2)?;
    Ok(())
}

fn set_filter_descriptors(
    d: &Bound<'_, PyDict>,
    parts: &MoleculeDescriptorParts<'_>,
) -> PyResult<()> {
    set_drug_filter_descriptors(d, parts)?;
    set_admet_descriptors(d, parts)
}

fn set_drug_filter_descriptors(
    d: &Bound<'_, PyDict>,
    parts: &MoleculeDescriptorParts<'_>,
) -> PyResult<()> {
    let m = parts.molecule;
    let rings = &parts.rings;
    d.set_item(
        "lipinski_passes",
        parts.mw <= 500.0 && parts.hbd <= 5 && rings.hba_count <= 10 && parts.logp <= 5.0,
    )?;
    d.set_item(
        "veber_passes",
        rings.rotatable_bond_count <= 10 && parts.tpsa <= 140.0,
    )?;
    d.set_item("egan_passes", parts.tpsa <= 131.6 && parts.logp <= 5.88)?;
    d.set_item("ghose_passes", chematic_chem::ghose_passes(m))?;
    d.set_item("reos_passes", chematic_chem::reos_passes(m))?;
    d.set_item("pains_passes", chematic_chem::pains_passes(m))?;
    d.set_item("ro3_passes", chematic_chem::ro3_passes(m))?;
    d.set_item("lead_like_passes", chematic_chem::lead_like_passes(m))?;
    d.set_item("pfizer_3_75_passes", chematic_chem::pfizer_3_75_passes(m))?;
    let pka_base = chematic_chem::pka_base(m).unwrap_or(0.0);
    d.set_item(
        "cns_mpo_score",
        chematic_chem::cns_mpo_from_parts(m, parts.logp, parts.tpsa, parts.mw, parts.hbd, pka_base),
    )?;
    d.set_item("mcf_passes", chematic_chem::mcf_passes(m))?;
    Ok(())
}

fn set_admet_descriptors(
    d: &Bound<'_, PyDict>,
    parts: &MoleculeDescriptorParts<'_>,
) -> PyResult<()> {
    let m = parts.molecule;
    let rings = &parts.rings;
    d.set_item(
        "bbb_score",
        chematic_chem::bbb_score_from_parts(parts.tpsa, parts.logp),
    )?;
    d.set_item(
        "bbb_passes",
        parts.tpsa < 90.0 && parts.mw < 400.0 && parts.hbd <= 3,
    )?;
    d.set_item(
        "caco2",
        chematic_chem::caco2_precomputed(parts.tpsa, parts.logp),
    )?;
    d.set_item(
        "herg_risk",
        chematic_chem::herg_risk_precomputed(m, parts.logp, parts.mw),
    )?;
    d.set_item(
        "cyp3a4_risk",
        chematic_chem::cyp3a4_precomputed(
            parts.mw,
            parts.logp,
            rings.num_aromatic_heterocycles,
            rings.hba_count,
        ),
    )?;
    let egg = chematic_chem::boiled_egg_from(parts.logp, parts.tpsa);
    d.set_item("gi_absorbed", egg.gi_absorbed)?;
    d.set_item("bbb_penetrant", egg.bbb_penetrant)?;
    d.set_item("pka_acid", chematic_chem::pka_acid(m))?;
    d.set_item("pka_base", chematic_chem::pka_base(m))?;
    Ok(())
}

fn set_vsa_descriptors(d: &Bound<'_, PyDict>, molecule: &chematic_core::Molecule) -> PyResult<()> {
    for (prefix, values) in [
        ("SlogP_VSA", chematic_chem::slogp_vsa(molecule)),
        ("SMR_VSA", chematic_chem::smr_vsa(molecule)),
        ("PEOE_VSA", chematic_chem::peoe_vsa(molecule)),
        ("EState_VSA", chematic_chem::estate_vsa(molecule)),
    ] {
        for (index, value) in values.into_iter().enumerate() {
            d.set_item(format!("{prefix}{}", index + 1), value)?;
        }
    }
    Ok(())
}

fn set_extended_descriptors(
    d: &Bound<'_, PyDict>,
    molecule: &chematic_core::Molecule,
) -> PyResult<()> {
    d.set_item(
        "num_valence_electrons",
        chematic_chem::num_valence_electrons(molecule),
    )?;
    d.set_item("hall_kier_alpha", chematic_chem::hall_kier_alpha(molecule))?;
    let bcut = chematic_chem::bcut2d(molecule);
    for (name, value) in [
        ("bcut2d_chghi", bcut.chghi),
        ("bcut2d_chglo", bcut.chglo),
        ("bcut2d_logphi", bcut.logphi),
        ("bcut2d_logplo", bcut.logplo),
        ("bcut2d_mrhi", bcut.mrhi),
        ("bcut2d_mrlo", bcut.mrlo),
        ("bcut2d_mwhi", bcut.mwhi),
        ("bcut2d_mwlo", bcut.mwlo),
    ] {
        d.set_item(name, value)?;
    }
    let carbon = chematic_chem::carbon_types(molecule);
    for (name, value) in [
        ("c1sp1", carbon.c1sp1),
        ("c2sp1", carbon.c2sp1),
        ("c1sp2", carbon.c1sp2),
        ("c2sp2", carbon.c2sp2),
        ("c3sp2", carbon.c3sp2),
        ("c1sp3", carbon.c1sp3),
        ("c2sp3", carbon.c2sp3),
        ("c3sp3", carbon.c3sp3),
    ] {
        d.set_item(name, value)?;
    }
    d.set_item("balaban_j", chematic_chem::balaban_j(molecule))?;
    d.set_item("num_amide_bonds", chematic_chem::num_amide_bonds(molecule))?;
    d.set_item("num_ester_bonds", chematic_chem::num_ester_bonds(molecule))?;
    set_element_counts(d, molecule)?;
    d.set_item("esol", chematic_chem::esol_solubility(molecule))?;
    d.set_item("logd_7_4", chematic_chem::logd_simple(molecule, 7.4))?;
    d.set_item("xlogp3", chematic_chem::xlogp3(molecule))?;
    d.set_item("drug_score", chematic_chem::drug_score(molecule))?;
    for (index, value) in chematic_chem::mqn(molecule).iter().enumerate() {
        d.set_item(format!("MQN{}", index + 1), value)?;
    }
    Ok(())
}

fn set_element_counts(d: &Bound<'_, PyDict>, molecule: &chematic_core::Molecule) -> PyResult<()> {
    d.set_item("num_carbons", chematic_chem::num_carbons(molecule))?;
    d.set_item("num_nitrogens", chematic_chem::num_nitrogens(molecule))?;
    d.set_item("num_oxygens", chematic_chem::num_oxygens(molecule))?;
    d.set_item("num_sulfurs", chematic_chem::num_sulfurs(molecule))?;
    d.set_item("num_phosphorus", chematic_chem::num_phosphorus(molecule))?;
    d.set_item("num_fluorines", chematic_chem::num_fluorines(molecule))?;
    d.set_item("num_chlorines", chematic_chem::num_chlorines(molecule))?;
    d.set_item("num_bromines", chematic_chem::num_bromines(molecule))?;
    d.set_item("num_iodines", chematic_chem::num_iodines(molecule))?;
    Ok(())
}
