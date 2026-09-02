//! 3D geometry, conformer, and force-field (MMFF94/UFF) bindings.

use crate::mol_io::coords_all_finite;
use crate::{
    MolHandle, WASM_MAX_ATOMS, WASM_MAX_INPUT_BYTES, WASM_MAX_JSON_STRING_BYTES, escape_json_string,
};
use wasm_bindgen::prelude::*;

/// Optimize molecular geometry using DREIDING force field.
///
/// Performs geometry minimization with DREIDING force field parameters.
/// Returns minimized coordinate PDB.
///
/// # Arguments
/// * `mol` - Molecule to optimize
///
/// # Returns
/// PDB format string with optimized coordinates
#[wasm_bindgen]
pub fn minimize_dreiding_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let min_coords = chematic_3d::minimize_dreiding(&mol.inner, coords);
    chematic_3d::write_pdb(&mol.inner, &min_coords)
}

/// Compute WHIM descriptors (Weighted Holistic Invariant Molecular) from 3D coordinates.
/// Returns JSON array of 22 values: 11 unit-weight descriptors followed by 11 mass-weight
/// descriptors. Each 11-element block is [λ₁, λ₂, λ₃, ν₁, ν₂, ν₃, T, A, V, K, D].
#[wasm_bindgen]
pub fn whim_descriptors_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let desc = chematic_3d::whim_descriptors(&mol.inner, &coords);
    let parts: Vec<String> = desc.iter().map(|&v| v.to_string()).collect();
    format!("[{}]", parts.join(","))
}

/// Compute GETAWAY descriptors (GEometry, Topology and Atom-Weights AssemblY) from 3D coords.
///
/// Returns a JSON array of **19** values:
/// - `[0..7]`  H[1..8]  — leverage autocorrelation at topological lags 1–8
/// - `[8..15]` R[1..8]  — H[k] normalised by pair count W_k
/// - `[16]` Hmax, `[17]` Hmean, `[18]` Htot — per-atom leverage statistics
///
/// Note: requires 3D coordinates (non-planar); for flat/2D structures the hat matrix
/// is degenerate and descriptors reflect squared centroid distances, not true leverage.
#[wasm_bindgen]
pub fn getaway_descriptors_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let desc = chematic_3d::getaway_descriptors(&mol.inner, &coords);
    let parts: Vec<String> = desc.iter().map(|&v| v.to_string()).collect();
    format!("[{}]", parts.join(","))
}

/// Compute combined WHIM + GETAWAY descriptors (**41** values total) as JSON array.
///
/// Returns WHIM[0..21] (22 values) followed by GETAWAY[0..18] (19 values) = 41 total.
/// Useful for ML pipelines requiring both shape and topologic features.
#[wasm_bindgen]
pub fn whim_getaway_combined_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let desc = chematic_3d::whim_getaway_combined(&mol.inner, &coords);
    let parts: Vec<String> = desc.iter().map(|&v| v.to_string()).collect();
    format!("[{}]", parts.join(","))
}

/// Gasteiger-Marsili PEOE partial charges as a JSON array of f64.
#[wasm_bindgen]
pub fn gasteiger_charges_json(mol: &MolHandle) -> String {
    let q = chematic_chem::gasteiger_charges(&mol.inner);
    let parts: Vec<String> = q.iter().map(|v| format!("{v:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// MMFF94 partial charges (BCI table, ±0.1e accuracy) as a JSON array of f64.
///
/// Uses Bond Charge Increment (BCI) model (Halgren 1996) for 25 common bond types.
/// Returns `[q0, q1, ..., qN]` — one value per heavy atom.
/// Total charge equals the sum of formal charges (charge conserved).
#[wasm_bindgen]
pub fn mmff94_charges_json(mol: &MolHandle) -> String {
    let q = chematic_chem::mmff94_charges(&mol.inner);
    let parts: Vec<String> = q.iter().map(|v| format!("{v:.6}")).collect();
    format!("[{}]", parts.join(","))
}

/// Compute MMFF94-style atom-typed partial charges (improved over element-pair BCI).
/// Returns JSON: {"charges":[f64,...]} or {"error":"..."}.
/// Uses atom-type classification (Csp3/Ccarbonyl/Ohydroxyl/Oester/Nar/NarH etc.)
/// for better accuracy (~±0.02e) vs element-pair BCI (~±0.05e).
#[wasm_bindgen]
pub fn mmff94_charges_typed_json(mol: &MolHandle) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    let q = chematic_chem::mmff94_charges_typed(&mol.inner);
    let vals: Vec<String> = q.iter().map(|x| format!("{x:.6}")).collect();
    format!(r#"{{"charges":[{}]}}"#, vals.join(","))
}

// ── MMFF94 Full Force Field (v0.2.10) ───────────────────────────────────────

/// Minimize geometry using MMFF94 steepest descent (Halgren 1996 full parameters).
/// Generates 3D coords internally if needed.
/// Returns JSON: {"energy":E,"rmsd":R,"converged":true,"iterations":N} or {"error":"..."}.
#[wasm_bindgen]
pub fn minimize_mmff94_json(mol: &MolHandle, max_iter: u32) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    let conf = chematic_3d::generate_coords(&mol.inner);
    let n = mol.inner.atom_count();
    let mut coords: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let p = conf.get(chematic_core::AtomIdx(i as u32));
            [p.x, p.y, p.z]
        })
        .collect();
    match chematic_ff::minimize_mmff94_full(&mol.inner, &mut coords, max_iter as usize) {
        Ok(r) => format!(
            r#"{{"energy":{:.4},"rmsd":{:.4},"converged":{},"iterations":{}}}"#,
            r.energy, r.rmsd, r.converged, r.iterations
        ),
        Err(e) => format!(r#"{{"error":"{}"}}"#, e),
    }
}

/// Minimize geometry using MMFF94 L-BFGS (faster convergence than steepest descent).
/// Returns JSON: {"energy":E,"rmsd":R,"converged":true,"iterations":N} or {"error":"..."}.
#[wasm_bindgen]
pub fn minimize_mmff94_lbfgs_json(mol: &MolHandle, max_iter: u32) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    let conf = chematic_3d::generate_coords(&mol.inner);
    let n = mol.inner.atom_count();
    let mut coords: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let p = conf.get(chematic_core::AtomIdx(i as u32));
            [p.x, p.y, p.z]
        })
        .collect();
    match chematic_ff::minimize_mmff94_lbfgs(&mol.inner, &mut coords, max_iter as usize) {
        Ok(r) => format!(
            r#"{{"energy":{:.4},"rmsd":{:.4},"converged":{},"iterations":{}}}"#,
            r.energy, r.rmsd, r.converged, r.iterations
        ),
        Err(e) => format!(r#"{{"error":"{}"}}"#, e),
    }
}

/// Computes energy on an internally generated conformer.
/// `MolHandle` stores topology only; coordinates previously read from PDB/XYZ
/// are not used by this function.
/// Use [`mmff94_energy_breakdown_from_coords_json`] for explicit coordinates.
///
/// Returns JSON: {"bond":B,"angle":A,"torsion":T,"vdw":V,"elec":E,"total":X} or {"error":"..."}.
#[wasm_bindgen]
pub fn mmff94_energy_breakdown_json(mol: &MolHandle) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    let conf = chematic_3d::generate_coords(&mol.inner);
    let n = mol.inner.atom_count();
    let coords: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let p = conf.get(chematic_core::AtomIdx(i as u32));
            [p.x, p.y, p.z]
        })
        .collect();
    match chematic_ff::mmff94_energy_breakdown(&mol.inner, &coords) {
        Ok(bd) => format!(
            r#"{{"bond":{:.4},"angle":{:.4},"torsion":{:.4},"vdw":{:.4},"elec":{:.4},"total":{:.4}}}"#,
            bd.bond, bd.angle, bd.torsion, bd.vdw, bd.electrostatic, bd.total
        ),
        Err(e) => format!(r#"{{"error":"{}"}}"#, e),
    }
}

/// Compute MMFF94 energy breakdown for EXPLICIT, caller-supplied 3D
/// coordinates -- unlike `mmff94_energy_breakdown_json`, this reads the
/// geometry the caller actually has (e.g. from `pdb_coords_json`,
/// `generate_3d_coords_json`, `generate_3d_etkdg_coords_json`, or an
/// externally computed conformer) instead of silently generating a fresh
/// rule-based one, matching the Python binding's
/// `mol.mmff94_energy_breakdown(coords)` contract (issue #90).
///
/// `coords_json` -- JSON array of `[x,y,z]` arrays (Å), one per heavy atom,
/// in the same atom order as `mol`.
///
/// Returns JSON `{"bond":B,"angle":A,"stretch_bend":S,"torsion":T,"oop":O,
/// "vdw":V,"electrostatic":E,"total":X}` at full `f64` round-trip precision
/// (not rounded -- this API exists specifically for oracle comparison
/// against the Python binding, where 4-decimal rounding would mask
/// sub-1e-4 discrepancies), or `{"error":"<msg>"}` on malformed JSON, a
/// non-finite coordinate, or a coordinate-count/atom-count mismatch. Never
/// falls back to a generated conformer.
#[wasm_bindgen]
pub fn mmff94_energy_breakdown_from_coords_json(mol: &MolHandle, coords_json: &str) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    if coords_json.len() > WASM_MAX_JSON_STRING_BYTES {
        return format!(
            r#"{{"error":"coords_json too large (max {} bytes)"}}"#,
            WASM_MAX_JSON_STRING_BYTES
        );
    }
    let coords: Vec<[f64; 3]> = match serde_json::from_str(coords_json) {
        Ok(v) => v,
        Err(e) => return format!(r#"{{"error":"coords parse error: {e}"}}"#),
    };
    if coords.len() != mol.inner.atom_count() {
        return format!(
            r#"{{"error":"coords length {} does not match atom count {}"}}"#,
            coords.len(),
            mol.inner.atom_count()
        );
    }
    if !coords_all_finite(&coords) {
        return r#"{"error":"coords contain a non-finite value (NaN or Infinity)"}"#.to_string();
    }
    match chematic_ff::mmff94_energy_breakdown(&mol.inner, &coords) {
        Ok(bd) => {
            let vals = [
                bd.bond,
                bd.angle,
                bd.stretch_bend,
                bd.torsion,
                bd.oop,
                bd.vdw,
                bd.electrostatic,
                bd.total,
            ];
            if vals.iter().any(|v| !v.is_finite()) {
                return r#"{"error":"computed energy is not finite (degenerate geometry, e.g. coincident atoms)"}"#
                    .to_string();
            }
            format!(
                r#"{{"bond":{},"angle":{},"stretch_bend":{},"torsion":{},"oop":{},"vdw":{},"electrostatic":{},"total":{}}}"#,
                bd.bond,
                bd.angle,
                bd.stretch_bend,
                bd.torsion,
                bd.oop,
                bd.vdw,
                bd.electrostatic,
                bd.total
            )
        }
        Err(e) => format!(r#"{{"error":"{}"}}"#, e),
    }
}

// torsion_scan_json removed to reduce WASM bundle size (force-field loop; rarely used in browser).

/// Compute MMFF94 partial charges using numeric atom types (Halgren 1996 eq. 15).
/// Returns JSON: {"charges":[-0.28,0.15,...]} or {"error":"..."}.
#[wasm_bindgen]
pub fn mmff94_partial_charges_json(mol: &MolHandle) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    match chematic_ff::mmff94_charges_numeric(&mol.inner) {
        Ok(q) => {
            let vals: Vec<String> = q.iter().map(|x| format!("{x:.4}")).collect();
            format!(r#"{{"charges":[{}]}}"#, vals.join(","))
        }
        Err(e) => format!(r#"{{"error":"{}"}}"#, e),
    }
}

/// Compute 3D pharmacophore fingerprint from generated 3D coordinates.
/// Returns simplified JSON with feature type counts (3D-aware version).
#[wasm_bindgen]
pub fn pharmacophore_fp_3d_summary(mol: &MolHandle) -> String {
    let _coords = chematic_3d::generate_coords(&mol.inner);
    let features = chematic_perception::detect_features(&mol.inner);
    let mut counts = [0; 6];
    for feat in features {
        let idx = match feat.ftype {
            chematic_perception::FeatureType::Donor => 0,
            chematic_perception::FeatureType::Acceptor => 1,
            chematic_perception::FeatureType::Aromatic => 2,
            chematic_perception::FeatureType::Hydrophobic => 3,
            chematic_perception::FeatureType::Positive => 4,
            chematic_perception::FeatureType::Negative => 5,
        };
        counts[idx] += 1;
    }
    format!(
        "{{\"Donor\":{},\"Acceptor\":{},\"Aromatic\":{},\"Hydrophobic\":{},\"Positive\":{},\"Negative\":{}}}",
        counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]
    )
}

/// AutoCorr3D descriptor (8 values: Euclidean distance bins 1-8 Å).
/// Requires 3D coordinates (generated automatically).
#[wasm_bindgen]
pub fn autocorr_3d_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let ac = chematic_3d::autocorr_3d(&mol.inner, &coords);
    format!(
        "[{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}]",
        ac[0], ac[1], ac[2], ac[3], ac[4], ac[5], ac[6], ac[7]
    )
}

/// Generate 3D coordinates as raw JSON array [[x,y,z], ...].
///
/// Unlike `generate_3d_pdb`, this returns coordinates that can be passed
/// to descriptor functions like `whim_descriptors_json` or `shape_descriptors_json`.
#[wasm_bindgen]
pub fn generate_3d_coords_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_and_minimize_dreiding(&mol.inner);
    let parts: Vec<String> = coords
        .points
        .iter()
        .map(|p| format!("[{:.4},{:.4},{:.4}]", p.x, p.y, p.z))
        .collect();
    format!("[{}]", parts.join(","))
}

/// Generate 3D coordinates using ETKDG as raw JSON array [[x,y,z], ...].
#[wasm_bindgen]
pub fn generate_3d_etkdg_coords_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords_etkdg(&mol.inner);
    let parts: Vec<String> = coords
        .points
        .iter()
        .map(|p| format!("[{:.4},{:.4},{:.4}]", p.x, p.y, p.z))
        .collect();
    format!("[{}]", parts.join(","))
}

/// Generate multiple conformers with RMSD-based pruning.
/// Returns JSON: `{"conformers": [[[x,y,z],...], ...], "count": int}`.
#[wasm_bindgen]
pub fn conformer_ensemble_json(mol: &MolHandle, n: u32, rmsd_threshold: f64) -> String {
    let smiles = chematic_smiles::canonical_smiles(&mol.inner);
    let fresh = match chematic_smiles::parse(&smiles) {
        Ok(m) => m,
        Err(_) => return r#"{"conformers":[],"count":0}"#.to_string(),
    };
    let config = chematic_3d::ConformerConfig {
        count: n as usize,
        rmsd_threshold,
        ..chematic_3d::ConformerConfig::default()
    };
    match chematic_3d::generate_conformer_ensemble_with_config(fresh, &config) {
        Ok(ensemble) => {
            let conf_jsons: Vec<String> = (0..ensemble.conformer_count())
                .filter_map(|i| ensemble.get_conformer(i))
                .map(|c| {
                    let pts: Vec<String> = c
                        .points
                        .iter()
                        .map(|p| format!("[{:.4},{:.4},{:.4}]", p.x, p.y, p.z))
                        .collect();
                    format!("[{}]", pts.join(","))
                })
                .collect();
            let count = conf_jsons.len();
            format!(
                r#"{{"conformers":[{}],"count":{count}}}"#,
                conf_jsons.join(",")
            )
        }
        Err(_) => r#"{"conformers":[],"count":0}"#.to_string(),
    }
}

/// 3D shape descriptors as a JSON object.
///
/// Keys: `pmi1`, `pmi2`, `pmi3`, `npr1`, `npr2`, `asphericity`, `eccentricity`,
/// `radiusOfGyration`, `planeOfBestFit`.  Non-finite values (e.g. single-atom
/// molecules where pmi3 = 0) are serialised as JSON `null`.
#[wasm_bindgen]
pub fn shape_descriptors_json(mol: &MolHandle) -> String {
    let coords = chematic_3d::generate_coords(&mol.inner);
    let m = &*mol.inner;

    fn fmt(v: f64) -> String {
        if v.is_finite() {
            format!("{v:.4}")
        } else {
            "null".to_string()
        }
    }

    format!(
        "{{\"pmi1\":{},\"pmi2\":{},\"pmi3\":{},\"npr1\":{},\"npr2\":{},\
\"asphericity\":{},\"eccentricity\":{},\"radiusOfGyration\":{},\"planeOfBestFit\":{}}}",
        fmt(chematic_3d::pmi1(m, &coords)),
        fmt(chematic_3d::pmi2(m, &coords)),
        fmt(chematic_3d::pmi3(m, &coords)),
        fmt(chematic_3d::npr1(m, &coords)),
        fmt(chematic_3d::npr2(m, &coords)),
        fmt(chematic_3d::asphericity(m, &coords)),
        fmt(chematic_3d::eccentricity(m, &coords)),
        fmt(chematic_3d::radius_of_gyration(m, &coords)),
        fmt(chematic_3d::plane_of_best_fit(m, &coords)),
    )
}

/// A conformer ensemble: one molecule geometry with multiple 3D coordinate sets.
///
/// Create with `new(smiles)`, then add conformers with `add_generated_conformer`
/// or `add_minimized_conformer`.  Retrieve coordinates as PDB strings via
/// `get_conformer_pdb(idx)`.  Compare conformers with `conformer_rmsd`.
#[wasm_bindgen]
pub struct ConformerHandle {
    inner: chematic_3d::ConformerEnsemble,
}

#[wasm_bindgen]
impl ConformerHandle {
    /// Create a new empty ensemble for the molecule given by `smiles`.
    ///
    /// Returns a JS error on SMILES parse failure.
    #[wasm_bindgen(constructor)]
    pub fn new(smiles: &str) -> Result<ConformerHandle, JsValue> {
        let mol = chematic_smiles::parse(smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
        Ok(ConformerHandle {
            inner: chematic_3d::ConformerEnsemble::new(mol),
        })
    }

    /// Number of conformers currently stored.
    pub fn conformer_count(&self) -> usize {
        self.inner.conformer_count()
    }

    /// The ensemble's molecule as a `MolHandle`.
    pub fn mol(&self) -> MolHandle {
        let smi = chematic_smiles::canonical_smiles(self.inner.mol());
        let mol = chematic_smiles::parse(&smi)
            .unwrap_or_else(|_| chematic_core::MoleculeBuilder::new().build());
        MolHandle {
            inner: std::rc::Rc::new(mol),
        }
    }

    /// Generate a new 3D conformer using distance-geometry and add it to the ensemble.
    ///
    /// Returns the index of the newly added conformer.
    pub fn add_generated_conformer(&mut self) -> usize {
        let coords = chematic_3d::generate_coords(self.inner.mol());
        self.inner.add_conformer(coords).unwrap_or(usize::MAX)
    }

    /// Generate a new 3D conformer, run force-field minimization, and add it.
    ///
    /// Returns the index of the newly added conformer.
    pub fn add_minimized_conformer(&mut self) -> usize {
        let coords = chematic_3d::generate_coords(self.inner.mol());
        let minimized = chematic_3d::minimize(self.inner.mol(), coords);
        self.inner.add_conformer(minimized).unwrap_or(usize::MAX)
    }

    /// Return conformer `idx` as a PDB string, or `null` if `idx` is out of range.
    pub fn get_conformer_pdb(&self, idx: usize) -> Option<String> {
        self.inner
            .get_conformer(idx)
            .map(|coords| chematic_3d::write_pdb(self.inner.mol(), coords))
    }

    /// Kabsch-aligned RMSD (Å) between conformers `a` and `b`.
    ///
    /// Returns `NaN` if either index is out of range.
    pub fn conformer_rmsd(&self, a: usize, b: usize) -> f64 {
        self.inner.conformer_rmsd(a, b).unwrap_or(f64::NAN)
    }

    /// Un-aligned (translation + rotation NOT removed) RMSD (Å) between conformers `a` and `b`.
    ///
    /// Returns `NaN` if either index is out of range.
    pub fn conformer_rmsd_no_align(&self, a: usize, b: usize) -> f64 {
        self.inner.conformer_rmsd_no_align(a, b).unwrap_or(f64::NAN)
    }

    /// Remove conformer `idx` and return `true`, or `false` if `idx` is out of range.
    pub fn remove_conformer(&mut self, idx: usize) -> bool {
        self.inner.remove_conformer(idx).is_some()
    }

    /// Cluster conformers by Kabsch-aligned RMSD and return a JSON object
    /// describing which conformers to keep.
    ///
    /// Uses greedy leader-linkage: conformers are visited in index order; each
    /// is compared against existing cluster representatives. If the RMSD to any
    /// representative is < `rms_threshold`, the conformer is discarded; otherwise
    /// it starts a new cluster and is kept.
    ///
    /// Returns `{"kept_indices":[0,3,7,...],"removed_count":5}` on success.
    pub fn cluster_conformers_json(&self, rms_threshold: f64) -> String {
        let kept = self.inner.cluster_conformers_by_rms(rms_threshold);
        let total = self.inner.conformer_count();
        let removed = total.saturating_sub(kept.len());
        let indices_str: Vec<String> = kept.iter().map(|i| i.to_string()).collect();
        format!(
            r#"{{"kept_indices":[{}],"removed_count":{}}}"#,
            indices_str.join(","),
            removed
        )
    }
}

// ---------------------------------------------------------------------------
// Mutable Molecule API
// ---------------------------------------------------------------------------

/// Compute structured depiction data using caller-supplied 2D coordinates.
///
/// `coords_json` — JSON array of `[x, y]` pairs, one per atom in order.
///
/// Returns the same JSON format as `depict_data_json`.
#[wasm_bindgen]
pub fn depict_data_with_coords_json(mol: &MolHandle, coords_json: &str) -> String {
    if coords_json.len() > WASM_MAX_JSON_STRING_BYTES {
        return "{\"error\":\"coords_json too large\"}".to_string();
    }
    if mol.atom_count() > WASM_MAX_ATOMS {
        return format!("{{\"error\":\"molecule too large (max {WASM_MAX_ATOMS} atoms)\"}}");
    }
    // Parse [[x,y],...] — simple format, no full JSON parser needed.
    let coords: Vec<(f64, f64)> = {
        let mut result = Vec::new();
        let inner = coords_json.trim().trim_matches(|c| c == '[' || c == ']');
        // Split on "]," pattern to get individual [x,y] pairs.
        for pair in inner.split("],[") {
            let clean = pair.trim().trim_matches(|c| c == '[' || c == ']');
            let nums: Vec<f64> = clean
                .split(',')
                .filter_map(|s| s.trim().parse::<f64>().ok())
                .collect();
            if nums.len() >= 2 {
                result.push((nums[0], nums[1]));
            }
        }
        result
    };

    let data = chematic_depict::depict_data_with_coords(&mol.inner, &coords);

    let atoms: Vec<String> = data
        .atoms
        .iter()
        .map(|a| {
            let label = match &a.label {
                Some(s) => format!("\"{}\"", escape_json_string(s)),
                None => "null".to_string(),
            };
            let charge_suffix = if a.charge != 0 {
                format!(",\"charge\":{}", a.charge)
            } else {
                String::new()
            };
            let label_prefix = format!("{label},");
            format!(
                r#"{{"idx":{},"element":"{}","x":{:.4},"y":{:.4},"label":{}"color":"{}"{}}}"#,
                a.idx.0,
                a.element.symbol(),
                a.pos.x,
                a.pos.y,
                label_prefix,
                escape_json_string(&a.color),
                charge_suffix,
            )
        })
        .collect();

    let bond_kind = |k: &chematic_depict::DepictBondKind| match k {
        chematic_depict::DepictBondKind::Single => "Single",
        chematic_depict::DepictBondKind::Double => "Double",
        chematic_depict::DepictBondKind::Triple => "Triple",
        chematic_depict::DepictBondKind::Aromatic => "Aromatic",
        chematic_depict::DepictBondKind::Up => "Up",
        chematic_depict::DepictBondKind::Down => "Down",
    };

    let bonds: Vec<String> = data
        .bonds
        .iter()
        .map(|b| {
            format!(
                r#"{{"idx":{},"atom1":{},"atom2":{},"kind":"{}"}}"#,
                b.idx.0,
                b.atom1.0,
                b.atom2.0,
                bond_kind(&b.kind)
            )
        })
        .collect();

    format!(
        r#"{{"atoms":[{}],"bonds":[{}]}}"#,
        atoms.join(","),
        bonds.join(",")
    )
}

// ---------------------------------------------------------------------------
// Sprint CC: MMP analysis (Matched Molecular Pairs)
// ---------------------------------------------------------------------------

/// Minimise a molecule's geometry using the Universal Force Field (UFF).
///
/// `coords_json` — JSON array of `[x,y,z]` arrays (Å), one per atom.
/// `max_iter` — maximum iterations (0 = default 500).
///
/// Returns JSON: `{"coords":[[x,y,z],...], "energy":float, "iterations":int, "converged":bool, "sound":bool}`
/// or `{"error":"<msg>"}` on failure. `sound` is all-finite coordinates and
/// no bond stretched past a sane covalent-bond length — independent of
/// `converged`, since steepest descent often reports `converged:false` on
/// geometries that are perfectly fine but simply haven't hit the tight
/// RMS-gradient threshold yet. Check `sound`, not just `converged`, before
/// trusting a result.
#[wasm_bindgen]
pub fn minimize_uff_json(smiles: &str, coords_json: &str, max_iter: u32) -> String {
    if smiles.len() > WASM_MAX_INPUT_BYTES || coords_json.len() > WASM_MAX_JSON_STRING_BYTES {
        return "{\"error\":\"input exceeds WASM size limits\"}".to_string();
    }
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => return format!("{{\"error\":\"{e}\"}}"),
    };
    if mol.atom_count() > WASM_MAX_ATOMS {
        return format!("{{\"error\":\"molecule too large (max {WASM_MAX_ATOMS} atoms)\"}}");
    }
    let coords_arr: Vec<[f64; 3]> = match serde_json::from_str(coords_json) {
        Ok(v) => v,
        Err(e) => return format!("{{\"error\":\"coords parse error: {e}\"}}"),
    };
    let types = chematic_ff::assign_uff_types(&mol);
    let iters = if max_iter == 0 {
        500
    } else {
        max_iter as usize
    };
    let result = chematic_ff::minimize_uff(&mol, &types, coords_arr, iters);
    let coords_out: Vec<[f64; 3]> = result.coords;
    let coords_str = coords_out
        .iter()
        .map(|c| format!("[{:.4},{:.4},{:.4}]", c[0], c[1], c[2]))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"coords\":[{coords_str}],\"energy\":{:.4},\"iterations\":{},\"converged\":{},\"sound\":{}}}",
        result.energy, result.iterations, result.converged, result.sound
    )
}

// ---------------------------------------------------------------------------
// InChI I/O
// ---------------------------------------------------------------------------

/// Get bond length in Ångströms between two atoms from a SMILES string.
/// Returns -1.0 if parsing fails or atom indices are out of range.
///
/// # Arguments
/// - `smiles`: SMILES string
/// - `a`: first atom index
/// - `b`: second atom index
///
/// # Example
/// ```javascript
/// const len = get_bond_length_json("CC", 0, 1);  // C-C single bond ≈ 1.54 Å
/// ```
#[wasm_bindgen]
pub fn get_bond_length_json(smiles: &str, a: u32, b: u32) -> f64 {
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return -1.0;
    }
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(_) => return -1.0,
    };
    let coords = chematic_3d::generate_coords(&mol);
    let a_idx = chematic_core::AtomIdx(a);
    let b_idx = chematic_core::AtomIdx(b);
    if a >= mol.atom_count() as u32 || b >= mol.atom_count() as u32 {
        return -1.0;
    }
    chematic_3d::get_bond_length(&coords, a_idx, b_idx)
}

/// Get dihedral angle A—B—C—D in degrees from a SMILES string.
/// Returns null (JSON null) if any atom index is out of range or atoms are collinear.
///
/// # Arguments
/// - `smiles`: SMILES string
/// - `a`, `b`, `c`, `d`: atom indices
///
/// # Example
/// ```javascript
/// const dihedral = get_dihedral_json("CCCC", 0, 1, 2, 3);  // A-B-C-D
/// ```
#[wasm_bindgen]
pub fn get_dihedral_json(smiles: &str, a: u32, b: u32, c: u32, d: u32) -> JsValue {
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return JsValue::NULL;
    }
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(_) => return JsValue::NULL,
    };
    if a >= mol.atom_count() as u32
        || b >= mol.atom_count() as u32
        || c >= mol.atom_count() as u32
        || d >= mol.atom_count() as u32
    {
        return JsValue::NULL;
    }
    let coords = chematic_3d::generate_coords(&mol);
    let a_idx = chematic_core::AtomIdx(a);
    let b_idx = chematic_core::AtomIdx(b);
    let c_idx = chematic_core::AtomIdx(c);
    let d_idx = chematic_core::AtomIdx(d);
    match chematic_3d::get_dihedral_deg(&coords, a_idx, b_idx, c_idx, d_idx) {
        Some(angle) => JsValue::from_f64(angle),
        None => JsValue::NULL,
    }
}

/// Set dihedral angle A—B—C—D and return PDB block with modified coordinates.
/// Rotates the D-side subtree around the B—C bond.
/// Returns a JS error if parsing fails or atom indices are invalid.
///
/// # Arguments
/// - `smiles`: SMILES string
/// - `a`, `b`, `c`, `d`: atom indices
/// - `angle_deg`: target dihedral angle in degrees
///
/// # Example
/// ```javascript
/// const pdbBlock = set_dihedral_json("CCCC", 0, 1, 2, 3, 120.0);
/// ```
#[wasm_bindgen]
pub fn set_dihedral_json(
    smiles: &str,
    a: u32,
    b: u32,
    c: u32,
    d: u32,
    angle_deg: f64,
) -> Result<String, String> {
    if smiles.len() > WASM_MAX_INPUT_BYTES {
        return Err("Input SMILES too long".to_string());
    }
    let mol = chematic_smiles::parse(smiles).map_err(|e| format!("Parse error: {}", e))?;
    if a >= mol.atom_count() as u32
        || b >= mol.atom_count() as u32
        || c >= mol.atom_count() as u32
        || d >= mol.atom_count() as u32
    {
        return Err("Atom index out of range".to_string());
    }
    let coords = chematic_3d::generate_coords(&mol);
    let a_idx = chematic_core::AtomIdx(a);
    let b_idx = chematic_core::AtomIdx(b);
    let c_idx = chematic_core::AtomIdx(c);
    let d_idx = chematic_core::AtomIdx(d);
    let angle_rad = angle_deg.to_radians();
    let new_coords =
        chematic_3d::set_dihedral(&coords, &mol, a_idx, b_idx, c_idx, d_idx, angle_rad);
    // Return PDB block with modified coordinates
    Ok(chematic_3d::write_pdb(&mol, &new_coords))
}

// ---------------------------------------------------------------------------
// random_smiles (SMILES augmentation)
// ---------------------------------------------------------------------------
