//! Descriptor/ADMET/pKa/QSAR bindings.

use crate::{MolHandle, escape_json_string};
use wasm_bindgen::prelude::*;

/// Per-atom EState values as a JSON array of f64.
///
/// Indices match `mol.atoms()` order.  Hydrogen atoms get 0.0.
#[wasm_bindgen]
pub fn estate_indices_json(mol: &MolHandle) -> String {
    let vals = chematic_chem::estate_indices(&mol.inner);
    let parts: Vec<String> = vals.iter().map(|v| format!("{v:.4}")).collect();
    format!("[{}]", parts.join(","))
}

// ---------------------------------------------------------------------------
// Tanimoto with topological path FP (Sprint P)
// ---------------------------------------------------------------------------

/// Detect pharmacophore features for virtual screening and lead optimization.
/// Returns JSON array of features: [{type, atom_idx, neighbor_count}, ...]
#[wasm_bindgen]
pub fn pharmacophore_features_json(mol: &MolHandle) -> String {
    let features = chematic_perception::detect_features(&mol.inner);
    let parts: Vec<String> = features
        .iter()
        .map(|f| {
            let ftype = match f.ftype {
                chematic_perception::FeatureType::Donor => "Donor",
                chematic_perception::FeatureType::Acceptor => "Acceptor",
                chematic_perception::FeatureType::Aromatic => "Aromatic",
                chematic_perception::FeatureType::Hydrophobic => "Hydrophobic",
                chematic_perception::FeatureType::Positive => "Positive",
                chematic_perception::FeatureType::Negative => "Negative",
            };
            format!(
                "{{\"type\":\"{}\",\"atom\":{},\"neighbors\":{}}}",
                ftype,
                f.atom.0,
                f.neighbors.len()
            )
        })
        .collect();
    format!("[{}]", parts.join(","))
}

/// Compute 2D pharmacophore fingerprint (2048 bits) as a JSON feature count summary.
/// Returns simplified JSON with feature type counts: {Donor, Acceptor, Aromatic, Hydrophobic, Positive, Negative}
#[wasm_bindgen]
pub fn pharmacophore_fp_2d_summary(mol: &MolHandle) -> String {
    let counts = chematic_fp::pharmacophore_feature_counts(&mol.inner);
    format!(
        "{{\"Donor\":{},\"Acceptor\":{},\"Aromatic\":{},\"Hydrophobic\":{},\"Positive\":{},\"Negative\":{}}}",
        counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]
    )
}

/// MQN descriptor (42 integer values: Molecular Quantum Numbers).
#[wasm_bindgen]
pub fn mqn_json(mol: &MolHandle) -> String {
    let desc = chematic_chem::mqn(&mol.inner);
    let values: Vec<String> = desc.iter().map(|&v| v.to_string()).collect();
    format!("[{}]", values.join(","))
}

/// AutoCorr2D descriptor (7 values: topological distance lags 1-7).
#[wasm_bindgen]
pub fn autocorr_2d_json(mol: &MolHandle) -> String {
    let ac = chematic_chem::autocorr_2d(&mol.inner);
    format!(
        "[{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}]",
        ac[0], ac[1], ac[2], ac[3], ac[4], ac[5], ac[6]
    )
}

/// XLogP3 partition coefficient (alternative to Crippen LogP).
/// Returns JSON: `{"xlogp3": float}`.
#[wasm_bindgen]
pub fn xlogp3_json(mol: &MolHandle) -> String {
    let v = chematic_chem::xlogp3(&mol.inner);
    format!(r#"{{"xlogp3":{v:.4}}}"#)
}

/// Per-atom XLogP3 contributions.
/// Returns JSON array of floats (one per heavy atom).
#[wasm_bindgen]
pub fn xlogp3_per_atom_json(mol: &MolHandle) -> String {
    let vals = chematic_chem::xlogp3_per_atom(&mol.inner);
    let parts: Vec<String> = vals.iter().map(|v| format!("{v:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// Synthetic Accessibility Score (1 = easy, 10 = hard).
#[wasm_bindgen]
pub fn sa_score(mol: &MolHandle) -> f64 {
    chematic_chem::sa_score(&mol.inner)
}

/// SlogP_VSA descriptors (12 bins) as a JSON array.
#[wasm_bindgen]
pub fn slogp_vsa_json(mol: &MolHandle) -> String {
    let v = chematic_chem::slogp_vsa(&mol.inner);
    let parts: Vec<String> = v.iter().map(|x| format!("{x:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// SMR_VSA descriptors (10 bins) as a JSON array.
#[wasm_bindgen]
pub fn smr_vsa_json(mol: &MolHandle) -> String {
    let v = chematic_chem::smr_vsa(&mol.inner);
    let parts: Vec<String> = v.iter().map(|x| format!("{x:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// PEOE_VSA descriptors (14 bins) as a JSON array.
#[wasm_bindgen]
pub fn peoe_vsa_json(mol: &MolHandle) -> String {
    let v = chematic_chem::peoe_vsa(&mol.inner);
    let parts: Vec<String> = v.iter().map(|x| format!("{x:.4}")).collect();
    format!("[{}]", parts.join(","))
}

// ---------------------------------------------------------------------------
// Sprint T: named functional groups, atom info, per-atom color
// ---------------------------------------------------------------------------

/// All scalar molecular descriptors as a single JSON object.
///
/// Keys use camelCase and match the individual `MolHandle` method names.
/// Drug-likeness rule outcomes are included as boolean fields.
#[wasm_bindgen]
pub fn get_descriptors_json(mol: &MolHandle) -> String {
    let m = &mol.inner;
    format!(
        concat!(
            "{{",
            "\"mw\":{mw:.4},\"exactMass\":{em:.6},\"tpsa\":{tpsa:.4},",
            "\"logP\":{logp:.4},\"molarRefractivity\":{mr:.4},",
            "\"hbd\":{hbd},\"hba\":{hba},\"rotatableBonds\":{rb},",
            "\"heavyAtomCount\":{hac},\"ringCount\":{rc},",
            "\"aromaticRingCount\":{arc},\"numHeteroatoms\":{nh},",
            "\"numStereocenters\":{nsc},\"numSpiroAtoms\":{nsp},",
            "\"numBridgeheadAtoms\":{nbh},\"fsp3\":{fsp3:.4},",
            "\"qed\":{qed:.4},\"saScore\":{sa:.4},",
            "\"formalChargeSum\":{fc},\"labuteASA\":{asa:.4},",
            "\"bertzCT\":{bertz:.4},\"wienerIndex\":{wi:.4},",
            "\"lipinskiPasses\":{lip},\"veberPasses\":{veb},",
            "\"eganPasses\":{egan},\"ghosePasses\":{ghose},",
            "\"reosPasses\":{reos},\"painsPasses\":{pains},",
            "\"bbbScore\":{bbb:.4},\"bbbPasses\":{bbp},",
            "\"caco2\":{caco:.4},\"hergRisk\":{herg:.4},",
            "\"cyp3a4Risk\":{cyp:.4},",
            "\"pkaAcid\":{acid},\"pkaBase\":{base}",
            "}}"
        ),
        mw = chematic_chem::molecular_weight(m),
        em = chematic_chem::exact_mass(m),
        tpsa = chematic_chem::tpsa(m),
        logp = chematic_chem::logp_crippen(m),
        mr = chematic_chem::molar_refractivity(m),
        hbd = chematic_chem::hbd_count(m),
        hba = chematic_chem::hba_count(m),
        rb = chematic_chem::rotatable_bond_count(m),
        hac = chematic_chem::heavy_atom_count(m),
        rc = chematic_chem::ring_count(m),
        arc = chematic_chem::aromatic_ring_count(m),
        nh = chematic_chem::num_heteroatoms(m),
        nsc = chematic_chem::num_stereocenters(m),
        nsp = chematic_chem::num_spiro_atoms(m),
        nbh = chematic_chem::num_bridgehead_atoms(m),
        fsp3 = chematic_chem::fsp3(m),
        qed = chematic_chem::qed(m),
        sa = chematic_chem::sa_score(m),
        fc = chematic_chem::formal_charge_sum(m),
        asa = chematic_chem::labute_asa(m),
        bertz = chematic_chem::bertz_ct(m),
        wi = chematic_chem::wiener_index(m),
        lip = chematic_chem::lipinski_passes(m),
        veb = chematic_chem::veber_passes(m),
        egan = chematic_chem::egan_passes(m),
        ghose = chematic_chem::ghose_passes(m),
        reos = chematic_chem::reos_passes(m),
        pains = chematic_chem::pains_passes(m),
        bbb = chematic_chem::bbb_score(m),
        bbp = chematic_chem::bbb_passes(m),
        caco = chematic_chem::caco2_permeability(m),
        herg = chematic_chem::herg_risk_score(m),
        cyp = chematic_chem::cyp3a4_inhibition_risk(m),
        acid = chematic_chem::pka_acid(m).map_or("null".to_string(), |v| format!("{v:.2}")),
        base = chematic_chem::pka_base(m).map_or("null".to_string(), |v| format!("{v:.2}")),
    )
}

// ---------------------------------------------------------------------------
// pKa + ADMET JSON functions
// ---------------------------------------------------------------------------

/// Predict pKa for all ionizable sites in a molecule.
///
/// Returns a JSON array: `[{"atom_idx":8,"pka":4.0,"type":"acid","group":"carboxylic_acid"},...]`
///
/// Returns `[]` if no ionizable sites are found, or `{"error":"..."}` on parse failure.
#[wasm_bindgen]
pub fn predict_pka_json(smiles: &str) -> String {
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => return format!(r#"{{"error":"{}"}}"#, escape_json_string(&e.to_string())),
    };
    let sites = chematic_chem::predict_pka(&mol);
    if sites.is_empty() {
        return "[]".to_string();
    }
    let items: Vec<String> = sites
        .iter()
        .map(|s| {
            let type_str = match s.site_type {
                chematic_chem::PkaSiteType::Acid => "acid",
                chematic_chem::PkaSiteType::Base => "base",
            };
            format!(
                r#"{{"atom_idx":{},"pka":{:.2},"type":"{}","group":"{}"}}"#,
                s.atom_idx.0, s.pka, type_str, s.group_name
            )
        })
        .collect();
    format!("[{}]", items.join(","))
}

/// Compute a full ADMET property profile for a molecule.
///
/// Returns a JSON object with fields:
/// `bbb_score`, `bbb_passes`, `caco2`, `herg_risk`, `cyp3a4_risk`,
/// `pka_acid` (null if absent), `pka_base` (null if absent),
/// `esol`, `logd74`, `mw`, `logp`, `tpsa`, `hbd`, `hba`, `rotatable_bonds`
///
/// Returns `{"error":"..."}` on parse failure.
#[wasm_bindgen]
pub fn admet_profile_json(smiles: &str) -> String {
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => return format!(r#"{{"error":"{}"}}"#, escape_json_string(&e.to_string())),
    };
    let p = chematic_chem::admet_profile(&mol);
    let acid_str = p.pka_acid.map_or("null".to_string(), |v| format!("{v:.2}"));
    let base_str = p.pka_base.map_or("null".to_string(), |v| format!("{v:.2}"));
    let egg = chematic_chem::boiled_egg(&mol);
    format!(
        concat!(
            "{{",
            "\"bbb_score\":{bbb:.4},\"bbb_passes\":{bbp},",
            "\"caco2\":{caco:.4},\"herg_risk\":{herg:.4},",
            "\"cyp3a4_risk\":{cyp:.4},",
            "\"pka_acid\":{acid},\"pka_base\":{base},",
            "\"esol\":{esol:.4},\"logd74\":{logd:.4},",
            "\"mw\":{mw:.4},\"logp\":{logp:.4},\"tpsa\":{tpsa:.4},",
            "\"hbd\":{hbd},\"hba\":{hba},\"rotatable_bonds\":{rb},",
            "\"gi_absorbed\":{gi},\"bbb_penetrant\":{bbbp}",
            "}}"
        ),
        bbb = p.bbb_score,
        bbp = p.bbb_passes,
        caco = p.caco2,
        herg = p.herg_risk,
        cyp = p.cyp3a4_risk,
        acid = acid_str,
        base = base_str,
        esol = p.esol,
        logd = p.logd74,
        mw = p.mw,
        logp = p.logp,
        tpsa = p.tpsa,
        hbd = p.hbd,
        hba = p.hba,
        rb = p.rotatable_bonds,
        gi = egg.gi_absorbed,
        bbbp = egg.bbb_penetrant,
    )
}

/// Predict GI absorption and BBB penetration using the BOILED-Egg method
/// (Daina & Zoete 2016).
///
/// Returns JSON: `{"gi_absorbed":bool,"bbb_penetrant":bool,"logp":f64,"tpsa":f64}`
#[wasm_bindgen]
pub fn boiled_egg_json(smiles: &str) -> String {
    let mol = match chematic_smiles::parse(smiles) {
        Ok(m) => m,
        Err(e) => return format!(r#"{{"error":"{}"}}"#, escape_json_string(&e.to_string())),
    };
    let e = chematic_chem::boiled_egg(&mol);
    format!(
        r#"{{"gi_absorbed":{},"bbb_penetrant":{},"logp":{:.4},"tpsa":{:.4}}}"#,
        e.gi_absorbed, e.bbb_penetrant, e.logp, e.tpsa
    )
}

// ---------------------------------------------------------------------------
// Sprint X: V3000 read, 3D minimization, SDF properties, highlighted grid
// ---------------------------------------------------------------------------

/// PAINS structural alert names matched by `mol` as a JSON array.
///
/// Returns `[]` when no alerts fire, or e.g. `["ene_six_het_A(483)"]`.
/// Use alongside `pains_passes()` to know *which* alerts triggered.
#[wasm_bindgen]
pub fn pains_matches_json(mol: &MolHandle) -> String {
    let names = chematic_chem::pains_matches(&mol.inner);
    let parts: Vec<String> = names
        .iter()
        .map(|n| format!("\"{}\"", n.replace('"', "\\\"")))
        .collect();
    format!("[{}]", parts.join(","))
}

/// CIP stereo assignments as a JSON array of `{atomIdx, cipCode}` objects.
///
/// `cipCode` is one of `"R"`, `"S"`, `"E"`, or `"Z"`.
/// Returns `[]` for molecules with no specified stereocenters.
#[wasm_bindgen]
pub fn cip_assignments_json(mol: &MolHandle) -> String {
    let assignment = chematic_chem::assign_cip(&mol.inner);
    cip_code_pairs_to_json(&assignment.assignments)
}

/// CIP stereo assignments via the accurate hierarchical-digraph engine, as a JSON
/// array of `{atomIdx, cipCode}` objects -- same shape as [`cip_assignments_json`],
/// but merges the accurate engine's tetrahedral R/S (~99.6% oracle-stable agreement,
/// see `docs/rfcs/cip_accurate_rfc.md`) with legacy's E/Z and allene answers (the accurate
/// engine computes neither). Atoms it can't resolve are omitted here -- see
/// [`cip_unresolved_json`] -- never a silently-guessed label. Returns `"null"` on an
/// internal engine error (budget-independent computations should not normally hit this).
#[wasm_bindgen]
pub fn cip_assignments_accurate_json(mol: &MolHandle) -> String {
    match chematic_chem::assign_cip_with_mode(&mol.inner, chematic_chem::CipMode::Accurate) {
        Ok(result) => cip_code_pairs_to_json(&result.assignments),
        Err(_) => "null".to_string(),
    }
}

/// Atoms the accurate CIP engine could not resolve a tetrahedral R/S for, as a JSON
/// array of `{atomIdx, reason}` objects. `reason` is `"tied"` (a genuine CIP-rule tie,
/// not a missing rule) or `"budgetExceeded"`. Always `[]` for the legacy engine (see
/// [`cip_assignments_json`]) -- it never reports "I don't know". Returns `"null"` on
/// an internal engine error.
#[wasm_bindgen]
pub fn cip_unresolved_json(mol: &MolHandle) -> String {
    match chematic_chem::assign_cip_with_mode(&mol.inner, chematic_chem::CipMode::Accurate) {
        Ok(result) => {
            let parts: Vec<String> = result
                .unresolved
                .iter()
                .map(|(idx, reason)| {
                    let reason_str = match reason {
                        chematic_chem::CipUnresolvedReason::Tied => "tied",
                        chematic_chem::CipUnresolvedReason::BudgetExceeded => "budgetExceeded",
                    };
                    format!("{{\"atomIdx\":{},\"reason\":\"{}\"}}", idx.0, reason_str)
                })
                .collect();
            format!("[{}]", parts.join(","))
        }
        Err(_) => "null".to_string(),
    }
}

fn cip_code_pairs_to_json(pairs: &[(chematic_core::AtomIdx, chematic_core::CipCode)]) -> String {
    let parts: Vec<String> = pairs
        .iter()
        .map(|(idx, code)| {
            let code_str = match code {
                chematic_core::CipCode::R => "R",
                chematic_core::CipCode::S => "S",
                chematic_core::CipCode::E => "E",
                chematic_core::CipCode::Z => "Z",
                chematic_core::CipCode::LowerR => "r",
                chematic_core::CipCode::LowerS => "s",
            };
            format!("{{\"atomIdx\":{},\"cipCode\":\"{}\"}}", idx.0, code_str)
        })
        .collect();
    format!("[{}]", parts.join(","))
}

/// Parse a ChemDraw XML (CDXML) string into a `MolHandle`.
///
/// Compute an HDF fingerprint and return it as a JSON array of float32 values.
///
/// Returns a unit-norm vector of length `dim` as a JSON number array.
/// Use cosine dot product for similarity: `a · b = sum(a[i]*b[i])`.
///
/// ```js
/// const fp = JSON.parse(hdf_json(mol));        // float[] of length 1024
/// const sim = fp.reduce((s, v, i) => s + v * fp2[i], 0);  // cosine similarity
/// ```
#[wasm_bindgen]
pub fn hdf_json(mol: &MolHandle, dim: usize, radius: usize, seed: u64) -> String {
    let config = chematic_fp::HdfConfig { dim, radius, seed };
    let fp = chematic_fp::hdf(&mol.inner, &config);
    serde_json::to_string(&fp.0).unwrap_or_default()
}

/// Per-atom Crippen LogP contributions as a JSON array of f64.
///
/// Index `i` corresponds to atom `i` in `mol.atoms()` order.
#[wasm_bindgen]
pub fn logp_per_atom_json(mol: &MolHandle) -> String {
    let vals = chematic_chem::logp_crippen_per_atom(&mol.inner);
    let parts: Vec<String> = vals.iter().map(|v| format!("{v:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// Per-atom molar refractivity contributions as a JSON array of f64.
#[wasm_bindgen]
pub fn mr_per_atom_json(mol: &MolHandle) -> String {
    let vals = chematic_chem::mr_per_atom(&mol.inner);
    let parts: Vec<String> = vals.iter().map(|v| format!("{v:.4}")).collect();
    format!("[{}]", parts.join(","))
}

/// Per-atom Labute approximate surface area contributions as a JSON array of f64.
///
/// Non-finite values (single-atom molecules etc.) are emitted as JSON `null`.
#[wasm_bindgen]
pub fn labute_asa_per_atom_json(mol: &MolHandle) -> String {
    let vals = chematic_chem::labute_asa_per_atom(&mol.inner);
    let parts: Vec<String> = vals
        .iter()
        .map(|v| {
            if v.is_finite() {
                format!("{v:.4}")
            } else {
                "null".to_string()
            }
        })
        .collect();
    format!("[{}]", parts.join(","))
}
