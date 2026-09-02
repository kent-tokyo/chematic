//! Fingerprint and similarity-search bindings: ECFP/FCFP/MACCS/MHFP, Tanimoto/Dice, nearest-neighbor search.

use crate::{
    MolHandle, WASM_MAX_ATOMS, enforce_wasm_input_len, enforce_wasm_molecule_size,
    parse_smiles_json_array,
};
use wasm_bindgen::prelude::*;

/// Tanimoto similarity between two molecules using ECFP4 fingerprints.
#[wasm_bindgen]
pub fn tanimoto_ecfp4(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::tanimoto_ecfp4(&a.inner, &b.inner)
}

/// Tanimoto similarity between two molecules using FCFP4 fingerprints (pharmacophore-based).
#[wasm_bindgen]
pub fn tanimoto_fcfp4(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::tanimoto_fcfp4(&a.inner, &b.inner)
}

/// Tanimoto similarity between two molecules using AtomPair fingerprints.
#[wasm_bindgen]
pub fn tanimoto_atom_pair(a: &MolHandle, b: &MolHandle) -> f64 {
    let fa = chematic_fp::atom_pair_fp(&a.inner);
    let fb = chematic_fp::atom_pair_fp(&b.inner);
    fa.tanimoto(&fb)
}

/// Tanimoto similarity between two molecules using Topological Torsion fingerprints.
#[wasm_bindgen]
pub fn tanimoto_torsion(a: &MolHandle, b: &MolHandle) -> f64 {
    let fa = chematic_fp::torsion_fp(&a.inner);
    let fb = chematic_fp::torsion_fp(&b.inner);
    fa.tanimoto(&fb)
}

/// Compute the ECFP4 fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
#[wasm_bindgen]
pub fn ecfp4_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::ecfp4(&mol.inner);
    // BitVec2048 is 2048 bits; extract them byte-by-byte via the public `get` method.
    (0..256usize)
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

/// Like `ecfp4_bitvec` but with explicit chirality control.
///
/// When `use_chirality=true`, tetrahedral stereochemistry is included in the
/// initial atom hash, making enantiomers have different fingerprints.
/// When `false` (default), chirality is ignored.
#[wasm_bindgen]
pub fn ecfp4_bitvec_with_chirality(mol: &MolHandle, use_chirality: bool) -> Vec<u8> {
    let config = chematic_fp::EcfpConfig {
        radius: 2,
        nbits: 2048,
        use_chirality,
        use_double_fold: false,
    };
    let fp = chematic_fp::ecfp(&mol.inner, &config);
    (0..256usize)
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

/// Like `ecfp6_bitvec` but with explicit chirality control.
///
/// When `use_chirality=true`, tetrahedral stereochemistry is included in the
/// initial atom hash, making enantiomers have different fingerprints.
/// When `false` (default), chirality is ignored.
#[wasm_bindgen]
pub fn ecfp6_bitvec_with_chirality(mol: &MolHandle, use_chirality: bool) -> Vec<u8> {
    let config = chematic_fp::EcfpConfig {
        radius: 3,
        nbits: 2048,
        use_chirality,
        use_double_fold: false,
    };
    let fp = chematic_fp::ecfp(&mol.inner, &config);
    (0..256usize)
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

/// Tanimoto similarity between two molecules using topological path fingerprints.
#[wasm_bindgen]
pub fn tanimoto_topo_path(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::tanimoto_topo_path(&a.inner, &b.inner)
}

// ---------------------------------------------------------------------------
// Sprint Q: IFG, VSA descriptors, Gasteiger charges, SA Score, Diversity
// ---------------------------------------------------------------------------

/// MinHash fingerprint (128 hashes) as JSON.
///
/// Returns `{"num_hashes":128,"hashes":[u64,...]}`.
/// Use `tanimoto_mhfp_smiles` for direct SMILES-to-SMILES similarity.
#[wasm_bindgen]
pub fn mhfp_hashes_json(mol: &MolHandle) -> String {
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return format!(
            r#"{{"error":"molecule too large (max {} atoms)"}}"#,
            WASM_MAX_ATOMS
        );
    }
    let fp = chematic_fp::mhfp(&mol.inner);
    let hs: Vec<String> = fp.hashes.iter().map(|h| h.to_string()).collect();
    format!(
        r#"{{"num_hashes":{},"hashes":[{}]}}"#,
        fp.num_hashes,
        hs.join(",")
    )
}

/// Tanimoto-like similarity between two SMILES via MHFP (MinHash Jaccard approximation).
#[wasm_bindgen]
pub fn tanimoto_mhfp_smiles(smi1: &str, smi2: &str) -> Result<f64, JsValue> {
    enforce_wasm_input_len("SMILES", smi1)?;
    enforce_wasm_input_len("SMILES", smi2)?;
    let m1 = chematic_smiles::parse(smi1).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let m2 = chematic_smiles::parse(smi2).map_err(|e| JsValue::from_str(&e.to_string()))?;
    enforce_wasm_molecule_size(&m1)?;
    enforce_wasm_molecule_size(&m2)?;
    Ok(chematic_fp::tanimoto_mhfp(&m1, &m2))
}

/// Tanimoto similarity between two molecules given only SMILES strings (ECFP4).
///
/// Returns a JS error on parse failure.
#[wasm_bindgen]
pub fn tanimoto_smiles(smiles1: &str, smiles2: &str) -> Result<f64, JsValue> {
    enforce_wasm_input_len("SMILES", smiles1)?;
    enforce_wasm_input_len("SMILES", smiles2)?;
    let m1 = chematic_smiles::parse(smiles1).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let m2 = chematic_smiles::parse(smiles2).map_err(|e| JsValue::from_str(&e.to_string()))?;
    enforce_wasm_molecule_size(&m1)?;
    enforce_wasm_molecule_size(&m2)?;
    Ok(chematic_fp::tanimoto_ecfp4(&m1, &m2))
}

/// MACCS 166-bit structural keys fingerprint as a byte array (21 bytes, LSB-first).
///
/// Bit `i` (0-indexed) corresponds to MACCS key `i+1`.
#[wasm_bindgen]
pub fn maccs_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::maccs(&mol.inner);
    // 166 bits → 21 bytes (bits 0–165; remaining bits in byte 20 are always 0).
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

/// Tanimoto similarity between `a` and `b` using MACCS 166-bit fingerprints.
#[wasm_bindgen]
pub fn tanimoto_maccs(a: &MolHandle, b: &MolHandle) -> f64 {
    let fa = chematic_fp::maccs(&a.inner);
    let fb = chematic_fp::maccs(&b.inner);
    fa.tanimoto(&fb)
}

/// ECFP6 (radius-3) fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
#[wasm_bindgen]
pub fn ecfp6_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::ecfp6(&mol.inner);
    (0..256usize)
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

/// Tanimoto similarity between `a` and `b` using ECFP6 fingerprints.
#[wasm_bindgen]
pub fn tanimoto_ecfp6(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::ecfp6(&a.inner).tanimoto(&chematic_fp::ecfp6(&b.inner))
}

/// Dice similarity between `a` and `b` using ECFP4 fingerprints.
#[wasm_bindgen]
pub fn dice_ecfp4(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::ecfp4(&a.inner).dice(&chematic_fp::ecfp4(&b.inner))
}

/// Dice similarity between `a` and `b` using MACCS 166-bit fingerprints.
#[wasm_bindgen]
pub fn dice_maccs(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::maccs(&a.inner).dice(&chematic_fp::maccs(&b.inner))
}

/// Select `n` maximally-diverse molecules (MaxMin algorithm, ECFP4 Tanimoto).
///
/// `smiles_json` — a JSON array of SMILES strings, e.g. `["CC","c1ccccc1","CCO"]`.
/// Returns a JSON array of 0-based indices into the input array.
/// Returns a JS error if any SMILES fails to parse (indices would otherwise shift).
#[wasm_bindgen]
pub fn maxmin_picks_ecfp4_json(smiles_json: &str, n: usize) -> Result<String, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let mols: Vec<chematic_core::Molecule> = smiles_list
        .iter()
        .map(|s| {
            let mol = chematic_smiles::parse(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect::<Result<_, _>>()?;
    let picks = chematic_chem::maxmin_picks(&mols, n, |a, b| {
        chematic_fp::ecfp4(a).tanimoto(&chematic_fp::ecfp4(b))
    });
    let parts: Vec<String> = picks.iter().map(|i| i.to_string()).collect();
    Ok(format!("[{}]", parts.join(",")))
}

/// Cluster molecules by structural similarity (Butina algorithm, ECFP4 Tanimoto).
///
/// `smiles_json` — a JSON array of SMILES strings.
/// `cutoff` — Tanimoto similarity threshold (0.0–1.0); molecules within this
///   distance of a cluster centre are assigned to that cluster.
/// Returns a JSON array of clusters, each cluster being an array of 0-based input indices.
/// Returns a JS error if any SMILES fails to parse.
#[wasm_bindgen]
pub fn butina_cluster_ecfp4_json(smiles_json: &str, cutoff: f64) -> Result<String, JsValue> {
    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let mols: Vec<chematic_core::Molecule> = smiles_list
        .iter()
        .map(|s| {
            let mol = chematic_smiles::parse(s).map_err(|e| JsValue::from_str(&e.to_string()))?;
            enforce_wasm_molecule_size(&mol)?;
            Ok::<_, JsValue>(mol)
        })
        .collect::<Result<_, _>>()?;
    let clusters = chematic_chem::butina_cluster(&mols, cutoff, |a, b| {
        chematic_fp::ecfp4(a).tanimoto(&chematic_fp::ecfp4(b))
    });
    let outer: Vec<String> = clusters
        .iter()
        .map(|cluster| {
            let inner: Vec<String> = cluster.iter().map(|i| i.to_string()).collect();
            format!("[{}]", inner.join(","))
        })
        .collect();
    Ok(format!("[{}]", outer.join(",")))
}

/// FCFP4 (pharmacophore, radius-2) fingerprint as a bit-packed byte vector (256 bytes).
#[wasm_bindgen]
pub fn fcfp4_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::fcfp4(&mol.inner);
    (0..256usize)
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

/// FCFP6 (pharmacophore, radius-3) fingerprint as a bit-packed byte vector (256 bytes).
#[wasm_bindgen]
pub fn fcfp6_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::fcfp6(&mol.inner);
    (0..256usize)
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

/// Dice similarity between `a` and `b` using ECFP6 fingerprints.
#[wasm_bindgen]
pub fn dice_ecfp6(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::ecfp6(&a.inner).dice(&chematic_fp::ecfp6(&b.inner))
}

/// AtomPair fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
#[wasm_bindgen]
pub fn atom_pair_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::atom_pair_fp(&mol.inner);
    (0..256usize)
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

/// Torsion fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
#[wasm_bindgen]
pub fn torsion_bitvec(mol: &MolHandle) -> Vec<u8> {
    let fp = chematic_fp::torsion_fp(&mol.inner);
    (0..256usize)
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

/// Tanimoto similarity between `a` and `b` using FCFP6 (radius-3 pharmacophore) fingerprints.
#[wasm_bindgen]
pub fn tanimoto_fcfp6(a: &MolHandle, b: &MolHandle) -> f64 {
    chematic_fp::fcfp6(&a.inner).tanimoto(&chematic_fp::fcfp6(&b.inner))
}

/// Compute a fingerprint bit-vector with configurable ECFP radius and bit width.
///
/// `radius` — Morgan radius (1 = ECFP2, 2 = ECFP4, 3 = ECFP6).
/// `nbits` — bit width; must be one of 256, 512, 1024, or 2048.
///   Returns a `Uint8Array` of `nbits/8` bytes.
///
/// The hash modulo is applied at fingerprint-generation time (`id % nbits`),
/// so no post-processing fold is needed.
/// Compute a custom ECFP (Extended Connectivity FingerPrint) with specified radius and bit count.
///
/// When `use_chirality=true`, tetrahedral stereochemistry is included in the initial
/// atom hash. When `false` (default), chirality is ignored.
#[wasm_bindgen]
pub fn ecfp_bitvec_custom(
    mol: &MolHandle,
    radius: u32,
    nbits: usize,
    use_chirality: bool,
) -> Vec<u8> {
    let nbits = match nbits {
        256 | 512 | 1024 | 2048 => nbits,
        _ => 2048,
    };
    let fp = chematic_fp::ecfp(
        &mol.inner,
        &chematic_fp::EcfpConfig {
            radius,
            nbits,
            use_chirality,
            use_double_fold: false,
        },
    );
    let byte_count = nbits / 8;
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

/// Find the k nearest neighbours of a query SMILES in a list of db SMILES.
///
/// `db_smiles_json`: JSON array of SMILES strings, e.g. `["CC","c1ccccc1"]`.
/// Returns JSON: `[{"index":0,"tanimoto":0.95},...]` sorted by descending Tanimoto.
/// Returns `"error:<msg>"` on parse failure.
#[wasm_bindgen]
pub fn nearest_neighbors_json(query_smiles: &str, db_smiles_json: &str, k: usize) -> String {
    if let Err(e) = enforce_wasm_input_len("query_smiles", query_smiles) {
        return format!(
            "error:{}",
            e.as_string()
                .unwrap_or_else(|| "input too large".to_string())
        );
    }
    let query = match chematic_smiles::parse(query_smiles) {
        Ok(m) => m,
        Err(e) => return format!("error:query parse failed: {e}"),
    };
    if let Err(e) = enforce_wasm_molecule_size(&query) {
        return format!(
            "error:{}",
            e.as_string()
                .unwrap_or_else(|| "molecule too large".to_string())
        );
    }

    let smiles_list = match parse_smiles_json_array(db_smiles_json) {
        Ok(values) => values,
        Err(e) => {
            return format!(
                "error:{}",
                e.as_string()
                    .unwrap_or_else(|| "db_smiles_json parse failed".to_string())
            );
        }
    };
    let mut db = Vec::with_capacity(smiles_list.len());
    let mut original_indices = Vec::with_capacity(smiles_list.len());
    for (idx, smiles) in smiles_list.iter().enumerate() {
        let mol = match chematic_smiles::parse(smiles) {
            Ok(mol) => mol,
            Err(e) => return format!("error:db parse failed at index {idx}: {e}"),
        };
        if let Err(e) = enforce_wasm_molecule_size(&mol) {
            return format!(
                "error:db molecule at index {idx}: {}",
                e.as_string()
                    .unwrap_or_else(|| "molecule too large".to_string())
            );
        }
        db.push(mol);
        original_indices.push(idx);
    }

    let results = chematic_fp::nearest_neighbors(&query, &db, k, chematic_fp::FpType::Ecfp4);
    let entries: Vec<String> = results
        .iter()
        .map(|(i, t)| {
            let original_idx = original_indices[*i];
            format!("{{\"index\":{original_idx},\"tanimoto\":{t:.6}}}")
        })
        .collect();
    format!("[{}]", entries.join(","))
}

/// Virtual screen a query SMILES against a database of SMILES using ECFP4 Tanimoto.
///
/// `db_smiles_json`: JSON array of SMILES strings (max 1024 via WASM_MAX_BATCH_ITEMS).
/// `k`: number of top hits to return; clamped to db size if larger.
///
/// Returns JSON: `{"results":[{"rank":1,"score":0.85,"smiles":"CCO","idx":42},...]}`.
/// Returns `"error:<msg>"` on any parse failure or oversized input.
#[wasm_bindgen]
pub fn virtual_screen_ecfp4_json(query_smi: &str, db_smiles_json: &str, k: u32) -> String {
    if let Err(e) = enforce_wasm_input_len("query_smi", query_smi) {
        return format!(
            "error:{}",
            e.as_string()
                .unwrap_or_else(|| "input too large".to_string())
        );
    }
    let query_mol = match chematic_smiles::parse(query_smi) {
        Ok(m) => m,
        Err(e) => return format!("error:query parse failed: {e}"),
    };
    if let Err(e) = enforce_wasm_molecule_size(&query_mol) {
        return format!(
            "error:{}",
            e.as_string()
                .unwrap_or_else(|| "molecule too large".to_string())
        );
    }
    let smiles_list = match parse_smiles_json_array(db_smiles_json) {
        Ok(v) => v,
        Err(e) => {
            return format!(
                "error:{}",
                e.as_string()
                    .unwrap_or_else(|| "db parse failed".to_string())
            );
        }
    };
    let mut db_fps = Vec::with_capacity(smiles_list.len());
    for (idx, smi) in smiles_list.iter().enumerate() {
        let mol = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(e) => return format!("error:db parse failed at index {idx}: {e}"),
        };
        if let Err(e) = enforce_wasm_molecule_size(&mol) {
            return format!(
                "error:db molecule at index {idx}: {}",
                e.as_string()
                    .unwrap_or_else(|| "molecule too large".to_string())
            );
        }
        db_fps.push(chematic_fp::ecfp4(&mol));
    }
    let query_fp = chematic_fp::ecfp4(&query_mol);
    let k = (k as usize).min(smiles_list.len());
    let hits = chematic_fp::top_k_similar(&query_fp, &db_fps, k);
    let entries: Vec<String> = hits
        .iter()
        .enumerate()
        .map(|(rank, (idx, score))| {
            let smi_escaped = smiles_list[*idx].replace('\\', "\\\\").replace('"', "\\\"");
            format!(
                r#"{{"rank":{},"score":{:.6},"smiles":"{}","idx":{}}}"#,
                rank + 1,
                score,
                smi_escaped,
                idx
            )
        })
        .collect();
    format!(r#"{{"results":[{}]}}"#, entries.join(","))
}

/// Compute ECFP4 Tanimoto similarity from one query SMILES to all db SMILES (dense output).
///
/// `db_smiles_json`: JSON array of SMILES strings (max 1024 via WASM_MAX_BATCH_ITEMS).
///
/// Returns a flat JSON array of f32 scores, one per db entry, e.g. `[0.12,0.0,0.85]`.
/// No zero-filtering: the length always equals the number of db entries.
/// Returns `"error:<msg>"` on parse failure or oversized input.
#[wasm_bindgen]
pub fn tanimoto_row_json(query_smi: &str, db_smiles_json: &str) -> String {
    if let Err(e) = enforce_wasm_input_len("query_smi", query_smi) {
        return format!(
            "error:{}",
            e.as_string()
                .unwrap_or_else(|| "input too large".to_string())
        );
    }
    let query_mol = match chematic_smiles::parse(query_smi) {
        Ok(m) => m,
        Err(e) => return format!("error:query parse failed: {e}"),
    };
    if let Err(e) = enforce_wasm_molecule_size(&query_mol) {
        return format!(
            "error:{}",
            e.as_string()
                .unwrap_or_else(|| "molecule too large".to_string())
        );
    }
    let smiles_list = match parse_smiles_json_array(db_smiles_json) {
        Ok(v) => v,
        Err(e) => {
            return format!(
                "error:{}",
                e.as_string()
                    .unwrap_or_else(|| "db parse failed".to_string())
            );
        }
    };
    let mut db_fps = Vec::with_capacity(smiles_list.len());
    for (idx, smi) in smiles_list.iter().enumerate() {
        let mol = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(e) => return format!("error:db parse failed at index {idx}: {e}"),
        };
        if let Err(e) = enforce_wasm_molecule_size(&mol) {
            return format!(
                "error:db molecule at index {idx}: {}",
                e.as_string()
                    .unwrap_or_else(|| "molecule too large".to_string())
            );
        }
        db_fps.push(chematic_fp::ecfp4(&mol));
    }
    let query_fp = chematic_fp::ecfp4(&query_mol);
    let scores = chematic_fp::tanimoto_slice(&query_fp, &db_fps);
    let parts: Vec<String> = scores.iter().map(|s| format!("{s:.6}")).collect();
    format!("[{}]", parts.join(","))
}

// ---------------------------------------------------------------------------
// RDKit-bit-exact Morgan/ECFP4 (promoted stable API)
// ---------------------------------------------------------------------------

/// RDKit-bit-exact ECFP4 (radius=2, 2048 bits, `useChirality=false`, `useBondTypes=true`,
/// RDKit's default atom invariant) as a bit-packed byte vector (256 bytes = 2048 bits).
///
/// Bit-for-bit identical to
/// `rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol)`
/// for every input this preprocessing handles. **Not** the same bits as `ecfp4_bitvec`
/// (that path uses chematic's own FNV-1a hash and is not RDKit-bit-compatible by
/// design -- the two are never silently interchanged).
///
/// Returns a JS error (its string carrying `RdkitMorganError`'s `Display` text, e.g.
/// `"rdkit-exact ecfp4: aromaticity: ..."`) if RDKit-parity aromaticity preprocessing
/// fails, never a silent fallback to `ecfp4_bitvec`'s Hückel-based engine -- the two
/// engines are not bit-compatible, so a silent substitution would look successful
/// while actually returning the wrong hash. See `docs/rfcs/ecfp4_bitexact_api_rfc.md`.
#[wasm_bindgen]
pub fn rdkit_ecfp4_bitvec(mol: &MolHandle) -> Result<Vec<u8>, JsValue> {
    let result = chematic_fp::rdkit_morgan_ecfp4_experimental(&mol.inner)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(bitvecn_to_bytes(&result.fingerprint.to_bitvecn()))
}

/// Same fingerprint as `rdkit_ecfp4_bitvec`, plus the raw (unfolded) data behind it, as
/// JSON: `{"fingerprint":[u8,...],"sparseCounts":{"rawId":count,...},
/// "rawBitInfo":{"rawId":[[atomIdx,radius],...],...},
/// "foldedBitInfo":{"bit":[[atomIdx,radius],...],...}}`.
///
/// Returns a JS error on the same preprocessing failures as `rdkit_ecfp4_bitvec`.
#[wasm_bindgen]
pub fn rdkit_ecfp4_detail_json(mol: &MolHandle) -> Result<String, JsValue> {
    let result = chematic_fp::rdkit_morgan_ecfp4_experimental(&mol.inner)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(rdkit_morgan_detail_json(
        &bitvecn_to_bytes(&result.fingerprint.to_bitvecn()),
        &sorted_pairs(&result.sparse_counts),
        &sorted_bit_info(&result.raw_bit_info),
        &sorted_bit_info(&result.folded_bit_info),
    ))
}

/// RDKit-bit-exact Morgan/ECFP fingerprint at a caller-chosen radius/bit-width, as a
/// bit-packed byte vector (`nbits / 8` bytes, LSB-first).
///
/// `radius` must be 0, 1, 2 (`rdkit_ecfp4_bitvec`'s ECFP4), or 3. `nbits` must be one
/// of 128, 256, 512, 1024, or 2048. Each of these 20 combinations is independently
/// re-verified against a live RDKit oracle (not assumed to generalize from
/// radius=2/2048 bits alone) -- see `validation/ecfp4_rdkit_stable_api_fixtures.json`.
/// An unsupported value returns a JS error rather than being silently coerced to the
/// nearest supported one.
#[wasm_bindgen]
pub fn rdkit_ecfp_config_bitvec(
    mol: &MolHandle,
    radius: u32,
    nbits: usize,
) -> Result<Vec<u8>, JsValue> {
    let config = wasm_rdkit_morgan_config(radius, nbits, false)?;
    let result = chematic_fp::rdkit_morgan_fingerprint(&mol.inner, &config)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(bitvecn_to_bytes(&result.fingerprint))
}

/// Same fingerprint as `rdkit_ecfp_config_bitvec`, plus the raw (unfolded) data -- see
/// `rdkit_ecfp4_detail_json` for the JSON shape (identical, generalized to this
/// function's `radius`/`nbits`).
#[wasm_bindgen]
pub fn rdkit_ecfp_config_detail_json(
    mol: &MolHandle,
    radius: u32,
    nbits: usize,
) -> Result<String, JsValue> {
    let config = wasm_rdkit_morgan_config(radius, nbits, false)?;
    let result = chematic_fp::rdkit_morgan_fingerprint(&mol.inner, &config)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(rdkit_morgan_detail_json(
        &bitvecn_to_bytes(&result.fingerprint),
        &sorted_pairs(&result.sparse_counts),
        &sorted_bit_info(&result.raw_bit_info),
        &sorted_bit_info(&result.folded_bit_info),
    ))
}

/// Chirality-enabled variant of [`rdkit_ecfp_config_bitvec`]. E/Z bond stereo
/// is not included by this API yet.
#[wasm_bindgen]
pub fn rdkit_ecfp_config_chiral_bitvec(
    mol: &MolHandle,
    radius: u32,
    nbits: usize,
) -> Result<Vec<u8>, JsValue> {
    let config = wasm_rdkit_morgan_config(radius, nbits, true)?;
    let result = chematic_fp::rdkit_morgan_fingerprint(&mol.inner, &config)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(bitvecn_to_bytes(&result.fingerprint))
}

/// Chirality-enabled variant of [`rdkit_ecfp_config_detail_json`].
#[wasm_bindgen]
pub fn rdkit_ecfp_config_chiral_detail_json(
    mol: &MolHandle,
    radius: u32,
    nbits: usize,
) -> Result<String, JsValue> {
    let config = wasm_rdkit_morgan_config(radius, nbits, true)?;
    let result = chematic_fp::rdkit_morgan_fingerprint(&mol.inner, &config)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(rdkit_morgan_detail_json(
        &bitvecn_to_bytes(&result.fingerprint),
        &sorted_pairs(&result.sparse_counts),
        &sorted_bit_info(&result.raw_bit_info),
        &sorted_bit_info(&result.folded_bit_info),
    ))
}

/// Maps the plain `radius`/`nbits` integers the WASM API accepts into chematic-fp's
/// closed `RdkitMorganConfig` enums. An unsupported value is an explicit JS error, not
/// a guessed/coerced conversion.
fn wasm_rdkit_morgan_config(
    radius: u32,
    nbits: usize,
    include_chirality: bool,
) -> Result<chematic_fp::RdkitMorganConfig, JsValue> {
    let radius = match radius {
        0 => chematic_fp::RdkitMorganRadius::R0,
        1 => chematic_fp::RdkitMorganRadius::R1,
        2 => chematic_fp::RdkitMorganRadius::R2,
        3 => chematic_fp::RdkitMorganRadius::R3,
        other => {
            return Err(JsValue::from_str(&format!(
                "unsupported radius {other} -- only 0, 1, 2, or 3 are supported"
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
            return Err(JsValue::from_str(&format!(
                "unsupported nbits {other} -- only 128, 256, 512, 1024, or 2048 are supported"
            )));
        }
    };
    Ok(chematic_fp::RdkitMorganConfig {
        radius,
        fp_size,
        include_chirality,
    })
}

fn bitvecn_to_bytes(fp: &chematic_fp::BitVecN) -> Vec<u8> {
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

/// Sorts a hash map's entries by key -- generic over the hasher (`S`) so this works on
/// both `std::collections::HashMap` and chematic-fp's `rustc_hash::FxHashMap` return
/// types without this crate needing its own `rustc-hash` dependency just to name the
/// type.
fn sorted_pairs<K: Ord + Copy, V: Copy, S>(
    map: &std::collections::HashMap<K, V, S>,
) -> Vec<(K, V)> {
    let mut v: Vec<(K, V)> = map.iter().map(|(&k, &v)| (k, v)).collect();
    v.sort_unstable_by_key(|(k, _)| *k);
    v
}

/// Same idea as [`sorted_pairs`], for the `{raw_id/bit -> [(atom_idx, radius), ...]}`
/// bit-info maps -- sorts by key, and sorts each value's `(atom_idx, radius)` pairs too
/// (insertion order isn't part of the API contract).
fn sorted_bit_info<K: Ord + Copy, S>(
    map: &std::collections::HashMap<K, Vec<(u32, u32)>, S>,
) -> Vec<(K, Vec<(u32, u32)>)> {
    let mut v: Vec<(K, Vec<(u32, u32)>)> = map
        .iter()
        .map(|(&k, pairs)| {
            let mut p = pairs.clone();
            p.sort_unstable();
            (k, p)
        })
        .collect();
    v.sort_unstable_by_key(|(k, _)| *k);
    v
}

fn rdkit_morgan_detail_json(
    fingerprint_bytes: &[u8],
    sparse_counts: &[(u32, u32)],
    raw_bit_info: &[(u32, Vec<(u32, u32)>)],
    folded_bit_info: &[(usize, Vec<(u32, u32)>)],
) -> String {
    let fp_json: Vec<String> = fingerprint_bytes.iter().map(|b| b.to_string()).collect();
    let sparse_json: Vec<String> = sparse_counts
        .iter()
        .map(|(k, v)| format!("\"{k}\":{v}"))
        .collect();
    let raw_json: Vec<String> = raw_bit_info
        .iter()
        .map(|(k, v)| format!("\"{k}\":{}", bit_info_pairs_json(v)))
        .collect();
    let folded_json: Vec<String> = folded_bit_info
        .iter()
        .map(|(k, v)| format!("\"{k}\":{}", bit_info_pairs_json(v)))
        .collect();

    format!(
        r#"{{"fingerprint":[{}],"sparseCounts":{{{}}},"rawBitInfo":{{{}}},"foldedBitInfo":{{{}}}}}"#,
        fp_json.join(","),
        sparse_json.join(","),
        raw_json.join(","),
        folded_json.join(","),
    )
}

/// Assumes `pairs` is already sorted (both [`sorted_pairs`]/[`sorted_bit_info`] callers
/// guarantee this).
fn bit_info_pairs_json(pairs: &[(u32, u32)]) -> String {
    let parts: Vec<String> = pairs.iter().map(|(a, r)| format!("[{a},{r}]")).collect();
    format!("[{}]", parts.join(","))
}

// ---------------------------------------------------------------------------
// MOL2 I/O
// ---------------------------------------------------------------------------
