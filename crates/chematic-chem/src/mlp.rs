//! ECFP4-based linear model for aqueous solubility prediction (logS).
//!
//! Uses a Ridge regression trained on the Delaney ESOL dataset (1128 molecules).
//! Weights are stored as binary files embedded at compile time via `include_bytes!`.
//!
//! ## Architecture
//!
//! Input: ECFP4 (2048 bits as f32)  →  Linear(2048 → 1)
//! Output: logS (log mol/L)
//!
//! ## Performance (Delaney ESOL dataset)
//!
//! - Train RMSE: 0.447  R²: 0.955
//! - 5-fold CV  R²: 0.626 ± 0.030
//!
//! ## Binary size
//!
//! - `mlp_w1.bin`: 2048 × 4 bytes = 8 192 bytes (~8 KB)
//! - `mlp_b1.bin`: 4 bytes
//!
//! ## How to retrain
//!
//! ```sh
//! python scripts/train_solubility_mlp.py --csv data/delaney.csv
//! ```
//!
//! ## Placeholder mode
//!
//! When the binary files are absent (feature `trained-solubility-mlp` disabled or
//! data files not found), `mlp_solubility` falls back to the Delaney ESOL linear
//! regression from `crate::esol`.

use chematic_core::Molecule;
#[cfg(feature = "trained-solubility-mlp")]
use chematic_fp::ecfp4;

// ---------------------------------------------------------------------------
// Embedded trained weights (compile-time include)
// ---------------------------------------------------------------------------

// Reject big-endian compilation immediately: the binary weight files are
// written in little-endian byte order and there is no runtime endian-swap.
#[cfg(all(feature = "trained-solubility-mlp", not(target_endian = "little")))]
compile_error!(
    "trained-solubility-mlp uses little-endian binary weight files; \
     big-endian compilation targets are not supported"
);

/// `true` when trained weight files are embedded.
pub const MLP_SOLUBILITY_TRAINED: bool = cfg!(feature = "trained-solubility-mlp");

// The weights are conditionally compiled when data files are present.
// If the files do not exist the fallback path (esol_solubility) is used.
#[cfg(feature = "trained-solubility-mlp")]
const W1_BYTES: &[u8] = include_bytes!("../../../data/mlp_w1.bin");
#[cfg(feature = "trained-solubility-mlp")]
const B1_BYTES: &[u8] = include_bytes!("../../../data/mlp_b1.bin");

// Cache parsed weights so that repeated calls to `mlp_solubility` (e.g. in a
// virtual-screening loop over 100 k molecules) pay the ~8 KB deserialization
// cost only once rather than on every invocation.
#[cfg(feature = "trained-solubility-mlp")]
static WEIGHTS: std::sync::OnceLock<(Vec<f32>, Vec<f32>)> = std::sync::OnceLock::new();

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Predict aqueous solubility (logS, log mol/L) using an ECFP4-based linear model.
///
/// When the `trained-solubility-mlp` Cargo feature is enabled and the binary weight
/// files (`data/mlp_w1.bin`, `data/mlp_b1.bin`) are present at compile time, this
/// runs a Ridge-regression forward pass on the molecule's ECFP4 fingerprint.
///
/// Otherwise falls back transparently to `esol_solubility` (Delaney 2004 formula).
///
/// Typical logS range: −10 (very insoluble) to 0 (fully miscible).
///
/// Performance with trained weights (Delaney ESOL, 5-fold CV): R² ≈ 0.63
pub fn mlp_solubility(mol: &Molecule) -> f64 {
    #[cfg(feature = "trained-solubility-mlp")]
    {
        let fp = ecfp4(mol);
        let (w, b) = WEIGHTS.get_or_init(|| (f32_from_bytes(W1_BYTES), f32_from_bytes(B1_BYTES)));

        if w.len() != 2048 || b.is_empty() {
            return crate::esol::esol_solubility(mol);
        }

        let mut acc = b[0];
        for (i, wi) in w.iter().enumerate().take(2048) {
            if fp.get(i) {
                acc += wi;
            }
        }
        acc as f64
    }
    #[cfg(not(feature = "trained-solubility-mlp"))]
    {
        crate::esol::esol_solubility(mol)
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Interpret a `&[u8]` byte slice as little-endian `f32` values.
///
/// Callers on big-endian targets are rejected at compile time via the
/// `compile_error!` above the feature guard.
#[cfg(feature = "trained-solubility-mlp")]
fn f32_from_bytes(bytes: &[u8]) -> Vec<f32> {
    bytes
        .as_chunks::<4>()
        .0
        .iter()
        .map(|c| f32::from_le_bytes(*c))
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn logs(smiles: &str) -> f64 {
        mlp_solubility(&parse(smiles).unwrap())
    }

    #[test]
    fn test_output_finite_and_plausible() {
        for smi in ["c1ccccc1", "CCO", "CC(=O)Oc1ccccc1C(=O)O", "CCCCCCCC"] {
            let s = logs(smi);
            assert!(s.is_finite(), "logS must be finite for {smi}");
            assert!(
                s > -15.0 && s < 5.0,
                "logS={s:.2} out of plausible range for {smi}"
            );
        }
    }

    #[test]
    fn test_water_more_soluble_than_octane() {
        // Water should be predicted more soluble than octane regardless of model mode.
        let water = logs("O");
        let octane = logs("CCCCCCCC");
        assert!(
            water > octane,
            "water ({water:.2}) should be more soluble than octane ({octane:.2})"
        );
    }

    #[test]
    fn test_placeholder_matches_esol() {
        // In placeholder mode (no trained-solubility-mlp feature), output == ESOL.
        if !MLP_SOLUBILITY_TRAINED {
            let mol = parse("c1ccccc1").unwrap();
            let esol = crate::esol::esol_solubility(&mol);
            let mlp = mlp_solubility(&mol);
            assert!(
                (mlp - esol).abs() < 1e-9,
                "placeholder must match ESOL: mlp={mlp:.4} esol={esol:.4}"
            );
        }
    }
}
