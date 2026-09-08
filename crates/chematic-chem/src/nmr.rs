//! Vendor-neutral NMR interchange contract; it does not parse or predict spectra.

use serde::{Deserialize, Serialize};
use serde_json::Value;

pub const NMR_SCHEMA_VERSION: u32 = 1;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct NmrSpectrum {
    pub schema_version: u32,
    pub nucleus: String,
    pub solvent: Option<String>,
    pub frequency_mhz: Option<f64>,
    pub temperature_k: Option<f64>,
    pub source: Option<String>,
    pub peaks: Vec<NmrPeak>,
    pub normalization: NmrNormalization,
    pub vendor_metadata: Value,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct NmrPeak {
    pub shift: f64,
    pub intensity: f64,
    pub assignment: Option<String>,
    pub multiplicity: Option<String>,
    pub coupling_hz: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum NmrNormalization {
    Raw,
    MaxOne,
    SumOne,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct NmrLimits {
    pub max_peaks: usize,
    pub max_serialized_bytes: usize,
}

impl Default for NmrLimits {
    fn default() -> Self {
        Self {
            max_peaks: 100_000,
            max_serialized_bytes: 1_000_000,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct NmrDiagnostic {
    pub code: String,
    pub path: String,
    pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct NmrValidation {
    pub schema_version: u32,
    pub valid: bool,
    pub diagnostics: Vec<NmrDiagnostic>,
}

impl NmrSpectrum {
    pub fn validate(&self, limits: &NmrLimits) -> NmrValidation {
        let mut diagnostics = Vec::new();
        let add =
            |diagnostics: &mut Vec<NmrDiagnostic>, code: &str, path: String, message: &str| {
                diagnostics.push(NmrDiagnostic {
                    code: code.into(),
                    path,
                    message: message.into(),
                });
            };
        if self.schema_version != NMR_SCHEMA_VERSION {
            add(
                &mut diagnostics,
                "schema_version",
                "schema_version".into(),
                "unsupported schema version",
            );
        }
        if self.nucleus.trim().is_empty() {
            add(
                &mut diagnostics,
                "missing_nucleus",
                "nucleus".into(),
                "nucleus is required",
            );
        }
        if self.peaks.len() > limits.max_peaks {
            add(
                &mut diagnostics,
                "peak_limit",
                "peaks".into(),
                "peak count exceeds limit",
            );
        }
        for (i, peak) in self.peaks.iter().enumerate() {
            if !peak.shift.is_finite() {
                add(
                    &mut diagnostics,
                    "non_finite_shift",
                    format!("peaks[{i}].shift"),
                    "chemical shift must be finite",
                );
            }
            if !peak.intensity.is_finite() || peak.intensity < 0.0 {
                add(
                    &mut diagnostics,
                    "invalid_intensity",
                    format!("peaks[{i}].intensity"),
                    "intensity must be finite and non-negative",
                );
            }
            for (j, coupling) in peak.coupling_hz.iter().enumerate() {
                if !coupling.is_finite() || *coupling < 0.0 {
                    add(
                        &mut diagnostics,
                        "invalid_coupling",
                        format!("peaks[{i}].coupling_hz[{j}]"),
                        "coupling must be finite and non-negative",
                    );
                }
            }
        }
        for (path, value) in [
            ("frequency_mhz", self.frequency_mhz),
            ("temperature_k", self.temperature_k),
        ] {
            if value.is_some_and(|v| !v.is_finite() || v < 0.0) {
                add(
                    &mut diagnostics,
                    "invalid_experiment_value",
                    path.into(),
                    "experiment value must be finite and non-negative",
                );
            }
        }
        NmrValidation {
            schema_version: NMR_SCHEMA_VERSION,
            valid: diagnostics.is_empty(),
            diagnostics,
        }
    }
}

pub fn validate_nmr_json(json: &str, limits: &NmrLimits) -> NmrValidation {
    if json.len() > limits.max_serialized_bytes {
        return NmrValidation {
            schema_version: NMR_SCHEMA_VERSION,
            valid: false,
            diagnostics: vec![NmrDiagnostic {
                code: "serialized_size_limit".into(),
                path: "$".into(),
                message: "serialized spectrum exceeds limit".into(),
            }],
        };
    }
    match serde_json::from_str::<NmrSpectrum>(json) {
        Ok(spectrum) => spectrum.validate(limits),
        Err(error) => NmrValidation {
            schema_version: NMR_SCHEMA_VERSION,
            valid: false,
            diagnostics: vec![NmrDiagnostic {
                code: "invalid_json".into(),
                path: "$".into(),
                message: error.to_string(),
            }],
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn clean() -> NmrSpectrum {
        NmrSpectrum {
            schema_version: 1,
            nucleus: "1H".into(),
            solvent: Some("CDCl3".into()),
            frequency_mhz: Some(400.0),
            temperature_k: None,
            source: Some("fixture".into()),
            peaks: vec![NmrPeak {
                shift: 7.26,
                intensity: 1.0,
                assignment: None,
                multiplicity: Some("s".into()),
                coupling_hz: vec![],
            }],
            normalization: NmrNormalization::MaxOne,
            vendor_metadata: serde_json::json!({"vendor":{"opaque":true}}),
        }
    }
    #[test]
    fn valid_round_trip_and_opaque_metadata() {
        let s = clean();
        let json = serde_json::to_string(&s).unwrap();
        assert!(validate_nmr_json(&json, &NmrLimits::default()).valid);
        assert_eq!(serde_json::from_str::<NmrSpectrum>(&json).unwrap(), s);
    }
    #[test]
    fn malformed_values_have_stable_paths() {
        let mut s = clean();
        s.peaks[0].intensity = f64::NAN;
        let r = s.validate(&NmrLimits::default());
        assert_eq!(r.diagnostics[0].path, "peaks[0].intensity");
    }
}
