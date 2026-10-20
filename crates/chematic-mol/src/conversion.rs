//! Checked conversion: MOL2/CJSON source retention and bounded graph projection.
//!
//! Reports account for source information before graph projection. Candidate
//! verification covers the declared graph fields only, not MOL2 chemistry in
//! opaque extensions or stereochemistry inferred from geometry.

use crate::{Mol2Error, Mol2ParseLimits, parse_mol2_source_record_with_limits};
use chematic_core::{AtomIdx, Molecule};
use serde::Serialize;
use std::collections::HashMap;

mod cjson_mol2;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ConversionPolicy {
    Strict,
    AllowLoss,
    PreserveSource,
}

impl std::str::FromStr for ConversionPolicy {
    type Err = ConversionError;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "strict" => Ok(Self::Strict),
            "allow_loss" => Ok(Self::AllowLoss),
            "preserve_source" => Ok(Self::PreserveSource),
            _ => Err(ConversionError::Unsupported(format!(
                "unknown conversion policy: {value}"
            ))),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum PreservationStatus {
    Preserved,
    Inferred,
    OpaqueRetained,
    Dropped,
    Unverified,
}

#[derive(Debug, Clone, Serialize)]
pub struct ConversionObservation {
    /// A stable field-group identifier, relative to source record zero.
    pub field: &'static str,
    pub status: PreservationStatus,
    /// `source_retention`, `projection`, `write`, or `candidate_verification`.
    pub stage: &'static str,
    pub reason: &'static str,
}

#[derive(Debug, Clone, Serialize)]
pub struct ConversionAtomMapping {
    pub source_id: u32,
    pub source_index: usize,
    pub candidate_index: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct ConversionBondMapping {
    pub source_id: u32,
    pub source_index: usize,
    pub candidate_index: usize,
}

/// JSON envelope shared by CLI/Python/WASM. Mapping indices refer to the
/// verified *candidate*, which strict policy may refuse to emit.
#[derive(Debug, Clone, Serialize)]
pub struct ConversionReport {
    pub schema: &'static str,
    pub input_format: &'static str,
    pub output_format: &'static str,
    pub policy: ConversionPolicy,
    /// `completed`, `completed_with_loss`, or `refused`.
    pub status: &'static str,
    pub source_record_index: usize,
    pub verification_scope: &'static str,
    pub candidate_graph_verified: bool,
    pub output: Option<String>,
    pub observations: Vec<ConversionObservation>,
    pub candidate_atom_mapping: Vec<ConversionAtomMapping>,
    pub candidate_bond_mapping: Vec<ConversionBondMapping>,
}

#[derive(Debug)]
pub enum ConversionError {
    Unsupported(String),
    Parse(Mol2Error),
    CjsonParse(crate::CjsonError),
    ResourceLimit(&'static str),
}

impl std::fmt::Display for ConversionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Unsupported(message) => f.write_str(message),
            Self::Parse(error) => error.fmt(f),
            Self::CjsonParse(error) => error.fmt(f),
            Self::ResourceLimit(resource) => write!(f, "conversion exceeds {resource} limit"),
        }
    }
}
impl std::error::Error for ConversionError {}

/// Bounded defaults; existing graph-only conversion APIs are unchanged.
#[derive(Debug, Clone, Copy)]
pub struct ConversionLimits {
    /// Shared input byte/atom/bond caps, also applied to CJSON input.
    pub mol2: Mol2ParseLimits,
    /// Applies to the molecule output and the serialized report separately.
    pub max_output_bytes: usize,
}
impl Default for ConversionLimits {
    fn default() -> Self {
        Self {
            mol2: Mol2ParseLimits {
                max_input_bytes: 16 * 1024 * 1024,
                max_atoms: 10_000,
                ..Default::default()
            },
            max_output_bytes: 16 * 1024 * 1024,
        }
    }
}

pub fn convert_with_report(
    text: &str,
    input_format: &str,
    output_format: &str,
    policy: ConversionPolicy,
) -> Result<ConversionReport, ConversionError> {
    convert_with_report_with_limits(
        text,
        input_format,
        output_format,
        policy,
        &ConversionLimits::default(),
    )
}

/// Return the shared JSON schema, including refusal reports. Invalid inputs,
/// unsupported routes/policies and resource limits are errors, not refusals.
/// The JSON size limit includes escaping and mapping/observation overhead.
pub fn convert_with_report_json(
    text: &str,
    input_format: &str,
    output_format: &str,
    policy: &str,
    limits: &ConversionLimits,
) -> Result<String, ConversionError> {
    let report = convert_with_report_with_limits(
        text,
        input_format,
        output_format,
        policy.parse()?,
        limits,
    )?;
    let json =
        serde_json::to_string(&report).expect("conversion report contains no fallible values");
    if json.len() > limits.max_output_bytes {
        return Err(ConversionError::ResourceLimit("report bytes"));
    }
    Ok(json)
}

pub fn convert_with_report_with_limits(
    text: &str,
    input_format: &str,
    output_format: &str,
    policy: ConversionPolicy,
    limits: &ConversionLimits,
) -> Result<ConversionReport, ConversionError> {
    let source = input_format
        .trim()
        .trim_start_matches('.')
        .to_ascii_lowercase();
    let target = output_format
        .trim()
        .trim_start_matches('.')
        .to_ascii_lowercase();
    let supported = match source.as_str() {
        "mol2" => matches!(target.as_str(), "mol2" | "smiles" | "smi" | "cjson"),
        "cjson" => matches!(target.as_str(), "cjson" | "smiles" | "smi" | "mol2"),
        _ => false,
    };
    if !supported {
        return Err(ConversionError::Unsupported(
            "checked conversion supports MOL2 → MOL2/SMILES/CJSON and CJSON → CJSON/SMILES/MOL2"
                .into(),
        ));
    }
    if policy == ConversionPolicy::PreserveSource && target != source {
        return Err(ConversionError::Unsupported(
            "preserve_source requires same-format conversion".into(),
        ));
    }
    if source == "cjson" {
        return convert_cjson(text, &target, policy, limits);
    }
    let record =
        parse_mol2_source_record_with_limits(text, &limits.mol2).map_err(ConversionError::Parse)?;
    let mut report = ConversionReport {
        schema: "chematic.conversion-report.v1",
        input_format: "mol2",
        output_format: match target.as_str() {
            "mol2" => "mol2",
            "cjson" => "cjson",
            _ => "smiles",
        },
        policy,
        status: "completed",
        source_record_index: 0,
        verification_scope: "exact_source_bytes",
        candidate_graph_verified: false,
        output: None,
        observations: Vec::new(),
        candidate_atom_mapping: Vec::new(),
        candidate_bond_mapping: Vec::new(),
    };
    let mut note = |field, status, stage, reason| {
        report.observations.push(ConversionObservation {
            field,
            status,
            stage,
            reason,
        })
    };
    if target == "mol2" {
        // Immutable record: all bytes, including unparsed columns and sections,
        // retain their source representation. No editable chemistry is claimed.
        let output = crate::write_mol2_source_record(&record);
        if output.len() > limits.max_output_bytes {
            return Err(ConversionError::ResourceLimit("output bytes"));
        }
        if output == text {
            note(
                "source_text",
                PreservationStatus::Preserved,
                "source_retention",
                "Exact input bytes retained after bounded structural validation; no edits or chemical interpretation promised.",
            );
        }
        note(
            "source_semantics",
            PreservationStatus::OpaqueRetained,
            "source_retention",
            "Source representations retained, including any unparsed status columns, formal-charge extensions and geometry stereo; their chemistry is not verified.",
        );
        report.output = Some(output);
        return Ok(report);
    }

    note(
        "source_text",
        PreservationStatus::Dropped,
        "write",
        "Candidate output omits MOL2 layout, comments, header flags and any unparsed optional columns.",
    );
    let partial_values: Option<Vec<f64>> =
        record.atoms().iter().map(|a| a.partial_charge).collect();
    let partial_method = record
        .charge_type()
        .filter(|method| !method.is_empty() && *method != "NO_CHARGES");
    let emit_partial = target == "cjson"
        && !record.atoms().is_empty()
        && partial_values.is_some()
        && partial_method.is_some();
    note(
        "molecule_metadata",
        PreservationStatus::Dropped,
        "write",
        if target == "cjson" {
            "Molecule type and other header flags are not emitted. Parsed name and any complete declared partial-charge dataset are checked separately."
        } else {
            "Molecule name, molecule type and charge-method header are not emitted."
        },
    );
    note(
        "atom_identity",
        PreservationStatus::Dropped,
        "write",
        "Source atom IDs and names are omitted; candidate correspondence is reported separately.",
    );
    note(
        "atom_types",
        PreservationStatus::Dropped,
        "projection",
        "Source MOL2 atom types are reduced to the bounded graph projection; the original typing model is not emitted.",
    );
    if record
        .atoms()
        .iter()
        .any(|atom| atom.substructure_id.is_some() || atom.substructure_name.is_some())
    {
        note(
            "substructures",
            PreservationStatus::Dropped,
            "projection",
            "Source atom substructure ownership and names are not emitted.",
        );
    }
    if record
        .atoms()
        .iter()
        .any(|atom| atom.partial_charge.is_some())
        && !emit_partial
    {
        note(
            "partial_charges",
            PreservationStatus::Dropped,
            "write",
            if target == "cjson" {
                "Partial charges require a non-NO_CHARGES method label and one value for every atom; missing values are not filled and undeclared charge columns are not labeled as a computed dataset."
            } else {
                "Source charge-column values are not formal charges; the candidate writer does not emit them as partial charges."
            },
        );
    }
    if target != "cjson" && !record.coordinates().is_empty() {
        note(
            "coordinates",
            PreservationStatus::Dropped,
            "write",
            "Plain SMILES has no source coordinates; geometry stereochemistry is not interpreted.",
        );
    }
    if !record.opaque_sections().is_empty() {
        note(
            "opaque_sections",
            PreservationStatus::Dropped,
            "projection",
            "Additional MOL2 sections are not interpreted and are not emitted.",
        );
    }
    if !record.bonds().is_empty() {
        note(
            "bond_records",
            PreservationStatus::Dropped,
            "projection",
            "Source bond IDs and MOL2 tokens such as am, du and nc are not retained as source metadata.",
        );
    }
    note(
        "formal_charge_interpretation",
        PreservationStatus::Unverified,
        "projection",
        "Explicit UNITY charge attributes are projected by source atom ID, with N.4 implying +1 only when no explicit charge exists. Other atom-type and extension rules are not interpreted. Preserving the projected charge is not verification of source chemistry.",
    );
    note(
        "stereochemistry",
        PreservationStatus::Unverified,
        "projection",
        "Coordinate-derived and extension stereochemistry are not interpreted or verified.",
    );
    let (candidate, order) = if target == "cjson" {
        // This writer emits atom arrays in source order and retains all three
        // coordinate components; no coordinates or hydrogens are generated.
        let mut json: serde_json::Value =
            serde_json::from_str(&crate::write_cjson(record.molecule(), record.coordinates()))
                .expect("CJSON writer emits valid JSON");
        if let Some(name) = record.name() {
            json["properties"] = serde_json::json!({"name": name});
        }
        if emit_partial {
            let mut datasets = serde_json::Map::new();
            datasets.insert(
                partial_method.unwrap().to_owned(),
                serde_json::json!(partial_values.as_ref().unwrap()),
            );
            json["partialCharges"] = serde_json::Value::Object(datasets);
        }
        (
            json.to_string(),
            (0..record.molecule().atom_count())
                .map(|i| AtomIdx(i as u32))
                .collect(),
        )
    } else {
        note(
            "implicit_hydrogens",
            PreservationStatus::Inferred,
            "write",
            "SMILES hydrogen counts follow the graph valence model, not an independently verified MOL2 hydrogen assignment.",
        );
        // Current SMILES writer uses recursive traversal. Bound this route
        // independently of reader limits until an iterative writer is available.
        if record.molecule().atom_count() > 512 {
            return Err(ConversionError::ResourceLimit(
                "SMILES projection atoms (512)",
            ));
        }
        chematic_smiles::write_with_atom_order(record.molecule())
    };
    if candidate.len() > limits.max_output_bytes {
        return Err(ConversionError::ResourceLimit("output bytes"));
    }
    let parsed = if target == "cjson" {
        crate::parse_cjson_with_limits(
            &candidate,
            &crate::CjsonParseLimits {
                max_input_bytes: limits.max_output_bytes,
                max_atoms: limits.mol2.max_atoms,
                max_bonds: limits.mol2.max_bonds,
                ..Default::default()
            },
        )
        .ok()
    } else {
        chematic_smiles::parse(&candidate)
            .ok()
            .map(|mol| (mol, Vec::new()))
    };
    let verified = parsed
        .as_ref()
        .and_then(|(mol, _)| verify_graph(record.molecule(), mol, &order));
    report.verification_scope = if target == "cjson" {
        "candidate_elements_aromaticity_projected_charges_bond_orders_coordinates_emitted_metadata"
    } else {
        "candidate_elements_aromaticity_projected_charges_bond_orders"
    };
    report.candidate_graph_verified = verified.is_some();
    report.observations.push(ConversionObservation {
        field: "projected_graph",
        status: if verified.is_some() { PreservationStatus::Preserved } else { PreservationStatus::Unverified },
        stage: "candidate_verification",
        reason: if verified.is_some() { "Candidate output reparsed; explicit atom counts, mapped elements/aromaticity/projected formal charges and mapped bond orders agree. Source chemical interpretation and inferred hydrogens are outside this check." } else { "Candidate graph could not be verified against the source projection; no correspondence is claimed." },
    });
    if target == "cjson" {
        let coords_verified = verified.is_some()
            && parsed.as_ref().is_some_and(|(_, coords)| {
                coords.len() == record.coordinates().len()
                    && order.iter().enumerate().all(|(candidate_index, source)| {
                        let before = record.coordinates()[source.0 as usize];
                        let after = coords[candidate_index];
                        before.0.to_bits() == after.0.to_bits()
                            && before.1.to_bits() == after.1.to_bits()
                            && before.2.to_bits() == after.2.to_bits()
                    })
            });
        report.observations.push(ConversionObservation {
            field: "coordinates",
            status: if coords_verified { PreservationStatus::Preserved } else { PreservationStatus::Unverified },
            stage: "candidate_verification",
            reason: if coords_verified { "Reparsed candidate coordinates agree bit for bit with source f64 x/y/z values through verified atom correspondence, including signed zero. Original decimal spelling and physical dimensionality are outside this check." } else { "Source-to-candidate graph correspondence or exact f64 coordinate equality could not be verified." },
        });
        // Reparse the emitted metadata independently of the writer's Value.
        // Partial charges are separate from formal charges and are never
        // computed here. The method is the source declaration only.
        let metadata: serde_json::Value =
            serde_json::from_str(&candidate).expect("CJSON candidate already parsed");
        if let Some(name) = record.name() {
            let equal = metadata
                .pointer("/properties/name")
                .and_then(serde_json::Value::as_str)
                == Some(name);
            report.observations.push(ConversionObservation {
                field: "molecule_name",
                status: if equal { PreservationStatus::Preserved } else { PreservationStatus::Unverified },
                stage: "candidate_verification",
                reason: "Reparsed properties.name equals the parsed (trimmed) source molecule name; original header whitespace is accounted for by source_text loss.",
            });
        }
        if emit_partial {
            let method = partial_method.unwrap();
            let datasets = metadata
                .get("partialCharges")
                .and_then(serde_json::Value::as_object);
            let values = datasets
                .and_then(|sets| sets.get(method))
                .and_then(serde_json::Value::as_array);
            let values_equal = verified.is_some()
                && values.is_some_and(|values| {
                    values.len() == record.atoms().len()
                        && order.iter().enumerate().all(|(candidate_index, source)| {
                            values[candidate_index].as_f64().is_some_and(|after| {
                                after.to_bits()
                                    == partial_values.as_ref().unwrap()[source.0 as usize].to_bits()
                            })
                        })
                });
            report.observations.push(ConversionObservation {
                field: "partial_charges",
                status: if values_equal { PreservationStatus::Preserved } else { PreservationStatus::Unverified },
                stage: "candidate_verification",
                reason: "Partial-charge preservation requires complete source values, reparsed bitwise equality and verified atom correspondence; no rounding into formal charges or filling missing values.",
            });
            report.observations.push(ConversionObservation {
                field: "charge_method",
                status: if datasets.is_some_and(|sets| sets.len() == 1 && sets.contains_key(method)) { PreservationStatus::Preserved } else { PreservationStatus::Unverified },
                stage: "candidate_verification",
                reason: "Source charge-method label retained verbatim as the partialCharges dataset key; the declared method is not recomputed or independently validated.",
            });
        }
        report.observations.push(ConversionObservation {
            field: "projected_formal_charges",
            status: if verified.is_some() { PreservationStatus::Preserved } else { PreservationStatus::Unverified },
            stage: "candidate_verification",
            reason: "CJSON formal charges are checked as part of mapped graph verification. MOL2 projection uses explicit UNITY charge attributes, then N.4 as a fallback; partial-charge preservation is checked separately. Source valence and other charge rules are not verified.",
        });
    }
    if let Some(bond_order) = verified {
        report.candidate_atom_mapping = order
            .iter()
            .enumerate()
            .map(|(candidate_index, source)| ConversionAtomMapping {
                source_id: record.atoms()[source.0 as usize].source_id,
                source_index: source.0 as usize,
                candidate_index,
            })
            .collect();
        report.candidate_bond_mapping = bond_order
            .iter()
            .enumerate()
            .map(|(source_index, &candidate_index)| ConversionBondMapping {
                source_id: record.bonds()[source_index].source_id,
                source_index,
                candidate_index,
            })
            .collect();
    }
    let has_loss = report.observations.iter().any(|observation| {
        matches!(
            observation.status,
            PreservationStatus::Dropped
                | PreservationStatus::Unverified
                | PreservationStatus::Inferred
        )
    });
    if (policy == ConversionPolicy::Strict && has_loss) || parsed.is_none() {
        report.status = "refused";
    } else {
        report.status = if has_loss {
            "completed_with_loss"
        } else {
            "completed"
        };
        report.output = Some(candidate);
    }
    Ok(report)
}

fn cjson_conversion_report(target: &str, policy: ConversionPolicy) -> ConversionReport {
    ConversionReport {
        schema: "chematic.conversion-report.v1",
        input_format: "cjson",
        output_format: match target {
            "cjson" => "cjson",
            "mol2" => "mol2",
            _ => "smiles",
        },
        policy,
        status: "completed",
        source_record_index: 0,
        verification_scope: "exact_source_bytes",
        candidate_graph_verified: false,
        output: None,
        observations: Vec::new(),
        candidate_atom_mapping: Vec::new(),
        candidate_bond_mapping: Vec::new(),
    }
}

fn observe(
    report: &mut ConversionReport,
    field: &'static str,
    status: PreservationStatus,
    stage: &'static str,
    reason: &'static str,
) {
    report.observations.push(ConversionObservation {
        field,
        status,
        stage,
        reason,
    });
}

fn retain_cjson_source(
    record: &crate::CjsonRecord,
    report: &mut ConversionReport,
    limits: &ConversionLimits,
) -> Result<(), ConversionError> {
    let output = crate::write_cjson_record(record);
    if output.len() > limits.max_output_bytes {
        return Err(ConversionError::ResourceLimit("output bytes"));
    }
    observe(
        report,
        "source_text",
        PreservationStatus::Preserved,
        "source_retention",
        "Exact UTF-8 input retained after bounded validation of the local CJSON profile; no edits promised.",
    );
    observe(
        report,
        "source_semantics",
        PreservationStatus::OpaqueRetained,
        "source_retention",
        "Unknown fields, decimal spellings and source representations retained verbatim; opaque chemistry and geometry stereo are not interpreted.",
    );
    report.output = Some(output);
    Ok(())
}

fn record_cjson_smiles_losses(record: &crate::CjsonRecord, report: &mut ConversionReport) {
    observe(
        report,
        "source_text",
        PreservationStatus::Dropped,
        "write",
        "CJSON layout, version declaration and array representation are not emitted in SMILES.",
    );
    let conditional = [
        (
            !record.opaque_fields().is_empty(),
            "opaque_fields",
            PreservationStatus::Dropped,
            "projection",
            "Fields listed by the CJSON record inspection API are retained in the source record only and omitted from SMILES; their semantics are not interpreted.",
        ),
        (
            record.name().is_some(),
            "molecule_name",
            PreservationStatus::Dropped,
            "write",
            "Plain SMILES omits source properties.name.",
        ),
        (
            record.coordinates_present(),
            "coordinates",
            PreservationStatus::Dropped,
            "write",
            "Plain SMILES has no source coordinates; geometry stereochemistry is not interpreted.",
        ),
        (
            !record.formal_charges_explicit(),
            "formal_charge_defaults",
            PreservationStatus::Inferred,
            "projection",
            "Absent formalCharges entries default to zero in the local CJSON graph projection.",
        ),
        (
            !record.isotopes_explicit(),
            "isotope_defaults",
            PreservationStatus::Inferred,
            "projection",
            "Absent isotope entries default to unspecified isotope in the local CJSON graph projection.",
        ),
        (
            !record.bond_orders_explicit() && record.molecule().bond_count() > 0,
            "bond_order_defaults",
            PreservationStatus::Inferred,
            "projection",
            "Absent bond orders default to single bonds in the local CJSON graph projection.",
        ),
    ];
    for (present, field, status, stage, reason) in conditional {
        if present {
            observe(report, field, status, stage, reason);
        }
    }
    if !record.partial_charge_datasets().is_empty() {
        observe(
            report,
            "partial_charges",
            PreservationStatus::Dropped,
            "write",
            "Plain SMILES omits source partial-charge datasets; partial charges are not converted into formal charges.",
        );
        observe(
            report,
            "charge_method",
            PreservationStatus::Dropped,
            "write",
            "Plain SMILES omits partial-charge method labels.",
        );
    }
    if record
        .molecule()
        .bonds()
        .any(|(_, bond)| bond.order == chematic_core::BondOrder::Aromatic)
    {
        observe(
            report,
            "aromaticity_interpretation",
            PreservationStatus::Unverified,
            "projection",
            "Orders 1.5 and 4 project to aromatic bonds; atom aromaticity is not perceived. Chemical aromaticity is outside candidate verification.",
        );
    }
    observe(
        report,
        "implicit_hydrogens",
        PreservationStatus::Inferred,
        "write",
        "SMILES hydrogen counts follow graph valence rules; source CJSON hydrogen assignments outside explicit atoms are not interpreted.",
    );
    observe(
        report,
        "stereochemistry",
        PreservationStatus::Unverified,
        "projection",
        "Coordinate-derived and opaque extension stereochemistry are not interpreted or verified.",
    );
}

fn finalize_cjson_smiles(
    record: &crate::CjsonRecord,
    mut report: ConversionReport,
    policy: ConversionPolicy,
    limits: &ConversionLimits,
) -> Result<ConversionReport, ConversionError> {
    if record.molecule().atom_count() > 512 {
        return Err(ConversionError::ResourceLimit(
            "SMILES projection atoms (512)",
        ));
    }
    let (candidate, order) = chematic_smiles::write_with_atom_order(record.molecule());
    if candidate.len() > limits.max_output_bytes {
        return Err(ConversionError::ResourceLimit("output bytes"));
    }
    let parsed = chematic_smiles::parse(&candidate).ok();
    let verified = parsed
        .as_ref()
        .and_then(|molecule| verify_graph(record.molecule(), molecule, &order));
    report.verification_scope =
        "candidate_elements_aromaticity_formal_charges_isotopes_bond_orders";
    report.candidate_graph_verified = verified.is_some();
    observe(
        &mut report,
        "projected_graph",
        if verified.is_some() {
            PreservationStatus::Preserved
        } else {
            PreservationStatus::Unverified
        },
        "candidate_verification",
        if verified.is_some() {
            "Candidate reparsed; atom counts, mapped elements/aromaticity/formal charges/isotopes and bond orders agree. CJSON source IDs in mappings are zero-based array indices, not opaque user metadata. Implicit hydrogens and source chemistry are outside this check."
        } else {
            "Candidate graph could not be verified against the source projection; no correspondence is claimed."
        },
    );
    if let Some(bonds) = verified {
        report.candidate_atom_mapping = order
            .iter()
            .enumerate()
            .map(|(candidate_index, source)| ConversionAtomMapping {
                source_id: source.0,
                source_index: source.0 as usize,
                candidate_index,
            })
            .collect();
        report.candidate_bond_mapping = bonds
            .iter()
            .enumerate()
            .map(|(source_index, &candidate_index)| ConversionBondMapping {
                source_id: source_index as u32,
                source_index,
                candidate_index,
            })
            .collect();
    }
    if policy == ConversionPolicy::Strict || parsed.is_none() {
        report.status = "refused";
    } else {
        report.status = "completed_with_loss";
        report.output = Some(candidate);
    }
    Ok(report)
}

fn convert_cjson(
    text: &str,
    target: &str,
    policy: ConversionPolicy,
    limits: &ConversionLimits,
) -> Result<ConversionReport, ConversionError> {
    let record = crate::parse_cjson_record_with_limits(
        text,
        &crate::CjsonParseLimits {
            max_input_bytes: limits.mol2.max_input_bytes,
            max_atoms: limits.mol2.max_atoms,
            max_bonds: limits.mol2.max_bonds,
            ..Default::default()
        },
    )
    .map_err(ConversionError::CjsonParse)?;
    let mut report = cjson_conversion_report(target, policy);
    if target == "cjson" {
        retain_cjson_source(&record, &mut report, limits)?;
        return Ok(report);
    }
    if target == "mol2" {
        return cjson_mol2::convert(&record, report, limits);
    }
    record_cjson_smiles_losses(&record, &mut report);
    finalize_cjson_smiles(&record, report, policy, limits)
}

/// Check through emitted atom correspondence, including branches, rings and
/// disconnected fragments. Equal formulas/canonical strings are insufficient.
fn verify_graph(source: &Molecule, candidate: &Molecule, order: &[AtomIdx]) -> Option<Vec<usize>> {
    if order.len() != source.atom_count()
        || candidate.atom_count() != source.atom_count()
        || candidate.bond_count() != source.bond_count()
    {
        return None;
    }
    let mut inverse = HashMap::with_capacity(order.len());
    for (output, &input) in order.iter().enumerate() {
        let a = source.atom_opt(input)?;
        let b = candidate.atom(AtomIdx(output as u32));
        if inverse.insert(input, AtomIdx(output as u32)).is_some()
            || a.element != b.element
            || a.aromatic != b.aromatic
            || a.charge != b.charge
            || a.isotope != b.isotope
        {
            return None;
        }
    }
    source
        .bonds()
        .map(|(_, bond)| {
            let (index, written) =
                candidate.bond_between(*inverse.get(&bond.atom1)?, *inverse.get(&bond.atom2)?)?;
            (bond.order == written.order).then_some(index.0 as usize)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::verify_graph;
    use chematic_core::AtomIdx;

    #[test]
    fn candidate_isotope_or_charge_change_prevents_verified_correspondence() {
        let source = chematic_smiles::parse("[13CH3][NH3+]").unwrap();
        let order = [AtomIdx(0), AtomIdx(1)];
        assert!(verify_graph(&source, &source, &order).is_some());
        for changed in ["[12CH3][NH3+]", "[13CH3]N"] {
            let candidate = chematic_smiles::parse(changed).unwrap();
            assert!(verify_graph(&source, &candidate, &order).is_none());
        }
    }
}
