//! Metadata-aware local CJSON → MOL2 adapter, separate from graph-only writers.
use super::*;
use crate::{CjsonRecord, Mol2SourceRecord};
use chematic_core::BondOrder;
use std::fmt::Write;

type PartialDataset<'a> = Option<(&'a String, &'a Vec<f64>)>;

fn note(
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

fn checked(report: &mut ConversionReport, field: &'static str, equal: bool, reason: &'static str) {
    note(
        report,
        field,
        if equal {
            PreservationStatus::Preserved
        } else {
            PreservationStatus::Unverified
        },
        "candidate_verification",
        reason,
    );
}

fn write_formal_charges(record: &CjsonRecord, candidate: &mut String) {
    // UNITY carries formal charge independently of the optional ATOM partial
    // charge. Open Babel 3.2.1 skips trailing UNITY, so call this before BOND.
    if record.molecule().atoms().any(|(_, atom)| atom.charge != 0) {
        candidate.push_str("@<TRIPOS>UNITY_ATOM_ATTR\n");
        for (idx, atom) in record
            .molecule()
            .atoms()
            .filter(|(_, atom)| atom.charge != 0)
        {
            writeln!(candidate, "{} 1\ncharge {}", idx.0 + 1, atom.charge).unwrap();
        }
    }
}

fn record_projection_contract(record: &CjsonRecord, report: &mut ConversionReport) {
    report.verification_scope = "candidate_elements_aromaticity_formal_charges_isotopes_bond_orders_coordinates_emitted_metadata";
    note(
        report,
        "source_text",
        PreservationStatus::Dropped,
        "write",
        "CJSON layout, version declaration and array representation are not emitted in MOL2.",
    );
    if !record.opaque_fields().is_empty() {
        note(
            report,
            "opaque_fields",
            PreservationStatus::Dropped,
            "projection",
            "Unknown CJSON fields are omitted from MOL2; opaque references and chemistry are not interpreted.",
        );
    }
    note(
        report,
        "mol2_metadata",
        PreservationStatus::Inferred,
        "write",
        "MOL2 uses generated sequential IDs, element-based atom names/types and SMALL molecule type; charged columns use placeholder substructure 1 LIG. Hybridization and force-field typing are not perceived.",
    );
    note(
        report,
        "stereochemistry",
        PreservationStatus::Unverified,
        "projection",
        "Coordinate-derived and opaque extension stereochemistry are not interpreted or verified.",
    );
    record_defaulted_fields(record, report);
    record_unrepresentable_fields(record, report);
}

fn record_defaulted_fields(record: &CjsonRecord, report: &mut ConversionReport) {
    let defaults = [
        (
            !record.formal_charges_explicit(),
            "formal_charge_defaults",
            "Absent CJSON formalCharges default to zero in the local graph projection.",
        ),
        (
            !record.isotopes_explicit(),
            "isotope_defaults",
            "Absent CJSON isotopes default to unspecified isotope in the local CJSON graph projection.",
        ),
        (
            !record.bond_orders_explicit() && record.molecule().bond_count() > 0,
            "bond_order_defaults",
            "Absent CJSON bond orders default to single bonds in the local CJSON graph projection.",
        ),
    ];
    for (present, field, reason) in defaults {
        if present {
            note(
                report,
                field,
                PreservationStatus::Inferred,
                "projection",
                reason,
            );
        }
    }
    if record
        .molecule()
        .bonds()
        .any(|(_, bond)| bond.order == BondOrder::Aromatic)
    {
        note(
            report,
            "aromaticity_interpretation",
            PreservationStatus::Unverified,
            "projection",
            "Local CJSON aromatic bonds emit ar tokens; atom aromaticity and chemical aromaticity are not perceived.",
        );
    }
}

fn record_unrepresentable_fields(record: &CjsonRecord, report: &mut ConversionReport) {
    if record
        .molecule()
        .atoms()
        .any(|(_, atom)| atom.isotope.is_some())
    {
        note(
            report,
            "isotopes",
            PreservationStatus::Dropped,
            "write",
            "This MOL2 adapter has no isotope representation; isotope loss also prevents full projected graph verification.",
        );
    }
}

fn projection_inputs<'a>(
    record: &'a CjsonRecord,
    report: &mut ConversionReport,
) -> Option<(&'a str, PartialDataset<'a>)> {
    // Never guess a charge dataset, invent geometry, trim names, or allow a
    // source string to inject a new MOL2 line/section. Same-format retention
    // remains available for every valid CJSON record rejected here.
    let sets = record.partial_charge_datasets();
    let partial = (sets.len() == 1).then(|| sets.iter().next()).flatten();
    let invalid_method = partial.is_some_and(|(method, _)| {
        method.is_empty()
            || method.eq_ignore_ascii_case("NO_CHARGES")
            || method.chars().any(|c| c.is_whitespace() || c.is_control())
            || method.to_ascii_uppercase().contains("@<TRIPOS>")
    });
    let invalid_name = record.has_conflicting_names()
        || record.name().is_some_and(|name| {
            name.is_empty()
                || name != name.trim()
                || name.chars().any(char::is_control)
                || name.to_ascii_uppercase().contains("@<TRIPOS>")
        });
    if sets.len() > 1 || invalid_method {
        note(
            report,
            "partial_charges",
            PreservationStatus::Unverified,
            "projection",
            "MOL2 requires one representable non-NO_CHARGES method token. Multiple datasets are not selected arbitrarily and invalid labels are not rewritten; no candidate emitted.",
        );
    }
    if invalid_name {
        note(
            report,
            "molecule_name",
            PreservationStatus::Unverified,
            "write",
            "Root name and legacy properties.name conflict, or the source name is empty, has boundary whitespace/control characters or a section marker. No title is selected or rewritten; this adapter refuses. Exact source retention remains available.",
        );
    }
    let invalid_geometry = !record.coordinates_present() || record.molecule().atom_count() == 0;
    if invalid_geometry {
        note(
            report,
            "coordinates",
            PreservationStatus::Unverified,
            "write",
            "MOL2 requires nonempty atoms and explicit 3D coordinates in this adapter; missing geometry is not filled or generated.",
        );
    }
    if sets.len() > 1 || invalid_method || invalid_name || invalid_geometry {
        report.status = "refused";
        return None;
    }
    if record.name().is_none() {
        note(
            report,
            "molecule_name",
            PreservationStatus::Inferred,
            "write",
            "Absent root name and legacy properties.name use the placeholder molecule name chematic.",
        );
    }
    Some((record.name().unwrap_or("chematic"), partial))
}

fn write_candidate(
    record: &CjsonRecord,
    name: &str,
    partial: PartialDataset<'_>,
    limits: &ConversionLimits,
) -> Result<String, ConversionError> {
    let method = partial.map_or("NO_CHARGES", |(method, _)| method.as_str());
    let mut candidate = format!(
        "@<TRIPOS>MOLECULE\n{name}\n{} {} 0 0 0\nSMALL\n{method}\n\n@<TRIPOS>ATOM\n",
        record.molecule().atom_count(),
        record.molecule().bond_count()
    );
    for (idx, atom) in record.molecule().atoms() {
        let i = idx.0 as usize;
        let symbol = atom.element.symbol();
        let atom_type = if symbol == "N" && atom.charge == 1 {
            "N.4"
        } else {
            symbol
        };
        let (x, y, z) = record.coordinates()[i];
        write!(
            candidate,
            "{} {symbol}{} {x} {y} {z} {atom_type}",
            i + 1,
            i + 1
        )
        .unwrap();
        if let Some((_, values)) = partial {
            write!(candidate, " 1 LIG {}", values[i]).unwrap();
        }
        candidate.push('\n');
    }
    write_formal_charges(record, &mut candidate);
    candidate.push_str("@<TRIPOS>BOND\n");
    for (idx, bond) in record.molecule().bonds() {
        let order = match bond.order {
            BondOrder::Single => "1",
            BondOrder::Double => "2",
            BondOrder::Triple => "3",
            BondOrder::Aromatic => "ar",
            _ => unreachable!("strict CJSON projection permits only these four orders"),
        };
        writeln!(
            candidate,
            "{} {} {} {order}",
            idx.0 + 1,
            bond.atom1.0 + 1,
            bond.atom2.0 + 1
        )
        .unwrap();
    }
    if candidate.len() > limits.max_output_bytes {
        return Err(ConversionError::ResourceLimit("output bytes"));
    }
    Ok(candidate)
}

fn verify_metadata(
    record: &CjsonRecord,
    partial: PartialDataset<'_>,
    parsed: Option<&Mol2SourceRecord>,
    graph_verified: bool,
    report: &mut ConversionReport,
) {
    checked(
        report,
        "formal_charges",
        graph_verified,
        "Explicit UNITY charge attributes (zero when absent) must reparse exactly through verified atom correspondence. Partial-charge values and method labels do not determine formal charge.",
    );
    checked(
        report,
        "coordinates",
        graph_verified
            && parsed.is_some_and(|after| {
                record.coordinates().len() == after.coordinates().len()
                    && record
                        .coordinates()
                        .iter()
                        .zip(after.coordinates())
                        .all(|(a, b)| {
                            a.0.to_bits() == b.0.to_bits()
                                && a.1.to_bits() == b.1.to_bits()
                                && a.2.to_bits() == b.2.to_bits()
                        })
            }),
        "Reparsed coordinates must agree bit for bit through verified atom correspondence, including signed zero and subnormals; original decimal spelling is outside this check.",
    );
    if let Some(name) = record.name() {
        checked(
            report,
            "molecule_name",
            parsed.is_some_and(|after| after.name() == Some(name)),
            "Reparsed MOL2 molecule name equals the unambiguous source root name or legacy properties.name verbatim; no trimming or placeholder substitution.",
        );
    }
    if let Some((method, values)) = partial {
        checked(
            report,
            "partial_charges",
            graph_verified
                && parsed.is_some_and(|after| {
                    values.len() == after.atoms().len()
                        && values.iter().zip(after.atoms()).all(|(a, b)| {
                            b.partial_charge.is_some_and(|b| a.to_bits() == b.to_bits())
                        })
                }),
            "Complete source partial charges must reparse bit for bit through verified atom correspondence. No calculation, rounding into formal charges or filling missing values.",
        );
        checked(
            report,
            "charge_method",
            parsed.is_some_and(|after| after.charge_type() == Some(method.as_str())),
            "Source dataset method label equals the reparsed MOL2 charge-method token; the declared method is not recomputed or independently validated.",
        );
    }
}

fn verify_candidate(
    record: &CjsonRecord,
    partial: PartialDataset<'_>,
    candidate: String,
    mut report: ConversionReport,
    limits: &ConversionLimits,
) -> ConversionReport {
    let parsed = parse_mol2_source_record_with_limits(
        &candidate,
        &Mol2ParseLimits {
            max_input_bytes: limits.max_output_bytes,
            ..limits.mol2
        },
    )
    .ok();
    let order: Vec<_> = (0..record.molecule().atom_count())
        .map(|i| AtomIdx(i as u32))
        .collect();
    let verified = parsed
        .as_ref()
        .and_then(|after| verify_graph(record.molecule(), after.molecule(), &order));
    report.candidate_graph_verified = verified.is_some();
    checked(
        &mut report,
        "projected_graph",
        verified.is_some(),
        "Candidate MOL2 reparsed; verification requires mapped counts, elements, aromatic flags, formal charges, isotopes and bond orders. Failed graph verification exposes no correspondence; source IDs are CJSON zero-based array indices.",
    );
    verify_metadata(
        record,
        partial,
        parsed.as_ref(),
        verified.is_some(),
        &mut report,
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
    if report.policy == ConversionPolicy::Strict || parsed.is_none() {
        report.status = "refused";
    } else {
        report.status = "completed_with_loss";
        report.output = Some(candidate);
    }
    report
}

pub(super) fn convert(
    record: &CjsonRecord,
    mut report: ConversionReport,
    limits: &ConversionLimits,
) -> Result<ConversionReport, ConversionError> {
    record_projection_contract(record, &mut report);
    let Some((name, partial)) = projection_inputs(record, &mut report) else {
        return Ok(report);
    };
    let candidate = write_candidate(record, name, partial, limits)?;
    Ok(verify_candidate(record, partial, candidate, report, limits))
}
