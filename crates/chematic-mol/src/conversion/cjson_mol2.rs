//! Metadata-aware local CJSON → MOL2 adapter, separate from graph-only writers.
use super::*;
use crate::CjsonRecord;
use chematic_core::BondOrder;
use std::fmt::Write;

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

pub(super) fn convert(
    record: &CjsonRecord,
    mut report: ConversionReport,
    limits: &ConversionLimits,
) -> Result<ConversionReport, ConversionError> {
    report.verification_scope = "candidate_elements_aromaticity_formal_charges_isotopes_bond_orders_coordinates_emitted_metadata";
    note(
        &mut report,
        "source_text",
        PreservationStatus::Dropped,
        "write",
        "CJSON layout, version declaration and array representation are not emitted in MOL2.",
    );
    if !record.opaque_fields().is_empty() {
        note(
            &mut report,
            "opaque_fields",
            PreservationStatus::Dropped,
            "projection",
            "Unknown CJSON fields are omitted from MOL2; opaque references and chemistry are not interpreted.",
        );
    }
    note(
        &mut report,
        "mol2_metadata",
        PreservationStatus::Inferred,
        "write",
        "MOL2 uses generated sequential IDs, element-based atom names/types and SMALL molecule type; charged columns use placeholder substructure 1 LIG. Hybridization and force-field typing are not perceived.",
    );
    note(
        &mut report,
        "stereochemistry",
        PreservationStatus::Unverified,
        "projection",
        "Coordinate-derived and opaque extension stereochemistry are not interpreted or verified.",
    );
    if !record.formal_charges_explicit() {
        note(
            &mut report,
            "formal_charge_defaults",
            PreservationStatus::Inferred,
            "projection",
            "Absent CJSON formalCharges default to zero in the local graph projection.",
        );
    }
    if !record.isotopes_explicit() {
        note(
            &mut report,
            "isotope_defaults",
            PreservationStatus::Inferred,
            "projection",
            "Absent CJSON isotopes default to unspecified isotope in the local graph projection.",
        );
    }
    if !record.bond_orders_explicit() && record.molecule().bond_count() > 0 {
        note(
            &mut report,
            "bond_order_defaults",
            PreservationStatus::Inferred,
            "projection",
            "Absent CJSON bond orders default to single bonds in the local graph projection.",
        );
    }
    if record
        .molecule()
        .bonds()
        .any(|(_, b)| b.order == BondOrder::Aromatic)
    {
        note(
            &mut report,
            "aromaticity_interpretation",
            PreservationStatus::Unverified,
            "projection",
            "Local CJSON aromatic bonds emit ar tokens; atom aromaticity and chemical aromaticity are not perceived.",
        );
    }
    if record.molecule().atoms().any(|(_, a)| a.isotope.is_some()) {
        note(
            &mut report,
            "isotopes",
            PreservationStatus::Dropped,
            "write",
            "This MOL2 adapter has no isotope representation; isotope loss also prevents full projected graph verification.",
        );
    }
    if record
        .molecule()
        .atoms()
        .any(|(_, a)| a.charge != 0 && !(a.element.symbol() == "N" && a.charge == 1))
    {
        note(
            &mut report,
            "formal_charges",
            PreservationStatus::Dropped,
            "write",
            "Only +1 nitrogen is represented by N.4 in this local adapter. Other nonzero formal charges are omitted; partial charges remain separate and are never substituted.",
        );
    }

    // Never guess a charge dataset, invent geometry, trim names, or allow a
    // source string to inject a new MOL2 line/section. Same-format retention
    // remains available for every valid CJSON record rejected here.
    let sets = record.partial_charge_datasets();
    let partial = if sets.len() == 1 {
        sets.iter().next()
    } else {
        None
    };
    let invalid_method = partial.is_some_and(|(method, _)| {
        method.is_empty()
            || method.eq_ignore_ascii_case("NO_CHARGES")
            || method.chars().any(|c| c.is_whitespace() || c.is_control())
            || method.to_ascii_uppercase().contains("@<TRIPOS>")
    });
    let invalid_name = record.name().is_some_and(|name| {
        name.is_empty()
            || name != name.trim()
            || name.chars().any(char::is_control)
            || name.to_ascii_uppercase().contains("@<TRIPOS>")
    });
    let mut refuse = false;
    if sets.len() > 1 || invalid_method {
        note(
            &mut report,
            "partial_charges",
            PreservationStatus::Unverified,
            "projection",
            "MOL2 requires one representable non-NO_CHARGES method token. Multiple datasets are not selected arbitrarily and invalid labels are not rewritten; no candidate emitted.",
        );
        refuse = true;
    }
    if invalid_name {
        note(
            &mut report,
            "molecule_name",
            PreservationStatus::Unverified,
            "write",
            "Source name is empty, has boundary whitespace/control characters or a section marker; this adapter refuses rather than trim or inject MOL2 syntax.",
        );
        refuse = true;
    }
    if !record.coordinates_present() || record.molecule().atom_count() == 0 {
        note(
            &mut report,
            "coordinates",
            PreservationStatus::Unverified,
            "write",
            "MOL2 requires nonempty atoms and explicit 3D coordinates in this adapter; missing geometry is not filled or generated.",
        );
        refuse = true;
    }
    if refuse {
        report.status = "refused";
        return Ok(report);
    }
    if record.name().is_none() {
        note(
            &mut report,
            "molecule_name",
            PreservationStatus::Inferred,
            "write",
            "Absent properties.name uses the placeholder molecule name chematic.",
        );
    }
    let name = record.name().unwrap_or("chematic");
    let method = partial.map_or("NO_CHARGES", |(method, _)| method.as_str());
    let mut candidate = format!(
        "@<TRIPOS>MOLECULE\n{name}\n{} {} 0 0 0\nSMALL\n{method}\n\n@<TRIPOS>ATOM\n",
        record.molecule().atom_count(),
        record.molecule().bond_count()
    );
    for (idx, atom) in record.molecule().atoms() {
        let i = idx.0 as usize;
        let sym = atom.element.symbol();
        let atom_type = if sym == "N" && atom.charge == 1 {
            "N.4"
        } else {
            sym
        };
        let (x, y, z) = record.coordinates()[i];
        // Shortest round-trip decimals retain every finite f64, including -0
        // and subnormals. No fixed-width coordinate/charge rounding.
        write!(
            candidate,
            "{} {sym}{} {x} {y} {z} {atom_type}",
            i + 1,
            i + 1
        )
        .unwrap();
        if let Some((_, values)) = partial {
            write!(candidate, " 1 LIG {}", values[i]).unwrap();
        }
        candidate.push('\n');
    }
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
    let after = parsed.as_ref();
    checked(
        &mut report,
        "coordinates",
        verified.is_some()
            && after.is_some_and(|after| {
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
            &mut report,
            "molecule_name",
            after.is_some_and(|after| after.name() == Some(name)),
            "Reparsed MOL2 molecule name equals source properties.name verbatim; no trimming or placeholder substitution.",
        );
    }
    if let Some((method, values)) = partial {
        checked(
            &mut report,
            "partial_charges",
            verified.is_some()
                && after.is_some_and(|after| {
                    values.len() == after.atoms().len()
                        && values.iter().zip(after.atoms()).all(|(a, b)| {
                            b.partial_charge.is_some_and(|b| a.to_bits() == b.to_bits())
                        })
                }),
            "Complete source partial charges must reparse bit for bit through verified atom correspondence. No calculation, rounding into formal charges or filling missing values.",
        );
        checked(
            &mut report,
            "charge_method",
            after.is_some_and(|after| after.charge_type() == Some(method.as_str())),
            "Source dataset method label equals the reparsed MOL2 charge-method token; the declared method is not recomputed or independently validated.",
        );
    }
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
    Ok(report)
}
