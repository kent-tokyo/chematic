use chematic_core::Element;

use super::{Occupancy, PoscarDocument, PoscarError};

pub(super) struct PoscarWritePlan {
    pub(super) group_order: Vec<Element>,
    pub(super) permutation: Vec<usize>,
    pub(super) counts: Vec<usize>,
}

pub(super) fn poscar_write_plan(doc: &PoscarDocument) -> Result<PoscarWritePlan, PoscarError> {
    let elements = writable_elements(doc)?;
    validate_sidecars(doc, elements.len())?;
    let group_order = first_appearance_order(&elements);
    let group_index = |element: Element| {
        group_order
            .iter()
            .position(|&group| group == element)
            .expect("group order is built from the same elements")
    };
    let mut permutation: Vec<usize> = (0..elements.len()).collect();
    permutation.sort_by_key(|&index| group_index(elements[index]));
    reject_opaque_reordering(doc, &permutation)?;
    let mut counts = vec![0; group_order.len()];
    for element in elements {
        counts[group_index(element)] += 1;
    }
    Ok(PoscarWritePlan {
        group_order,
        permutation,
        counts,
    })
}

fn writable_elements(doc: &PoscarDocument) -> Result<Vec<Element>, PoscarError> {
    doc.structure
        .sites()
        .iter()
        .enumerate()
        .map(|(index, site)| {
            if site.species.len() != 1 {
                return Err(PoscarError::Unwritable {
                    detail: format!(
                        "site {index} has {} species; POSCAR has no disorder/partial-occupancy representation (every site must be exactly one species)",
                        site.species.len()
                    ),
                });
            }
            let species = &site.species[0];
            if (species.occupancy.value() - 1.0).abs() > Occupancy::SUM_TOLERANCE {
                return Err(PoscarError::Unwritable {
                    detail: format!(
                        "site {index} has occupancy {} (POSCAR requires full occupancy 1.0)",
                        species.occupancy.value()
                    ),
                });
            }
            Ok(species.element)
        })
        .collect()
}

fn first_appearance_order(elements: &[Element]) -> Vec<Element> {
    let mut order = Vec::new();
    for &element in elements {
        if !order.contains(&element) {
            order.push(element);
        }
    }
    order
}

fn reject_opaque_reordering(
    doc: &PoscarDocument,
    permutation: &[usize],
) -> Result<(), PoscarError> {
    let reordered = permutation
        .iter()
        .enumerate()
        .any(|(index, &old)| index != old);
    if reordered && doc.predictor_corrector.is_some() {
        return Err(PoscarError::Unwritable {
            detail: "sites need reordering into contiguous species groups, but predictor_corrector holds opaque verbatim MD-restart data that cannot be safely permuted -- reorder sites (or clear predictor_corrector) before writing".to_string(),
        });
    }
    Ok(())
}

fn validate_sidecars(doc: &PoscarDocument, site_count: usize) -> Result<(), PoscarError> {
    validate_sidecar_len("velocities", doc.velocities.as_deref(), site_count)?;
    validate_sidecar_len(
        "selective_dynamics",
        doc.selective_dynamics.as_deref(),
        site_count,
    )?;
    if doc.predictor_corrector.is_some() && doc.velocities.is_none() {
        return Err(PoscarError::Unwritable {
            detail: "predictor_corrector is present but velocities is None -- VASP's CONTCAR layout always has ion velocities before predictor-corrector data".to_string(),
        });
    }
    Ok(())
}

fn validate_sidecar_len<T>(
    name: &str,
    values: Option<&[T]>,
    site_count: usize,
) -> Result<(), PoscarError> {
    if let Some(values) = values
        && values.len() != site_count
    {
        return Err(PoscarError::Unwritable {
            detail: format!(
                "{name} has {} entries but structure has {site_count} site(s)",
                values.len()
            ),
        });
    }
    Ok(())
}
