//! All-reference heavy-atom RMSD with explicit, conservative correspondence.
use crate::{
    Mol2Error, Mol2SourceReader, Mol2SourceRecord, Mol2StreamError, Mol2StreamLimits,
    parse_mol2_source_record_with_limits,
};
use chematic_core::Element;

#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct Mol2ReferenceRmsd {
    pub reference_index: usize,
    pub reference_name: String,
    pub probe_name: String,
    pub heavy_atom_count: usize,
    pub rmsd: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Mol2RmsdError {
    Probe(Mol2Error),
    Reference(Mol2StreamError),
    EmptyReference,
    NoHeavyAtoms,
    IdentityMismatch { reference_index: usize },
    NonFiniteComputation { reference_index: usize },
}

impl core::fmt::Display for Mol2RmsdError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Probe(error) => write!(f, "MOL2 probe: {error}"),
            Self::Reference(error) => write!(f, "MOL2 reference: {error}"),
            Self::EmptyReference => write!(f, "MOL2 reference contains no records"),
            Self::NoHeavyAtoms => write!(f, "MOL2 RMSD requires at least one heavy atom"),
            Self::IdentityMismatch { reference_index } => write!(
                f,
                "MOL2 reference {reference_index}: source atom identity/order or indexed graph differs from probe"
            ),
            Self::NonFiniteComputation { reference_index } => write!(
                f,
                "MOL2 reference {reference_index}: RMSD computation is not finite"
            ),
        }
    }
}

impl std::error::Error for Mol2RmsdError {}

fn same_identity(a: &Mol2SourceRecord, b: &Mol2SourceRecord) -> bool {
    if a.atoms().len() != b.atoms().len() || a.bonds().len() != b.bonds().len() {
        return false;
    }
    let atoms_match =
        a.atoms().iter().zip(b.atoms()).all(|(a, b)| {
            a.source_id == b.source_id && a.name == b.name && a.atom_type == b.atom_type
        });
    let graph_atoms_match =
        a.molecule()
            .atoms()
            .zip(b.molecule().atoms())
            .all(|((_, a), (_, b))| {
                a.element == b.element && a.charge == b.charge && a.isotope == b.isotope
            });
    // Bond row order and undirected endpoint spelling do not define atom identity.
    let bonds = |record: &Mol2SourceRecord| {
        let mut rows: Vec<_> = record
            .bonds()
            .iter()
            .map(|bond| {
                (
                    bond.atom_ids[0].min(bond.atom_ids[1]),
                    bond.atom_ids[0].max(bond.atom_ids[1]),
                    bond.bond_type.clone(),
                )
            })
            .collect();
        rows.sort_unstable();
        rows
    };
    atoms_match && graph_atoms_match && bonds(a) == bonds(b)
}

fn heavy_coords(record: &Mol2SourceRecord) -> Vec<[f64; 3]> {
    record
        .molecule()
        .atoms()
        .zip(record.coordinates())
        .filter_map(|((_, atom), &(x, y, z))| (atom.element != Element::H).then_some([x, y, z]))
        .collect()
}

/// Compare every reference record to exactly one MOL2 probe, in file order.
/// Coordinates are assumed to share their MOL2 unit (normally Å). Hydrogens
/// are excluded from the metric, but all atoms/bonds must retain the same
/// source IDs, atom names/types, row order, elements and formal-charge graph.
/// Partial charge values and opaque source metadata do not define identity.
/// No graph-isomorphism search, automatic atom remapping or reflection occurs.
/// `align = false` computes raw RMSD; `true` uses a proper Kabsch rotation and
/// translation over heavy atoms. Any invalid/mismatched record aborts the
/// whole result; records are never silently skipped. Limits cover the reference
/// stream and each record, including the probe. Values must compute finitely.
pub fn mol2_reference_rmsds(
    reference_text: &str,
    probe_text: &str,
    align: bool,
    limits: &Mol2StreamLimits,
) -> Result<Vec<Mol2ReferenceRmsd>, Mol2RmsdError> {
    let probe = parse_mol2_source_record_with_limits(probe_text, &limits.record)
        .map_err(Mol2RmsdError::Probe)?;
    let probe_coords = heavy_coords(&probe);
    if probe_coords.is_empty() {
        return Err(Mol2RmsdError::NoHeavyAtoms);
    }
    let mut results = Vec::new();
    for (reference_index, record) in
        Mol2SourceReader::with_limits(reference_text, *limits).enumerate()
    {
        let reference = record.map_err(Mol2RmsdError::Reference)?;
        if !same_identity(&reference, &probe) {
            return Err(Mol2RmsdError::IdentityMismatch { reference_index });
        }
        let coords = heavy_coords(&reference);
        let rmsd = if align {
            chematic_3d::align_coords(&coords, &probe_coords).rmsd
        } else {
            chematic_3d::rmsd_no_align(&coords, &probe_coords)
        };
        if !rmsd.is_finite() {
            return Err(Mol2RmsdError::NonFiniteComputation { reference_index });
        }
        results.push(Mol2ReferenceRmsd {
            reference_index,
            reference_name: reference.name().unwrap_or_default().into(),
            probe_name: probe.name().unwrap_or_default().into(),
            heavy_atom_count: probe_coords.len(),
            rmsd,
        });
    }
    if results.is_empty() {
        return Err(Mol2RmsdError::EmptyReference);
    }
    Ok(results)
}
