use crate::{error::MolParseError, mol2000::MolMetadata};
use chematic_core::{Coords3D, Molecule};

/// One group of SDF records that share the same order-sensitive molecular
/// graph and each carry a real 3D conformer.
pub struct ConformerEnsemble {
    /// The molecular graph shared by every member.
    pub mol: Molecule,
    /// Metadata from the first record in the group.
    pub metadata: MolMetadata,
    /// One conformer per member record, in file order.
    pub conformers: Vec<Coords3D>,
}

/// Parse an SDF string and group 3D records that share the same atom and bond
/// ordering into conformer ensembles.
///
/// Records without a 3D conformer are omitted. Grouping is independent of
/// record adjacency, but deliberately remains order-sensitive: repeated
/// conformers emitted by one embedding tool preserve atom and bond order.
pub fn read_sdf_conformer_ensembles(input: &str) -> Result<Vec<ConformerEnsemble>, MolParseError> {
    let mut ensembles: Vec<ConformerEnsemble> = Vec::new();

    for report in super::SdfRecordReader::new(input) {
        let report = report?;
        let Some(conformer) = report.conformer else {
            continue;
        };
        if let Some(existing) = ensembles
            .iter_mut()
            .find(|ensemble| same_graph_identity(&ensemble.mol, &report.mol))
        {
            existing.conformers.push(conformer);
        } else {
            ensembles.push(ConformerEnsemble {
                mol: report.mol,
                metadata: report.meta,
                conformers: vec![conformer],
            });
        }
    }

    Ok(ensembles)
}

fn same_graph_identity(a: &Molecule, b: &Molecule) -> bool {
    a.atom_count() == b.atom_count()
        && a.bond_count() == b.bond_count()
        && a.atoms().zip(b.atoms()).all(|((_, a), (_, b))| {
            a.element == b.element && a.charge == b.charge && a.isotope == b.isotope
        })
        && a.bonds()
            .zip(b.bonds())
            .all(|((_, a), (_, b))| a.atom1 == b.atom1 && a.atom2 == b.atom2 && a.order == b.order)
}
