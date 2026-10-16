//! Endpoint-aware inspection of declared 1-4 alkene stereo (RDKit #9719).
use super::RdkitEmbedError;
use chematic_core::Molecule;

/// Relative side of the two outer atoms of a declared alkene path.
/// This is a path relation, not a CIP E/Z label for the whole double bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlkenePathStereo {
    Cis,
    Trans,
}

/// Classify a real four-atom bonded path around its central double bond.
/// The complete path is normalized to the prepared bond's begin/end frame;
/// reversing all four indices therefore preserves the relation. Atom indices
/// use the input/prepared view order. To refer to an implicit H, first call
/// `add_hydrogens` and use its explicit atom index. No geometry or CIP labeling
/// is inferred. Undeclared/unspecified central stereo returns `None`; repeated,
/// out-of-range, disconnected or non-double-centered paths return an error.
/// This additive query changes no bounds matrices or embedding defaults.
pub fn rdkit_alkene_path_stereo(
    molecule: &Molecule,
    mut path: [usize; 4],
) -> Result<Option<AlkenePathStereo>, RdkitEmbedError> {
    let invalid = |reason: &str| RdkitEmbedError::Unsupported(format!("alkene path: {reason}"));
    if path.iter().any(|&a| a >= molecule.atom_count()) {
        return Err(invalid("atom index out of range"));
    }
    if (0..4).any(|i| (i + 1..4).any(|j| path[i] == path[j])) {
        return Err(invalid("four distinct atom indices are required"));
    }
    let view =
        chematic_smiles::rdkit_mol_view(molecule).map_err(|error| invalid(&error.to_string()))?;
    if path.iter().any(|&a| a >= view.num_atoms()) {
        return Err(invalid("prepared atom index out of range"));
    }
    if path
        .windows(2)
        .any(|pair| view.bond_between(pair[0], pair[1]).is_none())
    {
        return Err(invalid("consecutive atoms must be bonded"));
    }
    let b = view
        .bond_between(path[1], path[2])
        .expect("validated central edge");
    let bond = &view.bonds[b];
    if bond.bond_type != 2 {
        return Err(invalid("central bond must be a double bond"));
    }
    if bond.stereo <= 1 {
        return Ok(None);
    }
    if bond.stereo_atoms.len() != 2 {
        return Err(invalid("declared stereo lacks exactly two reference atoms"));
    }
    if path[1] == bond.end {
        path.reverse();
    }
    match super::bounds::atom_stereo(&view, b, path[0], path[3]) {
        2 => Ok(Some(AlkenePathStereo::Cis)),
        3 => Ok(Some(AlkenePathStereo::Trans)),
        _ => Err(invalid("unsupported declared stereo")),
    }
}
