use super::*;
pub(super) fn validate(mol: &Molecule, c: &RdkitLayeredConfig) -> Result<(), RdkitLayeredError> {
    if mol.atom_count() > 512 || mol.bond_count() > 1024 {
        return Err(error("graph_limit", "maximum 512 atoms and 1024 bonds"));
    }
    if c.min_path == 0 || c.min_path > c.max_path || c.max_path > 12 {
        return Err(error(
            "path_range",
            "require 1 <= min_path <= max_path <= 12",
        ));
    }
    if c.fp_size == 0 || c.fp_size > 65536 {
        return Err(error("fp_size", "fp_size must be 1..=65536"));
    }
    if c.max_visits == 0 || c.max_visits > 10_000_000 || c.max_paths == 0 || c.max_paths > 1_000_000
    {
        return Err(error("resource_limit", "invalid visits/path cap"));
    }
    if c.from_atoms.len() > 4096
        || c.from_atoms
            .iter()
            .any(|a| a.0 as usize >= mol.atom_count())
    {
        return Err(error("root_index", "invalid root list or atom index"));
    }
    // RDKit 2026.09.1's unrooted linear branch passes atom paths as bond indices;
    // the public wheel segfaults on CC/CCO. Never reproduce that invalid access.
    if !c.branched_paths && c.from_atoms.is_empty() {
        return Err(error(
            "unrooted_linear_unsupported",
            "branched_paths=false requires explicit from_atoms",
        ));
    }
    if let Some(counts) = &c.atom_counts
        && (counts.len() < mol.atom_count() || counts.len() > 65536)
    {
        return Err(error(
            "atom_counts",
            "counter list must cover atoms and contain <=65536 entries",
        ));
    }
    validate_mask(c)?;
    let mut edges = std::collections::BTreeSet::new();
    for (_, b) in mol.bonds() {
        if b.atom1 == b.atom2 || !edges.insert((b.atom1.min(b.atom2), b.atom1.max(b.atom2))) {
            return Err(error(
                "unsupported_graph",
                "loops and parallel edges are unsupported",
            ));
        }
        if !matches!(
            b.order,
            BondOrder::Single
                | BondOrder::Up
                | BondOrder::Down
                | BondOrder::Double
                | BondOrder::Triple
                | BondOrder::Quadruple
                | BondOrder::Aromatic
        ) {
            return Err(error(
                "unsupported_bond",
                "query, zero and dative bonds are unsupported",
            ));
        }
    }
    Ok(())
}
pub(super) fn same_indices(source: &Molecule, view: &Molecule) -> Result<(), RdkitLayeredError> {
    let same = source.atom_count() == view.atom_count()
        && source.bond_count() == view.bond_count()
        && source.atoms().zip(view.atoms()).all(|((_, a), (_, b))| {
            a.element == b.element && a.isotope == b.isotope && a.wildcard == b.wildcard
        })
        && source
            .bonds()
            .zip(view.bonds())
            .all(|((_, a), (_, b))| a.atom1 == b.atom1 && a.atom2 == b.atom2);
    if !same {
        return Err(error(
            "index_changing_view",
            "sanitized view changes atom/bond correspondence",
        ));
    }
    // A sanitization hook may introduce unsupported bond types.
    if view.bonds().any(|(_, b)| {
        matches!(
            b.order,
            BondOrder::Dative
                | BondOrder::Zero
                | BondOrder::QueryAny
                | BondOrder::QuerySingleOrDouble
                | BondOrder::QuerySingleOrAromatic
                | BondOrder::QueryDoubleOrAromatic
        )
    }) {
        return Err(error(
            "unsupported_bond",
            "sanitized view introduces unsupported bonds",
        ));
    }
    Ok(())
}

fn validate_mask(c: &RdkitLayeredConfig) -> Result<(), RdkitLayeredError> {
    if let Some(mask) = &c.set_only_bits {
        if mask.len() != c.fp_size.div_ceil(8) {
            return Err(error(
                "mask_size",
                "packed mask length does not match fp_size",
            ));
        }
        if !c.fp_size.is_multiple_of(8) && mask.last().is_some_and(|v| *v >> (c.fp_size % 8) != 0) {
            return Err(error("mask_padding", "unused high mask bits must be zero"));
        }
    }
    Ok(())
}
