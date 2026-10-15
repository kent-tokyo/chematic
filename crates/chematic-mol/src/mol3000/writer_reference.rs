use super::*;

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Serialise `mol` to a MOL V3000 (Extended Ctab) string.
///
/// `coords[i]` is the `(x, y)` position for atom `i`.  Atoms beyond
/// `coords.len()` receive `(0.0, 0.0, 0.0)`.
///
/// **Does not preserve `Chirality::SquarePlanar` stereo.** This is a
/// 2D-only writer with no z channel; MOL/CTfile has no other field for a
/// non-tetrahedral stereo tag either (see
/// `docs/rfcs/square_planar_mol_io_rfc.md`). A square-planar-tagged atom is
/// written with no indication anything was dropped. If `mol` may carry
/// `Chirality::SquarePlanar`, use
/// [`write_mol_v3000_with_conformer_checked`] instead, which fails closed
/// with a typed error rather than silently discarding the tag.
fn v3000_atom_symbol(mol: &Molecule, idx: AtomIdx, atom: &Atom) -> &'static str {
    if !atom.wildcard {
        return atom.element.symbol();
    }
    match mol.r_group_label(idx).and_then(|label| label.number()) {
        Some(_) => "R#",
        None if mol.r_group_label(idx).is_some() => "R",
        None => "*",
    }
}

fn append_v3000_rgroup_property(line: &mut String, mol: &Molecule, idx: AtomIdx) {
    if let Some(number) = mol.r_group_label(idx).and_then(|label| label.number()) {
        line.push_str(&format!(" RGROUPS=(1 {number})"));
    }
}

fn append_header(out: &mut String, mol: &Molecule, metadata: &MolMetadata, dimensional_code: &str) {
    out.push_str(&metadata.name);
    out.push('\n');
    out.push_str(dimensional_code);
    out.push_str(&metadata.comment);
    out.push('\n');
    out.push_str("  0  0  0  0  0  0  0  0  0  0999 V3000\n");
    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!(
        "M  V30 COUNTS {} {} {} 0 0\n",
        mol.atom_count(),
        mol.bond_count(),
        metadata.v3000_sgroups.len()
    ));
}

fn append_atom_properties(
    line: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
    idx: AtomIdx,
    atom: &Atom,
) {
    if atom.charge != 0 {
        line.push_str(&format!(" CHG={}", atom.charge));
    }
    if let Some(isotope) = atom.isotope {
        line.push_str(&format!(" MASS={isotope}"));
    }
    if let Some(hydrogens) = atom.hydrogen_count {
        line.push_str(&format!(" HCOUNT={hydrogens}"));
    }
    append_v3000_rgroup_property(line, mol, idx);
    if let Some((_, properties)) = metadata
        .v3000_atom_properties
        .iter()
        .find(|(atom_id, _)| *atom_id == idx.0 + 1)
    {
        line.push(' ');
        line.push_str(properties);
    }
}

fn append_atoms_2d(
    out: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) {
    out.push_str("M  V30 BEGIN ATOM\n");
    for (idx, atom) in mol.atoms() {
        let (x, y) = coords.get(idx.0 as usize).copied().unwrap_or((0.0, 0.0));
        let mut line = format!(
            "M  V30 {} {} ",
            idx.0 + 1,
            v3000_atom_symbol(mol, idx, atom)
        );
        crate::mol2000::push_fixed4(&mut line, x, 0);
        line.push(' ');
        crate::mol2000::push_fixed4(&mut line, y, 0);
        line.push_str(&format!(" 0.0000 {}", atom.atom_map.unwrap_or(0)));
        append_atom_properties(&mut line, mol, metadata, idx, atom);
        out.push_str(&line);
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
}

fn append_atoms_3d(out: &mut String, mol: &Molecule, metadata: &MolMetadata, conformer: &Coords3D) {
    out.push_str("M  V30 BEGIN ATOM\n");
    for (idx, atom) in mol.atoms() {
        let point = conformer
            .points
            .get(idx.0 as usize)
            .copied()
            .unwrap_or(Point3::zero());
        let mut line = format!(
            "M  V30 {} {} {:.4} {:.4} {:.4} {}",
            idx.0 + 1,
            v3000_atom_symbol(mol, idx, atom),
            point.x,
            point.y,
            point.z,
            atom.atom_map.unwrap_or(0)
        );
        append_atom_properties(&mut line, mol, metadata, idx, atom);
        out.push_str(&line);
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
}

fn v3000_bond_order(order: BondOrder) -> u8 {
    match order {
        BondOrder::Zero => 0,
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic | BondOrder::Quadruple => 4,
        BondOrder::QuerySingleOrDouble => 5,
        BondOrder::QuerySingleOrAromatic => 6,
        BondOrder::QueryDoubleOrAromatic => 7,
        BondOrder::QueryAny => 8,
        BondOrder::Dative => 9,
    }
}

fn bond_properties(metadata: &MolMetadata, bond_id: u32, cfg_written: bool) -> String {
    metadata
        .v3000_bond_properties
        .iter()
        .find(|(id, _)| *id == bond_id)
        .map(|(_, properties)| {
            if !cfg_written {
                return format!(" {properties}");
            }
            // A depiction may regenerate the explicit unknown marker.
            // Keep it exactly once while retaining all other properties.
            let properties = properties
                .split_whitespace()
                .filter(|token| !token.starts_with("CFG="))
                .collect::<Vec<_>>()
                .join(" ");
            if properties.is_empty() {
                String::new()
            } else {
                format!(" {properties}")
            }
        })
        .unwrap_or_default()
}

fn append_bonds_2d(
    out: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
    depiction: Option<&crate::stereo_depiction::StereoDepiction>,
) {
    out.push_str("M  V30 BEGIN BOND\n");
    let kekule = crate::mol2000::kekule_orders_for_writing(mol);
    for (bond_idx, bond) in mol.bonds() {
        let written = kekule
            .as_ref()
            .and_then(|orders| orders.get(&bond_idx).copied())
            .unwrap_or(bond.order);
        let atom1 = bond.atom1.0 + 1;
        let atom2 = bond.atom2.0 + 1;
        let (atom1, atom2, stereo) = match depiction {
            Some(depiction) => match depiction.wedges.get(&bond_idx) {
                Some(wedge) => {
                    let (start, end) = if wedge.start == bond.atom1 {
                        (atom1, atom2)
                    } else {
                        (atom2, atom1)
                    };
                    let cfg = if wedge.order == BondOrder::Up { 1 } else { 3 };
                    (start, end, format!(" CFG={cfg}"))
                }
                None if depiction.unexpressed_double_bonds.contains(&bond_idx)
                    || depiction.unspecified_double_bonds.contains(&bond_idx) =>
                {
                    (atom1, atom2, " CFG=2".to_string())
                }
                None => (atom1, atom2, String::new()),
            },
            None => {
                let cfg = match bond.order {
                    BondOrder::Up => " CFG=1",
                    BondOrder::Down => " CFG=3",
                    _ => "",
                };
                (atom1, atom2, cfg.to_string())
            }
        };
        let bond_id = bond_idx.0 + 1;
        out.push_str(&format!(
            "M  V30 {bond_id} {} {atom1} {atom2}{stereo}{}\n",
            v3000_bond_order(written),
            bond_properties(metadata, bond_id, !stereo.is_empty())
        ));
    }
    out.push_str("M  V30 END BOND\n");
}

fn append_bonds_3d(out: &mut String, mol: &Molecule, metadata: &MolMetadata) {
    out.push_str("M  V30 BEGIN BOND\n");
    let kekule = crate::mol2000::kekule_orders_for_writing(mol);
    for (bond_idx, bond) in mol.bonds() {
        let written = kekule
            .as_ref()
            .and_then(|orders| orders.get(&bond_idx).copied())
            .unwrap_or(bond.order);
        let bond_id = bond_idx.0 + 1;
        out.push_str(&format!(
            "M  V30 {bond_id} {} {} {}{}\n",
            v3000_bond_order(written),
            bond.atom1.0 + 1,
            bond.atom2.0 + 1,
            bond_properties(metadata, bond_id, false)
        ));
    }
    out.push_str("M  V30 END BOND\n");
}

fn append_sgroups(out: &mut String, metadata: &MolMetadata) {
    if metadata.v3000_sgroups.is_empty() {
        return;
    }
    out.push_str("M  V30 BEGIN SGROUP\n");
    for line in &metadata.v3000_sgroups {
        out.push_str("M  V30 ");
        out.push_str(line);
        out.push('\n');
    }
    out.push_str("M  V30 END SGROUP\n");
}

fn append_collections(out: &mut String, mol: &Molecule) {
    let groups = mol
        .stereo_groups()
        .iter()
        .enumerate()
        .filter(|(group_index, group)| {
            !group.atom_indices.is_empty()
                || !mol
                    .stereo_group_bond_indices(*group_index)
                    .unwrap_or_default()
                    .is_empty()
        })
        .collect::<Vec<_>>();
    if groups.is_empty() {
        return;
    }
    out.push_str("M  V30 BEGIN COLLECTION\n");
    for (group_index, group) in groups {
        let atom_key = match &group.kind {
            StereoGroupKind::Absolute => "MDLV30/STEABS".to_string(),
            StereoGroupKind::Or(number) => format!("MDLV30/STEREL{number}"),
            StereoGroupKind::And(number) => format!("MDLV30/STERAC{number}"),
        };
        let bond_key = match &group.kind {
            StereoGroupKind::Absolute => "MDLV30/STEBABS".to_string(),
            StereoGroupKind::Or(number) => format!("MDLV30/STEBREL{number}"),
            StereoGroupKind::And(number) => format!("MDLV30/STEBRAC{number}"),
        };
        let atoms = group
            .atom_indices
            .iter()
            .map(|atom| (atom.0 + 1).to_string())
            .collect::<Vec<_>>();
        if !atoms.is_empty() {
            out.push_str(&format!(
                "M  V30 {atom_key} ATOMS=({} {})\n",
                atoms.len(),
                atoms.join(" ")
            ));
        }
        let bonds = mol
            .stereo_group_bond_indices(group_index)
            .unwrap_or_default()
            .iter()
            .map(|bond| (bond.0 + 1).to_string())
            .collect::<Vec<_>>();
        if !bonds.is_empty() {
            out.push_str(&format!(
                "M  V30 {bond_key} BONDS=({} {})\n",
                bonds.len(),
                bonds.join(" ")
            ));
        }
    }
    out.push_str("M  V30 END COLLECTION\n");
}

fn finish_ctab(out: &mut String, mol: &Molecule, metadata: &MolMetadata) {
    append_sgroups(out, metadata);
    append_collections(out, mol);
    out.push_str("M  V30 END CTAB\nM  END\n");
}

pub fn write_mol_v3000(mol: &Molecule, metadata: &MolMetadata, coords: &[(f64, f64)]) -> String {
    write_mol_v3000_with_stereo_report(mol, metadata, coords).0
}

/// [`write_mol_v3000`] that also reports the stereo the block does not carry
/// (see [`crate::MolStereoLoss`]). An empty `coords` lays the molecule out.
pub fn write_mol_v3000_with_stereo_report(
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) -> (String, crate::MolStereoLoss) {
    let mut out = String::new();
    append_header(&mut out, mol, metadata, "  chematic\n");

    let depiction = crate::stereo_depiction::needs_stereo_depiction(mol)
        .then(|| crate::stereo_depiction::stereo_depiction(mol, coords));
    let coords = depiction.as_ref().map_or(coords, |value| &value.coords);

    append_atoms_2d(&mut out, mol, metadata, coords);
    append_bonds_2d(&mut out, mol, metadata, depiction.as_ref());
    finish_ctab(&mut out, mol, metadata);

    let loss = crate::MolStereoLoss::from_depiction(
        mol,
        depiction.as_ref(),
        crate::mol2000::MolFormat::V3000,
    );
    (out, loss)
}

/// Serialize `mol` to MOL V3000 (Extended Ctab) format using `conformer`'s
/// real 3D coordinates, stamping the header's line-2 dimensional code as
/// `3D` -- the V3000 counterpart of
/// [`crate::mol2000::write_mol_with_conformer`]. As with that function, no
/// wedge/hash `CFG` is ever emitted on a bond line here: a real 3D geometry
/// makes a 2D wedge symbol redundant at best, and round-tripping this output
/// back through [`read_mol_v3000_with_diagnostics`] must not manufacture a
/// fresh [`crate::mol2000::Stereo3DDiagnostic::WedgeVs3DParityConflict`] on
/// its own output. Enhanced stereo groups (`COLLECTION`/`STEABS`/`STEREL`/
/// `STERAC`) are unaffected -- they label which atoms form a stereo group,
/// not a direction, and remain meaningful for a 3D record.
///
/// **Does not validate `Chirality::SquarePlanar` stereo against
/// `conformer`** -- it writes whatever coordinates it is given, trusting the
/// caller. If `conformer`'s geometry doesn't actually match a declared
/// square-planar tag, this will silently write a self-inconsistent file.
/// Use [`write_mol_v3000_with_conformer_checked`] instead to fail closed on
/// that mismatch (or on a missing/flat conformer) rather than trusting it.
pub fn write_mol_v3000_with_conformer(
    mol: &Molecule,
    metadata: &MolMetadata,
    conformer: &Coords3D,
) -> String {
    let mut out = String::new();
    append_header(&mut out, mol, metadata, "  chematic          3D\n");
    append_atoms_3d(&mut out, mol, metadata, conformer);
    append_bonds_3d(&mut out, mol, metadata);
    finish_ctab(&mut out, mol, metadata);
    out
}

/// [`write_mol_v3000_with_conformer`], but fails closed with a typed
/// [`crate::mol2000::MolStereoWriteError`] instead of silently writing
/// coordinates that don't actually match a molecule's declared
/// square-planar stereo (or don't exist at all) -- the V3000 counterpart of
/// [`crate::mol2000::write_mol_with_conformer_checked`]. See
/// [`crate::mol2000::validate_square_planar_for_write`].
pub fn write_mol_v3000_with_conformer_checked(
    mol: &Molecule,
    metadata: &MolMetadata,
    conformer: &Coords3D,
) -> Result<String, crate::mol2000::MolStereoWriteError> {
    crate::mol2000::validate_square_planar_for_write(
        mol,
        Some(conformer),
        crate::mol2000::MolFormat::V3000,
    )?;
    Ok(write_mol_v3000_with_conformer(mol, metadata, conformer))
}

// Test-only access to the frozen pre-optimization property renderer.
pub(super) fn properties_for_test(
    metadata: &MolMetadata,
    bond_id: u32,
    cfg_written: bool,
) -> String {
    bond_properties(metadata, bond_id, cfg_written)
}
