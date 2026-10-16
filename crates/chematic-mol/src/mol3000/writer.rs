use super::*;
use std::fmt::Write as _;

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
        write!(line, " RGROUPS=(1 {number})").expect("writing to String cannot fail");
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
    writeln!(
        out,
        "M  V30 COUNTS {} {} {} 0 0",
        mol.atom_count(),
        mol.bond_count(),
        metadata.v3000_sgroups.len()
    )
    .expect("writing to String cannot fail");
}

fn append_atom_properties(
    line: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
    idx: AtomIdx,
    atom: &Atom,
) {
    if atom.charge != 0 {
        write!(line, " CHG={}", atom.charge).expect("writing to String cannot fail");
    }
    if let Some(isotope) = atom.isotope {
        write!(line, " MASS={isotope}").expect("writing to String cannot fail");
    }
    if let Some(hydrogens) = atom.hydrogen_count {
        write!(line, " HCOUNT={hydrogens}").expect("writing to String cannot fail");
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
        out.push_str("M  V30 ");
        crate::mol2000::push_right_aligned_u32(out, idx.0 + 1, 0);
        out.push(' ');
        out.push_str(v3000_atom_symbol(mol, idx, atom));
        out.push(' ');
        crate::mol2000::push_fixed4(out, x, 0);
        out.push(' ');
        crate::mol2000::push_fixed4(out, y, 0);
        if atom.atom_map.unwrap_or(0) == 0 {
            out.push_str(" 0.0000 0");
        } else {
            write!(out, " 0.0000 {}", atom.atom_map.unwrap_or(0))
                .expect("writing to String cannot fail");
        }
        append_atom_properties(out, mol, metadata, idx, atom);
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
        write!(
            out,
            "M  V30 {} {} {:.4} {:.4} {:.4} {}",
            idx.0 + 1,
            v3000_atom_symbol(mol, idx, atom),
            point.x,
            point.y,
            point.z,
            atom.atom_map.unwrap_or(0)
        )
        .expect("writing to String cannot fail");
        append_atom_properties(out, mol, metadata, idx, atom);
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

fn append_bond_properties(
    out: &mut String,
    metadata: &MolMetadata,
    bond_id: u32,
    cfg_written: bool,
) {
    if let Some((_, properties)) = metadata
        .v3000_bond_properties
        .iter()
        .find(|(id, _)| *id == bond_id)
    {
        if !cfg_written {
            out.push(' ');
            out.push_str(properties);
        } else {
            // Match the existing CFG deduplication, including its whitespace
            // normalization. Append surviving tokens without a temporary Vec.
            for token in properties
                .split_whitespace()
                .filter(|token| !token.starts_with("CFG="))
            {
                out.push(' ');
                out.push_str(token);
            }
        }
    }
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
                    let cfg = if wedge.order == BondOrder::Up {
                        " CFG=1"
                    } else {
                        " CFG=3"
                    };
                    (start, end, cfg)
                }
                None if depiction.unexpressed_double_bonds.contains(&bond_idx)
                    || depiction.unspecified_double_bonds.contains(&bond_idx) =>
                {
                    (atom1, atom2, " CFG=2")
                }
                None => (atom1, atom2, ""),
            },
            None => {
                let cfg = match bond.order {
                    BondOrder::Up => " CFG=1",
                    BondOrder::Down => " CFG=3",
                    _ => "",
                };
                (atom1, atom2, cfg)
            }
        };
        let bond_id = bond_idx.0 + 1;
        out.push_str("M  V30 ");
        crate::mol2000::push_right_aligned_u32(out, bond_id, 0);
        out.push(' ');
        out.push(char::from(b'0' + v3000_bond_order(written)));
        out.push(' ');
        crate::mol2000::push_right_aligned_u32(out, atom1, 0);
        out.push(' ');
        crate::mol2000::push_right_aligned_u32(out, atom2, 0);
        out.push_str(stereo);
        append_bond_properties(out, metadata, bond_id, !stereo.is_empty());
        out.push('\n');
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
        write!(
            out,
            "M  V30 {bond_id} {} {} {}",
            v3000_bond_order(written),
            bond.atom1.0 + 1,
            bond.atom2.0 + 1
        )
        .expect("writing to String cannot fail");
        append_bond_properties(out, metadata, bond_id, false);
        out.push('\n');
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

fn writer_buffer(mol: &Molecule, metadata: &MolMetadata) -> String {
    // Reserve ordinary atom/bond rows once. Large coordinates and opaque
    // metadata can still grow the buffer; no output or input limit changes.
    String::with_capacity(
        256 + metadata.name.len()
            + metadata.comment.len()
            + mol.atom_count() * 64
            + mol.bond_count() * 40,
    )
}

fn finish_ctab(out: &mut String, mol: &Molecule, metadata: &MolMetadata) {
    append_sgroups(out, metadata);
    append_collections(out, mol);
    out.push_str("M  V30 END CTAB\nM  END\n");
}

pub fn write_mol_v3000(mol: &Molecule, metadata: &MolMetadata, coords: &[(f64, f64)]) -> String {
    write_v3000_reporting(mol, metadata, coords, false).0
}

/// [`write_mol_v3000`] that also reports the stereo the block does not carry
/// (see [`crate::MolStereoLoss`]). An empty `coords` lays the molecule out.
pub fn write_mol_v3000_with_stereo_report(
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) -> (String, crate::MolStereoLoss) {
    write_v3000_reporting(mol, metadata, coords, true)
}

fn write_v3000_reporting(
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
    collect_loss: bool,
) -> (String, crate::MolStereoLoss) {
    let mut out = writer_buffer(mol, metadata);
    append_header(&mut out, mol, metadata, "  chematic\n");

    let depiction = crate::stereo_depiction::needs_stereo_depiction(mol)
        .then(|| crate::stereo_depiction::stereo_depiction(mol, coords));
    let coords = depiction.as_ref().map_or(coords, |value| &value.coords);

    append_atoms_2d(&mut out, mol, metadata, coords);
    append_bonds_2d(&mut out, mol, metadata, depiction.as_ref());
    finish_ctab(&mut out, mol, metadata);

    // Plain output discards this report. Keep depiction and serialization
    // identical, collecting losses only for the public reporting API.
    let loss = if collect_loss {
        crate::MolStereoLoss::from_depiction(
            mol,
            depiction.as_ref(),
            crate::mol2000::MolFormat::V3000,
        )
    } else {
        crate::MolStereoLoss::default()
    };
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
    let mut out = writer_buffer(mol, metadata);
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

#[cfg(test)]
mod write_tests {
    use super::*;
    use crate::mol3000::parse_mol_v3000;
    use chematic_core::Element;
    use chematic_core::{Atom, MoleculeBuilder};

    fn ethanol() -> Molecule {
        use chematic_core::BondOrder;
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        let c2 = b.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        let o = b.add_atom(Atom::new(Element::from_symbol("O").unwrap()));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, o, BondOrder::Single).unwrap();
        b.build()
    }

    #[test]
    fn write_v3000_roundtrip_ethanol() {
        let mol = ethanol();
        let meta = MolMetadata {
            name: "ethanol".into(),
            comment: String::new(),
            ..Default::default()
        };
        let v3k = write_mol_v3000(&mol, &meta, &[]);
        let (mol2, meta2) = parse_mol_v3000(&v3k).expect("round-trip parse");
        assert_eq!(mol.atom_count(), mol2.atom_count());
        assert_eq!(mol.bond_count(), mol2.bond_count());
        assert_eq!(meta2.name, "ethanol");
    }

    #[test]
    fn write_v3000_contains_v3000_tag() {
        let mol = ethanol();
        let meta = MolMetadata::default();
        let v3k = write_mol_v3000(&mol, &meta, &[]);
        assert!(v3k.contains("V3000"), "output should contain V3000 tag");
        assert!(
            v3k.contains("M  V30 BEGIN CTAB"),
            "should contain CTAB block"
        );
    }

    #[test]
    fn write_v3000_stereo_group_roundtrip() {
        use chematic_core::{AtomIdx, BondOrder, StereoGroup, StereoGroupKind};

        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        let c2 = b.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        let c3 = b.add_atom(Atom::new(Element::from_symbol("N").unwrap()));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        let mut mol = b.build();

        // Absolute group on atom 0; OR group on atoms 1 and 2; AND group on
        // atoms 0 and 2.
        mol.set_stereo_groups(vec![
            StereoGroup::new(StereoGroupKind::Absolute, vec![c1]),
            StereoGroup::new(StereoGroupKind::Or(1), vec![c2, c3]),
            StereoGroup::new(StereoGroupKind::And(2), vec![c1, c3]),
        ]);

        let meta = MolMetadata {
            name: "stereo_test".into(),
            comment: String::new(),
            ..Default::default()
        };
        let v3k = write_mol_v3000(&mol, &meta, &[]);

        // Verify COLLECTION block is present.
        assert!(
            v3k.contains("BEGIN COLLECTION"),
            "should have COLLECTION block"
        );
        assert!(v3k.contains("MDLV30/STEABS"), "should have STEABS entry");
        assert!(v3k.contains("MDLV30/STEREL1"), "should have STEREL1 entry");
        assert!(v3k.contains("MDLV30/STERAC2"), "should have STERAC2 entry");

        // Round-trip: parse back and verify stereo groups are preserved.
        let (mol2, _) = parse_mol_v3000(&v3k).expect("round-trip parse");
        assert_eq!(
            mol2.stereo_groups().len(),
            3,
            "should have 3 stereo groups after round-trip"
        );

        let abs_group = mol2
            .stereo_groups()
            .iter()
            .find(|g| g.kind == StereoGroupKind::Absolute)
            .expect("Absolute group should exist");
        assert_eq!(abs_group.atom_indices, vec![AtomIdx(0)]);

        let or_group = mol2
            .stereo_groups()
            .iter()
            .find(|g| g.kind == StereoGroupKind::Or(1))
            .expect("Or(1) group should exist");
        assert_eq!(or_group.atom_indices, vec![AtomIdx(1), AtomIdx(2)]);

        let and_group = mol2
            .stereo_groups()
            .iter()
            .find(|g| g.kind == StereoGroupKind::And(2))
            .expect("And(2) group should exist");
        assert_eq!(and_group.atom_indices, vec![AtomIdx(0), AtomIdx(2)]);
    }

    #[test]
    fn test_stereo_group_duplicate_atoms_deduplicated() {
        // RDKit PR #9258: duplicate atom indices in a StereoGroup must be silently
        // removed. In RDKit this caused heap-use-after-free; in chematic the
        // duplicate indices produce incorrect stereo group membership.
        //
        // Craft a minimal V3000 MOL block with STEABS ATOMS=(3 1 1 2)
        // — atom 1 appears twice, which should be collapsed to [AtomIdx(0)].
        let v3k = "\n\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\n\
M  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0 CFG=2\n\
M  V30 2 C 1 0 0 0 CFG=2\n\
M  V30 3 C 2 0 0 0\n\
M  V30 END ATOM\nM  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 2 1 2 3\n\
M  V30 END BOND\nM  V30 BEGIN COLLECTION\n\
M  V30 MDLV30/STEABS ATOMS=(3 1 1 2)\n\
M  V30 END COLLECTION\nM  V30 END CTAB\nM  END\n";

        let (mol, _) = parse_mol_v3000(v3k).expect("should parse despite duplicate atom");
        assert_eq!(mol.stereo_groups().len(), 1);
        let group = &mol.stereo_groups()[0];
        assert_eq!(
            group.kind,
            StereoGroupKind::Absolute,
            "group kind must be Absolute"
        );
        // Atom 1 should appear only ONCE (deduplicated), plus atom 2.
        assert_eq!(
            group.atom_indices.len(),
            2,
            "duplicate atom index must be removed: got {:?}",
            group.atom_indices
        );
        assert!(
            group.atom_indices.contains(&AtomIdx(0)),
            "atom 0 must be in group"
        );
        assert!(
            group.atom_indices.contains(&AtomIdx(1)),
            "atom 1 must be in group"
        );
    }

    #[test]
    fn v3000_sgroup_is_preserved_as_opaque_metadata() {
        // SGroup/polymer semantics are not represented in Molecule yet. Preserve
        // their logical lines opaquely so a core parse/write round trip does not
        // discard them or claim a typed interpretation that does not exist.
        let v3k = "\n\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 2 1 1 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0\n\
M  V30 2 C 1 0 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 END BOND\n\
M  V30 BEGIN SGROUP\n\
M  V30 1 SUP 0 ATOMS=(2 1 2)\n\
M  V30 END SGROUP\n\
M  V30 END CTAB\nM  END\n";

        let (mol, metadata) = parse_mol_v3000(v3k).expect("opaque SGroup should parse");
        assert_eq!(metadata.v3000_sgroups, vec!["1 SUP 0 ATOMS=(2 1 2)"]);
        let rewritten = write_mol_v3000(&mol, &metadata, &[]);
        assert!(rewritten.contains("M  V30 BEGIN SGROUP"));
        assert!(rewritten.contains("M  V30 1 SUP 0 ATOMS=(2 1 2)"));
        assert!(rewritten.find("BEGIN SGROUP").unwrap() < rewritten.find("END CTAB").unwrap());

        let with_collection = v3k.replace(
            "M  V30 END CTAB",
            "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/STEABS ATOMS=(1 1)\nM  V30 END COLLECTION\nM  V30 END CTAB",
        );
        let (_, collection_metadata) =
            parse_mol_v3000(&with_collection).expect("SGROUP and COLLECTION should compose");
        assert_eq!(collection_metadata.v3000_sgroups.len(), 1);
    }

    #[test]
    fn v3000_sgroup_has_strict_typed_view_without_losing_source_line() {
        let typed =
            parse_v3000_sgroup_line("1 COP 0 ATOMS=(2 1 2) LABEL=polymer BRKXYZ=(4 1 2 3 4)")
                .expect("valid SGROUP should have a typed view");
        assert_eq!(typed.id, 1);
        assert_eq!(typed.kind, V3000SGroupKind::Cop);
        assert_eq!(typed.parent_id, None);
        assert_eq!(typed.atom_ids, vec![1, 2]);
        assert_eq!(
            typed.attributes,
            vec![
                ("LABEL".into(), "polymer".into()),
                ("BRKXYZ".into(), "(4 1 2 3 4)".into())
            ]
        );
    }

    #[test]
    fn v3000_sgroup_typed_view_rejects_count_and_id_errors() {
        for line in [
            "0 SUP 0 ATOMS=(1 1)",
            "1 SUP 0 ATOMS=(2 1)",
            "1 SUP 0 ATOMS=(1 0)",
            "1 SUP 0 LABEL=missing_atoms",
        ] {
            assert!(parse_v3000_sgroup_line(line).is_err(), "must reject {line}");
        }
    }
}

#[cfg(test)]
#[path = "writer_reference.rs"]
mod reference;

#[cfg(test)]
#[path = "writer_speed_tests.rs"]
mod speed_tests;

#[cfg(test)]
#[path = "loss_collection_tests.rs"]
mod loss_collection_tests;
