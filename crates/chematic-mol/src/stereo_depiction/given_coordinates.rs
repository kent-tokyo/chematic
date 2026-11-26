use super::*;

/// Runs `draw` on `mol` with its stereo drawn on caller-supplied
/// coordinates, which are kept: `layout` is in depiction units (y down), for
/// example a file's own drawing scaled to [`chematic_depict::BOND_LEN`] per
/// bond. Each tetrahedral centre gets one wedge or hash chosen on these
/// coordinates and checked by re-reading it, as the MOL writers do for given
/// coordinates, so a centre whose wedge was on a hydrogen that
/// `remove_hydrogens` dropped is drawn with a wedge again. `mol`'s own wedges
/// and `/` `\` marks are replaced; a centre no wedge can express on these
/// coordinates is drawn without one. A molecule without stereo, or a layout
/// that does not place every atom, is drawn as given.
pub fn with_stereo_depiction_on<R>(
    mol: &Molecule,
    layout: &chematic_depict::Layout,
    draw: impl FnOnce(&Molecule, &chematic_depict::Layout) -> R,
) -> R {
    let n = mol.atom_count();
    let placed = layout.coords.len() >= n
        && layout.coords[..n]
            .iter()
            .all(|p| p.x.is_finite() && p.y.is_finite())
        && mol.bonds().all(|(_, bond)| {
            let a = layout.coords[bond.atom1.0 as usize];
            let b = layout.coords[bond.atom2.0 as usize];
            a.x != b.x || a.y != b.y
        });
    if !placed || !needs_stereo_depiction(mol) {
        return draw(mol, layout);
    }
    let scale = 1.5 / chematic_depict::layout::BOND_LEN;
    let coords: Vec<(f64, f64)> = layout.coords[..n]
        .iter()
        .map(|p| (p.x * scale, -p.y * scale))
        .collect();
    let depiction = stereo_depiction(mol, &coords);
    draw(&without_bond_marks(mol, &depiction.wedges), layout)
}

#[cfg(test)]
mod given_coordinates_tests {
    use super::*;
    use crate::mol2000::read_mol_with_diagnostics;
    use chematic_smiles::{canonical_smiles, parse};

    /// `drawn`'s tetrahedral stereo as a reader takes it from its wedges on
    /// `layout` alone (declared chirality cleared first).
    fn reread_wedges(drawn: &Molecule, layout: &chematic_depict::Layout) -> Molecule {
        let scale = 1.5 / chematic_depict::layout::BOND_LEN;
        let coords: Vec<_> = layout
            .coords
            .iter()
            .map(|p| (p.x * scale, -p.y * scale))
            .collect();
        let mut copy = drawn.clone();
        for i in 0..copy.atom_count() {
            let atom = AtomIdx(i as u32);
            copy.set_chirality(atom, chematic_core::Chirality::None);
            copy.clear_stereo_neighbor_order(atom);
        }
        chematic_perception::apply_ctab_local_parity_from_wedges_with_diagnostics(
            &mut copy, &coords,
        );
        copy
    }

    fn wedged_centres(mol: &Molecule) -> Vec<AtomIdx> {
        mol.bonds()
            .filter(|(_, b)| matches!(b.order, BondOrder::Up | BondOrder::Down))
            .map(|(_, b)| b.atom1)
            .collect()
    }

    /// L-alanine as a drawing program writes it: explicit alpha hydrogen,
    /// whose bond carries the centre's only wedge.
    const ALANINE_WEDGE_ON_H: &str = "L-alanine
  hand-written

  7  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  4  5  2  0
  4  6  1  0
  1  7  1  1
M  END
";

    #[test]
    fn file_drawing_without_its_hydrogens_gets_a_wedge_at_every_centre() {
        let read = read_mol_with_diagnostics(ALANINE_WEDGE_ON_H).unwrap();
        let expected = canonical_smiles(&parse("N[C@@H](C)C(=O)O").unwrap());
        let stripped = chematic_chem::hydrogen::remove_hydrogens(&read.mol);
        assert_eq!(canonical_smiles(&stripped), expected);
        // The heavy atoms keep their indices: the hydrogen is the last atom.
        let scale = 1.5 / chematic_depict::layout::BOND_LEN;
        let layout = chematic_depict::Layout {
            coords: read.coords[..stripped.atom_count()]
                .iter()
                .map(|&(x, y)| chematic_depict::Point::new(x / scale, -y / scale))
                .collect(),
        };
        // The bare drawing has lost the only wedge.
        assert!(wedged_centres(&stripped).is_empty());
        with_stereo_depiction_on(&stripped, &layout, |drawn, given| {
            assert_eq!(given.coords, layout.coords);
            assert_eq!(wedged_centres(drawn), vec![AtomIdx(0)]);
            assert_eq!(canonical_smiles(&reread_wedges(drawn, given)), expected);
        });
    }

    #[test]
    fn given_coordinates_are_kept_and_every_centre_reads_back() {
        for text in [
            "F[C@H](Cl)Br",
            "C[C@H]1CC[C@@H](O)CC1",
            "C[C@@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@@]43C)[C@@H]1CC[C@@H]2O",
            "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
            "C/C=C/[C@H](O)CC",
        ] {
            let mol = parse(text).unwrap();
            let expected = canonical_smiles(&mol);
            // A drawing on which every centre can be wedged: the stereo
            // layout, reused as the caller's coordinates.
            let (_, layout) = depiction_with_stereo(&mol);
            let centres = mol
                .atoms()
                .filter(|(_, a)| a.chirality.is_tetrahedral())
                .count();
            with_stereo_depiction_on(&mol, &layout, |drawn, given| {
                assert_eq!(given.coords, layout.coords, "{text}");
                assert_eq!(wedged_centres(drawn).len(), centres, "{text}");
                let reread = reread_wedges(drawn, given);
                let tetrahedral = |m: &Molecule| {
                    m.atoms()
                        .map(|(_, a)| a.chirality.is_tetrahedral())
                        .collect::<Vec<_>>()
                };
                assert_eq!(tetrahedral(&reread), tetrahedral(&mol), "{text}");
                let without_ez = |s: String| s.replace(['/', '\\'], "");
                assert_eq!(
                    without_ez(canonical_smiles(&reread)),
                    without_ez(expected.clone()),
                    "{text}"
                );
            });
        }
    }

    #[test]
    fn invalid_coordinates_preserve_the_original_molecule() {
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let mut layout = chematic_depict::compute_layout(&mol);
        layout.coords[0].x = f64::NAN;
        with_stereo_depiction_on(&mol, &layout, |drawn, _| {
            assert!(std::ptr::eq(drawn, &mol));
        });
        let mut layout = chematic_depict::compute_layout(&mol);
        let (_, bond) = mol.bonds().next().unwrap();
        layout.coords[bond.atom2.0 as usize] = layout.coords[bond.atom1.0 as usize];
        with_stereo_depiction_on(&mol, &layout, |drawn, _| {
            assert!(std::ptr::eq(drawn, &mol));
        });
    }

    #[test]
    fn molecule_without_stereo_or_layout_is_drawn_as_given() {
        let mol = parse("CC(O)CC").unwrap();
        let layout = chematic_depict::compute_layout(&mol);
        with_stereo_depiction_on(&mol, &layout, |drawn, given| {
            assert!(std::ptr::eq(drawn, &mol));
            assert_eq!(given.coords, layout.coords);
        });
        let chiral = parse("F[C@H](Cl)Br").unwrap();
        let empty = chematic_depict::Layout { coords: vec![] };
        with_stereo_depiction_on(&chiral, &empty, |drawn, _| {
            assert!(std::ptr::eq(drawn, &chiral));
        });
    }
}
