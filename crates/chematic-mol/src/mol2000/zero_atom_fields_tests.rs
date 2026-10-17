//! Byte preservation against the complete c84c68f2 atom-row writer.
use super::*;

fn compare_row(mol: &Molecule, idx: AtomIdx, coords: Option<(f64, f64)>) {
    let mut expected = String::from("境界:");
    let mut actual = expected.clone();
    frozen_atom_row(&mut expected, mol, idx, mol.atom(idx), coords);
    write_atom_row(&mut actual, mol, idx, mol.atom(idx), coords);
    assert_eq!(
        actual,
        expected,
        "atom={:?}, coords={coords:?}",
        mol.atom(idx)
    );
}

#[test]
fn atom_field_boundaries_preserve_zero_nonzero_and_overflow_bytes() {
    for element in [
        Element::C,
        Element::N,
        Element::O,
        Element::CL,
        Element::BR,
        Element::FE,
    ] {
        for isotope in [None, Some(0), Some(12), Some(13), Some(u16::MAX)] {
            for charge in [-4, -3, -1, 0, 1, 3, 4] {
                for atom_map in [
                    None,
                    Some(0),
                    Some(1),
                    Some(999),
                    Some(1000),
                    Some(u16::MAX),
                ] {
                    for hydrogen_count in [None, Some(0), Some(1), Some(4), Some(14), Some(15)] {
                        for wildcard in [false, true] {
                            let mut atom = Atom::new(element);
                            atom.isotope = isotope;
                            atom.charge = charge;
                            atom.atom_map = atom_map;
                            atom.hydrogen_count = hydrogen_count;
                            atom.wildcard = wildcard;
                            let mut builder = MoleculeBuilder::new();
                            let idx = builder.add_atom(atom);
                            let mol = builder.build();
                            for coords in [
                                None,
                                Some((0.0, -0.0)),
                                Some((-123.12345, 1e20)),
                                Some((f64::NAN, f64::INFINITY)),
                            ] {
                                compare_row(&mol, idx, coords);
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn every_atom_of_normal_and_unusual_corpora_keeps_its_complete_row() {
    let corpora = [
        include_str!("../../../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"),
        include_str!("../../../../validation/cosmolkit_comparison/corpora/stress-v1.smi"),
    ];
    let mut parsed = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let smiles = line.split_whitespace().next().unwrap();
            let mol = chematic_smiles::parse(smiles).unwrap();
            parsed += 1;
            for (idx, _) in mol.atoms() {
                for coords in [None, Some((0.0, -0.0)), Some((1.23445, -2.34565))] {
                    compare_row(&mol, idx, coords);
                }
            }
        }
    }
    assert_eq!(parsed, 14_072);
}

#[test]
fn r_group_symbols_keep_labels_and_exact_field_offsets() {
    for label in [
        None,
        Some(RGroupLabel::unnumbered()),
        RGroupLabel::numbered(1),
        RGroupLabel::numbered(9999),
    ] {
        let mut builder = MoleculeBuilder::new();
        let mut atom = Atom::new(Element::C);
        atom.wildcard = true;
        let idx = builder.add_atom(atom);
        let mut mol = builder.build();
        if let Some(label) = label {
            mol.set_r_group(idx, label);
        }
        compare_row(&mol, idx, None);
    }
}

// Complete baseline function, with its declaration renamed only.
fn frozen_atom_row(
    out: &mut String,
    mol: &Molecule,
    idx: AtomIdx,
    atom: &Atom,
    coords: Option<(f64, f64)>,
) {
    let symbol = mdl_atom_symbol(mol, idx, atom);
    let charge = encode_charge(atom.charge);
    let mass = encode_mass_difference(atom.element, atom.isotope).unwrap_or(0);
    let atom_map = atom.atom_map.unwrap_or(0);
    let valence = valence_field(mol, idx, atom);
    if let Some((x, y)) = coords {
        push_fixed4(out, x, 10);
        push_fixed4(out, y, 10);
        out.push_str("    0.0000 ");
    } else {
        out.push_str("    0.0000    0.0000    0.0000 ");
    }
    out.push_str(symbol);
    for _ in symbol.len()..3 {
        out.push(' ');
    }
    if mass == 0 && charge == 0 {
        out.push_str(" 0  0");
    } else {
        push_right_aligned_i16(out, mass, 2);
        push_right_aligned_u32(out, charge as u32, 3);
    }
    out.push_str("  0  0  0");
    push_right_aligned_u32(out, u32::from(valence), 3);
    out.push_str("  0  0  0");
    push_right_aligned_u32(out, atom_map as u32, 3);
    out.push_str("  0\n");
}
