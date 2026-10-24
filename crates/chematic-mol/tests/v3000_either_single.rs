//! Open Babel #2429: retain a wavy carrier and do not infer adjacent E/Z.
use chematic_core::{AtomIdx, BondOrder, Coords3D, Point3};
use chematic_mol::{
    read_mol_v3000_with_diagnostics, write_mol_v3000, write_mol_v3000_with_conformer,
};

const INPUT: &str = include_str!("../../../validation/openbabel-2429-2026-10-24/either-single.mol");

#[test]
fn parsed_unknown_alkene_remains_unknown_after_writing_real_coordinates() {
    let mol = read_mol_v3000_with_diagnostics(INPUT).unwrap().mol;
    let coords = [(-1., 1.), (0., 0.), (1., 0.), (2., -1.)];
    let text = write_mol_v3000(&mol, &Default::default(), &coords);
    assert!(text.contains("2 2 2 3 CFG=2"));
    let read = read_mol_v3000_with_diagnostics(&text).unwrap();
    assert!(chematic_chem::assign_ez_bonds(&read.mol).is_empty());
    let partial = write_mol_v3000(&mol, &Default::default(), &[(-1., 1.)]);
    let partial = read_mol_v3000_with_diagnostics(&partial).unwrap();
    assert_eq!(partial.coords, [(-1., 1.), (0., 0.), (0., 0.), (0., 0.)]);

    // A molecule that never carried an explicit unknown-stereo declaration
    // remains free to express the geometry supplied by its caller.
    let plain = chematic_smiles::parse("CC=CC").unwrap();
    let text = write_mol_v3000(&plain, &Default::default(), &coords);
    let read = read_mol_v3000_with_diagnostics(&text).unwrap();
    assert!(!chematic_chem::assign_ez_bonds(&read.mol).is_empty());

    for smiles in ["C/C=C/C", "C/C=C\\C"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let text = write_mol_v3000(&mol, &Default::default(), &[]);
        let read = read_mol_v3000_with_diagnostics(&text).unwrap();
        assert_eq!(
            chematic_smiles::canonical_smiles(&read.mol),
            chematic_smiles::canonical_smiles(&mol)
        );
    }
}

#[test]
fn wavy_single_does_not_turn_unknown_alkene_into_defined_ez() {
    let report = read_mol_v3000_with_diagnostics(INPUT).unwrap();
    let smiles = chematic_smiles::canonical_smiles(&report.mol);
    assert!(!smiles.contains(['/', '\\']), "{smiles}");
    assert!(
        report
            .metadata
            .v3000_bond_properties
            .iter()
            .any(|(id, p)| *id == 3 && p == "CFG=2")
    );
    let output = write_mol_v3000(&report.mol, &report.metadata, &report.coords);
    assert!(output.contains("M  V30 3 1 3 4 CFG=2"));
    let reread = read_mol_v3000_with_diagnostics(&output).unwrap();
    assert_eq!(chematic_smiles::canonical_smiles(&reread.mol), smiles);
    assert_eq!(
        write_mol_v3000(&reread.mol, &reread.metadata, &reread.coords),
        output
    );
}

#[test]
fn either_single_keeps_sparse_ids_orientation_continuations_and_other_attributes() {
    for (before, after, expected) in [
        ("70 1 30 40 CFG=2", "70 1 30 40 CFG=2", "3 1 3 4 CFG=2"),
        ("70 1 30 40 CFG=2", "70 1 40 30 CFG=2", "3 1 4 3 CFG=2"),
        (
            "70 1 30 40 CFG=2",
            "70 1 30 40 -\nM  V30 CFG=2 TOPO=1",
            "3 1 3 4 CFG=2 TOPO=1",
        ),
        ("50 1 10 20", "50 1 10 20 CFG=2", "1 1 1 2 CFG=2"),
    ] {
        let input = INPUT.replace(before, after);
        let report = read_mol_v3000_with_diagnostics(&input).unwrap();
        assert!(!chematic_smiles::canonical_smiles(&report.mol).contains(['/', '\\']));
        let output = write_mol_v3000(&report.mol, &report.metadata, &report.coords);
        assert!(output.contains(expected), "{output}");
    }
}

#[test]
fn unmarked_geometry_and_direct_double_either_keep_existing_behavior() {
    let ordinary = INPUT.replace(" CFG=2", "");
    let report = read_mol_v3000_with_diagnostics(&ordinary).unwrap();
    assert!(chematic_smiles::canonical_smiles(&report.mol).contains(['/', '\\']));
    let direct = ordinary.replace("60 2 20 30", "60 2 20 30 CFG=2");
    let report = read_mol_v3000_with_diagnostics(&direct).unwrap();
    assert!(!chematic_smiles::canonical_smiles(&report.mol).contains(['/', '\\']));
    assert!(
        write_mol_v3000(&report.mol, &report.metadata, &report.coords).contains("2 2 2 3 CFG=2")
    );
}

#[test]
fn three_d_writer_retains_the_explicit_either_carrier() {
    let input = INPUT
        .replace(" 2D", " 3D")
        .replace("40 C 2 -1 0", "40 C 2 -1 0.5");
    let report = read_mol_v3000_with_diagnostics(&input).unwrap();
    assert!(!chematic_smiles::canonical_smiles(&report.mol).contains(['/', '\\']));
    let conformer = Coords3D {
        points: vec![
            Point3::new(-1., 1., 0.),
            Point3::new(0., 0., 0.),
            Point3::new(1., 0., 0.),
            Point3::new(2., -1., 0.5),
        ],
    };
    let output = write_mol_v3000_with_conformer(&report.mol, &report.metadata, &conformer);
    assert!(output.contains("3 1 3 4 CFG=2"));
    let reread = read_mol_v3000_with_diagnostics(&output).unwrap();
    assert!(!chematic_smiles::canonical_smiles(&reread.mol).contains(['/', '\\']));
    assert_eq!(
        reread
            .mol
            .bond_between(AtomIdx(2), AtomIdx(3))
            .map(|(_, bond)| bond.order),
        Some(BondOrder::Single)
    );
}

#[test]
fn wavy_carrier_invalidates_both_adjacent_alkenes_but_not_remote_ones() {
    let expanded = INPUT
        .replace("COUNTS 4 3", "COUNTS 6 5")
        .replace(
            "M  V30 END ATOM",
            "M  V30 50 C 3 -1 0 0\nM  V30 60 C 4 -2 0 0\nM  V30 END ATOM",
        )
        .replace(
            "M  V30 END BOND",
            "M  V30 80 2 40 50\nM  V30 90 1 50 60\nM  V30 END BOND",
        );
    let unknown = read_mol_v3000_with_diagnostics(&expanded).unwrap();
    assert!(chematic_chem::assign_ez_bonds(&unknown.mol).is_empty());
    let specified = read_mol_v3000_with_diagnostics(&expanded.replace(" CFG=2", "")).unwrap();
    assert_eq!(chematic_chem::assign_ez_bonds(&specified.mol).len(), 2);

    let remote = INPUT
        .replace(" CFG=2", "")
        .replace("COUNTS 4 3", "COUNTS 6 4")
        .replace(
            "M  V30 END ATOM",
            "M  V30 50 C 3 3 0 0\nM  V30 60 C 4 3 0 0\nM  V30 END ATOM",
        )
        .replace(
            "M  V30 END BOND",
            "M  V30 80 1 50 60 CFG=2\nM  V30 END BOND",
        );
    let report = read_mol_v3000_with_diagnostics(&remote).unwrap();
    assert_eq!(chematic_chem::assign_ez_bonds(&report.mol).len(), 1);
}
