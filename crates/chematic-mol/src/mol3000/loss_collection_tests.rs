//! Frozen source-main writer; plain output and public loss fields must match.
use super::*;

fn frozen_reporting(
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) -> (String, crate::MolStereoLoss) {
    let mut out = writer_buffer(mol, metadata);
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

fn compare(mol: &Molecule, metadata: &MolMetadata, coords: &[(f64, f64)]) {
    let reference = mol.clone();
    let candidate = mol.clone();
    for pass in 0..2 {
        let expected = frozen_reporting(&reference, metadata, coords);
        let plain = write_mol_v3000(&candidate, metadata, coords);
        let reporting = write_mol_v3000_with_stereo_report(&candidate, metadata, coords);
        assert_eq!(plain, expected.0, "plain bytes, pass={pass}");
        assert_eq!(reporting.0, expected.0, "reporting bytes, pass={pass}");
        assert_eq!(reporting.1, expected.1, "all loss fields, pass={pass}");
    }
}

#[test]
fn plain_v3000_bytes_and_reported_losses_match_on_all_pinned_profiles() {
    use flate2::read::GzDecoder;
    use std::io::Read;
    let corpora: [(&[u8], usize); 3] = [
        (
            include_bytes!(
                "../../../../validation/benchmark_corpora/openbabel-v3000-speed-slice0-1000-v1.molblocks.gz"
            ),
            1000,
        ),
        (
            include_bytes!(
                "../../../../validation/benchmark_corpora/openbabel-v3000-speed-slice5000-1000-v1.molblocks.gz"
            ),
            1000,
        ),
        (
            include_bytes!(
                "../../../../validation/benchmark_corpora/openbabel-v3000-speed-chains64-v1.molblocks.gz"
            ),
            64,
        ),
    ];
    let mut rows = 0;
    for (corpus, expected) in corpora {
        let mut input = String::new();
        GzDecoder::new(corpus).read_to_string(&mut input).unwrap();
        let mut count = 0;
        for block in input
            .split("$$$$\n")
            .filter(|block| !block.trim().is_empty())
        {
            let report = crate::read_mol_v3000_with_diagnostics(block).unwrap();
            compare(&report.mol, &report.metadata, &report.coords);
            count += 1;
        }
        assert_eq!(count, expected);
        rows += count;
    }
    assert_eq!(rows, 2064);
}

#[test]
fn v2000_read_graphs_keep_v3000_plain_bytes_and_full_loss_reports() {
    use flate2::read::GzDecoder;
    use std::io::{BufReader, Cursor};
    let corpora: [&[u8]; 2] = [
        include_bytes!(
            "../../../../validation/benchmark_corpora/openbabel-sdf-speed-1000-v1.sdf.gz"
        ),
        include_bytes!(
            "../../../../validation/benchmark_corpora/openbabel-sdf-speed-slice5000-1000-v1.sdf.gz"
        ),
    ];
    let mut rows = 0;
    for corpus in corpora {
        for record in
            crate::sdf::SdfFileReader::new(BufReader::new(GzDecoder::new(Cursor::new(corpus))))
        {
            let record = record.unwrap();
            compare(&record.mol, &record.meta, &record.coords);
            rows += 1;
        }
    }
    assert_eq!(rows, 2000);
}

#[test]
fn reporting_api_keeps_losses_and_groups_with_empty_and_supplied_coordinates() {
    compare(
        &MoleculeBuilder::new().build(),
        &MolMetadata::default(),
        &[],
    );
    for smiles in [
        "CCO",
        "C[C@H](F)Cl",
        "F/C=C/F",
        "F/C=C\\F",
        "N->[Pt@SP1](<-N)(Cl)Cl",
    ] {
        let mut mol = chematic_smiles::parse(smiles).unwrap();
        mol.set_stereo_groups(vec![
            StereoGroup::new(StereoGroupKind::Absolute, vec![AtomIdx(0)]),
            StereoGroup::new(StereoGroupKind::Or(2), vec![AtomIdx(1)]),
        ]);
        let metadata = MolMetadata {
            name: "Unicode Δ\nname".into(),
            comment: "comment".into(),
            ..Default::default()
        };
        compare(&mol, &metadata, &[]);
        compare(&mol, &metadata, &vec![(0.0, 0.0); mol.atom_count()]);
    }
    let mol = chematic_smiles::parse("N->[Pt@SP1](<-N)(Cl)Cl").unwrap();
    let loss = write_mol_v3000_with_stereo_report(&mol, &MolMetadata::default(), &[]).1;
    assert!(!loss.non_tetrahedral_centres.is_empty());
    assert!(!loss.stereo_groups_dropped);
}
