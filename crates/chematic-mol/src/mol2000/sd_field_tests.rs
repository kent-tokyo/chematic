use super::*;
use std::collections::HashMap;
use std::fmt::Write as _;

// Frozen original formatting path; the same map is read without mutation so
// randomized map iteration order is identical in the two calls.
fn reference(out: &mut String, props: &HashMap<String, String>) {
    let additional = props
        .iter()
        .filter(|(key, _)| !key.starts_with('_'))
        .map(|(key, value)| key.len().saturating_add(value.len()).saturating_add(7))
        .sum::<usize>()
        .saturating_add(5);
    out.reserve(additional);
    for (key, value) in props {
        if !key.starts_with('_') {
            write!(out, "> <{key}>\n{value}\n\n").expect("writing to String cannot fail");
        }
    }
    out.push_str("$$$$\n");
}

#[test]
fn sd_field_bytes_match_original_for_hidden_unicode_and_literal_values() {
    let keys = [
        "",
        "_",
        "__private",
        "visible",
        "日本語",
        "雪💧",
        "a\n> <B>",
        "_秘密",
        "x\0y",
    ];
    let values = ["", "x", "line1\nline2", "\r\n", "{}\t \0", "値🍎"];
    let mut cases = 0;
    for mask in 0..1_usize << keys.len() {
        for rotation in 0..values.len() {
            let props: HashMap<_, _> = keys
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1 << i) != 0)
                .map(|(i, k)| {
                    (
                        k.to_string(),
                        values[(i + rotation) % values.len()].to_string(),
                    )
                })
                .collect();
            for prefix in ["", "M  END\n", "prior output\0💧"] {
                let mut a = prefix.to_string();
                let mut b = a.clone();
                append_sd_fields_and_delimiter(&mut a, &props);
                reference(&mut b, &props);
                assert_eq!(a.as_bytes(), b.as_bytes());
                cases += 1;
            }
        }
    }
    let props = HashMap::from([("large".to_string(), "値\n".repeat(65_536))]);
    let mut a = String::with_capacity(1);
    let mut b = String::with_capacity(1);
    append_sd_fields_and_delimiter(&mut a, &props);
    reference(&mut b, &props);
    assert_eq!(a, b);
    eprintln!(
        "SD-field exact byte comparisons: {}; large payload: {} bytes",
        cases + 1,
        a.len()
    );
}

fn reference_v2000_reporting(
    out: &mut String,
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
    either: &[chematic_core::BondIdx],
) -> MolStereoLoss {
    // A MOL block is dominated by fixed-width atom and bond rows.  Reserve
    // the common-size output up front so serialization does not repeatedly
    // grow and copy the String for every row (the old empty String was a
    // noticeable cost when writing many small records).
    out.clear();
    out.reserve(128 + mol.atom_count() * 80 + mol.bond_count() * 16);

    // Header lines 1–3
    out.push_str(&metadata.name);
    out.push('\n');
    out.push_str("  chematic\n");
    out.push_str(&metadata.comment);
    out.push('\n');

    // A molecule with stereo is written on coordinates that express it, with
    // one wedge or hash per tetrahedral centre (see `stereo_depiction`).
    let depiction = crate::stereo_depiction::needs_stereo_depiction(mol)
        .then(|| crate::stereo_depiction::stereo_depiction(mol, coords));
    let coords: &[(f64, f64)] = depiction.as_ref().map_or(coords, |d| &d.coords);

    // Counts line (line 4)
    let natoms = mol.atom_count();
    let nbonds = mol.bond_count();
    push_right_aligned_u32(out, natoms as u32, 3);
    push_right_aligned_u32(out, nbonds as u32, 3);
    out.push_str("  0  0  0  0  0  0  0  0999 V2000\n");

    write_atom_rows(out, mol, coords);
    write_bond_rows(out, mol, depiction.as_ref(), either);

    write_v2000_rgroups(out, mol);

    // Terminator
    out.push_str("M  END\n");
    MolStereoLoss::from_depiction(mol, depiction.as_ref(), MolFormat::V2000)
}

fn assert_loss_equal(a: &MolStereoLoss, b: &MolStereoLoss) {
    assert_eq!(a.centres, b.centres);
    assert_eq!(a.double_bonds, b.double_bonds);
    assert_eq!(a.non_tetrahedral_centres, b.non_tetrahedral_centres);
    assert_eq!(a.stereo_groups_dropped, b.stereo_groups_dropped);
}

#[test]
fn plain_bytes_and_reporting_losses_match_frozen_writer_on_pinned_records() {
    use flate2::read::GzDecoder;
    use std::io::{Cursor, Read};
    let corpora: [&[u8]; 2] = [
        include_bytes!(
            "../../../../validation/benchmark_corpora/openbabel-sdf-speed-1000-v1.sdf.gz"
        ),
        include_bytes!(
            "../../../../validation/benchmark_corpora/openbabel-sdf-speed-slice5000-1000-v1.sdf.gz"
        ),
    ];
    let mut comparisons = 0;
    for corpus in corpora {
        let mut input = Vec::new();
        GzDecoder::new(corpus).read_to_end(&mut input).unwrap();
        for record in crate::sdf::SdfFileReader::new(Cursor::new(&input)) {
            let record = record.unwrap();
            let either = if !record.coords.is_empty()
                && !crate::stereo_depiction::needs_stereo_depiction(&record.mol)
            {
                crate::stereo_depiction::undeclared_stereo_double_bonds(&record.mol)
            } else {
                Vec::new()
            };
            let mut old = String::new();
            let expected = reference_v2000_reporting(
                &mut old,
                &record.mol,
                &record.meta,
                &record.coords,
                &either,
            );
            let mut reporting = String::new();
            let observed = write_v2000_reporting(
                &mut reporting,
                &record.mol,
                &record.meta,
                &record.coords,
                &either,
                true,
            );
            let mut plain = String::new();
            write_v2000_reporting(
                &mut plain,
                &record.mol,
                &record.meta,
                &record.coords,
                &either,
                false,
            );
            assert_eq!(plain, old);
            assert_eq!(reporting, old);
            assert_loss_equal(&expected, &observed);
            let mut public_reference = String::new();
            let expected = reference_v2000_reporting(
                &mut public_reference,
                &record.mol,
                &record.meta,
                &record.coords,
                &[],
            );
            let (public_bytes, public_loss) =
                write_mol_with_stereo_report(&record.mol, &record.meta, &record.coords);
            assert_eq!(public_bytes, public_reference);
            assert_loss_equal(&expected, &public_loss);
            comparisons += 1;
        }
    }
    assert_eq!(comparisons, 2000);
    eprintln!("frozen writer byte/loss record comparisons: {comparisons}");
}
