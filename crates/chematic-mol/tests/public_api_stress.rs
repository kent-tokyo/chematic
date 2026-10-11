//! Cross-format panic-safety over the checked-in special-chemistry corpus.

use chematic_mol::{
    MolMetadata, parse_mol, parse_mol_v3000, read_mol_v3000_with_diagnostics,
    read_mol_with_diagnostics, write_mol_v3000_with_stereo_report, write_mol_with_stereo_report,
};
use chematic_smiles::parse;

const STRESS: &str = include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi");

#[test]
fn mol_writers_and_readers_are_panic_free_on_special_chemistry() {
    let metadata = MolMetadata::default();
    let mut exercised = 0usize;
    for line in STRESS.lines() {
        if exercised == 512 {
            break;
        }
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() == 0 || mol.atom_count() > 120 {
            continue;
        }
        exercised += 1;

        let (v2000, loss2000) = write_mol_with_stereo_report(&mol, &metadata, &[]);
        assert!(v2000.ends_with("M  END\n"));
        assert!(loss2000.centres.len() <= mol.atom_count());
        let _ = read_mol_with_diagnostics(&v2000);
        let _ = parse_mol(&v2000);

        let (v3000, loss3000) = write_mol_v3000_with_stereo_report(&mol, &metadata, &[]);
        assert!(v3000.contains("V3000"));
        assert!(loss3000.centres.len() <= mol.atom_count());
        let _ = read_mol_v3000_with_diagnostics(&v3000);
        let _ = parse_mol_v3000(&v3000);
    }
    assert_eq!(exercised, 512);
}
