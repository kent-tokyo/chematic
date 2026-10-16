//! Open Babel 3.2.1 SDF rewrites of RDKit 2025.09.3 source-gated 3D records.

use chematic_mol::{
    mol2000::{Stereo3DDiagnostic, read_mol_with_diagnostics},
    read_mol_v3000_with_diagnostics,
};
use chematic_smiles::rdkit_canonical_smiles;

#[test]
fn real_3d_geometry_resolves_conflicting_openbabel_wedges() {
    for (name, v2000, v3000, expected) in [
        (
            "isotope",
            include_str!("fixtures/openbabel-3d-conflict-isotope-v2000.mol"),
            include_str!("fixtures/openbabel-3d-conflict-isotope-v3000.mol"),
            "C[C@](C[13CH3])(c1ccc(O)cc1)[C@H](CCCCBr)c1ccc(O)cc1",
        ),
        (
            "cage",
            include_str!("fixtures/openbabel-3d-conflict-cage-v2000.mol"),
            include_str!("fixtures/openbabel-3d-conflict-cage-v3000.mol"),
            "C[C@@]1(CO)[C@H](O)CC[13C@@]2(C)[C@H]1CC[C@H]1C[C@@H]3C[C@@]12CC[C@]3(O)CO",
        ),
        (
            "ring",
            include_str!("fixtures/openbabel-3d-conflict-ring-v2000.mol"),
            include_str!("fixtures/openbabel-3d-conflict-ring-v3000.mol"),
            "CC1(C)C[C@H]2[C@@H]3OC(C)(C)O[C@@H]3CC[C@@]2(COS(N)(=O)=O)C1",
        ),
    ] {
        let report = read_mol_with_diagnostics(v2000).unwrap();
        assert!(
            report
                .stereo3d_diagnostics
                .iter()
                .any(|diagnostic| matches!(
                    diagnostic,
                    Stereo3DDiagnostic::WedgeVs3DParityConflict { .. }
                )),
            "{name}: conflicting Open Babel V2000 wedge must remain observable"
        );
        assert_eq!(
            rdkit_canonical_smiles(&report.mol).unwrap(),
            expected,
            "{name} V2000"
        );

        let report = read_mol_v3000_with_diagnostics(v3000).unwrap();
        assert_eq!(
            rdkit_canonical_smiles(&report.mol).unwrap(),
            expected,
            "{name} V3000"
        );
    }
}
