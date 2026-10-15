use std::io::Cursor;

use chematic_mol::{
    SdfFileReader, SdfReader, SdfRecordReader, read_mol_v3000_with_diagnostics,
    read_sdf_conformer_ensembles,
};
use chematic_smiles::rdkit_canonical_smiles;

const MIXED_SDF: &str =
    include_str!("../../../validation/cosmolkit_comparison/corpora/mixed-v2000-v3000-v1.sdf");
const V3000: &str = "M  END in title\n  chematic\nM  END in comment\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 O 1 0 0 0 CHG=-1\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";

#[test]
fn every_sdf_reader_recognizes_mixed_v2000_v3000_and_fields() {
    for input in [MIXED_SDF.to_owned(), MIXED_SDF.replace('\n', "\r\n")] {
        let simple: Vec<_> = SdfReader::new(&input).collect::<Result<_, _>>().unwrap();
        assert_eq!(
            simple
                .iter()
                .map(|(mol, _)| mol.atom_count())
                .collect::<Vec<_>>(),
            [1, 2]
        );
        let records: Vec<_> = SdfRecordReader::new(&input)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(records[1].mol.atom_count(), 2);
        assert_eq!(records[1].mol.bond_count(), 1);
        assert_eq!(records[1].mol.atoms().nth(1).unwrap().1.charge, -1);
        assert_eq!(records[1].meta.name, "M  END in title");
        assert_eq!(records[1].properties["id"], "new");
        assert_eq!(records[1].properties["note"], "M  END is data");

        for fast in [false, true] {
            let cursor = Cursor::new(input.as_bytes());
            let streamed: Vec<_> = if fast {
                SdfFileReader::fast(cursor).collect::<Result<_, _>>()
            } else {
                SdfFileReader::new(cursor).collect::<Result<_, _>>()
            }
            .unwrap();
            assert_eq!(
                streamed
                    .iter()
                    .map(|r| r.mol.atom_count())
                    .collect::<Vec<_>>(),
                [1, 2]
            );
            assert_eq!(streamed[1].properties["id"], "new");
            assert_eq!(streamed[1].properties["note"], "M  END is data");
        }
    }
}

#[test]
fn malformed_v3000_is_not_accepted_by_falling_back_to_v2000() {
    let malformed = V3000.replace("M  END\n", "");
    let input = format!("{malformed}$$$$\n");
    assert!(SdfReader::new(&input).next().unwrap().is_err());
    assert!(SdfRecordReader::new(&input).next().unwrap().is_err());
    assert!(
        SdfFileReader::new(Cursor::new(input.as_bytes()))
            .next()
            .unwrap()
            .is_err()
    );
    assert!(
        SdfFileReader::fast(Cursor::new(input.as_bytes()))
            .next()
            .unwrap()
            .is_err()
    );
}

#[test]
fn v3000_conformers_are_recognized_and_grouped() {
    let conformer = V3000.replace("M  V30 2 O 1 0 0 0", "M  V30 2 O 1 0 1 0");
    let input = format!("{conformer}$$$$\n{conformer}$$$$\n");
    let groups = read_sdf_conformer_ensembles(&input).unwrap();
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].conformers.len(), 2);
    assert_eq!(groups[0].mol.atom_count(), 2);
}

#[test]
fn v3000_tricoordinate_sulfoxide_wedge_retains_stereochemistry() {
    // RDKit 2025.09.3 MolToMolBlock(forceV3000=True) for
    // C[S@@+]([O-])c1ccccc1 after Compute2DCoords. Both RDKit and
    // COSMolKit 0.5.0 read this back with sulfur stereochemistry intact.
    let block = "sulfoxide\n     RDKit          2D\n\n\
  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 9 9 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 2.583333 -1.299038 0.000000 0\n\
M  V30 2 S 1.833333 0.000000 0.000000 0 CHG=1\n\
M  V30 3 O 2.583333 1.299038 0.000000 0 CHG=-1\n\
M  V30 4 C 0.333333 -0.000000 0.000000 0\n\
M  V30 5 C -0.416667 -1.299038 0.000000 0\n\
M  V30 6 C -1.916667 -1.299038 0.000000 0\n\
M  V30 7 C -2.666667 0.000000 0.000000 0\n\
M  V30 8 C -1.916667 1.299038 0.000000 0\n\
M  V30 9 C -0.416667 1.299038 0.000000 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 2 1 CFG=3\n\
M  V30 2 1 2 3\n\
M  V30 3 1 2 4\n\
M  V30 4 2 4 5\n\
M  V30 5 1 5 6\n\
M  V30 6 2 6 7\n\
M  V30 7 1 7 8\n\
M  V30 8 2 8 9\n\
M  V30 9 1 9 4\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n";
    let report = read_mol_v3000_with_diagnostics(block).unwrap();
    assert!(report.stereo_diagnostics.is_empty());
    assert_eq!(
        rdkit_canonical_smiles(&report.mol).unwrap(),
        "C[S@@+]([O-])c1ccccc1"
    );
}

#[test]
fn v3000_incoming_wedge_does_not_conflict_with_centres_own_wedge() {
    // Two adjacent stereocentres: atom 2 owns the 2->4 wedge and atom 4
    // owns the 4->5 hash. The 2->4 wedge is incoming at atom 4 and must not
    // enter its local parity calculation. RDKit 2025.09.3 reads this CTAB
    // as F[C@H](Cl)[C@@H](Br)I.
    let block = "adjacent stereocentres\n     RDKit          2D\n\n\
  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 6 5 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 F -1.500000 1.299038 0.000000 0\n\
M  V30 2 C -0.750000 -0.000000 0.000000 0\n\
M  V30 3 Cl -1.500000 -1.299038 0.000000 0\n\
M  V30 4 C 0.750000 0.000000 0.000000 0\n\
M  V30 5 Br 1.500000 -1.299038 0.000000 0\n\
M  V30 6 I 1.500000 1.299038 0.000000 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 2 1\n\
M  V30 2 1 2 3\n\
M  V30 3 1 2 4 CFG=1\n\
M  V30 4 1 4 5 CFG=3\n\
M  V30 5 1 4 6\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n";
    let report = read_mol_v3000_with_diagnostics(block).unwrap();
    assert!(report.stereo_diagnostics.is_empty());
    assert_eq!(
        rdkit_canonical_smiles(&report.mol).unwrap(),
        "F[C@H](Cl)[C@@H](Br)I"
    );
}

#[test]
fn v3000_total_valence_excludes_dative_bond_on_donor() {
    // RDKit 2025.09.3 V3000 output for C[CH](C)->N. VAL=3 means two
    // covalent C-C bonds plus one H; the dative bond spends no donor valence.
    let block = "dative donor\n     RDKit          2D\n\n\
  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 4 3 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C -1.299038 -0.750000 0.000000 0\n\
M  V30 2 C 0.000000 0.000000 0.000000 0 RAD=2 VAL=3\n\
M  V30 3 C 1.299038 -0.750000 0.000000 0\n\
M  V30 4 N 0.000000 1.500000 0.000000 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 2 1 2 3\n\
M  V30 3 9 2 4\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n";
    let report = read_mol_v3000_with_diagnostics(block).unwrap();
    assert_eq!(
        report.mol.atom(chematic_core::AtomIdx(1)).hydrogen_count,
        Some(1)
    );
    assert_eq!(rdkit_canonical_smiles(&report.mol).unwrap(), "C[CH](C)->N");
}

#[test]
fn v3000_ring_constrained_nitrogen_retains_lone_pair_stereochemistry() {
    // RDKit 2025.09.3 V3000 blocks: an asymmetric bicyclic bridgehead N and
    // a substituted aziridine N. Both have three bonds, zero H and a wedge.
    let cases = [
        (
            "bridgehead",
            "bridgehead\n     RDKit          2D\n\n\
  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 7 8 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 1.552539 -0.750000 0.000000 0\n\
M  V30 2 C 1.552539 0.750000 0.000000 0\n\
M  V30 3 N 0.125954 1.213525 0.000000 0\n\
M  V30 4 C -1.300631 0.750000 0.000000 0\n\
M  V30 5 C -1.300631 -0.750000 0.000000 0\n\
M  V30 6 C 0.125954 -1.213525 0.000000 0\n\
M  V30 7 C -0.755724 0.000000 0.000000 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2 CFG=3\n\
M  V30 2 1 3 2 CFG=3\n\
M  V30 3 1 3 4\n\
M  V30 4 1 4 5\n\
M  V30 5 1 6 5 CFG=3\n\
M  V30 6 1 6 7\n\
M  V30 7 1 6 1\n\
M  V30 8 1 7 3\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n",
            "C1C[N@]2CC[C@H]1C2",
        ),
        (
            "aziridine",
            "aziridine\n     RDKit          2D\n\n\
  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 5 5 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C -2.049038 -0.709808 0.000000 0\n\
M  V30 2 N -0.750000 0.040192 0.000000 0\n\
M  V30 3 C 0.750000 0.040192 0.000000 0\n\
M  V30 4 F 2.049038 -0.709808 0.000000 0\n\
M  V30 5 C 0.000000 1.339230 0.000000 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 2 1 CFG=1\n\
M  V30 2 1 2 3\n\
M  V30 3 1 3 4\n\
M  V30 4 1 3 5\n\
M  V30 5 1 5 2\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n",
            "C[N@@]1CC1F",
        ),
    ];
    for (name, block, expected) in cases {
        let report = read_mol_v3000_with_diagnostics(block).unwrap();
        assert!(report.stereo_diagnostics.is_empty(), "{name}");
        assert_eq!(
            rdkit_canonical_smiles(&report.mol).unwrap(),
            expected,
            "{name}"
        );
    }
}
