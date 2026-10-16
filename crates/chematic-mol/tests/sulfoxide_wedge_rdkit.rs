//! RDKit 2025.09.3 2D MOL fixtures for a three-coordinate sulfur center.

use chematic_mol::{mol2000::read_mol_with_diagnostics, read_mol_v3000_with_diagnostics};
use chematic_smiles::{canonical_smiles, parse, rdkit_canonical_smiles, write};

const V2000: &str = include_str!("fixtures/rdkit-sulfoxide-2d-v2000.mol");
const V3000: &str = include_str!("fixtures/rdkit-sulfoxide-2d-v3000.mol");
const OPENBABEL_SELENOETHER: &str = include_str!("fixtures/openbabel-nonstereo-selenide-v2000.mol");

#[test]
fn sulfoxide_wedge_stereo_survives_both_mol_formats_and_smiles_writers() {
    for (name, source, inverted) in [
        (
            "V2000",
            V2000.to_owned(),
            V2000.replace("  2  1  1  1\n", "  2  1  1  6\n"),
        ),
        (
            "V3000",
            V3000.to_owned(),
            V3000.replace("M  V30 1 1 2 1 CFG=1", "M  V30 1 1 2 1 CFG=3"),
        ),
    ] {
        for (block, expected) in [
            (source, "C[S@](=O)c1ccccc1"),
            (inverted, "C[S@@](=O)c1ccccc1"),
        ] {
            let mol = if name == "V2000" {
                read_mol_with_diagnostics(&block).unwrap().mol
            } else {
                read_mol_v3000_with_diagnostics(&block).unwrap().mol
            };
            for smiles in [write(&mol), canonical_smiles(&mol)] {
                let reparsed = parse(&smiles).unwrap();
                assert_eq!(
                    rdkit_canonical_smiles(&reparsed).unwrap(),
                    expected,
                    "{name} {smiles}"
                );
            }
        }
    }
}

#[test]
fn openbabel_wedge_on_nonstereogenic_selenium_does_not_create_a_chiral_center() {
    // Open Babel emits both atom parity and a wedge at this ring selenium,
    // although its Se=C bond does not make the atom a stereogenic sulfoxide.
    let mol = read_mol_with_diagnostics(OPENBABEL_SELENOETHER)
        .unwrap()
        .mol;
    for smiles in [write(&mol), canonical_smiles(&mol)] {
        assert_eq!(
            rdkit_canonical_smiles(&parse(&smiles).unwrap()).unwrap(),
            "C/C(CN(C)C)=C(/n1ccnc1)[Se]1=C(OCC(C)C)C=CC=C1",
            "{smiles}"
        );
    }
}
