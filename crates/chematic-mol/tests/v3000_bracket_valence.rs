//! CTfile valence and spin must preserve bracket atoms when another reader
//! would otherwise fill their vacant valence with implicit hydrogens.

use chematic_mol::{MolMetadata, write_mol_v3000};
use chematic_smiles::parse;

#[test]
fn isolated_neutral_atoms_and_noncarbon_radicals_have_ctfile_valence() {
    for (source, atom_line, expected) in [
        ("[H]", "M  V30 1 H ", " VAL=-1"),
        ("[Si]", "M  V30 1 Si ", " VAL=-1"),
        ("[Na]", "M  V30 1 Na ", " VAL=-1"),
        ("[Na+]", "M  V30 1 Na ", " VAL=-1"),
        ("[C+]", "M  V30 1 C ", " VAL=-1"),
        ("[O+2]", "M  V30 1 O ", " VAL=-1"),
        ("[N-2]", "M  V30 1 N ", " VAL=-1"),
        ("[SH]", "M  V30 1 S ", " RAD=2 VAL=1"),
        ("[PH2]", "M  V30 1 P ", " RAD=2 VAL=2"),
        ("C[S]", "M  V30 2 S ", " RAD=2 VAL=1"),
        ("C[Si]", "M  V30 2 Si ", " RAD=2 VAL=1"),
        ("C[S](=O)=O", "M  V30 2 S ", " RAD=2 VAL=5"),
        ("Cl[I]Cl", "M  V30 2 I ", " RAD=2 VAL=2"),
        ("C[N+]", "M  V30 2 N ", " RAD=2 VAL=1"),
        ("C=[N+]", "M  V30 2 N ", " RAD=3 VAL=2"),
        ("C#[N+]", "M  V30 2 N ", " RAD=2 VAL=3"),
        ("[NH3+]", "M  V30 1 N ", " RAD=2 VAL=3"),
    ] {
        let mol = parse(source).unwrap();
        let block = write_mol_v3000(&mol, &MolMetadata::default(), &[]);
        let line = block
            .lines()
            .find(|line| line.starts_with(atom_line))
            .expect("target atom line");
        assert!(line.contains(expected), "{source}: {line}");
    }
}

#[test]
fn ordinary_valence_and_aromatic_centers_are_not_marked_as_radicals() {
    for source in ["C", "[CH4]", "[NH4+]", "c1cc[nH]c1"] {
        let mol = parse(source).unwrap();
        let block = write_mol_v3000(&mol, &MolMetadata::default(), &[]);
        assert!(!block.contains(" VAL=-1"), "{source}");
        assert!(!block.contains(" RAD="), "{source}");
    }
}
