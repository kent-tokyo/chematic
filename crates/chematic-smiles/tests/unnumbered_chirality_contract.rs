use chematic_core::{AtomIdx, Chirality};
use chematic_smiles::parse;

#[test]
fn unnumbered_stereo_classes_are_not_treated_as_absent() {
    // RDKit 2025.09.3 recognizes these unnumbered class tags.
    for (smiles, expected) in [
        ("[C@TH](F)(Cl)(Br)I", Chirality::CounterClockwise),
        ("[P@TB](F)(Cl)(Br)(I)N", Chirality::TrigonalBipyramidal(0)),
        ("[Co@OH](F)(Cl)(Br)(I)(N)O", Chirality::Octahedral(0)),
    ] {
        let molecule = parse(smiles).unwrap();
        assert_eq!(molecule.atom(AtomIdx(0)).chirality, expected, "{smiles}");
    }
}
