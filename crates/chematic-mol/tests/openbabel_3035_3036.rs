//! Independent invariants from Open Babel issues #3035 and #3036.
//! The source SMILES define the configurations; Open Babel's output is not
//! used as the expected result.

use chematic_mol::{
    CoordinateDimension, GeometryRank, MolMetadata, read_mol_with_diagnostics, write_laid_out_mol,
};
use chematic_smiles::{canonical_smiles, parse};

#[test]
fn coordinate_free_stereo_is_drawn_before_mol_serialization() {
    // Open Babel #3035: a coordinate-free graph was labelled 2D while every
    // coordinate was zero, so a downstream InChI writer discarded stereo.
    let source = parse("C/C=C/[C@@H](N)O").unwrap();
    let (block, loss) = write_laid_out_mol(&source, &MolMetadata::default());
    assert!(loss.is_empty(), "{loss}: {block}");
    let parsed = read_mol_with_diagnostics(&block).unwrap();
    assert_eq!(canonical_smiles(&parsed.mol), canonical_smiles(&source));
    assert!(
        parsed.coords.iter().any(|&(x, y)| x != 0.0 || y != 0.0),
        "stereo-bearing MOL block must carry a real depiction: {block}"
    );
}

#[test]
fn openbabel_0d_mol_input_recovers_both_stereo_features() {
    // Open Babel 3.2.1 output for the #3035 example (with its all-zero atom
    // coordinates and a misleading 2D header). This checks the import path
    // independently of CheMatic's own writer.
    let block = include_str!("fixtures/openbabel-3035-0d-stereo.mol");
    let source = parse("C/C=C/[C@@H](N)O").unwrap();
    let parsed = read_mol_with_diagnostics(block).unwrap();
    assert_eq!(parsed.coordinate_dimension, CoordinateDimension::TwoD);
    assert_eq!(parsed.geometry_rank, GeometryRank::FlatZero);
    assert_eq!(canonical_smiles(&parsed.mol), canonical_smiles(&source));
}

#[test]
fn coupled_imine_stereo_round_trip_is_stable_across_repeated_writes() {
    // Open Babel #3036: pointer-ordered 0D MDL writer varies the exocyclic
    // C=N stereo between runs for this conjugated ring.
    let input = r"C/N=C/1\C(=C(\C1=N/C)O)O";
    let expected = canonical_smiles(&parse(input).unwrap());
    let opposite = canonical_smiles(&parse(r"C/N=C/1\C(=C(\C1=N\C)O)O").unwrap());
    assert_ne!(
        expected, opposite,
        "test must distinguish the wrong imine stereo"
    );
    let mut first = None;
    for _ in 0..40 {
        // Rebuild the graph each time so the writer cannot reuse cached
        // stereo decisions from a previous serialization.
        let source = parse(input).unwrap();
        let (block, loss) = write_laid_out_mol(&source, &MolMetadata::default());
        assert!(loss.is_empty(), "{loss}: {block}");
        if let Some(previous) = &first {
            assert_eq!(&block, previous);
        } else {
            first = Some(block.clone());
        }
        let parsed = read_mol_with_diagnostics(&block).unwrap();
        assert_eq!(canonical_smiles(&parsed.mol), expected, "{block}");
    }
}
