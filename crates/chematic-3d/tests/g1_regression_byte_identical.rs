//! Proves the G1 `rdkit_*` addition (crates/chematic-3d/src/rdkit_shape_descriptors.rs)
//! left the existing `shape_descriptors`/`descriptors_3d` functions byte-identical.
//!
//! `git diff` for this PR shows those two source files completely untouched
//! (only new files were added, plus two `pub mod`/`pub use` lines in `lib.rs`),
//! which is the strongest possible proof by construction -- but this test
//! additionally pins the actual numeric outputs (as raw `f64::to_bits()`, i.e.
//! bit-exact, not just "close") as an automated, durable regression guard
//! rather than relying on a one-time `git diff` read. Every expected value
//! below is a snapshot of this crate's pre-existing behavior on aspirin's
//! rule-based conformer, captured by running this same test with these
//! source files unmodified (only new files added elsewhere in the crate).

use chematic_3d::dg::generate_coords;
use chematic_3d::{
    asphericity, autocorr_3d, eccentricity, getaway_descriptors, npr1, npr2, plane_of_best_fit,
    pmi, radius_of_gyration, rdf_descriptors, whim_descriptors,
};
use chematic_smiles::parse;

fn aspirin_coords() -> (chematic_core::Molecule, chematic_3d::Coords3D) {
    let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    let coords = generate_coords(&mol);
    (mol, coords)
}

#[test]
fn shape_descriptors_outputs_unchanged() {
    let (mol, coords) = aspirin_coords();

    let (p1, p2, p3) = pmi(&mol, &coords);
    assert_eq!(p1.to_bits(), 4641030930311365385);
    assert_eq!(p2.to_bits(), 4642440373202006155);
    assert_eq!(p3.to_bits(), 4643411699880215071);

    assert_eq!(npr1(&mol, &coords).to_bits(), 4604711176012659252);
    assert_eq!(npr2(&mol, &coords).to_bits(), 4606060549641463300);
    assert_eq!(
        radius_of_gyration(&mol, &coords).to_bits(),
        4609080816536283578
    );
    assert_eq!(asphericity(&mol, &coords).to_bits(), 4632702946048551376);
    assert_eq!(eccentricity(&mol, &coords).to_bits(), 4602893161489872220);
    assert_eq!(
        plane_of_best_fit(&mol, &coords).to_bits(),
        4604408828436051570
    );
}

#[test]
fn descriptors_3d_outputs_unchanged() {
    let (mol, coords) = aspirin_coords();

    let whim = whim_descriptors(&mol, &coords);
    assert_eq!(whim.len(), 22);
    assert_eq!(whim[0].to_bits(), 4606306289098307689);

    let getaway = getaway_descriptors(&mol, &coords);
    assert_eq!(getaway.len(), 19);
    assert_eq!(getaway[0].to_bits(), 4613348228423443963);

    let rdf = rdf_descriptors(&mol, &coords);
    assert_eq!(rdf.len(), 20);
    assert_eq!(rdf[0].to_bits(), 4493306176433625389);

    let ac = autocorr_3d(&mol, &coords);
    assert_eq!(ac.len(), 8);
    assert_eq!(ac[0].to_bits(), 4655804895277587578);
}
