//! Example: build a CsCl-type periodic structure (one cation + one anion
//! per cubic cell, anion at the body center -- the same site arrangement
//! as CsCl, AlNi, or beta-brass; illustrated here with Na/Cl for
//! familiarity, not a claim about real NaCl, which is rock-salt-type and
//! needs 8 sites in its conventional cubic cell), inspect the lattice,
//! enumerate periodic neighbors, and expand a supercell.
//!
//! Run: `cargo run -p chematic-crystal --example basic_structure`

use chematic_core::Element;
use chematic_crystal::{
    CrystalError, FractionalCoord, Lattice, PeriodicSite, PeriodicStructure, SiteSpecies,
};

fn main() -> Result<(), CrystalError> {
    // Illustrative cubic cell, a = 5.64 Angstrom.
    let lattice = Lattice::cubic(5.64)?;
    println!(
        "lattice: volume={:.3} A^3, condition_indicator={:.3}",
        lattice.volume(),
        lattice.condition_indicator()
    );

    let sites = vec![
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::NA)],
            FractionalCoord::new([0.0, 0.0, 0.0]),
            Some("Na1".to_string()),
        )?,
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::CL)],
            FractionalCoord::new([0.5, 0.5, 0.5]),
            Some("Cl1".to_string()),
        )?,
    ];
    let structure = PeriodicStructure::new(lattice, sites)?;
    println!("structure: {} sites", structure.site_count());

    // Nearest cross-species distance (body diagonal / 2 ~= 4.885 A here),
    // via periodic minimum image.
    let neighbors = structure.neighbors_within(5.0)?;
    println!("neighbors within 5.0 A: {}", neighbors.len());
    for nb in neighbors.iter().take(3) {
        println!(
            "  center={} neighbor={} image={:?} distance={:.4}",
            nb.center_index, nb.neighbor_index, nb.image, nb.distance
        );
    }

    // 2x2x2 supercell: 8x the sites, 8x the volume, same local geometry.
    let supercell = structure.make_supercell([2, 2, 2])?;
    println!(
        "supercell [2,2,2]: {} sites, volume={:.3} A^3",
        supercell.site_count(),
        supercell.lattice().volume()
    );

    Ok(())
}
